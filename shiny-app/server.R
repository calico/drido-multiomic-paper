# server.R

library(shiny)
library(ggiraph)

source("server_fxns.R")

server <- function(input, output, session) {
  # ── Populate compound selector (server-side for large lists) ────────────────
  updateSelectizeInput(
    session,
    inputId = "selected_compounds",
    choices = compound_choices,
    server  = TRUE
  )

  # ── Compound Visualizer ─────────────────────────────────────────────────────

  # Filter all_data to selected compounds
  filtered_data <- reactive({
    req(input$selected_compounds)
    all_data %>% dplyr::filter(name_use %in% input$selected_compounds)
  })

  # Dynamically build UI rows based on number of selected compounds
  output$compound_plots_ui <- renderUI({
    compounds <- input$selected_compounds

    if (is.null(compounds) || length(compounds) == 0) {
      return(NULL)
    }

    n <- length(compounds)

    make_cell <- function(out_id, height = "320px") {
      div(class = "plotly-cell", girafeOutput(out_id, height = height))
    }

    # Build one row per plot type; columns = one per compound
    row_labels <- c("Diet and Age Effects", "Proportion of life lived trajectory", "Survival Correlation")
    row_subtitles <- c("", "LOESS spline shown with 80% confidence interval on all data", "Association between days of life remaining and compound abundance for each diet x timepoint group")
    row_heights <- c("320px", "320px", "640px")

    rows <- lapply(seq_along(row_labels), function(plot_idx) {
      sub <- row_subtitles[[plot_idx]]
      tagList(
        div(class = "plot-row-label", row_labels[[plot_idx]]),
        if (nchar(sub) > 0) div(class = "plot-row-subtitle", sub),
        div(
          class = "plotly-grid-row",
          lapply(seq_along(compounds), function(cmp_idx) {
            make_cell(paste0("cp_", plot_idx, "_", cmp_idx), height = row_heights[[plot_idx]])
          })
        )
      )
    })

    tagList(rows)
  })

  # Render each compound × plot combination dynamically
  observe({
    compounds <- input$selected_compounds
    req(compounds)

    plot_fns <- list(
      plot_compound_1,
      plot_compound_2,
      plot_compound_3
    )

    for (plot_idx in seq_along(plot_fns)) {
      for (cmp_idx in seq_along(compounds)) {
        local({
          pi <- plot_idx
          ci <- cmp_idx
          cmp <- compounds[[ci]]
          fn <- plot_fns[[pi]]
          oid <- paste0("cp_", pi, "_", ci)

          output[[oid]] <- renderGirafe({
            cmp_data <- filtered_data() %>% filter(name_use == cmp)
            show_outliers <- isTRUE(input$point_toggle == "Show Outliers")
            show_all_pts <- isTRUE(input$point_toggle == "Show All Data Points")
            if (pi == 1) {
              fn(cmp, cmp_data, show_outliers = show_outliers, show_all_points = show_all_pts)
            } else if (pi == 2) {
              fn(cmp, cmp_data, show_all_points = show_all_pts)
            } else {
              fn(cmp, cmp_data)
            }
          })
        })
      }
    }
  })

  # ── Regression Results ──────────────────────────────────────────────────────

  output$volcano_s4_subtitle <- renderUI({
    model_name <- unique(tbl_s4_diet_age$Model)[[1]]
    div(
      class = "plot-row-subtitle",
      paste0("Model estimates and p-values from ", model_name)
    )
  })

  output$volcano_s4_container <- renderUI({
    req(input$s4_terms)
    girafeOutput("volcano_s4", height = paste0(s4_n_platforms * 215, "px"))
  })

  output$volcano_s4 <- renderGirafe({
    req(input$s4_terms)
    plot_volcano_s4(tbl_s4_diet_age, selected_terms = input$s4_terms)
  })

  updateSelectizeInput(session, "s4_highlight_compounds",
    choices = s4_compound_choices, server = TRUE
  )

  output$s4_highlight_table <- renderUI({
    req(length(input$s4_highlight_compounds) > 0)

    recode_terms <- function(x) {
      dplyr::case_when(
        grepl("FastingFast20", x) ~ "Fasting-Diet20",
        grepl("FastingFast40", x) ~ "Fasting-Diet40",
        x == "(Intercept)" ~ "Intercept",
        x == "Diet1D" ~ "IF-1D",
        x == "Diet2D" ~ "IF-2D",
        x == "Diet20" ~ "CR-20",
        x == "Diet40" ~ "CR-40",
        TRUE ~ x
      )
    }

    tbl_filtered <- tbl_s4_diet_age %>%
      dplyr::filter(Compound %in% input$s4_highlight_compounds) %>%
      dplyr::mutate(Term = recode_terms(Term)) %>%
      dplyr::select(Compound, Platform, Term, Estimate, FDR) %>%
      dplyr::arrange(Compound, Platform, Term)

    if (nrow(tbl_filtered) == 0) {
      return(NULL)
    }

    compounds <- unique(tbl_filtered$Compound)

    cell_style <- "padding:3px 6px; border-bottom:1px solid #e8ede8;"
    header_style <- paste0(
      "text-align:left; border-bottom:2px solid #c8d8c8;",
      "padding:4px 6px; color:#4a644a; font-size:0.85rem;"
    )

    tagList(
      tags$hr(style = "border-top:1px solid #c8d8c8; margin:12px 0;"),
      lapply(compounds, function(cmp) {
        cmp_tbl <- tbl_filtered %>% dplyr::filter(Compound == cmp)
        tagList(
          tags$p(style = "font-weight:600; color:#2d4a2d; margin:10px 0 4px 0; font-size:0.95rem;", cmp),
          tags$table(
            style = "width:100%; border-collapse:collapse; font-size:0.88rem;",
            tags$thead(tags$tr(
              lapply(c("Platform", "Term", "Estimate", "FDR"), function(h) {
                tags$th(style = header_style, h)
              })
            )),
            tags$tbody(
              lapply(seq_len(nrow(cmp_tbl)), function(i) {
                row <- cmp_tbl[i, ]
                bg <- if (i %% 2 == 0) "#f0f4f0" else "#ffffff"
                tags$tr(
                  style = paste0("background:", bg, ";"),
                  tags$td(style = cell_style, row$Platform),
                  tags$td(style = cell_style, row$Term),
                  tags$td(style = cell_style, sprintf("%.3f", row$Estimate)),
                  tags$td(style = cell_style, formatC(row$FDR, format = "e", digits = 2))
                )
              })
            )
          )
        )
      })
    )
  })

  output$s4_corr_container <- renderUI({
    req(input$s4_corr_term_x, input$s4_corr_term_y)
    girafeOutput("s4_term_corr", height = "450px")
  })

  output$s4_term_corr <- renderGirafe({
    req(input$s4_corr_term_x, input$s4_corr_term_y)
    plot_s4_term_correlation(
      data                = tbl_s4_diet_age,
      term_x              = input$s4_corr_term_x,
      term_y              = input$s4_corr_term_y,
      show_platforms      = input$s4_show_platforms,
      highlight_compounds = if (length(input$s4_highlight_compounds) > 0) input$s4_highlight_compounds else character(0)
    )
  })

  # Combined top subtitle for Tab 3
  output$pll_tab_subtitle <- renderUI({
    aging_model <- unique(tbl_s6_pll_aging$Model[tbl_s6_pll_aging$Phase == "Aging"])[[1]]
    dying_model <- unique(tbl_s6_pll_aging$Model[tbl_s6_pll_aging$Phase == "Dying"])[[1]]
    div(
      class = "plot-row-subtitle", style = "margin-bottom: 10px;",
      tags$b("Aging Model:"),
      paste0(" Model estimates and p-values from ", aging_model, "."),
      tags$br(),
      tags$b("Dying Model:"),
      paste0(" Model estimates and p-values from ", dying_model, "."),
      tags$br(),
      tags$b("Aging vs. Dying PLL Comparison:"),
      " One point per compound per platform. Colored by platform if FDR < 0.01 on aging or dying slope; grey otherwise.",
      tags$br(),
      tags$b("Compound Selector:"),
      " Click points in the scatter or volcano plots to add compounds, or type to search below.",
      tags$br(),
      tags$b("Normalization and outlier settings"),
      " apply to all three plots. The outlier toggle removes points beyond Q1/Q3 ±6×IQR."
    )
  })

  # Aging
  output$volcano_s6_aging_container <- renderUI({
    req(input$s6_aging_terms)
    n_rows <- ceiling(length(input$s6_aging_terms) / 2)
    girafeOutput("volcano_s6_aging", height = paste0(n_rows * 250, "px"))
  })

  output$volcano_s6_aging <- renderGirafe({
    req(input$s6_aging_terms)
    d <- tbl_s6_pll_aging %>%
      filter(Phase == "Aging", Normalization == input$pll_norm_all)
    plot_volcano_s6(d,
      selected_terms = input$s6_aging_terms, phase = "Aging",
      show_outliers = isTRUE(input$pll_show_outliers)
    )
  })

  # Dying
  output$volcano_s6_dying_container <- renderUI({
    req(input$s6_dying_terms)
    n_rows <- ceiling(length(input$s6_dying_terms) / 2)
    girafeOutput("volcano_s6_dying", height = paste0(n_rows * 250, "px"))
  })

  output$volcano_s6_dying <- renderGirafe({
    req(input$s6_dying_terms)
    d <- tbl_s6_pll_aging %>%
      filter(Phase == "Dying", Normalization == input$pll_norm_all)
    plot_volcano_s6(d,
      selected_terms = input$s6_dying_terms, phase = "Dying",
      show_outliers = isTRUE(input$pll_show_outliers)
    )
  })

  # ── PLL Trajectories ─────────────────────────────────────────────────────────

  updateSelectizeInput(session, "pll_compounds",
    choices = pll_compound_choices, server = TRUE
  )

  pll_reset <- reactiveVal(0)

  output$pll_scatter <- renderGirafe({
    pll_reset()
    plot_pll_scatter(tbl_s6_pll_aging,
      normalization = input$pll_norm_all,
      show_outliers = isTRUE(input$pll_show_outliers)
    )
  })

  # Scatter / volcano clicks → additively update compound selector
  .add_pll_from_selection <- function(selected_safe) {
    new_sel <- pll_compound_choices[sapply(pll_compound_choices, .safe_id) %in% selected_safe]
    combined <- union(isolate(input$pll_compounds), new_sel)
    updateSelectizeInput(session, "pll_compounds", selected = combined)
  }

  observeEvent(input$pll_scatter_selected,
    {
      .add_pll_from_selection(input$pll_scatter_selected)
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  observeEvent(input$volcano_s6_aging_selected,
    {
      .add_pll_from_selection(input$volcano_s6_aging_selected)
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  observeEvent(input$volcano_s6_dying_selected,
    {
      .add_pll_from_selection(input$volcano_s6_dying_selected)
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  observeEvent(input$pll_clear, {
    updateSelectizeInput(session, "pll_compounds", selected = character(0))
    pll_reset(pll_reset() + 1)
  })

  # Trajectory plot UI (3 per row, fixed 1/3 width)
  output$pll_trajectories_ui <- renderUI({
    compounds <- input$pll_compounds
    if (is.null(compounds) || length(compounds) == 0) {
      return(NULL)
    }
    n <- length(compounds)
    row_idxs <- split(seq_len(n), ceiling(seq_len(n) / 3))
    tagList(
      div(class = "plot-row-label", "Trajectories"),
      lapply(row_idxs, function(idxs) {
        div(
          class = "plotly-grid-row",
          lapply(idxs, function(i) {
            div(
              class = "plotly-cell",
              style = "flex: 0 0 33.33%; max-width: 33.33%;",
              plotOutput(paste0("pll_traj_", i), height = "380px", width = "100%")
            )
          })
        )
      })
    )
  })

  # Individual trajectory renders
  observe({
    compounds <- input$pll_compounds
    req(compounds)
    for (i in seq_along(compounds)) {
      local({
        ci <- i
        cmp <- compounds[[ci]]
        oid <- paste0("pll_traj_", ci)
        output[[oid]] <- renderPlot(
          {
            platform <- tbl_s6_pll_aging$Platform[tbl_s6_pll_aging$Compound == cmp][1]
            plt_color <- unname(.pll_platform_colors[platform])
            if (is.na(plt_color)) plt_color <- "#4a5e50"
            plot_pll_trajectory(
              compound_name  = cmp,
              s6_data        = tbl_s6_pll_aging,
              raw_data       = all_data,
              normalization  = "Log2 Relative",
              death_shift    = isTRUE(input$pll_death_shift),
              platform_color = plt_color
            )
          },
          res = 96
        )
      })
    }
  })

  # ── QTL Results ─────────────────────────────────────────────────────────────

  updateSelectizeInput(
    session,
    inputId = "qtl_compounds",
    choices = qtl_compound_choices,
    server  = TRUE
  )

  output$manhattan_qtl_container <- renderUI({
    girafeOutput("manhattan_qtl", height = paste0(qtl_n_platforms * 300, "px"))
  })

  output$manhattan_qtl <- renderGirafe({
    selected <- input$qtl_compounds
    qtl_data <- if (length(selected) > 0) {
      tbl_s2_qtls %>% filter(name_use %in% selected)
    } else {
      tbl_s2_qtls
    }
    plot_manhattan_qtl(data = qtl_data)
  })

  # ── Data Download ────────────────────────────────────────────────────────────

  lapply(download_meta, function(meta) {
    output[[paste0("dl_", meta$id)]] <- downloadHandler(
      filename = function() sub("\\.Rds$", ".csv", meta$file),
      content = function(file) {
        readRDS(file.path("data", meta$file)) |> readr::write_csv(file)
      }
    )
  })
}
