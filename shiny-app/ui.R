# ui.R

library(shiny)
library(ggiraph)

# ── Custom CSS ────────────────────────────────────────────────────────────────
app_css <- "
  /* ── Fonts & Base ── */
  * { box-sizing: border-box; }

  body, .shiny-app, h1, h2, h3, h4, p, label, .selectize-input, button,
  .btn, .nav-tabs, .tab-content, .well {
    font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif !important;
  }

  body {
    background-color: #f4f5f2;
    color: #2b2e2b;
    margin: 0;
    padding: 0;
  }

  /* ── Header banner ── */
  #app-header {
    background: #4a5e50;
    padding: 36px 48px 28px 48px;
    border-bottom: 3px solid #1e331e;
    box-shadow: 0 2px 12px rgba(0,0,0,0.18);
  }

  #app-header h1 {
    color: #e8f0e8;
    font-size: 1.95rem;
    font-weight: 600;
    letter-spacing: 0.01em;
    margin: 0 0 12px 0;
  }

  #app-header .intro-text {
    color: #c4d4c4;
    font-size: 1.05rem;
    line-height: 1.65;
    max-width: 820px;
    margin: 0;
  }

  #app-header .pub-link {
    display: inline-block;
    margin-top: 10px;
    color: #a8c8a8;
    font-size: 0.98rem;
    text-decoration: underline;
    letter-spacing: 0.01em;
  }

  #app-header .pub-link:hover { color: #d0e8d0; }

  /* ── Main content area ── */
  #app-body {
    padding: 28px 40px 48px 40px;
    max-width: 1400px;
    margin: 0 auto;
  }

  /* ── Tabs ── */
  .nav-tabs {
    border-bottom: 2px solid #b5c9b5 !important;
    margin-bottom: 0 !important;
  }

  .nav-tabs > li > a {
    color: #4a644a !important;
    font-size: 1.02rem !important;
    font-weight: 500 !important;
    letter-spacing: 0.03em !important;
    padding: 10px 20px !important;
    border: 1px solid transparent !important;
    border-radius: 4px 4px 0 0 !important;
    background: transparent !important;
    transition: background 0.15s, color 0.15s;
  }

  .nav-tabs > li > a:hover {
    background: #dde8dd !important;
    color: #2d4a2d !important;
  }

  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    background: #ffffff !important;
    color: #2d4a2d !important;
    border-color: #b5c9b5 #b5c9b5 #ffffff !important;
    font-weight: 600 !important;
  }

  .tab-content {
    background: #ffffff;
    border: 1px solid #b5c9b5;
    border-top: none;
    border-radius: 0 0 6px 6px;
    padding: 28px 28px 32px 28px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
  }

  /* ── Selector panel ── */
  .selector-panel {
    background: #f0f4f0;
    border: 1px solid #c8d8c8;
    border-radius: 6px;
    padding: 18px 22px;
    margin-bottom: 24px;
  }

  .selector-panel label {
    font-weight: 600;
    color: #2d4a2d;
    font-size: 1.0rem;
    letter-spacing: 0.03em;
    text-transform: uppercase;
    margin-bottom: 8px;
    display: block;
  }

  .selectize-control .selectize-input {
    border: 1.5px solid #9ab89a !important;
    border-radius: 4px !important;
    background: #ffffff !important;
    font-size: 1.0rem !important;
    color: #2b2e2b !important;
    box-shadow: none !important;
    padding: 8px 12px !important;
  }

  .selectize-control .selectize-input.focus {
    border-color: #4e7c4e !important;
    box-shadow: 0 0 0 2px rgba(78,124,78,0.15) !important;
  }

  .selectize-dropdown {
    border: 1.5px solid #9ab89a !important;
    border-radius: 0 0 4px 4px !important;
  }

  .selectize-dropdown .option:hover,
  .selectize-dropdown .option.active {
    background: #e8f0e8 !important;
    color: #2d4a2d !important;
  }

  /* ── Platform checkbox — vertical stack ── */
  #s4_show_platforms .shiny-options-group {
    display: block !important;
    margin-left: 0 !important;
  }

  /* ── Inline radio buttons ── */
  .shiny-options-group {
    display: flex !important;
    flex-wrap: wrap;
    gap: 8px 32px;
    align-items: center;
    margin-left: 12px !important;
  }

  .radio-inline {
    margin-left: 0 !important;
    padding-left: 4px !important;
  }

  /* ── Plot row labels ── */
  .plot-row-label {
    font-size: 1.2rem;
    font-weight: 600;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: #5a7a5a;
    margin: 18px 0 8px 0;
    padding-left: 2px;
    border-left: 3px solid #4e7c4e;
    padding-left: 10px;
  }

  /* ── Plot row subtitles ── */
  .plot-row-subtitle {
    font-size: 1.05rem;
    color: #7a8c7a;
    margin: 2px 0 10px 13px;
    line-height: 1.6;
  }

  /* ── Download tab ── */
  .download-row {
    display: flex;
    align-items: center;
    gap: 24px;
    padding: 14px 0;
    border-bottom: 1px solid #e8ede8;
  }

  .download-row:last-child { border-bottom: none; }

  .download-desc {
    flex: 1;
    font-size: 1.0rem;
    color: #444c44;
    line-height: 1.55;
  }

  .download-header {
    font-size: 0.92rem;
    font-weight: 600;
    letter-spacing: 0.05em;
    text-transform: uppercase;
    color: #6a8a6a;
    margin-bottom: 16px;
    padding-bottom: 8px;
    border-bottom: 2px solid #c8d8c8;
  }

  .btn-download {
    background: #3d6b3d !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 4px !important;
    font-size: 0.95rem !important;
    font-weight: 500 !important;
    letter-spacing: 0.02em !important;
    padding: 8px 18px !important;
    min-width: 110px !important;
    transition: background 0.15s !important;
    white-space: nowrap;
  }

  .btn-download:hover {
    background: #2d4a2d !important;
  }

  /* ── Section headings inside tabs ── */
  .tab-section-heading {
    font-size: 1.1rem;
    font-weight: 600;
    color: #2d4a2d;
    letter-spacing: 0.02em;
    margin: 0 0 18px 0;
    padding-bottom: 8px;
    border-bottom: 1.5px solid #d0e0d0;
  }

  /* ── Responsive plotly container ── */
  .plotly-grid-row {
    display: flex;
    gap: 16px;
    align-items: flex-start;
    margin-bottom: 4px;
  }

  .plotly-grid-row .plotly-cell {
    flex: 1;
    min-width: 0;
  }
"


# ── UI definition ─────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML(app_css))),

  # ── Header ──────────────────────────────────────────────────────────────────
  div(
    id = "app-header",
    h1("A multiomic lifespan signature in genetically diverse, diet-restricted mice"),
    p(
      class = "intro-text",
      "Data explorer application for metabolomics, lipidomics, and proteomics data from the DRiDO study"
    )
    # tags$a(class = "pub-link",
    #   href  = "https://doi.org/PLACEHOLDER",
    #   target = "_blank",
    #   "\U0001F4C4 Link to Publication — DOI: 10.XXXX/PLACEHOLDER"
    # )
  ),

  # ── Main body ────────────────────────────────────────────────────────────────
  div(
    id = "app-body",
    tabsetPanel(
      id = "main_tabs", type = "tabs",

      # ── Tab 1: Compound Visualizer ──────────────────────────────────────────
      tabPanel(
        "Compound Visualizer",
        div(
          class = "selector-panel",
          selectizeInput(
            inputId = "selected_compounds",
            label = "Select Compound(s)",
            choices = NULL,
            multiple = TRUE,
            width = "50%",
            options = list(
              placeholder      = "Type to search and select one or more compounds…",
              maxOptions       = 2000,
              closeAfterSelect = FALSE
            )
          ),
          radioButtons("point_toggle",
            label = NULL,
            choices = c("None", "Show Outliers", "Show All Data Points"),
            selected = "None", inline = TRUE
          ),
          p(
            style = "color: #7a8c7a; font-size: 1.02rem; margin: 4px 0 0 0;",
            "Select one or more compounds above to view plots."
          )
        ),
        uiOutput("compound_plots_ui")
      ),

      # ── Tab 2: Regression Results ───────────────────────────────────────────
      tabPanel(
        "Regression Results",
        p(class = "tab-section-heading", "Volcano Plots"),
        div(class = "plot-row-label", "Diet & Age Model — Table S4"),
        uiOutput("volcano_s4_subtitle"),
        div(
          class = "selector-panel", style = "margin-bottom: 16px;",
          selectizeInput(
            inputId = "s4_terms",
            label = "Select Terms",
            choices = s4_term_choices,
            selected = intersect(
              c("AgeYear2", "AgeYear3", "IF-1D", "IF-2D", "CR-20", "CR-40"),
              s4_term_choices
            ),
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          )
        ),
        uiOutput("volcano_s4_container"),
        div(class = "plot-row-label", style = "margin-top: 28px;", "Term Correlation"),
        div(
          class = "selector-panel", style = "margin-bottom: 16px;",
          fluidRow(
            column(
              6,
              selectizeInput(
                inputId  = "s4_corr_term_x",
                label    = "X-Axis Term",
                choices  = s4_term_choices,
                selected = "Diet20",
                multiple = FALSE
              )
            ),
            column(
              6,
              selectizeInput(
                inputId  = "s4_corr_term_y",
                label    = "Y-Axis Term",
                choices  = s4_term_choices,
                selected = "Diet40",
                multiple = FALSE
              )
            )
          )
        ),
        fluidRow(
          column(
            3,
            div(
              class = "selector-panel",
              selectizeInput(
                inputId = "s4_highlight_compounds",
                label = "Show Compounds",
                choices = NULL,
                multiple = TRUE,
                width = "100%",
                options = list(
                  placeholder      = "Type to search compounds…",
                  maxOptions       = 2000,
                  closeAfterSelect = FALSE
                )
              ),
              uiOutput("s4_highlight_table"),
              tags$label(
                style = "font-weight: 600; color: #2d4a2d; font-size: 1.0rem;
                                   letter-spacing: 0.03em; text-transform: uppercase;
                                   margin-bottom: 8px; display: block;",
                "Platforms"
              ),
              checkboxGroupInput(
                inputId  = "s4_show_platforms",
                label    = NULL,
                choices  = c("Lipidomics", "Metabolomics", "Proteomics"),
                selected = c("Lipidomics", "Metabolomics", "Proteomics")
              )
            )
          ),
          column(
            9,
            uiOutput("s4_corr_container")
          )
        )
      ),

      # ── Tab 3: PLL Trajectories ─────────────────────────────────────────────
      tabPanel(
        "PLL Trajectories",
        uiOutput("pll_tab_subtitle"),
        div(
          class = "selector-panel", style = "padding: 14px 22px; margin-bottom: 16px;",
          radioButtons("pll_norm_all",
            label = "Normalization",
            choices = s6_normalization_choices,
            selected = "Z-score",
            inline = TRUE
          ),
          checkboxInput("pll_show_outliers",
            label = "Show outliers",
            value = FALSE
          )
        ),
        fluidRow(
          # Aging volcano
          column(
            4,
            div(class = "plot-row-label", "Aging Model"),
            div(
              class = "selector-panel", style = "margin-bottom: 16px;",
              selectizeInput(
                inputId  = "s6_aging_terms",
                label    = "Select Terms",
                choices  = s6_aging_term_choices,
                selected = intersect("PLL_Aging", s6_aging_term_choices),
                multiple = TRUE,
                options  = list(plugins = list("remove_button"))
              )
            ),
            uiOutput("volcano_s6_aging_container")
          ),
          # Dying volcano
          column(
            4,
            div(class = "plot-row-label", "Dying Model"),
            div(
              class = "selector-panel", style = "margin-bottom: 16px;",
              selectizeInput(
                inputId  = "s6_dying_terms",
                label    = "Select Terms",
                choices  = s6_dying_term_choices,
                selected = intersect("PLL_Dying", s6_dying_term_choices),
                multiple = TRUE,
                options  = list(plugins = list("remove_button"))
              )
            ),
            uiOutput("volcano_s6_dying_container")
          ),
          # PLL scatter
          column(
            4,
            div(class = "plot-row-label", "Aging vs. Dying PLL Comparison"),
            girafeOutput("pll_scatter", height = "560px")
          )
        ),

        # Full-width compound selector
        div(
          class = "selector-panel", style = "margin-top: 8px;",
          div(class = "plot-row-label", style = "margin-top: 0;", "Select Compound(s)"),
          selectizeInput(
            inputId = "pll_compounds",
            label = NULL,
            choices = NULL,
            multiple = TRUE,
            width = "100%",
            options = list(
              placeholder      = "Type to search compounds…",
              maxOptions       = 2000,
              closeAfterSelect = FALSE
            )
          ),
          div(
            style = "display: flex; align-items: center; gap: 24px; margin-top: 4px;",
            checkboxInput(
              inputId = "pll_death_shift",
              label   = "Align aging/dying slopes",
              value   = TRUE
            ),
            actionButton("pll_clear", "Clear All",
              style = "background: #3d6b3d; color: white; border: none; border-radius: 4px; padding: 8px 18px; font-size: 0.95rem; cursor: pointer;"
            )
          )
        ),
        uiOutput("pll_trajectories_ui")
      ),

      # ── Tab 4: QTL Results ──────────────────────────────────────────────────
      tabPanel(
        "QTL Results",
        p(class = "tab-section-heading", "Manhattan Plots"),
        div(
          class = "selector-panel",
          selectizeInput(
            inputId = "qtl_compounds",
            label = "Filter by Compound(s)",
            choices = NULL,
            multiple = TRUE,
            width = "50%",
            options = list(
              placeholder      = "Type to search and select compounds — leave empty to show all",
              maxOptions       = 2000,
              closeAfterSelect = FALSE
            )
          )
        ),
        div(class = "plot-row-label", "Molecular QTL"),
        div(
          class = "plot-row-subtitle",
          "Lipid, metabolite, and protein QTL detected across all timepoints. Click any point to open the QTL region in NCBI Genome Data Viewer (GRCm39).",
          tags$br(),
          "More genetic exploration of molecular and physiological trait QTL available through Jackson Labs: ",
          tags$a(
            href = "https://churchilllab.jax.org/qtlviewer/DRiDO",
            target = "_blank",
            "churchilllab.jax.org/qtlviewer/DRiDO"
          )
        ),
        uiOutput("manhattan_qtl_container")
      ),

      # ── Tab 5: Data Download ─────────────────────────────────────────────────
      tabPanel(
        "Data Download",
        p(class = "tab-section-heading", "Click to download file"),
        div(
          class = "download-header",
          "File description / Download"
        ),
        lapply(download_meta, function(tbl) {
          div(
            class = "download-row",
            div(class = "download-desc", tbl$desc),
            downloadButton(
              outputId = paste0("dl_", tbl$id),
              label    = tbl$label,
              class    = "btn-download"
            )
          )
        })
      )
    ) # end tabsetPanel
  ) # end app-body
)
