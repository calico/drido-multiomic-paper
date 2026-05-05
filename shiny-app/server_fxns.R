# server_fxns.R — Define plot-generating functions here.
# This file is sourced by server.R.
# Each function receives pre-loaded data and returns a girafe object.

# ── Shared aesthetics ─────────────────────────────────────────────────────────
diet_colors <- c("seashell4", "skyblue", "royalblue4", "orange", "firebrick")
names(diet_colors) <- c("AL", "IF-1D", "IF-2D", "CR-20", "CR-40")

theme_docr <- ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position    = "bottom",
    axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title         = ggplot2::element_text(size = 14),
    plot.subtitle      = ggplot2::element_text(size = 12),
    axis.title         = ggplot2::element_text(size = 12),
    axis.text          = ggplot2::element_text(size = 11),
    legend.text        = ggplot2::element_text(size = 11),
    legend.title       = ggplot2::element_text(size = 12),
    strip.text         = ggplot2::element_text(size = 11),
    panel.grid.major   = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
    strip.background   = ggplot2::element_rect(fill = "grey90", colour = "grey30", linewidth = 1),
    panel.border       = ggplot2::element_rect(linewidth = 1, color = "grey30")
  )

# ── Compound Visualizer ───────────────────────────────────────────────────────
# These five functions are called once per selected compound.
# `compound_name` is a single string from tbl_s1_annotations$name_use.
# `data` is all_data pre-filtered to rows where name_use == compound_name.

plot_compound_1 <- function(compound_name, data, show_outliers = FALSE, show_all_points = FALSE) {
  # Pre-compute box stats so the tooltip can show IQR summary per diet × Age
  box_tooltips <- data %>%
    dplyr::group_by(diet, Age) %>%
    dplyr::summarise(
      .med = median(trait_value, na.rm = TRUE),
      .q1 = quantile(trait_value, 0.25, na.rm = TRUE),
      .q3 = quantile(trait_value, 0.75, na.rm = TRUE),
      .iqr = IQR(trait_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(tooltip = sprintf(
      "Diet: %s<br>Median: %.2f<br>Q1: %.2f | Q3: %.2f<br>IQR: %.2f",
      diet, .med, .q1, .q3, .iqr
    ))

  data <- dplyr::left_join(data, dplyr::select(box_tooltips, diet, Age, tooltip),
    by = c("diet", "Age")
  )

  # Outlier points computed separately so they get mouse_id tooltips
  if (show_outliers) {
    outlier_pts <- data %>%
      dplyr::group_by(diet, Age) %>%
      dplyr::mutate(
        .q1  = quantile(trait_value, 0.25, na.rm = TRUE),
        .q3  = quantile(trait_value, 0.75, na.rm = TRUE),
        .iqr = .q3 - .q1
      ) %>%
      dplyr::filter(trait_value < .q1 - 1.5 * .iqr | trait_value > .q3 + 1.5 * .iqr) %>%
      dplyr::ungroup()
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = diet, color = diet, fill = diet, y = trait_value)) +
    ggiraph::geom_boxplot_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = diet),
      alpha = 0.5,
      outliers = FALSE,
      width = 0.7,
      linewidth = 1,
      position = ggplot2::position_dodge2(preserve = "single")
    )

  if (show_all_points) {
    p <- p + ggiraph::geom_jitter_interactive(
      ggplot2::aes(tooltip = mouse_id, data_id = mouse_id),
      width = 0.2, size = 1.2, alpha = 0.5
    )
  }

  if (show_outliers) {
    p <- p + ggiraph::geom_point_interactive(
      data = outlier_pts,
      ggplot2::aes(tooltip = mouse_id, data_id = mouse_id),
      size = 2
    )
  }

  p <- p +
    ggplot2::scale_color_manual(name = "Diet", values = diet_colors) +
    ggplot2::scale_fill_manual(name = "Diet", values = diet_colors) +
    theme_docr +
    ggplot2::scale_y_continuous(
      limits = \(x) c(min(x[1], -1), max(x[2], 1)),
      breaks = \(x) pretty(x, n = 3)
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(Age), scales = "free_x", space = "free_x") +
    ggplot2::labs(title = compound_name, y = expression("Log"[2] * " Abundance"), x = "") +
    ggplot2::theme(
      panel.spacing.x = ggplot2::unit(0.1, "cm"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  ggiraph::girafe(
    ggobj = p,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(1.05); transform-box:fill-box; transform-origin:center;")
    )
  )
}

plot_compound_2 <- function(compound_name, data, show_all_points = FALSE) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = PLL, y = trait_value, color = diet, fill = diet))

  if (show_all_points) {
    p <- p + ggiraph::geom_point_interactive(
      ggplot2::aes(tooltip = mouse_id, data_id = mouse_id),
      size = 1.2, alpha = 0.3
    )
  }

  p <- p +
    ggiraph::geom_smooth_interactive(
      ggplot2::aes(tooltip = diet, data_id = diet),
      method = "loess", level = 0.8, linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(name = "Diet", values = diet_colors) +
    ggplot2::scale_fill_manual(name = "Diet", values = diet_colors) +
    theme_docr +
    ggplot2::scale_y_continuous(
      limits = \(x) c(min(x[1], -1), max(x[2], 1)),
      breaks = \(x) pretty(x, n = 3)
    ) +
    ggplot2::labs(title = compound_name, y = expression("Log"[2] * " Abundance"), x = "Proportion of Life Lived")

  ggiraph::girafe(
    ggobj = p,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(1.05); transform-box:fill-box; transform-origin:center;")
    )
  )
}

plot_compound_3 <- function(compound_name, data) {
  # Per-facet correlation stats for annotation
  facet_stats <- data %>%
    dplyr::group_by(diet, Age) %>%
    dplyr::summarise(
      r = cor(days_remaining, trait_value, use = "complete.obs"),
      pval = tryCatch(cor.test(days_remaining, trait_value)$p.value, error = \(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r2    = r^2,
      p_lab = dplyr::if_else(pval < 0.001, "p<0.001", sprintf("p=%.3f", pval)),
      label = sprintf("R=%.2f, R\u00b2=%.2f, %s", r, r2, p_lab)
    )

  p <- ggplot2::ggplot(data, ggplot2::aes(x = days_remaining, y = trait_value, color = diet)) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(tooltip = mouse_id, data_id = mouse_id),
      size = 1.2, alpha = 0.4, show.legend = FALSE
    ) +
    ggiraph::geom_smooth_interactive(
      ggplot2::aes(tooltip = diet, data_id = diet, fill = diet),
      method = "lm", formula = "y ~ x", linewidth = 0.5, alpha = 0.8
    ) +
    ggplot2::geom_text(
      data = facet_stats,
      ggplot2::aes(x = -Inf, y = -Inf, label = label),
      hjust = -0.05, vjust = -0.8,
      size = 3.3, color = "grey30", inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(diet), cols = ggplot2::vars(Age), scales = "free") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.2, 0.05))) +
    ggplot2::scale_color_manual(name = "Diet", values = diet_colors) +
    ggplot2::scale_fill_manual(name = "Diet", values = diet_colors) +
    theme_docr +
    ggplot2::labs(
      title = compound_name,
      y     = expression("Log"[2] * " Abundance"),
      x     = "Days of Life Remaining"
    )

  ggiraph::girafe(
    ggobj = p,
    height_svg = 10,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(1.05); transform-box:fill-box; transform-origin:center;")
    )
  )
}

# ── Regression Results — Volcano Plots ───────────────────────────────────────
# Columns: Compound, Model, Term, Estimate, Est_SE, Pvalue, FDR

# Term color palette matching diet + age conventions
.diet_term_names <- c("IF-1D", "IF-2D", "CR-20", "CR-40")
.age_term_names <- c("AgeYear1", "AgeYear2", "AgeYear3")
.age_colors <- c(AgeYear1 = "#DCA8FB", AgeYear2 = "#C05CFA", AgeYear3 = "#5E2D7A")

.term_colors <- c(
  DietAL = unname(diet_colors["AL"]),
  setNames(unname(diet_colors[c("IF-1D", "IF-2D", "CR-20", "CR-40")]), .diet_term_names),
  .age_colors,
  PLL_Aging = "#1E88E5",
  PLL_Dying = "#E53935",
  `Fasting-Diet20` = "#555555",
  `Fasting-Diet40` = "#333333",
  BW_Loess = "#666666",
  Intercept = "#888888",
  No = "white"
)

.safe <- function(x) gsub("'", "&#39;", x)
.safe_id <- function(x) gsub("'", "_", x)

plot_volcano_s4 <- function(data, selected_terms = c("AgeYear2", "AgeYear3", "IF-1D", "IF-2D", "CR-20", "CR-40")) {
  plot_data <- data %>%
    dplyr::mutate(
      Term = dplyr::case_when(
        grepl("FastingFast20", Term) ~ "Fasting-Diet20",
        grepl("FastingFast40", Term) ~ "Fasting-Diet40",
        grepl("Intercept", Term) ~ "Intercept",
        Term == "Diet1D" ~ "IF-1D",
        Term == "Diet2D" ~ "IF-2D",
        Term == "Diet20" ~ "CR-20",
        Term == "Diet40" ~ "CR-40",
        TRUE ~ Term
      )
    ) %>%
    dplyr::filter(Term %in% selected_terms) %>%
    dplyr::mutate(
      Term = factor(Term, levels = selected_terms),
      neg_log10_p = -log10(Pvalue),
      # Significant points get their term color; non-significant get grey
      pt_color = dplyr::if_else(!is.na(FDR) & FDR < 0.05, as.character(Term), "n.s."),
      tooltip = sprintf(
        "%s<br>Estimate: %.3f<br>P-value: %.2e<br>FDR: %.2e",
        .safe(Compound), Estimate, Pvalue, FDR
      )
    )

  n_terms <- length(selected_terms)
  n_plat <- dplyr::n_distinct(plot_data$Platform)
  pt_size <- max(1.5, n_terms * 0.35)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Estimate, y = neg_log10_p)) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(color = pt_color, tooltip = tooltip, data_id = .safe_id(Compound)),
      size = pt_size, alpha = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::scale_color_manual(
      values   = c(.term_colors, "n.s." = "#d9d9d9"),
      na.value = "#777777",
      guide    = "none"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Platform), cols = ggplot2::vars(Term),
      scales = "free"
    ) +
    theme_docr +
    ggplot2::labs(x = "Estimate", y = expression("-log"[10] * "(p-value)"))

  ggiraph::girafe(
    ggobj = p,
    width_svg = n_terms * 2.5,
    height_svg = n_plat * 2.5,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(2.3); transform-box:fill-box; transform-origin:center; fill:#cc0000 !important; stroke:#cc0000 !important;")
    )
  )
}

plot_volcano_s6 <- function(data, selected_terms, phase, show_outliers = FALSE) {
  plot_data <- data %>%
    dplyr::mutate(Term = dplyr::case_when(
      Term == "(Intercept)" ~ "Intercept",
      Term == "Diet1D" ~ "IF-1D",
      Term == "Diet2D" ~ "IF-2D",
      Term == "Diet20" ~ "CR-20",
      Term == "Diet40" ~ "CR-40",
      TRUE ~ Term
    )) %>%
    dplyr::filter(Term %in% selected_terms) %>%
    dplyr::mutate(
      Term = factor(Term, levels = selected_terms),
      neg_log10_p = -log10(Pvalue),
      pt_color = dplyr::if_else(!is.na(FDR) & FDR < 0.05, as.character(Term), "n.s."),
      tooltip = sprintf(
        "%s (%s)<br>Estimate: %.3f<br>P-value: %.2e<br>FDR: %.2e",
        .safe(Compound), Platform, Estimate, Pvalue, FDR
      )
    )

  if (!show_outliers) {
    plot_data <- plot_data %>%
      dplyr::filter(
        Estimate >= quantile(Estimate, 0.25, na.rm = TRUE) - 6 * IQR(Estimate, na.rm = TRUE),
        Estimate <= quantile(Estimate, 0.75, na.rm = TRUE) + 6 * IQR(Estimate, na.rm = TRUE),
        neg_log10_p <= quantile(neg_log10_p, 0.75, na.rm = TRUE) + 6 * IQR(neg_log10_p, na.rm = TRUE)
      )
  }

  n_terms <- length(selected_terms)
  pt_size <- max(1.5, n_terms * 0.35)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Estimate, y = neg_log10_p)) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(color = pt_color, tooltip = tooltip, data_id = .safe_id(Compound)),
      size = pt_size, alpha = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::scale_color_manual(
      values   = c(.term_colors, "n.s." = "#d9d9d9"),
      na.value = "#777777",
      guide    = "none"
    ) +
    ggplot2::facet_wrap(~Term, ncol = 2, scales = "free") +
    theme_docr +
    ggplot2::labs(x = "Estimate", y = expression("-log"[10] * "(p-value)"))

  n_rows <- ceiling(n_terms / 2)
  ggiraph::girafe(
    ggobj = p,
    width_svg = min(n_terms, 2) * 2.5,
    height_svg = n_rows * 2.5,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(2.3); transform-box:fill-box; transform-origin:center; fill:#cc0000 !important; stroke:#cc0000 !important;"),
      ggiraph::opts_selection(type = "multiple", css = "fill:#ff8c00 !important; stroke:#ff8c00 !important;")
    )
  )
}


# ── PLL Trajectories ─────────────────────────────────────────────────────────

.pll_platform_colors <- c(
  Lipidomics    = "#39BC95",
  Metabolomics  = "#F77D20",
  Physiological = "#1E88E5",
  Proteomics    = "#938ED1",
  "n.s."        = "#d9d9d9"
)

plot_pll_scatter <- function(data, normalization = "Z-score",
                             show_outliers = FALSE, p_thresh = 0.01) {
  scatter_data <- data %>%
    dplyr::filter(
      Term %in% c("PLL_Aging", "PLL_Dying"),
      Normalization == normalization
    ) %>%
    dplyr::select(Platform, Compound, Term, Estimate, FDR) %>%
    tidyr::pivot_wider(
      names_from  = Term,
      values_from = c(Estimate, FDR),
      values_fn   = dplyr::first
    ) %>%
    dplyr::mutate(
      sig = (!is.na(FDR_PLL_Aging) & FDR_PLL_Aging < p_thresh) |
        (!is.na(FDR_PLL_Dying) & FDR_PLL_Dying < p_thresh),
      color_var = dplyr::if_else(sig, Platform, "n.s."),
      tooltip = sprintf(
        "%s (%s)<br>Aging slope: %.3f (FDR=%.2e)<br>Dying slope: %.3f (FDR=%.2e)",
        .safe(Compound), Platform,
        Estimate_PLL_Aging, FDR_PLL_Aging,
        Estimate_PLL_Dying, FDR_PLL_Dying
      )
    )

  if (!show_outliers) {
    .iqr_lo <- function(x) quantile(x, 0.25, na.rm = TRUE) - 6 * IQR(x, na.rm = TRUE)
    .iqr_hi <- function(x) quantile(x, 0.75, na.rm = TRUE) + 6 * IQR(x, na.rm = TRUE)
    scatter_data <- scatter_data %>%
      dplyr::filter(
        Estimate_PLL_Dying >= .iqr_lo(Estimate_PLL_Dying),
        Estimate_PLL_Dying <= .iqr_hi(Estimate_PLL_Dying),
        Estimate_PLL_Aging >= .iqr_lo(Estimate_PLL_Aging),
        Estimate_PLL_Aging <= .iqr_hi(Estimate_PLL_Aging)
      )
  }

  p <- ggplot2::ggplot(
    scatter_data,
    ggplot2::aes(x = Estimate_PLL_Dying, y = Estimate_PLL_Aging)
  ) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(color = color_var, tooltip = tooltip, data_id = .safe_id(Compound)),
      size = 1.5, alpha = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::scale_color_manual(
      name = "Platform",
      values = .pll_platform_colors,
      breaks = c("Lipidomics", "Metabolomics", "Physiological", "Proteomics"),
      na.value = "#d9d9d9"
    ) +
    theme_docr +
    ggplot2::theme(
      legend.position  = "top",
      legend.direction = "vertical",
      legend.title     = ggplot2::element_text(size = 11)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::labs(
      title = "PLL Aging vs. Dying Coefficients",
      x     = "Dying Slope (Post-0.85 PLL)",
      y     = "Aging Slope (Pre-0.85 PLL)"
    )

  ggiraph::girafe(
    ggobj = p,
    width_svg = 5,
    height_svg = 5.5,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(1.5); transform-box:fill-box; transform-origin:center;"),
      ggiraph::opts_selection(type = "multiple", css = "fill:#ff8c00; stroke:#ff8c00;")
    )
  )
}

plot_pll_trajectory <- function(compound_name, s6_data, raw_data,
                                normalization, death_shift = FALSE,
                                platform_color = "#4a5e50") {
  aging_rows <- s6_data %>%
    dplyr::filter(
      Phase == "Aging", Compound == compound_name,
      Normalization == normalization
    )
  dying_rows <- s6_data %>%
    dplyr::filter(
      Phase == "Dying", Compound == compound_name,
      Normalization == normalization
    )

  aging_coefs <- aging_rows %>%
    dplyr::select(Term, Estimate) %>%
    tibble::deframe()
  dying_coefs <- dying_rows %>%
    dplyr::select(Term, Estimate) %>%
    tibble::deframe()

  aging_slope_row <- aging_rows %>%
    dplyr::filter(Term == "PLL_Aging") %>%
    dplyr::slice(1)
  dying_slope_row <- dying_rows %>%
    dplyr::filter(Term == "PLL_Dying") %>%
    dplyr::slice(1)

  plot_df <- raw_data %>%
    dplyr::filter(name_use == compound_name, Age != "Year1") %>%
    tidyr::drop_na(trait_value, PLL)

  if (length(aging_coefs) == 0 && length(dying_coefs) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text",
        x = 0.5, y = 0.5,
        label = "No model data available", color = "#7a8c7a", size = 4
      ) +
      ggplot2::labs(title = compound_name) +
      ggplot2::theme_void())
  }

  has_raw <- nrow(plot_df) > 0
  pll_min <- if (has_raw) min(plot_df$PLL, na.rm = TRUE) else 0
  pll_max <- if (has_raw) max(plot_df$PLL, na.rm = TRUE) else 1
  bw_aging <- if ("BW_Loess" %in% names(aging_coefs)) aging_coefs[["BW_Loess"]] * 32 else 0
  bw_dying <- if ("BW_Loess" %in% names(dying_coefs)) dying_coefs[["BW_Loess"]] * 32 else 0

  line_aging <- if (pll_min < 0.845 && "PLL_Aging" %in% names(aging_coefs)) {
    tibble::tibble(x = seq(pll_min, min(0.845, pll_max), length.out = 100)) %>%
      dplyr::mutate(y = aging_coefs[["(Intercept)"]] + bw_aging +
        aging_coefs[["PLL_Aging"]] * x)
  } else {
    NULL
  }

  line_dying <- if (pll_max > 0.85 && "PLL_Dying" %in% names(dying_coefs)) {
    tibble::tibble(x = seq(max(0.85, pll_min), pll_max, length.out = 100)) %>%
      dplyr::mutate(y = dying_coefs[["(Intercept)"]] + bw_dying +
        dying_coefs[["PLL_Dying"]] * x)
  } else {
    NULL
  }

  if (death_shift && !is.null(line_aging) && !is.null(line_dying)) {
    shift <- dplyr::last(line_aging$y) - dplyr::first(line_dying$y)
    line_dying <- line_dying %>% dplyr::mutate(y = y + shift)
  }

  if (has_raw) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PLL, y = trait_value)) +
      ggplot2::geom_point(color = platform_color, size = 0.8, alpha = 0.3) +
      ggplot2::geom_smooth(
        method = "gam", formula = y ~ s(x), se = FALSE,
        color = platform_color, linewidth = 0.8
      )
  } else {
    all_y <- c(line_aging$y, line_dying$y)
    y_pad <- diff(range(all_y, na.rm = TRUE)) * 0.1
    p <- ggplot2::ggplot() +
      ggplot2::scale_x_continuous(limits = c(pll_min, pll_max)) +
      ggplot2::scale_y_continuous(limits = range(all_y, na.rm = TRUE) + c(-y_pad, y_pad))
  }

  p <- p +
    ggplot2::geom_vline(
      xintercept = 0.85, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    theme_docr +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(
      title = compound_name, x = "Proportion of Life Lived",
      y = expression("Log"[2] * " Abundance")
    )

  if (!is.null(line_aging)) {
    p <- p + ggplot2::geom_line(
      data = line_aging, ggplot2::aes(x = x, y = y),
      color = "black", linewidth = 0.7, alpha = 0.7,
      inherit.aes = FALSE
    )
  }
  if (!is.null(line_dying)) {
    p <- p + ggplot2::geom_line(
      data = line_dying, ggplot2::aes(x = x, y = y),
      color = "black", linewidth = 0.7, alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  ann_lines <- character(0)
  if (nrow(aging_slope_row) > 0) {
    ann_lines <- c(ann_lines, sprintf(
      "Aging slope: %.3f (FDR=%.2e)",
      aging_slope_row$Estimate, aging_slope_row$FDR
    ))
  }
  if (nrow(dying_slope_row) > 0) {
    ann_lines <- c(ann_lines, sprintf(
      "Dying slope:  %.3f (FDR=%.2e)",
      dying_slope_row$Estimate, dying_slope_row$FDR
    ))
  }

  if (length(ann_lines) > 0) {
    aging_est <- if (nrow(aging_slope_row) > 0) aging_slope_row$Estimate else 0
    ann_y <- if (aging_est >= 0) Inf else -Inf
    ann_vjust <- if (aging_est >= 0) 1.4 else -0.4
    p <- p + ggplot2::annotate(
      "text",
      x = -Inf, y = ann_y,
      label = paste(ann_lines, collapse = "\n"),
      hjust = -0.05, vjust = ann_vjust,
      size = 3, color = "grey30", family = "Helvetica"
    )
  }

  p
}

# ── Regression Results — Term Correlation Scatter ────────────────────────────

plot_s4_term_correlation <- function(data, term_x, term_y,
                                     show_platforms = c("Lipidomics", "Metabolomics", "Physiological", "Proteomics"),
                                     highlight_compounds = character(0),
                                     p_thresh = 0.05) {
  scatter_data <- data %>%
    dplyr::mutate(
      Term = dplyr::case_when(
        grepl("FastingFast20", Term) ~ "Fasting-Diet20",
        grepl("FastingFast40", Term) ~ "Fasting-Diet40",
        Term == "(Intercept)" ~ "Intercept",
        Term == "Diet1D" ~ "IF-1D",
        Term == "Diet2D" ~ "IF-2D",
        Term == "Diet20" ~ "CR-20",
        Term == "Diet40" ~ "CR-40",
        TRUE ~ Term
      )
    ) %>%
    dplyr::filter(Term %in% c(term_x, term_y)) %>%
    dplyr::select(Platform, Compound, Term, Estimate, FDR) %>%
    tidyr::pivot_wider(
      names_from  = Term,
      values_from = c(Estimate, FDR),
      values_fn   = dplyr::first
    )

  col_est_x <- paste0("Estimate_", term_x)
  col_est_y <- paste0("Estimate_", term_y)
  col_fdr_x <- paste0("FDR_", term_x)
  col_fdr_y <- paste0("FDR_", term_y)

  scatter_data <- scatter_data %>%
    dplyr::mutate(
      sig = (!is.na(.data[[col_fdr_x]]) & .data[[col_fdr_x]] < p_thresh) |
        (!is.na(.data[[col_fdr_y]]) & .data[[col_fdr_y]] < p_thresh),
      color_var = dplyr::if_else(sig, Platform, "n.s."),
      tooltip = sprintf(
        "%s (%s)<br>%s: %.3f (FDR=%.2e)<br>%s: %.3f (FDR=%.2e)",
        .safe(Compound), Platform,
        term_x, .data[[col_est_x]], .data[[col_fdr_x]],
        term_y, .data[[col_est_y]], .data[[col_fdr_y]]
      )
    )

  # Outlier filter (always on)
  .iqr_lo <- function(x) quantile(x, 0.25, na.rm = TRUE) - 5 * IQR(x, na.rm = TRUE)
  .iqr_hi <- function(x) quantile(x, 0.75, na.rm = TRUE) + 5 * IQR(x, na.rm = TRUE)
  scatter_data <- scatter_data %>%
    dplyr::filter(
      .data[[col_est_x]] >= .iqr_lo(.data[[col_est_x]]),
      .data[[col_est_x]] <= .iqr_hi(.data[[col_est_x]]),
      .data[[col_est_y]] >= .iqr_lo(.data[[col_est_y]]),
      .data[[col_est_y]] <= .iqr_hi(.data[[col_est_y]])
    )

  # Platform visibility filter
  if (length(show_platforms) > 0) {
    scatter_data <- scatter_data %>%
      dplyr::filter(Platform %in% show_platforms)
  }

  p <- ggplot2::ggplot(
    scatter_data,
    ggplot2::aes(x = .data[[col_est_x]], y = .data[[col_est_y]])
  ) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(color = color_var, tooltip = tooltip, data_id = .safe_id(Compound)),
      size = 1.5, alpha = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.4
    ) +
    ggplot2::scale_color_manual(
      name     = "Platform",
      values   = .pll_platform_colors,
      breaks   = c("Lipidomics", "Metabolomics", "Physiological", "Proteomics"),
      na.value = "#d9d9d9"
    ) +
    theme_docr +
    ggplot2::labs(
      x = paste0(term_x, " Estimate"),
      y = paste0(term_y, " Estimate")
    )

  # Highlight layer — drawn on top in bright red at 2× size
  if (length(highlight_compounds) > 0) {
    hi_data <- scatter_data %>% dplyr::filter(Compound %in% highlight_compounds)
    if (nrow(hi_data) > 0) {
      p <- p + ggiraph::geom_point_interactive(
        data = hi_data,
        ggplot2::aes(tooltip = tooltip, data_id = .safe_id(Compound)),
        color = "#E53935",
        size = 3,
        alpha = 1
      )
    }
  }

  ggiraph::girafe(
    ggobj = p,
    width_svg = 5,
    height_svg = 4.5,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(1.5); transform-box:fill-box; transform-origin:center;"),
      ggiraph::opts_selection(type = "multiple", css = "fill:#ff8c00; stroke:#ff8c00;")
    )
  )
}

# ── QTL Results — Manhattan Plots ────────────────────────────────────────────
# Data available: tbl_s2_qtls (loaded in global.R)
# Columns used: trait, chr, pos, lod, modality

plot_manhattan_qtl <- function(data = tbl_s2_qtls) {
  chr_levels <- c(as.character(1:19), "X")
  odd_chrs <- as.character(c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19))
  even_chrs <- c(as.character(c(2, 4, 6, 8, 10, 12, 14, 16, 18)), "X")

  qtl_sum <- data %>%
    dplyr::mutate(chr = factor(as.character(chr), levels = chr_levels)) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_pos = max(pos), .groups = "drop") %>%
    dplyr::arrange(chr) %>%
    dplyr::mutate(pos_add = dplyr::lag(cumsum(max_pos), default = 0)) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::select(chr, pos_add)

  plot_data <- data %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::inner_join(qtl_sum, by = "chr") %>%
    dplyr::mutate(
      chr = factor(chr, levels = chr_levels),
      pos_cum = pos + pos_add,
      Platform = stringr::str_to_title(modality),
      age_years = stringr::str_to_title(substr(trait, nchar(trait) - 4, nchar(trait))),
      nice_label = paste0(name_use, " - ", age_years),
      color_var = dplyr::case_when(
        Platform == "Lipidomics" & as.character(chr) %in% odd_chrs ~ "A",
        Platform == "Lipidomics" & as.character(chr) %in% even_chrs ~ "B",
        Platform == "Metabolomics" & as.character(chr) %in% odd_chrs ~ "C",
        Platform == "Metabolomics" & as.character(chr) %in% even_chrs ~ "D",
        Platform == "Proteomics" & as.character(chr) %in% odd_chrs ~ "E",
        Platform == "Proteomics" & as.character(chr) %in% even_chrs ~ "F",
        TRUE ~ NA_character_
      ),
      tooltip = sprintf(
        "Trait: %s<br>LOD: %.2f<br>Chr %s: %.1f Mb",
        .safe(nice_label), lod, as.character(chr), pos
      ),
      onclick = sprintf(
        'window.open("https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001635.27&chr=%s&from=%d&to=%d", "_blank")',
        as.character(chr),
        as.integer(dplyr::coalesce(ci_lo, pos - 1) * 1e6),
        as.integer(dplyr::coalesce(ci_hi, pos + 1) * 1e6)
      )
    )

  axis_set <- plot_data %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(center = mean(pos_cum), .groups = "drop") %>%
    dplyr::arrange(chr) %>%
    dplyr::filter(!(as.character(chr) %in% as.character(c(10, 12, 14, 16, 18))))

  my_colors <- c(
    "A" = "#39BC95", "B" = "#008059",
    "C" = "#F77D20", "D" = "#BB4100",
    "E" = "#938ED1", "F" = "#575295"
  )

  n_age <- dplyr::n_distinct(plot_data$age_years)
  n_plat <- dplyr::n_distinct(plot_data$Platform)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pos_cum, y = lod, color = color_var)) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(
        tooltip = tooltip,
        data_id = .safe_id(name_use),
        onclick = onclick
      ),
      size = 1.2, alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = my_colors, guide = "none", na.value = "grey70") +
    ggplot2::scale_x_continuous(labels = axis_set$chr, breaks = axis_set$center) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Platform), cols = ggplot2::vars(age_years),
      scales = "free_y"
    ) +
    theme_docr +
    ggplot2::labs(x = "Chromosome", y = "LOD")

  ggiraph::girafe(
    ggobj = p,
    width_svg = n_age * 4,
    height_svg = n_plat * 3,
    options = list(
      ggiraph::opts_hover(css = "transform:scale(2.5); transform-box:fill-box; transform-origin:center;")
    )
  )
}
