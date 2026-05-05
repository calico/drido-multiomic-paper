# For keeping features that have at least one compound significant in the dataset
# Assumes that padj and feature_id are in the dataset
keep_only_sig_features <- function(df,
                                   feat_col = "name_use",
                                   pval_col = "padj",
                                   p_cutoff = 0.05,
                                   n = 1) {
  df <- df %>%
    dplyr::mutate(Significant = ifelse(!!rlang::sym(pval_col) < p_cutoff, TRUE, FALSE)) %>%
    dplyr::group_by(!!rlang::sym(feat_col)) %>%
    dplyr::mutate(n_significant = sum(Significant)) %>% # Count TRUEs per group
    dplyr::filter(n_significant >= n) %>% # Filter based on the count
    dplyr::select(-n_significant) # Remove helper column

  return(df)
}


convert_string <- function(x) {
  x <- gsub("\\.M012.*", " (+)", x)
  x <- gsub("\\.M013.*", " (-)", x)
  x <- gsub("\\.M014.*", " (+)", x)
  x <- gsub("\\.M015.*", " (-)", x)
  return(x)
}

# Function to make unique names in a dataframe
make_unique_names <- function(df, col_name) {
  df$new_name <- convert_string(df[[col_name]])
  counts <- table(df$new_name)
  duplicates <- names(counts[counts > 1])

  for (dup in duplicates) {
    rows <- which(df$new_name == dup)
    df$new_name[rows] <- paste0(df$new_name[rows], " (", 1:length(rows), ")")
  }
  return(df)
}


# Function reads and cleans phenotype data
docr_read_phenotype_data <- function(phenotype_data_file_path_csv) {
  phenotype_data <- read.csv(phenotype_data_file_path_csv) %>%
    tidyr::pivot_longer(
      cols = grep("_", colnames(.), value = TRUE),
      values_to = "measurement"
    ) %>%
    dplyr::filter(!is.na(measurement)) %>%
    dplyr::rename(mouse_id = MouseID) %>%
    tidyr::separate(name,
      into = c("timepoint", "pheno_group", "phenotype"),
      sep = "_"
    ) %>%
    dplyr::group_by(pheno_group, timepoint, mouse_id) %>%
    dplyr::mutate(
      age_in_days = ifelse(any(phenotype == "AgeInDays"),
        measurement[phenotype == "AgeInDays"][1],
        NA
      ),
      days_remaining = SurvDays - age_in_days,
      bw_test = ifelse(any(phenotype == "BWTest"),
        measurement[phenotype == "BWTest"][1],
        NA
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      surv_days = SurvDays,
      diet_assignment = Diet
    ) %>%
    dplyr::filter(!phenotype %in% c("BWTest", "AgeInDays")) %>%
    dplyr::mutate(
      trait_id = paste0(pheno_group, "_", phenotype),
      diet = ifelse(age_in_days < 200, "AL", diet_assignment),
      diet = factor(diet, levels = c("AL", "1D", "2D", "20", "40")),
      diet_assignment = factor(diet_assignment, levels = c("AL", "1D", "2D", "20", "40")),
      PLL = as.numeric((surv_days - days_remaining) / surv_days),
      fasting = NA
    ) %>%
    dplyr::select(-pheno_group, -phenotype) %>%
    dplyr::group_by(trait_id) %>%
    dplyr::mutate(trait_value = scale(measurement, center = T, scale = F)) %>%
    dplyr::ungroup() %>%
    # Modifications for doing correlates against metabolites
    dplyr::mutate(
      age_years = dplyr::case_when(grepl("^Met_|^FACS_", trait_id) & timepoint == "Y2" ~ "year2",
        grepl("^Met_|^FACS_", trait_id) & timepoint == "Y3" ~ "year3",
        grepl("^Frailty_|^Grip_", trait_id) & timepoint == "Y2A" ~ "year2",
        grepl("^Frailty_|^Grip_", trait_id) & timepoint == "Y3A" ~ "year3",
        grepl("^Wheel_|^Void_|^Rota_|^PIXI_|^Echo_|^AS_|^CBC_", trait_id) & timepoint == "Y1" ~ "year2",
        grepl("^Wheel_|^Void_|^Rota_|^PIXI_|^Echo_|^AS_|^CBC_", trait_id) & timepoint == "Y2" ~ "year3",
        .default = NA
      ),
      modality = "physiological"
    )

  return(phenotype_data)
}


# Reads in metabolomics, proteomics, lipidomics, and phenotype data
# Removes bad features IDs from metabolomics
# Add good name IDs for plotting via name_conversion_key
# Add additional metadata
# Refactors data where necessary
docr_make_final_data <- function(metabolomics_data_filepath,
                                 proteomics_data_filepath,
                                 lipidomics_data_filepath,
                                 phenotype_data_filepath,
                                 name_conversion_key) {
  # Identified duplicate peaks for removal
  bad_duplicate_feature_ids <- docr_bad_duplicate_features()

  # Load plasma molecular data
  metab_data <- readRDS(metabolomics_data_filepath) %>%
    dplyr::filter(!(feature_id %in% bad_duplicate_feature_ids))
  lipid_data <- readRDS(lipidomics_data_filepath)
  prot_data <- readRDS(proteomics_data_filepath)

  ### Load and clean physiological data
  phenotype_data <- docr_read_phenotype_data(phenotype_data_filepath)

  ### Combine and clean molecular data
  all_molecular_data <- metab_data %>%
    dplyr::mutate(modality = "metabolomics") %>%
    dplyr::filter(!(feature_id %in% bad_duplicate_feature_ids)) %>%
    dplyr::bind_rows(prot_data %>%
      dplyr::mutate(modality = "proteomics")) %>%
    dplyr::bind_rows(lipid_data %>%
      dplyr::mutate(modality = "lipidomics")) %>%
    dplyr::mutate(
      days_remaining = as.numeric(delta_date),
      fasting = ifelse(diet %in% c("20") & weekday_collection == "Monday",
        "Fast20",
        ifelse(diet %in% c("40") & weekday_collection == "Monday",
          "Fast40",
          "No"
        )
      ),
      fasting = factor(fasting, levels = c("No", "Fast20", "Fast40"))
    ) %>%
    dplyr::mutate(
      trait_id = feature_id,
      bw_test = BW_Raw
    ) %>%
    dplyr::group_by(trait_id) %>%
    dplyr::mutate(trait_value = scale(norm_abundance, center = T, scale = F)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., name_conversion_use, by = c("trait_id" = "feature_id")) %>%
    # Metabolomics internal standards are removed here
    # Feature duplicates are removed here
    dplyr::filter(!is.na(name_use))

  ### Combine all data
  cols_use <- intersect(colnames(all_molecular_data), colnames(phenotype_data))
  data_use <- dplyr::bind_rows(
    all_molecular_data %>%
      dplyr::select(tidyselect::all_of(cols_use), "name_use"),
    phenotype_data %>%
      dplyr::select(tidyselect::all_of(cols_use))
  ) %>%
    dplyr::mutate(across(c("mouse_id", "diet", "fasting"), ~ factor(.x))) %>%
    dplyr::mutate(diet_fasting = paste0(diet, "-", fasting)) %>%
    dplyr::mutate(diet_fasting = stringr::str_trunc(diet_fasting, 7, ellipsis = "")) %>%
    dplyr::mutate(diet_fasting = factor(diet_fasting, levels = c("AL-No", "1D-No", "2D-No", "20-No", "40-No", "20-Fast", "40-Fast"))) %>%
    dplyr::mutate(Age = ifelse(age_years == "year1", "Year 1", ifelse(age_years == "year2", "Year 2", "Year 3"))) %>%
    dplyr::mutate(
      diet = factor(diet, levels = c("AL", "1D", "2D", "20", "40")),
      diet_assignment = factor(diet_assignment, levels = c("AL", "1D", "2D", "20", "40"))
    )

  invisible(gc())
  return(data_use)
}


# DOCR name clean helper function to remove trembl and swiss-prots
# that have the same gene identification, and adds ionization annotation
# if a single compound has been detected in multiple methods
docr_name_clean_helper <- function(name_conversion_use) {
  return_df <- name_conversion_use %>%
    dplyr::group_by(name_use) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::mutate(
      starts_with_sp = stringr::str_detect(feature_id, "^sp__"),
      starts_with_tr = stringr::str_detect(feature_id, "^tr__"),
      has_M012A = stringr::str_detect(feature_id, "M012A"),
      has_M013A = stringr::str_detect(feature_id, "M013A"),
      feature_id_length = stringr::str_length(feature_id)
    ) %>%
    # Prioritize sp__ and shorter IDs if there are duplicate trembl IDs
    dplyr::arrange(name_use, starts_with_sp, -feature_id_length) %>%
    # Remove tr__ if sp__ exists
    dplyr::filter(!(starts_with_tr & any(starts_with_sp))) %>%
    dplyr::group_by(name_use) %>%
    # Remove longer if both tr__ and sp__ exist
    dplyr::filter(!((starts_with_tr | starts_with_sp) & feature_id_length > min(feature_id_length))) %>%
    # If a lipid exists in muliple method, annotate metabolomics method
    dplyr::mutate(
      name_use = dplyr::case_when(
        has_M012A ~ paste0(name_use, " (+)"),
        has_M013A ~ paste0(name_use, " (-)"),
        TRUE ~ name_use
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-starts_with_sp, -starts_with_tr, -has_M012A, -has_M013A, -feature_id_length)
  return(return_df)
}


# Function cleans metabolomics and proteomics compound names
# Accepts a metabolomics dataframe lm_res containing feature_ids and modality column
# Accepts a proteimcs dataframe protein_names with uniprot IDs and gene names
docr_clean_compound_names <- function(lm_res,
                                      protein_names) {
  # Identified duplicate peaks for removal
  bad_duplicate_feature_ids <- docr_bad_duplicate_features()

  metab_name_conversion <- lm_res %>%
    # Remove IS, bad duplicate feature ids (manually annotations), and unknowns
    dplyr::filter(
      !grepl("-[d]\\d|15n|13c", tolower(feature_id)),
      !grepl("^unk ", feature_id),
      !(feature_id %in% bad_duplicate_feature_ids),
      modality == "metabolomics"
    ) %>%
    # Rename features that were given incorrect annotation
    dplyr::mutate(compoundName = dplyr::case_when(feature_id == "Tridecanedioic acid.M013A.7309" ~ "Undecanedicarboxylic acid",
      feature_id == "6-Hydroxynicotinic acid.M013A.2176" ~ "3-Hydroxypicolinic acid",
      feature_id == "N-acetylleucine.M012A.1418" ~ "Hexanoylglycine",
      feature_id == "Galactose.M013A.4469" ~ "Glucose",
      .default = compoundName
    )) %>%
    # Select feature_ids and assign unique names
    dplyr::group_by(feature_id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(feature_id, compoundName) %>%
    make_unique_names(., "feature_id") %>%
    dplyr::group_by(compoundName) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      name_use = ifelse(n == 1, compoundName, new_name),
      name_use = stringr::str_replace(name_use, "[+-][0-9.:;,]+$", ""),
      name_use = ifelse(stringr::str_count(name_use, "[A-Z]") >= 4,
        stringr::str_to_title(name_use),
        name_use
      )
    ) %>%
    dplyr::select(feature_id, name_use) %>%
    # Dealing with a few manually, some misannotated compounds here too
    dplyr::mutate(name_use = dplyr::case_when(grepl("INOSINE.M013A.8377", feature_id) ~ "Inosine (-)",
      grepl("Inosine.M012A.2988", feature_id) ~ "Inosine (+)",
      grepl("PYRIDOXAL+20.0:10", feature_id) ~ "Pyridoxal (+))",
      grepl("Pyridoxal.M012A.128", feature_id) ~ "Pyridoxal (-)",
      grepl("SARCOSINE+20.0:10.0:180", feature_id) ~ "Sarcosine (+)",
      grepl("SARCOSINE-20.0,5", feature_id) ~ "Sarcosine (-)",
      grepl("N-acetylleucine.M012A.1407", feature_id) ~ "N-acetylleucine (+)",
      grepl("Trigonelline \\(\\+\\)", name_use) ~ "N-methylnicotinate (+)",
      grepl("Trigonelline \\(-\\)", name_use) ~ "N-methylnicotinate (-)",
      grepl("N-methylglutamic acid.M012A.1193", feature_id) ~ "Aminoadipic acid (+)",
      grepl("Na-Acetylcitrulline; LC-tDDA; CE10.M013A.6140", feature_id) ~ "Acetylcitrulline (-)",
      grepl("5-Aminolevulinic acid.M012A.596", feature_id) ~ "N-acetylalanine (+)",
      grepl("5-AMINOLEVULINIC ACID.M013A.1756", feature_id) ~ "N-acetylalanine (-)",
      grepl("N-acetylalanine.M012A.605", feature_id) ~ "4-hydroxyproline (+)",
      grepl("N-acetylalanine.M013A.1757", feature_id) ~ "4-hydroxyproline (-)",
      .default = name_use
    ))

  # Make proteins unique and readable
  prot_name_conversion <- protein_names %>%
    dplyr::full_join(
      lm_res %>%
        dplyr::filter(modality == "proteomics") %>%
        dplyr::select(feature_id, compoundName) %>%
        dplyr::mutate(uniprot = sub(".*?__([^_]+)_.*", "\\1", feature_id)) %>%
        dplyr::distinct(),
      by = "uniprot",
      relationship = "many-to-many"
    ) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !(gene == "" & n > 1),
      !is.na(feature_id)
    ) %>%
    dplyr::select(feature_id, compoundName, gene) %>%
    dplyr::mutate(
      name_use = gene,
      name_type = ifelse(name_use != "", "gene", "protein"),
      name_use = ifelse(name_use == "" & grepl("^HVM", compoundName),
        paste0("Ighv-", compoundName),
        name_use
      ),
      name_use = ifelse(name_use == "" & grepl("^KV", compoundName),
        paste0("Igkv-", compoundName),
        name_use
      ),
      name_use = ifelse(name_use == "" & grepl("^LV", compoundName),
        paste0("Iglv-", compoundName),
        name_use
      ),
      name_use = ifelse(name_use == "", compoundName, name_use)
    )

  # Combine all modalities
  name_conversion_use <- metab_name_conversion %>%
    dplyr::select(name_use, feature_id) %>%
    dplyr::bind_rows(prot_name_conversion %>%
      dplyr::select(name_use, feature_id)) %>%
    dplyr::bind_rows(lm_res %>%
      dplyr::filter(modality == "lipidomics") %>%
      dplyr::select(name_use = compoundName, feature_id) %>%
      dplyr::mutate(name_use = ifelse(!grepl("^Ceramide|^HexCer|^SM", name_use),
        gsub("/", "_", name_use),
        name_use
      )) %>%
      dplyr::distinct())

  # After combining, some lipids and metabolites have the same name
  deduplicate_again <- name_conversion_use %>%
    docr_name_clean_helper()

  # Remove these from name_conversion_use, then add back corrected
  name_conversion_use <- name_conversion_use %>%
    dplyr::filter(!(name_use %in% deduplicate_again$name_use)) %>%
    dplyr::bind_rows(deduplicate_again)

  # Check AGAIN for some names that either did not get ion annotation, or
  # duplicates in proteomics due to capitalization differences
  fix_again <- name_conversion_use %>%
    dplyr::mutate(name_use = stringr::str_to_title(name_use)) %>%
    docr_name_clean_helper()

  # Remove these from name_conversion_use, add back corrected, fix
  # last protein to gene names needed for entrez ID conversion
  name_conversion_use <- name_conversion_use %>%
    dplyr::filter(
      !(stringr::str_to_title(name_use) %in% fix_again$name_use),
      !(feature_id %in% fix_again$feature_id)
    ) %>%
    dplyr::bind_rows(fix_again) %>%
    dplyr::mutate(name_use = dplyr::case_when(
      name_use == "C5" ~ "Hc",
      name_use == "Ca2" ~ "Car2",
      name_use == "C4bpa" ~ "C4bp",
      name_use == "Amy2" ~ "Amy2a5",
      name_use == "Fcn1" ~ "Fcna",
      name_use == "Ica" ~ "Inhca",
      name_use == "Igh-1a" ~ "Ighg2a",
      name_use == "Gpi" ~ "Gpi1",
      name_use == "Txn" ~ "Txn1",
      name_use == "Ca1" ~ "Car1",
      name_use == "Ca3" ~ "Car3",
      name_use == "Tf" ~ "Trf",
      TRUE ~ name_use
    ))

  return(name_conversion_use)
}


# Min max scale PLL trajectories so that all trajectories are between 0 and 1
docr_scale_trajectory <- function(df) {
  min_val <- min(df$fit)
  max_val <- max(df$fit)
  df$scaled_fit <- (df$fit - min_val) / (max_val - min_val)
  df$scaled_lower_ci <- (df$lower.95.ci - min_val) / (max_val - min_val)
  df$scaled_upper_ci <- (df$upper.95.ci - min_val) / (max_val - min_val)
  return(df)
}


# Generate a Scatter Plot with Conditional Labels
docr_plot_correlation_scatter <- function(data_corcor,
                                          x_col_name,
                                          y_col_name,
                                          label_col_name,
                                          label_threshold = 0.8,
                                          repel_params = list(
                                            box.padding = unit(0.6, "lines"),
                                            point.padding = unit(0.4, "lines"),
                                            max.overlaps = 20,
                                            min.segment.length = unit(0.3, "lines"),
                                            force = 15,
                                            size = 3
                                          )) {
  if (!is.data.frame(data_corcor) && !is.matrix(data_corcor)) {
    stop("'data_corcor' must be a data frame or matrix.")
  }
  if (!x_col_name %in% colnames(data_corcor)) {
    stop("Specified x_col_name '", x_col_name, "' not found in data_corcor columns.")
  }
  if (!y_col_name %in% colnames(data_corcor)) {
    stop("Specified y_col_name '", y_col_name, "' not found in data_corcor columns.")
  }

  gg <- ggplot(
    data = data_corcor,
    aes(x = .data[[x_col_name]], y = .data[[y_col_name]])
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", formula = y ~ x) +

    # Apply geom_text_repel only to a subset, using .data[[var_name]] in filter
    ggrepel::geom_label_repel(
      data = . %>% dplyr::filter(
        abs(.data[[x_col_name]]) > label_threshold | abs(.data[[y_col_name]]) > label_threshold
      ),
      aes(label = .data[[label_col_name]]),
      box.padding = unit(0.6, "lines"),
      point.padding = unit(0.4, "lines"),
      max.overlaps = 20,
      min.segment.length = unit(0.3, "lines"),
      force = 15,
      size = 3
    ) +
    theme_classic() +
    labs(
      title = paste("Correlation:", y_col_name, "vs", x_col_name),
      x = x_col_name,
      y = y_col_name
    )

  return(gg)
}


# Wrapper function for hclust to be used with clusGap
docr_hclust_wrapper <- function(x, k) {
  # Absolute distance matrix
  d <- as.dist(1 - abs(x))
  hc <- stats::hclust(d, method = "complete")
  clusters <- stats::cutree(hc, k = k)
  return(list(cluster = clusters))
}


# Function makes a dataframe to use with docr_ggplot_stats_label to print
# R2 and pvalue in faceted dataframe
docr_facet_stats <- function(df,
                             value_x,
                             value_y,
                             facet_1 = NULL,
                             facet_2 = NULL) {
  if (is.null(facet_1) && is.null(facet_2)) {
    lm_stats <- df
  } else if (!is.null(facet_1) && is.null(facet_2)) {
    lm_stats <- df %>%
      dplyr::group_by(!!rlang::sym(facet_1))
  } else if (!is.null(facet_1) && !is.null(facet_2)) {
    lm_stats <- df %>%
      dplyr::group_by(
        !!rlang::sym(facet_1),
        !!rlang::sym(facet_2)
      )
  } else {
    stop("Faceting error")
  }
  lm_stats <- lm_stats %>%
    dplyr::summarize(
      model = list(lm(as.formula(paste0(value_y, " ~ ", value_x)))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r.squared = purrr::map_dbl(model, ~ summary(.x)$adj.r.squared),
      r = purrr::map_dbl(model, ~ sqrt(summary(.x)$r.squared) * sign(coef(.x)[2])),
      f.statistic_p_value = purrr::map_dbl(model, ~ pf(summary(.x)$fstatistic[1],
        summary(.x)$fstatistic[2],
        summary(.x)$fstatistic[3],
        lower.tail = FALSE
      ))
    ) %>%
    dplyr::mutate(
      r.squared_label = paste("Adj.R² =", round(r.squared, 2)),
      r_label = paste("R =", round(r, 2)),
      p_value_label = paste("p =", format.pval(f.statistic_p_value, digits = 2)),
      combined_label = paste0(p_value_label, "\n", r_label, "\n", r.squared_label)
    )
  return(lm_stats)
}


# ggplot object for printing text in facets defined by facet_x and facet_y
# in docr_facet_stats
docr_ggplot_stats_label <- function(df,
                                    x = Inf,
                                    y = -Inf,
                                    hjust = 1.1,
                                    vjust = -0.4,
                                    size = 2.7,
                                    color = "red",
                                    lineheight = 0.8) {
  geom_return <- ggplot2::geom_text(
    data = df,
    aes(x = x, y = y, label = combined_label),
    hjust = hjust, vjust = vjust, size = size, color = color,
    lineheight = lineheight
  )
  return(geom_return)
}


#' Pretty Heatmaps with optional p-values from tidy data
#'
#' @description
#' This is a wrapper function for \code{pheatmap}. This function takes a tidy
#' dataframe as an input, uses a colorblind-friendly metabolomics standard color palette,
#' and centers the color palette on zero as white (ideal for centered data).
#' Value significance as stars or alternative text can be added by indicating which column
#' contains the desired p-values or text. This function accepts all \code{pheatmap}
#' arguments for customizing the figure or overriding defaults.
#'
#' @param heatmap_df tidy dataframe for creating pheatmap
#' @param x_axis_var column name of value on x-axis, as string
#' @param y_axis_var column name of value on y-axis, as string
#' @param value_var column name of value to represent in heatmap, as numeric
#' @param value_aggregate function for value aggregation if not unique, as string
#' @param display_var optional column name of values to print on heatmap, such as
#' pvalues
#' @param display_aggregate function for display aggregation if not unique, as string
#' @param pvalue_stars logical which indicates if \code{display_var} is a list of
#' pvalues to print as stars on the heatmap. Defaults to \code{FALSE}
#' @param sig_symbol if \code{pvalue_stars = TRUE}, symbol to use for significance
#' @param sig_numbers if \code{pvalue_stars = TRUE}, numbers to use for pvalue cutoffs
#' @param centered_data logical to indicate if color palette should be centered over
#' white. Defaults to \code{TRUE}
#' @param ... optional additional arguments to call from pheatmap formals
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @returns A pheatmap figure
#'
#' @export
docr_pheatmap <- function(heatmap_df,
                          x_axis_var = "model_term",
                          y_axis_var = "feature_id",
                          value_var = "coef",
                          value_aggregate = NULL,
                          display_var = NULL,
                          display_aggregate = NULL,
                          pvalue_stars = FALSE,
                          sig_symbol = "\u00D7",
                          sig_numbers = c(0.0001, 0.001, 0.01, 0.05),
                          centered_data = TRUE,
                          ...) {
  # Check dataframe
  checkmate::assertString(x_axis_var)
  checkmate::assertString(y_axis_var)
  checkmate::assertString(value_var)
  if (!(x_axis_var %in% colnames(heatmap_df))) {
    stop("\"x_axis_var\": ", x_axis_var, ", not present in heatmap_df dataframe")
  }
  if (!(y_axis_var %in% colnames(heatmap_df))) {
    stop("\"y_axis_var\": ", y_axis_var, ", not present in heatmap_df dataframe")
  }
  if (!(value_var %in% colnames(heatmap_df))) {
    stop("\"value_var\": ", value_var, ", not present in heatmap_df dataframe")
  }

  # Check pheatmap arguments
  dots <- list(...)
  method_args <- dots[intersect(names(dots), names(formals(pheatmap::pheatmap)))]
  unused_args <- setdiff(names(dots), names(formals(pheatmap::pheatmap)))
  if (length(unused_args) > 0) {
    warning(glue::glue("{length(unused_args)} passed arguments could not be used by pheatmap:\n{paste(unused_args, collapse = ', ')} "))
  }

  # Check if entries are unique
  unique_check <- docr_check_unique_entries(heatmap_df, x_axis_var, y_axis_var)

  ### CREATE HEATMAP VALUE MATRIX
  if (!unique_check) {
    if (is.null(value_aggregate)) {
      stop("Heatmap values are not unique, and no aggregation function supplied to \"value_aggregate\"")
    } else {
      message("Heatmap values are not unique and will be aggregated using function supplied in \"value_aggregate\"")
      value_aggregate <- eval(parse(text = value_aggregate))
    }
  }

  heatmap_df_wide <- heatmap_df %>%
    dplyr::select(tidyselect::all_of(c(x_axis_var, y_axis_var, value_var))) %>%
    tidyr::pivot_wider(
      names_from = x_axis_var,
      values_from = value_var,
      values_fn = value_aggregate
    ) %>%
    tibble::column_to_rownames(y_axis_var)

  # Column order can be set by specifying levels; drop levels to prevent errors
  heatmap_df <- heatmap_df %>%
    dplyr::mutate(`:=`(!!rlang::sym(x_axis_var), droplevels(as.factor(!!rlang::sym(x_axis_var)))))
  heatmap_df_wide <- heatmap_df_wide[, levels(heatmap_df[[x_axis_var]])]

  ### CREATE COLORS AND SCALES
  if (centered_data || "scale" %in% names(method_args) && method_args[["scale"]] %in% c("row", "column")) {
    if ("scale" %in% names(method_args)) {
      if (method_args[["scale"]] == "row") {
        temp_hm <- apply(heatmap_df_wide, 1, function(x) scale(x, center = TRUE, scale = T))
      } else if (method_args[["scale"]] == "column") {
        temp_hm <- apply(heatmap_df_wide, 2, function(x) scale(x, center = TRUE, scale = T))
      } else {
        warning(paste0("Issue with scale argument; not scaling colors"))
        temp_hm <- heatmap_df_wide
      }
    } else {
      temp_hm <- heatmap_df_wide
    }

    cols_and_breaks <- docr_make_centered_ramped_palette(temp_hm)
  } else {
    cols_and_breaks <- list(
      heatmap_colors = grDevices::colorRampPalette(docr_get_high_low_colors())(100),
      heatmap_breaks = NA
    )
  }

  # Identify maximum value on heatmap which will center heatmap around zero
  # This can be overridden by supplying breaks (length 100) in the function
  # mx <- max(abs(heatmap_df_wide), na.rm = TRUE)


  ### CREATE DISPLAY MATRIX
  if (!is.null(display_var) && display_var %in% colnames(heatmap_df)) {
    if (!unique_check) {
      if (is.null(display_aggregate)) {
        stop("Display values are not unique, and no aggregation function supplied to \"display_aggregate\"")
      } else {
        message("Display values are not unique and will be aggregated using function supplied in \"display_aggregate\"")
        display_aggregate <- eval(parse(text = display_aggregate))
      }
    }

    # Create pvalue stars if indicated
    if (pvalue_stars) {
      heatmap_df <- heatmap_df %>%
        docr_make_significant_pvalues(
          pval_var = display_var,
          sig_symbol = sig_symbol,
          sig_numbers = sig_numbers
        ) %>%
        dplyr::mutate(`:=`(!!rlang::sym(display_var), .data$significance))
    }

    disp_matrix <- heatmap_df %>%
      dplyr::select(tidyselect::all_of(c(x_axis_var, y_axis_var, display_var))) %>%
      tidyr::pivot_wider(
        names_from = x_axis_var,
        values_from = display_var,
        values_fn = display_aggregate
      ) %>%
      tibble::column_to_rownames(y_axis_var)

    # Set column names in order of provided levels
    disp_matrix <- disp_matrix[, levels(as.factor(heatmap_df[[x_axis_var]])), ]

    # If display_var is null, assign to disp_matrix to default for display_numbers (F)
  } else {
    disp_matrix <- FALSE
  }


  # Define default arguments from pheatmap
  default_args <- list(
    cluster_rows = TRUE,
    cluster_col = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    border_color = NA,
    fontsize = 9,
    fontsize_row = 5,
    color = cols_and_breaks$heatmap_colors,
    display_numbers = disp_matrix,
    # breaks = seq(-mx, mx, length.out = 100)
    breaks = cols_and_breaks$heatmap_breaks
  )

  # Only use default arguments that are not supplied as function arguments
  default_args_use <- default_args[setdiff(names(default_args), names(method_args))]
  method_args_use <- append(default_args_use, method_args)

  # Run pheatmap
  ph <- do.call(
    pheatmap::pheatmap,
    append(
      list(mat = heatmap_df_wide),
      method_args_use
    )
  )
  return(ph)
}


#' Make Pretty Heatmap Annotation Key
#'
#' @description
#' This function creates a the annotation dataframe for use with \code{annotation_col}
#' or \code{annotation_row} in \code{\link{docr_pheatmap}}. Given the same \code{heatmap_df}
#' and \code{x_axis_var} or \code{y_axis_var}, this function will make a compatible annotation
#' dataframe.
#'
#' @param heatmap_df tidy dataframe for creating pheatmap
#' @param axis_var column name of values to use for y-axis (for \code{pheatmap::annotation_row})
#' or x-axis (for \code{pheatmap::annotation_col}), as string
#' @param ... arguments to create annotation dataframe. Keys can be entered either
#' as named lists, or as a character vector. Function works by detecting string matches
#' on names, and assigning values to the string match. Values that do not
#' match are given an empty value " " that can be leveraged as a white color automatically
#' with \code{\link{docr_pheatmap_colors}}. See example.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table :=
#'
#' @returns A dataframe with group keys from \code{...} as factor columns
#'
#' @export
docr_pheatmap_annotations <- function(heatmap_df,
                                      axis_var = "model_term",
                                      ...) {
  # Check dataframe
  checkmate::assertString(axis_var)
  if (!(axis_var %in% colnames(heatmap_df))) {
    stop("\"axis_var\":", axis_var, ", not present in heatmap_df dataframe")
  }

  # Extract unique axis names
  annotation_col <- heatmap_df %>%
    dplyr::rename(a_var = axis_var) %>%
    dplyr::select(.data$a_var) %>%
    dplyr::distinct()

  dots <- list(...)
  for (n in names(dots)) {
    # Assign to blank on default
    annotation_col <- annotation_col %>%
      dplyr::mutate(!!rlang::sym(n) := " ")

    # If no names are provided, assign string as name
    n_options <- dots[[n]]
    if (is.null(names(n_options))) {
      names(n_options) <- n_options
    }

    # For each string option in a grouping key, assign to named value
    for (n_opt_name in names(n_options)) {
      n_opt_value <- n_options[[n_opt_name]]
      annotation_col <- annotation_col %>%
        dplyr::mutate(!!rlang::sym(n) := ifelse(grepl(n_opt_name, .data$a_var, fixed = T),
          n_opt_value,
          !!rlang::sym(n)
        ))
    }

    annotation_col <- annotation_col %>%
      dplyr::mutate(!!rlang::sym(n) := factor(!!rlang::sym(n), levels = c(unname(n_options), " ")))
  }

  annotation_col <- annotation_col %>%
    tibble::column_to_rownames("a_var")

  return(annotation_col)
}

#' Make Pretty Heatmap Annotation Colors
#'
#' @description
#' This function creates the annotation colors based off of the annotation dataframe
#' output from \code{\link{docr_pheatmap_annotations}} for use with
#' \code{\link{make_super_pheatmap}}. This function assigns white colors
#' to unmatched values in \code{\link{docr_pheatmap_annotations}} that are
#' assigned " ". The user can either provide desired color schemes as additional
#' arguments, or the function will automatically generate color schemes.
#' @param annotation_df dataframe output from \code{\link{docr_pheatmap_annotations}}
#' @param ... optional arguments to hardcode desired colors, as named lists where
#' item names match column names and variables in \code{annotation_df}.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A named list
#'
#' @export
docr_pheatmap_colors <- function(annotation_df,
                                 ...) {
  # Only make colors on annotation dataframe column names that are not provided
  annotation_color_list <- list(...)
  make_colors <- setdiff(colnames(annotation_df), names(annotation_color_list))

  # These will get added to a new list
  annotation_color_list_append <- list()
  if (length(make_colors) > 0) {
    possible_color_scales <- c("Blues", "Reds", "Greens", "Purples", "Oranges")

    # From preset color scales in RColorBrewer, make palettes
    for (i in 1:length(make_colors)) {
      mc <- make_colors[i]
      annotation_df[[mc]] <- factor(annotation_df[[mc]])
      mcx <- levels(annotation_df[[mc]])[!grepl("^\\s+$", levels(annotation_df[[mc]]))]

      # Index from back of palette for key groups with n < 3 -- this will select darker colors
      mcx_palette <- utils::tail(
        suppressWarnings(RColorBrewer::brewer.pal(length(mcx), name = possible_color_scales[i])),
        n = length(mcx)
      )
      names(mcx_palette) <- mcx
      annotation_color_list_append[[mc]] <- mcx_palette
    }
  }

  # Append to list of supplied colors
  annotation_color_list <- append(
    annotation_color_list,
    annotation_color_list_append
  )

  # For empty character values, assign white (looks clean and invisible on pheatmap)
  annotation_color_list <- lapply(
    annotation_color_list,
    function(x) {
      n <- c(x, ` ` = "white")
      return(n)
    }
  )

  return(annotation_color_list)
}


docr_make_significant_pvalues <- function(df,
                                          pval_var,
                                          star_var = "significance",
                                          sig_symbol = "×",
                                          sig_numbers = c(1e-04, 0.001, 0.01, 0.05)) {
  checkmate::assertString(pval_var)
  checkmate::assertString(star_var)
  checkmate::assertString(sig_symbol)
  if (!(pval_var %in% colnames(df))) {
    stop("\"pval_var\":", pval_var, ", not present in df")
  }
  if (!all(dplyr::between(df[[pval_var]], 0, 1))) {
    stop("All values in \"pval_var\" must be between 0 and 1")
  }
  if (!all(dplyr::between(sig_numbers, 0, 1))) {
    stop("All values in \"sig_numbers\" must be between 0 and 1")
  }
  ordered_sig_numbers <- sig_numbers[order(sig_numbers)]
  sig_symbols <- sapply(rev(order(ordered_sig_numbers)), function(x) {
    paste(rep(
      sig_symbol,
      x
    ), collapse = "")
  })
  df <- df %>%
    dplyr::rename(p_var = pval_var) %>%
    dplyr::mutate(`:=`(
      !!rlang::sym(star_var),
      NA
    ))
  for (i in 1:length(sig_symbols)) {
    df <- df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(`:=`(
        !!rlang::sym(star_var),
        dplyr::case_when(
          p_var < ordered_sig_numbers[i] &&
            is.na(!!rlang::sym(star_var)) ~ sig_symbols[i],
          .default = !!rlang::sym(star_var)
        )
      ))
  }
  df[is.na(df[[star_var]]), star_var] <- ""
  df <- df %>% dplyr::rename(`:=`(!!rlang::sym(pval_var), "p_var"))
  return(df)
}


docr_make_centered_ramped_palette <- function(numeric_df,
                                              break_count = 100,
                                              center_value = 0) {
  numeric_df <- sapply(numeric_df, as.numeric)
  df_max <- max(numeric_df, na.rm = T)
  df_min <- min(numeric_df, na.rm = T)
  df_range <- df_max - df_min
  percent_low <- abs(df_min - center_value) / df_range
  percent_high <- abs(df_max - center_value) / df_range
  heatmap_breaks <- c(seq(df_min, df_max, length.out = break_count))
  metab_colors <- docr_metabolomics_palette(center_white = T)
  colors_low <- (grDevices::colorRampPalette(metab_colors[c(
    "Low",
    "Mid"
  )]))(round(max(percent_high, percent_low) * break_count))
  colors_low <- utils::tail(colors_low, n = round(percent_low *
    break_count))
  colors_high <- (grDevices::colorRampPalette(metab_colors[c(
    "Mid",
    "High"
  )]))(round(max(percent_high, percent_low) * break_count))
  colors_high <- utils::head(colors_high, n = round(percent_high *
    break_count))
  heatmap_colors <- c(colors_low, colors_high)
  return(list(heatmap_colors = heatmap_colors, heatmap_breaks = heatmap_breaks))
}


# Before making heatmap, confirm unique
docr_check_unique_entries <- function(df, x_var, y_var) {
  unique_combinations <- df %>%
    dplyr::distinct(dplyr::across(c(x_var, y_var))) %>%
    nrow()

  if (unique_combinations == nrow(df)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


docr_metabolomics_palette <- function(center_white = FALSE) {
  if (center_white) {
    return(c(Low = "#21336B", Mid = "white", High = "#D63B24"))
  } else {
    return(c(Low = "#21336B", High = "#D63B24"))
  }
}


docr_find_elbow_smooth <- function(df, x_col, y_col, span = 0.5) {
  x <- df[[x_col]]
  y <- df[[y_col]]

  fit <- loess(y ~ x, span = span)

  x_seq <- seq(min(x), max(x), length.out = 500)
  y_hat <- predict(fit, newdata = data.frame(x = x_seq))

  # Remove NAs from edges before normalizing
  valid <- !is.na(y_hat)
  x_seq <- x_seq[valid]
  y_hat <- y_hat[valid]

  x_norm <- (x_seq - min(x_seq)) / (max(x_seq) - min(x_seq))
  y_norm <- (y_hat - min(y_hat)) / (max(y_hat) - min(y_hat))

  x1 <- x_norm[1]
  y1 <- y_norm[1]
  x2 <- x_norm[length(x_norm)]
  y2 <- y_norm[length(y_norm)]

  perp_dist <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)

  elbow_x <- x_seq[which.max(perp_dist)]

  return(elbow_x)
}


# Function to make models readable
make_names_pretty <- function(s) {
  s <- gsub("norm_abundance", "Compound", s)
  s <- gsub("age_years", "Age", s)
  s <- gsub("diet", "Diet", s)
  s <- gsub("bw_test", "Bodyweight", s)
  s <- gsub("BW_Loess", "Bodyweight", s)
  s <- gsub("surv_years", "Lifespan", s)
  s <- gsub("surv_days", "Lifespan", s)
  s <- gsub("surv_days.surv", "Lifespan", s)
  s <- gsub("days_remaining", "LifeRemaining", s)
  s <- gsub("generation_wave", "GenWave", s)
  s <- gsub("weekday_collection", "Weekday", s)
  s <- gsub("fasting", "Fasting", s)
  s <- gsub("mouse_id", "Mouse", s)
  s <- gsub("year1", "Year1", s)
  s <- gsub("year2", "Year2", s)
  s <- gsub("year3", "Year3", s)

  return(s)
}


# Bad IDs are either split peaks or in source fragments with bad identifications
docr_bad_duplicate_features <- function() {
  bad_duplicate_feature_ids <- c(
    "1-Aminocyclopropanecarboxylic acid.M012A.224", "1-Aminocyclopropanecarboxylic acid.M012A.223",
    "3-Amino-2-piperidone.M012A.333", "3-Amino-2-piperidone.M012A.348", "4-Acetamidobutanoic acid.M012A.864",
    "5-Valerolactone.M012A.187", "Glycerol 3-phosphate.M012A.1382",
    "Heptanoylcarnitine.M012A.3085", "Hexanoylcarnitine.M012A.2900",
    "LPC(16:0/0:0).M012A.5463", "LPC(18:1).M012A.5539", "LPE(18:1).M012A.5096",
    "LPE(18:1).M012A.5097", "Octanoylcarnitine.M012A.3268",
    "Urocanate.M012A.763", "Urocanate.M012A.765", "Valerylcarnitine.M012A.2706",
    "2-Ethyl-2-hydroxybutyric acid.M013A.1788", "Hypoxanthine.M013A.2117",
    "3,4-Dimethylbenzoic acid.M013A.2619", "Aminoadipic acid.M013A.3315",
    "Cholic acid.M013A.12262", "Dehydroascorbic acid.M013A.3968",
    "Dihydrouracil.M013A.1052", "Dihydrouracil.M012A.331",
    "Galactonic acid.M013A.5205", "Hist-ser.M013A.7219",
    "N-acetylphenylalanine.M013A.5607", "Ophthalmic acid.M013A.9093",
    "Propenoic acid.M013A.183", "Propenoic acid.M013A.184",
    "Taurocholic acid.M012A.5155.5156.5157", "3-Hydroxymethylglutarate.M013A.3420",
    "Aminoadipic acid.M012A.1195", "5-Hydroxyindoleacetic acid.M013A.4838",
    "5-Hydroxyindoleacetic acid.M012A.1757", "Indoleacrylic acid.M013A.4513",
    "Allothreonine.M012A.413", "Itaconic acid.M013A.1892", "Itaconic acid.M013A.1895",
    "Beta-alanine.M012A.115", "Indoleacetaldehyde.M012A.1139",
    "DEHYDROASCORBIC ACID.M013A.3964",
    "SN-GLYCERO-3-PHOSPHOCHOLINE.M013A.7940",
    "N-METHYL-L-GLUTAMATE+20.0:10.0:180.0:80.0:70.0:90.0:40.0:30.0:160.0:200.0:120.0:140.0:100.0:60.0:50.0.M012A.1194",
    "Cis-4-hydroxy-proline.M012A.588"
  )
  return(bad_duplicate_feature_ids)
}


# DO-CR diet colors
docr_get_diet_colors <- function() {
  diet_colors <- c("seashell4", "skyblue", "royalblue4", "orange", "firebrick")
  names(diet_colors) <- c("AL", "1D", "2D", "20", "40")
  return(diet_colors)
}

# DO-CR age colors
docr_get_age_colors <- function() {
  age_colors <- c("#DCA8FB", "#C05CFA", "#5E2D7A")
  names(age_colors) <- c("AgeYear1", "AgeYear2", "AgeYear3")
  return(age_colors)
}

# DOCR "pretty" diet terms
docr_get_diet_terms <- function() {
  diet_terms <- c("Diet1D", "Diet2D", "Diet20", "Diet40")
  return(diet_terms)
}

# DOCR "pretty" age terms
docr_get_age_terms <- function() {
  age_terms <- c("AgeYear2", "AgeYear3")
  return(age_terms)
}

# Get high and low colors for heatmaps and correlation plots
docr_get_high_low_colors <- function(center_white = TRUE) {
  if (center_white) {
    return(c(Low = "#21336B", Mid = "white", High = "#D63B24"))
  } else {
    return(c(Low = "#21336B", High = "#D63B24"))
  }
}

# ggplot theme for volcano plots or multi-facet plots
docr_ggplot_theme <- function(x_text_angle = 45,
                              x_text_hjust = 1,
                              legend_position = "none") {
  theme_return <- ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position,
      axis.text.x = ggplot2::element_text(
        angle = x_text_angle,
        hjust = x_text_hjust
      ),
      plot.title = element_text(size = 7),
      plot.subtitle = element_text(size = 6),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 6, hjust = 0)
    )
  return(theme_return)
}


# ggplot line for volcano plots and other scatter plots
docr_ggplot_line <- function(x = NULL,
                             y = NULL) {
  if (!is.null(x)) {
    line_return <- ggplot2::geom_vline(
      xintercept = x,
      linetype = "dashed",
      color = "grey50"
    )
  }
  if (!is.null(y)) {
    line_return <- ggplot2::geom_hline(
      yintercept = y,
      linetype = "dashed",
      color = "grey50"
    )
  }
  return(line_return)
}


# Easily save docr figures
docr_ggsave <- function(plot_object,
                        plot_file_name,
                        plot_height,
                        plot_width,
                        local_filepath = local_filepath,
                        figure_folder = "final_figures",
                        device_use = "pdf") {
  ggplot_filname <- file.path(local_filepath, figure_folder, paste0(plot_file_name, ".", device_use))

  ggplot2::ggsave(
    filename = ggplot_filname,
    plot = plot_object,
    height = plot_height,
    width = plot_width,
    units = "in",
    dpi = 450,
    device = device_use
  )
}
