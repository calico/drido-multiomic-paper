###
### Data cleaning and normalization for DO-CR metabolomics projects
### Johanna Fleischman, Calico Life Sciences 2024
###
###
### Install libraries and source functions
library(tidyverse)
library(claman)
library(romic)

# set
# local_filepath <- "~/drido-multiomic-paper"

source(file.path(local_filepath, "R/normalization_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

### Print plot toggle for viewing in IDE
### Defaults to FALSE for running script
### print_plot == FALSE stops correlation and PCA plots from printing
### print_plot == FALSE does NOT stop PCs from being generated to determine outliers
print_plot <- TRUE

print("Libraries and functions loaded")

### SQLite database files must be accessed from metabolomics workbench
qe_positive_file <- "DOCR_QE_Positive.mzrollDB"
qe_negative_file <- "DOCR_QE_Negative.mzrollDB"

# ## Process mzrollDB file ----
# ### Positive mode ----
mzroll_positive <- claman::process_mzroll(
  mzroll_db_path = file.path(local_filepath, qe_positive_file),
  only_identified = FALSE,
  validate = FALSE
)
#
# ### Negative mode ----
mzroll_negative <- claman::process_mzroll(
  mzroll_db_path = file.path(local_filepath, qe_negative_file),
  only_identified = FALSE,
  validate = FALSE
)

### Load pre-processed files
# mzroll_positive <- readRDS(file.path(local_filepath, "mzrollDB_full_positive_processed.Rds"))
# mzroll_negative <- readRDS(file.path(local_filepath, "mzrollDB_full_negative_processed.Rds"))

print("mzrollDBs processed")

### Clean files to extract metadata
mzroll_positive <- docr_clean_mzroll_samples(mzroll_positive)
mzroll_negative <- docr_clean_mzroll_samples(mzroll_negative)

invisible(gc())

print("mzrollDBs cleaned")

### Check discard original samples for which there have been re-runs
mzroll_positive <- docr_remove_reruns(mzroll_positive)
mzroll_negative <- docr_remove_reruns(mzroll_negative)

### Only keep "good" (any "g") annotated features
mzroll_positive <- romic::filter_tomic(mzroll_positive,
  filter_type = "quo",
  filter_table = "features",
  filter_value = rlang::quo(grepl("g", label))
)

mzroll_negative <- romic::filter_tomic(mzroll_negative,
  filter_type = "quo",
  filter_table = "features",
  filter_value = rlang::quo(grepl("g", label))
)

### Remove manually identified fragment ions and combine manually identified split peaks
positive_fragments_remove <- c("835", "1743", "581", "994", "692", "1386", "2374")
positive_peaks_combine <- list(
  c("5521", "5522"),
  c("1350", "1352"),
  c("1817", "1818"),
  c("329", "330"),
  c("5189", "5190", "5191"),
  c("1960", "1961"),
  c("5155", "5156", "5157"),
  c("5279", "5280", "5281"),
  c("812", "813")
)

negative_fragments_remove <- c("799", "1165", "2661", "3743", "6562")
negative_peaks_combine <- list(
  c("11394", "11395"),
  c("8378", "8379"),
  c("8851", "8852"),
  c("3909", "3910"),
  c("12312", "12313")
)

mzroll_positive <- docr_remove_redundant_fragments(mzroll_positive,
  fragments_to_remove = positive_fragments_remove,
  split_peaks = positive_peaks_combine
)
mzroll_negative <- docr_remove_redundant_fragments(mzroll_negative,
  fragments_to_remove = negative_fragments_remove,
  split_peaks = negative_peaks_combine
)


### Expand data, impute missing values, remove high missingness values
# Imputes around LOD half minimum, so set seed for replicability
mzroll_positive <- docr_impute_missingness(mzroll_positive,
  seed = 321
)
mzroll_negative <- docr_impute_missingness(mzroll_negative,
  seed = 321
)

invisible(gc())

print("Missing values imputed")

### Drift correct and center each feature
# Mass spec drift is correct by subtracting the predicted linear estimate on injection
# order for all samples from each feature value
mzroll_positive_norm <- mzroll_positive %>%
  romic::tomic_to("tidy_omic") %>%
  purrr::pluck(1) %>%
  tidyr::nest(data = -c(plate, groupId))

mzroll_positive_norm$data <- lapply(
  mzroll_positive_norm$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "log2_abundance",
      lm_norm_varname = "inj_order",
      sample_type_string = "Sample",
      new_varname = "samples_linear_NA"
    )
  }
)

mzroll_positive_norm$data <- lapply(
  mzroll_positive_norm$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "log2_abundance_halfmin",
      lm_norm_varname = "inj_order",
      sample_type_string = "Sample",
      new_varname = "samples_linear_halfmin"
    )
  }
)
invisible(gc())

mzroll_negative_norm <- mzroll_negative %>%
  romic::tomic_to("tidy_omic") %>%
  purrr::pluck(1) %>%
  tidyr::nest(data = -c(plate, groupId))

mzroll_negative_norm$data <- lapply(
  mzroll_negative_norm$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "log2_abundance",
      lm_norm_varname = "inj_order",
      sample_type_string = "Sample",
      new_varname = "samples_linear_NA"
    )
  }
)

mzroll_negative_norm$data <- lapply(
  mzroll_negative_norm$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "log2_abundance_halfmin",
      lm_norm_varname = "inj_order",
      sample_type_string = "Sample",
      new_varname = "samples_linear_halfmin"
    )
  }
)
invisible(gc())


print("Drift correction complete")


### Add animal biological metadata and create technical metadata
metadata <- readRDS(file.path(local_filepath, "DOCR_Sample_Metadata.Rds"))

mzroll_positive_norm <- mzroll_positive_norm %>%
  tidyr::unnest(data) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., metadata, by = c("mouse_id", "age_years")) %>%
  dplyr::group_by(sampleId) %>%
  dplyr::mutate(
    num_NAs = sum(!is.finite(log2_abundance)),
    sum_intensity = sum(log2_abundance, na.rm = T)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sampleId = as.character(sampleId),
    groupId = as.character(groupId)
  )

mzroll_negative_norm <- mzroll_negative_norm %>%
  tidyr::unnest(data) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., metadata, by = c("mouse_id", "age_years")) %>%
  dplyr::group_by(sampleId) %>%
  dplyr::mutate(
    num_NAs = sum(!is.finite(log2_abundance)),
    sum_intensity = sum(log2_abundance, na.rm = T)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sampleId = as.character(sampleId),
    groupId = as.character(groupId)
  )

### Plot in PCA space
for (ion_mode in c("Positive mode, ", "Negative mode, ")) {
  if (ion_mode == "Positive mode, ") {
    pca_use <- mzroll_positive_norm
  } else {
    pca_use <- mzroll_negative_norm
  }

  for (vcol in c(
    "log2_abundance_halfmin",
    "samples_linear_halfmin"
  )) {
    plot_title_temp <- paste0(ion_mode, vcol)

    sample_metadata_cols <- c(
      "age_years", "diet", "weekday_collection", "surv_days", "sample_type",
      "generation_wave", "plate", "inj_order", "well", "num_NAs", "sum_intensity"
    )

    if (print_plot) {
      ## All sample types
      m <- docr_extract_pca_data(pca_use,
        value_col = vcol,
        sample_metadata_cols = sample_metadata_cols
      )

      docr_plot_pca_data(m$pca_df,
        m$pca_var,
        color_vars = c("sample_type"),
        plot_title = vcol
      )

      docr_pca_correlation_plot(m$pca_df,
        m$pca_var,
        sample_metadata_cols = sample_metadata_cols,
        plot_title = plot_title_temp
      )
    }

    # Biological samples only
    m <- docr_extract_pca_data(
      pca_use %>%
        dplyr::filter(sample_type %in% "Sample"),
      value_col = vcol,
      sample_metadata_cols = c(sample_metadata_cols, "sample_name")
    )

    if (print_plot) {
      docr_plot_pca_data(m$pca_df,
        m$pca_var,
        color_vars = c("plate", "inj_order", "age_years", "diet"),
        plot_title = paste0(plot_title_temp, ", samples only")
      )

      docr_pca_correlation_plot(m$pca_df,
        m$pca_var,
        sample_metadata_cols = sample_metadata_cols,
        plot_title = paste0(plot_title_temp, ", samples only")
      )
    }

    if (vcol == "samples_linear_halfmin" && ion_mode == "Positive mode, ") {
      positive_outliers <- claman::check_outliers(m$pca_df)
      print(paste0(
        length(positive_outliers), " outliers identified in positive mode: ",
        paste(positive_outliers, collapse = ", ")
      ))
    } else if (vcol == "samples_linear_halfmin" && ion_mode == "Negative mode, ") {
      negative_outliers <- claman::check_outliers(m$pca_df)
      print(paste0(
        length(negative_outliers), " outliers identified in negative mode: ",
        paste(negative_outliers, collapse = ", ")
      ))
    } else {
      # Do nothing
    }
  }
}

mzroll_positive_norm <- mzroll_positive_norm %>%
  dplyr::filter(!sample_name %in% positive_outliers)

mzroll_negative_norm <- mzroll_negative_norm %>%
  dplyr::filter(!sample_name %in% negative_outliers)

print("Outliers removed")

### Sample-to-sample intensity normalization
mzroll_positive_samples <- mzroll_positive_norm %>%
  dplyr::filter(sample_type %in% "Sample") %>%
  docr_intensity_normalize(
    quant_var = "samples_linear_NA",
    med_var = "log2_abundance",
    norm_var = "samples_linear_NA_pqn"
  ) %>%
  docr_intensity_normalize(
    quant_var = "samples_linear_halfmin",
    med_var = "log2_abundance_halfmin",
    norm_var = "samples_linear_halfmin_pqn"
  )

mzroll_negative_samples <- mzroll_negative_norm %>%
  dplyr::filter(sample_type %in% "Sample") %>%
  docr_intensity_normalize(
    quant_var = "samples_linear_NA",
    med_var = "log2_abundance",
    norm_var = "samples_linear_NA_pqn"
  ) %>%
  docr_intensity_normalize(
    quant_var = "samples_linear_halfmin",
    med_var = "log2_abundance_halfmin",
    norm_var = "samples_linear_halfmin_pqn"
  )

invisible(gc())


### Well-drift correction
#### Positive mode
mzroll_positive_samples <- mzroll_positive_samples %>%
  tidyr::nest(data = -groupId)

mzroll_positive_samples$data <- lapply(
  mzroll_positive_samples$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "samples_linear_NA_pqn",
      lm_norm_varname = "well",
      new_varname = "samples_linear_NA_pqn_well"
    )
  }
)

mzroll_positive_samples$data <- lapply(
  mzroll_positive_samples$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "samples_linear_halfmin_pqn",
      lm_norm_varname = "well",
      new_varname = "samples_linear_halfmin_pqn_well"
    )
  }
)

mzroll_positive_samples <- mzroll_positive_samples %>%
  tidyr::unnest(data) %>%
  dplyr::ungroup()


#### Negative mode
mzroll_negative_samples <- mzroll_negative_samples %>%
  tidyr::nest(data = -groupId)

mzroll_negative_samples$data <- lapply(
  mzroll_negative_samples$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "samples_linear_NA_pqn",
      lm_norm_varname = "well",
      new_varname = "samples_linear_NA_pqn_well"
    )
  }
)

mzroll_negative_samples$data <- lapply(
  mzroll_negative_samples$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "samples_linear_halfmin_pqn",
      lm_norm_varname = "well",
      new_varname = "samples_linear_halfmin_pqn_well"
    )
  }
)

mzroll_negative_samples <- mzroll_negative_samples %>%
  tidyr::unnest(data) %>%
  dplyr::ungroup()

invisible(gc())


### PCA on all versions of normalized sample data
if (print_plot) {
  for (ion_mode in c("Positive mode, ", "Negative mode, ")) {
    if (ion_mode == "Positive mode, ") {
      pca_use <- mzroll_positive_samples
    } else {
      pca_use <- mzroll_negative_samples
    }

    for (vcol in c(
      "samples_linear_halfmin",
      "samples_linear_halfmin_pqn",
      "samples_linear_halfmin_pqn_well"
    )) {
      plot_title_temp <- paste0(ion_mode, vcol)

      sample_metadata_cols <- c(
        "age_years", "diet", "weekday_collection", "surv_days", "generation_wave",
        "plate", "inj_order", "well", "num_NAs", "sum_intensity"
      )

      m <- docr_extract_pca_data(pca_use,
        value_col = vcol,
        sample_metadata_cols = sample_metadata_cols
      )

      docr_plot_pca_data(m$pca_df,
        m$pca_var,
        color_vars = c("age_years"),
        plot_title = plot_title_temp
      )

      docr_pca_correlation_plot(m$pca_df,
        m$pca_var,
        sample_metadata_cols = sample_metadata_cols,
        plot_title = plot_title_temp
      )
    }
  }
}

### Combine methods
mzroll_all <- mzroll_positive_samples %>%
  dplyr::mutate(method = "M012A") %>%
  dplyr::bind_rows(mzroll_negative_samples %>%
    dplyr::mutate(method = "M013A")) %>%
  dplyr::mutate(
    norm_abundance = samples_linear_NA_pqn_well,
    feature_name = compoundName,
    feature_id = paste0(feature_name, ".", method, ".", groupId)
  ) %>%
  docr_clean_metadata() %>%
  dplyr::select(
    feature_id, compoundName, mz, rt, sample_id, mouse_id, age_years, diet, diet_assignment, method,
    date_birth, date_collection, date_exit, delta_date, PLL, BW_Raw, BW_Loess, is_ddm,
    weekday_collection, generation_wave, surv_days, surv_years, BW_Raw, BW_Loess, norm_abundance
  )

### Save file
saveRDS(object = mzroll_all, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Normalized-Metabolomics-Data.Rds")))

print("Files saved")
