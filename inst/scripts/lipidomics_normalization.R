###
### Data cleaning and normalization for DO-CR lipidomics
### Johanna Fleischman, Calico Life Sciences 2025
###
###
### Install libraries and source functions
suppressMessages(library(tidyverse))
library(claman)
library(romic)

# set
# local_filepath <- "~/drido-multiomic-paper"
# output_filepath <- "~/"

source(file.path(local_filepath, "R/normalization_functions.R"))
source(file.path(local_filepath, "R/dims_normalization.R"))

### Print plot toggle for viewing in IDE
### Defaults to FALSE for running script
### print_plot == FALSE stops correlation and PCA plots from printing
### print_plot == FALSE does NOT stop PCs from being generated to determine outliers
print_plot <- TRUE

print("Libraries and functions loaded")

### Load files
lipidomics <- readRDS(file.path(local_filepath, "inst/extdata/20250117-PreNorm-Lipidomics.Rds"))
metadata <- readRDS(file.path(local_filepath, "inst/extdata/DOCR_Sample_Metadata.Rds"))

print("Lipidomics data and metadata loaded")

### Extract/clean sample data and combine with metadata
lipidomics_tidy_omic_temp <- lipidomics %>%
  dplyr::mutate(
    log2_intensity_median = log2(intensity_median),
    sampleId = as.character(as.numeric(as.factor(name_unique))),
    groupId = compoundName
  )

lipid_positive_temp <- lipidomics_tidy_omic_temp %>%
  dplyr::filter(method == "M014B")

lipid_negative_temp <- lipidomics_tidy_omic_temp %>%
  dplyr::filter(method == "M015B")

### Convert to triple omic object for use with metabolomics functions
lipid_positive <- romic::create_triple_omic(
  measurement_df = lipid_positive_temp %>%
    dplyr::select(
      groupId, sampleId, centered_IS_norm_intensity, log2_centered_IS_norm_intensity,
      is_lowest_quant_class_CV, plate_BulkPool_CV, plate_biological_sample_CV,
      median_BulkPool_CV, intensity_median, log2_intensity_median, quant_type
    ) %>%
    dplyr::group_by(groupId, sampleId) %>%
    dplyr::slice_head(n = 1) %>% # Removes 2nd positive or negative control sample
    dplyr::ungroup() %>%
    dplyr::distinct(),
  feature_df = lipid_positive_temp %>%
    dplyr::select(
      groupId, lipidClass, compositionSummary, compoundName, adductName, quant_class,
      method, mode
    ) %>%
    dplyr::distinct(),
  sample_df = lipid_positive_temp %>%
    dplyr::select(
      sampleId, name_unique, sample_name, diet, mouse_id, timepoint,
      well, inj_order, run_num, plate, type, barcode, well_position,
      sample_type, order_num, injection_num
    ) %>%
    dplyr::distinct(),
  feature_pk = "groupId",
  sample_pk = "sampleId"
)

lipid_negative <- romic::create_triple_omic(
  measurement_df = lipid_negative_temp %>%
    dplyr::select(
      groupId, sampleId, centered_IS_norm_intensity, log2_centered_IS_norm_intensity,
      is_lowest_quant_class_CV, plate_BulkPool_CV, plate_biological_sample_CV,
      median_BulkPool_CV, intensity_median, log2_intensity_median, quant_type
    ) %>%
    dplyr::distinct(),
  feature_df = lipid_negative_temp %>%
    dplyr::select(
      groupId, lipidClass, compositionSummary, compoundName, adductName, quant_class,
      method, mode
    ) %>%
    dplyr::distinct(),
  sample_df = lipid_negative_temp %>%
    dplyr::select(
      sampleId, name_unique, sample_name, diet, mouse_id, timepoint,
      well, inj_order, run_num, plate, type, barcode, well_position,
      sample_type, order_num, injection_num
    ) %>%
    dplyr::distinct(),
  feature_pk = "groupId",
  sample_pk = "sampleId"
)

invisible(gc())
rm(list = ls()[grep("_temp", ls())])

### QC Visualization -----

### Expand data, impute missing values, remove high missingness values
# Imputes around LOD half minimum, so set seed for replicability
lipid_positive_temp <- docr_impute_missingness(lipid_positive,
  quant_var = "log2_centered_IS_norm_intensity",
  percent_missingness_threshold = 0.6,
  sample_type_filter = "biological_sample",
  sample_type_var = "type",
  seed = 321
)
lipid_negative_temp <- docr_impute_missingness(lipid_negative,
  quant_var = "log2_centered_IS_norm_intensity",
  percent_missingness_threshold = 0.6,
  sample_type_filter = "biological_sample",
  sample_type_var = "type",
  seed = 321
)

invisible(gc())

### Data cleaning for lipid data
lipid_positive_temp <- docr_clean_lipid_data(lipid_positive_temp) %>%
  dplyr::left_join(metadata, by = c("mouse_id", "age_years"))

lipid_negative_temp <- docr_clean_lipid_data(lipid_negative_temp) %>%
  dplyr::left_join(metadata, by = c("mouse_id", "age_years"))

### Plot in PCA space
for (ion_mode in c("Positive mode, ", "Negative mode, ")) {
  if (ion_mode == "Positive mode, ") {
    pca_use <- lipid_positive_temp
  } else {
    pca_use <- lipid_negative_temp
  }

  for (vcol in c("log2_abundance_halfmin")) {
    plot_title_temp <- paste0(ion_mode, vcol)

    sample_metadata_cols <- c(
      "age_years", "diet", "weekday_collection",
      "surv_days", "sample_type", "generation_wave", "plate",
      "inj_order", "well", "num_NAs", "sum_intensity"
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
        plot_title = plot_title_temp
      )
    }
  }
}


### Sample-only data generation -----

###' For lipidomics, outliers were handled on a mode-wise basis, due to
###' different missingness distributions. In positive mode, samples with
###' > 80% missingness are removed. In negative mode, samples with > 50%
###' missingness are removed. Samples with sum_intensity > 500 are removed.
###'
###' Outliers significant affect missingness imputation, which is why they
###' are removed at this point in the analysis

lipid_positive <- docr_remove_lipid_outliers(lipid_positive,
  per_NA_cutoff = 0.8
)
lipid_negative <- docr_remove_lipid_outliers(lipid_negative,
  per_NA_cutoff = 0.5
)

### Expand data, impute missing values, remove high missingness values
# Imputes around LOD half minimum of SAMPLES ONLY
lipid_positive_samples <- romic::filter_tomic(lipid_positive,
  filter_type = "category",
  filter_table = "samples",
  filter_variable = "type",
  filter_value = "biological_sample"
)
lipid_positive_samples <- docr_impute_missingness(lipid_positive_samples,
  quant_var = "log2_centered_IS_norm_intensity",
  percent_missingness_threshold = 0.6,
  sample_type_filter = "biological_sample",
  sample_type_var = "type",
  seed = 321
)

lipid_negative_samples <- romic::filter_tomic(lipid_negative,
  filter_type = "category",
  filter_table = "samples",
  filter_variable = "type",
  filter_value = "biological_sample"
)
lipid_negative_samples <- docr_impute_missingness(lipid_negative_samples,
  quant_var = "log2_centered_IS_norm_intensity",
  percent_missingness_threshold = 0.6,
  sample_type_filter = "biological_sample",
  sample_type_var = "type",
  seed = 321
)

invisible(gc())

### Data cleaning for lipid data
lipid_positive_samples <- docr_clean_lipid_data(lipid_positive_samples) %>%
  dplyr::left_join(metadata, by = c("mouse_id", "age_years"))

lipid_negative_samples <- docr_clean_lipid_data(lipid_negative_samples) %>%
  dplyr::left_join(metadata, by = c("mouse_id", "age_years"))

###' With removal of especially high and low intensity samples, the
###' Probabilistic Quotient Normalization (PQN) for global sample-to-sample
###' intensity differences is not performed on DIMS lipidomics data. PQN
###' normalization assumes chromatographic peak areas, which are not obtained
###' for direct infusion data; instead, total ion intensities would need to be
###' used. As total ion intensities can be affect by more than dilution
###' differences, and the lipidomics data is already normalized to plate-wise,
###' lipid-class internal standards and scaled to the plate-wise specific
###' compound median, no sample-to-sample intensity correction is performed.

### Well-drift correction
#### Positive mode
lipid_positive_samples <- lipid_positive_samples %>%
  tidyr::nest(data = -groupId)

lipid_positive_samples$data <- lapply(
  lipid_positive_samples$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "log2_abundance_halfmin",
      lm_norm_varname = "well",
      new_varname = "log2_abundance_halfmin_well"
    )
  }
)

lipid_positive_samples <- lipid_positive_samples %>%
  tidyr::unnest(data) %>%
  dplyr::ungroup()


#### Negative mode
lipid_negative_samples <- lipid_negative_samples %>%
  tidyr::nest(data = -groupId)

lipid_negative_samples$data <- lapply(
  lipid_negative_samples$data,
  function(x) {
    docr_linear_drift_normalize(
      mzroll_tidy_omic = x,
      quant_peak_varname = "log2_abundance_halfmin",
      lm_norm_varname = "well",
      new_varname = "log2_abundance_halfmin_well"
    )
  }
)

lipid_negative_samples <- lipid_negative_samples %>%
  tidyr::unnest(data) %>%
  dplyr::ungroup()

invisible(gc())


### PCA on all versions of normalized sample data
if (print_plot) {
  for (ion_mode in c("Positive mode, ", "Negative mode, ")) {
    if (ion_mode == "Positive mode, ") {
      pca_use <- lipid_positive_samples
    } else {
      pca_use <- lipid_negative_samples
    }

    for (vcol in c(
      "log2_abundance_halfmin",
      "log2_abundance_halfmin_well"
    )) {
      plot_title_temp <- paste0(ion_mode, vcol)

      sample_metadata_cols <- c(
        "age_years", "diet", "weekday_collection",
        "surv_days", "generation_wave", "plate",
        "inj_order", "well", "num_NAs", "sum_intensity"
      )

      m <- docr_extract_pca_data(pca_use,
        value_col = vcol,
        sample_metadata_cols = sample_metadata_cols
      )

      docr_plot_pca_data(m$pca_df,
        m$pca_var,
        color_vars = c("age_years"),
        pcs = c("PC1", "PC2"),
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
lipids_all <- lipid_positive_samples %>%
  dplyr::bind_rows(lipid_negative_samples) %>%
  dplyr::mutate(
    norm_abundance = log2_abundance_halfmin_well,
    feature_name = compoundName,
    feature_id = paste0(feature_name, ".", method)
  ) %>%
  docr_clean_metadata() %>%
  dplyr::mutate(is_imputed = ifelse(is.na(log2_centered_IS_norm_intensity),
    TRUE,
    FALSE
  )) %>%
  dplyr::select(
    feature_id, compoundName, sample_id, mouse_id, age_years, diet,
    diet_assignment, method, date_birth, date_collection, date_exit,
    delta_date, PLL, BW_Raw, BW_Loess, is_ddm, weekday_collection,
    generation_wave, surv_days, surv_years, BW_Raw, BW_Loess,
    is_imputed, norm_abundance
  )

### Save file
saveRDS(object = lipids_all, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Normalized-Lipidomics-Data.Rds")))

print("Files saved")
