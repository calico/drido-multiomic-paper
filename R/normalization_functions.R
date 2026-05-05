#' Clean mzrollDB samples dataframe and assign injection orders
#'
#' @description
#' This is function is used to clean the sample information in the DO-CR mzrollDB files
#'
#' @param mzroll an mzrollDB file
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns original mzrollDB file with updated \code{samples} object
#'
#' @export
docr_clean_mzroll_samples <- function(mzroll) {
  # Clean sample_name, correctly factorize other metadata
  suppressWarnings({
    mzroll$samples <- mzroll$samples %>%
      tidyr::separate(name,
        into = c(NA, NA, "sample_name", "plate", "position", "inj_order", NA, NA),
        sep = "_|\\.", remove = FALSE
      ) %>%
      dplyr::mutate(sample_type = dplyr::case_when(sample_name == "Blank" ~ "Blank",
        grepl("std", sample_name) ~ "Standard",
        grepl("NC", sample_name) ~ "NC",
        grepl("PC", sample_name) ~ "PC",
        sample_name == "BulkPool" ~ "BulkPool",
        sample_name == "BulkPool-DO" ~ "BulkPool",
        .default = "Sample"
      )) %>%
      dplyr::filter(plate != "p26A") %>% # Remove plate from separate experiment
      tidyr::separate(sample_name,
        into = c(NA, "diet", "mouse_id", "sample_age"),
        sep = "-", remove = FALSE
      ) %>%
      dplyr::mutate(
        inj_order = as.numeric(inj_order),
        sample_age = as.numeric(sample_age),
        diet = factor(diet, c("AL", "1D", "2D", "20", "40")),
        mouse_id = if_else(is.na(mouse_id), NA, paste("DO", diet, mouse_id, sep = "-")),
        age_years = as.factor(case_when(
          sample_age < 50 ~ "year1",
          sample_age > 50 & sample_age < 100 ~ "year2",
          sample_age > 100 ~ "year3"
        )),
        diet_assignment = diet,
        diet = ifelse(age_years == "year1", "AL", as.character(diet)),
        diet = factor(diet, levels = c("AL", "1D", "2D", "20", "40"))
      ) %>%
      dplyr::mutate(position = paste0(
        substr(position, 1, 1),
        stringr::str_pad(substr(position, 2, nchar(position)), 2, pad = "0", side = "left")
      )) %>%
      dplyr::mutate(well = as.numeric(as.factor(position))) %>%
      dplyr::mutate(plate_run = substr(plate, 4, 4)) %>%
      dplyr::arrange(inj_order) %>%
      dplyr::group_by(plate) %>%
      dplyr::mutate(inj_order = 1:n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(plate_order = paste(plate,
        stringr::str_pad(inj_order, 3, pad = "0", side = "left"),
        sep = "_"
      ))
  })

  # Update triple-omic object
  mzroll <- mzroll %>%
    romic::update_tomic(mzroll$samples)

  return(mzroll)
}


#' Clean lipid data for DO-CR
#'
#' @description
#' This is function is used to clean the sample information in the DO-CR
#' lipids dataframe
#'
#' @param mzroll_triple_omic a \code{romic} \code{triple_omic} object
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a cleaned \code{romic} \code{tidy_omic} object
#'
#' @export
docr_clean_lipid_data <- function(mzroll_triple_omic) {
  suppressWarnings({
    mzroll_triple_omic_return <- mzroll_triple_omic %>%
      romic::tomic_to("tidy_omic") %>%
      purrr::pluck(1) %>%
      tidyr::separate(sample_name,
        into = c(NA, "diet", "mouse_id", "sample_age"),
        sep = "-", remove = FALSE
      ) %>%
      dplyr::mutate(
        inj_order = as.numeric(inj_order),
        sample_age = as.numeric(sample_age),
        diet = factor(diet, c("AL", "1D", "2D", "20", "40")),
        mouse_id = if_else(is.na(mouse_id), NA, paste("DO", diet, mouse_id, sep = "-")),
        age_years = as.factor(case_when(
          sample_age < 50 ~ "year1",
          sample_age > 50 & sample_age < 100 ~ "year2",
          sample_age > 100 ~ "year3"
        )),
        diet_assignment = diet,
        diet = ifelse(age_years == "year1", "AL", as.character(diet)),
        diet = factor(diet, levels = c("AL", "1D", "2D", "20", "40"))
      ) %>%
      dplyr::select(
        groupId, lipidClass, compoundName, method, sampleId, mouse_id,
        age_years, diet, diet_assignment, well, inj_order, plate, type,
        sample_type, well_position, plate_biological_sample_CV,
        log2_centered_IS_norm_intensity,
        log2_abundance_halfmin,
        intensity_median,
        log2_intensity_median
      ) %>%
      dplyr::mutate(plate = as.factor(plate)) %>%
      dplyr::mutate(position = paste0(substr(well, 1, 1), str_pad(substr(well, 2, nchar(well)), 2, pad = "0", side = "left"))) %>%
      dplyr::mutate(well = as.numeric(as.factor(position))) %>%
      dplyr::mutate(inj_order = as.numeric(inj_order)) %>%
      dplyr::group_by(sampleId) %>%
      dplyr::mutate(
        num_NAs = sum(!is.finite(log2_centered_IS_norm_intensity)),
        per_NAs = num_NAs / dplyr::n(),
        sum_intensity = sum(log2_centered_IS_norm_intensity, na.rm = T)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        sampleId = as.character(sampleId),
        groupId = as.character(groupId)
      )
  })

  return(mzroll_triple_omic_return)
}


#' Remove fragment ions and combine split peaks
#'
#' @description
#' Remove manually identified fragment ions and combine manually identified split
#' peaks, and add technical metadata for each sample. Assumes sample columns as
#' \code{sampleId} and feature column as \code{groupId}. This function also calculates
#' the resulting \code{num_NAs} and \code{sum_intensity} meta value for each \code{sampleId}
#'
#' @param mzroll_triple_omic a \code{romic} \code{triple_omic} object
#' @param fragments_to_remove a character vector of \code{groupId} values to
#' remove from \code{mzroll_tidy_omic}
#' @param split_peaks a list of character vectors of \code{groupId} values to
#' combine
#' @param value_var column name as string of peak abundance value column to sum for
#' split peak groups
#'
#' @return \code{mzroll_tidy_omic} with fragments removed, split peaks removed, combined
#' peaks added, and technical metadata added
#'
#' @export
docr_remove_redundant_fragments <- function(mzroll_triple_omic,
                                            fragments_to_remove = c(),
                                            split_peaks = c(),
                                            value_var = "log2_abundance") {
  mzroll_tidy_omic <- mzroll_triple_omic %>%
    romic::tomic_to("tidy_omic")

  fixed_peaks_all <- data.frame()

  for (p in split_peaks) {
    fixed_peaks <- mzroll_tidy_omic %>%
      purrr::pluck(1) %>%
      dplyr::filter(groupId %in% p)

    if (length(unique(fixed_peaks$compoundName)) != 1) {
      print(paste0("For split peak group ", p, ", there is not exactly one unique compoundName; can't combine"))
      next
    }

    # Assign new groupId as groupId1.groupId2 etc
    new_groupId <- paste0(unique(fixed_peaks$groupId), collapse = ".")

    # Sum peak data
    fixed_peaks <- fixed_peaks %>%
      dplyr::group_by(sampleId) %>%
      dplyr::mutate(new_value_var = log2(sum(2^(!!rlang::sym(value_var)), na.rm = T))) %>%
      # We only want one entry per new groupId
      # NOTE that this extracts one data point for centered_log2_abundance
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      # Assign new feature_id and replace norm_abundance value
      dplyr::mutate(groupId = new_groupId) %>%
      dplyr::select(-!!rlang::sym(value_var)) %>%
      dplyr::rename(!!rlang::sym(value_var) := "new_value_var")

    fixed_peaks_all <- bind_rows(
      fixed_peaks_all,
      fixed_peaks
    )
  }

  ### Remove manually identified fragment ions and combine manually identified split peaks
  mzroll_tidy_omic_return <- mzroll_tidy_omic$data %>%
    dplyr::filter(
      !groupId %in% fragments_to_remove,
      (!groupId %in% unlist(split_peaks))
    ) %>%
    dplyr::bind_rows(fixed_peaks_all)

  ### Print messages to check that removal and additions were correct
  old_groupIds_removed <- paste0(setdiff(unique(mzroll_tidy_omic$data$groupId), unique(mzroll_tidy_omic_return$groupId)), collapse = ", ")
  new_groupIds_added <- paste0(setdiff(unique(mzroll_tidy_omic_return$groupId), unique(mzroll_tidy_omic$data$groupId)), collapse = ", ")

  print(paste0("groupIds removed: ", old_groupIds_removed))
  print(paste0("groupIds added: ", new_groupIds_added))

  ### Modify features dataframe to summarise peak group information
  mzroll_tidy_omic$data <- mzroll_tidy_omic_return
  mzroll_triple_omic_return <- mzroll_tidy_omic %>%
    romic::tomic_to("triple_omic")
  mzroll_triple_omic_return$features <- mzroll_triple_omic_return$features %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

  mzroll_triple_omic_return <- mzroll_triple_omic_return %>%
    romic::update_tomic(mzroll_triple_omic_return$features)

  return(mzroll_triple_omic_return)
}


#' Remove outliers from lipidomics data
#'
#' @description
#' Removes outliers based on % missingness and total summed intensity value
#'
#' @param mzroll a \code{romic} \code{triple_omic} object
#' @param sample_type_filter sample type or types on which to calculate missing values;
#' defaults to \code{"Sample"} (this prevents Blanks from being used to calculate missingness)
#' @param sample_type_var column name as string from which to filter \code{sample_type_filter};
#' defaults to \code{"sample_type"}
#' @param batch_var column name as string that define batches; defaults to \code{"plate"}
#' @param sum_intensitper_NA_cutoffy_cutoff upper bound cutoff value for per_NA_cutoff
#' @param sum_intensity_cutoff upper bound cutoff value for sum_intensity
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns mzroll_triple_omic a \code{romic} object with outliers removed
#'
#' @export
docr_remove_lipid_outliers <- function(mzroll,
                                       sample_type_filter = "Sample",
                                       sample_type_var = "sample_type",
                                       batch_var = "plate",
                                       per_NA_cutoff = 0.5,
                                       sum_intensity_cutoff = 500) {
  # Expand measurements so we can accurately count missing features
  expanded_measurements <- tidyr::expand_grid(
    groupId = mzroll$features$groupId,
    sampleId = mzroll$samples$sampleId
  ) %>%
    dplyr::full_join(mzroll$measurements,
      by = c("groupId", "sampleId")
    )

  # Identify outliers by percent missing or sum intensity cutoff
  outliers <- mzroll$samples %>%
    dplyr::left_join(expanded_measurements, ., by = "sampleId") %>%
    dplyr::group_by(sampleId) %>%
    dplyr::mutate(
      num_NAs = sum(!is.finite(log2_centered_IS_norm_intensity)),
      per_NAs = num_NAs / dplyr::n(),
      sum_intensity = sum(log2_centered_IS_norm_intensity, na.rm = T)
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!!rlang::sym(sample_type_var) == sample_type_filter) %>%
    dplyr::filter(per_NAs > per_NA_cutoff | sum_intensity > sum_intensity_cutoff) %>%
    dplyr::pull(sampleId)

  print(paste0(
    length(outliers),
    " outliers identified and removed: ",
    paste(outliers, collapse = ", ")
  ))

  # Filter out outliers with romic function
  mzroll_return <- romic::filter_tomic(mzroll,
    filter_type = "quo",
    filter_table = "samples",
    filter_value = rlang::quo(!sampleId %in% outliers)
  )

  return(mzroll_return)
}


#' Impute missing values and remove high missingness peaks
#'
#' @description
#' This function takes a \code{romic} \code{triple_omic} object and removes groupIds
#' that have a percent missing over the \code{percent_missing_threshold} value (default 60%)
#' in the provided sample type only (here, `sample_type == "Sample"`).
#'
#' This function imputes the missing values but taking the half-minimum of the groupId in each
#' `batch_var` (here, `plate`). This function assumes that \code{quant_var} is
#' in log2 space. If a groupId is completely missing on a plate, the half minimum
#' of the whole `mzroll` is used. A value is sampled around a distribution of the half-minimum.
#'
#' @param mzroll a \code{romic} \code{triple_omic} object
#' @param quant_var column on which to perform imputation and calculate missingness;
#' defaults to "log2_abundance" and assumes log2 space.
#' @param percent_missingness_threshold percent above which groupIds are removed;
#' defaults to \code{0.6} (60%)
#' @param sample_type_filter sample type or types on which to calculate missing values;
#' defaults to \code{"Sample"} (this prevents Blanks from being used to calculate missingness)
#' @param sample_type_var column name as string from which to filter \code{sample_type_filter};
#' defaults to \code{"sample_type"}
#' @param batch_var column name as string that define batches; defaults to \code{"plate"}
#' @param imp_sd standard deviation on which missing values are imputed; defaults to \code{0.15}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table :=
#'
#' @returns mzroll_tidy_omic a \code{romic} \code{tidy_omic} object with high missingness
#' groupIds removed and missing values imputed (dataframe is expanded)
#'
#' @export
docr_impute_missingness <- function(mzroll,
                                    quant_var = "log2_abundance",
                                    percent_missingness_threshold = 0.6,
                                    sample_type_filter = "Sample",
                                    sample_type_var = "sample_type",
                                    batch_var = "plate",
                                    imp_var = "log2_abundance_halfmin",
                                    imp_sd = 0.15,
                                    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  mzroll$measurements <- tidyr::expand_grid(
    groupId = mzroll$features$groupId,
    sampleId = mzroll$samples$sampleId
  ) %>%
    dplyr::full_join(mzroll$measurements,
      by = c("groupId", "sampleId")
    )

  # Find limit of detection as half-minimum values
  new_measurements <- mzroll$samples %>%
    dplyr::left_join(mzroll$measurements, ., by = "sampleId") %>%
    # Calculate the global half minimum and compute missingness / remove high missingness features
    dplyr::group_by(groupId) %>%
    dplyr::mutate(percent_missing = sum(!is.finite(!!rlang::sym(quant_var)) & !!rlang::sym(sample_type_var) %in% sample_type_filter) / sum(!!rlang::sym(sample_type_var) %in% sample_type_filter)) %>%
    dplyr::mutate(groupId_half_min = min(!!rlang::sym(quant_var), na.rm = TRUE) - 1) %>%
    dplyr::ungroup() %>%
    dplyr::filter(percent_missing < percent_missingness_threshold) %>%
    # For each plate, assign LOD to half-minimum on plate
    # If no value present on plate, assign LOD to half-minimum from experiment
    dplyr::group_by(groupId, !!rlang::sym(batch_var)) %>%
    dplyr::mutate(lod = ifelse(all(is.na(!!rlang::sym(quant_var))),
      groupId_half_min,
      min(!!rlang::sym(quant_var), na.rm = TRUE) - 1
    )) %>%
    dplyr::ungroup() %>%
    # Impute remaining NA values
    dplyr::rowwise() %>%
    dplyr::mutate(!!rlang::sym(imp_var) := ifelse(is.na(!!rlang::sym(quant_var)),
      stats::rnorm(1, mean = lod, sd = imp_sd),
      !!rlang::sym(quant_var)
    )) %>%
    dplyr::ungroup() %>%
    dplyr::select(groupId, sampleId, !!rlang::sym(imp_var))

  # Number of features removed
  remove_n <- length(unique(mzroll$measurements$groupId)) - length(unique(new_measurements$groupId))
  print(paste0(remove_n, " features removed for missingness above ", percent_missingness_threshold * 100, "%"))
  print("Old dimensions:")
  print(dim(mzroll$measurements))

  new_measurements <- new_measurements %>%
    dplyr::right_join(mzroll$measurements, ., by = c("groupId", "sampleId"))

  mzroll <- mzroll %>%
    romic::update_tomic(new_measurements)

  # Percent of dataset imputed
  print(paste0(
    round(sum(is.na(new_measurements[[quant_var]])) / nrow(new_measurements) * 100, 2),
    "% of remaining missing values are imputed"
  ))
  print("New dimensions:")
  print(dim(mzroll$measurements))

  return(mzroll)
}


#' Remove re-runs
#'
#' @description
#' This is function is used to remove samples that are re-run (issue with negative mode only).
#' This re-run function is assumes that the redundant re-runs on plate "A" are satisfied by plate "B"
#' (this is correct and based on a priori knowledge). Assumes that an individual samples is represented
#' by the \code{sampleId} column.
#'
#' @param mzroll_triple_omic a \code{romic} \code{triple_omic} object
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a \code{romic} \code{triple_omic} object with redundant sampleIds removed
#'
#' @export
docr_remove_reruns <- function(mzroll_triple_omic) {
  # Find redundant sample_names
  redundant_sample_names <- mzroll_triple_omic$samples %>%
    dplyr::filter(sample_type == "Sample") %>%
    dplyr::group_by(sampleId) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(sample_name)

  # Find sampleId corresponding to redundant sample_names for plate == "A"
  redundant_sampleIds <- mzroll_triple_omic$samples %>%
    dplyr::group_by(sampleId) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::filter(sample_name %in% redundant_sample_names) %>%
    dplyr::filter(plate_run == "A") %>%
    dplyr::pull(sampleId) %>%
    unique()

  print(paste0(length(redundant_sampleIds), " redundant sampleIds removed"))

  good_sampleIds <- setdiff(unique(mzroll_triple_omic$samples$sampleId), redundant_sampleIds)

  print(dim(mzroll_triple_omic$measurements))

  mzroll_triple_omic_return <- romic::filter_tomic(mzroll_triple_omic,
    filter_type = "category",
    filter_table = "samples",
    filter_variable = "sampleId",
    filter_value = good_sampleIds
  )

  print(dim(mzroll_triple_omic_return$measurements))

  return(mzroll_triple_omic_return)
}


#' Correct for batch-wise drift and center data
#'
#' @description
#' This function corrects for drift in each batch/plate. A linear model is fit on
#' all samples in a plate and the estimate on injection order is subtracted.
#' This removes the linear drift and centers the plate data for each feature on zero.
#' Designed to be run separately on each feature.
#'
#' @param mzroll_tidy_omic a \code{romic::tomic_to("tidy_omic")} data object
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns original \code{tidy_omic} object with a measurement column of drift
#' corrected and centered measurements
#'
#' @export
docr_linear_drift_normalize <- function(mzroll_tidy_omic,
                                        quant_peak_varname,
                                        lm_norm_varname,
                                        new_varname = "log2_abundance_norm",
                                        sample_type_string = NULL) {
  if (!is.null(sample_type_string)) {
    mzroll_tidy_omic_filtered <- mzroll_tidy_omic %>%
      dplyr::filter(sample_type == sample_type_string) %>%
      tidyr::drop_na(!!rlang::sym(quant_peak_varname))
  } else {
    mzroll_tidy_omic_filtered <- mzroll_tidy_omic %>%
      tidyr::drop_na(!!rlang::sym(quant_peak_varname))
  }

  tryCatch(
    {
      # If there are 2 or fewer valid values on the plate, we do not have enough information to batch correct
      # Remove data
      if (nrow(mzroll_tidy_omic_filtered) < 3) {
        mzroll_tidy_omic <- mzroll_tidy_omic %>%
          dplyr::mutate(!!rlang::sym(new_varname) := NA)

        # If there are 3 or 4 valid values, remove the median value of the sample filter
      } else if (nrow(mzroll_tidy_omic_filtered) < 5) {
        # Take median on sample values so that it is not biased by blanks and standards
        median_value <- median(mzroll_tidy_omic_filtered[[quant_peak_varname]], na.rm = TRUE)

        mzroll_tidy_omic <- mzroll_tidy_omic %>%
          dplyr::mutate(!!rlang::sym(new_varname) := !!rlang::sym(quant_peak_varname) - median_value)

        # If there are 5 or more valid values, fit an lm() line by lm_norm_varname and subtract predicted value
      } else {
        lfit <- lm(as.formula(paste0(quant_peak_varname, " ~ ", lm_norm_varname)),
          data = mzroll_tidy_omic_filtered
        )

        mzroll_tidy_omic[[new_varname]] <-
          mzroll_tidy_omic[[quant_peak_varname]] - predict(lfit, newdata = mzroll_tidy_omic)

        # Remove attributes on new column
        attributes(mzroll_tidy_omic[[new_varname]]) <- NULL
      }
    },
    error = function(cond) {
      message(cond)
    }
  ) # end tryCatch

  return(mzroll_tidy_omic)
}


#' Correct for sample-to-sample intensity differences
#'
#' @description
#' This function uses transforms dataframe and prepares values for
#' probabilistic quotient normalization methods via the \code{docr_pqn} function
#'
#' @param mzroll_tidy_omic a \code{romic} \code{tidy_omic} object
#' @param quant_var column name as string to normalize
#' @param med_var column name as string from which to extract median values to add back
#' @param norm_var column name as string to assign new PQN-normalized values
#' @param filter_on_features vector of features on which to filter calculations.
#' Defaults to NULL, which computes medians on all features.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns original \code{tidy_omic} object with additional column that is intensity-normalized
#'
#' @export
docr_intensity_normalize <- function(mzroll_tidy_omic,
                                     quant_var = "samples_linear_NA",
                                     med_var = "log2_abundance",
                                     norm_var = "samples_linear_NA_pqn",
                                     filter_on_features = NULL) {
  ## Prepare for Probabilistic Quotient Normalization
  mzroll_tidy_omic_prep <- mzroll_tidy_omic %>%
    dplyr::group_by(groupId) %>%
    # Calculate group medians
    dplyr::mutate(med_var_med = median(!!rlang::sym(med_var), na.rm = T)) %>%
    dplyr::ungroup() %>%
    # Add back medians and exponentiate values
    dplyr::mutate(quant_var_plus_med = !!rlang::sym(quant_var) + med_var_med) %>%
    dplyr::mutate(quant_var_plus_med_exp = 2^quant_var_plus_med)

  # Make matrices for ~f(x) docr_pqn() and perform calculation
  matrix_for_docr_pqn <- mzroll_tidy_omic_prep %>%
    dplyr::select(groupId, sampleId, quant_var_plus_med_exp) %>%
    tidyr::spread(key = "groupId", value = "quant_var_plus_med_exp") %>%
    column_to_rownames("sampleId") %>%
    as.matrix()

  # Remove median value to re-center and log2 transform
  use_features <- intersect(colnames(matrix_for_docr_pqn), filter_on_features)

  if (!is.null(filter_on_features) && length(use_features) > 0) {
    normalized_measurements <- docr_pqn(matrix_for_docr_pqn,
      filter_on_features = use_features
    )
  } else {
    normalized_measurements <- docr_pqn(matrix_for_docr_pqn)
  }

  normalized_measurements <- normalized_measurements %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sampleId") %>%
    tidyr::gather(2:ncol(.), key = "groupId", value = "quant_var_norm") %>%
    dplyr::mutate(quant_var_norm_log2 = log2(quant_var_norm)) %>%
    dplyr::left_join(mzroll_tidy_omic_prep, ., by = c("groupId", "sampleId")) %>%
    dplyr::mutate(quant_var_return = quant_var_norm_log2 - med_var_med) %>%
    dplyr::select(groupId, sampleId, !!rlang::sym(norm_var) := quant_var_return)

  return_df <- mzroll_tidy_omic %>%
    dplyr::left_join(normalized_measurements, by = c("groupId", "sampleId"))

  return(return_df)
}


#' Probabilistic Quotient Normalization
#'
#' @description
#' This function uses probabilistic quotient normalization methods to correct
#' for sample-to-sample intensity differences. This version was adapted and parsed
#' down for this paper from the more general \code{Rcpm::pqn} function. The original
#' paper can be found \href{https://pubs.acs.org/doi/epdf/10.1021/ac051632c}{here}
#'
#' @param feature_matrix a matrix with samples on rows and features on columns
#' @param filter_on_features vector of features on which to filter calculations.
#' Defaults to NULL, which computes medians on all features.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns original \code{feature_matrix} that is PQN normalized
#'
#' @export
docr_pqn <- function(feature_matrix,
                     filter_on_features = NULL) {
  feature_matrix_normed <- matrix(
    nrow = nrow(feature_matrix),
    ncol = ncol(feature_matrix)
  )
  colnames(feature_matrix_normed) <- colnames(feature_matrix)
  rownames(feature_matrix_normed) <- rownames(feature_matrix)

  if (is.null(filter_on_features)) {
    mX <- as.numeric(apply(feature_matrix, 2, function(x) median(x, na.rm = T)))
    for (i in 1:nrow(feature_matrix)) {
      feature_matrix_normed[i, ] <- as.numeric(feature_matrix[i, ] / median(as.numeric(feature_matrix[i, ] / mX), na.rm = T))
    }
  } else {
    filtered_feature_matrix <- feature_matrix[, (colnames(feature_matrix) %in% filter_on_features)]
    mX <- as.numeric(apply(filtered_feature_matrix, 2, function(x) median(x, na.rm = T)))
    for (i in 1:nrow(feature_matrix)) {
      feature_matrix_normed[i, ] <- as.numeric(feature_matrix[i, ] / median(as.numeric(filtered_feature_matrix[i, ] / mX), na.rm = T))
    }
  }

  return(feature_matrix_normed)
}


#' Clean DOCR metadata
#'
#' @description
#' This is function is used to clean the metadata for all omics data
#'
#' @param tidy_data a tidy dataframe with all DO-CR metadata columns
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns original \code{tidy_data} with columns cleans and refactored
#'
#' @export
docr_clean_metadata <- function(tidy_data) {
  ### Correct erroneous metadata
  tidy_data_return <- tidy_data %>%
    dplyr::mutate(
      mouse_id = case_when(
        mouse_id == "DO-2D-4077" & date_collection == "2019-01-28" ~ "DO-2D-4079",
        TRUE ~ mouse_id
      ),
      date_exit = case_when(
        mouse_id == "DO-40-2033" ~ as.Date("2017-10-17"),
        mouse_id == "DO-2D-4079" & date_exit == "2019-01-02" ~ as.Date("2019-12-23"),
        TRUE ~ date_exit
      ),
      surv_days = case_when(
        mouse_id == "DO-2D-4079" ~ 1189,
        TRUE ~ surv_days
      ),
      BW_Loess = case_when(
        mouse_id == "DO-2D-4079" & age_years == "year3" ~ 34.42,
        TRUE ~ BW_Loess
      ),
      BW_Raw = case_when(
        mouse_id == "DO-2D-4079" & age_years == "year3" ~ 41.12,
        TRUE ~ BW_Raw
      )
    ) %>%
    ### Clean remaining columns
    dplyr::mutate(
      sample_id = paste0(mouse_id, "-", age_years),
      mouse_id = factor(mouse_id),
      generation_wave = factor(generation_wave),
      weekday_collection = factor(weekday_collection, levels = c("TWT", "Monday")),
      age_years = factor(age_years, levels = c("year1", "year2", "year3")),
      diet = factor(diet, levels = c("AL", "1D", "2D", "20", "40")),
      diet_assignment = factor(diet_assignment, levels = c("AL", "1D", "2D", "20", "40")),
      date_exit = as.Date(as.character(date_exit), format = "%Y-%m-%d"),
      date_collection = as.Date(as.character(date_collection), format = "%Y-%m-%d"),
      date_birth = as.Date(as.character(date_birth), format = "%Y-%m-%d"),
      delta_date = date_exit - date_collection,
      surv_days = as.numeric(date_exit - date_birth),
      surv_years = surv_days / 365,
      is_ddm = ifelse(delta_date < 21, TRUE, FALSE),
      PLL = as.numeric(date_collection - date_birth) / as.numeric(date_exit - date_birth)
    )

  return(tidy_data_return)
}


#' Generate PCA data for visualization functions
#'
#' @param data_tbl dataframe containing, at minimum, a sample column, compound column, and
#' value column.
#' @param sample_col column name of samples, as a string. Default is `sampleId`.
#' @param compound_col column name of compounds, as a string. Default is `groupId`
#' @param value_col column name of values, as a string. Value must be numeric.
#' Default is `log2_abundance`
#' @param sample_metadata_cols column names of metadata in `data_tbl` to be appended to
#' PCA output data for ease of visualization downstream
#' @param scale_pca option to scale data within PCA function. Default is `FALSE`.
#' @param center_pca option to center data within PCA function. Default is `TRUE`.
#' @importFrom magrittr "%>%"
#' @importFrom rlang ".data"
#' @returns PCA data in a list, where `pca_df` is data, `pca_var` is the variance explained,
#' and `pca_loadings` is the rotation information
#' @export
docr_extract_pca_data <- function(data_tbl,
                                  sample_col = "sampleId",
                                  compound_col = "groupId",
                                  value_col = "log2_abundance",
                                  sample_metadata_cols = NULL,
                                  scale_pca = FALSE,
                                  center_pca = TRUE) {
  if (!sample_col %in% colnames(data_tbl)) {
    stop("sample_col must be a column in provided dataframe")
  }
  # require that the sample ID is a character since the conversion from
  # a PCA matrix will convert it to a character
  checkmate::assertCharacter(data_tbl[[sample_col]])
  if (!compound_col %in% colnames(data_tbl)) {
    stop("compound_col must be a column in provided dataframe")
  }
  if (!value_col %in% colnames(data_tbl)) {
    stop("value_col must be a column in provided dataframe")
  }
  checkmate::assertNumeric(data_tbl[[value_col]])

  # Create PCA dataframe
  pca_df <- data_tbl %>%
    dplyr::select(
      sid = rlang::sym(sample_col),
      gid = rlang::sym(compound_col),
      fid = rlang::sym(value_col)
    ) %>%
    stats::na.omit() %>%
    tidyr::spread(.data$gid, .data$fid) %>%
    tibble::column_to_rownames("sid")

  # Run PCA function
  pca_data <- stats::prcomp(
    x = pca_df,
    scale. = scale_pca,
    center = center_pca
  )

  # Extract PCA projections
  pca_datax <- pca_data$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = sample_col)

  # Extract explained variances
  var_explained <- summary(pca_data)$importance[2, ] %>%
    t() %>%
    as.data.frame()

  # Extract compound loadings
  pca_loadings <- pca_data$rotation %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = compound_col)

  # If provided, add metadata to PCA projection dataframe
  if (!is.null(sample_metadata_cols)) {
    # sample mdata columns should be present
    purrr::walk(
      sample_metadata_cols,
      ~ checkmate::assertChoice(.x, colnames(data_tbl))
    )

    sample_metadata_cols <- union(sample_col, sample_metadata_cols)
    sample_metadata_cols <- colnames(data_tbl)[colnames(data_tbl) %in% sample_metadata_cols]

    pca_return_data <- data_tbl %>%
      dplyr::group_by(!!rlang::sym(sample_col)) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::select(tidyselect::all_of(sample_metadata_cols)) %>%
      dplyr::left_join(pca_datax, by = sample_col)
  } else {
    pca_return_data <- pca_datax
  }

  # Return list object
  return_obj <- list(
    pca_df = pca_return_data,
    pca_var = var_explained,
    pca_loadings = pca_loadings
  )

  return(return_obj)
}


#' Create PCA scatter plot
#'
#' This function generates a PCA scatter plot.
#'
#' @param pca_df PCA projections as first list item
#' @param pca_var PCA variances as second list item
#' @param color_vars column names of metadata variables as columns in `pca_df`. For each column
#' name provided, a separate scatter plot will be generated. If none provided, will print an
#' unannotated scatter plot
#' @param add_ellipse when \code{TRUE}, ellipses are added to groups specified by
#' variables in the \code{color_var} columns
#' @param ellipse_level confidence interval to plot ellipses if \code{add_ellipse}
#' is \code{TRUE}
#' @param pcs names of PCs to be plotted, as strings. Default is `c("PC1","PC2")`
#' @param plot_title Title for plot, as string
#' @param print_plot default \code{TRUE}; change to \code{FALSE} to save ggplot
#' object without printing
#' @importFrom magrittr "%>%"
#' @returns A list of ggplot objects, for each variable indicated by `color_vars`
#' @export
docr_plot_pca_data <- function(pca_df,
                               pca_var,
                               color_vars = NULL,
                               add_ellipse = FALSE,
                               ellipse_level = 0.95,
                               pcs = c("PC1", "PC2"),
                               plot_title = "PCA Plot",
                               print_plot = TRUE) {
  if (!all(pcs %in% colnames(pca_var)) && all(pcs %in% colnames(pca_df))) {
    stop("pcs must be exactly 2 PC-identifying column names in pca_df and pca_var")
  }
  if (!is.null(color_vars)) {
    color_vars_use <- intersect(color_vars, colnames(pca_df))
    if (length(color_vars_use) == 0) {
      stop("Color variables provided but not present in pca_df data frame")
    }
  }

  return_list <- list()

  # make base ggplot theme
  base_theme <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(
      x = .data[[pcs[1]]],
      y = .data[[pcs[2]]]
    )
  ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste0(pcs[1], " (", round(pca_var[1, pcs[1]] * 100, 1), " %)"),
      y = paste0(pcs[2], " (", round(pca_var[1, pcs[2]] * 100, 1), " %)"),
      title = plot_title
    )

  # no colors
  if (is.null(color_vars)) {
    p <- base_theme + ggplot2::geom_point()
    return_list[[1]] <- p

    # with colors
  } else {
    for (color_var in color_vars_use) {
      p <- base_theme +
        ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]))

      # choose scale (continuous or categorical)
      if (is.numeric(pca_df[[color_var]])) {
        p <- p +
          ggplot2::scale_color_continuous(type = "viridis")
      } else {
        v <- rep(RColorBrewer::brewer.pal(8, "Dark2"),
          length.out = length(unique(pca_df[[color_var]]))
        )
        p <- p +
          ggplot2::scale_color_manual(values = v)
      }

      # add ellipse, if specified
      if (add_ellipse) {
        p <- p +
          ggplot2::stat_ellipse(
            mapping = ggplot2::aes(
              group = .data[[color_var]],
              color = .data[[color_var]]
            ),
            level = ellipse_level
          )
      }

      # save plot
      return_list[[color_var]] <- p
    }
  }

  # print plots
  if (print_plot) {
    lapply(return_list, print)
  }

  return(return_list)
}


#' Create PCA-metadata correlation plot
#'
#' This function generates a correlation plot between PC variables and metadata
#' variables.
#'
#' @param pca_df PCA projections as first list item
#' @param pca_var PCA variances as second list item
#' @param sample_col column name of samples, as a string. Default is `sampleId`.
#' @param sample_metadata_cols column names of metadata in `pca_df` on which to perform
#' correlation analyses against designated `pcs`
#' @param pcs names of PCs to be plotted, as strings. Default is `c("PC1","PC2","PC3","PC4","PC5")`
#' @param plot_title Title for plot, as string
#' @importFrom magrittr "%>%"
#' @param factor_encoding encoding for character or factor variables so they
#' can be converted to a numeric variable which can be correlated with PCs:
#' \describe{
#'   \item{numeric}{Convert characters to factors and then to integer values}
#'   \item{onehot}{If there are 3+ levels of a character or factor, 1-hot encode them as binary variables}
#' }
#' @param print_plot option to print the plot in addition to returning it. Default is `TRUE`
#'
#' @returns A single `ggcorrplot` plot as list object
#' @export
docr_pca_correlation_plot <- function(
  pca_df,
  pca_var,
  sample_col = "sampleId",
  sample_metadata_cols,
  pcs = c("PC1", "PC2", "PC3", "PC4", "PC5"),
  plot_title = "PCA Correlation Matrix",
  factor_encoding = "numeric",
  print_plot = TRUE
) {
  if (!sample_col %in% colnames(pca_df)) {
    stop("sample_col must be a column in provided pca_df")
  }
  if (!all(pcs %in% colnames(pca_var)) && all(pcs %in% colnames(pca_df))) {
    stop("pcs must be at least one PC-identifying column name in pca_df and pca_var")
  }
  # If sample column is included in sample_metadata_cols, remove for factor encoding
  if (sample_col %in% sample_metadata_cols) {
    sample_metadata_cols <- setdiff(sample_metadata_cols, sample_col)
  }

  purrr::walk(
    sample_metadata_cols,
    ~ checkmate::assertChoice(.x, colnames(pca_df))
  )

  checkmate::assertChoice(factor_encoding, c("numeric", "onehot"))
  checkmate::assertLogical(print_plot, len = 1)

  curr_pcax <- pca_df %>%
    dplyr::select(tidyselect::all_of(c(sample_col, pcs))) %>%
    stats::na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(sample_col)

  curr_meta <- pca_df %>%
    dplyr::select(tidyselect::any_of(c(sample_col, sample_metadata_cols))) %>%
    stats::na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(sample_col)

  # If rows are removed from na.omit(), make sure they are removed
  # from corresponding pca or metadata dataframe
  rownames_intersect <- intersect(rownames(curr_pcax), rownames(curr_meta))
  curr_pcax <- subset(curr_pcax, rownames(curr_pcax) %in% rownames_intersect)
  curr_meta <- subset(curr_meta, rownames(curr_meta) %in% rownames_intersect)

  for (var in sample_metadata_cols) {
    if (factor_encoding == "numeric") {
      curr_meta[[var]] <- as.numeric(as.factor(curr_meta[[var]]))
    } else {
      resultsvec_ <- curr_meta[[var]]
      if (class(resultsvec_) %in% c("character", "factor", "ordered")) {
        if (length(unique(resultsvec_)) > 2) {
          curr_meta <- curr_meta %>% dplyr::select(-!!rlang::sym(var))
          onehot_encoding <- stats::model.matrix(~ resultsvec_ + 0)
          colnames(onehot_encoding) <- gsub("resultsvec", var, colnames(onehot_encoding))
          curr_meta <- cbind(curr_meta, onehot_encoding)
        } else {
          curr_meta[[var]] <- as.numeric(as.factor(curr_meta[[var]]))
        }
      }
    }
  }

  # First, compute correlations between metadata and PCA data
  corrs <- stats::cor(curr_pcax, curr_meta, method = "spearman")

  # Make correlation plot between metadata and PCA data
  corr_plot <-
    ggcorrplot::ggcorrplot((corrs),
      lab = T,
      lab_size = 2.5,
      hc.order = F,
      method = "circle"
    ) +
    ggplot2::labs(title = plot_title) +
    ggplot2::scale_x_discrete(
      breaks = rownames(corrs),
      labels = paste0(rownames(corrs), " (", round(pca_var[1, rownames(corrs)] * 100, 1), " %)")
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  if (print_plot) {
    print(corr_plot)
  }

  return(list(corr_plot))
}
