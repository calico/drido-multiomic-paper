#' DO CR: process plate
#'
#' @export
docr_process_dims_plate <- function(
  plate_name,
  samples_file_path,
  lib_dir,
  is_lib_name,
  biological_lib_name,
  ms3_lib_name,
  bulkpool_compounds,
  is_search_params,
  bulkpool_params,
  biological_ms2_search_params,
  biological_ms3_search_params,
  rds_output_dir,
  mzrolldb_output_dir,
  is_ms3 = FALSE,
  save_mzrolldb_as_rds = FALSE
) {
  # files ######################################
  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  habc_lib_file <- paste(lib_dir, biological_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")
  is_ms3 <- !is.null(ms3_lib_name)
  method_name <- ifelse(is_ms3, "X0200-M014B", "X0200-M015B")
  plate_mode <- ifelse(is_ms3, "pos", "neg")

  # contains MS3, single measurement of post-partitioned MS1/MS2 compound ion with description
  rds_results_file <- paste0(rds_output_dir, "/", method_name, "_", plate_name, "_rows.rds")

  # Only MS1/MS2 compound ion. Contains MS1 and MS2 partitioned compound intensities,
  # along with ms1_intensity_is_nearest_scan_normalized and diagnostic_ms2_intensity
  col_quant_results_file <- paste0(rds_output_dir, "/", method_name, "_", plate_name, "_cols.rds")

  suffix <- ifelse(save_mzrolldb_as_rds == TRUE, ".rds", ".mzrollDB")
  mzrolldb_results_file <- paste0(mzrolldb_output_dir, "/", method_name, "_", plate_name, suffix)

  # samples ####################################
  samples_info <- docr_samples_df(samples_file_path, is_ms3, TRUE)
  samples_df <- samples_info$samples_df %>% dplyr::arrange(desc(type), order_num)
  habc4_ms2_ranges <- samples_info$ms2_ranges
  all_bulkpool_samples <- samples_info$samples_df %>% dplyr::filter(type == "BulkPool")

  # libraries ###################################
  lib_subset_string <- ifelse(is_ms3, "\\+$", "\\-$")
  is_lib <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl(lib_subset_string, adductName))
  is_lib_sliced <- mzkitcpp::DI_slice_library(habc4_ms2_ranges, is_lib)
  habc_lib <- mzkitcpp::import_msp_lipids_library(habc_lib_file)

  # good bulkpool samples ######################

  # MS2
  good_ms2_bulkpools <- habc4_good_samples_ms2(
    samples = all_bulkpool_samples$file,
    ms2_ranges = habc4_ms2_ranges,
    is_sliced_lib = is_lib_sliced,
    is_search_params = is_search_params,
    adducts_file = adducts_file,
    bulkpool_params = bulkpool_params
  )

  samples_df_ms2_filtered <- samples_df %>%
    dplyr::filter((type == "BulkPool" & sample_name %in% good_ms2_bulkpools) | type != "BulkPool")

  # MS3
  if (is_ms3) {
    good_ms3_bulkpools <- habc4_good_samples_ms3(
      samples = all_bulkpool_samples$file,
      is_search_params = is_search_params,
      search_params = biological_ms3_search_params,
      adducts_file = adducts_file,
      bulkpool_params = bulkpool_params
    )

    samples_df_ms3_filtered <- samples_df %>%
      dplyr::filter((type == "BulkPool" & sample_name %in% good_ms3_bulkpools) | type != "BulkPool")
  } else {
    habc4_bulkpool_ms3_good_samples <- c()
  }

  # MS2 searches ###############################

  is_plate_valid_ms2 <- TRUE
  if (is_plate_valid_ms2) {
    # IS library
    lib_subset_string <- ifelse(is_ms3, "\\+$", "\\-$")
    is_lib <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>%
      dplyr::filter(grepl(lib_subset_string, adductName))
    is_lib_sliced <- mzkitcpp::DI_slice_library(habc4_ms2_ranges, is_lib)

    # IS search
    habc4_IS_results <- mzkitcpp::DI_pipeline(
      samples = samples_df$file,
      ms2_ranges = habc4_ms2_ranges,
      is_sliced_lib = is_lib_sliced,
      is_search_params = is_search_params,
      sliced_lib = is_lib_sliced,
      search_params = is_search_params,
      adducts_file = adducts_file,
      debug = FALSE
    )

    # IS quant ions
    is_quant_table <- to_quant_table(habc4_IS_results$adduct_table, FALSE, TRUE) %>%
      dplyr::select(sample, compoundName, adductName, ms1_intensity, ms2_intensity)

    is_search_table_subset <- habc4_IS_results$search %>%
      dplyr::select(sample, compoundName, adductName, ms1_intensity_is_nearest_scan_normalized) %>%
      unique() %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized))

    is_adduct_table_subset <- habc4_IS_results$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::select(sample, compoundName, adductName, ms2_diagnostic_norm_intensity) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity)) %>%
      dplyr::rename(diagnostic_ms2_intensity = ms2_diagnostic_norm_intensity)

    is_quant_ions <- is_quant_table %>%
      dplyr::full_join(is_search_table_subset, by = c("compoundName", "adductName", "sample")) %>%
      dplyr::full_join(is_adduct_table_subset, by = c("compoundName", "adductName", "sample"))

    # MS2 library
    habc_lib_full <- mzkitcpp::import_msp_lipids_library(habc_lib_file)
    habc_lib_full_w_mc <- habc_lib_full %>%
      dplyr::mutate(monoCompoundName = gsub("_13C", "", compoundName)) %>%
      dplyr::inner_join(bulkpool_compounds, by = c("monoCompoundName" = "compoundName"))

    for (i in 1:nrow(habc_lib_full_w_mc)) {
      adduct_name <- habc_lib_full_w_mc$adductName[i]
      M1_adducts <- habc_lib_full_w_mc$M1_adducts[i]

      compound_name <- habc_lib_full_w_mc$compoundName[i]

      if (grepl(adduct_name, M1_adducts, fixed = TRUE)) {
        habc_lib_full_w_mc$compoundName[i] <- paste0(compound_name, "_M1")
      }
    }

    habc_lib_all <- habc_lib_full_w_mc %>%
      dplyr::filter(!grepl("_13C_M1", compoundName)) %>%
      dplyr::select(-monoCompoundName)

    habc_lib_sliced <- mzkitcpp::DI_slice_library(habc4_ms2_ranges, habc_lib_all)

    # MS2 biological search
    habc4_biological_search <- mzkitcpp::DI_pipeline(
      samples = samples_df_ms2_filtered$file,
      ms2_ranges = habc4_ms2_ranges,
      is_sliced_lib = is_lib_sliced,
      is_search_params = is_search_params,
      sliced_lib = habc_lib_sliced,
      search_params = biological_ms2_search_params,
      adducts_file = adducts_file,
      debug = FALSE
    )

    nearest_scan_IS <- habc4_biological_search$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized) &
        ms1_intensity_is_nearest_scan_normalized > 0 &
        ms1_intensity_is_nearest_scan_normalized < Inf)

    centered_nearest_scan_IS <- docr_dims_bulkpool_centered_quant_ion_table(nearest_scan_IS, "ms1_intensity_is_nearest_scan_normalized") %>%
      dplyr::select(lipidClass, compositionSummary, compoundName, adductName, sample, ms1_intensity_is_nearest_scan_normalized) %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized))

    diagnostic_ms2 <- habc4_biological_search$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity) &
        ms2_diagnostic_norm_intensity > 0 &
        ms2_diagnostic_norm_intensity < Inf)

    centered_diagnostic_ms2 <- docr_dims_bulkpool_centered_quant_ion_table(diagnostic_ms2, "ms2_diagnostic_norm_intensity") %>%
      dplyr::select(lipidClass, compositionSummary, compoundName, adductName, sample, ms2_diagnostic_norm_intensity) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity))

    # center after determination of preferred quant type to ensure accurate centering values
    habc4_quant_table <- to_quant_table(habc4_biological_search$adduct_table, FALSE, TRUE)

    # MS1
    habc4_quant_ions_ms1 <- habc4_quant_table %>%
      dplyr::select(-ms2_intensity, -ms2_quant_type) %>%
      dplyr::filter(!is.na(ms1_intensity) & ms1_intensity > 0 & ms1_intensity < Inf) %>%
      dplyr::mutate(ms1_quant_type = ifelse(grepl("^ms1", ms1_quant_type), ms1_quant_type, paste0("ms1_", ms1_quant_type))) %>%
      dplyr::rename(
        quant_type = ms1_quant_type,
        intensity = ms1_intensity
      ) %>%
      dplyr::mutate(quant_class = "ms1")

    # MS2
    habc4_quant_ions_ms2 <- habc4_quant_table %>%
      dplyr::select(-ms1_intensity, -ms1_quant_type) %>%
      dplyr::filter(!is.na(ms2_intensity) & ms2_intensity > 0 & ms2_intensity < Inf) %>%
      dplyr::mutate(ms2_quant_type = ifelse(grepl("^ms2", ms2_quant_type), ms2_quant_type, paste0("ms2_", ms2_quant_type))) %>%
      dplyr::rename(
        quant_type = ms2_quant_type,
        intensity = ms2_intensity
      ) %>%
      dplyr::mutate(quant_class = "ms2")

    # Center quant ions
    # Johanna modified function 2025-01-13
    habc4_centered_quant_ions_ms1 <- docr_dims_bulkpool_centered_quant_ion_table(habc4_quant_ions_ms1, "intensity")
    habc4_centered_quant_ions_ms2 <- docr_dims_bulkpool_centered_quant_ion_table(habc4_quant_ions_ms2, "intensity")

    di_ms1_quant_ions_condensed <- habc4_centered_quant_ions_ms1 %>%
      dplyr::rename(
        ms1_intensity = intensity,
        ms1_intensity_median = intensity_median
      ) %>%
      dplyr::select(sample, lipidClass, compositionSummary, compoundName, adductName, ms1_intensity, ms1_intensity_median) %>%
      dplyr::filter(!is.na(ms1_intensity))

    di_ms2_quant_ions_condensed <- habc4_centered_quant_ions_ms2 %>%
      dplyr::rename(
        ms2_intensity = intensity,
        ms2_intensity_median = intensity_median
      ) %>%
      dplyr::select(sample, lipidClass, compositionSummary, compoundName, adductName, ms2_intensity, ms2_intensity_median) %>%
      dplyr::filter(!is.na(ms2_intensity))

    # formatted for mzrollDB
    habc4_quant_ions <- di_ms1_quant_ions_condensed %>%
      dplyr::full_join(di_ms2_quant_ions_condensed, by = c("sample", "lipidClass", "compositionSummary", "compoundName", "adductName")) %>%
      dplyr::full_join(centered_nearest_scan_IS, by = c("sample", "lipidClass", "compositionSummary", "compoundName", "adductName")) %>%
      dplyr::full_join(centered_diagnostic_ms2, by = c("sample", "lipidClass", "compositionSummary", "compoundName", "adductName")) %>%
      dplyr::rename(diagnostic_ms2_intensity = ms2_diagnostic_norm_intensity)

    # formatted for stage 4 processing
    habc_centered_quant_ions <- rbind(habc4_centered_quant_ions_ms1, habc4_centered_quant_ions_ms2) %>%
      dplyr::inner_join(samples_df, by = c("sample" = "sample_name")) %>%
      dplyr::select(-file, -type_id) %>%
      unique()

    # label compounds for mzrollDB
    quant_labels <- biological_ms2_search_params$ms1IonList %>%
      dplyr::filter(ms1IsPreferredQuantIon) %>%
      dplyr::mutate(label = "l") %>%
      dplyr::select(lipidClass, adductName, label)

    compounds_13C_IS <- habc4_biological_search$search %>%
      dplyr::filter(grepl("_13C", compoundName) | grepl("_IS", compoundName)) %>%
      dplyr::mutate(label = "")

    compounds_no_13C_IS <- habc4_biological_search$search %>%
      dplyr::filter(!(grepl("_13C", compoundName) | grepl("_IS", compoundName))) %>%
      dplyr::left_join(quant_labels, by = c("lipidClass", "adductName")) %>%
      dplyr::mutate(label = ifelse(is.na(label), "", label))

    labeled_search_results <- rbind(compounds_13C_IS, compounds_no_13C_IS)
  } else {
    is_quant_ions <- NULL
    is_lib_sliced <- NULL
    habc4_IS_results <- NULL
    habc_lib_sliced <- NULL
    habc_centered_quant_ions <- NULL
    habc4_quant_ions <- NULL
    habc4_quant_table <- NULL
  }

  # MS3 searches ###################################
  is_plate_valid_ms3 <- is_ms3
  if (is_plate_valid_ms3) {
    # MS3 library
    habc_ms3_lib <- mzkitcpp::import_msp_lipids_library(file.path(lib_dir, ms3_lib_name))

    habc_ms3_lib_all_targets <- habc_ms3_lib %>%
      to_ms3_lib() %>%
      dplyr::select(prec_mzs) %>%
      unique()

    ms3_targets <- dplyr::inner_join(habc_ms3_lib_all_targets, samples_info$ms3_targets, by = c("prec_mzs"))

    habc_ms3_lib_filtered <- habc_ms3_lib %>%
      to_ms3_lib() %>%
      dplyr::filter(prec_mzs %in% ms3_targets$prec_mzs) %>%
      to_ms2_lib() %>%
      # for HABC, the TGs were added to the bulkpool_compounds list.
      # for DOCR, the TG library was re-made to include only compounds of interest.
      # dplyr::filter(compoundName %in% bulkpool_compounds$compoundName) %>%

      dplyr::mutate(ms2_intensity = 1) %>%
      dplyr::select(colnames(habc_ms3_lib))

    # MS3 biological search
    habc4_biological_ms3_search <- mzkitcpp::DI_pipeline_ms3_search(
      samples = samples_df_ms3_filtered$file,

      # package data
      is_lib = default_tg_is_ms3,
      is_search_params = biological_ms3_search_params,
      search_lib = habc_ms3_lib_filtered,
      search_params = biological_ms3_search_params,
      adducts_file = adducts_file,
      debug = F
    )

    # center results
    habc4_biological_ms3_median_column <- habc4_biological_ms3_search %>%
      dplyr::filter(!is.na(ms3_intensity_sum_norm) & ms3_intensity_sum_norm > 0 & ms3_intensity_sum_norm < Inf) %>%
      dplyr::select(sample, compoundName, adductName, ms3_intensity_sum_norm) %>%
      unique() %>%
      dplyr::group_by(compoundName, adductName) %>%
      dplyr::mutate(intensity_median = median(ms3_intensity_sum_norm, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(compoundName, adductName, intensity_median)

    habc4_biological_ms3_search_centered <- dplyr::inner_join(
      habc4_biological_ms3_search, habc4_biological_ms3_median_column,
      by = c("compoundName", "adductName")
    ) %>%
      dplyr::mutate(ms3_intensity_sum_norm = ms3_intensity_sum_norm / intensity_median) %>%
      unique()

    # quant table
    habc4_ms3_quant_table <- habc4_biological_ms3_search_centered %>%
      dplyr::inner_join(samples_df, by = c("sample" = "sample_name")) %>%
      dplyr::mutate(
        intensity = ms3_intensity_sum_norm,
        quant_type = "ms3_intensity_sum_norm",
        quant_class = "ms3"
      ) %>%
      dplyr::select(
        sample, lipidClass, compositionSummary, compoundName, adductName,
        intensity, intensity_median, quant_type, quant_class,
        method, mode, plate, order_num, type, barcode, well_position, injection_num
      ) %>%
      unique()
  } else {
    habc4_ms3_quant_table <- NULL
  }

  # color table ###################################
  color_table <- docr_type_color_table(samples_df, good_ms2_bulkpools)
  plate_name_vector <- rep(plate_name, nrow(samples_df))

  # save RDS results ###################################
  system(glue::glue("rm {old_rds_results_file} 2>&1", old_rds_results_file = rds_results_file))

  rds_results <- NULL
  if (!is.null(habc_centered_quant_ions) && !is.null(habc4_ms3_quant_table)) {
    rds_results <- rbind(habc_centered_quant_ions, habc4_ms3_quant_table)
  } else if (!is.null(habc_centered_quant_ions) && is.null(habc4_ms3_quant_table)) {
    rds_results <- habc_centered_quant_ions
  } else if (is.null(habc_centered_quant_ions) && !is.null(habc4_ms3_quant_table)) {
    rds_results <- habc4_ms3_quant_table
  }

  saveRDS(rds_results, file = rds_results_file)
  if (!is.null(habc_centered_quant_ions)) {
    saveRDS(habc4_quant_ions, file = col_quant_results_file)
  }

  # save mzrollDB results ###################################
  system(glue::glue("rm {old_mzroll_db_file} 2>&1", old_mzroll_db_file = mzrolldb_results_file))

  encoded_quantType <- DIMS_encoded_quant_type()

  if (save_mzrolldb_as_rds) {
    IS_search_results <- "NULL"
    IS_adduct_table <- "NULL"
    IS_quant_ions <- "NULL"
    if (!is.null(habc4_IS_results)) {
      IS_search_results <- habc4_IS_results$search
      IS_adduct_table <- habc4_IS_results$adduct_table
      IS_quant_ions <- is_quant_ions
    }

    HABC_search_results <- "NULL"
    HABC_adduct_table <- "NULL"
    HABC_quant_ions <- "NULL"
    if (!is.null(habc4_quant_table)) {
      HABC_search_results <- labeled_search_results
      HABC_adduct_table <- habc4_biological_search$adduct_table
      HABC_quant_ions <- habc4_quant_ions
    }

    MS3_search_results <- "NULL"
    if (!is.null(habc4_ms3_quant_table)) {
      MS3_search_results <- habc4_biological_ms3_search_centered
    } else { # Johanna added 2025-10-13
      samples_df_ms3_filtered <- NULL
      samples_df_ms3_filtered$file <- NULL
    }

    mzrolldb_results <-
      list(
        "ms2_ranges" = habc4_ms2_ranges,
        "adducts_file" = adducts_file,
        "color_table" = color_table,
        "encoded_quantType" = encoded_quantType,
        "IS_samples" = samples_df$file,
        "IS_library_name" = is_lib_file,
        "IS_search_lib_sliced" = is_lib_sliced,
        "IS_search_params" = is_search_params,
        "IS_search_results" = IS_search_results,
        "IS_adduct_table" = IS_adduct_table,
        "IS_quant_ions" = IS_quant_ions,
        "IS_set_name" = plate_name_vector,
        "HABC_samples" = samples_df_ms2_filtered$file,
        "HABC_library_name" = biological_lib_name,
        "HABC_search_lib_sliced" = habc_lib_sliced,
        "HABC_search_params" = biological_ms2_search_params,
        "HABC_search_results" = HABC_search_results,
        "HABC_adduct_table" = HABC_adduct_table,
        "HABC_quant_ions" = HABC_quant_ions,
        "HABC_set_name" = plate_name_vector,
        "MS3_samples" = samples_df_ms3_filtered$file,
        "MS3_library_name" = ms3_lib_name,
        "MS3_search_params" = biological_ms3_search_params,
        "MS3_search_results" = MS3_search_results,
        "MS3_set_name" = plate_name_vector
      )

    saveRDS(mzrolldb_results, file = mzrolldb_results_file)
  } else {
    if (!is.null(habc4_IS_results)) {
      # Add IS search results
      add_direct_infusion_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df$file,
        ms2_ranges = habc4_ms2_ranges,
        library_name = is_lib_file,
        search_lib_sliced = is_lib_sliced,
        search_params = is_search_params,
        adducts_file = adducts_file,
        di_search_results = habc4_IS_results$search,
        di_quant_ions = is_quant_ions,
        set_name = plate_name_vector,
        color_table = color_table
      )
    }

    if (!is.null(habc4_quant_table)) {
      # Add biological search
      add_direct_infusion_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df_ms2_filtered$file,
        ms2_ranges = habc4_ms2_ranges,
        library_name = biological_lib_name,
        search_lib_sliced = habc_lib_sliced,
        search_params = biological_ms2_search_params,
        adducts_file = adducts_file,
        di_search_results = labeled_search_results,
        di_quant_ions = habc4_quant_ions,
        set_name = plate_name_vector,
        color_table = color_table
      )
    }

    # Add MS3 search

    if (!is.null(habc4_ms3_quant_table)) {
      add_targeted_ms3_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df_ms3_filtered$file,
        ms3_search_results = habc4_biological_ms3_search_centered,
        library_name = ms3_lib_name,
        search_params = biological_ms3_search_params,
        set_name = plate_name_vector,
        color_table = color_table
      )
    }

    mzroll_db_con <- DBI::dbConnect(RSQLite::SQLite(), mzrolldb_results_file)
    ui_options <- tibble::tibble(
      key = c("quantType"),
      value = c(encoded_quantType)
    )
    DBI::dbAppendTable(mzroll_db_con, "ui", ui_options)
    DBI::dbDisconnect(mzroll_db_con)
  }

  return(invisible(0))
}


#' DO CR: Return table of quantitative measurements from collection of results.
#'
#' @description
#' combine stage 3 RDS result files into final table of quant results.
#'
#' @param stage_3_rds_dir collection of saved RDS results files.
#' @param bulkpool_params collection of params associated with bulkpool search.  Used here for quant ions.
#'
#' @export
docr_dims_formatting <- function(
  stage_3_rds_dir,
  bulkpool_params
) {
  # v4: quant ions explicltly come from bulkpool_params
  # autosummarization correction, when necessary
  # CVs should not be assessed until after summarization correction

  # quant ions ###############################
  quant_compound_ions <- rbind(bulkpool_params$ms1IonList, bulkpool_params$ms3IonList) %>%
    dplyr::filter(ms1IsPreferredQuantIon == TRUE) %>%
    dplyr::select(-ms1IsMajorIon, -ms1IsPreferredQuantIon)

  # import results ###############################
  row_stage3_results <- list.files(stage_3_rds_dir, pattern = "*_rows.rds$", full.names = TRUE)
  col_stage3_results <- list.files(stage_3_rds_dir, pattern = "*_cols.rds$", full.names = TRUE)

  stage3_quant_ions_cols <- NULL
  stage3_quant_ions <- NULL

  for (i in 1:length(row_stage3_results)) {
    # Extract quant data for each method:plate
    # Remove internal standards
    plate_quant_data <- readRDS(row_stage3_results[i]) %>%
      dplyr::inner_join(quant_compound_ions, by = c("lipidClass", "adductName")) %>%
      dplyr::filter(!grepl("_13C", compoundName) & !grepl("_IS", compoundName))

    col_plate_quant_data <- readRDS(col_stage3_results[i]) %>%
      dplyr::filter(!grepl("_13C", compoundName) & !grepl("_IS", compoundName))

    # Assign mode depending on adducts
    plate_quant_data <- plate_quant_data %>%
      dplyr::filter((grepl("\\+$", adductName) & mode == "pos") | (grepl("\\-$", adductName) & mode == "neg"))

    plate_mode <- plate_quant_data$mode[1]
    plate_name <- plate_quant_data$plate[1]

    # On first loop, assign quant data to stage3_quant_ions
    # After first loop, bind quant data together
    if (is.null(stage3_quant_ions)) {
      stage3_quant_ions <- plate_quant_data
      stage3_quant_ions_cols <- col_plate_quant_data
    } else {
      stage3_quant_ions <- rbind(stage3_quant_ions, plate_quant_data)
      stage3_quant_ions_cols <- rbind(stage3_quant_ions_cols, col_plate_quant_data)
    }
  } # End loop through all plates

  # Remove reruns ###################################
  all_samples <- stage3_quant_ions %>%
    dplyr::select(sample, method, type, plate, barcode) %>%
    unique() %>%
    # dplyr::mutate(plate_num = as.integer(stringr::str_extract(plate, "(?<=[pP])\\d+"))) # old
    dplyr::mutate(plate_num = as.integer(factor(plate,
      levels = unique(plate[order(
        substr(plate, 2, 3),
        substr(plate, 4, 4)
      )])
    ))) # New

  # Select barcode from latest run, determined as highest entry from plate_num
  biological_cleaned_samples <- all_samples %>%
    dplyr::filter(type == "biological_sample") %>%
    dplyr::arrange(method, barcode, desc(plate_num)) %>%
    dplyr::group_by(method, barcode) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # Select all non-biological samples
  pass_through_samples <- all_samples %>%
    dplyr::filter(type != "biological_sample")

  # Select on sample names from all non-biological, and most recent re-run samples
  rerun_adjusted_samples <- rbind(biological_cleaned_samples, pass_through_samples) %>%
    dplyr::select(sample)

  # Select quant data on data with re-runs removed
  stage3_quant_ions_rerun_adjusted <- stage3_quant_ions %>%
    dplyr::inner_join(rerun_adjusted_samples, by = c("sample"))
  stage3_quant_ions_cols_rerun_adjusted <- stage3_quant_ions_cols %>%
    dplyr::inner_join(rerun_adjusted_samples, by = c("sample"))

  # Summarization ###################################

  # Select on summarized compounds to summarize with grepl("^\\{", compoundName)
  stage3_quant_ions_to_summarize <- stage3_quant_ions_rerun_adjusted %>%
    dplyr::filter(grepl("^\\{", compoundName) & type == "biological_sample") %>%
    dplyr::select(compositionSummary, adductName) %>%
    unique()

  stage3_quant_ions_unaffected <- stage3_quant_ions_rerun_adjusted %>%
    dplyr::anti_join(stage3_quant_ions_to_summarize, by = c("compositionSummary", "adductName"))

  # stage3_quant_ions_to_alter <- stage3_quant_ions_rerun_adjusted %>%
  #  dplyr::inner_join(stage3_quant_ions_to_summarize, by = c("compositionSummary", "adductName"))

  stage3_quant_ions_cols_altered_subset <- stage3_quant_ions_cols_rerun_adjusted %>%
    dplyr::inner_join(stage3_quant_ions_to_summarize, by = c("compositionSummary", "adductName"))

  sample_info <- stage3_quant_ions_rerun_adjusted %>%
    dplyr::select(sample, method, mode, plate, order_num, type, barcode, well_position, injection_num) %>%
    unique()

  stage3_quant_ions_altered <- stage3_quant_ions_cols_altered_subset %>%
    dplyr::select(sample, lipidClass, compositionSummary, adductName, ms1_intensity_is_nearest_scan_normalized) %>%
    unique() %>%
    tidyr::pivot_longer(ms1_intensity_is_nearest_scan_normalized) %>%
    dplyr::select(-name) %>%
    dplyr::rename(intensity = value) %>%
    dplyr::mutate(quant_class = "ms1", quant_type = "ms1_intensity", compoundName = compositionSummary) %>%
    dplyr::inner_join(sample_info, by = c("sample" = "sample")) %>%
    dplyr::select(
      sample, lipidClass, compositionSummary, compoundName, adductName,
      intensity, quant_type, quant_class, method, mode, plate,
      order_num, type, barcode, well_position, injection_num
    )

  stage3_quant_ions_summarized_corrected <- rbind(stage3_quant_ions_unaffected, stage3_quant_ions_altered) %>%
    dplyr::filter(!grepl("^\\{", compoundName)) # remove any controls or bulkpool samples with summarized compounds.

  # add log2 intensity
  stage3_quant_ions_summarized_corrected <- stage3_quant_ions_summarized_corrected %>%
    dplyr::mutate(log2_intensity = log2(intensity))

  # Stage 3 Quant CVs ###############################
  plate_modes <- stage3_quant_ions_summarized_corrected %>%
    dplyr::select(plate, mode) %>%
    unique()

  stage3_quant_cvs <- NULL

  for (i in 1:nrow(plate_modes)) {
    ith_plate <- plate_modes$plate[i]
    ith_mode <- plate_modes$mode[i]

    plate_quant_data <- stage3_quant_ions_summarized_corrected %>%
      dplyr::filter(plate == ith_plate & mode == ith_mode)

    plate_quant_cvs <- habc_cv_comparison_v3(plate_quant_data, "biological_sample") %>%
      dplyr::mutate(mode = ith_mode, plate = ith_plate)

    if (is.null(stage3_quant_cvs)) {
      stage3_quant_cvs <- plate_quant_cvs
    } else {
      stage3_quant_cvs <- rbind(stage3_quant_cvs, plate_quant_cvs)
    }
  }

  # Median BulkPool CV ###############################
  BulkPool_CVs <- NULL

  ms_levels <- c("ms1", "ms2", "ms3")
  for (i in 1:length(ms_levels)) {
    ms_level <- ms_levels[i]

    ms_level_BulkPool_CVs <- stage3_quant_cvs %>%
      dplyr::filter(quant_class == ms_level) %>%
      dplyr::group_by(compoundName, adductName) %>%
      dplyr::mutate(median_BulkPool_CV = median(BulkPool_CV, na.rm = TRUE)) %>%
      dplyr::ungroup()

    if (is.null(BulkPool_CVs)) {
      BulkPool_CVs <- ms_level_BulkPool_CVs
    } else {
      BulkPool_CVs <- rbind(BulkPool_CVs, ms_level_BulkPool_CVs)
    }
  }

  stage3_median_bulkpool_CVs <- BulkPool_CVs %>%
    dplyr::select(compoundName, adductName, mode, quant_class, median_BulkPool_CV) %>%
    unique() %>%
    dplyr::arrange(compoundName, adductName, mode) %>%
    dplyr::group_by(compoundName) %>%
    dplyr::mutate(is_lowest_MS_CV = median_BulkPool_CV == min(median_BulkPool_CV, na.rm = TRUE)) %>%
    dplyr::ungroup()

  stage3_quant_cvs_plate <- stage3_quant_cvs %>%
    dplyr::rename(plate_BulkPool_CV = BulkPool_CV, plate_biological_sample_CV = bc_CV)

  # values have been centered, returning IS norm centered intensity and log2 IS norm centered intensity
  quantified_compounds <- stage3_quant_ions_summarized_corrected %>%
    dplyr::inner_join(stage3_median_bulkpool_CVs, by = c("compoundName", "adductName", "mode", "quant_class")) %>%
    dplyr::left_join(stage3_quant_cvs_plate, by = c("compoundName", "adductName", "mode", "plate", "quant_class")) %>%
    dplyr::rename(
      is_lowest_quant_class_CV = is_lowest_MS_CV,
      centered_IS_norm_intensity = intensity,
      log2_centered_IS_norm_intensity = log2_intensity
    ) %>%
    dplyr::select(
      sample, method, mode, plate, type, barcode, well_position, order_num, injection_num,
      lipidClass, compositionSummary, compoundName, adductName,
      centered_IS_norm_intensity, log2_centered_IS_norm_intensity,
      intensity_median, quant_class, quant_type, is_lowest_quant_class_CV,
      plate_BulkPool_CV, plate_biological_sample_CV, median_BulkPool_CV
    )

  return(quantified_compounds)
}


#' Bulkpool-centered adduct table
#'
#' @description
#' Center plate quant values after applying partitioning and other intensity transformations
#'
#' @param search_table
#'     Search table from HABC plate search
#'
#' @param samples_df
#'     HABC samples df as generated by \code{habc_samples_df()}
#'
#' @export
docr_dims_bulkpool_centered_quant_ion_table <- function(quant_ion_table, intensity_col_name) {
  quant_ion_medians <- quant_ion_table %>%
    dplyr::filter(grepl("BulkPool", sample)) %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(quant_col := !!sym(intensity_col_name)) %>%
    dplyr::mutate(intensity_median = median(quant_col, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(compoundName, adductName, intensity_median) %>%
    unique()

  centered_quant_ion_table <- quant_ion_table %>%
    dplyr::inner_join(quant_ion_medians, by = c("compoundName", "adductName")) %>%
    dplyr::mutate(!!sym(intensity_col_name) := (!!sym(intensity_col_name)) / intensity_median)

  return(centered_quant_ion_table)
}
