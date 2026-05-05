#' HABC4 Good Samples (MS2)
#'
#' @description
#' internal function to determine if samples are good or not.
#'
#' @param samples
#' vector of sample file paths
#'
#' @param ms2_ranges
#' isolation window information for samples
#'
#' @param is_sliced_lib
#' sliced IS library
#'
#' @param is_search_params
#' IS search parameters
#'
#' @param adducts_file
#' Necessary for MAVEN searches
#'
#' @param bulkpool_params
#' bulkpool params information
#'
#' @export
habc4_good_samples_ms2 <- function(samples,
                                   ms2_ranges,
                                   is_sliced_lib,
                                   is_search_params,
                                   adducts_file,
                                   bulkpool_params) {
  good_samples_ms2_results <- mzkitcpp::DI_pipeline(
    samples = samples,
    ms2_ranges = ms2_ranges,
    is_sliced_lib = is_sliced_lib,
    is_search_params = is_search_params,
    sliced_lib = is_sliced_lib,
    search_params = is_search_params,
    adducts_file = adducts_file,
    debug = FALSE
  )

  good_samples <- character(0)

  if (nrow(good_samples_ms2_results$adduct_table) > 0) {
    good_samples_ms2_results_filtered <- habc_good_bulkpools(
      good_samples_ms2_results$adduct_table,
      bulkpool_params$bulkpoolGoodMinFracIntensity,
      bulkpool_params$bulkpoolGoodMinFracDetected
    )

    good_samples <- good_samples_ms2_results_filtered$sample
  }

  return(good_samples)
}

#' HABC4 Good Samples (MS3)
#'
#' @description
#' internal function to determine if samples are good or not.
#'
#' @param samples
#' vector of sample file paths
#'
#' @param is_search_params
#' IS search parameters
#'
#' @param search_params
#' MS3 search parameters
#'
#' @param adducts_file
#' Necessary for MAVEN searches
#'
#' @param bulkpool_params
#' bulkpool params information
#'
#' @export
habc4_good_samples_ms3 <- function(samples,
                                   is_search_params,
                                   search_params,
                                   adducts_file,
                                   bulkpool_params) {
  good_samples_ms3_results <- mzkitcpp::DI_pipeline_ms3_search(
    samples = samples,
    is_lib = clamr::default_tg_is_ms3,
    is_search_params = is_search_params,
    search_lib = clamr::default_tg_is_ms3,
    search_params = search_params,
    adducts_file = adducts_file,
    debug = F
  )

  good_samples <- character(0)

  if (nrow(good_samples_ms3_results) > 0) {
    # filtering criteria: ms3 m/z matches and sum ms3 intensity
    good_samples_ms3_results_filtered <- good_samples_ms3_results %>%
      dplyr::select(sample, fragmentLabel, num_ms3_mz_matches, ms3_intensity_sum) %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(total_ms3_mz_matches = sum(num_ms3_mz_matches)) %>%
      dplyr::ungroup() %>%
      dplyr::select(sample, total_ms3_mz_matches, ms3_intensity_sum) %>%
      unique() %>%
      dplyr::filter(total_ms3_mz_matches >= bulkpool_params$bulkpoolGoodMs3MinTotalMatches &
        ms3_intensity_sum >= bulkpool_params$bulkpoolGoodMs3MinIntensitySum)

    good_samples <- good_samples_ms3_results_filtered$sample
  }

  return(good_samples)
}

#' HABC4 Default Quant Ions
#'
#' @description
#' Return default quant ions, formatted for writing to mzrollDB.
#' Note that this does not perform any adjustment to raw search results, it simply
#' reformats them into tables with the expected number and type of data columns.
#'
#' @param search_results
#' search results from mzkitcpp::DI_pipeline()
#'
#' @export
habc4_default_quant_ions <- function(search_results) {
  quant_table <- to_quant_table(search_results$adduct_table, FALSE, TRUE) %>%
    dplyr::select(sample, compoundName, adductName, ms1_intensity, ms2_intensity)

  search_table_subset <- search_results$search %>%
    dplyr::select(sample, compoundName, adductName, ms1_intensity_is_nearest_scan_normalized) %>%
    unique() %>%
    dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized))

  adduct_table_subset <- search_results$adduct_table %>%
    dplyr::filter(is_identified == TRUE) %>%
    dplyr::select(sample, compoundName, adductName, ms2_diagnostic_norm_intensity) %>%
    dplyr::filter(!is.na(ms2_diagnostic_norm_intensity)) %>%
    dplyr::rename(diagnostic_ms2_intensity = ms2_diagnostic_norm_intensity)

  is_quant_ions <- quant_table %>%
    dplyr::full_join(search_table_subset, by = c("compoundName", "adductName", "sample")) %>%
    dplyr::full_join(adduct_table_subset, by = c("compoundName", "adductName", "sample"))
}

#' HABC4 Compound Counts
#'
#' @description
#' When assessing presence/absence of a compound ion from a sample,
#' consider if the compound ion was found using a monoisotopic precursor mz,
#' or if the compound ion could only be found using the [M+1] precursor mz
#' (presumably as a result of ion coalescence.)
#'
#' @param compound_list
#' list of compounds saved from a stage 1 search
#'
#' @export
habc4_compound_counts <- function(compound_list) {
  compound_list_backup <- compound_list %>% dplyr::filter(is_13C_precursor == TRUE)

  compound_list_backup_summary <- compound_list_backup %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(num_samples_m_plus_one = n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(compoundName, adductName, num_samples_m_plus_one) %>%
    unique()

  compound_list_mono <- compound_list %>% dplyr::filter(is_13C_precursor == FALSE)

  compound_list_mono_summary <- compound_list_mono %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(num_samples_mono = n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(compoundName, adductName, num_samples_mono) %>%
    unique()

  compound_list_full <- compound_list_backup_summary %>%
    dplyr::full_join(., compound_list_mono_summary, by = c("compoundName", "adductName")) %>%
    dplyr::mutate(num_samples_mono = ifelse(is.na(num_samples_mono), 0, num_samples_mono)) %>%
    dplyr::mutate(num_samples_m_plus_one = ifelse(is.na(num_samples_m_plus_one), 0, num_samples_m_plus_one)) %>%
    dplyr::mutate(num_samples = num_samples_mono + num_samples_m_plus_one) %>%
    dplyr::arrange(-num_samples, compoundName, adductName)
}


#' HABC v4 stage 1: bulkpool samples analysis
#'
#' @param plate_name
#'        name of folder containing mzML files associated with a single HABC plate
#'
#' @param samples_file_path
#'        directory of HABC mzML files
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing IS info for nornmalization  (lib_dir/is_lib_name)
#'
#' @param habc_lib_name
#'        HABC library name (lib_dir/habc_lib_name)
#'
#' @param ms3_lib_name
#'        MS3 library name (lib_dir/ms3_lib_name)
#'
#' @param is_search_params
#'        search parameters for internal standard (IS) search.
#'
#' @param bulkpool_params
#'        Parameters specific to handling of bulkpool samples.
#'        Does not control how specific lipids are identified in bulkpool samples.
#'
#' @param whitelist_ms2_search_params
#'        Parameters associated with identification of lipids in bulkpool samples.
#'
#' @param whitelist_ms3_search_params
#'        Parameters associated with identification of lipids from targeted MS3 scans.
#'
#' @param stage_1_results_dir information about all bulkpool samples, and compounds identified therein
#' Three kinds of files are created in this directory:
#' "samples", "compounds", and "ctl".
#' Files are named as <plate>_<method>_<samples|compounds|ctl>.rds
#' Each file is a tibble, with information from the DIMS searches that can be later be
#' used to create a bulkpool library (in stage 2)
#'
#' @export
habc4_stage1_process_plate <- function(
  plate_name,
  samples_file_path,
  lib_dir,
  is_lib_name,
  habc_lib_name,
  ms3_lib_name,
  is_search_params,
  bulkpool_params,
  whitelist_ms2_search_params,
  whitelist_ms3_search_params,
  stage_1_results_dir
) {
  # files ######################################
  samples_files <- list.files(samples_file_path, pattern = ".*.mzML$", full.names = TRUE)
  if (length(samples_files) == 0) {
    cat(paste0("No samples found in samples file path '", samples_file_path, "'\n"))
    return(invisible(0))
  }

  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  habc_lib_file <- paste(lib_dir, habc_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")
  is_ms3 <- !is.null(ms3_lib_name)
  method_name <- ifelse(is_ms3, "X0158_M014A", "X0158_M015A")
  samples_output_file <- file.path(stage_1_results_dir, paste0(plate_name, "_", method_name, "_samples.rds"))
  compounds_output_file <- file.path(stage_1_results_dir, paste0(plate_name, "_", method_name, "_compounds.rds"))
  ctl_output_file <- file.path(stage_1_results_dir, paste0(plate_name, "_", method_name, "_ctl.rds"))

  # samples ####################################
  samples_info <- habc_samples_df(samples_file_path, is_ms3)
  samples_df <- samples_info$samples_df %>% dplyr::arrange(desc(type), order_num)
  habc4_ms2_ranges <- samples_info$ms2_ranges
  all_bulkpool_samples <- samples_df %>% dplyr::filter(type == "BulkPool")
  pos_ctl_samples <- samples_df %>% dplyr::filter(type == "posctl")
  neg_ctl_samples <- samples_df %>% dplyr::filter(type == "negctl")
  all_ctl_samples <- rbind(pos_ctl_samples, neg_ctl_samples)

  # libraries ###################################
  lib_subset_string <- ifelse(is_ms3, "\\+$", "\\-$")
  is_lib <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl(lib_subset_string, adductName))
  is_lib_sliced <- mzkitcpp::DI_slice_library(habc4_ms2_ranges, is_lib)
  habc_lib <- mzkitcpp::import_msp_lipids_library(habc_lib_file)

  # good bulkpool ms2 ###########################
  habc4_bulkpool_ms2_good_samples <- habc4_good_samples_ms2(
    samples = all_bulkpool_samples$file,
    ms2_ranges = habc4_ms2_ranges,
    is_sliced_lib = is_lib_sliced,
    is_search_params = is_search_params,
    adducts_file = adducts_file,
    bulkpool_params = bulkpool_params
  )

  habc4_posctl_ms2_good_samples <- habc4_good_samples_ms2(
    samples = pos_ctl_samples$file,
    ms2_ranges = habc4_ms2_ranges,
    is_sliced_lib = is_lib_sliced,
    is_search_params = is_search_params,
    adducts_file = adducts_file,
    bulkpool_params = bulkpool_params
  )

  # good bulkpool ms3 ###########################
  if (is_ms3) {
    habc4_bulkpool_ms3_good_samples <- habc4_good_samples_ms3(
      samples = all_bulkpool_samples$file,
      is_search_params = is_search_params,
      search_params = whitelist_ms3_search_params,
      adducts_file = adducts_file,
      bulkpool_params = bulkpool_params
    )

    habc4_posctl_ms3_good_samples <- habc4_good_samples_ms3(
      samples = pos_ctl_samples$file,
      is_search_params = is_search_params,
      search_params = whitelist_ms3_search_params,
      adducts_file = adducts_file,
      bulkpool_params = bulkpool_params
    )
  } else {
    habc4_bulkpool_ms3_good_samples <- c()
    habc4_posctl_ms3_good_samples <- c()
  }

  # bulkpool sample results #####################
  bulkpool_results_table <- all_bulkpool_samples %>%
    dplyr::mutate(
      is_good_ms2 = sample_name %in% habc4_bulkpool_ms2_good_samples,
      is_good_ms3 = sample_name %in% habc4_bulkpool_ms3_good_samples
    ) %>%
    dplyr::select(file, sample_name, method, mode, plate, well_position, is_good_ms2, is_good_ms3)

  posctl_results_table <- pos_ctl_samples %>%
    dplyr::mutate(
      is_good_ms2 = sample_name %in% habc4_posctl_ms2_good_samples,
      is_good_ms3 = sample_name %in% habc4_posctl_ms3_good_samples
    ) %>%
    dplyr::select(file, sample_name, method, mode, plate, well_position, is_good_ms2, is_good_ms3)

  # whitelist ms2 search ########################

  good_ms2_samples <- bulkpool_results_table %>% dplyr::filter(is_good_ms2)
  posctl_ms2_samples <- posctl_results_table %>% dplyr::filter(is_good_ms2)

  class_adduct_counts_table <- whitelist_ms2_search_params$ms1IonList %>%
    dplyr::group_by(lipidClass) %>%
    dplyr::summarize(num_adducts = n())

  ion_list <- whitelist_ms2_search_params$ms1IonList %>% dplyr::select(lipidClass, adductName)

  habc_whitelist_ms2_lib_ion_list <- habc_lib %>%
    dplyr::inner_join(., ion_list, by = c("lipidClass", "adductName")) %>%
    dplyr::filter(!grepl("_M1", compoundName))

  # [M+0]
  habc_whitelist_ms2_lib_mono <- habc_whitelist_ms2_lib_ion_list %>%
    dplyr::filter(!grepl("_13C$", compoundName))

  # [M+1]
  m_plus_one <- whitelist_ms2_search_params$ms1IsRequireMPlusOneByLipidClassAndAdduct %>%
    dplyr::filter(ms1IsRequireMPlusOne == TRUE) %>%
    dplyr::select(lipidClass, adductName)

  habc_whitelist_ms2_lib_m_plus_one <- habc_whitelist_ms2_lib_ion_list %>%
    dplyr::filter(grepl("[^13C]_13C$", compoundName)) %>%
    dplyr::inner_join(m_plus_one, by = c("lipidClass", "adductName"))

  # [M+2]
  m_plus_two <- whitelist_ms2_search_params$ms1IsRequireMPlusTwoByLipidClassAndAdduct %>%
    dplyr::filter(ms1IsRequireMPlusTwo == TRUE) %>%
    dplyr::select(lipidClass, adductName)

  habc_whitelist_ms2_lib_m_plus_two <- habc_whitelist_ms2_lib_ion_list %>%
    dplyr::filter(grepl("_13C_13C$", compoundName)) %>%
    dplyr::inner_join(m_plus_two, by = c("lipidClass", "adductName"))

  habc_whitelist_ms2_lib <- rbind(habc_whitelist_ms2_lib_mono, habc_whitelist_ms2_lib_m_plus_one, habc_whitelist_ms2_lib_m_plus_two)

  habc_whitelist_ms2_lib_sliced <- mzkitcpp::DI_slice_library(habc4_ms2_ranges, habc_whitelist_ms2_lib)

  whitelist_bulkpool_ms2_search <- mzkitcpp::DI_pipeline(
    samples = good_ms2_samples$file,
    ms2_ranges = habc4_ms2_ranges,
    is_sliced_lib = is_lib_sliced,
    is_search_params = is_search_params,
    sliced_lib = habc_whitelist_ms2_lib_sliced,
    search_params = whitelist_ms2_search_params,
    adducts_file = adducts_file,
    debug = FALSE
  )

  whitelist_ctl_ms2_assessment <- tibble::tibble(
    lipidClass = character(0),
    compoundName = character(0),
    adductName = character(0),
    plate = character(0),
    negctl = numeric(0),
    posctl = numeric(0),
    is_pass_ratio = logical(0)
  )

  # If no valid pos ctl samples, skip the pos ctl check
  if (nrow(posctl_ms2_samples) > 0) {
    ctl_samples <- c(posctl_ms2_samples$file, neg_ctl_samples$file)

    whitelist_ctl_ms2_search <- mzkitcpp::DI_pipeline(
      samples = ctl_samples,
      ms2_ranges = habc4_ms2_ranges,
      is_sliced_lib = is_lib_sliced,
      is_search_params = is_search_params,
      sliced_lib = habc_whitelist_ms2_lib_sliced,
      search_params = whitelist_ms2_search_params,
      adducts_file = adducts_file,
      debug = FALSE
    )

    # note that this assessment is currently done on the MS1 m/z intensity
    # (no partitioning)
    # This assumes that if an MS1 m/z is unreliable, MS2-type measurements will also be unreliable
    # It may be possible that in some cases, an MS1 m/z is unreadable while exclusively MS2-based
    # measurements are OK, e.g. diagnostic MS2 intensity and acyl chain measurements
    whitelist_ctl_ms2_assessment <- whitelist_ctl_ms2_search$adduct_table %>%
      dplyr::filter(is_identified) %>%
      dplyr::select(sample, lipidClass, compoundName, adductName, ms1_intensity_is_nearest_scan_normalized) %>%
      unique() %>%
      dplyr::inner_join(all_ctl_samples, by = c("sample" = "sample_name")) %>%
      dplyr::group_by(lipidClass, compoundName, adductName, type) %>%
      dplyr::mutate(median_intensity = median(ms1_intensity_is_nearest_scan_normalized)) %>%
      dplyr::ungroup() %>%
      dplyr::select(lipidClass, compoundName, adductName, type, median_intensity, plate) %>%
      unique() %>%
      tidyr::pivot_wider(names_from = type, values_from = median_intensity)

    if (!"posctl" %in% colnames(whitelist_ctl_ms2_assessment)) {
      whitelist_ctl_ms2_assessment <- whitelist_ctl_ms2_assessment %>%
        dplyr::mutate(posctl = NA)
    }

    if (!"negctl" %in% colnames(whitelist_ctl_ms2_assessment)) {
      whitelist_ctl_ms2_assessment <- whitelist_ctl_ms2_assessment %>%
        dplyr::mutate(negctl = NA)
    }

    whitelist_ctl_ms2_assessment <- whitelist_ctl_ms2_assessment %>%
      dplyr::mutate(is_pass_ratio = case_when(
        grepl("_IS", compoundName) ~ TRUE, # pass internal standards through
        is.na(negctl) & !is.na(posctl) ~ TRUE, # missing neg ctl values are treated like 0
        is.na(posctl) ~ FALSE, # missing posctl values should always fail
        posctl / negctl >= bulkpool_params$bulkpoolMinPosNegCtlRatio ~ TRUE, # ratio is sufficiently high
        TRUE == TRUE ~ FALSE # measurements in both posctl and negctl, ratio not sufficiently high
      )) %>%
      dplyr::mutate(quant_class = "ms2") %>%
      dplyr::select(lipidClass, compoundName, adductName, plate, quant_class, negctl, posctl, is_pass_ratio)
  }

  whitelist_bulkpool_ms2_compounds <- whitelist_bulkpool_ms2_search$search %>%
    dplyr::inner_join(all_bulkpool_samples, by = c("sample" = "sample_name")) %>%
    dplyr::rename(sample_name = sample) %>%
    dplyr::select(file, sample_name, method, mode, plate, well_position, lipidClass, compoundName, adductName, ms1_is_13C_precursor) %>%
    dplyr::rename(is_13C_precursor = ms1_is_13C_precursor) %>%
    unique() %>%
    dplyr::mutate(is_ms2_compound = TRUE, is_ms3_compound = FALSE) %>%
    dplyr::arrange(sample_name, compoundName, adductName)

  whitelist_bulkpool_ms3_compounds <- tibble::tibble(
    file = character(0),
    sample_name = character(0),
    method = character(0),
    mode = character(0),
    plate = character(0),
    well_position = character(0),
    lipidClass = character(0),
    compoundName = character(0),
    adductName = character(0),
    is_ms2_compound = logical(0),
    is_ms3_compound = logical(0)
  )

  whitelist_ctl_ms3_assessment <- tibble::tibble(
    lipidClass = character(0),
    compoundName = character(0),
    adductName = character(0),
    plate = character(0),
    negctl = numeric(0),
    posctl = numeric(0),
    is_pass_ratio = logical(0)
  )

  # whitelist ms3 search ###################################
  if (is_ms3) {
    good_ms3_samples <- bulkpool_results_table %>% dplyr::filter(is_good_ms3)

    habc_tg_lib <- mzkitcpp::import_msp_lipids_library(file.path(lib_dir, ms3_lib_name))

    habc_tg_lib_all_targets <- habc_tg_lib %>%
      to_ms3_lib() %>%
      dplyr::select(prec_mzs) %>%
      unique()

    ms3_targets <- dplyr::inner_join(habc_tg_lib_all_targets, samples_info$ms3_targets, by = c("prec_mzs"))

    habc_tg_lib_filtered <- habc_tg_lib %>%
      clamr::to_ms3_lib() %>%
      dplyr::filter(prec_mzs %in% ms3_targets$prec_mzs) %>%
      clamr::to_ms2_lib()

    whitelist_bulkpool_ms3_search <- mzkitcpp::DI_pipeline_ms3_search(
      samples = good_ms3_samples$file,
      is_lib = clamr::default_tg_is_ms3,
      is_search_params = whitelist_ms3_search_params,
      search_lib = habc_tg_lib_filtered,
      search_params = whitelist_ms3_search_params,
      adducts_file = adducts_file,
      debug = F
    )

    whitelist_bulkpool_ms3_compounds <- whitelist_bulkpool_ms3_search %>%
      dplyr::inner_join(all_bulkpool_samples, by = c("sample" = "sample_name")) %>%
      dplyr::rename(sample_name = sample) %>%
      dplyr::select(file, sample_name, method, mode, plate, well_position, lipidClass, compoundName, adductName) %>%
      unique() %>%
      dplyr::mutate(is_13C_precursor = FALSE, is_ms2_compound = FALSE, is_ms3_compound = TRUE) %>%
      dplyr::arrange(sample_name, compoundName, adductName)

    posctl_ms3_samples <- posctl_results_table %>% dplyr::filter(is_good_ms3)

    if (nrow(posctl_ms3_samples) > 0) {
      whitelist_ctl_ms3_search <- mzkitcpp::DI_pipeline_ms3_search(
        samples = ctl_samples,
        is_lib = clamr::default_tg_is_ms3,
        is_search_params = whitelist_ms3_search_params,
        search_lib = habc_tg_lib_filtered,
        search_params = whitelist_ms3_search_params,
        adducts_file = adducts_file,
        debug = F
      )

      whitelist_ctl_ms3_assessment <- whitelist_ctl_ms3_search %>%
        dplyr::select(sample, lipidClass, compoundName, adductName, ms3_intensity_sum_norm) %>%
        unique() %>%
        dplyr::inner_join(all_ctl_samples, by = c("sample" = "sample_name")) %>%
        dplyr::group_by(lipidClass, compoundName, adductName, type) %>%
        dplyr::mutate(median_intensity = median(ms3_intensity_sum_norm)) %>%
        dplyr::ungroup() %>%
        dplyr::select(lipidClass, compoundName, adductName, type, median_intensity, plate) %>%
        unique() %>%
        tidyr::pivot_wider(names_from = type, values_from = median_intensity)

      if (!"posctl" %in% colnames(whitelist_ctl_ms3_assessment)) {
        whitelist_ctl_ms3_assessment <- whitelist_ctl_ms3_assessment %>%
          dplyr::mutate(posctl = NA)
      }

      if (!"negctl" %in% colnames(whitelist_ctl_ms3_assessment)) {
        whitelist_ctl_ms3_assessment <- whitelist_ctl_ms3_assessment %>%
          dplyr::mutate(negctl = NA)
      }

      whitelist_ctl_ms3_assessment <- whitelist_ctl_ms3_assessment %>%
        dplyr::mutate(is_pass_ratio = case_when(
          grepl("_IS", compoundName) ~ TRUE, # pass internal standards through
          is.na(negctl) & !is.na(posctl) ~ TRUE, # missing neg ctl values are treated like 0
          is.na(posctl) ~ FALSE, # missing posctl values should always fail
          posctl / negctl >= bulkpool_params$bulkpoolMinPosNegCtlRatio ~ TRUE, # ratio is sufficiently high
          TRUE == TRUE ~ FALSE # measurements in both posctl and negctl, ratio not sufficiently high
        )) %>%
        dplyr::mutate(quant_class = "ms3") %>%
        dplyr::select(lipidClass, compoundName, adductName, plate, quant_class, negctl, posctl, is_pass_ratio)
    }
  }

  whitelist_bulkpool_compounds <- rbind(whitelist_bulkpool_ms2_compounds, whitelist_bulkpool_ms3_compounds)
  whitelist_ctl_assessment <- rbind(whitelist_ctl_ms2_assessment, whitelist_ctl_ms3_assessment)

  # save results ###################################
  saveRDS(bulkpool_results_table, samples_output_file)
  saveRDS(whitelist_bulkpool_compounds, compounds_output_file)
  saveRDS(whitelist_ctl_assessment, ctl_output_file)

  invisible(0)
}

#' HABC v4 stage 2: compound whitelist and plate validation
#'
#' @description
#' Examine bulkpool search results from all plates, determine valid plates and valid list of compounds
#' to be searched against biological samples.
#'
#' @param stage_1_results_dir information about all bulkpool samples, and compounds identified therein
#' Two kinds of files are created in this directory: "samples" and "compounds".
#' Files are named as <plate>_<method>_<samples|compounds>.rds
#' Each file is a tibble, with information from the DIMS searches that can be later be
#' used to create a bulkpool library (in stage 2)
#'
#' @param bulkpool_params
#'        Parameters specific to handling of bulkpool samples.
#'        Does not control how specific lipids are identified in bulkpool samples.
#'
#' @param plate_validation_file
#'        .rds File containing which plates are valid to use for subsequent pipeline stages
#'        tibble with 4 columns: plate_name, mode, is_ms2_good, is_ms3_good
#'
#' @param bulkpool_sample_validation_file
#'        .rds file containing "good" designation for all bulkpool samples
#'        saved as tibble
#'
#' @param compound_list_file
#'        .rds File containing list of compound names to retain in regular searches
#'
#' @export
habc4_stage2_build_library <- function(
  stage_1_results_dir,
  bulkpool_params,
  plate_validation_file,
  bulkpool_sample_validation_file,
  compound_list_file
) {
  # plate validation ###################################
  bulkpool_samples_files <- list.files(stage_1_results_dir, pattern = ".*_samples.rds$", full.names = TRUE)

  bulkpool_sample_data <- NULL
  for (i in 1:length(bulkpool_samples_files)) {
    plate_sample_data <- readRDS(bulkpool_samples_files[i])
    if (is.null(bulkpool_sample_data)) {
      bulkpool_sample_data <- plate_sample_data
    } else {
      bulkpool_sample_data <- rbind(bulkpool_sample_data, plate_sample_data)
    }
  }

  plate_validation_tibble <- bulkpool_sample_data %>%
    dplyr::group_by(plate, mode) %>%
    dplyr::mutate(
      num_good_ms2 = sum(is_good_ms2),
      num_good_ms3 = sum(is_good_ms3)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      is_plate_good_ms2 = num_good_ms2 >= bulkpool_params$bulkpoolMinGoodPerPlate,
      is_plate_good_ms3 = num_good_ms3 >= bulkpool_params$bulkpoolMinGoodPerPlate
    ) %>%
    dplyr::select(plate, mode, is_plate_good_ms2, is_plate_good_ms3) %>%
    unique()

  plate_counts <- plate_validation_tibble %>%
    dplyr::group_by(mode) %>%
    dplyr::mutate(
      num_good_ms2_plates = sum(is_plate_good_ms2),
      num_good_ms3_plates = sum(is_plate_good_ms3)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(mode, num_good_ms2_plates, num_good_ms3_plates) %>%
    unique()

  pos_mode_plate_counts <- plate_counts %>% dplyr::filter(mode == "pos")
  neg_mode_plate_counts <- plate_counts %>% dplyr::filter(mode == "neg")

  num_plates_pos_ms2 <- ifelse(nrow(pos_mode_plate_counts) == 1, pos_mode_plate_counts$num_good_ms2_plates[1], 0)
  num_plates_pos_ms3 <- ifelse(nrow(pos_mode_plate_counts) == 1, pos_mode_plate_counts$num_good_ms3_plates[1], 0)
  num_plates_neg_ms2 <- ifelse(nrow(neg_mode_plate_counts) == 1, neg_mode_plate_counts$num_good_ms2_plates[1], 0)

  mode_good_counts <- dplyr::inner_join(bulkpool_sample_data, plate_validation_tibble, by = c("plate", "mode")) %>%
    dplyr::group_by(mode) %>%
    dplyr::mutate(
      num_good_ms2_bulkpool = sum(is_plate_good_ms2 & is_good_ms2),
      num_good_ms3_bulkpool = sum(is_plate_good_ms3 & is_good_ms3)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(mode, num_good_ms2_bulkpool, num_good_ms3_bulkpool) %>%
    unique() %>%
    dplyr::arrange(mode)

  pos_mode_counts <- mode_good_counts %>% dplyr::filter(mode == "pos")
  neg_mode_counts <- mode_good_counts %>% dplyr::filter(mode == "neg")

  num_valid_pos_ms2 <- ifelse(nrow(pos_mode_counts) == 1, pos_mode_counts$num_good_ms2_bulkpool[1], 0)
  num_valid_pos_ms3 <- ifelse(nrow(pos_mode_counts) == 1, pos_mode_counts$num_good_ms3_bulkpool[1], 0)
  num_valid_neg_ms2 <- ifelse(nrow(neg_mode_counts) == 1, neg_mode_counts$num_good_ms2_bulkpool[1], 0)

  # build bulkpool and ctl compound lists ###################################
  bulkpool_compound_files <- list.files(stage_1_results_dir, pattern = ".*_compounds.rds$", full.names = TRUE)

  bulkpool_compound_data <- NULL
  for (i in 1:length(bulkpool_compound_files)) {
    plate_compound_data <- readRDS(bulkpool_compound_files[i])
    if (is.null(bulkpool_compound_data)) {
      bulkpool_compound_data <- plate_compound_data
    } else {
      bulkpool_compound_data <- rbind(bulkpool_compound_data, plate_compound_data)
    }
  }

  ctl_compound_files <- list.files(stage_1_results_dir, pattern = ".*_ctl.rds$", full.names = TRUE)

  # Issue 734: when no compounds are identified from the plate_ctl_compound_data,
  # assume that there were not good pos_ctl samples for that plate.

  num_plates_pos_ms2_ctl <- 0
  num_plates_neg_ms2_ctl <- 0
  num_plates_pos_ms3_ctl <- 0

  ctl_compound_data <- NULL

  for (i in 1:length(ctl_compound_files)) {
    ctl_compound_file <- ctl_compound_files[i]

    plate_ctl_compound_data <- readRDS(ctl_compound_file)

    if (nrow(plate_ctl_compound_data) > 0) {
      if (grepl("X0158_M015A", ctl_compound_file)) {
        num_plates_neg_ms2_ctl <- num_plates_neg_ms2_ctl + 1
      } else if (grepl("X0158_M014A", ctl_compound_file)) {
        plate_ctl_tgs <- plate_ctl_compound_data %>% dplyr::filter(lipidClass == "TG")
        plate_ctl_non_tgs <- plate_ctl_compound_data %>% dplyr::filter(lipidClass != "TG")

        if (nrow(plate_ctl_tgs) > 0) {
          num_plates_pos_ms3_ctl <- num_plates_pos_ms3_ctl + 1
        }

        if (nrow(plate_ctl_non_tgs) > 0) {
          num_plates_pos_ms2_ctl <- num_plates_pos_ms2_ctl + 1
        }
      }
    }

    if (is.null(ctl_compound_data)) {
      ctl_compound_data <- plate_ctl_compound_data
    } else {
      ctl_compound_data <- rbind(ctl_compound_data, plate_ctl_compound_data)
    }
  }

  # filtering ###################################

  all_compound_ions <- rbind(bulkpool_params$ms1IonList, bulkpool_params$ms3IonList)

  # add ctl data to bulkpool compound list
  ctl_ratio_filtered_compound_ions <- ctl_compound_data %>%
    dplyr::mutate(mode = ifelse(grepl("\\+$", adductName), "pos", "neg")) %>%
    dplyr::filter(is_pass_ratio == TRUE) %>%
    dplyr::mutate(
      is_ms2_compound = ifelse(quant_class == "ms3", FALSE, TRUE),
      is_ms3_compound = ifelse(quant_class == "ms3", TRUE, FALSE)
    ) %>%
    dplyr::group_by(mode, lipidClass, compoundName, adductName) %>%
    dplyr::mutate(num_plates = n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(lipidClass, compoundName, adductName, mode, is_ms2_compound, is_ms3_compound, num_plates) %>%
    unique() %>%
    dplyr::mutate(ctl_ratio = case_when(
      is_ms2_compound & mode == "pos" ~ num_plates / num_plates_pos_ms2_ctl,
      is_ms3_compound & mode == "pos" ~ num_plates / num_plates_pos_ms3_ctl,
      is_ms2_compound & mode == "neg" ~ num_plates / num_plates_neg_ms2_ctl,
      TRUE == TRUE ~ 0
    )) %>%
    dplyr::inner_join(all_compound_ions, by = c("lipidClass", "adductName")) %>%
    dplyr::mutate(is_apply_ctl_ratio = !grepl("13C", compoundName) & !grepl("IS", compoundName) & ms1IsMajorIon) %>%
    dplyr::filter(!is_apply_ctl_ratio | (is_apply_ctl_ratio & ctl_ratio >= bulkpool_params$bulkpoolCtlRatioMinPlateFraction)) %>%
    dplyr::select(compoundName, adductName)

  # Issue 737: Compound ions of the mode that are not associated with the quant mode can sometimes be
  # removed from consideration at an early stage (the "ONLY_QUANT_MODE" option)
  #
  # The old default behavior beforebulkpool_params$whitelistCompoundIonModePolicy was created is EITHER_MODE.
  #
  # Remove all ions from all species that are not associated with the quant mode
  if (bulkpool_params$whitelistCompoundIonModePolicy == "ONLY_QUANT_MODE") {
    quant_mode <- rbind(bulkpool_params$ms1IonList, bulkpool_params$ms3IonList) %>%
      dplyr::filter(ms1IsPreferredQuantIon == TRUE) %>%
      dplyr::mutate(quant_mode = ifelse(grepl("\\-$", adductName), "neg", "pos")) %>%
      dplyr::select(lipidClass, quant_mode) %>%
      unique()

    bulkpool_library_adducts <- rbind(bulkpool_params$ms1IonList, bulkpool_params$ms3IonList) %>%
      dplyr::inner_join(quant_mode, by = c("lipidClass")) %>%
      dplyr::filter((grepl("\\+$", adductName) & quant_mode == "pos") | (grepl("\\-$", adductName) & quant_mode == "neg")) %>%
      dplyr::select(lipidClass, adductName)

    bulkpool_compound_data_quant_mode_filtered <- bulkpool_compound_data %>%
      dplyr::inner_join(bulkpool_library_adducts, by = c("lipidClass", "adductName"))
  } else {
    bulkpool_compound_data_quant_mode_filtered <- bulkpool_compound_data
  }

  # filter out ctl-ratio compounds and invalid plates
  bulkpool_compound_data_plate_filtered <- bulkpool_compound_data_quant_mode_filtered %>%
    dplyr::inner_join(ctl_ratio_filtered_compound_ions, by = c("compoundName", "adductName")) %>%
    dplyr::inner_join(., plate_validation_tibble, by = c("plate", "mode")) %>%
    dplyr::filter((is_ms2_compound & is_plate_good_ms2) | (is_ms3_compound & is_plate_good_ms3)) %>%
    dplyr::select(-is_plate_good_ms2, -is_plate_good_ms3, -file)

  # filter based on bulkpool parameters
  bulkpool_adduct_filtered <- bulkpool_compound_data_plate_filtered %>%
    dplyr::inner_join(., all_compound_ions, by = c("lipidClass", "adductName")) %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(num_samples = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(presence_fraction = case_when(
      is_ms2_compound & mode == "neg" ~ num_samples / num_valid_neg_ms2,
      is_ms2_compound & mode == "pos" ~ num_samples / num_valid_pos_ms2,
      is_ms3_compound & mode == "pos" ~ num_samples / num_valid_pos_ms3,
      TRUE == TRUE ~ -1
    )) %>%
    # Issue 727: isotopic designations are stricter because of the possibility of precursor [M+1] fallback.
    dplyr::mutate(is_m_plus_zero = !grepl("_13C$", compoundName) & is_13C_precursor == FALSE) %>%
    dplyr::mutate(is_m_plus_one = !grepl("_13C_13C$", compoundName) & grepl("_13C$", compoundName) & is_13C_precursor == FALSE) %>%
    dplyr::mutate(is_m_plus_two = grepl("_13C_13C$", compoundName) & is_13C_precursor == FALSE) %>%
    dplyr::filter((ms1IsMajorIon == TRUE & presence_fraction >= bulkpool_params$bulkpoolMajorAdductMinFraction) |
      (ms1IsMajorIon == FALSE & presence_fraction >= bulkpool_params$bulkpoolMinorAdductMinFraction))

  # MS3 compounds have no isotope or multiple adduct constraints
  bulkpool_ms3_compounds <- bulkpool_adduct_filtered %>%
    dplyr::filter(is_ms3_compound == TRUE) %>%
    dplyr::select(compoundName) %>%
    unique()

  # Filter out compounds where _13C isotopes are required and not found
  isotope_constraints <- bulkpool_params$ms1IsRequireMPlusOneByLipidClassAndAdduct %>%
    dplyr::inner_join(bulkpool_params$ms1IsRequireMPlusTwoByLipidClassAndAdduct, by = c("lipidClass", "adductName"))

  # Issue 727: requiring a non-[M+0] isotope automatically adds requirement of also finding a true [M+0]
  bulkpool_isotope_filtered <- bulkpool_adduct_filtered %>%
    dplyr::filter(is_ms2_compound == TRUE) %>%
    dplyr::select(lipidClass, compoundName, adductName, is_m_plus_zero, is_m_plus_one, is_m_plus_two) %>%
    unique() %>%
    dplyr::mutate(monoCompoundName = gsub("_13C", "", compoundName)) %>%
    dplyr::group_by(monoCompoundName, adductName) %>%
    dplyr::mutate(
      is_has_m_plus_zero = any(is_m_plus_zero),
      is_has_m_plus_one = any(is_m_plus_one),
      is_has_m_plus_two = any(is_m_plus_two)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(isotope_constraints, by = c("lipidClass", "adductName")) %>%
    dplyr::select(
      lipidClass, monoCompoundName, adductName,
      is_has_m_plus_zero, is_has_m_plus_one, is_has_m_plus_two, ms1IsRequireMPlusOne, ms1IsRequireMPlusTwo
    ) %>%
    unique() %>%
    dplyr::rename(compoundName = monoCompoundName) %>%
    dplyr::filter(!ms1IsRequireMPlusOne | (ms1IsRequireMPlusOne & is_has_m_plus_one & is_has_m_plus_zero)) %>%
    dplyr::filter(!ms1IsRequireMPlusTwo | (ms1IsRequireMPlusTwo & is_has_m_plus_two & is_has_m_plus_zero)) %>%
    dplyr::select(lipidClass, compoundName, adductName) %>%
    unique()

  # Issues 737: filter out compounds where all required adduct forms are not found.
  # Depending on compound ion policy, This assessment is made within each mode separately
  # or across both modes simultaneously
  #
  # The old default behavior before bulkpool_params$whitelistCompoundIonModePolicy was created is EITHER_MODE.
  if (bulkpool_params$whitelistCompoundIonModePolicy == "EITHER_MODE" || bulkpool_params$whitelistCompoundIonModePolicy == "COMBINE_QUANT_MODE") {
    num_adducts_constraints <- bulkpool_params$ms1IonList %>%
      dplyr::mutate(mode = ifelse(grepl("\\+$", adductName), "pos", "neg")) %>%
      dplyr::group_by(lipidClass, mode) %>%
      dplyr::mutate(min_num_adducts = n()) %>%
      ungroup() %>%
      dplyr::select(lipidClass, mode, min_num_adducts) %>%
      unique()

    bulkpool_ms2_compound_ions_filtered <- bulkpool_isotope_filtered %>%
      dplyr::mutate(mode = ifelse(grepl("\\+$", adductName), "pos", "neg")) %>%
      dplyr::inner_join(num_adducts_constraints, by = c("lipidClass", "mode")) %>%
      dplyr::group_by(compoundName, mode) %>%
      dplyr::mutate(num_adducts = n()) %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter(num_adducts >= min_num_adducts) %>%
      dplyr::select(compoundName, adductName) %>%
      unique()
  } else {
    num_adducts_constraints <- bulkpool_params$ms1IonList %>%
      dplyr::group_by(lipidClass) %>%
      dplyr::mutate(min_num_adducts = n()) %>%
      ungroup() %>%
      dplyr::select(lipidClass, min_num_adducts) %>%
      unique()

    bulkpool_ms2_compound_ions_filtered <- bulkpool_isotope_filtered %>%
      dplyr::inner_join(num_adducts_constraints, by = c("lipidClass")) %>%
      dplyr::group_by(compoundName) %>%
      dplyr::mutate(num_adducts = n()) %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter(num_adducts >= min_num_adducts) %>%
      dplyr::select(compoundName, adductName) %>%
      unique()
  }

  # Add _M1 fallback precursor compounds ###################################
  bulkpool_compound_list_filtered <- bulkpool_compound_data_plate_filtered %>%
    dplyr::inner_join(bulkpool_ms2_compound_ions_filtered, by = c("compoundName", "adductName"))

  # add _M1 additions to the list, if the [M+1] fallback had to be triggered often enough.
  M1_assessments <- habc4_compound_counts(bulkpool_compound_list_filtered) %>%
    dplyr::filter(!grepl("_IS", compoundName) & !grepl("_13C", compoundName)) %>%
    dplyr::mutate(m_one_fraction = num_samples_m_plus_one / num_samples) %>%
    dplyr::filter(m_one_fraction >= bulkpool_params$bulkpoolM1MinFraction)

  M1_ions <- M1_assessments %>% dplyr::select(compoundName, adductName)

  ms2_compound_ions_no_M1 <- bulkpool_ms2_compound_ions_filtered %>%
    dplyr::anti_join(M1_ions, by = c("compoundName", "adductName")) %>%
    dplyr::mutate(adductName = "")

  ms3_compound_ions <- tibble::tibble(
    compoundName = bulkpool_ms3_compounds$compoundName,
    adductName = rep("", length(bulkpool_ms3_compounds$compoundName))
  )

  bulkpool_compounds <- rbind(ms2_compound_ions_no_M1, M1_ions, ms3_compound_ions) %>%
    dplyr::group_by(compoundName) %>%
    dplyr::mutate(M1_adducts = paste0(adductName, collapse = " ")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-adductName) %>%
    unique()

  # save results ###################################
  saveRDS(plate_validation_tibble, plate_validation_file)
  saveRDS(bulkpool_sample_data, bulkpool_sample_validation_file)
  saveRDS(bulkpool_compounds, compound_list_file)

  invisible(0)
}

#' HABC4 stage 3: bulkpool samples analysis
#'
#' @description
#' Stage 3 of pipeline - search biological samples using previously determined subset
#' of whitelisted compounds.
#'
#' @param plate_name
#'        name of folder containing mzML files associated with a single HABC plate
#'
#' @param plate_validation_data
#'        Table of plate validation results, indicating which plates are invalid and should not be analyzed.
#'
#' @param samples_file_path
#'        directory of HABC mzML files
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing IS info for nornmalization  (lib_dir/is_lib_name)
#'
#' @param habc_lib_name
#'        HABC library name (lib_dir/habc_lib_name)
#'
#' @param ms3_lib_name
#'        MS3 library name (lib_dir/ms3_lib_name)
#'
#' @param is_search_params
#'        search parameters for internal standard (IS) search.
#'
#' @param bulkpool_compounds
#'        List of all compounds cleared for search based on stage 1 and 2 analysis.
#'
#' @param biological_ms2_search_params
#'        Parameters associated with identification of lipids in biological samples.
#'
#' @param biological_ms3_search_params
#'        Parameters associated with identification of lipids from targeted MS3 scans in biological samples.
#'
#' @param stage_3_rds_dir
#'        Results of search, saved as table of compounds with centered quant measurements.
#'
#' @param stage3_mzrolldb_dir
#'        Plate-specific results saved into mzrollDB files
#'
#' @param save_mzrolldb_as_rds
#'        If true, save results as mzrollDB RDS, to be converted into mzrollDB later.
#'
#' @export
habc4_stage3_process_plate <- function(
  plate_name,
  plate_validation_data,
  samples_file_path,
  bulkpool_sample_data,
  lib_dir,
  is_lib_name,
  habc_lib_name,
  ms3_lib_name,
  bulkpool_compounds,
  is_search_params,
  biological_ms2_search_params,
  biological_ms3_search_params,
  stage_3_rds_dir,
  stage_3_mzrolldb_dir,
  save_mzrolldb_as_rds = FALSE
) {
  # files ######################################
  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  habc_lib_file <- paste(lib_dir, habc_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")
  is_ms3 <- !is.null(ms3_lib_name)
  method_name <- ifelse(is_ms3, "X0158_M014A", "X0158_M015A")
  plate_mode <- ifelse(is_ms3, "pos", "neg")

  # contains MS3, single measurement of post-partitioned MS1/MS2 compound ion with description
  rds_results_file <- paste0(stage_3_rds_dir, "/", method_name, "_", plate_name, "_rows.rds")

  # Only MS1/MS2 compound ion. Contains MS1 and MS2 partitioned compound intenstiies,
  # along with ms1_intensity_is_nearest_scan_normalized and diagnostic_ms2_intensity
  col_quant_results_file <- paste0(stage_3_rds_dir, "/", method_name, "_", plate_name, "_cols.rds")

  suffix <- ifelse(save_mzrolldb_as_rds == TRUE, ".rds", ".mzrollDB")
  mzrolldb_results_file <- paste0(stage_3_mzrolldb_dir, "/", method_name, "_", plate_name, suffix)

  # plate validation ############################
  plate_validation_data_plate <- plate_validation_data %>%
    dplyr::filter(plate == plate_name & mode == plate_mode)

  is_plate_valid_ms2 <- plate_validation_data_plate$is_plate_good_ms2[1]
  is_plate_valid_ms3 <- plate_validation_data_plate$is_plate_good_ms3[1]

  if (!is_plate_valid_ms2 && !is_plate_valid_ms3) {
    cat(paste0("plate ", plate_name, " ", plate_mode, " has insufficient valid ms2 or ms3 bulkpool samples.\n"))
    return(invisible(0))
  }

  # samples ####################################
  samples_info <- habc_samples_df(samples_file_path, is_ms3)
  samples_df <- samples_info$samples_df %>% dplyr::arrange(desc(type), order_num)
  habc4_ms2_ranges <- samples_info$ms2_ranges

  # good bulkpool samples ######################
  good_ms2_bulkpools <- bulkpool_sample_data %>% dplyr::filter(plate == plate_name, mode == plate_mode, is_good_ms2 == TRUE)
  good_ms3_bulkpools <- bulkpool_sample_data %>% dplyr::filter(plate == plate_name, mode == plate_mode, is_good_ms3 == TRUE)

  samples_df_ms2_filtered <- samples_df %>%
    dplyr::filter((type == "BulkPool" & sample_name %in% good_ms2_bulkpools$sample_name) | type != "BulkPool")

  samples_df_ms3_filtered <- samples_df %>%
    dplyr::filter((type == "BulkPool" & sample_name %in% good_ms3_bulkpools$sample_name) | type != "BulkPool")

  # MS2 searches ###############################
  if (is_plate_valid_ms2) {
    # IS library
    lib_subset_string <- ifelse(is_ms3, "\\+$", "\\-$")
    is_lib <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl(lib_subset_string, adductName))
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

    centered_nearest_scan_IS <- habc_bulkpool_centered_quant_ion_table(nearest_scan_IS, "ms1_intensity_is_nearest_scan_normalized") %>%
      dplyr::select(lipidClass, compositionSummary, compoundName, adductName, sample, ms1_intensity_is_nearest_scan_normalized) %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized))

    diagnostic_ms2 <- habc4_biological_search$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity) &
        ms2_diagnostic_norm_intensity > 0 &
        ms2_diagnostic_norm_intensity < Inf)

    centered_diagnostic_ms2 <- habc_bulkpool_centered_quant_ion_table(diagnostic_ms2, "ms2_diagnostic_norm_intensity") %>%
      dplyr::select(lipidClass, compositionSummary, compoundName, adductName, sample, ms2_diagnostic_norm_intensity) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity))

    # center after determination of preferred quant type to ensure accurate centering values
    habc4_quant_table <- to_quant_table(habc4_biological_search$adduct_table, FALSE, TRUE)

    # MS1
    habc4_quant_ions_ms1 <- habc4_quant_table %>%
      dplyr::select(-ms2_intensity, -ms2_quant_type) %>%
      dplyr::filter(!is.na(ms1_intensity) & ms1_intensity > 0 & ms1_intensity < Inf) %>%
      dplyr::mutate(ms1_quant_type = ifelse(grepl("^ms1", ms1_quant_type), ms1_quant_type, paste0("ms1_", ms1_quant_type))) %>%
      dplyr::rename(quant_type = ms1_quant_type, intensity = ms1_intensity) %>%
      dplyr::mutate(quant_class = "ms1")

    # MS2
    habc4_quant_ions_ms2 <- habc4_quant_table %>%
      dplyr::select(-ms1_intensity, -ms1_quant_type) %>%
      dplyr::filter(!is.na(ms2_intensity) & ms2_intensity > 0 & ms2_intensity < Inf) %>%
      dplyr::mutate(ms2_quant_type = ifelse(grepl("^ms2", ms2_quant_type), ms2_quant_type, paste0("ms2_", ms2_quant_type))) %>%
      dplyr::rename(quant_type = ms2_quant_type, intensity = ms2_intensity) %>%
      dplyr::mutate(quant_class = "ms2")

    # Center quant ions
    habc4_centered_quant_ions_ms1 <- habc_bulkpool_centered_quant_ion_table(habc4_quant_ions_ms1, "intensity")
    habc4_centered_quant_ions_ms2 <- habc_bulkpool_centered_quant_ion_table(habc4_quant_ions_ms2, "intensity")

    di_ms1_quant_ions_condensed <- habc4_centered_quant_ions_ms1 %>%
      dplyr::rename(ms1_intensity = intensity) %>%
      dplyr::select(sample, lipidClass, compositionSummary, compoundName, adductName, ms1_intensity) %>%
      dplyr::filter(!is.na(ms1_intensity))

    di_ms2_quant_ions_condensed <- habc4_centered_quant_ions_ms2 %>%
      dplyr::rename(ms2_intensity = intensity) %>%
      dplyr::select(sample, lipidClass, compositionSummary, compoundName, adductName, ms2_intensity) %>%
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
  if (is_plate_valid_ms3) {
    # MS3 library
    habc_ms3_lib <- mzkitcpp::import_msp_lipids_library(file.path(lib_dir, ms3_lib_name))

    habc_ms3_lib_all_targets <- habc_ms3_lib %>%
      to_ms3_lib() %>%
      dplyr::select(prec_mzs) %>%
      unique()

    ms3_targets <- dplyr::inner_join(habc_ms3_lib_all_targets, samples_info$ms3_targets, by = c("prec_mzs"))

    habc_ms3_lib_filtered <- habc_ms3_lib %>%
      clamr::to_ms3_lib() %>%
      dplyr::filter(prec_mzs %in% ms3_targets$prec_mzs) %>%
      clamr::to_ms2_lib() %>%
      dplyr::filter(compoundName %in% bulkpool_compounds$compoundName) %>%
      dplyr::mutate(ms2_intensity = 1) %>%
      dplyr::select(colnames(habc_ms3_lib))

    # MS3 biological search
    habc4_biological_ms3_search <- mzkitcpp::DI_pipeline_ms3_search(
      samples = samples_df_ms3_filtered$file,
      is_lib = clamr::default_tg_is_ms3,
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
        intensity, quant_type, quant_class,
        method, mode, plate, order_num, type, barcode, well_position, injection_num
      ) %>%
      unique()
  } else {
    habc4_ms3_quant_table <- NULL
  }

  # color table ###################################
  color_table <- habc_type_color_table(samples_df, good_ms2_bulkpools$sample_name)
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
      HABC_search_results <- habc4_biological_search$search
      HABC_adduct_table <- habc4_biological_search$adduct_table
      HABC_quant_ions <- habc4_quant_ions
    }

    MS3_search_results <- "NULL"
    if (!is.null(habc4_ms3_quant_table)) {
      MS3_search_results <- habc4_biological_ms3_search_centered
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
        "HABC_library_name" = habc_lib_name,
        "HABC_search_lib_sliced" = habc_lib_sliced,
        "HABC_search_params" = biological_ms2_search_params,
        "HABC_search_results" = HABC_search_results,
        "HABC_adduct_table" = HABC_adduct_table,
        "HABC_quant_ions" = HABC_quant_ions,
        "HABC_set_name" = plate_name_vector,
        "MS3_samples" = samples_df_ms3_filtered$file,
        "MS3_library_name" = ms3_lib_name,
        "MS3_search_params" = whitelist_ms3_search_params,
        "MS3_search_results" = MS3_search_results,
        "MS3_set_name" = plate_name_vector
      )

    saveRDS(mzrolldb_results, file = mzrolldb_results_file)
  } else {
    if (!is.null(habc4_IS_results)) {
      # Add IS search results
      clamr::add_direct_infusion_search_results(
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
      clamr::add_direct_infusion_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df_ms2_filtered$file,
        ms2_ranges = habc4_ms2_ranges,
        library_name = habc_lib_name,
        search_lib_sliced = habc_lib_sliced,
        search_params = biological_ms2_search_params,
        adducts_file = adducts_file,
        di_search_results = habc4_biological_search$search,
        di_quant_ions = habc4_quant_ions,
        set_name = plate_name_vector,
        color_table = color_table
      )
    }

    # Add MS3 search

    if (!is.null(habc4_ms3_quant_table)) {
      clamr::add_targeted_ms3_search_results(
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

  invisible(0)
}

#' HABC4: Return table of quantitative measurements from collection of results.
#'
#' @description
#' combine stage 3 RDS result files into final table of quant results.
#'
#' @param stage_3_rds_dir collection of saved RDS results files.
#' @param bulkpool_params collection of params associated with bulkpool search.  Used here for quant ions.
#'
#' @export
habc4_stage4_formatting <- function(
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
    # TESTING BLOCK
    # plate_quant_data <- readRDS(row_stage3_results[i])
    # col_plate_quant_data <- readRDS(col_stage3_results[i])

    plate_quant_data <- readRDS(row_stage3_results[i]) %>%
      dplyr::inner_join(quant_compound_ions, by = c("lipidClass", "adductName")) %>%
      dplyr::filter(!grepl("_13C", compoundName) & !grepl("_IS", compoundName))

    col_plate_quant_data <- readRDS(col_stage3_results[i]) %>%
      dplyr::filter(!grepl("_13C", compoundName) & !grepl("_IS", compoundName))

    plate_mode <- plate_quant_data$mode[1]
    plate_name <- plate_quant_data$plate[1]

    if (is.null(stage3_quant_ions)) {
      stage3_quant_ions <- plate_quant_data
      stage3_quant_ions_cols <- col_plate_quant_data
    } else {
      stage3_quant_ions <- rbind(stage3_quant_ions, plate_quant_data)
      stage3_quant_ions_cols <- rbind(stage3_quant_ions_cols, col_plate_quant_data)
    }
  }

  # Summarization ###################################

  stage3_quant_ions_to_summarize <- stage3_quant_ions %>%
    dplyr::filter(grepl("^\\{", compoundName) & type == "bc") %>%
    dplyr::select(compositionSummary, adductName) %>%
    unique()

  stage3_quant_ions_unaffected <- stage3_quant_ions %>%
    dplyr::anti_join(stage3_quant_ions_to_summarize, by = c("compositionSummary", "adductName"))

  stage3_quant_ions_to_alter <- stage3_quant_ions %>%
    dplyr::inner_join(stage3_quant_ions_to_summarize, by = c("compositionSummary", "adductName"))

  stage3_quant_ions_cols_altered_subset <- stage3_quant_ions_cols %>%
    dplyr::inner_join(stage3_quant_ions_to_summarize, by = c("compositionSummary", "adductName"))

  sample_info <- stage3_quant_ions %>%
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
      intensity, quant_type, quant_class, method, mode, plate, order_num, type, barcode, well_position, injection_num
    )

  stage3_quant_ions_summarized_corrected <- rbind(stage3_quant_ions_unaffected, stage3_quant_ions_altered) %>%
    dplyr::filter(!grepl("^\\{", compoundName)) # remove any controls or bulkpool samples with summarized compounds.

  # Median BulkPool CV ###############################
  plate_modes <- stage3_quant_ions_summarized_corrected %>%
    dplyr::select(plate, mode) %>%
    unique()

  stage3_quant_cvs <- NULL

  for (i in 1:nrow(plate_modes)) {
    ith_plate <- plate_modes$plate[i]
    ith_mode <- plate_modes$mode[i]

    plate_quant_data <- stage3_quant_ions_summarized_corrected %>%
      dplyr::filter(plate == ith_plate & mode == ith_mode)

    plate_quant_cvs <- habc_cv_comparison_v3(plate_quant_data) %>%
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

  # lowest CV designations ###############################
  stage3_median_bulkpool_CVs <- BulkPool_CVs %>%
    dplyr::select(compoundName, adductName, mode, quant_class, median_BulkPool_CV) %>%
    unique() %>%
    dplyr::arrange(compoundName, adductName, mode) %>%
    dplyr::group_by(compoundName) %>%
    dplyr::mutate(is_lowest_CV = median_BulkPool_CV == min(median_BulkPool_CV, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(compoundName, quant_class) %>%
    dplyr::mutate(is_lowest_MS_CV = median_BulkPool_CV == min(median_BulkPool_CV, na.rm = TRUE)) %>%
    dplyr::ungroup()

  quantified_compounds <- stage3_quant_ions_summarized_corrected %>%
    dplyr::inner_join(stage3_median_bulkpool_CVs, by = c("compoundName", "adductName", "mode", "quant_class"))
}
