load(file.path(here::here(), "data", "default_tg_is_ms3.rda"))

#' Print HABC Method info
#'
#' @description
#' Print HABC method names and mapping to positive and negative ionization mode.
#'
#' @export
habc_method_info <- function() {
  cat("HABC Method info:\n")
  cat("X0158_M014A (positive mode)\n")
  cat("X0158_M015A (negative mode)\n")
}

#' HABC Analysis Pipeline Version History
#'
#' @description
#' Detailed summary of version information
#'
#' @export
habc_version_history <- function() {
  cat("Version 1:
      \tDescription:
      \t\tInitial version of pipeline.
      \t\tUsed for initial analysis of 700 metadata-labled samples.
      \tCompletion Date:
      \t\tJanuary 11, 2021
      \tMass Spec Issues:
      \t\t#561 \"Port HABC process plates script into clamr\"
      \n")
  cat("Version 2:
      \tDescription:
      \t\tWhitelisting based on presence of common adduct form in bulkpool samples.
      \t\tApplies to negative mode only.
      \t\tWhitelisting parameters are strict.
      \t\tQuantitation is done on common adduct form.
      \t\tOrdinary search parameters are same as version 1 parameters.
      \t\t13C [M+1] isotopic forms are added based on whitelisting search results.
      \tCompletion Date:
      \t\tOctober 19, 2021
      \tMass Spec Issues:
      \t\t#700 \"HABC Search Version 2 (neg mode)\"
      \t\t#606 \"p58 Analysis\"
      \t\t#625 \"p58 Analysis 2\"
      \t\t#646 \"p58 Analysis 3 - Iterative DIMS Identification based on blacklisting\"
      \t\t#624 \"Utilize BulkPools to evaluate ms1_partition_fraction_SAF\"
      \n")
  cat("Version 3:
      \tDescription:
      \t\tWhitelisting/Bulkpool approach applied to positive mode data.
      \tCompletion Data:
      \t\tNot Yet Completed
      \tMass Spec Issues:
      \t\t#717 \"HABC Search Version 3 (positive mode)\"
      \t\t#716 \"Updates to positive mode lipids library based on DIMS IS data\"
      \t\t
      \n")

  invisible(0)
}

#' HABC Samples Information
#'
#' @description
#' Generate samples table containing sample names, filenames, and parsed metadata from the sample name.
#' Generate ms2 scans information and ms3 target lists based on samples content.
#' The ms2 and ms3 scans found in each sample are counted.
#' Samples that do not contain the maximal number of different ms2 m/z ranges and ms3 targets
#' are excluded.
#'
#' @param samples_file_path
#'       directory of HABC mzML files, or vector of HABC mzML files.
#'
#' @param is_has_ms3
#'        check for MS3 scans in mzML files, enumerate ms3 targets
#'
#' @return
#' \code{list("samples_df","ms2_ranges","ms3_targets")}
#'
#' \itemize{
#'    \item{\code{samples_df}}{: A table of samples with parsed metadata}
#'    \item{\code{ms2_ranges}}{: MS2 scan information (precursor m/z ranges)}
#'    \item{\code{ms3_targets}}{: vector of (MS1, MS2) precursor m/zs associated with MS3 scans.}
#' }
#'
#' @export
habc_samples_df <- function(samples_file_path, is_has_ms3 = FALSE) {
  if (length(samples_file_path) > 1) {
    # list of files - pass along, assuming full paths
    samples_files <- samples_file_path
  } else {
    # directory - take full paths of all mzML files in directory
    samples_files <- list.files(samples_file_path, pattern = ".*.mzML$", full.names = TRUE)
  }

  samples_df <- tibble::tibble(file = samples_files) %>%
    dplyr::mutate(sample_name = gsub("^.*\\/", "", file)) %>%
    # extract all components
    dplyr::mutate(components = stringr::str_split(sample_name, "_")) %>%
    dplyr::mutate(method = sapply(components, function(x) {
      x[2]
    })) %>%
    dplyr::mutate(mode = ifelse(method == "M014A", "pos", "neg")) %>%
    dplyr::mutate(plate = sapply(components, function(x) {
      x[3]
    })) %>%
    dplyr::mutate(order_num = as.numeric(sapply(components, function(x) {
      x[4]
    }))) %>%
    dplyr::mutate(type_id = sapply(components, function(x) {
      x[5]
    })) %>%
    dplyr::mutate(type = ifelse(grepl("bc[0-9]+", type_id), "bc", type_id)) %>%
    dplyr::mutate(barcode = stringr::str_extract(type_id, "(?<=bc)[0-9]+")) %>%
    dplyr::mutate(well_position = sapply(components, function(x) {
      x[6]
    })) %>%
    dplyr::mutate(injection_num = as.numeric(stringr::str_extract(sapply(components, function(x) {
      {
        x
      }[7]
    }), "[0-9]+(?=.mzML)"))) %>%
    # all components reformatted as new columns, remove list
    dplyr::select(-components) %>%
    # only retain latest injection when there are multiple injections
    dplyr::arrange(order_num, -injection_num) %>%
    dplyr::group_by(order_num) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  samples_df_scan_data <- tibble::tibble(
    file = samples_df$file,
    count_ms2_ranges = rep(0, nrow(samples_df)),
    count_ms3_targets = rep(0, nrow(samples_df))
  )

  best_ms2_ranges <- NULL
  best_ms3_targets <- NULL

  # Issue 580: Add extra guards, filters for ensuring that ms2_ranges and ms3_ranges are accurate.
  for (i in 1:nrow(samples_df_scan_data)) {
    ms2_ranges <- mzkitcpp::DI_ms2_ranges(samples_df_scan_data$file[i])

    if (is.null(best_ms2_ranges) || nrow(ms2_ranges) > nrow(best_ms2_ranges)) {
      best_ms2_ranges <- ms2_ranges
    }

    samples_df_scan_data$count_ms2_ranges[i] <- nrow(ms2_ranges)

    if (is_has_ms3) {
      ms3_targets <- mzkitcpp::DI_ms3_targets(samples_df_scan_data$file[i]) %>%
        dplyr::select(prec_mzs) %>%
        unique()

      num_ms3_ranges <- nrow(ms3_targets)

      if (is.null(best_ms3_targets) || nrow(ms3_targets) > nrow(best_ms3_targets)) {
        best_ms3_targets <- ms3_targets
      }

      samples_df_scan_data$count_ms3_targets[i] <- nrow(ms3_targets)
    }
  }

  num_best_ms3_targets <- 0
  if (!is.null(best_ms3_targets)) {
    num_best_ms3_targets <- nrow(best_ms3_targets)
  }

  samples_df_scan_data <- samples_df_scan_data %>%
    dplyr::filter(count_ms2_ranges == nrow(best_ms2_ranges) & count_ms3_targets == num_best_ms3_targets) %>%
    dplyr::filter(count_ms2_ranges > 0 & (count_ms3_targets > 0 | !is_has_ms3)) # case where no samples are valid

  samples_df <- samples_df %>%
    dplyr::filter(file %in% samples_df_scan_data$file)

  list(samples_df = samples_df, ms2_ranges = best_ms2_ranges, ms3_targets = best_ms3_targets)
}

#' Sample Type-based color table
#'
#' @description
#' Create a color table from a HABC samples data frame
#'
#' @param samples_df
#'     HABC samples df as generated by \code{habc_samples_df()}
#'
#' @param good_bulkpool_samples
#'     list of good bulkpool samples
#'
#' @export
habc_type_color_table <- function(samples_df, good_bulkpool_samples) {
  color_table <- tibble::tibble(
    sample_name = samples_df$sample_name,
    red = rep(0, length(samples_df$sample_name)),
    blue = rep(0, length(samples_df$sample_name)),
    green = rep(0, length(samples_df$sample_name))
  )

  for (i in 1:nrow(samples_df)) {
    type <- samples_df$type[i]
    sample <- samples_df$sample_name[i]

    if (type == "bc") {
      color_table$red[i] <- 0.9
      color_table$blue[i] <- 0.2
    } else if (type == "BulkPool") {
      if (sample %in% good_bulkpool_samples) {
        # good bulkpools - included in regular search
        color_table$blue[i] <- 0.8
        color_table$green[i] <- 0.5
      } else {
        # bad bulkpools - excluded from regular search
        color_table$red[i] <- 0.5
        color_table$blue[i] <- 0.5
        color_table$green[i] <- 0.5
      }
    } else if (type == "negctl" || type == "posctl") {
      color_table$green[i] <- 0.9
      color_table$blue[i] <- 0.2
    }
  }

  color_table
}

#' Bulkpool-centered adduct table
#'
#' @description
#' Center plate quant values based on median of quant measurements from BulkPool samples
#'
#' @param adduct_table
#'     Adduct table from HABC plate search
#'
#' @param samples_df
#'     HABC samples df as generated by \code{habc_samples_df()}
#'
#' @export
habc_bulkpool_centered_adduct_table <- function(adduct_table, samples_df) {
  adduct_table_ms1_median <- dplyr::left_join(adduct_table, samples_df, by = c("sample" = "sample_name")) %>%
    dplyr::filter(type == "BulkPool") %>%
    dplyr::group_by(compoundName, adductName, plate) %>%
    dplyr::mutate(ms1_intensity_median = median(ms1_intensity_is_nearest_scan_normalized, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(compoundName, adductName, ms1_intensity_median) %>%
    unique()

  adduct_table_ms2_medians <- dplyr::left_join(adduct_table, samples_df, by = c("sample" = "sample_name")) %>%
    dplyr::filter(type == "BulkPool") %>%
    dplyr::group_by(compoundName, adductName, plate) %>%
    dplyr::mutate(
      ms2_diagnostic_norm_intensity =
        ifelse(ms2_diagnostic_norm_intensity == -1, NA, ms2_diagnostic_norm_intensity),
      diagnostic_fragment_sum_is_normalized =
        ifelse(diagnostic_fragment_sum_is_normalized == -1, NA, diagnostic_fragment_sum_is_normalized),
      diagnostic_fragment_sum_SAF_is_normalized =
        ifelse(diagnostic_fragment_sum_SAF_is_normalized == -1, NA, diagnostic_fragment_sum_SAF_is_normalized),
      acyl_fragment_sum_is_normalized =
        ifelse(acyl_fragment_sum_is_normalized == -1, NA, acyl_fragment_sum_is_normalized),
      acyl_fragment_sum_SAF_is_normalized =
        ifelse(acyl_fragment_sum_SAF_is_normalized == -1, NA, acyl_fragment_sum_SAF_is_normalized)
    ) %>%
    dplyr::mutate(ms2_top_diagnostic_median = median(ms2_diagnostic_norm_intensity, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_diagnostic_sum_median = median(diagnostic_fragment_sum_is_normalized, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_diagnostic_sum_SAF_median = median(diagnostic_fragment_sum_SAF_is_normalized, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_acyl_sum_median = median(acyl_fragment_sum_is_normalized, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_acyl_sum_SAF_median = median(acyl_fragment_sum_SAF_is_normalized, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      compoundName, adductName,
      ms2_top_diagnostic_median, ms2_diagnostic_sum_median, ms2_acyl_sum_median, ms2_diagnostic_sum_SAF_median, ms2_acyl_sum_SAF_median
    ) %>%
    unique()

  adduct_table_centered <- adduct_table %>%
    dplyr::left_join(., adduct_table_ms1_median, by = c("compoundName", "adductName")) %>%
    dplyr::left_join(., adduct_table_ms2_medians, by = c("compoundName", "adductName")) %>%
    dplyr::mutate(
      ms1_intensity_is_nearest_scan_normalized = ms1_intensity_is_nearest_scan_normalized / ms1_intensity_median,
      diagnostic_fragment_sum_is_normalized =
        ifelse(diagnostic_fragment_sum_is_normalized == -1, NA, diagnostic_fragment_sum_is_normalized / ms2_diagnostic_sum_median),
      diagnostic_fragment_sum_SAF_is_normalized =
        ifelse(diagnostic_fragment_sum_SAF_is_normalized == -1, NA, diagnostic_fragment_sum_SAF_is_normalized / ms2_diagnostic_sum_SAF_median),
      acyl_fragment_sum_is_normalized =
        ifelse(acyl_fragment_sum_is_normalized == -1, NA, acyl_fragment_sum_is_normalized / ms2_acyl_sum_median),
      acyl_fragment_sum_SAF_is_normalized =
        ifelse(acyl_fragment_sum_SAF_is_normalized == -1, NA, acyl_fragment_sum_SAF_is_normalized / ms2_acyl_sum_SAF_median),
      ms2_diagnostic_norm_intensity =
        ifelse(ms2_diagnostic_norm_intensity == -1, NA, ms2_diagnostic_norm_intensity / ms2_top_diagnostic_median)
    )
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
habc_bulkpool_centered_quant_ion_table <- function(quant_ion_table, intensity_col_name) {
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
    dplyr::mutate(!!sym(intensity_col_name) := (!!sym(intensity_col_name)) / intensity_median) %>%
    dplyr::select(-intensity_median)
}

#' Bulkpool-centered adduct table
#'
#' @description
#' Center plate quant values based on median of quant measurements from BulkPool samples
#'
#' @param search_table
#'     Search table from HABC plate search
#'
#' @param samples_df
#'     HABC samples df as generated by \code{habc_samples_df()}
#'
#' @export
habc_bulkpool_centered_search_table <- function(search_table, samples_df) {
  search_table_ms1_median <- dplyr::left_join(search_table, samples_df, by = c("sample" = "sample_name")) %>%
    dplyr::filter(type == "BulkPool") %>%
    dplyr::group_by(compositionSummary, adductName, plate) %>%
    dplyr::mutate(ms1_intensity_median = median(ms1_intensity_is_nearest_scan_normalized, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      compositionSummary, adductName,
      ms1_intensity_median
    ) %>%
    unique()

  search_table_ms2_medians <- dplyr::left_join(search_table, samples_df, by = c("sample" = "sample_name")) %>%
    dplyr::filter(type == "BulkPool") %>%
    dplyr::group_by(compoundName, adductName, plate) %>%
    dplyr::mutate(
      diagnostic_fragment_sum_is_normalized =
        ifelse(diagnostic_fragment_sum_is_normalized == -1, NA, diagnostic_fragment_sum_is_normalized),
      diagnostic_fragment_sum_SAF_is_normalized =
        ifelse(diagnostic_fragment_sum_SAF_is_normalized == -1, NA, diagnostic_fragment_sum_SAF_is_normalized),
      acyl_fragment_sum_is_normalized =
        ifelse(acyl_fragment_sum_is_normalized == -1, NA, acyl_fragment_sum_is_normalized),
      acyl_fragment_sum_SAF_is_normalized =
        ifelse(acyl_fragment_sum_SAF_is_normalized == -1, NA, acyl_fragment_sum_SAF_is_normalized)
    ) %>%
    dplyr::mutate(ms2_diagnostic_sum_median = median(diagnostic_fragment_sum_is_normalized, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_diagnostic_sum_SAF_median = median(diagnostic_fragment_sum_SAF_is_normalized, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_acyl_sum_median = median(acyl_fragment_sum_is_normalized, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_acyl_sum_SAF_median = median(acyl_fragment_sum_SAF_is_normalized, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      compoundName, adductName,
      ms2_diagnostic_sum_median, ms2_acyl_sum_median, ms2_diagnostic_sum_SAF_median, ms2_acyl_sum_SAF_median
    ) %>%
    unique()

  search_table_centered <- search_table %>%
    dplyr::left_join(., search_table_ms1_median, by = c("compositionSummary", "adductName")) %>%
    dplyr::left_join(., search_table_ms2_medians, by = c("compoundName", "adductName")) %>%
    dplyr::mutate(
      diagnostic_fragment_sum_is_normalized =
        ifelse(diagnostic_fragment_sum_is_normalized == -1, NA, diagnostic_fragment_sum_is_normalized),
      diagnostic_fragment_sum_SAF_is_normalized =
        ifelse(diagnostic_fragment_sum_SAF_is_normalized == -1, NA, diagnostic_fragment_sum_SAF_is_normalized),
      acyl_fragment_sum_is_normalized =
        ifelse(acyl_fragment_sum_is_normalized == -1, NA, acyl_fragment_sum_is_normalized),
      acyl_fragment_sum_SAF_is_normalized =
        ifelse(acyl_fragment_sum_SAF_is_normalized == -1, NA, acyl_fragment_sum_SAF_is_normalized)
    ) %>%
    dplyr::mutate(
      ms1_intensity_is_nearest_scan_normalized = ms1_intensity_is_nearest_scan_normalized / ms1_intensity_median,
      diagnostic_fragment_sum_is_normalized = diagnostic_fragment_sum_is_normalized / ms2_diagnostic_sum_median,
      diagnostic_fragment_sum_SAF_is_normalized = diagnostic_fragment_sum_SAF_is_normalized / ms2_diagnostic_sum_SAF_median,
      acyl_fragment_sum_is_normalized = acyl_fragment_sum_is_normalized / ms2_acyl_sum_median,
      acyl_fragment_sum_SAF_is_normalized = acyl_fragment_sum_SAF_is_normalized / ms2_acyl_sum_SAF_median
    )
}

#' HABC IS search params
#'
#' @description
#' List of DIMS search params applied to IS searches
#'
#' @return
#' list of search params
#'
#' @export
habc_is_search_params <- function() {
  is_search_params <- list(
    "searchVersion" = "clamr HABC search",
    "consensusMinFractionMs2Scans" = 0.60, # 2/3 scans
    "ms2MinNumMatches" = 0,
    "ms2MinNumDiagnosticMatches" = 0,
    "spectralCompositionAlgorithm" = "AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION",
    "ms1PpmTolr" = 3,
    "ms2PpmTolr" = 10,
    "ms1ScanFilter" = "SIM",
    "consensusIntensityAgglomerationType" = "MEDIAN",
    "isReduceBySimpleParsimony" = TRUE,
    "consensusMs1PpmTolr" = 5,
    "ms1IsFindPrecursorIon" = TRUE,
    "ms1IsRequireMonoisotopic" = FALSE
  )
}

#' HABC MS3 search params
#'
#' @description
#' List of DIMS search params applied for MS3 searches (positive mode)
#'
#' @return
#' list of search params
#'
#' @export
habc_ms3_search_params <- function() {
  ms3_search_params <- list(
    "searchVersion" = "clamr HABC search",
    "ms1IsFindPrecursorIon" = FALSE,
    "ms1PpmTolr" = 5,
    "ms1ScanFilter" = "SIM",
    "ms3IsMs3Search" = TRUE,
    "ms3MinNumMatches" = 1,
    "ms3MinNumMs3MzMatches" = 1,
    "ms3AnalysisMs1PrecursorPpmTolr" = 20,
    "ms3PrecursorPpmTolr" = 20,
    "consensusIntensityAgglomerationType" = "MEDIAN",
    "ms3MinIntensity" = 100,
    "ms1IsRequireMonoisotopic" = FALSE,
    "isPreferSmallestScanMassWindow" = FALSE,
    "ms3MinFractionScans" = 0.5
  )
}


#' HABC positive mode search params
#'
#' @description
#' List of DIMS search params applied for negative mode MS1/MS2 searches
#'
#' @param version_num
#'     version number used for parameters
#'     By default, use most current version
#'
#' @return
#' list of search params
#'
#' @export
habc_pos_search_params <- function(version_num = 3) {
  ms2MinNumMatchesByLipidClassAndAdduct_pos <- tibble::tibble(
    "lipidClass" = c(
      "Alkyl_PC",
      "LPC",
      "SM",
      "Cholesterol-H2O",
      "CE",
      "DG"
    ),
    "adductName" = c(
      "[M+H]+",
      "[M+H]+",
      "[M+H]+",
      "*",
      "*",
      "*"
    ),
    "ms2MinNumMatches" = c(
      3,
      3,
      3,
      0,
      2,
      2
    )
  )

  if (version_num == 3) {
    ms2MinNumMatchesByLipidClassAndAdduct_pos$ms2MinNumMatches[1] <- 2 # Alkyl_PC [M+H]+
    ms2MinNumMatchesByLipidClassAndAdduct_pos$ms2MinNumMatches[2] <- 2 # LPC [M+H]+
    ms2MinNumMatchesByLipidClassAndAdduct_pos$ms2MinNumMatches[3] <- 1 # SM [M+H]+

    ms2MinNumMatchesByLipidClassAndAdduct_pos <- rbind(
      ms2MinNumMatchesByLipidClassAndAdduct_pos,
      tibble::tibble(
        lipidClass = c(
          "MG",
          "Carn",
          "LPE"
        ),
        adductName = c(
          "*",
          "*",
          "*"
        ),
        ms2MinNumMatches = c(
          0,
          0,
          2
        )
      )
    )
  }

  ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos <- tibble::tibble(
    "lipidClass" = c(
      "Alkyl_PC",
      "LPC",
      "SM",
      "Cholesterol-H2O",
      "CE",
      "DG"
    ),
    "adductName" = c(
      "[M+H]+",
      "[M+H]+",
      "[M+H]+",
      "*",
      "*",
      "*"
    ),
    "ms2MinNumDiagnosticMatches" = c(rep(0, 6))
  )

  if (version_num == 3) {
    ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos <- rbind(
      ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos,
      tibble::tibble(
        lipidClass = c(
          "MG",
          "Carn"
        ),
        adductName = c(
          "*",
          "*"
        ),
        ms2MinNumDiagnosticMatches = c(
          0,
          0
        )
      )
    )
  }

  pos_search_params <- list(
    "searchVersion" = "clamr HABC search",
    "consensusMinFractionMs2Scans" = 0.60, # 2/3 scans
    "ms2MinNumMatches" = 3,
    "ms2MinNumDiagnosticMatches" = 1,
    "ms2MinNumMatchesByLipidClassAndAdduct" = ms2MinNumMatchesByLipidClassAndAdduct_pos,
    "ms2MinNumDiagnosticMatchesByLipidClassAndAdduct" = ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos,
    "spectralCompositionAlgorithm" = "AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION",
    "ms1PpmTolr" = 3,
    "ms2PpmTolr" = 10,
    "ms1ScanFilter" = "SIM",
    "consensusIntensityAgglomerationType" = "MEDIAN",
    "isReduceBySimpleParsimony" = TRUE,
    "consensusMs1PpmTolr" = 5,
    "ms1IsFindPrecursorIon" = TRUE,
    "ms1IsRequireMonoisotopic" = FALSE
  )

  # Version 1 parameters up until this point
  if (version_num == 3) {
    # modifications for version 3
    pos_search_params["isReduceBySimpleParsimony"] <- FALSE
    pos_search_params["spectralCompositionAlgorithm"] <- "ALL_CANDIDATES"
  }

  pos_search_params
}

#' HABC whitelist search parameters
#'
#' @description
#' Initial search applied to bulk pool samples to generate a filtered library for HABC search.
#'
#' @param version_num
#'     version number used for parameters
#'     By default, use most current version
#'
#' @return
#' list of search params
#'
#' @export
habc_pos_whitelist_params <- function(version_num = 3) {
  whitelist_search_parameters <- habc_pos_search_params(version_num)

  ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos <- whitelist_search_parameters[["ms2MinNumDiagnosticMatchesByLipidClassAndAdduct"]]

  ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos_cleaned <- ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos %>%
    dplyr::filter(lipidClass != "DG" & lipidClass != "CE")

  ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos <- rbind(
    ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos_cleaned,
    tibble::tibble(
      lipidClass = c(
        "SM",
        "DG", "DG", "DG", "DG",
        "CE", "CE", "CE", "CE",
        "Ceramide",
        "HexCer",
        "LPE"
      ),
      adductName = c(
        "[M+Na]+",
        "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+",
        "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+",
        "[M+H]+",
        "[M+H]+",
        "[M+H]+"
      ),
      ms2MinNumDiagnosticMatches = c(
        2,
        0, 0, 0, 1,
        0, 0, 0, 1,
        2,
        1,
        1
      )
    )
  )

  whitelist_search_parameters[["ms2MinNumDiagnosticMatchesByLipidClassAndAdduct"]] <- ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_pos

  ms2sn1MinNumMatchesByLipidClassAndAdduct_pos <- tibble::tibble(
    lipidClass = c(
      "DG",
      "Ceramide",
      "HexCer",
      "LPE",
      "LPC",
      "Alkyl_PC"
    ),
    adductName = c(
      "[M+NH4]+",
      "[M+H]+",
      "[M+H]+",
      "[M+H]+",
      "[M+H]+",
      "[M+H]+"
    ),
    ms2sn1MinNumMatches = c(
      1,
      2,
      2,
      1,
      1,
      1
    )
  )

  ms2sn2MinNumMatchesByLipidClassAndAdduct_pos <- tibble::tibble(
    lipidClass = c(
      "SM",
      "DG",
      "Ceramide",
      "HexCer",
      "Alkyl_PC"
    ),
    adductName = c(
      "[M+Na]+",
      "[M+NH4]+",
      "[M+H]+",
      "[M+H]+",
      "[M+H]+"
    ),
    ms2sn2MinNumMatches = c(
      1,
      1,
      1,
      1,
      1
    )
  )

  ms2IsRequirePrecursorMatchByLipidClassAndAdduct_pos <- tibble::tibble(
    lipidClass = c(
      "Carn",
      "LPC",
      "Alkyl_PC"
    ),
    adductName = c(
      "[M+H]+",
      "[M+H]+",
      "[M+H]+"
    ),
    ms2IsRequirePrecursorMatch = c(
      TRUE,
      TRUE,
      TRUE
    )
  )

  whitelist_search_parameters <- c(
    whitelist_search_parameters,
    list(
      "ms2sn1MinNumMatchesByLipidClassAndAdduct" = ms2sn1MinNumMatchesByLipidClassAndAdduct_pos,
      "ms2sn2MinNumMatchesByLipidClassAndAdduct" = ms2sn2MinNumMatchesByLipidClassAndAdduct_pos,
      "ms2IsRequirePrecursorMatchByLipidClassAndAdduct" = ms2IsRequirePrecursorMatchByLipidClassAndAdduct_pos
    )
  )

  whitelist_search_parameters
}

#' HABC negative mode search params
#'
#' @description
#' List of DIMS search params applied for negative mode MS1/MS2 searches
#'
#' @param version_num
#'     version number used for parameters
#'     By default, use most current version
#'
#' @return
#' list of search params
#'
#' @export
habc_neg_search_params <- function(version_num = 2) {
  ms2MinNumMatchesByLipidClassAndAdduct_neg <- tibble::tibble(
    "lipidClass" = c("LPI", "PE", "PG", "FA", "LPC", "LPE", "LPG", "SM"),
    "adductName" = c(rep("*", 8)),
    "ms2MinNumMatches" = as.integer(c(3, 3, 3, 0, 2, 2, 2, 2))
  )

  ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_neg <- tibble::tibble(
    "lipidClass" = c("LPI", "PE", "PG", "FA", "LPC", "LPE", "LPG", "SM"),
    "adductName" = c(rep("*", 8)),
    "ms2MinNumDiagnosticMatches" = as.integer(c(rep(0, 8)))
  )

  neg_search_params <- list(
    "searchVersion" = "clamr HABC search",
    "consensusMinFractionMs2Scans" = 0.60, # 2/3 scans
    "ms2MinNumMatches" = 3,
    "ms2MinNumDiagnosticMatches" = 1,
    "ms2MinNumMatchesByLipidClassAndAdduct" = ms2MinNumMatchesByLipidClassAndAdduct_neg,
    "ms2MinNumDiagnosticMatchesByLipidClassAndAdduct" = ms2MinNumDiagnosticMatchesByLipidClassAndAdduct_neg,
    "spectralCompositionAlgorithm" = "AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION",
    "ms1PpmTolr" = 3,
    "ms2PpmTolr" = 10,
    "ms1ScanFilter" = "SIM",
    "consensusIntensityAgglomerationType" = "MEDIAN",
    "isReduceBySimpleParsimony" = TRUE,
    "consensusMs1PpmTolr" = 5,
    "ms1IsFindPrecursorIon" = TRUE,
    "ms1IsRequireMonoisotopic" = FALSE
  )

  # Version 1 parameters up until this point
  if (version_num == 2) {
    # modifications for version 2
    neg_search_params["isReduceBySimpleParsimony"] <- FALSE
    neg_search_params["spectralCompositionAlgorithm"] <- "ALL_CANDIDATES"
  }

  neg_search_params
}

#' HABC whitelist search parameters
#'
#' @description
#' Initial search applied to bulk pool samples to generate a filtered library for HABC search.
#'
#' @param version_num
#'     version number used for parameters
#'     By default, use most current version
#'
#' @return
#' list of search params
#'
#' @export
habc_neg_whitelist_params <- function(version_num = 2) {
  whitelist_search_parameters <- habc_neg_search_params(2)

  ms2sn1MinNumMatchesByLipidClassAndAdduct_neg <- tibble::tibble(
    "lipidClass" = c("Ceramide"),
    "adductName" = c("*"),
    "ms2sn1MinNumMatches" = c(2)
  )

  ms2sn2MinNumMatchesByLipidClassAndAdduct_neg <- tibble::tibble(
    "lipidClass" = c("Ceramide"),
    "adductName" = c("*"),
    "ms2sn2MinNumMatches" = c(2)
  )

  whitelist_search_parameters <- c(
    whitelist_search_parameters,
    list(
      "ms2sn1MinNumMatches" = 1,
      "ms2sn2MinNumMatches" = 1,
      "ms2sn1MinNumMatchesByLipidClassAndAdduct" = ms2sn1MinNumMatchesByLipidClassAndAdduct_neg,
      "ms2sn2MinNumMatchesByLipidClassAndAdduct" = ms2sn2MinNumMatchesByLipidClassAndAdduct_neg
    )
  )

  whitelist_search_parameters
}

#' HABC Exclude Bad Bulkpool samples
#'
#' @description
#' Return a filtered list of good bulkpool samples from an IS search results table
#'
#' @return
#' vector of bulkpool samples that are good
#'
#' @export
habc_good_bulkpools <- function(
  bulkpool_is_adduct_table,
  min_frac_intensity = 0.5,
  min_frac_detected = 0.5
) {
  # Determine intensity sum
  bulkpool_intensity_sum <- bulkpool_is_adduct_table %>%
    dplyr::filter(is_identified) %>%
    dplyr::select(sample, compoundName, adductName, ms1_scan_intensity) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(intensity_sum = sum(ms1_scan_intensity)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, intensity_sum) %>%
    unique()

  max_intensity <- max(bulkpool_intensity_sum$intensity_sum)

  bulkpool_TIC_stats <- bulkpool_intensity_sum %>%
    dplyr::mutate(frac_intensity = intensity_sum / max_intensity) %>%
    dplyr::mutate(frac_detected = -1)

  # determine presence/absence of each IS ion (ion must be found in at least one sample to be on the list)
  all_bulkpool_samples <- bulkpool_is_adduct_table$sample %>% unique()

  bulkpool_presence_absence <- bulkpool_is_adduct_table %>%
    dplyr::filter(is_identified) %>%
    dplyr::select(sample, compoundName, adductName, ms1_scan_intensity) %>%
    tidyr::pivot_wider(names_from = sample, values_from = ms1_scan_intensity)

  num_ions <- bulkpool_presence_absence %>% nrow()

  for (i in 1:length(all_bulkpool_samples)) {
    bulkpool_sample <- all_bulkpool_samples[i]

    num_detected <- length(which(!is.na(bulkpool_presence_absence[[bulkpool_sample]])))

    bulkpool_TIC_stats$frac_detected[i] <- num_detected / num_ions
  }

  # filter based on quality criteria
  bulkpool_filtered_TIC_stats <- bulkpool_TIC_stats %>%
    dplyr::filter(frac_intensity >= min_frac_intensity & frac_detected >= min_frac_detected)

  # return list of good samples
  bulkpool_filtered_TIC_stats
}


#' HABC CV comparison
#'
#' @description
#' Compare CV of BulkPool vs barcoded samples
#'
#' @param habc_quant_table
#'    Formatted quant output table from \code{habc_process_plate_neg_v2()}
#'
#' @return
#' Table with comparison of CVs for MS1 and MS2 quant types (if computable)
#'
#' @export
habc_cv_comparison <- function(habc_quant_table) {
  ms1_cv_comparison <- habc_quant_table %>%
    dplyr::filter(type %in% c("bc", "BulkPool")) %>%
    dplyr::filter(!is.na(ms1_intensity)) %>%
    dplyr::group_by(compoundName, adductName, type) %>%
    dplyr::mutate(cv = sd(ms1_intensity, na.rm = TRUE) / mean(ms1_intensity, na.rm = TRUE)) %>%
    dplyr::mutate(ms1_cv_quant_type = case_when(
      any(ms1_quant_type == "cross_class_acyl_SAF_partitioned_intensity") ~ "cross_class_acyl_SAF_partitioned_intensity",
      any(ms1_quant_type == "cross_class_acyl_partitioned_intensity") ~ "cross_class_acyl_partitioned_intensity",
      any(ms1_quant_type == "acyl_SAF_partitioned_intensity") ~ "acyl_SAF_partitioned_intensity",
      any(ms1_quant_type == "acyl_partitioned_intensity") ~ "acyl_partitioned_intensity",
      any(ms1_quant_type == "ms1_intensity") ~ "ms1_intensity",
      TRUE == TRUE ~ ""
    )) %>%
    dplyr::ungroup() %>%
    dplyr::select(compoundName, adductName, type, cv, ms1_cv_quant_type) %>%
    unique() %>%
    tidyr::pivot_wider(., names_from = type, values_from = cv) %>%
    dplyr::rename(bc_ms1 = bc, BulkPool_ms1 = BulkPool) %>%
    dplyr::mutate(bc_ms1 = ifelse(is.na(bc_ms1), -1, bc_ms1)) %>%
    dplyr::mutate(BulkPool_ms1 = ifelse(is.na(BulkPool_ms1), -1, BulkPool_ms1)) %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(
      bc_ms1 = max(bc_ms1),
      BulkPool_ms1 = max(BulkPool_ms1),
      ms1_cv_quant_type = case_when(
        any(ms1_cv_quant_type == "cross_class_acyl_SAF_partitioned_intensity") ~ "cross_class_acyl_SAF_partitioned_intensity",
        any(ms1_cv_quant_type == "cross_class_acyl_partitioned_intensity") ~ "cross_class_acyl_partitioned_intensity",
        any(ms1_cv_quant_type == "acyl_SAF_partitioned_intensity") ~ "acyl_SAF_partitioned_intensity",
        any(ms1_cv_quant_type == "acyl_partitioned_intensity") ~ "acyl_partitioned_intensity",
        any(ms1_cv_quant_type == "ms1_intensity") ~ "ms1_intensity",
        TRUE == TRUE ~ ""
      )
    ) %>%
    dplyr::ungroup() %>%
    unique()

  ms2_cv_comparison <- habc_quant_table %>%
    dplyr::filter(type %in% c("bc", "BulkPool")) %>%
    dplyr::group_by(compoundName, adductName, type) %>%
    dplyr::mutate(cv = sd(ms2_intensity, na.rm = TRUE) / mean(ms2_intensity, na.rm = TRUE)) %>%
    dplyr::mutate(ms2_cv_quant_type = case_when(
      any(ms2_quant_type == "acyl_fragment_sum_SAF_is_normalized") ~ "acyl_fragment_sum_SAF_is_normalized",
      any(ms2_quant_type == "acyl_fragment_sum_is_normalized") ~ "acyl_fragment_sum_is_normalized",
      any(ms2_quant_type == "diagnostic_ms2_intensity") ~ "diagnostic_ms2_intensity",
      TRUE == TRUE ~ ""
    )) %>%
    dplyr::ungroup() %>%
    dplyr::select(compoundName, adductName, type, cv, ms2_cv_quant_type) %>%
    unique() %>%
    tidyr::pivot_wider(., names_from = type, values_from = cv) %>%
    dplyr::rename(bc_ms2 = bc, BulkPool_ms2 = BulkPool) %>%
    dplyr::mutate(bc_ms2 = ifelse(is.na(bc_ms2), -1, bc_ms2)) %>%
    dplyr::mutate(BulkPool_ms2 = ifelse(is.na(BulkPool_ms2), -1, BulkPool_ms2)) %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(
      bc_ms2 = max(bc_ms2),
      BulkPool_ms2 = max(BulkPool_ms2),
      ms2_cv_quant_type = case_when(
        any(ms2_cv_quant_type == "acyl_fragment_sum_SAF_is_normalized") ~ "acyl_fragment_sum_SAF_is_normalized",
        any(ms2_cv_quant_type == "acyl_fragment_sum_is_normalized") ~ "acyl_fragment_sum_is_normalized",
        any(ms2_cv_quant_type == "diagnostic_ms2_intensity") ~ "diagnostic_ms2_intensity",
        TRUE == TRUE ~ ""
      )
    ) %>%
    dplyr::ungroup() %>%
    unique()

  cv_comparison <- dplyr::inner_join(ms1_cv_comparison, ms2_cv_comparison, by = c("compoundName", "adductName")) %>%
    # NA at this point means that there was only one barcoded sample, so no CV could be computed
    dplyr::mutate(
      bc_ms1 = ifelse(bc_ms1 == -1, NA, bc_ms1),
      bc_ms2 = ifelse(bc_ms2 == -1, NA, bc_ms2)
    )
}

#' HABC CV comparison
#'
#' @description
#' Compare CV of BulkPool vs barcoded samples
#' Updated version (separating on 'quant_class' as ms1, ms2, or ms3)
#'
#' @param habc_quant_table
#'    Formatted quant output table from \code{habc_process_plate_neg_v2()}
#' @param biological_type
#'    String code for biological sample. by default (for HABC), use \code{"bc"}
#'
#' @return
#' table with columns compoundName, adductName, quant_class, median_BulkPool_CV, median_bc_CV
#'
#' @export
habc_cv_comparison_v3 <- function(habc_quant_table_v3, biological_type = "bc") {
  ms_levels <- c("ms1", "ms2", "ms3")

  CV_results <- tibble::tibble(
    compoundName = character(0),
    adductName = character(0),
    quant_class = character(0),
    BulkPool_CV = numeric(0),
    bc_CV = numeric(0)
  )

  for (i in 1:length(ms_levels)) {
    ms_level <- ms_levels[i]

    quant_BulkPool <- habc_quant_table_v3 %>%
      dplyr::filter(quant_class == ms_level & type == "BulkPool") %>%
      dplyr::group_by(compoundName, adductName) %>%
      dplyr::mutate(BulkPool_CV = sd(intensity) / mean(intensity)) %>%
      dplyr::ungroup() %>%
      dplyr::select(compoundName, adductName, quant_class, BulkPool_CV) %>%
      unique()

    quant_bc <- habc_quant_table_v3 %>%
      dplyr::filter(quant_class == ms_level & type == biological_type) %>%
      dplyr::group_by(compoundName, adductName) %>%
      dplyr::mutate(bc_CV = sd(intensity) / mean(intensity)) %>%
      dplyr::ungroup() %>%
      dplyr::select(compoundName, adductName, quant_class, bc_CV) %>%
      unique()

    quant_combined <- quant_BulkPool %>%
      dplyr::full_join(quant_bc, by = c("compoundName", "adductName", "quant_class"))

    CV_results <- rbind(CV_results, quant_combined)
  }

  return(CV_results)
}

#' HABC process plate (positive mode data)
#'
#' @description
#' Process a HABC plate collected in positive mode and save results on file system
#'
#' @param plate_name
#'        name of folder containing mzML files associated with a single HABC plate
#'
#' @param samples_file_path
#'        directory of HABC mzML files
#'
#' @param rds_dir
#'        parent directory for saving .RDS files (rds_dir/plate_name/<file>)
#'
#' @param mzrolldb_dir
#'        parent directory for saving mzrollDB files (mzrolldb_dir/plate_name/<file>)
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing only post-lle IS compounds
#'        located at (lib_dir/is_lib_name)
#'
#' @param pre_lle_is_lib_name
#'        IS library containing only pre-lle IS compounds
#'        located at (lib_dir/pre_lle_is_lib_name)
#'
#' @param pos_lib_name
#'        HABC positive library name (lib_dir/pos_lib_name)
#'
#' @param tg_lib_name
#'        HABC TG library name (lib_dir/tg_lib_name)
#'
#' @param pos_search_params
#'        list of search params. Default behavior is habc::habc_pos_search_params()
#'
#' @export
habc_process_plate_pos_v1 <- function(
  plate_name,
  samples_file_path,
  rds_dir,
  mzrolldb_dir,
  lib_dir,
  is_lib_name,
  pre_lle_is_lib_name,
  pos_lib_name,
  tg_lib_name,
  pos_search_params = list()
) {
  # files ######################################
  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  is_pre_lle_lib_file <- paste(lib_dir, pre_lle_is_lib_name, sep = "/")
  habc_pos_lib_file <- paste(lib_dir, pos_lib_name, sep = "/")
  habc_tg_lib_file <- paste(lib_dir, tg_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")

  # samples ######################################
  samples_info <- habc_samples_df(samples_file_path, TRUE)
  samples_df_pos <- samples_info$samples_df
  if (nrow(samples_df_pos) == 0) {
    return(invisible(0))
  }
  ms2_ranges_pos <- samples_info$ms2_ranges
  ms3_targets <- samples_info$ms3_targets

  # libraries ####################################
  is_lib_pos <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl("\\+$", adductName))
  is_lib_pos_sliced <- mzkitcpp::DI_slice_library(ms2_ranges_pos, is_lib_pos)
  is_lib_pre_lle_pos <- mzkitcpp::import_msp_lipids_library(is_pre_lle_lib_file) %>% dplyr::filter(grepl("\\+$", adductName))
  is_lib_pre_lle_pos_sliced <- mzkitcpp::DI_slice_library(ms2_ranges_pos, is_lib_pre_lle_pos)
  habc_lib_pos <- mzkitcpp::import_msp_lipids_library(habc_pos_lib_file)
  habc_tg_lib <- mzkitcpp::import_msp_lipids_library(habc_tg_lib_file)

  # search parameters ############################
  is_search_params <- habc_is_search_params()

  if (length(pos_search_params) == 0) {
    pos_search_params <- habc_pos_search_params()
  }

  ms3_search_params <- habc_ms3_search_params()

  # pos IS Search ##################
  cat("Starting pos IS Search... \n")

  pos_raw_results_is <- mzkitcpp::DI_pipeline(
    samples = samples_df_pos$file,
    ms2_ranges = ms2_ranges_pos,
    is_sliced_lib = is_lib_pos_sliced,
    is_search_params = is_search_params,
    sliced_lib = is_lib_pos_sliced,
    search_params = is_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  pos_raw_results_is_w_metadata <- pos_raw_results_is$search %>%
    dplyr::inner_join(samples_df_pos, by = c("sample" = "sample_name"))

  cat("Finished pos IS Search.\n")

  # pos pre-lle Search ##################
  cat("Starting pos pre-lle IS Search... \n")

  pos_pre_lle_is_results_is <- mzkitcpp::DI_pipeline(
    samples = samples_df_pos$file,
    ms2_ranges = ms2_ranges_pos,
    is_sliced_lib = is_lib_pos_sliced,
    is_search_params = is_search_params,
    sliced_lib = is_lib_pre_lle_pos_sliced,
    search_params = is_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  pos_pre_lle_is_results_is_w_metadata <- pos_pre_lle_is_results_is$search %>%
    dplyr::inner_join(samples_df_pos, by = c("sample" = "sample_name"))

  cat("Finished pos pre-lle IS Search.\n")

  # pos Search #####################
  cat("Starting pos Search... \n")

  pos_is_lipid_class_adduct <- pos_raw_results_is$search %>%
    dplyr::select(lipidClass, adductName) %>%
    unique()

  habc_pos_lib_filtered <- habc_lib_pos %>%
    dplyr::inner_join(pos_is_lipid_class_adduct, by = c("lipidClass", "adductName"))

  habc_pos_lib_sliced <- mzkitcpp::DI_slice_library(ms2_ranges_pos, habc_pos_lib_filtered)

  pos_raw_results <- mzkitcpp::DI_pipeline(
    samples = samples_df_pos$file,
    ms2_ranges = ms2_ranges_pos,
    is_sliced_lib = is_lib_pos_sliced,
    is_search_params = is_search_params,
    sliced_lib = habc_pos_lib_sliced,
    search_params = pos_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  pos_raw_results_w_metadata <- pos_raw_results$search %>%
    dplyr::inner_join(samples_df_pos, by = c("sample" = "sample_name"))

  cat("Finished pos Search.\n")

  # TG Search ######################
  cat("Starting TG Search... \n")

  habc_tg_lib_all_targets <- habc_tg_lib %>%
    to_ms3_lib() %>%
    dplyr::select(prec_mzs) %>%
    unique()

  targets_both <- dplyr::inner_join(habc_tg_lib_all_targets, ms3_targets, by = c("prec_mzs"))

  habc_tg_lib_filtered <- habc_tg_lib %>%
    to_ms3_lib() %>%
    dplyr::filter(prec_mzs %in% targets_both$prec_mzs) %>%
    to_ms2_lib()

  tg_raw_results <- mzkitcpp::DI_pipeline_ms3_search(
    samples = samples_df_pos$file,
    is_lib = default_tg_is_ms3,
    is_search_params = ms3_search_params,
    search_lib = habc_tg_lib_filtered,
    search_params = ms3_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  tg_raw_results_w_metadata <- tg_raw_results %>%
    dplyr::inner_join(samples_df_pos, by = c("sample" = "sample_name"))

  cat("Finished TG Search.\n")

  # save rds results ###############

  results_summary <- list(
    "samples_df_pos" = samples_df_pos,
    "ms2_ranges_pos" = ms2_ranges_pos,
    "pos_is_results" = pos_raw_results_is_w_metadata,
    "pos_pre_lle_is_results" = pos_pre_lle_is_results_is_w_metadata,
    "pos_search_results" = pos_raw_results_w_metadata,
    "pos_adduct_table" = pos_raw_results$adduct_table,
    "tg_results" = tg_raw_results_w_metadata
  )

  rds_results_file <- paste0(rds_dir, "/X0158_M014A_", plate_name, ".rds")

  saveRDS(results_summary, file = rds_results_file)

  # save mzrollDB results ##########

  mzrolldb_results_file <- paste0(mzrolldb_dir, "/X0158_M014A_", plate_name, ".mzrollDB")

  system(glue::glue("rm {old_mzroll_db_file} 2>&1", old_mzroll_db_file = mzrolldb_results_file))

  add_direct_infusion_search_results(
    mzrolldb_results_file,
    samples_df_pos$file,
    ms2_ranges_pos,
    is_lib_file, # library_name
    is_lib_pos_sliced,
    is_search_params,
    adducts_file,
    pos_raw_results_is$search
  )

  add_direct_infusion_search_results(
    mzrolldb_results_file,
    samples_df_pos$file,
    ms2_ranges_pos,
    is_pre_lle_lib_file, # library_name
    is_lib_pre_lle_pos_sliced,
    is_search_params,
    adducts_file,
    pos_pre_lle_is_results_is$search
  )

  add_direct_infusion_search_results(
    mzrolldb_results_file,
    samples_df_pos$file,
    ms2_ranges_pos,
    habc_pos_lib_file, # library_name
    habc_pos_lib_sliced,
    pos_search_params,
    adducts_file,
    pos_raw_results$search
  )

  add_targeted_ms3_search_results(
    mzrolldb_results_file,
    samples_df_pos$file,
    tg_raw_results,
    habc_tg_lib_file,
    ms3_search_params
  )

  # end ############################
  invisible(0)
}

#' HABC process plate version 1 (negative mode data)
#'      Version 1 used late 2020-early 2021
#'
#' @description
#' Process a HABC plate collected in negative mode and save results on file system
#'
#' @param plate_name
#'        name of folder containing mzML files associated with a single HABC plate
#'
#' @param samples_file_path
#'        directory of HABC mzML files
#'
#' @param rds_dir
#'        parent directory for saving .RDS files (rds_dir/plate_name/<file>)
#'
#' @param mzrolldb_dir
#'        parent directory for saving mzrollDB files (mzrolldb_dir/plate_name/<file>)
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing only post-lle IS compounds
#'        located at (lib_dir/is_lib_name)
#'
#' @param pre_lle_is_lib_name
#'        IS library containing only pre-lle IS compounds
#'        located at (lib_dir/pre_lle_is_lib_name)
#'
#' @param neg_lib_name
#'        HABC positive library name (lib_dir/pos_lib_name)
#'
#' @param neg_search_params
#'        list of search params. default behavior is habc::habc_neg_search_params()
#'
#' @export
habc_process_plate_neg_v1 <- function(
  plate_name,
  samples_file_path,
  rds_dir,
  mzrolldb_dir,
  lib_dir,
  is_lib_name,
  pre_lle_is_lib_name,
  neg_lib_name,
  neg_search_params = list()
) {
  # files ######################################
  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  is_pre_lle_lib_file <- paste(lib_dir, pre_lle_is_lib_name, sep = "/")
  habc_neg_lib_file <- paste(lib_dir, neg_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")

  # samples ######################################
  samples_info <- habc_samples_df(samples_file_path, FALSE)
  samples_df_neg <- samples_info$samples_df
  if (nrow(samples_df_neg) == 0) {
    return(invisible(0))
  }
  ms2_ranges_neg <- samples_info$ms2_ranges

  # libraries ####################################
  is_lib_neg <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl("\\-$", adductName))
  is_lib_neg_sliced <- mzkitcpp::DI_slice_library(ms2_ranges_neg, is_lib_neg)
  is_lib_pre_lle_neg <- mzkitcpp::import_msp_lipids_library(is_pre_lle_lib_file) %>% dplyr::filter(grepl("\\-$", adductName))
  is_lib_pre_lle_neg_sliced <- mzkitcpp::DI_slice_library(ms2_ranges_neg, is_lib_pre_lle_neg)
  habc_lib_neg <- mzkitcpp::import_msp_lipids_library(habc_neg_lib_file)

  # search parameters ############################
  is_search_params <- habc_is_search_params()

  if (length(neg_search_params) == 0) {
    neg_search_params <- habc_neg_search_params()
  }

  # neg IS Search ##################
  cat("Starting neg IS Search... \n")

  neg_raw_results_is <- mzkitcpp::DI_pipeline(
    samples = samples_df_neg$file,
    ms2_ranges = ms2_ranges_neg,
    is_sliced_lib = is_lib_neg_sliced,
    is_search_params = is_search_params,
    sliced_lib = is_lib_neg_sliced,
    search_params = is_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  neg_raw_results_is_w_metadata <- neg_raw_results_is$search %>%
    dplyr::inner_join(samples_df_neg, by = c("sample" = "sample_name"))

  cat("Finished neg IS Search.\n")

  # neg pre-lle Search ##################
  cat("Starting neg pre-lle IS Search... \n")

  neg_pre_lle_is_results_is <- mzkitcpp::DI_pipeline(
    samples = samples_df_neg$file,
    ms2_ranges = ms2_ranges_neg,
    is_sliced_lib = is_lib_neg_sliced,
    is_search_params = is_search_params,
    sliced_lib = is_lib_pre_lle_neg_sliced,
    search_params = is_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  neg_pre_lle_is_results_is_w_metadata <- neg_pre_lle_is_results_is$search %>%
    dplyr::inner_join(samples_df_neg, by = c("sample" = "sample_name"))

  cat("Finished neg pre-lle IS Search.\n")

  # neg Search #####################
  cat("Starting neg Search... \n")

  neg_is_lipid_class_adduct <- neg_raw_results_is$search %>%
    dplyr::select(lipidClass, adductName) %>%
    unique()

  habc_neg_lib_filtered <- habc_lib_neg %>%
    dplyr::inner_join(neg_is_lipid_class_adduct, by = c("lipidClass", "adductName"))

  habc_neg_lib_sliced <- mzkitcpp::DI_slice_library(ms2_ranges_neg, habc_neg_lib_filtered)

  neg_raw_results <- mzkitcpp::DI_pipeline(
    samples = samples_df_neg$file,
    ms2_ranges = ms2_ranges_neg,
    is_sliced_lib = is_lib_neg_sliced,
    is_search_params = is_search_params,
    sliced_lib = habc_neg_lib_sliced,
    search_params = neg_search_params,
    adducts_file = adducts_file,
    debug = F
  )

  neg_raw_results_w_metadata <- neg_raw_results$search %>%
    dplyr::inner_join(samples_df_neg, by = c("sample" = "sample_name"))

  cat("Finished neg Search.\n")
  # save rds results ###############

  results_summary <- list(
    "samples_df_neg" = samples_df_neg,
    "ms2_ranges_neg" = ms2_ranges_neg,
    "neg_is_results" = neg_raw_results_is_w_metadata,
    "neg_pre_lle_is_results" = neg_pre_lle_is_results_is_w_metadata,
    "neg_search_results" = neg_raw_results_w_metadata,
    "neg_adduct_table" = neg_raw_results$adduct_table
  )

  rds_results_file <- paste0(rds_dir, "/X0158_M015A_", plate_name, ".rds")

  saveRDS(results_summary, file = rds_results_file)

  # save mzrollDB results ##########

  mzrolldb_results_file <- paste0(mzrolldb_dir, "/X0158_M015A_", plate_name, ".mzrollDB")

  system(glue::glue("rm {old_mzroll_db_file} 2>&1", old_mzroll_db_file = mzrolldb_results_file))

  add_direct_infusion_search_results(
    mzrolldb_results_file,
    samples_df_neg$file,
    ms2_ranges_neg,
    is_lib_file, # library_name
    is_lib_neg_sliced,
    is_search_params,
    adducts_file,
    neg_raw_results_is$search
  )

  add_direct_infusion_search_results(
    mzrolldb_results_file,
    samples_df_neg$file,
    ms2_ranges_neg,
    is_pre_lle_lib_file, # library_name
    is_lib_pre_lle_neg_sliced,
    is_search_params,
    adducts_file,
    neg_pre_lle_is_results_is$search
  )

  add_direct_infusion_search_results(
    mzrolldb_results_file,
    samples_df_neg$file,
    ms2_ranges_neg,
    habc_neg_lib_file, # library_name
    habc_neg_lib_sliced,
    neg_search_params,
    adducts_file,
    neg_raw_results$search
  )

  # end ############################
  invisible(0)
}

#' Retrieve Specific HABC RDS files
#'
#' @description
#' Retrieve specific RDS plates, and output in desired format.
#'
#' @param rds_dir
#'        Directory for saving .RDS files (rds_dir/<file>)
#'
#' @param rds_file_pattern
#'        Filter RDS files based on rds_file_pattern.
#'        Empty string will return all matches.
#'
#' @return
#' \code{
#' list("pos_is_results"=pos_is_results,
#' "pos_pre_lle_is_results"=pos_pre_lle_is_results,
#' "pos_search_results"=pos_search_results,
#' "pos_adduct_table"=pos_adduct_table,
#' "tg_results"=tg_results,
#' "neg_is_results"=neg_is_results,
#' "neg_pre_lle_is_results"=neg_pre_lle_is_results,
#' "neg_search_results"=neg_search_results,
#' "neg_adduct_table"=neg_adduct_table)
#' }
#'
#' \itemize{
#'    \item{\code{pos_is_results}}{: positive mode search results of post-LLE IS.}
#'    \item{\code{pos_pre_lle_is_results}}{: positive mode search results of pre-LLE IS.}
#'    \item{\code{pos_search_results}}{: positive mode search results of HABC compounds.}
#'    \item{\code{pos_adduct_table}}{: positive mode table of compound adduct MS1 scan intensities.}
#'    \item{\code{tg_results}}{: Positive mode MS3 search of TGs.}
#'    \item{\code{neg_is_results}}{: negativee mode search results of post-LLE IS.}
#'    \item{\code{neg_pre_lle_is_results}}{: negative mode search results of pre-LLE IS.}
#'    \item{\code{neg_search_results}}{: negative mode search results of HABC compounds.}
#'    \item{\code{neg_adduct_table}}{: negative mode table of compound adduct MS1 scan intensities.}
#' }
#'
#' @export
habc_get_rds_v1 <- function(rds_dir, rds_file_pattern = "") {
  # Read values from .rds, and concatenate into tables

  # pos
  pos_pattern <- paste0("^X0158_M014A.*", rds_file_pattern)
  pos_rds <- list.files(rds_dir, pattern = pos_pattern, full.names = TRUE)

  cat(paste0("Discovered ", length(pos_rds), " positive mode plates.\n"))

  if (length(pos_rds) > 0) {
    first_result <- readRDS(pos_rds[1])

    pos_search_results <- first_result$pos_search_results
    pos_adduct_table <- first_result$pos_adduct_table
    pos_is_results <- first_result$pos_is_results
    pos_pre_lle_is_results <- first_result$pos_pre_lle_is_results
    tg_results <- first_result$tg_results

    if (length(pos_rds) > 1) {
      for (x in 2:length(pos_rds)) {
        results <- readRDS(pos_rds[x])

        pos_search_results <- rbind(pos_search_results, results$pos_search_results)
        pos_adduct_table <- rbind(pos_adduct_table, results$pos_adduct_table)
        pos_is_results <- rbind(pos_is_results, results$pos_is_results)
        pos_pre_lle_is_results <- rbind(pos_pre_lle_is_results, results$pos_pre_lle_is_results)
        tg_results <- rbind(tg_results, results$tg_results)
      }
    }
  } else {
    pos_search_results <- tibble::tibble()
    pos_adduct_table <- tibble::tibble()
    pos_is_results <- tibble::tibble()
    pos_pre_lle_is_results <- tibble::tibble()
    tg_results <- tibble::tibble()
  }

  # neg
  neg_pattern <- paste0("^X0158_M015A.*", rds_file_pattern)
  neg_rds <- list.files(rds_dir, pattern = neg_pattern, full.names = TRUE)

  cat(paste0("Discovered ", length(neg_rds), " negative mode plates.\n"))

  if (length(neg_rds) > 0) {
    first_result <- readRDS(neg_rds[1])

    neg_is_results <- first_result$neg_is_results
    neg_pre_lle_is_results <- first_result$neg_pre_lle_is_results
    neg_search_results <- first_result$neg_search_results
    neg_adduct_table <- first_result$neg_adduct_table

    if (length(neg_rds) > 1) {
      for (x in 2:length(neg_rds)) {
        results <- readRDS(neg_rds[x])

        neg_is_results <- rbind(neg_is_results, results$neg_is_results)
        neg_pre_lle_is_results <- rbind(neg_pre_lle_is_results, results$neg_pre_lle_is_results)
        neg_search_results <- rbind(neg_search_results, results$neg_search_results)
        neg_adduct_table <- rbind(neg_adduct_table, results$neg_adduct_table)
      }
    }
  } else {
    neg_is_results <- tibble::tibble()
    neg_pre_lle_is_results <- tibble::tibble()
    neg_search_results <- tibble::tibble()
    neg_adduct_table <- tibble::tibble()
  }

  cat(paste0(
    "Returning large list(\n",
    "\tpos_is_results,\n",
    "\tpos_pre_lle_is_results,\n",
    "\tpos_search_results,\n",
    "\tpos_adduct_table,\n",
    "\ttg_results,\n",
    "\tneg_is_results,\n",
    "\tneg_pre_lle_is_results,\n",
    "\tneg_search_results,\n",
    "\tneg_adduct_table\n",
    ")\n"
  ))

  list(
    "pos_is_results" = pos_is_results,
    "pos_pre_lle_is_results" = pos_pre_lle_is_results,
    "pos_search_results" = pos_search_results,
    "pos_adduct_table" = pos_adduct_table,
    "tg_results" = tg_results,
    "neg_is_results" = neg_is_results,
    "neg_pre_lle_is_results" = neg_pre_lle_is_results,
    "neg_search_results" = neg_search_results,
    "neg_adduct_table" = neg_adduct_table
  )
}

#' HABC process all plates (positive and negative mode data)
#'
#' @description
#' Process all HABC plate data found on file system
#'
#' @param top_level_samples_dir
#'        HABC samples, ignorized into sub-folders by plate name.
#'
#'        Require this arrangement on file system:
#'
#'        top_level_samples_dir/
#'                     X0158_M014A/       (positive)
#'                              plate1/
#'                              plate2/
#'                              ...
#'
#'                      X0158_M015A/       (negative)
#'                              plate1/
#'                              plate2/
#'                              ...
#'
#' @param samples_file_path
#'        directory of HABC mzML files
#'
#' @param rds_dir
#'        Directory for saving .RDS files (rds_dir/<file>)
#'
#' @param mzrolldb_dir
#'        Directory for saving mzrollDB files (mzrolldb_dir/<file>)
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing only post-lle IS compounds
#'        located at (lib_dir/is_lib_name)
#'
#' @param pre_lle_is_lib_name
#'        IS library containing only pre-lle IS compounds
#'        located at (lib_dir/pre_lle_is_lib_name)
#'
#' @param neg_lib_name
#'        HABC positive library name (lib_dir/pos_lib_name)
#'
#' @param pos_lib_name
#'        HABC positive library name (lib_dir/pos_lib_name)
#'
#' @param tg_lib_name
#'        HABC TG library name (lib_dir/tg_lib_name)
#'
#' @param neg_search_params
#'        list of search params. default behavior is habc::habc_neg_search_params()
#'
#' @param pos_search_params
#'        list of search params. Default behavior is habc::habc_pos_search_params()
#'
#' @param reprocess_all_batches
#'        If TRUE, overwrite any existing data by reprocessing plates
#'        If FALSE, do not reprocess plates if that plate's results .rds file already exists.
#'
#' @return
#' \code{
#' list("pos_is_results"=pos_is_results,
#' "pos_pre_lle_is_results"=pos_pre_lle_is_results,
#' "pos_search_results"=pos_search_results,
#' "pos_adduct_table"=pos_adduct_table,
#' "tg_results"=tg_results,
#' "neg_is_results"=neg_is_results,
#' "neg_pre_lle_is_results"=neg_pre_lle_is_results,
#' "neg_search_results"=neg_search_results,
#' "neg_adduct_table"=neg_adduct_table)
#' }
#'
#' \itemize{
#'    \item{\code{pos_is_results}}{: positive mode search results of post-LLE IS.}
#'    \item{\code{pos_pre_lle_is_results}}{: positive mode search results of pre-LLE IS.}
#'    \item{\code{pos_search_results}}{: positive mode search results of HABC compounds.}
#'    \item{\code{pos_adduct_table}}{: positive mode table of compound adduct MS1 scan intensities.}
#'    \item{\code{tg_results}}{: Positive mode MS3 search of TGs.}
#'    \item{\code{neg_is_results}}{: negativee mode search results of post-LLE IS.}
#'    \item{\code{neg_pre_lle_is_results}}{: negative mode search results of pre-LLE IS.}
#'    \item{\code{neg_search_results}}{: negative mode search results of HABC compounds.}
#'    \item{\code{neg_adduct_table}}{: negative mode table of compound adduct MS1 scan intensities.}
#' }
#'
#' @export
habc_process_all_plates_v1 <- function(
  top_level_samples_dir,
  rds_dir,
  mzrolldb_dir,
  lib_dir,
  is_lib_name,
  pre_lle_is_lib_name,
  neg_lib_name,
  pos_lib_name,
  tg_lib_name,
  neg_search_params = list(),
  pos_search_params = list(),
  reprocess_all_batches = F
) {
  pos_dir <- paste(top_level_samples_dir, "X0158_M014A", sep = "/")
  neg_dir <- paste(top_level_samples_dir, "X0158_M015A", sep = "/")

  pos_plate_names <- list.files(pos_dir)
  neg_plate_names <- list.files(neg_dir)

  # check that .rds files are present in the appropriate directory, or not
  # only process those that do not have a corresponding .rds file
  if (!reprocess_all_batches) {
    missing_pos_plate_names <- character(0)

    if (length(pos_plate_names) > 0) {
      for (x in 1:length(pos_plate_names)) {
        plate_name <- pos_plate_names[x]

        rds_file <- paste0(rds_dir, "/X0158_M014A_", plate_name, ".rds")

        if (!file.exists(rds_file)) {
          missing_pos_plate_names <- c(missing_pos_plate_names, plate_name)
        }
      }
    }

    missing_neg_plate_names <- character(0)

    if (length(neg_plate_names) > 0) {
      for (x in 1:length(neg_plate_names)) {
        plate_name <- neg_plate_names[x]

        rds_file <- paste0(rds_dir, "/X0158_M015A_", plate_name, ".rds")

        if (!file.exists(rds_file)) {
          missing_neg_plate_names <- c(missing_neg_plate_names, plate_name)
        }
      }
    }

    pos_plate_names <- missing_pos_plate_names
    neg_plate_names <- missing_neg_plate_names
  }

  # Process positive plates
  if (length(pos_plate_names) >= 1) {
    cat("\nStarting processing positive mode for plates.\n\n")

    for (x in 1:length(pos_plate_names)) {
      plate_name <- pos_plate_names[x]

      cat(paste0("\nStarting processing positive mode plate ", x, "/", length(pos_plate_names), ": ", plate_name, "\n\n"))

      samples_file_path <- paste(pos_dir, plate_name, sep = "/")

      habc_process_plate_pos_v1(
        plate_name = plate_name,
        samples_file_path = samples_file_path,
        rds_dir = rds_dir,
        mzrolldb_dir = mzrolldb_dir,
        lib_dir = lib_dir,
        is_lib_name = is_lib_name,
        pre_lle_is_lib_name = pre_lle_is_lib_name,
        pos_lib_name = pos_lib_name,
        tg_lib_name = tg_lib_name,
        pos_search_params = pos_search_params
      )

      cat(paste0("\nFinished processing postive mode plate ", x, "/", length(pos_plate_names), ": ", plate_name, "\n\n"))
    }
    cat("\nFinished processing positive mode plates.\n\n")
  }

  # Process negative plates
  if (length(neg_plate_names) >= 1) {
    cat("\nStarting processing negative mode plates.\n\n")

    for (x in 1:length(neg_plate_names)) {
      plate_name <- neg_plate_names[x]

      cat(paste0("\nStarting processing negative mode plate ", x, "/", length(neg_plate_names), ": ", plate_name, "\n\n"))

      samples_file_path <- paste(neg_dir, plate_name, sep = "/")

      habc_process_plate_neg_v1(
        plate_name = plate_name,
        samples_file_path = samples_file_path,
        rds_dir = rds_dir,
        mzrolldb_dir = mzrolldb_dir,
        lib_dir = lib_dir,
        is_lib_name = is_lib_name,
        pre_lle_is_lib_name = pre_lle_is_lib_name,
        neg_lib_name = neg_lib_name,
        neg_search_params = neg_search_params
      )

      cat(paste0("\nFinished processing negative mode plate ", x, "/", length(neg_plate_names), ": ", plate_name, "\n\n"))
    }
    cat("\nFinished processing negative mode plates.\n\n")
  }

  # Issue 571
  habc_get_rds_v1(rds_dir = rds_dir)
}

#' HABC process plate version 1 (negative mode data)
#'      Version 2 introduced September 2021
#'
#' Function uses habc_neg_search_params() for simple search
#'
#' @description
#' Process a HABC plate collected in negative mode and save results on file system
#'
#' @param plate_name
#'        name of folder containing mzML files associated with a single HABC plate
#'
#' @param samples_file_path
#'        directory of HABC mzML files
#'
#' @param rds_dir
#'        parent directory for saving .RDS files (rds_dir/plate_name/<file>)
#'
#' @param mzrolldb_dir
#'        parent directory for saving mzrollDB files (mzrolldb_dir/plate_name/<file>)
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing only post-lle IS compounds
#'        located at (lib_dir/is_lib_name)
#'
#' @param neg_lib_name
#'        HABC positive library name (lib_dir/neg_lib_name)
#'
#' @param save_mzrolldb_as_rds
#'        Avoid saving results into an mzrolldb to avoid server SQLite issues.
#'        if TRUE, saves results as rds that can be reinflated and run again.
#'
#' @export
habc_process_plate_neg_v2 <- function(
  plate_name,
  samples_file_path,
  rds_dir,
  mzrolldb_dir,
  lib_dir,
  is_lib_name,
  neg_lib_name,
  save_mzrolldb_as_rds = FALSE
) {
  # files ######################################
  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  habc_neg_lib_file <- paste(lib_dir, neg_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")
  rds_results_file <- paste0(rds_dir, "/X0158_M015A_", plate_name, ".rds")

  suffix <- ifelse(save_mzrolldb_as_rds == TRUE, ".rds", ".mzrollDB")
  mzrolldb_results_file <- paste0(mzrolldb_dir, "/X0158_M015A_", plate_name, suffix)

  # samples ####################################
  samples_info <- habc_samples_df(samples_file_path, FALSE)
  samples_df <- samples_info$samples_df %>% dplyr::arrange(desc(type), order_num)
  if (nrow(samples_df) == 0) {
    cat(paste0("habc_process_plate_neg_v2():\nNo samples detected for plate '", plate_name, "'\nbased on samples files path '", samples_files_path, "'"))
    return(invisible(0))
  }
  habc2_ms2_ranges <- samples_info$ms2_ranges
  all_bulkpool_samples <- samples_df %>% dplyr::filter(type == "BulkPool")

  # bulkpool IS search #########################
  is_lib_file <- file.path(lib_dir, is_lib_name)
  is_lib_neg <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl("\\-$", adductName))
  is_lib_neg_sliced <- mzkitcpp::DI_slice_library(habc2_ms2_ranges, is_lib_neg)

  habc2_bulkpool_IS_results <- mzkitcpp::DI_pipeline(
    samples = all_bulkpool_samples$file,
    ms2_ranges = habc2_ms2_ranges,
    is_sliced_lib = is_lib_neg_sliced,
    is_search_params = habc_is_search_params(),
    sliced_lib = is_lib_neg_sliced,
    search_params = habc_is_search_params(),
    adducts_file = adducts_file,
    debug = FALSE
  )

  # filter bulkpool samples ####################
  good_bulkpools <- habc_good_bulkpools(habc2_bulkpool_IS_results$adduct_table)
  good_bulkpool_all_data <- samples_df %>%
    dplyr::filter(sample_name %in% good_bulkpools$sample)
  good_bulkpool_samples <- good_bulkpool_all_data$file

  if (length(good_bulkpool_samples) < 3) {
    cat(paste0("Only found ", length(good_bulkpool_samples), " good bulkpool samples!\nAt least 3 are required.\nExiting without processing plate."))
    return(invisible(0))
  }

  # retain all samples except for bad bulkpool samples
  samples_df_filtered <- samples_df %>% dplyr::filter(type != "BulkPool" | file %in% good_bulkpool_samples)

  # bulkpool whitelist search ##################
  whitelist_search_parameters <- habc_neg_whitelist_params()

  habc_neg_lib <- mzkitcpp::import_msp_lipids_library(habc_neg_lib_file)
  neg_lib_full_sliced <- mzkitcpp::DI_slice_library(habc2_ms2_ranges, habc_neg_lib)

  IS_quant_table <- dplyr::inner_join(habc2_bulkpool_IS_results$IS_quant, samples_df, by = c("sample" = "sample_name")) %>%
    dplyr::filter(type == "BulkPool" & rank == 1)

  IS_quant_table_class_adduct <- IS_quant_table %>%
    dplyr::select(lipidClass, adductName) %>%
    unique()

  lib_top_adducts <- dplyr::inner_join(habc_neg_lib, IS_quant_table_class_adduct, by = c("lipidClass", "adductName")) %>%
    dplyr::filter(!grepl("_13C$", compoundName))

  lib_top_adducts_sliced <- mzkitcpp::DI_slice_library(habc2_ms2_ranges, lib_top_adducts)

  habc2_bulkpool_search_results <- mzkitcpp::DI_pipeline(
    samples = good_bulkpool_samples,
    ms2_ranges = habc2_ms2_ranges,
    is_sliced_lib = is_lib_neg_sliced,
    is_search_params = habc_is_search_params(),
    sliced_lib = lib_top_adducts_sliced,
    search_params = whitelist_search_parameters,
    adducts_file = adducts_file,
    debug = F
  )

  # plate-specific library ################################
  bulkpool_compounds_w_frequency <- habc2_bulkpool_search_results$search %>%
    dplyr::select(compoundName, adductName, sample) %>%
    unique() %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(num_bulkpool_samples = n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-num_bulkpool_samples, compoundName, adductName) %>%
    # May be missing in one bulkpool sample per plate, but not more than that
    # This assumes that samples have been adequately QCed, and any bad bulkpool
    # samples have already been removed.
    dplyr::filter(num_bulkpool_samples >= length(good_bulkpool_samples) - 1) %>%
    dplyr::select(compoundName, adductName, num_bulkpool_samples) %>%
    unique()

  habc2_filtered_lib_mono <- habc_neg_lib %>% dplyr::filter(compoundName %in% bulkpool_compounds_w_frequency$compoundName)

  habc2_filtered_lib <- rbind(habc2_filtered_lib_mono, to_13C_lib(habc2_filtered_lib_mono))

  habc2_filtered_lib_sliced <- mzkitcpp::DI_slice_library(habc2_ms2_ranges, habc2_filtered_lib)

  # habc search ################################
  habc2_search_results <- mzkitcpp::DI_pipeline(
    samples = samples_df_filtered$file,
    ms2_ranges = habc2_ms2_ranges,
    is_sliced_lib = is_lib_neg_sliced,
    is_search_params = habc_is_search_params(),
    sliced_lib = habc2_filtered_lib_sliced,
    search_params = habc_neg_search_params(version_num = 2),
    adducts_file = adducts_file,
    debug = F
  )

  # center results ################################
  habc2_centered_adduct_table <- habc_bulkpool_centered_adduct_table(habc2_search_results$adduct_table, samples_df)
  habc2_centered_search_table <- habc_bulkpool_centered_search_table(habc2_search_results$search, samples_df)

  # quant table ################################
  habc2_quant_table <- to_quant_table(habc2_centered_adduct_table, FALSE, TRUE) %>%
    dplyr::inner_join(., bulkpool_compounds_w_frequency, by = c("compoundName", "adductName")) %>%
    dplyr::inner_join(samples_df, by = c("sample" = "sample_name")) %>%
    dplyr::select(-file, -type_id)

  # save results ###############################
  system(glue::glue("rm {old_rds_results_file} 2>&1", old_rds_results_file = rds_results_file))
  saveRDS(habc2_quant_table, file = rds_results_file)

  system(glue::glue("rm {old_mzroll_db_file} 2>&1", old_mzroll_db_file = mzrolldb_results_file))

  color_table <- habc_type_color_table(samples_df, good_bulkpool_all_data$sample_name)

  encoded_quantType <- DIMS_encoded_quant_type()

  if (save_mzrolldb_as_rds) {
    mzrolldb_results <-
      list(
        "ms2_ranges" = habc2_ms2_ranges,
        "adducts_file" = adducts_file,
        "color_table" = color_table,
        "encoded_quantType" = encoded_quantType,
        "IS_samples" = all_bulkpool_samples$file,
        "IS_library_name" = is_lib_file,
        "IS_search_lib_sliced" = is_lib_neg_sliced,
        "IS_search_params" = habc_is_search_params(),
        "IS_search_results" = habc2_bulkpool_IS_results$search,
        "IS_adduct_table" = habc2_bulkpool_IS_results$adduct_table,
        "IS_set_name" = rep(plate_name, nrow(all_bulkpool_samples)),
        "HABC_samples" = samples_df_filtered$file,
        "HABC_library_name" = neg_lib_name,
        "HABC_search_lib_sliced" = neg_lib_full_sliced,
        "HABC_search_params" = whitelist_search_parameters,
        "HABC_search_results" = habc2_centered_search_table,
        "HABC_adduct_table" = habc2_centered_adduct_table,
        "HABC_set_name" = samples_df_filtered$plate
      )
    saveRDS(mzrolldb_results, file = mzrolldb_results_file)
  } else {
    # Add IS search results
    add_direct_infusion_search_results(
      mzroll_db_path = mzrolldb_results_file,
      samples = all_bulkpool_samples$file,
      ms2_ranges = habc2_ms2_ranges,
      library_name = is_lib_file,
      search_lib_sliced = is_lib_neg_sliced,
      search_params = habc_is_search_params(),
      adducts_file = adducts_file,
      di_search_results = habc2_bulkpool_IS_results$search,
      di_adduct_table = habc2_bulkpool_IS_results$adduct_table,
      set_name = rep(plate_name, nrow(all_bulkpool_samples)),
      color_table = color_table
    )

    # Add HABC neg search
    add_direct_infusion_search_results(
      mzroll_db_path = mzrolldb_results_file,
      samples = samples_df_filtered$file,
      ms2_ranges = habc2_ms2_ranges,
      library_name = neg_lib_name,
      search_lib_sliced = neg_lib_full_sliced,
      search_params = whitelist_search_parameters,
      adducts_file = adducts_file,
      di_search_results = habc2_centered_search_table,
      di_adduct_table = habc2_centered_adduct_table,
      set_name = samples_df_filtered$plate,
      color_table = color_table
    )

    mzroll_db_con <- DBI::dbConnect(RSQLite::SQLite(), mzrolldb_results_file)
    ui_options <- tibble::tibble(
      key = c("quantType"),
      value = c(encoded_quantType)
    )
    DBI::dbAppendTable(mzroll_db_con, "ui", ui_options)
    DBI::dbDisconnect(mzroll_db_con)
  }

  # end ########################################
  invisible(0)
}

#' HABC process negative plates v2
#'
#' @description
#' Process all HABC plate data found on file system
#'
#' @param top_level_samples_dir
#'        HABC samples, organized into sub-folders by plate name.
#'
#'        Require this arrangement on file system:
#'
#'        top_level_samples_dir/
#'                      X0158_M015A/       (negative)
#'                              plate1/
#'                              plate2/
#'                              ...
#'
#' @param plates
#'        vector of plates to run
#'
#' @param rds_dir
#'        Directory for saving .RDS files (rds_dir/<file>)
#'
#' @param mzrolldb_dir
#'        Directory for saving mzrollDB files (mzrolldb_dir/<file>)
#'
#' @param lib_dir
#'        directory containing .msp spectral library files
#'
#' @param is_lib_name
#'        IS library containing only post-lle IS compounds
#'        located at (lib_dir/is_lib_name)
#'
#' @param neg_lib_name
#'        HABC positive library name (lib_dir/pos_lib_name)
#'
#' @param save_mzrolldb_as_rds
#'        Avoid saving results into an mzrolldb to avoid server SQLite issues.
#'        if TRUE, saves results as rds that can be reinflated and run again.
#'
#' @export
habc_process_neg_plates_v2 <- function(
  top_level_samples_dir,
  plates,
  rds_dir,
  mzrolldb_dir,
  lib_dir,
  is_lib_name,
  neg_lib_name,
  save_mzrolldb_as_rds = FALSE
) {
  for (i in 1:length(plates)) {
    plate_name <- plates[i]
    cat(paste0("Starting analysis of plate '", plate_name, "'\n"))

    samples_file_path <- file.path(top_level_samples_dir, "X0158_M015A", plate_name)

    cat(paste0("Samples for '", plate_name, "' retrieved from ", samples_file_path, "\n"))

    habc_process_plate_neg_v2(
      plate_name = plate_name,
      samples_file_path = samples_file_path,
      rds_dir = rds_dir,
      mzrolldb_dir = mzrolldb_dir,
      lib_dir = lib_dir,
      is_lib_name = is_lib_name,
      neg_lib_name = neg_lib_name,
      save_mzrolldb_as_rds = save_mzrolldb_as_rds
    )

    cat(paste0("Completed plate '", plate_name, "'.\n\n"))
  }

  cat(paste0("Successfully Completed All Plates.\n"))

  invisible(0)
}

#' HABC v3 stage 1: bulkpool samples analysis
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
#' Two kinds of files are created in this directory: "samples" and "compounds".
#' Files are named as <plate>_<method>_<samples|compounds>.rds
#' Each file is a tibble, with information from the DIMS searches that can be later be
#' used to create a bulkpool library (in stage 2)
#'
#' @export
habc_v3_stage1_process_plate <- function(
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
  is_lib_file <- paste(lib_dir, is_lib_name, sep = "/")
  habc_lib_file <- paste(lib_dir, habc_lib_name, sep = "/")
  adducts_file <- paste(lib_dir, "ADDUCTS.csv", sep = "/")
  is_ms3 <- !is.null(ms3_lib_name)
  method_name <- ifelse(is_ms3, "X0158_M014A", "X0158_M015A")
  samples_output_file <- file.path(stage_1_results_dir, paste0(plate_name, "_", method_name, "_samples.rds"))
  compounds_output_file <- file.path(stage_1_results_dir, paste0(plate_name, "_", method_name, "_compounds.rds"))

  # samples ####################################
  samples_info <- habc_samples_df(samples_file_path, is_ms3)
  samples_df <- samples_info$samples_df %>% dplyr::arrange(desc(type), order_num)
  habc3_ms2_ranges <- samples_info$ms2_ranges
  all_bulkpool_samples <- samples_df %>% dplyr::filter(type == "BulkPool")

  # libraries ###################################
  lib_subset_string <- ifelse(is_ms3, "\\+$", "\\-$")
  is_lib <- mzkitcpp::import_msp_lipids_library(is_lib_file) %>% dplyr::filter(grepl(lib_subset_string, adductName))
  is_lib_sliced <- mzkitcpp::DI_slice_library(habc3_ms2_ranges, is_lib)
  habc_lib <- mzkitcpp::import_msp_lipids_library(habc_lib_file)

  # good bulkpool ms2 ###################################
  habc3_bulkpool_ms2_results <- mzkitcpp::DI_pipeline(
    samples = all_bulkpool_samples$file,
    ms2_ranges = habc3_ms2_ranges,
    is_sliced_lib = is_lib_sliced,
    is_search_params = is_search_params,
    sliced_lib = is_lib_sliced,
    search_params = is_search_params,
    adducts_file = adducts_file,
    debug = FALSE
  )

  habc3_bulkpool_ms2_good <- habc_good_bulkpools(
    habc3_bulkpool_ms2_results$adduct_table,
    bulkpool_params$bulkpoolGoodMinFracIntensity,
    bulkpool_params$bulkpoolGoodMinFracDetected
  )

  habc3_bulkpool_ms2_good_samples <- habc3_bulkpool_ms2_good$sample

  # good bulkpool ms3 ###################################
  if (is_ms3) {
    habc3_bulkpool_ms3_results <- mzkitcpp::DI_pipeline_ms3_search(
      samples = all_bulkpool_samples$file,
      is_lib = default_tg_is_ms3,
      is_search_params = is_search_params,
      search_lib = default_tg_is_ms3,
      search_params = whitelist_ms3_search_params,
      adducts_file = adducts_file,
      debug = F
    )

    # filtering criteria: ms3 m/z matches and sum ms3 intensity
    habc3_bulkpool_ms3_good <- habc3_bulkpool_ms3_results %>%
      dplyr::select(sample, fragmentLabel, num_ms3_mz_matches, ms3_intensity_sum) %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(total_ms3_mz_matches = sum(num_ms3_mz_matches)) %>%
      dplyr::ungroup() %>%
      dplyr::select(sample, total_ms3_mz_matches, ms3_intensity_sum) %>%
      unique() %>%
      dplyr::filter(total_ms3_mz_matches >= bulkpool_params$bulkpoolGoodMs3MinTotalMatches &
        ms3_intensity_sum >= bulkpool_params$bulkpoolGoodMs3MinIntensitySum)

    habc3_bulkpool_ms3_good_samples <- habc3_bulkpool_ms3_good$sample
  } else {
    habc3_bulkpool_ms3_good_samples <- c()
  }

  # bulkpool sample results ###################################
  bulkpool_results_table <- all_bulkpool_samples %>%
    dplyr::mutate(
      is_good_ms2 = sample_name %in% habc3_bulkpool_ms2_good_samples,
      is_good_ms3 = sample_name %in% habc3_bulkpool_ms3_good_samples
    ) %>%
    dplyr::select(file, sample_name, method, mode, plate, well_position, is_good_ms2, is_good_ms3)

  # whitelist ms2 search ###################################

  good_ms2_samples <- bulkpool_results_table %>% dplyr::filter(is_good_ms2)

  class_adduct_counts_table <- whitelist_ms2_search_params$ms1IonList %>%
    dplyr::group_by(lipidClass) %>%
    dplyr::summarize(num_adducts = n())

  ion_list <- whitelist_ms2_search_params$ms1IonList %>% dplyr::select(lipidClass, adductName)

  habc_whitelist_ms2_lib_ion_list <- habc_lib %>%
    dplyr::inner_join(., ion_list, by = c("lipidClass", "adductName"))

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

  habc_whitelist_ms2_lib_sliced <- mzkitcpp::DI_slice_library(habc3_ms2_ranges, habc_whitelist_ms2_lib)

  whitelist_bulkpool_ms2_search <- mzkitcpp::DI_pipeline(
    samples = good_ms2_samples$file,
    ms2_ranges = habc3_ms2_ranges,
    is_sliced_lib = is_lib_sliced,
    is_search_params = is_search_params,
    sliced_lib = habc_whitelist_ms2_lib_sliced,
    search_params = whitelist_ms2_search_params,
    adducts_file = adducts_file,
    debug = FALSE
  )

  whitelist_bulkpool_ms2_compounds <- whitelist_bulkpool_ms2_search$search %>%
    dplyr::inner_join(all_bulkpool_samples, by = c("sample" = "sample_name")) %>%
    dplyr::rename(sample_name = sample) %>%
    dplyr::select(file, sample_name, method, mode, plate, well_position, lipidClass, compoundName, adductName) %>%
    unique() %>%
    dplyr::mutate(is_ms2_compound = TRUE, is_ms3_compound = FALSE) %>%
    dplyr::arrange(sample_name, compoundName, adductName)

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
      to_ms3_lib() %>%
      dplyr::filter(prec_mzs %in% ms3_targets$prec_mzs) %>%
      to_ms2_lib()

    whitelist_bulkpool_ms3_search <- mzkitcpp::DI_pipeline_ms3_search(
      samples = good_ms3_samples$file,
      is_lib = default_tg_is_ms3,
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
      dplyr::mutate(is_ms2_compound = FALSE, is_ms3_compound = TRUE) %>%
      dplyr::arrange(sample_name, compoundName, adductName)

    whitelist_bulkpool_compounds <- rbind(whitelist_bulkpool_ms2_compounds, whitelist_bulkpool_ms3_compounds)
  } else {
    whitelist_bulkpool_compounds <- whitelist_bulkpool_ms2_compounds
  }

  # save results ###################################
  saveRDS(bulkpool_results_table, samples_output_file)
  saveRDS(whitelist_bulkpool_compounds, compounds_output_file)

  invisible(0)
}

#' HABC v3 stage 2: compound whitelist and plate validation
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
#'        .rds File containing "good" designation for all bulkpool samples
#'        saved as tibble
#'
#' @param compound_list_file
#'        .rds File containing list of compound names to retain in regular searches
#'
#' @export
habc_v3_stage2_build_library <- function(
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

  num_valid_neg_ms2 <- mode_good_counts[[1, 2]] # (neg, ms2)
  num_valid_pos_ms2 <- mode_good_counts[[2, 2]] # (pos, ms2)
  num_valid_pos_ms3 <- mode_good_counts[[2, 3]] # (pos, ms3)

  # build compound list ###################################
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

  # filter out invalid plates
  bulkpool_compound_data_plate_filtered <- bulkpool_compound_data %>%
    dplyr::inner_join(., plate_validation_tibble, by = c("plate", "mode")) %>%
    dplyr::filter((is_ms2_compound & is_plate_good_ms2) | (is_ms3_compound & is_plate_good_ms3)) %>%
    dplyr::select(-is_plate_good_ms2, -is_plate_good_ms3, -file)

  # filter based on bulkpool parameters
  all_compound_ions <- rbind(bulkpool_params$ms1IonList, bulkpool_params$ms3IonList)

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
    dplyr::mutate(is_m_plus_two = grepl("_13C_13C$", compoundName)) %>%
    dplyr::mutate(is_m_plus_one = !is_m_plus_two & grepl("_13C$", compoundName)) %>%
    # TODO: should [M+1] and [M+2] have different ratios than major/minor adducts?

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

  bulkpool_isotope_filtered <- bulkpool_adduct_filtered %>%
    dplyr::filter(is_ms2_compound == TRUE) %>%
    dplyr::select(lipidClass, compoundName, adductName, is_m_plus_one, is_m_plus_two) %>%
    unique() %>%
    dplyr::mutate(monoCompoundName = gsub("_13C", "", compoundName)) %>%
    dplyr::group_by(monoCompoundName, adductName) %>%
    dplyr::mutate(is_has_m_plus_one = any(is_m_plus_one), is_has_m_plus_two = any(is_m_plus_two)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(isotope_constraints, by = c("lipidClass", "adductName")) %>%
    dplyr::select(
      lipidClass, monoCompoundName, adductName,
      is_has_m_plus_one, is_has_m_plus_two, ms1IsRequireMPlusOne, ms1IsRequireMPlusTwo
    ) %>%
    unique() %>%
    dplyr::rename(compoundName = monoCompoundName) %>%
    dplyr::filter(!ms1IsRequireMPlusOne | (ms1IsRequireMPlusOne & is_has_m_plus_one)) %>%
    dplyr::filter(!ms1IsRequireMPlusTwo | (ms1IsRequireMPlusTwo & is_has_m_plus_two)) %>%
    dplyr::select(lipidClass, compoundName, adductName) %>%
    unique()

  # filter out compounds where all required adduct forms are not found
  num_adducts_constraints <- bulkpool_params$ms1IonList %>%
    dplyr::mutate(mode = ifelse(grepl("\\+$", adductName), "pos", "neg")) %>%
    dplyr::group_by(lipidClass, mode) %>%
    dplyr::mutate(min_num_adducts = n()) %>%
    ungroup() %>%
    dplyr::select(lipidClass, mode, min_num_adducts) %>%
    unique()

  bulkpool_ms2_compounds <- bulkpool_isotope_filtered %>%
    dplyr::mutate(mode = ifelse(grepl("\\+$", adductName), "pos", "neg")) %>%
    dplyr::inner_join(num_adducts_constraints, by = c("lipidClass", "mode")) %>%
    dplyr::group_by(compoundName, mode) %>%
    dplyr::mutate(num_adducts = n()) %>%
    dplyr::ungroup() %>%
    unique() %>%
    dplyr::filter(num_adducts >= min_num_adducts) %>%
    dplyr::select(compoundName) %>%
    unique()

  bulkpool_compounds <- c(bulkpool_ms2_compounds$compoundName, bulkpool_ms3_compounds$compoundName)

  # save results ###################################
  saveRDS(plate_validation_tibble, plate_validation_file)
  saveRDS(bulkpool_sample_data, bulkpool_sample_validation_file)
  saveRDS(bulkpool_compounds, compound_list_file)

  invisible(0)
}

#' HABC v3 stage 3: bulkpool samples analysis
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
#' Two kinds of files are created in this directory: "samples" and "compounds".
#' Files are named as <plate>_<method>_<samples|compounds>.rds
#' Each file is a tibble, with information from the DIMS searches that can be later be
#' used to create a bulkpool library (in stage 2)
#'
#' @export
habc_v3_stage3_process_plate <- function(
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
  rds_results_file <- paste0(stage_3_rds_dir, "/", method_name, "_", plate_name, ".rds")
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
  habc3_ms2_ranges <- samples_info$ms2_ranges

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
    is_lib_sliced <- mzkitcpp::DI_slice_library(habc3_ms2_ranges, is_lib)

    # IS search
    habc3_IS_results <- mzkitcpp::DI_pipeline(
      samples = samples_df$file,
      ms2_ranges = habc3_ms2_ranges,
      is_sliced_lib = is_lib_sliced,
      is_search_params = is_search_params,
      sliced_lib = is_lib_sliced,
      search_params = is_search_params,
      adducts_file = adducts_file,
      debug = FALSE
    )

    # IS quant ions
    is_quant_table <- to_quant_table(habc3_IS_results$adduct_table, FALSE, TRUE) %>%
      dplyr::select(sample, compoundName, adductName, ms1_intensity, ms2_intensity)

    is_search_table_subset <- habc3_IS_results$search %>%
      dplyr::select(sample, compoundName, adductName, ms1_intensity_is_nearest_scan_normalized) %>%
      unique() %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized))

    is_adduct_table_subset <- habc3_IS_results$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::select(sample, compoundName, adductName, ms2_diagnostic_norm_intensity) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity)) %>%
      dplyr::rename(diagnostic_ms2_intensity = ms2_diagnostic_norm_intensity)

    is_quant_ions <- is_quant_table %>%
      dplyr::full_join(is_search_table_subset, by = c("compoundName", "adductName", "sample")) %>%
      dplyr::full_join(is_adduct_table_subset, by = c("compoundName", "adductName", "sample"))

    # MS2 library
    habc_lib <- mzkitcpp::import_msp_lipids_library(habc_lib_file) %>%
      dplyr::mutate(monoCompoundName = gsub("_13C", "", compoundName)) %>%
      dplyr::filter(monoCompoundName %in% bulkpool_compounds) %>%
      dplyr::select(-monoCompoundName)

    habc_lib_sliced <- mzkitcpp::DI_slice_library(habc3_ms2_ranges, habc_lib)

    # MS2 biological search
    habc3_biological_search <- mzkitcpp::DI_pipeline(
      samples = samples_df_ms2_filtered$file,
      ms2_ranges = habc3_ms2_ranges,
      is_sliced_lib = is_lib_sliced,
      is_search_params = is_search_params,
      sliced_lib = habc_lib_sliced,
      search_params = biological_ms2_search_params,
      adducts_file = adducts_file,
      debug = FALSE
    )

    nearest_scan_IS <- habc3_biological_search$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized) &
        ms1_intensity_is_nearest_scan_normalized > 0 &
        ms1_intensity_is_nearest_scan_normalized < Inf)

    centered_nearest_scan_IS <- habc_bulkpool_centered_quant_ion_table(nearest_scan_IS, "ms1_intensity_is_nearest_scan_normalized") %>%
      dplyr::select(compoundName, adductName, sample, ms1_intensity_is_nearest_scan_normalized) %>%
      dplyr::filter(!is.na(ms1_intensity_is_nearest_scan_normalized))

    diagnostic_ms2 <- habc3_biological_search$adduct_table %>%
      dplyr::filter(is_identified == TRUE) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity) &
        ms2_diagnostic_norm_intensity > 0 &
        ms2_diagnostic_norm_intensity < Inf)

    centered_diagnostic_ms2 <- habc_bulkpool_centered_quant_ion_table(diagnostic_ms2, "ms2_diagnostic_norm_intensity") %>%
      dplyr::select(compoundName, adductName, sample, ms2_diagnostic_norm_intensity) %>%
      dplyr::filter(!is.na(ms2_diagnostic_norm_intensity))

    # center after determination of preferred quant type to ensure accurate centering values
    habc3_quant_table <- to_quant_table(habc3_biological_search$adduct_table, FALSE, TRUE)

    # MS1
    habc3_quant_ions_ms1 <- habc3_quant_table %>%
      dplyr::select(-ms2_intensity, -ms2_quant_type) %>%
      dplyr::filter(!is.na(ms1_intensity) & ms1_intensity > 0 & ms1_intensity < Inf) %>%
      dplyr::mutate(ms1_quant_type = ifelse(grepl("^ms1", ms1_quant_type), ms1_quant_type, paste0("ms1_", ms1_quant_type))) %>%
      dplyr::rename(quant_type = ms1_quant_type, intensity = ms1_intensity) %>%
      dplyr::mutate(quant_class = "ms1")

    # MS2
    habc3_quant_ions_ms2 <- habc3_quant_table %>%
      dplyr::select(-ms1_intensity, -ms1_quant_type) %>%
      dplyr::filter(!is.na(ms2_intensity) & ms2_intensity > 0 & ms2_intensity < Inf) %>%
      dplyr::mutate(ms2_quant_type = ifelse(grepl("^ms2", ms2_quant_type), ms2_quant_type, paste0("ms2_", ms2_quant_type))) %>%
      dplyr::rename(quant_type = ms2_quant_type, intensity = ms2_intensity) %>%
      dplyr::mutate(quant_class = "ms2")

    # Center quant ions
    habc3_centered_quant_ions_ms1 <- habc_bulkpool_centered_quant_ion_table(habc3_quant_ions_ms1, "intensity")
    habc3_centered_quant_ions_ms2 <- habc_bulkpool_centered_quant_ion_table(habc3_quant_ions_ms2, "intensity")

    di_ms1_quant_ions_condensed <- habc3_centered_quant_ions_ms1 %>%
      dplyr::rename(ms1_intensity = intensity) %>%
      dplyr::select(sample, compoundName, adductName, ms1_intensity) %>%
      dplyr::filter(!is.na(ms1_intensity))

    di_ms2_quant_ions_condensed <- habc3_centered_quant_ions_ms2 %>%
      dplyr::rename(ms2_intensity = intensity) %>%
      dplyr::select(sample, compoundName, adductName, ms2_intensity) %>%
      dplyr::filter(!is.na(ms2_intensity))

    # formatted for mzrollDB
    habc3_quant_ions <- di_ms1_quant_ions_condensed %>%
      dplyr::full_join(di_ms2_quant_ions_condensed, by = c("sample", "compoundName", "adductName")) %>%
      dplyr::full_join(centered_nearest_scan_IS, by = c("sample", "compoundName", "adductName")) %>%
      dplyr::full_join(centered_diagnostic_ms2, by = c("sample", "compoundName", "adductName")) %>%
      dplyr::rename(diagnostic_ms2_intensity = ms2_diagnostic_norm_intensity)

    # formatted for stage 4 processing
    habc_centered_quant_ions <- rbind(habc3_centered_quant_ions_ms1, habc3_centered_quant_ions_ms2) %>%
      dplyr::inner_join(samples_df, by = c("sample" = "sample_name")) %>%
      dplyr::select(-file, -type_id) %>%
      unique()
  } else {
    is_quant_ions <- NULL
    is_lib_sliced <- NULL
    habc3_IS_results <- NULL
    habc_lib_sliced <- NULL
    habc_centered_quant_ions <- NULL
    habc3_quant_ions <- NULL
    habc3_quant_table <- NULL
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
      to_ms3_lib() %>%
      dplyr::filter(prec_mzs %in% ms3_targets$prec_mzs) %>%
      to_ms2_lib() %>%
      dplyr::filter(compoundName %in% bulkpool_compounds) %>%
      dplyr::mutate(ms2_intensity = 1) %>%
      dplyr::select(colnames(habc_ms3_lib))

    tg_IS_lib <- default_tg_is_ms3 %>%
      dplyr::mutate(ms2_intensity = 1) %>%
      dplyr::select(colnames(habc_ms3_lib))

    habc_ms3_lib_filtered_w_IS <- rbind(habc_ms3_lib_filtered, tg_IS_lib)

    # MS3 biological search
    habc3_biological_ms3_search <- mzkitcpp::DI_pipeline_ms3_search(
      samples = samples_df_ms3_filtered$file,
      is_lib = default_tg_is_ms3,
      is_search_params = biological_ms3_search_params,
      search_lib = habc_ms3_lib_filtered_w_IS,
      search_params = biological_ms3_search_params,
      adducts_file = adducts_file,
      debug = F
    )

    # center results
    habc3_biological_ms3_median_column <- habc3_biological_ms3_search %>%
      dplyr::filter(!is.na(ms3_intensity_sum_norm) & ms3_intensity_sum_norm > 0 & ms3_intensity_sum_norm < Inf) %>%
      dplyr::select(sample, compoundName, adductName, ms3_intensity_sum_norm) %>%
      unique() %>%
      dplyr::group_by(compoundName, adductName) %>%
      dplyr::mutate(intensity_median = median(ms3_intensity_sum_norm, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(compoundName, adductName, intensity_median)

    habc3_biological_ms3_search_centered <- dplyr::inner_join(
      habc3_biological_ms3_search, habc3_biological_ms3_median_column,
      by = c("compoundName", "adductName")
    ) %>%
      dplyr::mutate(ms3_intensity_sum_norm = ms3_intensity_sum_norm / intensity_median) %>%
      unique()

    # quant table
    habc3_ms3_quant_table <- habc3_biological_ms3_search_centered %>%
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
    habc3_ms3_quant_table <- NULL
  }

  # color table ###################################
  color_table <- habc_type_color_table(samples_df, good_ms2_bulkpools$sample_name)
  plate_name_vector <- rep(plate_name, nrow(samples_df))

  # save RDS results ###################################
  system(glue::glue("rm {old_rds_results_file} 2>&1", old_rds_results_file = rds_results_file))

  rds_results <- NULL
  if (!is.null(habc_centered_quant_ions) && !is.null(habc3_ms3_quant_table)) {
    rds_results <- rbind(habc_centered_quant_ions, habc3_ms3_quant_table)
  } else if (!is.null(habc_centered_quant_ions) && is.null(habc3_ms3_quant_table)) {
    rds_results <- habc_centered_quant_ions
  } else if (is.null(habc_centered_quant_ions) && !is.null(habc3_ms3_quant_table)) {
    rds_results <- habc3_ms3_quant_table
  }

  saveRDS(rds_results, file = rds_results_file)

  # save mzrollDB results ###################################
  system(glue::glue("rm {old_mzroll_db_file} 2>&1", old_mzroll_db_file = mzrolldb_results_file))

  encoded_quantType <- DIMS_encoded_quant_type()

  if (save_mzrolldb_as_rds) {
    IS_search_results <- "NULL"
    IS_adduct_table <- "NULL"
    IS_quant_ions <- "NULL"
    if (!is.null(habc3_IS_results)) {
      IS_search_results <- habc3_IS_results$search
      IS_adduct_table <- habc3_IS_results$adduct_table
      IS_quant_ions <- is_quant_ions
    }

    HABC_search_results <- "NULL"
    HABC_adduct_table <- "NULL"
    HABC_quant_ions <- "NULL"
    if (!is.null(habc3_quant_table)) {
      HABC_search_results <- habc3_biological_search$search
      HABC_adduct_table <- habc3_biological_search$adduct_table
      HABC_quant_ions <- habc3_quant_ions
    }

    MS3_search_results <- "NULL"
    if (!is.null(habc3_ms3_quant_table)) {
      MS3_search_results <- habc3_biological_ms3_search_centered
    }

    mzrolldb_results <-
      list(
        "ms2_ranges" = habc3_ms2_ranges,
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
    if (!is.null(habc3_IS_results)) {
      # Add IS search results
      add_direct_infusion_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df$file,
        ms2_ranges = habc3_ms2_ranges,
        library_name = is_lib_file,
        search_lib_sliced = is_lib_sliced,
        search_params = is_search_params,
        adducts_file = adducts_file,
        di_search_results = habc3_IS_results$search,
        di_quant_ions = is_quant_ions,
        set_name = plate_name_vector,
        color_table = color_table
      )
    }

    if (!is.null(habc3_quant_table)) {
      # Add biological search
      add_direct_infusion_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df_ms2_filtered$file,
        ms2_ranges = habc3_ms2_ranges,
        library_name = habc_lib_name,
        search_lib_sliced = habc_lib_sliced,
        search_params = biological_ms2_search_params,
        adducts_file = adducts_file,
        di_search_results = habc3_centered_search_table,
        di_quant_ions = habc3_quant_ions,
        set_name = plate_name_vector,
        color_table = color_table
      )
    }

    # Add MS3 search

    if (!is.null(habc3_ms3_quant_table)) {
      add_targeted_ms3_search_results(
        mzroll_db_path = mzrolldb_results_file,
        samples = samples_df_ms3_filtered$file,
        ms3_search_results = habc3_biological_ms3_search_centered,
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

#' Convert RDS to mzrolldb file
#'
#' @description
#' convert RDS files created using habc_v3_stage3_process_plate() into mzrolldB files.
#'
#' @param rds_results_file full path of RDS results file
#' @param mzrolldb_results_file full path of target mzrolldb file
#' @param labels table of (compoundName, adductName, label) column to add label information to search results
#'
#' @export
habc_v3_rds_to_mzrolldb <- function(
  rds_results_file,
  mzrolldb_results_file,
  labels = NULL
) {
  mzrolldb_results <- readRDS(rds_results_file)
  system(glue::glue("rm {old_mzroll_db_file} 2>&1", old_mzroll_db_file = mzrolldb_results_file))

  # Add IS search results
  if (class(mzrolldb_results$IS_search_results) != "character") {
    add_direct_infusion_search_results(
      mzroll_db_path = mzrolldb_results_file,
      samples = mzrolldb_results$IS_samples,
      ms2_ranges = mzrolldb_results$ms2_ranges,
      library_name = mzrolldb_results$IS_library_name,
      search_lib_sliced = mzrolldb_results$IS_search_lib_sliced,
      search_params = mzrolldb_results$IS_search_params,
      adducts_file = mzrolldb_results$adducts_file,
      di_search_results = mzrolldb_results$IS_search_results,
      di_quant_ions = mzrolldb_results$IS_quant_ions,
      set_name = mzrolldb_results$IS_set_name,
      color_table = mzrolldb_results$color_table
    )
  }

  # Add biological search
  if (class(mzrolldb_results$HABC_search_results) != "character") {
    if (!is.null(labels)) {
      labeled_search_table <- mzrolldb_results$HABC_search_results %>%
        dplyr::left_join(labels, by = c("compoundName", "adductName")) %>%
        dplyr::mutate(label = ifelse(is.na(label), "", label))
    } else {
      labeled_search_table <- mzrolldb_results$HABC_search_results
    }

    add_direct_infusion_search_results(
      mzroll_db_path = mzrolldb_results_file,
      samples = mzrolldb_results$HABC_samples,
      ms2_ranges = mzrolldb_results$ms2_ranges,
      library_name = mzrolldb_results$HABC_library_name,
      search_lib_sliced = mzrolldb_results$HABC_search_lib_sliced,
      search_params = mzrolldb_results$HABC_search_params,
      adducts_file = mzrolldb_results$adducts_file,
      di_search_results = labeled_search_table,
      di_quant_ions = mzrolldb_results$HABC_quant_ions,
      set_name = mzrolldb_results$HABC_set_name,
      color_table = mzrolldb_results$color_table
    )
  }

  # Add MS3 search
  if (class(mzrolldb_results$MS3_search_results) != "character") {
    add_targeted_ms3_search_results(
      mzroll_db_path = mzrolldb_results_file,
      samples = mzrolldb_results$MS3_samples,
      ms3_search_results = mzrolldb_results$MS3_search_results,
      library_name = mzrolldb_results$MS3_library_name,
      search_params = mzrolldb_results$MS3_search_params,
      set_name = mzrolldb_results$MS3_set_name,
      color_table = mzrolldb_results$color_table
    )
  }

  mzroll_db_con <- DBI::dbConnect(RSQLite::SQLite(), mzrolldb_results_file)
  ui_options <- tibble::tibble(
    key = c("quantType"),
    value = c(mzrolldb_results$encoded_quantType)
  )
  DBI::dbAppendTable(mzroll_db_con, "ui", ui_options)
  DBI::dbDisconnect(mzroll_db_con)

  invisible(0)
}

#' Return table of quantitative measurements from collection of results.
#'
#' @description
#' combine stage 3 RDS result files into final table of quant results.
#'
#' @param stage_3_rds_dir collection of saved RDS results files.
#' @param bulkpool_params collection of params associated with bulkpool search.  Used here for quant ions.
#'
#' @export
habc_v3_stage4_formatting <- function(
  stage_3_rds_dir,
  bulkpool_params
) {
  # v3:  major ions are quant ions

  # quant ions ###############################
  quant_compound_ions <- rbind(bulkpool_params$ms1IonList, bulkpool_params$ms3IonList) %>%
    dplyr::filter(ms1IsMajorIon == TRUE) %>%
    dplyr::select(-ms1IsMajorIon)

  # import results ###############################
  all_stage3_results <- list.files(stage_3_rds_dir, pattern = "*.rds$", full.names = TRUE)

  stage3_quant_ions <- NULL
  stage3_quant_cvs <- NULL
  for (i in 1:length(all_stage3_results)) {
    plate_quant_data <- readRDS(all_stage3_results[i]) %>%
      dplyr::inner_join(quant_compound_ions, by = c("lipidClass", "adductName")) %>%
      dplyr::filter(!grepl("_13C", compoundName) & !grepl("_IS", compoundName))

    plate_mode <- plate_quant_data$mode[1]
    plate_name <- plate_quant_data$plate[1]

    plate_quant_cvs <- habc_cv_comparison_v3(plate_quant_data) %>%
      dplyr::mutate(mode = plate_mode, plate = plate_name)

    if (is.null(stage3_quant_ions)) {
      stage3_quant_ions <- plate_quant_data
      stage3_quant_cvs <- plate_quant_cvs
    } else {
      stage3_quant_ions <- rbind(stage3_quant_ions, plate_quant_data)
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
      dplyr::mutate(median_BulkPool_CV = median(BulkPool_CV)) %>%
      dplyr::ungroup()

    if (is.null(BulkPool_CVs)) {
      BulkPool_CVs <- ms_level_BulkPool_CVs
    } else {
      BulkPool_CVs <- rbind(BulkPool_CVs, ms_level_BulkPool_CVs)
    }
  }

  # Median BulkPool CV ###############################
  stage3_median_bulkpool_CVs <- BulkPool_CVs %>%
    dplyr::select(compoundName, adductName, mode, quant_class, median_BulkPool_CV) %>%
    unique() %>%
    dplyr::arrange(compoundName, adductName, mode) %>%
    dplyr::group_by(compoundName) %>%
    dplyr::mutate(is_lowest_CV = median_BulkPool_CV == min(median_BulkPool_CV)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(compoundName, mode) %>%
    dplyr::mutate(is_lowest_mode_CV = median_BulkPool_CV == min(median_BulkPool_CV)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(compoundName, quant_class) %>%
    dplyr::mutate(is_lowest_MS_CV = median_BulkPool_CV == min(median_BulkPool_CV)) %>%
    dplyr::ungroup()

  quantified_compounds <- stage3_quant_ions %>%
    dplyr::inner_join(stage3_median_bulkpool_CVs, by = c("compoundName", "adductName", "mode", "quant_class"))
}
