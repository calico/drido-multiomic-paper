###
### Data cleaning and normalization for DO-CR lipidomics
### Requires mzkitcpp
### Johanna Fleischman, Calico Life Sciences 2025
###
###
### Install libraries and source functions
suppressMessages(library(tidyverse))

# set
# local_filepath <- "~/drido-multiomic-paper"

source(file.path(local_filepath, "R/dims_normalization.R"))
source(file.path(local_filepath, "R/habc_lipids.R"))
source(file.path(local_filepath, "R/habc4_lipids.R"))
source(file.path(local_filepath, "R/docr_lipids.R"))
source(file.path(local_filepath, "R/direct_infusion.R"))

# Inputs for all processes
top_level_dir <- file.path(local_filepath, "data")
load(file.path(top_level_dir, "default_tg_is_ms3.rda"))
bulkpool_compounds <- readRDS(file.path(top_level_dir, "20230213-bulkpool-compounds.rds"))
is_search_params <- DIMS_read_search_params(file.path(top_level_dir, "habc_v4_IS.json"))
bulkpool_params <- DIMS_read_search_params(file.path(top_level_dir, "habc_v4_bulkpool_20250116.json"))
is_lib_name <- "20221020-Calico-Lipids-IS-dims_norm.msp"
rds_output_dir <- file.path(top_level_dir, "rds_output_2") # make
mzrolldb_output_dir <- file.path(top_level_dir, "mzrolldb_output_dir_2") # make
save_mzrolldb_as_rds <- TRUE
save_directory <- "~/"

# Inputs for positive mode
# positive_sample_plates <- "~/*.mzML" # mzML files
biological_lib_name <- "do_cr_pos_reduced_20230211.msp"
biological_ms2_search_params <- DIMS_read_search_params(
  file.path(top_level_dir, "habc_v4_biological_pos.json")
)
biological_ms3_search_params <- DIMS_read_search_params(
  file.path(top_level_dir, "habc_v4_biological_ms3_20250116.json")
)
is_ms3 <- TRUE
ms3_lib_name <- "20230213-DOCR-tg.msp"

for (p in list.files(positive_sample_plates)) {
  print("=======================================")
  print(paste0("Processing Positive Mode Plate ", p))

  plate_name <- sub("^(.)", "\\L\\1", p, perl = TRUE)
  samples_file_path <- file.path(positive_sample_plates, p)
  print("")

  docr_process_dims_plate(
    plate_name,
    samples_file_path,
    top_level_dir,
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
    is_ms3,
    save_mzrolldb_as_rds,
    default_tg_is_ms3 = default_tg_is_ms3
  )

  print("")
}

# Inputs for negative mode
negative_sample_plates <- "~/" # mzML files available on online repository
biological_lib_name <- "do_cr_neg_reduced_20230211.msp"
biological_ms2_search_params <- DIMS_read_search_params(
  file.path(top_level_dir, "habc_v4_biological_neg.json")
)
biological_ms3_search_params <- NULL
is_ms3 <- FALSE
ms3_lib_name <- NULL

for (p in list.files(negative_sample_plates)[nchar(list.files(negative_sample_plates)) == 4]) {
  print("=======================================")
  print(paste0("Processing Negative Mode Plate ", p))
  print("")

  plate_name <- sub("^(.)", "\\L\\1", p, perl = TRUE)
  samples_file_path <- file.path(negative_sample_plates, p)

  docr_process_dims_plate(
    plate_name,
    samples_file_path,
    top_level_dir,
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
    is_ms3,
    save_mzrolldb_as_rds
  )

  print("")
}

### Clean data
### Select on best ion quant
### Remove re-runs
bulkpool_params_file <- file.path(top_level_dir, "habc_v4_bulkpool_20250116.json")
stage_3_rds_dir <- file.path(top_level_dir, "rds_output")

bulkpool_params <- DIMS_read_search_params(bulkpool_params_file)

combined_results <- docr_dims_formatting(
  stage_3_rds_dir,
  bulkpool_params
)

combined_results_file <- file.path(top_level_dir, "combined_results.rds")
# saveRDS(combined_results, file = combined_results_file)
# combined_results <- readRDS(file.path(top_level_dir, "combined_results.rds"))

### Clean names
combined_results <- combined_results %>%
  dplyr::mutate(sample = gsub("X0200_M014B_|X0200_M015B_|_20230422161324|\\.mzML", "", sample)) %>%
  dplyr::rename(name_unique = sample) %>%
  tidyr::separate(name_unique,
    into = c("sample_name", "rm_plate", "well", "inj_order", "run_num"),
    remove = FALSE, sep = "_"
  ) %>%
  tidyr::separate(sample_name,
    into = c("rm_exp", "diet", "mouse_id", "timepoint"),
    remove = FALSE, fill = "right", sep = "-"
  ) %>%
  dplyr::mutate(sample_type = ifelse(grepl("Blank", sample_name), "Blank",
    ifelse(grepl("BulkPool", sample_name), "BulkPool",
      ifelse(grepl("Std", sample_name), "Standard",
        ifelse(grepl("PosControl", sample_name), "PosCtrl",
          ifelse(grepl("NegControl", sample_name), "NegCtrl", "Sample")
        )
      )
    )
  )) %>%
  ### Remove extra samples from adjacent fasting experiment
  dplyr::mutate(fasting_exp = ifelse(grepl("Fast|AdLib|B6BulkPool|B6PosControl", name_unique), TRUE, FALSE)) %>%
  dplyr::filter(!fasting_exp) %>%
  dplyr::select(-tidyr::starts_with("rm_"), -fasting_exp)

### Make finalized dataframe
### Remove TGs that are found in < 75% of biological & bulk samples
### Save ms1 quants for remaining lipids
good_tgs <- combined_results %>%
  dplyr::filter(
    lipidClass == "TG",
    sample_type %in% c("Sample", "BulkPool")
  ) %>%
  dplyr::group_by(compoundName) %>%
  dplyr::mutate(sample_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fraction_found = sample_count / max(sample_count)) %>%
  dplyr::filter(fraction_found > 0.75) %>%
  dplyr::pull(compoundName) %>%
  unique()

results_for_reporting <- combined_results %>%
  dplyr::filter(quant_class == "ms1") %>%
  dplyr::filter(lipidClass != "TG") %>%
  dplyr::bind_rows(combined_results %>%
    dplyr::filter(compoundName %in% good_tgs))

saveRDS(results_for_reporting, file = file.path(save_directory, paste0(gsub("-", "", Sys.Date()), "-PreNorm-Lipidomics.Rds")))

print("Combined file saved")
