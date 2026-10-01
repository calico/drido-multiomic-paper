# Install and load necessary packages
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
med_analysis_folder <- "feature_mediation_baseline_dr" # output folder from mediation-test-all-features.R
output_filepath <- "~/workspace/docr_data/MS1553"

source(file.path(local_filepath, "R/statistics_functions.R"))

future::plan(future::multisession, workers = future::availableCores() - 1)

## Determine mediators
all_med_test <- furrr::future_map_dfr(
  .x = list.files(file.path(output_filepath, med_analysis_folder)),
  .f = ~ {
    temp <- readRDS(file.path(output_filepath, med_analysis_folder, .x))
    if (nrow(temp) == 0) {
      return(NULL)
    }
    return(temp)
  }
)

future::plan(future::sequential)

all_med_test_adj <- all_med_test %>%
  dplyr::mutate(
    modality =
      dplyr::case_when(
        grepl("\\.M012|\\.M013", mediation_var) ~ "metabolomics",
        grepl("\\.M014|\\.M015", mediation_var) ~ "lipidomics",
        TRUE ~ "proteomics"
      )
  ) %>%
  dplyr::filter(!is.na(sobel_p)) %>%
  fdr_multi(
    pval_var = "sobel_p",
    nest_vars = c("model_term", "outcome_var", "modality", "pll_cutoff", "covar_set"),
    padj_var = "sobel_p_adj"
  )

saveRDS(all_med_test_adj, file = file.path(local_filepath,
                                           "inst/extdata",
                                           paste0(gsub("-", "", Sys.Date()), "-Mediation-Sobel-Test-Summary-AllFeatures-DRCutoff.Rds")))



