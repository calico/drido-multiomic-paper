# Install and load necessary packages
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
med_analysis_folder <- "feature_mediation_baseline" # output folder from mediation-test-all-features.R
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

    # Extract metadata from filename: {cutoff}_{covar_set}_{trait}-Intervention.Rds
    parts <- stringr::str_match(.x, "^(with_trim|without_trim)_(with_bw|without_bw)_(.+)-Intervention\\.Rds$")
    if (!is.na(parts[1, 1])) {
      temp$pll_cutoff <- parts[1, 2]
      temp$covar_set <- parts[1, 3]
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

saveRDS(all_med_test_adj, file = file.path(output_filepath, paste0(gsub("-", "", Sys.Date()), "-Mediation-Sobel-Test-ALL-Mol-Ints.Rds")))
