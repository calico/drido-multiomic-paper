# Install and load necessary packages
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"
med_analysis_folder <- "diet_mediation_20261006" # output folder from mediation-test-diet-final.R

source(file.path(local_filepath, "R/statistics_functions.R"))

future::plan(future::multisession, workers = future::availableCores() - 1)

## Determine mediators
all_med_test <- furrr::future_map_dfr(
  .x = list.files(file.path(output_filepath, med_analysis_folder)),
  .f = ~ {
    temp <- readRDS(file.path(output_filepath, med_analysis_folder, .x))
    if (nrow(temp) == 0) {
      return(NULL)
    } else {
      return(temp)
    }
  }
)

future::plan(future::sequential)

all_med_test_adj <- all_med_test %>%
  dplyr::mutate(
    feature_id = mediation_var,
    modality =
      dplyr::case_when(
        grepl("\\.M012|\\.M013", feature_id) ~ "metabolomics",
        grepl("\\.M014|\\.M015", feature_id) ~ "lipidomics",
        TRUE ~ "proteomics"
      )
  ) %>%
  dplyr::select(-intervention_var, -mediation_var) %>%
  fdr_multi(
    pval_var = "acme_control_p",
    nest_vars = c("model_term", "outcome_var", "modality", "pll_cutoff", "covar_set"),
    padj_var = "acme_control_p_adj"
  ) %>%
  fdr_multi(
    pval_var = "acme_treated_p",
    nest_vars = c("model_term", "outcome_var", "modality", "pll_cutoff", "covar_set"),
    padj_var = "acme_treated_p_adj"
  ) %>%
  fdr_multi(
    pval_var = "acme_avg_p",
    nest_vars = c("model_term", "outcome_var", "modality", "pll_cutoff", "covar_set"),
    padj_var = "acme_avg_p_adj"
  ) %>%
  fdr_multi(
    pval_var = "total_effect_p",
    nest_vars = c("model_term", "outcome_var", "modality", "pll_cutoff", "covar_set"),
    padj_var = "total_effect_p_adj"
  )

saveRDS(all_med_test_adj, file = file.path(local_filepath,
                                           "inst/extdata",
                                           paste0(gsub("-", "", Sys.Date()), "-Mediation-CounterfactualTest-Summary-Diet.Rds")))
