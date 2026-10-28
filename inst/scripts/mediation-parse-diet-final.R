# Install and load necessary packages
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"

# local_filepath <- "~/GitHub/drido-multiomic-paper"
# output_filepath <- "~/GitHub/drido-multiomic-paper/inst/extdata/"

med_analysis_folder <- "diet_mediation_20261028" # output folder from mediation-test-diet-final.R

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

all_med_test_clean <- all_med_test %>%
  dplyr::mutate(
    feature_id = mediation_var,
    modality =
      dplyr::case_when(
        grepl("\\.M012|\\.M013", feature_id) ~ "metabolomics",
        grepl("\\.M014|\\.M015", feature_id) ~ "lipidomics",
        TRUE ~ "proteomics"
      )
  ) %>%
  dplyr::select(-mediation_var)

p_cols <- c("acme_control_p", "acme_treated_p", "acme_avg_p", "acme_diff_p", "total_effect_p")

all_med_test_adj <- purrr::reduce(
  p_cols,
  function(df, p) {
    fdr_multi(
      df,
      pval_var = p,
      nest_vars = c("model_term", "outcome_var", "modality", "pll_cutoff", "covar_set"),
      padj_var = paste0(p, "_adj")
    )
  },
  .init = all_med_test_clean
)

# Diet and body weight interventions saved separately; FDR is already
# within intervention via model_term ("diet<level>" vs "bw_test")
saveRDS(all_med_test_adj %>%
          dplyr::filter(intervention_var == "diet") %>%
          dplyr::select(-intervention_var),
        file = file.path(local_filepath,
                         "inst/extdata",
                         paste0(gsub("-", "", Sys.Date()), "-Mediation-CounterfactualTest-Summary-Diet.Rds")))

saveRDS(all_med_test_adj %>%
          dplyr::filter(intervention_var == "bw_test") %>%
          dplyr::select(-intervention_var),
        file = file.path(local_filepath,
                         "inst/extdata",
                         paste0(gsub("-", "", Sys.Date()), "-Mediation-CounterfactualTest-Summary-BW.Rds")))
