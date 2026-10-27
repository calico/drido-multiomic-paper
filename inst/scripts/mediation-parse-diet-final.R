# Install and load necessary packages
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"

# local_filepath <- "~/GitHub/drido-multiomic-paper"
# output_filepath <- "~/GitHub/drido-multiomic-paper/inst/extdata/"

med_analysis_folder <- "diet_mediation_20261027" # output folder from mediation-test-diet-final.R

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

# Normal p-value using the CI half facing zero; avoids the bootstrap p = 0 floor (< 2/n_sims)
ci_pvalue <- function(est, lo, hi) {
  se <- ifelse(est > 0, est - lo, hi - est) / 1.96
  2 * pnorm(-abs(est) / se)
}

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
  dplyr::select(-intervention_var, -mediation_var) %>%
  dplyr::mutate(
    acme_control_p_norm = ci_pvalue(acme_control, acme_control_lo, acme_control_hi),
    acme_treated_p_norm = ci_pvalue(acme_treated, acme_treated_lo, acme_treated_hi),
    acme_avg_p_norm = ci_pvalue(acme_avg, acme_avg_lo, acme_avg_hi),
    total_effect_p_norm = ci_pvalue(total_effect, total_effect_lo, total_effect_hi)
  )

p_cols <- c("acme_control_p", "acme_treated_p", "acme_avg_p", "total_effect_p",
            "acme_control_p_norm", "acme_treated_p_norm", "acme_avg_p_norm", "total_effect_p_norm")

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

saveRDS(all_med_test_adj, file = file.path(local_filepath,
                                           "inst/extdata",
                                           paste0(gsub("-", "", Sys.Date()), "-Mediation-CounterfactualTest-Summary-Diet.Rds")))
