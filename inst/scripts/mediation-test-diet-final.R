###
### Mediation Analyses for DO-CR project - Diet Interventions (Final)
###
### Counterfactual mediation (natural direct/indirect effects, via the
### `mediation` package) in place of the Sobel/product-of-coefficients
### approach used in mediation-test-diet-baseline.R and
### mediation-test-diet-agingPLL.R. The outcome model now includes a
### diet x mediator interaction term -- per Reviewer 2 comment 4.1, this is
### the one mediation design in the paper (diet -> feature -> lifespan) that
### is temporally ordered, and a counterfactual formulation allowing
### exposure-mediator interaction strengthens it. See
### docr_counterfactual_mediation() in R/statistics_functions.R for the
### model-fitting details.
###
### Johanna Fleischman, Calico Life Sciences 2026
###
###
### Install libraries and source functions
suppressMessages(library(tidyverse))
suppressMessages(library(mediation))
library(future)
library(furrr)

# ### Assign file paths
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"

# local_filepath <- "~/GitHub/drido-multiomic-paper"
# output_filepath <- "~/GitHub/drido-multiomic-paper/inst/extdata/"

source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))
source(file.path(local_filepath, "R/normalization_functions.R"))

lipidomics_data <- "inst/extdata/20250128-Normalized-Lipidomics-Data.Rds"
metabolomics_data <- "inst/extdata/20250128-Normalized-Metabolomics-Data.Rds"
proteomics_data <- "inst/extdata/20250129-Normalized-Proteomics-Data.Rds"
phenotype_data <- "inst/extdata/DOCR_Phenotype_Data.csv"
name_conv_file <- "inst/supp_tables/Table_S1_CompoundAnnotations.csv"
name_conversion_use <- read.csv(file.path(local_filepath, name_conv_file)) %>%
  dplyr::select(feature_id, name_use)

# import data
data_use <- docr_make_final_data(
  metabolomics_data_filepath = file.path(local_filepath, metabolomics_data),
  proteomics_data_filepath = file.path(local_filepath, proteomics_data),
  lipidomics_data_filepath = file.path(local_filepath, lipidomics_data),
  phenotype_data_filepath = file.path(local_filepath, phenotype_data),
  name_conversion_key = name_conversion_use
)

print("Functions and files loaded")

# get baseline data
data_use <- data_use %>%
  dplyr::filter(modality != "physiological") %>%
  dplyr::group_by(trait_id) %>%
  dplyr::mutate(trait_value = as.numeric(scale(trait_value, center = T, scale = F))) %>%
  dplyr::ungroup()

baseline <- data_use %>%
  dplyr::group_by(mouse_id, trait_id) %>%
  dplyr::arrange(desc(days_remaining)) %>%
  dplyr::slice_head(n = 1)

data_use2 <- data_use %>%
  dplyr::anti_join(baseline,
                   by = c("mouse_id", "trait_id", "days_remaining")) %>%
  dplyr::left_join(baseline %>%
                     dplyr::select(mouse_id, trait_id, baseline_value = trait_value),
                   by = c("mouse_id", "trait_id")) %>%
  dplyr::filter(age_years == "year2")

invisible(gc())

data_wide <- data_use2 %>%
  dplyr::select(-diet_fasting, -Age, -age_years, -diet_assignment,
                -modality, -name_use) %>%
  tidyr::pivot_wider(values_from = c("baseline_value", "trait_value"),
                     names_from = "trait_id")

print("Data combined")

# Create safe name mapping for formula-unfriendly column names
all_features <- unique(data_use2$trait_id)
name_key <- data.frame(
  trait_id = all_features,
  safe = paste0("V", seq_along(all_features))
)

trait_col_map <- setNames(
  c(paste0("tv_", name_key$safe), paste0("bl_", name_key$safe)),
  c(paste0("trait_value_", name_key$trait_id), paste0("baseline_value_", name_key$trait_id))
)
colnames(data_wide) <- ifelse(
  colnames(data_wide) %in% names(trait_col_map),
  trait_col_map[colnames(data_wide)],
  colnames(data_wide)
)

rm(data_use, data_use2, baseline)
invisible(gc())

# Versions: with and without bodyweight, with and without end-life (PLL > 0.86) cutoff
pll_cutoffs <- list(with_trim = 0.86,
                    without_trim = 1)
covars <- list(with_bw = c("generation_wave", "fasting", "diet", "bw_test"),
               without_bw = c("generation_wave", "fasting", "diet"))

# Diet levels contrasted against ad lib
diet_control <- "AL"
diet_treat_values <- c("1D", "2D", "20", "40")

# Number of quasi-Bayesian draws mediate() uses per ACME/ADE CI
# Tested at higher simulation for the traits that mattered
n_sims <- 500

### Start loop ----
cat("------ Starting counterfactual mediation --------\n")

for (co in names(pll_cutoffs)) {
  for (cv_name in names(covars)) {
    dir.create(
      file.path(output_filepath, "diet_mediation_20261006"),
      showWarnings = FALSE, recursive = TRUE
    )
  }
}

future::plan(future::multisession, workers = future::availableCores() - 1)

for (co in names(pll_cutoffs)) {

  data_filtered <- data_wide %>%
    dplyr::filter(PLL <= pll_cutoffs[[co]])

  for (cv_name in names(covars)) {
    cv_use <- covars[[cv_name]]

    new_filepath <- file.path(
      output_filepath, "diet_mediation_20261006",
      paste0(co, "_", cv_name, "_diet-Intervention.Rds")
    )

    if (file.exists(new_filepath)) next

    final_mediation_df <- furrr::future_map_dfr(seq_len(nrow(name_key)), function(j) {
      mediator_var <- paste0("tv_", name_key$safe[j])
      mediator_bl <- paste0("bl_", name_key$safe[j])
      if (!(mediator_var %in% colnames(data_filtered))) return(data.frame())

      docr_counterfactual_mediation(
        data_use = data_filtered,
        outcome_var = "surv_days",
        intervention_var = "diet",
        mediation_var = mediator_var,
        co_vars_a = c(cv_use, mediator_bl),
        co_vars_direct = c(cv_use, mediator_bl),
        control_value = diet_control,
        treat_values = diet_treat_values,
        n_sims = n_sims
      ) %>%
        dplyr::mutate(mediation_var = name_key$trait_id[j])
    })

    if (nrow(final_mediation_df) > 0) {
      final_mediation_df <- final_mediation_df %>%
        dplyr::mutate(
          intervention_var = "diet",
          pll_cutoff = co,
          pll_cutoff_value = pll_cutoffs[[co]],
          covar_set = cv_name
        )
    }

    saveRDS(final_mediation_df, new_filepath)
    cat(paste0("\n[", co, "/", cv_name, "] Saved: diet intervention"))
  }
}

future::plan(future::sequential)

print("Data saved")
