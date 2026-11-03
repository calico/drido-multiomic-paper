###
### Mediation Analyses for DO-CR project - Diet -> Bodyweight -> Compound (Year 2)
###
### Counterfactual mediation (natural direct/indirect effects, via the
### `mediation` package, quasi-Bayesian draws) of whether diet's effect on
### each Year 2 compound runs through bodyweight. Same setup as
### mediation-test-diet-final.R, with bodyweight (bw_test) as the mediator
### and the compound as the outcome. Adjusts for baseline (Year 1)
### bodyweight and the baseline compound value. See
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

# Baseline (Year 1, pre-diet) bodyweight per mouse
bw_baseline <- data_use %>%
  dplyr::filter(age_years == "year1") %>%
  dplyr::summarise(bl_bw = dplyr::first(bw_test), .by = mouse_id)

data_use2 <- data_use %>%
  dplyr::anti_join(baseline,
                   by = c("mouse_id", "trait_id", "days_remaining")) %>%
  dplyr::left_join(baseline %>%
                     dplyr::select(mouse_id, trait_id, baseline_value = trait_value),
                   by = c("mouse_id", "trait_id")) %>%
  dplyr::filter(age_years == "year2")

invisible(gc())

data_wide <- data_use2 %>%
  dplyr::select(-diet_fasting, -Age, -diet_assignment,
                -modality, -name_use) %>%
  tidyr::pivot_wider(values_from = c("baseline_value", "trait_value"),
                     names_from = "trait_id") %>%
  dplyr::left_join(bw_baseline, by = "mouse_id")

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

rm(data_use, data_use2, baseline, bw_baseline)
invisible(gc())

# Versions: with and without end-life (PLL > 0.86) cutoff
pll_cutoffs <- list(with_trim = 0.86,
                    without_trim = 1)
covars <- c("generation_wave", "fasting", "bl_bw")

# Diet levels contrasted against ad lib
diet_control <- "AL"
diet_treat_values <- c("1D", "2D", "20", "40")

# Quasi-Bayesian draws (no bootstrap)
n_sims <- 2000
boot_use <- FALSE

med_folder <- "diet_bw_compound_mediation_20261102"

### Start loop ----
cat("------ Starting counterfactual mediation --------\n")

dir.create(file.path(output_filepath, med_folder), showWarnings = FALSE, recursive = TRUE)

future::plan(future::multisession, workers = future::availableCores() - 1)

for (co in names(pll_cutoffs)) {
  data_filtered <- data_wide %>%
    dplyr::filter(PLL <= pll_cutoffs[[co]])

  new_filepath <- file.path(output_filepath, med_folder, paste0(co, "_diet-bw-compound.Rds"))

  if (file.exists(new_filepath)) next

  final_mediation_df <- furrr::future_map_dfr(seq_len(nrow(name_key)), function(j) {
    outcome_col <- paste0("tv_", name_key$safe[j])
    outcome_bl <- paste0("bl_", name_key$safe[j])
    if (!(outcome_col %in% colnames(data_filtered))) return(data.frame())

    docr_counterfactual_mediation(
      data_use = data_filtered,
      outcome_var = outcome_col,
      intervention_var = "diet",
      mediation_var = "bw_test",
      co_vars_a = c(covars, outcome_bl),
      control_value = diet_control,
      treat_values = diet_treat_values,
      boot_use = boot_use,
      n_sims = n_sims
    ) %>%
      dplyr::mutate(feature_id = name_key$trait_id[j])
  }, .options = furrr::furrr_options(seed = TRUE))

  if (nrow(final_mediation_df) > 0) {
    final_mediation_df <- final_mediation_df %>%
      dplyr::mutate(pll_cutoff = co, pll_cutoff_value = pll_cutoffs[[co]])
  }

  saveRDS(final_mediation_df, new_filepath)
  cat(paste0("\n[", co, "] Saved: diet -> bw -> compound"))
}

future::plan(future::sequential)

### Combine, FDR-correct (as in mediation-parse-diet-final.R), and save ----
all_med_test <- purrr::map_dfr(list.files(file.path(output_filepath, med_folder), full.names = TRUE), readRDS) %>%
  dplyr::mutate(modality = dplyr::case_when(
    grepl("\\.M012|\\.M013", feature_id) ~ "metabolomics",
    grepl("\\.M014|\\.M015", feature_id) ~ "lipidomics",
    TRUE ~ "proteomics"
  ))

p_cols <- c("acme_control_p", "acme_treated_p", "acme_avg_p", "acme_diff_p", "total_effect_p")

all_med_test_adj <- purrr::reduce(
  p_cols,
  function(df, p) {
    fdr_multi(df, pval_var = p, nest_vars = c("model_term", "modality", "pll_cutoff"), padj_var = paste0(p, "_adj"))
  },
  .init = all_med_test
)

saveRDS(all_med_test_adj,
        file = file.path(local_filepath,
                         "inst/extdata",
                         paste0(gsub("-", "", Sys.Date()), "-Mediation-CounterfactualTest-Summary-DietBWCompound.Rds")))

print("Data saved")
