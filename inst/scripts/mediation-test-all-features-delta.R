###
### Mediation Analyses for DO-CR project - Molecular-Molecular Hub Interactions
### Delta (change from baseline) approach
### Johanna Fleischman, Calico Life Sciences 2026
###
###
### Install libraries and source functions
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# ### Assign file paths
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"

# local_filepath <- "~/GitHub/drido-multiomic-paper"
# output_filepath <- "~/Desktop"

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

# get baseline data and compute delta (change from baseline)
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
  dplyr::filter(age_years == "year2") %>%
  dplyr::mutate(delta_value = trait_value - baseline_value)

invisible(gc())

# Pivot wide using delta values instead of separate tv/bl columns
data_wide <- data_use2 %>%
  dplyr::select(-diet_fasting, -Age, -age_years, -diet_assignment,
                -modality, -name_use, -trait_value, -baseline_value) %>%
  tidyr::pivot_wider(values_from = "delta_value",
                     names_from = "trait_id")

print("Data combined")

# Create safe name mapping for formula-unfriendly column names
all_features <- unique(data_use2$trait_id)
name_key <- data.frame(
  trait_id = all_features,
  safe = paste0("V", seq_along(all_features))
)

trait_col_map <- setNames(
  paste0("d_", name_key$safe),
  name_key$trait_id
)
colnames(data_wide) <- ifelse(
  colnames(data_wide) %in% names(trait_col_map),
  trait_col_map[colnames(data_wide)],
  colnames(data_wide)
)

rm(data_use, data_use2, baseline)
invisible(gc())

# Versions: with and without bodyweight, with and without end-life cut off
pll_cutoffs <- list(with_trim = 0.86,
                    without_trim = 1)
cv_use <- c("generation_wave", "fasting", "diet")

### Start loop ----
cat("------ Starting mediation (delta approach) --------\n")

future::plan(future::multisession, workers = future::availableCores() - 1)

for (co in names(pll_cutoffs)) {

  data_filtered <- data_wide %>%
    dplyr::filter(PLL <= pll_cutoffs[[co]])

  for (i in seq_len(nrow(name_key))) {

    afh <- name_key$trait_id[i]
    afh_safe <- name_key$safe[i]
    afh_name <- gsub("/", "-", afh)

    new_filepath <- file.path(
      output_filepath, "feature_mediation_delta",
      paste0(co, "_", afh_name, "-Intervention.Rds")
    )

    if (file.exists(new_filepath)) next

    intervention_var <- paste0("d_", afh_safe)
    if (!(intervention_var %in% colnames(data_filtered))) next

    other_idx <- setdiff(seq_len(nrow(name_key)), i)

    final_mediation_df <- furrr::future_map_dfr(other_idx, function(j) {
      mediator_var <- paste0("d_", name_key$safe[j])
      if (!(mediator_var %in% colnames(data_filtered))) return(data.frame())

      docr_sobel_mediation(
        data_use = data_filtered,
        outcome_var = "surv_days",
        intervention_var = intervention_var,
        mediation_var = mediator_var,
        co_vars_a = cv_use,
        co_vars_total = cv_use,
        co_vars_direct = cv_use
      ) %>%
        dplyr::mutate(mediation_var = name_key$trait_id[j])
    })

    if (nrow(final_mediation_df) > 0) {
      final_mediation_df <- final_mediation_df %>%
        dplyr::mutate(
          intervention_var = afh,
          pll_cutoff = co,
          pll_cutoff_value = pll_cutoffs[[co]]
        )
    }

    saveRDS(final_mediation_df, new_filepath)
    cat(paste0("\n[", co, "] Saved: ", afh_name))
  }
}

future::plan(future::sequential)

print("Data saved")
