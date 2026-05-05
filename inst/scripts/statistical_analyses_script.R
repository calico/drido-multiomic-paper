###
### Statistical modeling for DO-CR project
### Johanna Fleischman, Calico Life Sciences 2024
###
###
### Install libraries and source functions
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
# local_filepath <- "~/drido-multiomic-paper"

source(file.path(local_filepath, "R/normalization_functions.R"))
source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

print("Libraries loaded")

### Load files
metabolomics <- readRDS(file.path(local_filepath, "20250128-Normalized-Metabolomics-Data.Rds")) %>%
  dplyr::mutate(modality = "metabolomics")
docr_check_factors(metabolomics)

lipidomics <- readRDS(file.path(local_filepath, "20250128-Normalized-Lipidomics-Data.Rds")) %>%
  dplyr::mutate(modality = "lipidomics")
docr_check_factors(lipidomics)

proteomics <- readRDS(file.path(local_filepath, "20250129-Normalized-Proteomics-Data.Rds")) %>%
  dplyr::mutate(modality = "proteomics")
docr_check_factors(proteomics)

phenotype_data <- docr_read_phenotype_data(file.path(local_filepath, "DOCR_Phenotype_Data.csv"))

print("Data loaded and clean")

all <- dplyr::bind_rows(
  metabolomics,
  lipidomics,
  proteomics
) %>%
  dplyr::mutate(
    days_remaining = as.numeric(delta_date),
    fasting = ifelse(diet %in% c("20") & weekday_collection == "Monday",
      "Fast20",
      ifelse(diet %in% c("40") & weekday_collection == "Monday",
        "Fast40",
        "No"
      )
    ),
    fasting = factor(fasting, levels = c("No", "Fast20", "Fast40"))
  )

rm(metabolomics, lipidomics, proteomics)

### Linear Mixed Modeling -----
#### Age and Diet ----

# Parallelize computing
future::plan(future::multisession, workers = future::availableCores() - 1)

n_samples <- length(unique(all$sample_id))
lm_formulas <- c(
  "norm_abundance ~ age_years + diet + BW_Loess + fasting + (1|mouse_id)" # ,
  #  "norm_abundance ~ age_years + diet + BW_Loess + (1|mouse_id)",
  #  "norm_abundance ~ age_years + diet + fasting + (1|mouse_id)",
  #  "norm_abundance ~ age_years + diet + (1|mouse_id)"
)

mzroll_all_da <- all %>%
  dplyr::select(-c(
    "mz", "rt", "date_birth", "date_collection",
    "diet_assignment", "BW_Raw", "surv_days", "surv_years",
    "date_exit", "delta_date", "is_ddm", "weekday_collection",
    "is_imputed", "PA.Gene.Symbol", "days_remaining"
  )) %>%
  tidyr::nest(data = -c(feature_id, compoundName, method, modality)) %>%
  dplyr::mutate(da_data = furrr::future_map2(
    .x = data,
    .y = feature_id,
    .f = ~ lmer_multi_formula(
      data_use = .x,
      feature_id_test = .y,
      forms_to_test = lm_formulas,
      n_samples = n_samples
    ),
    .options = furrr::furrr_options(
      globals = c("lmer_multi_formula", "lm_formulas", "n_samples")
    )
  )) %>%
  dplyr::select(-data, -feature_id) %>%
  tidyr::unnest(da_data) %>%
  fdr_multi(
    pval_var = "p_value_coef",
    nest_vars = c("model", "model_term", "modality")
  )
future::plan(future::sequential)

saveRDS(object = mzroll_all_da, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Linear-MM-Age-Diet-Results.Rds")))

print("Linear mixed modeling results saved")


#### PLL Molecular ----

future::plan("multisession", workers = future::availableCores() - 1)

# PLL pre-85%, post year 1; 85% is the median greatest inflection point on
# non-linear PLL trajectories via GAM modeling
lm_formulas <- c(
  "norm_abundance ~ PLL + diet + BW_Loess + fasting + (1|mouse_id)"
)
pll_cutoff <- 0.85

# Reduction for furrr
all_reduced <- all %>%
  dplyr::select(-c(
    "mz", "rt", "date_birth", "date_collection",
    "diet_assignment", "BW_Raw", "surv_days", "surv_years",
    "date_exit", "delta_date", "is_ddm", "weekday_collection",
    "is_imputed", "PA.Gene.Symbol", "days_remaining"
  ))

rm(all)
gc()

mzroll_all_da <- all_reduced %>%
  dplyr::filter(
    age_years != "year1",
    PLL <= pll_cutoff
  ) %>%
  # For scaled version only
  dplyr::group_by(feature_id) %>%
  dplyr::mutate(norm_abundance = as.numeric(
    scale(norm_abundance, center = TRUE, scale = TRUE)
  )) %>%
  dplyr::ungroup() %>%
  tidyr::nest(data = -c(feature_id, compoundName, method, modality)) %>%
  dplyr::mutate(da_data = furrr::future_map2(
    .x = data,
    .y = feature_id,
    .f = ~ lmer_multi_formula(
      data_use = .x,
      feature_id_test = .y,
      forms_to_test = lm_formulas,
      n_samples = 150
    )
  )) %>%
  dplyr::select(-data, -feature_id) %>%
  tidyr::unnest(da_data) %>%
  fdr_multi(
    pval_var = "p_value_coef",
    nest_vars = c("model", "model_term", "modality")
  )

saveRDS(object = mzroll_all_da, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Linear-Scaled-MM-PLL-Aging-Results.Rds")))


# PLL post-85%
mzroll_all_da <- all_reduced %>%
  dplyr::filter(
    age_years != "year1",
    PLL > pll_cutoff
  ) %>%
  # For scaled version only
  dplyr::group_by(feature_id) %>%
  dplyr::mutate(norm_abundance = as.numeric(
    scale(norm_abundance, center = TRUE, scale = TRUE)
  )) %>%
  dplyr::ungroup() %>%
  # Resume stats
  tidyr::nest(data = -c(feature_id, compoundName, method, modality)) %>%
  dplyr::mutate(da_data = furrr::future_map2(
    .x = data,
    .y = feature_id,
    .f = ~ lmer_multi_formula(
      data_use = .x,
      feature_id_test = .y,
      forms_to_test = lm_formulas,
      n_samples = 150
    )
  )) %>%
  dplyr::select(-data, -feature_id) %>%
  tidyr::unnest(da_data) %>%
  fdr_multi(
    pval_var = "p_value_coef",
    nest_vars = c("model", "model_term", "modality")
  )

future::plan(future::sequential)

saveRDS(object = mzroll_all_da, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Linear-Scaled-MM-PLL-Dying-Results.Rds")))

#### PLL Phenotypic ----

future::plan("multisession", workers = future::availableCores() - 1)

# PLL pre-85%, post year 1; 85% is the median greatest inflection point on
# non-linear PLL trajectories via GAM modeling
pll_cutoff <- 0.85

# So I can run the same model on all physiological traits without erroring
phenotype_data$bw_test[is.na(phenotype_data$bw_test) & phenotype_data$trait_id %in% c("BW_BW", "BW_Delta", "BW_PhenoDelta")] <- 0

mzroll_all_da <- phenotype_data %>%
  dplyr::rename(
    generation_wave = "Generation", # For error handling
    #  norm_abundance = "trait_value",
    feature_id = "trait_id"
  ) %>% # For consistency with molecular data
  dplyr::filter(
    age_in_days > 200, # remove pre-CR intervention data
    PLL <= pll_cutoff
  ) %>%
  # Z-score physiological traits so that they are be compared to molecular traits
  dplyr::group_by(feature_id) %>%
  dplyr::mutate(norm_abundance = as.numeric(
    scale(trait_value, center = TRUE, scale = TRUE)
  )) %>%
  dplyr::ungroup() %>%
  # Resume stats
  tidyr::nest(data = -c(feature_id, modality)) %>%
  dplyr::mutate(da_data = furrr::future_map2(
    .x = data,
    .y = feature_id,
    .f = ~ lmer_multi_formula(
      data_use = .x,
      feature_id_test = .y,
      forms_to_test = "norm_abundance ~ PLL + diet + bw_test + (1|mouse_id)",
      n_samples = 150
    )
  )) %>%
  dplyr::select(-data, -feature_id) %>%
  tidyr::unnest(da_data) %>%
  fdr_multi(
    pval_var = "p_value_coef",
    nest_vars = c("model", "model_term")
  )


saveRDS(object = mzroll_all_da, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Linear-MM-PLL-Aging-Phenotype-Results.Rds")))


# PLL post-85%
mzroll_all_da <- phenotype_data %>%
  dplyr::rename(
    generation_wave = "Generation", # For error handling
    #  norm_abundance = "trait_value",
    feature_id = "trait_id"
  ) %>% # For consistency with molecular data
  dplyr::filter(
    age_in_days > 200, # remove pre-CR intervention data
    PLL > pll_cutoff
  ) %>%
  # Z-score physiological traits so that they are be compared to molecular traits
  dplyr::group_by(feature_id) %>%
  dplyr::mutate(norm_abundance = as.numeric(
    scale(trait_value, center = TRUE, scale = TRUE)
  )) %>%
  dplyr::ungroup() %>%
  # Resume stats
  tidyr::nest(data = -c(feature_id, modality)) %>%
  dplyr::mutate(da_data = furrr::future_map2(
    .x = data,
    .y = feature_id,
    .f = ~ lmer_multi_formula(
      data_use = .x,
      feature_id_test = .y,
      forms_to_test = "norm_abundance ~ PLL + diet + bw_test + (1|mouse_id)",
      n_samples = 150
    )
  )) %>%
  dplyr::select(-data, -feature_id) %>%
  tidyr::unnest(da_data) %>%
  fdr_multi(
    pval_var = "p_value_coef",
    nest_vars = c("model", "model_term", "modality")
  )

future::plan(future::sequential)

saveRDS(object = mzroll_all_da, file = file.path(local_filepath, paste0(gsub("-", "", Sys.Date()), "-Linear-MM-PLL-Dying-Phenotype-Results.Rds")))

print("Linear mixed modeling results saved for both PLL models")

future::plan(future::sequential)
