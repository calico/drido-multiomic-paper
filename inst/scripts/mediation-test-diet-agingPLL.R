###
### Mediation Analyses for DO-CR project
### Johanna Fleischman, Calico Life Sciences 2024
###
###
### Install libraries and source functions
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# ### Assign file paths
# local_filepath <- "~/drido-multiomic-paper"
# output_filepath <- "~/"

source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))
source(file.path(local_filepath, "R/normalization_functions.R"))

# import data
lipidomics <- readRDS(file.path(local_filepath, "inst/extdata/20250128-Normalized-Lipidomics-Data.Rds"))
metabolomics <- readRDS(file.path(local_filepath, "inst/extdata/20250128-Normalized-Metabolomics-Data.Rds"))
proteomics <- readRDS(file.path(local_filepath, "inst/extdata/20250129-Normalized-Proteomics-Data.Rds"))
good_ids <- read.csv(file.path(local_filepath, "inst/extdata/Table_S1_CompoundAnnotations.csv"))

### Combine all data ----
all_da <- metabolomics %>%
  dplyr::mutate(modality = "metabolomics") %>%
  dplyr::bind_rows(lipidomics %>%
    dplyr::mutate(modality = "lipidomics")) %>%
  dplyr::bind_rows(proteomics %>%
    dplyr::mutate(modality = "proteomics")) %>%
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
  ) %>%
  dplyr::mutate(
    trait_id = feature_id,
    bw_test = BW_Raw
  ) %>%
  dplyr::group_by(trait_id) %>%
  dplyr::mutate(trait_value = scale(norm_abundance, center = T, scale = F)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    age_years == "year2", # Remove pre-CR-intervention data
    feature_id %in% good_ids$feature_id,
    PLL <= 0.85
  ) %>%
  tidyr::nest(data = -c(feature_id, age_years)) %>%
  purrr::transpose()

rm(list = c("lipidomics", "metabolomics", "proteomics"))

cat("Data loaded")

cat("------ Starting mediation --------")

### Mediation Analyses -----

# Parallelize computing
future::plan(future::multisession, workers = 10)

n_samples <- 1000

print("Starting modeling")

furrr::future_walk(all_da, function(i) {
  current_data <- i$data
  feat <- gsub("/", "-", i$feature_id)
  age <- i$age_years

  all_da_temp <- docr_lm_sobel_mediation_test(
    data_use = current_data,
    outcome_var = "surv_years",
    intervention_var = "diet",
    mediation_var = "norm_abundance",
    co_vars = c("generation_wave", "fasting", "bw_test"),
    n_samples = n_samples
  ) %>%
    dplyr::mutate(
      feature_id = i$feature_id,
      age_years = age
    )

  saveRDS(all_da_temp, file = file.path(output_filepath, "med-analysis-diet-lifespan-20260310", paste0(feat, ".", age, ".Rds")))
})

future::plan(future::sequential)

print("Data saved")
