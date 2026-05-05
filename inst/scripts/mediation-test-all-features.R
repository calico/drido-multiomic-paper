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
all_molecular_data <- metabolomics %>%
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
  )

rm(list = c("lipidomics", "metabolomics", "proteomics"))

cat("Data loaded")

cat("------ Starting mediation --------")

all_features <- unique(all_molecular_data$trait_id)

for (afh in all_features) {
  ### Combine and run models ----
  all_da <- all_molecular_data %>%
    dplyr::filter(trait_id != afh) %>% # remove mediator as exposure
    dplyr::left_join(
      all_molecular_data %>%
        dplyr::filter(trait_id == afh) %>%
        dplyr::select(sample_id, PROTEIN_MEDIATOR = norm_abundance),
      by = "sample_id"
    ) %>%
    dplyr::rename(!!rlang::sym(afh) := "PROTEIN_MEDIATOR") %>%
    tidyr::nest(data = -c(feature_id)) %>%
    purrr::transpose()

  invisible(gc())

  cat(paste0("Data combined for mediator: ", afh))

  ### Mediation Analyses -----

  print(paste0("Starting modeling: ", afh))

  afh_name <- gsub("/", "-", afh)
  new_filepath <- file.path(
    output_filepath, "feature_mediation", paste0(afh_name, "-Mediator.Rds")
  )

  if (!file.exists(new_filepath)) {
    # Parallelize computing
    future::plan(future::multisession, workers = future::availableCores() - 1)

    final_mediation_df <- furrr::future_map_dfr(all_da, function(i) {
      current_data <- i$data
      feat <- gsub("/", "-", i$feature_id)
      afh_name <- gsub("/", "-", afh)

      all_da_temp <- docr_lm_sobel_mediation_test(
        data_use = current_data,
        outcome_var = "surv_years",
        intervention_var = afh,
        mediation_var = "norm_abundance",
        co_vars = c("generation_wave", "fasting", "diet", "bw_test"),
        n_samples = 1000 # will not run if less than 100 sample points
      ) %>%
        dplyr::mutate(
          mediation_var = i$feature_id
        )
      return(all_da_temp)
    })
    future::plan(future::sequential)
    saveRDS(final_mediation_df, new_filepath)
  }
}

print("Data saved")
