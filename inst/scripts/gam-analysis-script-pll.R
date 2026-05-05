# Install and load necessary packages
suppressMessages(library(mgcv))
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
# local_filepath <- "~/drido-multiomic-paper"
# gam_output_filepath <- "~/"

source(file.path(local_filepath, "R/normalization_functions.R"))
source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

print("Libraries loaded")
save_files <- FALSE # Change to true before running, to save outputs

### Load and clean molecular data ----
metabolomics <- readRDS(file.path(local_filepath, "inst/extdata/20250128-Normalized-Metabolomics-Data.Rds"))
docr_check_factors(metabolomics)

lipidomics <- readRDS(file.path(local_filepath, "inst/extdata/20250128-Normalized-Lipidomics-Data.Rds"))
docr_check_factors(lipidomics)

proteomics <- readRDS(file.path(local_filepath, "inst/extdata/20250129-Normalized-Proteomics-Data.Rds"))
docr_check_factors(proteomics)

### Extract names of features that were run through the LMM ---
### Only run GAMS for features run through LLM - removes low-n features
mol_traits <- readRDS(file.path(local_filepath, "inst/extdata/20250417-Linear-MM-Age-Diet-Results.Rds")) %>%
  dplyr::select(feature_id) %>%
  dplyr::distinct() %>%
  dplyr::pull()

### Combine molecular data ----
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
    age_years != "year1", # Remove pre-CR-intervention data
    trait_id %in% mol_traits
  ) # Only run traits with high enough n

### Load and clean physiological data ----
phenotype_data <- docr_read_phenotype_data(file.path(local_filepath, "/inst/extdata/DOCR_Phenotype_Data.csv"))
phenotype_data <- phenotype_data %>%
  dplyr::rename(
    generation_wave = "Generation", # For error handling
    feature_id = "trait_id"
  ) %>% # For consistency with molecular data
  dplyr::filter(age_in_days > 200) %>% # Remove pre-intervention data
  dplyr::group_by(feature_id) %>%
  dplyr::mutate(norm_abundance = as.numeric(
    scale(trait_value, center = TRUE, scale = TRUE)
  )) %>%
  dplyr::ungroup()

print("Data loaded")

### Combine and run models ----
cols_use <- intersect(colnames(all_molecular_data), colnames(phenotype_data))
data_use <- dplyr::bind_rows(
  all_molecular_data %>%
    dplyr::select(tidyselect::all_of(cols_use)),
  phenotype_data %>%
    dplyr::select(tidyselect::all_of(cols_use))
) %>%
  dplyr::mutate(across(c("mouse_id", "diet", "fasting"), ~ factor(.x))) %>%
  dplyr::mutate(diet_group = ifelse(diet %in% c("AL", "1D", "2D"),
    "AL-1D-2D",
    "20-40"
  )) %>%
  dplyr::mutate(diet_group = factor(diet_group, levels = c("AL-1D-2D", "20-40")))

rm(list = c(
  "lipidomics", "proteomics", "metabolomics",
  "phenotype_data", "all_molecular_data"
))
invisible(gc())

print("Data combined; starting models")

k_use <- -1 # mgcv::gam() default

future::plan(future::multisession, workers = 6)

furrr::future_walk(
  .x = unique(data_use$trait_id),
  .progress = TRUE,
  function(current_phenotype = .x) {
    tryCatch(
      {
        phenotype_subset <- data_use %>%
          dplyr::filter(trait_id == current_phenotype)

        if (all(is.na(phenotype_subset$fasting))) {
          model1 <- mgcv::gam(
            trait_value ~
              s(PLL,
                by = diet,
                k = k_use
              ) +
              diet +
              s(mouse_id, bs = "re"),
            data = phenotype_subset,
            method = "REML"
          )

          model2 <- mgcv::gam(
            trait_value ~
              s(PLL,
                by = diet_group,
                k = k_use
              ) +
              diet +
              s(mouse_id, bs = "re"),
            data = phenotype_subset,
            method = "REML"
          )

          model3 <- mgcv::gam(
            trait_value ~
              s(PLL,
                k = k_use
              ) +
              diet +
              s(mouse_id, bs = "re"),
            data = phenotype_subset,
            method = "REML"
          )
        } else {
          model1 <- mgcv::gam(
            trait_value ~
              s(PLL,
                by = diet,
                k = k_use
              ) +
              fasting +
              diet +
              s(mouse_id, bs = "re"),
            data = phenotype_subset,
            method = "REML"
          )

          model2 <- mgcv::gam(
            trait_value ~
              s(PLL,
                by = diet_group,
                k = k_use
              ) +
              fasting +
              diet +
              s(mouse_id, bs = "re"),
            data = phenotype_subset,
            method = "REML"
          )

          model3 <- mgcv::gam(
            trait_value ~
              s(PLL,
                k = k_use
              ) +
              fasting +
              diet +
              s(mouse_id, bs = "re"),
            data = phenotype_subset,
            method = "REML"
          )
        }

        # Save models
        model_pll_return <- list()
        model_pll_return[[current_phenotype]] <- list(
          by_single_diet = model1,
          by_diet_group = model2,
          all = model3
        )

        saveRDS(model_pll_return,
          file = paste0(file.path(gam_output_filepath, "gam_models_PLL_5", gsub("\\/|:|\\\\", "-", current_phenotype)), ".Rds")
        )

        return(invisible())
      },
      error = function(e) {},
      warning = function(w) {}
    ) # end tryCatch
  }
) # End future walk

future::plan(sequential)
