# Install and load necessary packages
suppressMessages(library(mgcv))
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"

source(file.path(local_filepath, "R/normalization_functions.R"))
source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

print("Libraries loaded")
save_files <- TRUE # Change to true before running, to save outputs

### Load and clean molecular data ----
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
  dplyr::group_by(trait_id) %>%
  dplyr::mutate(trait_value = as.numeric(scale(trait_value, center = T, scale = F))) %>%
  dplyr::ungroup()

baseline <- data_use %>%
  dplyr::group_by(mouse_id, trait_id) %>%
  dplyr::arrange(desc(days_remaining)) %>%
  dplyr::slice_head(n = 1)

non_baseline <- data_use %>%
  dplyr::anti_join(baseline,
                   by = c("mouse_id", "trait_id", "days_remaining")) %>%
  dplyr::left_join(baseline %>% 
                     dplyr::select(mouse_id, trait_id, baseline_value = trait_value),
                   by = c("mouse_id", "trait_id"))

invisible(gc())

print("Data combined; starting models")

future::plan(future::multisession, workers = 11)

furrr::future_walk(
  .x = unique(non_baseline$trait_id),
  .progress = TRUE,
  function(current_phenotype = .x) {
    tryCatch(
      {
        phenotype_subset <- non_baseline %>%
          dplyr::filter(trait_id == current_phenotype)
    
        has_fasting <- !all(is.na(phenotype_subset$fasting))

        if (has_fasting) {
          model_nob <- mgcv::gam(
            trait_value ~ s(PLL) + diet + fasting + s(mouse_id, bs = "re"),
            data = phenotype_subset, method = "REML"
          )
          model_b <- mgcv::gam(
            trait_value ~ s(PLL) + diet + fasting + baseline_value + s(mouse_id, bs = "re"),
            data = phenotype_subset, method = "REML"
          )
        } else {
          model_nob <- mgcv::gam(
            trait_value ~ s(PLL) + diet + s(mouse_id, bs = "re"),
            data = phenotype_subset, method = "REML"
          )
          model_b <- mgcv::gam(
            trait_value ~ s(PLL) + diet + baseline_value + s(mouse_id, bs = "re"),
            data = phenotype_subset, method = "REML"
          )
        }

        # Save models
        model_pll_return <- list()
        model_pll_return[[current_phenotype]] <- list(
          with_baseline = model_b,
          without_baseline = model_nob
          
        )

        saveRDS(model_pll_return,
          file = paste0(file.path(gam_output_filepath, "gam_models_PLL_6", gsub("\\/|:|\\\\", "-", current_phenotype)), ".Rds")
        )

        return(invisible())
      },
      error = function(e) {},
      warning = function(w) {}
    ) # end tryCatch
  }
) # End future walk

future::plan(sequential)
