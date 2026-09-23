# Install and load necessary packages
suppressMessages(library(mgcv))
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
local_filepath <- "~/workspace/drido-multiomic-paper"
gam_output_filepath <- "~/workspace/docr_data/MS1553/gam_models_PLL_6"

source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

name_conv_file <- "inst/supp_tables/Table_S1_CompoundAnnotations.csv"

### Process Models ### -------

# Make new data
# mouse_id added in docr_gam_predict() fxn
new_data <- data.frame(
  PLL = seq(0.3, 1, by = 0.01), # predict every 1% PLL
  fasting = "No",
  baseline_value = 0
)

# Initialize data
all_summary_data <- tibble()
all_prediction_data <- tibble()
all_residual_data <- list()

print("Reading in files")
# Collect results from models
future::plan(future::multisession, workers = future::availableCores() - 3)
results <- furrr::future_map(
  .x = list.files(file.path(gam_output_filepath), full.names = TRUE),
  ~ docr_gam_process(
    mt = .x,
    new_data = new_data
  ),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)
future::plan(future::sequential)
print("Files read")

# Combine all results
all_summary_data <- purrr::map_dfr(results, ~ .x$summary)
all_prediction_data <- purrr::map_dfr(results, ~ .x$prediction)

# Process residual data
for (res in results) {
  if (is.null(res$residual)) next        
  trait <- names(res$residual)[1]
  all_residual_data[[trait]] <- res$residual[[trait]]
}

# Save the results
saveRDS(
  list(all_summary_data, all_prediction_data, all_residual_data),
  file = file.path(gam_output_filepath, paste0(gsub("-", "", Sys.Date()), "-GAM-Models-Parsed-Baseline.Rds"))
)
print("Summary RDS saved")