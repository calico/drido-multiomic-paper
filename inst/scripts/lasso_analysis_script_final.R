library(glmnet)
library(survival)
library(dplyr)
library(purrr)
library(furrr)
library(limma)
library(ggplot2)
library(future)

# args = commandArgs(trailingOnly = TRUE)
#
# if (length(args) < 4) {
#   stop("Need at least 4 arguments: local_filepath, output_filepath, n_bootstrap", call. = FALSE)
# }
#
# local_filepath <- args[1]
# output_filepath <- args[2]
# n_bootstrap <- as.numeric(args[3])


## arguments
model_alpha <- 1
n_bootstrap <- 1000

# ### Assign file paths
local_filepath <- "~/Github/drido-multiomic-paper"
output_filepath <- "~/Github/drido-multiomic-paper/inst/extdata"

source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

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

# filter data for regularization
# molecular traits only; year2 only; pre-terminal decline only
dat <- data_use %>%
  as.data.frame() %>%
  dplyr::filter(
    !is.na(age_years),
    age_years == "year2",
    modality != "physiological"
  ) %>%
  dplyr::mutate(age_days = surv_days - days_remaining)

# get compound data
X_metabo_df_raw <- dat %>%
  dplyr::select(mouse_id, trait_id, trait_value) %>%
  tidyr::pivot_wider(names_from = trait_id, values_from = trait_value)

# get covariate and lifespan data
X_control_df_raw <- dat %>%
  dplyr::distinct(mouse_id, fasting, bw_test, diet, age_days, surv_days)

# make final data frame; mouse_ids in same order
final_data <- dplyr::inner_join(X_control_df_raw,
  X_metabo_df_raw,
  by = "mouse_id"
) %>%
  dplyr::arrange(mouse_id) %>%
  dplyr::distinct(mouse_id, .keep_all = TRUE) %>%
  tidyr::drop_na(fasting, bw_test, diet)

# make test and train sets
train_frac <- 0.8
X_final <- docr_elastic_train_test_parse(final_data,
                                         train_frac = train_frac,
                                         seed = 212)

print(paste0("Data parsed for training fraction = ", train_frac))

# Initial fit to determine lambda.1se
penalty_factor <- rep(1, ncol(X_final$X_metabo_matrix_train))
set.seed(212)
cv_fit_initial <- docr_elastic_model_fit(
  i = 1,
  alpha_value = model_alpha,
  penalty_vec = penalty_factor,
  X_data = X_final$X_metabo_matrix_train,
  Y_data = X_final$surv_obj_train
)
lambda_1se <- cv_fit_initial$lambda.1se
print(paste0("lambda.1se from initial fit: ", lambda_1se))

# bootstrap result
workers <- min(future::availableCores() - 1, 32)

print(paste0("Starting bootstrapping for alpha = ", model_alpha))
print(paste0("Bootstrapping ", n_bootstrap, " times"))
print(paste0("Running with ", workers, " cores"))
print(paste0("Using fixed lambda.1se = ", lambda_1se))

future::plan(future::multisession, workers = workers)
coef_list <- furrr::future_map(
  .x = 1:n_bootstrap,
  .f = ~ docr_elastic_model_fit(
    i = .x,
    alpha_value = model_alpha,
    penalty_vec = penalty_factor,
    X_data = X_final$X_metabo_matrix_train,
    Y_data = X_final$surv_obj_train,
    bootstrap = TRUE,
    fixed_lambda = lambda_1se
  ),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

coef_bootstrap <- do.call(rbind, coef_list)
colnames(coef_bootstrap) <- colnames(X_final$X_metabo_matrix_train)

# Refitted CoxPh on features selected > 60% of bootstrap iterations
selection_freq <- apply(coef_bootstrap, 2, function(col) mean(col != 0))
selected_features <- names(selection_freq)[selection_freq > 0.60]
print(paste0("Features selected > 60% of the time: ", length(selected_features)))

coxph_summary <- data.frame(variable = character(), coxph_pvalue = numeric())
if (length(selected_features) > 0) {
  refitted_cox <- survival::coxph(y ~ ., data = data.frame(
    y = X_final$surv_obj_train,
    x = X_final$X_metabo_matrix_train[, selected_features]
  ))
  coxph_summary <- summary(refitted_cox)$coefficients %>%
    as.data.frame() %>%
    dplyr::mutate(variable = selected_features) %>%
    dplyr::mutate(coxph_pvalue = `Pr(>|z|)`) %>%
    as.data.frame() %>%
    dplyr::select(variable, coxph_pvalue)
}

bootstrap_results <- data.frame(
  variable = colnames(X_final$X_metabo_matrix_train),
  coef_estimate = coef(cv_fit_initial, s = "lambda.1se") %>% as.numeric(),
  ci_lower = apply(coef_bootstrap, 2, function(col) {
    quantile(x = col, probs = 0.025)
  }),
  ci_upper = apply(coef_bootstrap, 2, function(col) {
    quantile(x = col, probs = 0.975)
  }),
  se_coefs = apply(coef_bootstrap, 2, sd),
  n_selection = apply(coef_bootstrap, 2, function(col) {
    length(col[col != 0])
  })
) %>%
  dplyr::mutate(
    per_selection = n_selection / nrow(coef_bootstrap),
    z_scores = coef_estimate / se_coefs,
    zscore_pvalue = 2 * pnorm(-abs(z_scores)),
    ci_significant = (ci_lower > 0 & ci_upper > 0) | (ci_lower < 0 & ci_upper < 0)
  ) %>%
  dplyr::left_join(coxph_summary, by = "variable") %>%
  dplyr::arrange(zscore_pvalue)

bootstrap_file_name1 <- file.path(output_filepath, paste0(
  gsub("-", "", Sys.Date()),
  "_LassoCalcs_Alpha", model_alpha,
  "_Nboot", n_bootstrap, "ALL_PLL.Rds"
))
bootstrap_file_name2 <- file.path(output_filepath, paste0(
  gsub("-", "", Sys.Date()),
  "_LassoCoefs_Alpha", model_alpha,
  "_Nboot", n_bootstrap, "ALL_PLL.Rds"
))
saveRDS(bootstrap_results, bootstrap_file_name1)
saveRDS(coef_bootstrap, bootstrap_file_name2)

print(paste0("Bootstrap calculations saved as ", bootstrap_file_name1))
print(paste0("Bootstrap coefficients saved as ", bootstrap_file_name2))
