library(glmnet)
library(survival)
library(dplyr)
library(purrr)
library(furrr)
library(limma)
library(future)

## arguments
model_alpha <- 1

# Stability selection: pool bootstrap draws across many independent
# train/test splits (seeds), not just one, due to high variability by split.
# Less bootstraps per seed and then taking an average of all seeds
# Leakage determined to be less problematic than train/test seed-bias
# Changed original (arbitrary) seed to starting seed = 1 to reduce bias
seed_use <- 1
n_seeds <- 20
seeds_use <- seed_use + seq_len(n_seeds) - 1
n_bootstrap <- 200

# ### Assign file paths
local_filepath <- "~/workspace/drido-multiomic-paper"
output_filepath <- "~/workspace/docr_data/MS1553"

# local_filepath <- "~/GitHub/drido-multiomic-paper"
# output_filepath <- "~/GitHub/drido-multiomic-paper/inst/extdata"

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

# Sensitivity configurations -- the four core lambda.min configs from
# lasso_analysis_script_final_2.R (with/without BW x with/without trim).
# lambda.1se is a separate validation axis (see that script's _1se configs)
# and isn't crossed with multi-seed stability selection here.
sensitivity_configs <- list(
  list(label = "All_Year2_wBW", filter_expr = TRUE, regress_bw = TRUE),
  list(label = "All_Year2_noBW", filter_expr = TRUE, regress_bw = FALSE),
  list(label = "Trim_Year2_wBW", filter_expr = quote(PLL <= 0.86), regress_bw = TRUE),
  list(label = "Trim_Year2_noBW", filter_expr = quote(PLL <= 0.86), regress_bw = FALSE)
)

# Parallel processing
workers <- future::availableCores() - 1
future::plan(future::multisession, workers = workers)

for (config in sensitivity_configs) {

  config_label <- config$label
  print(paste0("===== Running configuration: ", config_label, " ====="))

  # filter data for regularization
  # molecular traits only; year2 only -- identical across seeds, built once
  dat <- data_use %>%
    as.data.frame() %>%
    dplyr::filter(
      !is.na(age_years),
      age_years == "year2",
      modality != "physiological"
    ) %>%
    dplyr::filter(eval(config$filter_expr)) %>%
    dplyr::mutate(age_days = surv_days - days_remaining)

  print(paste0("Filtered to ", dplyr::n_distinct(dat$mouse_id), " mice"))

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

  train_frac <- 0.8

  # Accumulators across seeds
  coef_bootstrap_list <- vector("list", n_seeds)
  coef_estimate_list <- vector("list", n_seeds)
  seed_summary_list <- vector("list", n_seeds)

  for (i in seq_along(seeds_use)) {
    seed_i <- seeds_use[i]
    print(paste0("--- ", config_label, ": seed ", seed_i, " (", i, "/", n_seeds, ") ---"))

    # make test and train sets
    X_final <- docr_elastic_train_test_parse(final_data,
                                             train_frac = train_frac,
                                             seed = seed_i,
                                             regress_bw = config$regress_bw)

    # Initial fit to determine lambda.min
    penalty_factor <- rep(1, ncol(X_final$X_metabo_matrix_train))
    set.seed(seed_i)
    cv_fit_initial <- docr_elastic_model_fit(
      i = 1,
      alpha_value = model_alpha,
      penalty_vec = penalty_factor,
      X_data = X_final$X_metabo_matrix_train,
      Y_data = X_final$surv_obj_train
    )

    fixed_lambda <- cv_fit_initial$lambda.min
    print(paste0("lambda.min from seed ", seed_i, ": ", fixed_lambda))

    # bootstrap result
    set.seed(seed_i)
    coef_list <- furrr::future_map(
      .x = 1:n_bootstrap,
      .f = ~ docr_elastic_model_fit(
        i = .x,
        alpha_value = model_alpha,
        penalty_vec = penalty_factor,
        X_data = X_final$X_metabo_matrix_train,
        Y_data = X_final$surv_obj_train,
        bootstrap = TRUE,
        fixed_lambda = fixed_lambda
      ),
      .options = furrr::furrr_options(seed = TRUE)
    )

    coef_bootstrap_seed <- do.call(rbind, coef_list)
    colnames(coef_bootstrap_seed) <- colnames(X_final$X_metabo_matrix_train)

    # selection frequency calculation for the seed
    selection_freq_seed <- apply(coef_bootstrap_seed, 2, function(col) mean(col != 0))

    seed_summary_list[[i]] <- data.frame(
      variable = colnames(X_final$X_metabo_matrix_train),
      seed = seed_i,
      coef_estimate_seed = coef(cv_fit_initial, s = "lambda.min") %>% as.numeric(),
      per_selection_seed = selection_freq_seed,
      lambda_min_seed = fixed_lambda
    )

    coef_bootstrap_list[[i]] <- coef_bootstrap_seed
    coef_estimate_list[[i]] <- coef(cv_fit_initial, s = "lambda.min") %>% as.numeric()
  }

  # Pool bootstrap coefficients across all seeds
  # Selection frequency is computed on the full pooled set
  coef_bootstrap <- do.call(rbind, coef_bootstrap_list)
  seed_summary <- dplyr::bind_rows(seed_summary_list)

  print(paste0(
    "Pooled ", nrow(coef_bootstrap), " bootstrap draws across ", n_seeds,
    " seeds (", n_bootstrap, " each) for ", config_label
  ))

  # Coxph p-value refit on final feature set determined from the pooled set
  selection_freq <- apply(coef_bootstrap, 2, function(col) mean(col != 0))
  selected_features <- names(selection_freq)[selection_freq > 0.60]
  print(paste0("Features selected > 60% of pooled iterations: ", length(selected_features)))

  # selection-stability diagnostic: of the n_seeds individual seeds, what
  # fraction also independently cleared 60% selection on their own?
  seed_robustness <- seed_summary %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(frac_seeds_selected = mean(per_selection_seed > 0.60), .groups = "drop")

  # coxPH p-values for the final (pooled-average) selected feature set,
  # refit once on the full dataset rather than per-seed.
  # docr_elastic_train_test_parse(train_frac = 1) reuses the same cleaning
  # (removeBatchEffect) pipeline and returns every mouse in "train"
  coxph_summary <- data.frame(variable = character(), coxph_pvalue = numeric())
  if (length(selected_features) > 0) {
    X_full <- docr_elastic_train_test_parse(final_data,
                                            train_frac = 1,
                                            seed = seed_use,
                                            regress_bw = config$regress_bw)
    refitted_cox <- survival::coxph(y ~ ., data = data.frame(
      y = X_full$surv_obj_train,
      x = X_full$X_metabo_matrix_train[, selected_features, drop = FALSE]
    ))
    coxph_summary <- summary(refitted_cox)$coefficients %>%
      as.data.frame() %>%
      dplyr::mutate(variable = selected_features, coxph_pvalue = `Pr(>|z|)`) %>%
      as.data.frame() %>%
      dplyr::select(variable, coxph_pvalue)
  }

  bootstrap_results <- data.frame(
    variable = colnames(coef_bootstrap),
    coef_estimate = rowMeans(do.call(cbind, coef_estimate_list)),
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
    dplyr::left_join(seed_robustness, by = "variable") %>%
    dplyr::left_join(coxph_summary, by = "variable") %>%
    dplyr::arrange(zscore_pvalue)

  bootstrap_file_name1 <- file.path(output_filepath, paste0(
    gsub("-", "", Sys.Date()),
    "_LassoCalcs_StabilitySelection_Alpha", model_alpha,
    "_Nseeds", n_seeds, "_Nboot", n_bootstrap, "_", config_label, ".Rds"
  ))
  bootstrap_file_name2 <- file.path(output_filepath, paste0(
    gsub("-", "", Sys.Date()),
    "_LassoCoefs_StabilitySelection_Alpha", model_alpha,
    "_Nseeds", n_seeds, "_Nboot", n_bootstrap, "_", config_label, ".Rds"
  ))
  seed_summary_file_name <- file.path(output_filepath, paste0(
    gsub("-", "", Sys.Date()),
    "_LassoSeedSummary_StabilitySelection_Alpha", model_alpha,
    "_Nseeds", n_seeds, "_Nboot", n_bootstrap, "_", config_label, ".Rds"
  ))
  saveRDS(bootstrap_results, bootstrap_file_name1)
  saveRDS(coef_bootstrap, bootstrap_file_name2)
  saveRDS(seed_summary, seed_summary_file_name)

  print(paste0("Pooled bootstrap calculations saved as ", bootstrap_file_name1))
  print(paste0("Pooled bootstrap coefficients saved as ", bootstrap_file_name2))
  print(paste0("Per-seed summary saved as ", seed_summary_file_name))
}

future::plan(future::sequential)

print("Complete")
