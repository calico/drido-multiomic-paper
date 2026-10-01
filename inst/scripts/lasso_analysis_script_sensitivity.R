library(glmnet)
library(survival)
library(dplyr)
library(purrr)
library(furrr)
library(limma)
library(ggplot2)
library(future)

## arguments
model_alpha <- 1
n_bootstrap <- 1000

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

# Sensitivity configurations:
# 1) All year2 data (no end-life trim)
# 2) Year2 with days_remaining > 100
sensitivity_configs <- list(
  list(label = "ALL_Year2", filter_expr = TRUE, lambda_use = "lambda.min", regress_bw = TRUE),
  list(label = "ALL_Year2_noBW", filter_expr = TRUE, lambda_use = "lambda.min", regress_bw = FALSE),
  list(label = "ALL_Year2_1se", filter_expr = TRUE, lambda_use = "lambda.1se", regress_bw = TRUE),
  list(label = "DR100", filter_expr = quote(days_remaining > 100), lambda_use = "lambda.min", regress_bw = TRUE),
  list(label = "DR120", filter_expr = quote(days_remaining > 120), lambda_use = "lambda.min", regress_bw = TRUE)
)

# parallelize
workers <- future::availableCores() - 1
future::plan(future::multisession, workers = workers)

for (config in sensitivity_configs) {

  config_label <- config$label
  print(paste0("===== Running configuration: ", config_label, " ====="))

  # filter data for regularization
  # molecular traits only; year2 only
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

  # make test and train sets
  train_frac <- 0.8
  X_final <- docr_elastic_train_test_parse(final_data,
                                           train_frac = train_frac,
                                           seed = 212,
                                           regress_bw = config$regress_bw)

  print(paste0("Data parsed for training fraction = ", train_frac))

  # Initial fit to determine lambda.min
  penalty_factor <- rep(1, ncol(X_final$X_metabo_matrix_train))
  set.seed(212)
  cv_fit_initial <- docr_elastic_model_fit(
    i = 1,
    alpha_value = model_alpha,
    penalty_vec = penalty_factor,
    X_data = X_final$X_metabo_matrix_train,
    Y_data = X_final$surv_obj_train
  )

  lambda_min <- cv_fit_initial$lambda.min
  lambda_1se <- cv_fit_initial$lambda.1se
  lambda_use <- config$lambda_use
  fixed_lambda <- cv_fit_initial[[lambda_use]]
  print(paste0(lambda_use, " from initial fit: ", fixed_lambda))

  # bootstrap result
  print(paste0("Starting bootstrapping for alpha = ", model_alpha))
  print(paste0("Bootstrapping ", n_bootstrap, " times"))
  print(paste0("Running with ", workers, " cores"))
  print(paste0("Using fixed ", lambda_use, " = ", fixed_lambda))

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
      x = X_final$X_metabo_matrix_train[, selected_features, drop = FALSE]
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
    coef_estimate = coef(cv_fit_initial, s = lambda_use) %>% as.numeric(),
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
    "_Nboot", n_bootstrap, "_", config_label, ".Rds"
  ))
  bootstrap_file_name2 <- file.path(output_filepath, paste0(
    gsub("-", "", Sys.Date()),
    "_LassoCoefs_Alpha", model_alpha,
    "_Nboot", n_bootstrap, "_", config_label, ".Rds"
  ))
  saveRDS(bootstrap_results, bootstrap_file_name1)
  saveRDS(coef_bootstrap, bootstrap_file_name2)

  print(paste0("Bootstrap calculations saved as ", bootstrap_file_name1))
  print(paste0("Bootstrap coefficients saved as ", bootstrap_file_name2))

  # ---- Diagnostic Plots ----
  plot_list <- list()

  # CV curve (base R — wrap with recordPlot)
  plot(cv_fit_initial)
  title(paste0(config_label, " - Min CVM = ", round(min(cv_fit_initial$cvm), 3)), line = 2.5)
  plot_list[["cv_curve"]] <- recordPlot()

  # Regularization paths
  plot_list[["reg_paths"]] <- docr_regularization_paths(
    cv_fit_initial, name_conversion_use, ylim = c(-1, 1),
    plot_title = paste0("Regularization Paths - ", config_label))

  # C-index across lambda
  all_train_pred_ci <- predict(cv_fit_initial,
    newx = X_final$X_metabo_matrix_train,
    s = cv_fit_initial$lambda, type = "response") %>% as.data.frame()
  colnames(all_train_pred_ci) <- paste0("lambda", cv_fit_initial$lambda)
  c_indices <- apply(all_train_pred_ci, 2, function(p) {
    glmnet::Cindex(pred = p, y = X_final$surv_obj_train)
  })
  c_at_min <- c_indices[[paste0("lambda", lambda_min)]]

  plot_list[["c_index"]] <- ggplot(data.frame(log_lambda = log(cv_fit_initial$lambda),
                          c_index = as.numeric(c_indices)),
               aes(x = log_lambda, y = c_index)) +
    geom_line() +
    geom_vline(xintercept = log(lambda_min), linetype = "dashed", color = "#2166AC") +
    geom_vline(xintercept = log(lambda_1se), linetype = "dotted", color = "red") +
    annotate("text", x = log(lambda_min), y = min(c_indices),
             label = "lambda*'.min'", parse = TRUE,
             vjust = 1.5, hjust = -0.1, size = 3, color = "#2166AC") +
    annotate("text", x = log(lambda_1se), y = min(c_indices),
             label = "lambda*'.1se'", parse = TRUE,
             vjust = 3, hjust = -0.1, size = 3, color = "red") +
    labs(x = expression(log(lambda)), y = "C-index",
         title = bquote(.(config_label) ~ "- C-index (C at" ~ lambda*.min ~ "=" ~ .(round(c_at_min, 3)) * ")")) +
    theme_classic()

  # Train and test predictions by diet
  diet_names <- c("All Diets", "AL", "IF-1D", "IF-2D", "CR-20", "CR-40")
  diet_colors <- c("grey20", docr_get_diet_colors())
  names(diet_colors) <- diet_names

  for (set_label in c("Train", "Test")) {
    if (set_label == "Train") {
      X_mat <- X_final$X_metabo_matrix_train
      control_df <- X_final$X_control_df_train
    } else {
      X_mat <- X_final$X_metabo_matrix_test
      control_df <- X_final$X_control_df_test
    }

    pred_df <- predict(cv_fit_initial, newx = X_mat, s = fixed_lambda) %>%
      as.data.frame() %>%
      setNames("predicted_log_hr") %>%
      tibble::rownames_to_column("mouse_id") %>%
      dplyr::right_join(control_df %>% as.data.frame(), by = "mouse_id") %>%
      dplyr::mutate(
        days_remaining = surv_days - age_days,
        diet = ifelse(grepl("20|40", diet), paste0("CR-", diet),
                      ifelse(grepl("1D|2D", diet), paste0("IF-", diet), "AL"))
      )
    pred_diet <- dplyr::bind_rows(pred_df, pred_df %>% dplyr::mutate(diet = "All Diets")) %>%
      dplyr::mutate(diet = factor(diet, levels = diet_names))

    temp_label <- docr_facet_stats(pred_diet,
      value_x = "days_remaining", value_y = "predicted_log_hr", facet_1 = "diet")

    plot_list[[paste0("pred_", tolower(set_label))]] <-
      ggplot(pred_diet, aes(x = days_remaining, y = predicted_log_hr, color = diet)) +
      facet_wrap(~diet, scales = "free_x", nrow = 1) +
      geom_point(size = 0.4, alpha = 0.5) +
      scale_color_manual(name = "Diet", values = diet_colors) +
      geom_smooth(method = "lm", show.legend = FALSE) +
      docr_ggplot_theme() +
      docr_ggplot_stats_label(temp_label,
        y = Inf, vjust = 1.3, size = 5 / ggplot2::.pt, color = "black") +
      labs(y = "log(HR)", x = "Days of Life Remaining",
           title = paste0(config_label, " - ", set_label,
                          " Set at ", lambda_use, " (n=", nrow(control_df), ")"))
  }

  # Bootstrap coefficient distributions
  per_select_thresh <- 0.60
  lasso_coef_clean <- coef_bootstrap %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    dplyr::mutate(feature_id = gsub("^`|`$", "", feature_id)) %>%
    dplyr::inner_join(name_conversion_use, by = "feature_id") %>%
    tidyr::pivot_longer(cols = -c("feature_id", "name_use"),
                        names_to = "model_id", values_to = "coef") %>%
    dplyr::group_by(name_use, feature_id) %>%
    dplyr::summarise(per_selection = sum(coef != 0) / n_bootstrap, .groups = "drop")

  selected_names <- lasso_coef_clean %>%
    dplyr::filter(per_selection > per_select_thresh) %>%
    dplyr::arrange(desc(per_selection)) %>%
    dplyr::pull(name_use)

  per_selection_red <- 0.80

  tryCatch({
    if (length(selected_names) > 0) {
      lasso_results_clean <- bootstrap_results %>%
        dplyr::mutate(feature_id = gsub("^`|`$", "", variable)) %>%
        dplyr::inner_join(name_conversion_use, by = "feature_id")

      box_data <- coef_bootstrap %>%
        t() %>% as.data.frame() %>%
        tibble::rownames_to_column("feature_id") %>%
        dplyr::mutate(feature_id = gsub("^`|`$", "", feature_id)) %>%
        dplyr::inner_join(name_conversion_use, by = "feature_id") %>%
        tidyr::pivot_longer(cols = -c("feature_id", "name_use"),
                            names_to = "model_id", values_to = "coef") %>%
        dplyr::left_join(lasso_results_clean,
                          by = c("name_use", "feature_id")) %>%
        dplyr::filter(name_use %in% selected_names) %>%
        dplyr::mutate(
          name_use = factor(name_use, levels = rev(selected_names)),
          group = ifelse(coxph_pvalue < 0.05 & per_selection > per_selection_red, "Sig", "Insig")
        )

      plot_list[["coef_dist"]] <- ggplot(box_data, aes(y = name_use, x = coef, color = group)) +
        geom_boxplot(outliers = FALSE, fill = NA,
                     outlier.alpha = 0.2, outlier.fill = NA, outlier.shape = 21) +
        geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
        scale_color_manual(values = c("Sig" = "#009E73", "Insig" = "grey80")) +
        labs(x = paste0("Bootstrap Coefficients (n=", n_bootstrap, ")"), y = "",
             title = paste0(config_label, " - Coefficient Distributions (>",
                            per_select_thresh * 100, "% selected)")) +
        docr_ggplot_theme(legend_position = "none")

      # Selection frequency bars
      bar_data <- lasso_results_clean %>%
        dplyr::filter(name_use %in% selected_names) %>%
        dplyr::mutate(
          name_use = factor(name_use, levels = rev(selected_names)),
          group = ifelse(coxph_pvalue < 0.05 & per_selection > per_selection_red, "Sig", "Insig")
        )

      plot_list[["sel_freq"]] <- ggplot(bar_data, aes(y = name_use, x = per_selection * 100,
                                 color = group, fill = group)) +
        geom_col(width = 0.7) +
        geom_vline(xintercept = per_selection_red * 100, linetype = "dashed", color = "grey50") +
        scale_color_manual(values = c("Sig" = "#009E73", "Insig" = "grey80")) +
        scale_fill_manual(values = c("Sig" = "#009E73", "Insig" = "grey80")) +
        scale_x_continuous(breaks = c(0, per_selection_red * 100, 100)) +
        labs(x = "% Selection", y = "",
             title = paste0(config_label, " - Selection Frequency")) +
        docr_ggplot_theme(legend_position = "none") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
  }, error = function(e) message("Bootstrap coefficient plots skipped: ", e$message))

  # Save plot list
  plot_file_name <- file.path(output_filepath, paste0(
    gsub("-", "", Sys.Date()),
    "_LassoDiagnosticPlots_Alpha", model_alpha, "_Nboot", n_bootstrap, "_", config_label, ".Rds"
  ))
  saveRDS(plot_list, plot_file_name)
  print(paste0("Diagnostic plots saved as ", plot_file_name))
}

future::plan(future::sequential)

print("Complete")
