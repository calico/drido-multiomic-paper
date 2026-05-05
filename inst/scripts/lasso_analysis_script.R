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
#   stop("Need at least 4 arguments: local_filepath, output_filepath, n_bootstrap, model_alpha", call. = FALSE)
# }
#
# local_filepath <- args[1]
# output_filepath <- args[2]
# n_bootstrap <- as.numeric(args[3])
# model_alpha <- as.numeric(args[4])


## arguments
# n_bootstrap <- 5000
# model_alpha <- 1

# ### Assign file paths
# local_filepath <- "~/drido-multiomic-paper"
# output_filepath <- "~/"

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
    PLL < 0.85,
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
X_final <- docr_elastic_train_test_parse(final_data, train_frac = train_frac)

print(paste0("Data parsed for training fraction = ", train_frac))

# define alphas
# alphas <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
alphas <- c(1)

# remove anything with low variance
penalty_factor <- rep(1, ncol(X_final$X_metabo_matrix_train))

# run elastic net models
workers <- future::availableCores() - 1
print(paste0("Running initial cv.glmnet model with ", workers, " cores"))

set.seed(21)
future::plan(future::multisession(), workers = workers)
results_list <- furrr::future_map(
  .x = alphas,
  .f = ~ docr_elastic_model_fit(
    alpha_value = .x,
    penalty_vec = penalty_factor,
    X_data = X_final$X_metabo_matrix_train,
    Y_data = X_final$surv_obj_train
  ),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

# analyze results
names(results_list) <- paste0("alpha", alphas)

diagnostic_analysis <- TRUE

if (diagnostic_analysis) {
  # diagnostic plots
  for (alph in alphas) {
    plot(results_list[[paste0("alpha", alph)]]$glmnet.fit,
      xvar = "lambda",
      main = paste0("Alpha = ", alph)
    )
  }
  for (alph in alphas) {
    plot(results_list[[paste0("alpha", alph)]],
      main = paste0("Min CVM for Alpha ", alph, " = ", min(results_list[[paste0("alpha", alph)]]$cvm))
    )
  }

  # predictions
  for (alph in alphas) {
    fit.name <- paste0("alpha", alph)
    cv_model <- results_list[[fit.name]]

    # make predictions on all lambda values
    all_train_predict <- predict(cv_model,
      newx = X_final$X_metabo_matrix_train,
      s = cv_model$lambda,
      type = "response"
    ) %>%
      as.data.frame()
    colnames(all_train_predict) <- paste0("lambda", cv_model$lambda)
    c_indices <- apply(all_train_predict, 2, function(pred_vec) {
      glmnet::Cindex(pred = pred_vec, y = X_final$surv_obj_train)
    })

    # isolate lambda values
    best_lambda_min <- cv_model$lambda.min
    best_lambda_1se <- cv_model$lambda.1se
    best_lambda_ci <- gsub("lambda", "", names(c_indices)[which.max(c_indices)]) %>% as.numeric()

    # select predictions from viable lamdba values
    all_train_predict <- all_train_predict %>%
      dplyr::select(c(
        paste0("lambda", best_lambda_min),
        paste0("lambda", best_lambda_1se),
        paste0("lambda", best_lambda_ci)
      ))
    colnames(all_train_predict) <- c("lambda.min", "lambda.1se", "lambda.ci")
    all_train_predict <- all_train_predict %>%
      tibble::rownames_to_column("mouse_id") %>%
      tidyr::gather(-mouse_id, key = "lambda", value = "predicted_log_hr") %>%
      dplyr::right_join(
        X_final$X_control_df_train %>%
          as.data.frame(),
        by = "mouse_id"
      )

    # make predictions on test data
    all_test_predict <- predict(cv_model,
      newx = X_final$X_metabo_matrix_test,
      s = c(best_lambda_min, best_lambda_1se, best_lambda_ci)
    ) %>%
      as.data.frame()

    colnames(all_test_predict) <- c("lambda.min", "lambda.1se", "lambda.ci")

    all_test_predict <- all_test_predict %>%
      tibble::rownames_to_column("mouse_id") %>%
      tidyr::gather(-mouse_id, key = "lambda", value = "predicted_log_hr") %>%
      dplyr::right_join(
        X_final$X_control_df_test %>%
          as.data.frame(),
        by = "mouse_id"
      ) %>%
      dplyr::mutate(
        hazard_ratio = exp(predicted_log_hr),
        survival_score = -predicted_log_hr,
        se = (surv_days - survival_score)^2,
        mse = mean(se, na.rm = TRUE)
      ) %>%
      dplyr::mutate(days_remaining = surv_days - age_days)

    g <- ggplot(
      all_test_predict,
      aes(x = surv_days, y = survival_score)
    ) +
      facet_grid(cols = vars(diet), rows = vars(lambda), scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )) +
      labs(
        y = "-(Hazard Ratio)",
        title = paste0("Alpha = ", alph)
      ) +
      geom_point() +
      geom_smooth(method = "lm")
    print(g)


    # make predictions on train data (supplement)
    all_train_predict <- predict(cv_model,
      newx = X_final$X_metabo_matrix_train,
      s = c(best_lambda_min, best_lambda_1se, best_lambda_ci)
    ) %>%
      as.data.frame()

    colnames(all_train_predict) <- c("lambda.min", "lambda.1se", "lambda.ci")

    all_train_predict <- all_train_predict %>%
      tibble::rownames_to_column("mouse_id") %>%
      tidyr::gather(-mouse_id, key = "lambda", value = "predicted_log_hr") %>%
      dplyr::right_join(
        X_final$X_control_df_train %>%
          as.data.frame(),
        by = "mouse_id"
      ) %>%
      dplyr::mutate(
        hazard_ratio = exp(predicted_log_hr),
        survival_score = -predicted_log_hr,
        se = (surv_days - survival_score)^2,
        mse = mean(se, na.rm = TRUE)
      ) %>%
      dplyr::mutate(days_remaining = surv_days - age_days)

    g <- ggplot(
      all_train_predict,
      aes(x = surv_days, y = survival_score)
    ) +
      facet_grid(cols = vars(diet), rows = vars(lambda), scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )) +
      labs(
        y = "-(Hazard Ratio)",
        title = paste0("Alpha = ", alph)
      ) +
      geom_point() +
      geom_smooth(method = "lm")
    print(g)
  }

  diet_names <- c("All Diets", "AL", "IF-1D", "IF-2D", "CR-20", "CR-40")
  diet_colors <- c("grey20", docr_get_diet_colors())
  names(diet_colors) <- diet_names

  # Fig 4D - Test set
  test_data_plot <- dplyr::bind_rows(
    all_test_predict %>%
      dplyr::filter(lambda == "lambda.min") %>%
      dplyr::mutate(diet = ifelse(grepl("20|40", diet), paste0("CR-", diet),
        ifelse(grepl("1D|2D", diet), paste0("IF-", diet),
          "AL"
        )
      )),
    all_test_predict %>%
      dplyr::filter(lambda == "lambda.min") %>%
      dplyr::mutate(diet = "All Diets")
  ) %>%
    dplyr::mutate(diet = factor(diet, levels = diet_names))

  temp_label <- docr_facet_stats(test_data_plot,
    value_x = "days_remaining",
    value_y = "predicted_log_hr",
    facet_1 = "diet"
  )

  n_test <- X_final$X_control_df_test$mouse_id %>%
    unique() %>%
    length()

  g <- ggplot(
    test_data_plot,
    aes(x = days_remaining, y = predicted_log_hr, color = diet)
  ) +
    facet_wrap(~diet, scales = "free_x", nrow = 1) +
    theme_bw() +
    labs(
      y = "log(HR)",
      x = "Days of Life Remaining",
      title = paste0("Predictions on held out test set (20%, n = ", n_test, ")")
    ) +
    geom_point(size = 0.4, alpha = 0.5) +
    scale_color_manual(name = "Diet", values = diet_colors) +
    geom_smooth(method = "lm", show.legend = FALSE) +
    docr_ggplot_theme() +
    docr_ggplot_stats_label(temp_label,
      y = Inf,
      vjust = 1.3,
      size = 5 / ggplot2::.pt,
      color = "black"
    ) +
    theme(
      strip.text = element_text(margin = margin(0, 2, 2, 2)),
      panel.spacing.y = unit(3, "pt")
    ) +
    coord_cartesian(ylim = c(-0.29, 1))
  print(g)


  # docr_ggsave(plot_object = g,
  #             "Fig_4D_Test_Predictions",
  #             plot_width = 5,
  #             plot_height = 1.55,
  #             local_filepath = local_filepath)


  # Supplement Figure - Train Set
  temp_label <- docr_facet_stats(
    all_train_predict %>%
      dplyr::filter(lambda == "lambda.min"),
    value_x = "days_remaining",
    value_y = "predicted_log_hr",
    facet_1 = "diet"
  )

  n_test <- X_final$X_control_df_train$mouse_id %>%
    unique() %>%
    length()

  g <- ggplot(
    all_train_predict %>%
      dplyr::filter(lambda == "lambda.min"),
    aes(x = days_remaining, y = predicted_log_hr, color = diet)
  ) +
    facet_wrap(~diet, scales = "free_x", ncol = 5) +
    theme_bw() +
    labs(
      y = "log(HR)",
      x = "Days of Life Remaining",
      title = paste0("Predictions on training set (80%, n = ", n_test, ")")
    ) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(name = "Diet", values = docr_get_diet_colors()) +
    geom_smooth(method = "lm", show.legend = FALSE) +
    docr_ggplot_theme() +
    docr_ggplot_stats_label(temp_label,
      y = Inf,
      vjust = 1.3,
      size = 5 / ggplot2::.pt,
      color = "black"
    )
  print(g)

  # docr_ggsave(plot_object = g,
  #             "Train_Predictions",
  #             plot_width = 4.7,
  #             plot_height = 1.7,
  #             local_filepath = local_filepath)
} else {
  print("Skipping diagnostics")
}

# bootstrap result
workers <- future::availableCores() - 1

print(paste0("Starting bootstrapping for alpha = ", model_alpha))
print(paste0("Bootstrapping ", n_bootstrap, " times"))
print(paste0("Running with ", workers, " cores"))

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
    bootstrap_lamdba = "lambda.min"
  ),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

coef_bootstrap <- do.call(rbind, coef_list)


# Refitted CoxPh on lambda min for p-values -----
coefs_temp <- coef(results_list[["alpha1"]], s = "lambda.min")
selected_features <- rownames(coefs_temp)[which(coefs_temp != 0)]
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
  coef_estimate = coef(results_list[["alpha1"]], s = "lambda.min") %>% as.numeric(),
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

colnames(coef_bootstrap) <- colnames(X_final$X_metabo_matrix_train)

bootstrap_file_name1 <- file.path(output_filepath, paste0(
  gsub("-", "", Sys.Date()),
  "_LassoCalcs_Alpha", model_alpha,
  "_Nboot", n_bootstrap, ".Rds"
))
bootstrap_file_name2 <- file.path(output_filepath, paste0(
  gsub("-", "", Sys.Date()),
  "_LassoCoefs_Alpha", model_alpha,
  "_Nboot", n_bootstrap, ".Rds"
))
saveRDS(bootstrap_results, bootstrap_file_name1)
saveRDS(coef_bootstrap, bootstrap_file_name2)

print(paste0("Bootstrap calculations saved as ", bootstrap_file_name1))
print(paste0("Bootstrap coefficients saved as ", bootstrap_file_name2))
