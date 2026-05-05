#' Linear mixed modeling
#'
#' @description
#' This function trains a linear mixed model through \code{lmer} as specified
#'
#' @param data_use a dataframe with \code{value_var} numeric values
#' corresponding to the feature identified by \code{feature_id_test},
#' containing columns associated with covariates in \code{forms_to_test}
#' @param feature_id_test name of feature being tested. This is not a string
#' for filtering, but is used for print messages during model running. Helpful
#' when model fitting features in parallel
#' @param forms_to_test formulates to test, as character vector. Only include
#' covariates to the right of the tilde
#' @param value_var column names of dependent values, as a string. Defaults to
#' \code{norm_abundance} and is pasted to the left of the tilde in the formula
#' @param n_samples number of unique samples in experiment. Defaults to
#' \code{NULL} to ignore parameter.
#' @param n_sample_threshold percent of \code{n_samples} below which model
#' fitting will skip. Defaults to 0.1 for 10%.
#' @param model_type type of model to run, one of \code{"mixed"} for
#' \code{lmer()} fit, or \code{"simple"} for \code{lm()} fit
#'
#' @importFrom magrittr %>%
#'
#' @returns a dataframe containing the feature_id, model_term, model, coef,
#' se_coef, p_value_coef, t_value, and AIC of the model fit
#'
#' @export
lmer_multi_formula <- function(data_use,
                               feature_id_test,
                               forms_to_test,
                               value_var = "norm_abundance",
                               n_samples = NULL,
                               n_sample_threshold = 0.1,
                               model_type = "mixed",
                               center_data = TRUE) {
  checkmate::assertNumeric(data_use[[value_var]])

  # Initialized empty dataframe
  da <- data.frame()

  # Drop NA values and center, if assigned
  dfx <- data_use %>%
    tidyr::drop_na(!!rlang::sym(value_var))

  if (center_data) {
    dfx <- dfx %>%
      dplyr::mutate(!!rlang::sym(value_var) := as.numeric(
        scale(!!rlang::sym(value_var),
          scale = FALSE,
          center = TRUE
        )
      ))
  }

  # If there is no data or if data is only present in < n_sample_threshold
  # percent of n_samples, skip modeling
  if (nrow(dfx) == 0) {
    print(paste0(
      "For feature_id ", feature_id_test,
      ", skipping modeling; no data points"
    ))
    return(da)
  }
  if (!is.null(n_samples)) {
    if (nrow(dfx) < n_samples * n_sample_threshold) {
      print(paste0(
        "For feature_id ", feature_id_test,
        paste0(
          ", skipping modeling; data present in < ",
          as.character(n_sample_threshold * 100), " % of samples"
        )
      ))
      return(da)
    }
  }

  # Run analysis through all formulas
  for (f_use in forms_to_test) {
    # f_use <- f
    #   # Properly format equation
    #   if (!grepl("~", f)) {
    #     f_use <- paste0("norm_abundance ~ ", f)
    #   } else {
    #     f_use <- paste0("norm_abundance ", f)
    #   }

    if (model_type == "mixed") {
      # If only one data point per mouse id, change random intercept to
      # generation_wave
      # Will skip if these columns are not in dataframe (for DO-CR)
      if (all(c("mouse_id", "generation_wave") %in% colnames(dfx))) {
        if (length(unique(dfx$mouse_id)) == nrow(dfx)) {
          f_use <- gsub("mouse_id", "generation_wave", f_use)
          print(paste0(
            "For feature_id ", feature_id_test, ", formula ",
            f_use, ", random intercept changed to generation_wave"
          ))
        }
      }

      # Run linear mixed model on data subset, for each formula
      model1 <- suppressMessages(lmerTest::lmer(formula = f_use, data = dfx))

      # Extract model summary data, including coefficient
      # P-values derived from the lmerTest library
      dfc <- as.data.frame(coef(summary(model1)))
      colnames(dfc) <- c("coef", "se_coef", "df", "t_value", "p_value_coef")
    } else if (model_type == "simple") {
      # Run simple linear model on data subset, for each formula
      model1 <- stats::lm(formula = f_use, data = dfx)

      # Extract model summary data, including coefficient
      dfc <- as.data.frame(coef(summary(model1)))
      colnames(dfc) <- c("coef", "se_coef", "t_value", "p_value_coef")
    } else {
      print("'model_type' must be one of 'mixed' or 'simple'; skipping modeling")
      return(da)
    }

    # Add other information that we want in DA dataframe
    dfc <- dfc %>%
      tibble::rownames_to_column(var = "model_term") %>%
      dplyr::mutate(
        model_term = factor(model_term),
        feature_id = .env$feature_id_test,
        model = .env$f_use,
        AIC = AIC(.env$model1)
      ) %>%
      dplyr::select(
        feature_id, model_term, model,
        coef, se_coef, p_value_coef, t_value, AIC
      )

    # Bind with results from other models
    da <- da %>%
      dplyr::bind_rows(dfc)
  }

  return(da)
}


#' ANOVA modeling for variance partitioning
#'
#' @description
#' This function runs \code{anova()} on one or more formulas to determine the
#' percent variance explained by each variable.
#'
#' @param data_use a dataframe with \code{value_var} numeric values
#' corresponding to the feature identified by \code{feature_id_test},
#' containing columns associated with covariates in \code{forms_to_test}
#' @param feature_id_test name of feature being tested. This is not a string
#' for filtering, but is used for print messages during model running.
#' Helpful when model fitting features in parallel
#' @param forms_to_test formulates to test, as character vector. Only
#' include covariates to the right of the tilde
#' @param value_var column names of dependent values, as a string. Defaults to
#' \code{norm_abundance} and is pasted to the left of the tilde in the formula
#' @param n_min minimum number of unique samples required to run ANOVA
#' @param model_type type of model to run, one of \code{"mixed"} for
#' \code{lmer()} fit, or \code{"simple"} for \code{lm()} fit
#' @param per_exp_var variable on which to calculate the percent of variance
#' explained, one of \code{"Sum Sq"} for net sum of squares, or
#' \code{"Mean Sq"} for mean sum of squares for each variable (sum of squares
#' divided by the degrees of freedom)
#'
#' @importFrom magrittr %>%
#'
#' @returns a dataframe containing the model terms and percent variance
#' explained from \code{anova()}
#'
#' @export
anova_multi_formula <- function(data_use,
                                feature_id_test,
                                forms_to_test,
                                value_var = "norm_abundance",
                                n_min = 100,
                                model_type = "simple",
                                per_exp_var = "Sum Sq") {
  checkmate::assertNumeric(data_use[[value_var]])

  # Initialized empty dataframe
  da <- data.frame()

  # Drop NA values and assign `norm_abundance` to the value variable
  dfx <- data_use %>%
    tidyr::drop_na(!!rlang::sym(value_var)) %>%
    dplyr::mutate(norm_abundance = scale(!!rlang::sym(value_var), scale = F, center = T))

  # If there is no data or if data is only present in < n_min samples, skip modeling
  if (nrow(dfx) == 0 || nrow(dfx) < n_min) {
    print(paste0("For feature_id ", feature_id_test, ", skipping modeling; not enough data points"))
    return(da)
  }

  # DOCR specific check
  if ("diet" %in% colnames(dfx)) {
    if (length(unique(dfx$diet)) < 5) {
      print(paste0("For feature_id ", feature_id_test, ", skipping modeling; < 5 diets represented"))
      return(da)
    }
  }

  # Run Variance analysis
  for (f in forms_to_test) {
    # Properly format equation
    if (!grepl("~", f)) {
      f_use <- paste0("norm_abundance ~ ", f)
    } else {
      f_use <- paste0("norm_abundance ", f)
    }

    # Run linear model on data subset, and run anova to determine percent
    # of variance explained. The variance can be decomposed on the linear mixed
    # effect model as well, but the residual variance cannot be extracted
    if (model_type == "simple") {
      anova_model <- stats::anova(stats::lm(formula = f_use, data = dfx))
    } else if (model_type == "mixed") {
      anova_model <- stats::anova(lme4::lmer(formula = f_use, data = dfx))
    } else {
      stop("model_type must be one of `simple` for lm() linear model, or `mixed` for lmer() mixed effects model")
    }

    ss <- anova_model[[per_exp_var]]
    percent_explained <- anova_model %>%
      dplyr::bind_cols(PerExp = ss / sum(ss) * 100) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("model_term") %>%
      dplyr::mutate(
        model = .env$f,
        feature_id = .env$feature_id_test
      )

    # Bind with results from other models
    da <- da %>%
      dplyr::bind_rows(percent_explained)
  }

  return(da)
}


#' Calculate qvalues
#'
#' @param term_data dataframe of linear regression containing p.value
#'
#' @param pval_var column name of pvalues in \code{term_data}, as a string.
#' Defaults to \code{"p.value"}
#'
#' @returns dataframe of linear regression with qvalues added
#'
#' @export
fdr_test <- function(term_data,
                     pval_var = "p.value") {
  if (nrow(term_data) == 0) {
    warning("No data; returning empty dataframe")
    return(data.frame())
  }
  if (!all(dplyr::between(term_data[[pval_var]], 0, 1), na.rm = TRUE)) {
    stop("All values in \"pval_var\" must be between 0 and 1")
  }

  p_values <- term_data[[pval_var]]
  q_values <- try(qvalue::qvalue(p_values)$qvalues, silent = TRUE)

  if ("try-error" %in% class(q_values)) {
    # if qvalue fails this is probably because there are no p-values greater
    # than 0.95 (the highest lambda value)
    # if so add a single p-value of 1 to try to combat the problem

    q_values <- try(qvalue::qvalue(c(p_values, 1))$qvalues, silent = TRUE)

    # If q-value STILL won't calculate, perform BH
    if ("try-error" %in% class(q_values)) {
      q_values <- qvalue::qvalue(p_values, pi0 = 1)$qvalues

      # If p_values didn't need to be handled by BH, remove the last
      # value, since it was corrected by adding a 1 to the end
    } else {
      q_values <- q_values[-length(q_values)]
    }
  }

  term_data <- term_data %>%
    dplyr::mutate(qvalue = q_values)
  return(term_data)
}

#' Calculate q-values on multi-model or multi-term dataframes
#'
#' @param term_data dataframe containing pvalues columns requiring
#' adjustment
#'
#' @param pval_var column name of pvalues to adjust, as a string
#'
#' @param nest_vars column name(s) of nesting variables, if there are
#' multiple models, model terms, outcomes, etc. contained in \code{term_data}
#'
#' @param padj_var column name for column with adjusted pvalues
#'
#' @importFrom magrittr %>%
#'
#' @returns dataframe of linear regression with qvalues added
#'
#' @export
fdr_multi <- function(term_data,
                      pval_var = "pvalue",
                      nest_vars = NULL,
                      padj_var = "padj") {
  if (!(pval_var %in% colnames(term_data))) {
    stop("\"pval_var\":", pval_var, ", not present in df")
  }
  if (!all(nest_vars %in% colnames(term_data))) {
    stop(
      "Not all \"nest_vars\":", paste(nest_vars, collapse = ", "),
      ", present as column names in df"
    )
  }

  term_data_return <- term_data %>%
    tidyr::nest(., data = -all_of(nest_vars)) %>%
    dplyr::mutate(adjusted_data = purrr::map(
      .x = data,
      .f = ~ if (nrow(.x) > 2) {
        fdr_test(
          term_data = .x,
          pval_var = pval_var
        )
      } else {
        .x %>%
          dplyr::mutate(qvalue = NA)
      }
    )) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(adjusted_data, keep_empty = TRUE) %>%
    dplyr::rename(`:=`(!!rlang::sym(padj_var), qvalue))

  return(term_data_return)
}


#' Check DOCR factors for statistical analyses
#'
#' @description
#' This is function checks that all required columns are present in a dataframe
#' for statistical analyses, and that factors are correctly ordered.
#'
#' @param df a dataframe with columns corresponding to necessary variables
#'
#' @returns \code{invisible()}
#'
#' @export
docr_check_factors <- function(df) {
  missing_columns <- setdiff(
    c(
      "weekday_collection", "generation_wave",
      "surv_years", "norm_abundance", "mouse_id",
      "BW_Loess", "is_ddm", "age_years", "diet"
    ),
    colnames(df)
  )

  if (length(missing_columns) > 0) {
    stop(paste0(
      "Columns: ", paste(missing_columns, collapse = ", "),
      " are required but not present in dataframe"
    ))
  } else {
    checkmate::assertFactor(df$generation_wave)
    checkmate::assertFactor(df$mouse_id)

    checkmate::assertFactor(df$diet)
    checkmate::assertTRUE(identical(
      levels(df$diet),
      c("AL", "1D", "2D", "20", "40")
    ))

    checkmate::assertFactor(df$age_years)
    checkmate::assertTRUE(identical(
      levels(df$age_years),
      c("year1", "year2", "year3")
    ))

    checkmate::assertFactor(df$weekday_collection)
    checkmate::assertTRUE(identical(
      levels(df$weekday_collection),
      c("TWT", "Monday")
    ))

    checkmate::assertNumeric(df$surv_years)
    checkmate::assertNumeric(df$BW_Loess)
    checkmate::assertNumeric(df$norm_abundance)
    checkmate::assertLogical(df$is_ddm)
  }

  print("All columns are present and correctly factored")
  return(invisible())
}


#' Extract summary statistics from DOCR mediation test
#'
#' @description
#' This is function runs an \code{lm()} for all mediation tests
#' and extract summary statistics and AIC, tryCatch-wrapped for
#' error handling
#'
#' @param df a dataframe with columns corresponding to necessary variables
#' @param form a n \code{lm()} formula
#' @param mt the model type as string, for identification and sorting
#'
#' @returns a dataframe
#'
#' @export
docr_lm_med_test <- function(df,
                             form,
                             mt) {
  return_data <- tryCatch(
    expr = {
      model1 <- lm(as.formula(form), data = df)
      model1_data <- as.data.frame(summary(model1)$coef)
      colnames(model1_data) <- c("coef", "se_coef", "t-stat", "p_value_coef")
      model1_data$model_term <- row.names(model1_data)
      model1_data <- model1_data %>%
        dplyr::mutate(
          model = form,
          model_type = mt,
          df = df.residual(model1),
          AIC = AIC(model1)
        )

      return(model1_data)
    }, error = function(e) {
      return(data.frame())
    }
  )

  return(return_data)
}


#' Systemically run all DOCR mediation tests
#'
#' @description
#' This is function wraps \code{docr_lm_med_test()} for all mediation tests,
#' organizes dataframe, and computes mediation test
#'
#' @param data_use a dataframe with columns corresponding to necessary variables
#' @param outcome_var column name for outcome variable of interest, as as string
#' @param intervention_var column name for intervention variable of interest, as as string
#' @param mediation_var column name for mediation variable of interest, as as string
#' @param co_vars column name(s) for covariate variables to adjust for, as a string
#' @param n_samples number of samples in experiment
#' @param n_sample_percent minimum percent of \code{n_samples} that are non-missing
#' that are required to run model
#'
#' @returns a dataframe of all outcome results
#'
#' @export
docr_lm_sobel_mediation_test <- function(data_use,
                                         outcome_var,
                                         intervention_var,
                                         mediation_var,
                                         co_vars,
                                         n_samples,
                                         n_sample_percent = 0.1) {
  # Input validation
  if (!is.data.frame(data_use)) {
    stop("data_use must be a data.frame")
  }

  if (!all(c(outcome_var, intervention_var, mediation_var, co_vars) %in% names(data_use))) {
    stop("All specified variables must exist in data_use")
  }

  # Clean data for mediation variable
  data_use <- data_use %>%
    tidyr::drop_na(rlang::sym(mediation_var))

  # Early return if insufficient data for mediation variable
  min_sample_size <- n_samples * n_sample_percent
  if (nrow(data_use) == 0 || nrow(data_use) < min_sample_size) {
    return(data.frame())
  }

  all_outcome_results <- data.frame()

  for (ov in outcome_var) {
    data_current <- data_use %>%
      tidyr::drop_na(rlang::sym(ov))

    # Skip if insufficient data for this outcome variable
    if (nrow(data_current) == 0 || nrow(data_current) < min_sample_size) {
      next
    }

    # Filter covariates to only include those with variation
    co_vars_use <- co_vars[sapply(co_vars, function(cv) {
      length(unique(data_current[[cv]])) > 1
    })]

    # Create covariate string
    cov_string <- if (length(co_vars_use) > 0) {
      paste(" +", paste(co_vars_use, collapse = " + "))
    } else {
      ""
    }

    # Model formulas
    f2 <- paste0(mediation_var, " ~ ", intervention_var, cov_string)
    f3 <- paste0(ov, " ~ ", intervention_var, cov_string)
    f4 <- paste0(ov, " ~ ", intervention_var, " + ", mediation_var, cov_string)

    # Fit models with error handling
    models <- list(
      m2 = tryCatch(docr_lm_med_test(df = data_current, form = f2, mt = "Intervention>Mediator"),
        error = function(e) data.frame()
      ),
      m3 = tryCatch(docr_lm_med_test(df = data_current, form = f3, mt = "Intervention>Outcome"),
        error = function(e) data.frame()
      ),
      m4 = tryCatch(docr_lm_med_test(df = data_current, form = f4, mt = "Intervention+Mediator>Outcome"),
        error = function(e) data.frame()
      )
    )

    # Check if any models failed
    if (any(sapply(models, function(x) nrow(x) == 0))) {
      next
    }

    # Determine intervention variable terms
    intervention_vars <- if (is.character(data_current[[intervention_var]]) || is.factor(data_current[[intervention_var]])) {
      paste0(intervention_var, levels(as.factor(data_current[[intervention_var]])))
    } else {
      intervention_var
    }

    # Extract coefficients for Sobel test
    sobel_results <- docr_extract_sobel_med_effects(
      models,
      intervention_vars,
      mediation_var
    )

    # Add metadata to results
    if (nrow(sobel_results) > 0) {
      sobel_results$outcome_var <- ov
      sobel_results$mediation_var <- mediation_var
      sobel_results$intervention_var <- intervention_var
      sobel_results$n_obs <- nrow(data_current)

      all_outcome_results <- dplyr::bind_rows(
        all_outcome_results,
        sobel_results
      )
    }

    # Explicit garbage collection
    rm(data_current, models)
    invisible(gc())
  }

  return(all_outcome_results)
}

# Helper function to extract all mediation effects
docr_extract_sobel_med_effects <- function(models,
                                           intervention_vars,
                                           mediation_var) {
  m2 <- models$m2 # Intervention -> Mediator
  m3 <- models$m3 # Intervention -> Outcome (total effect)
  m4 <- models$m4 # Intervention + Mediator -> Outcome

  # Extract path a coefficients (intervention -> mediator)
  path_a <- m2 %>%
    dplyr::filter(model_term %in% intervention_vars) %>%
    dplyr::select(model_term, coef, se_coef, p_value_coef) %>%
    dplyr::rename(a_coeff = coef, a_se = se_coef, a_p = p_value_coef)

  # Extract path b coefficient (mediator -> outcome, controlling for intervention)
  path_b <- m4 %>%
    dplyr::filter(model_term == mediation_var) %>%
    dplyr::select(coef, se_coef, p_value_coef) %>%
    dplyr::rename(b_coeff = coef, b_se = se_coef, b_p = p_value_coef)

  # Extract total effect (c path: intervention -> outcome without mediator)
  total_effect <- m3 %>%
    dplyr::filter(model_term %in% intervention_vars) %>%
    dplyr::select(model_term, coef, se_coef, p_value_coef) %>%
    dplyr::rename(total_effect = coef, total_se = se_coef, total_p = p_value_coef)

  # Extract direct effect (c' path: intervention -> outcome controlling for mediator)
  direct_effect <- m4 %>%
    dplyr::filter(model_term %in% intervention_vars) %>%
    dplyr::select(model_term, coef, se_coef, p_value_coef) %>%
    dplyr::rename(direct_effect = coef, direct_se = se_coef, direct_p = p_value_coef)

  # Check if we have the required coefficients
  if (nrow(path_a) == 0 || nrow(path_b) == 0 || nrow(total_effect) == 0 || nrow(direct_effect) == 0) {
    return(data.frame())
  }

  # Merge all effects by intervention term
  all_effects <- path_a %>%
    dplyr::left_join(total_effect, by = "model_term") %>%
    dplyr::left_join(direct_effect, by = "model_term") %>%
    dplyr::mutate(
      # Add path b coefficients (same for all intervention terms)
      b_coeff = path_b$b_coeff[1],
      b_se = path_b$b_se[1],
      b_p = path_b$b_p[1],
      # Calculate indirect effect and Sobel test
      indirect_effect = a_coeff * b_coeff,
      sobel_se = sqrt(a_se^2 * b_coeff^2 + b_se^2 * a_coeff^2),
      sobel_Z_stat = indirect_effect / sobel_se,
      sobel_p = 2 * pnorm(-abs(sobel_Z_stat)),
      # Calculate proportion mediated (indirect/total)
      prop_mediated = ifelse(total_effect != 0, indirect_effect / total_effect, NA)
    ) %>%
    dplyr::select(
      model_term,
      # Path coefficients
      a_coeff, a_se, a_p, # Intervention -> Mediator
      b_coeff, b_se, b_p, # Mediator -> Outcome
      # Effects
      total_effect, total_se, total_p, # Total effect (c)
      direct_effect, direct_se, direct_p, # Direct effect (c')
      indirect_effect, sobel_se, sobel_p, # Indirect effect (a*b)
      sobel_Z_stat, # Sobel Z statistic
      prop_mediated
    ) # Proportion mediated

  return(all_effects)
}


# Extract summary stats and predictions from GAM model
docr_gam_predict <- function(gam_model,
                             smooth_term,
                             new_data) {
  # Add a single real RE to the model for correct error estimation
  # Also minimizes warnings
  temp_re <- as.character(unique(gam_model$model$mouse_id))[1]
  new_data <- new_data %>%
    dplyr::mutate(mouse_id = as.factor(temp_re))

  # Error from RE is included, but estimate is zeroed via exclude
  predictions <- mgcv::predict.gam(
    gam_model,
    newdata = new_data,
    se.fit = TRUE,
    exclude = "s(mouse_id)"
  )

  first_derivative <- gratia::derivatives(
    gam_model,
    select = smooth_term,
    type = "central",
    order = 1L,
    data = new_data
  ) %>%
    dplyr::select(
      deriv.1st = .derivative,
      se.deriv.1st = .se,
      upper.95.ci.derive.1st = .upper_ci,
      lower.95.ci.derive.1st = .lower_ci,
      PLL
    )

  second_derivative <- gratia::derivatives(
    gam_model,
    select = smooth_term,
    type = "central",
    order = 2L,
    data = new_data
  ) %>%
    dplyr::select(
      deriv.2nd = .derivative,
      se.deriv.2nd = .se,
      upper.95.ci.derive.2nd = .upper_ci,
      lower.95.ci.derive.2nd = .lower_ci,
      PLL
    ) %>%
    dplyr::mutate(
      inflection.pt = !is.na(lag(deriv.2nd)) & abs(sign(deriv.2nd) - sign(lag(deriv.2nd))) == 2
    )

  # Calculate confidence intervals and errors
  # Calculating 95% with 1.96 multiplier
  new_data_return <- new_data %>%
    dplyr::select(-mouse_id) %>%
    dplyr::mutate(
      smooth_term = smooth_term,
      fit = predictions$fit,
      se.fit = predictions$se.fit,
      upper.95.ci = fit + (1.96 * se.fit),
      lower.95.ci = fit - (1.96 * se.fit)
    ) %>%
    dplyr::full_join(first_derivative,
      by = "PLL",
      relationship = "many-to-many"
    ) %>%
    dplyr::full_join(second_derivative,
      by = "PLL",
      relationship = "many-to-many"
    )

  # Residuals
  resids <- residuals(gam_model)

  # Summary stats
  sts <- data.frame(
    pvalue = summary(gam_model)$s.table[smooth_term, "p-value"],
    adj.rsquared = summary(gam_model)$r.sq,
    deviation.exp = summary(gam_model)$dev.expl,
    edf = summary(gam_model)$s.table[smooth_term, "edf"],
    AIC = gam_model$aic,
    smooth_term = smooth_term
  )

  return(list(
    new_data = new_data_return,
    residuals = resids,
    summary_stats = sts
  ))
}


docr_gam_process <- function(mt,
                             new_data) {
  # Initialize data
  all_summary_data <- tibble()
  all_prediction_data <- tibble()
  all_residual_data <- list()

  rds <- tryCatch(readRDS(mt), error = function(e) NULL)
  if (is.null(rds)) {
    return(NULL)
  }

  # Extract information from model file
  trait_id <- names(rds)

  # Extract information from each model file
  for (gm in names(rds[[1]])) {
    gam_model <- rds[[1]][[gm]]

    # Extract information from each smooth_grep smooth term
    terms_use <- grep("s\\(PLL\\)",
      rownames(summary(gam_model)$s.table),
      value = TRUE
    )

    for (tu in terms_use) {
      # If on diet_group before diet to avoid errors
      if (grepl("diet_group", tu)) {
        new_data_temp <- data.frame(
          diet = strsplit(gsub("s\\(PLL\\):diet_group", "", tu), "-")[[1]]
        ) %>%
          tidyr::expand_grid(new_data) %>%
          dplyr::mutate(diet_group = gsub(
            "s\\(PLL\\):diet_group",
            "", tu
          ))
      } else if (grepl("diet", tu)) {
        new_data_temp <- new_data %>%
          dplyr::mutate(diet = gsub("s\\(PLL\\):diet", "", tu))
      } else {
        new_data_temp <- new_data %>%
          dplyr::mutate(diet = "AL")
      }

      # Extract prediction data
      res1 <- docr_gam_predict(
        gam_model = gam_model,
        smooth_term = tu,
        new_data = new_data_temp
      )

      # Combine summary data
      all_summary_data <- all_summary_data %>%
        dplyr::bind_rows(res1$summary_stats) %>%
        dplyr::mutate(trait = trait_id)

      # Combine prediction data
      all_prediction_data <- all_prediction_data %>%
        dplyr::bind_rows(res1$new_data) %>%
        dplyr::mutate(trait = trait_id)

      # Create residual data list
      all_residual_data[[tu]] <- res1$residuals
    } # End term for-loop
  } # End model for-loop

  # Return results as a list
  # Assumes models will not share coefficients!

  return(list(
    summary = all_summary_data,
    prediction = all_prediction_data,
    residual = setNames(list(all_residual_data), trait_id)
  ))
}


docr_elastic_train_test_parse <- function(final_data,
                                          train_frac = 0.8) {
  checkmate::assertSubset(c(
    "mouse_id", "fasting", "bw_test", "diet",
    "age_days", "surv_days"
  ), colnames(final_data))

  # make metabolite matrix
  X_metabo_matrix_initial <- final_data %>%
    dplyr::select(-fasting, -bw_test, -diet, -age_days, -surv_days) %>%
    tibble::column_to_rownames("mouse_id") %>%
    dplyr::select(where(~ mean(is.na(.)) <= 0.20)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), scale)) %>%
    glmnet::makeX(na.impute = TRUE)

  # clear covariates
  X_metabo_matrix_clean <- limma::removeBatchEffect(
    x = t(X_metabo_matrix_initial),
    batch = final_data %>% dplyr::pull(fasting),
    batch2 = final_data %>% dplyr::pull(diet),
    covariates = final_data %>% dplyr::pull(bw_test) %>% scale()
  ) %>%
    t()

  # make test and train mouse_id sets
  set.seed(212)
  mouse_id_train <- final_data %>%
    dplyr::group_by(diet, fasting) %>%
    dplyr::mutate(
      surv_quartile = ntile(surv_days, 4),
      mouse_id = as.character(mouse_id)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(fasting, diet, surv_quartile) %>%
    dplyr::slice_sample(prop = train_frac) %>%
    dplyr::ungroup() %>%
    dplyr::pull(mouse_id)
  mouse_id_test <- final_data %>%
    dplyr::mutate(mouse_id = as.character(mouse_id)) %>%
    dplyr::filter(!mouse_id %in% mouse_id_train) %>%
    dplyr::pull(mouse_id)

  # make control data frames
  X_control_df_train <- final_data %>%
    dplyr::select(mouse_id, fasting, bw_test, diet, age_days, surv_days) %>%
    dplyr::filter(mouse_id %in% mouse_id_train)
  X_control_df_test <- final_data %>%
    dplyr::select(mouse_id, fasting, bw_test, diet, age_days, surv_days) %>%
    dplyr::filter(mouse_id %in% mouse_id_test)

  # order rows
  train_idx <- match(X_control_df_train$mouse_id, rownames(X_metabo_matrix_clean))
  test_idx <- match(X_control_df_test$mouse_id, rownames(X_metabo_matrix_clean))
  X_metabo_matrix_train <- X_metabo_matrix_clean[train_idx, ]
  X_metabo_matrix_test <- X_metabo_matrix_clean[test_idx, ]

  # make survival objects
  surv_obj_train <- survival::Surv(
    time = X_control_df_train$age_days,
    time2 = X_control_df_train$surv_days,
    event = rep(1, nrow(X_control_df_train))
  )
  surv_obj_test <- survival::Surv(
    time = X_control_df_test$age_days,
    time2 = X_control_df_test$surv_days,
    event = rep(1, nrow(X_control_df_test))
  )

  return(list(
    X_control_df_train = X_control_df_train,
    X_control_df_test = X_control_df_test,
    X_metabo_matrix_train = X_metabo_matrix_train,
    X_metabo_matrix_test = X_metabo_matrix_test,
    surv_obj_train = surv_obj_train,
    surv_obj_test = surv_obj_test
  ))
}


docr_elastic_train_test_parse_diet <- function(final_data,
                                               train_frac = 0.8) {
  checkmate::assertSubset(c(
    "mouse_id", "fasting", "bw_test", "dietAL",
    "diet1D", "diet2D", "diet20", "diet40",
    "age_days", "surv_days"
  ), colnames(final_data))

  # make metabolite matrix
  X_metabo_matrix_initial <- final_data %>%
    dplyr::select(
      -dietAL, -diet1D, -diet2D, -diet20, -diet40,
      -fasting, -bw_test, -age_days, -surv_days
    ) %>%
    tibble::column_to_rownames("mouse_id") %>%
    dplyr::select(where(~ mean(is.na(.)) <= 0.20)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), scale)) %>%
    glmnet::makeX(na.impute = TRUE)

  # clear covariates
  X_metabo_matrix_clean <- limma::removeBatchEffect(
    x = t(X_metabo_matrix_initial),
    batch = final_data %>% dplyr::pull(fasting),
    covariates = final_data %>% dplyr::pull(bw_test) %>% scale()
  ) %>%
    t() %>%
    dplyr::bind_cols(final_data %>%
      dplyr::select(dietAL, diet1D, diet2D, diet20, diet40), .)
  rownames(X_metabo_matrix_clean) <- rownames(X_metabo_matrix_initial)

  # make test and train mouse_id sets
  set.seed(212)
  mouse_id_train <- final_data %>%
    tidyr::pivot_longer(
      cols = starts_with("diet", ignore.case = FALSE),
      names_to = "diet",
      values_to = "Is_Present"
    ) %>%
    dplyr::filter(Is_Present == 1) %>%
    dplyr::group_by(diet, fasting) %>%
    dplyr::mutate(
      surv_quartile = ntile(surv_days, 4),
      mouse_id = as.character(mouse_id)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(fasting, diet, surv_quartile) %>%
    dplyr::slice_sample(prop = train_frac) %>%
    dplyr::ungroup() %>%
    dplyr::pull(mouse_id)
  mouse_id_test <- final_data %>%
    dplyr::mutate(mouse_id = as.character(mouse_id)) %>%
    dplyr::filter(!mouse_id %in% mouse_id_train) %>%
    dplyr::pull(mouse_id)

  # make control data frames
  X_control_df_train <- final_data %>%
    dplyr::select(mouse_id, fasting, bw_test, age_days, surv_days) %>%
    dplyr::filter(mouse_id %in% mouse_id_train)
  X_control_df_test <- final_data %>%
    dplyr::select(mouse_id, fasting, bw_test, age_days, surv_days) %>%
    dplyr::filter(mouse_id %in% mouse_id_test)

  # order rows
  train_idx <- match(X_control_df_train$mouse_id, rownames(X_metabo_matrix_clean))
  test_idx <- match(X_control_df_test$mouse_id, rownames(X_metabo_matrix_clean))
  X_metabo_matrix_train <- X_metabo_matrix_clean[train_idx, ]
  X_metabo_matrix_test <- X_metabo_matrix_clean[test_idx, ]

  # make survival objects
  surv_obj_train <- survival::Surv(
    time = X_control_df_train$age_days,
    time2 = X_control_df_train$surv_days,
    event = rep(1, nrow(X_control_df_train))
  )
  surv_obj_test <- survival::Surv(
    time = X_control_df_test$age_days,
    time2 = X_control_df_test$surv_days,
    event = rep(1, nrow(X_control_df_test))
  )

  return(list(
    X_control_df_train = X_control_df_train,
    X_control_df_test = X_control_df_test,
    X_metabo_matrix_train = X_metabo_matrix_train,
    X_metabo_matrix_test = X_metabo_matrix_test,
    surv_obj_train = surv_obj_train,
    surv_obj_test = surv_obj_test
  ))
}


# define elastic net model
docr_elastic_model_fit <- function(i,
                                   X_data,
                                   Y_data,
                                   alpha_value,
                                   penalty_vec,
                                   bootstrap = FALSE,
                                   bootstrap_lamdba = c(
                                     "lambda.min",
                                     "lambda.1se",
                                     "lambda.ci"
                                   ),
                                   maxit = 1e+06,
                                   nlambda = 100,
                                   type_measure = c("default", "C") # default is partial likelihood
) {
  if (length(type_measure) > 1) {
    type_measure <- type_measure[1]
  }

  if (bootstrap) {
    boot_idx <- sample(nrow(X_data), replace = TRUE)
    X_data <- X_data[boot_idx, ]
    Y_data <- Y_data[boot_idx, ]
  }

  cv_model <- cv.glmnet(
    x = X_data,
    y = Y_data,
    family = "cox",
    alpha = alpha_value,
    penalty.factor = penalty_vec,
    maxit = maxit,
    nlambda = nlambda,
    type.measure = type_measure
  )

  if (bootstrap) {
    if (length(bootstrap_lamdba) > 1) {
      bootstrap_lamdba <- bootstrap_lamdba[1]
    }
    if (bootstrap_lamdba == "lambda.ci") {
      all_predictions <- predict(cv_model,
        newx = X_data,
        s = cv_model$lambda,
        type = "response"
      )
      c_indices <- apply(all_predictions, 2, function(pred_vec) {
        glmnet::Cindex(pred = pred_vec, y = Y_data)
      })
      bootstrap_lamdba <- cv_model$lambda[which.max(c_indices)]
    }
    return(as.numeric(coef(cv_model, s = bootstrap_lamdba)))
  } else {
    return(cv_model)
  }
}
