# Install and load necessary packages
suppressMessages(library(mgcv))
suppressMessages(library(tidyverse))
library(future)
library(furrr)

# set
# local_filepath <- "~/drido-multiomic-paper"
# gam_output_filepath <- "~/"

source(file.path(local_filepath, "R/statistics_functions.R"))
source(file.path(local_filepath, "R/figure_functions.R"))

name_conv_file <- "inst/supp_tables/Table_S1_CompoundAnnotations.csv"

### Process Models ### -------

# Make new data
# mouse_id added in docr_gam_predict() fxn
new_data <- data.frame(
  PLL = seq(0.3, 1, by = 0.01), # predict every 1% PLL
  fasting = "No"
)

# Initialize data
all_summary_data <- tibble()
all_prediction_data <- tibble()
all_residual_data <- list()

# Collect results from models
future::plan(future::multisession, workers = future::availableCores() - 1)
results <- furrr::future_map(
  .x = list.files(file.path(gam_output_filepath, "gam_models_PLL_5"), full.names = TRUE),
  ~ docr_gam_process(
    mt = .x,
    new_data = new_data
  ),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)
future::plan(future::sequential)

# Combine all results
all_summary_data <- purrr::map_dfr(results, ~ .x$summary)
all_prediction_data <- purrr::map_dfr(results, ~ .x$prediction)

# Process residual data
for (res in results) {
  trait <- names(res$residual)[1]
  all_residual_data[[trait]] <- res$residual[[trait]]
}

# Save the results
# saveRDS(
#   list(all_summary_data, all_prediction_data, all_residual_data),
#   file = file.path(gam_output_filepath, "20250310-gam-models-parsed-PLL-k-1.Rds")
# )

### Cluster ### -----

# pll_gam_fit_temp <- readRDS(file.path(gam_output_filepath, "20250310-gam-models-parsed-PLL-k-1-all.Rds"))
# all_summary_data <- pll_gam_fit_temp[[1]]
# all_prediction_data <- pll_gam_fit_temp[[2]]

## Extract good compound IDs
name_conversion_use <- read.csv(file.path(local_filepath, name_conv_file)) %>%
  dplyr::select(modality, feature_id, name_use)

# Clean summary df
all_summary_data <- all_summary_data %>%
  dplyr::filter(smooth_term == "s(PLL)") %>%
  dplyr::group_by(smooth_term) %>%
  dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  dplyr::ungroup()

# Clean prediction df
all_prediction_data <- all_prediction_data %>%
  dplyr::filter(smooth_term == "s(PLL)") %>%
  dplyr::left_join(name_conversion_use, by = c("trait" = "feature_id")) %>%
  dplyr::filter(!(is.na(name_use) & !is.na(modality)))

# identify significant traits
fdr_cutoff <- 0.05

# use maximum R2 of metabolomics internal standard as minimum R2 threshold
r2_cutoff <- all_summary_data %>%
  dplyr::filter(
    padj < fdr_cutoff,
    grepl("13C|[Dd]\\d\\.", trait)
  ) %>%
  {
    max(.$adj.rsquared, na.rm = TRUE)
  }

# filter summary data to good traits; keep physiological traits
all_summary_data <- all_summary_data %>%
  dplyr::left_join(name_conversion_use, by = c("trait" = "feature_id")) %>%
  dplyr::filter(!(is.na(name_use) & !is.na(modality)))

# filter summary by R2 and adj.pvalue
significant_all_traits <- all_summary_data %>%
  dplyr::filter(
    padj < fdr_cutoff,
    adj.rsquared > r2_cutoff
  ) %>%
  dplyr::pull(trait)

# add significance column and save
all_summary_data <- all_summary_data %>%
  dplyr::mutate(significant = trait %in% significant_all_traits)
all_prediction_data <- all_prediction_data %>%
  dplyr::mutate(significant = trait %in% significant_all_traits)

# To find max on single slope
top_PLL_inflection <- all_prediction_data %>%
  dplyr::filter(significant) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(is_max_value = ifelse(abs(deriv.2nd) == max(abs(deriv.2nd)), TRUE, FALSE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is_max_value) %>%
  dplyr::select(trait, extrema_trait = PLL)

predict_matrix_fit <- all_prediction_data %>%
  dplyr::filter(significant) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(fit.scale = scale(fit, scale = TRUE, center = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::select(PLL, trait, fit.scale) %>%
  tidyr::spread(PLL, fit.scale) %>%
  tibble::column_to_rownames("trait")

# Run prcomp for PCA and extract matrix
p <- prcomp(predict_matrix_fit)
pca_df <- p$x %>%
  as.data.frame()

# Extract rotation data
pca_rot <- p$rotation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PLL") %>%
  dplyr::select(1:7) %>%
  tidyr::gather(2:7, key = "PC", value = "loading")

# Create annotation dataframe
pca_var <- summary(p)$importance[2, ]
annotation_df <- data.frame(
  PC = names(pca_var[1:6]),
  annotation_value = paste0(round(pca_var[1:6], 3) * 100, "% of var")
) %>%
  dplyr::mutate(
    PLL = -Inf,
    loading = Inf
  )

# plot for QC
qc_plots <- FALSE
if (qc_plots) {
  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point() +
    theme_classic() +
    ggtitle("PCA Projections")

  ggplot(pca_rot) +
    geom_point(aes(
      x = PLL,
      y = loading
    )) +
    facet_wrap(~PC) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    geom_text(
      data = annotation_df,
      aes(x = PLL, y = loading, label = annotation_value),
      hjust = -0.5, vjust = 1.5, size = 4, color = "red"
    ) +
    ggtitle("PCA Loadings")
}

# Only use PCs that describe any variance
important_pcs <- pca_var[pca_var > 0]
important_pc_df <- pca_df %>%
  dplyr::select(1:length(important_pcs))

cat(paste0("PC1 to PC", length(important_pcs), " account for 100% of the variability in the dataset"))
cat("\nThe corresponding values (percent variability described) are used as weights in k-means clustering of principle components:\n\n")
cat(important_pcs)

# Weighted Kmeans clustering
kmeans_flexclust_wrapper <- function(x, k) {
  kcca_result <- suppressMessages(
    flexclust::kcca(x, k,
      family = kccaFamily("kmeans"),
      weights = important_pcs
    )
  )
  result <- list(cluster = flexclust::clusters(kcca_result))
  return(result)
}


# Run gap stat to get cluster number
set.seed(2017)
gap_stat <- suppressMessages(
  cluster::clusGap(important_pc_df,
    FUN = kmeans_flexclust_wrapper,
    K.max = 20,
    B = 50
  )
)

gap_stat_fig <- factoextra::fviz_gap_stat(gap_stat)
print(gap_stat_fig)
opt_k <- gap_stat_fig$layers[[4]]$data$xintercept

# Run weighted KMeans clustering with PCA weights
set.seed(2017)
clusters_row <- kmeans_flexclust_wrapper(
  x = important_pc_df,
  k = opt_k
)$cluster %>%
  as.data.frame() %>%
  dplyr::rename(cluster = ".") %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>%
  tibble::rownames_to_column("trait") %>%
  dplyr::left_join(top_PLL_inflection, by = "trait") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(
    extrema_mean = mean(extrema_trait, na.rm = TRUE),
    extrema_se = sd(extrema_trait) / sqrt(dplyr::n())
  ) %>%
  dplyr::ungroup()

# change numbers for visualization
clusters_change <- clusters_row %>%
  dplyr::mutate(
    cluster = as.numeric(cluster),
    cluster = dplyr::case_when(cluster == 1 ~ 1,
      cluster == 2 ~ 7,
      cluster == 3 ~ 10,
      cluster == 4 ~ 9,
      cluster == 5 ~ 3,
      cluster == 6 ~ 8,
      cluster == 7 ~ 4,
      cluster == 9 ~ 6,
      cluster == 8 ~ 5,
      cluster == 10 ~ 2,
      .default = cluster
    )
  )

# add cluster information to prediction and summary df
all_summary_data <- all_summary_data %>%
  dplyr::left_join(clusters_change, by = "trait") %>%
  dplyr::mutate(
    modality = ifelse(is.na(modality), "physiological", modality),
    cluster = as.factor(cluster)
  ) %>%
  dplyr::mutate(name_use = ifelse(modality == "physiological", trait, name_use))
all_prediction_data <- all_prediction_data %>%
  dplyr::left_join(clusters_change, by = "trait") %>%
  dplyr::mutate(
    modality = ifelse(is.na(modality), "physiological", modality),
    cluster = as.factor(cluster)
  ) %>%
  dplyr::mutate(name_use = ifelse(modality == "physiological", trait, name_use))


## Save
# saveRDS(all_summary_data, file.path(local_filepath, "20251008-PLL-GAM-Summary-Data.Rds"))
# saveRDS(all_prediction_data, file.path(local_filepath, "20251008-PLL-GAM-Fit-Data.Rds"))
