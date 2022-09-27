# cross_val_visualisation.R
# Visualisation of the cross-validation results for the 5 datasets from Breckels
# et al.
#
# Output: two-layered boxplot of the accuracy and quadratic loss for MDI vs TAGM
#

library(magrittr)
library(tibble)
library(ggplot2)
library(patchwork)

mdiHelpR::setMyTheme()

# Where the output objects live
data_dir <- "../TAGMMDI/Validation/"

# The data files
data_files <- list.files(path = data_dir, pattern = "*75.rds", full.names = T)

# The names of the datasets
short_names <- c("Dunkley", "E14TG2aS1", "Groen", "HEK293T2011", "Tan")

# Models used
models <- c("MDI", "TAGM")

# Where to save the output to
save_dir <- data_dir
plot_name <- paste0(save_dir, "model_performance.png")

cv_data <- lapply(data_files, readRDS) %>%
  set_names(short_names)

# These are iterated over
n_datasets <- length(data_files)
n_models <- length(cv_data[[1]])
n_chains <- length(cv_data[[1]][[1]]$cmlist)

# Break the results of interest out of the lists and into a single data frame
for (ii in seq(1, n_datasets)) {
  .d <- short_names[ii]
  .d_results <- cv_data[[ii]]

  for (jj in seq(1, n_models)) {
    .mod <- models[jj]
    .curr <- .d_results[[jj]]
    .loss <- unlist(.curr$quadloss)

    for (kk in seq(1, n_chains)) {
      .tab <- .curr$cmlist[[kk]]

      .tp <- sum(diag(.tab))
      .n <- sum(.tab)
      .acc <- .tp / .n

      if (kk == 1) {
        .accuracy <- .acc
      } else {
        .accuracy <- c(.accuracy, .acc)
      }
    }

    .df <- data.frame(
      "Dataset" = .d,
      "Model" = .mod,
      "Quad_loss" = .loss,
      "Accuracy" = .accuracy
    )

    if (jj == 1 & ii == 1) {
      results_df <- .df
    } else {
      results_df <- rbind(results_df, .df)
    }
  }
}

# The model accuracy
p_acc <- results_df %>%
  ggplot(aes(x = Model, y = Accuracy)) +
  geom_boxplot(fill = "gold") +
  facet_grid(~Dataset, nrow = 1) +
  ylim(0, 1)

# The quadratic loss
p_quad <- results_df %>%
  ggplot(aes(x = Model, y = Quad_loss)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Dataset, nrow = 1) +
  labs(y = "Quadratic loss")

p_patch <- p_acc / p_quad
p_saved <- p_patch + plot_annotation(tag_levels = "A")

ggsave(plot_name, plot = p_saved, width = 9, height = 6)
