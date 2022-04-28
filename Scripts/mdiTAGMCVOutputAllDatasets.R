#!/usr/bin/env Rscript
# 
# Title: MDI & TAGM Cross validation outputs
# Description: This script takes the .RDS object outputted by 
# ``mdiTagmCVSingleLoop.R`` and plots the complete log-likelihood for each chain
# of MDI and the mixture model, and the performance of each method under 
# different performance scores.
# Output: Model fit plot, model performance plot.
# 
# Example call: 
# Rscript mdiTAGMCVOutputs.R --data_dir "./tan2009r1" --save_dir "./"

suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(optparse))

# === Functions ================================================================

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    
    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--data_dir"),
      type = "character",
      default = "./CV_output/",
      help = "Directory to read input from [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
    )
  )
  
  
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

readFile <- function(filename) {
  out <- tryCatch(
    {
      readRDS(filename)
    },
    error=function(cond) {
      message(paste("File does not seem to exist:", filename))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("Filename caused a warning:", filename))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )  
  
  out
}


# === Cross-validation =========================================================

cat("\n\n=== SETUP =========================================================\n")

# ggplot2 theme
setMyTheme()

# Pass arguments from the command line
args <- input_arguments()

cat("\nRead in data.")

data_dir <- args$data_dir
save_dir <- args$save_dir

# The plot filenames
performance_save_name <- paste0(save_dir, "performance.png")
model_fit_save_name <- paste0(save_dir, "complete_likelihood.png")
phi_plot_save_name <- paste0(save_dir, "phis.png")
accuracy_save_name <- paste0(save_dir, "accuracy.png")

datasets <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
n_datasets <- length(datasets)

cat("\nDatasets found:", datasets, "\n")

map <- list("Callus", "Mouse", "Root", "Human", "Fly")
names(map) <- datasets

for(ii in seq(1, n_datasets)) {
  .d <- datasets[ii]
  dataset_name <- map[[.d]]
  
  .curr_dir <- paste0(data_dir, "/", .d)
  filenames <- list.files(.curr_dir,
                          pattern = "*_nChains_5_testSize_", 
                          full.names = TRUE)
  knntl_filenames <- list.files(.curr_dir, 
                                pattern = "*knnTL_numberWeights_5_seed_",
                                full.names = TRUE
                                )
  n_folds <- length(filenames)
  
  knn_folds <- mdi_folds <- tagm_folds <- vector("list", n_folds)
  
  for (jj in seq(n_folds)) {
    .f <- filenames[jj]
    .df <- readRDS(.f)
    
    .knn_f <- knntl_filenames[jj]
    .knn_df <- readFile(.knn_f)
    transfer_learning_output_missing <- is.na(.knn_df)
    
    mdi_folds[[jj]] <- .mdi <- .df$mdi_fold
    tagm_folds[[jj]] <- .tagm <- .df$tagm_fold
    
    
    # Recover the model performance scores for MDI and TAGM
    .perf <- data.frame(
      "Model" = "MDI",
      "Dataset" = dataset_name,
      "Brier score" = .mdi$quadloss,
      "Accuracy" = .mdi$prediction_scores$accuracy,
      "Macro F1" = .mdi$prediction_scores$macro_f1,
      "Weighted F1" = .mdi$prediction_scores$weighted_f1
    )
    
    if ((jj == 1) & (ii == 1)) {
      performance_df <- .perf
    } else {
      performance_df <- rbind(performance_df, .perf)
    }
    
    .perf <- data.frame(
      "Model" = "TAGM",
      "Dataset" = dataset_name,
      "Brier score" = .tagm$quadloss,
      "Accuracy" = .tagm$prediction_scores$accuracy,
      "Macro F1" = .tagm$prediction_scores$macro_f1,
      "Weighted F1" = .tagm$prediction_scores$weighted_f1
    )
    
    performance_df <- rbind(performance_df, .perf)
    
    if(! transfer_learning_output_missing) {
    knn_folds[[jj]] <- .knn <- .knn_df$knn_fold
    .perf <- data.frame(
      "Model" = "KNN_TL",
      "Dataset" = dataset_name,
      "Brier score" = .knn$quadloss,
      "Accuracy" = .knn$prediction_scores$accuracy,
      "Macro F1" = .knn$prediction_scores$macro_f1,
      "Weighted F1" = .knn$prediction_scores$weighted_f1
    )
    performance_df <- rbind(performance_df, .perf)
    
    }
    
    # Recover the model fit
    .mdi$likelihood_df$Model <- "MDI"
    .tagm$likelihood_df$Model <- "TAGM"
    
    .l_df <- rbind(.mdi$likelihood_df, .tagm$likelihood_df) %>%
      mutate(Dataset = dataset_name) %>% 
      pivot_longer(-c(Iteration, Fold, Model, Dataset), names_to = "Chain", values_to = "Complete_log_likelihood")
    
    .phi_df <- .mdi$phi_df  %>% 
      mutate(Fold = .mdi$likelihood_df$Fold[1], Dataset = dataset_name) %>% 
      pivot_longer(-c(Iteration, Fold, Dataset), names_to = "Chain", values_to = "Phi")
    
    if ((jj == 1) & (ii == 1)) {
      likelihood_df <- .l_df
      phi_df <- .phi_df
    } else {
      likelihood_df <- rbind(likelihood_df, .l_df)
      phi_df <- rbind(phi_df, .phi_df)
    }
  }
}


cat("\n\n=== PLOTTING ======================================================\n")

cat("\nPlot sampled distribution for phi parameter of MDI.")

p_phi <- phi_df %>% 
  ggplot(aes(x = Phi, fill = Chain)) +
  geom_density(alpha = 0.3) +
  facet_grid(Dataset ~ Fold, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  labs(y = "Complete log-likelihood")

cat("\nPlot model fit using the complete log-likelihood.")

p_likelihood <- likelihood_df %>%
  ggplot(aes(x = Iteration, y = Complete_log_likelihood, colour = Chain)) +
  geom_line() +
  facet_grid(Model + Dataset ~ Fold, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  labs(y = "Complete log-likelihood")

cat("\nPlot model performance across 4 scores.")

p_performance <- performance_df %>%
  pivot_longer(-c(Dataset, Model), names_to = "Score", values_to = "Value") %>%
  ggplot(aes(x = Model, y = Value)) +
  geom_boxplot() +
  facet_grid(Score~Dataset, scales = "free_y")

p_accuracy <- performance_df %>%
  select(Model, Dataset, Brier.score, Accuracy) %>% 
  pivot_longer(-c(Dataset, Model), names_to = "Score", values_to = "Value") %>%
  ggplot(aes(x = Model, y = Value)) +
  geom_boxplot() +
  facet_grid(Score~Dataset, scales = "free_y")

ggsave(phi_plot_save_name, plot = p_phi, height = 5, width = 10)
ggsave(model_fit_save_name, plot = p_likelihood, height = 5, width = 10)
ggsave(accuracy_save_name, plot = p_accuracy, height = 5, width = 5)
ggsave(performance_save_name, plot = p_performance, height = 5, width = 5)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")
