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
      default = "./",
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

# === Cross-validation =========================================================

cat("\n\n=== SETUP =========================================================\n")

# ggplot2 theme
setMyTheme()

# Pass arguments from the command line
args <- input_arguments()

cat("\nRead in data.")

data_dir <- args$data_dir # "~/Documents/PhD/TAGMMDI/CV_output/tan2009r1/"
save_dir <- args$save_dir

# The plot filenames
performance_save_name <- paste0(save_dir, "performance.png")
model_fit_save_name <- paste0(save_dir, "model_fit.png")

fold_outputs <- list.files(data_dir, pattern = "*.rds", full.names = TRUE)
n_folds <- length(fold_outputs)

mdi_folds <- tagm_folds <- vector("list", n_folds)

for (ii in seq(n_folds)) {
  .f <- fold_outputs[ii]
  .d <- readRDS(.f)

  mdi_folds[[ii]] <- .mdi <- .d$mdi_fold
  tagm_folds[[ii]] <- .tagm <- .d$tagm_fold

  # Recover the model performance scores for MDI and TAGM
  .perf <- data.frame(
    "Model" = "MDI",
    "Brier score" = .mdi$quadloss,
    "Accuracy" = .mdi$prediction_scores$accuracy,
    "Macro F1" = .mdi$prediction_scores$macro_f1,
    "Weighted F1" = .mdi$prediction_scores$weighted_f1
  )

  if (ii == 1) {
    performance_df <- .perf
  } else {
    performance_df <- rbind(performance_df, .perf)
  }

  .perf <- data.frame(
    "Model" = "TAGM",
    "Brier score" = .tagm$quadloss,
    "Accuracy" = .tagm$prediction_scores$accuracy,
    "Macro F1" = .tagm$prediction_scores$macro_f1,
    "Weighted F1" = .tagm$prediction_scores$weighted_f1
  )

  performance_df <- rbind(performance_df, .perf)

  # Recover the model fit
  .mdi$likelihood_df$Model <- "MDI"
  .tagm$likelihood_df$Model <- "TAGM"

  .l_df <- rbind(.mdi$likelihood_df, .tagm$likelihood_df) %>%
    pivot_longer(-c(Iteration, Fold, Model), names_to = "Chain", values_to = "Complete_log_likelihood")

  if (ii == 1) {
    likelihood_df <- .l_df
  } else {
    likelihood_df <- rbind(likelihood_df, .l_df)
  }
}

cat("\n\n=== PLOTTING ======================================================\n")

cat("\nPlot model fit using the complete log-likelihood.")

p_likelihood <- likelihood_df %>%
  ggplot(aes(x = Iteration, y = Complete_log_likelihood, colour = Chain)) +
  geom_line() +
  facet_grid(Model ~ Fold, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  labs(y = "Complete log-likelihood")

cat("\nPlot model performance across 4 scores.")

p_performance <- performance_df %>%
  pivot_longer(-Model, names_to = "Score", values_to = "Value") %>%
  ggplot(aes(x = Model, y = Value)) +
  geom_boxplot() +
  facet_wrap(~Score, scales = "free_y")

ggsave(model_fit_save_name, plot = p_likelihood, height = 5, width = 10)
ggsave(performance_save_name, plot = p_performance, height = 5, width = 5)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")
