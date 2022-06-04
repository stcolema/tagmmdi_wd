library(tibble)
library(tidyr)
library(ggplot2)
library(mdiHelpR)
library(tagmReDraft)
library(magrittr)
library(ggthemes)
library(pRolocdata)
library(optparse)



setMyTheme()
sim_col_pal <- simColPal()

data_dir <- "./Results/CV_output/"
fracs <- c(10, 30, 50)
datasets <- c("E14TG2aS1", "groen2014r1", "dunkley2006", "HEK293T2011")
result_df <- NULL

result_csv_file <- paste0(data_dir, "result_df.csv")
result_plot_file <- paste0(data_dir, "result_plot.png")

for (dataset in datasets) {
  for (test_frac in fracs) {
    curr_dir <- paste0(data_dir, "test_", test_frac, "/", dataset)

    bayesian_pattern <- paste0("*_nChains_5_testSize_", test_frac, "_K_50_*")
    knn_tl_pattern <- paste0("*_knnTL_numberWeights_5_*")

    bayesian_files <- list.files(curr_dir, pattern = bayesian_pattern, full.names = TRUE)
    knn_tl_files <- list.files(curr_dir, pattern = knn_tl_pattern, full.names = TRUE)

    n_files <- length(bayesian_files)
    n_files_2 <- length(knn_tl_files)
    if (n_files != n_files_2) {
      stop("Different numbers of files between transfer learner and Bayesian methods.")
    }

    for (ii in seq(1, n_files)) {
      bayes_f <- bayesian_files[ii]
      knn_tl_f <- knn_tl_files[ii]

      .x <- readRDS(bayes_f)

      .mdi_df <- data.frame(
        "Model" = "MDI",
        "Accuracy" = .x$mdi_fold$prediction_scores$accuracy,
        "Macro_F1" = .x$mdi_fold$prediction_scores$macro_f1,
        "Weighted_F1" = .x$mdi_fold$prediction_scores$weighted_f1,
        "Brier_loss" = .x$mdi_fold$quadloss,
        "Seed" = .x$seed
      )

      .tagm_df <- data.frame(
        "Model" = "TAGM",
        "Accuracy" = .x$tagm_fold$prediction_scores$accuracy,
        "Macro_F1" = .x$tagm_fold$prediction_scores$macro_f1,
        "Weighted_F1" = .x$tagm_fold$prediction_scores$weighted_f1,
        "Brier_loss" = .x$tagm_fold$quadloss,
        "Seed" = .x$seed
      )

      .x <- readRDS(knn_tl_f)


      .knn_df <- data.frame(
        "Model" = "KNN_TL",
        "Accuracy" = .x$knn_fold$prediction_scores$accuracy,
        "Macro_F1" = .x$knn_fold$prediction_scores$macro_f1,
        "Weighted_F1" = .x$knn_fold$prediction_scores$weighted_f1,
        "Brier_loss" = .x$knn_fold$quadloss,
        "Seed" = .x$seed
      )

      .cur_df <- rbind(.mdi_df, .tagm_df, .knn_df)
      .cur_df$Dataset <- dataset
      .cur_df$Test_frac <- test_frac
      if (is.null(result_df)) {
        result_df <- .cur_df
      } else {
        result_df <- rbind(result_df, .cur_df)
      }
    }
  }
}

result_df$Dataset <- factor(result_df$Dataset, 
  labels = c("Callus", "Root", "Human", "Mouse"), 
  levels = c("dunkley2006", "groen2014r1", "HEK293T2011", "E14TG2aS1")
)

result_df$Model <- factor(result_df$Model, 
  labels = c("MDI", "TAGM", "TL"),
  levels = c("TAGM", "MDI", "KNN_TL")
)

long_result_df <- result_df %>%
  pivot_longer(-c(Model, Seed, Dataset, Test_frac),
    names_to = "Score",
    values_to = "Value"
  ) 

long_result_df$Score <- factor(long_result_df$Score,
  labels = c("Accuracy", "Macro F1", "Weighted F1", "Brier loss"),
  levels = c("Accuracy", "Macro_F1", "Weighted_F1", "Brier_loss")
)

long_result_df[["Test fraction"]] <- long_result_df$Test_frac * 0.01

# Use the colorblind palette from ggthemes, skipping black
my_palette <- ggthemes::colorblind_pal()(4)[2:4]

p_out <- long_result_df %>% 
  ggplot(aes(x = Model, y = Value)) +
  geom_boxplot(aes(fill = Model)) +
  facet_grid(Score ~ Dataset + `Test fraction`, scales = "free", labeller = label_both) +
  scale_fill_manual(values = my_palette)

write.csv(result_df, file = result_csv_file)
ggsave(filename = result_plot_file, plot = p_out, width = 16, height = 8)
