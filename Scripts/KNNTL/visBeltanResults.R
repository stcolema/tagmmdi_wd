library(ggplot2)
library(tidyr)
mdiHelpR::setMyTheme()

my_dir <- "~/Desktop/TestMDIcall/Beltran/"

model_runs_files <- list.files(my_dir, full.names = TRUE)
n_files <- length(model_runs_files)
n_models <- n_files * 2
L <- 5
my_df <- data.frame(
  "Model" = c(rep("MDI", L * n_files), rep("TAGM", L * n_files)),
  "Accuracy" = 0,
  "Macro_F1" = 0,
  "Weighted_F1" = 0,
  "Brier_loss" = 0,
  "Modality" = 0
)

entry <- 0
for (ii in seq(1, n_files)) {
  run <- readRDS(model_runs_files[ii])
  N_test <- length(run$test.idx)
  pred_scores <- run$mdi_fold$prediction_scores

  for (ll in seq(1, L)) {
    entry <- entry + 1
    my_df$Model[entry] <- "MDI"
    my_df$Accuracy[entry] <- pred_scores[[ll]]$accuracy
    my_df$Macro_F1[entry] <- pred_scores[[ll]]$macro_f1
    my_df$Weighted_F1[entry] <- pred_scores[[ll]]$weighted_f1
    my_df$Brier_loss[entry] <- pred_scores[[ll]]$quadloss
    my_df$Modality[entry] <- ll
  }
}



for (ii in seq(1, n_files)) {
  run <- readRDS(model_runs_files[ii])
  N_test <- length(run$test.idx)


  for (ll in seq(1, L)) {
    entry <- entry + 1
    my_df$Model[entry] <- "TAGM"
    pred_scores <- run$tagm_fold[[ll]]$prediction_scores
    my_df$Accuracy[entry] <- pred_scores[[1]]$accuracy
    my_df$Macro_F1[entry] <- pred_scores[[1]]$macro_f1
    my_df$Weighted_F1[entry] <- pred_scores[[1]]$weighted_f1
    my_df$Brier_loss[entry] <- pred_scores[[1]]$quadloss
    my_df$Modality[entry] <- ll
  }

  # my_df$Model[n_files + ii] <- "TAGM"
  # my_df$Accuracy[n_files + ii] <- pred_scores$accuracy
  # my_df$Macro_F1[n_files + ii] <- pred_scores$macro_f1
  # my_df$Weighted_F1[n_files + ii] <- pred_scores$weighted_f1
  # my_df$Brier_loss[n_files + ii] <- pred_scores$brier_loss
}


my_df$Modality <- factor(my_df$Modality, labels = c(
  "HCMV24",
  "HCMV48",
  "HCMV72",
  "HCMV96",
  "HCMV120"
))

my_df |>
  ggplot(aes(x = Model, y = Weighted_F1, fill = Model)) +
  geom_boxplot() +
  facet_wrap(~Modality) +
  ggthemes::scale_fill_colorblind()

my_df |>
  ggplot(aes(x = Model, y = Brier_loss, fill = Model)) +
  geom_boxplot() +
  facet_wrap(~Modality) +
  ggthemes::scale_fill_colorblind()

my_df |>
  pivot_longer(-c(Model, Modality), names_to = "Score", values_to = "Value") |>
  ggplot(aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot() +
  facet_grid(Score ~ Modality, scales = "free_y") +
  ggthemes::scale_fill_colorblind()

write.csv(my_df, "./KNN_TL_comp/Beltran_perf.csv")
