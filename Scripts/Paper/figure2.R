suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(ggthemes))
suppressMessages(library(optparse))

setMyTheme()
sim_col_pal <- simColPal()

col_pal <- c(
  "#a7bed3",
  "#adf7b6",
  "#c6e2e9",
  "#f1ffc4",
  "#ffcaaf",
  "#dab894",
  "#70d6ff",
  "#ff70a6",
  "#ff9770",
  "#ffd670",
  "#e9ff70",
  "#f08080",
  "#fbc4ab",
  "#dfb2f4",
  "#f5e960",
  "#f5e5f0",
  "#55d6c2",
  "#ffffff",
  "#84dcc6",
  "#a5ffd6",
  "#79addc",
  "#ffc09f",
  "#ffee93",
  "#fcf5c7",
  "#ffa69e",
  "#ff686b",
  "#808080"
)


data("HEK293T2011")
data("groen2014r1")
data("dunkley2006")
data("E14TG2aS1")


data_dir <- "./KNN_TL_comp/"
fracs <- c(30)
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
      N_test <- length(.x$test.idx)
      .mdi_df <- data.frame(
        "Model" = "MDI",
        "Accuracy" = .x$mdi_fold$prediction_scores$accuracy,
        "Macro_F1" = .x$mdi_fold$prediction_scores$macro_f1,
        "Weighted_F1" = .x$mdi_fold$prediction_scores$weighted_f1,
        "Brier_loss" = .x$mdi_fold$quadloss / N_test,
        "Seed" = .x$seed
      )

      .tagm_df <- data.frame(
        "Model" = "TAGM",
        "Accuracy" = .x$tagm_fold$prediction_scores$accuracy,
        "Macro_F1" = .x$tagm_fold$prediction_scores$macro_f1,
        "Weighted_F1" = .x$tagm_fold$prediction_scores$weighted_f1,
        "Brier_loss" = .x$tagm_fold$quadloss / N_test,
        "Seed" = .x$seed
      )

      .x <- readRDS(knn_tl_f)
      N_test <- length(.x$test.idx)

      .knn_df <- data.frame(
        "Model" = "KNN_TL",
        "Accuracy" = .x$knn_fold$prediction_scores$accuracy,
        "Macro_F1" = .x$knn_fold$prediction_scores$macro_f1,
        "Weighted_F1" = .x$knn_fold$prediction_scores$weighted_f1,
        "Brier_loss" = .x$knn_fold$quadloss / N_test,
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
  labels = c("MDI", "TAGM", "KNN TL"),
  levels = c("MDI", "TAGM", "KNN_TL")
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

p_performance <- long_result_df %>%
  dplyr::filter(Score != "Weighted F1") |>
  ggplot(aes(x = Model, y = Value)) +
  geom_boxplot(aes(fill = Model)) +
  facet_grid(Score ~ Dataset, scales = "free") +
  scale_fill_manual(values = my_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # +
# theme(legend.position = "bottom")

# write.csv(result_df, file = result_csv_file)
# ggsave(filename = result_plot_file, plot = p_out, width = 16, height = 7)

# === Human dataset ============================================================

human_marker_proteins_df <- pRoloc:::subsetAsDataFrame(HEK293T2011, "markers", train = TRUE)
human_lopit <- human_marker_proteins_df |>
  dplyr::select(-markers)

human_pcs <- prcomp(human_lopit, center = TRUE, scale. = TRUE)

human_df <- as.data.frame(human_pcs$x[, c(1, 2)]) |>
  mutate(Marker = human_marker_proteins_df$markers, Dataset = "Human")

human_col_pal <- c(
  "#a7bed3",
  "#70d6ff",
  "#c6e2e9",
  "#ff9770",
  "#dfb2f4",
  "#f5e960",
  "#ffffff",
  "#79addc",
  "#fcf5c7",
  "#ffa69e",
  "#ff686b",
  "#808080"
)

names(human_col_pal) <- levels(human_df$Marker)

p_human <- human_df |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Marker), color = "#808080", shape = 21) +
  facet_wrap(~Dataset) +
  scale_fill_manual(values = human_col_pal) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# === Root dataset =============================================================

root_marker_proteins_df <- pRoloc:::subsetAsDataFrame(groen2014r1, "markers", train = TRUE)
root_lopit <- root_marker_proteins_df |>
  dplyr::select(-markers)

root_pcs <- prcomp(root_lopit, center = TRUE, scale. = TRUE)

root_df <- as.data.frame(root_pcs$x[, c(1, 2)]) |>
  mutate(Marker = root_marker_proteins_df$markers, Dataset = "Root")


root_col_pal <- c(
  "#c6e2e9",
  "#f5e960",
  "#ffffff",
  "#79addc",
  "#ff686b",
  "#808080"
)

names(root_col_pal) <- levels(root_df$Marker)

p_root <- root_df |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Marker), color = "#808080", shape = 21) +
  facet_wrap(~Dataset) +
  scale_fill_manual(values = root_col_pal) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# === Callus dataset ===========================================================

callus_marker_proteins_df <- pRoloc:::subsetAsDataFrame(dunkley2006, "markers", train = TRUE)
callus_lopit <- callus_marker_proteins_df |>
  dplyr::select(-markers)

callus_pcs <- prcomp(callus_lopit, center = TRUE, scale. = TRUE)

callus_df <- as.data.frame(callus_pcs$x[, c(1, 2)]) |>
  mutate(Marker = callus_marker_proteins_df$markers, Dataset = "Callus")


callus_col_pal <- c(
  "#70d6ff",
  "#c6e2e9",
  "#ff9770",
  "#f5e960",
  "#ffffff",
  "#79addc",
  "#fcf5c7",
  "#ff686b",
  "#808080"
)

names(callus_col_pal) <- levels(callus_df$Marker)

p_callus <- callus_df |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Marker), color = "#808080", shape = 21) +
  facet_wrap(~Dataset) +
  scale_fill_manual(values = callus_col_pal) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# === Mouse dataset ============================================================

mouse_marker_proteins_df <- pRoloc:::subsetAsDataFrame(E14TG2aS1, "markers", train = TRUE)
mouse_lopit <- mouse_marker_proteins_df |>
  dplyr::select(-markers)

mouse_pcs <- prcomp(mouse_lopit, center = TRUE, scale. = TRUE)

mouse_df <- as.data.frame(mouse_pcs$x[, c(1, 2)]) |>
  mutate(Marker = mouse_marker_proteins_df$markers, Dataset = "Mouse")

mouse_col_pal <- c(
  "#a7bed3",
  "#70d6ff",
  "#ff9770",
  "#dfb2f4",
  "#f5e960",
  "#ffffff",
  "#79addc",
  "#fcf5c7",
  "#ffa69e",
  "#ff686b",
  "#808080"
)

names(mouse_col_pal) <- levels(mouse_df$Marker)


p_mouse <- mouse_df |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Marker), color = "#808080", shape = 21) +
  facet_wrap(~Dataset) +
  scale_fill_manual(values = mouse_col_pal) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# === Visualisation ============================================================

plot_df <- rbind(callus_df, root_df, human_df, mouse_df)

p_data <- plot_df |>
  ggplot(aes(x = PC1, y = PC2, color = Marker)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Dataset) +
  theme(legend.position = "bottom")

p_data_alt <- (p_callus + p_root) / (p_human + p_mouse) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "bottom")

layout <- c("
  AB
  AB
  CD
  CD
  EE
  EE
  EE
")

data_layout <- c("
  AB
  CD
")

p_data <-  p_callus + p_root + p_human + p_mouse +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = data_layout) 

p_final <- p_callus + p_root + p_human + p_mouse + p_performance +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout) +
  theme(legend.position = "bottom")

ggsave("Plots/fig2ValidationStudy.pdf",
  plot = p_final,
  device = "pdf",
  height = 18.0,
  width = 12.0
)

ggsave("Plots/Fig2/fig2ValidationStudyData.pdf",
  plot = p_data,
  device = "pdf",
  height = 10.0,
  width = 12.0
)


ggsave("Plots/Fig2/fig2ACallus.pdf",
  plot = p_callus,
  device = "pdf",
  height = 4.5,
  width = 6.0
)
ggsave("Plots/Fig2/fig2BRoot.pdf",
  plot = p_root,
  device = "pdf",
  height = 4.5,
  width = 6.0
)
ggsave("Plots/Fig2/fig2CHuman.pdf",
  plot = p_human,
  device = "pdf",
  height = 4.5,
  width = 6.0
)
ggsave("Plots/Fig2/fig2DMouse.pdf",
  plot = p_mouse,
  device = "pdf",
  height = 4.5,
  width = 6.0
)
ggsave("Plots/Fig2/fig2EPerformance.pdf",
  plot = p_performance,
  device = "pdf",
  height = 9.0,
  width = 8.0
)
