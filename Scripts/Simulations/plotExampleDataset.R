
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

mdiHelpR::setMyTheme()

data_dir <- "./Simulations/Data/"
scenarios <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
n_scn <- length(scenarios)
jj <- 1
for (ii in seq(1, n_scn)) {
  scn <- scenarios[ii]
  scn_input_dir <- paste0(data_dir, scn)
  input_files <- list.files(scn_input_dir, pattern = "*.rds", full.names = TRUE) |>
    stringr::str_sort(numeric = TRUE)
  n_files <- length(input_files)
  # for (jj in seq(1, n_files)) {
  f_in <- input_files[jj]
  x <- readRDS(f_in)

  row_order <- findOrder(x$Data$View_1)

  view_1_gg_df <- x$Data$View_1 |>
    prepDataForggHeatmap(row_order = row_order) |>
    mutate(View = 1)
  view_2_gg_df <- x$Data$View_2 |>
    prepDataForggHeatmap(row_order = NULL) |>
    mutate(View = 2)
  view_3_gg_df <- x$Data$View_3 |>
    prepDataForggHeatmap(row_order = NULL) |>
    mutate(View = 3)

  gg_df <- rbind(view_1_gg_df, view_2_gg_df, view_3_gg_df)
  gg_df$View <- factor(gg_df$View)
  gg_df$Scenario <- scn

  labels <- do.call(cbind, x$Clustering)

  labels_gg_df <- labels |>
    prepDataForggHeatmap(row_order = row_order, col_order = c(1, 2, 3)) |>
    mutate(View = "Labels")

  labels_gg_df$Scenario <- scn

  if (ii == 1) {
    big_gg_df <- gg_df
    big_labels_gg_df <- labels_gg_df
  } else {
    big_gg_df <- rbind(big_gg_df, gg_df)
    big_labels_gg_df <- rbind(big_labels_gg_df, labels_gg_df)
  }
}

p_data <- big_gg_df |>
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() +
  facet_grid(Scenario ~ View, scales = "free_x", labeller = label_both) +
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900") +
  labs(y = NULL, x = "Feature") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank()
  )

p_labels <- big_labels_gg_df |>
  ggplot(aes(x = x, y = y, fill = factor(Entry))) +
  geom_tile() +
  facet_grid(Scenario ~ View, scales = "free_x") +
  ggthemes::scale_fill_colorblind() +
  labs(fill = "Label", x = "View", y = "Row number") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  )

p_gg <- p_labels + p_data +
  plot_layout(widths = c(1, 6), guides = "collect")
# plot_annotation(title = gg_title)

ggsave("./Simulations/Simulation1AllScenarios.png", p_gg, height = 6, width = 8)
