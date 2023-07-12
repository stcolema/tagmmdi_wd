
suppressMessages(library(mdir))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

mdiHelpR::setMyTheme()

prepDataForggHeatmap <- function(X, row_order = NULL, col_order = NULL) {
  N <- nrow(X)
  P <- ncol(X)
  
  X_not_matrix <- ! is.matrix(X)
  if(X_not_matrix) {
    stop("X should be a matrix please.")
  }
  
  X_not_DF <- !is.data.frame(X)
  
  cluster_rows <- is.null(row_order)
  cluster_cols <- is.null(col_order)
  
  no_row_ordering <- isFALSE(row_order)
  no_col_ordering <- isFALSE(col_order)
  
  if (X_not_DF) {
    X <- data.frame(X)
  }
  Y <- X
  
  Y$y <- seq(1, N, by = 1)
  # Y$y <- seq(N, 1, by = -1)
  if(no_row_ordering) {
    row_order <- seq(1, N)
  }
  if (cluster_rows) {
    row_order <-  findOrder(X)
  }
  
  Y$y <-match(Y$y, row_order)
  feature_name_order <- colnames(Y)[-ncol(Y)]
  Z <- tidyr::pivot_longer(Y, -tidyr::any_of("y"), names_to = "Feature", values_to = "Entry")
  # Z$x <- as.numeric(stringr::str_extract(Z$Feature, "[:digit:]+"))
  Z$x <- match(Z$Feature, feature_name_order)
  
  if(no_col_ordering) {
    col_order <- Z$x
  }
  if (cluster_cols) {
    col_order <- findOrder(t(X))
  }
  
  Z$x <- match(Z$x, col_order)
  Z$Feature <- NULL
  Z
}


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

big_gg_df$Dataset <- big_gg_df$View
big_labels_gg_df$Dataset <- big_labels_gg_df$View
big_labels_gg_df$View <- big_gg_df$View <- NULL

p_data <- big_gg_df |>
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() +
  facet_grid(Scenario ~ Dataset, scales = "free_x", labeller = label_both) +
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
  facet_grid(Scenario ~ Dataset, scales = "free_x") +
  ggthemes::scale_fill_colorblind() +
  labs(fill = "Label", x = "Dataset", y = "Row number") +
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


data_dir <- "./Simulations/Data/"
scenarios <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
n_scn <- length(scenarios)
for (ii in seq(1, n_scn)) {
  scn <- scenarios[ii]
  scn_input_dir <- paste0(data_dir, scn)
  input_files <- list.files(scn_input_dir, pattern = "*.rds", full.names = TRUE) |>
    stringr::str_sort(numeric = TRUE)
  n_files <- length(input_files)
  for (jj in seq(1, n_files)) {
    f_in <- input_files[jj]
    x <- readRDS(f_in)
    ari_12 <- mcclust::arandi(x$Clustering$View_1, x$Clustering$View_2)
    ari_13 <- mcclust::arandi(x$Clustering$View_1, x$Clustering$View_3)
    ari_23 <- mcclust::arandi(x$Clustering$View_2, x$Clustering$View_3)
    
    .entry <- data.frame("ARI_12" = ari_12, "ARI_13" = ari_13, "ARI_23" = ari_23, Scenario = scn)
    if ((ii == 1) & (jj == 1)) {
      ari_df <- .entry
    } else {
      ari_df <- rbind(ari_df, .entry)
    }
  }
}

long_ari_df <- ari_df |>
  pivot_longer(-Scenario, names_to = "Views", values_to = "ARI")

long_ari_df$Views <- factor(long_ari_df$Views, labels = c("(1, 2)", "(1, 3)", "(2, 3)"))
long_ari_df$Datasets <- long_ari_df$Views


p_view_comp <- long_ari_df |>
  ggplot(aes(x = Datasets, y = ARI, fill = Datasets)) +
  geom_boxplot() +
  ggthemes::scale_fill_colorblind() +
  ylim(c(0, 1))

p_view_comp + p_gg + plot_annotation(tag_levels = "A")
