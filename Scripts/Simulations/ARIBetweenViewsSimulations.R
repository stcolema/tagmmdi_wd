
library(ggplot2)
library(tagmReDraft)

mdiHelpR::setMyTheme()

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

p_view_comp <- long_ari_df |>
  ggplot(aes(x = Views, y = ARI, fill = Views)) +
  geom_boxplot() +
  ggthemes::scale_fill_colorblind() +
  ylim(c(0, 1))

ggsave("./Simulations/ComparisonTrueClusterings.png", p_view_comp, height = 6, width = 6)
