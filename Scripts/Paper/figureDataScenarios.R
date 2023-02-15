library(pheatmap)
library(ggplot2)
library(patchwork)
library(magrittr)
library(mdiHelpR)

mdiHelpR::setMyTheme()
set.seed(1)

N <- 10
P <- 20
K <- c(2, 4)

row_names <- paste0("Gene ", seq(1, N))
col_names <- paste0("Measurement ", seq(1, P))

X2 <- X1 <- matrix(rnorm(N * P), N, P) |> 
  set_rownames(row_names) |> 
  set_colnames(col_names)

gen_data_1 <- generateSimulationDataset(K[1], N, P)
X1 <- gen_data_1$data
labels1 <- gen_data_1$cluster_IDs

gen_data_2 <- generateSimulationDataset(K[2], N, P)
X2 <- gen_data_2$data
labels1 <- gen_data_2$cluster_IDs


labels1 <- sample(seq(1, K[1]), N, replace = TRUE)
labels2 <- sample(seq(1, K[2]), N, replace = TRUE)

col_pal <- dataColPal()
my_breaks <- defineDataBreaks(X1, col_pal)
my_breaks[1:50] <- my_breaks[1:50] - 1e5
my_breaks[52:101] <- my_breaks[52:101] + 1e5
col_pal[50:52] <- col_pal[1]

# Create the annotation data.frame for the rows
anno_row <- data.frame(Localisation = rep(factor(paste("Organelle", 1)), N)) |> 
  magrittr::set_rownames(rownames(X1))

anno_colours <- list("Organelle" = "red")

# my_breaks[1:50] <- 

pheatmap(X1,
         annotation_row = anno_row, 
         annotation_colors = anno_colours,
         legend = FALSE,
         color = col_pal,
         breaks = my_breaks, 
         border_color = NA,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         gaps_row = 5, 
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         filename = "~/Desktop/data_cartoon_base.png"
)

annotatedHeatmap(X1, labels1, cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = 5, col_pal = col_pal, my_breaks = my_breaks)
annotatedHeatmap(X2, labels2, cluster_rows = FALSE, cluster_cols = FALSE)
