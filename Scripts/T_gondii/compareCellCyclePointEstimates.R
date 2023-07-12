
library(pRolocdata)
library(ggplot2)
library(pheatmap)
library(mdiHelpR)
library(magrittr)

mdiHelpR::setMyTheme()

mdi_obj <- readRDS("./T_gondii/ConsensusClustering/CC_R_15000_K_125_full.rds")
mix_obj <- readRDS("./T_gondii/ConsensusClustering/CC_R_15000_K_300_cell_cycle_mix.rds")


mdi_salso <- mdi_obj$cell_cycle_predictions$salso
mdi_vi <- mdi_obj$cell_cycle_predictions$vi

mix_salso <- mix_obj$cell_cycle_predictions$salso
mix_vi <- mix_obj$cell_cycle_predictions$vi

# === Input data ===============================================================

inputdata_dir <- "./T_gondii/Prepared_data/" #
lopit_file <- paste0(inputdata_dir, "LOPITreduced.csv")
data(Barylyuk2020ToxoLopit)

lopit_data <- read.csv(lopit_file, row.names = 1, header = T)

data_modelled <- readRDS(paste0(inputdata_dir, "TGondiiMDI_K_125_input.rds"))
datasets <- c("LOPIT", "Cell_cycle", "RNA-seq")


cols_used <- c("Description", "markers", "tagm.mcmc.allocation", "tagm.mcmc.probability")

proteins_modelled <- row.names(data_modelled$data_modelled[[1]])
tagm_comparison <- fData(Barylyuk2020ToxoLopit)[proteins_modelled, cols_used]

label_to_organelle <- data.frame(
  "Organelle" = levels(tagm_comparison$markers)[-27],
  "Label" = seq(1, 26)
)

marker_inds <- which(tagm_comparison$markers != "unknown")
marker_proteins <- tagm_comparison$markers[marker_inds]
vi_predictions <- mdi_vi[marker_inds]
salso_predictions <- mdi_salso[marker_inds]

mcclust::arandi(marker_proteins, vi_predictions)
mcclust::arandi(marker_proteins, salso_predictions)
mcclust::arandi(vi_predictions, salso_predictions)
mcclust::arandi(mdi_vi, mdi_salso)
mcclust::arandi(mdi_vi, mix_vi)
mcclust::arandi(mdi_salso, mix_salso)

vi_df <- as.data.frame(table(mdi_vi))
salso_df <- as.data.frame(table(mdi_salso))

mix_vi_df <- as.data.frame(table(mix_vi))
mix_salso_df <- as.data.frame(table(mix_salso))

mix_vi_df$Method <- vi_df$Method <- "vi"
mix_salso_df$Method <- salso_df$Method <- "salso"

salso_df$Model <- vi_df$Model <- "MDI"
mix_salso_df$Model <- mix_vi_df$Model <- "Mixture model"

comp_df <- rbind(vi_df[, -1], salso_df[, -1], mix_vi_df[, -1], mix_salso_df[, -1])

comp_df |>
  ggplot(aes(x = Method, y = Freq, fill = Model)) +
  geom_boxplot() +
  labs(y = "Number of members in a cluster") +
  ggthemes::scale_fill_colorblind()

ggsave("./T_gondii/Analysis/ComparisonSALSOVIMembership.png")

p_occ_mdi <- (1 * (mdi_obj$N_k[[2]] != 0)) |>
  pheatmap(cluster_cols = TRUE, cluster_rows = FALSE, main = "MDI", color = simColPal())

p_occ_mix <- (1 * (mix_obj$N_k[[1]] != 0)) |>
  pheatmap(
    cluster_cols = TRUE, cluster_rows = FALSE, main = "Mixture model", color = simColPal()
  )

rowSums((1 * (mdi_obj$N_k[[2]] != 0))) |> mean()
rowSums((1 * (mix_obj$N_k[[1]] != 0))) |> mean()

table
mdi_pred <- mdi_vi
mix_pred <- mix_vi

cl_14 <- which(mdi_pred == 20)

mix_pred_in_cl_14 <- mix_pred[cl_14]
mix_pred_cl_14 <- mdi_pred[cl_14]

data_modelled$data_modelled[[2]][cl_14, ] |>
  annotatedHeatmap(mix_pred_in_cl_14,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    main = NA
  )

# Create the annotation data.frame for the rows
anno_row <- data.frame("Mixture.model" = factor(mix_pred_in_cl_14))
row.names(anno_row) <- row.names(data_modelled$data_modelled[[2]][cl_14, ])

# The number of cololurs to use
K <- length(unique(mix_pred_in_cl_14))

# Create the annotation colours
ann_colours <- list("Mixture.model" = viridis::viridis(K))
names(ann_colours$Mixture.model) <- sort(unique(mix_pred_in_cl_14))

col_pal <- dataColPal()
breaks <- defineDataBreaks(data_modelled$data_modelled[[2]][cl_14, ], dataColPal(), mid_point = 0)
pheatmap::pheatmap(data_modelled$data_modelled[[2]][cl_14, ],
  color = col_pal,
  breaks = breaks,
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  filename = "./T_gondii/Analysis/CellCycleMDICluster20AnnotatedByMixtureModel.png"
)

pred_cell_cycle <- data.frame("MDI" = c(mdi_pred), "Mixture.model" = c(mix_pred))

overlap_mat <- pred_cell_cycle |>
  table() |>
  apply(1, function(x) {
    x / sum(x)
  })

overlap_mat |>
  pheatmap(
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    # main = "Fraction of mixture model clusters in MDI clusters",
    color = simColPal()
    # filename = "./T_gondii/Analysis/cellCycleMDIMixtureOverlap.png"
  )

overlap_mat |> apply(2, quantile)
