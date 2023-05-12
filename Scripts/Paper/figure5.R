#!/usr/bin/Rscript
#
# Example call:
#  Rscript processTGondiiCCoutput.R --data_dir ~/rds/hpc-work/tagmmdi/T_gondii/Data/ --save_dir ~/rds/hpc-work/tagmmdi/T_gondii/ConsensusClustering/ --model_output_dir ~/rds/hpc-work/tagmmdi/T_gondii/ConsensusClustering/ --R 5000 --K 125


suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(mdir))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))
suppressMessages(library(ggforce))

# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--data_dir"),
      type = "character",
      help = "Directory here the data that were modelled are stored.",
      metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--save_dir"),
      type = "character",
      help = "Directory to save the outputs of this file to.",
      metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--model_output_dir"),
      type = "character",
      help = "Directory where the model outputs are saved.",
      metavar = "character"
    ),

    # Number of MCMC iterations
    optparse::make_option(c("-R", "--R"),
      type = "numeric",
      default = 5000,
      help = "Number of iterations to run in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-b", "--burn"),
      type = "numeric",
      default = NULL,
      help = paste(
        "Number of iterations to burn of the warm-up period in each MCMC chain",
        "[default= 0.2 * R]"
      ),
      metavar = "numeric"
    ),
    optparse::make_option(c("-K", "--K"),
      type = "numeric",
      default = 125,
      help = "Number of components modelled in the unsupervised views.",
      metavar = "numeric"
    )
  )

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

#' @title Creates an annoted heatmap of simulated data
#' @description Plots the output of ``generateBatchData`` as an annotated
#' heatmap
#' @param X Output of ``generateBatchData``.
#' @param row_order Order for rows to be represented in heatmap along the
#' y-axis.
#' @param col_order Order for columns to be represented in heatmap along the
#' x-axis.
#' @param include_dendogram Flag indicating if the dendogram associated with the
#'  ordering should be included in the final patchwork plot.
#' @returns A long data.frame containing columns `Feature` (the x-axis position
#' of the entry for geom_tile()), `Item` (the y-axis position of the entry for
#' geom_tile()), and `Entry` (value in similarity  matrix).
heatmapAnnotatedByClassAndBatch <- function(X, include_dendogram = TRUE, lead = "Coexpression") {
  labels <- factor(X$labels)
  markers <- factor(X$markers)

  X_co <- X$coexpression_data
  X_ab <- X$abundance_data

  row_order <- NA
  if (lead == "Coexpression") {
    row_order <- findOrder(X_co)
    data_levels <- c("Coexpression", "Abundance")
  }
  if (lead == "Abundance") {
    row_order <- findOrder(X_ab)
    data_levels <- c("Abundance", "Coexpression")
  }

  no_order_applied <- FALSE
  length_of_one <- length(row_order) == 1
  if (length_of_one) {
    no_order_applied <- is.na(row_order)
  }
  if (no_order_applied) {
    stop("`lead` must be one of 'Coexpression' or 'Abundance'.")
  }
  col_order <- FALSE

  protein_names <- factor(row.names(X_co), levels = row.names(X_co)[row_order])

  annotation_row_1 <- data.frame(Cluster = labels, ID = protein_names)
  annotation_row_2 <- data.frame(Marker = markers, ID = protein_names)

  n_colors <- length(unique(labels))
  if (n_colors < 8) {
    large_quant_pal <- ggthemes::colorblind_pal()(n_colors + 1)[seq(2, n_colors + 1)]
  } else {
    large_quant_pal <- Polychrome::green.armytage.colors(n_colors)
    # large_quant_pal <- ggsci::pal_ucscgb()(length(unique(labels)))
  }

  names(large_quant_pal) <- NULL
  p_anno_1 <- annotation_row_1 |>
    ggplot(aes(x = NA, y = ID, fill = Cluster)) +
    geom_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    labs(fill = "Cluster", x = "Cluster") +
    scale_fill_manual(values = large_quant_pal)


  p_anno_2 <- annotation_row_2 |>
    ggplot(aes(x = NA, y = ID, fill = Marker)) +
    geom_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    labs(fill = "Cluster", x = "Cluster") +
    scale_fill_manual(values = c("#FFFFFF", "#146EB4"))
  # ggthemes::scale_fill_colorblind()

  # Order data based on clusters and sample types
  coexpression_mat <- tagmReDraft::prepDataForggHeatmap(X_co,
    row_order = row_order,
    col_order = col_order
  )

  abundance_mat <- tagmReDraft::prepDataForggHeatmap(X_ab,
    row_order = row_order,
    col_order = col_order
  )

  coexpression_mat$y <- factor(coexpression_mat$y)
  abundance_mat$y <- factor(abundance_mat$y)

  coexpression_mat$x <- factor(coexpression_mat$x, labels = colnames(X_co))
  abundance_mat$x <- factor(abundance_mat$x, labels = colnames(X_ab))

  abundance_mat$Type <- "Abundance"
  coexpression_mat$Type <- "Coexpression"
  p_df <- rbind(coexpression_mat, abundance_mat)

  p_df$Type <- factor(p_df$Type, levels = data_levels)

  # Color pallette and color breaks
  col_pal <- dataColPal()
  breaks <- defineDataBreaks(X$coexpression_data, col_pal = col_pal, mid_point = 0.0)

  # Range in data
  limits <- c(min(c(X_co, X_ab)), max(c(X_co, X_ab)))

  # Data heatmap
  p_facet <- p_df |>
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = Entry)) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      # axis.title.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(0, 0, 0, 0), # removes margins
      legend.key.height = unit(1, "cm")
    ) +
    labs(x = "Time (hours)", y = NA) +
    scale_fill_gradient2(
      low = col_pal[1],
      mid = "white", # col_pal[50],
      high = col_pal[100],
      limits = limits
    ) +
    facet_wrap(~Type)

  # Layout of patchwork
  design <- "ABCCCCC"

  # Combine plots
  plot <- p_anno_2 + p_anno_1 + p_facet +
    plot_layout(
      design = design,
      guides = "collect", # Specify layout, collect legends

      # Adjust widths and heights to align plots.
      # When annotation plot is larger, it might not fit into its column/row.
      # Then you need to make column/row larger.

      # Relative widths and heights of each column and row:
      # Currently, the width of the first column is 15 % and the height of
      # first two rows are 30 % the size of others

      # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal,
      # but then the gaps between plots are larger.
      widths = c(0.075, 0.075, 1.0)
    )
  if (include_dendogram) {
    # Hierarchical clustering
    taxa_hclust <- hclust(dist(X_co), method = "complete")

    # Creates a phylogenetic tree
    taxa_tree <- ape::as.phylo(taxa_hclust)

    # Plot taxa tree
    taxa_tree <- ggtree::ggtree(taxa_tree) +
      theme(plot.margin = ggplot2::margin(0, 0, 0, 0)) # removes margins

    # Get order of taxa in plot
    taxa_ordered <- ggtree::get_taxa_name(taxa_tree)

    # Create layout
    design <- c(
      "ABCDDDDDD"
    )

    # Combine plots
    plot <- taxa_tree + p_anno_1 + p_anno_2 + p_facet +
      plot_layout(
        design = design, guides = "collect", # Specify layout, collect legends

        # Adjust widths and heights to align plots.
        # When annotation plot is larger, it might not fit into its column/row.
        # Then you need to make column/row larger.

        # Relative widths and heights of each column and row:
        # Currently, the width of the first column is 15 % and the height of
        # first two rows are 30 % the size of others

        # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal,
        # but then the gaps between plots are larger.
        widths = c(0.15, 0.15, 0.15, 1.0),
        heights = c(1.0, 1.0, 1.0, 1.0)
      )
  }
  plot
}

ggPheatmap <- function(X, markers, col_pal, row_order = NULL, include_dendogram = TRUE) {
  # markers <- factor(marrel_col_palkers)
  if (is.null(row_order)) {
    row_order <- findOrder(X)
  }
  data_levels <- c("Coexpression")
  protein_names <- factor(row.names(X), levels = row.names(X)[row_order])
  annotation_row <- data.frame(Marker = markers, ID = protein_names)

  rel_col_pal <- col_pal[names(col_pal) %in% unique(markers)]
  p_anno <- annotation_row |>
    ggplot(aes(x = NA, y = ID, fill = Marker)) +
    geom_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    labs(fill = "Cluster", x = "Cluster") +
    scale_fill_manual(values = rel_col_pal)

  # Order data based on clusters and sample types
  coexpression_mat <- MDIr::prepDataForggHeatmap(X, row_order = row_order, col_order = FALSE)
  coexpression_mat$Type <- "Coexpression"
  coexpression_mat$Feature <- factor(coexpression_mat$Feature,
    levels = levels(coexpression_mat$Feature),
    labels = paste0(seq(0, 12))
  )

  # Color pallette and color breaks
  col_pal <- dataColPal()
  breaks <- defineDataBreaks(X, col_pal = col_pal, mid_point = 0.0)

  # Range in data
  limits <- c(min(X), max(X))

  # Data heatmap
  p_facet <- coexpression_mat |>
    ggplot() +
    geom_tile(aes(x = Feature, y = Item, fill = Entry)) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      # axis.title.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(0, 0, 0, 0), # removes margins
      legend.key.height = unit(1, "cm")
    ) +
    labs(x = "Time (hours)", y = NA) +
    scale_fill_gradient2(
      low = col_pal[1],
      mid = "white", # col_pal[50],
      high = col_pal[100],
      limits = limits
    )

  # Layout of patchwork
  design <- "A"

  # Layout of patchwork
  design <- "ACCCC"

  # Combine plots
  plot <- p_anno + p_facet +
    plot_layout(
      design = design,
      guides = "collect", # Specify layout, collect legends

      # Adjust widths and heights to align plots.
      # When annotation plot is larger, it might not fit into its column/row.
      # Then you need to make column/row larger.

      # Relative widths and heights of each column and row:
      # Currently, the width of the first column is 15 % and the height of
      # first two rows are 30 % the size of others

      # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal,
      # but then the gaps between plots are larger.
      widths = c(0.075, 1.0)
    )

  # Combine plots
  # plot <- p_facet

  if (include_dendogram) {
    # Hierarchical clustering
    taxa_hclust <- hclust(dist(X), method = "complete")

    # Creates a phylogenetic tree
    taxa_tree <- ape::as.phylo(taxa_hclust)

    # Plot taxa tree
    taxa_tree <- ggtree::ggtree(taxa_tree) +
      theme(plot.margin = ggplot2::margin(0, 0, 0, 0)) # removes margins

    # Get order of taxa in plot
    taxa_ordered <- ggtree::get_taxa_name(taxa_tree)

    # Create layout
    design <- c(
      "ABBBBBB"
    )

    # Combine plots
    plot <- taxa_tree + p_facet +
      plot_layout(
        design = design, guides = "collect", # Specify layout, collect legends

        # Adjust widths and heights to align plots.
        # When annotation plot is larger, it might not fit into its column/row.
        # Then you need to make column/row larger.

        # Relative widths and heights of each column and row:
        # Currently, the width of the first column is 15 % and the height of
        # first two rows are 30 % the size of others

        # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal,
        # but then the gaps between plots are larger.
        widths = c(0.15, 1.0),
        heights = c(1.0, 1.0)
      )
  }
  plot
}


# === MAIN =====================================================================

cat("\n=== Running ``processTGondiiCCoutput.R`` ===============================")
setMyTheme()
set.seed(1)

args <- input_arguments()

# Directories for input and output respectively
inputdata_dir <- "./T_gondii/Prepared_data/" #
save_dir <- "./T_gondii/Analysis/" #
model_output_dir <- "./T_gondii/ConsensusClustering/" #

# === READ IN MODEL OUTPUT =====================================================
cat("\n=== Reading in files ===================================================")

mdi_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_125_full.rds"

D <- 15000
W <- 150
mdi_mod <- readRDS(mdi_file)
# mix_mod <- readRDS(mix_file)
V <- 2

vi_cl <- mcclust.ext::minVI(mdi_mod$cm[[2]], max.k = 50)
vi_cb <- mcclust.ext::credibleball(vi_cl$cl, mdi_mod$allocations[[2]])

# maxpear_cl <- mcclust::maxpear(mdi_mod$cm[[2]], max.k = 50)

pred_cl <- mdi_mod$pred
# pred_cl[[2]] <- maxpear_cl$cl
pred_cl[[2]] <- vi_cl$cl

prob_cl <- mdi_mod$prob
fused_genes_1 <- which(colMeans(mdi_mod$allocations[[1]] == mdi_mod$allocations[[2]]) > 0.5)

plotting <- TRUE
plot_height <- 6
plot_width <- 8

# === Input data ===============================================================

microarray_file <- paste0(inputdata_dir, "cellCycleNormalised.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
lopit_file <- paste0(inputdata_dir, "LOPITreduced.csv")
data(Barylyuk2020ToxoLopit)


microarray_data <- read.csv(microarray_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 212)
)

lopit_data <- read.csv(lopit_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 255)
)

data_modelled <- readRDS(paste0(inputdata_dir, "TGondiiMDI_K_125_input.rds"))
datasets <- c("LOPIT", "Cell_cycle")

# Proteins in final analysis
proteins_modelled <- row.names(data_modelled$data_modelled[[1]])

# Relevant columns from the pRolocdata object
cols_used <- c("Description", "markers", "tagm.mcmc.allocation", "tagm.mcmc.probability")

# Relevant data for comparing allocations
tagm_comparison <- fData(Barylyuk2020ToxoLopit)[proteins_modelled, cols_used]

# Map between symbolic labels and localisation name
label_to_organelle <- data.frame(
  "Organelle" = levels(tagm_comparison$markers)[-27],
  "Label" = seq(1, 26)
)


marker_labels <- label_to_organelle$Organelle[data_modelled$initial_labels[, 1]]
marker_labels[data_modelled$fixed[, 1] == 0] <- NA

# The MDI predictions and probability of allocation
mdi_predictions <- label_to_organelle$Organelle[pred_cl[[1]]]
mdi_probabilities <- prob_cl[[1]]

tagm_comparison <- tagm_comparison[proteins_modelled, ]
tagm_comparison$mdi.mcmc.allocation <- mdi_predictions
tagm_comparison$mdi.mcmc.probability <- mdi_probabilities
tagm_predictions <- tagm_comparison$tagm.mcmc.allocation

# === T-SNE ====================================================================

markers <- label_to_organelle$Organelle[data_modelled$initial_labels[, 1]]
markers[data_modelled$fixed[, 1] == 0] <- "unknown"

my_tsne <- Rtsne::Rtsne(lopit_data[, seq(1, 30)])

marker_levels <- c(
  "apical 1",
  "apical 2",
  "micronemes",
  "rhoptries 1",
  "rhoptries 2",
  "dense granules",
  "IMC",
  "tubulin cytoskeleton",
  "PM - peripheral 1",
  "endomembrane vesicles",
  "PM - peripheral 2",
  "PM - integral",
  "Golgi",
  "ER",
  "ER 2",
  "apicoplast",
  "mitochondrion - membranes",
  "mitochondrion - soluble",
  "nucleus - chromatin",
  "nucleus - non-chromatin",
  "nucleolus",
  "40S ribosome",
  "60S ribosome",
  "cytosol",
  "19S proteasome",
  "20S proteasome",
  "unknown"
)
marker_labels <- c(
  "apical 1",
  "apical 2",
  "micronemes",
  "rhoptries 1",
  "rhoptries 2",
  "dense granules",
  "IMC",
  "tubulin cytoskeleton",
  "PM - peripheral 1",
  "endomembrane vesicles",
  "PM - peripheral 2",
  "PM - integral",
  "Golgi",
  "ER 1",
  "ER 2",
  "apicoplast",
  "mitochondrion - membranes",
  "mitochondrion - soluble",
  "nucleus - chromatin",
  "nucleus - non-chromatin",
  "nucleolus",
  "40S ribosome",
  "60S ribosome",
  "cytosol",
  "19S proteasome",
  "20S proteasome",
  "all other proteins"
)

markers <- factor(markers, levels = marker_levels, labels = marker_labels)
predicted_localisation <- factor(tagm_comparison$mdi.mcmc.allocation,
  levels = marker_levels[-27],
  labels = marker_labels[-27]
)

localisation_probability <- tagm_comparison$mdi.mcmc.probability

names(predicted_localisation) <- row.names(tagm_comparison)

plot_df <- data.frame(
  "tSNE_1" = my_tsne$Y[, 1],
  "tSNE_2" = my_tsne$Y[, 2],
  markers = markers,
  predicted_localisation = predicted_localisation,
  Probability = localisation_probability
) |>
  mutate(Marker_protein = markers != "all other proteins")
row.names(plot_df) <- row.names(tagm_comparison)

col_pal <- c(pals::alphabet())
names(col_pal) <- marker_labels[-27] # c(names(pals::alphabet()), "grey50")


# plot_df[row.names(golgi_changes_and_markers), ] |>
#   dplyr::filter(tSNE_1 > -25) |>
#   ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation, shape = Marker_protein)) +
#   geom_point(aes(alpha = Probability), size = 3) +
#   labs(x = "tSNE 1", y = "tSNE 2", color = "Prediction") +
#   theme(legend.position = "bottom") +
#   scale_color_manual(values = col_pal) +
#   guides(color = guide_legend(ncol = 4))

# p_tsne <- plot_df |>
#   ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
#   geom_point(aes(alpha = Probability)) +
#   labs(x = "tSNE 1", y = "tSNE 2", color = "Markers") +
#   theme(legend.position = "bottom") +
#   scale_color_manual(values = col_pal) +
#   guides(color = guide_legend(ncol = 4)) +
#   facet_zoom(xy = predicted_localisation %in% secretory_organelles)

microarray_mat <- as.matrix(microarray_data)

# microarray_data |>
#   rownames_to_column("Gene") |>
#   mutate(Marker = markers, Predicted_localisation = predicted_localisation, Marker_protein = markers != "all other proteins") |>
#   filter(Gene %in% row.names(golgi_changes_and_markers)) |>
#   pivot_longer(-c(Marker, Predicted_localisation, Marker_protein, Gene), names_to = "Timepoint", values_to = "Expression") |>
#   ggplot(aes(x = Timepoint, y = Expression, group = Gene, color = Marker_protein)) +
#   geom_line() +
#   facet_wrap(~Predicted_localisation)

# exprs(Barylyuk2020ToxoLopit)[proteins_modelled, ] |>
#   as.data.frame() |>
#   rownames_to_column("Protein") |>
#   mutate(Marker = markers, Predicted_localisation = predicted_localisation, Marker_protein = markers != "all other proteins") |>
#   filter(Protein %in% row.names(golgi_changes_and_markers)) |>
#   pivot_longer(-c(Marker, Predicted_localisation, Marker_protein, Protein), names_to = "Timepoint", values_to = "Expression") |>
#   ggplot(aes(x = Timepoint, y = Expression, group = Protein, color = Marker_protein)) +
#   geom_line() +
#   facet_wrap(~Predicted_localisation)
# geom_smooth(method = "loess", se = FALSE)

# === Dense granules ===========================================================

dense_granule_organelle <- c(
  "dense granules"
)

dg_subset <- predicted_localisation %in% dense_granule_organelle
dg_proteins <- row.names(plot_df)[dg_subset]
dg_gene_expression_data <- microarray_mat[dg_proteins, ]
dg_predictions <- predicted_localisation[dg_subset]

p_tsne <- plot_df |>
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
  geom_point(alpha = 0.8) +
  labs(x = "tSNE 1", y = "tSNE 2", color = "Markers") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = col_pal) +
  guides(color = guide_legend(ncol = 4)) # +
  # facet_zoom(xy = predicted_localisation %in% dense_granule_organelle)

p_tsne_dg <- plot_df |>
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
  geom_point(alpha = 0.8) +
  labs(x = "tSNE 1", y = "tSNE 2", color = "Markers") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = col_pal) +
  guides(color = guide_legend(ncol = 4)) +
  coord_cartesian(xlim = c(-0, 20), ylim = c(-12, 17.5)) 
# facet_zoom(xy = predicted_localisation %in% dense_granule_organelle)

p_tsne_dg2 <- plot_df |>
  mutate(GE_cluster = factor(pred_cl[[2]])) |> 
  filter(predicted_localisation == "dense granules") |> 
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = GE_cluster)) +
  geom_point(alpha = 0.8) +
  labs(x = "tSNE 1", y = "tSNE 2", color = "Markers") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = col_pal) +
  guides(color = guide_legend(ncol = 4)) +
  coord_cartesian(xlim = c(-0, 20), ylim = c(-12, 17.5)) 
# facet_zoom(xy = predicted_localisation %in% dense_granule_organelle)


row_order <- c()
for (organelle in dense_granule_organelle) {
  curr_proteins <- which(dg_predictions == organelle)
  organelle_gene_expression_data <- dg_gene_expression_data[curr_proteins, ]
  ordered_indices <- findOrder(organelle_gene_expression_data)
  row_order <- c(row_order, curr_proteins[ordered_indices])
}

# row_order <- order(dg_predictions)
dg_clusters_ge <- factor(pred_cl[[2]][dg_subset])
ge_col_pal <- Polychrome::green.armytage.colors(length(unique(dg_clusters_ge)))
names(ge_col_pal) <- levels(dg_clusters_ge)
p_ge <- ggPheatmap(dg_gene_expression_data,
  dg_clusters_ge,
  row_order = row_order,
  col_pal = ge_col_pal,
  include_dendogram = FALSE
) +
  labs(fill = "Gene\nexpression")

kept_clusters <- levels(dg_clusters_ge)[table(dg_clusters_ge) > 10]
kept_genes <- which(dg_clusters_ge %in% kept_clusters)

# annotatedHeatmap(dg_gene_expression_data[kept_genes, ], dg_clusters_ge[kept_genes])

# p_tsne
# p_ge

# group_1 <- seq(1, 15)
# group_2 <- seq(71, 78)
# group_3 <- seq(110, 134)
# group_4 <- c(seq(16, 70), seq(79, 109))
# 
# n_group_1 <- length(group_1)
# n_group_2 <- length(group_2)
# n_group_3 <- length(group_3)
# n_group_4 <- length(group_4)

# group_ordering <- c(
#   group_1,
#   group_2,
#   group_3,
#   group_4
# )
# 
# merged_clusters <- c(
#   rep(1, n_group_1),
#   rep(2, n_group_2),
#   rep(3, n_group_3),
#   rep(4, n_group_4)
# )
# 
# dg_gene_expression_data[row_order, ][group_1, ] |> pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = FALSE)
# dg_gene_expression_data[row_order, ][group_2, ] |> pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = FALSE)
# dg_gene_expression_data[row_order, ][group_3, ] |> pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = FALSE)
# dg_gene_expression_data[row_order, ][group_4, ] |> pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = FALSE)
# 
# clusters_considered <- as.numeric(names(which(table(pred_cl[[2]][dg_subset]) > 2)))
# final_dg_gene_set <- pred_cl[[2]][dg_subset][row_order] %in% clusters_considered
# 
# annotatedHeatmap(dg_gene_expression_data[row_order, ][final_dg_gene_set, ], pred_cl[[2]][dg_subset][row_order][final_dg_gene_set], cluster_cols = FALSE, cluster_rows = FALSE)
# pheatmap::pheatmap(dg_gene_expression_data[row_order, ][group_ordering, ], cluster_cols = FALSE, cluster_rows = FALSE)
# annotatedHeatmap(dg_gene_expression_data[row_order, ][group_ordering, ], merged_clusters, cluster_cols = FALSE, cluster_rows = FALSE)
# 
# lg_hg_clusters <- c(merged_clusters == 4) * 1 + 1
# lg_hg_clusters <- c((merged_clusters == 3 * 1) + ((merged_clusters == 4) * 2) + 1)
# annotatedHeatmap(dg_gene_expression_data[row_order, ][group_ordering, ], lg_hg_clusters, cluster_cols = FALSE, cluster_rows = FALSE)

# === Biochemical properties ===================================================

protein_biochemical_properties <- read.csv("./T_gondii/TGondiiGolgiProteinsAttributes.csv", row.names = 1)

tagm_dg_subset <- tagm_predictions %in% dense_granule_organelle
tagm_dg_proteins <- row.names(plot_df)[tagm_dg_subset]
all_dg_proteins <- c(dg_proteins, tagm_dg_proteins) |> unique()

rel_protein_biochemical_properties <- protein_biochemical_properties[dg_proteins, ]
dg_tSNE_df <- plot_df[dg_proteins, ]

dg_df <- cbind(dg_tSNE_df, rel_protein_biochemical_properties[, -c(2, 4)])
dg_df$GE_cluster <- dg_clusters_ge
# |>
# pivot_longer(c(mdi.mcmc.allocation, tagm.mcmc.allocation), names_to = "Model", values_to = "Allocation")

reduced_dg_df <- dg_df |> 
  filter(GE_cluster %in% kept_clusters)

# reduced_dg_df <- dg_df[row_order, ][group_ordering, ]
# reduced_dg_df$Merged_cluster <- factor(lg_hg_clusters)

reduced_dg_df$pI <- as.numeric(reduced_dg_df$pI)
dg_df$Fixed <- dg_df$markers != "all other proteins"
p_pI <- reduced_dg_df |>
  ggplot(aes(x = GE_cluster, y = log(1 + pI), fill = GE_cluster)) +
  # facet_wrap(~ Fixed) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual(values = ge_col_pal)

p_dNdS <- reduced_dg_df |>
  ggplot(aes(x = GE_cluster, y = log(1 + dNdS), fill = GE_cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = ge_col_pal) +
  labs(x = "Inferred cluster") +
  theme(legend.position = "None")

# dg_df$Num.TMDs.TMHMM |> table()
# dg_df$TMD.within.first.60.AA <- my_df$TMD.within.first.60.AA != "FALSE"

# p_num_tmds_TMHMM <- dg_df |>
#   filter(!Num.TMDs.TMHMM %in% c("apical 2", "ER 1", "ER 2", "micronemes", "rhoptries 2", "unknown")) |>
#   mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
#   ggplot(aes(x = predicted_localisation, y = Num.TMDs.TMHMM, fill = predicted_localisation)) +
#   facet_wrap(~Fixed) +
#   geom_boxplot() +
#   scale_fill_manual(values = col_pal)

# p_conservation_score <- dg_df |>
#   ggplot(aes(x = predicted_localisation, y = Conservation.score, fill = predicted_localisation)) +
#   # facet_grid(~Model) +
#   geom_boxplot() +
#   scale_fill_manual(values = col_pal)
# 
# p_tmd_first_60_aa <- dg_df |>
#   ggplot(aes(y = TMD.within.first.60.AA, group = predicted_localisation, fill = predicted_localisation)) +
#   stat_count() +
#   scale_fill_manual(values = col_pal)

# === Save plots ==============================================================
# p_patch <- p_dNdS / p_conservation_score / p_num_tmds_TMHMM / p_tmd_first_60_aa + plot_layout(guides = "collect")

layout <- c("
  BC
  DC
")

p_full <- (p_tsne_dg + p_ge + p_dNdS) +
  plot_layout(design = layout) # & #, guides = "collect") &
# theme(legend.position = "bottom")
# p_ge

ggsave("Plots/fig5DenseGranules.png", plot = p_full, height = 10.0, width = 16.0)

ggsave("Plots/Fig5/fig5AtSNE.png", plot = p_tsne_dg, height = 5.0, width = 8.0)
ggsave("Plots/Fig5/fig5BGE.png", plot = p_ge, height = 10.0, width = 8.0)
ggsave("Plots/Fig5/fig5CdNdS.png", plot = p_dNdS, height = 5.0, width = 8.0)

ggsave("Plots/Fig5/fig5AtSNE.pdf", plot = p_tsne_dg, height = 5.0, width = 8.0)
ggsave("Plots/Fig5/fig5BGE.pdf", plot = p_ge, height = 10.0, width = 8.0)
ggsave("Plots/Fig5/fig5CdNdS.pdf", plot = p_dNdS, height = 5.0, width = 8.0)
