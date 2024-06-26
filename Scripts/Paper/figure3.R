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
suppressMessages(library(png))

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
    labs(fill = "Marker", x = "Marker") +
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

ggPheatmap <- function(X, markers, col_pal, include_dendogram = TRUE) {
  markers <- factor(markers)
  row_order <- findOrder(X)
  data_levels <- c("Coexpression")
  protein_names <- factor(row.names(X), levels = row.names(X)[row_order])
  annotation_row <- data.frame(Marker = markers, ID = protein_names)
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
    labs(fill = "Marker", x = "Marker") +
    scale_fill_manual(values = col_pal)

  # Order data based on clusters and sample types
  coexpression_mat <- mdir::prepDataForggHeatmap(X, row_order = row_order, col_order = FALSE)
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

  # Combine plots
  plot <- p_facet

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

min_vi_cl <- mcclust.ext::minVI(mdi_mod$cm[[2]], max.k = 50)


pred_cl <- mdi_mod$pred
pred_cl[[2]] <- min_vi_cl$cl
prob_cl <- mdi_mod$prob
fused_genes_1 <- which(colMeans(mdi_mod$allocations[[1]] == mdi_mod$allocations[[2]]) > 0.5)

plotting <- TRUE
plot_height <- 6
plot_width <- 8

# === Input data ===============================================================

microarray_file <- paste0(inputdata_dir, "cellCycleNormalised.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
rna_seq_file <- paste0(inputdata_dir, "rnaSeqMacrophage.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_RNA-Seq.txt")
lopit_file <- paste0(inputdata_dir, "LOPITreduced.csv")
data(Barylyuk2020ToxoLopit)


microarray_data <- read.csv(microarray_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 212)
)

rna_seq_data <- read.csv(rna_seq_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 255)
)

lopit_data <- read.csv(lopit_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 255)
)

data_modelled <- readRDS(paste0(inputdata_dir, "TGondiiMDI_K_125_input.rds"))
datasets <- c("LOPIT", "Cell_cycle", "RNA-seq")

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
plot_df <- data.frame("tSNE_1" = my_tsne$Y[, 1], "tSNE_2" = my_tsne$Y[, 2], organelle = markers)


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

# col_pal <- c(pals::alphabet(), "#808080")
names(col_pal) <- marker_labels # c(names(pals::alphabet()), "grey50")

plot_df$Alpha <- 0.8 * (plot_df$organelle != "all other proteins") + 0.6 * (plot_df$organelle == "all other proteins")
p_tsne <- plot_df |>
  ggplot(aes(x = tSNE_1, y = tSNE_2, fill = organelle), color = "#808080", ) +
  geom_point(aes(alpha = Alpha), shape = 21, size = 2.25) +
  labs(x = "tSNE 1", y = "tSNE 2", fill = "Markers") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = col_pal) +
  guides(color = guide_legend(ncol = 4)) +
  scale_alpha(guide = "none") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

microarray_mat <- as.matrix(microarray_data)
# colnames(microarray_mat) <- paste0("Hour ", seq(0, 12))
p_ge <- ggPheatmap(microarray_mat, markers, col_pal, include_dendogram = FALSE) +
  labs(fill = "Gene\nexpression")

# p_tsne

# p_tsne + p_ge + plot_annotation(tag_levels = "A")

behnke_cell_cycle_diagram_filename <- "/home/stephen/Documents/TAGM_MDI_paper/Figure/1/cellCycleBehnkePlot.png"
behnke_cell_cycle_diagram <- readPNG(behnke_cell_cycle_diagram_filename)

design <- "
 A##BBB
 A#CCC#
"

p_patch <- (p_tsne
# + plot_spacer()
+ p_ge
  # + plot_spacer()
  + grid::rasterGrob(behnke_cell_cycle_diagram)
) + plot_layout(
  design = design,
  widths = c(1.25, 0.1, 0.04, 1.0, 0.2, 0.0005),
  heights = c(1.0, 1.0)
) + plot_annotation(tag_levels = "A")


ggsave("Plots/Fig3/Fig3TGondiiData.pdf",
  plot = p_patch,
  device = "pdf",
  height = 12.0,
  width = 16.0
)

ggsave("Plots/Fig3/fig3ACaseTGondiiDatatSNE.pdf",
  plot = p_tsne,
  device = "pdf",
  height = 12.0,
  width = 8.0
)
ggsave("Plots/Fig3/fig3BCaseTGondiiDataGE.pdf",
  plot = p_ge,
  device = "pdf",
  height = 9.0,
  width = 8.0
)
# lopit_data |>
#   dplyr::mutate(Organelle = markers) |>
#   dplyr::filter(Fixed == 1) |>
#   group_by(Organelle) |>
#   summarise(N = n())
#   pivot_longer(-c("Protein", "Label", "Fixed", "Organelle")) |>
#   ggplot(aes(x = name, y = value, color = Organelle, group = Protein)) +
#   geom_line() +
#   facet_wrap(~Organelle)
