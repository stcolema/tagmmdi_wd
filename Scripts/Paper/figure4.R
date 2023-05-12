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

ggPheatmap <- function(X, markers, col_pal, row_order = NULL, include_dendogram = TRUE) {
  # markers <- factor(marrel_col_palkers)
  if (is.null(row_order)) {
    row_order <- findOrder(X)
  }
  data_levels <- c("Coexpression")
  protein_names <- factor(row.names(X), levels = row.names(X)[row_order])
  annotation_row <- data.frame(Localisation = markers, ID = protein_names)

  rel_col_pal <- col_pal[names(col_pal) %in% unique(markers)]
  p_anno <- annotation_row |>
    ggplot(aes(x = NA, y = ID, fill = Localisation)) +
    geom_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    labs(fill = "Localisation", x = "Localisation") +
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

maxpear_cl <- mcclust::maxpear(mdi_mod$cm[[2]], max.k = 50)

pred_cl <- mdi_mod$pred
pred_cl[[2]] <- maxpear_cl$cl
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

# === Golgi disagreements from TAGM ============================================
#
# lopit_disagreement <- tagm_comparison$mdi.mcmc.allocation != tagm_comparison$tagm.mcmc.allocation
# tagm_uncertain <- tagm_comparison$tagm.mcmc.probability < 0.7
# mdi_certain <- tagm_comparison$mdi.mcmc.probability > 0.7
# non_golgi <- tagm_comparison$mdi.mcmc.allocation != "Golgi"
# disagreeing_uncertain_inds <- which(lopit_disagreement & tagm_uncertain & non_golgi & mdi_certain)
#
# diagreement_uncertain_table <- tagm_comparison[disagreeing_uncertain_inds, ]
#
# diagreement_uncertain_table$tagm.mcmc.probability <- diagreement_uncertain_table$tagm.mcmc.probability |>
#   round(digits = 3)
# diagreement_uncertain_table$mdi.mcmc.probability <- diagreement_uncertain_table$mdi.mcmc.probability |>
#   round(digits = 3)
#
# pm_localisation_options <- c("PM - integral", "PM - peripheral 1", "PM - peripheral 2")
#
# golgi_changes_to_pm <- tagm_comparison |>
#   dplyr::filter(tagm.mcmc.allocation == "Golgi", mdi.mcmc.allocation %in% pm_localisation_options)
#
# pm_marker_proteins <- tagm_comparison |>
#   dplyr::filter(markers %in% pm_localisation_options)
#
# golgi_marker_proteins <- tagm_comparison |>
#   dplyr::filter(markers == "Golgi")
#
# golgi_changes_and_pm_markers <- rbind(pm_marker_proteins, golgi_changes_to_pm)
# golgi_changes_and_markers <- rbind(pm_marker_proteins, golgi_changes_to_pm, golgi_marker_proteins)
#
# protein_biochemical_properties <- read.csv("./T_gondii/TGondiiGolgiProteinsAttributes.csv", row.names = 1)
#
# rel_protein_biochemical_properties <- protein_biochemical_properties[row.names(golgi_changes_and_markers), ]
#
# my_df <- cbind(golgi_changes_and_markers, rel_protein_biochemical_properties[, -c(2, 4, 5, 6)]) |>
#   pivot_longer(c(mdi.mcmc.allocation, tagm.mcmc.allocation), names_to = "Model", values_to = "Allocation")
#
# my_df$pI <- as.numeric(my_df$pI)
# my_df$Fixed <- my_df$markers != "unknown"
# my_df |>
#   ggplot(aes(x = Allocation, y = pI, group = Allocation)) +
#   facet_grid(Model ~ Fixed) +
#   geom_boxplot()
#
# p_dNdS <- my_df |>
#   ggplot(aes(x = Allocation, y = dNdS, fill = Allocation)) +
#   facet_grid(~Model) +
#   geom_boxplot() +
#   scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])
#
#
# my_df$Num.TMDs.TMHMM |> table()
# my_df$TMD.within.first.60.AA <- my_df$TMD.within.first.60.AA != "FALSE"
#
# my_df |>
#   filter(!Num.TMDs.TMHMM %in% c("Golgi", "PM - integral", "PM - peripheral 1", "PM - peripheral 2", "unknown")) |>
#   mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
#   ggplot(aes(x = Allocation, y = Num.TMDs.TMHMM, group = Allocation)) +
#   facet_grid(Model ~ Fixed) +
#   geom_boxplot()
#
# my_df |>
#   filter(!Num.TMDs.TMHMM %in% c("Golgi", "PM - integral", "PM - peripheral 1", "PM - peripheral 2", "unknown")) |>
#   mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
#   ggplot(aes(x = Allocation, y = Num.TMDs.TMHMM, group = Allocation)) +
#   facet_grid(~Model) +
#   geom_boxplot()
#
# p_conservation_score <- my_df |>
#   ggplot(aes(x = Allocation, y = Conservation.score, group = Allocation, fill = Allocation)) +
#   facet_grid(~Model) +
#   geom_boxplot() +
#   scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])
#
# p_num_tmds_TMHMM <- my_df |>
#   filter(!Num.TMDs.TMHMM %in% c("Golgi", "PM - integral", "PM - peripheral 1", "PM - peripheral 2", "unknown")) |>
#   mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
#   ggplot(aes(x = Allocation, y = Num.TMDs.TMHMM, group = Allocation, fill = Allocation)) +
#   facet_grid(~Model) +
#   geom_boxplot() +
#   scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])
#
# p_tmd_first_60_aa <- my_df |>
#   ggplot(aes(y = TMD.within.first.60.AA, group = Allocation, fill = Allocation)) +
#   facet_grid(~Model) +
#   stat_count() +
#   scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])
#
# p_patch <- p_dNdS / p_conservation_score / p_num_tmds_TMHMM / p_tmd_first_60_aa + plot_layout(guides = "collect")


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


# plot_df[row.names(golgi_changes_and_markers),] |>
#   dplyr::filter(tSNE_1 > -25) |>
#   ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation, shape = Marker_protein)) +
#   geom_point(aes(alpha = Probability), size = 3) +
#   labs(x = "tSNE 1", y = "tSNE 2", color = "Prediction") +
#   theme(legend.position = "bottom") +
#   scale_color_manual(values = col_pal) +
#   guides(color = guide_legend(ncol = 4))
#
# p_tsne <- plot_df |>
#   ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
#   geom_point(aes(alpha = Probability)) +
#   labs(x = "tSNE 1", y = "tSNE 2", color = "Markers") +
#   theme(legend.position = "bottom") +
#   scale_color_manual(values = col_pal) +
#   guides(color = guide_legend(ncol = 4)) +
#   facet_zoom(xy = predicted_localisation %in% secretory_organelles)
#
microarray_mat <- as.matrix(microarray_data)
#
# microarray_data |>
#   rownames_to_column("Gene") |>
#   mutate(Marker = markers, Predicted_localisation = predicted_localisation, Marker_protein = markers != "all other proteins") |>
#   filter(Gene %in% row.names(golgi_changes_and_markers)) |>
#   pivot_longer(-c(Marker, Predicted_localisation, Marker_protein, Gene), names_to = "Timepoint", values_to = "Expression") |>
#   ggplot(aes(x = Timepoint, y = Expression, group = Gene, color = Marker_protein)) +
#   geom_line() +
#   facet_wrap(~Predicted_localisation)
#
# exprs(Barylyuk2020ToxoLopit)[proteins_modelled,] |>
#   as.data.frame() |>
#   rownames_to_column("Protein") |>
#   mutate(Marker = markers, Predicted_localisation = predicted_localisation, Marker_protein = markers != "all other proteins") |>
#   filter(Protein %in% row.names(golgi_changes_and_markers)) |>
#   pivot_longer(-c(Marker, Predicted_localisation, Marker_protein, Protein), names_to = "Timepoint", values_to = "Expression") |>
#   ggplot(aes(x = Timepoint, y = Expression, group = Protein, color = Marker_protein)) +
#   geom_line() +
#   facet_wrap(~Predicted_localisation)
#   # geom_smooth(method = "loess", se = FALSE)
#
# golgi_reallocations_subset <- match(row.names(golgi_changes_and_markers), row.names(microarray_mat))
# golgi_reallocated_ge_data <- microarray_mat[golgi_reallocations_subset, ]
# golgi_new_predictions <- predicted_localisation[golgi_reallocations_subset]
#
#
# p_ge <- ggPheatmap(golgi_reallocated_ge_data, golgi_new_predictions, col_pal,
#                    include_dendogram = FALSE
# )

# === Secretory organelles =====================================================

secretory_organelles <- c(
  "apical 1",
  "apical 2",
  "micronemes",
  "rhoptries 1",
  "rhoptries 2",
  "ER 1",
  "ER 2"
)

secretory_subset <- predicted_localisation %in% secretory_organelles
secretory_proteins <- row.names(plot_df)[secretory_subset]

secretory_gene_expression_data <- microarray_mat[secretory_subset, ]
secretory_predictions <- predicted_localisation[secretory_subset]


p_tsne_full <- plot_df |>
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
  geom_point(alpha = 0.4) +
  labs(x = "tSNE 1", y = "tSNE 2", color = "Localisation") +
  scale_color_manual(values = col_pal, drop = FALSE) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4))
# facet_zoom(xy = predicted_localisation %in% secretory_organelles)

p_tsne_seccretory <- plot_df |>
  dplyr::filter(predicted_localisation %in% secretory_organelles) |>
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
  geom_point(alpha = 0.4) +
  labs(x = "tSNE 1", y = "tSNE 2", color = "Localisation") +
  scale_color_manual(values = col_pal, drop = FALSE, guide = "none") +
  theme(legend.position = "bottom") # +
# guides(color = guide_legend(ncol = 4))

row_order <- c()
for (organelle in secretory_organelles) {
  curr_proteins <- which(secretory_predictions == organelle)
  organelle_gene_expression_data <- secretory_gene_expression_data[curr_proteins, ]
  ordered_indices <- findOrder(organelle_gene_expression_data)
  row_order <- c(row_order, curr_proteins[ordered_indices])
}

# row_order <- order(secretory_predictions)

p_ge <- ggPheatmap(secretory_gene_expression_data,
  secretory_predictions,
  row_order = row_order,
  col_pal = col_pal,
  include_dendogram = FALSE
) + 
  labs(fill = "Gene\nexpression")

# outlier_rhoptry_proteins <- c("TGME49_294630",
#                               "TGME49_261440",
#                               "TGME49_270320",
#                               "TGME49_246550",
#                               "TGME49_215775",
#                               "TGME49_315160",
#                               "TGME49_209985",
#                               "TGME49_210095"
#                               )
#
# rhoptries_proteins <- which(secretory_predictions %in% c("rhoptries 1", "rhoptries 2"))
# annotatedHeatmap(secretory_gene_expression_data[rhoptries_proteins, ][outlier_rhoptry_proteins, ],
#   secretory_predictions[rhoptries_proteins][outlier_rhoptry_proteins],
#   cluster_cols = FALSE)
# findOrder(secretory_gene_expression_data[rhoptries_proteins, ])
#
# plot_df[outlier_rhoptry_proteins, ] |>
#   ggplot(aes(x = tSNE_1, y = tSNE_2, color = predicted_localisation)) +
#   geom_point(alpha = 0.4) +
#   labs(x = "tSNE 1", y = "tSNE 2", color = "Predicted localisation") +
#   theme(legend.position = "bottom") +
#   scale_color_manual(values = col_pal) +
#   guides(color = guide_legend(ncol = 4)) #+
#
# outlier_micronemes <- plot_df |>
#   filter(tSNE_1 > 17, tSNE_1 < 19, tSNE_2 < -22, tSNE_2 > -25) |>
#   row.names()
#
# outlier_micronemes_inds <- which(names(secretory_predictions) %in% outlier_micronemes)
# annotatedHeatmap(secretory_gene_expression_data[outlier_micronemes_inds, ],
#                  secretory_predictions[outlier_micronemes_inds],
#                  include_rownames = FALSE, cluster_cols = FALSE)

# colnames(microarray_mat) <- paste0("Hour ", seq(0, 12))
# p_ge <- ggPheatmap(microarray_mat, predicted_localisation, col_pal, include_dendogram = FALSE)

# === Biochemical properties ===================================================

protein_biochemical_properties <- read.csv("./T_gondii/TGondiiGolgiProteinsAttributes.csv", row.names = 1)

tagm_secretory_subset <- tagm_predictions %in% secretory_organelles
tagm_secretory_proteins <- row.names(plot_df)[tagm_secretory_subset]
all_secretory_proteins <- c(secretory_proteins, tagm_secretory_proteins) |> unique()

rel_protein_biochemical_properties <- protein_biochemical_properties[all_secretory_proteins, ]
secretory_tSNE_df <- plot_df[all_secretory_proteins, ]

secretory_df <- cbind(secretory_tSNE_df, rel_protein_biochemical_properties[, -c(2, 4)]) |>
  pivot_longer(c(predicted_localisation, tagm.mcmc.allocation), names_to = "Prediction_origin", values_to = "Localisation") |>
  dplyr::mutate(Model = ifelse(Prediction_origin == "predicted_localisation", "MDI", "TAGM"))

secretory_df$pI <- as.numeric(secretory_df$pI)
secretory_df$Fixed <- secretory_df$markers != "all other proteins"
p_pI <- secretory_df |>
  dplyr::filter(Localisation %in% secretory_organelles) |>
  ggplot(aes(x = Localisation, y = log(1 + pI), fill = Model)) +
  # facet_wrap(~Model) +
  geom_boxplot() 
  # scale_fill_manual(values = col_pal)

# x12 <- secretory_df |> 
#   dplyr::filter(Localisation == "ER 2", Model == "MDI") |>
#   dplyr::select(dNdS) 
# hist(x12$dNdS)

p_dNdS <- secretory_df |>
  dplyr::filter(Localisation %in% secretory_organelles) |>
  ggplot(aes(x = Localisation, y = log(1 + dNdS), fill = Model)) +
  geom_boxplot(alpha = 0.5) +
  # facet_wrap(~ Localisation) +
  # scale_fill_manual(values = col_pal) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "bottom")

secretory_df$Num.TMDs.TMHMM |> table()
secretory_df$TMD.within.first.60.AA <- my_df$TMD.within.first.60.AA != "FALSE"

p_num_tmds_TMHMM <- secretory_df |>
  dplyr::filter(Localisation %in% secretory_organelles) |>
  filter(!Num.TMDs.TMHMM %in% c("apical 2", "ER 1", "ER 2", "micronemes", "rhoptries 2", "unknown")) |>
  mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
  ggplot(aes(x = Localisation, y = Num.TMDs.TMHMM, fill = Localisation)) +
  facet_wrap(~Model) +
  geom_boxplot() +
  scale_fill_manual(values = col_pal)

p_conservation_score <- secretory_df |>
  dplyr::filter(Localisation %in% secretory_organelles) |>
  ggplot(aes(x = Localisation, y = Conservation.score, fill = Localisation)) +
  facet_grid(~Model) +
  geom_boxplot() +
  scale_fill_manual(values = col_pal)

p_tmd_first_60_aa <- secretory_df |>
  dplyr::filter(Localisation %in% secretory_organelles) |>
  ggplot(aes(y = TMD.within.first.60.AA, group = Localisation, fill = Localisation)) +
  stat_count() +
  facet_grid(~Model) +
  scale_fill_manual(values = col_pal)

p_patch <- p_dNdS / p_conservation_score / p_num_tmds_TMHMM / p_tmd_first_60_aa + plot_layout(guides = "collect")

# === Save plots ===============================================================
layout <- c("
  AB
  CD
")
p_full <- (p_tsne_full + p_tsne_seccretory + p_ge + p_dNdS) +
  plot_layout(design = layout) # & #, guides = "collect") &
# theme(legend.position = "bottom")
# p_ge

# ggsave("Plots/Fig4/fig4SecretoryOrganelles.png", plot = p_full, height = 12.0, width = 16.0)
# ggsave("Plots/Fig4/fig4CSecretoryOrganellesGE.png", plot = p_ge, height = 6.0, width = 8.0)
# ggsave("Plots/Fig4/fig4ASecretoryOrganellestSNEFull.png", plot = p_tsne_full, height = 6.0, width = 8.5)
# ggsave("Plots/Fig4/fig4BSecretoryOrganellestSNESecretory.png", plot = p_tsne_seccretory, height = 6.0, width = 8.0)

ggsave("Plots/Fig4/fig4CSecretoryOrganellesGE.pdf", plot = p_ge, height = 6.0, width = 8.0)
ggsave("Plots/Fig4/fig4ASecretoryOrganellestSNEFull.pdf", plot = p_tsne_full, height = 6.0, width = 8.5)
ggsave("Plots/Fig4/fig4BSecretoryOrganellestSNESecretory.pdf", plot = p_tsne_seccretory, height = 6.0, width = 8.0)
ggsave("Plots/Fig4/fig4DSecretoryOrganellesdNdS.pdf", plot = p_dNdS, height = 6.0, width = 8.0)

# # === Additional EDA ===========================================================
# 
# outlier_rhoptry_proteins
# 
# rhoptry_inds <- which(predicted_localisation %in% c("rhoptries 1", "rhoptries 2"))
# rhoptry_proteins <- row.names(plot_df)[rhoptry_inds]
# 
# rhoptry_df <- secretory_df[rhoptry_proteins, ] |>
#   mutate(Outlier = rhoptry_proteins %in% outlier_rhoptry_proteins)
# 
# rhoptry_df |>
#   ggplot(aes(x = predicted_localisation, y = Conservation.score, fill = predicted_localisation)) +
#   facet_grid(~Outlier) +
#   geom_boxplot() +
#   scale_fill_manual(values = col_pal)
# 
# 
# rhoptry_df |>
#   ggplot(aes(x = predicted_localisation, y = log(pI), fill = predicted_localisation)) +
#   facet_grid(~Outlier) +
#   geom_boxplot() +
#   scale_fill_manual(values = col_pal)
# 
# rhoptry_df |>
#   ggplot(aes(x = predicted_localisation, y = log(dNdS), fill = predicted_localisation)) +
#   facet_grid(~Outlier) +
#   geom_boxplot() +
#   geom_jitter() +
#   scale_fill_manual(values = col_pal)
