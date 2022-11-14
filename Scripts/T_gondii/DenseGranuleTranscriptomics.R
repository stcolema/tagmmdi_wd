#!/usr/bin/Rscript
#
# Example call:
#  Rscript processTGondiiCCoutput.R --data_dir ~/rds/hpc-work/tagmmdi/T_gondii/Data/ --save_dir ~/rds/hpc-work/tagmmdi/T_gondii/ConsensusClustering/ --model_output_dir ~/rds/hpc-work/tagmmdi/T_gondii/ConsensusClustering/ --R 5000 --K 125


suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))
suppressMessages(library(mdiHelpR))

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
  if(lead == "Coexpression") { 
    row_order <- findOrder(X_co)
    data_levels <- c("Coexpression", "Abundance")
  } 
  if(lead == "Abundance") {
    row_order <- findOrder(X_ab)
    data_levels <- c("Abundance", "Coexpression")
  }
  
  no_order_applied <- FALSE
  length_of_one <- length(row_order) == 1
  if(length_of_one) {
    no_order_applied <- is.na(row_order)
  }
  if(no_order_applied) {
    stop("`lead` must be one of 'Coexpression' or 'Abundance'.")
  }
  col_order <- FALSE
  
  protein_names <- factor(row.names(X_co), levels = row.names(X_co)[row_order])
  
  annotation_row_1 <- data.frame(Cluster = labels, ID = protein_names)
  annotation_row_2 <- data.frame(Marker = markers, ID = protein_names)
  
  n_colors <- length(unique(labels))
  if(n_colors < 8) {
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


cat("\n=== Running ``processTGondiiCCoutput.R`` ===============================")
setMyTheme()

args <- input_arguments()

# Directories for input and output respectively
inputdata_dir <- "./T_gondii/Prepared_data/" #
save_dir <- "./T_gondii/Analysis/" #
model_output_dir <- "./T_gondii/ConsensusClustering/" #

# === READ IN MODEL OUTPUT =====================================================
cat("\n=== Reading in files ===================================================")

mdi_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_125_full.rds"
mix_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_300_cell_cycle_mix.rds"

D <- 15000
W <- 150
mdi_mod <- readRDS(mdi_file)
mix_mod <- readRDS(mix_file)
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

# === LOPIT ===================================================================

# Compare MDI and TAGM allocations

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

# Check which types of disagreements emerge
tagm_comparison[, c(3, 5)] |> unique()
tagm_comparison[fused_genes_1, c(3, 5)] |> unique()
tagm_comparison[
  fused_genes_1[which(tagm_comparison[fused_genes_1, ]$mdi.mcmc.allocation != "Golgi")],
  c(3, 5)
]



# === Investigate cell cycle splitting/merging ======================================

model_output <- readRDS("./T_gondii/Output/processedModelOutputs.rds")

p_occ_mdi <- (1 * (model_output$MDI$N_k[[2]] != 0)) |>
  pheatmap(cluster_cols = TRUE, cluster_rows = FALSE, main = "MDI", silent = TRUE, color = simColPal())

p_occ_mix <- (1 * (model_output$Cell_cycle$N_k[[1]] != 0)) |>
  pheatmap(
    cluster_cols = TRUE, cluster_rows = FALSE, main = "Mixture model", silent = TRUE,
    color = simColPal()
  )

rowSums((1 * (model_output$MDI$N_k[[2]] != 0))) |> mean()
rowSums((1 * (model_output$Cell_cycle$N_k[[1]] != 0))) |> mean()

combinePheatmaps(list(p_occ_mdi$gtable, p_occ_mix$gtable), save_name = "~/Desktop/cellCycleOccupiedComponents.png", main = "Cell cycle data occupied components")

model_output$Cell_cycle$pred[[1]] |> table()
model_output$MDI$pred[[2]] |> table()

cl_14 <- which(pred_cl[[2]] == 14)

pred_cl[[1]][cl_14]
pred_cl[[2]][cl_14]

data_modelled$data_modelled[[2]] |>
  annotatedHeatmap()

cl_14_order <- order(pred_cl[[1]][cl_14])

cl_mix_in_14 <- model_output$Cell_cycle$pred[[1]] %in% c(4, 5, 7, 24, 32, 37, 57, 76, 94, 96, 102, 105, 109)

pred_cell_cycle <- data.frame("MDI" = c(model_output$MDI$pred[[2]]), "Mixture.model" = c(model_output$Cell_cycle$pred[[1]]))

# === DG clusters ====

dg_indices <- which(mdi_predictions == "dense granules")

dg_clusters <- as.numeric(names(which(table(pred_cl[[2]][dg_indices]) > 10)))

dg_cluster_members <- pred_cl[[2]] %in% dg_clusters

head(microarray_data)

dg_cell_cycle_data <- microarray_data[dg_cluster_members, ]
dg_cell_cycle_clusters <- pred_cl[[2]][dg_cluster_members]

# Create the annotation data.frame for the rows
anno_row <- data.frame(Cluster = factor(paste("Cluster", dg_cell_cycle_clusters)),
                       Localisation = mdi_predictions[dg_cluster_members])

cluster_IDs <- anno_row$Cluster

row.names(anno_row) <- row.names(dg_cell_cycle_data)

# The number of cololurs to use
K <- length(unique(cluster_IDs))

ann_colors <- list(Cluster = viridis::viridis(length(unique(dg_cell_cycle_clusters))),
                   Localisation = viridis::viridis(length(unique(mdi_predictions[dg_cluster_members]))))

names(ann_colors$Cluster) <- sort(unique(cluster_IDs))
names(ann_colors$Localisation) <- sort(unique(mdi_predictions[dg_cluster_members]))

col_pal <- dataColPal()
my_breaks <- defineDataBreaks(microarray_data[dg_cluster_members, ], col_pal, mid_point = 0)
pheatmap::pheatmap(microarray_data[dg_cluster_members, ], 
                   color = col_pal, 
                   breaks = my_breaks, 
                   annotation_row = anno_row,
                   annotation_colors = ann_colors,
                   show_rownames = FALSE, 
                   cluster_cols = FALSE
                   )

annotatedHeatmap(microarray_data[dg_indices, ], pred_cl[[2]][dg_indices],
                 cluster_cols = FALSE, 
                 main = "Dense granule protein transcription profiles",
                 show_rownames = FALSE #,
                 # filename = "~/Desktop/DGtranscriptomics.png"
)

dg_markers <- which(tagm_comparison$markers == "dense granules")

annotatedHeatmap(microarray_data[dg_markers, ], pred_cl[[2]][dg_markers],
  cluster_cols = FALSE, 
  main = "Dense granule marker protein transcription profiles",
  show_rownames = FALSE #,
  # filename = "~/Desktop/DGMarkersTranscriptomics.png"
)

raw_cell_cycle_data <- read.csv("./T_gondii/Original_data/ToxoDB_TgME49_Protein-coding_DNA_microarray.txt", sep = "\t", row.names = 1)

unnormalised_cell_cycle_data <- raw_cell_cycle_data[row.names(tagm_comparison), 4:16]

for(jj in seq(1, ncol(unnormalised_cell_cycle_data))) {
  unnormalised_cell_cycle_data[ , jj] <- as.numeric(unnormalised_cell_cycle_data[, jj])
}
unnormalised_cell_cycle_data <- scale(unnormalised_cell_cycle_data)

colnames(unnormalised_cell_cycle_data) <- colnames(microarray_data)

col_pal
dg_breaks <- defineDataBreaks(unnormalised_cell_cycle_data[dg_indices, ], col_pal)
annotatedHeatmap(unnormalised_cell_cycle_data[dg_markers, ], pred_cl[[2]][dg_markers],
                 my_breaks = dg_breaks,
                 cluster_cols = FALSE, 
                 main = "Dense granule marker protein absolute transcription profiles",
                 show_rownames = FALSE,
                 filename = "~/Desktop/DGMarkersAbsoluteTranscriptomics.png"
)

annotatedHeatmap(unnormalised_cell_cycle_data[dg_indices, ], pred_cl[[2]][dg_indices],
                 my_breaks = dg_breaks,
                 cluster_cols = FALSE, 
                 main = "Dense granule protein absolute transcription profiles",
                 show_rownames = FALSE,
                 filename = "~/Desktop/DGAbsoluteTranscriptomics.png"
)

dg_markers <- tagm_comparison$markers[dg_indices] == "dense granules"

heatmap_input <- list(
  coexpression_data = as.matrix(microarray_data[dg_indices, ]),
  abundance_data = as.matrix(unnormalised_cell_cycle_data[dg_indices, ]),
  labels = pred_cl[[2]][dg_indices],
  markers = dg_markers
)

p1 <- heatmapAnnotatedByClassAndBatch(heatmap_input)
p2 <- heatmapAnnotatedByClassAndBatch(heatmap_input, lead = "Abundance")

dg_larger_clusters <- as.numeric(names(which(table(pred_cl[[2]][dg_indices]) > 5)))
dg_large_cluster_indices <- which((mdi_predictions == "dense granules") & (pred_cl[[2]] %in% dg_larger_clusters))
dg_large_cluster_markers <- tagm_comparison$markers[dg_large_cluster_indices] == "dense granules"

heatmap_input_large_clusters <- list(
  coexpression_data = as.matrix(microarray_data[dg_large_cluster_indices, ]),
  abundance_data = as.matrix(unnormalised_cell_cycle_data[dg_large_cluster_indices, ]),
  labels = pred_cl[[2]][dg_large_cluster_indices],
  markers = dg_large_cluster_markers
)

p3 <- heatmapAnnotatedByClassAndBatch(heatmap_input3)
p4 <- heatmapAnnotatedByClassAndBatch(heatmap_input3, lead = "Abundance")

p1 + plot_annotation(title = "Dense granule proteins", subtitle = "Cell cycle transcriptomics data")
ggsave("~/Desktop/DG_proteins_cell_cycle_all.png", height = 8, width = 12)
p2 + plot_annotation(title = "Dense granule proteins", subtitle = "Cell cycle transcriptomics data")
ggsave("~/Desktop/DG_proteins_cell_cycle_all_ordered_by_abundance.png", height = 8, width = 12)
p3 + plot_annotation(title = "Subset of dense granule proteins", subtitle = "Cell cycle transcriptomics data")
ggsave("~/Desktop/DG_proteins_cell_cycle_reduced.png", height = 8, width = 12)
p4 + plot_annotation(title = "Subset of dense granule proteins", subtitle = "Cell cycle transcriptomics data")
ggsave("~/Desktop/DG_proteins_cell_cycle_reduced_ordered_by_abundance.png", height = 8, width = 12)


annoatedThreeLevels <- function(X, organelles, labels, markers, ...) {
  
  # Make sure the appropriate objects are factors
  organelles <- factor(organelles)
  labels <- factor(labels)
  markers <- factor(markers)
  
  # Create the annotation data.frame for the rows
  annotation_row <- data.frame(Localisation = organelles,
    Markers = markers,
    Cluster = labels
  )
  row.names(annotation_row) <- row.names(X)
  
  # The number of cololurs to use
  n_clusters <- length(levels(labels))
  n_organelles <- length(levels(organelles))
  
  annotation_colors <- list(
    Cluster = viridis::viridis(n_clusters),
    Markers = c("#FFFFFF", "#146EB4"),
    Localisation = viridis::viridis(n_organelles)
  )
  
  names(annotation_colors$Cluster) <- levels(labels)
  names(annotation_colors$Localisation) <- levels(organelles)
  names(annotation_colors$Markers) <- levels(markers)
  
  col_pal <- dataColPal()
  my_breaks <- defineDataBreaks(X, col_pal, mid_point = 0)
  pheatmap::pheatmap(X, 
                     color = col_pal, 
                     breaks = my_breaks, 
                     annotation_row = annotation_row,
                     annotation_colors = annotation_colors,
                     show_rownames = FALSE, 
                     cluster_cols = FALSE,
                     ...
  )
  
}
rhoptry_proteins <- which(mdi_predictions %in% c("rhoptries 1", "rhoptries 2"))
apical_proteins <- which(mdi_predictions %in% c("apical 1", "apical 2"))
apical_rhoptry_indices <- c(rhoptry_proteins, apical_proteins)

localisations <- mdi_predictions[apical_rhoptry_indices]
cluster_IDs <- pred_cl[[2]][apical_rhoptry_indices]
apical_rhoptry_markers <- tagm_comparison$markers[apical_rhoptry_indices] != "unknown"
apical_rhoptry_coexpression_data <- microarray_data[apical_rhoptry_indices, ]

annoatedThreeLevels(apical_rhoptry_coexpression_data, 
  localisations,
  cluster_IDs,
  apical_rhoptry_markers,
  main = "Apical and rhoptry protein co-expression data",
  filename = "~/Desktop/ApicalRhoprtryCoexpressionData.png"
)
# dev.off()

annotatedHeatmap(microarray_data[rhoptry_proteins, ], pred_cl[[2]][rhoptry_proteins])
annotatedHeatmap(microarray_data[apical_proteins, ], pred_cl[[2]][apical_proteins])
