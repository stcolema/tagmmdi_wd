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


prepCMssForGGplot <- function(cms, model_description_df,
                              matrix_setting_order = 1,
                              use_common_ordering = TRUE) {
  not_list <- !is.list(cms)
  if (not_list) {
    stop("`similarity_matrices` must be a list of matrices.")
  }
  
  n_matrices <- length(cms)
  n_models <- nrow(model_description_df)
  if (n_matrices != n_models) {
    .err <- paste(
      "Number of consensus matrices and number of models described",
      "in ``model_description_df`` are not matching."
    )
    stop(.err)
  }
  row_order <- col_order <- findOrder(cms[[matrix_setting_order]])
  
  depths <- model_description_df$Depth
  widths <- model_description_df$Width
  
  for (ii in seq(1, n_matrices)) {
    first_iteration <- ii == 1
    
    d <- depths[ii]
    w <- widths[ii]
    .df <- prepDataForggHeatmap(cms[[ii]], row_order, col_order)
    .df$Depth <- d
    .df$Width <- w
    
    if (first_iteration) {
      cm_df <- .df
    } else {
      cm_df <- rbind(cm_df, .df)
    }
  }
  cm_df$Depth <- factor(cm_df$Depth)
  cm_df$Width <- factor(cm_df$Width)
  cm_df
}

ccCalcAllocProbs <- function(allocation_probabilities, view,
                             burn = 0,
                             method = "median") {
  .alloc <- allocation_probabilities[[view]] # mcmc_samples$allocation_probabilities[[view]]
  
  probs <- NULL
  
  if (method == "median") {
    probs <- apply(.alloc, c(1, 2), median)
  }
  if (method == "mean") {
    probs <- rowSums(.alloc, dims = 2) / dim(.alloc)[3]
  }
  if (length(probs) == 1) {
    if (is.null(probs)) {
      stop("``method`` must be one of 'mean' or 'median'")
    }
  }
  probs
}

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

cat("\n=== Running ``processTGondiiCCoutput.R`` ===============================")

t0 <- Sys.time()
setMyTheme()

args <- input_arguments()

# Directories for input and output respectively
inputdata_dir <- "./T_gondii/Prepared_data/" #
save_dir <- "./T_gondii/ConsensusClustering/" #
model_output_dir <- "./T_gondii/ConsensusClustering/" #

plot_dir <- paste0(save_dir, "Plots/")
# dir.create(plot_dir, showWarnings = FALSE)

R <- 12000 #
K <- 125

cc_file <- paste0(save_dir, "CC_R_", R, "_K_", K, ".rds")
cc_file <- "~/Desktop/CC_R_12000_K_125_V_2.rds"
cc_out <- readRDS(cc_file)

cat("\n=== Reading in files ===================================================")

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

# microarray_data <- data_modelled$data_modelled[[1]] |> 
#   as_tibble()
# rna_seq_data <- data_modelled$data_modelled[[2]] |> 
#   as_tibble(rownames = "Gene")

prediction_mat <- do.call(cbind, cc_out$predicted_partitions)

# .alloc <- cc_out$allocations
# cc_out$allocations <- array(NA, c(50, 3643, 3))
# cc_out$allocations[, , 1] <- .alloc[[1]]
# cc_out$allocations[, , 2] <- .alloc[[2]]
# cc_out$allocations[, , 3] <- .alloc[[3]]

D <- nrow(cc_out$allocations)
fused_gene_probs_12 <- colMeans(cc_out$allocations[[1]] == cc_out$allocations[[2]] )
fused_genes_12 <- which(fused_gene_probs_12 > 0.5)

fused_gene_probs_13 <- colMeans(cc_out$allocations[[1]] == cc_out$allocations[[3]])
fused_genes_13 <- which(fused_gene_probs_13 > 0.5)

fused_gene_probs_23 <- colMeans(cc_out$allocations[[2]] == cc_out$allocations[[3]] )
fused_genes_23 <- which(fused_gene_probs_23 > 0.5)

prediction_mat[fused_genes_12, ]
c("nucelolus", "cytosol", "dense granules", "nucleus chromatim")

data("Barylyuk2020ToxoLopit")
cols_used <- c("Description", "markers", "tagm.mcmc.allocation","tagm.mcmc.probability")

proteins_modelled <- row.names(data_modelled$data_modelled[[1]])
tagm_comparison <- fData(Barylyuk2020ToxoLopit)[proteins_modelled, cols_used]

label_to_organelle <- data.frame("Organelle" = levels(tagm_comparison$markers)[-27],
                                 "Label"= seq(1, 26)
)

mdi_predictions <- label_to_organelle$Organelle[cc_out$predicted_partitions[[1]]]
mdi_probabilities <- cc_out$classification_probability

tagm_comparison <- tagm_comparison[proteins_modelled,  ]
# proteins_modelled <- which(row.names(tagm_comparison) %in% row.names(lopit_data))
# tagm_comparison <- tagm_comparison[proteins_modelled, ]
# tagm_comparison <- tagm_comparison[sort(row.names(tagm_comparison)), ]
tagm_comparison$mdi.mcmc.allocation <- mdi_predictions
tagm_comparison$mdi.mcmc.probability <- mdi_probabilities
tagm_comparison |> head()
lopit_data |> head()

tagm_comparison[fused_genes_12, ] |> head()

lopit_data$MDI_prediction <- mdi_predictions
lopit_data$TAGM_prediction <- tagm_comparison$tagm.mcmc.allocation
lopit_data |> head()

p_mdi <- lopit_data |> 
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |> 
  mutate(Fraction = factor(Fraction)) |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line() +
  facet_wrap(~MDI_prediction)

p_tagm <- lopit_data |> 
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |> 
  mutate(Fraction = factor(Fraction)) |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line() +
  facet_wrap(~TAGM_prediction)

lopit_plot_data <- lopit_data |> 
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |> 
  mutate(Fraction = factor(Fraction)) |> 
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |> 
  mutate(Model = case_when(Model == "MDI_prediction" ~ "MDI",
                           Model == "TAGM_prediction" ~ "TAGM")) 

lopit_plot_data |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 1)


lopit_plot_data |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 2)

lopit_plot_data |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 3)

lopit_plot_data |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 4)

lopit_plot_data |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 5)

lopit_plot_data |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 6)

lopit_data |> 
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |> 
  mutate(Fraction = factor(Fraction)) |> 
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) + 
  geom_line() +
  facet_grid(Organelle ~ Model)
  

lopit_data |> 
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |> 
  mutate(Fraction = factor(Fraction)) |> 
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, colour = Model)) + 
  geom_line(alpha = 0.3) +
  facet_wrap(~ Organelle) + 
  ggthemes::scale_color_colorblind()


tagm_comparison[which(tagm_comparison$tagm.mcmc.allocation != tagm_comparison$mdi.mcmc.allocation), c(1, 3, 5)] |> 

tagm_comparison[fused_genes_12, ]

p_tagm_mdi_contrast <- lopit_data |> 
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |> 
  mutate(Fraction = factor(Fraction)) |> 
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |> 
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, colour = Model)) + 
  geom_line(alpha = 0.3) +
  facet_wrap(~ Organelle) + 
  ggthemes::scale_color_colorblind()

ggsave("./tagm_mdi_contrast.png", p_tagm_mdi_contrast)

microarray_data |> 
  mutate(Predicted_label = factor(cc_out$predicted_partitions[[2]]),
         Gene = row.names(microarray_data)) |> 
  pivot_longer(-c(Predicted_label, Gene), names_to = "Time", values_to = "Expression") |> 
  mutate(Time = factor(Time, labels = seq(0, 12))) |> 
  # filter(Predicted_label %in% c(1:6, 9)) |> 
  ggplot(aes(x = Time, y = Expression, group = Gene)) + 
  geom_line(alpha = 0.3) +
  # ggthemes::scale_color_colorblind() + 
  facet_wrap(~Predicted_label)

rna_seq_data |>
  mutate(Predicted_label = factor(cc_out$predicted_partitions[[3]]),
         Gene = row.names(rna_seq_data)) |>
  pivot_longer(-c(Predicted_label, Gene), names_to = "Experiment", values_to = "Expression") |>
  # mutate(Time = factor(Time, labels = seq(0, 12))) |>
  # filter(Predicted_label %in% c(1:6, 9)) |>
  ggplot(aes(x = Experiment, y = Expression, group = Gene)) +
  geom_line(aes(color = Predicted_label)) +
  # ggthemes::scale_color_colorblind() +
  facet_wrap(~Predicted_label)

annotatedHeatmap(rna_seq_data, cc_out$predicted_partitions[[2]], 
                 show_colnames = FALSE, show_rownames = FALSE,
                 main = "RNA-seq data annotated by predicted cluster" # ,
                 # filename = "./rnaseq_predicted_clustering.png"
                 )


clusters_bigger_50 <- which(table(cc_out$predicted_partitions[[1]]) > 50)
genes_in_bigger_cl <- which(cc_out$predicted_partitions[[1]] %in% clusters_bigger_50)

annotatedHeatmap(microarray_data[genes_in_bigger_cl, ], 
                 cc_out$predicted_partitions[[1]][genes_in_bigger_cl], 
                      show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE,
                 main = "Cell cycle data annotated by predicted cluster (clusters with membership exceeding 50)",
                 filename = "./cell_cycle_predicted_clustering.png")

marker_genes <- which(tagm_comparison$markers != "unknown")
annotatedHeatmap(microarray_data[marker_genes, ], 
                 cc_out$predicted_partitions[[1]][marker_genes], 
                 show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE,
                 main = "Marker genes in cell cycle data")

# Create the annotation data.frame for the rows
anno_row <- data.frame(Cluster = factor(paste("Cluster", cc_out$predicted_partitions[[1]][marker_genes])),
                       Organelle = tagm_comparison$mdi.mcmc.allocation[marker_genes]) %>%
  magrittr::set_rownames(rownames(microarray_data)[marker_genes])

# The number of cololurs to use
# K <- 26 # length(unique(cluster_IDs))

# Create the annotation colours
ann_colours <- list(Cluster = viridis::viridis(19), Organelle = viridis::viridis(26))
names(ann_colours[[1]]) <- paste("Cluster", unique(cc_out$predicted_partitions[[1]][marker_genes]))
names(ann_colours[[2]]) <- unique(tagm_comparison$mdi.mcmc.allocation[marker_genes])
pheatmap::pheatmap(microarray_data[marker_genes, ],
                   color = dataColPal(),
                   breaks = defineDataBreaks(microarray_data[marker_genes, ], dataColPal(), mid_point = 0),
                   annotation_row = anno_row,
                   annotation_colors = ann_colours,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster",
                   filename = "./markger_genes_cell_cycle_data_comp.png"
)

organelles_of_interest <- c("micronemes", 
                            "apical 1",
                            "apical 2",
                            "rhoptries 1",
                            "rhoptries 2",
                            "dense granules")
proteins_in_organelles_of_interest <- which(tagm_comparison$mdi.mcmc.allocation %in% organelles_of_interest)
cell_cycle_clusters_plotted <- which(table(cc_out$predicted_partitions[[1]][proteins_in_organelles_of_interest]) > 15)
genes_of_interest <- which(cc_out$predicted_partitions[[1]] %in% cell_cycle_clusters_plotted)

genes_plotted <- proteins_in_organelles_of_interest[which(proteins_in_organelles_of_interest %in% genes_of_interest)]

table(tagm_comparison$mdi.mcmc.allocation[proteins_in_organelles_of_interest])


pheatmap::pheatmap(microarray_data[proteins_in_organelles_of_interest, ],
                   color = dataColPal(),
                   breaks = defineDataBreaks(microarray_data[proteins_in_organelles_of_interest, ], dataColPal(), mid_point = 0),
                   annotation_row = anno_row,
                   annotation_colors = ann_colours,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster",
                   filename = "./markger_genes_cell_cycle_data_comp.png"
)


annotatedHeatmap(microarray_data[genes_plotted, ],
                 tagm_comparison$mdi.mcmc.allocation[genes_plotted],
                 show_colnames = FALSE, show_rownames = FALSE)



# Create the annotation data.frame for the rows
anno_row <- data.frame(Cluster = factor(paste("Cluster", cc_out$predicted_partitions[[1]][genes_plotted])),
                       Organelle = tagm_comparison$mdi.mcmc.allocation[genes_plotted]) %>%
  magrittr::set_rownames(rownames(microarray_data)[genes_plotted])

# The number of cololurs to use
# K <- 26 # length(unique(cluster_IDs))

# Create the annotation colours
ann_colours <- list(Cluster = ggthemes::colorblind_pal()(7), Organelle = ggthemes::colorblind_pal()(6)) # viridis::viridis(6))
names(ann_colours[[1]]) <- paste("Cluster", unique(cc_out$predicted_partitions[[1]][genes_plotted]))
names(ann_colours[[2]]) <- unique(tagm_comparison$mdi.mcmc.allocation[genes_plotted])
pheatmap::pheatmap(microarray_data[genes_plotted, ],
                   color = dataColPal(),
                   breaks = defineDataBreaks(microarray_data[genes_plotted, ], dataColPal(), mid_point = 0),
                   annotation_row = anno_row,
                   annotation_colors = ann_colours,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster"
                   # filename = "./markger_genes_cell_cycle_data_comp.png"
)

tagm_comparison$Description[cc_out$predicted_partitions[[1]] == 3 & (tagm_comparison$Description != "hypothetical protein")]
hypo_prots <- which(tagm_comparison$Description == "hypothetical protein")

hypo_prots_in_genes_plotted <- genes_plotted[genes_plotted %in% hypo_prots]

ann_colours <- list(Cluster = ggthemes::colorblind_pal()(7), Organelle = ggthemes::colorblind_pal()(6)) # viridis::viridis(6))
names(ann_colours[[1]]) <- paste("Cluster", unique(cc_out$predicted_partitions[[1]][hypo_prots_in_genes_plotted]))
names(ann_colours[[2]]) <- unique(tagm_comparison$mdi.mcmc.allocation[hypo_prots_in_genes_plotted])
pheatmap::pheatmap(microarray_data[hypo_prots_in_genes_plotted, ],
                   color = dataColPal(),
                   breaks = defineDataBreaks(microarray_data[hypo_prots_in_genes_plotted, ], dataColPal(), mid_point = 0),
                   annotation_row = anno_row,
                   annotation_colors = ann_colours,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster"
                   # filename = "./markger_genes_cell_cycle_data_comp.png"
)

row_order <- findOrder(microarray_data)
micro_df <- prepDataForggHeatmap(as.matrix(microarray_data), row_order = row_order, col_order = seq(1, ncol(microarray_data)))
lopit_df <- prepDataForggHeatmap(scale(exprs(Barylyuk2020ToxoLopit)),
                                 row_order = row_order, 
                                 col_order = seq(1, ncol(exprs(Barylyuk2020ToxoLopit))))
rnaseq_df <- prepDataForggHeatmap(as.matrix(rna_seq_data), row_order = row_order)

lopit_df$Dataset <- "LOPIT"
micro_df$Dataset <- "Cell cycle"
rnaseq_df$Dataset <- "RNA seq"
heat_df <- rbind(micro_df, lopit_df, rnaseq_df)

heat_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() +
  facet_wrap(~Dataset, scales = "free") + 
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")

row_order <- findOrder(microarray_data[marker_genes,])
micro_df <- prepDataForggHeatmap(as.matrix(microarray_data[marker_genes,]), row_order = row_order, col_order = seq(1, ncol(microarray_data)))
lopit_df <- prepDataForggHeatmap(scale(exprs(Barylyuk2020ToxoLopit)[marker_genes, ]),
                                 row_order = row_order, 
                                 col_order = seq(1, ncol(exprs(Barylyuk2020ToxoLopit))))
rnaseq_df <- prepDataForggHeatmap(as.matrix(rna_seq_data)[marker_genes, ], row_order = row_order)

lopit_df$Dataset <- "LOPIT"
micro_df$Dataset <- "Cell cycle"
rnaseq_df$Dataset <- "RNA seq"
heat_df <- rbind(micro_df, lopit_df, rnaseq_df)

heat_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() +
  facet_wrap(~Dataset, scales = "free") + 
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")

microarray_data <- fread("./T_gondii/Original_data/ToxoDB_TgME49_Protein-coding_DNA_microarray.txt",
            na.strings = "N/A",
            strip.white = T,
            header = T,
            select = seq(1, 212)
)

rna_seq_data <- fread("./T_gondii/Original_data/ToxoDB_TgME49_Protein-coding_RNA-Seq.txt",
                      na.strings = "N/A",
                      strip.white = T,
                      header = T,
                      select = seq(1, 255)
)


m_white_cell_cycle
m_white_cell_cycle_normalised

tagm_comparison |> head()

marker_genes <- row.names(tagm_comparison)[tagm_comparison$markers != "unknown"]

pheatmap::pheatmap(m_white_cell_cycle_normalised[match(marker_genes, row.names(m_white_cell_cycle)),], show_colnames = FALSE)

m_white_cell_cycle_normalised_markers <- m_white_cell_cycle_normalised[match(marker_genes, row.names(m_white_cell_cycle)),]

m_white_cell_cycle_normalised_ordered <- m_white_cell_cycle_normalised[match(row.names(tagm_comparison), row.names(m_white_cell_cycle)),]

row_order <- findOrder(m_white_cell_cycle_normalised_ordered[ ,-1])

micro_df <- prepDataForggHeatmap(m_white_cell_cycle_normalised_ordered[, -1], row_order = row_order, col_order = seq(1, ncol(microarray_data)))
lopit_df <- prepDataForggHeatmap(scale(as.matrix(final_protein_df[, seq(1, 30)])),
                                 row_order = row_order, 
                                 col_order = seq(1, 30))
rnaseq_df <- prepDataForggHeatmap(as.matrix(normalised_rnaseq_macrophages_infected_by_T_gondii)[, ], row_order = row_order)

lopit_df$Dataset <- "LOPIT"
micro_df$Dataset <- "Cell cycle"
rnaseq_df$Dataset <- "RNA seq"
heat_df <- rbind(micro_df, lopit_df , rnaseq_df)

heat_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() +
  facet_wrap(~Dataset, scales = "free") + 
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")


normalised_rnaseq_macrophages_infected_by_T_gondii |> head()
