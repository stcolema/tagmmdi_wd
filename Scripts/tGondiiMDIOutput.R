
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

setMyTheme()

#' @title Prepare MS Object
#' @description Prepares a mass spectrometry experiment (as stored in
#' the Bioconductor package ``pRolocdata``) for modelling, extracting the
#' numerical data, the classes and the indicator matrix for which labels are
#' observed.
#' @param MS_object A mass spectrometry experiment such as ``tan2009r1`` from
#' ``pRolocdata``.
#' @return A list of ``X``, the fracitonation data from a LOPIT experiment,
#'  ``fixed``, the matrix indicating which labels are observed,
#'  ``initial_labels``, the matrix of the initial labels to be input into
#'  ``runMCMCChains`` (note that ``"unknown"`` organelles are represented
#'  arbitrarily with a ``1`` as these will be sampled again within the wrapper
#'  of ``callMDI``) and ``class_key`` which maps the numeric representation of
#'  organelles back to the original naming.
#' @export
prepareMSObject <- function(MS_object) {
  
  # Extract the LOPIT data and the organelles
  X <- Biobase::exprs(MS_object)
  organelles <- fData(MS_object)[, "markers"]
  
  # Create a data frame of the classes present and their associated number;\
  # this can be used to map the numeric representation of the classes back to
  # an organelle
  organelles_present <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(
    Organelle = organelles_present,
    Key = 1:length(organelles_present)
  )
  
  # Number of components modelled
  K <- length(organelles_present)
  
  # Number of views modelled
  V <- 1
  
  # Number of samples modelled
  N <- nrow(X)
  
  # Prepare initial labels
  initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)
  
  # Fix training points, allow test points to move component
  fix_vec <- (organelles != "unknown") * 1
  fixed[, 1] <- fix_vec
  
  # Assign initial labels
  initial_labels[, 1] <- class_key$Key[match(organelles, class_key$Organelle)]
  
  # Any unknown labels are given an arbitrary value which will be reset in the
  # model call function.
  initial_labels[is.na(initial_labels)] <- 1
  
  data_modelled <- list(
    X
  )
  
  # Return the prepared objects
  list(
    X = X,
    fixed = fixed,
    initial_labels = initial_labels,
    class_key = class_key
  )
}


library(data.table)
library(tibble)
library(tagmReDraft)
library(dplyr)
library(tidyr)
library(ggplot2)

plot_dir <- "./T_gondii/Output/Plots/"
dir.create(plot_dir, showWarnings = FALSE)


# Directories for input and output respectively
inputdata_dir <- "./T_gondii/Prepared_data/" # args$data_dir
save_dir <- "./" # args$save_dir
data_dir <- "./T_gondii/Output/"

files <- list.files(data_dir, 
                    pattern="(TGondiiMDI).*\\_K_125_R_65000.rds$",
                    full.names = TRUE)

plotting <- FALSE
plot_height = 6
plot_width = 8

n_files <- length(files)
mcmc_output <- vector("list", n_files)
burn <- 40000

for(ii in seq(1, n_files)) {
  .f <- files[ii]
  .x <- readRDS(.f)[[1]]
  mcmc_output[[ii]] <- .mcmc <- processMCMCChain(.x, 
                                        burn = burn, 
                                        point_estimate_method = "median")
  
  mcmc_output[[ii]]$Chain <- ii
  
  applied_burn_in <- .mcmc$thin * floor(.mcmc$burn / .mcmc$thin)
  iterations <- seq(applied_burn_in, .mcmc$R, .mcmc$thin)
  
  .phi_df <- .mcmc$phis %>% 
    data.frame() %>% 
    set_colnames(c("Phi_12", "Phi_13", "Phi_23")) %>% 
    mutate(Chain = ii, Iteration = iterations)

  if(ii == 1) {
    phi_df <- .phi_df
  } else {
    phi_df <- rbind(phi_df, .phi_df)
  }
}

phi_df$Chain <- factor(phi_df$Chain)

if(plotting) {
phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  # filter(Value > ) %>% 
  ggplot(aes(x = Value, fill = factor(Chain))) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  ggthemes::scale_fill_colorblind() + 
  labs(title = "Sampled phis across chains",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_across_chains_full.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Variable %in% c("Phi_12", "Phi_23") ) %>%
  ggplot(aes(x = Iteration, y = Value, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  ggthemes::scale_fill_colorblind() + 
  labs(title = "Sampled phis across chains",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_across_chains_trace.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Value < 50 ) %>%
  ggplot(aes(x = Value, fill = factor(Chain))) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  ggthemes::scale_fill_colorblind() + 
  labs(title = "Sampled phis across chains (less than 50)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")


ggsave(paste0(plot_dir, "sampled_phis_across_chains_lt_50.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Value > 1) %>% 
  ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  geom_density() +
  facet_grid(Variable~Chain, scales = "free") + 
  labs(title = "Sampled phis across chains (greater than 1)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_across_chains_gt_1_grid.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Chain == 1) %>% 
  ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  geom_density() +
  facet_wrap(~Variable, ncol = 1, scales = "free") + 
  labs(title = "Sampled phis (chain 1)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_chain_1.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Chain == 2, Value > 1.0) %>% 
  ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  geom_density() +
  facet_wrap(~Variable, ncol = 1, scales = "free") + 
  labs(title = "Sampled phis (chain 2, values less than 1 dropped)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_chain_2_gt_1.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Chain == 3) %>% 
  ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  geom_density() +
  facet_wrap(~Variable, ncol = 1, scales = "free") + 
  labs(title = "Sampled phis (chain 3)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_chain_3.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Chain == 4) %>% 
  ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  geom_density() +
  facet_wrap(~Variable, ncol = 1, scales = "free") + 
  labs(title = "Sampled phis (chain 4)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_chain_4.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Chain == 5) %>% 
  ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  geom_density() +
  facet_wrap(~Variable, ncol = 1, scales = "free") + 
  labs(title = "Sampled phis (chain 5)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phis_chain_5.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Variable == "Phi_12") %>% 
  ggplot(aes(x = Value, y = stat(ndensity), group = Chain, fill = Chain)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  ggthemes::scale_color_colorblind() + 
  labs(title = "Sampled phi_12 (all chains)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phi_12_all_chains.png"),
       height = plot_height,
       width = plot_width)


phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Variable == "Phi_13", Value < 80) %>% 
  ggplot(aes(x = Value, y = stat(ndensity), group = Chain, fill = Chain)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  ggthemes::scale_color_colorblind() + 
  labs(title = "Sampled phi_13 (all chains)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phi_13_all_chains.png"),
       height = plot_height,
       width = plot_width)

phi_df %>% 
  pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>% 
  filter(Variable == "Phi_23", Value < 80) %>% 
  ggplot(aes(x = Value, y = stat(ndensity), group = Chain, fill = Chain)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  ggthemes::scale_color_colorblind() + 
  labs(title = "Sampled phi_23 (all chains)",
       caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data")

ggsave(paste0(plot_dir, "sampled_phi_23_all_chains.png"),
       height = plot_height,
       width = plot_width)
}

# mcmc_output <- mcmc_output[[1]] # x1[[1]]
# 
# R <- mcmc_output$R
# thin <- mcmc_output$thin
# iter <- seq(1, R + 1, thin) - 1
# phi_df <- mcmc_output$phis %>% 
#   as.data.frame() %>% 
#   magrittr::set_colnames(c("Phi_12", "Phi_13", "Phi_23")) %>% 
#   mutate(Iteration = iter) %>% 
#   pivot_longer(-Iteration, values_to = "Sampled_value", names_to = "Phi")
# 
# phi_df %>% 
#   filter(Iteration > 5000) %>% 
#   ggplot(aes(x = Sampled_value)) +
#   geom_histogram() +
#   facet_wrap(~Phi, scales = "free_x")
# 
# phi_df %>% 
#   filter(Iteration > 5000, Sampled_value < 2000) %>% 
#   ggplot(aes(x = Sampled_value)) +
#   geom_histogram() +
#   facet_wrap(~Phi)
# 
# eff_r <- R / thin

predictions <- predictFromMultipleChains(mcmc_output, burn = burn, chains_already_processed = TRUE)


allocations <- predictions$allocations
# str(allocations)

fusion_probs_12 <- colSums(allocations[[1]] == allocations[[2]]) / nrow(allocations[[1]])
fusion_probs_13 <- colSums(allocations[[1]] == allocations[[3]]) / nrow(allocations[[1]])
fusion_probs_23 <- colSums(allocations[[2]] == allocations[[3]]) / nrow(allocations[[1]])

which(fusion_probs_12 > 0.5)
which(fusion_probs_13 > 0.5)
which(fusion_probs_23 > 0.5)

fused_genes_12 <- which(fusion_probs_12 > 0.5)
fused_genes_13 <- which(fusion_probs_13 > 0.5)
fused_genes_23 <- which(fusion_probs_23 > 0.5)

fused_genes <- list(
  "microarray-rnaseq" = fused_genes_12,
  "microarray-lopit" = fused_genes_13,
  "rnaseq-lopit" = fused_genes_23
)



microarray_file <- paste0(inputdata_dir, "cellCycleNormalised.csv")  # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
rna_seq_file <- paste0(inputdata_dir, "rnaSeqMacrophage.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_RNA-Seq.txt")
lopit_file <- paste0(inputdata_dir, "LOPITreduced.csv")
data(Barylyuk2020ToxoLopit)

microarray_data <- read.csv(microarray_file, row.names = 1, 
                         # na.strings = "N/A",
                         # strip.white = T,
                         header = T
                         # select = seq(1, 212)
)

rna_seq_data <- read.csv(rna_seq_file, row.names = 1, 
                      # na.strings = "N/A",
                      # strip.white = T,
                      header = T
                      # select = seq(1, 255)
)

lopit_data <- read.csv(lopit_file, row.names = 1, 
                         # na.strings = "N/A",
                         # strip.white = T,
                         header = T
                         # select = seq(1, 255)
)

data_modelled <- readRDS(paste0(data_dir, "TGondiiMDI_K_125_input.rds"))

# list(
#   cell_cycle_data,
#   normalised_cell_cycle_data,
#   protein_mat
# )

point_estimate <- predictions$pred

data_df <- lopit_data
data_df$Predicted_label <- point_estimate[[3]]

protein_lst <- prepareMSObject(Barylyuk2020ToxoLopit)
measurements <- colnames(Barylyuk2020ToxoLopit)

long_data_df <- data_df[fused_genes_13, ] %>%
  pivot_longer(-c(
    Protein,
    Label,
    Fixed,
    Predicted_label
  ), names_to = "Fraction", values_to = "Value") %>%
  mutate(Fraction = factor(Fraction, levels = measurements, ordered = TRUE),
         Label = factor(Label, 
                        levels = protein_lst$class_key$Key, 
                        labels = protein_lst$class_key$Organelle),
         Predicted_label = factor(Predicted_label, 
                        levels = protein_lst$class_key$Key, 
                        labels = protein_lst$class_key$Organelle))


long_data_df %>% 
  ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
  geom_line() +
  facet_grid(Fixed~Predicted_label)

# fused_microarray_data <- microarray_mat[fused_genes, c(1:27, 124:140)]
fused_microarray_data <- microarray_data[fused_genes_13, ]
fused_normalised_cell_cycle_data <- normalised_cell_cycle_data[fused_genes_12, ]

pheatmap::pheatmap(fused_microarray_data[order(data_df$Predicted_label[fused_genes_13]), ], 
                   cluster_rows = F, 
                   cluster_cols = F,
                   show_rownames = FALSE, show_colnames = FALSE)

annotation_row <- factor(data_df$Predicted_label, 
                         levels = protein_lst$class_key$Key, 
                         labels = protein_lst$class_key$Organelle) %>% 
  set_names(row.names(data_df))
annotatedHeatmap(fused_microarray_data[order(data_df$Predicted_label[fused_genes_13]), ], 
                 cluster_IDs = annotation_row[fused_genes_13][order(data_df$Predicted_label[fused_genes_13])], 
                 show_colnames = FALSE, 
                 show_rownames = FALSE,
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE,
                 main = "Fused genes across cell cycle microarray and LOPIT data",
                 filename = paste0(save_dir, "cellCycleDataFusedGenes.png"))

pheatmap::pheatmap(fused_cell_cycle_data[order(data_df$Label[fused_genes_13]), ], 
                   cluster_rows = F, 
                   cluster_cols = F,
                   show_rownames = FALSE, show_colnames = FALSE)



pheatmap::pheatmap(fused_normalised_cell_cycle_data[order(data_df$Label[fused_genes_12]), ], 
                   cluster_rows = F, 
                   cluster_cols = F,
                   show_rownames = FALSE, show_colnames = FALSE)


protein_mat <- as.matrix(final_protein_df[, -c(31:33)])

# data_modelled <- list(
#   microarray_mat[, c(1:27, 124:140)],
#   rna_mat,
#   protein_mat
# )

saveRDS(fused_genes, file = paste0(save_dir, "fusedGenes.rds"))
saveRDS(predictions, file = paste0(save_dir, "MDIpredictions.rds"))

