

suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(tagmReDraft))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))


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
      default = 50000,
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


setMyTheme()

args <- input_arguments()

# Directories for input and output respectively
inputdata_dir <- args$data_dir # "./T_gondii/Prepared_data/"
save_dir <- args$save_dir # "./" #
model_output_dir <- args$model_output_dir # "./T_gondii/Output/"

plot_dir <- paste0(save_dir, "Plots/")
dir.create(plot_dir, showWarnings = FALSE)

R <- args$R # 25000
K <- args$K # 125

burn <- args$burn
if(is.null(burn)) {
  burn <- floor(0.2 * R)
}

pattern <- paste0("(TGondiiMDI).*\\_K_", K, "_R_", R, ".rds$")
# pattern <- paste0("(TGondii_RNAseq).*\\_K_", K, "_R_", R, "_type_G.rds$")
files <- list.files(model_output_dir,
  pattern = pattern, # "(TGondiiMDI).*\\_K_125_R_25000.rds$",
  full.names = TRUE
)


plotting <- TRUE
plot_height <- 6
plot_width <- 8

n_files <- length(files)
mcmc_output <- vector("list", n_files)

cell_cycle_psms <- list()
rna_seq_psms <- list()
lopit_psms <- list()

# === Processing and input data ================================================

cat("\n# === Processing and input data =========================================")

for (ii in seq(1, n_files)) {
  .f <- files[ii]
  .x <- readRDS(.f)[[1]]
  mcmc_output[[ii]] <- .mcmc <- processMCMCChain(.x,
    burn = burn,
    point_estimate_method = "median"
  )

  mcmc_output[[ii]]$Chain <- ii

  applied_burn_in <- .mcmc$thin * floor(.mcmc$burn / .mcmc$thin)
  iterations <- seq(applied_burn_in + .mcmc$thin, .mcmc$R, .mcmc$thin)

  .phi_df <- .mcmc$phis %>%
    data.frame() %>%
    set_colnames(c("Phi_12", "Phi_13", "Phi_23")) %>%
    mutate(Chain = ii, Iteration = iterations)

  .alpha_df <- .mcmc$mass |>
    data.frame() %>%
    set_colnames(c("alpha_1", "alpha_2", "alpha_3")) %>%
    mutate(Chain = ii, Iteration = iterations)

  .evidence_df <- .mcmc$evidence |>
    data.frame() %>%
    set_colnames(c("Evidence")) %>%
    mutate(Chain = ii, Iteration = iterations)


  cell_cycle_psms[[ii]] <- mdiHelpR::makePSM(.mcmc$allocations[, , 1])
  rna_seq_psms[[ii]] <- mdiHelpR::makePSM(.mcmc$allocations[, , 2])
  lopit_psms[[ii]] <- mdiHelpR::makePSM(.mcmc$allocations[, , 3])

  if (ii == 1) {
    phi_df <- .phi_df
    alpha_df <- .alpha_df
    evidence_df <- .evidence_df
  } else {
    phi_df <- rbind(phi_df, .phi_df)
    alpha_df <- rbind(alpha_df, .alpha_df)
    evidence_df <- rbind(evidence_df, .evidence_df)
  }
}

phi_df$Chain <- factor(phi_df$Chain)
alpha_df$Chain <- factor(alpha_df$Chain)
evidence_df$Chain <- factor(evidence_df$Chain)

predictions <- predictFromMultipleChains(mcmc_output, burn = burn, chains_already_processed = TRUE)


allocations <- predictions$allocations
# str(allocations)

fusion_probs_12 <- colSums(allocations[[1]] == allocations[[2]]) / nrow(allocations[[1]])
fusion_probs_13 <- colSums(allocations[[1]] == allocations[[3]]) / nrow(allocations[[1]])
fusion_probs_23 <- colSums(allocations[[2]] == allocations[[3]]) / nrow(allocations[[1]])

fused_genes_12 <- which(fusion_probs_12 > 0.5)
fused_genes_13 <- which(fusion_probs_13 > 0.5)
fused_genes_23 <- which(fusion_probs_23 > 0.5)

fused_genes <- list(
  "microarray-rnaseq" = fused_genes_12,
  "microarray-lopit" = fused_genes_13,
  "rnaseq-lopit" = fused_genes_23
)



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

# === Plotting =================================================================

cat("\n# === Plotting ==========================================================")

cell_cycle_psm_df <- prepSimilarityMatricesForGGplot(cell_cycle_psms, ignore_checks = TRUE)
rna_seq_psm_df <- prepSimilarityMatricesForGGplot(rna_seq_psms, ignore_checks = TRUE)
lopit_psm_df <- prepSimilarityMatricesForGGplot(lopit_psms, ignore_checks = TRUE)

if (plotting) {
  cat("\nEvidence plots.")
  p_evidence_density <- evidence_df %>%
    ggplot(aes(x = Evidence, fill = Chain)) +
    geom_density(alpha = 0.3) +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Marginal likelihood across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )

  ggsave(paste0(plot_dir, "marginal_llikelihood_density.png"),
    plot = p_evidence_density,
    height = plot_height,
    width = plot_width
  )


  p_evidence_trace <- evidence_df %>%
    ggplot(aes(x = Iteration, y = Evidence, colour = Chain)) +
    geom_line() +
    ggthemes::scale_color_colorblind() +
    labs(
      title = "Marginal likelihood across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
   )
  
  ggsave(paste0(plot_dir, "marginal_llikelihood_trace.png"),
    plot = p_evidence_trace,
    height = plot_height,
    width = plot_width
  )

  cat("\nPhi plots.")

  p_phis <- phi_df %>%
    pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
    # filter(Value > ) %>%
    ggplot(aes(x = Value, fill = factor(Chain))) +
    geom_density(alpha = 0.3) +
    facet_wrap(~Variable, ncol = 1, scales = "free") +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Sampled phis across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )

  ggsave(paste0(plot_dir, "sampled_phis_density.png"),
    plot = p_phis,
    height = plot_height,
    width = plot_width
  )

  p_phis_trace <- phi_df %>%
    pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
    filter(Variable %in% c("Phi_12", "Phi_23")) %>%
    ggplot(aes(x = Iteration, y = Value, colour = Chain)) +
    geom_line() +
    facet_wrap(~Variable, ncol = 1, scales = "free") +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Sampled phis across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )

  ggsave(paste0(plot_dir, "sampled_phis_trace.png"),
    plot = p_phis_trace,
    height = plot_height,
    width = plot_width
  )

  cat("\nMass parameter.")
    p_alphas <- alpha_df %>%
    pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
    ggplot(aes(x = Value, fill = factor(Chain))) +
    geom_density(alpha = 0.3) +
    facet_wrap(~Variable, ncol = 1, scales = "free") +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Sampled mass parameter across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )


  p_alphas_trace <- alpha_df %>%
    pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
    ggplot(aes(x = Iteration, y = Value, colour = Chain)) +
    geom_line() +
    facet_wrap(~Variable, ncol = 1, scales = "free") +
    ggthemes::scale_color_colorblind() +
    labs(
      title = "Sampled mass parameter across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )

  ggsave(paste0(plot_dir, "sampled_alpha_density.png"),
    plot = p_alphas,
    height = plot_height,
    width = plot_width
  )

  ggsave(paste0(plot_dir, "sampled_alpha_trace.png"),
    plot = p_alphas_trace,
    height = plot_height,
    width = plot_width
  )


  cat("\nPSMs.")

  p_cc_psm <- cell_cycle_psm_df |>
    ggplot(aes(x = x, y = y, fill = Entry)) +
    geom_tile() +
    facet_wrap(~Chain) +
    scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
    labs(title = "Cell cycle PSMs across chains")

  p_rna_psm <- rna_seq_psm_df |>
    ggplot(aes(x = x, y = y, fill = Entry)) +
    geom_tile() +
    facet_wrap(~Chain) +
    scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
    labs(title = "CRNA-seq PSMs across chains")

  p_lopit_psm <- lopit_psm_df |>
    ggplot(aes(x = x, y = y, fill = Entry)) +
    geom_tile() +
    facet_wrap(~Chain) +
    scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
    labs(title = "LOPIT PSMs across chains")

  ggsave(paste0(plot_dir, "cell_cycle_psms.png"),
    plot = p_cc_psm,
    height = plot_height,
    width = plot_width
  )

  ggsave(paste0(plot_dir, "rna_psms.png"),
    plot = p_rna_psm,
    height = plot_height,
    width = plot_width
  )

  ggsave(paste0(plot_dir, "lopit_psms.png"),
    plot = p_lopit_psm,
    height = plot_height,
    width = plot_width
  )
  
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Value < 50) %>%
  #   ggplot(aes(x = Value, fill = factor(Chain))) +
  #   geom_density(alpha = 0.3) +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   ggthemes::scale_fill_colorblind() +
  #   labs(
  #     title = "Sampled phis across chains (less than 50)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  #
  # ggsave(paste0(plot_dir, "sampled_phis_across_chains_lt_50.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Value > 1) %>%
  #   ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  #   geom_density() +
  #   facet_grid(Variable ~ Chain, scales = "free") +
  #   labs(
  #     title = "Sampled phis across chains (greater than 1)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phis_across_chains_gt_1_grid.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Chain == 1) %>%
  #   ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  #   geom_density() +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   labs(
  #     title = "Sampled phis (chain 1)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phis_chain_1.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Chain == 2, Value > 1.0) %>%
  #   ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  #   geom_density() +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   labs(
  #     title = "Sampled phis (chain 2, values less than 1 dropped)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phis_chain_2_gt_1.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Chain == 3) %>%
  #   ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  #   geom_density() +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   labs(
  #     title = "Sampled phis (chain 3)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phis_chain_3.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Chain == 4) %>%
  #   ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  #   geom_density() +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   labs(
  #     title = "Sampled phis (chain 4)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phis_chain_4.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Chain == 5) %>%
  #   ggplot(aes(x = Value, group = Variable, fill = Variable)) +
  #   geom_density() +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   labs(
  #     title = "Sampled phis (chain 5)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phis_chain_5.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Variable == "Phi_12") %>%
  #   ggplot(aes(x = Value, y = stat(ndensity), group = Chain, fill = Chain)) +
  #   geom_density(alpha = 0.3) +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   ggthemes::scale_color_colorblind() +
  #   labs(
  #     title = "Sampled phi_12 (all chains)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phi_12_all_chains.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Variable == "Phi_13", Value < 80) %>%
  #   ggplot(aes(x = Value, y = stat(ndensity), group = Chain, fill = Chain)) +
  #   geom_density(alpha = 0.3) +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   ggthemes::scale_color_colorblind() +
  #   labs(
  #     title = "Sampled phi_13 (all chains)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phi_13_all_chains.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
  #
  # phi_df %>%
  #   pivot_longer(-c(Chain, Iteration), names_to = "Variable", values_to = "Value") %>%
  #   filter(Variable == "Phi_23", Value < 80) %>%
  #   ggplot(aes(x = Value, y = stat(ndensity), group = Chain, fill = Chain)) +
  #   geom_density(alpha = 0.3) +
  #   facet_wrap(~Variable, ncol = 1, scales = "free") +
  #   ggthemes::scale_color_colorblind() +
  #   labs(
  #     title = "Sampled phi_23 (all chains)",
  #     caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
  #   )
  #
  # ggsave(paste0(plot_dir, "sampled_phi_23_all_chains.png"),
  #   height = plot_height,
  #   width = plot_width
  # )
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

# list(
#   cell_cycle_data,
#   normalised_cell_cycle_data,
#   protein_mat
# )

# === Cell cycle data ==========================================================

cat("\n# === Investigate cell cycle data =======================================")

microarray_data |> pheatmap::pheatmap(show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE)

cell_cycle_psm <- predictions$allocations[[1]] |> makePSM()
cell_cycle_pred <- mcclust::maxpear(cell_cycle_psm, max.k = 125)

cell_cycle_psm_1 <- mcmc_output[[1]]$allocations[, , 1] |> makePSM()
cell_cycle_psm_2 <- mcmc_output[[2]]$allocations[, , 1] |> makePSM()
cell_cycle_psm_3 <- mcmc_output[[3]]$allocations[, , 1] |> makePSM()

point_estimate <- predictions$pred
point_estimate[[1]] <- cell_cycle_pred$cl

cell_cycle_pred_1 <- mcclust::maxpear(cell_cycle_psm_1, max.k = 30)$cl
cell_cycle_pred_2 <- mcclust::maxpear(cell_cycle_psm_2, max.k = 30)$cl
cell_cycle_pred_3 <- mcclust::maxpear(cell_cycle_psm_3, max.k = 30)$cl

used_inds1 <- which(cell_cycle_pred_1 < 7)
used_inds2 <- which(cell_cycle_pred_2 < 8)

annotatedHeatmap(microarray_data, cell_cycle_pred$cl, show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE, main = "Cell cycle data annotated by cluster")
annotatedHeatmap(microarray_data[used_inds1, ], cell_cycle_pred_1[used_inds1], show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE, main = "Cell cycle data annotated by cluster")
annotatedHeatmap(microarray_data[used_inds2, ], cell_cycle_pred_2[used_inds2], show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE, main = "Cell cycle data annotated by cluster")
annotatedHeatmap(microarray_data, cell_cycle_pred_3, show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE, main = "Cell cycle data annotated by cluster")


# Create the annotation data.frame for the rows
anno_row <- data.frame(
  Chain_1 = factor(paste("Cluster", cell_cycle_pred_1)),
  Chain_2 = factor(paste("Cluster", cell_cycle_pred_2))
) %>%
  magrittr::set_rownames(rownames(microarray_data))

# The number of cololurs to use
K <- c(length(unique(cell_cycle_pred_1)), length(unique(cell_cycle_pred_2)))

# Create the annotation colours
ann_colours <- list(Chain_1 = viridis::viridis(K[1]), Chain_2 = viridis::viridis(K[2]))

names(ann_colours[[1]]) <- paste("Cluster", sort(unique(cell_cycle_pred_1)))
names(ann_colours[[2]]) <- paste("Cluster", sort(unique(cell_cycle_pred_2)))

col_pal <- dataColPal()
my_breaks <- defineDataBreaks(microarray_data, col_pal)
pheatmap::pheatmap(microarray_data,
  color = col_pal,
  breaks = my_breaks,
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_colnames = FALSE, show_rownames = FALSE
)


annotatedHeatmap(microarray_data, cell_cycle_pred$cl, show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE, main = "Cell cycle data annotated by cluster")
# pheatmap::pheatmap(cell_cycle_psm, breaks = defineBreaks(simColPal(), lb = 0), color = simColPal())

timepoints <- paste0(seq(0, 12), "HR")
long_cell_cycle_df <- microarray_data |>
  set_colnames(timepoints) |>
  mutate(Predicted_label = factor(cell_cycle_pred$cl), Gene = row.names(microarray_data)) |>
  pivot_longer(-c(
    Gene,
    Predicted_label
  ), names_to = "Timepoint", values_to = "Value") |>
  mutate(Timepoint = factor(Timepoint, levels = timepoints, ordered = TRUE))

long_cell_cycle_df %>%
  ggplot(aes(x = Timepoint, y = Value, colour = Predicted_label, group = Gene)) +
  geom_line() +
  facet_wrap(~Predicted_label)

# === LOPIT data ===============================================================

data_df <- lopit_data
data_df$Predicted_label <- point_estimate[[3]]

protein_lst <- prepareMSObject(Barylyuk2020ToxoLopit)
measurements <- colnames(Barylyuk2020ToxoLopit)

long_data_df <- data_df |> # [fused_genes_23, ] %>%
  filter(Predicted_label %in% c(24, 25)) |>
  pivot_longer(-c(
    Protein,
    Label,
    Fixed,
    Predicted_label
  ), names_to = "Fraction", values_to = "Value") %>%
  mutate(
    Fraction = factor(Fraction, levels = measurements, ordered = TRUE),
    Label = factor(Label,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    ),
    Predicted_label = factor(Predicted_label,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    )
  )


long_data_df %>%
  # filter(Fixed == 0) |>
  ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
  geom_line() +
  # facet_wrap(~Predicted_label)
  facet_grid(Fixed ~ Predicted_label)

# === Merged proteins ==========================================================


tagm_pred <- fData(Barylyuk2020ToxoLopit)$tagm.mcmc.allocation.pred
pred_order <- match(row.names(lopit_data), row.names(fData(Barylyuk2020ToxoLopit)))

mdi_pred <- factor(data_df$Predicted_label,
  levels = c(protein_lst$class_key$Key, 27),
  labels = c(protein_lst$class_key$Organelle, "unknown")
)

names(mdi_pred) <- data_df$Protein

tagm_pred <- tagm_pred[pred_order]
fixed <- data_df$Fixed
test_inds <- which(fixed == 1)

# Quick check that things are not completely wrong
all(mdi_pred[test_inds] == tagm_pred[test_inds])
disagreeing_pred <- which(mdi_pred != tagm_pred)

mdi_pred[disagreeing_pred] |> table()
tagm_pred[disagreeing_pred] |> table()

novel_rhoptries_1 <- which(mdi_pred[disagreeing_pred] == "rhoptries 1")
novel_rhoptries_2 <- which(mdi_pred[disagreeing_pred] == "rhoptries 2")

tagm_pred[disagreeing_pred][novel_rhoptries_2]
tagm_pred[disagreeing_pred][novel_rhoptries_1]

long_data_df <- data_df |> # [fused_genes_23, ] %>%
  filter(Predicted_label %in% c(24, 25)) |>
  pivot_longer(-c(
    Protein,
    Label,
    Fixed,
    Predicted_label
  ), names_to = "Fraction", values_to = "Value") %>%
  mutate(
    Fraction = factor(Fraction, levels = measurements, ordered = TRUE),
    Label = factor(Label,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    ),
    Predicted_label = factor(Predicted_label,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    )
  )


long_data_df %>%
  # filter(Fixed == 0) |>
  ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
  geom_line() +
  # facet_wrap(~Predicted_label)
  facet_grid(Fixed ~ Predicted_label)

# fused_microarray_data <- microarray_mat[fused_genes, c(1:27, 124:140)]
fused_microarray_data <- microarray_data[fused_genes_13, ]
fused_rna_seq_data <- rna_seq_data[fused_genes_23, ]

pheatmap::pheatmap(fused_microarray_data[order(data_df$Predicted_label[fused_genes_13]), ],
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = FALSE, show_colnames = FALSE
)

annotation_row <- factor(data_df$Predicted_label,
  levels = protein_lst$class_key$Key,
  labels = protein_lst$class_key$Organelle
) %>%
  set_names(row.names(data_df))

annotatedHeatmap(fused_microarray_data[order(data_df$Predicted_label[fused_genes_13]), ],
  cluster_IDs = annotation_row[fused_genes_13][order(data_df$Predicted_label[fused_genes_13])],
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Fused genes across cell cycle microarray and LOPIT data",
  filename = paste0(save_dir, "cellCycleDataFusedGenes.png")
)

pheatmap::pheatmap(fused_rna_seq_data[order(data_df$Label[fused_genes_23]), ],
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = FALSE, show_colnames = FALSE
)



pheatmap::pheatmap(fused_normalised_cell_cycle_data[order(data_df$Label[fused_genes_12]), ],
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = FALSE, show_colnames = FALSE
)


protein_mat <- as.matrix(final_protein_df[, -c(31:33)])

# data_modelled <- list(
#   microarray_mat[, c(1:27, 124:140)],
#   rna_mat,
#   protein_mat
# )

saveRDS(fused_genes, file = paste0(save_dir, "fusedGenes.rds"))
saveRDS(predictions, file = paste0(save_dir, "MDIpredictions.rds"))
