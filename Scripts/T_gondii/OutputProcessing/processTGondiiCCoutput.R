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
    optparse::make_option(c("-r", "--R"),
      type = "numeric",
      default = 15000L,
      help = "Number of iterations to run in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-t", "--thin"),
      type = "numeric",
      default = 1000L,
      help = "Number of iterations to thin each MCMC chain by [default= %default]",
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
      default = 125L,
      help = "Number of components modelled in the unsupervised views. [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-v", "--V"),
      type = "numeric",
      default = 3L,
      help = "Number of views modelled.",
      metavar = "numeric"
    ),
    optparse::make_option(c("-p", "--plotting"),
      type = "logical",
      default = TRUE,
      help = "Flag indicating if plotting should be performed. [default= %default]",
      metavar = "logical"
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
inputdata_dir <- args$data_dir # "./T_gondii/Prepared_data/" #
save_dir <- args$save_dir # "./" #
model_output_dir <- args$model_output_dir # "./T_gondii/ConsensusClustering/" #

plot_dir <- paste0(save_dir, "Plots/")
dir.create(plot_dir, showWarnings = FALSE)

R <- args$R # 5000 #
K <- args$K # 125
thin <- args$thin

V <- args$V
pattern <- paste0("(TGondiiMDI).*\\_K_", K, "_R_", R, "_V_", V, ".rds$")
# pattern <- paste0("(TGondii_RNAseq).*\\_K_", K, "_R_", R, "_type_G.rds$")

output_file <- paste0(save_dir, "CC_R_", R, "_K_", K, "_V_", V, ".rds")

cat("\n=== Reading in files ===================================================")

files <- list.files(model_output_dir,
  pattern = pattern, # "(TGondiiMDI).*\\_K_125_R_25000.rds$",
  full.names = TRUE
)

plotting <- args$plotting
plot_height <- 6
plot_width <- 8

n_files <- length(files)
mcmc_output <- vector("list", n_files)
burn <- 0 # floor(0.2 * R)

cell_cycle_psms <- list()
rna_seq_psms <- list()
lopit_psms <- list()

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

# === Processing and input data ================================================
cat("\n# === Processing and input data =======================================")

N <- nrow(microarray_data)
# V <- length(data_modelled$data_modelled)
view_inds <- seq(1, V)
lopit_ind <- 1 # which(data_modelled$types == "TAGM")

D_considered <- seq(3000, 15000, 3000) # seq(thin, R, thin * 2)
# D_considered <- c(1000, 3000, 5000, 8000, 12000)
D <- max(D_considered)
number_depths <- length(D_considered)

W_considered <- c(25, 50, 75, 100, 125, 150)
# W_considered <- c(10, 20, 40, 60, 80, 100)
# W_considered <- c(floor(n_files / 4), floor(n_files / 2), n_files)
W <- max(W_considered)
number_chains <- length(W_considered)

models <- expand.grid(D_considered, W_considered)
colnames(models) <- c("Depth", "Width")
n_models <- nrow(models)

iterations <- seq(0, D, thin)

allocation_list <- consensus_chains <- vector("list", W)

consensus_matrix <- matrix(0, nrow = N, ncol = N)

consensus_matrices <- vector("list", n_models)

allocations <- allocation_probs <- vector("list", length(D_considered))

for (ii in seq(1, number_depths)) {
  allocations[[ii]] <- array(0, dim = c(W, N, V))
  allocation_probs[[ii]] <- vector("list", V)
  for (v in view_inds) {
    K_v <- data_modelled$K[v]
    allocation_probs[[ii]][[v]] <- array(0, dim = c(N, K_v, W))
  }
}

cat("\nNumber of components in each view:\n", data_modelled$K)

# account for saving of 0th iterations
samples_extracted <- (D_considered / thin) + 1

n_phis <- choose(V, 2)
phi_names <- c("Phi_12", "Phi_13", "Phi_23")[seq(1, n_phis)]

for (ii in seq(1, n_files)) {
  .f <- files[ii]
  .mcmc <- readRDS(.f)[[1]]

  for (jj in seq(1, number_depths)) {
    curr_d <- D_considered[jj]
    sample_used <- samples_extracted[jj]

    .prob <- .mcmc$allocation_probabilities
    .alloc <- .mcmc$allocations[sample_used, , ]
    .ll <- .mcmc$complete_likelihood[sample_used]
    .phis <- .mcmc$phis[sample_used, ]
    .mass <- .mcmc$mass[sample_used, ]
    .evidence <- .mcmc$evidence[sample_used - 1]

    .fit_entry <- data.frame(
      "Complete_log_likelihood" = .ll,
      "Evidence" = .evidence,
      "Depth" = curr_d,
      "Chain" = ii
    )

    .phi_df <- .phis %>%
      t() %>%
      data.frame() %>%
      set_colnames(phi_names) %>%
      mutate(Chain = ii, Depth = curr_d)

    .mass_entry <- as.data.frame(t(.mass)) %>%
      set_colnames(c(paste0("Mass_", view_inds))) %>%
      mutate("Depth" = curr_d, "Chain" = ii)

    for (v in view_inds) {
      allocations[[jj]][ii, , v] <- .alloc[, v]
      allocation_probs[[jj]][[v]][, , ii] <- .prob[[v]][, , sample_used]
    }

    if (ii == 1) {
      # allocations[[jj]] <- .alloc
      fit_df <- .fit_entry
      # lkl_df <- .lkl_entry
      mass_df <- .mass_entry
      phi_df <- .phi_df
    } else {

      # allocations[[jj]] <- rbind(allocations[[jj]], .alloc)
      # lkl_df <- rbind(lkl_df, .lkl_entry)
      fit_df <- rbind(fit_df, .fit_entry)
      mass_df <- rbind(mass_df, .mass_entry)
      phi_df <- rbind(phi_df, .phi_df)
    }
  }
}

fit_df$Chain <- factor(fit_df$Chain)
phi_df$Chain <- factor(phi_df$Chain)
mass_df$Chain <- factor(mass_df$Chain)

fit_df$Depth <- factor(fit_df$Depth)
phi_df$Depth <- factor(phi_df$Depth)
mass_df$Depth <- factor(mass_df$Depth)

cat("\nMaking consensus matrices.")

# Now the assessing of consensus clustering
cms <- vector("list", V)
predicted_partitions <- vector("list", V)
classification_probability <- vector("list", n_models)
mean_absolute_difference <- vector("list", V)
mean_absolute_difference_df <- NULL
for (v in view_inds) {
  cms[[v]] <- vector("list", n_models)
  predicted_partitions[[v]] <- vector("list", n_models)
}

# Make CMs and compare across depths
for (ii in seq(1, number_chains)) {
  curr_w <- W_considered[ii]
  chains_used <- seq(1, curr_w)

  for (jj in seq(1, number_depths)) {
    curr_d <- D_considered[jj]

    curr_cm_index <- which((models$Depth == curr_d) & (models$Width == curr_w)) # (ii - 1) * number_depths + jj
    for (v in view_inds) {
      .cm <- createSimilarityMat(matrix(allocations[[jj]][chains_used, , v], ncol = N))
      row.names(.cm) <- colnames(.cm) <- row.names(microarray_data)
      cms[[v]][[curr_cm_index]] <- .cm

      if (v != lopit_ind) {
        predicted_partitions[[v]][[curr_cm_index]] <- matrix(allocations[[jj]][chains_used, , v], ncol = N)
        # predicted_partitions[[v]][[curr_cm_index]] <- mcclust::maxpear(.cm)$cl
      } else {
        .alloc_prob <- ccCalcAllocProbs(allocation_probs[[jj]], view = lopit_ind)
        classification_probability[[curr_cm_index]] <- apply(.alloc_prob, 1, max)
        predicted_partitions[[v]][[curr_cm_index]] <- apply(.alloc_prob, 1, which.max)
      }

      if (jj > 1) {
        comparison_d <- D_considered[jj - 1]
        comparison_cm_index <- which((models$Depth == comparison_d) & (models$Width == curr_w))

        # mean_absolute_difference[[v]] <-
        mad <- mean(abs(cms[[v]][[curr_cm_index]] - cms[[v]][[comparison_cm_index]]))
        mad_entry <- data.frame(
          "Depth" = curr_d,
          "Width" = curr_w,
          "Mean_absolute_difference" = mad,
          "Depth_compared" = comparison_d,
          "Width_compared" = curr_w,
          "View" = v,
          "Dataset" = datasets[v],
          "Quantity_varied" = "Depth"
        )

        if (is.null(mean_absolute_difference_df)) {
          mean_absolute_difference_df <- mad_entry
        } else {
          mean_absolute_difference_df <- rbind(mean_absolute_difference_df, mad_entry)
        }
      }
    }
  }
}

# Compare across widths
for (jj in seq(1, number_depths)) {
  curr_d <- D_considered[jj]
  for (ii in seq(2, number_chains)) {
    curr_w <- W_considered[ii]
    chains_used <- seq(1, curr_w)
    curr_cm_index <- which((models$Depth == curr_d & models$Width == curr_w)) # (ii - 1) * number_depths + jj
    comparison_w <- W_considered[ii - 1]
    comparison_cm_index <- which((models$Depth == curr_d) & (models$Width == comparison_w))
    for (v in view_inds) {

      # mean_absolute_difference[[v]] <-
      mad <- mean(abs(cms[[v]][[curr_cm_index]] - cms[[v]][[comparison_cm_index]]))
      mad_entry <- data.frame(
        "Depth" = curr_d,
        "Width" = curr_w,
        "Mean_absolute_difference" = mad,
        "Depth_compared" = curr_d,
        "Width_compared" = comparison_w,
        "View" = v,
        "Dataset" = datasets[v],
        "Quantity_varied" = "Width"
      )

      mean_absolute_difference_df <- rbind(mean_absolute_difference_df, mad_entry)
    }
  }
}

# === Plotting =================================================================
cat("\n === Plotting =========================================================")

cell_cycle_cm_df <- prepCMssForGGplot(cms[[2]], models, matrix_setting_order = n_models)
if (V == 3) {
  rna_seq_cm_df <- prepCMssForGGplot(cms[[3]], models, matrix_setting_order = n_models)
}
lopit_cm_df <- prepCMssForGGplot(cms[[1]], models, matrix_setting_order = n_models)

cat("\nCMs prepared.")

if (plotting) {
  p_phis <- phi_df %>%
    pivot_longer(-c(Chain, Depth), names_to = "Variable", values_to = "Value") %>%
    # filter(Value > ) %>%
    ggplot(aes(x = Value, fill = Depth)) +
    geom_density(alpha = 0.3) +
    facet_grid(Depth ~ Variable) + # , scales = "free") +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Sampled phis across chains",
      caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )

  ggsave(paste0(plot_dir, "sampled_phis_across_chain_depths.png"),
    plot = p_phis,
    height = plot_height,
    width = plot_width
  )

  p_alphas <- mass_df %>%
    pivot_longer(-c(Chain, Depth), names_to = "Variable", values_to = "Value") %>%
    ggplot(aes(x = Value, fill = Depth)) +
    geom_density(alpha = 0.3) +
    facet_grid(Depth ~ Variable) + # , ncol = 1, scales = "free") +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Sampled mass parameter across different chain depths"
      # caption = "Dataset 1 is the cell cycle data, 2 is the RNA-seq data and 3 is the LOPIT data"
    )

  ggsave(paste0(plot_dir, "sampled_alphas_across_depths.png"),
    plot = p_alphas,
    height = plot_height,
    width = plot_width
  )

  p_cc_psm <- cell_cycle_cm_df %>%
    ggplot(aes(x = x, y = y, fill = Entry)) +
    geom_tile() +
    facet_grid(Depth ~ Width, labeller = label_both) +
    scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
    labs(title = "Cell cycle CMs")

  if(V == 3) {
    p_rna_psm <- rna_seq_cm_df %>%
      ggplot(aes(x = x, y = y, fill = Entry)) +
      geom_tile() +
      facet_grid(Depth ~ Width, labeller = label_both) +
      scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
      labs(title = "RNA-seq CMs")

    ggsave(paste0(plot_dir, "rna_cms.png"),
      plot = p_rna_psm,
      height = plot_height,
      width = plot_width
    )

  }

  p_lopit_psm <- lopit_cm_df %>%
    ggplot(aes(x = x, y = y, fill = Entry)) +
    geom_tile() +
    facet_grid(Depth ~ Width, labeller = label_both) +
    scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
    labs(title = "LOPIT CMs")

  ggsave(paste0(plot_dir, "cell_cycle_cms.png"),
    plot = p_cc_psm,
    height = plot_height,
    width = plot_width
  )

  ggsave(paste0(plot_dir, "lopit_cms.png"),
    plot = p_lopit_psm,
    height = plot_height,
    width = plot_width
  )

  p_evidence <- fit_df %>%
    ggplot(aes(x = Evidence, fill = Depth)) +
    geom_density() +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Evidence for different chain depths"
    )

  p_cll <- fit_df %>%
    ggplot(aes(x = Complete_log_likelihood, fill = Depth)) +
    geom_density() +
    ggthemes::scale_fill_colorblind() +
    labs(
      title = "Complete log ikelihood for different chain depths"
    )

  ggsave(paste0(plot_dir, "evidence_densitiy.png"),
    plot = p_evidence,
    height = plot_height,
    width = plot_width
  )

  ggsave(paste0(plot_dir, "complete_log_likelihood_density.png"),
    plot = p_cll,
    height = plot_height,
    width = plot_width
  )

  p_mad_depth <- mean_absolute_difference_df |>
    dplyr::filter(Quantity_varied == "Depth") |>
    ggplot(aes(x = Depth, y = Mean_absolute_difference, colour = factor(Width))) +
    geom_line() +
    facet_wrap(~Dataset) +
    ggthemes::scale_color_colorblind() +
    labs(title = "Mean absoulte difference for varying depths")


  p_mad_width <- mean_absolute_difference_df |>
    dplyr::filter(Quantity_varied == "Width") |>
    ggplot(aes(x = Width, y = Mean_absolute_difference, colour = factor(Depth))) +
    geom_line() +
    facet_wrap(~Dataset) +
    ggthemes::scale_color_colorblind() +
    labs(title = "Mean absoulte difference for varying widths")

  ggsave(paste0(plot_dir, "mean_abs_diff_depth.png"),
    plot = p_mad_depth,
    height = plot_height,
    width = plot_width
  )

  ggsave(paste0(plot_dir, "mean_abs_diff_width.png"),
    plot = p_mad_width,
    height = plot_height,
    width = plot_width
  )
}

cat("\nPlotting complete.\nSave output.")

output <- list()
output$classification_probability <- classification_probability[[n_models]]
output$predicted_partitions <- vector("list", V)
output$cms <- vector("list", V)
output$allocations <- vector("list", V)
for (v in view_inds) {
  output$predicted_partitions[[v]] <- predicted_partitions[[v]][[n_models]]
  output$cms[[v]] <- cms[[v]][[n_models]]
  output$allocations[[v]] <- allocations[[number_depths]][chains_used, , v]
}
output$fit <- fit_df
output$mass <- mass_df
output$phis <- phi_df

saveRDS(output, file = output_file)

cat("\n# === SCRIPT COMPLETED ================================================")
t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))
