

suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(batchmix))
suppressMessages(library(magrittr))
suppressMessages(library(mcclust))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(RcppParallel))
suppressMessages(library(optparse))


#' @title multi-class F1
#' @description Calculate the accuracy, the per-class F1 score and the macro-
#' and weighted-F1 scores.
#' @param pred Factor. The predicted classes.
#' @param truth Factor wit the same levels as ``pred``. The true classification.
#' @return A list containing ``accuracy``, the prediction accuracy, ``f1``, a
#' vector of the f1 score in each class, ``macro_f1``, the macro f1 score across
#' classes, and ``weighted_f1``, the f1 score weighted by class membership.
multiClassF1 <- function(pred, truth) {
  mismatch_in_lengths <- length(pred) != length(truth)
  if (mismatch_in_lengths) {
    stop("Prediction vector and truth must be of the same length.")
  }

  # Confusion matrix for current fold
  conf <- caret::confusionMatrix(
    data = pred,
    reference = truth
  )$table

  # conf <- t( get.confusion_matrix(
  #   truth = truth,
  #   predicted = pred
  # ) )

  N <- length(pred)
  n_levels <- ncol(conf)
  seq_levels <- seq(1, n_levels)

  f1 <- rep(0, n_levels)

  for (ii in seq_levels) {
    precision <- conf[ii, ii] / sum(conf[ii, ])
    recall <- conf[ii, ii] / sum(conf[, ii])
    f1[ii] <- (2 * precision * recall) / (precision + recall)
    if (is.na(f1[ii])) {
      f1[ii] <- 0
    }
  }

  macro_f1 <- mean(f1)
  weighted_f1 <- sum(f1 * colSums(conf)) / N
  accuracy <- sum(pred == truth, na.rm = TRUE) / N
  list(
    "accuracy" = accuracy,
    "f1" = f1,
    "macro_f1" = macro_f1,
    "weighted_f1" = weighted_f1
  )
}



makePhiDF <- function(mcmc) {
  n_phis <- choose(mcmc$V, 2)
  phi_names <- c()
  for (v in seq(1, mcmc$V - 1)) {
    for (u in seq(v + 1, mcmc$V)) {
      .name <- paste0("Phi_", v, u)
      phi_names <- c(phi_names, .name)
    }
  }

  phi_df <- mcmc$phis |>
    data.frame() |>
    magrittr::set_colnames(phi_names) |>
    dplyr::mutate(Iteration = seq(mcmc$burn + mcmc$thin, mcmc$R, mcmc$thin)) |>
    tidyr::pivot_longer(-Iteration, values_to = "Sampled_value", names_to = "Parameter")

  phi_df
}

input_arguments <- function() {
  option_list <- list(

    # Number of MCMC iterations
    optparse::make_option(c("-R", "--R"),
      type = "numeric",
      default = 20000L,
      help = "Number of iterations to run in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-t", "--thin"),
      type = "numeric",
      default = 125L,
      help = "Thinning factor in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-b", "--burn"),
      type = "numeric",
      default = 4000L,
      help = "Thinning factor in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-i", "--index"),
      type = "numeric",
      default = 1L,
      help = "Index of simulation to analyse [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--n_chains"),
      type = "numeric",
      default = 1L,
      help = "Number of MCMC chains to run [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-K", "--K"),
      type = "numeric",
      default = 50,
      help = paste(
        "Number of components modelled in each dataset. If a dataset is",
        "semi-supervised then the number of unique labels is modelled, if",
        "unsupervised we default to 50."
      ),
      metavar = "numeric"
    ),
    optparse::make_option(c("--test_frac"),
      type = "numeric",
      default = 0.8,
      help = "Fraction of labels unobserved in dataset 1.",
      metavar = "numeric"
    ),
    optparse::make_option(c("--scn"),
      type = "character",
      default = "Gaussian",
      help = "Simulation scenario [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./Simulations/Output/",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(c("--data_dir"),
      type = "character",
      default = "./Simulations/Data/",
      help = "File to read input from from",
      metavar = "character"
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

cat("\n# === SIMULATION DATA GENERATION =====================================\n")

cat("\nRead in arguments.")
set.seed(1)
mdiHelpR::setMyTheme()

# Pass the inputs from the command line
args <- input_arguments()

# Dataset index used defining this fold
index <- args$index

# Scenario
scn <- args$scn

# Directories for input and output respectively
data_dir <- paste0(args$data_dir, scn, "/")
save_dir <- paste0(args$save_dir, scn, "/")

# Fraction of labels to consider unobserved in dataset 1
test_frac <- args$test_frac

# MCMC sampler arguments
R <- args$R
thin <- args$thin
burn <- args$burn

# Number of chains modelled
n_chains <- args$n_chains

# The number of components modelled
K <- args$K

# Number of clusters modelled in the unsupervised views
n_clust_unsupervised <- 50

# Input and output files
sim_data_file <- paste0(data_dir, "dataset_", index, ".rds")
test_perc <- test_frac * 100
save_file <- paste0(save_dir, "dataset_", index, "_testFrac_", test_perc, "_nChains_", n_chains, "_K_", K, ".png")

cat("\n Arguments passed. Reading in data.")

# Read in data
sim_data_lst <- readRDS(sim_data_file)
sim_data <- sim_data_lst$Data
sim_cl <- sim_data_lst$Clustering
fixed_view_1 <- sim_data_lst$Fixed

N <- nrow(sim_data$View_1)
V <- length(sim_data)

# Modelling inputs
types <- c("G", "G", "G")
K_max <- c(K, n_clust_unsupervised, n_clust_unsupervised)

fixed <- initial_labels <- matrix(0, nrow = N, 3)
fixed[, 1] <- fixed_view_1

test_inds <- which(fixed_view_1 == 0)

initial_labels[, 1] <- sim_cl$View_1
X <- list(
  sim_data$View_1,
  sim_data$View_2,
  sim_data$View_3
)

cat("\n# === MODELLING ======================================================\n")

cat("\nSemi-supervised MDI.")
mcmc_semi <- tagmReDraft::runMCMCChains(X, n_chains,
  R = R,
  thin = thin,
  types = types,
  K = K_max,
  fixed = fixed,
  initial_labels = initial_labels
)

cat("\nUnsupervised MDI.")
mcmc_un <- tagmReDraft::runMCMCChains(X, n_chains,
  R = R,
  thin = thin,
  types = types,
  K = rep(n_clust_unsupervised, V),
  fixed = matrix(0, N, V),
  initial_labels = initial_labels
)

cat("\n# === MODELLING COMPLETE =============================================\n")

cat("\nProcess output.")
new_semi <- predictFromMultipleChains(mcmc_semi, burn)
new_un <- predictFromMultipleChains(mcmc_un, burn)
psms_un <- psms_semi <- list()
K_pred_un <- K_pred_semi <- c()
for (v in seq(1, V)) {
  psms_semi[[v]] <- createSimilarityMat(new_semi$allocations[[v]])
  psms_un[[v]] <- createSimilarityMat(new_un$allocations[[v]])

  .cl_semi <- maxpear(psms_semi[[v]])$cl
  .cl_un <- maxpear(psms_un[[v]])$cl

  K_pred_semi[v] <- .k_semi <- length(unique(.cl_semi))
  K_pred_un[v] <- .k_un <- length(unique(.cl_un))

  k_true <- length(unique(sim_cl[[v]]))

  if (.k_semi < k_true) {
    K_pred_semi[v] <- k_true
  }
  if (.k_un < k_true) {
    K_pred_un[v] <- k_true
  }

  new_semi$pred[[v]] <- factor(.cl_semi, levels = seq(1, K_pred_semi[v]))
  new_un$pred[[v]] <- factor(.cl_un, levels = seq(1, K_pred_un[v]))
}

new_semi$psms <- psms_semi
new_un$psms <- psms_un

# multiClassF1(new_semi$pred[[1]], factor(sim_cl$View_1, levels = seq(1, K_pred_semi[1])))
# multiClassF1(new_un$pred[[1]], factor(sim_cl$View_1, levels = seq(1, K_pred_un[1])))
#
# multiClassF1(new_semi$pred[[2]], factor(sim_cl$View_2, levels = seq(1, K_pred_semi[2])))
# multiClassF1(new_un$pred[[2]], factor(sim_cl$View_2, levels = seq(1, K_pred_un[2])))
#
# multiClassF1(new_semi$pred[[3]], factor(sim_cl$View_3, levels = seq(1, K_pred_semi[3])))
# multiClassF1(new_un$pred[[3]], factor(sim_cl$View_3, levels = seq(1, K_pred_un[3])))

cat("\nCompare to truth using the adjusted rand index.")

ari_semi_1 <- mcclust::arandi(new_semi$pred[[1]][test_inds], sim_cl$View_1[test_inds])
ari_un_1 <- mcclust::arandi(new_un$pred[[1]][test_inds], sim_cl$View_1[test_inds])

ari_semi_2 <- mcclust::arandi(new_semi$pred[[2]], sim_cl$View_2)
ari_un_2 <- mcclust::arandi(new_un$pred[[2]], sim_cl$View_2)

ari_semi_3 <- mcclust::arandi(new_semi$pred[[3]], sim_cl$View_3)
ari_un_3 <- mcclust::arandi(new_un$pred[[3]], sim_cl$View_3)

results_df <- data.frame(
  "Scenario" = rep(scn, V),
  "Index" = rep(index, V),
  "Semi-supservised" = c(ari_semi_1, ari_semi_2, ari_semi_3),
  "Unsupservised" = c(ari_un_1, ari_un_2, ari_un_3),
  "Difference" = c(
    ari_semi_1 - ari_un_1,
    ari_semi_2 - ari_un_2, ari_semi_3 - ari_un_3
  )
)

cat("\nSave results\n.")
knitr::kable(results_df, digits = 3)

out_lst <- list(
  "MCMC" = list(
    "Semisupservised" = new_semi,
    "Unsupservised" = new_un
  ),
  "ARI" = results_df
)

saveRDS(out_lst, save_file)

cat("\n# === SCRIPT COMPLETE ================================================\n")
