
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

generateGaussianView <- function(cluster_means, std_devs, N, P, labels,
                                 row_names = paste0("Person_", 1:n),
                                 col_names = paste0("Gene_", 1:p)) {
  gen_data <- matrix(0, N, P)

  for (ii in seq(1, P)) {
    reordered_cluster_means <- sample(cluster_means)
    reordered_std_devs <- sample(std_devs)

    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (jj in seq(1, N)) {
      .mu <- reordered_cluster_means[labels[jj]]
      .sd <- reordered_std_devs[labels[jj]]

      # print(labels[jj])
      # print(.mu)
      # print(.sd)

      gen_data[jj, ii] <- rnorm(1,
        mean = .mu,
        sd = .sd
      )
    }
  }
  gen_data
}

generateCategoricalView <- function(probs, N, P, labels,
                                    row_names = paste0("Person_", 1:n),
                                    col_names = paste0("Gene_", 1:p)) {
  gen_data <- matrix(0, N, P)

  for (ii in seq(1, P)) {
    reordered_probs <- sample(probs)

    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (jj in seq(1, N)) {
      .p <- reordered_probs[labels[jj]]

      gen_data[jj, ii] <- sample(c(0, 1), 1, prob = c(1 - .p, .p))
    }
  }
  gen_data
}

generateViewGivenStructure <- function(generating_structure, frac_present, P, K, class_weights, type, params) {
  N <- length(generating_structure)
  labels_transferred <- sample(c(0, 1), N, replace = TRUE, prob = c(1 - frac_present, frac_present))
  labels_transferred_ind <- which(labels_transferred == 1)
  N_transferred <- sum(labels_transferred)
  N_new <- N - N_transferred
  identical_view <- N_new == N
  K_ind <- seq(1, K)
  labels <- rep(0, N)
  labels[labels_transferred_ind] <- generating_structure[labels_transferred_ind]
  if (!identical_view) {
    new_labels <- sample(K_ind, N - N_transferred, replace = TRUE, prob = class_weights)
    labels[-labels_transferred_ind] <- new_labels
  }
  gaussian <- type == "G"
  categorical <- type == "C"

  if (gaussian) {
    means <- params$means
    std_devs <- params$std_devs
    gen_data <- generateGaussianView(means, std_devs, N, P, labels)
  }
  if (categorical) {
    probs <- params$probs
    gen_data <- generateCategoricalView(probs, N, P, labels)
  }

  row.names(gen_data) <- names(generating_structure)

  list(
    gen_data = gen_data,
    labels = labels
  )
}


generateMDIDataLabels <- function(N, V, K, gammas, phis) {
  
  N_inds <- seq(1, N)
  V_inds <- seq(1, V)
  
  labels <- matrix(0, N, V)
  
  for(n in seq(1, N)) {
    for(v in seq(1, V)) {
      V_comp_inds <- V_inds[-v]
      comps <- seq(1, K[v])
      weights <- gammas[, v]
      if(v > 1) {
        for(w in seq(1, v - 1)) {
          label_in_view_w <- labels[n, w]
          weights[label_in_view_w] <- weights[label_in_view_w] * (1.0 + phis[v, w])
        }
      }
      weights <- weights / sum(weights)
      labels[n, v] <- sample(comps, prob = weights, size = 1)
    }
  }
  labels
}


input_arguments <- function() {
  option_list <- list(

    # Number of MCMC iterations
    optparse::make_option(c("-n", "--n_datasets"),
      type = "numeric",
      default = 100L,
      help = "Number of datasets to generate in each scenario [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-s", "--seed"),
      type = "numeric",
      default = 1,
      help = "Random seed that defines the data generation [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--test_frac"),
      type = "numeric",
      default = 0.8,
      help = "Fraction of labels unobserved in dataset 1.",
      metavar = "numeric"
    ),
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./Simulations/Data/",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

cat("\n# === SIMULATION DATA GENERATION ======================================")

cat("\nRead in arguments.")

# Pass the inputs from the command line
args <- input_arguments()

# Random seed used defining this fold
seed <- args$seed
set.seed(seed)

# Directories for input and output respectively
save_dir <- args$save_dir

# Fraction of labels to consider unobserved in dataset 1
test_frac <- args$test_frac

sample_file <- paste0(save_dir, "mcmcMDItestdataTestFrac", test_frac, "Seed", seed, ".rds")
plot_file <- paste0(save_dir, "mcmcMDItestdataTestFrac", test_frac, "Seed", seed, ".png")

N <- 200
P <- c(30, 25, 20)
K <- c(7, 9, 11)
frac_present <- c(1, 0.5, 0.4)
frac_known <- 0.3

# Scenarios are defined by changing how view 1 is misspecified
scenarios <- c("Gaussian", "MVT", "MixGaussians", "MixMVTs")
mix_scenarios <- c("MixGaussians", "MixMVTs")
n_scenario <- length(scenarios)
n_datasets <- args$n_datasets

for (n in seq(1, n_datasets)) {
  P_1 <- P[1]
  K_1 <- K[1]

  # Parameters for view 1
  group_means <- rnorm(K_1, sd = 1.75)
  group_sds <- rgamma(K_1, 4, 2)
  group_dfs <- c(5, 10, 5, 5, 10, 5, 20)
  group_weights <- c(0.25, 0.2, 0.15, rep(0.4 / 4, 4))

  # If each cluster is generated from a mixture these come into play
  B <- 3
  batch_shift <- rnorm(B, sd = 0.4)
  batch_scale <- 1 + rgamma(B, 10, 20)
  batch_weights <- c(0.4, 0.4, 0.2)

  # Parameters for view 2
  frac_present_2 <- frac_present[2]
  P_2 <- P[2]
  K_2 <- K[2]
  class_weights_2 <- rgamma(K_2, 10)
  class_weights_2 <- class_weights_2 / sum(class_weights_2)
  params_2 <- list(
    means = rnorm(K_2, sd = 0.85),
    std_devs = rgamma(K_2, 2, 2)
  )

  # Parameters for view 3
  frac_present_3 <- frac_present[3]
  P_3 <- P[3]
  K_3 <- K[3]
  class_weights_3 <- rgamma(K_3, 10, 10)
  class_weights_3 <- class_weights_3 / sum(class_weights_3)
  params_3 <- list(
    means = rnorm(K_3, sd = 1.05),
    std_devs = rgamma(K_3, 2, 2)
  )
  
  fixed <- sample(c(0, 1), N, replace = TRUE, prob = c(1 - frac_known, frac_known))

  for (scn in scenarios) {
    if (scn == "Gaussian") {
      view_1 <- generateBatchData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = rep(0, B),
        batch_scale = rep(1, B),
        group_weights = group_weights,
        batch_weights = rep(1 / B, B),
        type = "MVN"
      )
    }
    if (scn == "MVT") {
      view_1 <- generateBatchData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = rep(0, B),
        batch_scale = rep(1, B),
        group_weights = group_weights,
        batch_weights = rep(1 / B, B),
        group_dfs = group_dfs,
        type = "MVT"
      )
    }
    if (scn == "MixGaussians") {
      view_1 <- generateBatchData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = batch_shift,
        batch_scale = batch_scale,
        group_weights = group_weights,
        batch_weights = batch_weights,
        type = "MVN"
      )
    }
    if (scn == "MixMVTs") {
      view_1 <- generateBatchData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = batch_shift,
        batch_scale = batch_scale,
        group_weights = group_weights,
        batch_weights = batch_weights,
        group_dfs = group_dfs,
        type = "MVT"
      )
    }

    view_2 <- generateViewGivenStructure(view_1$group_IDs,
      frac_present = frac_present_2,
      P = P_2,
      K = K_2,
      type = "G",
      class_weights = class_weights_2,
      params = params_2
    )

    view_3 <- generateViewGivenStructure(view_1$group_IDs,
      frac_present = frac_present_3,
      P = P_3,
      K = K_3,
      type = "G",
      class_weights = class_weights_3,
      params = params_3
    )
    
    sim_data <- list(
      "Data" = list(
        "View_1" = view_1$observed_data,
        "View_2" = view_2$gen_data,
        "View_3" = view_3$gen_data
      ), 
      "Clustering" = list(
        "View_1" = view_1$group_IDs,
        "View_2" = view_2$labels,
        "View_3" = view_3$labels
        ),
      "Fixed" = fixed,
      "Frac_known" = frac_known
      )
    
    file_name <- paste0(save_dir, scn, "/dataset_", n, ".rds")
    saveRDS(sim_data, file_name)
  }
}

cat("\n === SCRIPT COMPLETE ==================================================")
