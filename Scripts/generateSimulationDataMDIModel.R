
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


generateLogPoissonData <- function(lambdas, N, P, labels) {
  
  classes <- unique(labels)
  data <- matrix(0, N, P)
  for(p in seq(1, P)) {
    lambdas_p <- sample(lambdas)
    for(k in classes) {
      class_inds <- which(labels == k)
      n_k <- length(class_inds)
      data[class_inds, p] <- log(1.0 + rpois(n_k, lambdas_p[k])) + rnorm(n_k)
    }
  }
  row.names(data) <- paste0("Person_", seq(1, N))
  data
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

generateMVTView <- function(cluster_means, std_devs, dfs, N, P, labels,
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
      # Draw a point from a standard normal
      x <- stats::rnorm(1)
      chi_draw <- stats::rchisq(1, dfs[k])
      
      gen_data[jj, ii] <- x * .sd * sqrt(dfs[k] / chi_draw) + .mu
      
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
  labels_transferred <- rep(1, N) # sample(c(0, 1), N, replace = TRUE, prob = c(1 - frac_present, frac_present))
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

  for (n in seq(1, N)) {
    for (v in seq(1, V)) {
      V_comp_inds <- V_inds[-v]
      comps <- seq(1, K[v])
      weights <- gammas[comps, v]
      if (v > 1) {
        for (w in seq(1, v - 1)) {
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


generateSimData <- function(N,
                            P,
                            group_means,
                            group_std_devs,
                            batch_shift,
                            batch_scale,
                            batch_weights,
                            group_IDs,
                            type = "MVN",
                            group_dfs = NULL,
                            frac_known = 0.2,
                            permute_variables = TRUE,
                            scale_data = FALSE) {
  mvn_generated <- type == "MVN"
  mvt_generated <- type == "MVT"

  # The number of batches and groups to generate
  B <- length(batch_weights)
  K <- length(group_means)

  # # Allow for varying groups within batches
  # varying_group_within_batch <- is.matrix(group_weights)

  # The batch labels for the N points
  batch_IDs <- sample(seq(1, B), N, replace = TRUE, prob = batch_weights)

  # Fixed labels
  fixed <- sample(seq(0, 1), N,
    replace = TRUE,
    prob = c(1 - frac_known, frac_known)
  )

  # The data matrices
  observed_data <- matrix(nrow = N, ncol = P)

  # Iterate over each of the columns permuting the means associated with each
  # label.
  for (p in seq(1, P))
  {

    # To provide different information in each column, randomly sample the
    # parameters with each group and batch
    reordered_group_means <- group_means
    reordered_group_std_devs <- group_std_devs

    reordered_batch_shift <- batch_shift
    reordered_batch_scale <- batch_scale

    if (permute_variables) {
      # To provide different information in each column, randomly sample the
      # parameters with each group and batch
      reordered_group_means <- sample(group_means)
      reordered_group_std_devs <- sample(group_std_devs)

      reordered_batch_shift <- sample(batch_shift)
      reordered_batch_scale <- sample(batch_scale)
    }

    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (n in seq(1, N)) {

      # Draw a point from a standard normal
      x <- stats::rnorm(1)
      k <- group_IDs[n]
      b <- batch_IDs[n]

      # For ease of reading the following lines, create group and batch parameters
      .mu <- reordered_group_means[k]
      .sd <- reordered_group_std_devs[k]
      .m <- reordered_batch_shift[b]
      .s <- reordered_batch_scale[b]

      if (mvn_generated) {
        # Adjust to the batched group distribution
        observed_data[n, p] <- x * .sd * .s + .mu + .m
      }

      if (mvt_generated) {
        chi_draw <- stats::rchisq(1, group_dfs[k])

        # Adjust to the batched group distribution
        observed_data[n, p] <- x * .sd * .s * sqrt(group_dfs[k] / chi_draw) + .mu + .m
      }
    }
  }

  row.names(observed_data) <- paste0("Item_", seq(1, N))
  if (scale_data) {
    observed_data <- scale(observed_data)
  }

  # Return the data, the data without batch effects, the allocation labels and
  # the batch labels.
  list(
    observed_data = observed_data,
    group_IDs = group_IDs,
    batch_IDs = batch_IDs,
    fixed = fixed
  )
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
      default = "./Simulations/DataPhiModel/",
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
test_frac_str <- paste0(test_frac * 100, "%")
sample_file <- paste0(save_dir, "mcmcMDISimPhiModelTestFrac", test_frac_str, "Seed", seed, ".rds")
plot_file <- paste0(save_dir, "mcmcMDISimPhiModelTestFrac", test_frac_str, "Seed", seed, ".png")

N <- 200
P <- c(30, 25, 20)
K <- c(7, 9, 11)
K_max <- max(K)
V <- length(K)

weight_shape <- 10
weight_rate <- 2

frac_present <- c(1, 1, 1) # 0.5, 0.4)
frac_known <- 0.3

# gammas <- matrix(0, K_max, V)
# for(v in seq(1, V)) {
#   gammas[seq(1, K[v]), v] <- rgamma(K[v], weight_shape, weight_rate)
# }

# for(v in seq(1, V - 1)) {
#   for(w in seq(v + 1, V)) {
#     phis[v, w] <- rgamma(1, )
#   }
# }


phis <- matrix(c(0, 8, 16, 8, 0, 1, 16, 1, 0), V, V)


# Scenarios are defined by changing how view 1 is misspecified
scenarios <- c("Gaussian", "MVT", "MixGaussians", "MixMVTs")
mix_scenarios <- c("MixGaussians", "MixMVTs")
n_scenario <- length(scenarios)
n_datasets <- args$n_datasets

for (n in seq(1, n_datasets)) {
  gammas <- matrix(0, K_max, V)
  for (v in seq(1, V)) {
    gammas[seq(1, K[v]), v] <- rgamma(K[v], weight_shape, weight_rate)
  }

  labels <- generateMDIDataLabels(N, V, K, gammas, phis)

  P_1 <- P[1]
  K_1 <- K[1]

  # Parameters for view 1
  group_means <- rnorm(K_1, sd = 1.75)
  group_sds <- rgamma(K_1, 4, 2)
  group_dfs <- c(5, 10, 5, 5, 10, 5, 20)
  # group_weights <- c(0.25, 0.2, 0.15, rep(0.4 / 4, 4))

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
      view_1 <- generateSimData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = rep(0, B),
        batch_scale = rep(1, B),
        # group_weights = group_weights,
        batch_weights = rep(1 / B, B),
        group_IDs = labels[, 1],
        type = "MVN"
      )
    }
    if (scn == "MVT") {
      view_1 <- generateSimData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = rep(0, B),
        batch_scale = rep(1, B),
        # group_weights = group_weights,
        batch_weights = rep(1 / B, B),
        group_IDs = labels[, 1],
        group_dfs = group_dfs,
        type = "MVT"
      )
    }
    if (scn == "MixGaussians") {
      view_1 <- generateSimData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = rep(0, B),
        batch_scale = rep(1, B),
        group_IDs = labels[, 1],
        batch_weights = batch_weights,
        type = "MVN"
      )
    }
    if (scn == "MixMVTs") {
      view_1 <- generateSimData(N, P_1,
        group_means = group_means,
        group_std_devs = group_sds,
        batch_shift = rep(0, B),
        batch_scale = rep(1, B),
        group_IDs = labels[, 1],
        batch_weights = batch_weights,
        group_dfs = group_dfs,
        type = "MVT"
      )
    }
    
    if(scn == "LogPoisson") {
      view_1 <- generateLogPoissonData(N, P, lambdas, labels[, 1])
    }
    
    names(view_1$group_IDs) <- row.names(view_1$observed_data)
    
    if(scn == "ViewsG") {
      cluster_means <- params_2$means
      std_devs <- params_2$std_devs
      view_2 <- generateGaussianView(cluster_means, std_devs, N, P, labels[, 2])
      
      cluster_means <- params_3$means
      std_devs <- params_3$std_devs
      view_3 <- generateGaussianView(cluster_means, std_devs, N, P, labels[, 3])
    }
    
    if(scn == "ViewsMVT") {
      cluster_means <- params_2$means
      std_devs <- params_2$std_devs
      dfs <- params_2$dfs
      
      view_2 <- generateMVTView(cluster_means, std_devs, dfs, N, P, labels[, 2])
      
      cluster_means <- params_3$means
      std_devs <- params_3$std_devs
      dfs <- params_3$dfs
      view_3 <- generateMVTView(cluster_means, std_devs, dfs, N, P, labels[, 3])
    }

    if(scn == "ViewsLogPoisson") {
      rates <- params_2$rates
      view_2 <- generateLogPoissonData(lambdas, N, P, labels[, 2])
      
      rates <- params_3$rates
      view_3 <- generateLogPoissonData(lambdas, N, P, labels[, 3])
    }
    
    sim_data <- list(
      "Data" = list(
        "View_1" = view_1$observed_data,
        "View_2" = view_2,
        "View_3" = view_3
      ),
      "Clustering" = list(
        "View_1" = view_1$group_IDs,
        "View_2" = labels[, 2],
        "View_3" = labels[, 3]
      ),
      "Fixed" = fixed,
      "Frac_known" = frac_known
    )

    plot_titles <- paste0(scn, " dataset ", n, ": view ", seq(1, V))
    plot_names <- paste0(save_dir, scn, "/dataset_", n, "_view_", seq(1, V), ".png")
    annotatedHeatmap(view_1$observed_data, view_1$group_IDs, 
                     main = plot_titles[1], 
                     show_rownames = FALSE, 
                     filename = plot_names[1])
    annotatedHeatmap(view_2$gen_data, view_2$labels,
                     main = plot_titles[2],
                     show_rownames = FALSE, 
                     filename = plot_names[2])
    annotatedHeatmap(view_3$gen_data, view_3$labels,
                     main = plot_titles[3], 
                     show_rownames = FALSE, 
                     filename = plot_names[3]
                     )
    
    file_name <- paste0(save_dir, scn, "/dataset_", n, ".rds")
    saveRDS(sim_data, file_name)
  }
}

cat("\n === SCRIPT COMPLETE ==================================================")
