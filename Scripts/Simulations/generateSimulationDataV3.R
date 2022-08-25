
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(batchmix))
suppressMessages(library(magrittr))
suppressMessages(library(mcclust))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(RcppParallel))
suppressMessages(library(optparse))


generateLogPoissonData <- function(lambdas, N, P, labels,
                                   row_names = paste0("Person_", seq(1, N)),
                                   col_names = paste0("Gene_", seq(1, P))) {
  classes <- unique(labels)
  gen_data <- matrix(0, N, P)
  for (p in seq(1, P)) {
    lambdas_p <- sample(lambdas)
    for (k in classes) {
      class_inds <- which(labels == k)
      n_k <- length(class_inds)
      gen_data[class_inds, p] <- log(1.0 + rpois(n_k, lambdas_p[k])) + rnorm(n_k)
    }
  }
  row.names(gen_data) <- row_names
  colnames(gen_data) <- col_names
  gen_data <- scale(gen_data)
  gen_data
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
                                 row_names = paste0("Person_", seq(1, N)),
                                 col_names = paste0("Gene_", seq(1, P))) {
  gen_data <- matrix(0, N, P)

  for (ii in seq(1, P)) {
    reordered_cluster_means <- sample(cluster_means)
    reordered_std_devs <- sample(std_devs)

    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (jj in seq(1, N)) {
      k <- labels[jj]
      .mu <- reordered_cluster_means[k]
      .sd <- reordered_std_devs[k]

      # print(labels[jj])
      # print(.mu)
      # print(.sd)

      gen_data[jj, ii] <- rnorm(1,
        mean = .mu,
        sd = .sd
      )
    }
  }
  row.names(gen_data) <- row_names
  colnames(gen_data) <- col_names

  gen_data <- scale(gen_data)
  gen_data
}

generateMVTView <- function(cluster_means, std_devs, dfs, N, P, labels,
                            row_names = paste0("Person_", seq(1, N)),
                            col_names = paste0("Gene_", seq(1, P))) {
  gen_data <- matrix(0, N, P)

  for (ii in seq(1, P)) {
    reordered_cluster_means <- sample(cluster_means)
    reordered_std_devs <- sample(std_devs)

    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (jj in seq(1, N)) {
      k <- labels[jj]
      .mu <- reordered_cluster_means[k]
      .sd <- reordered_std_devs[k]

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
  row.names(gen_data) <- row_names
  colnames(gen_data) <- col_names

  gen_data <- scale(gen_data)
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
    ),
    optparse::make_option(c("--plot_data"),
      type = "logical",
      default = FALSE,
      help = "Instruction to make heatmaps of generated data [default= %default]",
      metavar = "logical"
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

plot_data <- args$plot_data

# Fraction of labels to consider unobserved in dataset 1
test_frac <- args$test_frac
test_frac_str <- paste0(test_frac * 100, "%")
sample_file <- paste0(save_dir, "mcmcMDISimPhiModelTestFrac", test_frac_str, "Seed", seed, ".rds")
plot_file <- paste0(save_dir, "mcmcMDISimPhiModelTestFrac", test_frac_str, "Seed", seed, ".png")

N <- 200
P <- c(30, 25, 20)
K <- c(6, 7, 8)
K_max <- max(K)
V <- length(K)

weight_shape <- 10
weight_rate <- 2

frac_present <- c(1, 1, 1) # 0.5, 0.4)
frac_known <- 0.3

phis <- matrix(c(0, 15, 10, 15, 0, 5, 10, 5, 0), V, V)

# Scenarios are defined by changing how view 1 is misspecified
scenarios <- c("Gaussian", "MVT", "LogPoisson")
n_scenario <- length(scenarios)
n_datasets <- args$n_datasets

gammas <- matrix(c(
  4.0, 7.0, 7.0, 5.5, 3.0, 5.5, 0.0, 0.0,
  5.5, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 0.0,
  6.0, 6.0, 5.0, 3.0, 5.5, 7.0, 4.0, 4.0
), ncol = V)


# gammas <- matrix(c(
#   4.0, 7.0, 7.0, 5.5, 3.0, 5.5, 6.0, 0.0, 0.0, 0.0, 0.0,
#   5.5, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 6.0, 5.5, 0.0, 0.0,
#   6.0, 6.0, 5.0, 3.0, 5.5, 7.0, 4.0, 4.0, 4.0, 4.0, 4.0
# ), ncol = V)

for (n in seq(1, n_datasets)) {
  # gammas <- matrix(0, K_max, V)
  # for (v in seq(1, V)) {
  #   gammas[seq(1, K[v]), v] <- rgamma(K[v], weight_shape, weight_rate)
  # }

  labels <- generateMDIDataLabels(N, V, K, gammas, phis)

  P_1 <- P[1]
  K_1 <- K[1]

  # Parameters for view 2
  P_2 <- P[2]
  K_2 <- K[2]


  # Parameters for view 3
  P_3 <- P[3]
  K_3 <- K[3]

  fixed <- sample(c(0, 1), N, replace = TRUE, prob = c(1 - frac_known, frac_known))

  for (scn in scenarios) {
    if (scn %in% c("Gaussian", "MVT")) {
      # Parameters for view 1
      group_means <- c(-1.50, -1.25, 1.25, 0.50, -1.00, 1.00, 1.75) #  rnorm(K_1, sd = 1.75)
      group_sds <- c(1.50, 1.25, 1.75, 2.25, 2.75, 1.50, 3.0) # rgamma(K_1, 4, 2)
      group_dfs <- c(5, 10, 5, 5, 10, 5, 10)

      params_2 <- list(
        means = c(-1.00, -0.50, -0.75, 0.00, 0.25, 0.50, 0.75, 0.50, 0.00), # rnorm(K_2, sd = 0.85),
        std_devs = c(0.5, 1.0, 0.5, 1.0, 0.5, 0.25, 1.25, 0.5, 1.0), # rgamma(K_2, 2, 2),
        dfs = sample(c(5, 10), size = K_2, replace = TRUE)
      )
      params_3 <- list(
        means = c(-1.00, -0.50, -1.25, -0.50, -0.25, -0.25, 0.25, 1.25, 1.25, 1.00, 0.50), # rnorm(K_3, sd = 1.05),
        std_devs = c(0.50, 0.75, 1.50, 1.00, 0.50, 0.75, 0.75, 0.50, 1.50, 1.00, 0.50), # rgamma(K_3, 2, 2),
        dfs = sample(c(5, 10), size = K_3, replace = TRUE)
      )
    }

    if (scn == "LogPoisson") {
      rates <- c(5, 15, 25, 35, 30, 15, 40) # rgamma(K_1, 2.5, 0.1) + 1
      params_2 <- list(rates = c(4, 13, 21, 26, 55, 43, 31, 28, 15)) #  rgamma(K_2, 2.5, 0.1) + 1)
      params_3 <- list(rates = c(4, 8, 44, 32, 16, 28, 55, 24, 12, 18, 28)) # rgamma(K_3, 2.5, 0.1) + 1)
    }

    if (scn == "Gaussian") {
      view_1 <- generateGaussianView(
        group_means,
        group_sds,
        N,
        P_1,
        labels[, 1]
      )

      view_2 <- generateGaussianView(
        params_2$means,
        params_2$std_devs,
        N,
        P_2,
        labels[, 2]
      )
      view_3 <- generateGaussianView(
        params_3$means,
        params_3$std_devs,
        N,
        P_3,
        labels[, 3]
      )
    }
    if (scn == "MVT") {
      view_1 <- generateMVTView(
        group_means,
        group_sds,
        group_dfs,
        N,
        P_1,
        labels[, 1]
      )

      view_2 <- generateMVTView(
        params_2$means,
        params_2$std_devs,
        params_2$dfs,
        N,
        P_2,
        labels[, 2]
      )
      view_3 <- generateMVTView(
        params_3$means,
        params_3$std_devs,
        params_3$dfs,
        N,
        P_3,
        labels[, 3]
      )
    }
    if (scn == "LogPoisson") {
      view_1 <- generateLogPoissonData(rates, N, P_1, labels[, 1])
      view_2 <- generateLogPoissonData(params_2$rates, N, P_2, labels[, 2])
      view_3 <- generateLogPoissonData(params_3$rates, N, P_3, labels[, 3])
    }

    row_order <- findOrder(view_1)

    view_1_gg_df <- view_1 |>
      prepDataForggHeatmap(row_order = row_order) |>
      mutate(View = 1)
    view_2_gg_df <- view_2 |>
      prepDataForggHeatmap(row_order = row_order) |>
      mutate(View = 2)
    view_3_gg_df <- view_3 |>
      prepDataForggHeatmap(row_order = row_order) |>
      mutate(View = 3)

    gg_df <- rbind(view_1_gg_df, view_2_gg_df, view_3_gg_df)

    labels_gg_df <- labels |>
      prepDataForggHeatmap(row_order = row_order, col_order = c(1, 2, 3)) |>
      mutate(View = "Labels")

    sim_data <- list(
      "Data" = list(
        "View_1" = view_1,
        "View_2" = view_2,
        "View_3" = view_3
      ),
      "Clustering" = list(
        "View_1" = labels[, 1],
        "View_2" = labels[, 2],
        "View_3" = labels[, 3]
      ),
      "Fixed" = fixed,
      "Frac_known" = frac_known
    )

    plot_titles <- paste0(scn, " dataset ", n, ": view ", seq(1, V))
    plot_names <- paste0(save_dir, scn, "/dataset_", n, "_view_", seq(1, V), ".png")
    gg_plot_name <- paste0(save_dir, scn, "/dataset_", n, ".png")
    gg_title <- paste0(scn, " dataset ", n, " annotated by true labels")

    if (plot_data) {
      p_data <- gg_df |>
        ggplot(aes(x = x, y = y, fill = Entry)) +
        geom_tile() +
        facet_wrap(~View, scales = "free_x", nrow = 1, labeller = label_both) +
        scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900") +
        labs(y = NULL, x = "Feature") +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank()
        )

      p_labels <- labels_gg_df |>
        ggplot(aes(x = x, y = y, fill = factor(Entry))) +
        geom_tile() +
        facet_wrap(~View, scales = "free_x", nrow = 1) +
        ggthemes::scale_fill_colorblind() +
        labs(fill = "Label", x = "View", y = "Row number") +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank()
        )

      p_gg <- p_labels + p_data +
        plot_layout(widths = c(1, 6), guides = "collect") +
        plot_annotation(title = gg_title)

      ggsave(gg_plot_name, p_gg, height = 8, width = 12)

      annotatedHeatmap(view_1, labels[, 1],
        main = plot_titles[1],
        show_rownames = FALSE,
        filename = plot_names[1]
      )

      annotatedHeatmap(view_2, labels[, 2],
        main = plot_titles[2],
        show_rownames = FALSE,
        filename = plot_names[2]
      )
      annotatedHeatmap(view_3, labels[, 3],
        main = plot_titles[3],
        show_rownames = FALSE,
        filename = plot_names[3]
      )
    }
    file_name <- paste0(save_dir, scn, "/dataset_", n, ".rds")
    saveRDS(sim_data, file_name)
  }
}

cat("\n === SCRIPT COMPLETE ==================================================")
