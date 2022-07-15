
library(tagmReDraft)
library(mdiHelpR)
library(batchmix)

setMyTheme()
set.seed(1)

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
  if(! identical_view) {
    new_labels <- sample(K_ind, N - N_transferred, replace = TRUE, prob = class_weights)
    labels[-labels_transferred_ind] <- new_labels
  }
  gaussian <- type == "G"
  categorical <- type == "C"

  if(gaussian) {
    means <- params$means
    std_devs <- params$std_devs
    gen_data <- generateGaussianView(means, std_devs, N, P, labels)
  }
  if(categorical) {
    probs <- params$probs
    gen_data <- generateCategoricalView(probs, N, P, labels)
  }
  
  row.names(gen_data) <- names(generating_structure)
  
  list(
    gen_data = gen_data,
    labels = labels
  )
}

N <- 200
P <- 1000
K <- 7
B <- 3

group_means <- rnorm(K)
group_sds <- rgamma(K, 4, 2)
group_dfs <- c(4, 10, 20, 50, 10, 20, 20)

batch_shift <- rnorm(B, sd = 0.5)
batch_scale <- 1 + rgamma(B, 10, 20)

group_weights <- c(0.3, 0.2, 0.1, rep(0.4 / 4, 4))
batch_weights <- c(0.4, 0.4, 0.2)

gen_data <- generateBatchData(N, P, 
                  group_means = group_means, 
                  group_std_devs = group_sds,
                  batch_shift = batch_shift, 
                  batch_scale = batch_scale,
                  group_weights = group_weights, 
                  batch_weights = batch_weights, 
                  group_dfs = group_dfs)

names(gen_data$group_IDs) <- row.names(gen_data$observed_data)

# frac_present <- 0.2
# P_2 <- 25
# K_2 <- 15
# class_weights_2 <- rgamma(K_2, 10)
# class_weights_2 <- class_weights_2 / sum(class_weights_2)
# params_2 <- list(means = rnorm(K_2, sd = 1),
#                  std_devs = rgamma(K_2, 50, 50)
#                  )
# gen_view_2 <- generateViewGivenStructure(gen_data$group_IDs,
#                                          frac_present = frac_present, 
#                                          P = P_2, 
#                                          K = K_2, 
#                                          type = "G", 
#                                          class_weights = class_weights_2,
#                                          params = params_2)
# 
# frac_present <- 0.2
# P_3 <- 20
# K_3 <- 10
# class_weights_3 <- rgamma(K_3, 10, 10)
# class_weights_3 <- class_weights_3 / sum(class_weights_3)
# probs <- rbeta(K_3, 1, 1)
# params_3 <- list(probs = probs)
# 
# gen_view_3 <- generateViewGivenStructure(gen_data$group_IDs,
#                                          frac_present = frac_present, 
#                                          P = P_3, 
#                                          K = K_3, 
#                                          type = "C", 
#                                          class_weights = class_weights_3,
#                                          params = params_3)


annotatedHeatmap(gen_data$observed_data, gen_data$group_IDs, show_rownames = FALSE)
annotatedHeatmap(gen_data$corrected_data, gen_data$group_IDs, show_rownames = FALSE)

annotatedHeatmap(gen_view_2$gen_data, gen_view_2$labels, show_rownames = FALSE)
annotatedHeatmap(gen_view_3$gen_data, gen_view_3$labels, show_rownames = FALSE, col_pal = simColPal(), defineBreaks(simColPal(), lb = 0))


R <- 1000
thin <- 25
type <- "G"
K_max <- 25
fixed <- gen_data$fixed
initial_labels <- gen_data$group_IDs

mcmc_semi <- callMixtureModel(gen_data$observed_data,
                 R = R,
                 thin = thin,
                 type = type, 
                 K = K_max,
                 fixed = fixed, 
                 initial_labels = initial_labels)

mcmc_un <- callMixtureModel(gen_data$observed_data,
                         R = R,
                         thin = thin,
                         type = type, 
                         K = K_max,
                         fixed = rep(0, N), 
                         initial_labels = initial_labels)

psm1 <- mcmc_semi$allocations |> createSimilarityMat()
row.names(psm1) <- row.names(gen_data$observed_data)

psm2 <- mcmc_un$allocations |> createSimilarityMat()
row.names(psm2) <- row.names(gen_data$observed_data)

annotatedHeatmap(psm1, gen_data$group_IDs, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(psm2, gen_data$group_IDs, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
