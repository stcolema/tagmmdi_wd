
library(tagmReDraft)
library(ggplot2)
library(tidyr)
library(magrittr)
library(mdiHelpR)
library(pheatmap)

annotatedData <- function(X, v, filename) {
  K_used <- 1
  for(ii in seq(1, n_files)) {
    .pred <- data.frame(X = factor(paste0("Cluster ", mcmc_out[[ii]]$predicted_clustering[[v]])))
    K_used <- max(K_used, length(unique(.pred$X)))
    # K_used[ii] <- length(unique(.pred$X))
    colnames(.pred) <- paste0("Chain ", ii)
    if(ii == 1) {
      anno_row = .pred
    } else {
      anno_row = cbind(anno_row, .pred)
    }
  }
  
  row.names(anno_row) <- row.names(X)
  
  ann_colours <- lapply(seq(1, n_files), function(x) ggthemes::colorblind_pal()(K_used))
  names(ann_colours) <- paste("Chain", seq(1, n_files))
  for(ii in seq(1, n_files)) {
    names(ann_colours[[ii]]) <- paste0("Cluster ", seq(1, K_used))
  }
  
  # Create the heatmap
  pheatmap(X,
           color = dataColPal(),
           # breaks = my_breaks,
           annotation_row = anno_row,
           annotation_colors = ann_colours,
           main = paste0("View ", v, " annotated by predicted clusterings"),
           filename = filename
  )
}

findFusedGenes <- function(x, y, threshold = 0.5) {
  colMeans(x == y) > threshold
}

whichGenesFused <- function(x, y, threshold = 0.5) {
  which(findFusedGenes(x, y, threshold))
}

findFusedGenesAcrossAllViews <- function(mcmc, threshold = 0.5) {
  V <- mcmc$V
  N <- mcmc$N
  n_pairs <- choose(V, 2)
  # fused_genes <- list()
  fused_genes <- matrix(0, nrow = N, ncol = n_pairs)
  count <- 0
  for(ii in seq(1, V - 1)) {
    for(jj in seq(ii + 1, V)) {
      count <- count + 1
      fused_genes[, count] <- findFusedGenes(mcmc$allocations[ , , ii], mcmc$allocations[ , , jj], threshold)
      # fused_genes[[count]] <- findFusedGenes(mcmc$allocations[ , , ii], mcmc$allocations[ , , jj], threshold)
    }
  }
  fused_genes
}

out_dir <- "./MDITestData/Output/"

mditest1 <- read.csv("./MDITestData/MDItestdata1.csv", row.names = 1)
mditest2 <- read.csv("./MDITestData/MDItestdata2.csv", row.names = 1)
mditest3 <- read.csv("./MDITestData/MDItestdata3.csv", row.names = 1)
mditest4 <- read.csv("./MDITestData/MDItestdata4.csv", row.names = 1)
gene_names <- row.names(mditest1)

files <- list.files(out_dir, pattern = "mcmcMDItestdataTestFrac08Seed*", full.names = TRUE)
n_files <- length(files)

mcmc_out <- list()
burn <- 2e4

psms <- vector("list", 4)

for(ii in seq(1, n_files)) {
  f <- files[ii]
  .mcmc <- readRDS(f)[[1]]
  # .mcmc$V <- 
  mcmc_out[[ii]] <- processMCMCChain(.mcmc, burn)
  V <- mcmc_out[[ii]]$V
  n_pairs <- choose(V, 2)
  mcmc_out[[ii]]$psm <- list()
  mcmc_out[[ii]]$predicted_clustering <- list()
  for(v in seq(1, V)) {
    mcmc_out[[ii]]$psm[[v]] <- psms[[v]][[ii]] <- psm <- mdiHelpR::makePSM(mcmc_out[[ii]]$allocations[, , v]) |> 
      set_rownames(gene_names) |> 
      set_colnames(gene_names)
    mcmc_out[[ii]]$predicted_clustering[[v]] <- pred_cl <- mcclust::maxpear(psm)$cl
  }
  
  .phis <- data.frame(mcmc_out[[ii]]$phis)
  colnames(.phis) <- c("Phi_12", "Phi_13", "Phi_14", "Phi_23", "Phi_24", "Phi_34")
  .phis$Chain <- ii
  if(ii == 1) {
    phi_df <- .phis
  } else {
    phi_df <- rbind(phi_df, .phis)
  }
}

phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  dplyr::filter(Chain == 9, Sampled_value < 125) |> 
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Parameter~Chain, scales = "free_y")


phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Parameter~Chain, scales = "free_y")

phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30) |>  #, Parameter != "Phi_12") |>
  dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Chain ~ Parameter, scales = "free_y")

V <- mcmc_out[[1]]$V
n_phis <- choose(V, 2)

phi_medians <- lapply(mcmc_out, function(x) apply(x$phis, 2, median)) |> 
  unlist() |> 
  matrix(nrow = n_files, byrow = TRUE)

phi_means <- lapply(mcmc_out, function(x) colMeans(x$phis)) |> 
  unlist() |> 
  matrix(nrow = n_files, byrow = TRUE)

phis <- lapply(mcmc_out, function(x) x$phis) |> 
  unlist() |> 
  matrix(ncol = n_phis, byrow = FALSE)

boxplot(phis)
boxplot(phi_means)
boxplot(phi_medians)

fused_genes_1 <- findFusedGenesAcrossAllViews(mcmc_out[[1]])
fused_genes_2 <- findFusedGenesAcrossAllViews(mcmc_out[[2]])
fused_genes_3 <- findFusedGenesAcrossAllViews(mcmc_out[[3]])

mcmc_out[[1]]$predicted_clustering

mditest1 <- read.csv("./MDITestData/MDItestdata1.csv", row.names = 1)
mditest2 <- read.csv("./MDITestData/MDItestdata2.csv", row.names = 1)
mditest3 <- read.csv("./MDITestData/MDItestdata3.csv", row.names = 1)
mditest4 <- read.csv("./MDITestData/MDItestdata4.csv", row.names = 1)

K_used <- 1
for(ii in seq(1, n_files)) {
  .pred <- data.frame(X = factor(paste0("Cluster ", mcmc_out[[ii]]$predicted_clustering[[1]])))
  K_used <- max(K_used, length(unique(.pred$X)))
  # K_used[ii] <- length(unique(.pred$X))
  colnames(.pred) <- paste0("Chain ", ii)
  if(ii == 1) {
    anno_row = .pred
  } else {
    anno_row = cbind(anno_row, .pred)
  }
}

row.names(anno_row) <- row.names(mditest1)

ann_colours <- lapply(seq(1, n_files), function(x) ggthemes::colorblind_pal()(K_used))
names(ann_colours) <- paste("Chain", seq(1, n_files))
for(ii in seq(1, n_files)) {
  names(ann_colours[[ii]]) <- paste0("Cluster ", seq(1, K_used))
}

# Create the heatmap
pheatmap(mditest1,
                         color = dataColPal(),
                         # breaks = my_breaks,
                         annotation_row = anno_row,
                         annotation_colors = ann_colours,
                         main = "View 1 annotated by predicted clusterings")


annotatedData(mditest1, 1, filename = "MDItestdata1_all_chains.png")
annotatedData(mditest2, 2, filename = "MDItestdata2_all_chains.png")
annotatedData(mditest3, 3, filename = "MDItestdata3_all_chains.png")
annotatedData(mditest4, 4, filename = "MDItestdata4_all_chains.png")

annotatedHeatmap(mditest1, mcmc_out[[1]]$predicted_clustering[[1]], main = "View 1 annotated by chain 1", filename = "MDItestdata1_ch1.png")
annotatedHeatmap(mditest1, mcmc_out[[2]]$predicted_clustering[[1]], main = "View 2 annotated by chain 1", filename = "MDItestdata2_ch1.png")
annotatedHeatmap(mditest1, mcmc_out[[3]]$predicted_clustering[[1]], main = "View 3 annotated by chain 1", filename = "MDItestdata3_ch1.png")

annotatedHeatmap(mditest2, mcmc_out[[1]]$predicted_clustering[[2]], main = "View 1 annotated by chain 2", filename = "MDItestdata1_ch2.png")
annotatedHeatmap(mditest2, mcmc_out[[2]]$predicted_clustering[[2]], main = "View 2 annotated by chain 2", filename = "MDItestdata2_ch2.png")
annotatedHeatmap(mditest2, mcmc_out[[3]]$predicted_clustering[[2]], main = "View 3 annotated by chain 2", filename = "MDItestdata3_ch2.png")

annotatedHeatmap(mditest3, mcmc_out[[1]]$predicted_clustering[[3]], main = "View 1 annotated by chain 3", filename = "MDItestdata1_ch3.png")
annotatedHeatmap(mditest3, mcmc_out[[2]]$predicted_clustering[[3]], main = "View 2 annotated by chain 3", filename = "MDItestdata2_ch3.png")
annotatedHeatmap(mditest3, mcmc_out[[3]]$predicted_clustering[[3]], main = "View 3 annotated by chain 3", filename = "MDItestdata3_ch3.png")

annotatedHeatmap(mditest4, mcmc_out[[1]]$predicted_clustering[[4]])
annotatedHeatmap(mditest4, mcmc_out[[2]]$predicted_clustering[[4]])
annotatedHeatmap(mditest4, mcmc_out[[3]]$predicted_clustering[[4]])

annotatedHeatmap(mcmc_out[[1]]$psm[[4]], mcmc_out[[1]]$predicted_clustering[[4]], col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[2]]$psm[[4]], mcmc_out[[2]]$predicted_clustering[[4]], col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[3]]$psm[[4]], mcmc_out[[3]]$predicted_clustering[[4]], col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))

psms_v1 <- prepSimilarityMatricesForGGplot(psms[[1]])
psms_v2 <- prepSimilarityMatricesForGGplot(psms[[2]])
psms_v3 <- prepSimilarityMatricesForGGplot(psms[[3]])
psms_v4 <- prepSimilarityMatricesForGGplot(psms[[4]])

psms_v1$View <- 1
psms_v2$View <- 2
psms_v3$View <- 3
psms_v4$View <- 4

psm_df <- rbind(psms_v1, psms_v2, psms_v3, psms_v4)

mdiHelpR::setMyTheme()
psm_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() + 
  facet_grid(View~Chain, labeller = label_both) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#FDF9DC", colour = "#FDF9DC")
  )

psms_v1 |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() + 
  facet_wrap(~Chain) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")

psms_v1 |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() + 
  facet_wrap(~Chain) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")

mcclust::arandi(mcmc_out[[1]]$predicted_clustering[[1]], mcmc_out[[1]]$predicted_clustering[[2]])
mcclust::arandi(mcmc_out[[2]]$predicted_clustering[[1]], mcmc_out[[2]]$predicted_clustering[[2]])
mcclust::arandi(mcmc_out[[3]]$predicted_clustering[[1]], mcmc_out[[3]]$predicted_clustering[[2]])

mcclust::arandi(mcmc_out[[1]]$predicted_clustering[[1]], mcmc_out[[1]]$predicted_clustering[[3]])
mcclust::arandi(mcmc_out[[2]]$predicted_clustering[[1]], mcmc_out[[2]]$predicted_clustering[[3]])
mcclust::arandi(mcmc_out[[3]]$predicted_clustering[[1]], mcmc_out[[3]]$predicted_clustering[[3]])

mcclust::arandi(mcmc_out[[1]]$predicted_clustering[[1]], mcmc_out[[1]]$predicted_clustering[[4]])
mcclust::arandi(mcmc_out[[2]]$predicted_clustering[[1]], mcmc_out[[2]]$predicted_clustering[[4]])
mcclust::arandi(mcmc_out[[3]]$predicted_clustering[[1]], mcmc_out[[3]]$predicted_clustering[[4]])

mcclust::arandi(mcmc_out[[1]]$predicted_clustering[[2]], mcmc_out[[1]]$predicted_clustering[[3]])
mcclust::arandi(mcmc_out[[2]]$predicted_clustering[[2]], mcmc_out[[2]]$predicted_clustering[[3]])
mcclust::arandi(mcmc_out[[3]]$predicted_clustering[[2]], mcmc_out[[3]]$predicted_clustering[[3]])

mcclust::arandi(mcmc_out[[1]]$predicted_clustering[[2]], mcmc_out[[1]]$predicted_clustering[[4]])
mcclust::arandi(mcmc_out[[2]]$predicted_clustering[[2]], mcmc_out[[2]]$predicted_clustering[[4]])
mcclust::arandi(mcmc_out[[3]]$predicted_clustering[[2]], mcmc_out[[3]]$predicted_clustering[[4]])

mcclust::arandi(mcmc_out[[1]]$predicted_clustering[[3]], mcmc_out[[1]]$predicted_clustering[[4]])
mcclust::arandi(mcmc_out[[2]]$predicted_clustering[[3]], mcmc_out[[2]]$predicted_clustering[[4]])
mcclust::arandi(mcmc_out[[3]]$predicted_clustering[[3]], mcmc_out[[3]]$predicted_clustering[[4]])

ari_lst <- list()
phi_lst <- list()
for(ii in seq(1, n_files)) {
  count <- 0
  ari_lst[[ii]] <- matrix(1, nrow = V, ncol = V)
  phi_lst[[ii]] <- matrix(NA, nrow = V, ncol = V)
  for(vv in seq(1, V - 1)) {
    for(uu in seq(vv + 1, V)) {
      ari_lst[[ii]][uu, vv] <- ari_lst[[ii]][vv, uu] <- mcclust::arandi(mcmc_out[[ii]]$predicted_clustering[[uu]], mcmc_out[[ii]]$predicted_clustering[[vv]]) 
      count <- count + 1
      phi_lst[[ii]][uu, vv] <- phi_lst[[ii]][vv, uu] <- mean(mcmc_out[[ii]]$phis[, count])
    }
  }
}

ari_lst[[1]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[2]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[3]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

phi_lst[[1]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[2]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[3]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
