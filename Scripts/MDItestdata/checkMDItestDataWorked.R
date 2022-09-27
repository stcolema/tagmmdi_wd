
library(tagmReDraft)
library(ggplot2)
library(tidyr)
library(magrittr)
library(mdiHelpR)
library(pheatmap)

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

out_dir <- "./MDITestData/Output/Unsupervised//"

labels <- c(
  rep(1, 15),
  rep(2, 35),
  rep(3, 6),
  rep(4, 8),
  rep(5, 11),
  rep(6, 8),
  rep(7, 17)
)

mditest1 <- read.csv("./MDITestData/MDItestdata1.csv", row.names = 1)
mditest2 <- read.csv("./MDITestData/MDItestdata2.csv", row.names = 1)
mditest3 <- read.csv("./MDITestData/MDItestdata3.csv", row.names = 1)
mditest4 <- read.csv("./MDITestData/MDItestdata4.csv", row.names = 1)
gene_names <- row.names(mditest1)


row_order1 <- findOrder(mditest1)
row_order2 <- findOrder(mditest2)
row_order3 <- findOrder(mditest3)
row_order4 <- findOrder(mditest4)

labels2a <- c(
  rep(1, 15),
  rep(4, 8),
  rep(3, 6),
  rep(5, 11),
  rep(6, 8),
  rep(7, 17),
  rep(2, 35)
)

labels2 <- labels2a[match(row.names(mditest2), row.names(mditest2[row_order2, ]))]
labels3 <- labels2a[match(row.names(mditest3), row.names(mditest3[row_order3, ]))]
labels4 <- labels2a[match(row.names(mditest4), row.names(mditest4[row_order4, ]))]

files <- list.files(out_dir, pattern = "*rds", full.names = TRUE)
n_files <- length(files)

mcmc_out <- list()
burn <- 15000 # 1e4

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
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Parameter~Chain, scales = "free_y")

phi_df |> 
  # filter(Chain == c(6, 8)) |> 
  # filter(Chain == c(2, 3, 4, 6, 8, 9)) |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  # xlim(c(0, 50)) +
  facet_grid(Chain ~ Parameter, scales = "free")

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

annotatedHeatmap(mcmc_out[[1]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[2]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[3]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[4]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[5]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[6]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[7]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[8]]$psm[[1]], labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))

annotatedHeatmap(mcmc_out[[1]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[2]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[3]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[4]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[5]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[6]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[7]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[8]]$psm[[2]], labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))

annotatedHeatmap(mcmc_out[[1]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[2]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[3]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[4]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[5]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[6]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[7]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[8]]$psm[[3]], labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))

annotatedHeatmap(mcmc_out[[1]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[2]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[3]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[4]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[5]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[6]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[7]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(mcmc_out[[8]]$psm[[4]], labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0))



psms_v1 <- prepSimilarityMatricesForGGplot(psms[[1]], matrix_setting_order = 8)
psms_v2 <- prepSimilarityMatricesForGGplot(psms[[2]], matrix_setting_order = 8)
psms_v3 <- prepSimilarityMatricesForGGplot(psms[[3]], matrix_setting_order = 8)
psms_v4 <- prepSimilarityMatricesForGGplot(psms[[4]], matrix_setting_order = 8)

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
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")

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

