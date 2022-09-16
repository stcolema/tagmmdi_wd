
library(tagmReDraft)
library(ggplot2)
library(tidyr)
library(magrittr)
library(mdiHelpR)
library(pheatmap)
library(patchwork)

mdiHelpR::setMyTheme()

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

out_dir <- "./MDITestData/Output/Unsupervised/"

mditest1 <- read.csv("./MDITestData/MDItestdata1.csv", row.names = 1)
mditest2 <- read.csv("./MDITestData/MDItestdata2.csv", row.names = 1)
mditest3 <- read.csv("./MDITestData/MDItestdata3.csv", row.names = 1)
mditest4 <- read.csv("./MDITestData/MDItestdata4.csv", row.names = 1)
gene_names <- row.names(mditest1)


row_order1 <- findOrder(mditest1)
row_order2 <- findOrder(mditest2)
row_order3 <- findOrder(mditest3)
row_order4 <- findOrder(mditest4)


labels <- c(
  rep(1, 15),
  rep(2, 35),
  rep(3, 6),
  rep(4, 8),
  rep(5, 11),
  rep(6, 8),
  rep(7, 17)
)

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
burn <- 15000

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

across_chain <- predictFromMultipleChains(mcmc_out, burn = 0, chains_already_processed = T)
across_chain$psm <- list()
across_chain$psm[[1]] <- mdiHelpR::createSimilarityMat(across_chain$allocations[[1]]) |> 
  set_rownames(gene_names)
across_chain$psm[[2]] <- mdiHelpR::createSimilarityMat(across_chain$allocations[[2]]) |> 
  set_rownames(gene_names)
across_chain$psm[[3]] <- mdiHelpR::createSimilarityMat(across_chain$allocations[[3]]) |> 
  set_rownames(gene_names)
across_chain$psm[[4]] <- mdiHelpR::createSimilarityMat(across_chain$allocations[[4]]) |> 
  set_rownames(gene_names)

annotatedHeatmap(across_chain$psm[[1]], labels, 
                 main = "Across chain PSM view 1 annotated by true labels",
                 col_pal = simColPal(), 
                 my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(across_chain$psm[[2]], labels2,
                 main = "Across chain PSM view 2 annotated by true labels",
                 col_pal = simColPal(), 
                 my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(across_chain$psm[[3]], labels3, 
                 main = "Across chain PSM view 3 annotated by true labels",
                 col_pal = simColPal(), 
                 my_breaks = defineBreaks(simColPal(), lb = 0))
annotatedHeatmap(across_chain$psm[[4]], labels4, 
                 main = "Across chain PSM view 4 annotated by true labels",
                 col_pal = simColPal(), 
                 my_breaks = defineBreaks(simColPal(), lb = 0))

across_chain$phis |> boxplot(main = "Sampled phis across 8 chains")

across_chain$phis |> 
  as.data.frame() |> 
  pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Chain == 9, Sampled_value < 125) |> 
  ggplot(aes(x = Sampled_value, y = Parameter, fill = Parameter)) +
  geom_boxplot()
# +
  # facet_wrap(~Parameter, scales = "free_y")

across_chain_psm_df <- tagmReDraft::prepSimilarityMatricesForGGplot(across_chain$psm)

phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Chain == 9, Sampled_value < 125) |> 
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
  # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Chain ~ Parameter, scales = "free")

# phi_12 <- phi_df |> 
#   pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
#   dplyr::filter(Parameter == "Phi_12") |>
#   # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
#   ggplot(aes(x = Sampled_value, group = Chain,)) +
#   geom_density(alpha = 0.3, fill = "#000000") +
#   facet_wrap(~Chain, ncol = 1, scales = "free_y") + 
#   labs(subtitle = expression(phi[12]))
# 
# phi_13 <- phi_df |> 
#   pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
#   dplyr::filter(Parameter == "Phi_13") |>
#   # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
#   ggplot(aes(x = Sampled_value, group = Chain,)) +
#   geom_density(alpha = 0.3, fill = "#E69F00") +
#   facet_wrap(~Chain, ncol = 1, scales = "free_y") + 
#   labs(subtitle = expression(phi[13]))
# 
# phi_14 <- phi_df |> 
#   pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
#   dplyr::filter(Parameter == "Phi_14") |>
#   # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
#   ggplot(aes(x = Sampled_value, group = Chain,)) +
#   geom_density(alpha = 0.3, fill = "#56B4E9") +
#   facet_wrap(~Chain, ncol = 1, scales = "free_y") + 
#   labs(subtitle = expression(phi[14]))
# 
# phi_23 <- phi_df |> 
#   pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
#   dplyr::filter(Parameter == "Phi_23") |>
#   # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
#   ggplot(aes(x = Sampled_value, group = Chain,)) +
#   geom_density(alpha = 0.3, fill = "#009E73") +
#   facet_wrap(~Chain, ncol = 1, scales = "free_y") + 
#   labs(subtitle = expression(phi[23]))
# 
# phi_24 <- phi_df |> 
#   pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
#   dplyr::filter(Parameter == "Phi_24") |>
#   # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
#   ggplot(aes(x = Sampled_value, group = Chain,)) +
#   geom_density(alpha = 0.3, fill = "#F0E442") +
#   facet_wrap(~Chain, ncol = 1, scales = "free_y") + 
#   labs(subtitle = expression(phi[24]))
# 
# phi_34 <- phi_df |> 
#   pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
#   dplyr::filter(Parameter == "Phi_34") |>
#   # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
#   ggplot(aes(x = Sampled_value, group = Chain,)) +
#   geom_density(alpha = 0.3, fill = "#0072B2") +
#   facet_wrap(~Chain, ncol = 1, scales = "free_y") + 
#   labs(subtitle = expression(phi[34]))
# 
# phi_12 + phi_13 + phi_14 + phi_23 + phi_24 + phi_34 +
#   plot_annotation(title = "Sampled phi values across chains")
# ggsave("./Sampled_phi_values_across_chains_semisupervised.png", width = 16, height = 11)


phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30, Parameter != "Phi_12") |>
  # dplyr::filter(Chain %in% c(1, 5, 7)) |>
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Parameter~Chain, scales = "free_y")

phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30) |>  #, Parameter != "Phi_12") |>
  # dplyr::filter(Chain %in% c(1, 5, 7)) |>
  ggplot(aes(x = Sampled_value, fill = Parameter)) +
  geom_density() +
  facet_grid(Chain ~ Parameter, scales = "free_y")

phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30) |>  #, Parameter != "Phi_12") |>
  # dplyr::filter(Chain %in% c(1, 5, 7)) |>
  ggplot(aes(x = Sampled_value, y = Parameter, fill = Parameter)) +
  geom_boxplot() +
  facet_wrap(~ Chain, labeller = label_both)

ggsave("./MDITestData/Output/Unsupervised/all_phis.png", width = 9, height = 6)

phi_df |> 
  pivot_longer(-Chain, names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30) |>  #, Parameter != "Phi_12") |>
  # dplyr::filter(Chain %in% c(1, 5, 7)) |>
  ggplot(aes(x = Sampled_value, y = Parameter, fill = Parameter)) +
  geom_boxplot() +
  facet_wrap(~ Chain, labeller = label_both) +
  coord_cartesian(xlim=c(0, 2.5))

# ggsave("../Phi_naive_calc.png")
# 
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

ggsave("./MDITestData/Output/Unsupervised/all_psms.png", width = 9, height = 6)

across_chain_psm_df$View <- across_chain_psm_df$Chain
across_chain_psm_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() + 
  facet_wrap(~View, labeller = label_both) + 
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
ari_lst[[4]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[5]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[6]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[7]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[8]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
ari_lst[[9]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

phi_lst[[1]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[2]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[3]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[4]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[5]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[6]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[7]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[8]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)
phi_lst[[9]] |> pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

