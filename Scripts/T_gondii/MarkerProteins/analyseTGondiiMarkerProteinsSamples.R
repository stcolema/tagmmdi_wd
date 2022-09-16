

library(tidyverse)
library(tagmReDraft)
library(magrittr)
library(mdiHelpR)

set.seed(1)
setMyTheme()

my_dir <- "./T_gondii/MarkersOnlyOutput/"
files <- list.files(my_dir, full.names = TRUE, pattern = "90000.rds$")
n_files <- length(files)

burn <- 20000

mcmc_chains <- list()
for(ii in seq(1, n_files)) {
  f <- files[ii]
  x <- readRDS(f)[[1]]
  mcmc_chains[[ii]] <- processMCMCChain(x, burn, construct_psm = TRUE)
  mcmc_chains[[ii]]$Chain <- ii
}

# new_mcmc <- tagmReDraft::processMCMCChains(x, burn, point_estimate_method = "median", construct_psm = TRUE)

mcmc_chains[[1]]$phis |> boxplot()
mcmc_chains[[2]]$phis |> boxplot()
mcmc_chains[[3]]$phis |> boxplot()

fused_12_chain_1 <- which(colMeans(mcmc_chains[[1]]$allocations[, , 1] == mcmc_chains[[1]]$allocations[, , 2]) > 0.5)
fused_13_chain_1 <- which(colMeans(mcmc_chains[[1]]$allocations[, , 1] == mcmc_chains[[1]]$allocations[, , 3]) > 0.5)
fused_23_chain_1 <- which(colMeans(mcmc_chains[[1]]$allocations[, , 2] == mcmc_chains[[1]]$allocations[, , 3]) > 0.5)

fused_12_chain_2 <- which(colMeans(mcmc_chains[[2]]$allocations[, , 1] == mcmc_chains[[2]]$allocations[, , 2]) > 0.5)
fused_13_chain_2 <- which(colMeans(mcmc_chains[[2]]$allocations[, , 1] == mcmc_chains[[2]]$allocations[, , 3]) > 0.5)
fused_23_chain_2 <- which(colMeans(mcmc_chains[[2]]$allocations[, , 2] == mcmc_chains[[2]]$allocations[, , 3]) > 0.5)

fused_12_chain_3 <- which(colMeans(mcmc_chains[[3]]$allocations[, , 1] == mcmc_chains[[3]]$allocations[, , 2]) > 0.5)
fused_13_chain_3 <- which(colMeans(mcmc_chains[[3]]$allocations[, , 1] == mcmc_chains[[3]]$allocations[, , 3]) > 0.5)
fused_23_chain_3 <- which(colMeans(mcmc_chains[[3]]$allocations[, , 2] == mcmc_chains[[3]]$allocations[, , 3]) > 0.5)

fused_123 <- fused_12[fused_12 %in% fused_13]

psms_v1 <- mcmc_chains |> lapply(function(x) {x$psm[[1]]}) |> 
  prepSimilarityMatricesForGGplot()
psms_v2 <- mcmc_chains |> lapply(function(x) {x$psm[[2]]}) |> 
  prepSimilarityMatricesForGGplot()
psms_v3 <- mcmc_chains |> lapply(function(x) {x$psm[[3]]}) |> 
  prepSimilarityMatricesForGGplot()

psms_v1$View <- 1
psms_v2$View <- 2
psms_v3$View <- 3

psm_df <- rbind(psms_v1, psms_v2, psms_v3)

psms_v1 <- NULL
psms_v2 <- NULL
psms_v3 <- NULL

psm_df |> 
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() +
  facet_grid(View~Chain) +
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")

inputs <- readRDS("./T_gondii/PreparedData/TGondiiMDI_K_125_input_marker_proteins.rds")

proteins <- inputs$data_modelled[[1]] |> row.names()


fused_cl <- mcmc_chains[[1]]$pred[[1]][fused_123]
names(fused_cl) <- proteins[fused_123]
library(pRolocdata)
data("Barylyuk2020ToxoLopit")

fData(Barylyuk2020ToxoLopit[proteins[fused_123], ])$"markers" |> 
  table()

microneme_cl <- fused_cl[fused_cl == 15]
apicoplast_cl <- fused_cl[fused_cl == 7]

fused_cl2 <- mcmc_chains[[1]]$pred[[2]][fused_123]
names(fused_cl2) <- proteins[fused_123]
fused_cl3 <- mcmc_chains[[1]]$pred[[3]][fused_123]
names(fused_cl3) <- proteins[fused_123]

fused_cl2[fused_cl == 15]
fused_cl2[fused_cl == 7]

fused_cl3[fused_cl == 15]
fused_cl3[fused_cl == 7]

inputs$data_modelled[[2]][names(fused_cl[fused_cl == 15]), ] |> 
  pheatmap::pheatmap(show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE)
inputs$data_modelled[[3]][names(fused_cl[fused_cl == 15]), ] |>
  pheatmap::pheatmap(show_colnames = FALSE, show_rownames = FALSE)

inputs$data_modelled[[2]][names(fused_cl[fused_cl == 7]), ] |> 
  pheatmap::pheatmap(show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE)
inputs$data_modelled[[3]][names(fused_cl[fused_cl == 7]), ] |> 
  pheatmap::pheatmap(show_colnames = T, show_rownames = FALSE)


inputs$data_modelled[[3]][names(fused_cl[fused_cl == 7]), findOrder(t(inputs$data_modelled[[3]][names(fused_cl[fused_cl == 7]), ] ))] |> colnames()
