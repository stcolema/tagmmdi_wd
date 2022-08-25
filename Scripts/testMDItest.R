
library(tagmReDraft)
library(mdiHelpR)
library(magrittr)
library(ggplot2)
library(tidyr)
mdiHelpR::setMyTheme()
set.seed(1)

labels <- c(
  rep(1, 15),
  rep(2, 35),
  rep(3, 6),
  rep(4, 8),
  rep(5, 11),
  rep(6, 8),
  rep(7, 17)
)

mditest1 <- read.csv("/home/stephen/Documents/PhD/tagmmdi_wd/MDITestData/MDItestdata1.csv", row.names = 1)
mditest2 <- read.csv("/home/stephen/Documents/PhD/tagmmdi_wd/MDITestData/MDItestdata2.csv", row.names = 1)
mditest3 <- read.csv("/home/stephen/Documents/PhD/tagmmdi_wd/MDITestData/MDItestdata3.csv", row.names = 1)
mditest4 <- read.csv("/home/stephen/Documents/PhD/tagmmdi_wd/MDITestData/MDItestdata4.csv", row.names = 1)

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

annotatedHeatmap(mditest1, labels, cluster_rows = T)
annotatedHeatmap(mditest2, labels2, cluster_rows = T)
annotatedHeatmap(mditest3, labels3, cluster_rows = T)
annotatedHeatmap(mditest4, labels4, cluster_rows = T)

data_modelled <- list(
  as.matrix(mditest1),
  as.matrix(mditest2),
  as.matrix(mditest3),
  as.matrix(mditest4)
)

N <- nrow(mditest1)
V <- length(data_modelled)

initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)
initial_labels_1 <- labels # generateInitialUnsupervisedLabels(100, 1, 25)
initial_labels_2 <- labels2 # generateInitialUnsupervisedLabels(100, 1, 25)
initial_labels_3 <- labels3 # generateInitialUnsupervisedLabels(100, 1, 25)
initial_labels_4 <- labels4 # generateInitialUnsupervisedLabels(100, 1, 25)

initial_labels[, 1] <- initial_labels_1
initial_labels[, 2] <- initial_labels_2
initial_labels[, 3] <- initial_labels_3
initial_labels[, 4] <- initial_labels_4

fixed <- matrix(0, N, V)

type <- rep("GP", V)
R <- 5000
thin <- 25
burn <- 1500
K_max <- 35
K <- rep(K_max, V)

# gp_mix_mod <- tagmReDraft::callMixtureModel(data_modelled[[1]],
#                                R = 5000,
#                                thin = 25,
#                                type = "GP",
#                                K = K[1]
# )
# 
# psm_mix_mod <- gp_mix_mod$allocations |>
#   createSimilarityMat() |>
#   set_rownames(row.names(mditest1))
# 
# annotatedHeatmap(psm_mix_mod, labels, my_breaks = defineBreaks(simColPal(), lb = 0), col_pal = simColPal())

gp_mix <- tagmReDraft::callMDI(
  data_modelled,
  R = R,
  thin = thin,
  type = type,
  K = K,
  fixed = fixed,
  initial_labels = initial_labels, 
  initial_labels_as_intended = FALSE
)

gp_mix$N_k

# gp_mix <- tagmReDraft::callMixtureModel(as.matrix(mditest1),
#                                R = 10000, 
#                                thin = 25, 
#                                type = "GP",
#                                K = 25,
#                                fixed = rep(0, 100),
#                                initial_labels = initial_labels_1
#                                
# )

new_gp <- processMCMCChain(gp_mix, burn)

new_gp$phis |> boxplot()
new_gp$mass |> boxplot()

new_gp$weights[, , 1] |> boxplot(main = "Component weights view 1")
new_gp$weights[, , 2] |> boxplot(main = "Component weights view 2")
new_gp$weights[, , 3] |> boxplot(main = "Component weights view 3")
new_gp$weights[, , 4] |> boxplot(main = "Component weights view 4")

psm1 <- new_gp$allocations[, , 1] |> createSimilarityMat()
psm2 <- new_gp$allocations[, , 2] |> createSimilarityMat()
psm3 <- new_gp$allocations[, , 3] |> createSimilarityMat()
psm4 <- new_gp$allocations[, , 4] |> createSimilarityMat()

row.names(psm1) <- row.names(psm2) <- row.names(psm3) <- row.names(psm4) <- row.names(mditest1)

annotatedHeatmap(psm1, labels, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0), main = "Synthetic Data 1: annotated by true labels")
annotatedHeatmap(psm2, labels2, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0), main = "Synthetic Data 2: annotated by true labels")
annotatedHeatmap(psm3, labels3, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0), main = "Synthetic Data 3: annotated by true labels")
annotatedHeatmap(psm4, labels4, col_pal = simColPal(), my_breaks = defineBreaks(simColPal(), lb = 0), main = "Synthetic Data 4: annotated by true labels")

phis <- data.frame(new_gp$phis)
colnames(phis) <- c("Phi_12", "Phi_13", "Phi_14", "Phi_23", "Phi_24", "Phi_34")

phis |> 
  pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_value") |> 
  # dplyr::filter(Sampled_value < 30) |>  #, Parameter != "Phi_12") |>
  # dplyr::filter(Chain %in% c(1, 5, 7)) |>
  ggplot(aes(x = Sampled_value, y = Parameter, fill = Parameter)) +
  geom_boxplot() +
  ggthemes::scale_fill_colorblind()

gamma_names <- paste0("Gamma_", seq(1, K_max))
for(v in seq(1, V)) {
  w_df <- new_gp$weights[, , v] |> 
    set_colnames(gamma_names) |> 
    as.data.frame()
  w_df$Iteration <- seq(0, R, thin)
  w_df <- w_df |> 
    pivot_longer(-Iteration, names_to = "Parameter", values_to = "Sampled_value")
  w_df$View <- v
  if(v == 1){
    weight_df <- w_df
  } else {
    weight_df <- rbind(weight_df, w_df)
  }
}

weight_df$Parameter <- factor(weight_df$Parameter, levels = gamma_names)

weight_df |> 
  ggplot(aes(x = Sampled_value, y = Parameter, fill = Parameter)) +
  geom_boxplot() +
  # ggthemes::scale_fill_colorblind() +
  facet_wrap(~View) 
