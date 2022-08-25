
library(tagmReDraft)
library(mdiHelpR)
library(ggplot2)
library(tidyr)
library(patchwork)

setMyTheme()

out_dir <- "Simulations/OutputPhiModel/"

scenarios <- list.dirs(out_dir, full.names = F, recursive = F)
scn_dirs <- paste0(out_dir, scenarios)[c(1,2, 5)]

files <- lapply(scn_dirs, list.files, ".rds$", F, T) |> unlist()

k_6_files <- lapply(scn_dirs, list.files, "*K_6_*", F, T) |> unlist()

n_files <- length(files)

V <- 3
phi_names <- c("Phi_12", "Phi_13", "Phi_23")

for(ii in seq(1, n_files)) {
  f <- files[ii]
  .x <- readRDS(f)
  .ari <- .x$ARI
  scn <- .x$ARI$Scenario |> unique()
  dataset_index <- .x$ARI$Index |> unique()
  
  # sum(colMeans(.x$MCMC$Semisupservised$allocations[[1]] == .x$MCMC$Semisupservised$allocations[[2]]) > 0.5)
  # sum(colMeans(.x$MCMC$Semisupservised$allocations[[1]] == .x$MCMC$Semisupservised$allocations[[3]]) > 0.5)
  # sum(colMeans(.x$MCMC$Semisupservised$allocations[[2]] == .x$MCMC$Semisupservised$allocations[[3]]) > 0.5)
  # colMeans(.x$MCMC$Semisupservised$phis)
  # 
  # sum(colMeans(.x$MCMC$Unsupservised$allocations[[1]] == .x$MCMC$Unsupservised$allocations[[2]]) > 0.5)
  # sum(colMeans(.x$MCMC$Unsupservised$allocations[[1]] == .x$MCMC$Unsupservised$allocations[[3]]) > 0.5)
  # sum(colMeans(.x$MCMC$Unsupservised$allocations[[2]] == .x$MCMC$Unsupservised$allocations[[3]]) > 0.5)
  # colMeans(.x$MCMC$Unsupservised$phis)
  
  .phi_semi <- .x$MCMC$Semisupservised$phis |> 
    as.data.frame() |> 
    magrittr::set_colnames(phi_names) |> 
    pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |> 
    dplyr::mutate(Method = "Semisupervised (Overfitted)", Scenario = scn, Dataset_index = dataset_index)
  
  .phi_un <- .x$MCMC$Unsupservised$phis |> 
    as.data.frame() |> 
    magrittr::set_colnames(phi_names) |> 
    pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |> 
    dplyr::mutate(Method = "Unsupervised", Scenario = scn, Dataset_index = dataset_index)
  
  .phi_df <- rbind(.phi_un, .phi_semi)
  
  
  if(f %in% k_6_files) {
    .phi_df <- .x$MCMC$Semisupservised$phis |> 
      as.data.frame() |> 
      magrittr::set_colnames(phi_names) |> 
      pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |> 
      dplyr::mutate(Method = "Semisupervised", Scenario = scn, Dataset_index = dataset_index)
    
    .ari$Unsupservised <- NA
    .ari$Semi.supservised.overfitted <- NA
    .ari$Mixture.model.overfitted <- NA
  } else {
    .ari$Semi.supservised.overfitted <- .ari$Semi.supservised
    .ari$Mixture.model.overfitted <- .ari$Mixture_model
    .ari$Semi.supservised <- NA
    .ari$Mixture_model <- NA
  }
  
  # for(v in seq(1, V)) {
  # .x$MCMC$Semisupservised$weights[[v]] |> 
  #   as.data.frame() |> 
  #   magrittr::set_colnames(phi_names) |> 
  #   pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |> 
  #   dplyr::mutate(Method = "Unsupervised", Scenario = scn, Dataset_index = dataset_index)
  # }
  
  if(ii == 1) {
    ari_df <- .ari
    phi_df <- .phi_df
  } else {
    ari_df <- rbind(ari_df, .ari)
    phi_df <- rbind(phi_df, .phi_df)
  }
}

ari_df$View <- factor(ari_df$View)

p_diff <- ari_df |> 
  pivot_longer(c(Difference_unsupervised, Difference_mixture_model), names_prefix = "Difference_", values_to = "Difference", names_to = "Model") |> 
  dplyr::mutate(Model = stringr::str_to_sentence(Model)) |> 
  ggplot(aes(x = Difference, y = View, fill = Model)) +
  geom_boxplot() +
  facet_wrap(~ Scenario, nrow = 1) +
  geom_vline(xintercept = 0.0, lty = 2, color = "red") +
  scale_fill_manual(values = c("#E69F00", "#009E73")) +
  labs(x = "Difference between point estimates")
  
p_model <- ari_df |>
  pivot_longer(c(Semi.supservised, Unsupservised, Mixture_model, Semi.supservised.overfitted, Mixture.model.overfitted), names_to = "Model", values_to = "ARI") |> 
  ggplot(aes(x = ARI, fill = Model, y = View)) +
  geom_boxplot() +
  facet_wrap(~ Scenario, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  labs(x = "ARI between inferred point clusterings and truth")
  # ggthemes::scale_fill_colorblind()

p_model / p_diff  +
  plot_layout(guides = 'collect')

intercept_df <- data.frame(Parameter = phi_names, Value = c(15, 10, 5))

phi_df |> 
  ggplot(aes(x = Sampled_values, y = Scenario, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~ Parameter) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  coord_cartesian(xlim=c(0, 40)) +
  geom_vline(aes(xintercept = Value),data=intercept_df, lty = 2, colour = "#21677e")

# +
  # labs(x = "ARI between inferred point clusterings and truth")
# ggthemes::scale_fill_colorblind()

# psm_mix <- .x$MCMC$Mixture_model$psms |> prepSimilarityMatricesForGGplot(use_common_ordering = FALSE)
# psm_mix |>
#   ggplot(aes(x = x, y = y, fill = Entry)) +
#   geom_tile() +
#   facet_wrap(~Chain) +
#   scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
# 
# 
# psm_mix <- .x$MCMC$Mixture_model$psms |> prepSimilarityMatricesForGGplot()
# psm_mix <- .x$MCMC$Mixture_model$psms |> prepSimilarityMatricesForGGplot()
