library(tidyr)
library(ggplot2)
mdiHelpR::setMyTheme()

main_dir <- "./Simulations/Output"
scn_dirs <- list.dirs(main_dir, full.names = FALSE, recursive = FALSE)
n_scn <- length(scn_dirs)
files <- list()
for(ii in seq(1, n_scn)) {
  scn <- scn_dirs[ii]
  .curr_dir <- paste0(main_dir, "/", scn)
  files[[scn]] <- .files <- list.files(.curr_dir, pattern = "*testFrac_80_nChains_10_K_50_kUnsupservised_50.png", full.names = TRUE)
  n_sims <- length(.files)
  for(jj in seq(1, n_sims)) {
    .f <- .files[jj]
    .x <- readRDS(.f)
    
    semi_phis <- .x$MCMC$Semisupservised$phis |> 
      as.data.frame() |> 
      magrittr::set_colnames(c("Phi_12", "Phi_13", "Phi_23")) |> 
      dplyr::mutate(Scenario = scn, Method = "Semisupervised")
    un_phis <- .x$MCMC$Unsupservised$phis |> 
      as.data.frame() |> 
      magrittr::set_colnames(c("Phi_12", "Phi_13", "Phi_23")) |> 
      dplyr::mutate(Scenario = scn, Method = "Unsupervised")
    
    phi_df <- rbind(semi_phis, un_phis)
    
    if(ii == 1 && jj == 1) {
      ari_df <- .x$ARI
      phis_df <- phi_df
    } else {
      ari_df <- rbind(ari_df, .x$ARI)
      phis_df <- rbind(phis_df, phi_df)
    }
  }
}

long_ari_df <- ari_df |> 
  pivot_longer(c(Semi.supservised, Unsupservised, Difference), names_to = "Method", values_to = "ARI") |> 
  dplyr::mutate(Scenario = factor(Scenario, levels = c("Gaussian", "MVT", "MixGaussians", "MixMVTs")))

p_perf <- long_ari_df |> 
  dplyr::filter(Method != "Difference") |> 
  ggplot(aes(x = Method, y = ARI)) + 
  geom_boxplot(aes(fill = Method)) +
  facet_grid(View ~ Scenario, labeller = label_both) +
  ggthemes::scale_fill_colorblind()

ggsave("./Simulations/PerfAcrossScenariosAndViews.png", plot = p_perf, width = 11, height = 6.5)

long_ari_df |> 
  dplyr::filter(Method == "Difference") |> 
  ggplot(aes(x = Scenario, y = ARI)) + 
  geom_boxplot(aes(fill = Scenario)) +
  facet_wrap(~View, labeller = label_both) +
  ggthemes::scale_fill_colorblind() +
  labs(title = "Difference in performance between semi-supervised and unsupervised MDI across views")

phis_df |> head()

phis_df |> 
  pivot_longer(c(Phi_12, Phi_13, Phi_23), names_to = "Parameter", values_to = "Sampled_value") |> 
  dplyr::filter(Sampled_value < 100) |> 
  ggplot(aes(x = Parameter, y = Sampled_value)) +
  geom_boxplot(aes(fill = Method)) +
  facet_wrap(~Scenario) +
  ggthemes::scale_fill_colorblind()

psm_df <- .x$MCMC$Unsupservised$psms |> 
  tagmReDraft::prepSimilarityMatricesForGGplot()

psm_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() + 
  facet_wrap(~Chain) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
