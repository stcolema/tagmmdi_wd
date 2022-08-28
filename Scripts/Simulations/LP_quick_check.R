
library(tidyverse)
library(tagmReDraft)
library(magrittr)
library(mdiHelpR)

set.seed(1)
setMyTheme()

my_dir <- "./Simulations/Output/"
my_dirs <- list.dirs(my_dir, recursive = F)[c(1, 2 ,5)]

files <- lapply(my_dirs, list.files, full.names = TRUE) |> unlist()
# files <- list.files(my_dir, full.names = TRUE)
n_files <- length(files)
# outputs <- lapply(files, readRDS)


phi_names <- c("Phi_12", "Phi_13", "Phi_23")
phi_df <- NULL
for(ii in seq(1, n_files)) {
  f <- files[ii]
  output <- readRDS(f)
  models <- names(output$MCMC)
  n_models <- length(output$MCMC)
  n_chains <- length(output$MCMC$MDI_unsupervised)
  mdi_models <- models[seq(1, 3)]
  scn <- output$ARI$Scenario |> unique()
  if(ii == 1) {
    perf_df <- output$ARI
  } else {
    perf_df <- rbind(perf_df, output$ARI)
  }
  
  for(jj in seq(1, n_models)) {
    model <- names(output$MCMC)[jj]
    if(model %in% mdi_models) {
      for(kk in seq(1, n_chains)) {
        R <- output$MCMC[[jj]][[kk]]$R
        thin <- output$MCMC[[jj]][[kk]]$thin
        burn <- output$MCMC[[jj]][[kk]]$burn
        iter <- seq(burn + thin, R, thin)
        .phis <- output$MCMC[[jj]][[kk]]$phis |> 
          as.data.frame() |> 
          set_colnames(phi_names) |> 
          mutate(Iteration = iter) |> 
          pivot_longer(-Iteration, names_to = "Parameter", values_to = "Sampled_value") |> 
          mutate(Scenario = scn, Index = ii, Chain = kk, Model = model)
        if(is.null(phi_df)) {
          phi_df <- .phis
        } else {
          phi_df <- rbind(phi_df, .phis)
        }
      } 
    }
  }
}


perf_df$View <- factor(perf_df$View)
perf_df$Chain <- factor(perf_df$Chain)
perf_df$Index <- factor(perf_df$Index)

phi_df$Chain <- factor(phi_df$Chain)
phi_df$Index <- factor(phi_df$Index)

perf_df |> 
  ggplot(aes(x = Model, y = ARI.test.labels, fill = Model)) +
  geom_boxplot() +
  facet_grid(Scenario~View) +
  ggthemes::scale_fill_colorblind()

perf_df |> 
  ggplot(aes(x = Model, y = ARI.all.labels, fill = Model)) +
  geom_boxplot() +
  facet_grid(Scenario~View) +
  ggthemes::scale_fill_colorblind()

phi_df |> 
  ggplot(aes(y = Parameter, x = Sampled_value, fill = Model)) +
  geom_boxplot() +
  facet_grid(Scenario~Model) +
  ggthemes::scale_fill_colorblind()

