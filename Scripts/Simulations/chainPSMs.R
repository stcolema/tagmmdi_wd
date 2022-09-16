library(tagmReDraft)

comparePSMsAcrossChains <- function(mcmc,
                                    matrix_setting_order = 1,
                                    use_common_ordering = TRUE,
                                    ignore_checks = TRUE) {
  n_chains <- length(mcmc)
  V <- mcmc[[1]]$V
  
  chain_inds <- seq(1, n_chains)
  view_inds <- seq(1, V)
  
  psm_df <- NULL
  first_iteration <- is.null(psm_df)
  psm_lst <- list()
  for(v in view_inds) {
    psm_lst[[v]] <- list()
    for(ii in chain_inds) {
      psm_lst[[v]][[ii]] <- mcmc[[ii]]$psm[[v]]
    }
  }
  for(v in view_inds) {
    .psm_df <- prepSimilarityMatricesForGGplot(psm_lst[[v]], 
      matrix_setting_order = matrix_setting_order,
      use_common_ordering = use_common_ordering,
      ignore_checks = ignore_checks
    )
    .psm_df$View <- v
    if(first_iteration) {
      psm_df <- .psm_df
      first_iteration <- FALSE
    } else {
      psm_df <- rbind(psm_df, .psm_df)
    }
  }
  psm_df$Chain <- factor(psm_df$Chain, levels = chain_inds)
  psm_df$View <- factor(psm_df$View, levels = view_inds)
  psm_df
}


data_dir <- "./Simulations/Output/"
scenarios <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
n_scn <- length(scenarios)
V <- 3
scn_dirs <- paste0(data_dir, scenarios)

for(ii in seq(1, n_scn)) {
  scn <- scn_dirs[ii]
  files <- list.files(scn, full.names = TRUE)
  n_files <- length(files)
  for(jj in seq(1, n_files)) {
    f <- files[jj]
    x <- readRDS(f)
    mcmc <- x$MCMC
    methods <- names(mcmc)
    n_methods <- length(methods)
    for(kk in seq(1, n_methods)) {
      m <- methods[kk]
      mixture_model <- m %in% methods[c(4, 5)]
      if(mixture_model) {
        for(v in seq(1, V)) {
          psm_df <- comparePSMsAcrossChains(mcmc[[m]][[v]], matrix_setting_order = 6)
        }
        
      } else {
        psm_df <- comparePSMsAcrossChains(mcmc[[m]], matrix_setting_order = 5)
      }
      
      psm_df |>
        ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = Entry)) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(View~Chain) +
        ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
    }
  }
}

chains_used <- data.frame(
  Model = c(
    "MDI_unsupervised",
    "MDI_semisupervised_overfitted",
    "MDI_semisupervised_k_known", 
    "Mixture_model_overfitted",
    "Mixture_model_k_known",
    "MDI_unsupervised",
    "MDI_semisupervised_overfitted",
    "MDI_semisupervised_k_known",
    "Mixture_model_overfitted",
    "Mixture_model_k_known",
    "MDI_unsupervised",
    "MDI_semisupervised_overfitted",
    "MDI_semisupervised_k_known",
    "Mixture_model_overfitted",
    "Mixture_model_k_known",
    "MDI_unsupervised",
    "MDI_semisupervised_overfitted",
    "MDI_semisupervised_k_known",
    "Mixture_model_overfitted",
    "Mixture_model_k_known",
    "MDI_unsupervised",
    "MDI_semisupervised_overfitted",
    "MDI_semisupervised_k_known",
    "Mixture_model_overfitted",
    "Mixture_model_k_known",
  ),
  Chain = c(6, 1, 1, 2, 2, 3, 4, 5, 2, 2, 3, 5, 6, 4, 10, 1, 7, 7, 4, 4, 6, 3, 5, 6, 6),
  Seed = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5),
  Scenario = c(
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian",
    "Gaussian"
  )
)
