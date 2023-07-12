suppressMessages(library(mdir))
suppressMessages(library(mdiHelpR))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))
suppressMessages(library(optparse))

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--model_output_dir"),
      type = "character",
      default = "./Simulations/Output/",
      help = "Directory where the model outputs are saved.",
      metavar = "character"
    ),
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./Simulations/Output/",
      help = "Directory where the plots will be saved.",
      metavar = "character"
    )
  )

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


mdiHelpR::setMyTheme()

args <- input_arguments()

out_dir <- args$model_output_dir
save_dir <- args$save_dir
selection_df <- read.csv(paste0(out_dir, "OutputChainsUsed.csv"))

scenarios <- list.dirs(out_dir, full.names = F, recursive = F)
scn_dirs <- paste0(out_dir, scenarios)

files <- lapply(scn_dirs, list.files, ".rds$", F, T) |> unlist()

k_5_files <- lapply(scn_dirs, list.files, "*K_5_*", F, T) |> unlist()

n_files <- length(files)

V <- 3
phi_names <- c("Phi_12", "Phi_13", "Phi_23")

matching_model_names <- selection_df$Model |>
  stringr::str_replace_all("_", ".") |>
  stringr::str_replace_all("semisupervised", "semi-supervised") |>
  stringr::str_replace_all(".k.known", "")

for (ii in seq(1, n_files)) {
  f <- files[ii]
  .x <- readRDS(f)
  .ari <- .x$ARI

  scn <- .x$ARI$Scenario |> unique()
  dataset_index <- .x$ARI$Index |> unique()

  models <- .ari$Model |> unique()
  n_model <- length(models)
  for (jj in seq(1, n_model)) {
    mod <- models[jj]

    rel_selection_df_ind <- which(
      (matching_model_names == mod) &
        (selection_df$Seed == dataset_index) &
        (selection_df$Scenario == scn)
    )

    chain_used <- selection_df$Chain_used[rel_selection_df_ind]

    if (jj %in% c(1, 2, 3)) {
      chain_used <- unique(chain_used)
      rel_phis <- .x$MCMC[[jj]][[chain_used]]$phis |>
        as.data.frame()
      rel_phis$Model <- mod
      rel_phis$Index <- dataset_index
      rel_phis$Scenario <- scn
      rel_ari <- .ari |>
        dplyr::filter(Model == mod, Scenario == scn, Index == dataset_index, Chain == chain_used)
    } else {
      for (v in seq(1, V)) {
        if (v == 1) {
          rel_ari <- .ari |>
            dplyr::filter(Model == mod, Scenario == scn, Index == dataset_index, Chain == chain_used[v])
        } else {
          .entry <- .ari |>
            dplyr::filter(Model == mod, Scenario == scn, Index == dataset_index, Chain == chain_used[v])
          rel_ari <- rbind(rel_ari, .entry)
        }
      }
    }


    if ((ii == 1) & (jj == 1)) {
      ari_model_selection <- rel_ari
      phi_df <- rel_phis
    } else {
      ari_model_selection <- rbind(ari_model_selection, rel_ari)
      if (jj %in% c(1, 2, 3)) {
        phi_df <- rbind(phi_df, rel_phis)
      }
    }
  }
}

# sum(colMeans(.x$MCMC$Semisupservised$allocations[[1]] == .x$MCMC$Semisupservised$allocations[[2]]) > 0.5)
# sum(colMeans(.x$MCMC$Semisupservised$allocations[[1]] == .x$MCMC$Semisupservised$allocations[[3]]) > 0.5)
# sum(colMeans(.x$MCMC$Semisupservised$allocations[[2]] == .x$MCMC$Semisupservised$allocations[[3]]) > 0.5)
# colMeans(.x$MCMC$Semisupservised$phis)
#
# sum(colMeans(.x$MCMC$Unsupservised$allocations[[1]] == .x$MCMC$Unsupservised$allocations[[2]]) > 0.5)
# sum(colMeans(.x$MCMC$Unsupservised$allocations[[1]] == .x$MCMC$Unsupservised$allocations[[3]]) > 0.5)
# sum(colMeans(.x$MCMC$Unsupservised$allocations[[2]] == .x$MCMC$Unsupservised$allocations[[3]]) > 0.5)
# colMeans(.x$MCMC$Unsupservised$phis)
#
#   .phi_semi <- .x$MCMC$Semisupservised$phis |>
#     as.data.frame() |>
#     magrittr::set_colnames(phi_names) |>
#     pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |>
#     dplyr::mutate(Method = "Semisupervised (Overfitted)", Scenario = scn, Dataset_index = dataset_index)
#
#   .phi_un <- .x$MCMC$Unsupservised$phis |>
#     as.data.frame() |>
#     magrittr::set_colnames(phi_names) |>
#     pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |>
#     dplyr::mutate(Method = "Unsupervised", Scenario = scn, Dataset_index = dataset_index)
#
#   .phi_df <- rbind(.phi_un, .phi_semi)
#
#
#   if(f %in% k_6_files) {
#     .phi_df <- .x$MCMC$Semisupservised$phis |>
#       as.data.frame() |>
#       magrittr::set_colnames(phi_names) |>
#       pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |>
#       dplyr::mutate(Method = "Semisupervised", Scenario = scn, Dataset_index = dataset_index)
#
#     .ari$Unsupservised <- NA
#     .ari$Semi.supservised.overfitted <- NA
#     .ari$Mixture.model.overfitted <- NA
#   } else {
#     .ari$Semi.supservised.overfitted <- .ari$Semi.supservised
#     .ari$Mixture.model.overfitted <- .ari$Mixture_model
#     .ari$Semi.supservised <- NA
#     .ari$Mixture_model <- NA
#   }
#
#   # for(v in seq(1, V)) {
#   # .x$MCMC$Semisupservised$weights[[v]] |>
#   #   as.data.frame() |>
#   #   magrittr::set_colnames(phi_names) |>
#   #   pivot_longer(everything(), names_to = "Parameter", values_to = "Sampled_values") |>
#   #   dplyr::mutate(Method = "Unsupervised", Scenario = scn, Dataset_index = dataset_index)
#   # }
#
#   if(ii == 1) {
#     ari_df <- .ari
#     phi_df <- .phi_df
#   } else {
#     ari_df <- rbind(ari_df, .ari)
#     phi_df <- rbind(phi_df, .phi_df)
#   }
# }

# ari_df$View <- factor(ari_df$View)
ari_model_selection$View <- factor(ari_model_selection$View)
ari_model_selection$Scenario <- factor(ari_model_selection$Scenario, levels = c("Gaussian", "MVT", "LogPoisson"), labels = c("Gaussian", "MVT", "Log-Poisson"))
ari_model_selection$Model <- factor(ari_model_selection$Model,
  levels = c(
    "Mixture.model",
    "Mixture.model.overfitted",
    "MDI.unsupervised",
    "MDI.semi-supervised.overfitted",
    "MDI.semi-supervised"
  ),
  labels = c(
    "Mixture model",
    "Mixture model (overfitted)",
    "MDI (unsupervised)",
    "MDI (semi-supervised; overfitted)",
    "MDI (semi-supervised)"
  )
)


phi_df$Scenario <- factor(phi_df$Scenario, levels = c("Gaussian", "MVT", "LogPoisson"), labels = c("Gaussian", "MVT", "Log-Poisson"))
phi_df$Model <- factor(phi_df$Model,
  levels = c(
    "MDI.unsupervised",
    "MDI.semi-supervised.overfitted",
    "MDI.semi-supervised"
  ),
  labels = c(
    "MDI (unsupervised)",
    "MDI (semi-supervised; overfitted)",
    "MDI (semi-supervised)"
  )
)
colnames(phi_df)[c(1, 2, 3)] <- c("phi_12", "phi_13", "phi_23")

phi_df |>
  pivot_longer(c(phi_12, phi_13, phi_23), names_to = "Parameter", values_to = "Sampled_value") |>
  ggplot(aes(x = Parameter, y = Sampled_value, fill = Model)) +
  geom_boxplot() +
  ggthemes::scale_fill_colorblind()

ari_model_selection$Dataset <- ari_model_selection$View
p_ari <- ari_model_selection |> ggplot(aes(y = Model, x = ARI.test.labels, fill = Model)) +
  geom_boxplot() +
  facet_grid(Scenario ~ Dataset, labeller = label_both) +
  ggthemes::scale_fill_colorblind() +
  labs(x = "ARI between inferred and true labels", y = NULL) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
# + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(paste0(save_dir, "PerfAcrossScenariosAndViews.png"),
  plot = p_ari,
  height = 7,
  width = 10
)

ggsave(paste0(save_dir, "PerfAcrossScenariosAndViews.pdf"),
  plot = p_ari,
  device = "pdf",
  height = 7,
  width = 10
)


# p_diff <- ari_df |>
#   pivot_longer(c(Difference_unsupervised, Difference_mixture_model), names_prefix = "Difference_", values_to = "Difference", names_to = "Model") |>
#   dplyr::mutate(Model = stringr::str_to_sentence(Model)) |>
#   ggplot(aes(x = Difference, y = View, fill = Model)) +
#   geom_boxplot() +
#   facet_wrap(~ Scenario, nrow = 1) +
#   geom_vline(xintercept = 0.0, lty = 2, color = "red") +
#   scale_fill_manual(values = c("#E69F00", "#009E73")) +
#   labs(x = "Difference between point estimates")
#
# p_model <- ari_df |>
#   pivot_longer(c(Semi.supservised, Unsupservised, Mixture_model, Semi.supservised.overfitted, Mixture.model.overfitted), names_to = "Model", values_to = "ARI") |>
#   ggplot(aes(x = ARI, fill = Model, y = View)) +
#   geom_boxplot() +
#   facet_wrap(~ Scenario, nrow = 1) +
#   scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
#   labs(x = "ARI between inferred point clusterings and truth")
#   # ggthemes::scale_fill_colorblind()
#
# p_model / p_diff  +
#   plot_layout(guides = 'collect')
#
# intercept_df <- data.frame(Parameter = phi_names, Value = c(15, 10, 5))
#
# phi_df |>
#   ggplot(aes(x = Sampled_values, y = Scenario, fill = Method)) +
#   geom_boxplot() +
#   facet_wrap(~ Parameter) +
#   scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
#   coord_cartesian(xlim=c(0, 40)) +
#   geom_vline(aes(xintercept = Value),data=intercept_df, lty = 2, colour = "#21677e")
#
# # +
#   # labs(x = "ARI between inferred point clusterings and truth")
# # ggthemes::scale_fill_colorblind()
#
# # psm_mix <- .x$MCMC$Mixture_model$psms |> prepSimilarityMatricesForGGplot(use_common_ordering = FALSE)
# # psm_mix |>
# #   ggplot(aes(x = x, y = y, fill = Entry)) +
# #   geom_tile() +
# #   facet_wrap(~Chain) +
# #   scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
# #
# #
# # psm_mix <- .x$MCMC$Mixture_model$psms |> prepSimilarityMatricesForGGplot()
# # psm_mix <- .x$MCMC$Mixture_model$psms |> prepSimilarityMatricesForGGplot()
