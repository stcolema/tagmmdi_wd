
library(ggplot2)
library(tagmReDraft)

mdiHelpR::setMyTheme()

data_dir <- "./Simulations/Data/"
output_dir <- "./Simulations/Output/"
scenarios <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
n_scn <- length(scenarios)
V <- 3

plotting_psms <- FALSE
selection_df <- NULL

for (ii in seq(1, n_scn)) {
  scn <- scenarios[ii]
  scn_input_dir <- paste0(data_dir, scn)
  scn_output_dir <- paste0(output_dir, scn)

  input_files <- list.files(scn_input_dir, pattern = "*.rds", full.names = TRUE) |>
    stringr::str_sort(numeric = TRUE)
  output_files <- list.files(scn_output_dir, pattern = "*.rds", full.names = TRUE) |>
    stringr::str_sort(numeric = TRUE)

  n_files <- length(input_files)
  for (jj in seq(1, n_files)) {
    f_in <- input_files[jj]
    f_out <- paste0(scn_output_dir, "/dataset_", jj, "_testFrac_80_nChains_10_K_5_kUnsupservised_50.rds")
    # f_out <- output_files[jj]
    output_found <- file.exists(f_out)

    if (output_found) {
      x <- readRDS(f_in)
      truth <- list(x$Clustering$View_2, x$Clustering$View_3)

      mod_lst <- readRDS(f_out)
      models <- names(mod_lst$MCMC)
      n_models <- length(models)
      for (kk in seq(1, n_models)) {
        mod <- mod_lst$MCMC[[kk]]
        n_chains <- length(mod)
        chain_setting_order <- c(1, 1, 1)
        best_ari <- c(0, 0, 0)
        chain_setting_order_ari <- list("View_2" = c(0, 0), "View_3" = c(0, 0))
        mixture_model <- kk %in% c(4, 5)
        for (ll in seq(1, n_chains)) {
          ari_comp <- c(0, 0)
          if (!mixture_model) {
            pred_cl <- mod[[ll]]$pred

            for (v in seq(2, V)) {
              ari <- mcclust::arandi(pred_cl[[v]], truth[[v - 1]])
              if (ari > best_ari[v]) {
                chain_setting_order[v] <- ll
                best_ari[v] <- ari

                for (mm in seq(1, 2)) {
                  ari_comp[mm] <- mcclust::arandi(pred_cl[[mm + 1]], truth[[mm]])
                }
                chain_setting_order_ari[[v - 1]] <- ari_comp
              }
            }
          } else {
            for (v in seq(2, V)) {
              pred_cl <- mod[[ll]][[v]]$pred[[1]]
              ari <- mcclust::arandi(pred_cl, truth[[v - 1]])
              if (ari > best_ari[v]) {
                chain_setting_order[v] <- ll
                best_ari[v] <- ari

                for (mm in seq(1, 2)) {
                  ari_comp[mm] <- mcclust::arandi(
                    mod[[ll]][[mm + 1]]$pred[[1]],
                    truth[[mm]]
                  )
                }
                chain_setting_order_ari[[v - 1]] <- ari_comp
              }
            }
          }
        }
        chain_setting_order_view_used <- lapply(chain_setting_order_ari, mean) |>
          unlist() |>
          which.max()
        chain_setting_order_used <- chain_setting_order[chain_setting_order_view_used + 1]

        if (plotting_psms) {
          if (mixture_model) {
            for (v in seq(1, V)) {
              .df <- comparePSMsAcrossChains(mod[[v]], matrix_setting_order = chain_setting_order_used)
              .df$View <- v
              if (v == 1) {
                psm_df <- .df
              } else {
                psm_df <- rbind(psm_df, .df)
              }
            }
          } else {
            psm_df <- comparePSMsAcrossChains(mod, matrix_setting_order = chain_setting_order_used)
          }
        }

        # p_psm <- psm_df |>
        #   ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = Entry)) +
        #   ggplot2::geom_tile() +
        #   ggplot2::facet_grid(View ~ Chain) +
        #   ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
        # save_name <- paste0("~/Desktop/TAGMMDI_SIMS/", scn, "/", models[kk], "seed", jj, ".png")
        # ggsave(save_name, p_psm, height = 8, width = 12)

        entry <- data.frame(Scenario = scn, Model = models[kk], Seed = jj, Chain_used = chain_setting_order_used)
        if (is.null(selection_df)) {
          selection_df <- entry
        } else {
          selection_df <- rbind(selection_df, entry)
        }
      }
    }
  }
}

write.csv(selection_df, "./Simulations/Output/OutputChainsUsed.csv")
