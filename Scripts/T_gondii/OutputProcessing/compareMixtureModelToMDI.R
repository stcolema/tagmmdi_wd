
suppressPackageStartupMessages(library(tagmReDraft))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(mdiHelpR))
suppressMessages(library(optparse))


ariAcrossChains <- function(mcmc_lst) {
  n_chains <- length(mcmc_lst)

  V <- mcmc_lst[[1]]$V
  R <- mcmc_lst[[1]]$R
  burn <- mcmc_lst[[1]]$burn
  thin <- mcmc_lst[[1]]$thin
  no_burn_in <- burn == 0
  
  first_iter <- burn + thin * (1 - no_burn_in)
  iteration <- seq(first_iter, R, thin)
  ari_lst <- list()
  for(ii in seq(1, n_chains)) {
    ari_lst[[ii]] <- ariAcrossViews(mcmc_lst[[ii]])
  }
  for(ii in seq(1, n_chains)) {
    .df <- as.data.frame(ari_lst[[ii]])
    colnames(.df) <- paste0("View", seq(1, V))
    .df$Chain <- mcmc_lst[[ii]]$Chain
    .df$Iteration <- iteration
    if(ii == 1) {
      ari_df <- .df
    } else {
      ari_df <- rbind(ari_df, .df)
    }
  }
  ari_df$Chain <- factor(ari_df$Chain)
  ari_df
}

ariAcrossViews <- function(mcmc) {
  V <- mcmc$V
  is_list <- is.list(mcmc$allocations)
  
  R <- mcmc$R
  burn <- mcmc$burn
  no_burn_in <- burn == 0
  
  n_samples <- (mcmc$R - mcmc$burn) / mcmc$thin + no_burn_in
  ari_mat <- matrix(0, nrow = n_samples, ncol = V)
  for(v in seq(1, V)) {
    if(is_list) {
      alloc <- mcmc$allocations[[v]]
    } else {
      alloc <- mcmc$allocations[, , v]
    }
    
    ari_mat[, v] <- apply(alloc, 1, function(x) {
      mcclust::arandi(x, alloc[n_samples, ])
    })
    
  }
  ari_mat
}

setMyTheme()

# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--file_dir"),
      type = "character",
      default = "/home/sdc56/rds/hpc-work/tagmmdi/T_gondii/Output/",
      help = "Directory containing model output [default= %default]",
      metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "/home/sdc56/rds/hpc-work/tagmmdi/T_gondii/Output/",
      help = "Directory to save processed model output to [default= %default]",
      metavar = "character"
    ),

    optparse::make_option(c("--burn"),
      type = "numeric",
      default = 10000L,
      help = "Number of MCMC iterations to discard [default= %default]",
      metavar = "numeric"
    ),

    optparse::make_option(c("--plot_psm"),
      type = "logical",
      default = TRUE,
      help = "Flag to plot PSMs [default= %default]",
      metavar = "logical"
    ),

   optparse::make_option(c("--make_predictions"),
     type = "logical",
      default = TRUE,
      help = "Flag to combine chains and make predictions [default= %default]",
      metavar = "logical"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

cat("\n# === Processing T. gondii model output =============================\n")

t0 <- Sys.time()

args <- input_arguments()
file_dir <- args$file_dir
save_dir <- args$save_dir
burn <- args$burn

plot_psm <- args$plot_psm
make_predictions <- args$make_predictions

cell_cycle_pattern <- "TGondii_CellCycle_seed_*"
lopit_file_pattern <- "TGondii_TAGPM_seed_*"
mdi_file_pattern <- "*_K_125_R_30000_V_2.rds"
mdi_v3_file_pattern <- "*_K_125_R_30000_V_3.rds"

cell_cycle_files <- list.files(file_dir, pattern = cell_cycle_pattern, full.names = TRUE)
lopit_files <- list.files(file_dir, pattern = lopit_file_pattern, full.names = TRUE)
mdi_files <- list.files(file_dir, pattern = mdi_file_pattern, full.names = TRUE)
mdi_v3_files <- list.files(file_dir, pattern = mdi_v3_file_pattern, full.names = TRUE)

files <- list(
  "Cell_cycle" = cell_cycle_files,
  "LOPIT" = lopit_files,
  "MDI" = mdi_files,
  "MDI_v3" = mdi_v3_files
)

model_names <- names(files)
plot_filenames <- paste0(save_dir, model_names, "PSMcomparison.png")
plot_width <- 12

main_output_filename <- paste0(save_dir, "processedModelOutputs.rds")

n_files <- length(mdi_files)
n_models <- length(files)

models <- list()
ari_lst <- list()
cat("\nRead in files.\n")

for(ii in seq(1, n_models)) {
  mod <- model_names[ii]
  cat("\n\n", mod, "\n")
  mod_files <- files[[ii]]
  x <- lapply(mod_files, function(y) readRDS(y)[[1]])
  for(jj in seq(1, n_files)) {
    x[[jj]]$Chain <- jj
  }

  # Data.frame comparing ARI between ith sample and last 
  ari_df <- ariAcrossChains(x)

  long_ari_df <- tidyr::pivot_longer(ari_df, -c(Chain, Iteration), 
    names_to = "View", 
    values_to = "ARI"
  )

  p_ari <- ggplot(long_ari_df, aes(x = Iteration, y = ARI, group = Chain)) +
    geom_line() +
    facet_wrap(~View) +
    labs(title = "ARI between current and final sampled partitions")

  ggsave(paste0(save_dir, mod, "ARI.png"), p_ari)

  # Apply a burn in and construct PSMs
  x <- tagmReDraft::processMCMCChains(x, burn, construct_psm = plot_psm)

  if(make_predictions) {
  # Combine samples across chains
  models[[ii]] <- predictFromMultipleChains(x, burn = burn, 
    chains_already_processed = TRUE
  )
  }

  if(plot_psm) {
  # Make the PSMs ready for ggplot2 and plot
  psm_df <- tagmReDraft::comparePSMsAcrossChains(x)

  # Set the view names
  if(mod == "MDI") {
    plot_height <- 10
    psm_df$View <- factor(psm_df$View, labels = c("LOPIT", "Cell-cycle microarray"))
  } else {
    plot_height <- 6
    if(mod == "Cell_cycle") {
      psm_df$View <- factor(psm_df$View, labels = c("Cell-cycle microarray"))
    }
    if(mod == "LOPIT") {
      psm_df$View <- factor(psm_df$View, labels = c("LOPIT"))
    }
  }

     p_psm <- psm_df |> 
    ggplot(aes(x = x, y = y, fill = Entry)) +
    geom_tile() +
    facet_grid(View ~ Chain, labeller = label_both) +
    scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
    labs(x = "Item", y = "Item", fill = "Coclustering\nproportion") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.title.y = element_text(size = 10.5),
      axis.title.x = element_text(size = 10.5),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14),
      strip.text.x = element_text(size = 10.5),
      legend.text = element_text(size = 10.5)
    )
  
  cat("\nSave PSM comparison across chains:\n", plot_filenames[ii])
  ggsave(plot_filenames[ii], p_psm, height = plot_height, width = plot_width)
  }
}

if(make_predictions) {
  cat("\nSave model predictions:\n", main_output_filename)

  names(models) <- names(files)
  saveRDS(models, paste0(save_dir, "processedModelOutputs.rds"))
}

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")

t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))


