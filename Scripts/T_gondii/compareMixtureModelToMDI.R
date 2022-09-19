
suppressPackageStartupMessages(library(tagmReDraft))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(mdiHelpR))
suppressMessages(library(optparse))

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
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

t0 <- Sys.time()

args <- input_arguments()
file_dir <- args$file_dir
save_dir <- args$save_dir
burn <- args$burn

cell_cycle_pattern <- "TGondii_CellCycle_seed_*"
lopit_file_pattern <- "TGondii_TAGPM_seed_*"
mdi_file_pattern <- "*_K_125_R_30000_V_2.rds"

cell_cycle_files <- list.files(file_dir, pattern = cell_cycle_pattern, full.names = TRUE)
lopit_files <- list.files(file_dir, pattern = lopit_file_pattern, full.names = TRUE)
mdi_files <- list.files(file_dir, pattern = mdi_file_pattern, full.names = TRUE)

files <- list(
  "Cell_cycle" = cell_cycle_files,
  "LOPIT" = lopit_files,
  "MDI" = mdi_files
)

n_files <- length(mdi_files)
n_models <- length(files)

models <- list()

cat("\nRead in files.\n")

for(ii in seq(1, n_models)) {
  cat("\n", names(files)[ii], "\n")
  mod_files <- files[[ii]]
  x <- lapply(mod_files, function(y) readRDS(y)[[1]])
  for(jj in seq(1, n_files)) {
    x[[jj]]$Chain <- jj
  }
  models[[ii]] <- predictFromMultipleChains(x, burn = burn)
}

cat("\nSave model predictions.")

names(models) <- names(files)
saveRDS(models, paste0(save_dir, "processedModelOutputs.rds"))

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")

t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))


