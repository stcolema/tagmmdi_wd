#!/usr/bin/Rscript
# Name: callTgondiiRNAseq.R
# Summary: Runs a mixutre model on the output RNA-seq dataset output by
# ``tGondiiMDIInputs.R``.
#
# Example: Rscript callTgondiiCellCycle.R --R 25000 --thin 50 --seed 1
#   --K 125 --save_dir "./" --data_file "./TGondiiMDI_seed_1_K_125.rds"
#
# Author: Stephen Coleman
suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(tagmReDraft))
suppressMessages(library(optparse))


# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Number of MCMC iterations
    optparse::make_option(c("-R", "--R"),
      type = "numeric",
      default = 10000,
      help = "Number of iterations to run in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-t", "--thin"),
      type = "numeric",
      default = 50,
      help = "Thinning factor in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-s", "--seed"),
      type = "numeric",
      default = 1,
      help = "Random seed that defines the test/train partition [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--n_chains"),
      type = "numeric",
      default = 1L,
      help = "Number of MCMC chains to run [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-K", "--K"),
      type = "numeric",
      default = NULL,
      help = paste(
        "Number of components modelled in each dataset. If a dataset is",
        "semi-supervised then the number of unique labels is modelled, if",
        "unsupervised we default to 100."
      ),
      metavar = "numeric"
    ),
    optparse::make_option(c("--type"),
      type = "character",
      default = "MVN",
      help = "Choice of density in the mixture model. [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(c("--data_file"),
      type = "character",
      help = "File to read input from from",
      metavar = "character"
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === T. gondii analysis =======================================================

cat("\n\n=== PASS INPUTS ===================================================\n")

t0 <- Sys.time()

# Pass the inputs from the command line
args <- input_arguments()

# Directories for input and output respectively
data_file <- args$data_file
save_dir <- args$save_dir

# MCMC sampler arguments
R <- args$R
thin <- args$thin

# Number of chains modelled
n_chains <- args$n_chains

# The number of components modelled
K <- args$K

# Number of clusters modelled in the categorical dataset
n_clust_unsupervised <- K
if (is.null(K)) {
  K <- 100
}

types <- args$type

# Random seed used defining this fold
seed <- args$seed

# Output saved to
save_file <- paste0(
  save_dir,
  "TGondii_RNAseq_",
  "seed_",
  seed,
  "_K_",
  n_clust_unsupervised,
  "_R_",
  R,
  ".rds"
)

# random seed
set.seed(seed)

mcmc_input <- readRDS(data_file)


data_modelled <- list(mcmc_input$data_modelled[[2]])
initial_labels <- mcmc_input$initial_labels[, 2, drop = FALSE]
fixed <- mcmc_input$fixed[, 2, drop = FALSE]

cat("\n\n=== MODELLING =====================================================\n")

mcmc_output <- runMCMCChains(data_modelled, n_chains,
  R = R,
  thin = thin,
  initial_labels = initial_labels,
  fixed = fixed,
  K = K,
  types = types
)

saveRDS(mcmc_output, file = save_file)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")

t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))
