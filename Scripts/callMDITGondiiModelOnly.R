#!/usr/bin/Rscript 
# Summary: Runs MDI using the output of ``tGondiiMDIInputs.R``, a microarray 
# cell cycle dataset, some RNA-seq data and the LOPIT data. the first two are 
# unsupervised modelled using MVN densities locally. The LOPIT data is modelled
# using a TAGM model and is semi-supservised.
# 
# Example: Rscript callMDITGondiiModelOnly.R --R 25000 --thin 50 --seed 1 
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
                            "unsupervised we default to 50."
                          ),
                          metavar = "numeric"
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

# ggplot2 theme
setMyTheme()

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

# Random seed used defining this fold
seed <- args$seed

# Output saved to 
save_file <- paste0(
  save_dir,
  "TGondiiMDI_",
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

# Number of clusters modelled in the unsupervised views
n_clust_unsupervised <- K

data_modelled <- mcmc_input$data_modelled
initial_labels <- mcmc_input$initial_labels
fixed <- mcmc_input$fixed
K <- mcmc_input$K

# The number of components is a little awkward as it is set in the 
# semi-supservised view but should be changeable elsewhere, so this slightly 
# hack-y solution handles that.
unsupervised_changed <- (K[1] != n_clust_unsupervised) & (K[2] != n_clust_unsupervised)

if (unsupervised_changed) {
  K[1] <- K[2] <- n_clust_unsupervised
}

types <- mcmc_input$types

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
