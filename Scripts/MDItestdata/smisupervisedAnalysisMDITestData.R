
suppressMessages(library(optparse))
library(tagmReDraft)
library(tidyr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(RcppParallel)

input_arguments <- function() {
  option_list <- list(
    
    # Number of MCMC iterations
    optparse::make_option(c("-R", "--R"),
                          type = "numeric",
                          default = 50000,
                          help = "Number of iterations to run in each MCMC chain [default= %default]",
                          metavar = "numeric"
    ),
    optparse::make_option(c("-t", "--thin"),
                          type = "numeric",
                          default = 500,
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
    optparse::make_option(c("--n_views"),
                          type = "numeric",
                          default = 5L,
                          help = "Number of datasets modelled [default= %default]",
                          metavar = "numeric"
    ),
    
    optparse::make_option(c("-K", "--K"),
                          type = "numeric",
                          default = 50,
                          help = paste(
                            "Number of components modelled in each dataset. If a dataset is",
                            "semi-supervised then the number of unique labels is modelled, if",
                            "unsupervised we default to 50."
                          ),
                          metavar = "numeric"
    ),
    
    optparse::make_option(c("--test_frac"),
                          type = "numeric",
                          default = 0.8,
                          help = "Fraction of labels unobserved in dataset 1.",
                          metavar = "numeric"
    ),
    
    optparse::make_option(c("--save_dir"),
                          type = "character",
                          default = "./",
                          help = "Directory to save output to [default= %default]",
                          metavar = "character"
    ),
    optparse::make_option(c("--data_dir"),
                          type = "character",
                          default = "~/rds/hpc-work/tagmmdi/Data/MDITest/",
                          help = "File to read input from from",
                          metavar = "character"
    )
  )
  
  
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

setThreadOptions(numThreads = "auto", stackSize = "auto")

cat("\n# === MDITESTDATA ANALYSIS ============================================")

cat("\nRead in data.")

# Pass the inputs from the command line
args <- input_arguments()

# Random seed used defining this fold
seed <- args$seed
set.seed(seed)

# Directories for input and output respectively
data_dir <- args$data_dir
save_dir <- args$save_dir

# Fraction of labels to consider unobserved in dataset 1
test_frac <- args$test_frac

sample_file <- paste0(save_dir, "mcmcMDItestdataTestFrac", test_frac, "Seed", seed, ".rds")
plot_file <- paste0(save_dir, "mcmcMDItestdataTestFrac", test_frac, "Seed", seed, ".png")

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
  n_clust_unsupervised <- 100
}


# Where deos the MDI test data live
# data_dir <- "../tagmmdi_wd/MDITestData/"
# data_dir <- "~/rds/hpc-work/tagmmdi/Data/MDITest/"

data_files <- list.files(data_dir, pattern = "MDItestdata", full.names = TRUE)
V <- n_files <- length(data_files)
view_inds <- seq(1, V)
test_data <- list()
for(ii in seq(1, n_files)) {
  f <- data_files[ii]
  test_data[[ii]] <- as.matrix(read.csv(f, row.names = 1))
}

test_1_true_labels <- c(
  rep(1, 15),
  rep(2, 35),
  rep(3, 6),
  rep(4, 8),
  rep(5, 11),
  rep(6, 8),
  rep(7, 17)
)

N <- nrow(test_data[[1]])
R <- args$R # 50000
thin <- args$thin # 500
n_chains <- args$n_chains # 1
K_max <- rep(K, V)
types <- rep("GP", V)
fixed <- initial_labels <- matrix(0, nrow = N, ncol = V)
for(v in view_inds) {
  initial_labels[, v] <- generateInitialUnsupervisedLabels(N, 1, K_max[v])
}

fixed[, 1] <- sample(c(0, 1), N, replace = TRUE, prob = c(test_frac, 1 - test_frac))
initial_labels[, 1] <- generateInitialSemiSupervisedLabels(test_1_true_labels, fixed[, 1])

views_used <- args$n_views # 5


cat("\n# === Modelling =======================================================")
cat("\nModel", views_used, "views and", n_chains, "using MDI.")

t0 <- Sys.time()
mcmc_out <- runMCMCChains(test_data[1:views_used], 
                          n_chains = n_chains, 
                          R = R, 
                          thin = thin, 
                          initial_labels = initial_labels[, 1:views_used], 
                          fixed = fixed[, 1:views_used], 
                          types = types[1:views_used], 
                          K = K_max[1:views_used],
                          initial_labels_as_intended = TRUE)

# mcmc_out |> str()
t1 <- Sys.time()
t_diff <- round(t1 - t0, 2)
unit <- attr(t_diff, "units")
cat("\nMCMC finished.\nTime taken:", t_diff, unit ,"\nSave samples.")

saveRDS(mcmc_out, file = sample_file) # "./mcmcMDItestdata.rds")

cat("\n# === Plotting ========================================================")

if(views_used == 5) {
  phi_names <- c("Phi_12", "Phi_13", "Phi_14", "Phi_15",
                 "Phi_23", "Phi_24", "Phi_25",
                 "Phi_34", "Phi_35",
                 "Phi_45"
  )
}
if(views_used == 3) {
  phi_names <- c("Phi_12", "Phi_13", "Phi_23")
}
if(views_used == 4) {
  phi_names <- c("Phi_12", "Phi_13", "Phi_14", "Phi_23", "Phi_24", "Phi_34")
}
if(views_used == 6) {
  phi_names <- c("Phi_12", "Phi_13", "Phi_14", "Phi_15", "Phi_16",
                 "Phi_23", "Phi_24", "Phi_25", "Phi_26",
                 "Phi_34", "Phi_35", "Phi_36",
                 "Phi_45", "Phi_46",
                 "Phi_56"
  )
}

p_phis <- mcmc_out[[1]]$phis |> 
  as.data.frame() |> 
  set_colnames(phi_names) |> 
  mutate(Iteration = seq(0, R, thin)) |> 
  pivot_longer(-Iteration, names_to = "Parameter", values_to = "Sampled_value") |> 
  ggplot(aes(x = Iteration, y = Sampled_value, colour = Parameter)) +
  geom_line() +
  facet_wrap(~Parameter, scales = "free")


ggsave(plot_file, plot = p_phis)

cat("\n# === SCRIPT COMPLETE =================================================")
