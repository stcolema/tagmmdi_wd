#!/usr/bin/env Rscript
# 
# Title: MDI / TAGM Cross validation (single fold)
# Description: This script performs a single fold of cross-validation for data
# from ``pRolocdata``. The user inputs the datasets, the number of iterations to
# run, the thinning factor within each MCMC chain, the number of samples to 
# discard as part of the warm up, the random seed that defines this fold,
# the fraction of the data held as a test set, the number of components to model
# in the unsupervised auxiliary data used in the unsupervised view of MDI,
# the number of chains to run for each model, the minimum number of entries 
# needed in a column in the auxiliary dataset to have that column included in 
# the modelled data, and the directory to save the output to, within which a new 
# directory with the same name as the main dataset will be created.
# Output: An .RDS object that holds a list of the model outputs from the fold 
# and some measures of model performance and fit. 
# 
# Example call: 
# Rscript mdiTagmCVSingleLoop.R --datasets "tan2009r1 tan2009r1goCC" --R 10000
# --thin 50 --burn 1000 --seed 1 --test_size 0.75 --K 100 --n_chains 5 
# --save_dir "./" --categorical_column_threshold 5

suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

# === Functions ================================================================

#' @title multi-class F1
#' @description Calculate the accuracy, the per-class F1 score and the macro- 
#' and weighted-F1 scores.
#' @param pred Factor. The predicted classes.
#' @param truth Factor wit the same levels as ``pred``. The true classification.
#' @return A list containing ``accuracy``, the prediction accuracy, ``f1``, a 
#' vector of the f1 score in each class, ``macro_f1``, the macro f1 score across
#' classes, and ``weighted_f1``, the f1 score weighted by class membership.
multiClassF1 <- function(pred, truth) {

  mismatch_in_lengths <- length(pred) != length(truth)
  if(mismatch_in_lengths)
    stop("Prediction vector and truth must be of the same length.")
  
  # Confusion matrix for current fold
  conf <- caret::confusionMatrix(
    data = pred,
    reference = truth
  )$table
  
  N <- length(pred)
  n_levels <- nrow(conf)
  seq_levels <- seq(1, n_levels)
  
  f1 <- rep(0, n_levels)
  
  for(ii in seq_levels) {
    precision <- conf[ii, ii] / sum(conf[ii, ])
    recall <- conf[ii, ii] / sum(conf[ , ii])
    f1[ii] <- (2 * precision * recall) / (precision + recall)
    if(is.na(f1[ii])) {
      f1[ii] <- 0
    }
  }
  
  macro_f1 <- mean(f1)
  weighted_f1 <- sum(f1 * colSums(conf)) / N
  accuracy <- sum(pred == truth) / N
  list(
    "accuracy" = accuracy,
    "f1" = f1,
    "macro_f1" = macro_f1,
    "weighted_f1" = weighted_f1
  )
}

prepInputsForModelRun <- function(MS_object, test.idx,
                                  MS_cat_object = NULL,
                                  # test_size = 0.2,
                                  n_clust_cat = 50,
                                  categorical_column_threshold = 0
) {
  
  # cat("\nPreparing inputs.")
  
  marker.data <- pRoloc::markerMSnSet(MS_object)
  
  # cat("\nMarker data obtained.")
  
  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)
  K <- length(pRoloc::getMarkerClasses(MS_object))
  
  # Create a data frame of the classes present and their associated number
  classes_pres <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(Class = classes_pres, Key = 1:length(classes_pres))
  
  # Number of views modelled
  V <- 1
  
  if (!is.null(MS_cat_object)) {
    # cat("\nAccessing categorical data.\n")
    
    marker.data.cat <- pRoloc::markerMSnSet(MS_cat_object)
    K <- c(K, n_clust_cat)
    V <- 2
    
    # cat("\nCategorical data accessed.")
  }
  
  # cat("\nSplitting data into training and test sets.")
  
  ## 'unseen' test set
  .test <- MSnbase::MSnSet(
    MSnbase::exprs(marker.data)[test.idx, ],
    MSnbase::fData(marker.data[test.idx, ]),
    MSnbase::pData(marker.data)
  )
  
  ## 'seen' training set
  .train <- MSnbase::MSnSet(
    MSnbase::exprs(marker.data)[-test.idx, ],
    MSnbase::fData(marker.data[-test.idx, ]),
    MSnbase::pData(marker.data)
  )
  
  # save true marker labels
  test.markers <- MSnbase::fData(.test)$markers
  test.labels <- match(test.markers, class_key$Class)
  
  # training labels
  train.markers <- MSnbase::fData(.train)$markers
  train.labels <- match(train.markers, class_key$Class)
  
  # cat("\nCombine trainig and test sets.")
  
  # create new combined MSnset
  main_data <- BiocGenerics::combine(.train, .test)
  N <- nrow(main_data)
  
  # Set levels of markers cateogries
  levels(MSnbase::fData(main_data)$markers) <- c(
    levels(
      MSnbase::fData(main_data)$markers
    ),
    "unknown"
  )
  
  # hide marker labels
  MSnbase::fData(main_data)[rownames(.test), "markers"] <- "unknown"
  
  initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)
  
  # Fix training points, allow test points to move component
  fix_vec <- (MSnbase::fData(main_data)[, "markers"] != "unknown") * 1
  
  fixed[, 1] <- fix_vec
  initial_labels[, 1] <- c(train.labels, test.labels)
  
  types <- c("TAGM")
  
  data_modelled <- list(
    MSnbase::exprs(main_data)
  )
  
  if (!is.null(MS_cat_object)) {
    # cat("\nEnsure auxiliary data has the same ordering as the main dataset.")
    
    # create new combined MSnset
    cat_data <- BiocGenerics::combine(
      marker.data.cat[-test.idx, ],
      marker.data.cat[test.idx, ]
    )
    
    
    # Set levels of markers cateogries
    levels(MSnbase::fData(cat_data)$markers) <- c(
      levels(
        MSnbase::fData(cat_data)$markers
      ),
      "unknown"
    )
    
    # hide marker labels
    MSnbase::fData(cat_data)[rownames(.test), "markers"] <- "unknown"
   
    cat("\nRemoving uninformative columns from the GO dataset.")
    cat("\nThreshold:", categorical_column_threshold)
    cat("\nInitial number of columns:", ncol(cat_data))
    informative_terms <- colSums(exprs(cat_data)) > categorical_column_threshold
    if (any(!informative_terms)) {
      cat_data <- cat_data[, informative_terms]
    }
    cat("\nNumber of columns after reduction:", ncol(cat_data))
    
    types <- c("TAGM", "C")
    
    if(any(fData(cat_data)[, "markers"] != fData(main_data)[, "markers"])) {
      stop("\n\nERROR: Main and auxiliary data have different markers.\n")
    }
    
    if(any( row.names(cat_data) != row.names(main_data)) ) {
      stop("\n\nERROR: Main and auxiliary data have differing row names.\n")
    }
    
    # cat("\nList of datasets prepared for model call.")
    
    data_modelled <- list(
      MSnbase::exprs(main_data),
      MSnbase::exprs(cat_data)
    )
  }
  
  # cat("\nOutput list of various inputs to the MDI model.")
  
  list(
    
    # Inputs for MDI/mixture model
    data_modelled = data_modelled,
    types = types,
    fixed = fixed,
    initial_labels = initial_labels,
    K = K,
    V = V,
    
    # Objects used in assessing performance
    test.markers = test.markers,
    test.labels = test.labels,
    
    train.markers = train.markers,
    train.labels = train.labels,
    
    classes_pres = classes_pres
  )
}

#' @title CV single fold
#' @description Performs a single fold of cross validation for either MDI or a 
#' mixture model using MS objects available from ``pRolocdata``.
cvSingleFold <- function(MS_object, test.idx,
                          MS_cat_object = NULL,
                          # test_size = 0.2,
                          num_iter = 1000,
                          burn = 0,
                          thin = 25,
                          n_clust_cat = 50,
                          n_chains = 4,
                          categorical_column_threshold = 0,
                          ...) {
  
  # Flag indicating if we are using MDI or a mixture model
  integrative_analysis <- ! is.null(MS_cat_object)

  # Used in the model fit matrices
  eff_r <- floor(num_iter / thin) + 1
  eff_burn <- floor(burn / thin) + 1
  n_recorded <- eff_r - eff_burn
  recorded_iterations <- seq(burn + thin, num_iter, thin)
  recorded_iteration_indices <- seq(eff_burn + 1, eff_r, 1)

  cat("\nPreparing model inputs.")
  
  # Convert the input data into the appropriate format for the model call
  model_inputs <- prepInputsForModelRun(MS_object, test.idx,
    MS_cat_object = MS_cat_object,
    # test_size = test_size,
    n_clust_cat = n_clust_cat,
    categorical_column_threshold = categorical_column_threshold
  )
  
  cat("\n\nModel inputs prepared.")
  
  data_modelled <- model_inputs$data_modelled
  initial_labels <- model_inputs$initial_labels
  fixed <- model_inputs$fixed
  types <- model_inputs$types
  K <- model_inputs$K
  V <- model_inputs$V
  
  cat("\nFit a maximum of N/2 components in the unsupervised data.")
  N <- nrow(data_modelled[[1]])


  cat("\nK:", K)
  cat("\nN:", N)
  cat("\nN/2:", floor(N/2))
 
  if(V == 2) { 
    # Set a limit on the number of components modelled
    if(K[2] > floor(N / 2)) {
      cat("\nReducing K.")
      K[2] <- floor(N / 2)
      cat("\nK now set to", K[2], "in the GO data.")
    }
  }
  
  # Used in validation steps
  test.markers <- model_inputs$test.markers
  test.labels <- model_inputs$test.labels
  classes_pres <- model_inputs$classes_pres
  
  # Create a data frame of the classes present and their associated number
  class_key <- data.frame(Class = classes_pres, Key = 1:length(classes_pres))
  
  cat("\nModel call.")
  
  # MDI
  mcmc_chains <- runMCMCChains(
    data_modelled,
    n_chains,
    num_iter,
    thin,
    initial_labels,
    fixed,
    types,
    K
  )
  
  cat("\nModel has been run.\nMake predictions.")

  # Check convergence
  likelihood_mat <- phi_mat <- matrix(0, 
    nrow = n_recorded, 
    ncol = n_chains
  )

  for(jj in seq(1, n_chains)) {
    likelihood_mat[, jj] <- mcmc_chains[[jj]]$complete_likelihood[recorded_iteration_indices]
    
    if(integrative_analysis) {
      phi_mat[, jj] <- mcmc_chains[[jj]]$phis[recorded_iteration_indices]
    }
    
  }
  
  # Add the column names to the likelihood matrix
  colnames(likelihood_mat) <- paste0("Chain ", seq(1, n_chains))
  colnames(phi_mat) <- paste0("Chain ", seq(1, n_chains))

  # Use all of the chains for the point estimates
  combined_chains <- predictFromMultipleChains(mcmc_chains, burn)
  allocation_matrix <- combined_chains$allocation_probability[[1]]
    
  # predicted classes
  predicted_classes <- combined_chains$pred[[1]]
  predicted_organelles <- class_key$Class[predicted_classes]

  # True allocation for test data
  reference <- factor(test.markers, levels = classes_pres)
  
  test_indices <- which(fixed[, 1] == 0)
  
  comparison <- factor(predicted_organelles[test_indices],
    levels = classes_pres
  )
  
  cat("\nCompute model performance scores.")
  
  # Calculate the accuracy and the per-class, macro and weighted F1 scores.
  prediction_scores <- multiClassF1(comparison, reference)

  cat("\nCalculate the Brier score.")
  
  # Create allocation matrices for truth, filled initially with 0's
  allocmatrix <- matrix(0,
    nrow = N_test,
    ncol = length(classes_pres)
  )

  # The numbers associated with the classes (i.e. numerical representation of
  # the classes)
  class_numerics <- seq(1, length(unique(test.markers)))

  # create allocation matrix
  for (j in seq_along(test.idx)) {
    # The class the current individual belongs to
    alloc <- as.numeric(test.markers, class_numerics)[j]

    # Enter this in the allocation matrix
    allocmatrix[j, alloc] <- 1
  }

  # Compute quadratic loss
  quadloss <- sum((allocmatrix - allocation_matrix[test_indices, ])^2)
  
  # Convert hte likelihood matrix to a data frame and add the iterations as a 
  # variable
  likelihood_df <- likelihood_mat %>% 
    data.frame() %>% 
    mutate(Iteration = recorded_iteration_indices * thin)
  
  # Declare an object to return for the mixture model
  phi_df <- NULL
  if(integrative_analysis) {
    phi_df <- phi_mat %>% 
      data.frame() %>% 
      mutate(Iteration = recorded_iteration_indices * thin)
  }
  
  # Return the chains, the point estimates, the true organelles, the quadratic 
  # loss score and the prediction scores
  list(mcmc_output = mcmc_chains,
       point_esimates = combined_chains,
       truth = reference,
       quadloss = quadloss,
       prediction_scores = prediction_scores,
       likelihood_df = likelihood_df,
       class_key = class_key,
       phi_df = phi_df # This is non-empty only if MDI is run
  )
       
}

# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--datasets"),
      type = "character",
      help = "Which dataset should be analysed. E.g., 'tan2009r1 tan2009r1goCC'.",
      metavar = "character"
    ),

    # Number of MCMC iterations
    optparse::make_option(c("-R", "--R"),
      type = "numeric",
      default = 5000,
      help = "Number of iterations to run in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-t", "--thin"),
      type = "numeric",
      default = 25,
      help = "Thinning factor in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-b", "--burn"),
      type = "numeric",
      default = 1000,
      help = paste(
        "Number of iterations to burn of the warm-up period in each MCMC chain",
        "[default= %default]"
      ),
      metavar = "numeric"
    ),
    optparse::make_option(c("-s", "--seed"),
      type = "numeric",
      default = 1,
      help = "Random seed that defines the test/train partition [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--test_size"),
      type = "numeric",
      default = 0.2,
      help = "Fraction of observed labels used as a test set in each fold [default= %default]",
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
    optparse::make_option(c("--n_chains"),
      type = "numeric",
      default = 4,
      help = "Number of MCMC chains run in each fold [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(c("-c", "--categorical_column_threshold"),
      type = "numeric",
      default = 0L,
      help = "Required percentage representation in category [default= %default]",
      metavar = "numeric"
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === Cross-validation =========================================================

cat("\n\n=== PASS INPUTS ===================================================\n")

t0 <- Sys.time()

# ggplot2 theme
setMyTheme()

# Pass the inputs from the command line
args <- input_arguments()

# Random seed used defining this fold
seed <- args$seed

# MCMC sampler arguments
num_iter <- args$R
thin <- args$thin
burn <- args$burn

# Fraction of data used for validation in fold
test_size <- args$test_size

# Number of chains modelled
n_chains <- args$n_chains

# The number of components modelled
K <- args$K

# The minimum number of entries required for a column of the GO term data
# to be included
categorical_column_threshold <- args$categorical_column_threshold

# Number of clusters modelled in the categorical dataset
n_clust_cat <- K
if(is.null(K))
  n_clust_cat <- 75

# random seed
set.seed(seed)

# Convert the single string input into multiple possible strings
datasets <- unlist(stringr::str_split(args$datasets, " "))

test_size_str <- paste0(test_size * 100)

# Use the passed directory appened by the dataset. Create this directory if 
# required.
save_dir <- paste0(args$save_dir, datasets[1], "/")
dir.create(save_dir, showWarnings = FALSE)

# Details of the save files
save_name <- paste0(
  save_dir,
  paste(datasets, collapse = "_"),
  "_R_",
  num_iter,
  "_seed_",
  seed,
  "_nChains_",
  n_chains,
  "_testSize_",
  test_size_str
)

# What will the saved object be called
list_save_name <- paste0(save_name, ".rds")

cat("\nLoading data.")

# Load the relevant data
data(list = datasets[1])
data(list = datasets[2])

# Convert to the correct format to be passed to the validation function
d1 <- eval(parse(text = datasets[1]))
d2 <- eval(parse(text = datasets[2]))

cat("\nData loaded.")

# Number of samples modelled
N <- nrow(d1)

# Possibly reduce the dimensionality of the categorical dataset depending on 
# the number of entries in each column

# Use the same test indices across methods
marker.data <- pRoloc::markerMSnSet(d1)
X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)

# get sizes
.size <- ceiling(table(MSnbase::fData(marker.data)$markers) * test_size)

# strata needs size to be ordered as they appear in data
.size <- .size[unique(MSnbase::fData(marker.data)$markers)]

# get strata indices
test.idx <- sampling::strata(X, "markers",
  size = .size,
  method = "srswor"
)$ID_unit

N_test <- length(test.idx)

cat("\n\n=== BEGIN MAIN FUNCTION ===========================================\n")

cat("\nMDI validation fold.\n")

# Do a fold of the cross validation for MDI and the single-view mixture model
mdi_fold <- cvSingleFold(MS_object = d1, 
  test.idx = test.idx,
  MS_cat_object = d2,
  num_iter = num_iter,
  burn = burn,
  thin = thin,
  n_clust_cat = n_clust_cat,
  n_chains = n_chains,
  categorical_column_threshold = categorical_column_threshold
)

cat("\n\nMDI has run. Now run a TAGM model on the main dataset.\n")

tagm_fold <- cvSingleFold(MS_object = d1, 
  test.idx = test.idx,
  num_iter = num_iter,
  burn = burn,
  thin = thin,
  n_clust_cat = n_clust_cat,
  n_chains = n_chains
)

cat("\n\nTAGM has run.\nSaving outputs.")

# Record the seed used in the likelihood / model fit data frame
mdi_fold$likelihood_df$Fold <- seed
tagm_fold$likelihood_df$Fold <- seed

# Combine the two objects into a single list and save it
output <- list(
  mdi_fold = mdi_fold, 
  tagm_fold = tagm_fold, 
  test.idx = test.idx,
  seed = seed
)

saveRDS(output, file = list_save_name)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")

t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))
