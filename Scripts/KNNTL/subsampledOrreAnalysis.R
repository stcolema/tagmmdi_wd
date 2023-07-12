#!/usr/bin/env Rscript
#
# Title: KNN Transfer Learner Cross Validation
# Description: Performs a single loop of cross-validation of data from
# ``pRolocdata`` for the knn-transfer learner of Breckels et al. Takes an input
# from ``mdiTagmCVSingleLoop.R``, using this to ensure the same test/training
# data is used in this model as MDI and the mixture model.
# Output:
#
# Example call:
# Rscript knnTLCVSingleLoop.R --datasets 'tan2009r1 tan2009r1goCC' --seed 1
# --test_indices ./test_50/tan2009r1/tan2009r1_seed_1_testSize_50_trainingAdjustedForTL_TRUE.rds
# --save_dir ./test_50/ --categorical_column_threshold 3 --number_weights 5
# --number_weights_sampled 20000


suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(mdir))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))
suppressPackageStartupMessages(library(MSnbase))

subsetToKnownMarkers <- function(MSnSet_obj) {
  .X <- pRoloc:::subsetAsDataFrame(MSnSet_obj, "markers", train = TRUE)
  drop_na(.X)
}

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
  if (mismatch_in_lengths) {
    stop("Prediction vector and truth must be of the same length.")
  }

  # Confusion matrix for current fold
  conf <- caret::confusionMatrix(
    data = pred,
    reference = truth
  )$table

  N <- length(pred)
  n_levels <- nrow(conf)
  seq_levels <- seq(1, n_levels)

  f1 <- rep(0, n_levels)

  for (ii in seq_levels) {
    precision <- conf[ii, ii] / sum(conf[ii, ])
    recall <- conf[ii, ii] / sum(conf[, ii])
    f1[ii] <- (2 * precision * recall) / (precision + recall)
    if (is.na(f1[ii])) {
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


# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Datasets to analyse

    # Random seed (i.e. fold number)
    optparse::make_option(c("-s", "--seed"),
      type = "numeric",
      default = 1,
      help = "Random seed that defines the test/train partition [default= %default]",
      metavar = "numeric"
    ),

    # Save the output to this directory
    optparse::make_option(c("-d", "--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
    ),

    # Save the output to this directory
    optparse::make_option(c("--test_size"),
      type = "numeric",
      default = 0.7,
      help = "Fraction of data used as test set [default= %default]",
      metavar = "numeric"
    ),

    # Save the output to this directory
    optparse::make_option(c("-r", "--R"),
      type = "integer",
      default = 10000,
      help = "Number of iterations to run in each chain [default= %default]",
      metavar = "integer"
    ),
    optparse::make_option(c("-t", "--thin"),
      type = "integer",
      default = 50,
      help = "MCMC thinning factor [default= %default]",
      metavar = "integer"
    ),
    optparse::make_option(c("-b", "--burn"),
      type = "integer",
      default = 2500,
      help = "MCMC burn in [default= %default]",
      metavar = "integer"
    ),
    optparse::make_option(c("-n", "--n_chains"),
      type = "integer",
      default = 4,
      help = "MCMC burn in [default= %default]",
      metavar = "integer"
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

calcPerformanceMetrics <- function(combined_chains,
                                   class_key,
                                   classes_pres,
                                   test_markers,
                                   test_indices) {
  prediction_scores <- list()
  L <- length(combined_chains$allocation_probability)
  N_test <- length(test_indices)
  for (ll in seq(1, L)) {
    allocation_matrix <- combined_chains$allocation_probability[[ll]]

    # predicted classes
    predicted_classes <- combined_chains$pred[[ll]]
    predicted_organelles <- factor(predicted_classes, labels = levels(test_markers)) # class_key$Class[predicted_classes]

    # True allocation for test data
    reference <- test_markers

    comparison <- predicted_organelles[test_indices] # factor(predicted_organelles[test_indices],
    #   levels = classes_pres
    # )

    cat("\nCompute model performance scores.")

    # Calculate the accuracy and the per-class, macro and weighted F1 scores.
    prediction_scores[[ll]] <- multiClassF1(comparison, reference)

    cat("\nCalculate the Brier score.")

    # Create allocation matrices for truth, filled initially with 0's
    allocmatrix <- matrix(0,
      nrow = N_test,
      ncol = length(classes_pres)
    )

    # The numbers associated with the classes (i.e. numerical representation of
    # the classes)
    class_numerics <- seq(1, length(unique(test_markers)))

    # create allocation matrix
    for (j in seq_along(test_indices)) {
      # The class the current individual belongs to
      alloc <- as.numeric(test_markers, class_numerics)[j]

      # Enter this in the allocation matrix
      allocmatrix[j, alloc] <- 1
    }

    # Compute quadratic loss
    prediction_scores[[ll]]$quadloss <- sum((allocmatrix - allocation_matrix[test_indices, ])^2) / N_test
  }

  list(
    truth = reference,
    prediction_scores = prediction_scores
  )
}

fitModel <- function(mcmc_chains,
                     R,
                     thin,
                     burn,
                     test_markers,
                     test_indices,
                     class_key,
                     classes_pres,
                     MDI_run = TRUE) {
  # Used in the model fit matrices
  eff_r <- floor(R / thin) + 1
  eff_burn <- floor(burn / thin) + 1
  n_recorded <- eff_r - eff_burn
  recorded_iterations <- seq(burn + thin, R, thin)
  recorded_iteration_indices <- seq(eff_burn + 1, eff_r, 1)

  n_chains <- length(mcmc_chains)

  # Check convergence
  likelihood_mat <- phi_mat <- matrix(0,
    nrow = n_recorded,
    ncol = n_chains
  )

  for (jj in seq(1, n_chains)) {
    likelihood_mat[, jj] <- mcmc_chains[[jj]]$complete_likelihood[recorded_iteration_indices]

    if (MDI_run) {
      phi_mat[, jj] <- mcmc_chains[[jj]]$phis[recorded_iteration_indices]
    }
  }

  # Add the column names to the likelihood matrix
  colnames(likelihood_mat) <- paste0("Chain ", seq(1, n_chains))
  colnames(phi_mat) <- paste0("Chain ", seq(1, n_chains))

  # Use all of the chains for the point estimates
  combined_chains <- predictFromMultipleChains(mcmc_chains, burn)
  perf_lst <- calcPerformanceMetrics(
    combined_chains,
    class_key,
    classes_pres,
    test_markers,
    test_indices
  )

  # Convert the likelihood matrix to a data frame and add the iterations as a
  # variable
  likelihood_df <- likelihood_mat %>%
    data.frame() %>%
    mutate(Iteration = recorded_iteration_indices * thin)

  # Declare an object to return for the mixture model
  phi_df <- NULL
  if (MDI_run) {
    phi_df <- phi_mat %>%
      data.frame() %>%
      mutate(Iteration = recorded_iteration_indices * thin)
  }

  # Return the chains, the point estimates, the true organelles, the quadratic
  # loss score and the prediction scores
  list(
    mcmc_output = mcmc_chains,
    point_esimates = combined_chains,
    truth = perf_lst$reference,
    prediction_scores = perf_lst$prediction_scores,
    likelihood_df = likelihood_df,
    class_key = class_key,
    phi_df = phi_df # This is non-empty only if MDI is run
  )
}

# === Cross-validation =========================================================

cat("\n\n# === ORRE DATA =======================================================\n")

t0 <- Sys.time()

# Pass the inputs from the command line
args <- input_arguments()

# Random seed used defining this fold
seed <- args$seed

# MCMC sampler arguments
R <- num_iter <- args$R
thin <- args$thin
burn <- args$burn

# Fraction of data used for validation in fold
test_size <- args$test_size
test_size_str <- paste0(test_size * 100)

# Number of chains modelled
n_chains <- args$n_chains

save_dir <- args$save_dir

# Use the passed directory appened by the dataset. Create this directory if
# required.
save_dir <- paste0(args$save_dir, "Orre/")
dir.create(save_dir, showWarnings = FALSE)

# Details of the save files
list_save_name <- paste0(
  save_dir,
  "Orre_R_",
  num_iter,
  "_seed_",
  seed,
  "_nChains_",
  n_chains,
  "_testSize_",
  test_size_str,
  ".rds"
)

cat("\n\n# === LOAD DATA =======================================================")

set.seed(seed)


data("orre2019a431")
data("orre2019h322")
data("orre2019hcc827")
# data("orre2019hcc827gef")
# data("orre2019hcc827rep1")
data("orre2019mcf7")
data("orre2019u251")

# Reduce to the marker proteins
orre2019a431_markers <- subsetToKnownMarkers(orre2019a431) # "markers", train = TRUE)
orre2019h322_markers <- subsetToKnownMarkers(orre2019h322)
orre2019hcc827_markers <- subsetToKnownMarkers(orre2019hcc827)
# orre2019hcc827gef_markers <- subsetToKnownMarkers(orre2019hcc827gef)
# orre2019hcc827rep1_markers <- subsetToKnownMarkers(orre2019hcc827rep1)
orre2019mcf7_markers <- subsetToKnownMarkers(orre2019mcf7)
orre2019u251_markers <- subsetToKnownMarkers(orre2019u251)

# The names of these proteins
orre2019a431_proteins <- row.names(orre2019a431_markers)
orre2019h322_proteins <- row.names(orre2019h322_markers)
orre2019hcc827_proteins <- row.names(orre2019hcc827_markers)
# orre2019hcc827gef_proteins <- row.names(orre2019hcc827gef_markers)
# orre2019hcc827rep1_proteins <- row.names(orre2019hcc827rep1_markers)
orre2019mcf7_proteins <- row.names(orre2019mcf7_markers)
orre2019u251_proteins <- row.names(orre2019u251_markers)

# Model only the common proteins
proteins_captured <- Reduce(
  intersect,
  list(
    orre2019a431_proteins,
    orre2019h322_proteins,
    orre2019hcc827_proteins,
    # orre2019hcc827gef_proteins,
    # orre2019hcc827rep1_proteins,
    orre2019mcf7_proteins,
    orre2019u251_proteins
  )
)

n_measurements <- ncol(orre2019a431_markers) - 1
measurements <- seq(1, n_measurements)
marker_col <- n_measurements + 1
N <- length(proteins_captured)

# === Generate test indices ====================================================

# Adjust the trest indices to ensure all organelles are represented in the
# training data as required by the transfer learner
all_classes_represented_in_training <- TRUE

# Use the same test indices across methods
marker.data <- orre2019a431_markers[proteins_captured, marker_col]
X <- orre2019a431_markers[proteins_captured, ]

# get sizes
# .size <- ceiling(table(MSnbase::fData(marker.data)$markers) * test_size)
organelle_representation <- table(X$markers)

sampled_dataset_size <- 1100 / N
data_used <- ceiling(organelle_representation * sampled_dataset_size)
data_used <- data_used[unique(marker.data)]

# get strata indices
sampled_inds <- sampling::strata(X, "markers",
  size = data_used,
  method = "srswor"
)$ID_unit

X <- X[sampled_inds, ]
updated_organelle_representation <- table(X$markers)

.size <- ceiling(updated_organelle_representation * test_size)

if (all_classes_represented_in_training) {
  someOrganellesAdjusted <- FALSE
  organelles_adjusted <- c()
  for (ii in seq(1, length(updated_organelle_representation))) {
    no_training_members <- .size[ii] == updated_organelle_representation[ii]
    if (no_training_members) {
      someOrganellesAdjusted <- TRUE
      org <- names(updated_organelle_representation)[ii]
      organelles_adjusted <- c(organelles_adjusted, org)
      .size[ii] <- .size[ii] - 1
    }
  }
  if (someOrganellesAdjusted) {
    increased_training_size_str <- paste0(organelles_adjusted, collapse = ", ")
    cat(
      "\nNo training data for",
      increased_training_size_str,
      "at original test size (",
      test_size,
      ")."
    )
  }
}

# strata needs size to be ordered as they appear in data
.size <- .size[unique(marker.data)]

# get strata indices
test_inds <- sampling::strata(X, "markers",
  size = .size,
  method = "srswor"
)$ID_unit

# Sample test indices
N_test <- length(test_inds)

# === ===

proteins_captured <- proteins_captured[sampled_inds]

# The expression data
orre2019a431_exprs <- orre2019a431_markers[proteins_captured, measurements]
orre2019h322_exprs <- orre2019h322_markers[proteins_captured, measurements]
orre2019hcc827_exprs <- orre2019hcc827_markers[proteins_captured, measurements]
# orre2019hcc827gef_exprs <- orre2019hcc827gef_markers[proteins_captured, measurements]
# orre2019hcc827rep1_exprs <- orre2019hcc827rep1_markers[proteins_captured, measurements]
orre2019mcf7_exprs <- orre2019mcf7_markers[proteins_captured, measurements]
orre2019u251_exprs <- orre2019u251_markers[proteins_captured, measurements]

# Dataset dimensions
N <- length(proteins_captured)
L <- 5

# The ground truth to compare to
true_allocations <- matrix(c(
  orre2019a431_markers[proteins_captured, marker_col],
  orre2019h322_markers[proteins_captured, marker_col],
  orre2019hcc827_markers[proteins_captured, marker_col],
  # orre2019hcc827gef_markers[proteins_captured, marker_col],
  # orre2019hcc827rep1_markers[proteins_captured, marker_col],
  orre2019mcf7_markers[proteins_captured, marker_col],
  orre2019u251_markers[proteins_captured, marker_col]
), N, L)

train_inds <- seq(1, N)[-test_inds]

# === Model inputs =============================================================

# Create the matrix describing which proteins have known allocations
fixed <- matrix(0, N, L)
fixed[train_inds, ] <- 1

# The number of components to model
K <- c(
  length(pRoloc::getMarkerClasses(orre2019a431)),
  length(pRoloc::getMarkerClasses(orre2019h322)),
  length(pRoloc::getMarkerClasses(orre2019hcc827)),
  # length(pRoloc::getMarkerClasses(orre2019hcc827gef)),
  # length(pRoloc::getMarkerClasses(orre2019hcc827rep1)),
  length(pRoloc::getMarkerClasses(orre2019mcf7)),
  length(pRoloc::getMarkerClasses(orre2019u251))
)


# Create a data frame of the classes present and their associated number
classes_pres <- pRoloc::getMarkerClasses(orre2019a431)
classes_numeric_representation <- seq(1, length(classes_pres))
class_key <- data.frame(Class = classes_pres, Key = classes_numeric_representation)

test_markers <- orre2019a431_markers[proteins_captured, marker_col][test_inds]

numeric_truth <- class_key$Key[match(c(true_allocations), class_key$Class)] |>
  matrix(N, L)

# Initial classification for MCMC input
initial_labels <- matrix(1, N, L)
for (l in seq(1, L)) {
  initial_labels[test_inds, l] <- sample(classes_numeric_representation, size = N_test, replace = TRUE)
  initial_labels[train_inds, l] <- numeric_truth[train_inds, l]
}

# Place the expression data in a list to pass to MDI
X <- list(
  as.matrix(orre2019a431_exprs),
  as.matrix(orre2019h322_exprs),
  as.matrix(orre2019hcc827_exprs),
  # as.matrix(orre2019hcc827gef_exprs),
  # as.matrix(orre2019hcc827rep1_exprs),
  as.matrix(orre2019mcf7_exprs),
  as.matrix(orre2019u251_exprs)
)


# Do a fold of the cross validation for MDI and the single-view mixture model
types <- rep("TAGM", L)

# Used in the model fit matrices
eff_r <- floor(num_iter / thin) + 1
eff_burn <- floor(burn / thin) + 1
n_recorded <- eff_r - eff_burn
recorded_iterations <- seq(burn + thin, num_iter, thin)
recorded_iteration_indices <- seq(eff_burn + 1, eff_r, 1)

# === Fit models ===============================================================

cat("\n\n=== BEGIN MAIN FUNCTION ===========================================\n")

cat("\nMDI validation fold.\n")

mdi_chains <- mdir::runMCMCChains(X, n_chains, R, thin,
  types = types,
  K = K,
  fixed = fixed,
  initial_labels = initial_labels
)

mdi_perf_lst <- fitModel(mdi_chains,
  R,
  thin,
  burn,
  test_markers,
  test_inds,
  class_key,
  classes_pres,
  MDI_run = TRUE
)

mdi_perf_lst$Seed <- seed
mdi_perf_lst$likelihood_df$Fold <- seed

cat("\n\nMDI has run. Now run a TAGM model on the main dataset.\n")

tagm_perf_lst <- vector("list", L)
for (ll in seq(1, L)) {
  tagm_chains <- mdir::runMCMCChains(X[ll], n_chains, R, thin,
    types = types[l],
    K = K[l],
    fixed = fixed[, l, drop = FALSE],
    initial_labels = initial_labels[, l, drop = FALSE]
  )

  tagm_perf_lst[[ll]] <- fitModel(tagm_chains,
    R,
    thin,
    burn,
    test_markers,
    test_inds,
    class_key,
    classes_pres,
    MDI_run = FALSE
  )
  tagm_perf_lst[[ll]]$Seed <- seed
  tagm_perf_lst[[ll]]$likelihood_df$Fold <- seed
}

cat("\n\nTAGM has run.")

t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))

# Combine the two objects into a single list and save it
output <- list(
  mdi_fold = mdi_perf_lst,
  tagm_fold = tagm_perf_lst,
  test.idx = test_inds,
  seed = seed,
  time_taken = time_taken
)

cat("\nSaving output to:\n", list_save_name)
saveRDS(output, file = list_save_name)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")
