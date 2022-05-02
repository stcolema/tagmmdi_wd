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
# Rscript mdiTagmCVSingleLoop.R --datasets "tan2009r1 tan2009r1goCC" --seed 1
# --test_indices "./tan2009r1_tan2009r1goCC_R_15000_seed_1_nChains_5_testSize_75.rds"
# --save_dir "./" --categorical_column_threshold 5 --adjust_for_TL TRUE

suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressPackageStartupMessages(library(MSnbase))

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
    optparse::make_option(c("--datasets"),
      type = "character",
      help = "Which dataset should be analysed. E.g., 'tan2009r1 tan2009r1goCC'.",
      metavar = "character"
    ),

    # File containing the data regarding the test/train split of the data
    optparse::make_option(c("-t", "--test_indices"),
      type = "character",
      help = paste(
        "File containing the data regarding the test/train split of the data.",
        "Output of ``mdiTagmCVSingleLoop.R``"
      ),
      metavar = "character"
    ),

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
    optparse::make_option(c("-c", "--categorical_column_threshold"),
      type = "numeric",
      default = 0L,
      help = "Required percentage representation in category [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--number_weights_sampled"),
      type = "numeric",
      default = 500L,
      help = "Number of transfer weights actually considered, randomly sampled from the grid generated using the pRoloc functions [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-n", "--number_weights"),
      type = "numeric",
      default = 4L,
      help = "Number of transfer weights considered [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("--adjust_for_TL"),
      type = "logical",
      default = FALSE,
      help = paste0(
        "Adjust the test indices to ensure all classes are ",
        "represented in the training data (as required by the Transfer learner)",
        " [default= %default]"
      ),
      metavar = "logical"
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


weightCombinations <- function(n_classes, n_weights, n_used = NULL) {
  weights <- seq(0, 1, by = 1 / (n_weights - 1))
  total_combinations <- n_weights**n_classes
  if (n_used > total_combinations) {
    warning_message <- paste(
      "Number of requested combinations exceeds the number of possible unique",
      "combinations.\nReducing to the largest possible amount:",
      total_combinations
    )
    message(warning_message)
    n_used <- total_combinations
  }
  combinations_used <- sample(seq(1, total_combinations), size = n_used)

  combinations_matrix <- matrix(0, nrow = n_used, ncol = n_classes)
  for (ii in seq(n_classes)) {
    combinations_matrix[, ii] <- sample(weights, size = n_used, replace = TRUE)
  }
  duplicated_rows <- duplicated.matrix(combinations_matrix)
  if (any(duplicated_rows)) {
    combinations_matrix <- combinations_matrix[-duplicated_rows, ]
  }
  combinations_matrix
}

prepInputsForTransferLearner <- function(MS_object,
                                         MS_cat_object,
                                         test.idx,
                                         categorical_column_threshold = 0) {

  # cat("\nPreparing inputs.")

  marker.data <- pRoloc::markerMSnSet(MS_object)

  # cat("\nMarker data obtained.")

  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)
  K <- length(pRoloc::getMarkerClasses(MS_object))

  # Create a data frame of the classes present and their associated number
  classes_pres <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(Class = classes_pres, Key = 1:length(classes_pres))

  # cat("\nAccessing categorical data.\n")

  marker.data.cat <- pRoloc::markerMSnSet(MS_cat_object)
  GO_df <- pRoloc:::subsetAsDataFrame(marker.data.cat, "markers", train = TRUE)

  cat("\nSplitting data into training and test sets.")
  cat("\nsize of test set:", (length(test.idx)))
  cat("\nSize of full dataset:", nrow(MSnbase::exprs(marker.data)))

  ## 'unseen' test set
  .test <- MSnbase::MSnSet(
    MSnbase::exprs(marker.data)[test.idx, ],
    MSnbase::fData(marker.data)[test.idx, ],
    MSnbase::pData(marker.data)
  )

  ## 'seen' training set
  .train <- MSnbase::MSnSet(
    MSnbase::exprs(marker.data)[-test.idx, ],
    MSnbase::fData(marker.data)[-test.idx, ],
    MSnbase::pData(marker.data)
  )

  # Samples in each split
  N_test <- nrow(.test)
  N_train <- nrow(.train)
  # test_indices_in_new_data <- seq(N_train + 1, N_test + N_train, 1)

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
  
  # Record internally which proteins are testing/training data
  fData(main_data)$test.protein <- c(rep(0, N_train), rep(1, N_test))

  # Set levels of markers cateogries
  levels(MSnbase::fData(main_data)$markers) <- c(
    levels(
      MSnbase::fData(main_data)$markers
    ),
    "unknown"
  )

  # hide marker labels in the column used for prediction, record truth too
  MSnbase::fData(main_data)$true.markers <- MSnbase::fData(main_data)$markers
  MSnbase::fData(main_data)[row.names(.test), "markers"] <- "unknown"

  # cat("\nEnsure auxiliary data has the same ordering as the main dataset.")

  # create new combined MSnset
  auxiliary_data <- BiocGenerics::combine(
    marker.data.cat[-test.idx, ],
    marker.data.cat[test.idx, ]
  )

  cat("\nRemoving uninformative columns from the GO dataset.")
  cat("\nThreshold:", categorical_column_threshold)
  cat("\nInitial number of columns:", ncol(auxiliary_data))
  informative_terms <- colSums(exprs(auxiliary_data)) > categorical_column_threshold
  if (any(!informative_terms)) {
    auxiliary_data <- auxiliary_data[, informative_terms]
  }
  cat("\nNumber of columns after reduction:", ncol(auxiliary_data))


  # Set levels of markers cateogries
  levels(MSnbase::fData(auxiliary_data)$markers) <- c(
    levels(
      MSnbase::fData(auxiliary_data)$markers
    ),
    "unknown"
  )

  # hide marker labels
  MSnbase::fData(auxiliary_data)[row.names(.test), "markers"] <- "unknown"

  if (any(fData(auxiliary_data)[, "markers"] != fData(main_data)[, "markers"])) {
    stop("\n\nERROR: Main and auxiliary data have different markers.\n")
  }

  if (any(row.names(auxiliary_data) != row.names(main_data))) {
    stop("\n\nERROR: Main and auxiliary data have differing row names.\n")
  }

  # cat("\nList of datasets prepared for model call.")

  data_modelled <- list(
    main_data = main_data,
    auxiliary_data = auxiliary_data
  )

  # cat("\nOutput list of various inputs to the MDI model.")

  list(

    # Inputs for KNN Transfer learner model
    data_modelled = data_modelled,

    # Objects used in assessing performance
    test.markers = test.markers,
    test.labels = test.labels,
    classes_pres = classes_pres,
    class_key = class_key
  )
}


knnSingleFold <- function(MS_object,
                          MS_cat_object,
                          test.idx,
                          number_weights,
                          seed,
                          number_weights_sampled = NULL,
                          categorical_column_threshold = 0) {
  inputs <- prepInputsForTransferLearner(MS_object, MS_cat_object, test.idx,
    categorical_column_threshold = categorical_column_threshold
  )

  d1 <- inputs$data_modelled$main_data
  d2 <- inputs$data_modelled$auxiliary_data

  # Map between organelle as a string and numeric representation
  class_key <- inputs$class_key
  classes_pres <- inputs$classes_pres

  # Used in assessing performance
  # test.markers <- inputs$test.markers

  # Define the weights to be explored
  f_data_col <- "markers"
  m <- sort(unique(fData(MS_object)[["markers"]])) # $markers.tl)

  # if(is.null(m)) {
  #   f_data_col <- "markers.tl"
  #   m <- sort(unique(fData(MS_object)[[f_data_col]])) # $markers)
  # }

  m <- m[m != "unknown"]

  cat("\nNumber of classes:", length(m))

  # As for the HEK dataset theta explodes in size (and memory) and gets killed
  # on the HPC, we use this function which does not guarantee that the number of
  # weights sampled is actually the desired amount as there can be reptition and
  # we reduce to the unique rows.
  cat("\nSampling weight combinations.")
  th <- weightCombinations(length(m), number_weights,
    n_used = number_weights_sampled
  )

  # n_combinations <- number_weights**length(m)
  # subsetting_of_weights_intended <- ! is.null(number_weights_sampled)
  # if(subsetting_of_weights_intended) {
  #   cat("\nSampling weight combinations.")
  #   number_weights_sampled <- min(number_weights_sampled, n_combinations)
  #   subset_of_weights_considered <- number_weights_sampled != n_combinations
  #   if(subset_of_weights_considered) {
  #     th_used <- sample(seq(1, n_combinations),
  #       size = number_weights_sampled,
  #       replace = FALSE
  #     )
  #     th <- thetas(length(m), length.out = number_weights, verbose = FALSE)[th_used, ]
  #   } else {
  #     th <- thetas(length(m), length.out = number_weights, verbose = FALSE)
  #   }
  # } else {
  #   th <- thetas(length(m), length.out = number_weights, verbose = FALSE)
  # }

  t0 <- Sys.time()

  cat("\nFinding choice of k for main dataset.")

  # find the best choice of k for the knn part of the transfer learner
  kopt <- knnOptimisation(d1,
    fcol = "markers",
    times = 100,
    k = seq(3, 20, 2),
    verbose = FALSE,
    seed = seed
  )

  best_k_main <- getParams(kopt)
  cat("\nChoice of k for main dataset:", best_k_main)

  cat("\n\nFinding choice of k for auxiliary dataset.")

  kopt <- knnOptimisation(d2,
    fcol = "markers",
    times = 100,
    k = seq(3, 20, 2),
    verbose = FALSE,
    seed = seed
  )

  best_k_aux <- getParams(kopt)
  cat("\nChoice of k for auxiliary dataset:", best_k_aux)

  cat("\n\nFinding best transfer weights for transfer learning algorithm.")

  # Find the best transfer weights
  topt <- knntlOptimisation(
    primary = d1,
    auxiliary = d2,
    th = th,
    k = c(best_k_main, best_k_aux),
    fcol = "markers",
    times = 50,

    # We only use a single thread on the HPC
    BPPARAM = BiocParallel::SerialParam(),
    seed = seed
  )

  bw <- getParams(topt)

  cat("\nBest transfer weights found.\n", bw)
  cat("\n\nPerforming classification.")

  ## Applying best *theta* weights {#sec:thclass}
  # Perform the final prediction
  d1 <- knntlClassification(d1, d2,
    bestTheta = bw,
    k = c(best_k_main, best_k_aux),
    fcol = "markers" # "markers.tl"
  )

  t1 <- Sys.time()
  time_taken <- t1 - t0

  d1 <- getPredictions(d1, fcol = "knntl")

  cat("\nPrediction complete.\n")

  predicted_organelle <- fData(d1)$knntl
  predicted_class <- class_key$Key[match(predicted_organelle, class_key$Class)]

  true_test_ids <- which(fData(d1)$test.protein == 1)
  N_test <- length(true_test_ids)
  true_markers <- fData(d1)$true.markers[true_test_ids]
  predicted_markers <- fData(d1)$knntl[true_test_ids]
  
  # True allocation for test data
  reference <- factor(true_markers, levels = classes_pres)
  comparison <- factor(predicted_markers, levels = classes_pres)

  cat("\nCompute model performance scores.")

  # Calculate the accuracy and the per-class, macro and weighted F1 scores.
  prediction_scores <- multiClassF1(comparison, reference)

  cat("\nCalculate the Brier score.")

  # Create allocation matrices for truth, filled initially with 0's
  true_allocation_matrix <- knn_allocation_matrix <- matrix(0,
    nrow = N_test,
    ncol = length(classes_pres)
  )

  # The numbers associated with the classes (i.e. numerical representation of
  # the classes)
  class_numerics <- seq(1, length(unique(true_markers)))

  # create allocation matrix
  for (j in seq_along(true_test_ids)) {
    # The class the current individual belongs to
    true_alloc <- as.numeric(true_markers[j], class_numerics)
    knn_alloc <- as.numeric(predicted_markers[j], class_numerics)
    # alloc <- as.numeric(test.markers, class_numerics)[j]

    # Enter this in the allocation matrix
    true_allocation_matrix[j, true_alloc] <- 1
    knn_allocation_matrix[j, knn_alloc] <- 1
  }

  # Compute quadratic loss
  quadloss <- sum((true_allocation_matrix - knn_allocation_matrix)^2)

  # Return the chains, the point estimates, the true organelles, the quadratic
  # loss score and the prediction scores
  list(
    MS_object = d1,
    predicted_organelle = predicted_organelle,
    predicted_class = predicted_class,
    truth = reference,
    quadloss = quadloss,
    prediction_scores = prediction_scores,
    class_key = class_key,
    time_taken = time_taken
  )
}

# setStockcol(paste0(getStockcol(), "80"))
# ptsze <- exp(fData(tan2009r1)$knntl.scores) - 1
# plot2D(tan2009r1, fcol = "knntl", cex = ptsze)
# setStockcol(NULL)
# addLegend(tan2009r1,
#   where = "topright",
#   fcol = "markers.tl",
#   bty = "n", cex = .7
# )


# === Cross-validation =========================================================

cat("\n\n=== PASS INPUTS ===================================================\n")

t0 <- Sys.time()

# ggplot2 theme
setMyTheme()

# Pass the inputs from the command line
args <- input_arguments()

# Random seed used defining this fold
seed <- args$seed

# Fraction of data used for validation in fold
test_indices_obj_filename <- args$test_indices

# The minimum number of entries required for a column of the GO term data
# to be included
categorical_column_threshold <- args$categorical_column_threshold

# Number of weights considered for each class in the transfer elarning algorithm
number_weights <- args$number_weights

# For computational reasons we consider a random sample of the
# generated weight combinations
number_weights_sampled <- args$number_weights_sampled

# Adjust the trest indices to ensure all organelles are represented in the
# training data as required by the transfer learner
adjust_training_for_transfer_learner <- args$adjust_for_TL
adjust_for_TL <- adjust_training_for_transfer_learner #

true_string_passed <- (adjust_training_for_transfer_learner == "true" ||
  adjust_training_for_transfer_learner == "TRUE"
)

false_string_passed <- (adjust_training_for_transfer_learner == "false" ||
  adjust_training_for_transfer_learner == "FALSE"
)

if (true_string_passed) {
  adjust_for_TL <- TRUE
}
if (false_string_passed) {
  adjust_for_TL <- FALSE
}

# random seed
set.seed(seed)

# Convert the single string input into multiple possible strings
datasets <- unlist(stringr::str_split(args$datasets, " "))

# Read in the object containing the test indices for the train/test split
test_indices_obj <- readRDS(test_indices_obj_filename)
test.idx <- test_indices_obj$test.idx

# Use the passed directory appened by the dataset. Create this directory if
# required.
save_dir <- paste0(args$save_dir, datasets[1], "/")
dir.create(save_dir, showWarnings = FALSE)

# Details of the save files
save_name <- paste0(
  save_dir,
  paste(datasets, collapse = "_"),
  "_knnTL",
  "_numberWeights_",
  number_weights,
  "_seed_",
  seed,
  "_trainingAdjustedForTL_",
  adjust_for_TL
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


# marker.data <- pRoloc::markerMSnSet(d1)
# cat("\nMarker data obtained.")
# X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)

# Number of samples modelled
# N <- nrow(X)

# d2 <- d2[, colSums(exprs(d2)) > categorical_column_threshold]

cat("\n\n=== BEGIN MAIN FUNCTION ===========================================\n")

cat("\nKNN validation fold.\n")

# Do a fold of the cross validation for MDI and the single-view mixture model
knn_fold <- knnSingleFold(
  MS_object = d1,
  MS_cat_object = d2,
  test.idx = test.idx,
  number_weights = number_weights,
  seed = seed,
  number_weights_sampled = number_weights_sampled,
  categorical_column_threshold = categorical_column_threshold
)

cat("\n\nKNN transfer learner has run.\n")

cat("\nSaving outputs.")

# Combine the two objects into a single list and save it
output <- list(
  knn_fold = knn_fold,
  test.idx = test.idx,
  seed = seed
)

saveRDS(output, file = list_save_name)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")

t1 <- Sys.time()
time_taken <- t1 - t0

cat("\n\nTIME TAKEN:", round(time_taken, 2), attr(time_taken, "units"))
