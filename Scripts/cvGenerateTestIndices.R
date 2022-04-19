#!/usr/bin/env Rscript
#
# Title: Generate test indices for Cross validation (single fold)
# Description: This script generates the indices of  the proteins used as test
# data in a CV fold for a given dataset and seed. If instructed, it will ensure 
# that all classes are represented in the training data.
#
# Example call:
# Rscript cvGenerateTestIndices.R --data "tan2009r1" --seed 1
# --test_size 0.75 --save_dir "./" --adjust_for_TL TRUE

suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

# === Functions ================================================================


# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--dataset"),
      type = "character",
      help = "Which dataset should be analysed. E.g., 'tan2009r1'.",
      default = "tan2009r1",
      metavar = "character"
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
    optparse::make_option(c("--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save output to [default= %default]",
      metavar = "character"
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

# === Cross-validation =========================================================

cat("\n\n=== PASS INPUTS ===================================================\n")

t0 <- Sys.time()

# Pass the inputs from the command line
args <- input_arguments()

# Random seed used defining this fold
seed <- args$seed

# Fraction of data used for validation in fold
test_size <- args$test_size

# Adjust the trest indices to ensure all organelles are represented in the
# training data as required by the transfer learner
adjust_training_for_transfer_learner <- args$adjust_for_TL

is_true <- (adjust_training_for_transfer_learner == "true" | adjust_training_for_transfer_learner == "TRUE")
is_false <- (adjust_training_for_transfer_learner == "false" | adjust_training_for_transfer_learner == "FALSE")

if(is_false) {
  adjust_training_for_transfer_learner <- FALSE
}
if(is_true) {
  adjust_training_for_transfer_learner <- TRUE
}

# random seed
set.seed(seed)

# Convert the single string input into multiple possible strings
dataset <- args$dataset

test_size_str <- paste0(test_size * 100)

# Use the passed directory appened by the dataset. Create this directory if
# required.
save_dir <- paste0(args$save_dir, dataset, "/")
dir.create(save_dir, showWarnings = FALSE)

# Details of the save files
save_name <- paste0(
  save_dir,
  dataset,
  "_seed_",
  seed,
  "_testSize_",
  test_size_str,
  "_trainingAdjustedForTL_",
  adjust_training_for_transfer_learner
)

# What will the saved object be called
list_save_name <- paste0(save_name, ".rds")

cat("\nLoading data.")

# Load the relevant data
data(list = dataset)

# Convert to the correct format to be passed to the validation function
d1 <- eval(parse(text = dataset))

cat("\nData loaded.")

# Number of samples modelled
N <- nrow(d1)

# Possibly reduce the dimensionality of the categorical dataset depending on
# the number of entries in each column

# === Generate test indices ====================================================

# Use the same test indices across methods
marker.data <- pRoloc::markerMSnSet(d1)
X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)

# get sizes
# .size <- ceiling(table(MSnbase::fData(marker.data)$markers) * test_size)
organelle_representation <- table(MSnbase::fData(marker.data)$markers)
.size <- ceiling(organelle_representation * test_size)

if (adjust_training_for_transfer_learner) {
  someOrganellesAdjusted <- FALSE
  organelles_adjusted <- c()
  for (ii in seq(1, length(organelle_representation))) {
    no_training_members <- .size[ii] == organelle_representation[ii]
    if (no_training_members) {
      someOrganellesAdjusted <- TRUE
      org <- names(organelle_representation)[ii]
      organelles_adjusted <- c(organelles_adjusted, org)
      .size[ii] <- .size[ii] - 1
    }
  }
  if(someOrganellesAdjusted) {
    increased_training_size_str <- paste0(organelles_adjusted, collapse = ", ")
    cat(
      "\nNo training data for",
      increased_training_size_str,
      "at original test size (",
      test_size,
      ").\nCorrecting for the sake of the KNN transfer learner."
    )
  }
}

# strata needs size to be ordered as they appear in data
.size <- .size[unique(MSnbase::fData(marker.data)$markers)]

# get strata indices
test.idx <- sampling::strata(X, "markers",
  size = .size,
  method = "srswor"
)$ID_unit


# Combine the two objects into a single list and save it
output <- list(test.idx = test.idx, seed = seed)

cat("\nSaving output to:\n", list_save_name)
saveRDS(output, file = list_save_name)

cat("\n\n=== SCRIPT COMPLETE ===============================================\n")
