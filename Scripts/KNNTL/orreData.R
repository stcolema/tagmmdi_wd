
library(mdir)
library(pRoloc)
library(pRolocdata)
library(ggplot2)
library(tidyr)

subsetToKnownMarkers <- function(MSnSet_obj) {
  .X <- pRoloc:::subsetAsDataFrame(MSnSet_obj, "markers", train = TRUE)
  drop_na(.X)
}

set.seed(1)

data("orre2019a431")
data("orre2019h322")
data("orre2019hcc827")
data("orre2019hcc827gef")
data("orre2019hcc827rep1")
data("orre2019mcf7")
data("orre2019u251")

L <- 7

orre2019a431_markers <- subsetToKnownMarkers(orre2019a431) # "markers", train = TRUE)
orre2019h322_markers <- subsetToKnownMarkers(orre2019h322)
orre2019hcc827_markers <- subsetToKnownMarkers(orre2019hcc827)
orre2019hcc827gef_markers <- subsetToKnownMarkers(orre2019hcc827gef)
orre2019hcc827rep1_markers <- subsetToKnownMarkers(orre2019hcc827rep1)
orre2019mcf7_markers <- subsetToKnownMarkers(orre2019mcf7)
orre2019u251_markers <- subsetToKnownMarkers(orre2019u251)

orre2019a431_proteins <- row.names(orre2019a431_markers)
orre2019h322_proteins <- row.names(orre2019h322_markers)
orre2019hcc827_proteins <- row.names(orre2019hcc827_markers)
orre2019hcc827gef_proteins <- row.names(orre2019hcc827gef_markers)
orre2019hcc827rep1_proteins <- row.names(orre2019hcc827rep1_markers)
orre2019mcf7_proteins <- row.names(orre2019mcf7_markers)
orre2019u251_proteins <- row.names(orre2019u251_markers)

proteins_captured <- Reduce(
  intersect,
  list(
    orre2019a431_proteins,
    orre2019h322_proteins,
    orre2019hcc827_proteins,
    orre2019hcc827gef_proteins,
    orre2019hcc827rep1_proteins,
    orre2019mcf7_proteins,
    orre2019u251_proteins
  )
)

n_measurements <- ncol(orre2019a431_markers) - 1
measurements <- seq(1, n_measurements)
marker_col <- n_measurements + 1

orre2019a431_exprs <- orre2019a431_markers[proteins_captured, measurements]
orre2019h322_exprs <- orre2019h322_markers[proteins_captured, measurements]
orre2019hcc827_exprs <- orre2019hcc827_markers[proteins_captured, measurements]
orre2019hcc827gef_exprs <- orre2019hcc827gef_markers[proteins_captured, measurements]
orre2019hcc827rep1_exprs <- orre2019hcc827rep1_markers[proteins_captured, measurements]
orre2019mcf7_exprs <- orre2019mcf7_markers[proteins_captured, measurements]
orre2019u251_exprs <- orre2019u251_markers[proteins_captured, measurements]

N <- length(proteins_captured)

true_allocations <- matrix(c(
  orre2019a431_markers[proteins_captured, marker_col],
  orre2019h322_markers[proteins_captured, marker_col],
  orre2019hcc827_markers[proteins_captured, marker_col],
  orre2019hcc827gef_markers[proteins_captured, marker_col],
  orre2019hcc827rep1_markers[proteins_captured, marker_col],
  orre2019mcf7_markers[proteins_captured, marker_col],
  orre2019u251_markers[proteins_captured, marker_col]
), N, L)

N_test <- floor(0.75 * N)

test_inds <- sample(seq(1, N), replace = FALSE, size = N_test)
train_inds <- seq(1, N)[-test_inds]

fixed <- matrix(0, N, L)
fixed[train_inds, ] <- 1

K <- c(
  length(pRoloc::getMarkerClasses(orre2019a431)),
  length(pRoloc::getMarkerClasses(orre2019h322)),
  length(pRoloc::getMarkerClasses(orre2019hcc827)),
  length(pRoloc::getMarkerClasses(orre2019hcc827gef)),
  length(pRoloc::getMarkerClasses(orre2019hcc827rep1)),
  length(pRoloc::getMarkerClasses(orre2019mcf7)),
  length(pRoloc::getMarkerClasses(orre2019u251))
)


# Create a data frame of the classes present and their associated number
classes_pres <- pRoloc::getMarkerClasses(orre2019a431)
classes_numeric_representation <- seq(1, length(classes_pres))
class_key <- data.frame(Class = classes_pres, Key = classes_numeric_representation)

numeric_truth <- class_key$Key[match(c(true_allocations), class_key$Class)] |>
  matrix(N, L)


initial_labels <- matrix(1, N, L)

for (l in seq(1, L)) {
  initial_labels[test_inds, l] <- sample(classes_numeric_representation, size = N_test, replace = TRUE)
  initial_labels[train_inds, l] <- numeric_truth[train_inds, l]
}

X <- list(
  as.matrix(orre2019a431_exprs),
  as.matrix(orre2019h322_exprs),
  as.matrix(orre2019hcc827_exprs),
  as.matrix(orre2019hcc827gef_exprs),
  as.matrix(orre2019hcc827rep1_exprs),
  as.matrix(orre2019mcf7_exprs),
  as.matrix(orre2019u251_exprs)
)

R <- 10
burn <- 2
thin <- 1
types <- rep("TAGM", L)

samples <- mdir::callMDI(X, R, thin,
                         types = types,
                         K = K,
                         fixed = fixed,
                         initial_labels = initial_labels
)

processed_samples <- mdir::processMCMCChain(samples, burn, point_estimate_method = "mean")

lapply(processed_samples$fusion_probabilities, mean)

processed_samples$phi[, 1] |> hist()
processed_samples$phi[, 2] |> hist()
processed_samples$phi[, 3] |> hist()
processed_samples$phi[, 4] |> hist()
processed_samples$phi[, 5] |> hist()
processed_samples$phi[, 6] |> hist()

pred <- processed_samples$pred |> unlist() |> matrix(N, L)

res_gtagm <- mdi_acc <- colMeans(pred[test_inds, ] == numeric_truth[test_inds, ])

tagm_acc <- rep(0, L)
for(l in seq(1, L)) {
  tagm_smaples <- callMDI(X[l], R, thin,
                          type = types[l],
                          K = K[l],
                          fixed = fixed[, l, drop = FALSE],
                          initial_labels = initial_labels[, l, drop = FALSE])
  
  proc_tagm_smaples <- processMCMCChain(tagm_smaples, burn, "mean")
  tagm_pred <- proc_tagm_smaples$pred |> unlist() |> matrix(N, 1)
  tagm_acc[l] <- colMeans(tagm_pred[test_inds, , drop = FALSE] == numeric_truth[test_inds, l, drop = FALSE])
}

cc_samples <- mdir::runMCMCChainsSavingToFile(X, 1000, 10, 10, 
                                              dir_path = "/home/stephen/Desktop/Test/", types = types, K, initial_labels, fixed)
