
library(tagmReDraft)

library(data.table)
library(pheatmap)
library(pRolocdata)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(mdiHelpR)
setMyTheme()

#' @title Prepare MS Object
#' @description Prepares a mass spectrometry experiment (as stored in
#' the Bioconductor package ``pRolocdata``) for modelling, extracting the
#' numerical data, the classes and the indicator matrix for which labels are
#' observed.
#' @param MS_object A mass spectrometry experiment such as ``tan2009r1`` from
#' ``pRolocdata``.
#' @return A list of ``X``, the fracitonation data from a LOPIT experiment,
#'  ``fixed``, the matrix indicating which labels are observed,
#'  ``initial_labels``, the matrix of the initial labels to be input into
#'  ``runMCMCChains`` (note that ``"unknown"`` organelles are represented
#'  arbitrarily with a ``1`` as these will be sampled again within the wrapper
#'  of ``callMDI``) and ``class_key`` which maps the numeric representation of
#'  organelles back to the original naming.
#' @export
prepareMSObject <- function(MS_object) {

  # Extract the LOPIT data and the organelles
  X <- Biobase::exprs(MS_object)
  organelles <- fData(MS_object)[, "markers"]

  # Create a data frame of the classes present and their associated number;\
  # this can be used to map the numeric representation of the classes back to
  # an organelle
  organelles_present <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(
    Organelle = organelles_present,
    Key = 1:length(organelles_present)
  )

  # Number of components modelled
  K <- length(organelles_present)

  # Number of views modelled
  V <- 1

  # Number of samples modelled
  N <- nrow(X)

  # Prepare initial labels
  initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)

  # Fix training points, allow test points to move component
  fix_vec <- (organelles != "unknown") * 1
  fixed[, 1] <- fix_vec

  # Assign initial labels
  initial_labels[, 1] <- class_key$Key[match(organelles, class_key$Organelle)]

  # Any unknown labels are given an arbitrary value which will be reset in the
  # model call function.
  initial_labels[is.na(initial_labels)] <- 1

  data_modelled <- list(
    X
  )

  # Return the prepared objects
  list(
    X = X,
    fixed = fixed,
    initial_labels = initial_labels,
    class_key = class_key
  )
}

data(Barylyuk2020ToxoLopit)

microarray_data <- fread("./T_gondii/ToxoDB_TgME49_Protein-coding_DNA_microarray.txt",
  na.strings = "N/A",
  strip.white = T,
  header = T,
  select = seq(1, 212)
)

rna_seq_data <- fread("./T_gondii/ToxoDB_TgME49_Protein-coding_RNA-Seq.txt",
  na.strings = "N/A",
  strip.white = T,
  header = T,
  select = seq(1, 255)
)

mismatching_order <- any(microarray_data[, 1] != rna_seq_data[, 1])

if (mismatching_order) {
  stop("Rownames not matching.")
}


# Check that we have no NAs in the data
nas_in_rna_seq <- any(apply(rna_seq_data, 2, function(x) {
  any(is.na(x))
}))

columns_containing_nas_in_microarray <- apply(microarray_data, 2, function(x) {
  sum(is.na(x))
})

rows_containing_nas_in_microarray <- apply(microarray_data, 1, function(x) {
  sum(is.na(x))
})

table(columns_containing_nas_in_microarray)
table(rows_containing_nas_in_microarray)

which_rows_containing_nas_in_microarray <- which(apply(microarray_data, 1, function(x) {
  any(is.na(x))
}))

cleaned_microarray_data <- microarray_data[-which_rows_containing_nas_in_microarray, ]

remaining_genes <- cleaned_microarray_data[[1]]
rna_genes <- rna_seq_data[[1]]

kept_rows <- match(remaining_genes, rna_genes)
cleaned_rna_seq_data <- rna_seq_data[kept_rows, ]

rna_seq_numeric_columns <- seq(4, ncol(cleaned_rna_seq_data))
microarray_numeric_columns <- seq(4, ncol(cleaned_microarray_data))

rna_mat <- log(as.matrix(cleaned_rna_seq_data[, ..rna_seq_numeric_columns]) + 1)
microarray_mat <- as.matrix(cleaned_microarray_data[, ..microarray_numeric_columns])

# pheatmap(rna_mat,
# show_colnames = F,
# show_rownames = F,
# main = "Bulk RNA-seq data"
# )
#
# pheatmap(microarray_mat, show_colnames = F, show_rownames = F, cluster_cols = F, main = "Cell cycle microarray data")

protein_lst <- prepareMSObject(Barylyuk2020ToxoLopit)
measurements <- colnames(protein_lst$X)
proteins_represented <- row.names(protein_lst$X)

protein_df <- data.frame(protein_lst$X) %>%
  mutate(
    Protein = row.names(protein_lst$X),
    Label = factor(protein_lst$initial_labels[, 1]),
    Fixed = protein_lst$fixed[, 1]
  )

long_protein_df <- protein_df %>%
  pivot_longer(-c(
    Protein,
    Label,
    Fixed
  ), names_to = "Fraction", values_to = "Value") %>%
  mutate(Fraction = factor(Fraction, levels = measurements, ordered = TRUE))

index_protein_in_rna_data <- match(remaining_genes, proteins_represented)
index_protein_in_rna_data <- index_protein_in_rna_data[!is.na(index_protein_in_rna_data)]
proteins_in_rna_data <- proteins_represented[index_protein_in_rna_data]

index_genes_in_protein_data <- match(proteins_in_rna_data, remaining_genes)
index_genes_in_protein_data <- index_genes_in_protein_data[!is.na(index_genes_in_protein_data)]
genes_in_protein_data <- remaining_genes[index_genes_in_protein_data]

mismatch_in_items <- any(genes_in_protein_data != proteins_in_rna_data)
if (mismatch_in_items) {
  stop("Still mismatch in items analysed.")
}

retained_rows <- match(genes_in_protein_data, cleaned_microarray_data$`Gene ID`)

final_microarray_data <- cleaned_microarray_data[retained_rows, ]
final_rna_seq_data <- cleaned_rna_seq_data[retained_rows, ]

final_protein_df <- protein_df[match(proteins_in_rna_data, protein_df$Protein), ]

all_items_matching <- (all(final_protein_df$Protein == final_rna_seq_data$`Gene ID`) &
  all(final_protein_df$Protein == final_microarray_data$`Gene ID`)
)

if (!all_items_matching) {
  stop("There's a mismatch in the order / membership of genes represented.")
}

rna_mat <- log(as.matrix(final_rna_seq_data[, ..rna_seq_numeric_columns]) + 1)
microarray_mat <- as.matrix(final_microarray_data[, ..microarray_numeric_columns])

# pheatmap(rna_mat,
# show_colnames = F,
# show_rownames = F,
# main = "Bulk RNA-seq data"
# )
#
# pheatmap(microarray_mat, show_colnames = F, show_rownames = F, cluster_cols = F, main = "Cell cycle microarray data")
#
# long_final_protein_df<- final_protein_df %>%
# pivot_longer(-c(Protein,
# Label,
# Fixed), names_to = "Fraction", values_to = "Value") %>%
# mutate(Fraction = factor(Fraction, levels = measurements, ordered = TRUE))
#
#
# long_final_protein_df %>%
# filter(Fixed == 1, Label %in% seq(1, 9)) %>%
# ggplot(aes(x = Fraction, y= Value, group = Protein)) +
# geom_line() +
# facet_wrap(~Label)
#
# long_final_protein_df %>%
# filter(Fixed == 1, Label %in% seq(10, 18)) %>%
# ggplot(aes(x = Fraction, y= Value, group = Protein)) +
# geom_line() +
# facet_wrap(~Label)
#
# long_final_protein_df %>%
# filter(Fixed == 1, Label %in% seq(19, 26)) %>%
# ggplot(aes(x = Fraction, y= Value, group = Protein)) +
# geom_line() +
# facet_wrap(~Label)
#
# pheatmap(microarray_mat[, c(1:27, 124:140)],
# show_colnames = F,
# show_rownames = F,
# cluster_cols = F,
# main = "Cell cycle microarray data; reduced to two experiments"
# )

protein_mat <- as.matrix(final_protein_df[, -c(31:33)])

data_modelled <- list(
  microarray_mat[, c(1:27, 124:140)],
  rna_mat,
  protein_mat
)

N <- nrow(rna_mat)
V <- length(data_modelled)

initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)

fixed[, 3] <- final_protein_df$Fixed
initial_labels[, 3] <- final_protein_df$Label

types <- c("MVN", "MVN", "TAGM")

K <- c(100, 100, 26)

R <- 5
thin <- 1
burn <- 2

mcmc_output <- callMDI(data_modelled,
  R = R,
  thin = thin,
  initial_labels = initial_labels,
  fixed = fixed,
  K = K,
  types = types
)

