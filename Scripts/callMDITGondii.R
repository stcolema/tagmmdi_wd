
suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(pheatmap))

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


# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--data_dir"),
      type = "character",
      help = "Path to the directory containing the data.",
      metavar = "character"
    ),

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
    optparse::make_option(c("-b", "--burn"),
      type = "numeric",
      default = 1000,
      help = paste(
        "Number of iterations to burn of the warm-up period in each MCMC chain [default= %default]",
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
    optparse::make_option(c("--n_chains"),
      type = "numeric",
      default = 4L,
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
    )
  )


  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

my_heatmap <- function(X, annotation_labels, cluster_rows = FALSE, cluster_cols = FALSE, ...) {
  organelles_present <- sort(unique(annotation_labels$Organelle))

  K <- length(organelles_present)
  if (any(is.na(organelles_present))) {
    K <- K - 1
  }

  ann_colours <- list("Organelle" = viridis::viridis(K))
  names(ann_colours$Organelle) <- factor(sort(levels(annotation_labels$Organelle)))

  col_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)

  my_breaks <- defineDataBreaks(X, col_pal, mid_point = 0)
  if (min(X) >= 0.0) {
    my_breaks <- defineDataBreaks(X, col_pal, mid_point = median(X))
  }

  ordering <- c() # seq(1, nrow(X))
  for (org in organelles_present) {
    org_indices <- which(annotation_labels$Organelle == org)
    if (length(org_indices) > 1) {
      org_row_order <- mdiHelpR::findOrder(X[org_indices, ])
      org_indices <- org_indices[org_row_order]
    }
    ordering <- c(ordering, org_indices)
  }

  X <- X[ordering, ]
  annotation_labels <- annotation_labels[ordering, , drop = FALSE]

  ph <- pheatmap(X,
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    color = col_pal,
    breaks = my_breaks,
    annotation_row = annotation_labels,
    annotation_colors = ann_colours,
    ...
  )
  ph
}

# === T. gondii analysis =======================================================

cat("\n\n=== PASS INPUTS ===================================================\n")

t0 <- Sys.time()

# ggplot2 theme
setMyTheme()

# Pass the inputs from the command line
args <- input_arguments()

# Directories for input and output respectively
data_dir <- args$data_dir
save_dir <- args$save_dir

# Random seed used defining this fold
seed <- args$seed

# MCMC sampler arguments
R <- args$R
thin <- args$thin
burn <- args$burn

# Number of chains modelled
n_chains <- args$n_chains

# The number of components modelled
K <- args$K

# Number of clusters modelled in the categorical dataset
n_clust_unsupervised <- K
if (is.null(K)) {
  n_clust_unsupervised <- 100
}

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

microarray_file <- paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
rna_seq_file <- paste0(data_dir, "ToxoDB_TgME49_Protein-coding_RNA-Seq.txt")

data(Barylyuk2020ToxoLopit)

microarray_data <- fread(microarray_file,
  na.strings = "N/A",
  strip.white = T,
  header = T,
  select = seq(1, 212)
)

rna_seq_data <- fread(rna_seq_file,
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

# table(columns_containing_nas_in_microarray)
# table(rows_containing_nas_in_microarray)

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
protein_mat <- as.matrix(final_protein_df[, -c(31:33)])

gene_ids <- row.names(final_protein_df)
row.names(microarray_mat) <- row.names(rna_mat) <- gene_ids
fixed <- which(final_protein_df$Fixed == 1)

normaliseForCoexpression <- function(X) {
  na.omit(t(scale(t(X))))
}

rnaseq_macrophages_infected_by_T_gondii_inds <- seq(89, 142)
rnaseq_macrophages_infected_by_T_gondii <- rna_mat[, rnaseq_macrophages_infected_by_T_gondii_inds]
nonunique_rnaseq_macrophages_infected_by_T_gondii <- rnaseq_macrophages_infected_by_T_gondii[, seq(26, 54)]
unique_rnaseq_macrophages_infected_by_T_gondii <- rnaseq_macrophages_infected_by_T_gondii[, seq(1, 25)]

normalised_rnaseq_macrophages_infected_by_T_gondii <- normaliseForCoexpression(rnaseq_macrophages_infected_by_T_gondii)
normalised_nonuinque_rnaseq_macrophages_infected_by_T_gondii <- normaliseForCoexpression(nonunique_rnaseq_macrophages_infected_by_T_gondii)
normalised_unique_rnaseq_macrophages_infected_by_T_gondii <- normaliseForCoexpression(unique_rnaseq_macrophages_infected_by_T_gondii)



m_white_cell_cycle_inds <- seq(1, 14) # 27)
m_white_cell_cycle <- microarray_mat[, m_white_cell_cycle_inds]
m_white_cell_cycle_normalised <- m_white_cell_cycle %>%
  t() %>%
  scale() %>%
  t() %>%
  na.omit()

write.csv(m_white_cell_cycle_normalised, paste0(save_dir, "cellCycleNormalised.csv"))
write.csv(normalised_rnaseq_macrophages_infected_by_T_gondii, paste0(save_dir, "rnaSeqMacrophage.csv"))
write.csv(final_protein_df, paste0(save_dir, "LOPITreduced.csv"))

data_modelled <- list(
  m_white_cell_cycle_normalised,
  normalised_rnaseq_macrophages_infected_by_T_gondii,
  protein_mat
)

saveRDS(data_modelled, file = paste0(save_dir, "dataModelled.rds"))

N <- nrow(rna_mat)
V <- length(data_modelled)

initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)

fixed[, 3] <- final_protein_df$Fixed
initial_labels[, 3] <- final_protein_df$Label

# types <- c("MVN", "G", "TAGM")
types <- c("MVN", "MVN", "TAGM")

K <- c(
  n_clust_unsupervised,
  n_clust_unsupervised,
  length(pRoloc::getMarkerClasses(Barylyuk2020ToxoLopit))
)

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
