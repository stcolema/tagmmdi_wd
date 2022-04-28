
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
  if(min(X) >= 0.0) 
    my_breaks <- defineDataBreaks(X, col_pal, mid_point = median(X))
  
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

gene_ids <- row.names(final_protein_df)
row.names(microarray_mat) <- row.names(rna_mat) <- gene_ids
col_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
# my_breaks <- defineDataBreaks(m_white_cell_cycle_normalised, col_pal, mid_point = 0)

fixed <- which(final_protein_df$Fixed == 1)

# data_modelled <- list(
#   microarray_mat[, c(1:27, 124:140)],
#   rna_mat,
#   protein_mat
# )

rnaseq_exp1_inds <- 1:4
rnaseq_exp2_inds <- 5:16
rnaseq_exp3_inds <- 17:32
rnaseq_exp4_inds <- 33:36
rnaseq_exp5_inds <- 37:48
rnaseq_exp6_inds <- 49:72
rnaseq_exp7_inds <- 73:76
rnaseq_exp8_inds <- 77:88
rnaseq_exp9_inds <- 89:142
rnaseq_exp10_inds <- 143:148
rnaseq_exp11_inds <- 149:152
rnaseq_exp12_inds <- 153:160
rnaseq_exp13_inds <- 161:184
rnaseq_exp14_inds <- 185:212
rnaseq_exp15_inds <- 213:220
rnaseq_exp16_inds <- 221:236
rnaseq_exp17_inds <- 237:240
rnaseq_exp18_inds <- 241:252

rnaseq_exp1 <- rna_mat[, rnaseq_exp1_inds]
rnaseq_exp2 <- rna_mat[, rnaseq_exp2_inds]
rnaseq_exp3 <- rna_mat[, rnaseq_exp3_inds]
rnaseq_exp4 <- rna_mat[, rnaseq_exp4_inds]
rnaseq_exp5 <- rna_mat[, rnaseq_exp5_inds]
rnaseq_exp6 <- rna_mat[, rnaseq_exp6_inds]
rnaseq_exp7 <- rna_mat[, rnaseq_exp7_inds]
rnaseq_exp8 <- rna_mat[, rnaseq_exp8_inds]
rnaseq_exp9 <- rna_mat[, rnaseq_exp9_inds]
rnaseq_exp10 <- rna_mat[, rnaseq_exp10_inds]
rnaseq_exp11 <- rna_mat[, rnaseq_exp11_inds]
rnaseq_exp12 <- rna_mat[, rnaseq_exp12_inds]
rnaseq_exp13 <- rna_mat[, rnaseq_exp13_inds]
rnaseq_exp14 <- rna_mat[, rnaseq_exp14_inds]
rnaseq_exp15 <- rna_mat[, rnaseq_exp15_inds]
rnaseq_exp16 <- rna_mat[, rnaseq_exp16_inds]
rnaseq_exp17 <- rna_mat[, rnaseq_exp17_inds]
rnaseq_exp18 <- rna_mat[, rnaseq_exp18_inds]

pheatmap(rnaseq_exp18, show_colnames = FALSE)

normaliseForCoexpression <- function(X) {
  na.omit(t(scale(t(X))))
}

normalised_rnaseq_exp1 <- normaliseForCoexpression(rnaseq_exp1)
normalised_rnaseq_exp2 <- normaliseForCoexpression(rnaseq_exp2)
normalised_rnaseq_exp3 <- normaliseForCoexpression(rnaseq_exp3)
normalised_rnaseq_exp4 <- normaliseForCoexpression(rnaseq_exp4)
normalised_rnaseq_exp5 <- normaliseForCoexpression(rnaseq_exp5)
normalised_rnaseq_exp6 <- normaliseForCoexpression(rnaseq_exp6)
normalised_rnaseq_exp7 <- normaliseForCoexpression(rnaseq_exp7)
normalised_rnaseq_exp8 <- normaliseForCoexpression(rnaseq_exp8)
normalised_rnaseq_exp9 <- normaliseForCoexpression(rnaseq_exp9)
normalised_rnaseq_exp10 <- normaliseForCoexpression(rnaseq_exp10)
normalised_rnaseq_exp11 <- normaliseForCoexpression(rnaseq_exp11)
normalised_rnaseq_exp12 <- normaliseForCoexpression(rnaseq_exp12)
normalised_rnaseq_exp13 <- normaliseForCoexpression(rnaseq_exp13)
normalised_rnaseq_exp14 <- normaliseForCoexpression(rnaseq_exp14)
normalised_rnaseq_exp15 <- normaliseForCoexpression(rnaseq_exp15)
normalised_rnaseq_exp16 <- normaliseForCoexpression(rnaseq_exp16)
normalised_rnaseq_exp17 <- normaliseForCoexpression(rnaseq_exp17)
normalised_rnaseq_exp18 <- normaliseForCoexpression(rnaseq_exp18)

my_heatmap(rnaseq_exp1[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp2[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp3[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp4[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp5[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp6[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp7[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp8[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp9[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp10[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp11[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp12[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp13[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp14[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp15[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp16[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp17[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)
my_heatmap(rnaseq_exp18[fixed, ], annotation_labels[fixed, , drop = FALSE], cluster_cols = TRUE)

my_heatmap(normalised_rnaseq_exp1[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp2[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp3[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp4[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp5[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp6[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp7[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp8[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp9[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp10[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp11[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp12[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp13[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp14[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp15[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp16[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp17[fixed, ], annotation_labels[fixed, , drop = FALSE])
my_heatmap(normalised_rnaseq_exp18[fixed, ], annotation_labels[fixed, , drop = FALSE])


colnames(normalised_rnaseq_exp9)

my_heatmap(normalised_rnaseq_exp9[fixed, 1:25], annotation_labels[fixed, , drop = FALSE], cluster_rows = TRUE, cluster_cols = TRUE)
my_heatmap(normalised_rnaseq_exp9[fixed, -seq(1, 25)], annotation_labels[fixed, , drop = FALSE], cluster_rows = TRUE, cluster_cols = TRUE)
my_heatmap(normalised_rnaseq_exp9, annotation_labels, cluster_rows = TRUE, cluster_cols = TRUE)


my_heatmap(normalised_rnaseq_exp18[fixed, ],
           annotation_labels[fixed, , drop = FALSE],
           # main = "Behnke et al. (2010)"
           # filename = "./T_gondii/Plots/behnke_cell_cycle_heatmap_fixed.png"
)


m_white_cell_cycle_inds <- seq(1, 14) # 27)
ME49_Bradyzoite_Differentiation_inds <- RH_delta_HXGPRT_delta_UPRT_strain_alkaline_bradyzoite_inds <- seq(132, 140)
Tachyzoite_transcriptome_during_invasion_inds <- TgRH_Tachy_invade_Marray_inds <- seq(161, 163)
Bradyzoite_Differentiation_72_hrs_inds <- seq(150, 153)

m_white_cell_cycle <- microarray_mat[, m_white_cell_cycle_inds]
ME49_Bradyzoite_Differentiation <- microarray_mat[, ME49_Bradyzoite_Differentiation_inds]
Tachyzoite_transcriptome_during_invasion <- microarray_mat[, Tachyzoite_transcriptome_during_invasion_inds]
Bradyzoite_Differentiation_72_hrs <- microarray_mat[, Bradyzoite_Differentiation_72_hrs_inds]

m_white_cell_cycle_normalised <- m_white_cell_cycle %>%
  t() %>%
  scale() %>%
  t() %>%
  na.omit()

ME49_Bradyzoite_Differentiation_normalised <- ME49_Bradyzoite_Differentiation %>%
  t() %>%
  scale() %>%
  t() %>%
  na.omit()

Tachyzoite_transcriptome_during_invasion_normalised <- Tachyzoite_transcriptome_during_invasion %>%
  t() %>%
  scale() %>%
  t() %>%
  na.omit()

dropped_ind <- is.na(t(scale(t(Tachyzoite_transcriptome_during_invasion)))) %>%
  apply(1, any) %>%
  which()

Bradyzoite_Differentiation_72_hrs_normalised <- Bradyzoite_Differentiation_72_hrs %>%
  t() %>%
  scale() %>%
  t() %>%
  na.omit()

dim(m_white_cell_cycle_normalised)
dim(ME49_Bradyzoite_Differentiation_normalised)
dim(Tachyzoite_transcriptome_during_invasion_normalised)
dim(Bradyzoite_Differentiation_72_hrs_normalised)

Bradyzoite_Differentiation_72_hrs_normalised[match(row.names(Bradyzoite_Differentiation_72_hrs_normalised), row.names(Tachyzoite_transcriptome_during_invasion_normalised)), ]

# protein_lst$class_key$Organelle[final_protein_df$Label]



# annotatedHeatmap(m_white_cell_cycle_normalised, annotation_labels)


K <- length(unique(annotation_labels$Organelle)) - 1

ann_colours <- list("Organelle" = viridis::viridis(K))
names(ann_colours$Organelle) <- factor(sort(levels(annotation_labels$Organelle)))

# ann_colours$Organelle[K] <-

nrow(Tachyzoite_transcriptome_during_invasion_normalised)

my_heatmap(m_white_cell_cycle_normalised[fixed, ],
  annotation_labels[fixed, , drop = FALSE],
  main = "Behnke et al. (2010)",
  filename = "./T_gondii/Plots/behnke_cell_cycle_heatmap_fixed.png"
)

my_heatmap(ME49_Bradyzoite_Differentiation_normalised[fixed, ],
  annotation_labels[fixed, , drop = FALSE],
  main = "Bradyzoite inducing conditions",
  filename = "./T_gondii/Plots/brad_conditions_heatmap_fixed.png"
)


my_heatmap(Bradyzoite_Differentiation_72_hrs_normalised[fixed, ],
  annotation_labels[fixed, , drop = FALSE],
  main = "Bradyzoite Differentiation (days 0 - 3)",
  filename = "./T_gondii/Plots/brad_diff_heatmap_fixed.png"
)

Tachyzoite_labels <- annotation_labels[-dropped_ind, , drop = FALSE]
Tachyzoite_fixed <- fixed[-dropped_ind]

my_heatmap(Tachyzoite_transcriptome_during_invasion_normalised[Tachyzoite_fixed, ],
  Tachyzoite_labels[Tachyzoite_fixed, , drop = FALSE],
  main = "Gaji et al. (2011)",
  filename = "./T_gondii/Plots/gaji_invasion_heatmap_fixed.png"
)


m_white_cell_cycle_normalised %>%
  pheatmap(
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    color = col_pal,
    breaks = my_breaks,
    annotation_row = annotation_labels,
    annotation_colors = ann_colours,
    main = "Behnke et al. (2010)",
    filename = "./T_gondii/Plots/behnke_cell_cycle_heatmap.png"
  )

ME49_Bradyzoite_Differentiation_normalised %>%
  pheatmap(
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    color = col_pal,
    breaks = my_breaks,
    annotation_row = annotation_labels,
    annotation_colors = ann_colours,
    main = "Bradyzoite inducing conditions",
    filename = "./T_gondii/Plots/brad_conditions_heatmap.png"
  )

Tachyzoite_transcriptome_during_invasion_normalised %>%
  pheatmap(
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    color = col_pal,
    breaks = my_breaks,
    annotation_row = annotation_labels,
    annotation_colors = ann_colours,
    main = "Gaji et al. (2011)",
    filename = "./T_gondii/Plots/gaji_invasion_heatmap.png"
  )

Bradyzoite_Differentiation_72_hrs_normalised %>%
  pheatmap(
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    color = col_pal,
    breaks = my_breaks,
    annotation_row = annotation_labels,
    annotation_colors = ann_colours,
    main = "Bradyzoite Differentiation (days 0 - 3)",
    filename = "./T_gondii/Plots/brad_diff_heatmap.png",
    width = 6,
    height = 8
  )

# pheatmap(m_white_cell_cycle_normalised)

# m_white_spline_cell_cycle <- seq(28, 87)
# Pru_dHXGPRT_strain_CO2_starvation_bradyzoite <- seq(88, 97)
# Pru_VEF_strain_CO2_starvation_bradyzoite <- seq(98, 102)
# Pru_dHXGPRT_strain_alkaline_bradyzoite <- seq(103, 111)
# Pru_dHXGPRT_strain_CO2_starvation_bradyzoite_hourly <- seq(112, 122)
# Pru_dHXGPRT_strain_CO2_sodium_nitroprusside_bradyzoite <- seq(123, 131)


library(pheatmap)
pheatmap(rna_mat, cluster_cols = FALSE, show_colnames = FALSE)
pheatmap(rna_mat[, 1:4], cluster_cols = FALSE, show_colnames = FALSE)
colnames(rna_mat)[1:10]

pheatmap(microarray_mat[, TgRH_Tachy_invade_Marray], show_colnames = F, cluster_cols = F)

TgRH_Tachy_invade_Marray_data <- microarray_mat[, TgRH_Tachy_invade_Marray]
normalised_TgRH_Tachy_invade_Marray_data <- t(scale(t(TgRH_Tachy_invade_Marray_data)))
normalised_TgRH_Tachy_invade_Marray_data <- na.omit(normalised_TgRH_Tachy_invade_Marray_data)

Bradyzoite_Differentiation_72_hrs_data <- microarray_mat[, Bradyzoite_Differentiation_72_hrs]
normalised_Bradyzoite_Differentiation_72_hrs_data <- Bradyzoite_Differentiation_72_hrs_data %>%
  t() %>%
  scale() %>%
  t()

Bradyzoite_Differentiation_72_hrs_data %>%
  pheatmap(
    show_colnames = F,
    cluster_cols = F
  )

normalised_Bradyzoite_Differentiation_72_hrs_data %>%
  pheatmap(
    show_colnames = F,
    cluster_cols = F
  )

pheatmap(TgRH_Tachy_invade_Marray_data,
  show_colnames = F,
  cluster_cols = F,
  main = ""
)

pheatmap(normalised_TgRH_Tachy_invade_Marray_data,
  show_colnames = F,
  cluster_cols = F,
)


cell_cycle_data <- microarray_mat[, m_white_cell_cycle]
normalised_cell_cycle_data <- t(scale(t(cell_cycle_data)))

data_modelled <- list(
  cell_cycle_data,
  normalised_cell_cycle_data,
  protein_mat
)

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
