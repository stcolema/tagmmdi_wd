#!/usr/bin/Rscript
# Summary: Creates 3 .csv files of the data modelled and an .rds file of all the
# inputs for the MDI model in their appropriate format.
# 
# Example: Rscript tGondiiMDIInputs.R --K 125 --save_dir "./T_gondii/Output/" 
#   --data_dir "./T_gondii/Original_data/"
# 
# Author: Stephen Coleman
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

normaliseForCoexpression <- function(X) {
  na.omit(t(scale(t(X))))
}


setMyTheme()

# #' @title Prepare MS Object
# #' @description Prepares a mass spectrometry experiment (as stored in
# #' the Bioconductor package ``pRolocdata``) for modelling, extracting the
# #' numerical data, the classes and the indicator matrix for which labels are
# #' observed.
# #' @param MS_object A mass spectrometry experiment such as ``tan2009r1`` from
# #' ``pRolocdata``.
# #' @return A list of ``X``, the fracitonation data from a LOPIT experiment,
# #'  ``fixed``, the matrix indicating which labels are observed,
# #'  ``initial_labels``, the matrix of the initial labels to be input into
# #'  ``runMCMCChains`` (note that ``"unknown"`` organelles are represented
# #'  arbitrarily with a ``1`` as these will be sampled again within the wrapper
# #'  of ``callMDI``) and ``class_key`` which maps the numeric representation of
# #'  organelles back to the original naming.
# #' @export
# prepareMSObject <- function(MS_object) {
#   
#   # Extract the LOPIT data and the organelles
#   X <- Biobase::exprs(MS_object)
#   organelles <- fData(MS_object)[, "markers"]
#   
#   # Create a data frame of the classes present and their associated number;\
#   # this can be used to map the numeric representation of the classes back to
#   # an organelle
#   organelles_present <- pRoloc::getMarkerClasses(MS_object)
#   class_key <- data.frame(
#     Organelle = organelles_present,
#     Key = 1:length(organelles_present)
#   )
#   
#   # Number of components modelled
#   K <- length(organelles_present)
#   
#   # Number of views modelled
#   V <- 1
#   
#   # Number of samples modelled
#   N <- nrow(X)
#   
#   # Prepare initial labels
#   initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)
#   
#   # Fix training points, allow test points to move component
#   fix_vec <- (organelles != "unknown") * 1
#   fixed[, 1] <- fix_vec
#   
#   # Assign initial labels
#   initial_labels[, 1] <- class_key$Key[match(organelles, class_key$Organelle)]
#   
#   # Any unknown labels are given an arbitrary value which will be reset in the
#   # model call function.
#   initial_labels[is.na(initial_labels)] <- 1
#   
#   data_modelled <- list(
#     X
#   )
#   
#   # Return the prepared objects
#   list(
#     X = X,
#     fixed = fixed,
#     initial_labels = initial_labels,
#     class_key = class_key
#   )
# }


# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    
    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--data_dir"),
                          type = "character",
                          help = "Path to the directory containing the data.",
                          metavar = "character"
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
data_dir <- args$data_dir # "./T_gondii/Original_data/"
save_dir <- args$save_dir # "./T_gondii/Prepared_data/"

# The number of components modelled
K <- args$K

# Number of clusters modelled in the categorical dataset
n_clust_unsupervised <- K
if (is.null(K)) {
  n_clust_unsupervised <- 125
}

save_file <- paste0(
  save_dir,
  "TGondiiMDI",
  "_K_",
  n_clust_unsupervised,
  "_input_marker_proteins",
  ".rds"
)

# Read in the data
microarray_file <- paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
rna_seq_file <- paste0(data_dir, "ToxoDB_TgME49_Protein-coding_RNA-Seq.txt")
data(Barylyuk2020ToxoLopit)

# Use the same test indices across methods
marker.data <- pRoloc::markerMSnSet(Barylyuk2020ToxoLopit)

marker_proteins <- row.names(fData(marker.data))

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

kept_rows <- which(rna_genes %in% remaining_genes) # match(remaining_genes, rna_genes)
cleaned_rna_seq_data <- rna_seq_data[kept_rows, ]

rna_seq_numeric_columns <- seq(4, ncol(cleaned_rna_seq_data))
microarray_numeric_columns <- seq(4, ncol(cleaned_microarray_data))

rna_mat <- log(as.matrix(cleaned_rna_seq_data[, ..rna_seq_numeric_columns]) + 1)
microarray_mat <- as.matrix(cleaned_microarray_data[, ..microarray_numeric_columns])

row.names(rna_mat) <- cleaned_rna_seq_data[[1]]
row.names(microarray_mat) <- cleaned_microarray_data[[1]]

proteins_present_ind <- which(marker_proteins %in% cleaned_rna_seq_data[[1]])
proteins_present <- marker_proteins[proteins_present_ind]

protein_lst <- prepareMSObject(Barylyuk2020ToxoLopit[proteins_present, ], order_by_protein_name = FALSE)
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

rna_reducced_mat <- rna_mat[proteins_present, ]
microarray_reduced_mat <- microarray_mat[proteins_present, ]

rnaseq_macrophages_infected_by_T_gondii_inds <- seq(85, 142)
rnaseq_macrophages_infected_by_T_gondii <- rna_reducced_mat[, rnaseq_macrophages_infected_by_T_gondii_inds]
nonunique_rnaseq_macrophages_infected_by_T_gondii <- rnaseq_macrophages_infected_by_T_gondii[, seq(26, 54)]
unique_rnaseq_macrophages_infected_by_T_gondii <- rnaseq_macrophages_infected_by_T_gondii[, seq(1, 25)]

normalised_rnaseq_macrophages_infected_by_T_gondii <- normaliseForCoexpression(rnaseq_macrophages_infected_by_T_gondii)
normalised_nonuinque_rnaseq_macrophages_infected_by_T_gondii <- normaliseForCoexpression(nonunique_rnaseq_macrophages_infected_by_T_gondii)
normalised_unique_rnaseq_macrophages_infected_by_T_gondii <- normaliseForCoexpression(unique_rnaseq_macrophages_infected_by_T_gondii)


# Do not use the aynchronous reading
## m_white_cell_cycle_inds <- seq(1, 14) # 27)
m_white_cell_cycle_inds <- seq(2, 14) # 27)
m_white_cell_cycle <- microarray_reduced_mat[, m_white_cell_cycle_inds]
m_white_cell_cycle_normalised <- m_white_cell_cycle %>%
  t() %>%
  scale() %>%
  t() %>%
  na.omit()

final_protein_df <-  data.frame(protein_lst$X) %>%
  mutate(
    Protein = row.names(protein_lst$X),
    Label = factor(protein_lst$initial_labels[, 1]),
    Fixed = protein_lst$fixed[, 1]
  )

protein_mat <- protein_lst$X

shortened_col_names <- normalised_rnaseq_macrophages_infected_by_T_gondii |> 
  colnames() |>
  stringr::str_remove("Murine macrophages infected by 29 different strains of T. gondii - ") |> 
  stringr::str_remove(" \\(Tg 29strains inMouse RNA-Seq\\)") |> 
  stringr::str_remove("infected ")

colnames(normalised_rnaseq_macrophages_infected_by_T_gondii) <- shortened_col_names

write.csv(m_white_cell_cycle_normalised, paste0(save_dir, "markerProteinsCellCycleNormalised.csv"))
write.csv(normalised_rnaseq_macrophages_infected_by_T_gondii, paste0(save_dir, "markerProteinsRNASeqMacrophage.csv"))
write.csv(final_protein_df, paste0(save_dir, "markerProteinsLOPIT.csv"))

data_modelled <- list(
  protein_mat,
  m_white_cell_cycle_normalised,
  normalised_rnaseq_macrophages_infected_by_T_gondii
)

N <- nrow(m_white_cell_cycle_normalised)
V <- length(data_modelled)

initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)

fixed[, 1] <- final_protein_df$Fixed
initial_labels[, 1] <- final_protein_df$Label

# types <- c("MVN", "G", "TAGM")
types <- c("TAGPM", "GP", "G")

K <- c(
  length(pRoloc::getMarkerClasses(Barylyuk2020ToxoLopit)),
  n_clust_unsupervised,
  n_clust_unsupervised
  
)

train_inds <- which( final_protein_df$Fixed == 1)
row_order <- findOrder(protein_mat[ train_inds, ])
p_gg <- prepDataForggHeatmap(protein_mat[train_inds, ], row_order = row_order)
rna_gg <- prepDataForggHeatmap(normalised_rnaseq_macrophages_infected_by_T_gondii[train_inds, ], row_order = row_order)
microarray_gg <- prepDataForggHeatmap(m_white_cell_cycle_normalised[train_inds, ], row_order = row_order)

p_gg$Dataset <- "hyperLOPIT"
rna_gg$Dataset <- "Murine macrophages"
microarray_gg$Dataset <- "Cell-cycle"

gg_df <- rbind(p_gg, microarray_gg, rna_gg)
p_heatmap_comparison <- gg_df |> 
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() + 
  facet_wrap(~Dataset, scales = "free_x") + 
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")

cat("\n\n=== INPUT PREPARED ================================================\n")

mcmc_input <- list(
  data_modelled = data_modelled, 
  initial_labels = initial_labels,
  fixed = fixed,
  K = K,
  types = types
)

saveRDS(mcmc_input, file = save_file)

cat("\n\n=== INPUT SAVED ===================================================\n")
R <- 5000
thin <- 100
n_chains <- 3
burn <- 1250
K <- c(26, 75, 75)

# saveRDS(mdi_mcmc, file = "~/Desktop/mdi3chainsR5000.rds")

mdi_mcmc <- runMCMCChains(data_modelled, n_chains, R, thin, 
                          initial_labels = initial_labels, 
                          fixed = fixed, 
                          types = types, 
                          K = K
)

processed_mcmc <- predictFromMultipleChains(mdi_mcmc, burn, construct_psm = TRUE)

mdi_mcmc[[1]]$phis |> boxplot()

processed_mcmc$phis |> boxplot()

psm1 <- mdi_mcmc[[1]]$allocations[, , 1] |> createSimilarityMat()
psm2 <- mdi_mcmc[[1]]$allocations[, , 2] |> createSimilarityMat()
psm3 <- mdi_mcmc[[1]]$allocations[, , 3] |> createSimilarityMat()

pheatmap::pheatmap(processed_mcmc$psm[[2]])
pheatmap::pheatmap(processed_mcmc$psm[[3]])

processed_mcmc$allocations |> str()

which(colMeans(processed_mcmc$allocations[[1]] == processed_mcmc$allocations[[2]]) > 0.5)
which(colMeans(processed_mcmc$allocations[[1]] == processed_mcmc$allocations[[3]]) > 0.5)
which(colMeans(processed_mcmc$allocations[[2]] == processed_mcmc$allocations[[3]]) > 0.5)

phi_names <- paste0("Phi_", c("12", "13", "23"))

for(ii in seq(1, n_chains)) {
  
  .phi_df <- mdi_mcmc[[ii]]$phi |> 
    as.data.frame() |> 
    set_colnames(phi_names) |> 
    mutate(Iteration = seq(0, R, thin), Chain = ii) |> 
    pivot_longer(-c(Iteration, Chain), names_to = "Parameter", values_to = "Sampled_value")
  
  .mass_df <- mdi_mcmc[[ii]]$mass |> 
    as.data.frame() |> 
    set_colnames(paste0("Mass_", seq(1, V))) |> 
    mutate(Iteration = seq(0, R, thin), Chain = ii) |> 
    pivot_longer(-c(Iteration, Chain), names_to = "Parameter", values_to = "Sampled_value")
  if(ii == 1) {
    mass_df <- .mass_df
    phi_df <- .phi_df
  }  else {
    mass_df <- rbind(mass_df, .mass_df)
    phi_df <- rbind(phi_df, .phi_df)
  }
}

phi_df$Chain <- factor(phi_df$Chain)
mass_df$Chain <- factor(mass_df$Chain)

mass_df |> 
  ggplot(aes(x = Iteration, y = Sampled_value, colour = Chain)) +
  geom_line() +
  facet_wrap(~Parameter) +
  ggthemes::scale_color_colorblind()

phi_df |> 
  ggplot(aes(x = Iteration, y = Sampled_value, colour = Chain)) +
  geom_line() +
  facet_wrap(~Parameter) +
  ggthemes::scale_color_colorblind()

cbind(mdi_mcmc[[1]]$evidence,
      mdi_mcmc[[2]]$evidence,
      mdi_mcmc[[3]]$evidence
) |> 
  as.data.frame() |> 
  set_colnames(paste0("Chain_", seq(1, n_chains))) |> 
  mutate(Iteration = seq(thin, R, thin)) |> 
  pivot_longer(-Iteration, names_to = "Chain", values_to = "Evidence") |> 
  ggplot(aes(x = Iteration, y = Evidence, colour = Chain)) +
  geom_line() +
  ggthemes::scale_color_colorblind()

chain_used <- processMCMCChain(mdi_mcmc[[3]], 4000)

chain_used$allocations

which(colMeans(chain_used$allocations[, ,1] == chain_used$allocations[, ,2]) > 0.5)
which(colMeans(chain_used$allocations[, ,1] == chain_used$allocations[, , 3]) > 0.5)
which(colMeans(chain_used$allocations[, ,2] == chain_used$allocations[, , 3]) > 0.5)
