


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
suppressMessages(library(tibble))
suppressMessages(library(ggforce))
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


plot_dir <- "./T_gondii/Output/Plots/"
dir.create(plot_dir, showWarnings = FALSE)

# Directories for input data (to model), model output and where to save results
inputdata_dir <- "./T_gondii/Prepared_data/" # args$data_dir
microarray_file <- paste0(inputdata_dir, "cellCycleNormalised.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
rna_seq_file <- paste0(inputdata_dir, "rnaSeqMacrophage.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_RNA-Seq.txt")
lopit_file <- paste0(inputdata_dir, "LOPITreduced.csv")

save_dir <- "./" # args$save_dir
data_dir <- "./T_gondii/Output/"

files <- list.files(data_dir,
  # pattern = "TGondii_TAGM_seed_*",
  pattern = "(TGondii_TAGM_seed_).*\\_R_65000.rds$", full.names = TRUE
)


mdi_output <- readRDS("./MDIpredictions.rds")
fused_genes <- readRDS("./fusedGenes.rds")

plotting <- FALSE
plot_height <- 6
plot_width <- 8

n_files <- length(files)
mcmc_output <- vector("list", n_files)
burn <- 40000

for (ii in seq(1, n_files)) {
  .f <- files[ii]
  .x <- readRDS(.f)[[1]]
  mcmc_output[[ii]] <- .mcmc <- processMCMCChain(.x,
    burn = burn,
    point_estimate_method = "median"
  )

  mcmc_output[[ii]]$Chain <- ii
  iterations <- seq(.mcmc$burn, .mcmc$R, .mcmc$thin)


  .lkl_df <- .mcmc$complete_likelihood %>%
    data.frame() %>%
    set_colnames(c("Complete_log_likelihood")) %>%
    mutate(Chain = ii, Iteration = iterations)

  if (ii == 1) {
    lkl_df <- .lkl_df
  } else {
    lkl_df <- rbind(lkl_df, .lkl_df)
  }

  #
  # .phi_df <- .mcmc$phis %>%
  #   data.frame() %>%
  #   set_colnames(c("Phi_12", "Phi_13", "Phi_23")) %>%
  #   mutate(Chain = ii, Iteration = iterations)
  #
  # if(ii == 1) {
  #   phi_df <- .phi_df
  # } else {
  #   phi_df <- rbind(phi_df, .phi_df)
  # }
}


lkl_df$Chain <- factor(lkl_df$Chain)

# Find a point estimate from the chains
predictions <- predictFromMultipleChains(mcmc_output, burn = burn, chains_already_processed = TRUE)
allocations <- predictions$allocations

# Read in the data used in the original analysis
data(Barylyuk2020ToxoLopit)

microarray_data <- read.csv(microarray_file,
  row.names = 1,
  header = T
)

rna_seq_data <- read.csv(rna_seq_file,
  row.names = 1,
  header = T
)

lopit_data <- read.csv(lopit_file,
  row.names = 1,
  header = T
)

# The point estimate classification/clustering
point_estimate <- predictions$pred

protein_lst <- prepareMSObject(Barylyuk2020ToxoLopit)
measurements <- colnames(Barylyuk2020ToxoLopit)
fused_proteins <- data_df$Protein[fused_genes$`microarray-lopit`]

# This data frame will hold the LOPIT data and additional data
data_df <- lopit_data %>%
  mutate(
    TAGM = point_estimate[[1]],
    MDI = mdi_output$pred[[3]],
    TAGM_prob = predictions$prob[[1]],
    MDI_prob = mdi_output$prob[[3]],
    Fused = Protein %in% fused_proteins,
    Disagree = TAGM != MDI
  )

disagreed_proteins <- data_df$Protein[data_df$Disagree]

# data_df$Fused[fused_genes$`microarray-lopit`] <- TRUE

# data_df$TAGM <- point_estimate[[1]]
# data_df$MDI <- mdi_output$pred[[3]]

# data_df$Fused <- FALSE
# data_df$Fused[fused_genes$`microarray-lopit`] <- TRUE



long_data_df <- data_df %>%
  pivot_longer(-c(
    Protein,
    Label,
    Fixed,
    TAGM,
    MDI,
    TAGM_prob,
    MDI_prob,
    Fused,
    Disagree
  ), names_to = "Fraction", values_to = "Value") %>%
  mutate(
    Fraction = factor(Fraction, levels = measurements, ordered = TRUE),
    Label = factor(Label,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    ),
    TAGM = factor(TAGM,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    ),
    MDI = factor(MDI,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    )
  )

if (plotting) {
  long_data_df %>%
    # filter(Protein %in% fused_proteins) %>%
    pivot_longer(c(TAGM, MDI), names_to = "Model", values_to = "Predicted_label") %>%
    ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
    geom_line() +
    facet_grid_paginate(Model ~ Predicted_label, ncol = 5, nrow = 2, page = 1)
  # facet_grid(Model~Predicted_label)

  long_data_df %>%
    # mutate(Disagree = TAGM != MDI) %>%
    filter(Fused, Disagree) %>%
    pivot_longer(c(TAGM, MDI), names_to = "Model", values_to = "Predicted_label") %>%
    ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
    geom_line() +
    facet_grid(Model ~ Predicted_label)

  long_data_df %>%
    # mutate(Disagree = TAGM != MDI) %>%
    filter(Disagree) %>%
    pivot_longer(c(TAGM, MDI), names_to = "Model", values_to = "Predicted_label") %>%
    ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
    geom_line() +
    facet_grid(Model ~ Predicted_label)
}


long_data_df %>%
  # mutate(Disagree = TAGM != MDI) %>%
  filter(Fused, Disagree) %>%
  pivot_longer(c(TAGM, MDI), names_to = "Model", values_to = "Predicted_label") %>%
  ggplot(aes(x = Fraction, y = Value, colour = Predicted_label, group = Protein)) +
  geom_line() +
  facet_grid(Model ~ Predicted_label)

tagm_classification_posterior <- predictions$allocation_probability[[1]] %>%
  set_rownames(data_df$Protein) %>%
  set_colnames(protein_lst$class_key$Organelle)
mdi_classification_posterior <- mdi_output$allocation_probability[[3]] %>%
  set_rownames(data_df$Protein) %>%
  set_colnames(protein_lst$class_key$Organelle)

tagm_prob_df <- tagm_classification_posterior %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "Organelle", values_to = "TAGM")


mdi_prob_df <- mdi_classification_posterior %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "Organelle", values_to = "MDI")

prob_df <- full_join(mdi_prob_df, tagm_prob_df, by = c("Protein", "Organelle"))

prob_df %>%
  filter(
    # Protein %in% fused_proteins,
    Protein %in% disagreed_proteins[sample(1:496, size = 10)]
  ) %>% # ,
  # Organelle %in% c("cytosol", "Golgi", "nucleolus", "nucleus - chromatin", "nucleus - non-chromatin")) %>%
  pivot_longer(c(MDI, TAGM), names_to = "Model", values_to = "Probability") %>%
  ggplot(aes(x = Organelle, y = Probability, group = Protein)) +
  geom_line(aes(colour = Protein)) +
  facet_wrap(~Model, nrow = 2) +
  theme(legend.position = "none")


disagreeing_df <- data_df
# %>%
# filter(Protein %in% disagreed_proteins)

lopit_exprs <- disagreeing_df[, 1:30] %>% as.matrix()

# 19S proteasome
# apical 2
# PM - peripheral 1
# rhoptries 2
# tubulin cytoskeleton
empty_organelles <- c(1, 6, 25, 26)

anno_row <- disagreeing_df[, 34:35] %>%
  mutate(
    TAGM = factor(TAGM,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    ),
    MDI = factor(MDI,
      levels = protein_lst$class_key$Key,
      labels = protein_lst$class_key$Organelle
    )
  )

col_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
my_breaks <- defineDataBreaks(lopit_exprs, col_pal)

K <- length(protein_lst$class_key$Organelle)
ann_colours <- viridis::viridis(K) %>%
  magrittr::set_names(protein_lst$class_key$Organelle)

ann_colours <- list("MDI" = ann_colours, "TAGM" = ann_colours)

pheatmap::pheatmap(lopit_exprs,
  color = col_pal,
  breaks = my_breaks,
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = FALSE
)

write.csv(data_df, "./lopitDataWithTAGMandMDIpredictions.csv")
