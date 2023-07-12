


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
  pattern = "(TGondii_CellCycle_seed_).*\\_K_125_R_65000.rds$", full.names = TRUE
)


plotting <- FALSE
plot_width <- 8
plot_height <- 6

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
saveRDS(predictions, file = "./MVNMMCellCyclePredictions.rds")

allocations <- predictions$allocations[[1]]

rm(mcmc_output, predictions)
psm <- mdiHelpR::createSimilarityMat(allocations)
mm_predicted <- mcclust.ext::minVI(psm, max.k = 100)
mvn_predicted_clustering <- mm_predicted$cl
# pheatmap(psm)

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

names(mvn_predicted_clustering) <- row.names(psm) <- row.names(microarray_data)

# annotatedHeatmap(psm, mvn_predicted_clustering, col_pal = simColPal(),
#                  my_breaks = defineBreaks(col_pal = simColPal(), lb = 0, ub = 1))

# Read in the output from MDI.
mdi_output <- readRDS("./MDIpredictions.rds")
fused_genes <- readRDS("./fusedGenes.rds")

# lopit_df <- read.csv("./lopitDataWithTAGMandMDIpredictions.csv")

# The point estimate classification/clustering
# point_estimate <- predictions$pred


mdi_psm <- mdiHelpR::createSimilarityMat(mdi_output$allocations[[1]])
mdi_predicted <- mcclust.ext::minVI(mdi_psm, max.k = 100)
mdi_predicted_clustering <- mdi_predicted$cl

# This data frame will hold the LOPIT data and additional data
data_df <- microarray_data %>%
  set_colnames(c("Asynchronous", paste0(0:12, "HR")))
measurements <- colnames(data_df)

data_df$MVN <- mvn_predicted_clustering
data_df$MDI <- mdi_predicted_clustering
data_df$Gene <- row.names(microarray_data)

protein_lst <- prepareMSObject(Barylyuk2020ToxoLopit)

fused_proteins <- data_df$Gene[fused_genes$`microarray-lopit`]

long_data_df <- data_df %>%
  pivot_longer(-c(
    Gene,
    MVN,
    MDI
  ), names_to = "Timepoint", values_to = "Value") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = measurements, ordered = TRUE),
    Disagree = MVN != MDI
    # MVN = factor(MVN,
    #               levels = protein_lst$class_key$Key,
    #               labels = protein_lst$class_key$Organelle
    # ),
    # MDI = factor(MDI,
    #              levels = protein_lst$class_key$Key,
    #              labels = protein_lst$class_key$Organelle
    # )
  )




long_data_df %>%
  # filter(Protein %in% fused_proteins) %>%
  pivot_longer(c(MVN, MDI), names_to = "Model", values_to = "Predicted_label") %>%
  ggplot(aes(x = Timepoint, y = Value, colour = Predicted_label, group = Gene)) +
  geom_line() +
  facet_grid(Model ~ Predicted_label)
# facet_grid(Model~Predicted_label)

long_data_df %>%
  # mutate(Disagree = MVN != MDI) %>%
  filter(Gene %in% fused_proteins, Disagree) %>%
  pivot_longer(c(MVN, MDI), names_to = "Model", values_to = "Predicted_label") %>%
  ggplot(aes(x = Timepoint, y = Value, colour = Predicted_label, group = Gene)) +
  geom_line() +
  facet_grid(Model ~ Predicted_label)

long_data_df %>%
  # mutate(Disagree = TAGM != MDI) %>%
  filter(Disagree) %>%
  pivot_longer(c(MVN, MDI), names_to = "Model", values_to = "Predicted_label") %>%
  ggplot(aes(x = Timepoint, y = Value, colour = Predicted_label, group = Gene)) +
  geom_line() +
  facet_grid(Model ~ Predicted_label)

col_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
my_breaks <- defineDataBreaks(microarray_data, col_pal)

K <- max(c(mvn_predicted_clustering, mdi_predicted_clustering))

anno_row <- data_df[, 15:16] %>%
  mutate(
    MVN = factor(MVN,
      levels = seq(K),
      labels = paste("Cluster", seq(K))
    ),
    MDI = factor(MDI,
      levels = seq(K),
      labels = paste("Cluster", seq(K))
    )
  )



ann_colours <- viridis::viridis(K) %>%
  magrittr::set_names(paste("Cluster", seq(K)))

ann_colours <- list("MDI" = ann_colours, "MVN" = ann_colours)


pheatmap::pheatmap(microarray_data,
  color = col_pal,
  breaks = my_breaks,
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = TRUE
)

mvn_order <- findOrder(psm)
mdiHelpR::compareHeatmapPSMandData(
  psm[mvn_order, mvn_order], 
  microarray_data[mvn_order, ], 
  save_name = "MVNcellCyclePSM_mvn_ordered.png", 
  main = "MVN psm and cell cycle data (ordered by MVN)"
)

mdiHelpR::compareHeatmapPSMandData(
  mdi_psm[mvn_order, mvn_order], 
  microarray_data[mvn_order, ], 
  save_name = "MDIcellCyclePSM_mvn_ordered.png", 
  main = "MDI psm and cell cycle data (ordered by MVN)"
)


mdi_order <- findOrder(mdi_psm)
mdiHelpR::compareHeatmapPSMandData(
  psm[mdi_order, mdi_order], 
  microarray_data[mdi_order, ], 
  save_name = "MVNcellCyclePSM_mdi_ordered.png", 
  main = "MVN psm and cell cycle data (ordered by MDI)"
)

mdiHelpR::compareHeatmapPSMandData(
  mdi_psm[mdi_order, mdi_order], 
  microarray_data[mdi_order, ], 
  save_name = "MDIcellCyclePSM_mdi_ordered.png", 
  main = "MDI psm and cell cycle data (ordered by MDI)"
)

data_order <- findOrder(microarray_data)
mdiHelpR::compareHeatmapPSMandData(
  psm[data_order, data_order], 
  microarray_data[data_order, ], 
  save_name = "MVNcellCyclePSM_data_ordered.png", 
  main = "MVN psm and cell cycle data (ordered by data)"
)

mdiHelpR::compareHeatmapPSMandData(
  mdi_psm[data_order, data_order], 
  microarray_data[data_order, ], 
  save_name = "MDIcellCyclePSM_data_ordered.png", 
  main = "MDI psm and cell cycle data (ordered by data)"
)
