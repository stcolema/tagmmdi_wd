#!/usr/bin/Rscript
#
# Example call:
#  Rscript processTGondiiCCoutput.R --data_dir ~/rds/hpc-work/tagmmdi/T_gondii/Data/ --save_dir ~/rds/hpc-work/tagmmdi/T_gondii/ConsensusClustering/ --model_output_dir ~/rds/hpc-work/tagmmdi/T_gondii/ConsensusClustering/ --R 5000 --K 125


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
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))


# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--data_dir"),
      type = "character",
      help = "Directory here the data that were modelled are stored.",
      metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--save_dir"),
      type = "character",
      help = "Directory to save the outputs of this file to.",
      metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--model_output_dir"),
      type = "character",
      help = "Directory where the model outputs are saved.",
      metavar = "character"
    ),

    # Number of MCMC iterations
    optparse::make_option(c("-R", "--R"),
      type = "numeric",
      default = 5000,
      help = "Number of iterations to run in each MCMC chain [default= %default]",
      metavar = "numeric"
    ),
    optparse::make_option(c("-b", "--burn"),
      type = "numeric",
      default = NULL,
      help = paste(
        "Number of iterations to burn of the warm-up period in each MCMC chain",
        "[default= 0.2 * R]"
      ),
      metavar = "numeric"
    ),
    optparse::make_option(c("-K", "--K"),
      type = "numeric",
      default = 125,
      help = "Number of components modelled in the unsupervised views.",
      metavar = "numeric"
    )
  )

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

cat("\n=== Running ``processTGondiiCCoutput.R`` ===============================")
setMyTheme()

args <- input_arguments()

# Directories for input and output respectively
inputdata_dir <- "./T_gondii/Prepared_data/" #
save_dir <- "./T_gondii/Analysis/" #
model_output_dir <- "./T_gondii/ConsensusClustering/" #

# === READ IN MODEL OUTPUT =====================================================
cat("\n=== Reading in files ===================================================")

mdi_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_125_full.rds"
mix_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_300_cell_cycle_mix.rds"

D <- 15000
W <- 150
mdi_mod <- readRDS(mdi_file)
mix_mod <- readRDS(mix_file)
V <- 2

pred_cl <- mdi_mod$pred
prob_cl <- mdi_mod$prob
fused_genes_1 <- which(colMeans(mdi_mod$allocations[[1]] == mdi_mod$allocations[[2]]) > 0.5)

plotting <- TRUE
plot_height <- 6
plot_width <- 8

# === Input data ===============================================================

microarray_file <- paste0(inputdata_dir, "cellCycleNormalised.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_DNA_microarray.txt")
rna_seq_file <- paste0(inputdata_dir, "rnaSeqMacrophage.csv") # paste0(data_dir, "ToxoDB_TgME49_Protein-coding_RNA-Seq.txt")
lopit_file <- paste0(inputdata_dir, "LOPITreduced.csv")
data(Barylyuk2020ToxoLopit)

microarray_data <- read.csv(microarray_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 212)
)

rna_seq_data <- read.csv(rna_seq_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 255)
)

lopit_data <- read.csv(lopit_file,
  row.names = 1,
  # na.strings = "N/A",
  # strip.white = T,
  header = T
  # select = seq(1, 255)
)

data_modelled <- readRDS(paste0(inputdata_dir, "TGondiiMDI_K_125_input.rds"))
datasets <- c("LOPIT", "Cell_cycle", "RNA-seq")

# === LOPIT ===================================================================

# Compare MDI and TAGM allocations

# Proteins in final analysis
proteins_modelled <- row.names(data_modelled$data_modelled[[1]])

# Relevant columns from the pRolocdata object
cols_used <- c("Description", "markers", "tagm.mcmc.allocation", "tagm.mcmc.probability")

# Relevant data for comparing allocations
tagm_comparison <- fData(Barylyuk2020ToxoLopit)[proteins_modelled, cols_used]

# Map between symbolic labels and localisation name
label_to_organelle <- data.frame(
  "Organelle" = levels(tagm_comparison$markers)[-27],
  "Label" = seq(1, 26)
)


marker_labels <- label_to_organelle$Organelle[data_modelled$initial_labels[, 1]]
marker_labels[data_modelled$fixed[, 1] == 0] <- NA

# The MDI predictions and probability of allocation
mdi_predictions <- label_to_organelle$Organelle[pred_cl[[1]]]
mdi_probabilities <- prob_cl[[1]]

tagm_comparison <- tagm_comparison[proteins_modelled, ]
tagm_comparison$mdi.mcmc.allocation <- mdi_predictions
tagm_comparison$mdi.mcmc.probability <- mdi_probabilities

# Check which types of disagreements emerge
tagm_comparison[, c(3, 5)] |> unique()
tagm_comparison[fused_genes_1, c(3, 5)] |> unique()
tagm_comparison[
  fused_genes_1[which(tagm_comparison[fused_genes_1, ]$mdi.mcmc.allocation != "Golgi")],
  c(3, 5)
]

lopit_disagreement <- tagm_comparison$mdi.mcmc.allocation != tagm_comparison$tagm.mcmc.allocation
tagm_uncertain <- tagm_comparison$tagm.mcmc.probability < 0.7
mdi_certain <- tagm_comparison$mdi.mcmc.probability > 0.7
non_golgi <- tagm_comparison$mdi.mcmc.allocation != "Golgi"
disagreeing_uncertain_inds <- which(lopit_disagreement & tagm_uncertain & non_golgi & mdi_certain)

diagreement_uncertain_table <- tagm_comparison[disagreeing_uncertain_inds, ]

diagreement_uncertain_table$tagm.mcmc.probability <- diagreement_uncertain_table$tagm.mcmc.probability |>
  round(digits = 3)
diagreement_uncertain_table$mdi.mcmc.probability <- diagreement_uncertain_table$mdi.mcmc.probability |>
  round(digits = 3)

write.csv(diagreement_uncertain_table, "~/Desktop/DisagreeingLocalisationsNonGolgiUncertain.csv")

pm_localisation_options <- c("PM - integral", "PM - peripheral 1", "PM - peripheral 2")

golgi_changes_to_pm <- tagm_comparison |>
  dplyr::filter(tagm.mcmc.allocation == "Golgi", mdi.mcmc.allocation %in% pm_localisation_options)

pm_marker_proteins <- tagm_comparison |>
  dplyr::filter(markers %in% pm_localisation_options)

golgi_marker_proteins <- tagm_comparison |>
  dplyr::filter(markers == "Golgi")

golgi_changes_and_pm_markers <- rbind(pm_marker_proteins, golgi_changes_to_pm)
golgi_changes_and_markers <- rbind(pm_marker_proteins, golgi_changes_to_pm, golgi_marker_proteins)

protein_biochemical_properties <- read.csv("./T_gondii/TGondiiGolgiProteinsAttributes.csv", row.names = 1)

rel_protein_biochemical_properties <- protein_biochemical_properties[row.names(golgi_changes_and_markers), ]

my_df <- cbind(golgi_changes_and_markers, rel_protein_biochemical_properties[, -c(2, 4, 5, 6)]) |>
  pivot_longer(c(mdi.mcmc.allocation, tagm.mcmc.allocation), names_to = "Model", values_to = "Allocation")

my_df$pI <- as.numeric(my_df$pI)
my_df$Fixed <- my_df$markers != "unknown"
my_df |>
  ggplot(aes(x = Allocation, y = pI, group = Allocation)) +
  facet_grid(Model ~ Fixed) +
  geom_boxplot()

p_dNdS <- my_df |>
  ggplot(aes(x = Allocation, y = dNdS, fill = Allocation)) +
  facet_grid(~Model) +
  geom_boxplot() +
  scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])


my_df$Num.TMDs.TMHMM |> table()
my_df$TMD.within.first.60.AA <- my_df$TMD.within.first.60.AA != "FALSE"

my_df |>
  filter(!Num.TMDs.TMHMM %in% c("Golgi", "PM - integral", "PM - peripheral 1", "PM - peripheral 2", "unknown")) |>
  mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
  ggplot(aes(x = Allocation, y = Num.TMDs.TMHMM, group = Allocation)) +
  facet_grid(Model ~ Fixed) +
  geom_boxplot()

my_df |>
  filter(!Num.TMDs.TMHMM %in% c("Golgi", "PM - integral", "PM - peripheral 1", "PM - peripheral 2", "unknown")) |>
  mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
  ggplot(aes(x = Allocation, y = Num.TMDs.TMHMM, group = Allocation)) +
  facet_grid(~Model) +
  geom_boxplot()

p_conservation_score <- my_df |>
  ggplot(aes(x = Allocation, y = Conservation.score, group = Allocation, fill = Allocation)) +
  facet_grid(~Model) +
  geom_boxplot() +
  scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])

p_num_tmds_TMHMM <- my_df |>
  filter(!Num.TMDs.TMHMM %in% c("Golgi", "PM - integral", "PM - peripheral 1", "PM - peripheral 2", "unknown")) |>
  mutate(Num.TMDs.TMHMM = as.numeric(Num.TMDs.TMHMM)) |>
  ggplot(aes(x = Allocation, y = Num.TMDs.TMHMM, group = Allocation, fill = Allocation)) +
  facet_grid(~Model) +
  geom_boxplot() +
  scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])

p_tmd_first_60_aa <- my_df |>
  ggplot(aes(y = TMD.within.first.60.AA, group = Allocation, fill = Allocation)) +
  facet_grid(~Model) +
  stat_count() +
  scale_fill_manual(values = ggthemes::colorblind_pal()(5)[-1])

library(patchwork)
p_patch <- p_dNdS / p_conservation_score / p_num_tmds_TMHMM / p_tmd_first_60_aa + plot_layout(guides = "collect")

ggsave("~/Desktop/golgi_new_allocations_biochemical__properties.png", plot = p_patch, height = 8, width = 7)

my_df |>
  ggplot(aes(x = Allocation, y = dNdS, group = Allocation)) +
  facet_wrap(~Model) +
  geom_boxplot()

my_df |>
  ggplot(aes(x = tagm.mcmc.allocation, y = pI, group = tagm.mcmc.allocation)) +
  geom_boxplot()
