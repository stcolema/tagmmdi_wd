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

mdi_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_125_V_2_all_other_objects.rds"
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

# === Hyper parameters =========================================================

amplitude_inds <- seq(1, mix_mod$K)
length_scale_inds <- seq(mix_mod$K + 1, 2 * mix_mod$K)
noise_inds <- seq(2 * mix_mod$K + 1, 3 * mix_mod$K)

mix_amp_acc_df <- data.frame("Parameter" = "Amplitude", "Acceptance_rate" = c(mix_mod$acceptance_count[[1]][, amplitude_inds]))
mix_length_acc_df <- data.frame("Parameter" = "Length-scale", "Acceptance_rate" = c(mix_mod$acceptance_count[[1]][, length_scale_inds]))
mix_noise_acc_df <- data.frame("Parameter" = "Noise", "Acceptance_rate" = c(mix_mod$acceptance_count[[1]][, noise_inds]))

amplitude_inds <- seq(1, mdi_mod$K[2])
length_scale_inds <- seq(mdi_mod$K[2] + 1, 2 * mdi_mod$K[2])
noise_inds <- seq(2 * mdi_mod$K[2] + 1, 3 * mdi_mod$K[2])

mdi_amp_acc_df <- data.frame("Parameter" = "Amplitude", "Acceptance_rate" = c(mdi_mod$acceptance_count[[2]][, amplitude_inds]))
mdi_length_acc_df <- data.frame("Parameter" = "Length-scale", "Acceptance_rate" = c(mdi_mod$acceptance_count[[2]][, length_scale_inds]))
mdi_noise_acc_df <- data.frame("Parameter" = "Noise", "Acceptance_rate" = c(mdi_mod$acceptance_count[[2]][, noise_inds]))


mix_acc_df <- rbind(mix_amp_acc_df, mix_length_acc_df, mix_noise_acc_df)
mdi_acc_df <- rbind(mdi_amp_acc_df, mdi_length_acc_df, mdi_noise_acc_df)

mdi_acc_df |>
  ggplot(aes(x = Parameter, y = Acceptance_rate)) +
  geom_boxplot(fill = "gold") +
  labs(y = "Acceptance rate across chains")

mix_acc_df |>
  ggplot(aes(x = Parameter, y = Acceptance_rate)) +
  geom_boxplot(fill = "gold") +
  labs(y = "Acceptance rate across chains")

mix_mod$hypers[[1]]$length |>
  log() |>
  boxplot()
mix_mod$hypers[[1]]$noise |>
  log() |>
  boxplot()
mix_mod$hypers[[1]]$amplitude |>
  log() |>
  boxplot()

filled_inds <- seq(7, 24)

mdi_mod$hypers[[2]]$length[, filled_inds] |>
  log() |>
  boxplot()
mdi_mod$hypers[[2]]$noise[, filled_inds] |>
  log() |>
  boxplot()
mdi_mod$hypers[[2]]$amplitude[, filled_inds] |>
  log() |>
  boxplot()

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

write.csv(golgi_changes_and_pm_markers, "~/Desktop/DisagreeingChangesFromGolgi.csv")


disagreeing_inds <- which(lopit_disagreement & non_golgi)
diagreement_table <- tagm_comparison[disagreeing_inds, ]

diagreement_table$tagm.mcmc.probability <- diagreement_table$tagm.mcmc.probability |>
  round(digits = 3)
diagreement_table$mdi.mcmc.probability <- diagreement_table$mdi.mcmc.probability |>
  round(digits = 3)

write.csv(diagreement_table, "~/Desktop/DisagreeingLocalisationsNonGolgi.csv")

# Compare predictions in the meassurement data
lopit_data$MDI_prediction <- mdi_predictions
lopit_data$TAGM_prediction <- tagm_comparison$tagm.mcmc.allocation

lopit_data$MDI_probability <- mdi_probabilities
lopit_data$TAGM_probability <- tagm_comparison$tagm.mcmc.probability

# Plot MDI and TAGM predicted clusters' profiles
p_mdi <- lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, MDI_probability, TAGM_prediction, TAGM_probability),
    values_to = "Measurement",
    names_to = "Fraction"
  ) |>
  mutate(Fraction = factor(Fraction), Localisation_known = factor(Fixed, levels = c(1, 0), labels = c("True", "False"))) |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, color = Localisation_known)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~MDI_prediction, ncol = 5) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))

p_tagm <- lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, MDI_probability, TAGM_prediction, TAGM_probability),
    values_to = "Measurement",
    names_to = "Fraction"
  ) |>
  mutate(Fraction = factor(Fraction), Localisation_known = factor(Fixed, levels = c(1, 0), labels = c("True", "False"))) |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, color = Localisation_known)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~TAGM_prediction) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))


p_mdi / p_tagm


golgi_inds <- which(mdi_predictions == "Golgi")
pm_peripheral_2_inds <- which(mdi_predictions == "PM - peripheral 2")
nuc_chromatin_inds <- which(mdi_predictions == "nucleus - chromatin")
pheatmap(mdi_mod$cms[[1]][nuc_chromatin_inds, nuc_chromatin_inds],
  color = simColPal(), breaks = defineBreaks(simColPal(), lb = 0),
  show_rownames = FALSE, show_colnames = FALSE
)

annotatedHeatmap(mdi_mod$cms[[1]][c(pm_peripheral_2_inds, golgi_inds), c(pm_peripheral_2_inds, golgi_inds)],
  mdi_predictions[c(pm_peripheral_2_inds, golgi_inds)],
  col_pal = simColPal(),
  my_breaks = defineBreaks(simColPal(), lb = 0),
  show_rownames = FALSE, show_colnames = FALSE
)


# Create the annotation data.frame for the rows
anno_row <- data.frame(Localisation = factor(mdi_predictions[c(pm_peripheral_2_inds, golgi_inds)]))
row.names(anno_row) <- rownames(mdi_mod$cms[[1]][c(pm_peripheral_2_inds, golgi_inds), ])

ann_colors <- list(Localisation = ggthemes::colorblind_pal()(2))
names(ann_colors$Localisation) <- unique(mdi_predictions[c(pm_peripheral_2_inds, golgi_inds)])

pheatmap::pheatmap(mdi_mod$cms[[1]][c(pm_peripheral_2_inds, golgi_inds), c(pm_peripheral_2_inds, golgi_inds)],
  color = simColPal(),
  annotation_row = anno_row,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  filename = paste0(save_dir, "MDI_Goldi_PM2_cm.png")
)

ggsave(paste0(save_dir, "MDI_LOPIT_localisations_fixed.png"), plot = p_mdi, height = 10, width = 18)

lopit_data |>
  mutate(MDI_prob = mdi_probabilities) |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction, MDI_prob), values_to = "Measurement", names_to = "Fraction") |>
  mutate(Fraction = factor(Fraction), Localisation_known = factor(Fixed, levels = c(1, 0), labels = c("True", "False"))) |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, color = Localisation_known)) +
  geom_line(aes(alpha = MDI_prob)) +
  facet_wrap(~MDI_prediction) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

lopit_plot_data <- lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |>
  mutate(Fraction = factor(Fraction), Localisation_known = factor(Fixed, levels = c(1, 0), labels = c("True", "False"))) |>
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |>
  mutate(Model = case_when(
    Model == "MDI_prediction" ~ "MDI",
    Model == "TAGM_prediction" ~ "TAGM"
  ))

# These organelles are a lot noisier in MDI vs TAGM
misbehaving_organelles <- c(
  "Golgi",
  "nucleolus"
)

p1 <- lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, color = Localisation_known)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Model ~ Organelle, nrow = 2, ncol = 13, page = 1) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# ggsave("~/Desktop/mdi_tagm_comp_for_OC.png", height = 8, width = 32)

p2 <- lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, color = Localisation_known)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Model ~ Organelle, nrow = 2, ncol = 13, page = 2) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p1 / p2 +
  plot_layout(guides = "collect")

ggsave("T_gondii/Analysis/comparison_tagm_mdi_localisation.png", height = 16, width = 24)

plotting <- TRUE
n_rows <- 10
lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = n_rows, page = 1)

if (plotting) {
  ggsave(paste0("./T_gondii/Analysis/MDI_TAGM_comp_first_", n_rows, "_organelles.png"),
    height = 12,
    width = 8
  )
}

lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = n_rows, page = 2)

if (plotting) {
  ggsave(paste0("./T_gondii/Analysis/MDI_TAGM_comp_second_", n_rows, "_organelles.png"),
    height = 12,
    width = 8
  )
}

# MDI puts a lot into GOLGI
lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = n_rows, page = 3)

if (plotting) {
  ggsave(paste0("./T_gondii/Analysis/MDI_TAGM_comp_third_", n_rows, "_organelles.png"),
    height = 12,
    width = 8
  )
}

# TAGM puts a lot more into the nulceus related organelles
lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = n_rows, page = 4)

lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = n_rows, page = 5)

lopit_plot_data |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line(alpha = 0.3) +
  ggforce::facet_grid_paginate(Organelle ~ Model, ncol = 2, nrow = 5, page = 6)

lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |>
  mutate(Fraction = factor(Fraction)) |>
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein)) +
  geom_line() +
  facet_grid(Model ~ Organelle)

lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |>
  mutate(Fraction = factor(Fraction)) |>
  pivot_longer(c(TAGM_prediction, MDI_prediction), names_to = "Model", values_to = "Organelle") |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, colour = Model)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~Organelle) +
  ggthemes::scale_color_colorblind()

tagm_comparison[which(tagm_comparison$tagm.mcmc.allocation != tagm_comparison$mdi.mcmc.allocation), c(1, 3, 5)]

tagm_comparison[fused_genes_1, c(2, 3, 5)]

p_tagm_mdi_contrast <- lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |>
  mutate(Fraction = factor(Fraction)) |>
  pivot_longer(c(TAGM_prediction, MDI_prediction), names_to = "Model", values_to = "Organelle") |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, colour = Model)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~Organelle) +
  ggthemes::scale_color_colorblind()

ggsave("./T_gondii/Analysis/tagm_mdi_contrast.png", p_tagm_mdi_contrast, height = 12, width = 10)

p_mdi_tagm_contrast <- lopit_data |>
  pivot_longer(-c(Label, Fixed, Protein, MDI_prediction, TAGM_prediction), values_to = "Measurement", names_to = "Fraction") |>
  mutate(Fraction = factor(Fraction)) |>
  pivot_longer(c(MDI_prediction, TAGM_prediction), names_to = "Model", values_to = "Organelle") |>
  ggplot(aes(x = Fraction, y = Measurement, group = Protein, colour = Model)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~Organelle) +
  ggthemes::scale_color_colorblind()

ggsave("./T_gondii/Analysis/mdi_tagm_contrast.png", p_mdi_tagm_contrast, height = 12, width = 10)

microarray_data |>
  mutate(
    Predicted_label = factor(pred_cl[[2]]),
    Gene = row.names(microarray_data)
  ) |>
  pivot_longer(-c(Predicted_label, Gene), names_to = "Time", values_to = "Expression") |>
  mutate(Time = factor(Time, labels = seq(0, 12))) |>
  # filter(Predicted_label %in% c(1:6, 9)) |>
  ggplot(aes(x = Time, y = Expression, group = Gene)) +
  geom_line(alpha = 0.3) +
  # ggthemes::scale_color_colorblind() +
  facet_wrap(~Predicted_label)

ggsave("./T_gondii/Analysis/mdi_cell_cycle_cluster_profiles.png", height = 12, width = 16)


annotatedHeatmap(microarray_data, pred_cl[[2]],
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Cell-cycle data annotated by predicted cluster",
  filename = "./T_gondii/Analysis/cell_cycle_predicted_clustering.png"
)


annotatedHeatmap(microarray_data[fused_genes_1, ], pred_cl[[2]][fused_genes_1],
  show_colnames = TRUE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  main = "Fused genes in cell-cycle data annotated by predicted cluster",
  filename = "./T_gondii/Analysis/cell_cycle_fused_genes_predicted_clustering.png"
)


annotation_row <- data.frame(
  "Cell-cycle.cluster" = factor(c(pred_cl[[2]])[fused_genes_1]),
  "Predicted.localisation" = mdi_predictions[fused_genes_1]
)
row.names(annotation_row) <- row.names(microarray_data)[fused_genes_1]

K_1 <- mdi_predictions[fused_genes_1] |>
  unique() |>
  length()
K_2 <- pred_cl[[2]][fused_genes_1] |>
  unique() |>
  length()

annotation_colors <- list(
  "Cell-cycle.cluster" = RColorBrewer::brewer.pal(K_2, "Set3"),
  "Predicted.localisation" = RColorBrewer::brewer.pal(K_1, "Set3")
)

names(annotation_colors$`Cell-cycle.cluster`) <- sort(unique(pred_cl[[2]][fused_genes_1]))
names(annotation_colors$`Predicted.localisation`) <- sort(unique(mdi_predictions[fused_genes_1]))

pheatmap::pheatmap(microarray_data[fused_genes_1, ],
  color = dataColPal(),
  breaks = defineDataBreaks(microarray_data, dataColPal()),
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  main = "Fused  genes in cell-cycle data\nannotated by predicted cluster and localisation"
  # filename = "~/Desktop/mdiCellCycleFusedGenes.png"
)

fused_microarray_data <- microarray_data[fused_genes_1, ]
fused_lopit_data <- data_modelled$data_modelled[[1]][fused_genes_1, ]
row_order <- order(mdi_predictions[fused_genes_1])

ph_lopit <- pheatmap(fused_lopit_data[row_order, ],
  color = dataColPal(),
  breaks = defineDataBreaks(fused_lopit_data, dataColPal()),
  annotation_row = annotation_row[row_order, ],
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE
  # main = "Fused  genes in cell-cycle data\nannotated by predicted cluster and localisation"
  # filename = "~/Desktop/mdiCellCycleFusedGenes.png"
)


ph_cell_cycle <- pheatmap(fused_microarray_data[row_order, ],
  color = dataColPal(),
  breaks = defineDataBreaks(fused_lopit_data, dataColPal()),
  annotation_row = annotation_row[row_order, ],
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE
  # main = "Fused  genes in cell-cycle data\nannotated by predicted cluster and localisation"
  # filename = "~/Desktop/mdiCellCycleFusedGenes.png"
)

combinePheatmaps(list(ph_lopit$gtable, ph_cell_cycle$gtable), save_name = "~/Desktop/CombinedData.png", main = NULL)

RColorBrewer::brewer.pal(10, "Paired")

# Create the annotation colours
ann_colours <- list(Cluster = viridis::viridis(K) %>%
  magrittr::set_names(paste("Cluster", sort(unique(cluster_IDs)))))

marker_genes <- which(tagm_comparison$markers != "unknown")
annotatedHeatmap(microarray_data[marker_genes, ],
  pred_cl[[1]][[2]][marker_genes],
  show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE,
  main = "Marker genes in cell cycle data"
)

# Create the annotation data.frame for the rows
anno_row <- data.frame(
  Cluster = factor(paste("Cluster", pred_cl[[1]][[2]][marker_genes])),
  Organelle = tagm_comparison$mdi.mcmc.allocation[marker_genes]
) %>%
  magrittr::set_rownames(rownames(microarray_data)[marker_genes])

# The number of cololurs to use
# K <- 26 # length(unique(cluster_IDs))

# Create the annotation colours
ann_colours <- list(Cluster = viridis::viridis(19), Organelle = viridis::viridis(26))
names(ann_colours[[1]]) <- paste("Cluster", unique(pred_cl[[1]][[2]][marker_genes]))
names(ann_colours[[2]]) <- unique(tagm_comparison$mdi.mcmc.allocation[marker_genes])
pheatmap::pheatmap(microarray_data[marker_genes, ],
  color = dataColPal(),
  breaks = defineDataBreaks(microarray_data[marker_genes, ], dataColPal(), mid_point = 0),
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster",
  filename = "./marker_genes_cell_cycle_data_comp.png"
)

organelles_of_interest <- c(
  "micronemes",
  "apical 1",
  "apical 2",
  "rhoptries 1",
  "rhoptries 2",
  "dense granules"
)
proteins_in_organelles_of_interest <- which(tagm_comparison$mdi.mcmc.allocation %in% organelles_of_interest)
cell_cycle_clusters_plotted <- which(table(cc_out$predicted_partitions[[1]][proteins_in_organelles_of_interest]) > 15)
genes_of_interest <- which(cc_out$predicted_partitions[[1]] %in% cell_cycle_clusters_plotted)

genes_plotted <- proteins_in_organelles_of_interest[which(proteins_in_organelles_of_interest %in% genes_of_interest)]

table(tagm_comparison$mdi.mcmc.allocation[proteins_in_organelles_of_interest])


pheatmap::pheatmap(microarray_data[proteins_in_organelles_of_interest, ],
  color = dataColPal(),
  breaks = defineDataBreaks(microarray_data[proteins_in_organelles_of_interest, ], dataColPal(), mid_point = 0),
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster",
  filename = "./markger_genes_cell_cycle_data_comp.png"
)


annotatedHeatmap(microarray_data[genes_plotted, ],
  tagm_comparison$mdi.mcmc.allocation[genes_plotted],
  show_colnames = FALSE, show_rownames = FALSE
)



# Create the annotation data.frame for the rows
anno_row <- data.frame(
  Cluster = factor(paste("Cluster", cc_out$predicted_partitions[[1]][genes_plotted])),
  Organelle = tagm_comparison$mdi.mcmc.allocation[genes_plotted]
) %>%
  magrittr::set_rownames(rownames(microarray_data)[genes_plotted])

# The number of cololurs to use
# K <- 26 # length(unique(cluster_IDs))

# Create the annotation colours
ann_colours <- list(Cluster = ggthemes::colorblind_pal()(7), Organelle = ggthemes::colorblind_pal()(6)) # viridis::viridis(6))
names(ann_colours[[1]]) <- paste("Cluster", unique(cc_out$predicted_partitions[[1]][genes_plotted]))
names(ann_colours[[2]]) <- unique(tagm_comparison$mdi.mcmc.allocation[genes_plotted])
pheatmap::pheatmap(microarray_data[genes_plotted, ],
  color = dataColPal(),
  breaks = defineDataBreaks(microarray_data[genes_plotted, ], dataColPal(), mid_point = 0),
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster"
  # filename = "./markger_genes_cell_cycle_data_comp.png"
)

tagm_comparison$Description[cc_out$predicted_partitions[[1]] == 3 & (tagm_comparison$Description != "hypothetical protein")]
hypo_prots <- which(tagm_comparison$Description == "hypothetical protein")

hypo_prots_in_genes_plotted <- genes_plotted[genes_plotted %in% hypo_prots]

ann_colours <- list(Cluster = ggthemes::colorblind_pal()(7), Organelle = ggthemes::colorblind_pal()(6)) # viridis::viridis(6))
names(ann_colours[[1]]) <- paste("Cluster", unique(cc_out$predicted_partitions[[1]][hypo_prots_in_genes_plotted]))
names(ann_colours[[2]]) <- unique(tagm_comparison$mdi.mcmc.allocation[hypo_prots_in_genes_plotted])
pheatmap::pheatmap(microarray_data[hypo_prots_in_genes_plotted, ],
  color = dataColPal(),
  breaks = defineDataBreaks(microarray_data[hypo_prots_in_genes_plotted, ], dataColPal(), mid_point = 0),
  annotation_row = anno_row,
  annotation_colors = ann_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Genes of marker proteins in cell cycle data\nannotated by organelle and inferred cluster"
  # filename = "./markger_genes_cell_cycle_data_comp.png"
)

row_order <- findOrder(microarray_data)
micro_df <- prepDataForggHeatmap(as.matrix(microarray_data), row_order = row_order, col_order = seq(1, ncol(microarray_data)))
lopit_df <- prepDataForggHeatmap(scale(exprs(Barylyuk2020ToxoLopit)),
  row_order = row_order,
  col_order = seq(1, ncol(exprs(Barylyuk2020ToxoLopit)))
)
rnaseq_df <- prepDataForggHeatmap(as.matrix(rna_seq_data), row_order = row_order)

lopit_df$Dataset <- "LOPIT"
micro_df$Dataset <- "Cell cycle"
rnaseq_df$Dataset <- "RNA seq"
heat_df <- rbind(micro_df, lopit_df, rnaseq_df)

heat_df |>
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() +
  facet_wrap(~Dataset, scales = "free") +
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")

row_order <- findOrder(microarray_data[marker_genes, ])
micro_df <- prepDataForggHeatmap(as.matrix(microarray_data[marker_genes, ]), row_order = row_order, col_order = seq(1, ncol(microarray_data)))
lopit_df <- prepDataForggHeatmap(scale(exprs(Barylyuk2020ToxoLopit)[marker_genes, ]),
  row_order = row_order,
  col_order = seq(1, ncol(exprs(Barylyuk2020ToxoLopit)))
)
rnaseq_df <- prepDataForggHeatmap(as.matrix(rna_seq_data)[marker_genes, ], row_order = row_order)

lopit_df$Dataset <- "LOPIT"
micro_df$Dataset <- "Cell cycle"
rnaseq_df$Dataset <- "RNA seq"
heat_df <- rbind(micro_df, lopit_df, rnaseq_df)

heat_df |>
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() +
  facet_wrap(~Dataset, scales = "free") +
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")

microarray_data <- fread("./T_gondii/Original_data/ToxoDB_TgME49_Protein-coding_DNA_microarray.txt",
  na.strings = "N/A",
  strip.white = T,
  header = T,
  select = seq(1, 212)
)

rna_seq_data <- fread("./T_gondii/Original_data/ToxoDB_TgME49_Protein-coding_RNA-Seq.txt",
  na.strings = "N/A",
  strip.white = T,
  header = T,
  select = seq(1, 255)
)


m_white_cell_cycle
m_white_cell_cycle_normalised

tagm_comparison |> head()

marker_genes <- row.names(tagm_comparison)[tagm_comparison$markers != "unknown"]

pheatmap::pheatmap(m_white_cell_cycle_normalised[match(marker_genes, row.names(m_white_cell_cycle)), ], show_colnames = FALSE)

m_white_cell_cycle_normalised_markers <- m_white_cell_cycle_normalised[match(marker_genes, row.names(m_white_cell_cycle)), ]

m_white_cell_cycle_normalised_ordered <- m_white_cell_cycle_normalised[match(row.names(tagm_comparison), row.names(m_white_cell_cycle)), ]

row_order <- findOrder(m_white_cell_cycle_normalised_ordered[, -1])

micro_df <- prepDataForggHeatmap(m_white_cell_cycle_normalised_ordered[, -1], row_order = row_order, col_order = seq(1, ncol(microarray_data)))
lopit_df <- prepDataForggHeatmap(scale(as.matrix(final_protein_df[, seq(1, 30)])),
  row_order = row_order,
  col_order = seq(1, 30)
)
rnaseq_df <- prepDataForggHeatmap(as.matrix(normalised_rnaseq_macrophages_infected_by_T_gondii)[, ], row_order = row_order)

lopit_df$Dataset <- "LOPIT"
micro_df$Dataset <- "Cell cycle"
rnaseq_df$Dataset <- "RNA seq"
heat_df <- rbind(micro_df, lopit_df, rnaseq_df)

heat_df |>
  ggplot(aes(x = x, y = y, fill = Entry)) +
  geom_tile() +
  facet_wrap(~Dataset, scales = "free") +
  scale_fill_gradient2(mid = "#FFFFFF", low = "#146EB4", high = "#FF9900")


normalised_rnaseq_macrophages_infected_by_T_gondii |> head()

# === GO term over-representation ========================================================

data_modelled$data_modelled[[2]][fused_genes_1, ]

pred_cl[[1]][fused_genes_1]
pred_cl[[2]][fused_genes_1]

cc_pred <- mdi_mod$cell_cycle_predictions[[2]][fused_genes_1]

pred_df <- data.frame("LOPIT" = mdi_predictions, "Cell-cycle" = mdi_mod$cell_cycle_predictions$vi) # c(pred_cl[[2]]))
row.names(pred_df) <- row.names(data_modelled$data_modelled[[1]])

pred_df[fused_genes_1, ] |> table() # |> write.csv("~/Desktop/LOPIT_CellCycle_label_map.csv")
print(xtable::xtable(table(pred_df[fused_genes_1, ])))

dense_granule_cell_cycle <- which((pred_cl[[2]] == 3) & colMeans(mdi_mod$allocations[[1]] == mdi_mod$allocations[[2]]) > 0.5)
dense_granule_cell_cycle <- which((pred_cl[[2]] == 3) & (model_output$s[[1]] > 0.5))
data_modelled$data_modelled[[2]][dense_granule_cell_cycle, ] |>
  pheatmap(
    cluster_cols = FALSE,
    show_colnames = FALSE,
    main = "Dense granule localised proteins cell-cycle data",
    # filename = "./T_gondii/Analysis/cell_cycle_dense_granules.png"
  )

rhoptries_1_cell_cycle <- which((pred_cl[[2]] == 14) & colMeans(mdi_mod$allocations[[1]] == mdi_mod$allocations[[2]]) > 0.5)
# rhoptries_1_cell_cycle <- which((pred_cl[[2]] == 11) & (model_output$MDI$fusion_probabilities[[1]] > 0.5))
data_modelled$data_modelled[[2]][rhoptries_1_cell_cycle, ] |>
  pheatmap(
    cluster_cols = FALSE,
    show_colnames = FALSE,
    main = "Rhoptries 1 localised proteins cell-cycle data",
    filename = "./T_gondii/Analysis/cell_cycle_rhoptries_1.png"
  )

row.names(data_modelled$data_modelled[[1]])[dense_granule_cell_cycle]

fused_predictions <- pred_df[fused_genes_1, ]

fused_predictions
fixed_inds <- lopit_data$Fixed
names(fixed_inds) <- row.names(lopit_data)

fused_predictions$Fixed <- fixed_inds[fused_genes_1]

unfused_predictions <- pred_df[-fused_genes_1, ]
unfused_predictions$Fixed <- fixed_inds[-fused_genes_1]

# write.csv(fused_predictions[order(fused_predictions$Cell.cycle), ], file = "./T_gondii/Analysis/tGondiiFusedClusters.csv")
# write.csv(unfused_predictions[order(unfused_predictions$Cell.cycle), ], file = "./T_gondii/Analysis/tGondiiUnfusedClusters.csv")


# === Visualise results =======================================================


marker_labels <- label_to_organelle$Organelle[data_modelled$initial_labels[, 1]]
marker_labels[data_modelled$fixed[, 1] == 0] <- NA

anno_row <- data.frame(Localisation = factor(marker_labels, levels = sort(unique(marker_labels))))
row.names(anno_row) <- row.names(data_modelled$data_modelled[[1]])

ann_colours <- list(Localisation = viridis::viridis(26))
names(ann_colours$Localisation) <- sort(unique(label_to_organelle$Organelle))

pheatmap::pheatmap(data_modelled$data_modelled[[1]],
                   color = dataColPal(),
                   breaks = defineDataBreaks(data_modelled$data_modelled[[1]], dataColPal(), mid_point = 0),
                   annotation_row = anno_row,
                   annotation_colors = ann_colours,
                   show_rownames = FALSE,
                   cluster_cols = FALSE,
                   filename = "./T_gondii/Original_data/LOPIT_heatmap.png"
)

colnames(data_modelled$data_modelled[[2]]) <- paste0(seq(0, 12), " HR")
pheatmap::pheatmap(data_modelled$data_modelled[[2]],
                   color = dataColPal(),
                   breaks = defineDataBreaks(data_modelled$data_modelled[[2]], dataColPal(), mid_point = 0),
                   # annotation_row = anno_row,
                   # annotation_colors = ann_colours,
                   show_rownames = FALSE,
                   cluster_cols = FALSE,
                   filename = "./T_gondii/Original_data/CellCycle_heatmap.png"
)


# === Investigate cell cycle splitting/merging ======================================


p_occ_mdi <- (1 * (model_output$MDI$N_k[[2]] != 0)) |>
  pheatmap(cluster_cols = TRUE, cluster_rows = FALSE, main = "MDI", silent = TRUE, color = simColPal())

p_occ_mix <- (1 * (model_output$Cell_cycle$N_k[[1]] != 0)) |>
  pheatmap(
    cluster_cols = TRUE, cluster_rows = FALSE, main = "Mixture model", silent = TRUE,
    color = simColPal()
  )

rowSums((1 * (model_output$MDI$N_k[[2]] != 0))) |> mean()
rowSums((1 * (model_output$Cell_cycle$N_k[[1]] != 0))) |> mean()

combinePheatmaps(list(p_occ_mdi$gtable, p_occ_mix$gtable), save_name = "~/Desktop/cellCycleOccupiedComponents.png", main = "Cell cycle data occupied components")

model_output$Cell_cycle$pred[[1]] |> table()
model_output$MDI$pred[[2]] |> table()

cl_14 <- which(pred_cl[[2]] == 14)

pred_cl[[1]][cl_14]
pred_cl[[2]][cl_14]

data_modelled$data_modelled[[2]] |>
  annotatedHeatmap()

cl_14_order <- order(pred_cl[[1]][cl_14])

cl_mix_in_14 <- model_output$Cell_cycle$pred[[1]] %in% c(4, 5, 7, 24, 32, 37, 57, 76, 94, 96, 102, 105, 109)

pred_cell_cycle <- data.frame("MDI" = c(model_output$MDI$pred[[2]]), "Mixture.model" = c(model_output$Cell_cycle$pred[[1]]))

pred_cell_cycle |>
  table() |>
  apply(1, function(x) {
    x / sum(x)
  }) |>
  pheatmap(
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Fraction of mixture model clusters in MDI clusters",
    color = simColPal() # ,
    # filename = "~/Desktop/cellCycleMDIMixtureOverlap.png"
  )

annotatedHeatmap(microarray_data[cl_14[cl_14_order], ], model_output$Cell_cycle$pred[[1]][cl_14][cl_14_order],
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Cell-cycle data\nCluster 14 from MDI annotated by predicted cluster in mixture model" # ,
  # filename = "~/Desktop/cellCycleMDIclust14AnnotatedMixLabels.png"
)

annotatedHeatmap(microarray_data[cl_mix_in_14, ], model_output$MDI$pred[[2]][cl_mix_in_14],
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Cell-cycle data annotated by predicted cluster" # ,
  # filename = "./rnaseq_predicted_clustering.png"
)

annotatedHeatmap(microarray_data[cl_mix_in_14, ], model_output$Cell_cycle$pred[[1]][cl_mix_in_14],
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Cell-cycle data annotated by predicted cluster" # ,
  # filename = "./rnaseq_predicted_clustering.png"
)
