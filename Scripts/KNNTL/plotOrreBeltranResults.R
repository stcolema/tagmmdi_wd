# This script takes the outputs of the Orre and Beltran analyses and plots them

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggthemes)
library(patchwork)

devices <- c("pdf")
save_dir <- "./Plots/ValidationModernData/"

orre_filename <- "orrePerf."
beltran_filename <- "beltranPerf."
patch_filename <- "orreBeltranPatch."

orre_df <- read.csv("./KNN_TL_comp/Orre_perf.csv", row.names = 1) |>
  pivot_longer(-c(Model, Modality), names_to = "Score", values_to = "Value")

beltran_df <- read.csv("./KNN_TL_comp/Beltran_perf.csv", row.names = 1) |>
  pivot_longer(-c(Model, Modality), names_to = "Score", values_to = "Value")

orre_df$Score <- factor(orre_df$Score,
  levels = c("Accuracy", "Macro_F1", "Weighted_F1", "Brier_loss"),
  labels = c("Accuracy", "Macro F1", "Weighted F1", "Brier loss")
)

orre_df$Modality <- factor(orre_df$Modality,
  levels = c("Epidermoid carcinoma (A431)", "Lung cancer (NCI-H322)", "Lung cancer (HCC-827)", "Breast cancer (MCF7)", "Glioblastoma (U251)"),
  labels = c("A431", "NCI-H322", "HCC-827", "MCF7", "U251")
)



beltran_df$Score <- factor(beltran_df$Score,
  levels = c("Accuracy", "Macro_F1", "Weighted_F1", "Brier_loss"),
  labels = c("Accuracy", "Macro F1", "Weighted F1", "Brier loss")
)

beltran_df$Modality <- factor(beltran_df$Modality,
  levels = c("HCMV24", "HCMV48", "HCMV72", "HCMV96", "HCMV120"),
  labels = paste0(c(24, 48, 72, 96, 120), " hours\npost infection")
)

col_pal <- ggthemes::colorblind_pal()(4)[c(2, 3)]
names(col_pal) <- c("MDI", "TAGM")

p_beltran <- beltran_df |>
  dplyr::filter(Score %in% c("Accuracy", "Macro F1", "Brier loss")) |>
  ggplot(aes(y = Value, x = Model, fill = Model)) +
  geom_boxplot() +
  facet_grid(Score ~ Modality, scales = "free") + # , nrow = 1) +
  scale_fill_manual(values = col_pal)


p_orre <- orre_df |>
  dplyr::filter(Score %in% c("Accuracy", "Macro F1", "Brier loss")) |>
  ggplot(aes(y = Value, x = Model, fill = Model)) +
  geom_boxplot() +
  facet_grid(Score ~ Modality, scales = "free") + # , nrow = 1) +
  scale_fill_manual(values = col_pal)

p_patch <- p_beltran / p_orre +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")

for (dev in devices) {
  filename <- paste0(save_dir, orre_filename, dev)
  ggsave(filename, p_orre, device = dev)

  filename <- paste0(save_dir, beltran_filename, dev)
  ggsave(filename, p_beltran, device = dev)

  filename <- paste0(save_dir, patch_filename, dev)
  ggsave(filename, p_patch, device = dev)
}
