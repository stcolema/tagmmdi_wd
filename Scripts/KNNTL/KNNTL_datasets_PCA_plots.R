
suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

setMyTheme()

data("HEK293T2011")
data("groen2014r1")
data("dunkley2006")
data("E14TG2aS1")

# === Human dataset ============================================================

human_marker_proteins_df <- pRoloc:::subsetAsDataFrame(HEK293T2011, "markers", train = TRUE)
human_lopit <- human_marker_proteins_df |> 
  dplyr::select(-markers)

human_pcs <- prcomp(human_lopit, center = TRUE, scale. = TRUE)

human_df <- as.data.frame(human_pcs$x[ , c(1, 2)]) |> 
  mutate(Marker = human_marker_proteins_df$markers, Dataset = "Human")

p_human <- human_df |> 
  ggplot(aes(x = PC1, y = PC2, color = Marker)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Dataset)

# === Root dataset =============================================================

root_marker_proteins_df <- pRoloc:::subsetAsDataFrame(groen2014r1, "markers", train = TRUE)
root_lopit <- root_marker_proteins_df |> 
  dplyr::select(-markers)

root_pcs <- prcomp(root_lopit, center = TRUE, scale. = TRUE)

root_df <- as.data.frame(root_pcs$x[ , c(1, 2)]) |> 
  mutate(Marker = root_marker_proteins_df$markers, Dataset = "Root")

p_root <- root_df |> 
  ggplot(aes(x = PC1, y = PC2, color = Marker)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Dataset)

# === Callus dataset ===========================================================

callus_marker_proteins_df <- pRoloc:::subsetAsDataFrame(dunkley2006, "markers", train = TRUE)
callus_lopit <- callus_marker_proteins_df |> 
  dplyr::select(-markers)

callus_pcs <- prcomp(callus_lopit, center = TRUE, scale. = TRUE)

callus_df <- as.data.frame(callus_pcs$x[ , c(1, 2)]) |> 
  mutate(Marker = callus_marker_proteins_df$markers, Dataset = "Callus")

p_callus <- callus_df |> 
  ggplot(aes(x = PC1, y = PC2, color = Marker)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Dataset)

# === Mouse dataset ============================================================

mouse_marker_proteins_df <- pRoloc:::subsetAsDataFrame(E14TG2aS1, "markers", train = TRUE)
mouse_lopit <- mouse_marker_proteins_df |> 
  dplyr::select(-markers)

mouse_pcs <- prcomp(mouse_lopit, center = TRUE, scale. = TRUE)

mouse_df <- as.data.frame(mouse_pcs$x[ , c(1, 2)]) |> 
  mutate(Marker = mouse_marker_proteins_df$markers, Dataset = "Mouse")

p_mouse <- mouse_df |> 
  ggplot(aes(x = PC1, y = PC2, color = Marker)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Dataset)

# === Visualisation ============================================================

plot_df <- rbind(callus_df, root_df, human_df, mouse_df)

plot_df |> 
  ggplot(aes(x = PC1, y = PC2, color = Marker)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Dataset)

p_patch <- (p_callus + p_root) / (p_human + p_mouse) + 
  plot_annotation(tag_levels = 'A')

ggsave("./KNN_TL_comp/DatasetsPCPlots.png", height = 6, width = 10)
