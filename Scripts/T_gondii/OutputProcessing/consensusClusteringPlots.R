
library(ggplot2)
library(patchwork)

mdiHelpR::setMyTheme()

x <- readRDS("./T_gondii/ConsensusClustering/CC_R_15000_K_125_V_2.rds")

x$stability$Dataset[x$stability$Dataset == "Cell_cycle"] <- "Cell cycle"

p1 <- x$stability |> 
  dplyr::filter(Quantity_varied == "Depth") |>
  ggplot(aes(x = Depth, y = Mean_absolute_difference, colour = factor(Width))) +
  geom_line() +
  facet_wrap(~Dataset) +
  ggthemes::scale_color_colorblind() +
  labs(y = "Mean absolute difference", color = "Width") +
  ylim(c(0, 0.025))

p2 <- x$stability |> 
  dplyr::filter(Quantity_varied == "Width") |>
  ggplot(aes(x = Width, y = Mean_absolute_difference, colour = factor(Depth))) +
  geom_line() +
  facet_wrap(~Dataset) +
  ggthemes::scale_color_colorblind() +
  labs(y = "Mean absolute difference", color = "Depth") +
  ylim(c(0, 0.025))

p3 <- x$phis |> 
  ggplot(aes(x = Phi_12, fill = Depth)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Depth, ncol = 1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = expression(phi[ "(1, 2)"]))


layout <- "
AACC
BBCC
"
p_patch <- p1 + p2 + p3 + 
  plot_layout(design = layout)

ggsave("./T_gondii/Analysis/consensus_stability.png", 
  plot = p_patch,
  height = 8, 
  width = 12
)

x$phis |> 
  dplyr::filter(Depth == 15000) |> 
  ggplot(aes(x = Phi_12, y = Depth)) +
  geom_boxplot() +
  ggthemes::scale_fill_colorblind()

x$stability |> 
   dplyr::filter(Quantity_varied == "Depth", Dataset == "LOPIT") |>
   ggplot(aes(x = Depth, y = Adjusted_Rand_index, colour = factor(Width))) +
   geom_line() +
   # facet_wrap(~Dataset) +
   ggthemes::scale_color_colorblind() +
   labs(y = "Adjusted Rand index", color = "Width") +
   ylim(c(0.5, 1.0))

x$stability |> 
  dplyr::filter(Quantity_varied == "Width") |>
  ggplot(aes(x = Width, y = Adjusted_Rand_index, colour = factor(Depth))) +
  geom_line() +
  facet_wrap(~Dataset) +
  ggthemes::scale_color_colorblind() +
  labs(y = "Adjusted Rand index", color = "Width") +
  ylim(c(0.5, 1.0)) 



z <- readRDS("./T_gondii/ConsensusClustering/CC_R_15000_K_300_cell_cycle_mix.rds")

cm_lst <- list()
cm_lst[[1]] <- x$cms[[1]]
cm_lst[[2]] <- x$cms[[2]]
cm_lst[[3]] <- z$cm[[1]]

row.names(cm_lst[[3]]) <- colnames(cm_lst[[3]]) <- row.names(cm_lst[[1]])

lopit_cm_df <- cm_lst[1] |> 
  tagmReDraft::prepSimilarityMatricesForGGplot(matrix_setting_order = 1)

cellcycle_cm_df <- cm_lst[c(2, 3)] |> 
  tagmReDraft::prepSimilarityMatricesForGGplot(matrix_setting_order = 1)

lopit_cm_df$Case <- "MDI: LOPIT"
cellcycle_cm_df$Case <- "MDI: Cell cycle"
cellcycle_cm_df$Case[cellcycle_cm_df$Chain == 2] <- "Mixture model: Cell cycle"
cm_df <- rbind(lopit_cm_df, cellcycle_cm_df)

cm_df$Case <- factor(cm_df$Case, levels = c("MDI: LOPIT", "MDI: Cell cycle", "Mixture model: Cell cycle"))

data.table::fwrite(cm_df, "./T_gondii/ConsensusClustering/CMsComparison.csv")

p_cm <- cm_df |> 
  ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = Entry)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~Case) +
  ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 10.5),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text.x = element_text(size = 10.5),
    legend.text = element_text(size = 10.5)
  ) +
  labs(fill = "Coclustering\nproportion", x = "", y = "")

ggsave(p_cm, filename = "./CM_comparison_t_gondii.png", height = 12, width = 30)

