# Script plotting the sampled phis for a single chain from the dunkley data

library(tidyr)
library(dplyr)
library(mdir)

cc_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_125_full.rds"

plot_dir <- "./Plots/"
plot_name <- "TGondiiSampledPhis."
devices <- c("png", "pdf")

x <- readRDS(cc_file)
phi_df <- data.frame("Sampled_phi" = x$phis)

p1 <- phi_df |>
  ggplot(aes(x = Sampled_phi)) +
  geom_density(alpha = 0.3, fill = "#000000") +
  ggthemes::scale_fill_colorblind() +
  labs(x = expression(paste("Sampled ", phi)), y = "Density")

## Number of merged proteins/genes
# sum(colMeans(x$allocations[[1]] == x$allocations[[2]]) > 0.5)

for (dev in devices) {
  ggsave(paste0(plot_dir, "/", plot_name, dev),
         plot = p1,
         device = dev,
         height = 6,
         width = 5
  )
}

# p1 <- x$mdi_fold$phi_df |>
#   pivot_longer(-Iteration, names_to = "Chain", values_to = "Sampled_phi", names_prefix = "Chain.") |>
#   ggplot(aes(x = Sampled_phi, fill = Chain)) +
#   geom_density(alpha = 0.3) +
#   ggthemes::scale_fill_colorblind() +
#   labs(x = expression(paste("Sampled ", phi)), y = "Density")
