# Script plotting the sampled phis for a single chain from the dunkley data

library(tidyr)
library(dplyr)
library(mdir)

n_seeds <- 10
file_used <- 1
output_dir <- "./KNN_TL_comp/test_30/"
datasets <- c("dunkley2006", "E14TG2aS1", "groen2014r1", "HEK293T2011")
dataset_names <- c("Callus", "Root", "Mouse", "Human")

plot_dir <- "./Plots"
plot_name <- "validationStudyPhis."
devices <- c("png", "pdf")

gen_filename <- paste0("_R_15000_seed_", seq(1, n_seeds), "_nChains_5_testSize_30_K_50_trainingAdjustedForTL_TRUE.rds")

my_df <- NULL
# dat <- datasets[1]
for (dat in datasets) {
  filenames <- paste0(output_dir, dat, "/", dat, "_", dat, "goCC", gen_filename)

  x <- readRDS(filenames[file_used])

  phis <- x$mdi_fold$phi_df |>
    pivot_longer(-Iteration, names_to = "Chain", values_to = "Sampled_phi", names_prefix = "Chain.") |>
    mutate(Dataset = dat)

  if (is.null(my_df)) {
    my_df <- phis
  } else {
    my_df <- rbind(my_df, phis)
  }
}

my_df$Dataset <- factor(my_df$Dataset, levels = datasets, labels = dataset_names)

p1 <- my_df |>
  ggplot(aes(x = Sampled_phi, fill = Chain)) +
  geom_density(alpha = 0.3) +
  ggthemes::scale_fill_colorblind() +
  labs(x = expression(paste("Sampled ", phi)), y = "Density") +
  facet_wrap(~Dataset)

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
