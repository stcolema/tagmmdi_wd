# Script plotting the sampled phis for a single chain from the dunkley data

library(tidyr)
library(dplyr)
library(mdir)

cc_file <- "./T_gondii/ConsensusClustering/CC_R_15000_K_125_full.rds"

plot_dir <- "./Plots/"
plot_name <- "TGondiiLOPITCM."
devices <- c("png", "pdf")

x <- readRDS(cc_file)

sim_col <- mdiHelpR::simColPal()
cm_breaks <- mdiHelpR::defineBreaks(sim_col, lb = 0, ub = 1)

pheatmap::pheatmap(x$cm[[1]], 
                   color = sim_col, 
                   breaks = cm_breaks, 
                   filename = paste0(plot_dir, plot_name, devices[2]), 
                   silent = TRUE, 
                   show_rownames = FALSE, 
                   show_colnames = FALSE
)
