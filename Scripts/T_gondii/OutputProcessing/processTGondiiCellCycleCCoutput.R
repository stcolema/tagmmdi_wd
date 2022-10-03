
library(tagmReDraft)

file_path <- "./"
pattern <- "TGondii_CellCycle_seed_*"
files <- list.files(file_path, pattern = pattern)
n_files <- length(files)

mcmc_lst <- list()
for(ii in seq(1, n_files)) {
  f <- files[ii]
  .x <- readRDS(f)[[1]]
  .x$Chain <- ii
  mcmc_lst[[ii]] <- .x
}

cc_out <- compileConsensusClustering(mcmc_lst)

cc_out$pred <- salso::salso(cc_out$allocations[[1]])

file_name <- paste0(file_path, "CC_R_15000_K_300_cell_cycle_mix.rds")
saveRDS(cc_out, file_name)
