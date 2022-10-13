
#!/usr/bin/env Rscript
#
# Title: KNN Transfer Learner Timing plots
# Description: Plot the time taken for the KNN-TL and Bayesian models on each dataset
#
# Example call:
# Rscript knnTLCVSingleLoop.R --datasets 'tan2009r1 tan2009r1goCC' --seed 1 
# --test_indices ./test_50/tan2009r1/tan2009r1_seed_1_testSize_50_trainingAdjustedForTL_TRUE.rds
# --save_dir ./test_50/ --categorical_column_threshold 3 --number_weights 5
# --number_weights_sampled 20000


suppressMessages(library(pRolocdata))
suppressMessages(library(pRoloc))
suppressMessages(library(ggplot2))
suppressMessages(library(tagmReDraft))
suppressMessages(library(mdiHelpR))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(patchwork))
suppressPackageStartupMessages(library(MSnbase))

mdiHelpR::setMyTheme()


# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    
    # Save the output to this directory
    optparse::make_option(c("-d", "--data_dir"),
                          type = "character",
                          default = "./CV_output/",
                          help = "Directory to read model output from [default= %default]",
                          metavar = "character"
    ),
    
    # Save the output to this directory
    optparse::make_option(c("-d", "--save_dir"),
                          type = "character",
                          default = "./",
                          help = "Directory to save output to [default= %default]",
                          metavar = "character"
    )
  )
  
  
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

args <- input_arguments()

data_dir <- args$data_dir
save_dir <- args$save_dir

test_splits <- c(10, 30, 50)
test_dirs <- paste0(data_dir, "test_", test_splits, "/")
n_splits <- length(test_splits)

datasets <- c("dunkley2006",
              "E14TG2aS1",
              "groen2014r1",
              "HEK293T2011"
)
n_datasets <- length(datasets)

time_df <- NULL
for(ii in seq(1, n_splits)) {
  for(jj in seq(1, n_datasets)) {
    
    test_frac <- test_splits[ii]
    dataset <- datasets[jj]
    curr_dir <- paste0(test_dirs[ii], dataset)
    
    bayesian_pattern <- paste0("*_nChains_5_testSize_", test_frac, "*")
    bayesian_files <- list.files(curr_dir, pattern = bayesian_pattern, full.names = TRUE)
    n_bayesian_files <- length(bayesian_files)
    for(kk in seq(1, n_bayesian_files)) {
      f <- bayesian_files[kk]
      x <- readRDS(f)
      mdi_mcmc <- x$mdi_fold$mcmc_output
      n_chains <- length(mdi_mcmc)
      t <- 0
      for(ll in seq(1, n_chains)) {
        t_ll <- mdi_mcmc[[ll]]$Time
        t <- t + t_ll
      }
      entry <- data.frame(Model = "MDI", Dataset = dataset, Time = t, Test_frac = test_frac)
      
      if(is.null(time_df)) {
        time_df <- entry
      } else {
        time_df <- rbind(time_df, entry)
      }
      
      tagm <- x$tagm_fold$mcmc_output
      n_chains <- length(tagm)
      t <- 0
      for(ll in seq(1, n_chains)) {
        t_ll <- tagm[[ll]]$Time
        t <- t + t_ll
      }
      entry <- data.frame(Model = "TAGM", Dataset = dataset, Time = t, Test_frac = test_frac)
      time_df <- rbind(time_df, entry)
    }
    
    knn_tl_pattern <- paste0("*_knnTL_numberWeights_5_*")
    knn_tl_files <- list.files(curr_dir, pattern = knn_tl_pattern, full.names = TRUE)
    n_knn_tl_files <- length(knn_tl_files)
    
    for(kk in seq(1, n_knn_tl_files)) {
      f <- knn_tl_files[kk]
      x <- readRDS(f)
      t <- x$knn_fold$time_taken
      entry <- data.frame(Model = "KNN-TL", Dataset = dataset, Time = t, Test_frac = test_frac)
      time_df <- rbind(time_df, entry)
    }
  }
}

write.csv(time_df, "./KNN_TL_time_results.csv")
# time_df <- read.csv("./KNN_TL_comp/KNN_TL_time_results.csv", row.names = 1)

time_df$Dataset <- factor(time_df$Dataset, 
  levels = c("dunkley2006", "groen2014r1", "HEK293T2011", "E14TG2aS1"),
  labels = c("Callus", "Root", "Human", "Mouse")
)

time_df$Model <- factor(time_df$Model, levels = c("MDI", "TAGM", "KNN-TL"))

time_df$`Test fraction` <- time_df$Test_frac * 0.01

p_time <- time_df |> 
  ggplot(aes(x = Model, y = Time, fill = Model)) +
  geom_boxplot() +
  facet_grid(`Test fraction` ~ Dataset, labeller = label_both) +
  ggthemes::scale_fill_colorblind() +
  labs(y = "Time [seconds]")

p_log_time <- time_df |> 
  ggplot(aes(x = Model, y = log(Time, 10), fill = Model)) +
  geom_boxplot() +
  facet_grid(`Test fraction` ~ Dataset, labeller = label_both) +
  ggthemes::scale_fill_colorblind() +
  labs(y = expression(paste(log[10], " (Time) [", log[10], "(seconds)]", sep = "")))

p_patch <- p_time + p_log_time + plot_layout(guides = 'collect')

ggsave("./KNN_TL_comp/log_time_plot.png", p_log_time)
ggsave("./KNN_TL_comp/time_plot.png", p_time)
ggsave("./KNN_TL_comp/time_plot_patchwork.png", p_patch, height = 8, width = 16)
