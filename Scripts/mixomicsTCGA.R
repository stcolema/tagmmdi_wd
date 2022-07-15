
library(klic)
library(abind)
library(tagmReDraft)
library(mixOmics)
library(mdiHelpR)
library(magrittr)


#' @title multi-class F1
#' @description Calculate the accuracy, the per-class F1 score and the macro-
#' and weighted-F1 scores.
#' @param pred Factor. The predicted classes.
#' @param truth Factor wit the same levels as ``pred``. The true classification.
#' @return A list containing ``accuracy``, the prediction accuracy, ``f1``, a
#' vector of the f1 score in each class, ``macro_f1``, the macro f1 score across
#' classes, and ``weighted_f1``, the f1 score weighted by class membership.
multiClassF1 <- function(pred, truth) {
  mismatch_in_lengths <- length(pred) != length(truth)
  if (mismatch_in_lengths) {
    stop("Prediction vector and truth must be of the same length.")
  }
  
  # Confusion matrix for current fold
  conf <- caret::confusionMatrix(
    data = pred,
    reference = truth
  )$table
  
  conf <- t( get.confusion_matrix(
    truth = truth, 
    predicted = pred
  ) )
  
  N <- length(pred)
  n_levels <- ncol(conf)
  seq_levels <- seq(1, n_levels)
  
  f1 <- rep(0, n_levels)
  
  for (ii in seq_levels) {
    precision <- conf[ii, ii] / sum(conf[ii, ])
    recall <- conf[ii, ii] / sum(conf[, ii])
    f1[ii] <- (2 * precision * recall) / (precision + recall)
    if (is.na(f1[ii])) {
      f1[ii] <- 0
    }
  }
  
  macro_f1 <- mean(f1)
  weighted_f1 <- sum(f1 * colSums(conf)) / N
  accuracy <- sum(pred == truth, na.rm = TRUE) / N
  list(
    "accuracy" = accuracy,
    "f1" = f1,
    "macro_f1" = macro_f1,
    "weighted_f1" = weighted_f1
  )
}

data("breast.TCGA")

summary(breast.TCGA)

N_train <- breast.TCGA$data.train$mirna |> nrow()
N_test <- breast.TCGA$data.test$mirna |> nrow()
N <- N_train + N_test

train_inds <- seq(1, N_train)
test_inds <- N_train + seq(1, N_test)

data_modelled <- list(
  "mirna" = scale(rbind(breast.TCGA$data.train$mirna,
        breast.TCGA$data.test$mirna
  )),
  "mrna" = scale(rbind(breast.TCGA$data.train$mrna,
        breast.TCGA$data.test$mrna
  )),
  "subtype" = c(breast.TCGA$data.train$subtype,
        breast.TCGA$data.test$subtype
  )
)

annotatedHeatmap(scale(data_modelled$mrna), data_modelled$subtype)
annotatedHeatmap(t(scale(t(data_modelled$mrna))), data_modelled$subtype)

annotatedHeatmap(scale(data_modelled$mirna), data_modelled$subtype)
annotatedHeatmap(t(scale(t(data_modelled$mirna))), data_modelled$subtype)

initial_labels <- fixed <- matrix(0, nrow = N, ncol = 2)
initial_labels[train_inds, 1] <- as.numeric(breast.TCGA$data.train$subtype) 
initial_labels[train_inds, 2] <- as.numeric(breast.TCGA$data.train$subtype) 
fixed[train_inds, ] <- 1

R <- 10000
thin <- 50
types <- c("G", "G")
mcmc_out <- callMDI(data_modelled[1:2], R = R, thin = thin, initial_labels = initial_labels, fixed = fixed, types = types, K = c(3, 3))

burn <- 0.2 * R
new_out <- processMCMCChain(mcmc_out, burn)

truth <- data_modelled$subtype[test_inds]
pred_v1 <- factor(new_out$pred[[1]], labels = c("Basal", "Her2", "LumA"))
pred_v2 <- factor(new_out$pred[[2]], labels = c("Basal", "Her2", "LumA"))

multiClassF1(pred_v1[test_inds], truth)
multiClassF1(pred_v2[test_inds], truth)

psm_v1 <- new_out$allocations[, , 1] |> createSimilarityMat()
psm_v2 <- new_out$allocations[, , 2] |> createSimilarityMat()

annotatedHeatmap(data_modelled$mirna, new_out$pred[[1]])
annotatedHeatmap(data_modelled$mirna, new_out$pred[[2]])

row_names <- data_modelled$mirna |> row.names()
tumour_types <- levels(data_modelled$subtype)

# Create the annotation data.frame for the rows
anno_row <- data.frame("Truth" = data_modelled$subtype,
                       "MIRNA" = pred_v1,
                       "MRNA" = pred_v2) %>%
  magrittr::set_rownames(row_names)

# Create the annotation colours
ann_colours <- list("Truth" = ggthemes::colorblind_pal()(3),
                    "MIRNA" = ggthemes::colorblind_pal()(3),
                    "MRNA" = ggthemes::colorblind_pal()(3)
)

names(ann_colours$Truth) <- tumour_types
names(ann_colours$MIRNA) <- tumour_types
names(ann_colours$MRNA) <- tumour_types

# Create the heatmap
pheatmap::pheatmap(data_modelled$mirna,
                         color = dataColPal(),
                         # breaks = my_breaks,
                         annotation_row = anno_row[, c(1, 2)],
                         annotation_colors = ann_colours[c(1, 2)]
)

pheatmap::pheatmap(data_modelled$mrna,
                   color = dataColPal(),
                   # breaks = my_breaks,
                   annotation_row = anno_row[, c(1, 3)],
                   annotation_colors = ann_colours[c(1, 3)]
)




# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
summary(Y)

list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)
circosPlot(MyResult.diablo, cutoff=0.7)
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")

# prepare test set data: here one block (proteins) is missing
X.test <- list(mRNA = breast.TCGA$data.test$mrna, 
               miRNA = breast.TCGA$data.test$mirna)

Mypredict.diablo <- predict(MyResult.diablo, newdata = X.test)

diablo_pred <- factor(Mypredict.diablo$MajorityVote$centroids.dist[,2])

confusion.mat <- get.confusion_matrix(
  truth = breast.TCGA$data.test$subtype, 
  predicted = Mypredict.diablo$MajorityVote$centroids.dist[,2])
knitr::kable(confusion.mat)

multiClassF1(diablo_pred, breast.TCGA$data.test$subtype)
multiClassF1(pred_v1[test_inds], truth)
multiClassF1(pred_v2[test_inds], truth)

psms <- list(psm_v1, psm_v2)

psm_df <- prepSimilarityMatricesForGGplot(psms)
psm_df |> 
  ggplot(aes(x = x, y= y, fill = Entry)) + 
  geom_tile() + 
  facet_grid(~Chain) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#FDF9DC", colour = "#FDF9DC")
  )

prob_global <- (new_out$allocation_probability[[1]] + new_out$allocation_probability[[2]]) / 2
pred_global <- factor(apply(prob_global, 1, which.max), labels = c("Basal", "Her2", "LumA"))

new_out$phis |> plot()

multiClassF1(pred_v1[test_inds], truth)
multiClassF1(pred_v2[test_inds], truth)
multiClassF1(pred_global[test_inds], truth)
multiClassF1(diablo_pred, breast.TCGA$data.test$subtype)
