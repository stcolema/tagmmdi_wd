
prepCMssForGGplot <- function(cms, model_description_df,
        matrix_setting_order = 1,
        use_common_ordering = TRUE) {
  not_list <- !is.list(cms)
  if (not_list) {
    stop("`similarity_matrices` must be a list of matrices.")
  }
  
  n_matrices <- length(cms)
  n_models <- length(model_description_df)
  if(n_matrices != n_models) {
    .err <- paste("Number of consensus matrices and number of models described",
      "in ``model_description_df`` are not matching."
    )
    stop(.err)
  }
  row_order <- col_order <- findOrder(cms[[matrix_setting_order]])
  
  depths <- model_description_df$Depth
  widths <- model_description_df$Width
  
  for (ii in seq(1, n_matrices)) {
    first_iteration <- ii == 1
    
    d <- depths[ii]
    w <- widths[ii]
    .df <- prepDataForggHeatmap(cms[[ii]], row_order, col_order)
    .df$Depth <- d
    .df$Width <- w
    
    if (first_iteration) {
      cm_df <- .df
    } else {
      cm_df <- rbind(cm_df, .df)
    }
  }
  cm_df$Depth <- factor(cm_df$Depth)
  cm_df$Width <- factor(cm_df$Width)
  cm_df
}