# Calculate certain dist matrix
calculate_distmatrix <- function(countMatrixFile, signatureFile) {
    # Download library
    suppressMessages(library("amap"))

    # Download data from files
    signature <- readLines(signatureFile)
    cm <- read.table(countMatrixFile, header=T, sep=",")
    cm <- subset(cm, !duplicated(cm$X))
    cm <- data.frame(cm[, -1], row.names=cm[, 1])

    # Calculate dist matrix
    cm_subset <- cm[signature, ]
    cat("SUBSET SHAPE: ")
    cat(dim(cm_subset))
    cat("\n")
    distMatrix <- Dist(t(cm_subset), method="kendall", nbproc=6)

    return(distMatrix)
}

heatmap_v <- function(countMatrixFile, tool_name, condition, sorting_condition, out) {
    signatureFile <- file.path(out, paste0(tool_name, "_sig.txt"))
    if (file.info(signatureFile)$size == 0) {
        cat(paste0('There is no gene in ', tool_name, ' for that signature\n'))
        return()
    }
    signature <- readLines(signatureFile)
    if (length(signature) == 1) {
        cat(paste0('There is not enough genes (only 1) in ', tool_name, ' signature\n'))
        return()
    }
    cm <- read.table(countMatrixFile, header=T, sep=",")
    cm <- subset(cm, !duplicated(cm$X))
    cm <- data.frame(cm[, -1], row.names=cm[, 1])
    cm_subset <- cm[signature, ]
    #cm_subset <- filter_if(cm_subset, is.numeric, all_vars((.) != 0))
    cm_subset <- cm_subset[is.finite(rowSums(cm_subset)),]
    my_sample_col <- data.frame(condition)
    row.names(my_sample_col) <- colnames(cm_subset)
    cm_subset <- cm_subset[, row.names(my_sample_col)]
    cal_z_score <- function(x){
      if(max(x) - min(x) != 0) {
	    (x - min(x)) / (max(x) - min(x))
      } else {
	    0
      }

    }

    data_subset_norm <- t(apply(cm_subset, 1, cal_z_score))

    library('pheatmap')
    pheatmap(data_subset_norm, annotation_col = my_sample_col, annotation_names_col = FALSE, annotation_legend = TRUE,
            annotation_names_row = FALSE, main = paste0(tool_name, ' DE analysis', sorting_condition), drop_levels = FALSE, show_rownames = F, show_colnames = T,
            cluster_cols = T, cutree_rows = 2, cutree_cols = 4)

}
