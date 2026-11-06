# Calculate diff expression analysis
deseq2_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(DESeq2))

    # Prepare data
    condition <- coldata$condition
    contrast <- c("condition", unique(coldata[c("condition")])$condition[1], unique(coldata[c("condition")])$condition[2])
    dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)), colData=coldata, design = ~condition)
    dds$condition <- relevel(dds$condition, unique(coldata[c("condition")])$condition[2])
    # Calculate diff expression
    dds <- DESeq(dds)
    res <- results(dds)
    res <- lfcShrink(dds, contrast=contrast,
                     res=res, type='normal')
    res <- results(dds, contrast=contrast)
    print(res)
    resultsNames(dds)
    # Summary for results
    summary(res)
    # Return results of diff expression analysis
    return(res)
}

# Visualize function
deseq2_v <- function(deseq2_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(DESeq2))
    suppressMessages(library(EnhancedVolcano))

    vis_res <- EnhancedVolcano(deseq2_res,
            lab = rownames(deseq2_res),
            x = 'log2FoldChange',
            y = 'pvalue',
            pCutoff = 0.05,
            FCcutoff = 1,
            title = "DESeq2 results")
    return(vis_res)
}

# Make a signature of top-n genes
deseq2_top <- function(results, n, out) {
    # Import libraries
    suppressMessages(library("DESeq2"))
    suppressMessages(library("biomaRt"))

    # Filter results by logFC > 0.25 or logFC < -0.25
    filtered_results <- results[!is.na(results$log2FoldChange), ]
    filtered_results <- results[!is.na(results$pvalue), ]
    filtered_results_info <- as.data.frame(do.call(cbind, filtered_results@listData))
    row.names(filtered_results_info) <- row.names(filtered_results)
    write.csv(filtered_results_info, file = file.path(out, "DESeq_info.csv"))

    filtered_results <- filtered_results[abs(filtered_results$log2FoldChange) > 0.25, ]

    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$pvalue), ], n)
    tmp <- gsub("\\..*","",row.names(top))


    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and p-value
deseq2_filtered <- function(results, out) {
  suppressMessages(library("DESeq2"))

  # filtering results by log2FC >= 1 and p-value < 0.05
  filtered_results <- results[!is.na(results$log2FoldChange), ]
  filtered_results <- results[!is.na(results$padj), ]
  filtered_results_info <- as.data.frame(do.call(cbind, filtered_results@listData))
  row.names(filtered_results_info) <- row.names(filtered_results)
  write.csv(filtered_results_info, file = file.path(out, "DESeq_info.csv"))

  filtered_results <- filtered_results[abs(filtered_results$log2FoldChange) >= 1, ]
  filtered_results <- filtered_results[abs(filtered_results$padj) < 0.05, ]
  filtered_genes <- gsub("\\..*","",row.names(filtered_results))


  return(filtered_genes)
}




