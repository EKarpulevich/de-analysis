# Calculate diff expression analysis
voom_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(edgeR))
    suppressMessages(library(limma))

    # Prepare data
    condition <- coldata$condition
    d <- DGEList(counts)
    # Calculate diff expression
    design <- model.matrix(~0 + condition, data = coldata)
    y <- voom(d, design)
    fit <- lmFit(y, design)
    head(coef(fit))
    # Make contrast
    contrast <- c(unique(coldata[c("condition")])$condition[1], unique(coldata[c("condition")])$condition[2])
    contrast[1] <- paste("condition", contrast[1], sep="")
    contrast[2] <- paste("condition", contrast[2], sep="")
    contrast <- paste(contrast[2], " - ", contrast[1])
    contrasts <- makeContrasts(contrasts = contrast, levels = colnames(coef(fit)))
    contrasts
    # Fit
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2, trend=TRUE)
    res <- topTable(fit2, sort.by = "P", n = Inf)
    head(res, 20)
    # Return results of diff expression analysis
    return(res)
}

# Visualize function
voom_v <- function(voom_res, out) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(edgeR))
    suppressMessages(library(limma))
    suppressMessages(library(EnhancedVolcano))

    #Plot
    vis_res <- EnhancedVolcano(voom_res,
            lab = rownames(voom_res),
            x = 'logFC',
            y = 'P.Value',
            pCutoff = 0.05,
            FCcutoff = 1,
            title = "limma-voom results")
    return(vis_res)
}

# Make a signature of top-n genes
voom_top <- function(results, n, out) {

    # Import libraries
    suppressMessages(library("edgeR"))
    suppressMessages(library("biomaRt"))

    # Filter results by logFC > 1 or logFC < -1
    filtered_results <- results[!is.na(results$logFC), ]
    filtered_results <- filtered_results[abs(filtered_results$logFC) > 0.25, ]
    write.csv(filtered_results, file = file.path(out, "voom_info.csv"))

    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$P.Value), ], n)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and p-value
voom_filtered <- function(results, out) {
    suppressMessages(library("edgeR"))
    # filtering results by log2FC >= 1 and p-value < 0.05
    filtered_results <- results[!is.na(results$logFC), ]
    filtered_results <- results[!is.na(results$P.Value), ]
    write.csv(filtered_results, file = file.path(out, "voom_info.csv"))

    filtered_results <- filtered_results[abs(filtered_results$logFC) >= 1, ]
    filtered_results <- filtered_results[abs(filtered_results$P.Value) < 0.05, ]
    filtered_genes <- gsub("\\..*","",row.names(filtered_results))

    return(filtered_genes)
}
