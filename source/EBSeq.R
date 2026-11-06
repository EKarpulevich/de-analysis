# Calculate diff expression analysis
ebseq_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(EBSeq))

    # Prepare data
    condition <- coldata$condition
    # Calculate diff expression
    Sizes = MedianNorm(as.matrix(counts))
    EBOut = EBTest(Data = as.matrix(counts),
                    Conditions = as.factor(condition),
                    sizeFactors = Sizes, maxround = 5)
    Out = GetDEResults(EBOut)
    res <- Out$PPMat
    # Summary for results
    summary(res)
    # Return results of diff expression analysis
    return(EBOut)
}

# Visualize function
ebseq_v <- function(ebseq_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(EBSeq))
    suppressMessages(library(EnhancedVolcano))
    res <- GetDEResults(ebseq_res)$PPMat
    GeneFC <- PostFC(ebseq_res, SmallNum = 0.05)    # Create dataframe with groups
    x = -log2(GeneFC$PostFC)
    y = ebseq_res$PPMat[,"PPDE"]

    xcond = 1   # Change group colors
    ycond = 0.05
    group <- ifelse((x < -xcond | x > xcond) & y > ycond , "DE", ifelse(y > ycond , "ADE2", ifelse((x < -xcond | x > xcond), "ADE1", "NDE")))
    df <- data.frame(x = x, y = y, group = factor(group))
    colors <- c("red", "blue", "green", "gray")

      # Make plot
    plot(df$x, df$y, col = colors[df$group], pch = 16,
        xlim =c(-5,5),
        ylim=c(0,1),
        main= "EBSeq results",
        xlab="EBSeq Log2 Fold Change", ylab="EBSeq PPDE")
    DEgplots <- data.frame(x = ifelse((x < -xcond | x > xcond) & y > ycond, x, 0), y = ifelse((x < -xcond | x > xcond) & y > ycond , y, 0),
        name = ifelse((x < -xcond | x > xcond) & y > ycond , rownames(res), ""))
    text(DEgplots$x, DEgplots$y,
         labels = DEgplots$name,
         cex = 0.6, pos = 4, col = "black")

}

# Make a signature of top-n genes
ebseq_top <- function(results, n, out) {
    # Import libraries
    suppressMessages(library("EBSeq"))
    suppressMessages(library("biomaRt"))

    # Extract results of analysis
    ppmat = GetDEResults(results)$PPMat
    GeneFC <- PostFC(results, SmallNum = 0.01)
    postfc <- as.data.frame(GeneFC$PostFC)
    colnames(postfc) <- c("PostFC")

    # Extract top-n differentially expressed genes ordered by PostFC
    top <- head(postfc[order(postfc$PostFC, decreasing = TRUE), , drop = FALSE], n)
    tmp <- gsub("\\..*","",row.names(top))

    GeneFC <- PostFC(results, SmallNum = 0.00)    # Create dataframe with groups
    filtered_info <- data.frame(results$PPMat[,"PPDE"], -log2(GeneFC$PostFC))
    colnames(filtered_info) <- c("PPDE", "log2PostFC")
    write.csv(filtered_info, file = file.path(out, "EBSeq_info.csv"))
    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)


}

# Filter genes by PPDE
ebseq_filtered <- function(results, out) {
    suppressMessages(library("EBSeq"))

    # filtering results by PPDE >= 0.95
    ppde <- results$PPMat[,"PPDE"]
    ppde <- ppde[ppde >= 0.95]
    ppde <- as.data.frame(ppde)
    colnames(ppde) <- c("PPDE")
    filtered_genes <- gsub("\\..*","",row.names(ppde))
    print(str(results))

    GeneFC <- PostFC(results, SmallNum = 0.00)    # Create dataframe with groups
    filtered_info <- data.frame(results$PPMat[,"PPDE"], -log2(GeneFC$PostFC))
    colnames(filtered_info) <- c("PPDE", "log2PostFC")

    write.csv(filtered_info, file = file.path(out, "EBSeq_info.csv"))
    return(filtered_genes)
}
