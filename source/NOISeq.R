# Calculate diff expression analysis
noiseq_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(NOISeq))

    # Prepare data
    factors <- as.data.frame(coldata$condition)
    colnames(factors) <- c("condition")
    noiseq_data <- readData(data = counts, factors = factors)
    head(noiseq_data)
    # Calculate diff expression
    noiseq_res <- noiseqbio(noiseq_data, factor = "condition")

    # Return results of diff expression analysis
    return(noiseq_res)
}

# Visualize function
noiseq_v <- function(noiseq_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(NOISeq))
    suppressMessages(library(EnhancedVolcano))

    # Built in graphic function
     DE.plot(noiseq_res, q = 0.9, graphic = "MD")


}


# Make a signature of top-n genes
noiseq_top <- function(results, n, out) {
    # Import libraries
    suppressMessages(library("NOISeq"))
    suppressMessages(library("biomaRt"))

    # Filter results by q threshold
    filtered_results <- degenes(results, q = 0.9, M = NULL)

    # Extract top-n differentially expressed genes ordered by prob
    top <- head(filtered_results[order(filtered_results$prob, decreasing = TRUE), ], n)
    tmp <- gsub("\\..*","",row.names(top))
    write.csv(degenes(results, q = 0.0, M = NULL), file = file.path(out, "NOISeq_info.csv"))
    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and probability
noiseq_filtered <- function(results, out) {
  suppressMessages(library("NOISeq"))

  # Filter results by q threshold
  filtered_results <- degenes(results, q = 0.9, M = NULL)
  filtered_genes <- gsub("\\..*","",row.names(filtered_results))
  write.csv(degenes(results, q = 0.0, M = NULL), file = file.path(out, "NOISeq_info.csv"))
  return(filtered_genes)
}
