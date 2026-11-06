
# Calculate hobotnica score for one tool
calculate_hobotnica <- function(distmat, condition) {
    source("hobotnica/R/Hobotnica.R")
    return(Hobotnica(distmat, condition))
}

# Calculate hobotnica scores for all tools and write it in the files
use_hobotnica <- function(deseq_distmat, ebseq_distmat, edger_distmat, voom_distmat, noiseq_distmat, condition, out) {
    results <- c()
    results <- c(results, paste0("DESeq ", calculate_hobotnica(deseq_distmat, condition)))
    results <- c(results, paste0("EBSeq ", calculate_hobotnica(ebseq_distmat, condition)))
    results <- c(results, paste0("edgeR ", calculate_hobotnica(edger_distmat, condition)))
    results <- c(results, paste0("voom ", calculate_hobotnica(voom_distmat, condition)))
    results <- c(results, paste0("NOISeq ", calculate_hobotnica(noiseq_distmat, condition)))
    writeLines(results, file.path(out, "hobotnica_scores.txt"))
    return(results)
}
