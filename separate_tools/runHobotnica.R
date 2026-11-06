# Logger config
library(logging)
basicConfig()
loginfo('Start')

# Load data
cm_args <- commandArgs(trailingOnly = TRUE)
count <- cm_args[1]
cols <- cm_args[2]

if (!file.exists(count)) {
    stop("File ", count, " does not exists!")
}
if (!file.exists(cols)) {
    stop("File ", cols, " does not exists!")
}

if (length(cm_args) <= 2) {
    stop("Enter output directory name")
}
out <- cm_args[3]
if (!file.exists(out)) {
    dir.create(out)
}
save_flag <- 0
arg_counter <- 4
n <- 30
while (length(cm_args)  >= arg_counter) {
    option <- cm_args[arg_counter]
    if (option == "-s") {
        save_flag <- 1
        arg_counter <- arg_counter + 1
    }
    if (option == "-n") {
        if (length(cm_args)  > arg_counter) {
            arg_counter <- arg_counter + 1
            n <- as.numeric(cm_args[arg_counter])
        } else {
            loginfo("-n option should have a number of genes in top after itself. Default number is 30, run with it")
        }
        arg_counter <- arg_counter + 1
    }
}

loginfo('Start table reading')
# Read count matrix table and correct its form
counts <- read.table(file=count, header = TRUE, sep = ",", dec = ".")
nrow(counts)
counts <- subset(counts, !duplicated(counts$X))
nrow(counts)
gene_names <- counts$X
counts <- counts[,-c(1)]
rownames(counts) <- gene_names

loginfo('Start annotation reading')
# Read annotation and make table with useful form
coldata <- read.table(file=cols, sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]

loginfo('Start comparison of 5 differential analysis tools with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))
loginfo(paste0('Output will be written in ', out))


source("source/signatures_utils.R")
source("source/calculate_distmatrix_utils.R")
loginfo('Visualize signature crossing')
# Draw a Venn diagram
draw_venn_diag(out)
loginfo('Visualized')

loginfo('Save results of distance matrices making')
deseq_distmat <- read.table(file=file.path(out,"DESeq_sig.txt.distmatrix"), header = TRUE, sep = ",")
sample_names <- deseq_distmat$X
deseq_distmat <- deseq_distmat[,-c(1)]
rownames(deseq_distmat) <- sample_names
ebseq_distmat <- read.table(file=file.path(out,"EBSeq_sig.txt.distmatrix"), header = TRUE, sep = ",")
sample_names <- ebseq_distmat$X
ebseq_distmat <- ebseq_distmat[,-c(1)]
rownames(ebseq_distmat) <- sample_names
edger_distmat <- read.table(file=file.path(out,"edgeR_sig.txt.distmatrix"), header = TRUE, sep = ",")
sample_names <- edger_distmat$X
edger_distmat <- edger_distmat[,-c(1)]
rownames(edger_distmat) <- sample_names
voom_distmat <- read.table(file=file.path(out,"voom_sig.txt.distmatrix"), header = TRUE, sep = ",")
sample_names <- voom_distmat$X
voom_distmat <- voom_distmat[,-c(1)]
rownames(voom_distmat) <- sample_names
noiseq_distmat <- read.table(file=file.path(out,"NOISeq_sig.txt.distmatrix"), header = TRUE, sep = ",")
sample_names <- noiseq_distmat$X
noiseq_distmat <- noiseq_distmat[,-c(1)]
rownames(noiseq_distmat) <- sample_names
loginfo('Results of distance matrices making are made')

loginfo('Ready for Hobotnica using')
loginfo('Start scoring Hobotnica')
source("source/hobotnica_using.R")
use_hobotnica(deseq_distmat, ebseq_distmat, edger_distmat, voom_distmat, noiseq_distmat, coldata$condition, out)
loginfo('Hobotnica scores are ready')

h_results <- read.table(file=file.path(out, "hobotnica_scores.txt"), sep = " ", dec = ".")
colnames(h_results) <- c("names", "scores")
h_results <- h_results[order(h_results$scores, decreasing = TRUE), ]
loginfo(paste0('According to Hobotnica best tool is ', h_results$names[1]))

loginfo('Save top genes crossing with best tool\'s top genes')
# Save top genes crossing
best_crossing(h_results, out)
loginfo('Saved')
loginfo('Save top genes crossing')
# Save top genes crossing
crossing(out)
loginfo('Saved')
loginfo('The end')
