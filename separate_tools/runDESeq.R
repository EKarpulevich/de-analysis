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

loginfo('Rundifferential analysis tool and score with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))
loginfo(paste0('Output will be written in ', out))

# Run tools
# Import functions from modules
source("source/DESeq.R")
loginfo('Start DESeq analysis')
deseq2_res <- deseq2_f(counts, coldata)
loginfo('DESeq analysis is completed')

# Visualize results of differential expression
loginfo('Visualize results of differential expression into one pdf')
deseq_de <- deseq2_v(deseq2_res, out)
loginfo('Vizualize results of differential expression into png')
png(file.path(out, "DESeq_de.png"))
deseq_de
dev.off()
loginfo('Visualization results of differential expression is done')

loginfo('Make signatures of differential expression analysis')
source("source/signatures_utils.R")
writeLines(deseq2_filtered(deseq2_res, out), file.path(out, "DESeq_sig.txt"))
loginfo('Signatures of differential expression analysis are made')

pdf(file.path(out, "heatmaps.pdf"))
# Visualize using heat maps
condition <- coldata$condition
loginfo('Calculate heat maps for tools')
source("source/calculate_distmatrix_utils.R")
loginfo('Calculate DESeq heat map')
deseq_hm <- heatmap_v (count, "DESeq", condition, ', |log2FC| > 1, p < 0.05', out)
loginfo('DESeq heat map is done')
dev.off()

loginfo('Vizualize heat map')
png(file.path(out, "DESeq_heatmap.png"))
deseq_hm
dev.off()
loginfo('Heat map visualization is done')


loginfo('Calculate DESeq distance matrix')
deseq_distmat <- calculate_distmatrix (count, file.path(out, "DESeq_sig.txt"))
loginfo('DESeq distance matrix is done')

if (save_flag == 1) {
     # Save results of distance matrices
     loginfo('Save results of distance matrices making')
     write.csv(as.data.frame(as.matrix(deseq_distmat)), file=file.path(out,"DESeq_sig.txt.distmatrix"))
     loginfo('Results of distance matrices making are made')
}

loginfo('The end')
