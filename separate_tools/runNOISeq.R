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

loginfo('Run differential analysis tool and score with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))
loginfo(paste0('Output will be written in ', out))

# Run tools
# Import functions from modules
source("source/NOISeq.R")
loginfo('Start NOISeq analysis')
noiseq_res <- noiseq_f(counts, coldata)
loginfo('NOISeq analysis is completed')


# Visualize results of differential expression
loginfo('Visualize results of differential expression into one pdf')
png(file.path(out, "NOISeq_de.png"))
noiseq_v(noiseq_res, out)
dev.off()
loginfo('Visualization results of differential expression is done')

loginfo('Make signatures of differential expression analysis')
source("source/signatures_utils.R")
writeLines(noiseq_filtered(noiseq_res, out), file.path(out, "NOISeq_sig.txt"))
loginfo('Signatures of differential expression analysis are made')

pdf(file.path(out, "heatmaps.pdf"))
# Visualize using heat maps
condition <- coldata$condition
loginfo('Calculate heat maps for tools')
source("source/calculate_distmatrix_utils.R")
loginfo('Calculate NOISeq heat map')
noiseq_hm <- heatmap_v (count, "NOISeq", condition, ', q > 0.9', out)
loginfo('NOISeq heat map is done')
loginfo('Heat maps for tools are made')

loginfo('Vizualize heat map')
png(file.path(out, "NOISeq_heatmap.png"))
noiseq_hm
dev.off()
loginfo('Heat map visualization is done')

loginfo('Calculate NOISeq distance matrix')
noiseq_distmat <- calculate_distmatrix (count, file.path(out, "NOISeq_sig.txt"))
loginfo('NOISeq distance matrix is done')

if (save_flag == 1) {
     # Save results of distance matrices
     loginfo('Save results of distance matrices making')
     write.csv(as.data.frame(as.matrix(noiseq_distmat)), file=file.path(out,"NOISeq_sig.txt.distmatrix"))
     loginfo('Results of distance matrices making are made')
}

loginfo('The end')
