library(biomaRt)
library(r2r)
library(sjmisc)

#https://github.com/grimbough/biomaRt/issues/31
httr::set_config(httr::config(ssl_verifypeer = FALSE))

# Load data
cm_args <- commandArgs(trailingOnly = TRUE)
filename <- cm_args[1]
if (!file.exists(filename)) {
    stop("File ", filename, " does not exists!")
}

if (length(cm_args) < 2) {
    stop("Enter output filename name")
}
output_filename <- cm_args[2]

if (length(cm_args) <= 2) {
    stop("Enter data type name")
}
data_type <- cm_args[3]

# Read count matrix table and correct its form
counts <- read.table(file=filename, header = TRUE, sep = ",", dec = ".")
sample_names <- colnames(counts)
nrow(counts)
gene_names <- counts$X
head(gene_names)
gene_names <- sub("\\.[^.]*$", "", gene_names)
counts$X <- gene_names
head(gene_names)

data_to_mgi_df <- data.frame()
if (data_type == "mouse") {
    ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "asia")
    transcript = getBM(attributes=c("ensembl_transcript_id","external_gene_name"),
                 filters="ensembl_transcript_id",
                 values=gene_names,
                 mart=ensembl)
    head(transcript)
    colnames(transcript) <- c("data_symbol", "final_symbol")
    data_to_mgi_df <- rbind(transcript)
    transcript <- NULL
} else if (data_type == "human") {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")
    transcript = getBM(attributes=c("ensembl_transcript_id","external_gene_name"),
                 filters="ensembl_transcript_id",
                 values=gene_names,
                 mart=ensembl)
    head(transcript)
    colnames(transcript) <- c("data_symbol", "final_symbol")
    data_to_mgi_df <- rbind(transcript)
    transcript <- NULL
} else if (data_type == "rat") {
    ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", mirror = "asia")
    transcript = getBM(attributes=c("ensembl_transcript_id","external_gene_name"),
                 filters="ensembl_transcript_id",
                 values=gene_names,
                 mart=ensembl)
    head(transcript)
    colnames(transcript) <- c("data_symbol", "final_symbol")
    data_to_mgi_df <- rbind(transcript)
    transcript <- NULL
} else {
    print("We don't have this species yet to transform it")
}
head(data_to_mgi_df)

data_to_mgi_dict <- hashmap()
for(i in 1:nrow(data_to_mgi_df)) {
  if (!(is.na(data_to_mgi_df[i, "data_symbol"]) || is_empty(data_to_mgi_df[i, "data_symbol"])) ) {
      data_to_mgi_dict[[data_to_mgi_df[i, "data_symbol"]]]<- data_to_mgi_df[i, "final_symbol"]
  }
}
data_to_mgi_df <- NULL

counts_dict <- hashmap()
for(i in 1:nrow(counts)) {
  if (is.na(counts[i, "X"]) || is_empty(counts[i, "X"])) {
        next
  }
  if (has_key(data_to_mgi_dict, counts[i, "X"])) {
    gene_name <- data_to_mgi_dict[[counts[i, "X"]]]
  } else {
    gene_name <- counts[i, "X"]
  }
  if (is.na(gene_name) || is_empty(gene_name)) {
    gene_name <- counts[i, "X"]
  }
  if (has_key(counts_dict, gene_name)) {
    for (j in 1:ncol(counts_dict[[gene_name]])) {
        counts_dict[[gene_name]][1, j] <- counts_dict[[gene_name]][1, j] + counts[i, 1 + j]
    }
  } else {
    counts_dict[[gene_name]] <- counts[i, 2:ncol(counts)]
  }
}
counts <- NULL
write.row <- function(x, file = output_filename, append = TRUE, quote = TRUE, sep=",",
                      eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE,
                      qmethod = c("escape", "double"), fileEncoding = ""){
  write.table(as.matrix(t(x)), file, append , quote, sep,
              eol, na, dec, row.names, col.names,
              qmethod, fileEncoding)
}

write.table(as.matrix(t(sample_names)), output_filename, append = FALSE, row.names = FALSE, col.names = FALSE, sep=",")
for(gene_name in keys(counts_dict)) {
    write.row(c(c(gene_name), counts_dict[[gene_name]]))
}
