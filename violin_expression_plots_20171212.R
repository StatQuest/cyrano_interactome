library(ggplot2)

library.size.scale.factor <- 100000

##########################
#
#
#  First deal with the protein stuff...
#
#
##########################

counts <- read.delim(file="rna_seq_20131014_htseq_counts_gencode.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

head(counts)
dim(counts)

protein.lengths <- read.delim(file="gencode_protein_gene_lengths.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

head(protein.lengths)
dim(protein.lengths)

protein.counts <- counts[counts$Gene %in% protein.lengths$ID,]
head(protein.counts)
dim(protein.counts)

## To convert to TPM....
## 1) divide by gene length (in kb)
## 2) divide by sequencing depth (in M)

protein.counts <- merge(protein.lengths[,c("ID", "Length")], protein.counts, by.x="ID", by.y="Gene")
head(protein.counts)
## scale length by kb...
protein.counts$Length <- protein.counts$Length / 1000
head(protein.counts)


## verify this worked correctly...
protein.lengths[protein.lengths$ID == "ENSMUSG00000000001.4",]
counts[counts$Gene == "ENSMUSG00000000001.4",]

## now normalize each gene by length...
protein.counts.matrix <- protein.counts[,3:ncol(protein.counts)]
rownames(protein.counts.matrix) <- protein.counts$ID
head(protein.counts.matrix)

protein.rpk.matrix <- protein.counts.matrix / protein.counts$Length
head(protein.rpk.matrix)

## now normalize each sample by library size...
protein.lib.sizes <- colSums(protein.rpk.matrix)
summary(protein.lib.sizes)
protein.lib.sizes <- protein.lib.sizes / library.size.scale.factor
summary(protein.lib.sizes)

protein.tpm.matrix <- protein.rpk.matrix
for (i in 1:length(protein.lib.sizes)) {
  protein.tpm.matrix[,i] <- protein.tpm.matrix[,i] / protein.lib.sizes[i]
}

head(protein.tpm.matrix)
dim(protein.tpm.matrix)
summary(log10(protein.tpm.matrix))


################################
##
## now load in the lncRNA data
##
################################

lncrna.counts <- read.delim(file="rna_seq_20131014_htseq_counts_gencode_lncrna.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

head(lncrna.counts)
dim(lncrna.counts)


lncrna.lengths <- read.delim(file="gencode_lncrna_gene_lengths.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

head(lncrna.lengths)


## To convert to TPM....
## 1) divide by gene length (in kb)
## 2) divide by sequencing depth (in M)

lncrna.counts <- merge(lncrna.lengths[,c("ID", "Length")], lncrna.counts, by.x="ID", by.y="Gene")
head(lncrna.counts)
## scale length by kb...
lncrna.counts$Length <- lncrna.counts$Length / 1000
head(lncrna.counts)


## verify this worked correctly...
lncrna.lengths[lncrna.lengths$ID == "ENSMUSG00000000031.15",]
lncrna.counts[lncrna.counts$ID == "ENSMUSG00000000031.15",]

## now normalize each gene by length...
lncrna.counts.matrix <- lncrna.counts[,3:ncol(lncrna.counts)]
rownames(lncrna.counts.matrix) <- lncrna.counts$ID
head(lncrna.counts.matrix)

lncrna.rpk.matrix <- lncrna.counts.matrix / lncrna.counts$Length
head(lncrna.rpk.matrix)

## now normalize each sample by library size...
lncrna.lib.sizes <- colSums(lncrna.rpk.matrix)
summary(lncrna.lib.sizes)
lncrna.lib.sizes <- lncrna.lib.sizes / library.size.scale.factor
summary(lncrna.lib.sizes)

lncrna.tpm.matrix <- lncrna.rpk.matrix
for (i in 1:length(lncrna.lib.sizes)) {
  lncrna.tpm.matrix[,i] <- lncrna.tpm.matrix[,i] / lncrna.lib.sizes[i]
}

head(lncrna.tpm.matrix)
dim(lncrna.tpm.matrix)
summary(log10(lncrna.tpm.matrix))


##################################
#
#
#  Draw violin plots comparing lncrRNA to protein...
#
#
#################################

protein.tpm.matrix[1:4, 1:3]
lncrna.tpm.matrix[1:5, 1:3]

protein.output <- merge(protein.lengths[, c("ID", "Gene", "Length")],
  protein.tpm.matrix, by.x="ID", by.y="row.names")
head(protein.output)

# write.table(protein.output, file="gencode_protein_genes_tpm.txt", row.names=FALSE, sep="\t", quote=FALSE)

lncrna.output <- merge(lncrna.lengths[, c("ID", "Gene", "Length")],
  lncrna.tpm.matrix, by.x="ID", by.y="row.names")
head(lncrna.output)
# write.table(lncrna.output, file="gencode_lncrna_genes_tpm.txt", row.names=FALSE, sep="\t", quote=FALSE)

colnames(protein.tpm.matrix)
colnames(lncrna.tpm.matrix)

## the column names are the same...
colnames(protein.tpm.matrix) == colnames(lncrna.tpm.matrix)

for (i in 1:ncol(protein.tpm.matrix)) {
  
  # i<-1
  
  sample.name <-  colnames(protein.tpm.matrix)[i]
  
  full.data <- data.frame(
    tpm=c(protein.tpm.matrix[,sample.name],
      lncrna.tpm.matrix[,sample.name]),
    type=c(rep("protein", times=nrow(protein.tpm.matrix)),
      rep("lncRNA", times=nrow(lncrna.tpm.matrix))))
  
  full.data <- full.data[full.data$tpm > 1,]
  
  output.file.name <- paste0(sample.name, "_tpm_plots.pdf")
  print(output.file.name)
  
  
  ggplot(data=full.data, aes(x=type, y=log10(tpm))) +
    geom_violin(aes(fill=type)) +
    geom_boxplot(aes(fill=type), outlier.shape=NA, width=0.25, alpha=0.5) +
    theme_bw() +
    xlab("") +
    scale_fill_discrete("") +
    ggtitle(sample.name)
  ggsave(output.file.name)
}
