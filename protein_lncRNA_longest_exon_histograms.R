library(ggplot2)

bin.width <- 250
the.lwd <- 0.5

#################################
##
## first, human.
##
#################################

gene.lengths <- read.delim(file="gencode.human.v27.longest_exon.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

gene.lengths$type <- "Protein Coding"
colnames(gene.lengths)[2] <- "Length"
head(gene.lengths)

summary(gene.lengths$Length)
gene.lengths[gene.lengths$Length == max(gene.lengths$Length),]

lncrna.lengths <- read.delim(file="gencode.human.v27.long_noncoding_RNAs.longest_exon.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

lncrna.lengths$type <- "lncRNA"
colnames(lncrna.lengths)[2] <- "Length"
head(lncrna.lengths)

summary(lncrna.lengths$Length)
lncrna.lengths[lncrna.lengths$Length == max(lncrna.lengths$Length),]

## OIP5-AS1 is the name for Cyrano lnc-RNA in the human genome
lncrna.lengths[lncrna.lengths$Gene == "OIP5-AS1",]$Length

## remove lncRNAs from gene.lengths
dim(gene.lengths)
gene.lengths <- gene.lengths[!(gene.lengths$Gene %in% lncrna.lengths$Gene),]
dim(gene.lengths)

cyrano.length <- lncrna.lengths[lncrna.lengths$Gene == "OIP5-AS1",]$Length
cyrano.length

full.data <- rbind(gene.lengths, lncrna.lengths)
head(full.data)
tail(full.data)

dim(full.data)
## remove genes with length = 0
full.data <- full.data[full.data$Length > 0,]
dim(full.data)

cyrano.pvalue <- 1 - (nrow(full.data[full.data$type == "lncRNA" & 
    full.data$Length <= 
    full.data[full.data$Gene == "OIP5-AS1",]$Length, ]) / 
    nrow(full.data[full.data$type == "lncRNA",]))

cyrano.pvalue.full <- 1 - (nrow(full.data[full.data$Length <= 
    full.data[full.data$Gene == "OIP5-AS1",]$Length, ]) / nrow(full.data))


summary(full.data$Length)
summary(full.data[full.data$type == "Protein Coding",]$Length)
summary(full.data[full.data$type == "lncRNA",]$Length)

print(paste0("OIP5-AS1 longest exon length: ", full.data[full.data$Gene == "OIP5-AS1",]$Length))

print(paste0("OIP5-AS1 vs human lncRNA longest exon length p-value: ", round(cyrano.pvalue, digits=4)))
print(paste0("OIP5-AS1 vs human lncRNA + protein longest exon length p-value: ", round(cyrano.pvalue.full, digits=4)))



# full.data$type <- factor(full.data$type, levels=c("Protein Coding", "lncRNA")) 

ggplot(data=full.data, 
  aes(x=Length, ..density.., color=type)) +
  # geom_freqpoly(binwidth=bin.width, lwd=the.lwd) +
  geom_freqpoly(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  xlim(c(0, 10000)) +
  ylab("Density") +
  xlab("Largest Exon Size (nt)") +
  scale_color_discrete(name="")
ggsave(file="human_protein_lncrna_longest_exon_freqpoly.pdf", width=6, height=4)

ggplot(data=full.data, 
  aes(x=Length, ..density.., fill=type)) +
  geom_histogram(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  xlim(c(0, 10000)) +
  ylab("Density") +
  xlab("Largest Exon Size (nt)") +
  scale_fill_discrete(name="") 
ggsave(file="human_protein_lncrna_longest_exon_histogram.pdf", width=6, height=4)


#################################
##
## second, mouse.
##
#################################

gene.lengths <- read.delim(file="gencode.mouse.vM15.longest_exon.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

gene.lengths$type <- "Protein Coding"
colnames(gene.lengths)[2] <- "Length"
head(gene.lengths)

summary(gene.lengths$Length)
gene.lengths[gene.lengths$Length == max(gene.lengths$Length),]

lncrna.lengths <- read.delim(file="gencode.mouse.vM15.long_noncoding_RNAs.longest_exon.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

lncrna.lengths$type <- "lncRNA"
colnames(lncrna.lengths)[2] <- "Length"
head(lncrna.lengths)

summary(lncrna.lengths$Length)
lncrna.lengths[lncrna.lengths$Length == max(lncrna.lengths$Length),]

## remove lncRNAs from gene.lengths
dim(gene.lengths)
gene.lengths <- gene.lengths[!(gene.lengths$Gene %in% lncrna.lengths$Gene),]
dim(gene.lengths)

## 1700020I14Rik is the name for Cyrano lnc-RNA in the mouse genome
lncrna.lengths[lncrna.lengths$Gene == "1700020I14Rik",]$Length

cyrano.length <- lncrna.lengths[lncrna.lengths$Gene == "1700020I14Rik",]$Length

full.data <- rbind(gene.lengths, lncrna.lengths)
head(full.data)
tail(full.data)

dim(full.data)
## remove genes with length = 0
full.data <- full.data[full.data$Length > 0,]
dim(full.data)

cyrano.pvalue <- 1 - (nrow(full.data[full.data$type == "lncRNA" & 
    full.data$Length <= 
    full.data[full.data$Gene == "1700020I14Rik",]$Length, ]) / 
    nrow(full.data[full.data$type == "lncRNA",]))

cyrano.pvalue.full <- 1 - (nrow(full.data[full.data$Length <= 
    full.data[full.data$Gene == "1700020I14Rik",]$Length, ]) / nrow(full.data))


summary(full.data$Length)
summary(full.data[full.data$type == "Protein Coding",]$Length)
summary(full.data[full.data$type == "lncRNA",]$Length)

print(paste0("1700020I14Rik longest exon length: ", full.data[full.data$Gene == "1700020I14Rik",]$Length))

print(paste0("1700020I14Rik vs mouse lncRNA longest exon length p-value: ", round(cyrano.pvalue, digits=4)))
print(paste0("1700020I14Rik vs mouse lncRNA + protein longest exon length p-value: ", round(cyrano.pvalue.full, digits=4)))


ggplot(data=full.data, 
  aes(x=Length, ..density.., color=type)) +
  geom_freqpoly(binwidth=bin.width, lwd=the.lwd) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  xlim(c(0, 10000)) +
  ylab("Density") +
  xlab("Largest Exon Size (nt)") +
  scale_color_discrete(name="")
ggsave(file="mouse_protein_lncrna_longest_exon_freqpoly.pdf", width=6, height=4)

ggplot(data=full.data, 
  aes(x=Length, ..density.., fill=type)) +
  geom_histogram(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  xlim(c(0, 10000)) +
  ylab("Density") +
  xlab("Largest Exon Size (nt)") +
  scale_fill_discrete(name="") 
ggsave(file="mouse_protein_lncrna_longest_exon_histogram.pdf", width=6, height=4)

