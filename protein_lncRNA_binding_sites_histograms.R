library(ggplot2)

bin.width <- 1
the.lwd <- 0.5

#################################
##
## first, human.
##
#################################

protein.data <- read.delim(file="gencode.human.protein.clipb.overlaps.bed", sep="\t", stringsAsFactors=FALSE, header=FALSE)

head(protein.data)
colnames(protein.data) <- c("chr", "start", "stop", "gene", "blah", "strand", "length", "count")
protein.data$type <- "Protein Coding"
head(protein.data)

summary(protein.data$count)

lncrna.data <- read.delim(file="gencode.human.lncRNA_binding.clipdb.overlaps.bed", sep="\t", stringsAsFactors=FALSE, header=FALSE)

head(lncrna.data)
colnames(lncrna.data) <- c("chr", "start", "stop", "gene", "blah", "strand", "length", "count")
lncrna.data$type <- "lncRNA"
head(lncrna.data)

summary(lncrna.data$count)

lncrna.data[lncrna.data$gene == "OIP5-AS1",]$count
lncrna.data[lncrna.data$gene == "OIP5-AS1",]


## remove lncRNAs from protein.data
dim(protein.data)
protein.data <- protein.data[!(protein.data$gene %in% lncrna.data$gene),]
dim(protein.data)

full.data <- rbind(protein.data, lncrna.data)
head(full.data)
tail(full.data)

dim(full.data)
## remove genes with length = 0
full.data <- full.data[full.data$length > 0,]
dim(full.data)

## scale the counts per kb of gene
full.data$countPerKb <- full.data$count / (full.data$length / 1000)
head(full.data)
tail(full.data)

summary(full.data$countPerKb)
summary(full.data[full.data$type == "Protein Coding",]$countPerKb)
summary(full.data[full.data$type == "lncRNA",]$countPerKb)

full.data[full.data$gene == "OIP5-AS1",]

## calculate p-values...
## p-value for just lncRNAs...
cyrano.pvalue <- 1 - (nrow(full.data[full.data$type == "lncRNA" & 
    full.data$countPerKb <= 
      full.data[full.data$gene == "OIP5-AS1",]$countPerKb, ]) / 
  nrow(full.data[full.data$type == "lncRNA",]))

cyrano.pvalue.full <- 1 - (nrow(full.data[full.data$countPerKb <= 
      full.data[full.data$gene == "OIP5-AS1",]$countPerKb, ]) / nrow(full.data))


print(paste0("OIP5-AS1 clipdb count: ", full.data[full.data$gene == "OIP5-AS1",]$count))
print(paste0("OIP5-AS1 length: ", full.data[full.data$gene == "OIP5-AS1",]$length))
print(paste0("OIP5-AS1 countPerKb: ", round(full.data[full.data$gene == "OIP5-AS1",]$countPerKb, digits=4)))

print(paste0("OIP5-AS1 vs human lncRNA clipdb binding sites p-value: ", round(cyrano.pvalue, digits=4)))

print(paste0("OIP5-AS1 vs human lncRNA + protein clipdb binding sites p-value: ", round(cyrano.pvalue.full, digits=4)))


# full.data$type <- factor(full.data$type, levels=c("Protein Coding", "lncRNA")) 

ggplot(data=full.data, 
  aes(x=countPerKb, ..density.., color=type)) +
  # geom_freqpoly(binwidth=bin.width, lwd=the.lwd) +
  geom_freqpoly(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  xlim(c(0, 20)) +
  ylab("Density") +
  xlab("Binding sites per Kb") +
  scale_color_discrete(name="")
ggsave(file="human_protein_lncrna_clipdb_overlap.pdf", width=6, height=4)

ggplot(data=full.data, 
  aes(x=countPerKb, ..density.., fill=type)) +
  geom_histogram(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  ylim(c(0, 0.5)) +
  xlim(c(0, 20)) +
  ylab("Density") +
  xlab("Binding sites per Kb") +
  scale_fill_discrete(name="") 
ggsave(file="human_protein_lncrna_clipdb_overlap_histogram.pdf", width=6, height=4)

#####
#
# Now for eclip...
#
######

protein.data <- read.delim(file="gencode.human.protein.eclip.overlaps.bed", sep="\t", stringsAsFactors=FALSE, header=FALSE)

head(protein.data)
colnames(protein.data) <- c("chr", "start", "stop", "gene", "blah", "strand", "length", "count")
protein.data$type <- "Protein Coding"
head(protein.data)

summary(protein.data$count)

lncrna.data <- read.delim(file="gencode.human.lncrna.eclip.overlaps.bed", sep="\t", stringsAsFactors=FALSE, header=FALSE)

head(lncrna.data)
colnames(lncrna.data) <- c("chr", "start", "stop", "gene", "blah", "strand", "length", "count")
lncrna.data$type <- "lncRNA"
head(lncrna.data)

summary(lncrna.data$count)

lncrna.data[lncrna.data$gene == "OIP5-AS1",]$count
lncrna.data[lncrna.data$gene == "OIP5-AS1",]


## remove lncRNAs from protein.data
dim(protein.data)
protein.data <- protein.data[!(protein.data$gene %in% lncrna.data$gene),]
dim(protein.data)

full.data <- rbind(protein.data, lncrna.data)
head(full.data)
tail(full.data)

dim(full.data)
## remove genes with length = 0
full.data <- full.data[full.data$length > 0,]
dim(full.data)

## scale the counts per kb of gene
full.data$countPerKb <- full.data$count / (full.data$length / 1000)
head(full.data)
tail(full.data)

summary(full.data$countPerKb)
summary(full.data[full.data$type == "Protein Coding",]$countPerKb)
summary(full.data[full.data$type == "lncRNA",]$countPerKb)

full.data[full.data$gene == "OIP5-AS1",]

## calculate p-values...
## p-value for just lncRNAs...
cyrano.pvalue <- 1 - (nrow(full.data[full.data$type == "lncRNA" & 
    full.data$countPerKb <= 
    full.data[full.data$gene == "OIP5-AS1",]$countPerKb, ]) / 
    nrow(full.data[full.data$type == "lncRNA",]))

cyrano.pvalue.full <- 1 - (nrow(full.data[full.data$countPerKb <= 
    full.data[full.data$gene == "OIP5-AS1",]$countPerKb, ]) / nrow(full.data))


print(paste0("OIP5-AS1 eclip count: ", full.data[full.data$gene == "OIP5-AS1",]$count))
print(paste0("OIP5-AS1 length: ", full.data[full.data$gene == "OIP5-AS1",]$length))
print(paste0("OIP5-AS1 countPerKb: ", round(full.data[full.data$gene == "OIP5-AS1",]$countPerKb, digits=4)))

print(paste0("OIP5-AS1 vs human lncRNA eclip binding sites p-value: ", round(cyrano.pvalue, digits=4)))

print(paste0("OIP5-AS1 vs human lncRNA + protein eclip binding sites p-value: ", round(cyrano.pvalue.full, digits=4)))

# full.data$type <- factor(full.data$type, levels=c("Protein Coding", "lncRNA")) 

ggplot(data=full.data, 
  aes(x=countPerKb, ..density.., color=type)) +
  # geom_freqpoly(binwidth=bin.width, lwd=the.lwd) +
  geom_freqpoly(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  xlim(c(0, 75)) +
  ylab("Density") +
  xlab("Binding sites per Kb") +
  scale_color_discrete(name="")
ggsave(file="human_protein_lncrna_eclip_overlap.pdf", width=6, height=4)

ggplot(data=full.data, 
  aes(x=countPerKb, ..density.., fill=type)) +
  geom_histogram(binwidth=bin.width) +
  # geom_vline(xintercept=cyrano.length, color="blue") +
  theme_bw() +
  ylim(c(0, 0.5)) +
  xlim(c(0, 75)) +
  ylab("Density") +
  xlab("Binding sites per Kb") +
  scale_fill_discrete(name="") 
ggsave(file="human_protein_lncrna_eclip_overlap_histogram.pdf", width=6, height=4)
