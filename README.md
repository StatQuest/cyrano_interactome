The following commands will re-create the graphs in:

"A Multimodal Long Noncoding RNA-Centered Network to Support Self-Renewal and
Pluripotency"

By Smith, Starmer and Magnuson

## NOTE: These scripts depend on...

Python3:
HTSeq

R:
ggplot2

bedtools: http://bedtools.readthedocs.io/en/latest/index.html

####################
#
# Figure1, A, B, C and D
#
####################

#####
Calculate the number of exons and the lengths of the protein coding genes
and long non-coding RNAs...
#####

## Human

./transcriptLengthsByExon.py -p -g gencode.human.v27.basic.annotation.gtf > gencode.human.v27.max_transcript_lengths.txt

./transcriptLengthsByExon.py -g gencode.human.v27.long_noncoding_RNAs.gtf > gencode.human.v27.long_noncoding_RNAs.max_transcript_lengths.txt

## Mouse

./transcriptLengthsByExon.py -p -g gencode.mouse.vM15.basic.annotation.gtf > gencode.mouse.vM15.max_transcript_lengths.txt

./transcriptLengthsByExon.py -g gencode.mouse.vM15.long_noncoding_RNAs.gtf > gencode.mouse.vM15.long_noncoding_RNAs.max_transcript_lengths.txt

#####
Verify the lengths of Cyrano (aka OIP5-AS1 in human and 1700020I14Rik in mouse)
#####

## Human:OIP5-AS1	ENST00000500949.6	8829	4

grep OIP5-AS1 gencode.human.v27.long_noncoding_RNAs.max_transcript_lengths.txt

## Mouse: 1700020I14Rik	ENSMUST00000153581.1	8300	3

grep 1700020I14Rik gencode.mouse.vM15.long_noncoding_RNAs.max_transcript_lengths.txt

#####
Now create gene length and number of exon histograms and p-values with R
#####

Rscript protein_lncRNA_histigorams.R





####################
#
# Figure1, E and F
#
####################

#####
Now determine the longest exon for all protein coding genes and long, non-coding
RNAs.
#####

## Human

./longestExon.py -p -g gencode.human.v27.basic.annotation.gtf > gencode.human.v27.longest_exon.txt

./longestExon.py -g gencode.human.v27.long_noncoding_RNAs.gtf > gencode.human.v27.long_noncoding_RNAs.longest_exon.txt

## Mouse

./longestExon.py -p -g gencode.mouse.vM15.basic.annotation.gtf > gencode.mouse.vM15.longest_exon.txt

./longestExon.py -g gencode.mouse.vM15.long_noncoding_RNAs.gtf > gencode.mouse.vM15.long_noncoding_RNAs.longest_exon.txt

#####
Verify the lengths of Cyrano (aka OIP5-AS1 in human and 1700020I14Rik in mouse)
#####

## Human:OIP5-AS1	8494

grep OIP5-AS1 gencode.human.v27.long_noncoding_RNAs.longest_exon.txt

## Mouse:1700020I14Rik	7915

grep 1700020I14Rik gencode.mouse.vM15.long_noncoding_RNAs.longest_exon.txt 


#####
Now create largest exon histograms and p-values with R
#####

Rscript protein_lncRNA_longest_exon_histograms.R




####################
#
# Figure3, A
#
####################

#####
Calculate the total length of exons for protein coding genes and
long non-coding RNAs... NOTE: We are using an older version of the
gencode annotations here because we need it to calculate TPM for an
older RNA-seq experiment...
#####

./geneLengthsByExon.py -p -g gencode.vM13.annotation.gtf > gencode_protein_gene_lengths.txt

./geneLengthsByExon.py -g gencode.vM13.long_noncoding_RNAs.gtf > gencode_lncrna_gene_lengths.txt


#####
Now create violin plots with R
#####

Rscript violin_expression_plots_20171212.R




####################
#
# Figure2, D and E
#
####################

#####
Now download CLIPdb and eCLIP data sets (these give the locations of
protein binding sites) for HUMAN

NOTE: to simplify this process, rather than downloading the orignal files,
I have included the "liftover" output that has the correct hg38 coordinates.
#####

# CLIPdb data

#bigBedToBed http://lulab.life.tsinghua.edu.cn/postar/browse/bigbed/human_clipdb_par2.table.bed.sorted.bb hg19_clipdb.bed -udcDir=.

# NOTE: These are hg19 coordinates, and we need hg38 coordinates. So we convert
# these coordinates using LiftOver:
#
# https://genome.ucsc.edu/cgi-bin/hgLiftOver
#
# To download the new, hg38, coordinates, click on the "View Conversions" link
# and save the file as: hg19_liftover_hg38_clipdb.bed

# eCLIP data

#bigBedToBed http://lulab.life.tsinghua.edu.cn/postar/browse/bigbed/human_eclip.table.bed.sorted.bb hg19_eclip.bed -udcDir=.

# NOTE: These are hg19 coordinates and need to be converted, via LiftOver, to
# hg38. Save the transformed file as: hg19_liftover_hg38_eclip.bed


#####
Now get the longest transcript for all protein coding genes and put the
coordinates in a bed file...
#####

./longestTranscript.py -b -p -g gencode.human.v27.basic.annotation.gtf > gencode.human.v27.longest_transcripts.bed 

./longestTranscript.py -b -g gencode.human.v27.long_noncoding_RNAs.gtf > gencode.human.v27.long_noncoding_RNA_longest_transcript.bed

#####
Now use bedtools to find the intersections between the binding sites and the
transcripts
#####

# CLIPdb

bedtools intersect -a gencode.human.v27.longest_transcripts.bed -b hg19_liftover_hg38_clipdb.bed -c > gencode.human.protein.clipb.overlaps.bed

bedtools intersect -a gencode.human.v27.long_noncoding_RNA_longest_transcript.bed -b hg19_liftover_hg38_clipdb.bed -c > gencode.human.lncRNA_binding.clipdb.overlaps.bed

# eCLIP

bedtools intersect -a gencode.human.v27.longest_transcripts.bed -b hg19_liftover_hg38_eclip.bed -c > gencode.human.protein.eclip.overlaps.bed

bedtools intersect -a gencode.human.v27.long_noncoding_RNA_longest_transcript.bed -b hg19_liftover_hg38_eclip.bed -c > gencode.human.lncrna.eclip.overlaps.bed


#####
Now draw histograms and generate p-values with R
#####

Rscript protein_lncRNA_binding_sites_histograms.R 
