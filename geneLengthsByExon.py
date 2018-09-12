#!/usr/bin/python

from __future__ import print_function
import sys
import argparse

import HTSeq

parser = argparse.ArgumentParser(description='Gene lengths from a GTF file.')

parser.add_argument('-p', '--protein', help='only examine protein coding genes', required=False, action="store_true")
parser.add_argument('-g', '--gtf', help='a GTF file with genome features', required=True)

args=parser.parse_args()

if args.protein:
    print("\nFiltering out all genes by protein coding genes\n", file=sys.stderr)

#
# main()
#            
def main():

    gtf_file = HTSeq.GFF_Reader(args.gtf)
    
    transcript_lengths = {}
    transcript_exon_counts = {}
    gene_name_to_transcripts = {}
    gene_name_to_gene_id = {}
    
    for feature in gtf_file:
        if feature.type == "exon":
            
            gene_type = feature.attr["gene_type"]
            if args.protein:
                if gene_type != "protein_coding":
                    continue
                
            gene_name = feature.attr["gene_name"]
            gene_id = feature.attr["gene_id"]
            transcript_name = feature.attr["transcript_id"]
            exon_length = feature.iv.length
            
            if transcript_name in transcript_lengths:
                transcript_lengths[transcript_name] += exon_length
                transcript_exon_counts[transcript_name] += 1
            else:
                transcript_lengths[transcript_name] = exon_length
                transcript_exon_counts[transcript_name] = 1
                if gene_name in gene_name_to_transcripts:
                    gene_name_to_transcripts[gene_name] += "\t" + transcript_name
                else:
                    gene_name_to_gene_id[gene_name] = gene_id
                    gene_name_to_transcripts[gene_name] = transcript_name

                    
    print("Gene", "ID", "Transcript", "Length", "Exons", sep="\t", file=sys.stdout)
                    
    for gene_name in gene_name_to_transcripts:
#        print("Gene Name:", gene_name, file=sys.stderr)

        transcript_names = gene_name_to_transcripts[gene_name].split("\t")
        longest_transcript_name = ""
        longest_transcript_length = 0
        longest_transcript_exon_count = 0
        
        for transcript_name in transcript_names:
#            print("\t", transcript_name, sep="", file=sys.stderr)
            transcript_length = transcript_lengths[transcript_name]
#            print("\tLength:", transcript_length, file=sys.stderr)
            transcript_exon_count = transcript_exon_counts[transcript_name]
#            print("\tNum Exons:", transcript_exon_count, file=sys.stderr)
            if transcript_length > longest_transcript_length:
                longest_transcript_name = transcript_name
                longest_transcript_length = transcript_length
                longest_transcript_exon_count = transcript_exon_count

        print(gene_name, gene_name_to_gene_id[gene_name], longest_transcript_name, longest_transcript_length, longest_transcript_exon_count, sep="\t", file=sys.stdout)
#        print("Longest transcript:", file=sys.stderr)
#        print("\t", longest_transcript_name, sep="", file=sys.stderr)
#        print("\t", longest_transcript_length, sep="", file=sys.stderr)
#        print("\t", longest_transcript_exon_count, sep="", file=sys.stderr)
                   
#
# end: main()
#            


if __name__ == '__main__': 
    main()
