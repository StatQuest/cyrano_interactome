#!/usr/bin/python

from __future__ import print_function
import sys
import argparse

import HTSeq

parser = argparse.ArgumentParser(description='Extra gene lengths from a GTF file.')

parser.add_argument('-p', '--protein', help='only examine protein coding genes', required=False, action="store_true")
parser.add_argument('-g', '--gtf', help='a GTF file with genome features', required=True)
parser.add_argument('-b', '--bed', help='output in BED format', required=False, action="store_true")

args=parser.parse_args()

if args.protein:
    print("\nFiltering out everything but protein coding genes\n", file=sys.stderr)

if args.bed:
    print("\nOutput will be in BED format\n", file=sys.stderr)

#
# main()
#            
def main():

    gtf_file = HTSeq.GFF_Reader(args.gtf)
    transcript_lengths = {}
    starts = {}
    stops = {}
    chrs  = {}
    strands = {}
    
    for feature in gtf_file:
        if feature.type == "transcript":
            #gene_id = feature.attr["gene_id"]
            gene_type = feature.attr["gene_type"]
            if args.protein:
                if gene_type != "protein_coding":
                    continue
            gene_name = feature.attr["gene_name"]
            if gene_name in transcript_lengths:
                old_length = transcript_lengths[gene_name]
                if feature.iv.length > old_length:
                    transcript_lengths[gene_name] = feature.iv.length
                    starts[gene_name]  = feature.iv.start
                    stops[gene_name]   = feature.iv.end
                    chrs[gene_name]    = feature.iv.chrom
                    strands[gene_name] = feature.iv.strand
            else:
                transcript_lengths[gene_name] = feature.iv.length
                starts[gene_name]  = feature.iv.start
                stops[gene_name]   = feature.iv.end
                chrs[gene_name]    = feature.iv.chrom
                strands[gene_name] = feature.iv.strand

    if not args.bed:
        print("gene_name\tLength", sep="", file=sys.stdout)

        
    for gene_name in transcript_lengths:

        
        
        if args.bed:
            print(chrs[gene_name], starts[gene_name], stops[gene_name], gene_name, "1000", strands[gene_name], transcript_lengths[gene_name], sep="\t", file=sys.stdout)
            #print(chrs[gene_name], starts[gene_name], stops[gene_name], gene_name, "1000", strands[gene_name], transcript_lengths[gene_name][0], sep="\t", file=sys.stdout)
            
        else:
            print(gene_name, "\t", transcript_lengths[gene_name], sep="", file=sys.stdout)
#
# end: main()
#            


if __name__ == '__main__': 
    main()
