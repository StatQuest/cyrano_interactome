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
    print("Filtering out all GTF entries except those for protein coding", file=sys.stderr)

#
# main()
#            
def main():

    gtf_file = HTSeq.GFF_Reader(args.gtf)

    gene_name_to_longest_exon = {}
        
    for feature in gtf_file:
        if feature.type == "exon":
            
            gene_type = feature.attr["gene_type"]
            if args.protein:
                if gene_type != "protein_coding":
                    continue
                
            gene_name = feature.attr["gene_name"]
            exon_length = feature.iv.length
            
            if gene_name in gene_name_to_longest_exon:
                if exon_length > gene_name_to_longest_exon[gene_name]:
                    gene_name_to_longest_exon[gene_name] = exon_length
            else:
                gene_name_to_longest_exon[gene_name] = exon_length

                    
    print("Gene", "LongestExonLength", sep="\t", file=sys.stdout)
    for gene_name in gene_name_to_longest_exon:
        print(gene_name, gene_name_to_longest_exon[gene_name], sep="\t", file=sys.stdout)
                   
#
# end: main()
#            


if __name__ == '__main__': 
    main()
