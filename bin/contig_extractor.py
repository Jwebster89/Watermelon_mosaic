#!/usr/bin/env python3

import os,sys,argparse
from Bio import SeqIO

def fasta_parse(fasta,substring):
    subset=[]
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if substring in seq_record.id :
            subset.append(seq_record)
    SeqIO.write(subset, f'extracted_contig_{fasta}', 'fasta')

def main():
    fasta_parse(sys.argv[1],sys.argv[2])


if __name__ == '__main__':
    main()
