#!/usr/bin/env python3

import os,sys,argparse
from Bio import SeqIO

def fasta_parse(fasta,min, max):
    subset=[]
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if len(seq_record.seq) > int(min) and len(seq_record.seq) < int(max):
            subset.append(seq_record)
    SeqIO.write(subset, f'subset_{min}-{max}.fasta', 'fasta')

def main():
    fasta_parse(sys.argv[1],sys.argv[2], sys.argv[3])


if __name__ == '__main__':
    main()
