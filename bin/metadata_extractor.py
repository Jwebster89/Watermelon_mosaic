#!/usr/bin/env python3
import os, sys
from Bio import SeqIO

def read_genbank(genbank):
    count=0
    for seq_record in SeqIO.parse(genbank, "gb"):
        if count >0:
            for feature in seq_record.features:
                if feature.type == "source":
                    if "country" in feature.qualifiers:
                        print(seq_record.id,feature.qualifiers["country"])
        count=count+1

read_genbank(sys.argv[1])