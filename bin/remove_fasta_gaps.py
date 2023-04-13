#!/usr/bin/env python3

import sys
from Bio import SeqIO

# Get input and output file names from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Replace "-" symbols with "N's" in each sequence record
with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
		
	for seq_record in SeqIO.parse(in_handle, "fasta"):
		seq_record.seq = seq_record.seq.replace("-", "N")
		SeqIO.write(seq_record, out_handle, "fasta")
