#!/usr/bin/env python3

import sys
import os
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline

# Define the input and output files
if os.path.islink(sys.argv[1]):
    multi_seq_fasta = os.readlink(sys.argv[1])
else:
    multi_seq_fasta = sys.argv[1]

label = sys.argv[2]

# Load the multi-sequence FASTA and the reference sequence
sequences = list(SeqIO.parse(multi_seq_fasta, "fasta"))

# Find the length of the longest sequence
max_length = max(len(seq) for seq in sequences)

# Loop through the sequences and add "n" symbols as needed
for seq in sequences:
    if len(seq) < max_length:
        diff = max_length - len(seq)
        n_before = diff // 2
        n_after = diff - n_before
        seq.seq = "n" * n_before + seq.seq + "n" * n_after

# Write the updated sequences to a new FASTA file
SeqIO.write(sequences, "samelength.fasta", "fasta")

# import same length FASTA
alignment = AlignIO.read("samelength.fasta", "fasta")

# perform a muscle alignment with the modified alignment and reference
align_out_file = f"{label}-aligned-centroids.fasta"
muscle_cline = MuscleCommandline(input=alignment, out=align_out_file, maxiters=1, quiet=True)
muscle_cline()

# Replace "/" characters in the deflines with underscores
for record in AlignIO.read(align_out_file, "fasta"):
    record.id = record.id.replace("/", "_")
    record.description = record.description.replace("/", "_")
