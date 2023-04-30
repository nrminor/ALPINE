#!/usr/bin/env python3

import sys
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline

# Define the input and output files
multi_seq_fasta = sys.argv[1]
ref_seq_fasta = sys.argv[2]
label = sys.argv[3]
aligned_output = str(label) + "-centroids-with-ref.fasta"

# Load the pre-aligned multi-sequence alignment and the reference sequence
alignment = AlignIO.read(multi_seq_fasta, "fasta")
ref_seq = SeqIO.read(ref_seq_fasta, "fasta")

# Determine the length of the sequences in the pre-aligned multi-sequence FASTA
seq_length = len(alignment[0])

# Determine the length difference between the pre-aligned multi-sequence FASTA and the reference sequence
ref_length_diff = seq_length - len(ref_seq)

# Add half the difference of 'n' characters to the start and end of the reference sequence
ref_seq_start = "n" * (ref_length_diff // 2)
ref_seq_end = "n" * (ref_length_diff - (ref_length_diff // 2))
ref_seq.seq = ref_seq_start + ref_seq.seq + ref_seq_end

# Add the reference sequence to the alignment
alignment.append(ref_seq)

# Re-align the sequences using MUSCLE
muscle_cline = MuscleCommandline(input=alignment, out=aligned_output)
muscle_cline()

# Replace "/" characters in the deflines with underscores
for record in AlignIO.read(aligned_output, "fasta"):
    record.id = record.id.replace("/", "_")
    record.description = record.description.replace("/", "_")

# Write the aligned sequences to a new file
with open(aligned_output, "w") as handle:
    SeqIO.write(AlignIO.read(aligned_output, "fasta"), handle, "fasta")
