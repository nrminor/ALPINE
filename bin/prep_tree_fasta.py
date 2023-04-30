#!/usr/bin/env python3

import sys
import os
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline

# parse command line arguments
if len(sys.argv) != 4:
    sys.exit("Usage: python3 prep_tree_fasta.py <input_alignment.fasta> <reference.fasta> <label>")

# Define the input and output files
multi_seq_fasta = sys.argv[1]
ref_seq_fasta = sys.argv[2]
label = sys.argv[3]

# Load the pre-aligned multi-sequence alignment and the reference sequence
alignment = AlignIO.read(multi_seq_fasta, "fasta")
ref_seq = SeqIO.read(ref_seq_fasta, "fasta")

# determine the length difference between the alignment and reference sequences
align_len = len(alignment[0].seq)
ref_len = len(ref_seq.seq)
diff_len = align_len - ref_len

# add n's to the beginning and end of the reference sequence
half_diff_len = int(diff_len/2)
ref_seq.seq = 'n'*half_diff_len + ref_seq.seq + 'n'*half_diff_len

# write the modified reference sequence to a file
tmp_out_file = f"{label}-tmp.fasta"
with open(tmp_out_file, "w") as f:
    SeqIO.write(ref_seq, f, "fasta")

# write the modified alignment to a file
with open(tmp_out_file, "a") as f:
    SeqIO.write(alignment, f, "fasta")

# perform a muscle alignment with the modified alignment and reference
align_out_file = f"{label}-centroids-with-ref.fasta"
muscle_cline = MuscleCommandline(input=tmp_out_file, out=align_out_file, maxiters=1, quiet=True)
muscle_cline()

# Replace "/" characters in the deflines with underscores
for record in AlignIO.read(align_out_file, "fasta"):
    record.id = record.id.replace("/", "_")
    record.description = record.description.replace("/", "_")

# Remove the temporary file
os.remove(tmp_out_file)