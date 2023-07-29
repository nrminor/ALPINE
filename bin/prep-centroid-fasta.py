#!/usr/bin/env python3

import sys
import os
# from numba import njit
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline

# Define the input and output files
if os.path.islink(sys.argv[1]):
    MULTI_SEQ_FASTA = os.readlink(sys.argv[1])
else:
    MULTI_SEQ_FASTA = sys.argv[1]

LABEL = sys.argv[2]

# define main script function
def main(input_fasta: str, input_tag: str):

    """Parse sequences and align them to each other so that they
    are all the same length, as required for calling distance 
    matrices downstream."""

    # Load the multi-sequence FASTA and the reference sequence
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

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

    # perform a muscle alignment with the modified alignment and reference
    align_out_file = f"{input_tag}-aligned-centroids.fasta"
    muscle_cline = MuscleCommandline(input="samelength.fasta",
                                    out=align_out_file,
                                    maxiters=1, quiet=True)
    muscle_cline()

    # Replace "/" characters in the deflines with underscores
    for record in AlignIO.read(align_out_file, "fasta"):
        record.id = record.id.replace("/", "_")
        record.description = record.description.replace("/", "_")

# run main function if not invoked as a module
if __name__ == "__main__":
    main(MULTI_SEQ_FASTA, LABEL)
