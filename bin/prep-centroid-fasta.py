#!/usr/bin/env python3

"""
CLEAN THE ALIGNED CENTROID FASTA FOR DOWNSTREAM USAGE
-----------------------------------------------------

Vsearch --cluster_fast can provide a multi-sequence alignment
as one of its outputs, which we need to call a distance matrix.
But that alignment has a few quirks. It will include cluster
consensus sequences, which we are not interested in (at this
stage, at least). It also adds a pesky "*" to one or more of
the aligned sequences, which is not explained in the vsearch
documentation. We do away with the entire consensus sequences,
and clean the "*" symbols out of any remaining records here.
"""

import argparse
from Bio import SeqIO

def main(input_path: str, label: str, count: int):
    """
    This function parses a FASTA file and removes all records with 
    "consensus" in their defline and then removes an asterisk 
    ("*") symbols in the deflines of the remaining records.

    Args:
    fasta_path: The path to the FASTA file.
    label: A label or prefix for the output file.
    count: An expected count to sanity check with.

    Returns:
    None
    """

    output_filename = f"{label}-centroids.fasta"

    with open(input_path, "r", encoding="utf-8") as infile, open(output_filename, "w", encoding="utf-8") as outfile:
        records = SeqIO.parse(infile, "fasta")
        filtered_records = [record for record in records if "consensus" not in record.description]

        assert len(filtered_records) == count

        for record in filtered_records:
            record.description = record.description.replace("*", "")
            record.id = record.id.replace("*", "")
            SeqIO.write(record, outfile, "fasta")
# end main function def

if __name__ == "__main__":

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_path", help="The path to the FASTA file.")
    parser.add_argument("label", help="Label to use when naming the output file.")
    parser.add_argument("count", type=int, help="An expected count of the records in the FASTA.")
    args = parser.parse_args()

    # run main
    main(args.fasta_path, args.label, args.count)
