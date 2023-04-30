#!/usr/bin/env python3

import sys
import fcntl
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count

# Define a function to process a single sequence record
def process_seq_record(seq_record):
	seq_str = str(seq_record.seq).replace("-", "N")
	seq_record.seq = Seq(seq_str)
	return seq_record

# Define a function to write a sequence record to a file with a file lock
def write_seq_record(record, file_handle):
	fcntl.flock(file_handle, fcntl.LOCK_EX)
	SeqIO.write(record, file_handle, "fasta")
	fcntl.flock(file_handle, fcntl.LOCK_UN)

# Get input and output file names from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Define the number of processes to use (default to number of CPU cores)
num_processes = int(sys.argv[3])

# Open the output file with a file lock
with open(output_file, "w") as out_handle:
	fcntl.flock(out_handle, fcntl.LOCK_EX)

	# Read in the input FASTA file
	with open(input_file, "r") as in_handle:
		# Process the sequence records in parallel
		with Pool(processes=num_processes) as pool:
			for processed_seq_record in pool.imap(process_seq_record, SeqIO.parse(in_handle, "fasta")):
				# Write the record to the output file with a file lock
				write_seq_record(processed_seq_record, out_handle)

	fcntl.flock(out_handle, fcntl.LOCK_UN)
