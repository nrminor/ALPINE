#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count, Lock

# Define a function to process a single sequence record
def process_seq_record(seq_record, lock, out_handle):
	seq_str = str(seq_record.seq).replace("-", "N")
	seq_record.seq = Seq(seq_str)
	
	# Acquire the lock to write to the output file
	lock.acquire()
	SeqIO.write(seq_record, out_handle, "fasta")
	lock.release()

# Get input and output file names from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Define the number of processes to use (default to number of CPU cores)
num_processes = int(sys.argv[3])

# Read in the input FASTA file
with open(input_file, "r") as in_handle:
	seq_records = list(SeqIO.parse(in_handle, "fasta"))

# Process the sequence records in parallel
lock = Lock()
with open(output_file, "w") as out_handle:
	with Pool(processes=num_processes) as pool:
		for seq_record in seq_records:
			pool.apply_async(process_seq_record, (seq_record, lock, out_handle))

		# Wait for all the processes to finish before closing the output file
		pool.close()
		pool.join()

