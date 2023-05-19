#!/usr/bin/env julia

# saving command line arguments supplied by nextflow
fasta_path = ARGS[1]
output_filename = ARGS[2]

# Check if the input FASTA file is a symlink, and if it is, follow the symlink
if islink(fasta_path)
    fasta_path = readlink(fasta_path)
end

# create the output file
touch(output_filename)

# iterate through each line in the FASTA to minimuze memory usage
open(fasta_path) do infile
    open(output_filename, "w") do outfile
        # Read, replace, and write one line at a time
        for line in eachline(infile)
            if startswith(line, ">")
                println(outfile, line)
            else
                new_line = replace(line, '-' => 'N')
                println(outfile, new_line)
            end
        end
    end
end