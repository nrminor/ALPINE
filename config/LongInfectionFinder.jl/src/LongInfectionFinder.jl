module LongInfectionFinder

# load dependencies
using DelimitedFiles, DataFrames, CSV, FastaIO, FileIO, Dates, BioSequences, Distances, Statistics, Pipe, CodecZstd, CodecZlib
import Base.Threads

export filter_tsv_by_geo, filter_by_geo, replace_gaps, filter_by_n, lookup_date, separate_by_month, set_to_uppercase, weight_by_cluster_size, prep_for_clustering

### FUNCTION(S) TO FILTER GENBANK METADATA TO A PARTICULAR GEOGRAPHY STRING ###
### ----------------------------------------------------------------------- ###
function filter_tsv_by_geo(input_tsv::String,geography::String)

    # Read in the TSV file with metadata
    metadata_df = CSV.read(input_tsv, DataFrame, delim="\t")

    # Double check the column name for geographic locations
    if "Geographic location" in names(metadata_df)
        idx = findfirst(isequal("Geographic location"), names(metadata_df))
        rename!(metadata_df, names(metadata_df)[idx] => "Geographic Location")
    end

    # filter metadata based on desired geography
    filtered = metadata_df[[contains(string(value), geography) for value in metadata_df[!,"Geographic Location"]], :]

    # Writing filtered metadata
    CSV.write("filtered-to-geography.tsv", filtered, delim="\t")

    # separating out accessions
    accessions = filtered[!,"Accession"]

    # Writing accessions to a text file for use by seqtk subseq or subseq.rs
    writedlm("accessions.txt", accessions, "\n")
    
end

function filter_by_geo(input_tsv::String,fasta_path::String,geography::String)

    # Read in the TSV file with metadata
    metadata_df = CSV.read(input_tsv, DataFrame, delim="\t")

    # Double check the column name for geographic locations
    if "Geographic location" in names(metadata_df)
        idx = findfirst(isequal("Geographic location"), names(metadata_df))
        rename!(metadata_df, names(metadata_df)[idx] => "Geographic Location")
    end

    # filter metadata based on desired geography
    filtered_meta = metadata_df[[contains(string(value), geography) for value in metadata_df[!,"Geographic Location"]], :]

    # Writing filtered metadata
    CSV.write("filtered-to-geography.tsv", filtered_meta, delim="\t")

    # separating out accessions
    accessions = Set(filtered_meta[!,"Accession"])

    # create a file lock to prevent threads from corrupting the output file
    u = ReentrantLock()

    # create the output file
    filtered_seqs = "filtered-to-geography.fasta.zst"

    # use accessions list to filter FASTA records into a ZSTD
    # compressed FASTA
    open(ZstdCompressorStream, filtered_seqs, "w") do outstream
        FastaWriter(outstream) do fa
            FastaReader(fasta_path) do fr
                # @sync for (name, seq) in fr
                #     Threads.@spawn begin
                #         accession = split(name, " ")[1]
                #         if accession in accessions
                #             Threads.lock(u) do
                #                 writeentry(fa, accession, seq)
                #             end
                #         end
                #     end
                # end
                for (name, seq) in fr
                    accession = split(name, " ")[1]
                    if accession in accessions
                        writeentry(fa, accession, seq)
                    end
                end
            end
        end
    end

    return filtered_seqs

end

### ----------------------------------------------------------------------- ###



### FUNCTION(S) TO REPLACE "-" SYMBOLS WITH "N" CHARACTERS ###
### ------------------------------------------------------ ###
function replace_gaps(fasta_path::String, output_filename::String)

    # create the output file
    touch(output_filename)

    # iterate through each line in the FASTA to minimuze memory usage
    open(ZstdDecompressorStream, fasta_path, "r") do infile
        open(GzipCompressorStream, output_filename, "w") do outfile
            # Read, replace, and write one line at a time
            for line in eachline(infile)
                if startswith(line, ">")
                    accession = split(line, ' ')[1]
                    println(outfile, accession)
                else
                    new_line = replace(line, '-' => 'N')
                    println(outfile, new_line)
                end
            end
        end
    end
end

### ------------------------------------------------------ ###



### FUNCTION(S) TO FILTER A FASTA BY THE NUMBER OF MASKED BASES ("N") IT CONTAINS ###
### ----------------------------------------------------------------------------- ###
function filter_by_n(input_fasta_path::String, output_filename::String)

    touch(output_filename)

    # create a file lock to prevent threads from corrupting the output file
    # u = ReentrantLock()

    # creating a loop that goes through sequence records and writes out
    # any sequences that have less than the minimum N count
    FastaWriter(output_filename , "a") do fa
        FastaReader(input_fasta_path) do fr
            # @sync for (name, seq) in fr
                # Threads.@spawn begin
            max_n_count = floor(length(seq) * 0.1)
            n_count = count("N", convert(String, seq))
            if n_count < max_n_count
                        # Threads.lock(u) do
                writeentry(fa, name, seq)
                        # end
            end
                # end
            # end
        end
    end

end

### ----------------------------------------------------------------------------- ###



### FUNCTION(S) TO SEPARATE FASTA INTO ONE FASTA FOR EACH MONTH IN GENBANK METADATA ### 
### ------------------------------------------------------------------------------- ###
# function that accesses the collection date for each record name
function lookup_date(record_name::String, lookup::Dict)

    # Look up the collection date for the accession number
    date = get(lookup, record_name, "")
    try
        return Date(date)
    catch
        return nothing
    end
end

# function that pulls in metadata to separate each sequence into a FASTA for the month it was collected in
function separate_by_month(input_fasta::String, accession_to_date::Dict)

    # creating a dictionary of FastaWriters for each year_month
    open_writers = Dict{String, FastaWriter}()

    # create a finalizer to close all writers when script ends
    atexit(() -> foreach(close, values(open_writers)))

    # process each sequence in the input FASTA file
    FastaReader(input_fasta) do fr
        for (name, seq) in fr
            record_date = lookup_date(name,accession_to_date)

            if record_date === nothing
                @warn "Skipping record with missing or unparseable date:" name
                continue
            end

            year_month = Dates.format(record_date, "yyyy-mm")

            # get the writer for this year_month, or create a new one if it doesn't exist yet
            writer = get!(open_writers, year_month) do
                output_filename = string(year_month, ".fasta")
                touch(output_filename)
                FastaWriter(output_filename, "a")
            end

            writeentry(writer, name, seq)
        end
    end
end

### ------------------------------------------------------------------------------- ###



### FUNCTION(S) TO COMPUTE WEIGHTED DISTANCE MATRICES FOR GENBANK FASTA FILES ###
### ------------------------------------------------------------------------- ###
# replace lowercase n symbols with uppercase Ns
function set_to_uppercase(fasta_filename::String, temp_filename::String)
    open(fasta_filename) do infile
        open(temp_filename, "w") do outfile
            # Read, replace, and write one line at a time
            for line in eachline(infile)
                if startswith(line, ">")
                    println(outfile, line)
                else
                    line = uppercase(line)
                    println(outfile, line)
                end
            end
        end
    end
end

# determine how many sequences are in each cluster to compute weighted Distances
function weight_by_cluster_size(seq_name::String, cluster_table::DataFrame)

    # filter down to centroid rows only
    centroids = filter(:1 => x -> x == "C", cluster_table)

    # filter down to hits only
    hits = filter(:1 => x -> x != "C", cluster_table)
    if nrow(hits) == 0
        month_total = 1
    else
        month_total = nrow(hits)
    end

    # find the number of rows for the accession in question
    cluster_size = nrow(filter(:9 => x -> x == seq_name, centroids))

    # return the weight for this centroid's distance
    weight = cluster_size / month_total
    return weight

end

# compute cluster-size-weighted weighted distance matrices
function distance_matrix(temp_filename::String, cluster_table::DataFrame, count::Int, majority_centroid::String, yearmonth::String)

    # Collect both names and sequences
    seqs = [seq for (name, seq) in FastaReader(temp_filename)]
    seq_names = [name for (name, seq) in FastaReader(temp_filename)]

    # Convert the sequences to BioSequence objects, replacing 'N' characters
    filtered_seqs = [replace(seq, 'N' => '-') for seq in seqs]
    seq_vectors = [LongSequence{DNAAlphabet{4}}(seq) for seq in filtered_seqs]

    # Compute the Hamming distance matrix
    dist_matrix = pairwise(Hamming(), seq_vectors, seq_vectors)

    # Convert the distance matrix to a DataFrame
    dist_df = DataFrame(dist_matrix, :auto)
    rename!(dist_df, seq_names)

    # weight distance estimates by relative cluster size
    for seq in names(dist_df)
        col_weight = weight_by_cluster_size(seq, cluster_table)
        dist_df[!, seq] = dist_df[!, seq] * col_weight
    end

    # Add a column for the sequence names
    dist_df[!, :Sequence_Name] = seq_names

    # Move the Sequence_Name column to the front
    select!(dist_df, :Sequence_Name, :)

    if count > 2

        # Write the distance matrix to a CSV file
        CSV.write("$yearmonth-dist-matrix.csv", dist_df)
    
    else

        # constrain down to a simple one by one with the non-majority cluster 
        # as the only row
        filter!(:Sequence_Name => !=(majority_centroid), dist_df)
        select!(dist_df, Not(Symbol(dist_df[1,:1])))

        # Write the distance matrix to a CSV file
        CSV.write("$yearmonth-dist-matrix.csv", dist_df)

    end

end

### ------------------------------------------------------------------------- ###



### FUNCTION THAT IS A COMBINATION OF THE ABOVE FUNCTIONS TO COMPLETE ALL DESIRED OPERATIONS ### 
### ---------------------------------------------------------------------------------------- ###
function prep_for_clustering(ncbi_metadata::String, desired_geography::String, ncbi_fasta::String)

    # ------------------------------------------------------------------------
    # Note that the custom functions below are defined and precompiled in the
    # LongInfectionFinder.jl module created for this project. These functions
    # are best run in the context of the project Docker container, the Dockerfile 
    # for which is provided and built in the config subdirectory of the project
    # root directory.
    # ------------------------------------------------------------------------
    
    # filter metadata and fasta to determine which accessions will be retained downstream
    fasta = filter_by_geo(ncbi_metadata,ncbi_fasta,desired_geography)

    # Quickly go through the FASTA and replace all "-" characters with "N" characters
    replace_gaps(fasta, "no-gaps.fasta.gz")

    # Now that "-" characters have been substituted for "N" characters, go through all
    # the sequences from the geographic region of interest, and remove any records 
    # with sequences whose bases are >10% N
    filter_by_n("no-gaps.fasta.gz", "filtered-by-n.fasta.gz")

    # Finally, we figure out which months each sequence comes from, and separate them 
    # out into one FASTA per month (or more specifically, per year-month)
    metadata_df = CSV.read("filtered-to-geography.tsv", DataFrame, delim="\t")
    accession_to_date = Dict(zip(metadata_df[!,"Accession"], metadata_df[!,"Isolate Collection date"]))
    separate_by_month("filtered-by-n.fasta.gz", accession_to_date)

end

### ---------------------------------------------------------------------------------------- ###

end # LongInfectionFinder module