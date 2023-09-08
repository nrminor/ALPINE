module ALPINE

# load dependencies
using DelimitedFiles, DataFrames, CSV, Arrow, FastaIO, FileIO, Dates, BioSequences, Distances, Statistics, Pipe, CodecZstd, CodecZlib, FLoops, Missings, RCall, Scratch
import Base.Threads

export filter_metadata_by_geo, filter_by_geo, replace_gaps, filter_by_n, date_accessions, lookup_date, separate_by_month, distance_matrix, set_to_uppercase, weight_by_cluster_size, date_pango_calls, create_rarity_lookup, assign_anachron, find_anachron_seqs, find_double_candidates

### FUNCTION(S) TO FILTER GENBANK METADATA TO A PARTICULAR GEOGRAPHY STRING ###
### ----------------------------------------------------------------------- ###

"""
Filter the Apache Arrow representation of the metadata, supplied by the argument
`input_table::String`, to only those rows that contain the string in argument
`geography::String`. Then, write the output to a new filtered arrow IPC file as
well as a text file listing one filter-passing accession per line.
"""
function filter_metadata_by_geo(input_table::String, geography::String)

    # Read in the metadata file as a dataframe
    metadata_df = DataFrame(Arrow.Table(input_table))

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert "Geographic Location" in names(metadata_df)
    @assert "Isolate Collection date" in names(metadata_df)
    @assert "Accession" in names(metadata_df)

    # filter metadata based on desired geography
    filtered = @pipe metadata_df |>
    filter(Symbol("Geographic Location") => value -> contains(string(value), geography),_)

    # Writing filtered metadata
    Arrow.write("filtered-to-geography.arrow", filtered)

    # separating out accessions
    accessions = filtered[!,"Accession"]

    # Writing accessions to a text file for use by seqtk subseq or seqkit grep
    writedlm("accessions.txt", accessions, "\n")
    
end

"""
Filter the Apache Arrow representation of the metadata, supplied by the argument
`input_table::String`, to only those rows that contain the string in argument
`geography::String` and contain a collection date that is after `min_date::String`. 
Then, write the output to a new filtered arrow IPC file as well as a text file 
listing one filter-passing accession per line.
"""
function filter_metadata_by_geo(input_table::String, geography::String, min_date::String)

    # Read in the metadata file as a dataframe
    metadata_df = DataFrame(Arrow.Table(input_table))

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert "Geographic Location" in names(metadata_df)
    @assert "Isolate Collection date" in names(metadata_df)
    @assert "Accession" in names(metadata_df)

    # filter metadata based on desired geography
    filtered = @pipe metadata_df |>
    filter(Symbol("Geographic Location") => value -> contains(string(value), geography),_)

    # filter metadata to desired minimum range
    @pipe filtered |>
    filter!(Symbol("Isolate Collection date") => date -> date > Dates.Date(min_date),_)

    # Writing filtered metadata
    Arrow.write("filtered-to-geography.arrow", filtered)

    # separating out accessions
    accessions = filtered[!,"Accession"]

    # Writing accessions to a text file for use by seqtk subseq or seqkit grep
    writedlm("accessions.txt", accessions, "\n")
    
end

"""
Filter the Apache Arrow representation of the metadata, supplied by the argument
`input_table::String`, to only those rows that contain the string in argument
`geography::String` and contain a collection date that is after `min_date::String`
and before `max_date::String`. Then, write the output to a new filtered arrow IPC
file as well as a text file listing one filter-passing accession per line.
"""
function filter_metadata_by_geo(input_table::String, geography::String, min_date::String, max_date::String)

    # Read in the metadata file as a dataframe
    metadata_df = DataFrame(Arrow.Table(input_table))

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert "Geographic Location" in names(metadata_df)
    @assert "Isolate Collection date" in names(metadata_df)
    @assert "Accession" in names(metadata_df)

    # filter metadata based on desired geography
    filtered = @pipe metadata_df |>
    filter(Symbol("Geographic Location") => value -> contains(string(value), geography),_)

    # filter metadata to desired date range
    @pipe filtered |>
    filter!(Symbol("Isolate Collection date") => date -> date > Dates.Date(min_date),_) |>
    filter!(Symbol("Isolate Collection date") => date -> date < Dates.Date(max_date),_)

    # Writing filtered metadata
    Arrow.write("filtered-to-geography.arrow", filtered)

    # separating out accessions
    accessions = filtered[!,"Accession"]

    # Writing accessions to a text file for use by seqtk subseq or seqkit grep
    writedlm("accessions.txt", accessions, "\n")
    
end

### ----------------------------------------------------------------------- ###



### FUNCTION(S) TO REPLACE "-" SYMBOLS WITH "N" CHARACTERS ###
### ------------------------------------------------------ ###

"""
Replace all "-" symbols in the nucleotide sequence of each input FASTA 
record (supplied by the argument `fasta_path::String`) with the masked base
character "N." Then, write the updated FASTA with the file name provided by
the argument `output_filename::String`.
"""
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

"""
Read the input FASTA from the file path provided by `input_fasta_path::String` and
remove any records with more than 10% of their bases masked with the "N" character.
This is potentially crucial quality control to ensure that consensus sequences 
processed downstream are not the result of un-rigorous bioinformatic processing.
The filtered FASTA is then written with the file name provided by the argument 
`output_filename::String`.
"""
function filter_by_n(input_fasta_path::String, max_ambiguity::Float64, output_filename::String)

    touch(output_filename)

    # create a file lock to prevent threads from corrupting the output file
    # u = ReentrantLock()

    # creating a loop that goes through sequence records and writes out
    # any sequences that have less than the minimum N count
    FastaWriter(output_filename , "a") do fa
        FastaReader(input_fasta_path) do fr
            # @sync for (name, seq) in fr
                # Threads.@spawn begin
            for (name, seq) in fr
                max_n_count = floor(length(seq) * max_ambiguity)
                n_count = count("N", convert(String, seq))
                if n_count < max_n_count
                            # Threads.lock(u) do
                    writeentry(fa, name, seq)
                            # end
                end
                # end
            end
        end
    end
end

### ----------------------------------------------------------------------------- ###



### FUNCTION(S) TO SEPARATE FASTA INTO ONE FASTA FOR EACH MONTH IN GENBANK METADATA ### 
### ------------------------------------------------------------------------------- ###

"""
Read an input metadata file in Apache Arrow format, supplied by the argument 
`input_path::String`, and create a more compact dictionary where accession strings
are the keys, and collection dates are the values, and return that dictionary for
usage downstream.
"""
function date_accessions(input_path::String)

    # read in metadata
    metadata_df = DataFrame(Arrow.Table(input_path))

    # create lookup of dates and accessions
    @assert "Accession" in names(metadata_df)
    @assert "Isolate Collection date" in names(metadata_df)
    accession_to_date = Dict(zip(metadata_df[!,"Accession"], metadata_df[!,"Isolate Collection date"]))

    return(accession_to_date)

end

"""
Return the date associated with an accession of interest provided by the argument
`record_name::String`. The lookup itself is provided with the argument `lookup::Dict`.
"""
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
"""
Read an input FASTA, supplied by the argument `input_fasta::String` and use the dictionary
supplied by the argument `accession_to_date::Dict` to separate each input FASTA record
into a new FASTA according to its collection month. Each FASTA record will be output into
a FASTA for the month it was collected in, such that records will be separated into one
FASTA per month present in the provided dataset.
"""
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
"""
The `vsearch --cluster_fast` algorithm produces some unusual patterns in its
outputs, including lowercase "n" symbols for masked bases. To avoid confusing
algorithms downstream, this function opens the input FASTA provided by the
argument `fasta_filename::String`, sets all its characters to uppercase, and
writes them to a temporary FASTA file with a name supplied by the argument
`temp_filename::String`.
"""
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

# determine how many sequences are in each cluster to compute weighted distances
"""
For the sequence identified by argument `seq_name::String` in the dataframe `dist_df::DataFrame`
and the `vsearch --cluster_fast` output table provided by argument `cluster_table::DataFrame`,
use the size of that sequence's cluster to "weight" its distances from other clusters. This 
weight, adjusted by `stringency::String`, will penalize large clusters and advantage small
clusters
"""
function weight_by_cluster_size(seq_name::String, stringency::String, dist_df::DataFrame, cluster_table::DataFrame)

    # filter down to centroid rows only
    centroids = filter(:1 => x -> x == "C", cluster_table)

    # make sure centroids are in the same order as the distance matrix
    centroids = centroids[indexin(names(dist_df), centroids.Column9),:]
    @assert centroids[:,9] == names(dist_df)

    # extract a vector of cluster sizes that are in the same order as the 
    # sequences in the distance matrix
    all_sizes = centroids[:,3]

    # filter down to hits only
    hits = filter(:1 => x -> x != "C", cluster_table)
    month_total = nrow(hits) == 0 ? 1 : nrow(hits)

    # find the number of sequences for the accession in question
    cluster_freq = filter(:9 => x -> x == seq_name, centroids)[:,3][1] / month_total

    # return the weight for this centroid's distance
    weights = stringency == "strict" ? ((all_sizes .* -log(cluster_freq)) ./ month_total) : ((all_sizes .* (1 - cluster_freq)) ./ month_total)
    @assert length(weights) == nrow(dist_df)
    return weights

end

# compute cluster-size-weighted weighted distance matrices

"""
Compute a pairwise Hamming distance matrix between all the sequences in the FASTA provided
by argument `temp_filename::String`. Use the information in argument `cluster_table::DataFrame`
and `stringency::String` to weight the distances by sequence ID, and use arguments `count::Int`,
`majority_centroid::String`, and `yearmonth::String` to structure and name output CSV files.
"""
function distance_matrix(temp_filename::String, cluster_table::DataFrame, yearmonth::String, stringency::String)

    # Collect both names and sequences
    seqs = [seq for (_, seq) in FastaReader(temp_filename)]
    seq_names = [name for (name, _) in FastaReader(temp_filename)]

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
        col_weights = weight_by_cluster_size(seq, stringency, dist_df, cluster_table)
        dist_df[!, seq] = dist_df[!, seq] .* col_weights
    end

    # Add a column for the sequence names
    dist_df[!, :Sequence_Name] = seq_names

    # Move the Sequence_Name column to the front
    select!(dist_df, :Sequence_Name, :)

    # Write the distance matrix to a CSV file
    CSV.write("$yearmonth-dist-matrix.csv", dist_df)

end

### ------------------------------------------------------------------------- ###



### FUNCTIONS TO IDENTIFY ANACHRONISTIC SEQUENCES WITH THE OUTBREAK.INFO API ###
### ------------------------------------------------------------------------ ###

"""
Take a the pangolin output report supplied by argument `report_path::String` and 
associated metadata supplied by the argument `metadata::DataFrame`, join the
metadata collection dates onto the pangolin report dataframe, and return it
for downstream usage.
"""
function date_pango_calls(report_path::String, metadata::DataFrame)

    # Read in the necessary pangolin report columns
    pango_report = CSV.read(report_path, DataFrame)

    @assert "taxon" in names(pango_report)
    @assert "Accession" in names(metadata)

    pango_with_dates = @pipe pango_report |>
        leftjoin(metadata, _, on = :taxon => :Accession) |>
        select(-:Accession) |>
        sort(:taxon)
    
    return pango_with_dates

end

"""
Use the [outbreak.info R package](https://outbreak-info.github.io/R-outbreak-info/) to access the outbreak.info
API, which we use to access prevalence estimates for the lineages provided by argument
`pango_report::DataFrame` in the United States. This is possible thanks to the Julia package
`RCall.jl`, which makes it easy to call R libraries and functions and convert their output
objects into Julia objects. Then, use the prevalence estimates to compute a "rarity date" for
each lineage. After that rarity date, any additional collections of that lineage may be considered
anachronistic, and therefore likely stem from prolonged infections or animal spillbacks.
"""
function create_rarity_lookup(pango_report::DataFrame)

    # create vector of unique lineages to iterate through
    lineages = unique(pango_report[:,"lineage"])

    # load R outbreakinfo package
    @rlibrary outbreakinfo

    # create empty rare dates vector
    rarity_lookup = Dict{String, Int}()

    # create rarity lookup vector
    for lineage in lineages

        # convert lineage into R object
        rlineage = robject(lineage)

        # call the outbreak.info API for prevalence estimates
        prevalence_table = rcopy(R"getPrevalence(pangolin_lineage = rlineage, location = 'United States')")

        @assert "proportion" in names(prevalence_table)
        @assert "date" in names(prevalence_table)

        # get rid of zeroes, which will skew our methods below. We also get
        # rid of 1.0's, which are likely erroneous.
        @pipe prevalence_table |>
            filter!(:proportion => prop -> prop > 0.0, _) |>
            filter!(:proportion => prop -> prop < 1.0, _)
        
        # deal with potential empty tables
        if nrow(prevalence_table) == 0
            rare_date = passmissing(Dates.Date(missing))
            my_dict[lineage] = rare_date
            break
        end

        # define a rarity level and then find the date when that level is reached
        rare_date = @pipe quantile(prevalence_table.proportion, 0.05) |>
            prevalence_table.proportion[argmin(abs.(prevalence_table.proportion .- _))] |>
            maximum(prevalence_table[prevalence_table.proportion .== _, :date])

        # make sure we are using a rarity date that is after the crest
        if rare_date <= crest_date
            @pipe prevalence_table |>
                filter!(:date => date -> date > crest_date, _)
            rare_date = @pipe post_crest.proportion[argmin(abs.(post_crest.proportion .- q))] |>
                maximum(post_crest[post_crest.proportion .== _, :date])
        end

        # make sure a lineage isn't too recent
        if (Dates.today() - crest_date) <= 30
            rare_date = passmissing(Dates.Date(missing))
            my_dict[lineage] = rare_date
            break
        end

        # if the lineage has not yet crested, by definition it cannot be
        # anachronistic
        if nrow(filter(:date => dates -> date > crest_date, prevalence_table)) == 0
            rare_date = passmissing(Dates.Date(missing))
            my_dict[lineage] = rare_date
            break
        end

        # push any rarity dates that passed the above conditions
        my_dict[lineage] = rare_date

    end

    return rarity_lookup

end

"""
Take the pangolin report, with dates, provided by argument `report_path::String`,
and use the lineages in that report and the collection dates from argument
`metadata_path::String` to assign "anachronicities" to each row. Here, we define
anachronicity as the number of days past the date at which a SARS-CoV-2 lineage
can be considered rare.
"""
function assign_anachron(metadata_path::String, report_path::String)

    # read in the metadata and sort by accession
    metadata = Arrow.Table(metadata_path) |> DataFrame
    @assert "Accession" in metadata.columns
    new_metadata = sort(metadata, :Accession)

    # date pango date_pango_calls
    pango_report = date_pango_calls(report_path, new_metadata)

    # make sure the accessions are the same and in the sample
    # order between the metadata and the pango report
    @assert new_metadata.Accession == pango_report.taxon
    @assert "Isolate Collection date" in names(pango_report)

    # create rare-date lookup with the outbreak.info API
    rarity_lookup = create_rarity_lookup(pango_report)

    # define anachronicities
    anachronicities = []
    for (lineage, date) in zip(pango_report.lineage, pango_report."Isolate Collection date")

        # compute anachronicity as an integer
        anachronicity = (date - get(rarity_lookup, lineage, Dates.today())).value

        # make sure anachronicity isn't missing
        if ismissing(anachronicity) 
            return 0
        end

        # append in place to the list anachronicities
        push!(anachronicities, anachronicity)
        
    end

    # make sure enough anachronicities were computed
    @assert length(anachronicities) == nrow(new_metadata)

    # add anachronicities as a column in the metadata
    new_metadata[:, "Anachronicity (days)"] = anachronicities

    # filter to anachronicities greater than 0
    filter!(Symbol("Anachronicity (days)") => days -> days > 0, new_metadata)

    return new_metadata

end

"""
Take the outputs of `assign_anachron()` and filter a FASTA file
to only those sequences with a positive "anachronicity", where 
anachronicity is defined as the number of days past the date at
which a SARS-CoV-2 lineage can be considered rare.
"""
function find_anachron_seqs(metadata::DataFrame, fasta_path::String)

    # separate out accessions into a list
    accessions = metadata.Accession

    # iterate through the input fasta and write anachronistics
    # to the output
    output_handle = "anachronistic_seq_candidates.fasta"
    touch(output_handle)
    FastaWriter(output_handle , "a") do fasta_append
        FastaReader(fasta_path) do fasta_read
            for (name, seq) in fasta_read
                if name in accessions
                    writeentry(fasta_append, name, seq)
                end
            end
        end
    end

end

### ------------------------------------------------------------------------ ###




### FUNCTION THAT FINDS THE OVERLAP BETWEEN HIGH-DISTANCE AND ANACHRONISTIC SEQUENCES ###
### --------------------------------------------------------------------------------- ###

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method compares between the output of a highly evolved detector and a single
anachronistic reporter, and then filters a FASTA to any identifiable double candidates.
"""
function find_double_candidates(metadata1::DataFrame, metadata2::DataFrame, seqs::String)

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadmetadata2[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata2_filtered = filter(row -> row[1] in common_accessions, metadata2)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata2_filtered, :Accession)
    
    if occursin("Anachronicity", names(metadata2_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata2_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim='\t')

    else

        @assert occursin("Anachronicity", names(metadata1_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata2_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata2_filtered, delim='\t')

    end

    # Open FASTA files and perform the filtering operation
    FastaReader(seqs) do fr
        FastaWriter("double_candidates.fasta", "w") do fw
            for (name, seq) in fr
                if name in common_accessions
                    writeentry(fw, name, seq)
                end
            end
        end
        
    end

end

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method looks for double candidates among a single highly evolved sequence
detector and two anachronistic sequence detectors, and then filters a FASTA to any
identifiable double candidates.
"""
function find_double_candidates(metadata1::DataFrame, metadata2::DataFrame, metadata3::DataFrame, seqs::String)

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert ncol(metadata3) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)
    @assert "Accession" in names(metadata3)

    # join the metadata based on anachronistic accessions
    if nrow(metadata2) >= nrow(metadata3)
        metadata_joined = leftjoin(metadata2, metadata3, on = :Accession)
    else
        metadata_joined = leftjoin(metadata3, metadata2, on = :Accession)
    end

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadmetadata_joined[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata_joined_filtered = filter(row -> row[1] in common_accessions, metadata_joined)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata_joined_filtered, :Accession)

    if occursin("Anachronicity", names(metadata_joined_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata_joined_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim='\t')

    else

        @assert occursin("Anachronicity", names(metadata1_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata_joined_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata_joined_filtered, delim='\t')

    end

    # Open FASTA files and perform the filtering operation
    FastaReader(seqs) do fr
        FastaWriter("double_candidates.fasta", "w") do fw
            for (name, seq) in fr
                if name in common_accessions
                    writeentry(fw, name, seq)
                end
            end
        end
        
    end

end

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method looks for double candidates from three anachronistic sequence detectors.
"""
function find_double_candidates(metadata1::DataFrame, metadata2::DataFrame, metadata3::DataFrame)

    # use assertions to catch runtime errors more informatively and make
    # assumptions about the data explicit
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert ncol(metadata3) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)
    @assert "Accession" in names(metadata3)

    # join the metadata based on anachronistic accessions
    if nrow(metadata2) >= nrow(metadata3)
        metadata_joined = leftjoin(metadata2, metadata3, on = :Accession)
    else
        metadata_joined = leftjoin(metadata3, metadata2, on = :Accession)
    end

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadmetadata_joined[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata_joined_filtered = filter(row -> row[1] in common_accessions, metadata_joined)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata_joined_filtered, :Accession)

    if "Anachronicity" in names(metadata_joined_filtered)

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata_joined_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim='\t')

    else

        @assert occursin("Anachronicity", names(metadata1_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata_joined_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata_joined_filtered, delim='\t')

    end

end

"""
Here, "double candidates" are defined as sequences that are both highly evolved,
with a large number of nucleotide substitutions relative to co-occurring sequences,
and are anachronistic, occurring long after a lineage is circulating and prevalent.
This method looks for double candidates from two anachronistic sequence detectors.
"""
function find_double_candidates(metadata1::DataFrame, metadata2::DataFrame)

    # use assertions to catch runtime errors more informatively
    @assert ncol(metadata1) > 0
    @assert ncol(metadata2) > 0
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)
    @assert "Accession" in names(metadata1)
    @assert "Accession" in names(metadata2)

    # Get the intersecting accessions
    common_accessions = intersect(metadata1[!, :Accession], metadmetadata2[!, :Accession])

    # Filter metadata to only common accessions
    metadata1_filtered = filter(row -> row[1] in common_accessions, metadata1)
    metadata2_filtered = filter(row -> row[1] in common_accessions, metadata2)

    # Sort both dataframes by Accession
    sort!(metadata1_filtered, :Accession)
    sort!(metadata2_filtered, :Accession)

    if "Anachronicity" in names(metadata2_filtered)

        # Add the last column from metadata2 to metadata1
        metadata1_filtered[!, :Anachronicity] = metadata2_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata1_filtered, delim='\t')

    else

        @assert occursin("Anachronicity", names(metadata1_filtered)[end])

        # Add the last column from metadata2 to metadata1
        metadata2_filtered[!, :Anachronicity] = metadata1_filtered[!, end]

        # Write the combined metadata to a TSV file
        CSV.write("double_candidate_metadata.tsv", metadata2_filtered, delim='\t')

    end

end

"""
Function `estimate_prevalence` uses the outputs from `seqkit stats` to
estimate the proportion of sequences flagged at the end of ALPINE in the
sequences used as the workflow's input.
"""
function estimate_prevalence(early_stats::DataFrame, late_stats::DataFrame)

    # use an assertion to make sure the necessary column names are present
    @assert "num_seqs" in names(early_stats)
    @assert "num_seqs" in names(late_stats)

    # make sure that candidates were found at all before pulling out
    # the final sample size
    candidate_count = nrow(late_stats) == 0 ? 0 : unique(late_stats.num_seqs)

    # pull out the input sample size
    sample_size = unique(early_stats.num_seqs)

    # compute the prevalence estimate
    prevalence = (candidate_count / sample_size) * 100

    return (prevalence, sample_size)

end

### ---------------------------------------------------------------------------------------- ###

end # LongInfectionFinder module