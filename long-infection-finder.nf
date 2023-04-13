#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// handling the case where no geography or date filters are provided
if( params.geography.isEmpty() || params.min_date.isEmpty() || params.max_date.isEmpty() ){
	params.ncbi_results = params.results + "/GenBank"
} else {
	params.ncbi_results = params.results + "/GenBank_" + params.geography + "_" + params.min_date + "_to_" + params.max_date
}

// creating results subfolders for the three orthogonal anachronistic
// sequence search methods
params.high_distance_candidates = params.ncbi_results + "/high_distance_cluster"
params.anachronistic_candidates = params.ncbi_results + "/anachronistic_candidates"
params.metadata_candidates = params.ncbi_results + "/metadata_candidates"

// --------------------------------------------------------------- //



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// Before anything else, make sure pangolin is up to date
	UPDATE_PANGO_CONTAINER ( )
	
	println "Pangolin updated to version:"
	UPDATE_PANGO_CONTAINER.out.cue.view()
	

	// Download the latest Pangolin designation dates from GitHub
	// GET_DESIGNATION_DATES ( )


	// NCBI/GENBANK BRANCH:
	/* Here we use three orthogonal methods for identifying prolonged
	infection candidates:
		1 - Creating a distance matrix of nucleotide differences between
		consensus sequences, run a clustering algorithm on those distances,
		and identify candidates in the highest-distance cluster.
		2 - Reclassifying sequences using with the most up-to-date version
		of Pangolin and comparing each isolate's collection date with 
		lineage prevalence estimates from outbreak.info
		3 - The fastest option: trusting pango lineages in the metadata 
		as being mostly up-to-date and comparing those lineage's collection 
		dates with outbreak.info prevalence estimates.
	We suggest users run all of these methods (see nextflow.config) and 
	cross reference the results from each.
	*/
	DOWNLOAD_SC2_PACKAGE ( )

	PREP_FOR_CLUSTERING (
		DOWNLOAD_SC2_PACKAGE.out.fasta
	)

	CLUSTER_BY_DISTANCE (
		PREP_FOR_CLUSTERING.out.fasta,
		DOWNLOAD_SC2_PACKAGE.out.metadata
	)

	COLLATE_CLUSTER_METADATA (
		CLUSTER_BY_DISTANCE.out,
		DOWNLOAD_SC2_PACKAGE.out.metadata
	)

	FILTER_NCBI_METADATA (
		DOWNLOAD_SC2_PACKAGE.out.metadata
	)

	// PULL_NCBI_SEQUENCES ( 
	// 	FILTER_NCBI_METADATA.out
	// 		.splitCsv ( header: true, sep: "\t" )
	// 		.map {row -> tuple(row.accession, row.date, row.location, row.pango)}
	// )

	// ALIGN_NCBI_SEQUENCES ()

	// SEQUENCE_DISTANCE_MATRIX ()
	
	HIGH_THROUGHPUT_PANGOLIN ( 
		UPDATE_PANGO_CONTAINER.out.cue,
		DOWNLOAD_SC2_PACKAGE.out.fasta
			.splitFasta( by: 5000, file: true )
	)

	CONCAT_PANGOLIN_REPORTS (
		HIGH_THROUGHPUT_PANGOLIN.out.collect()
	)

	FIND_CANDIDATE_LINEAGES_BY_DATE ( 
		GET_DESIGNATION_DATES.out,
		HIGH_THROUGHPUT_PANGOLIN.out
	)

	// FIND_CANDIDATE_LINEAGES_BY_DATE ( 
	// 	GET_DESIGNATION_DATES.out,
	// 	HIGH_THROUGHPUT_PANGOLIN.out
	// )

	// CONCAT_DATE_BASED_CANDIDATES (
	// 	FIND_CANDIDATE_LINEAGES_BY_DATE.out.collect()
	// )

	SEARCH_NCBI_METADATA ( 
		FILTER_NCBI_METADATA.out,
		GET_DESIGNATION_DATES.out
	)

	PULL_NCBI_METADATA_CANDIDATES (
		SEARCH_NCBI_METADATA.out
			.splitCsv( header: true, sep: "\t" )
			.map { row -> tuple( row.accession, row.infection_duration ) }
			.filter  { it[1].toInteger() >= params.duration_of_interest }
	)

	CONCAT_NCBI_METADATA_CANDIDATES (
		PULL_NCBI_METADATA_CANDIDATES.out.collect()
	)
	

}
// --------------------------------------------------------------- //



// PROCESS SPECIFICATIONS
// --------------------------------------------------------------- //

process UPDATE_PANGO_CONTAINER {
	
	// This process builds a new docker image with the latest available pangolin version
	
	output:
	env(version), emit: cue
	
	when:
	(workflow.profile == 'standard' || workflow.profile == 'docker' || workflow.profile == 'singularity') && (params.search_gisaid_seqs == true || params.search_genbank_seqs == true || !params.local_input_path.isEmpty() )
	
	script:
	"""
	pangolin --update --update-data
	version=`pangolin --version | sed 's/pangolin//g' | xargs`
	"""
}

process GET_DESIGNATION_DATES {
	
	// This process downloads a table of pangolin lineage designation dates
	// from Cornelius Roemer's GitHub. These dates represent when each lineage was
	// added to pangolin, after which point sequences could be classified as such 
	
	publishDir params.resources, mode: 'copy', overwrite: true
	
	output:
	path "*.csv"
	
	script:
	"""
	curl -fsSL https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv > lineage_designation_dates.csv
	"""
}


// NCBI/GENBANK PROCESSES:
process DOWNLOAD_SC2_PACKAGE {

	/*
	Here we download two, large files from NCBI: the FASTA of all 
	SARS-CoV-2 consensus sequences in GenBank, and a tab-delimited
	table of metadata for all those sequences. Depending on the 
	settings specified in nextflow.config, various processing
	will be performed on these files downstream.
	*/

	cpus 3

	output:
	path "*.fasta", emit: fasta
	path "*.tsv", emit: metadata

	script:
	"""

	datasets download virus genome taxon SARS-CoV-2 \
	--complete-only \
	--filename ${params.date}.zip && \
	unzip ${params.date}.zip

	mv ncbi_dataset/data/genomic.fna ./genbank_sc2_${params.date}.fasta

	mv ncbi_dataset/data/data_report.jsonl ./genbank_sc2_${params.date}.jsonl && \
	dataformat tsv virus-genome genbank_sc2_${params.date}.jsonl | \
	genbank_sc2_${params.date}.tsv

	rm -rf ${params.date}/
	"""

}

process PREP_FOR_CLUSTERING {

	/*
	In this process, we replace all dashes in the NCBI FASTA with
	"N's", effectively converting them into masked bases instead of
	gaps. We also sort the FASTA by sequence length in descending order,
	which will make the clustering algorithm run faster.
	*/

	input:
	path metadata

	output:

	when:
	params.make_distance_matrix == true

	script:
	"""
	"""
	
}

process CLUSTER_BY_DISTANCE {

	/*
	Here we run a version of Robert Edgar's UCLUST algorithm that 
	is part of the open-source VSEARCH package. This step may take 
	a large amount of time and RAM, and thus may be better suited 
	for a high-RAM cluster environment.
	*/

	cpus ${params.max_cpus}

	input:
	path fasta

	output:
	
	script:
	"""
	"""
	
}

process COLLATE_CLUSTER_METADATA {

	/*
	In this process, the clustering results are combined with metadata
	to provide as much information as possible with the candidate
	evolutionarily advanced sequences. We recommend users visually 
	inspect these candidates and compare them with candidates from
	the two additional methods used in this pipeline.
	*/

	publishDir params.high_distance_candidates, mode: 'copy'

	input:
	path clusters
	path metadata
}

process FILTER_NCBI_METADATA {

	/* 
	In parallel with the distance matrix method, this pipeline also pans
	the GenBank metadata for anachronistic sequences, which may have come
	from evolutionarily advanced virus lineages in prolonged infections.
	NCBI's pre-classified pango lineages, which come with the metadata,
	are what make this method possible. In this process, we tee off this
	method by filtering NCBI metadata for all rows with correctly formatted
	dates and pango lineages.
	*/

	input:
	path tsv

	output:
	path "*.tsv"

	script:
	"""
	filter_ncbi_metadata.R ${tsv} \
	${params.min_date} ${params.max_date} \
	${params.geography}
	"""
}

// process PULL_NCBI_SEQUENCES {
	
// 	tag "${accession}"
	
// 	cpus 1
	
// 	input:
// 	tuple val(accession), val(date), val(loc), val(pango)
	
// 	output:
// 	tuple val(accession), val(date), path("*.fasta")
	
// 	when:
// 	params.search_genbank_seqs == true

// 	script:
// 	"""
// 	datasets download virus genome accession "${accession}" \
// 	--exclude-cds --exclude-protein
// 	unzip ncbi_dataset.zip
// 	mv ncbi_dataset/data/genomic.fna ./"${accession}".fasta
// 	rm -rf ncbi_dataset/
// 	"""

// }

// process ALIGN_NCBI_SEQUENCES {

// 	cpus ${params.max_cpus}

// 	input:
// 	path fasta

// 	output:
// 	path "*.fasta"

// 	script:
// 	"""
// 	muscle -align seqs.fa -output aln.afa
// 	"""

// }

// process SEQUENCE_DISTANCE_MATRIX {}

process HIGH_THROUGHPUT_PANGOLIN {
	
	cpus 4
	errorStrategy 'retry'
	maxRetries 1
	
	input:
	each cue
	path fasta
	
	output:
	path "*.csv"

	script:
	"""
	pangolin \
	--skip-scorpio --skip-designation-cache \
	--threads ${task.cpus} \
	--outfile lineage_report.csv \
	"${fasta}"
	"""

}

process CONCAT_PANGOLIN_REPORTS {

	input:
	path lineage_report, stageAs: 'report??.csv'

	output:
	path csv

	script:
	"""
	"""
}

process FIND_CANDIDATE_LINEAGES_BY_DATE {

	publishDir params.ncbi_results, mode: 'copy'

	cpus ${params.max_cpus}

	input:
	path csv

	output:
	path "*putative_long_infections_ncbi*.csv"

	script:
	"""
	compare_lineage_prevalences.R ${csv} ${task.cpus}
	"""
}

// process FIND_CANDIDATE_LINEAGES_BY_DATE {

// 	publishDir params.ncbi_results, mode: 'copy'

// 	input:
// 	each path(lineage_dates)
// 	tuple path(lineage_csv), val(accession), val(date)

// 	output:
// 	path "*putative_long_infections_ncbi*.csv"

// 	script:
// 	"""
// 	ncbi_long_infection_finder.R ${lineage_dates} ${lineage_csv} ${date} ${params.days_of_infection}
// 	"""
// }

// process CONCAT_DATE_BASED_CANDIDATES {

// 	publishDir params.results, mode: 'copy'

// 	input:
// 	path file_list, stageAs: 'infections??.csv'

// 	output:
// 	path "*.csv"

// 	script:
// 	"""
// 	concat_long_infections.R ${params.days_of_infection}
// 	"""
// }

process SEARCH_NCBI_METADATA {

	publishDir params.ncbi_results, mode: 'copy'

	cpus 7

	input:
	path metadata
	path lineage_dates

	output:
	path "*.tsv"

	when:
	params.search_genbank_metadata == true
		
	script:
	"""
	search_ncbi_metadata.R ${metadata} ${lineage_dates} ${params.days_of_infection} ${task.cpus}
	"""

}

process PULL_NCBI_METADATA_CANDIDATES {

	tag "${accession}"

	cpus 1
	time { 1.minute * task.attempt }
	errorStrategy 'retry'
	maxRetries 4

	input:
	tuple val(accession), val(duration)

	output:
	path "*.fasta"

	script:
	"""
	datasets download virus genome accession "${accession}" && \
	unzip ncbi_dataset.zip && \
	mv ncbi_dataset/data/genomic.fna ./"${accession}".fasta && \
	rm -rf ncbi_dataset/
	"""

}

process CONCAT_NCBI_METADATA_CANDIDATES {

	publishDir params.ncbi_results, mode: 'copy'

	input:
	path fasta_list

	output:
	path "*.fasta.xz"

	script:
	"""
	find . -name "*.fasta" > fasta.list && \
	for i in `cat fasta.list`;
	do
		cat \$i >> ncbi_long_infection_candidates_${params.date}.fasta
	done && \
	xz -9 ncbi_long_infection_candidates_${params.date}.fasta
	"""
}


// --------------------------------------------------------------- //
