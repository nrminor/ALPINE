#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	UPDATE_PANGO_CONTAINER ( )
	
	if ( update_pango == true ){
		println "This workflow will use the following Pangolin version:"
		UPDATE_PANGO_CONTAINER.out.cue.view()
	}

	// Data setup steps
	DOWNLOAD_NCBI_PACKAGE ( )

	DOWNLOAD_REFSEQ ( )

	GET_DESIGNATION_DATES ( )

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

	// Distance matrix clustering steps
	REMOVE_FASTA_GAPS ( 
		DOWNLOAD_NCBI_PACKAGE.out.fasta
	)

	APPEND_DATES (
		DOWNLOAD_NCBI_PACKAGE.out.metadata,
		REMOVE_FASTA_GAPS.out
	)

	SEPARATE_BY_MONTH (
		APPEND_DATES.out
	)

	CLUSTER_BY_DISTANCE (
		PREP_FOR_CLUSTERING.out.flatten(),
		DOWNLOAD_NCBI_PACKAGE.out.metadata
	)

	PREP_CENTROID_FASTAS (
		CLUSTER_BY_DISTANCE.out.centroid_fasta
	)

	BUILD_CENTROID_TREE (
		PREP_CENTROID_FASTAS.out,
		DOWNLOAD_REFSEQ.out
	)

	GENERATE_CLUSTER_REPORT (
		CLUSTER_BY_DISTANCE.out.cluster_table.collect(),
		CLUSTER_BY_DISTANCE.out.cluster_fastas.collect(),
		BUILD_CENTROID_TREE.out.collect(),
		DOWNLOAD_NCBI_PACKAGE.out.metadata
	)
	
	// Steps for re-running pangolin and comparing dates
	HIGH_THROUGHPUT_PANGOLIN ( 
		UPDATE_PANGO_CONTAINER.out.cue,
		DOWNLOAD_NCBI_PACKAGE.out.fasta
			.splitFasta( by: 5000, file: true )
	)

	CONCAT_PANGOLIN_REPORTS (
		HIGH_THROUGHPUT_PANGOLIN.out.collect()
	)

	FIND_CANDIDATE_LINEAGES_BY_DATE (
		HIGH_THROUGHPUT_PANGOLIN.out
	)

	// Steps for inspecting NCBI metadata
	FILTER_NCBI_METADATA (
		DOWNLOAD_NCBI_PACKAGE.out.metadata
	)

	SEARCH_NCBI_METADATA ( 
		FILTER_NCBI_METADATA.out,
		GET_DESIGNATION_DATES.out
	)

	COLLATE_NCBI_METADATA_CANDIDATES (
		SEARCH_NCBI_METADATA.out
	)
	

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// Defining number of cpus to use base on execution environment
if ( workflow.profile == "chtc" ){
	params.max_cpus = executor.cpus
} else {
	params.max_cpus = params.max_local_cpus
}

// handling the case where no geography or date filters are provided
if( params.geography.isEmpty() || params.min_date.isEmpty() || params.max_date.isEmpty() ){
	params.ncbi_results = params.results + "/GenBank"
} else {
	params.ncbi_results = params.results + "/GenBank_" + params.geography + "_" + params.min_date + "_to_" + params.max_date
}

// creating results subfolders for the three orthogonal anachronistic
// sequence search methods
params.clustering_results = params.ncbi_results + "/all_clustering_results"
params.high_distance_candidates = params.ncbi_results + "/high_distance_cluster"
params.anachronistic_candidates = params.ncbi_results + "/anachronistic_candidates"
params.metadata_candidates = params.ncbi_results + "/metadata_candidates"

// --------------------------------------------------------------- //



// PROCESS SPECIFICATIONS
// --------------------------------------------------------------- //

process UPDATE_PANGO_CONTAINER {
	
	// This process builds a new docker image with the latest available pangolin version
	
	output:
	env(version), emit: cue
	
	when:
	workflow.profile == 'standard' || workflow.profile == 'docker' || workflow.profile == 'singularity'
	
	script:
	"""
	pangolin --update --update-data
	version=`pangolin --version | sed 's/pangolin//g' | xargs`
	"""
}

process DOWNLOAD_NCBI_PACKAGE {

	/*
	Here we download two, large files from NCBI: the FASTA of all 
	pathogen consensus sequences in GenBank, and a tab-delimited
	table of metadata for all those sequences. Depending on the 
	settings specified in nextflow.config, various processing
	will be performed on these files downstream.
	*/

	cpus 3

	output:
	path "*.fasta", emit: fasta
	path "*.tsv", emit: metadata

	script:
	if ( params.geograpy != "" )
		"""

		datasets download virus genome taxon ${params.pathogen} \
		--complete-only \
		--filename ${params.date}.zip \
		--geo-location ${params.geography} && \
		unzip ${params.date}.zip

		mv ncbi_dataset/data/genomic.fna ./genbank_${params.date}.fasta

		mv ncbi_dataset/data/data_report.jsonl ./genbank_${params.date}.jsonl && \
		dataformat tsv virus-genome \
		--force \
		--inputfile genbank_${params.date}.jsonl \
		| genbank_${params.date}.tsv

		rm -rf ${params.date}/
		"""
	else
		"""

		datasets download virus genome taxon ${params.pathogen} \
		--complete-only \
		--filename ${params.date}.zip && \
		unzip ${params.date}.zip

		mv ncbi_dataset/data/genomic.fna ./genbank_${params.date}.fasta

		mv ncbi_dataset/data/data_report.jsonl ./genbank_${params.date}.jsonl && \
		dataformat tsv virus-genome \
		--force \
		--inputfile genbank_${params.date}.jsonl \
		| genbank_${params.date}.tsv

		rm -rf ${params.date}/
		"""	

}

process DOWNLOAD_REFSEQ {

	publishDir params.resources, mode: 'copy'

	output:
	path "refseq.fasta"

	script:
	"""
	datasets download virus genome taxon ${params.pathogen} \
	--refseq && \
	unzip ncbi_dataset.zip
	mv ncbi_dataset/data/genomic.fna ./refseq.fasta
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
	-fsSL https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv > lineage_designation_dates.csv
	"""
}

process REMOVE_FASTA_GAPS {

	/*
	In this process, we replace all dashes in the NCBI FASTA with
	"N's", effectively converting them into masked bases instead of
	gaps.
	*/

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.fasta"

	when:
	params.make_distance_matrix == true

	script:
	"""
	remove_fasta_gaps.py ${fasta} no_gaps.fasta ${params.max_cpus}
	"""
	
}

process APPEND_DATES {

	/*
	This process appends collection dates onto the ends of each
	FASTA defline by cross referencing the GenBank accessions in 
	the FASTA with GenBank accessions in the NCBI metadata.
	*/

	input:
	path fasta
	path metadata

	output:
	path "*.fasta"

	script:
	"""
	append_dates.py ${metadata} ${fasta} dated_seqs.fasta
	"""

}

process SEPARATE_BY_MONTH {

	/*
	In this process, we get ready to parallelize by splitting the big
	NCBI FASTA into one for each year-month combination. This will 
	allow us to identify particularly advanced viruses in each 
	month instead of simply flagging the newest variants as the 
	most evolved (which, almost by definition, they are!)
	*/

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.fasta"

	script:
	"""
	separate_by_year-month.py ${fasta} ${task.cpus}
	"""

}

process CLUSTER_BY_DISTANCE {

	/*
	Here we run a version of Robert Edgar's UCLUST algorithm that 
	is part of the open-source VSEARCH package. This step may take 
	a large amount of time and RAM, and thus may be better suited 
	for a high-RAM cluster environment.
	*/

	tag "${yearmonth}"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	cpus ${params.max_cpus}

	input:
	path fasta
	each path(metadata)

	output:
	path "*.uc", emit: cluster_table
	tuple path("*centroids.fasta"), val(yearmonth), emit: centroid_fasta
	path "*-cluster", emit: cluster_fastas
	
	script:
	yearmonth = fasta.getSimpleName()
	"""
	vsearch --cluster_fast ${fasta} \
	--id ${params.id_threshold} \
	--centroids ${yearmonth}-centroids.fasta \
	--uc ${yearmonth}-clusters.uc \
	--clusters ${yearmonth}-cluster \
	--threads ${task.cpus}
	"""
	
}

process PREP_CENTROID_FASTAS {

	/*
	In this process, we align centroid sequences to Wuhan-1, which
	will be used as an outgroup, and replace "/" symbols with underscores
	to keep iqTree happy.
	*/

	input:
	tuple path(fasta), val(yearmonth)
	path refseq

	output:
	tuple path("*.fasta"), val(yearmonth)

	script:
	"""
	prep_tree_fasta.py ${fasta} ${refseq} ${yearmonth}
	"""

}

process BUILD_CENTROID_TREE {

	/*
	Next, we take the centroid sequences aligned with Wuhan-1
	and build a tree with Wuhan-1 as the outgroup. The goal is
	to identify which centroid has the longest branch length 
	for each year-month. The sequences that cluster with that
	long-branch centroid will be classified as a evolutionarily
	advanced in the final report.
	*/

	tag "${yearmonth}"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	input:
	tuple path(fasta), val(yearmonth)
	path refseq

	output:
	path "*.treefile"

	script:
	"""
	ref_id=$(grep "^>" ${refseq} | head -n 1 | cut -d' ' -f1 | sed 's/^>//')
	iqtree -s ${fasta} -o \${ref_id} -pre ${yearmonth} -m MFP -bb 1000 -nt ${task.cpus}
	"""
}

process GENERATE_CLUSTER_REPORT {

	/*
	In this process, the clustering results are combined with metadata
	to provide as much information as possible with the candidate
	evolutionarily advanced sequences. We recommend users visually 
	inspect these candidates and compare them with candidates from
	the two additional methods used in this pipeline.
	*/

	publishDir params.high_distance_candidates, mode: 'copy'

	input:
	path cluster_tables
	path cluster_fastas
	path cluster_trees
	path metadata

	output:

	script:
	"""
	generate_cluster_report.R ${metadata}
	"""

}

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

	publishDir params.anachronistic_candidates, mode: 'copy'

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

	when:
	params.inspect_ncbi_metadata == true && update_pango == false

	script:
	"""
	filter_ncbi_metadata.R ${tsv} \
	${params.min_date} ${params.max_date} \
	${params.geography}
	"""
}

process SEARCH_NCBI_METADATA {

	publishDir params.ncbi_results, mode: 'copy'

	cpus 7

	input:
	path metadata
	path lineage_dates

	output:
	path "*.tsv"

	when:
	params.inspect_ncbi_metadata == true
		
	script:
	"""
	search_ncbi_metadata.R ${metadata} ${lineage_dates} ${params.days_of_infection} ${task.cpus}
	"""

}

process COLLATE_NCBI_METADATA_CANDIDATES {

	input:
	tuple val(accession), val(duration)

	output:
	path "*.fasta"

	script:
	"""
	"""

}


// --------------------------------------------------------------- //
