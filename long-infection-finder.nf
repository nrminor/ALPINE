#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	if ( params.update_pango == true ){

		UPDATE_PANGO_CONTAINER ( )

		println "This workflow will use the following Pangolin version:"
		UPDATE_PANGO_CONTAINER.out.cue.view()

	}

	// Data setup steps
	DOWNLOAD_REFSEQ ( )
	
	if ( params.compare_lineage_dates == true ){

		GET_DESIGNATION_DATES ( )

	}

	if ( params.fasta_path == "" || params.metadata_path == "" ) {

		DOWNLOAD_NCBI_PACKAGE ( )

		UNZIP_NCBI_METADATA (
			DOWNLOAD_NCBI_PACKAGE.out.zip_archive
		)

		UNZIP_NCBI_FASTA (
			DOWNLOAD_NCBI_PACKAGE.out.zip_archive
		)

		FILTER_TSV_TO_GEOGRAPHY (
			UNZIP_NCBI_METADATA.out
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			UNZIP_NCBI_FASTA.out
				.splitFasta( by: 5000, file: "genbank-${params.pathogen}.fasta" ),
			FILTER_TSV_TO_GEOGRAPHY.out.accessions
		)

	} else {

		ch_local_fasta = Channel
			.fromPath( params.fasta_path )

		ch_local_metadata = Channel
			.fromPath( params.metadata_path )

		FILTER_TSV_TO_GEOGRAPHY (
			ch_local_metadata
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			ch_local_fasta
				.splitFasta( by: 5000, file: "genbank-${params.pathogen}.fasta" ),
			FILTER_TSV_TO_GEOGRAPHY.out.accessions
		)

	}

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
		FILTER_SEQS_TO_GEOGRAPHY.out
	)

	FILTER_BY_MASKED_BASES (
		REMOVE_FASTA_GAPS.out
	)

	APPEND_DATES (
		FILTER_TSV_TO_GEOGRAPHY.out.metadata,
		FILTER_BY_MASKED_BASES.out
	)

	SEPARATE_BY_MONTH (
		APPEND_DATES.out
			.collectFile( name: "${params.pathogen}_prepped.fasta", newLine: true )
	)

	CLUSTER_BY_DISTANCE (
		SEPARATE_BY_MONTH.out.flatten()
	)

	PREP_CENTROID_FASTAS (
		CLUSTER_BY_DISTANCE.out.centroid_fasta,
		DOWNLOAD_REFSEQ.out
	)

	BUILD_CENTROID_TREE (
		PREP_CENTROID_FASTAS.out,
		PREP_CENTROID_FASTAS.out
			.map { fasta, yearmonth -> file(fasta).countFasta() },
		DOWNLOAD_REFSEQ.out
	)

	// MDS_PLOT (
	// 	CLUSTER_BY_DISTANCE.out.cluster_fastas
	// )

	// PLOT_TREE (
	// 	BUILD_CENTROID_TREE.out
	// )

	GENERATE_CLUSTER_REPORT (
		CLUSTER_BY_DISTANCE.out.cluster_table.collect(),
		CLUSTER_BY_DISTANCE.out.cluster_fastas.collect(),
		BUILD_CENTROID_TREE.out.collect(),
		FILTER_TSV_TO_GEOGRAPHY.out.metadata
	)

	RUN_META_CLUSTER (
		GENERATE_CLUSTER_REPORT.out.high_dist_seqs
	)

	// META_CLUSTER_REPORT (
	// 	RUN_META_CLUSTER.out.cluster_fastas
	// 		.flatten(),
	// 	RUN_META_CLUSTER.out.cluster_fastas
	// 		.flatten()
	// 		.map { fasta, yearmonth -> file(fasta).countFasta() }
	// )
	
	// Steps for re-running pangolin and comparing dates
	// HIGH_THROUGHPUT_PANGOLIN ( 
	// 	UPDATE_PANGO_CONTAINER.out.cue,
	// 	FILTER_SEQS_TO_GEOGRAPHY.out
	// 		.splitFasta( by: 5000, file: true )
	// )

	// CONCAT_PANGOLIN_REPORTS (
	// 	HIGH_THROUGHPUT_PANGOLIN.out.collect()
	// )

	// FIND_CANDIDATE_LINEAGES_BY_DATE (
	// 	HIGH_THROUGHPUT_PANGOLIN.out
	// )

	// Steps for inspecting NCBI metadata
	// FILTER_NCBI_METADATA (
	// 	FILTER_TSV_TO_GEOGRAPHY.out.metadata
	// )

	// SEARCH_NCBI_METADATA ( 
	// 	FILTER_NCBI_METADATA.out,
	// 	GET_DESIGNATION_DATES.out
	// )

	// COLLATE_NCBI_METADATA_CANDIDATES (
	// 	SEARCH_NCBI_METADATA.out
	// )
	

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// Defining number of cpus to use base on execution environment
if ( workflow.profile == "chtc" ){
	params.max_cpus = executor.cpus
} else {
	params.max_cpus = params.max_shared_cpus
}

// specifying whether to run in low disk mode
if( params.low_disk_mode == true ) {
    params.publishMode = 'symlink'
}
else {
    params.publishMode = 'copy'
}

// Making a date-stamped results folder
params.dated_results = params.results + "/" + params.date

// handling the case where no geography or date filters are provided
if( params.geography == "" || params.min_date == "" || params.max_date == "" ){
	params.ncbi_results = params.dated_results + "/GenBank"
} else {
	params.ncbi_results = params.dated_results + "/GenBank_" + params.geography + "_" + params.min_date + "_to_" + params.max_date
}

// creating results subfolders for the three orthogonal anachronistic
// sequence search methods
params.clustering_results = params.ncbi_results + "/all_clustering_results"
params.repeat_lineages = params.clustering_results + "/repeat_lineages"
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

process DOWNLOAD_REFSEQ {

	/*
	Here the NCBI RefSeq for the selected pathogen, if 
	available, is downloaded for downstream usage.
	*/

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
	curl -fsSL https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv > lineage_designation_dates.csv
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

	tag "${params.pathogen}"
	publishDir params.dated_results, mode: params.publishMode

	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 5

	output:
	path "*.zip", emit: zip_archive

	script:
	"""
	datasets download virus genome taxon ${params.pathogen} --complete-only
	"""	

}

process UNZIP_NCBI_METADATA {

	/*
	Here the pathogen metadata in TSV format is extracted from the 
	lightweight NCBI zip archive.
	*/

	publishDir params.dated_results, mode: params.publishMode

	cpus 3

	input:
	path zip

	output:
	path "*.tsv"

	script:
	"""
	unzip -p ${zip} ncbi_dataset/data/data_report.jsonl \
	| dataformat tsv virus-genome --force > genbank_metadata.tsv
	"""

}

process UNZIP_NCBI_FASTA {

	/*
	Here the pathogen sequences in FASTA format are extracted 
	from the lightweight NCBI zip archive.
	*/

	publishDir params.dated_results, mode: params.publishMode

	cpus 3

	input:
	path zip

	output:
	path "*.fasta"

	script:
	"""
	unzip -p ${zip} ncbi_dataset/data/genomic.fna | cat > genbank_sequences.fasta
	"""

}

process FILTER_TSV_TO_GEOGRAPHY {

	/*
	The geography filter in the NCBI Datasets command line
	interface appears to be broken. Until this is fixed,
	this process uses a simple julia script to filter
	the metadata down to a geography of interest.
	*/

	publishDir params.dated_results, mode: params.publishMode

	cpus 8

	input:
	path metadata

	output:
	path "*.tsv", emit: metadata
	path "*.txt", emit: accessions

	script:
	"""
	filter-to-geography.jl ${metadata} ${params.geography}
	"""

}

process FILTER_SEQS_TO_GEOGRAPHY {

	/*
	This process takes the aceessions list from 
	FILTER_TSV_TO_GEOGRAPHY and filters down the full FASTA
	to those accessions, ensuring that both the metadata and
	the sequences reflect the same geography filtering.

	In the future, this process will be parallelized to split
	up the large genbank FASTA into many more computationally
	bite-sized pieces.
	*/

	publishDir params.dated_results, mode: params.publishMode

	cpus 8

	input:
	each path(fasta)
	path accessions

	output:
	path "*.fasta"

	script:
	"""
	seqtk subseq ${fasta} ${accessions} > filtered_to_geography.fasta 
	"""

}

process REMOVE_FASTA_GAPS {

	/*
	In this process, we replace all dashes in the NCBI FASTA with
	"N's", effectively converting them into masked bases instead of
	gaps.
	*/

	cpus 8

	input:
	path fasta

	output:
	path "*.fasta"

	when:
	params.make_distance_matrix == true

	script:
	"""
	remove_fasta_gaps.py ${fasta} no_gaps.fasta ${task.cpus}
	"""
	
}

process FILTER_BY_MASKED_BASES {

	/*
	This process makes sure that none of the sequences in the
	Genbank FASTA have fewer than 10% of their bases masked,
	i.e., are "N" instead of a defined base.
	*/

	input:
	path fasta

	output:
	path "*.fasta"

	script:
	"""
	filter-by-n-count.jl `realpath ${fasta}`
	"""

}

process APPEND_DATES {

	/*
	This process appends collection dates onto the ends of each
	FASTA defline by cross referencing the GenBank accessions in 
	the FASTA with GenBank accessions in the NCBI metadata.
	*/

	input:
	path metadata
	path fasta

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

	cpus 8

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

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.uc", emit: cluster_table
	path "*centroids.fasta", emit: centroid_fasta
	path "*-cluster*", emit: cluster_fastas
	
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

	tag "${yearmonth}"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	input:
	each path(fasta)
	path refseq

	output:
	path "${yearmonth}-centroids-with-ref.fasta"

	script:
	yearmonth = fasta.getSimpleName.replace("-centroids", "")
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
	each path(fasta)
	val seq_count
	path refseq

	output:
	path "*.treefile"

	when:
	seq_count.toInteger() > 3

	script:
	yearmonth = fasta.getSimpleName.replace("-centroids-with-ref", "")
	"""
	ref_id=\$(grep "^>" ${refseq} | head -n 1 | cut -d' ' -f1 | sed 's/^>//')
	iqtree -s ${fasta} -o \${ref_id} -pre ${yearmonth} -m MFP -bb 1000 -nt ${task.cpus}
	"""

}

process MDS_PLOT {

	/*
	*/

	tag "${yearmonth}"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	input:
	path cluster_seqs

	output:
	path "*.pdf"

	script:
	"""
	"""

}

process PLOT_TREE {

	/*
	*/

	tag "${yearmonth}"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	input:
	path treefile

	output:
	path "*.pdf"

	script:
	"""
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
	path "*.tsv", emit: metadata
	path "*.fasta", emit: high_dist_seqs

	script:
	"""
	generate_cluster_report.R ${metadata}
	"""

}

process RUN_META_CLUSTER {

	/*
	*/

	publishDir params.repeat_lineages, mode: 'copy'

	cpus params.max_cpus

	input:
	path fastas

	output:
	path "*.uc", emit: cluster_table
	tuple path("centroids.fasta"), val(yearmonth), emit: centroid_fasta
	path "cluster*", emit: cluster_fastas

	script:
	"""
	vsearch --cluster_fast ${fasta} \
	--id 0.9999 \
	--centroids centroids.fasta \
	--uc clusters.uc \
	--clusters cluster \
	--threads ${task.cpus}
	"""
}

process META_CLUSTER_REPORT {

	/*
	*/

	publishDir params.repeat_lineages, mode: 'copy'

	cpus params.max_cpus

	input:
	path fasta
	path seq_count

	output:

	when:
	seq_count.toInteger() > 1

	script:
	"""
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

	cpus params.max_cpus

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
	params.inspect_ncbi_metadata == true && params.update_pango == false

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
