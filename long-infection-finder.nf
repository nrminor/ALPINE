#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	ch_gisaid_token = Channel
		.fromPath ( params.gisaid_token )

	GET_DESIGNATION_DATES ( )

	if ( params.compare_lineage_dates == true ){

		UPDATE_PANGO_CONTAINER ( )

		println "This workflow will use the following Pangolin version:"
		UPDATE_PANGO_CONTAINER.out.cue.view()

	}

	// Data setup steps
	DOWNLOAD_REFSEQ ( )

	if ( params.fasta_path == "" || params.metadata_path == "" ) {

		DOWNLOAD_NCBI_PACKAGE ( )

		UNZIP_NCBI_METADATA (
			DOWNLOAD_NCBI_PACKAGE.out.zip_archive
		)

		EXTRACT_NCBI_FASTA (
			DOWNLOAD_NCBI_PACKAGE.out.zip_archive
		)

		VALIDATE_METADATA (
			UNZIP_NCBI_METADATA.out
		)

		// STORE_METADATA_WITH_ARROW (
		// 	VALIDATE_METADATA.out
		// )

		FILTER_META_TO_GEOGRAPHY (
			VALIDATE_METADATA.out
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			EXTRACT_NCBI_FASTA.out,
			FILTER_META_TO_GEOGRAPHY.out.accessions
		)

		REMOVE_FASTA_GAPS ( 
			FILTER_SEQS_TO_GEOGRAPHY.out.fasta
		)

	} else {

		if ( workflow.profile == 'standard' || workflow.profile == 'docker' ){

			ch_local_fasta = Channel
				.fromPath( params.fasta_path )

		} else {

			ch_local_fasta = Channel
				.fromPath( params.fasta_path )
				.splitFasta( by: 5000, file: "genbank-${params.pathogen}.fasta" )

		}

		ch_local_metadata = Channel
			.fromPath( params.metadata_path )

		VALIDATE_METADATA (
			ch_local_metadata
		)

		// STORE_METADATA_WITH_ARROW (
		// 	VALIDATE_METADATA.out
		// )

		FILTER_META_TO_GEOGRAPHY (
			VALIDATE_METADATA.out
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			ch_local_fasta,
			FILTER_META_TO_GEOGRAPHY.out.accessions
		)

		REMOVE_FASTA_GAPS ( 
			FILTER_SEQS_TO_GEOGRAPHY.out.fasta
		)

	}

	// ADVANCED VIRUS &/OR LONG INFECTION FINDER:
	/* This workflow uses three orthogonal methods for identifying prolonged
	infection candidates:
		1 - Running a clustering algorithm on the sequences from each month to reduce
		the dimensionality of the data, creating a distance matrix of nucleotide 
		differences between clusters, and then flagging the sequences from the 
		highest-distance clusters.
		2 - Reclassifying sequences using with the most up-to-date version
		of Pangolin and comparing each isolate's collection date with 
		lineage prevalence estimates from outbreak.info
		3 - The fastest option: trusting pango lineages in the metadata 
		as being mostly up-to-date and comparing those lineages' collection 
		dates with outbreak.info prevalence estimates.
	We suggest users run all of these methods (see nextflow.config) and 
	cross reference the results from each.
	*/

	// Distance matrix clustering steps
	FILTER_BY_MASKED_BASES (
		REMOVE_FASTA_GAPS.out
	)

	SEPARATE_BY_MONTH (
		FILTER_BY_MASKED_BASES.out
			.collectFile( name: "${params.pathogen}_prepped.fasta", newLine: true ),
		FILTER_META_TO_GEOGRAPHY.out.metadata
	)

	CLUSTER_BY_IDENTITY (
		SEPARATE_BY_MONTH.out.flatten()
	)

	PREP_CENTROID_FASTAS (
		CLUSTER_BY_IDENTITY.out.centroid_fasta
			.filter { it[1].toInteger() > 1 }
			.map { fasta, count -> fasta }
	)

	COMPUTE_DISTANCE_MATRIX (
		CLUSTER_BY_IDENTITY.out.cluster_table,
		PREP_CENTROID_FASTAS.out
			.filter { it[1].toInteger() > 1 }
			.map { fasta, count -> fasta }
	)

	MULTIDIMENSIONAL_SCALING (
		CLUSTER_BY_IDENTITY.out.cluster_table,
		PREP_CENTROID_FASTAS.out
			.filter { it[1].toInteger() > 2 }
			.map { fasta, count -> fasta }
	)

	SUMMARIZE_CANDIDATES (
		CLUSTER_BY_IDENTITY.out.cluster_table.collect(),
		CLUSTER_BY_IDENTITY.out.cluster_fastas.collect(),
		COMPUTE_DISTANCE_MATRIX.out.collect(),
		FILTER_META_TO_GEOGRAPHY.out.metadata
	)

	RUN_META_CLUSTER (
		SUMMARIZE_CANDIDATES.out.high_dist_seqs
			.filter { it[1].toInteger() > 2 }
			.map { fasta, count -> fasta }
	)

	META_CLUSTER_REPORT (
		RUN_META_CLUSTER.out.cluster_table,
		RUN_META_CLUSTER.out.cluster_fastas,
		SUMMARIZE_CANDIDATES.out.metadata,
		RUN_META_CLUSTER.out.whether_repeats
	)
	
	// Steps for re-running pangolin and comparing dates
	CLASSIFY_SC2_WITH_PANGOLIN (
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta
	)

	FIND_CANDIDATE_LINEAGES_BY_DATE (
		CLASSIFY_SC2_WITH_PANGOLIN.out
			.collectFile(name: 'new_pango_calls.csv', newLine: true),
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta,
		FILTER_META_TO_GEOGRAPHY.out.metadata,
		ch_gisaid_token
	)

	// Steps for inspecting NCBI metadata
	SEARCH_NCBI_METADATA ( 
		FILTER_META_TO_GEOGRAPHY.out.metadata,
		GET_DESIGNATION_DATES.out
	)

	// Steps for analyzing and visualizing the results from the 
	// approaches above
	FIND_DOUBLE_CANDIDATES (
		SUMMARIZE_CANDIDATES.out.metadata,
		FIND_CANDIDATE_LINEAGES_BY_DATE.out.metadata,
		FIND_CANDIDATE_LINEAGES_BY_DATE.out.sequences
			.filter { it[1].toInteger() > 2 }
			.map { fasta, count -> fasta }
	)

	// PREP_FOR_ESCAPE_CALC ()

	// RUN_ESCAPE_CALC ()

	// BUILD_NEXTSTRAIN_TREE ()

	// PREP_FOR_USHER ()

	// PLACE_WITH_USHER ()
	

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

if ( params.pathogen == "SARS-CoV-2" || params.pathogen == "sars-cov-2" ){
	pathogen = "2697049"
}

// Making a date-stamped results folder
params.dated_results = params.results + "/" + params.date

// handling the case where no geography or date filters are provided
if( params.geography == "" ){
	params.ncbi_results = params.dated_results + "/GenBank"
} else {
	params.ncbi_results = params.dated_results + "/GenBank_" + params.geography
}

// creating results subfolders for the three orthogonal anachronistic
// sequence search methods
params.clustering_results = params.ncbi_results + "/all_clustering_results"
params.high_distance_candidates = params.ncbi_results + "/high_distance_clusters"
params.repeat_lineages = params.high_distance_candidates + "/repeat_lineages"
params.anachronistic_candidates = params.ncbi_results + "/anachronistic_candidates"
params.metadata_candidates = params.ncbi_results + "/metadata_candidates"
params.double_candidates = params.ncbi_results + "/double_candidates"

// --------------------------------------------------------------- //



// PROCESS SPECIFICATIONS
// --------------------------------------------------------------- //

process UPDATE_PANGO_CONTAINER {
	
	// This process builds a new docker image with the latest available pangolin version
	
	output:
	env(version), emit: cue
	
	when:
	(params.pathogen == "SARS-CoV-2" || params.pathogen == "sars-cov-2") && params.compare_lineage_dates == true
	
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

	label "lif_container"
	tag "${params.pathogen}"
	publishDir params.resources, mode: 'copy'

	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 5

	output:
	path "*.fasta", emit: ref_fasta
	env ref_id, emit: ref_id

	script:
	"""
	datasets download virus genome taxon ${pathogen} \
	--refseq && \
	unzip ncbi_dataset.zip
	mv ncbi_dataset/data/genomic.fna ./${params.pathogen}_refseq.fasta
	ref_id=\$(grep "^>" ${params.pathogen}_refseq.fasta | head -n 1 | cut -d' ' -f1 | sed 's/^>//')
	"""
}

process GET_DESIGNATION_DATES {
	
	// This process downloads a table of pangolin lineage designation dates
	// from Cornelius Roemer's GitHub. These dates represent when each lineage was
	// added to pangolin, after which point sequences could be classified as such 

	label "lif_container"
	publishDir params.resources, mode: 'copy', overwrite: true
	
	output:
	path "*.csv"

	when:
	params.pathogen == "SARS-CoV-2"
	
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
	publishDir params.ncbi_results, mode: params.publishMode

	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 5

	output:
	path "*.zip", emit: zip_archive

	script:
	"""
	datasets download virus genome taxon ${pathogen} --complete-only
	"""	

}

process UNZIP_NCBI_METADATA {

	/*
	Here the pathogen metadata in TSV format is extracted from the 
	lightweight NCBI zip archive.
	*/

	tag "${params.pathogen}"

	publishDir params.ncbi_results, pattern: "*.tsv", mode: params.publishMode

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

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

process EXTRACT_NCBI_FASTA {

	/*
	Here the pathogen sequences in FASTA format are extracted 
	from the lightweight NCBI zip archive.
	*/

	tag "${params.pathogen}"

	label "lif_container"
	publishDir params.ncbi_results, mode: params.publishMode

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 3

	input:
	path zip

	output:
	path "*.fasta.zst"

	script:
	"""
	unzip -p ${zip} ncbi_dataset/data/genomic.fna | zstd -10 -o genbank_sequences.fasta.zst
	"""

}

process VALIDATE_METADATA {

	/*
	This step checks the typing and column header names 
	for the input metadata, ensuring in particular that 
	GISAID metadata are compatible with the scripts used
	in this workflow.
	*/

	tag "${params.pathogen}"

	label "lif_container"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 8

	input:
	path metadata

	output:
	path "*.arrow"

	script:
	"""
	JULIA_NUM_THREADS=8 && \
	validate-metadata.jl ${metadata}
	"""
}

// process STORE_METADATA_WITH_ARROW {

// 	/*
// 	This workflow uses the Apache Arrow in-memory representation
// 	of metadata. This both improves the speed of the computations
// 	it runs on large metadata and also the speed read/write, 
// 	which is the predominant bottleneck for most processes 
// 	in this workflow.

// 	To convert the metadata to an arrow representation, we
// 	use a suite of tools written in Rust by Dominik Moritz 
// 	called arrow-tools, which is available at:
// 	https://github.com/domoritz/arrow-tools
// 	*/

// 	tag "${params.pathogen}"

// 	label "lif_container"

// 	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
// 	maxRetries 2

// 	input:
// 	path metadata

// 	output:
// 	path "*.arrow"

// 	shell:
// 	'''
// 	if grep -q "GC-Content" <(cut -f 1 !{metadata}); then
// 		csv2arrow --schema-file !{params.gisaid_schema} --header true --delimiter $'\t' \
// 		!{metadata} full_database.arrow
// 	else
// 		csv2arrow --schema-file !{params.genbank_schema} --header true --delimiter $'\t' \
// 		!{metadata} full_database.arrow
// 	fi
// 	'''

// }

process FILTER_META_TO_GEOGRAPHY {

	/*
	The geography filter in the NCBI Datasets command line
	interface appears to be broken. Until this is fixed,
	this process uses a simple julia script to filter
	the metadata down to a geography of interest.
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "lif_container"
	publishDir params.ncbi_results, mode: params.publishMode, pattern: "*.tsv"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path metadata

	output:
	path "*.arrow", emit: metadata
	path "*.txt", emit: accessions

	when:
	params.download_only == false

	script:
	"""
	filter-to-geography.jl ${metadata} ${params.geography}
	"""

}

process FILTER_SEQS_TO_GEOGRAPHY {

	/*
	This process takes the aceessions list from 
	FILTER_META_TO_GEOGRAPHY and filters down the full FASTA
	to those accessions, ensuring that both the metadata and
	the sequences reflect the same geography filtering.
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "lif_container"
	publishDir params.ncbi_results, mode: params.publishMode

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 8

	input:
	each path(fasta)
	path accessions

	output:
	path "*.fasta.zst", emit: fasta
	env count, emit: count

	script:
	"""
	seqkit grep -j ${task.cpus} -f ${accessions} ${fasta} -o filtered-to-geography.fasta.zst
	count=\$(cat ${accessions} | wc -l)
	"""

}

process REMOVE_FASTA_GAPS {

	/*
	In this process, we replace all dashes in the NCBI FASTA with
	"N's", effectively converting them into masked bases instead of
	gaps. As is the case throughout this pipeline, this process uses
	a fasta julia script to scan through each sequence.
	*/
	
	label "lif_container"

	tag "${params.pathogen}, ${params.geography}"
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	path fasta

	output:
	path "*.fasta.gz"

	when:
	params.make_distance_matrix == true

	script:
	"""
	remove-fasta-gaps.jl ${fasta}
	"""
	
}

process FILTER_BY_MASKED_BASES {

	/*
	This process makes sure that none of the sequences in the
	Genbank FASTA have fewer than 10% of their bases masked,
	i.e., are "N" instead of a defined base.
	*/

	label "lif_container"

	tag "${params.pathogen}, ${params.geography}"
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.fasta.gz"

	script:
	"""
	filter-by-n-count.jl ${fasta}
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

	label "lif_container"

	tag "${params.pathogen}, ${params.geography}"
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta
	path metadata

	output:
	path "*.fasta"

	script:
	"""
	separate-by-yearmonth.jl ${fasta} ${metadata}
	"""

}

process CLUSTER_BY_IDENTITY {

	/*
	Here we run a version of Robert Edgar's UCLUST algorithm that 
	is part of the open-source VSEARCH package. This step may take 
	a large amount of time and RAM, and thus may be better suited 
	for a high-RAM cluster environment.
	*/

	tag "${yearmonth}"
	label "lif_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.uc", emit: cluster_table
	tuple path("*centroids.fasta"), env(count), emit: centroid_fasta
	path "*-cluster-seqs*", emit: cluster_fastas
	
	script:
	yearmonth = fasta.getSimpleName()
	"""
	vsearch --cluster_fast ${fasta} \
	--id ${params.id_threshold} \
	--centroids ${yearmonth}-centroids.fasta \
	--uc ${yearmonth}-clusters.uc \
	--clusters ${yearmonth}-cluster-seqs \
	--msaout ${yearmonth}-centroids-msa.fasta \
	--threads ${task.cpus} && \
	count=\$(grep -c "^>" ${yearmonth}-centroids.fasta)
	"""
	
}

process PREP_CENTROID_FASTAS {

	/*
	In this process, we align centroid sequences to Wuhan-1, which
	will be used as an outgroup, and replace "/" symbols with underscores
	to keep iqTree happy.
	*/

	tag "${yearmonth}"
	label "lif_container"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta

	output:
	tuple path("*-aligned-centroids.fasta"), env(count)

	script:
	yearmonth = file(fasta.toString()).getSimpleName().replace("-centroids", "")
	"""
	count=\$(grep -c "^>" ${fasta})
	prep-centroid-fasta.py ${fasta} ${yearmonth}
	"""

}

process COMPUTE_DISTANCE_MATRIX {

	/*
	In parallel to clustering each month's sequences by nucleotide 
	identity, the workflow will also compute a simple nucleotide
	distance matrix, which will be used to assign distances for each
	sequence in each month.
	*/

	tag "${yearmonth}"
	label "lif_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	each path(cluster_table)
	path fasta

	output:
	path "*-dist-matrix.csv"

	when:
	file(fasta.toString()).getSimpleName().contains(file(cluster_table.toString()).getSimpleName().replace("-clusters", ""))

	script:
	yearmonth = file(cluster_table.toString()).getSimpleName().replace("-clusters", "")
	"""
	compute-distance-matrix.jl \
	${fasta} \
	${cluster_table} \
	${yearmonth} \
	${params.strictness_mode}
	"""

}

process MULTIDIMENSIONAL_SCALING {

	/*
	This plot runs multi-dimensional scaling to produce a "bee-swarm"
	plot of the sequences in each cluster. This plot will visualize 
	the cluster's distances relative to each other, and make it more
	intuitive to discern when a cluster is exceptionally distant, 
	and therefore evolutionarily advanced.
	*/

	tag "${yearmonth}"
	label "lif_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path cluster_table
	each path(centroid_fasta)

	output:
	path "*"

	when:
	file(cluster_table.toString()).getSimpleName().contains(file(centroid_fasta.toString()).getSimpleName().replace("-aligned-centroids", ""))

	script:
	yearmonth = file(cluster_table.toString()).getSimpleName().replace("-clusters", "")
	"""
	plot-mds.R ${yearmonth} ${cluster_table} ${centroid_fasta}
	"""

}

process SUMMARIZE_CANDIDATES {

	/*
	In this process, the clustering results are combined with metadata
	to provide as much information as possible with the candidate
	evolutionarily advanced sequences. We recommend users visually 
	inspect these candidates and compare them with candidates from
	the two additional methods used in this pipeline.
	*/

	label "lif_container"
	publishDir params.high_distance_candidates, mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path cluster_tables
	path cluster_fastas
	path distance_matrices
	path metadata

	output:
	path "*.tsv", emit: metadata
	tuple path("*.fasta"), env(count), emit: high_dist_seqs
	path "*.pdf", emit: plots

	script:
	"""
	report-high-dist-candidates.R ${metadata} ${params.strictness_mode} && \
	count=\$(grep -c "^>" high_distance_candidates.fasta)
	"""

}

process RUN_META_CLUSTER {

	/*
	This process runs clustering on the previous clustering
	results, but with a tighter identity threshold. In doing 
	so, the workflow can identify multiple sequences that
	are potentially from the same prolonged infection.
	*/

	label "lif_container"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.uc", emit: cluster_table
	path "*meta-centroids.fasta", emit: centroid_fasta
	path "*meta-cluster-seqs*", emit: cluster_fastas
	env repeats, emit: whether_repeats

	script:
	"""
	vsearch --cluster_fast ${fasta} \
	--id 0.9999 \
	--centroids meta-centroids.fasta \
	--uc meta-clusters.uc \
	--clusters meta-cluster-seqs \
	--threads ${task.cpus} && \
	if grep -q "H" <(cut -f 1 meta-clusters.uc); then
		repeats="true"
	else
		repeats="false"
	fi
	"""
}

process META_CLUSTER_REPORT {

	/*
	In the final step of the distance matrix branch of the 
	workflow, this process generates a CSV report highlighting
	any promising results from the meta-clustering step. If
	all high-distance sequences are in their own clusters,
	indicating that they do not stem from prolonged infections
	sampled multiple times, no report will be generated.
	*/

	label "lif_container"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	publishDir params.repeat_lineages, mode: 'copy'

	input:
	path cluster_table
	path fastas
	path metadata
	val whether_repeats

	output:
	path "*"

	when:
	whether_repeats == "true"

	script:
	"""
	report-repeat-lineages.py ${cluster_table} ${metadata}
	"""

}

process CLASSIFY_SC2_WITH_PANGOLIN {

	/*
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "lif_container"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus
	
	input:
	path fasta
	
	output:
	path "*.csv"

	when:
	params.compare_lineage_dates == true

	script:
	"""
	zstd -d `realpath ${fasta}` -o ./decompressed_genbank.fasta && \
	pangolin \
	--skip-scorpio --skip-designation-cache \
	--threads ${task.cpus} \
	--outfile lineage_report.csv \
	decompressed_genbank.fasta && \
	rm -f decompressed_genbank.fasta
	"""

}

process FIND_CANDIDATE_LINEAGES_BY_DATE {

	/*
	*/
	
	tag "${params.pathogen}, ${params.geography}"
	label "lif_container"
	publishDir params.anachronistic_candidates, mode: 'copy'
	
	errorStrategy { task.attempt < 3 ? {sleep(Math.pow(2, task.attempt) * 2000 as long); return 'retry'} : 'ignore' }
	maxRetries 2

	input:
	path lineages
	path fasta
	path metadata
	path token

	output:
	path "*.tsv", emit: metadata
	tuple path("*.fasta"), env(sample_size), emit: sequences
	path "*.pdf"

	script:
	"""
	zstd -d `realpath ${fasta}` -o ./decompressed_genbank.fasta && \
	compare-lineage-prevalences.R \
	${lineages} \
	decompressed_genbank.fasta \
	${metadata} && \
	rm -f decompressed_genbank.fasta && \
	sample_size=\$(grep -c "^>" anachronistic_seq_candidates.fasta)
	"""
}

process SEARCH_NCBI_METADATA {

	/* 
	In parallel with the distance matrix method, this pipeline also pans
	the GenBank metadata for anachronistic sequences, which may have come
	from evolutionarily advanced virus lineages in prolonged infections.
	NCBI's pre-classified pango lineages, which come with the metadata,
	are what make this method possible. In this process, we tee off this
	method by filtering NCBI metadata for all rows with correctly formatted
	dates and pango lineages.
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "lif_container"
	publishDir params.metadata_candidates, mode: 'copy'
	
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path metadata
	path lineage_dates

	output:
	path "*.tsv"

	when:
	params.inspect_ncbi_metadata == true
		
	script:
	"""
	search-ncbi-metadata.py ${metadata} ${lineage_dates} ${params.days_of_infection} ${task.cpus}
	"""

}

process FIND_DOUBLE_CANDIDATES {

	/*
	Check the outputs from the lineage and date-based approach with 
	the outputs of the distance matrix approach to find the "Venn
	overlap" between the two, where two lines of evidence suggest
	that a sequence comes from a prolonged infection within one host.
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "lif_container"
	publishDir params.double_candidates, mode: 'copy'
	
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path high_dist_metadata
	path anachronistic_meta
	path anachronistic_seqs

	output:
	path "*.tsv"
	path "*.fasta"

	when:
	(params.compare_lineage_dates == true && params.make_distance_matrix == true) || (params.make_distance_matrix == true && params.inspect_ncbi_metadata == true)

	script:
	"""
	find-double-candidates.jl ${high_dist_metadata} \
	${anachronistic_meta} \
	${anachronistic_seqs}
	"""

}

// process PREP_FOR_ESCAPE_CALC {}

// process RUN_ESCAPE_CALC {}

// process BUILD_NEXTSTRAIN_TREE {

// 	script:
// 	"""
// 	augur tree \
//     --alignment ${fasta} \
//     --output augur_output/wi_sequences_tree_raw.nwk
// 	"""
// }


// --------------------------------------------------------------- //
