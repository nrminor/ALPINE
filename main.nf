#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// prints to the screen and to the log
log.info	"""
			_______  ___      _______  ___   __    _  _______ 
			|   _   ||   |    |       ||   | |  |  | ||       |
			|  |_|  ||   |    |    _  ||   | |   |_| ||    ___|
			|       ||   |    |   |_| ||   | |       ||   |___ 
			|       ||   |___ |    ___||   | |  _    ||    ___|
			|   _   ||       ||   |    |   | | | |   ||   |___ 
			|__| |__||_______||___|    |___| |_|  |__||_______|

			ALPINE: Anachronistic Lineage and Persistent INfection Explorer
			(version v0.2.0-alpha)
			===============================================================
			fastq dir          : ${params.fasta_path}
			metadata dir       : ${params.metadata_path}
			results dir        : ${params.results}

			Filtering parameters:
			-----------------------------------
			Minimum date       : ${params.min_date}
			Maximum date       : ${params.max_date}
			Geography          : ${params.geography}
			Pathogen           : ${params.pathogen}
			ID threshold       : ${params.id_threshold}
			Ambiguity limit    : ${params.max_ambiguity}
			Stringency         : ${params.strictness_mode}
			Infection duration : ${params.days_of_infection}

			Run settings:
			-----------------------------------
			Download only?     : ${params.download_only}
			Search metadata    : ${params.search_metadata_dates}
			Pango-reclassify   : ${params.reclassify_sc2_lineages}
			Distance matrix    : ${params.make_distance_matrix}
			Debug mode         : ${params.debugmode}
			Cleanup mode       : ${params.cleanup}
			Workflow profile   : ${workflow.profile}
			Run date           : ${params.date}

			Available resources:
			-----------------------------------
			Available CPUs     : ${params.available_cpus}
			Max shared CPUs    : ${params.max_shared_cpus}
			
			"""
			.stripIndent()


// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	ch_gisaid_token = Channel
		.fromPath ( params.gisaid_token )
	
	ch_still_schemas = Channel
		.fromPath( "${params.resources}/*.schema" )
		.collect()
	
	ch_prql_query = Channel
		.fromPath( params.prql_query )

	// Data setup steps
	GET_DESIGNATION_DATES ( )

	DOWNLOAD_REFSEQ ( )

	if ( params.fasta_path == "" || params.metadata_path == "" ) {

		DOWNLOAD_NCBI_PACKAGE ( )

		UNZIP_NCBI_METADATA (
			DOWNLOAD_NCBI_PACKAGE.out.zip_archive,
			ch_prql_query
		)

		EXTRACT_NCBI_FASTA (
			DOWNLOAD_NCBI_PACKAGE.out.zip_archive
		)

		VALIDATE_METADATA (
			UNZIP_NCBI_METADATA.out
		)

		FILTER_META_TO_GEOGRAPHY (
			VALIDATE_METADATA.out
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			EXTRACT_NCBI_FASTA.out,
			FILTER_META_TO_GEOGRAPHY.out.accessions
		)

		EARLY_STATS (
			FILTER_SEQS_TO_GEOGRAPHY.out.fasta
		)

	} else {


		ch_local_fasta = Channel
			.fromPath( params.fasta_path )

		ch_local_metadata = Channel
			.fromPath( params.metadata_path )
			.map { metadata -> tuple( file(metadata), file(metadata) ) }

		VALIDATE_METADATA (
			ch_local_metadata,
			ch_still_schemas
		)

		CHECK_COLUMN_NAMES (
			VALIDATE_METADATA.out
		)

		VALIDATE_SEQUENCES (
			ch_local_fasta
		)

		FILTER_META_TO_GEOGRAPHY (
			CHECK_COLUMN_NAMES.out
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			VALIDATE_SEQUENCES.out,
			FILTER_META_TO_GEOGRAPHY.out.accessions
		)

		EARLY_STATS (
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
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta,
		DOWNLOAD_REFSEQ.out.ref_fasta
	)

	SEPARATE_BY_MONTH (
		FILTER_BY_MASKED_BASES.out
			.collectFile( name: "${params.pathogen}_prepped.fasta", newLine: true ),
		FILTER_META_TO_GEOGRAPHY.out.metadata
	)

	CLUSTER_BY_IDENTITY (
		SEPARATE_BY_MONTH.out.flatten()
	)

	COMPUTE_DISTANCE_MATRIX (
		CLUSTER_BY_IDENTITY.out.cluster_table,
		CLUSTER_BY_IDENTITY.out.centroid_fasta
			.filter { it[2].toInteger() > 2 }
			.map { fasta, yearmonth, count -> file(fasta) }
	)

	// MULTIDIMENSIONAL_SCALING (
	// 	CLUSTER_BY_IDENTITY.out.cluster_table,
	// 	PREP_CENTROID_FASTAS.out
	// 		.filter { it[1].toInteger() > 2 }
	// 		.map { fasta, count -> fasta }
	// )

	REPORT_HIGHDIST_CANDIDATES (
		CLUSTER_BY_IDENTITY.out.cluster_table.collect(),
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta,
		COMPUTE_DISTANCE_MATRIX.out.collect(),
		FILTER_META_TO_GEOGRAPHY.out.metadata
	)

	RUN_META_CLUSTER (
		REPORT_HIGHDIST_CANDIDATES.out.high_dist_seqs
	)

	META_CLUSTER_REPORT (
		RUN_META_CLUSTER.out.cluster_table,
		RUN_META_CLUSTER.out.cluster_fastas,
		REPORT_HIGHDIST_CANDIDATES.out.metadata,
		RUN_META_CLUSTER.out.whether_repeats
	)
	
	// Steps for re-running pangolin and comparing dates
	RECLASSIFY_SC2_WITH_PANGOLIN (
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta
	)

	FIND_CANDIDATE_LINEAGES_BY_DATE (
		RECLASSIFY_SC2_WITH_PANGOLIN.out
			.collectFile(name: 'new_pango_calls.csv', newLine: true),
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta,
		FILTER_META_TO_GEOGRAPHY.out.metadata,
		ch_gisaid_token
	)

	// Steps for inspecting NCBI metadata
	SEARCH_NCBI_METADATA ( 
		FILTER_META_TO_GEOGRAPHY.out.metadata,
		FILTER_SEQS_TO_GEOGRAPHY.out.fasta,
		GET_DESIGNATION_DATES.out
	)

	if ( params.make_distance_matrix == true && params.search_metadata_dates == true && params.reclassify_sc2_lineages == false ){

		FIND_DOUBLE_CANDIDATES (
			REPORT_HIGHDIST_CANDIDATES.out.metadata
				.mix(
					REPORT_HIGHDIST_CANDIDATES.out.high_dist_seqs,
					SEARCH_NCBI_METADATA.out.metadata
				).collect()
		)

		LATE_STATS (
			FIND_DOUBLE_CANDIDATES.out.fasta
		)

		COMPUTE_PREVALENCE_ESTIMATE (
			EARLY_STATS.out,
			LATE_STATS.out
		)

		COMPUTE_PREVALENCE_ESTIMATE.out.view()

		SUMMARIZE_RUN_RESULTS (
			LATE_STATS.out
		)

	} else if ( params.make_distance_matrix == true && params.search_metadata_dates == false && params.reclassify_sc2_lineages == true ){

		FIND_DOUBLE_CANDIDATES (
			REPORT_HIGHDIST_CANDIDATES.out.metadata
				.mix(
					REPORT_HIGHDIST_CANDIDATES.out.high_dist_seqs,
					FIND_CANDIDATE_LINEAGES_BY_DATE.out.metadata,
					FIND_CANDIDATE_LINEAGES_BY_DATE.out.sequences
						.filter { it[1].toInteger() > 2 }
						.map { fasta, count -> fasta }
				).collect()
		)

		LATE_STATS (
			FIND_DOUBLE_CANDIDATES.out.fasta
		)

		COMPUTE_PREVALENCE_ESTIMATE (
			EARLY_STATS.out,
			LATE_STATS.out
		)

		COMPUTE_PREVALENCE_ESTIMATE.out.view()

		SUMMARIZE_RUN_RESULTS (
			LATE_STATS.out
		)

	} else if ( params.make_distance_matrix == true && params.search_metadata_dates == true && params.reclassify_sc2_lineages == true ){

		FIND_DOUBLE_CANDIDATES (
			REPORT_HIGHDIST_CANDIDATES.out.metadata
				.mix(
					REPORT_HIGHDIST_CANDIDATES.out.high_dist_seqs,
					SEARCH_NCBI_METADATA.out.metadata,
					FIND_CANDIDATE_LINEAGES_BY_DATE.out.metadata,
					FIND_CANDIDATE_LINEAGES_BY_DATE.out.sequences
						.filter { it[1].toInteger() > 2 }
						.map { fasta, count -> fasta }
				).collect()
		)

		LATE_STATS (
			FIND_DOUBLE_CANDIDATES.out.fasta
		)

		COMPUTE_PREVALENCE_ESTIMATE (
			EARLY_STATS.out,
			LATE_STATS.out
		)

		COMPUTE_PREVALENCE_ESTIMATE.out.view()

		SUMMARIZE_RUN_RESULTS (
			LATE_STATS.out
		)

	} else if ( params.make_distance_matrix == false && params.search_metadata_dates == true && params.reclassify_sc2_lineages == true ){

		FIND_DOUBLE_CANDIDATES (
			SEARCH_NCBI_METADATA.out.metadata
				.mix(
					FIND_CANDIDATE_LINEAGES_BY_DATE.out.metadata,
					FIND_CANDIDATE_LINEAGES_BY_DATE.out.sequences
						.filter { it[1].toInteger() > 2 }
						.map { fasta, count -> fasta }
				).collect()
		)

		LATE_STATS (
			FIND_DOUBLE_CANDIDATES.out.fasta
		)

		COMPUTE_PREVALENCE_ESTIMATE (
			EARLY_STATS.out,
			LATE_STATS.out
		)

		COMPUTE_PREVALENCE_ESTIMATE.out.view()

		SUMMARIZE_RUN_RESULTS (
			LATE_STATS.out
		)

	}

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// Using debugmode setting to decide how to handle errors
if ( params.debugmode == true ){
	params.errorMode = 'terminate'
} else {
	params.errorMode = 'ignore'
}

// Defining number of cpus to use based on execution environment
if ( workflow.configFiles.contains("chtc_hpc.config") ){
	params.max_cpus = 20
} else {
	params.max_cpus = params.max_shared_cpus
}

// calling pathogen taxID if pulling SARS-CoV-2
if ( params.pathogen == "SARS-CoV-2" || params.pathogen == "sars-cov-2" ){
	pathogen = "2697049"
}

// Making a date-stamped results folder
params.dated_results = params.results + "/" + params.date

// handling the case where no geography or date filters are provided
if( params.geography == "" ){
	if ( params.fasta_path == "" ){
		params.results_subdir = params.dated_results + "/GenBank"
	} else {
		params.results_subdir = params.dated_results + "/LocalDataset"
	}
} else {
	if ( params.fasta_path == "" ){
		params.results_subdir = params.dated_results + "/GenBank_" + params.geography.replace(" ", "_")
	} else {
		params.results_subdir = params.dated_results + "/LocalDataset_" + params.geography.replace(" ", "_")
	}
}

// creating results subfolders for the three orthogonal anachronistic
// sequence search methods
params.clustering_results = params.results_subdir + "/all_clustering_results"
params.high_distance_candidates = params.results_subdir + "/high_distance_clusters"
params.repeat_lineages = params.high_distance_candidates + "/repeat_lineages"
params.anachronistic_candidates = params.results_subdir + "/anachronistic_candidates"
params.metadata_candidates = params.results_subdir + "/metadata_candidates"
params.double_candidates = params.results_subdir + "/double_candidates"

// --------------------------------------------------------------- //



// PROCESS SPECIFICATIONS
// --------------------------------------------------------------- //

process DOWNLOAD_REFSEQ {

	/*
	Here the NCBI RefSeq for the selected pathogen, if 
	available, is downloaded for downstream usage.
	*/

	tag "${params.pathogen}"
	publishDir params.resources, mode: 'copy'

	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 5

	output:
	path "${params.pathogen}_refseq.fasta", emit: ref_fasta

	script:
	"""
	datasets download virus genome taxon ${pathogen} \
	--refseq && \
	unzip ncbi_dataset.zip
	mv ncbi_dataset/data/genomic.fna ./${params.pathogen}_refseq.fasta
	"""
}

process GET_DESIGNATION_DATES {
	
	/*
	This process downloads a table of pangolin lineage designation dates
	from Cornelius Roemer's GitHub. These dates represent when each lineage was
	added to pangolin, after which point sequences could be classified as such 
	*/

	label "alpine_container"
	publishDir params.resources, mode: 'copy', overwrite: true
	
	output:
	path "lineage_designation_dates.csv"

	when:
	params.pathogen == "SARS-CoV-2"
	
	script:
	"""
	curl -fsSL ${params.lineage_dates} > lineage_designation_dates.csv
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
	// publishDir params.results_subdir, mode: 'copy'

	errorStrategy { sleep(Math.pow(2, task.attempt) * 1000 as long); return 'retry' }
	maxRetries 3

	output:
	path "ncbi_dataset.zip", emit: zip_archive

	script:
	"""
	datasets download virus genome taxon ${pathogen}
	"""	

}

process UNZIP_NCBI_METADATA {

	/*
	Here the pathogen metadata in TSV format is extracted from the 
	lightweight NCBI zip archive.
	*/
	
	label "alpine_container"
	tag "${params.pathogen}"
	publishDir params.results_subdir, pattern: "*.tsv", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	path zip
	path query

	output:
	tuple path("genbank_metadata.csv"), path("genbank_metadata.arrow")

	script:
	"""
	prqlc compile ${query} > query.sql && \
	unzip -p ${zip} ncbi_dataset/data/data_report.jsonl \
	| dataformat tsv virus-genome --force \
	| qsv fmt --delimiter "\t" --out-delimiter , -o genbank_metadata.csv

	qsv index genbank_metadata.csv
	qsv validate --invalid invalid_accessions.tsv --jobs 8 genbank_metadata.csv > validation_report.txt

	qsv sqlp genbank_metadata.csv query.sql \
	--try-parsedates --date-format '%Y-%m-%d' --ignore-errors \
	--format arrow --compression 'zstd' \
	--output genbank_metadata.arrow
	"""

}

process EXTRACT_NCBI_FASTA {

	/*
	Here the pathogen sequences in FASTA format are extracted 
	from the lightweight NCBI zip archive.
	*/

	tag "${params.pathogen}"

	label "alpine_container"
	publishDir params.results_subdir, mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 3

	input:
	path zip

	output:
	path "genbank_sequences.fasta.zst"

	script:
	"""
	unzip -p ${zip} ncbi_dataset/data/genomic.fna \
	| zstd -10 -o genbank_sequences.fasta.zst
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

	label "alpine_container"

	input:
	tuple path(raw_text), path(arrow)
	each path(still_schemas)

	output:
	tuple path(arrow), env(db)

	when:
	params.download_only == false

	script:
	"""
	if head -n 1 ${metadata} | grep -q "Virus name"; then
		db="GISAID"
		still validate gisaid.schema ${raw_text}
	elif head -n 1 ${metadata} | grep -q ^"Virus name"; then
		db="Genbank"
		still validate genbank.schema ${raw_text}
	else
		echo "Database source not recognized and thus column names and types cannot be ascertained."
		exit 1
	fi
	"""
}

process CHECK_COLUMN_NAMES {

	/* */

	tag "${params.pathogen}"

	label "alpine_container"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	input:
	tuple path(metadata), val(db)

	output:
	path "validated.arrow"

	script:
	"""
	check-columns.py ${metadata}
	"""

}

process VALIDATE_SEQUENCES {

	/*
	This step quickly runs through the fasta to make sure that
	the accession identifier for each record is accessible when
	the defline is parsed by space (" "). This is mostly because
	GISAID delimits its deflines with a pipe symbol ("|"). This
	step does not run when data are coming directly from NCBI.
	*/

	tag "${params.pathogen}"

	label "alpine_container"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "validated.fasta.zst"

	when:
	params.download_only == false

	script:
	"""
	seqkit replace \
	-j ${task.cpus} --f-by-name --keep-untouch \
	--pattern "\\|" --replacement " " \
	${fasta} -o validated.fasta.zst
	"""

}

process FILTER_META_TO_GEOGRAPHY {

	/*
	The geography filter in the NCBI Datasets command line
	interface appears to be broken. Until this is fixed,
	this process uses a simple julia script to filter
	the metadata down to a geography of interest.
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "alpine_container"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	input:
	path metadata

	output:
	path "filtered-to-geography.arrow", emit: metadata
	path "accessions.txt", emit: accessions

	script:
	"""
	filter-to-geography.py \
	--metadata ${metadata} \
	--max_date ${params.max_date} \
	--min_date ${params.min_date} \
	--geography ${params.geography}
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
	label "alpine_container"
	publishDir params.results_subdir, mode: 'copy'

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus

	input:
	path fasta
	path accessions

	output:
	path "filtered-to-geography.fasta.zst", emit: fasta

	script:
	"""
	seqkit grep -j ${task.cpus} -f ${accessions} ${fasta} -o filtered-to-geography.fasta.zst
	"""

}

process EARLY_STATS {

	/*
	This process computes a variety of sequence statistics for
	the sequences flagged as double candidates, with the most
	important being the number of sequences. This number will
	be used to compute the prevalence of these high-distance,
	anachronistic sequences among the input sequences. The
	input sequence database is surely a biased sample, but we
	attempt to use it here as a proxy for the prevalence of
	these pathogens in the general population.
	*/
	
	tag "${params.pathogen}, ${params.geography}"
	label "alpine_container"
	publishDir params.results_subdir, mode: 'copy'

	errorStrategy { task.attempt < 1 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "early_stats.tsv"

	script:
	"""
	seqkit stats -j ${task.cpus} --tabular ${fasta} > early_stats.tsv
	"""

}

process FILTER_BY_MASKED_BASES {

	/*
	This process makes sure that none of the sequences in the
	Genbank FASTA have fewer than 10% of their bases masked,
	i.e., are "N" instead of a defined base.
	*/

	label "alpine_container"

	tag "${params.pathogen}, ${params.geography}"
	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	input:
	path fasta
	path reference

	output:
	path "filtered-by-n.fasta.gz"

	when:
	params.make_distance_matrix == true

	script:
	"""
	alpine filter-by-n \
	--fasta ${fasta} \
	--ambiguity ${params.max_ambiguity} \
	--reference ${reference}
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

	label "alpine_container"

	tag "${params.pathogen}, ${params.geography}"
	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	input:
	path fasta
	path metadata

	output:
	path "*.fasta"

	script:
	"""
	alpine separate-by-month --fasta ${fasta} --metadata ${metadata}
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
	label "alpine_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy', pattern: "*.uc"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy', pattern: "*-cluster-seqs*"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "${yearmonth}-clusters.uc", emit: cluster_table
	tuple path("${yearmonth}-centroids.fasta"), val(yearmonth), env(count), emit: centroid_fasta
	path "${yearmonth}-cluster-seqs*", emit: cluster_fastas
	
	script:
	yearmonth = fasta.getSimpleName()
	"""
	vsearch --cluster_fast ${fasta} \
	--id ${params.id_threshold} \
	--centroids ${yearmonth}-centroids.fasta \
	--uc ${yearmonth}-clusters.uc \
	--clusters ${yearmonth}-cluster-seqs \
	--threads ${task.cpus} && \
	count=\$(grep -c "^>" ${yearmonth}-centroids.fasta)
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
	label "alpine_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy 'retry' // { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus 1

	input:
	each path(cluster_table)
	path fasta

	output:
	path "${yearmonth}-dist-matrix.csv"

	when:
	file(fasta.toString()).getSimpleName().contains(file(cluster_table.toString()).getSimpleName().replace("-clusters", ""))

	script:
	yearmonth = file(cluster_table.toString()).getSimpleName().replace("-clusters", "")
	"""
	ln -s /opt/.julia/alpine.so alpine.so
	compute-distance-matrix.jl ${fasta} ${cluster_table} ${yearmonth} ${params.strictness_mode}
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
	label "alpine_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

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

process REPORT_HIGHDIST_CANDIDATES {

	/*
	In this process, the clustering results are combined with metadata
	to provide as much information as possible with the candidate
	evolutionarily advanced sequences. We recommend users visually 
	inspect these candidates and compare them with candidates from
	the two additional methods used in this pipeline.
	*/

	label "alpine_container"
	publishDir params.high_distance_candidates, mode: 'copy'

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	input:
	path cluster_tables
	path fasta
	path distance_matrices
	path metadata

	output:
	path "high_distance_candidates.tsv", emit: metadata
	path "high_distance_candidates.fasta", emit: high_dist_seqs
	// path "*.pdf", emit: plots

	script:
	"""
	report-high-dist-candidates.py \
	--sequences ${fasta} \
	--stringency ${params.strictness_mode}
	"""

}

process RUN_META_CLUSTER {

	/*
	This process runs clustering on the previous clustering
	results, but with a tighter identity threshold. In doing 
	so, the workflow can identify multiple sequences that
	are potentially from the same prolonged infection.
	*/

	label "alpine_container"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

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
	--id 0.9995 \
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

	label "alpine_container"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	publishDir params.repeat_lineages, mode: 'copy'

	input:
	path cluster_table
	path fastas
	path metadata
	val whether_repeats

	output:
	path "repeat-lineage*"

	when:
	whether_repeats == "true"

	script:
	"""
	report-repeat-lineages.py ${cluster_table} ${metadata}
	"""

}

process RECLASSIFY_SC2_WITH_PANGOLIN {

	/*
	Here, the workflow runs the filtered Genbank FASTA through
	the Pangolin tool with its fastest settings. These updated 
	pango lineage classifications are then used in conjunction
	with the Outbreak.info API to identify anachronistic 
	sequences.
	*/

	tag "${params.pathogen}, ${params.geography}"
	label "alpine_container"

	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus
	
	input:
	path fasta
	
	output:
	path "*.csv"

	when:
	params.reclassify_sc2_lineages == true && params.pathogen == "SARS-CoV-2"

	script:
	"""
	pangolin --update --update-data && \
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
	Here the workflow uses a GISAID authentication token,
	which must be provided by the user, to access lineage prevalence
	estimates via the Outbreak.info API. Lineage collection dates
	that are far past what is called the "rarity date" in the script
	are flagged as candidate anachronistics.

	Note that the Outbreak.info is currently quite finicky, and
	may boot out the script for making too many API calls. This
	process will dynamically retry with longer and longer wait times
	to address this, but in some cases the API may cause this step to
	fail. In such cases, we recommend 
	*/
	
	tag "${params.pathogen}, ${params.geography}"
	label "alpine_container"
	publishDir params.anachronistic_candidates, mode: 'copy'
	
	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path lineages
	path fasta
	path metadata
	path token

	output:
	path "*.tsv", emit: metadata
	path "*.fasta", emit: sequences
	path "*.pdf"

	script:
	"""
	compare-lineage-prevalences.R \
	${lineages} \
	${metadata} \
	${fasta}
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
	label "alpine_container"
	publishDir params.metadata_candidates, mode: 'copy'
	
	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus

	input:
	path metadata
	path fasta
	path lineage_dates

	output:
	path "anachronistic_metadata_only_candidates.tsv", emit: metadata
	path "anachronistic_metadata_only_candidates.fasta", emit: fasta

	when:
	params.search_metadata_dates == true && params.pathogen == "SARS-CoV-2"
		
	script:
	"""
	search-ncbi-metadata.py ${metadata} ${lineage_dates} ${params.days_of_infection} ${task.cpus} && \
	cut -f 1 anachronistic_metadata_only_candidates.tsv \
	| tail -n +2 > anachronistic_accessions.txt && \
	seqkit grep -j ${task.cpus} -nri -f anachronistic_accessions.txt ${fasta} \
	-o anachronistic_metadata_only_candidates.fasta
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
	label "alpine_container"
	publishDir params.double_candidates, mode: 'copy'
	
	errorStrategy { task.attempt < 2 ? 'retry' : params.errorMode }
	maxRetries 1

	input:
	path collected_files

	output:
	path "double_candidate*.fasta", emit: fasta
	path "double_candidate*.tsv"

	when:
	collected_files.toString().count(".tsv") > 1

	script:
	"""
	ln -s /opt/.julia/alpine.so alpine.so
	find-double-candidates.jl
	"""

}

process LATE_STATS {

	/*
	This process computes the same statistics post-run as the
	process `EARLY_STATS` above. Refer to the docstring there
	for more explanation and caveats.
	*/
	
	tag "${params.pathogen}, ${params.geography}"
	label "alpine_container"
	publishDir params.results_subdir, mode: 'copy'

	errorStrategy { task.attempt < 1 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "late_stats.tsv"

	script:
	"""
	seqkit stats -j ${task.cpus} --tabular ${fasta} > late_stats.tsv
	"""

}

process COMPUTE_PREVALENCE_ESTIMATE {

	/*
	Here the workflow compares the sequence statistics
	from above to estimate and print out a prevalence estimate for
	the double-candidates identified above.
	*/
	
	tag "${params.pathogen}, ${params.geography}"
	label "alpine_container"
	publishDir params.results_subdir, mode: 'copy'

	errorStrategy { task.attempt < 1 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus 1

	input:
	path early_stats
	path late_stats

	output:
	stdout

	script:
	"""
	ln -s /opt/.julia/alpine.so alpine.so
	estimate-prevalence.jl ${early_stats} ${late_stats}
	"""

}

process SUMMARIZE_RUN_RESULTS {

	/*
	As a final step, the workflow walks through the geographies
	searched on this run-date and summarizes their findings in
	a readable spreadsheet.
	*/
	
	tag "${params.pathogen}, ${params.date}"
	label "alpine_container"
	publishDir params.dated_results, mode: 'copy'

	errorStrategy { task.attempt < 1 ? 'retry' : params.errorMode }
	maxRetries 1

	cpus 1

	input:
	path tsv

	output:
	path "*.tsv"

	script:
	"""
	summarize-run-results.R ${params.dated_results}
	"""

}

// --------------------------------------------------------------- //
