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
	// DOWNLOAD_REFSEQ ( )

	println "This run will search NCBI GenBank data for the pathogen ${params.pathogen}"
	
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
			UNZIP_NCBI_FASTA.out,
			FILTER_TSV_TO_GEOGRAPHY.out.accessions
		)

		if ( workflow.profile == 'standard' || workflow.profile == 'docker' ){

			REMOVE_FASTA_GAPS ( 
				FILTER_SEQS_TO_GEOGRAPHY.out.fasta
			)

		} else {

			REMOVE_FASTA_GAPS ( 
				FILTER_SEQS_TO_GEOGRAPHY.out.fasta
					.splitFasta( by: 5000, file: "genbank-${params.pathogen}.fasta" )
			)
				
		}

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

		FILTER_TSV_TO_GEOGRAPHY (
			ch_local_metadata
		)

		FILTER_SEQS_TO_GEOGRAPHY (
			ch_local_fasta,
			FILTER_TSV_TO_GEOGRAPHY.out.accessions
		)

		REMOVE_FASTA_GAPS ( 
			FILTER_SEQS_TO_GEOGRAPHY.out.fasta
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
	FILTER_BY_MASKED_BASES (
		REMOVE_FASTA_GAPS.out
	)

	SEPARATE_BY_MONTH (
		FILTER_BY_MASKED_BASES.out
			.collectFile( name: "${params.pathogen}_prepped.fasta", newLine: true ),
		FILTER_TSV_TO_GEOGRAPHY.out.metadata
	)

	CLUSTER_BY_IDENTITY (
		SEPARATE_BY_MONTH.out.flatten()
	)

	COUNT_FASTA_RECORDS (
		CLUSTER_BY_IDENTITY.out.centroid_fasta
	)

	PREP_CENTROID_FASTAS (
		COUNT_FASTA_RECORDS.out
			.filter { it[1].toInteger() > 1 }
			.map { fasta, count -> fasta }
	)

	FIND_MAJORITY_CLUSTER (
		CLUSTER_BY_IDENTITY.out.cluster_table,
		PREP_CENTROID_FASTAS.out
	)

	COMPUTE_DISTANCE_MATRIX (
		FIND_MAJORITY_CLUSTER.out
	)

	// BUILD_CENTROID_TREE (
	// 	PREP_CENTROID_FASTAS.out
	// 		.filter { it[1].toInteger() > 3 }
	// 		.map { fasta, count -> fasta }
	// )

	MDS_PLOT (
		CLUSTER_BY_IDENTITY.out.cluster_table,
		PREP_CENTROID_FASTAS.out
			.filter { it[1].toInteger() > 2 }
			.map { fasta, count -> fasta }
	)

	// PLOT_TREE (
	// 	BUILD_CENTROID_TREE.out,
	// 	DOWNLOAD_REFSEQ.out.ref_id
	// )

	GENERATE_CLUSTER_REPORT (
		CLUSTER_BY_IDENTITY.out.cluster_table.collect(),
		CLUSTER_BY_IDENTITY.out.cluster_fastas.collect(),
		COMPUTE_DISTANCE_MATRIX.out.collect(),
		FILTER_TSV_TO_GEOGRAPHY.out.metadata
	)

	RUN_META_CLUSTER (
		GENERATE_CLUSTER_REPORT.out.high_dist_seqs
			.filter { it[1].toInteger() > 2 }
			.map { fasta, count -> fasta }
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
	// 	FILTER_SEQS_TO_GEOGRAPHY.out.fasta
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

if ( params.pathogen == "SARS-CoV-2" || params.pathogen == "sars-cov-2" ){
	pathogen = "2697049"
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
params.high_distance_candidates = params.ncbi_results + "/high_distance_clusters"
params.repeat_lineages = params.high_distance_candidates + "/repeat_lineages"
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
	(workflow.profile == 'standard' || workflow.profile == 'docker' || workflow.profile == 'singularity') && (params.pathogen == "SARS-CoV-2" || params.pathogen == "sars-cov-2")
	
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

	publishDir params.ncbi_results, mode: params.publishMode

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

process UNZIP_NCBI_FASTA {

	/*
	Here the pathogen sequences in FASTA format are extracted 
	from the lightweight NCBI zip archive.
	*/

	publishDir params.ncbi_results, mode: params.publishMode

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

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

	label "lif_container"
	publishDir params.ncbi_results, mode: params.publishMode, pattern: "*.tsv"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	path metadata

	output:
	path "*.tsv", emit: metadata
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
	FILTER_TSV_TO_GEOGRAPHY and filters down the full FASTA
	to those accessions, ensuring that both the metadata and
	the sequences reflect the same geography filtering.

	This process is parallelized, meaning that it splits up
	the large genbank FASTA into many more computationally
	bite-sized pieces of ~5,000 sequences.
	*/

	label "lif_container"
	publishDir params.ncbi_results, mode: params.publishMode

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	each path(fasta)
	path accessions

	output:
	path "*.fasta", emit: fasta
	env count, emit: count

	script:
	"""
	subseq_rs ${fasta} ${accessions} && \
	count=\$(grep -c "^>" filtered-to-geography.fasta)
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

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 1

	input:
	path fasta

	output:
	path "*.fasta"

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

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.fasta"

	script:
	"""
	JULIA_NUM_THREADS=${task.cpus} \
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
	JULIA_NUM_THREADS=${task.cpus} \
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
	path "*centroids.fasta", emit: centroid_fasta
	path "*-cluster-seqs*", emit: cluster_fastas
	
	script:
	yearmonth = fasta.getSimpleName()
	"""
	vsearch --cluster_fast ${fasta} \
	--id ${params.id_threshold} \
	--centroids ${yearmonth}-centroids.fasta \
	--uc ${yearmonth}-clusters.uc \
	--clusters ${yearmonth}-cluster-seqs \
	--threads ${task.cpus}
	"""
	
}

process COUNT_FASTA_RECORDS {

	/*
	This process counts the number of records in each year-month's
	centroid FASTA so it can filter low-cluster months prior to
	iqTree.
	*/

	tag "${yearmonth}"
	label "lif_container"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path fasta

	output:
	tuple path(fasta), env(count)

	shell:
	yearmonth = file(fasta.toString()).getSimpleName().replace("-centroids", "")
	'''
	count=$(grep -c "^>" !{fasta})
	echo "There are" ${count} "records in the FASTA" !{fasta}
	'''

}

process PREP_CENTROID_FASTAS {

	/*
	In this process, we align centroid sequences to Wuhan-1, which
	will be used as an outgroup, and replace "/" symbols with underscores
	to keep iqTree happy.
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
	tuple path("*-aligned-centroids.fasta"), env(count)

	script:
	yearmonth = file(fasta.toString()).getSimpleName().replace("-centroids", "")
	"""
	count=\$(grep -c "^>" ${fasta})
	prep-centroid-fasta.py ${fasta} ${yearmonth}
	"""

}

process FIND_MAJORITY_CLUSTER {

	/*
	This process replaces the previous "prep_centroid_fasta" approach that brought in
	the Wuhan-1 sequence with an approach that uses the biggest cluster, representing
	the majority of sequences from each month, as a sort of "reference." If there is
	just one additional cluster, the downstream distance matrix step will create
	a single distance and assign it to the small cluster as a potential 
	evolutionarily advanced lineage.
	*/

	tag "${yearmonth}"
	label "lif_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	each path(cluster_table)
	tuple path(fasta), val(count)

	output:
	tuple path(fasta), path(cluster_table), val(yearmonth), env(majority_cluster), env(majority_centroid), env(cluster_count)

	when:
	file(fasta.toString()).getSimpleName().contains(file(cluster_table.toString()).getSimpleName().replace("-clusters", ""))

	script:
	yearmonth = file(cluster_table.toString()).getSimpleName().replace("-clusters", "")
	"""
	find-majority-cluster.py ${cluster_table}
	cluster_count=`cat cluster_count.txt`
	majority_cluster=`cat majority_cluster.txt`
	majority_centroid=`cat majority_centroid.txt`
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
	tuple path(fasta), path(cluster_table), val(yearmonth), val(majority_cluster), val(majority_centroid), val(cluster_count)

	output:
	path "*-dist-matrix.csv"

	script:
	"""
	compute-distance-matrix.jl \
	${fasta} \
	${cluster_table} \
	${yearmonth} ${cluster_count} ${majority_centroid}
	"""

}

process BUILD_CENTROID_TREE {

	/*
	Next, we take the centroid sequences aligned with Wuhan-1
	and build a tree with Wuhan-1 as the outgroup. The goal is
	to identify which centroid has the longest branch length 
	for each year-month. The sequences that cluster with that
	long-branch centroid will be classified as a evolutionarily
	advanced.
	*/

	tag "${yearmonth}"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path fasta

	output:
	path "*.treefile"

	script:
	yearmonth = file(fasta.toString()).getSimpleName().replace("-aligned-centroids", "")
	"""
	iqtree -s ${fasta} -pre ${yearmonth} -m MFP -bb 1000 -nt ${task.cpus}
	"""

}

process MDS_PLOT {

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

process PLOT_TREE {

	/*
	This process builds a simple phylogeny plot to show how the 
	centroid sequences from each month are related to each other
	and to the pathogen RefSeq.
	*/

	tag "${yearmonth}"
	label "lif_container"
	publishDir "${params.clustering_results}/${yearmonth}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path treefile
	val ref_id

	output:
	path "*.pdf"

	script:
	yearmonth = file(treefile.toString()).getSimpleName().replace(".treefile", "")
	"""
	plot-tree.R ${treefile} ${yearmonth} ${ref_id}
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
	generate-cluster-report.R ${metadata} && \
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
	publishDir params.repeat_lineages, mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus params.max_cpus

	input:
	path fasta

	output:
	path "*.uc", emit: cluster_table
	path "*meta-centroids.fasta", emit: centroid_fasta
	path "*meta-cluster-seqs*", emit: cluster_fastas

	script:
	"""
	vsearch --cluster_fast ${fasta} \
	--id 0.9999 \
	--centroids meta-centroids.fasta \
	--uc meta-clusters.uc \
	--clusters meta-cluster-seqs \
	--threads ${task.cpus}
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
