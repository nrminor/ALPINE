#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
if( params.local_database == "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data"){
	params.local_input_path = params.local_database + "/DHO*/gisaid/*.fasta"
} else if( params.local_database.isEmpty() ){
	params.local_input_path = launchDir + "/*.fasta"
} else {
	params.local_input_path = params.local_database + "/**/*.fasta"
}

if( params.local_metadata.isEmpty() ){
	params.local_metadata = launchDir
}

if( params.gisaid_seq_dir.isEmpty() ){
	params.gisaid_seq_dir = launchDir
}

if( params.gisaid_metadata_dir.isEmpty() ){
	params.gisaid_metadata_dir = launchDir
}

params.ncbi_results = params.results + "/GenBank_" + params.geography + "_" + params.min_date + "_to_" + params.max_date
params.gisaid_results = params.results + "/GISAID_" + params.geography + "_" + params.min_date + "_to_" + params.max_date
params.local_results = params.results + "/local_database_" + params.geography + "_" + params.min_date + "_to_" + params.max_date
// --------------------------------------------------------------- //



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels for consensus sequences 
	ch_local_seqs = Channel
		.fromPath( params.local_input_path )
		.filter { !it.toString().contains(" copy") }
		.map { fasta -> tuple( file(fasta), fasta.getParent(), fasta.getSimpleName() ) }
	
	ch_gisaid_seqs = Channel
		.fromPath( "${params.gisaid_seq_dir}/*sequences*.fasta" )
		.splitFasta( by: 1 )
	
	ch_gisaid_metadata = Channel
		.fromPath( "${params.gisaid_metadata_dir}/*metadata.tsv" )
	
	

	// Before anything else, make sure pangolin is up to date
	UPDATE_PANGO_CONTAINER ( )
	
	println "Pangolin updated to version:"
	UPDATE_PANGO_CONTAINER.out.cue.view()
	

	// Download the latest Pangolin designation dates from GitHub
	GET_DESIGNATION_DATES ( )


	// NCBI/GENBANK BRANCH:
	// Pull NCBI metadata and, if params.search_genbank_seqs = true,
	// download sequences as well and use those to find long
	// infection candidates
	PULL_NCBI_METADATA ( )

	FILTER_NCBI_METADATA (
		PULL_NCBI_METADATA.out
	)

	PULL_NCBI_SEQUENCES ( 
		FILTER_NCBI_METADATA.out
			.splitCsv ( header: true, sep: "\t" )
			.map {row -> tuple(row.accession, row.date, row.location, row.pango)}
	)
	
	RECLASSIFY_NCBI_SEQUENCES ( 
		UPDATE_PANGO_CONTAINER.out.cue,
		PULL_NCBI_SEQUENCES.out
	)

	FIND_NCBI_LONG_INFECTIONS ( 
		GET_DESIGNATION_DATES.out,
		RECLASSIFY_NCBI_SEQUENCES.out
	)

	CONCAT_NCBI_LONG_INFECTIONS (
		FIND_NCBI_LONG_INFECTIONS.out.collect()
	)

	SEARCH_NCBI_METADATA ( 
		FILTER_NCBI_METADATA.out,
		GET_DESIGNATION_DATES.out
	)

	PULL_NCBI_CANDIDATES (
		SEARCH_NCBI_METADATA.out
			.splitCsv( header: true, sep: "\t" )
			.map { row -> row.accession }
	)

	CONCAT_NCBI_CANDIDATES (
		PULL_NCBI_CANDIDATES.out.collect()
	)


	// LOCAL DATABASE BRANCH:
	// Search for long infection candidates in FASTAs from a local
	// database. 
	// RECLASSIFY_LOCAL_SEQS (
	// 	UPDATE_PANGO_CONTAINER.out.cue,
	// 	ch_local_seqs
	// )

	// FIND_LOCAL_LONG_INFECTIONS (
	// 	GET_DESIGNATION_DATES.out,
	// 	RECLASSIFY_LOCAL_SEQS.out,
	// )
	
	// CONCAT_LOCAL_LONG_INFECTIONS (
	// 	FIND_LOCAL_LONG_INFECTIONS.out.collect()
	// )

	// SEARCH_DHOLAB_METADATA ( 
	// 	GET_DESIGNATION_DATES.out
	// )


	// GISAID BRANCH:
	// This branch of the workflow filters the full GISAID metadata 
	// down to the date range and geography specified in nextflow.config.
	// It then reclassify GISAID EpiCov FASTA sequences with pangolin, 
	// if they have been made available with params.gisaid_seq_dir
	FILTER_GISAID_METADATA (
		ch_gisaid_metadata
	)

	DATE_FASTAS (
		ch_gisaid_seqs
	)

	RECLASSIFY_GISAID_SEQS (
		UPDATE_PANGO_CONTAINER.out.cue,
		DATE_FASTAS.out
			.splitCsv( header: false )
			.map { row -> tuple( file(row[0]), row[1], row[2] ) }
	)

	FIND_GISAID_LONG_INFECTIONS ( 
		GET_DESIGNATION_DATES.out,
		RECLASSIFY_GISAID_SEQS.out
	)

	CONCAT_GISAID_LONG_INFECTIONS (
		FIND_GISAID_LONG_INFECTIONS.out.collect()
	)
	
	SEARCH_GISAID_METADATA (
		FILTER_GISAID_METADATA.out,
		GET_DESIGNATION_DATES.out
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
process PULL_NCBI_METADATA {

	publishDir params.ncbi_results, mode: 'copy'

	output:
	path "*.tsv"

	script:
	"""
	datasets summary virus genome taxon sars-cov-2 \
	--as-json-lines | dataformat tsv virus-genome \
	--fields accession,isolate-lineage,isolate-lineage-source,geo-location,geo-region,isolate-collection-date,virus-pangolin,virus-name,sra-accs,submitter-affiliation,submitter-country,submitter-names,update-date \
	> sarscov2-metadata.tsv
	"""

}

process FILTER_NCBI_METADATA {

	input:
	path tsv

	output:
	path "*.tsv"

	script:
	"""
	filter_ncbi_metadata.R ${tsv} ${params.min_date} ${params.max_date} ${params.geography}
	"""
}

process PULL_NCBI_SEQUENCES {
	
	tag "${accession}"
	
	cpus 1
	
	input:
	tuple val(accession), val(date), val(loc), val(pango)
	
	output:
	tuple val(accession), val(date), path("*.fasta")
	
	when:
	params.search_genbank_seqs == true

	script:
	"""
	datasets download virus genome accession "${accession}" \
	--exclude-cds --exclude-protein
	unzip ncbi_dataset.zip
	mv ncbi_dataset/data/genomic.fna ./"${accession}".fasta
	rm -rf ncbi_dataset/
	"""

}

process RECLASSIFY_NCBI_SEQUENCES {
	
	tag "${accession}"
	
	cpus 1
	time { 5.minutes * task.attempt }
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	each cue
	tuple val(accession), val(date), path(fasta)
	
	output:
	tuple path("*.csv"), val(accession), val(date)

	script:
	"""
	pangolin \
	--threads ${task.cpus} \
	--outfile ${accession}_lineages_updated_${date}.csv \
	"${fasta}"
	"""

}

process FIND_NCBI_LONG_INFECTIONS {

	publishDir params.ncbi_results, mode: 'copy'

	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(accession), val(date)

	output:
	path "*putative_long_infections_ncbi*.csv"

	script:
	"""
	ncbi_long_infection_finder.R ${lineage_dates} ${lineage_csv} ${date} ${params.days_of_infection}
	"""
}

process CONCAT_NCBI_LONG_INFECTIONS {

	publishDir params.results, mode: 'copy'

	input:
	path file_list, stageAs: 'infections??.csv'

	output:
	path "*.csv"

	script:
	"""
	concat_long_infections.R ${params.days_of_infection}
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
	params.search_genbank_metadata == true
		
	script:
	"""
	search_ncbi_metadata.R ${metadata} ${lineage_dates} ${params.days_of_infection} ${task.cpus}
	"""

}

process PULL_NCBI_CANDIDATES {

	tag "${accession}"

	errorStrategy 'retry'
	maxRetries 4

	input:
	val accession

	output:
	path "*.fasta"

	script:
	"""
	datasets download virus genome accession "${accession}" && \
	unzip ncbi_dataset.zip
	mv ncbi_dataset/data/genomic.fna ./"${accession}".fasta
	rm -rf ncbi_dataset/
	"""

}

process CONCAT_NCBI_CANDIDATES {

	publishDir params.ncbi_results, mode: 'copy'

	input:
	path fasta_list

	output:
	path "*.fasta"

	script:
	"""
	find . -name "*.fasta" > fasta.list

	for i in `cat fasta.list`;
	do
		cat \$i >> ncbi_long_infection_candidates.fasta
	done
	"""
}


// LOCAL DATABASE PROCESSES
process RECLASSIFY_LOCAL_SEQS {
	
	// This process identifies lineages for all samples in every local
	// sequencing run. It pulls FASTAs directly from the supplied path to do so.
	
	tag "${experiment_number}"
	// publishDir params.lineage_reports, pattern: '*.csv', mode: 'copy'
	
	cpus 1
	time { 5.minutes * task.attempt }
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	each cue
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	path("*.csv"), val(run_name), val(run_dir), val(experiment_number), env(experiment_date)
	
	script:
	experiment_number = 'DHO_' + parentdir.toString().split("DHO_")[1].replaceAll('/gisaid','')
	"""
	experiment_date=`date -r ${fasta} "+%Y-%m-%d"`
	
	pangolin \
	--threads ${task.cpus} \
	--outfile ${experiment_number}_lineages_updated_${params.date}.csv \
	"${fasta}"
	"""
	
}

process FIND_LOCAL_LONG_INFECTIONS {
	
	// This process calls a script that identifies all samples that classify as "old"
	// lineages, i.e. lineages that were first designated in pangolin long before a 
	// sequence classified as it in a sequencing run. By default, this amount of time is
	// eight months, but this is subject to change in the future. In such cases, where
	// a lineage appears long after it arose and subsided, it is more likely that the 
	// infection individual has sustained a prolonged infection since that lineage 
	// was prevalent, and less likely that the old lineage re-appeared despite competition
	// from newer, more fit lineages.
	
	tag "${experiment_number}"
	
	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(run_name), val(run_dir), val(experiment_number), val(experiment_date)
	
	output:
	path "*putative_long_infections*.csv"
	
	script:
	run = run_name.toString().replaceAll(" copy", "")
	
	"""
	local_long_infection_finder.R ${run} ${experiment_number} ${experiment_date} ${lineage_csv} ${lineage_dates} ${params.days_of_infection}
	"""
}

process CONCAT_LOCAL_LONG_INFECTIONS {
	
	// This step simply concatenates any putative prolonged infection samples into on
	// table, which is then exported into the results directory.
	
	publishDir params.results, mode: 'copy'
	
	input:
	path file_list, stageAs: 'infections??.csv'
	
	output:
	path "*.csv"
	
	script:
	"""
	concat_long_infections.R "${params.days_of_infection}"
	"""
	
}

process SEARCH_DHOLAB_METADATA {

	publishDir params.local_results, mode: 'copy'

	input:
	path lineage_dates

	output:
	path "*.csv"

	when:
	params.search_metadata == true
		
	script:
	"""
	search_local_metadata.R
	"""
	
}


// GISAID PROCESSES
process FILTER_GISAID_METADATA {

	input:
	path tsv

	output:
	path "*.tsv"

	when:
	params.search_gisaid_metadata == true

	script:
	"""
	cut -f 1,5,6,7,14,15 ${tsv} > gisaid_metadata_reduced.tsv && \
	filter_gisaid_metadata.R gisaid_metadata_reduced.tsv \
	${params.min_date} ${params.max_date} ${params.geography}
	"""

}

process DATE_FASTAS {

	input:
	path fasta

	output:
	path "*.csv"

	script:
	"""
	date_fastas.R ${fasta}
	"""
	
}

process RECLASSIFY_GISAID_SEQS {
	
	cpus 1
	time { 5.minutes * task.attempt }
	errorStrategy 'retry'
	maxRetries 4

	input:
	each cue
	tuple path(fasta), val(accession), val(date)

	output:
	tuple path("*.csv"), val(accession), val(date)

	when:
	params.search_gisaid_seqs == true

	script:
	"""
	pangolin \
	--threads ${task.cpus} \
	--outfile ${accession}_lineages_updated_${params.date}.csv \
	"${fasta}"
	"""

}

process FIND_GISAID_LONG_INFECTIONS {

	publishDir params.gisaid_results, mode: 'copy'

	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(accession), val(date)

	output:
	path "*putative_long_infections_gisaid*.tsv"

	script:
	"""
	gisaid_long_infection_finder.R ${lineage_dates} ${lineage_csv} ${date} ${params.days_of_infection}
	"""

}

process CONCAT_GISAID_LONG_INFECTIONS {

	input:
	path file_list

	output:

	script:
	"""
	concat_long_infections.R ${params.days_of_infection}
	"""
}

process SEARCH_GISAID_METADATA {

	publishDir params.gisaid_results, mode: 'copy'

	cpus 7

	input:
	path metadata
	path lineage_dates

	output:
	path "*.tsv"

	when:
	params.search_gisaid_metadata == true
		
	script:
	"""
	search_gisaid_metadata.R ${metadata} ${lineage_dates} ${params.days_of_infection} ${task.cpus}
	"""

}

// --------------------------------------------------------------- //
