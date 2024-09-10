#!/usr/bin/env nextflow

/*
Test nextflow script for stylo workflow
*/

/*
PROCESSES
*/
process READFILTERING {
	tag "Nanoq on $sample_id"
	publishDir "${params.outdir}/${sample_id}/reads/", mode: 'copy'
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}_nanoq.fastq.gz"), path("nanoq.log")
	
	script:
	"""
	nanoq -i $reads -l ${params.nanoq_length} -o ${sample_id}_nanoq.fastq.gz -vvv -H >> nanoq.log 2>&1
	"""
}

process DOWNSAMPLE {
	tag "Rasusa on $sample_id"
	publishDir(
		path: "${params.outdir}/${sample_id}/reads/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith('.fastq.gz')) "$filename"
					else if (filename.endsWith('.log')) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}_nanoq_rasusa.fastq.gz"), path("rasusa.log")
	
	script:
	"""
	rasusa -g ${params.rasusa_genome_size} -c ${params.rasusa_coverage} -i ${reads} -o ${sample_id}_nanoq_rasusa.fastq.gz >> rasusa.log 2>&1
	"""
}

process ASSEMBLE {
	tag "Flye on $sample_id"
	publishDir (
		path: "${params.outdir}/${sample_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches 'flye') "$filename"
					else null
		}		
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("flye/assembly.fasta"), path(reads), path("flye")
	
	script:
	if (params.flye_minoverlap == '' )
		"""
		flye ${params.flye_read_type} $reads -g ${params.flye_genome_size} -o flye --threads ${params.flye_threads}
		"""
	else
		"""
		flye ${params.flye_read_type} $reads -g ${params.flye_genome_size} -o flye --threads ${params.flye_threads} --min-overlap ${params.flye_minoverlap}
		"""
}

process ROTATE {
	tag "Circlator on $sample_id"
	publishDir(
		path: "${params.outdir}/${sample_id}/flye",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith '.fasta') "$filename"
					else null
		}		
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads)
	
	output:
	tuple val(sample_id), path("${params.circ_prefix}.fasta"), path(reads)
	
	script:
	"""
	circlator fixstart ${assembly} ${params.circ_prefix}
	"""
}

process POLISH {
	tag "Medaka on $sample_id"
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads)
	
	output:
	tuple val(sample_id), path("${params.medaka_outdir}/consensus.fasta"), path(reads), path("${params.medaka_outdir}")
	
	script:
	"""
	medaka_consensus -i ${reads} -d ${assembly} -o ${params.medaka_outdir} -m ${params.medaka_model}
	"""
}

process RENAME {
	tag "Rename final assembly for $sample_id"
	publishDir(
		path: "${params.outdir}/${sample_id}/medaka/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith '.fasta') "${sample_id}.consensus.fasta"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads)
	
	output:
	tuple val(sample_id), path(assembly), path(reads), stdout
	
	script:
	"""
	basename ${assembly}
	"""
}

process FORMATREADS {
	tag "Reformatting reads for $sample_id"
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads)
	
	output:
	tuple val(sample_id), path(assembly), path(reads), path("reads.fasta")
	
	script:
	"""
	seqtk seq -a ${reads} > reads.fasta
	"""
}

process SOCRU {
	tag "Socru on $sample_id"
	publishDir (
		path: "${params.outdir}/${sample_id}/socru/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith '.txt') "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	containerOptions {
		'-B $PWD:$PWD'
		'--cleanenv'
		'--no-home'
	}	
	
	input:
	tuple val(sample_id), path(assembly), path(reads), val(genus), val(species)
	
	output:
	tuple val(sample_id), path(assembly), path(reads), val(genus), val(species), path('operon_directions.txt'), path("${params.socru_output}.txt"), path("${params.socru_blastoutput}.txt"), optional:true
	
	script:
	if (genus == 'Salmonella')
		"""
		socru Salmonella_enterica ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Escherichia')
		"""
		socru Escherichia_coli ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Shigella')
		"""
		socru Escherichia_coli ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Campylobacter')
		"""
		socru Campylobacter_sp. ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Campylobacter' && species == 'jejuni')
		"""
		socru Campylobacter_jejuni ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Campylobacter' && species == 'coli')
		"""
		socru Campylobacter_coli ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Campylobacter' && species == 'upsaliensis')
		"""
		socru Campylobacter_sp. ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Vibrio')
		"""
		socru Vibrio_sp. ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Vibrio' && species == 'cholerae')
		"""
		socru Vibrio_sp. ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""
	else if (genus == 'Vibrio' && species == 'parahaemolyticus')
		"""
		socru Vibrio_parahaemolyticus ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""	
	else if (genus == 'Yersiniae')
		"""
		socru Yersinia_enterocolitica ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""	
}

process ASSEMBLYQC {
	tag "Busco on $sample_id"
	publishDir (
		path: "${params.outdir}/${sample_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches "${params.busco_output}") "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads), val(genus), val(species)
	
	output:
	tuple val(sample_id), path(assembly), path(reads), val(genus), val(species), path("${params.busco_output}")
	
	script:
	"""
	busco -m ${params.busco_mode} -i ${assembly} -o ${params.busco_output} --auto-lineage-prok
	"""
}

/*
WORKFLOW
*/
workflow {
	Channel //set channel for cat'd fastq files
		.fromPath(params.reads, checkIfExists: true)
		.map{ file -> tuple(file.SimpleName, file) }
		//.view()
		.set{rawreads_ch}
	Channel //set channel for sampleinfo to get genus&species information
		.fromPath(params.sampleinfo, checkIfExists: true)
		.splitCsv(sep: '\t', header:true, strip:true)
		.map{ row -> tuple(row.WGSID, row.GENUS, row.SPECIES) }
		//.view()
		.set{import_ch}

	nanoq_ch = READFILTERING(rawreads_ch)
		.map{ outTuple -> outTuple[0,1] }
		//.view()
	rasusa_ch = DOWNSAMPLE(nanoq_ch)
		.map{ outTuple -> outTuple[0,1] }
		//.view()
	assembly_ch = ASSEMBLE(rasusa_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	circlator_ch = ROTATE(assembly_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	medaka_ch = POLISH(circlator_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	rename_ch = RENAME(medaka_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		.join(import_ch, by: 0)
		//.view()
	busco_ch = ASSEMBLYQC(rename_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4] }
		//.view()
	socru_ch = SOCRU(busco_ch)
}
