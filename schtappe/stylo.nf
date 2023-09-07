#!/usr/bin/env nextflow

/*
Test nextflow script for stylo workflow
*/

/*
PROCESSES
*/

process READFILTERING {
	/*
	- READFILTERING uses Nanoq to filter basecalled ONT reads based on READ LENGTH, default value in config file is set to 1000bp
	- Inputs (per sample): sample name and path to reads from rawreads channel
	- Output: filtered reads (sample_id_nanoq.fastq.gz) and log of nanoq (redirected stdout)
	*/
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
	/*
	- DOWNSAMPLE uses rasusa to randomly subsample basecalled ONT reads
	- Input: sampleid, nanoq-filtered reads
	- Output: sampleid, rasusa&nanoq filtered reads, rasusa log file (redirected stdout)
	*/
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
	/*
	- ASSEMBLE uses flye to assemble long-reads
	- Input: sampleid, path to filtered&downsampled reads
	- Output: sampleid, flye assembly & output directory
	*/
	
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
	if (params.flye_asm_coverage)
	"""
	flye ${params.flye_read_type} $reads -g ${params.flye_genome_size} -o flye --threads ${params.flye_threads} --asm-coverage ${params.flye_asm_coverage}
	"""
	else
	"""
	flye ${params.flye_read_type} $reads -g ${params.flye_genome_size} -o flye --threads ${params.flye_threads}
	"""
}

process HYBRID {
	/*
	- HYBRID uses Unicycler to assemble using both short reads and long reads in conjunction
	- Input: sampleid, path to paired end reads, path to long reads
	- Output: sampleid, hybrid assembly
	*/
	tag "Unicycler on $sample_id"
	publishDir (
		path: "${params.outdir}/${sample_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches "${params.unicycler_output}") "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${params.unicycler_output}/assembly.fasta"), path(reads), path("${params.unicycler_output}")
	
	script:
	"""
	unicycler --unpaired ${reads} -l ${reads} -o ${params.unicycler_output} -t ${params.unicycler_threads} --min_fasta_length ${params.unicycler_min_fasta_length} --mode ${params.unicycler_mode} --keep ${params.unicycler_keep} --verbosity ${params.unicycler_verbosity}
	"""
}

process ROTATE {
	/*
	ROTATE uses circlator to change the starting point of the chromosomal contig to a default point
	- Input: sampleid, flye assembly
	- Output: sampleid, rotated flye assembly
	*/
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
	/*
	POLISH uses medaka to polish the flye assembly
	- Input: sampleid, flye assembly
	- Output: sampleid, polished assembly, medaka output directory
	*/
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
	/*
	RENAME changes the file name of the polished assembly using the sampleid that's been carried through the entire workflow
	- Input: sampleid, polished assembly
	- Output: renamed assembly
	*/
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
	/*
	FORMATEREADS uses seqtk to convert the reads from fastq to fasta (staramr takes fasta only)
	- Input: long reads
	- Output: sampleid, reformatted assembly
	*/
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

process PLASMIDCHECK {
	/*
	PLASMIDCHECK uses starmar to run staramr on both the reads and assembly a check for known, major plasmid contigs + resistance determinants
	- Input: sampleid, path to reads, path to assembly
	- Output: sampleid, staramr output for assembly, staramr output for reads
	*/
	tag "Staramr on ${sample_id}"
	publishDir (
		path: "${params.outdir}/${sample_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches("staramr_assembly")) "$filename"
					else if (filename.matches("staramr_reads")) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads), path(reads_fasta)
	
	output:
	tuple val(sample_id), path(assembly), path(reads), path("staramr_reads"), path("staramr_assembly")
	
	script:
	"""
	staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
	staramr search --pid-threshold 90 --percent-length-overlap-resfinder 50 --no-exclude-genes --database databases -o "staramr_reads" ${reads_fasta}
	staramr search --pid-threshold 90 --percent-length-overlap-resfinder 50 --no-exclude-genes --database databases -o "staramr_assembly" ${assembly}
	"""
}

process SOCRU {
	/*
	SOCRU uses socru to check the assembly for structural biological improbability, returning a color: green, yellow, amber, red - as inditcators for 'assembly health' 
	- Input: sampleid, assembly, reads, genus, and species
	- OUtput: sampleid, socru output files
	*/
	tag "Socru on $sample_id"
	publishDir (
		path: "${params.outdir}/${sample_id}/socru/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith '.txt') "$filename"
					else if (filename.endsWith '.pdf') "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), path(assembly), path(reads), val(genus), val(species)
	
	output:
	tuple val(sample_id), path(assembly), path(reads), val(genus), val(species), path('operon_directions.txt'), path("${params.socru_output}.txt"), path("${params.socru_blastoutput}.txt"), path('genome_structure.pdf')
	
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
	else if (genus == 'Campylobacter' && species == 'sp')
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
	else if (genus == 'Vibrio' && species == 'sp')
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
	else if (genus == 'Staphylococcus' && species == 'aureus')
		"""
		socru Staphylococcus_aureus ${assembly} --output_file "${params.socru_output}.txt" --top_blast_hits "${params.socru_blastoutput}.txt"
		"""	
}

process ASSEMBLYQC {
	/*
	ASSEMBLYQC uses BUSCO as a final assembly metrics + quality check
	- Input: sampleid, final assembly
	- OUtput: sampleid, busco output directory
	*/
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
	if (params.unicycler) {
		assembly_ch = HYBRID(rasusa_ch)
		.map{ outTuple -> outTuple[0,1,2]}
		//.view()
	}
	else {
		assembly_ch = ASSEMBLE(rasusa_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	}
	circlator_ch = ROTATE(assembly_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	medaka_ch = POLISH(circlator_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	rename_ch = RENAME(medaka_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		//.view()
	seqtk_ch = FORMATREADS(rename_ch)
		.map{ outTuple -> outTuple[0,1,2,3] }
		//.view()
 	staramr_ch = PLASMIDCHECK(seqtk_ch)
		.map{ outTuple -> outTuple[0,1,2] }
		.join(import_ch, by: 0)
		//.view()
	socru_ch = SOCRU(staramr_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4] }
		//.view()
	busco_ch = ASSEMBLYQC(socru_ch)
}
