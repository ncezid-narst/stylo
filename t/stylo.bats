#!/usr/bin/env bats

#https://github.com/bats-core/bats-core

thisDir=$BATS_TEST_DIRNAME
nfDir=$(realpath "$thisDir/../schtappe")
configDir=$(realpath "$thisDir/../config")

function note(){
  echo "# $@" >&3 2>&3
  stdbuf -i0 -o0 -e0 echo -ne "" >&3
}

#note "DEBUG: resetting bats tmp dir"; export BATS_SUITE_TMPDIR="./tmp"

@test "Environment" {
	note "Checking nextflow and singularity installation"
	note "============="
	note "which nextflow"
	which nextflow
	note "which singularity"
	which singularity
}

@test "Stylo - Full workflow (12min)" {
	note "Downloading test dataset and assembling"
	note "============="
	note "Check and create stylo.test directory"
	if [ -e ${thisDir}/stylo.test ]; then
		rm -rf ${thisDir}/stylo.test
	fi
	if [ -e ${thisDir}/.nextflow ]; then
		rm -rf ${thisDir}/.nextflow
	fi
	if [ -e ${thisDir}/.nextflow.log ]; then
		rm -f ${thisDir}/.nextflow.log
	fi
	note "Run singularity container of sratoolkit"
	mkdir -pv ${thisDir}/stylo.test/
	cd ${thisDir}/stylo.test/
	singularity exec -B $PWD:/data --no-home --cleanenv docker://pegi3s/sratoolkit:3.0.1 fasterq-dump --outdir fastq_pass --outfile C347.fastq.gz --progress --concatenate-reads SRR22859768
	note "Create samplesheet, add headers, add metadata"
	touch ${thisDir}/stylo.test/sampleinfo.txt
	printf 'BARCODE\tWGSID\tGENUS\tSPECIES' > ${thisDir}/stylo.test/sampleinfo.txt
	printf '\nbarcode01\tC347\tStaphylococcus\taureus' >> ${thisDir}/stylo.test/sampleinfo.txt
	note "Run stylo"
	nextflow run ${nfDir}/stylo.nf -c ${configDir}/stylo.config -profile standard --rasusa_genome_size 2.9MB --rasusa_coverage 20 --flye_genome_size 2.9m --medaka_model r941_min_sup_g507
}

@test "Stylo - Output" {
	note "stylo"
	[ -d ${thisDir}/stylo.test/stylo/ ]
}

@test "Stylo - Process outputs" {
	note "============="
	note "reads"
	[ -d ${thisDir}/stylo.test/stylo/C347/reads ]
	note "flye"
	[ -d ${thisDir}/stylo.test/stylo/C347/flye ]
	note "medaka"
	[ -d ${thisDir}/stylo.test/stylo/C347/medaka ]
	note "staramr_reads"
	[ -d ${thisDir}/stylo.test/stylo/C347/staramr_reads ]
	note "staramr_assembly"
	[ -d ${thisDir}/stylo.test/stylo/C347/staramr_assembly ]
	note "socru"
	[ -d ${thisDir}/stylo.test/stylo/C347/socru ]
	note "busco"
	[ -d ${thisDir}/stylo.test/stylo/C347/busco ]
}

@test "READFILTERING/DOWNSAMPLE - Output" {
	note "C347_nanoq.fastq.gz"
	[ -f ${thisDir}/stylo.test/stylo/C347/reads/C347_nanoq.fastq.gz ]
	note "C347_nanoq_rasusa.fastq.gz"
	[ -f ${thisDir}/stylo.test/stylo/C347/reads/C347_nanoq_rasusa.fastq.gz ]
	note "nanoq.log"
	[ -f ${thisDir}/stylo.test/stylo/C347/reads/nanoq.log ]
	note "rasusa.log"
	[ -f ${thisDir}/stylo.test/stylo/C347/reads/rasusa.log ]
}

@test "ASSEMBLE/ROTATE - Output" {
	note "Flye subdirectories"
	[ -d ${thisDir}/stylo.test/stylo/C347/flye/00-assembly ]
	[ -d ${thisDir}/stylo.test/stylo/C347/flye/10-consensus ]
	[ -d ${thisDir}/stylo.test/stylo/C347/flye/20-repeat ]
	[ -d ${thisDir}/stylo.test/stylo/C347/flye/30-contigger ]
	[ -d ${thisDir}/stylo.test/stylo/C347/flye/40-polishing ]
	note "assembly.fasta"
	[ -f ${thisDir}/stylo.test/stylo/C347/flye/assembly.fasta ]
	note "assembly_info.txt"
	[ -f ${thisDir}/stylo.test/stylo/C347/flye/assembly_info.txt ]
	note "flye.circ.fasta"
	[ -f ${thisDir}/stylo.test/stylo/C347/flye/flye.circ.fasta ]
}

@test "POLISH - Output" {
	note "C347.consensus.fasta"
	[ -f ${thisDir}/stylo.test/stylo/C347/medaka/C347.consensus.fasta ]
}

@test "PLASMIDCHECK - Output" {
	note "reads - plasmidfinder.tsv"
	[ -f ${thisDir}/stylo.test/stylo/C347/staramr_reads/plasmidfinder.tsv ]
	note "assembly - plasmidfinder.tsv"
	[ -f ${thisDir}/stylo.test/stylo/C347/staramr_assembly/plasmidfinder.tsv ]
}

@test "SOCRU - Output" {
	note "socru_output.txt"
	[ -f ${thisDir}/stylo.test/stylo/C347/socru/socru_output.txt ]
	note "GREEN"
	result=$(cat ${thisDir}/stylo.test/stylo/C347/socru/socru_output.txt | awk '{print $2 }')
	[ $result == "GREEN" ]
}

@test "ASSEMBLYQC - Output" {
	note "Busco subdirectories"
	[ -d ${thisDir}/stylo.test/stylo/C347/busco/auto_lineage ]
	[ -d ${thisDir}/stylo.test/stylo/C347/busco/logs ]
	[ -d ${thisDir}/stylo.test/stylo/C347/busco/prodigal_output ]
	[ -d ${thisDir}/stylo.test/stylo/C347/busco/run_bacillales_odb10 ]
	[ -d ${thisDir}/stylo.test/stylo/C347/busco/run_bacteria_odb10 ]
	note "Busco summary files"
	[ -f ${thisDir}/stylo.test/stylo/C347/busco/short_summary.generic.bacteria_odb10.busco.json ]
	[ -f ${thisDir}/stylo.test/stylo/C347/busco/short_summary.generic.bacteria_odb10.busco.txt ]
	[ -f ${thisDir}/stylo.test/stylo/C347/busco/short_summary.specific.bacillales_odb10.busco.json ]
	[ -f ${thisDir}/stylo.test/stylo/C347/busco/short_summary.specific.bacillales_odb10.busco.txt ]
}
