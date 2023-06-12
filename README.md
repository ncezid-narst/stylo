# Schtappe

Welcome to Schtappe, what is going to be a suite of Nextflow workflows created to address a variety of CDC/NCEZID/DFWED/EDLB/NARST Nanopore sequencing bioinformatics needs. Although Schtappe is intended for CDC users, the portable nature of Nextflow should allow for external users to easily adapt these workflows to their specific environments and/or needs. To facilitate that, all of the processes use Singularity containers.

As of 4/17/2023, the workflows currently available include:

 * Stylo - Nanopore assembly workflow from basecalled reads to polished assembly plus assembly QC, metrics, and plasmid replicon detection

See [Future Plans](#future-plans) for information on updates to existing and in-development workflows.

## Install
Navigate to your home directory and git clone the repository.
```bash
$ git clone https://github.com/ncezid-narst/Schtappe.git
```
You will need to install Nextflow if you don't already have it: https://www.nextflow.io/docs/latest/getstarted.html

You will also need to install Singularity if you don't already have it: https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html

If you're working on CDC servers, run `module load nextflow/XX.XX.X` to load Nextflow, Singularity, and Java modules.

As of 4/17/2023, this workflow was developed on Nextflow version 22.10.6.

## Workflows

### Stylo

#### Overview:
* READFILTERING - Filters reads based on read-length using Nanoq
* DOWNSAMPLE - Randomly downsamples read set based on organism genome size and desired coverage using Rasusa
* ASSEMBLE - Generates long-read assembly using Flye
* HYBRID - Generates hybrid assembly using Unicycler (alternative option to using Flye)
* ROTATE - Changes start position of contigs using Circlator
* POLISH - Creates consensus calls using Medaka
* RENAME - Renames polished assembly for convenience
* FORMATREADS - Preps reads to be run through staramr using Seqtk 
* PLASMIDCHECK - Detects plasmid replicons in reads and assembly using Staramr
* SOCRU - Does Socru (*shrug) 
* ASSEMBLYQC - Generates assembly quality metrics using BUSCO

![image](https://user-images.githubusercontent.com/112518496/232568302-f03f21ca-0918-402c-a3e5-2bd499e1a135.png)

#### Parameters:
Parameters for each process can be changed in `stylo.config` under the first bracketed section `params`. Check out [Resources](#resources) for links to each process's main github page to learn more about process-specific parameters.

Prior to running stylo, make sure the INITIAL PARAMETERS are set accurately - the default settings are as follows:
```java
	//Initial parameters
	reads = 'fastq_pass/**.fastq.gz'
	sampleinfo = 'sampleinfo.txt'
	outdir = 'stylo'
	unicycler = false
```

`reads`: Stylo will look for any and all fastq.gz files under a directory and assume _each_ one is a unique sample. Prior to running stylo, you should concatenate, rename, and compress (if they're not already) your reads. Anything prior to the file extension will be used as the Sample ID. For example:
```bash
fastq_pass/
├── 01_2014C-3598
│   └── 01_2014C-3598_all.fastq.gz
├── 02_2014C-3599
│   └── 02_2014C-3599_all.fastq.gz
├── 03_2014C-3857
│   └── 03_2014C-3857_all.fastq.gz
```

`sampleinfo`: Tab-delimited text file with sample information. For example:
```bash
BARCODE WGSID   GENUS   SPECIES
barcode01       01_2014C-3598_all    Salmonella      enterica
barcode02       02_2014C-3599_all    Salmonella      enterica
barcode03       03_2014C-3857_all    Salmonella      enterica
```
* BARCODE: Standard barcode output from MinKNOW e.g. barcode01-barcode96
* WGSID: Sample ID. **Must match filename of concatenated reads**.
* GENUS: Sample genus, used in SOCRU
* SPECIES: Sample species, used in SOCRU

`outdir`: Name of Stylo output directory. Default name is set to `stylo`.
`unicycler`: Option to use Unicycler instead of Flye as the assembler. Default is set to `false`.

You can see how parameters are used in the next section **Usage**.

NOTE: Support for hybrid assemblies using short-reads hasn't been added yet. This option was added as an experiment to test how well k-mer based assemblers perform with ONT's v14/r10.4.1 chemistry.

#### Processes:
Directives for each process can be changed in `stylo.config` under the second bracketed section `process`. This is where you can update the containers used in each process. Check out [Resources](#resources) to see a full list of all the containers and the tools' githubs.

#### Profiles:
Configuration settings for each profile can be changed in `stylo.config` under the third bracketed section `profiles`. This is where you can update or create profiles that will dictate where and how each process is run. By default, there are two main profiles and three auxiliary profiles:

* `standard`: Will execute stylo using the 'local' executor, running processes on the computer where Nextflow is launched. 
* `sge`: Will execute stylo using a Sun Grid Engine cluster, running processes on the HPC (qsub).
* `short`: Auxiliary profile to change the sge default queue to short queue
* `gpu`: Auxiliary profile to chage the sge default queue to gpu queue
* `highmem`: Auxiliary profile to change the sge default queue to highmem queue

You can see how profiles are used in the next section **Usage**.

NOTE: The default profile settings were mostly pulled from recommendations made by CDC Scicomp in their Nextflow training called 'Reproducible, scalable, and shareable analysis workflows with Nextflow'. There is a good chance you will have to create/modify your own profile to run stylo using your institution's computing environment. Check out [Resources](#resources) to learn more about creating profiles.

#### Usage
Once you've made the necessary changes to the configuration file to run the workflow on your computing environment and have set up inital parameters, you can run stylo just as you would any nextflow workflow:
```bash
nextflow run /path/to/nanoporeWorkflow/schtappe/stylo.nf -c /path/to/nanoporeWorkflow/config/stylo.config
```
Nextflow is picky about single-hyphen flags vs. double-hyphen flags. Single-hyphens affect the nextflow command while double-hyphens affect the parameters in the configuration file. For example, to change the initial parameters without directly editing `stylo.config`:
```bash
nextflow run /path/to/nanoporeWorkflow/schtappe/stylo.nf -c /path/to/nanoporeWorkflow/config/stylo.config \
  --reads path/to/your/reads/**.fastq.gz \
  --sampleinfo yoursampleinfofile.txt \
  --outdir youroutputdirectory \
  --unicycler true
```

By default, nextflow will run locally. If you want to specify a profile, use the `-profile` flag. For example, to qsub stylo's processes:
```bash
nextflow run /.../nanoporeWorkflow/schtappe/stylo.nf -c /.../nanoporeWorkflow/config/stylo.config -profile sge
```

You can change the queue by adding the auxiliary profile name, separated by a comma:
```bash
nextflow run /.../nanoporeWorkflow/schtappe/stylo.nf -c /.../nanoporeWorkflow/config/stylo.config -profile sge,highmem
```
Run `nextflow help` or `nextflow run -help` for more information on nextflow flags.

NOTE: Nextflow applies the same parameters to each sample being processed. This means you'll want to run stylo on read sets all of the same organism or at least the same genome size and all have been generated using the same chemistry and guppy basecaller version (affects flye_read_type and medaka_model) This could change in the future by adding more fields to the sampleinfo sheet, but for now it is what it is.

#### Output
Here's what stylo output looks like per sample(directories only):
```bash
stylo/
└── PNUSAS002131
    ├── busco
    │   ├── auto_lineage
    │   │   ├── run_archaea_odb10
    │   │   │   ├── busco_sequences
    │   │   │   │   ├── fragmented_busco_sequences
    │   │   │   │   ├── multi_copy_busco_sequences
    │   │   │   │   └── single_copy_busco_sequences
    │   │   │   └── hmmer_output
    │   │   └── run_bacteria_odb10
    │   │       ├── busco_sequences
    │   │       │   ├── fragmented_busco_sequences
    │   │       │   ├── multi_copy_busco_sequences
    │   │       │   └── single_copy_busco_sequences
    │   │       ├── hmmer_output
    │   │       └── placement_files
    │   ├── logs
    │   ├── prodigal_output
    │   │   └── predicted_genes
    │   │       └── tmp
    │   ├── run_bacteria_odb10
    │   │   ├── busco_sequences
    │   │   │   ├── fragmented_busco_sequences
    │   │   │   ├── multi_copy_busco_sequences
    │   │   │   └── single_copy_busco_sequences
    │   │   ├── hmmer_output
    │   │   └── placement_files
    │   └── run_enterobacterales_odb10
    │       ├── busco_sequences
    │       │   ├── fragmented_busco_sequences
    │       │   ├── multi_copy_busco_sequences
    │       │   └── single_copy_busco_sequences
    │       └── hmmer_output
    ├── flye
    │   ├── 00-assembly
    │   ├── 10-consensus
    │   ├── 20-repeat
    │   ├── 30-contigger
    │   └── 40-polishing
    ├── medaka
    ├── reads
    ├── socru
    ├── staramr_assembly
    │   └── hits
    └── staramr_reads
        └── hits
```

## Future plans
Small bug/typo/formatting issues aside, updates to stylo will likely include: 
* Support for short-reads to run through Unicycler
* More param options for each process (ideally I guess all of them?)
* A barcode renaming script that will handle the read concatenating necessary to run the workflow
* Test set of sample data

Workflows in development:
* Messer: Plasmid-recovery workflow to individually fix and assemble plasmid contigs
* Geteilt: Assembly-based antibiotic-resistance screening workflow

## Resources
### Containers:
* Nanoq:
  * https://hub.docker.com/r/jimmyliu1326/nanoq
  * https://github.com/esteinig/nanoq
* Rasusa: 
  * https://hub.docker.com/r/staphb/rasusa
  * https://github.com/mbhall88/rasusa
* Flye: 
  * https://hub.docker.com/r/staphb/flye
  * https://github.com/fenderglass/Flye
* Unicycler: 
  * https://hub.docker.com/r/staphb/unicycler
  * https://github.com/rrwick/Unicycler
* Circlator: 
  * https://hub.docker.com/r/staphb/circlator
  * https://github.com/sanger-pathogens/circlator
* Medaka: 
  * https://hub.docker.com/r/ontresearch/medaka
  * https://github.com/nanoporetech/medaka
* Seqtk: 
  * https://hub.docker.com/r/staphb/seqtk
  * https://github.com/lh3/seqtk
* Staramr:
  * https://hub.docker.com/r/staphb/staramr
  * https://github.com/phac-nml/staramr
* Socru:
  * https://hub.docker.com/r/quadraminstitute/socru
  * https://github.com/quadram-institute-bioscience/socru 
* BUSCO:
  * https://hub.docker.com/r/ezlabgva/busco
  * https://busco.ezlab.org/

### Nextflow:
If you're unfamiliar with Nextflow or would like to just learn more, consider doing these free trainings found here: https://training.nextflow.io/
The Nextflow documentation is super helpful as well, especially to learn more about what process directives and profile configurations you can include in your local copy of `stylo.config`. An entire section is dedicated to just containers which should help troubleshoot any issues with Singularity or assist in using Docker instead. Nextflow documentation can be found here: https://www.nextflow.io/docs/latest/index.html

If you're a CDC user, you should register for Scicomp's Nextflow tutorials here: https://info.biotech.cdc.gov/info/training-portal-updated-tracks/
Or you can access the material directly here: https://training.biotech.cdc.gov/nextflow/

## Thanks
Special thanks to: 
* Jess Chen, Lee Katz, and Curtis Kapsak for their enthusiastic guidance and wealth of data science expertise
* Kaitlin Tagg and Hattie Webb for their encouragement and passion for plasmids that started this whole venture
* StaPH-B who are an indispensable resource for bioinformatics in public health
* Seqera whose workshop is a fantastic introduction to writing in Nextflow
* OpenAI - ChatGPT is an amazing tool for looking up documentation, especially when you're a bad coder

## Etymology
The names 'Schtappe', 'Stylo', 'Messer', and 'Geteilt' are all spells from a wonderful book series called Ascendance of a Bookworm by Miya Kazuki. Obviously, they're all German words. The meaning of the spells all somewhat capture their respective purposes. 
