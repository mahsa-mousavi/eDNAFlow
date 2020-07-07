# eDNAFlow
## About the workflow
eDNAFlow, is a fully automated pipeline that employs a number of state-of-the-art applications to process eDNA data from raw sequences (single-end or paired-end) to generation of curated and non-curated zero-radius operational taxonomic units (ZOTUs) and their abundance tables. This pipeline is based on Nextflow and Singularity which enable a scalable, portable and reproducible workflow using software containers on a local computer, clouds and high-performance computing (HPC) clusters. We also present an in-house Python script to assign taxonomy to ZOTUs based on user specified thresholds for assigning Lowest Common Ancestor (LCA).

## Setup the pipeline
To run the pipeline, first Nextflow and Singularity have to be installed or made available for loading as modules (e.g. in the case of running it on an HPC cluster) on your system. This pipeline was built and tested with versions 19.10 and 3.5.2 of Nextflow and Singularity, respectively. For easy installation of Nextflow, follow the instructions at [nextflow installation](https://www.nextflow.io/docs/latest/getstarted.html). To install Singularity version 3.5.2, follow the instructions at [singularity installation](https://sylabs.io/guides/3.5/admin-guide/installation.html).

Download eDNAFlow scripts (including eDNAFlow.nf, eDNAFlow.config and lulu.R) to the directory where your unzipped sequencing and Multiplex identifier (MID) tag (here defined as “barcode”) files exist. Alternatively, for easy accessibility and testing, the user can also clone the Git repository so that all the scripts and test data are downloaded and in one folder. To clone the repository to your directory, run this command: 

`git clone https://github.com/mahsa-mousavi/eDNAFlow.git` 

## Download database
One of the mandatory parameters to run eDNAFlow is to provide a path to a local GenBank nucleotide (nt) and/or your custom database. To download the NCBI nucleotide database locally, follow the steps below.

1) Download the official [BLAST+ container](https://github.com/ncbi/blast_plus_docs#show-blast-databases-available-for-download-from-ncbi) with Singularity using the below command (tested on Ubuntu 18.04):

`singularity pull --dir  directoryName docker://ncbi/blast:2.10.0`

*directoryName* is the path to the directory where you want to keep the container image

2) Make a folder where you want to keep the database and from there run the following command: 

`singularity run directoryName/blast_2.10.0.sif update_blastdb.pl --decompress nt`

\* Please be aware step 2 will take a long time and will need a large space available on the disk due to the size of GenBank nucleotide database. We are unable to provide an estimate of these values.  


## Test run Basic usage
The pipeline scripts, along with two small test datasets (one single-end and one paired-end) and a small test database are available to download from Github page.
We suggest running the pipeline first on the test datasets provided, to ensure it is setup properly. 

**Basic usage:**
To run the pipeline on your local machine on the single-end test dataset use the following command:

`nextflow run eDNAFlow.nf -profile local --reads 'test_30000reads.fastq' --barcode 'se_bc*' --blast_db 'Path2TestBlastDataset/fakeDB.fasta'`

## Description of run parametrs
The following parameters can be adjusted on the command line:

`-profile`: depending on the execution platform you would like to run the pipeline on this parameter must be set to `-profile cluster` for running on HPC; `-profile cloud` for running on cloud, and `-profile local` for running on local machine   

`--reads 'read.fastq'`: provide the name of your raw fastq file; no need to specify this option if reads are paired-end as they get identified automatically by default; reads must be unzipped

`--barcode 'bc.tab'`: your barcode file name; barcode file format must match [OBITools requirement](https://pythonhosted.org/OBITools/scripts/ngsfilter.html); if multiple barcode files exist (e.g. bc_1.txt; bc_2.txt) it can be specified like this: bc*.txt

`--blast_db 'path2/LocalGenbankDatabase'`: 

`custom_db ''`: 

`--minQuality '20'`: the minimum Phred quality score to apply for quality control of raw sequences; Default is 20; must be an integer 

`--minAlignLeng '12'`: the minimum alignment length for merging read1 and read2 in case of paired-end sequences; Default is 12; must be an integer

`--minLen '50'`: the minimum length allowed for sequences; Default is 50; must be an integer

`--minsize '8'`: the minimum abundance; input sequences with lower abundances are removed; Default is 8; to check how adjusting this option affects the results check out [Usearch documentation](https://drive5.com/usearch/manual/cmd_unoise3.html)

`--maxTarSeq '10'`: a blast parameter; the maximum number of target sequences for hits per query to be returned by Blast; Default is 10

`--perc_identity '95'`: a blast parameter; percentage of identical matches; Default is 95

`--evalue '1e-3'`: a blast parameter; expected value for saving blast hits; Default is 1e-3

`--qcov '100'`: a blast parameter; the percent of the query that has to form an alignment against the reference to be retained; Higher values prevent spurious alignments of only a short portion of the query to a reference; Default is 100

`--lulu 'lulu.R'`: an R script to run post-clustering curation with default setting of LULU; this file has been provided and must be present in the same directory as other scripts; by default eDNAFlow will be looking for this file

`--mode 'usearch32'`: by default eDNAFlow uses free version of usearch (i.e. usearch 32 bit version); if you have access to 64bit version it can be set via changing mode as `--mode 'usearch64'`

`--usearch64 'Path2/Usearch64/executable'`: if mode is set to usearch64, then this option has to be specified; the full path must point to usearch64 executable 
