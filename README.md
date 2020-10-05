# eDNAFlow
## About the workflow
eDNAFlow, is a fully automated pipeline that employs a number of state-of-the-art applications to process eDNA data from raw sequences (single-end or paired-end) to generation of curated and non-curated zero-radius operational taxonomic units (ZOTUs) and their abundance tables. This pipeline is based on Nextflow and Singularity which enable a scalable, portable and reproducible workflow using software containers on a local computer, clouds and high-performance computing (HPC) clusters. We also present an in-house Python script to assign taxonomy to ZOTUs based on user specified thresholds for assigning Lowest Common Ancestor (LCA).

## Setup and test the pipeline
To run the pipeline, first Nextflow and Singularity have to be installed or made available for loading as modules (e.g. in the case of running it on an HPC cluster) on your system. This pipeline was built and tested with versions 19.10 and 3.5.2 of Nextflow and Singularity, respectively. We strongly suggest that you first try the pipeline on your local machine using the test dataset provided. 

We are providing scripts that will install Nextflow and Singularity on your local machine if they don't already exist and will run the pipeline on the test dataset. The scripts have been successfully tested on Ubuntu 16, 18 and 20. Alternatively, for manual installation of Nextflow, follow the instructions at [nextflow installation](https://www.nextflow.io/docs/latest/getstarted.html). To install Singularity version 3.5.2, follow the instructions at [singularity installation](https://sylabs.io/guides/3.5/admin-guide/installation.html). If working on HPC, you may need to contact your HPC helpdesk.  

**Follow the steps below:**

1- Clone the Git repository so that all the scripts and test data are downloaded and in one folder. To clone the repository to your directory, run this command: 

`git clone https://github.com/mahsa-mousavi/eDNAFlow.git` 

2- Next, in your terminal go to the "install" directory which is located inside the "eDNAFlow" directory (e.g. `cd eDNAFlow/install`) 

3- Once inside the install directory run: `bash install_and_se_testRun.sh` to try the pipeline on single-end test data or run `bash install_and_pe_testRun.sh` for testing it on paired-end test data. As every step gets completed you will see a ✔ next to the relavant step or ✘ if it fails.        

4- If all goes well you can now find all the results inside folder "testData2_Play"  

**What should you expect in result folders:**

In folders 01 to 07 you will find soft link(s) to the final result files of each step as explained below. To check any intermediate files or debug any particular step check the relevant folder inside work directory (check out the files starting with .command* for debugging, checking the log, etc if necessary). 

`01_a_quality_Filtering_YourSequenceFileName`: The filtered fastq file (adapterRemoval package) 

`01_b_assigned_dmux_YourSequenceFileName_yourBarcodeFileName`: Demultiplexed file for each barcode file (OBITools package)

`02_Length_filtered_YourSequenceFileName`: Length filtered demultiplexed file (OBITools package)

`03_splitSamples_YourSequenceFileName`: Splitted files per sample (OBITools package)

`04_relabel_Cat_YourSequenceFileName`: Count of filtered demultiplexed reads in each sample (i.e. CountOfSeq.txt), each sample, concatenated fastq and fasta files

`05_Uniques_ZOTUs_YourSequenceFileName`: Unique file, Final ZOTU fasta file and Final ZOTU table (this table is uncurated) (USEARCH package)

`06_blast_YourSequenceFileName`: Blast result, match file and table for generating curated result (BLAST package)  

`07_lulu_YourSequenceFileName`: LULU results including map file and curated ZOTU table (LULU package)

`work`: Holds all the results, intermediate files, ...  

`.nextflow`: Nextflow generated folder holding history info 

`.nextflow.log`: Nextflow generated log file which can be used for debugging and checking what each number in *work directory* maps to    

## Running eDNAFlow on your data
Make sure eDNAFlow scripts (including eDNAFlow.nf, nextflow.config and lulu.R) are in the same directory where your unzipped sequencing and Multiplex identifier (MID) tag (here defined as “barcode”) files exist. 

## Download database
One of the mandatory parameters to run eDNAFlow is to provide a path to a local GenBank nucleotide (nt) and/or your custom database. To download the NCBI nucleotide database locally, follow the steps below.

1) Download the official [BLAST+ container](https://github.com/ncbi/blast_plus_docs#show-blast-databases-available-for-download-from-ncbi) with Singularity using the below command (tested on Ubuntu 18.04):

`singularity pull --dir  directoryName docker://ncbi/blast:2.10.0`

**directoryName** is the path to the directory where you want to keep the container image

2) Make a folder where you want to keep the database and from there run the following command: 

`singularity run directoryName/blast_2.10.0.sif update_blastdb.pl --decompress nt`

\* Please be aware step 2 will take some time and will need a large space available on the disk due to the size of GenBank nucleotide database. For us it took under 2 hours on the NCBI default 1 core setting (~10MB per second), and was done a lot faster using an HPC data transfer node (hpc-data.pawsey.org.au) or copyq (with 16 cores) on Zeus, at almost 100MB per second.


## Basic command usage

Example of basic command to run the pipeline on your local machine on single-end/paired-end data with multiple barcode files using blast and/or custom database:

`nextflow run eDNAFlow.nf -profile local --reads 'file.fastq' --barcode 'bc_*.txt' --blast_db 'path2/LocalGenbankDatabase/nt' [OPTIONS]` 

`nextflow run eDNAFlow.nf -profile local --barcode 'pe_bc*'  --blast_db 'Path2TestBlastDataset/file.fasta' --custom_db 'path2/customDatabase/myDb' [OPTIONS]` 

## Description of run options
The following parameters can be adjusted on the command line:

**Mandatory parameters**

`-profile`: depending on the execution platform you would like to run the pipeline on this parameter must be set to `-profile cluster` for running on HPC; `-profile cloud` for running on cloud, and `-profile local` for running on local machine   

`--reads 'read.fastq'`: provide the name of your raw fastq file; **You should NOT specify this option if reads are paired-end as they get identified automatically by default**; reads must be unzipped

`--barcode 'bc.tab'`: your barcode file name; barcode file format must match [OBITools requirement](https://pythonhosted.org/OBITools/scripts/ngsfilter.html); if multiple barcode files exist (e.g. bc_1.txt; bc_2.txt) it can be specified like this: bc*.txt

* At least one of the below databases must be specified. 

`--blast_db 'path2/LocalGenbankDatabase/nt'`: the path to where nt databse is stored 

`--custom_db 'path2/customDatabase/myDb'`: the path to where custom database is stored

**Optional parameters**

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

## Adjusting configuration file
All run parameters listed above, as well as settings of execution platforms are prearranged in the nextflow.config file. As described above, all run parameters should be adjusted according to user needs at the command line when running the script and **should not be changed** in the config file itself.

However, if the user wants to run the script on HPC or cloud there are couple of lines that have to be adjusted in configuration file itself:

1) line 63 for HPC or line 137 for cloud: `runOptions = '-B path'` where path is the working directory/volume on your HPC/cloud to add Singularity bind paths to your Nextflow call

2) line 64 for HPC: `cacheDir = 'path'` where path is the directory where remote Singularity images are going to be stored


**NOTE**

eDNAFlow was tested on an HPC with slurm resource manager, but it can be easily adjusted to other executors (e.g. pbs, moab, etc) by setting line 76 (e.g. executor = 'pbs') and adjusting cpus, queue, time, memory and clusterOptions parameters in accordance with the HPC resource available. For the complete list of the executors supported by nextflow please see [nextflow documentation about executors](https://www.nextflow.io/docs/latest/executor.html).  

## LCA (Lowest Common Ancestor) script for assigning taxonomy

To run the LCA script (written with Python version 3) you need to download the two Python scripts (i.e. working_function.py and runAssign_collapsedTaxonomy.py), and have them in the same folder with a copy or a soft-link of your blast result as well as your ZOTU table. 

**NOTE 1:**

This script is not limited to eDNAFlow created files. It can be used for assigning taxonomy of any OTU, ZOTU, and or ASV files, as long as the blastn of the fasta file is performed using the following format:

-outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

***If you have used eDNAFlow to generate your files, then the above has already been taken care of.***

**NOTE 2:**

If you want to use the curated ZOTU table, first you need to make some little changes in that file: 

1) Remove all occurrence of ײ in the curated file. 
2) In the first line add #ID followed by a tab.

**Run the script:**

On the terminal type: python3 runAssign_collapsedTaxonomy.py

Provided you have wget* installed on your system the script starts downloading the taxonomy dump file from NCBI and will put it in a folder with the date of that day. It will then ask users for the input file names (i.e. OTU/ZOTU/ASV table file and blast result file), filtering thresholds and what the output file name should be.

The filtering applied in this script is based on a set of user specified thresholds, including query coverage (qCov), percentage identity (% identity) and the difference (Diff) between % identities of two hits when their qCov is equal. Setting qCov and % identity thresholds ensures that only BLAST hits >= to those thresholds will progress to the Diff comparison step. Setting Diff means that if the absolute value for the difference between % identity of hit1 and hit2 is > Diff, then a species level taxonomy will be returned, otherwise taxonomy of that ZOTU will be dropped to the lowest common ancestor. This script produces two files, a file in which the taxonomy is assigned (the final result), and an intermediate file which will give the user an idea of why some ZOTUs may have been assigned to the lowest common ancestor 

 

*Mac users can easily install wget following instruction [here](https://www.cyberciti.biz/faq/howto-install-wget-om-mac-os-x-mountain-lion-mavericks-snow-leopard/)
