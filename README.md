# eDNAFlow
## About the workflow
eDNAFlow, is a fully automated pipeline that employs a number of state-of-the-art applications to process eDNA data from raw sequences (single-end or paired-end) to generation of curated and non-curated zero-radius operational taxonomic units (ZOTUs) and their abundance tables. This pipeline is based on Nextflow and Singularity which enable a scalable, portable and reproducible workflow using software containers on a local computer, clouds and high-performance computing (HPC) clusters. We also present an in-house Python script to assign taxonomy to ZOTUs based on user specified thresholds for assigning Lowest Common Ancestor (LCA).

## Setup the pipeline
To run the pipeline, first Nextflow and Singularity have to be installed or made available for loading as modules (e.g. in the case of running it on an HPC cluster) on your system. This pipeline was built and tested with versions 19.10 and 3.5.2 of Nextflow and Singularity, respectively. For easy installation of Nextflow, follow the instructions at [nextflow installation](https://www.nextflow.io/docs/latest/getstarted.html). To install Singularity version 3.5.2, follow the instructions at [singularity installation](https://sylabs.io/guides/3.5/admin-guide/installation.html).

Download eDNAFlow scripts (including eDNAFlow.nf, eDNAFlow.config and lulu.R) to the directory where your unzipped sequencing and Multiplex identifier (MID) tag (here defined as “barcode”) files exist. Alternatively, for easy accessibility and testing, the user can also clone the Git repository so that all the scripts and test data are downloaded and in one folder. To clone the repository to your directory, run this command: `git clone https://github.com/mahsa-mousavi/eDNAFlow.git` 

## Download database
One of the mandatory parameters to run eDNAFlow is to provide a path to a local GenBank nucleotide (nt) and/or your custom database. To download the NCBI nucleotide database locally, follow the steps below.

1) Download the official [BLAST+ container](https://github.com/ncbi/blast_plus_docs#show-blast-databases-available-for-download-from-ncbi) with Singularity using the below command (tested on Ubuntu 18.04):

`singularity pull --dir  directoryName docker://ncbi/blast:2.10.0`

*directoryName* is the path to the directory where you want to keep the container image

2) Make a folder where you want to keep the database and from there run the following command: 

`singularity run directoryName/blast_2.10.0.sif update_blastdb.pl --decompress nt`

\* Please be aware step 2 will take a long time and will need a large space available on the disk due to the size of GenBank nucleotide database. We are unable to provide an estimate of these values.  
