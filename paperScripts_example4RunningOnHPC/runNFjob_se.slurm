#!/bin/bash -l

# This script was used to run eDNAFlow.nf on Pawsey Supercomuter Zeus on partition longq ; Pawsey uses slurm workload manager
# General conditions are set by what follows #SBATCH lines  
# The script was run on a single-end dataset

#SBATCH --account=pawsey0159
#SBATCH --partition=longq
#SBATCH --job-name=nxf-se
#SBATCH --time=4-00:00:00

unset SBATCH_EXPORT

module load singularity
module load nextflow

# for description of flags check the README file 
nextflow run eDNAFlow.nf -profile zeus  --reads 'cook_georg_goc1_goc2_mind_plc_sach.fastq' --barcode '*.txt' --minsize '4' --minLen '50' --perc_identity '90' --maxTarSeq '10' --blast_db '/group/data/blast_v5/nt' 
