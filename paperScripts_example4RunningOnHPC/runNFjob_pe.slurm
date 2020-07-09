#!/bin/bash -l

# This script was used to run eDNAFlow.nf on Pawsey Supercomuter Zeus on partition longq ; Pawsey uses slurm workload manager
# General conditions are set by what follows #SBATCH lines  
# The script was run on a paire-end dataset

#SBATCH --account=pawsey0159
#SBATCH --partition=longq
#SBATCH --job-name=nxf-pip01
#SBATCH --time=4-00:00:00
#SBATCH --no-requeue

unset SBATCH_EXPORT

module load singularity
module load nextflow

# for description of flags check the README file
nextflow run eDNAFlow.nf -profile cluster --barcode '*.tab' --minAlignLeng '12' --minsize '2' --minLen '100' --perc_identity '95' --maxTarSeq '10' --blast_db '/group/data/blast/nt' 
