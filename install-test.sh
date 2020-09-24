#!/bin/bash

# USAGE: bash install-test.sh 
# This script will set up the singularity and nextflow if don't already exist, and will run the workflow on the test dataset  

#  Clone the eDNAFlow repository 
git clone https://github.com/mahsa-mousavi/eDNAFlow.git

# install nextflow if doesn't exist
if ! [ -x "$(command -v nextflow)" ]; then
	echo "Installing Nextflow..."
	./eDNAFlow/install/install-nextflow.sh
else
	echo "Nextflow is already installed"
fi

# install singularity if doesn't exist 
if ! [ -x "$(command -v singularity)" ]; then
	echo "Installing Singularity..."
	./eDNAFlow/install/install-singularity.sh
else
	echo "Singualrity is already installed"
fi

# set up the path for test data analysis
EDNA=$(pwd)/eDNAFlow
DATA=$EDNA/testData2_Play
cd $EDNA

# define binding path and container directory in config file for singularity based on user path 
cp -n nextflow.config nextflow.config.orig
head -n 157 ./nextflow.config.orig > .nextflow.config.tmp
echo "runOptions = '-B $EDNA'" >> .nextflow.config.tmp
echo "cacheDir = \"$EDNA/cache\"" >> .nextflow.config.tmp
tail -n +158 ./nextflow.config.orig >> .nextflow.config.tmp
mv .nextflow.config.tmp nextflow.config

# run eDNAFlow on local machine for single-end test data
nextflow run ./eDNAFlow.nf -profile local --reads "$DATA/test_30000reads.fastq" --barcode "$DATA/se_bc*" --blast_db "$DATA/fake_db/fakeDB.fasta"

