#!/bin/bash


# USAGE: bash install_and_pe_testRun.sh
# This script will set up the singularity and nextflow if don't already exist, and will run the workflow on the test dataset


# find the path of install script directory inside eDNAFlow
PWD=$(dirname $(realpath $0))


# install nextflow if doesn't exist
if ! [ -x "$(command -v nextflow)" ]; then
        echo "Installing Nextflow..."
        $PWD/install-nextflow.sh
else
        echo "Nextflow is already installed"
fi

# install singularity if doesn't exist
if ! [ -x "$(command -v singularity)" ]; then
        echo "Installing Singularity..."
        $PWD/install-singularity.sh
else
        echo "Singualrity is already installed"
fi

# find the path of main eDNAFlow script directory
EDNA=$(dirname $PWD)
DATA=$EDNA/testData2_Play
cd $EDNA

# define binding path and container directory in config file for singularity based on user path
cp -n nextflow.config nextflow.config.orig
head -n 157 ./nextflow.config.orig > .nextflow.config.tmp
echo "runOptions = '-B $EDNA'" >> .nextflow.config.tmp
echo "cacheDir = \"$EDNA/cache\"" >> .nextflow.config.tmp
tail -n +158 ./nextflow.config.orig >> .nextflow.config.tmp
mv .nextflow.config.tmp nextflow.config


# run eDNAFlow on local machine for paired-end test data
cd $DATA
nextflow -c $EDNA/nextflow.config run $EDNA/eDNAFlow.nf -profile local --barcode "$DATA/pe_bc*"  --blast_db "$DATA/fake_db/fakeDB.fasta" --minAlignLeng '15' --minsize '3' --qcov '98' --lulu "$EDNA/lulu.R"
