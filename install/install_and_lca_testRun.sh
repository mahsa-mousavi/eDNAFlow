#!/bin/bash

# USAGE: bash install_and_se_testRun.sh
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
mv conf/initialTest.config conf/initialTest.config.bak 2>/dev/null
echo "singularity {" >> conf/initialTest.config
echo "runOptions = '-B $EDNA'" >> conf/initialTest.config
echo "cacheDir = \"$EDNA/cache\"" >> conf/initialTest.config
echo "}" >> conf/initialTest.config

# run eDNAFlow on local machine for single-end test data
cd $DATA
nextflow -c $EDNA/nextflow.config run $EDNA/eDNAFlow.nf --test --taxonomyAssignment --zotuTable "$DATA/curated_zotuTable_2test_LCAscript.tab" --blastFile "$DATA/blastResult_2test_LCAscript.tab" --lca_qcov '100' --lca_pid '98' --lca_diff '1' --lca_output "$EDNA/lca_result"

