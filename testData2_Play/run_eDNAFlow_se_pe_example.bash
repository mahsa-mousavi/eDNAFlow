#!/usr/bin/bash

# Example of running the eDNAFlow script for single-end test data using default parameters 
nextflow run eDNAFlow.nf -profile local --barcode 'se_bc*' --reads 'test_30000reads.fastq' --blast_db '$HOME/Desktop/test-nextFlow-local/fake_db/fakeDB.fasta' --custom_db '$HOME/Desktop/test-nextFlow-local/custom_db/cus_fakeDB.fasta'

# Example of running the eDNAFlow script for paired-end test data using some user specified parameters 
nextflow run eDNAFlow.nf -profile local --barcode 'pe_bc*'  --blast_db '$HOME/Desktop/test-nextFlow-local/fake_db/fakeDB.fasta' --minAlignLeng '15' --minsize '3' --qcov '98' 

