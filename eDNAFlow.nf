#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""

    Usage:

    Examples of basic commands to run the pipeline on local machine:

    For single-end run:
    nextflow run eDNAFlow.nf --reads 'file.fastq' --barcode 'bc_*.txt' --blast_db 'path2/LocalGenbankDatabase/nt' [OPTIONS] 
    
    For paired-end run:
    nextflow run eDNAFlow.nf --barcode 'pe_bc*' --blast_db 'Path2TestBlastDataset/file.fasta' --custom_db 'path2/customDatabase/myDb' [OPTIONS]
    
    For running LCA taxonomy assignment script
    nextflow run eDNAFlow.nf --taxonomyAssignment --zotuTable "path2/curatedOruncurated_ZotuTable_file" --blastFile "path2/blastResult_file" --lca_output "my_lca_result" [OPTIONS]
    
    Mandatory arguments if sequences are NOT demultiplexed
      --reads [file]                  Input data name (must be surrounded with quotes);
                                      do NOT specify this option if reads are paired-end as they get identified automatically by default;
                                      reads must be unzipped
      --barcode [file]                Barcode file name; barcode file format must match OBITools requirement;
                                      if multiple barcode files exist (e.g. bc_1.txt; bc_2.txt) it can be specified like this: 'bc*.txt'  

    At least one of the below databases must be specified.
      --blast_db [dir]                Absolute path to local nt databse 
      --custom_db [dir]               Absolute path to custom database
    
    Mandatory arguments if sequences are demultiplexed
      --skipDemux [bool]
      --demuxedInput [file]           A fasta file holding all the demultiplexed sequences;
                                      Format of the sample identifier must match USEARCH requirements
    
    At least one of the below databases must be specified.
      --blast_db [dir]                Absolute path to local nt databse 
      --custom_db [dir]               Absolute path to custom database

    Mandatory and optional arguments for running LCA taxonomy assignment script
      --taxonomyAssignment [bool]
      --zotuTable [dir]               Path to Zotu table file
      --blastFile [dir]               Path to blast result file

      Optional
      --lca_qcov    [num]             Percent of query coverage; Default is 100
      --lca_pid     [num]             Percent of identity; Default is 97
      --lca_diff    [num]             Diff, Default is 1
      --lca_output  [string]          Output file name
    
    Skipping                          Skip any of the mentioned steps
      --skipDemux [bool]              If this option is set, then --demuxedInput [file] must be provided
      --skipFastqc [bool]

    Parameters to run eDNAFlow on Cloud/HPC
      -profile [string]                 Currently can choose between "nimbus" (can be used if user has access to more memory i.e. cloud or HPC)
                                      and "zeus" (it's specific to users who have access to ZEUS - a high-throughput HPC cluster at the Pawsey Supercomputing Centre)
                                      e.g. -profile nimbus
      --bindDir [dir]                   If you run eDNAFlow on Cloud or HPC, you will need to specify this option
                                      On HPC, it usually will be /scratch or /group. On Cloud, it could be your mounted volume.
                                      e.g. --bindDir "/scratch"
    
    General optional parameters
      --help                            Show this help message
      --publish_dir_mode [string]       Choose between symlink (Default), copy, link
      --singularityDir [dir]            Directory where singularity images will be stored
    
    Demultiplexing
      --onlyDemux [bool]
      

    Quality filtering / Merging
      --minQuality [num]              The minimum Phred quality score to apply for quality control of raw sequences; Default is 20
      --minAlignLeng [num]            The minimum alignment length for merging read1 and read2; Default is 12
      --minLen [num]                  The minimum length allowed for sequences; Default is 50

    ZOTU formation
      --minsize [num]                 The minimum abundance; input sequences with lower abundances are removed; Default is 8

    Blast parameters
      --maxTarSeq [num]               the maximum number of target sequences for hits per query to be returned by Blast; Default is 10
      --perc_identity [num]           percentage of identical matches; Default is 95
      --evalue [num]                  expected value for saving blast hits; Default is 1e-3
      --qcov [num]                    the percent of the query that has to form an alignment against the reference to be retained;
                                      Higher values prevent alignments of only a short portion of the query to a reference; Default is 100

    Choice of USEARCH32 vs USEARCH64 
      --mode [str]                   Default is usearch32; for running with 64 version the mode has to be set to --mode 'usearch64'
                                     and below option has to be specified as well
      --usearch64 [dir]              Full path to where usearch64 bit version is stored locally

    LULU
      --lulu [file]                  An R script to run post-clustering curation with default settings of LULU;
                                     This file has been provided and must be present in the same directory as other scripts;
                                     by default eDNAFlow will be looking for this file in the same directory where eDNAFlow.nf is. 

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// setting up variables 

fastQC_input_ch = Channel.fromFilePairs(params.reads, size: -1)
reads_ch      = Channel.fromFilePairs(params.reads, size: -1)
bc_ch         = Channel.fromPath(params.barcode)
minQuality    = params.minQuality
minAlignLeng  = params.minAlignLeng
minLen        = params.minLen
minsize       = params.minsize
maxTarSeq     = params.maxTarSeq
perc_identity = params.perc_identity
evalue        = params.evalue
qcov          = params.qcov
lulu          = file(params.lulu)
mode	      = params.mode
usearch64     = params.usearch64


// Setting up channels & inputs for taxonomy assignment step

params.taxonomyAssignment = false

if (params.taxonomyAssignment) {
	reads_ch = Channel.empty()
  fastQC_input_ch = Channel.empty()
	lca_script = file(params.lca_script)
	zotuTable_ch = Channel.fromPath(params.zotuTable)
	blastFile_ch = Channel.fromPath(params.blastFile)
	lca_qcov = params.lca_qcov
	lca_pid = params.lca_pid
	lca_diff = params.lca_diff
	lca_output = params.lca_output
} else {
	zotuTable_ch = Channel.empty()
	blastFile_ch = Channel.empty()
	lca_script = Channel.empty()
}


// Running only demultiplexing or skipping 

params.onlyDemux = false

params.skipDemux = false

params.demuxedInput = "*.fasta"

params.skipFastqc = false


if (params.skipDemux) {
	reads_ch = Channel.empty()
	Channel.fromPath(params.demuxedInput).set{ dmuxed_relabeled_input_ch }
	name_ch = Channel.value(params.demuxedInput)
}



/*
 * Process 0: FastQC to check the initial quality of raw reads
 */

process '00_fastQC' {
        label 'fastqc'

  publishDir "00_fastQC_${sample_id}", mode: params.publish_dir_mode

  input:
      tuple val(sample_id), path(read) from fastQC_input_ch

  output:
      tuple val(sample_id), path('*_fastqc.{zip,html}') into res_fastQC_ch

  when:
  !params.skipFastqc

  script:
  
  if( read instanceof Path ) {
  """
  fastqc -q ${read}
  """
  } else {
  """
  fastqc -q ${read[0]} ${read[1]}
  """
  }
}

 
/*
 * Process 1: Quality filtering of raw reads
 */

process '01_a_quality_Filtering' {
	label 'adapterRemoval'
 
  publishDir "01_a_quality_Filtering_${sample_id}", mode: params.publish_dir_mode

  input:
      tuple val(sample_id), path(read) from reads_ch

  output:
      tuple val(sample_id), path('*_QF.fastq') into QF_ch 
      
  script:

  if( read instanceof Path ) {   
     """
     AdapterRemoval --threads ${task.cpus} --file1 ${read} \
                    --trimns --trimqualities \
                    --minquality ${minQuality} \
                    --basename ${sample_id}

     mv ${sample_id}.truncated ${sample_id}_QF.fastq

    """
  }

// if reads are paired-end then merge 
  else {  
     """
     AdapterRemoval --threads ${task.cpus} --file1 ${read[0]} --file2 ${read[1]} \
                    --collapse --trimns --trimqualities \
                    --minquality $minQuality \
                    --minalignmentlength $minAlignLeng \
                    --basename ${sample_id}

    mv ${sample_id}.collapsed ${sample_id}_QF.fastq  
    """
  }
}

QF_ch.into { QFres_4dmux_ch; QFres_4fastQC_ch}
barcode_reads_mixed = QFres_4dmux_ch.combine(bc_ch)
                           .view()



process '01_b_fastQC' {
        label 'fastqc'

  publishDir "01_b_fastQC_${sample_id}", mode: params.publish_dir_mode

  input:
      tuple val(sample_id), path("${sample_id}_QF.fastq") from QFres_4fastQC_ch

  output:
      path("*_fastqc.{zip,html}") into last_fastQC_ch

  when:
  !params.skipFastqc

  script:
  """
  fastqc -q ${"${sample_id}_QF.fastq"}
  """
}
   


/*
 * Process 2: Assigning reads to samples/demultiplexing
 */


process '02_assigned_dmux' {
  label 'obitools'

  publishDir "02_assigned_dmux_${sample_id}_${barcode.baseName}", mode: params.publish_dir_mode

  input:
      tuple val(sample_id), path(read), path(barcode) from barcode_reads_mixed
    

  output:
      tuple val(sample_id), val("${barcode.baseName}"), path("*_Dmux.fastq") into dmux_ch

  script:
  """
  ngsfilter -t ${barcode} -u "orphan.fastq" ${read} > "${sample_id}_${barcode.baseName}_QF_Dmux.fastq"
  """
}

grouped_demux = dmux_ch.groupTuple()

/*
 * Process 3: Cat multiple file after demux, Filtering sequences shorter than min length and single barcoded reads
 */

process '03_Length_filtered' {
  label 'obitools'

  publishDir "03_Length_filtered_${sample_id}", mode: params.publish_dir_mode

  input: 
  tuple val(sample_id), val(barcode_files), path(fastq_files) from grouped_demux
  output:
  tuple val(sample_id), path('*_QF_Dmux_minLF.fastq') into lenFilt_ch


  script:
  """
  cat ${fastq_files} > "${sample_id}_QF_Dmux.fastq" 
  obigrep -l $minLen -p 'forward_tag is not None and reverse_tag is not None' ${sample_id}_QF_Dmux.fastq > ${sample_id}_QF_Dmux_minLF.fastq
 
  """

}


/*
 * Process 4: Split the file based on samples
 */


process '04_splitSamples_relabel_Cat' {
  label 'obitools'

  publishDir "04_splitSamples_${sample_id}", mode: params.publish_dir_mode

  input:
  tuple val(sample_id), path(fastqs) from lenFilt_ch
  output:
  tuple val(sample_id), path("$sample_id/*.fastq") into split_ch
  script:
  """
  mkdir ${sample_id}
  obisplit -t sample -u "noSampleID.fastq" $fastqs
  mv *.fastq ${sample_id}
  mv ${sample_id}/$fastqs ..
  mv ${sample_id}/noSampleID.fastq  noSampleID.fastq.ignore
  """
}


/*
 * Process 5: Relabel file for usearch
 */


process '05_relabel_Cat' {

   label 'usearch'
   
   publishDir "${task.process}_${sample_id}", mode: params.publish_dir_mode

   input:
   tuple val(sample_id), path(fastqs) from split_ch

   output:
   tuple val(sample_id), path("*.relabeled.fastq"), path("CountOfSeq.txt"), path("*_relabeled4Usearch.fastq") into addition_ch 
   tuple val(sample_id), path("*_upper.fasta") into relabel_ch


   script:
   if(mode == 'usearch32')
	   """
	   for files in ${fastqs}
	   do
	   label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
	   usearch -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq 
	   done
	 
	   for files in *.relabeled.fastq
	   do
	   name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
	   echo \${name} >> CountOfSeq.txt
	   grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
	   done 

	   cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"
	  
	   usearch -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

	   usearch -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}.fasta
	   
	   awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta
	 
	   """
   else if(mode == 'usearch64')
   	   """
	   for files in ${fastqs}
           do
           label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
           $usearch64 -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq
           done
  
           for files in *.relabeled.fastq
           do
           name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
           echo \${name} >> CountOfSeq.txt
           grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
           done

           cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"
   
           $usearch64 -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

           $usearch64 -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}.fasta
    
           awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta
	   """
 }

// redirecting channel choice depending on whether or not --skipDemux parameter is set
demux_channel = (params.skipDemux ? name_ch.combine(dmuxed_relabeled_input_ch) : relabel_ch)

/*
 * Process 6: Dereplication, ZOTUs creation, ZOTU table creation
 */


process '06_Uniques_ZOTUs' {
  label 'usearch'
  echo true

  publishDir "${task.process}_${sample_id}", mode: params.publish_dir_mode
  input:
    tuple val(sample_id), path(upper_fasta) from demux_channel 

  output:
    tuple val(sample_id), path("${sample_id}_Unq.fasta"), path("${sample_id}_zotus.fasta"), path("zotuTable.txt") into zotu_ch

  when:
  !params.onlyDemux

  script:
  if(mode == 'usearch32')
	   """
	   usearch -fastx_uniques ${upper_fasta} -sizeout -fastaout "${sample_id}_Unq.fasta"

	   usearch -unoise3 "${sample_id}_Unq.fasta"  -zotus "${sample_id}_zotus.fasta" -tabbedout "${sample_id}_Unq_unoise3.txt" -minsize $minsize
	    
	   usearch -otutab ${upper_fasta} -zotus ${sample_id}_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
	   """
   else if(mode == 'usearch64')
	   """
           $usearch64 -fastx_uniques ${upper_fasta} -sizeout -fastaout "${sample_id}_Unq.fasta"

           $usearch64 -unoise3 "${sample_id}_Unq.fasta"  -zotus "${sample_id}_zotus.fasta" -tabbedout "${sample_id}_Unq_unoise3.txt" -minsize $minsize

           $usearch64 -otutab ${upper_fasta} -zotus ${sample_id}_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
           """
}


/*
 * Process 7: Blast
 */

process '07_blast' {
  label 'blast'

  publishDir "${task.process}_${sample_id}", mode: params.publish_dir_mode
  input:
     tuple val(sample_id), path(a), path(zotus_fasta), path(zotuTable) from zotu_ch
  output:
     tuple val(sample_id), path("${sample_id}_blast_Result.tab"), path(zotuTable), path("match_list.txt") into blast_ch
  script:
      """
      blastn -db "${params.blast_db} ${params.custom_db}" \
             -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
             -perc_identity $perc_identity -evalue $evalue \
             -best_hit_score_edge 0.05 -best_hit_overhang 0.25 \
             -qcov_hsp_perc $qcov -max_target_seqs $maxTarSeq \
             -query ${zotus_fasta} -out ${sample_id}_blast_Result.tab

      makeblastdb -in ${zotus_fasta} -parse_seqids -dbtype nucl -out ${sample_id}_zotus

      blastn -db ${sample_id}_zotus \
             -outfmt "6 qseqid sseqid pident" \
             -out match_list.txt -qcov_hsp_perc 80 \
             -perc_identity 84 -query ${zotus_fasta}

      """
}

/*
 * Process 8: LULU curation
 */


process '08_lulu' {
  label 'lulu'
  publishDir "${task.process}_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(a), path(zotuTable), path(match_list) from blast_ch
    file script from lulu

  output:
    tuple path("curated_zotuTable.tab"), path("lulu_zotu_map.tab") into lulu_ch

  script:
    """
    Rscript $lulu
    """
}



/*
 * Process 9: Taxonomy assignment with LCA
 */

process '09_taxonomyAssigned' {
  label 'lca_python3'
  publishDir "${task.process}_${lca_output}_qCov${lca_qcov}_id${lca_pid}_diff${lca_diff}"
  
  input:
    tuple path(table), path(blastRes) from zotuTable_ch.combine(blastFile_ch)
    file lcaScript from lca_script

  output:
    tuple path("interMediate_res.tab"), path("${lca_output}_qCov${lca_qcov}_id${lca_pid}_diff${lca_diff}.tab") into last_ch


  script:
    """
    python3 $lca_script ${table} ${blastRes} ${lca_qcov} ${lca_pid} ${lca_diff} ${lca_output}_qCov${lca_qcov}_id${lca_pid}_diff${lca_diff}.tab
    """
}  
  	
