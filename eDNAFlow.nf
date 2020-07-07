#!/usr/bin/env nextflow


// setting up variables 

reads_ch = Channel.fromFilePairs(params.reads, size: -1)
bc_ch = Channel.fromPath(params.barcode)
minQuality = params.minQuality
minAlignLeng = params.minAlignLeng
minLen = params.minLen
minsize = params.minsize
maxTarSeq = params.maxTarSeq
perc_identity = params.perc_identity
evalue = params.evalue
qcov = params.qcov
lulu = file(params.lulu)
mode = params.mode
usearch64 = params.usearch64


/*
 * Process 1: Quality filtering of raw reads
 */

process '01_a_quality_Filtering' {
	label 'adapterRemoval'
 
  publishDir "01_a_quality_Filtering_${sample_id}"

  input:
      tuple val(sample_id), path(read) from reads_ch

  output:
      tuple val(sample_id), path('*_QF.fastq') into QC_ch 
      
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

barcode_reads_mixed = QC_ch.combine(bc_ch)
                           .view()
/*
 * Process 2: Assigning reads to samples/demultiplexing
 */


process '01_b_assigned_dmux' {
  label 'obitools'

  publishDir "01_b_assigned_dmux_${sample_id}_${barcode.baseName}"

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

process '02_Length_filtered' {
  label 'obitools'

  publishDir "02_Length_filtered_${sample_id}"

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


process '03_splitSamples_relabel_Cat' {
  label 'obitools'

  publishDir "03_splitSamples_${sample_id}"

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
 * Process 5: Relable file for usearch
 */


process '04_relabel_Cat' {

   label 'usearch'
   
   publishDir "${task.process}_${sample_id}"

   input:
   tuple val(sample_id), path(fastqs) from split_ch

   output:
   tuple val(sample_id), path("*.relabeled.fastq"), path("CountOfSeq.txt"), path("*_relabeled4Usearch.fastq"), path("*_upper.fasta") into relabel_ch
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

/*
 * Process 6: Dereplication, ZOTUs creation, ZOTU table creation
 */


process '05_Uniques_ZOTUs_' {
  label 'usearch'

  publishDir "${task.process}_${sample_id}"
  input:
    tuple val(sample_id), path(a), path(b), path(c), path(upper_fasta) from relabel_ch
  output:
    tuple val(sample_id), path("${sample_id}_Unq.fasta"), path("${sample_id}_zotus.fasta"), path("zotuTable.txt") into zotu_ch

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

process '06_blast' {
  label 'blast'

  publishDir "${task.process}_${sample_id}"
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


process '07_lulu' {
  label 'lulu'
  publishDir "${task.process}_${sample_id}"

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
