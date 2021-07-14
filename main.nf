#! /usr/bin/env nextflow

def helpMessage() {
  log.info """
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run main.nf --inputdir "inputDirName" --outdir "outDirName"

  Arguments:
  --input_seq_dir           The direcoty containing the input fastq or fastq.gz files (default is the current direcoty)
  --input_kraken2_dir       The direcoty containing the kraken2 output (default is the current direcoty)
  --outdir                  The directory to place processing results (default is ./outdir)
  --paired_end              the fastq filename pattern (default *_R{1,2}*.fastq.gz)
  --single_end              the fastq filename pattern (default *.fastq.gz)
  --ref_genome              the fasta format reference genome for mapping (default tbdb.fasta )
  --mq                      mapping quality cutoff to filter the uniquelly mapped reads (default 10)
  --type                    paired-end or single-end (default: paired-end)
  --help                    This usage statement.
  """.stripIndent()
}
//--host                     The host reference fasta file (defualt is ./ref/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz)
seq_status = "raw"

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

kraken2_output_ch = Channel
  .fromPath(params.input_kraken2_dir + '/' + params.kraken2_output)
  .ifEmpty { exit 1, "Cannot find  kraken2 report matching: ${params.kraken2_output}"}
  .map { file -> tuple(file.simpleName, file) }
  

if (params.type == 'single-end') {
    reads_ch = Channel
                        .fromPath(params.input_seq_dir + '/' + params.single_end)
                        .map { file -> tuple(file.simpleName, file) }
                
} else {
  Channel
      .fromFilePairs(params.input_seq_dir + '/' + params.paired_end, size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.paired_end}"}
      .set {reads_ch}
}


/*
Step 1
*/
process mtbc_namelist {
  
  publishDir "${params.outdir}", mode: "copy" 

  
  input:
  set val(name), file(kraken2_output_file) from kraken2_output_ch
  
  output:
  file "*_mtbc_accession.txt" into mtbc_namelist_ch

  script:
  
  """
  taxonFilter.pl $kraken2_output_file  > ${name}_mtbc_accession.txt
  """
}

/*
Step 2
*/
process mtbc_reads {
  
  publishDir "${params.outdir}", mode: "copy" 

  
  input:
  file(namelist) from mtbc_namelist_ch
  set val(name), file(reads) from reads_ch
  
  output:
  set val(name), file("*.mtbc.fastq") into mtbc_reads_ch

  script:
  prefix = name  
  if (params.type == 'single-end') {

    """
    seqtk subseq ${reads} ${namelist} > ${prefix}.mtbc.fastq
    """
  }
  else{
    """
    seqtk subseq ${reads[0]} ${namelist} > ${prefix}_R1.mtbc.fastq
    seqtk subseq ${reads[1]} ${namelist} > ${prefix}_R2.mtbc.fastq
    """

  }
}

process indexReference {
    
    publishDir "${params.outdir}/bwa_index", mode: "copy"

    input:
      path genome from params.ref_genome
    output:
      file '*' into bwa_index

    script:
    """
    bwa index ${genome}
    """
}

process mapping_and_filter {

  label "small_mem"

  publishDir "${params.outdir}", mode: "copy" 

  input:
  set val(name), file(reads) from mtbc_reads_ch
  path index from bwa_index.first()
  path genome from params.ref_genome
  
  output:
  set val(name), file("*.sorted.bam") into bam_ch
  file "*.sorted.bam.bai" into bam_index_ch
 
  script:

   if (params.type == 'single-end') {
    """
     bwa mem -t ${task.cpus} \
     -c 100 -R '@RG\\tID:ERR036228\\tSM:ERR036228\\tPL:illumina' -M -T 50 \
     ${genome} ${reads} \
     | samtools view -@ ${task.cpus} -Sb1 -q ${params.mq} - \
     | samtools sort -@ ${task.cpus} -o ${name}.sorted.bam -

     samtools index -@ ${task.cpus} ${name}.sorted.bam
    """
  }

   else {
    """
     bwa mem -t ${task.cpus} \
     -c 100 -R '@RG\\tID:ERR036228\\tSM:ERR036228\\tPL:illumina' -M -T 50 \
     ${genome} ${reads[0]} ${reads[1]} \
     | samtools view -@ ${task.cpus} -Sb1 -q ${params.mq} - \
     | samtools sort -@ ${task.cpus} -o ${name}.sorted.bam -

     samtools index -@ ${task.cpus} ${name}.sorted.bam
    """
  }
}

process call_tbprofiler_on_bam {
  
  label "small_mem"

  publishDir "${params.outdir}", mode: "copy" 

  input:
  set val(name), file(bam) from bam_ch
  
  output:
  file("*") into tbprofiler_output_ch

 //samtools sort -m 10G -@ ${task.cpus} -o ${name}.sorted.bam -) 3>&1 1>&2 2>&3 \
  script:
  
    """
    tb-profiler profile --threads ${task.cpus} --bam ${bam} --prefix ${name} --dir tbprofiler \
      --csv --txt  --call_whole_genome --verbose 2 
    """
  

}
/*
process call_tbprofiler_on_reads {
  
  label "small_mem"

  publishDir "${params.outdir}", mode: "copy" 

  input:
  set val(name), file(reads) from mtbc_reads_ch
  
  output:
  file("*") into tbprofiler_output_ch

  script:
  
    """
    tb-profiler  profile -t ${task.cpus} \
      -1 ${reads[0]} -2 ${reads[1]} \
      --prefix ${name} \
      --calling_params "-q ${params.mq}"  \
      --call_whole_genome \
      --dir tbprofiler --csv --txt 
    """
  

}

*/


workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
