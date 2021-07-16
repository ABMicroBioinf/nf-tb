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
  --min_samples_locus       Minimum of samples required to be present at a locus (default=1)
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
process get_mtbc_reads {
  
  publishDir "${params.outdir}/mtbc", mode: "copy" 

  input:
  set val(name), file(kraken2_output_file) from kraken2_output_ch
  set val(name), file(reads) from reads_ch
  
  output:
  set val(name), file("*.mtbc.fastq") into mtbc_reads_ch

  script:
  
   if (params.type == 'single-end') {

    """
    taxonFilter.pl $kraken2_output_file  > ${name}_mtbc_accession.txt
    seqtk subseq ${reads} ${name}_mtbc_accession.txt > ${name}.mtbc.fastq
    """
  }
  else{
    """
    taxonFilter.pl $kraken2_output_file  > ${name}_mtbc_accession.txt
    seqtk subseq ${reads[0]} ${name}_mtbc_accession.txt > ${name}_R1.mtbc.fastq
    seqtk subseq ${reads[1]} ${name}_mtbc_accession.txt > ${name}_R2.mtbc.fastq
    """

  }
}

process indexRefGenome {
    //publishDir "${params.outdir}/bwa_index", mode: "copy"
    input:
      path genome from params.ref_genome
    output:
      file '*' into bwa_index

    script:
    """
    bwa index ${genome}
    """
}

process mapping_mtbc_and_filter {

  label "small_mem"

  publishDir "${params.outdir}/bam", mode: "copy" 

  input:
  set val(name), file(reads) from mtbc_reads_ch
  path index from bwa_index
  path genome from params.ref_genome
  
  output:
  set val(name), file("*.sort.bam") into bam_to_tbprofiling
  set val(name), file("*.sort.bam") into bam_to_snippy
  file "*" into mapped_ch
 
  script:

   if (params.type == 'single-end') {
    """
     bwa mem -t ${task.cpus} \
     -c 100 -R '@RG ID:${name} SM:${name} PL:illumina' -M -T 50 \
     ${genome} ${reads} \
     | samtools view -@ ${task.cpus} -Sb1 -q ${params.mq} - \
     | samtools sort -@ ${task.cpus} -o ${name}.sort.bam -

     samtools index -@ ${task.cpus} ${name}.sort.bam
    """
  }

   else {
    """
     bwa mem -t ${task.cpus} \
     -c 100 -R '@RG\\tID:${name}\\tSM:${name}\\tPL:illumina' -M -T 50 \
     ${genome} ${reads[0]} ${reads[1]} \
     | samtools view -@ ${task.cpus} -Sb1 -q ${params.mq} - \
     | samtools sort -@ ${task.cpus} -o ${name}.sort.bam -

     samtools index -@ ${task.cpus} ${name}.sort.bam
    """
  }
}

process run_snippy {
  
  label "small_mem"
  publishDir "${params.outdir}/snippy", mode: 'copy'
  
  input:
  set val(name), file(bam) from bam_to_snippy
  path genome from params.ref_genome

  output:
  file "*.whole.filter.vcf.gz" into build_tree_ch
  file "*.whole.filter.vcf.gz.csi" into csi_ch 

  script:
  
    """
    echo ${bam}
    snippy --outdir ./ --prefix ${name}.whole --ref ${genome} --bam ${bam} --rgid ${name} --minfrac 0.1 --cpus ${task.cpus} --force
    tb_variant_filter  -I  --indel_window_size 5 --region_filter pe_ppe,uvp  -P --min_percentage_alt 90  -D --min_depth 30 ${name}.whole.vcf ${name}.whole.filter.vcf
    bgzip -c ${name}.whole.filter.vcf > ${name}.whole.filter.vcf.gz
    bcftools index ${name}.whole.filter.vcf.gz
    """

}

process tbprofiling {
  
  label "small_mem"
  publishDir "${params.outdir}/tbprofiler", mode: 'copy'
  
  input:
  set val(name), file(bam) from bam_to_tbprofiling
  
  output:
  file "vcf/*" into vcf_ch
  file "results/*" into results_ch

  script:
  
    """
    echo ${bam}

    tb-profiler profile --threads ${task.cpus} --bam ${bam} \
      --prefix ${name} --csv --txt --verbose 2 
    
    """

}

Channel.fromPath("${params.outdir}/tbprofiler/results", type: 'dir')
  .ifEmpty { exit 1, "Cannot find  directory: ${params.outdir}tbprofiler/results"}
  .set{ results_dir_ch}
 

process tbprofiler_collate {
  label "small_mem"

  publishDir "${params.outdir}/tbprofiler/collate", mode: "copy" 

  input:
  file dir from results_dir_ch
  path result_file from results_ch.collect()
  
  output:
  file "tbprofiler_collate*" into collate_results


  script:
  """
  
  echo $dir
  tb-profiler collate -d $dir -p tbprofiler_collate 
  """
 
}

process build_tree {
  label "small_mem"

  publishDir "${params.outdir}/tree", mode: "copy" 

  input:
  //file vcf_whole from whole_vcf_ch.collect()
  file  vcf_whole from build_tree_ch.collect()
  file csi from  csi_ch.collect()
  output:
  file("*") into tree_out_ch

  script:
  """
  echo ${vcf_whole}
  bcftools merge -Oz -o merged.whole.vcf *.vcf.gz
  vcf2phylip.py -i merged.whole.vcf -m ${params.min_samples_locus}
  raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -s merged.whole.min${params.min_samples_locus}.phy -n whole.tree -T ${task.cpus}  -x 0123 -N 100 
  """
  //
  
  //https://digital.csic.es/bitstream/10261/181837/2/2019_J%20Infect%20Dis_suppl_supplementary_material.pdf
}

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
