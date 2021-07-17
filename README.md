# nf-tb
ABMicroBioinf/nf-tb is a bioinformatics pipeline designed to analyze Mycobacterium tuberculosis genomes. It takes quality controlled reads and kraken2 output report as input. The pipeline aligns filtered MTBC reads to the H37Rv reference using bowtie2, BWA or minimap2 and then calls variants using bcftools. These variants are then compared to a drug-resistance database. The drug resistance classes and calls were done using tb-profiler. The pipeline also uses Snippy to do rapid SNP calling over the Mycobacterium tuberculosis whole genome. Snippy finds both substitutions (snps) and insertions/deletions (indels) over the whole genome.The identified variants were furthur filtered to exclude 1) the snps in PE/PPE, UVP, and repeat / insertion sequence sites; 2) the minimum alternate allele percentage to accept < 90.0; 3) the minimum read depth (coverage) < 30. The filtered vcf format snps are used to build phylogenetic trees. 

# Third-party software
This pipeline are depending on a number of the third-party software. Please install the following 3rd party dependencies and make sure they are on your system path
* nextflow
* seqtk
* bwa
* samtools
* bcftools
* snippy
* tb-profiler
* tb_variant_filter
* vcf2phylip
* RaxML

# Quick Start
```
 nextflow -bg run /data/deve/workflows/nf-tb/main.nf --input_seq_dir ../nf-qc/outdir/qc/seqs --input_kraken2_dir ../nf-qc/outdir/classifier --outdir outdir --min_samples_locus 1 
```
