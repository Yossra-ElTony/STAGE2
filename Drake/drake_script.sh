#!/bin/bash
#perform QC on raw reads

echo **************PERFORMING QC*****************
mkdir fastqc

fastqc *.gz -o fastqc

#Aggregate the QC resylts together
mkdir multiqc

multiqc fastqc -o multiqc

echo ******************TRIMMING READS***************

#Trim raw reads
mkdir trimmed

fastp -i ERR8774458_1.fastq.gz -o trimmed/R1.fastq.gz -I ERR8774458_2.fastq.gz -O trimmed/R2.fastq.gz

echo *****************MAPPING**************************
#mapping

mkdir mapping
# first we index the ref_genome with bwa index


bwa index Reference.fasta 

#then we do alignment of our forward and reverse to ref_genome with bwa mem

bwa mem  Reference.fasta trimmed/R1.fastq.gz trimmed/R2.fastq.gz > mapping/mappped.sam

#covert sam file to bam format using sam tools


samtools view -b mapping/mappped.sam > mapping/mappped.bam

samtools view -b -F 0xc mapping/mappped.bam > mapping/filtered_mappped.bam
#then sort the bam file using sam tools

samtools sort -n mapping/filtered_mappped.bam > mapping/sorted_map.bam

samtools fixmate -m mapping/sorted_map.bam mapping/sorted_fixmate_map.bam

samtools sort mapping/sorted_fixmate_map.bam -o mapping/sortedp_fixmate_map.bam

samtools markdup -r mapping/sortedp_fixmate_map.bam mapping/dedup.bam

samtools index mapping/dedup.bam

echo ******************VARIANT CALLING***********************
#Perform variant calling
mkdir freebayes_variantcalls
mkdir bcftools_variantcalls

freebayes -f Reference.fasta -b mapping/dedup.bam --vcf freebayes_variantcalls/bayescalls.vcf

#zip the vcf file
bgzip freebayes_variantcalls/bayescalls.vcf

#index the zipped vcf file
bcftools index freebayes_variantcalls/bayescalls.vcf.gz

#view variants _snps
bcftools view -v snps freebayes_variantcalls/bayescalls.vcf.gz| grep -v -c '^#'

#view variants _indels

bcftools view -v indels freebayes_variantcalls/bayescalls.vcf.gz| grep -v -c '^#'







