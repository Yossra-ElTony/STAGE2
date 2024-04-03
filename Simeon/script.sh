#!/bin/bash
# create a parent folder purposely for variant calling
mkdir variant_calling
cd variant_calling

# create sub-folders to contain the reads and the reference genome
mkdir read_data
mkdir genome_data

# move into the read_data directory and downloads the read1 and read2
cd read_data
wget -O read1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1
wget -O read2.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1

# move into genome_data and downloads the reference genome
cd ../genome_data
wget -O reference.fasta https://zenodo.org/records/10886725/files/Reference.fasta?download=1

# print parts of the files to have a fair idea of they look like
head reference.fasta
head ../read_data/read1.fastq.gz
head ../read_data/read2.fastq.gz

# change into the parent directory to install fastqc
cd ..
sudo apt-get install fastqc

# move into the read directory and run fastqc to do quality control
cd read_data
ls
fastqc *.gz

# install multiqc to aggregate the two files from fastqc
sudo apt-get install multiqc
ls
multiqc *.zip

# install fastp and run to do trimming on the reads
sudo apt-get install fastp
fastp -i read1.fastq.gz -I read2.fastq.gz -o trimmed_read1.fastq.gz -O trimmed_read2.fastq.gz -h fastp_report.html 
ls

# install bwa to do alignment or mapping
sudo apt-get install bwa

# index the reference genome to aid in location while mapping
cd ../genome_data
bwa index reference.fasta

# make a new directory bam_data to hold all .bam files
cd ..
mkdir bam_data

# install samtools to aid mapping and print all the contents in the various files
sudo apt-get install samtools
samtools --help
ls genome_data
ls read_data

# run bwa and pipe to samtools to finish alignment
bwa mem -t 8 genome_data/reference.fasta read_data/trimmed_read1.fastq.gz read_data/trimmed_read2.fastq.gz | samtools view -h -b -o bam_data/aligned.bam
ls bam_data
cd bam_data

# use samtools flagstat to output the mapping statistics 
samtools flagstat aligned.bam>mapping_statistics.txt
cat mapping_statistics.txt

# use samtools -F for flagging or filtering all the reads that were not mapped
samtools view --help
samtools view -b -F 4 aligned.bam > filtered_aligned.bam

# sort and index the filtered reads to make variant calling fast and easy
samtools sort filtered_aligned.bam > sorted_filtered_aligned.bam
samtools index sorted_filtered_aligned.bam

# install bcftools to do variant calling
sudo apt-get install bcftools
bcftools
bcftools call --help

# move to the parent directory and creates a new folder vcf_data store files associated with variant calling
cd ..
mkdir vcf_data

# print the file contents of the two folders involved in variant calling
ls bam_data
ls genome_data

# run bcftools mpileup to pileup all the alleles (similar genes at the same location) and output the bcf file
bcftools mpileup -O b -o vcf_data/pileup_file.bcf -f genome_data/reference.fasta bam_data/sorted_filtered_aligned.bam

# call the variants and output them in a .vcf format
bcftools call --ploidy 1 -mv -O v -o vcf_data/variant_calls.vcf vcf_data/pileup_file.bcf
ls vcf_data

# output the variant file
less vcf_data/variant_calls.vcf
