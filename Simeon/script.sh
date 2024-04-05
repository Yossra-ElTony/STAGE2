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
# head reference.fasta
# head ../read_data/read1.fastq.gz
# head ../read_data/read2.fastq.gz

# change into the parent directory to install fastqc
cd ..
sudo apt-get install fastqc

# move into the read directory and run fastqc to do quality control
cd read_data
# ls
fastqc *.gz

# install multiqc to aggregate the two files from fastqc
sudo apt-get install multiqc
# ls
multiqc *.zip

# install fastp and run to do trimming on the reads
sudo apt-get install fastp
fastp -i read1.fastq.gz -I read2.fastq.gz -o trimmed_read1.fastq.gz -O trimmed_read2.fastq.gz -h fastp_report.html 
# ls

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
# samtools --help
# ls genome_data
# ls read_data

# run bwa and pipe to samtools to finish alignment
bwa mem -t 8 genome_data/reference.fasta read_data/trimmed_read1.fastq.gz read_data/trimmed_read2.fastq.gz | samtools view -h -b -o bam_data/aligned.bam
# ls bam_data
cd bam_data

# use samtools flagstat to output the mapping statistics 
samtools flagstat aligned.bam>mapping_statistics.txt
# cat mapping_statistics.txt

# use samtools -F for flagging or filtering all the reads that were not mapped
# samtools view --help
samtools view -b -F 4 aligned.bam > filtered_aligned.bam

# sort and index the filtered reads to make variant calling fast and easy
samtools sort filtered_aligned.bam > sorted_filtered_aligned.bam
samtools index sorted_filtered_aligned.bam

# install bcftools to do variant calling
sudo apt-get install bcftools

# move to the parent directory and creates a new folder vcf_data store files associated with variant calling
cd ..
mkdir vcf_data

# print the file contents of the two folders involved in variant calling
# ls bam_data
# ls genome_data

# run bcftools mpileup to pileup all the alleles (similar genes at the same location) and output the bcf file
bcftools mpileup -O b -o vcf_data/pileup_file.bcf -f genome_data/reference.fasta bam_data/sorted_filtered_aligned.bam

# call the variants and output them in a .vcf format
bcftools call --ploidy 1 -mv -O v -o vcf_data/variant_calls.vcf vcf_data/pileup_file.bcf
ls vcf_data

# output the variant file
# less vcf_data/variant_calls.vcf



# start the loop on the other dataset

# make the parent directory and the read folder
mkdir variant_calling1 && cd $_
mkdir read_data
cd read_data

# download all reads into read data folder
wget -O acbarrie_read1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
wget -O acbarrie_read2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz

wget -O alsen_read1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
wget -O alsen_read2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz

wget -O baxter_read1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz
wget -O baxter_read2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz

wget -O chara_read1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
wget -O chara_read2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz

wget -O drysdale_read1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
wget -O drysdale_read2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz

# create the reference (genome) folder and download the reference genome into it
cd ..
mkdir genome_data
cd genome_data

wget -O reference.fasta https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

cd ..
cd read_data
# install and run fastqc to do quality control
sudo apt-get install fastqc -y
for file in *.fastq.gz; do fastqc $file; done

# install and run multiqc to aggregate all the files from fastqc
sudo apt-get install multiqc -y
multiqc *.zip

# install fastp and run to do trimming on the reads
sudo apt-get install fastp -y

# do a for loop to run through each of the reads and trim
samples=("acbarrie" "alsen" "baxter" "chara" "drysdale")
for sample in "${samples[@]}"; do
        fastp -i "${sample}_read1.fastq.gz" -I "${sample}_read2.fastq.gz" -o "${sample}_trimmed_read1.fastq.gz" -O "${sample}_trimmed_read2.fastq.gz" -h "${sample}_fastp_report.html"
done

# install bwa for indexing and mapping
sudo apt-get install bwa -y

# index the reference genome to aid in location while mapping
cd ../genome_data
bwa index reference.fasta

# make a new directory bam_data to hold all .bam files
cd ..
mkdir bam_data

# install samtools to aid mapping and print all the contents in the various files
sudo apt-get install samtools -y

# do a for loop to run through each of the reads and perform alignment with bwa and pipe to samtools
# samples=("acbarrie" "alsen" "baxter" "chara" "drysdale") --not necessary do call the array again
for sample in "${samples[@]}"; do
        bwa mem -t 8 genome_data/reference.fasta read_data/"${sample}_trimmed_read1.fastq.gz" read_data/"${sample}_trimmed_read2.fastq.gz" | samtools view -h -b -o bam_data/"${sample}_aligned.bam"
done

# use samtools flagstat to output the mapping statistics
cd bam_data
for sample in "${samples[@]}"; do
        samtools flagstat "${sample}_aligned.bam" > "${sample}_mapping_statistics.txt"
done

# use samtools -F for flagging or filtering all the reads that were not mapped
for sample in "${samples[@]}"; do
        samtools view -b -F 4 "${sample}_aligned.bam" > "${sample}_filtered_aligned.bam"
done

# sort and index the filtered reads to make variant calling fast and easy
for sample in "${samples[@]}"; do
        samtools sort "${sample}_filtered_aligned.bam" > "${sample}_sorted_filtered_aligned.bam"
        samtools index "${sample}_sorted_filtered_aligned.bam"
done

# install bcftools to do variant calling
sudo apt-get install bcftools -y

# move to the parent directory and creates a new folder vcf_data store files associated with variant calling
cd ..
mkdir vcf_data

# run bcftools mpileup to pileup all the alleles (similar genes at the same location) and output the bcf file
for sample in "${samples[@]}"; do
        bcftools mpileup -O b -o vcf_data/"${sample}_pileup_file.bcf" -f genome_data/reference.fasta bam_data/"${sample}_sorted_filtered_aligned.bam"
        # call the variants and output them in a .vcf format
        bcftools call -mv -O v -o vcf_data/"${sample}_variant_calls.vcf" vcf_data/"${sample}_pileup_file.bcf"
done
