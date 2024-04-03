#!/bin/bash

#create directory of stage 2 task
mkdir STAGE2

# go inside the directory
cd STAGE2/

#running pipeling on original dataset
mkdir original-dataset

#go inside the directory
cd original-dataset/

#download the original dataset
echo "Downoading R1 and R2 of original dataset"
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz

#create directory of reference of original dataset
mkdir ref


# go inside the directory of reference
cd ref

#ownload the reference of original dataset
echo "Downloading the ref for original datasets"
wget https://zenodo.org/records/10886725/files/Reference.fasta

cd ..

#conduct quality control by tool; fastqc
echo "Conducting Qulaity control by tool; fastqc"
for f in *.fastq.gz; do fastqc -t 1 -f fastq -noextract $f; done

#conduct trimming by tool; fastp
echo "Conducting Trimming by tool; fastp"
fastp -i *_1.fastq.gz -I *_2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq --qualified_quality_phred 20 --html report.html --json report.json

#conduct quality control by tool fastqc after trimming
echo "Conducting Qulaity control after trimming"
fastqc trimmed_*.fastq

#index the reference of original data
echo "Indexing the reference of original datasets"
cd ref
bwa index -a bwtsw Reference.fasta
cd ..

#conduct the mapping by tool; bwa
echo "Conducting the mapping by tool; bwa; convert to BAM, sort and index"
bwa mem ref/Reference.fasta trimmed_R1.fastq trimmed_R2.fastq > ERR8774458.sam

#Convert the SAM file into a BAM file
samtools view -hbo ERR8774458.bam ERR8774458.sam

#Sort the BAM file by position in genome
samtools sort ERR8774458.bam -o ERR8774458.sorted.bam

#Index the sorted BAM file
samtools index ERR8774458.sorted.bam

#Conduct varient calling by tool; bcftools
echo "Conducting the Varient calling by tool; bcftools"
bcftools mpileup -Ou -f ref/Reference.fasta ERR8774458.sorted.bam | bcftools call -Ov -mv > ERR8774458.vcf

cd ..
#==========================================


# More datasets
mkdir bigger-analysis
cd bigger-analysis

#1-download data

echo "Downoading R1 and R2, ACBarrie"
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz

echo "Downoading R1 and R2, Alsen"
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz

echo "Downoading R1 and R2, Baxter"
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz

echo "Downoading R1 and R2, Chara"
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz

echo "Downoading R1 and R2, Drysdale"
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz

#download Ref
mkdir ref
cd ref
echo "Downloading the ref for additional datasets"
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

cd ..


#2-Qulaity control
echo "Conducting Qulaity control"

mkdir qc-output

for f in *.fastq.gz; do fastqc -t 1 -f fastq -noextract $f -o qc-output; done



#3-Trimming
echo "Conducting Trimming"

mkdir trimming-output

SAMPLES=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

for SAMPLE in "${SAMPLES[@]}"; do fastp -i "$PWD/${SAMPLE}_R1.fastq.gz" -I "$PWD/${SAMPLE}_R2.fastq.gz" -o "trimming-output/${SAMPLE}_trimmed_R1.fastq.gz" -O "trimming-output/${SAMPLE}_trimmed_R2.fastq.gz" --html "trimming-output/${SAMPLE}_fastp.html"; done




#4-Qulaity control after trimming
echo "Conducting Qulaity control after trimming" 

mkdir qc-after-trimming-output

for f in trimming-output/*.fastq.gz; do fastqc -t 1 -f fastq -noextract $f -o qc-after-trimming-output; done




#5-Alignment
echo "Indexing the reference of additional datasets"
cd ref
bwa index -a bwtsw reference.fasta
cd ..

echo "Conducting the mapping"
mkdir alignment-output

SAMPLES=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

for SAMPLE in "${SAMPLES[@]}"; do
    bwa mem ref/reference.fasta trimming-output/${SAMPLE}_trimmed_R1.fastq.gz trimming-output/${SAMPLE}_trimmed_R2.fastq.gz > alignment-output/${SAMPLE}.sam;
done



#Convert the SAM file into a BAM file, sort and index BAM file
echo "Converting the SAM into BAM"

SAMPLES=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

for SAMPLE in "${SAMPLES[@]}"; do samtools view -hbo alignment-output/${SAMPLE}.bam alignment-output/${SAMPLE}.sam; samtools sort alignment-output/${SAMPLE}.bam -o alignment-output/${SAMPLE}.sorted.bam; samtools index alignment-output/${SAMPLE}.sorted.bam;
done






#6-Varient calling
echo "Conducting the Varient calling"
mkdir  vc-output

SAMPLES=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

for SAMPLE in "${SAMPLES[@]}"; do
bcftools mpileup -Ou -f ref/reference.fasta alignment-output/${SAMPLE}.sorted.bam | bcftools call -Ov -mv > vc-output/${SAMPLE}.vcf; done




































































