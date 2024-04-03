#####################################################
#!/bin/bash

# change working directory to the directory where the data is stored
#1
mkdir /home/nuero/ngs_pip
cd /home/nuero/ngs_pip
mkdir data_set_1
cd data_set_1
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget https://zenodo.org/records/10886725/files/Reference.fasta


# STEP 1: Run fastqc
fastqc ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz

mkdir quality_controlled_1
cd 
cd ngs_pip
cd data_set_1
mv ERR8774458_1_fastqc.zip ERR8774458_2_fastqc.zip quality_controlled_1 
cd quality_controlled_1
multiqc .

cd 
cd ngs_pip
cd data_set_1
# Step 2: trimming with fastp
fastp -i  ERR8774458_1.fastq.gz -I ERR8774458_2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq

mkdir trimmed
mv trimmed_R1.fastq trimmed_R2.fastq trimmed
cd trimmed 
fastqc trimmed_R1.fastq trimmed_R2.fastq
multiqc .
cd 
cd ngs_pip
cd data_set_1

# Step 2: Mapping with BWA
cp Reference.fasta trimmed
cd trimmed

#indexing the reference file
bwa index -a bwtsw Reference.fasta
bwa mem Reference.fasta trimmed_R1.fastq trimmed_R2.fastq > output_data_set_1.sam

#converting to bam file
samtools view -hbo  output_data_set_1.bam  output_data_set_1.sam
#sorting the bam file
samtools sort  output_data_set_1.bam -o  output_data_sorted_set_1.bam
#Indexing the sorted bam file
samtools index  output_data_sorted_set_1.bam
#Variant calling
bcftools mpileup -Ou -f Reference.fasta output_data_sorted_set_1.bam | bcftools  call -Ov -mv > variant_data_set_1.vcf  




#####################################################
#Data set 2

cd /home/nuero/ngs_pip
mkdir data_set_2
cd data_set_2

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta


# STEP 1: Run fastqc
fastqc ACBarrie_R1.fastq.gz ACBarrie_R2.fastq.gz

mkdir quality_controlled_2
cd 
cd ngs_pip
cd data_set_2
mv ACBarrie_R1_fastqc.zip ACBarrie_R2_fastqc.zip quality_controlled_2
cd quality_controlled_2
multiqc .

cd 
cd ngs_pip
cd data_set_2


# Step 2: trimming with fastp
fastp -i ACBarrie_R1.fastq.gz -I ACBarrie_R2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq

mkdir trimmed
mv trimmed_R1.fastq trimmed_R2.fastq trimmed
cd trimmed
fastqc trimmed_R1.fastq trimmed_R2.fastq
multiqc .
cd

cd ngs_pip
cd data_set_2
# Step 2: Mapping with BWA
cp reference.fasta trimmed
cd trimmed


#indexing the reference file
bwa index -a bwtsw reference.fasta
bwa mem reference.fasta trimmed_R1.fastq trimmed_R2.fastq > output_data_set_2.sam

#converting to bam file
samtools view -hbo  output_data_set_2.bam  output_data_set_2.sam
#sorting the bam file
samtools sort  output_data_set_2.bam -o  output_data_sorted_set_2.bam
#Indexing the sorted bam file
samtools index  output_data_sorted_set_2.bam
#Variant calling
bcftools mpileup -Ou -f reference.fasta output_data_sorted_set_2.bam | bcftools  call -Ov -mv > variant_data_set_2.vcf



#####################################################
#Data set 3
cd /home/nuero/ngs_pip
mkdir data_set_3
cd data_set_3

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta 


# STEP 1: Run fastqc
fastqc Alsen_R1.fastq.gz Alsen_R2.fastq.gz



mkdir quality_controlled_3
cd 
cd ngs_pip
cd data_set_3
mv  Alsen_R1_fastqc.zip Alsen_R2_fastqc.zip quality_controlled_3
cd quality_controlled_3
multiqc .

cd 
cd ngs_pip
cd data_set_3

# Step 2: trimming with fastp
fastp -i Alsen_R1.fastq.gz  -I Alsen_R2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq

mkdir trimmed
mv trimmed_R1.fastq trimmed_R2.fastq trimmed
cd trimmed
fastqc trimmed_R1.fastq trimmed_R2.fastq
multiqc .
cd

cd ngs_pip
cd data_set_3
# Step 2: Mapping with BWA
cp reference.fasta trimmed
cd trimmed

#indexing the reference file
bwa index -a bwtsw reference.fasta
bwa mem reference.fasta trimmed_R1.fastq trimmed_R2.fastq > output_data_set_3.sam

#converting to bam file
samtools view -hbo  output_data_set_3.bam  output_data_set_3.sam
#sorting the bam file
samtools sort  output_data_set_3.bam -o  output_data_sorted_set_3.bam
#Indexing the sorted bam file
samtools index  output_data_sorted_set_3.bam
#Variant calling
bcftools mpileup -Ou -f reference.fasta output_data_sorted_set_3.bam | bcftools  call -Ov -mv > variant_data_set_3.vcf


#####################################################
#Data set 4
cd /home/nuero/ngs_pip
mkdir data_set_4
cd data_set_4

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta  

# STEP 1: Run fastqc
fastqc Baxter_R1.fastq.gz Baxter_R2.fastq.gz

mkdir quality_controlled_4
cd 
cd ngs_pip
cd data_set_4
mv Baxter_R1_fastqc.zip Baxter_R2_fastqc.zip quality_controlled_4
cd quality_controlled_4
multiqc .

cd 
cd ngs_pip
cd data_set_4

# Step 2: trimming with fastp
fastp -i Baxter_R1.fastq.gz  -I  Baxter_R2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq

mkdir trimmed
mv trimmed_R1.fastq trimmed_R2.fastq trimmed
cd trimmed
fastqc trimmed_R1.fastq trimmed_R2.fastq
multiqc .
cd

cd ngs_pip
cd data_set_4
# Step 2: Mapping with BWA
cp reference.fasta trimmed
cd trimmed

#indexing the reference file
bwa index -a bwtsw reference.fasta
bwa mem reference.fasta trimmed_R1.fastq trimmed_R2.fastq > output_data_set_4.sam

#converting to bam file
samtools view -hbo  output_data_set_4.bam  output_data_set_4.sam
#sorting the bam file
samtools sort  output_data_set_4.bam -o  output_data_sorted_set_4.bam
#Indexing the sorted bam file
samtools index  output_data_sorted_set_4.bam
#Variant calling
bcftools mpileup -Ou -f reference.fasta output_data_sorted_set_4.bam | bcftools  call -Ov -mv > variant_data_set_4.vcf


#####################################################
#Data set 5
cd /home/nuero/ngs_pip
mkdir data_set_5
cd data_set_5

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

# STEP 1: Run fastqc
fastqc  Chara_R1.fastq.gz Chara_R2.fastq.gz

mkdir quality_controlled_5
cd
cd ngs_pip
cd data_set_5
mv  Chara_R1_fastqc.zip Chara_R2_fastqc.zip  quality_controlled_5
cd quality_controlled_5
multiqc .


cd 
cd ngs_pip
cd data_set_5
# Step 2: trimming with fastp
fastp -i   Chara_R1.fastq.gz -I Chara_R2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq

mkdir trimmed
mv trimmed_R1.fastq trimmed_R2.fastq trimmed
cd trimmed
fastqc trimmed_R1.fastq trimmed_R2.fastq
multiqc .
cd

cd ngs_pip
cd data_set_5
# Step 2: Mapping with BWA
cp reference.fasta trimmed
cd trimmed

#indexing the reference file
bwa index -a bwtsw reference.fasta
bwa mem reference.fasta trimmed_R1.fastq trimmed_R2.fastq > output_data_set_5.sam

#converting to bam file
samtools view -hbo  output_data_set_5.bam  output_data_set_5.sam
#sorting the bam file
samtools sort  output_data_set_5.bam -o  output_data_sorted_set_5.bam
#Indexing the sorted bam file
samtools index  output_data_sorted_set_5.bam
#Variant calling
bcftools mpileup -Ou -f reference.fasta output_data_sorted_set_5.bam | bcftools  call -Ov -mv > variant_data_set_5.vcf

#####################################################
#Data set 6
cd /home/nuero/ngs_pip
mkdir data_set_6
cd data_set_6

wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

# STEP 1: Run fastqc
fastqc  Drysdale_R1.fastq.gz Drysdale_R2.fastq.gz

mkdir quality_controlled_6
cd
cd ngs_pip
cd data_set_6
mv   Drysdale_R1_fastqc.zip Drysdale_R2_fastqc.zip  quality_controlled_6
cd quality_controlled_6
multiqc .

cd 
cd ngs_pip
cd data_set_6
# Step 2: trimming with fastp
fastp -i  Drysdale_R1.fastq.gz -I Drysdale_R2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq

mkdir trimmed
mv trimmed_R1.fastq trimmed_R2.fastq trimmed
cd trimmed

fastqc trimmed_R1.fastq trimmed_R2.fastq
multiqc .
cd


cd ngs_pip
cd data_set_6
# Step 2: Mapping with BWA
cp reference.fasta trimmed
cd trimmed

#indexing the reference file
bwa index -a bwtsw reference.fasta
bwa mem reference.fasta trimmed_R1.fastq trimmed_R2.fastq > output_data_set_6.sam

#converting to bam file
samtools view -hbo  output_data_set_6.bam  output_data_set_6.sam
#sorting the bam file
samtools sort  output_data_set_6.bam -o  output_data_sorted_set_6.bam
#Indexing the sorted bam file
samtools index  output_data_sorted_set_6.bam
#Variant calling
bcftools mpileup -Ou -f reference.fasta output_data_sorted_set_6.bam | bcftools  call -Ov -mv > variant_data_set_6.vcf


#####################################################