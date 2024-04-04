#!/bin/bash

# Downloading the first file
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1

# Renaming the downloaded files
mv 'ERR8774458_1.fastq.gz?download=1' forwardstr.fastq.gz
mv 'ERR8774458_2.fastq.gz?download=1' reversestr.fastq.gz

# Downloading the first reference genome
mkdir -p ref
cd ref
wget https://zenodo.org/records/10886725/files/Reference.fasta?download=1 -O Reference.fasta
cd ..

# Processing the first file
for f in *.fastq.gz; do fastqc -t 1 -f fastq -noextract "$f"; done

mkdir -p qc_output
mv *.html qc_output
mv *.zip qc_output

cd qc_output
multiqc -z -o . .

cd ..

# Trimming
fastp -i forwardstr.fastq.gz -I reversestr.fastq.gz -o trimmed_forward.fastq.gz -O trimmed_reverse.fastq.gz --detect_adapter_for_pe

# Indexing the first reference genome
bwa index ref/Reference.fasta

# Processing the first file
bwa mem ref/Reference.fasta trimmed_forward.fastq.gz trimmed_reverse.fastq.gz > aligned_reads.sam
samtools view -bS aligned_reads.sam > aligned_reads.bam
samtools sort aligned_reads.bam -o sorted_aligned_reads.bam
samtools index sorted_aligned_reads.bam
bcftools mpileup -Ou -f ref/Reference.fasta sorted_aligned_reads.bam | bcftools call -mv -Ov -o variants.vcf

# Downloading and indexing the second reference genome
mkdir -p ref2
cd ref2
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta
bwa index reference.fasta
cd ..

# Now processing the other five  files
declare -a files=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

for file in "${files[@]}"; do
    # Downloading the files
    wget "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/${file}_R1.fastq.gz"
    wget "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/${file}_R2.fastq.gz"

    # Trimming
    fastp -i "${file}_R1.fastq.gz" -I "${file}_R2.fastq.gz" -o "${file}_trimmed_R1.fastq.gz" -O "${file}_trimmed_R2.fastq.gz" --detect_adapter_for_pe

    # Aligning and variant calling using the second reference genome
    bwa mem ref2/reference.fasta "${file}_trimmed_R1.fastq.gz" "${file}_trimmed_R2.fastq.gz" > "${file}_aligned_reads.sam"
    samtools view -bS "${file}_aligned_reads.sam" > "${file}_aligned_reads.bam"
    samtools sort "${file}_aligned_reads.bam" -o "${file}_sorted_aligned_reads.bam"
    samtools index "${file}_sorted_aligned_reads.bam"
    bcftools mpileup -Ou -f ref2/reference.fasta "${file}_sorted_aligned_reads.bam" | bcftools call -mv -Ov -o "${file}_variants.vcf"
done
echo "All processing completed."

