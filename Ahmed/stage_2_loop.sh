# Function to analyze a dataset
analyze_dataset() {
    local dataset_name="$1"
    local url_R1="$2"
    local url_R2="$3"
    local reference_url="$4"
    
    echo "Analyzing dataset: $dataset_name"

    # Create directory for the dataset
    mkdir -p "$dataset_name" && cd "$dataset_name"

    # Download reference genome
    wget -O Reference.fasta "$reference_url"

    # Download dataset files
    wget -O "${dataset_name}_1.fastq.gz" "$url_R1"
    wget -O "${dataset_name}_2.fastq.gz" "$url_R2"

    # Quality control
    mkdir -p quality_control && \
    for f in *.fastq.gz; do fastqc -t 1 -f fastq -noextract "$f" -o quality_control; done

    multiqc -z -o quality_control .

    # Trimming 
    sickle pe -f "${dataset_name}_1.fastq.gz" -r "${dataset_name}_2.fastq.gz" -t sanger -q 20 -g -o "${dataset_name}_1.trimmed.fastq.gz" -p "${dataset_name}_2.trimmed.fastq.gz" -s "${dataset_name}_trimms.fastq.gz"

    # Reference genome indexing
    bwa index -a bwtsw Reference.fasta

    # Generating SAM file
    bwa mem -t 8 Reference.fasta "${dataset_name}_1.trimmed.fastq.gz" "${dataset_name}_2.trimmed.fastq.gz" > "${dataset_name}.sam"

    # Converting SAM into BAM
    samtools view -S -b "${dataset_name}.sam" > "${dataset_name}.bam"

    # Checking size difference
    du -sh "${dataset_name}.sam"
    du -sh "${dataset_name}.bam"

    # Sorting the created BAM file
    samtools sort -o "${dataset_name}_sorted.bam" "${dataset_name}.bam"

    # Indexing Bam file
    samtools index "${dataset_name}_sorted.bam"

    # Mapping statistics
    samtools flagstat "${dataset_name}_sorted.bam" > "${dataset_name}_mappingstats.txt"
    cat "${dataset_name}_mappingstats.txt"

    # Generate coverage information
    bcftools mpileup -O b -o "${dataset_name}.bcf" -f Reference.fasta --threads 8 "${dataset_name}_sorted.bam"

    # Variant calling
    bcftools mpileup -Ou -f Reference.fasta "${dataset_name}_sorted.bam" | bcftools call -Ov -mv > "${dataset_name}_variants.vcf"

    cd ..
}

# Define the datasets
datasets=(
    "ERR8774458 https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1 https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1 https://zenodo.org/records/10886725/files/Reference.fasta?download=1"
    "ACBarrie https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"
    "Alsen https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"
    "Baxter https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"
    "Chara https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"
    "Drysdale https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"
)

# Iterate over the datasets and analyze each one
for dataset_info in "${datasets[@]}"; do
    read -r -a dataset <<< "$dataset_info"
    analyze_dataset "${dataset[@]}"
done

