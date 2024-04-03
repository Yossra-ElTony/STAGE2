#!/bin/bash

# Install necessary tools using conda
conda install -y -c bioconda bwa samtools bcftools fastqc multiqc

# Install Sickle using apt-get
sudo apt-get update
sudo apt-get install -y sickle

echo "Setup completed."

