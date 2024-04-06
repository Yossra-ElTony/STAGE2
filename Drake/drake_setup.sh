#!/bin/bash

# Update Conda and its repositories
conda update -n base -c defaults conda

# Install necessary tools using Conda
conda install -y bcftools bgzip bwa fastqc freebayes multiqc samtools

echo "Setup complete."
