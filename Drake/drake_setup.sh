#!/bin/bash

# Check if Conda is installed
if ! command -v conda &> /dev/null
then
    echo "Conda not found. Installing Miniconda..."
    # Download and install Miniconda 
    curl -o ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash ~/miniconda.sh -b -p $HOME/miniconda
    rm ~/miniconda.sh
    # Activate Conda (might need to change based on your shell, here we assume bash)
    source $HOME/miniconda/bin/activate
    echo 'source $HOME/miniconda/bin/activate' >> ~/.bashrc  # Add to .bashrc for persistence
fi

# Update Conda and its repositories
conda update -n base -c defaults conda

# Install necessary tools using Conda
conda install -y bcftools bgzip bwa conda fastqc freebayes multiqc samtools

echo "Setup complete."
