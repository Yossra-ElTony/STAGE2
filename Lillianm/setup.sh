#!/bin/bash

# Function to check if a package is installed
package_installed() {
  dpkg -l "$1" &> /dev/null
}

# Function to install a package if it's not already installed
install_package() {
  local package_name="$1"
  if package_installed "$package_name"; then
    echo "$package_name is already installed."
  else
    sudo apt-get install -y "$package_name"
  fi
}

# Install FastQC
install_package "fastqc"

# Install Fastp
install_package "fastp"

# Install BWA
install_package "bwa"

# Install Samtools
install_package "samtools"

# Install Bcftools
install_package "bcftools"

echo "Installation complete."

# Make this script executable
chmod +x "$0"

