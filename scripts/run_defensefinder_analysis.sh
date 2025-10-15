#!/bin/bash
"""
Defense Finder Analysis Script
==============================

This script runs DefenseFinder on all FASTA files to identify bacterial defense systems.
Generates {contig}_defense_finder_systems.tsv files for downstream analysis.

Requirements:
- Conda/Miniconda installed
- FASTA files in the current directory
- Internet connection for package installation and database updates

Date: June 2025
"""

set -e  # Exit on any error

# Configuration
ENV_NAME="defensefinder_updated"
PYTHON_VERSION="3.10"

echo "DefenseFinder Analysis Script"
echo "============================="
echo ""

# Step 1: Create and activate conda environment
echo "Setting up conda environment: $ENV_NAME"
conda create -n "$ENV_NAME" python="$PYTHON_VERSION" click=8.0.3 -c conda-forge -y
source ~/miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

# Step 2: Install DefenseFinder
echo "Installing DefenseFinder..."
pip install mdmparis-defense-finder

# Step 3: Update DefenseFinder database
echo "Updating DefenseFinder database..."
defense-finder update

# Step 4: Run analysis on all FASTA files
echo "Running DefenseFinder analysis on all FASTA files..."
echo "This may take several hours depending on dataset size..."

for i in *.fasta; do
    echo "Processing: $i"
    defense-finder run "$i"
done

# Step 5: Generate summary
echo "Generating results summary..."
wc -l *_defense_finder_systems.tsv > total-system-number.tsv

echo ""
echo "Analysis completed successfully!"
echo "Results:"
echo "- Individual TSV files: {contig}_defense_finder_systems.tsv"
echo "- Summary counts: total-system-number.tsv"
echo ""