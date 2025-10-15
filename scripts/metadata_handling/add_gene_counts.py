#!/usr/bin/env python3
"""
Annotate metadata CSV files with total gene counts from Bakta JSON outputs.

This script:
  - Scans all bakta result JSON files in scripts/bakta_results/,
    counts CDS features per file, and maps contig name to gene count.
  - For each metadata CSV matching 'metadata*.csv', adds (or updates) a
    'gene_count' column with the count for the Contig value.

Usage:
  python3 scripts/add_gene_counts.py
  (Modifies all metadata*.csv files in place)
"""
import csv
import glob
import json
import os

def main():
    # Build mapping from contig to gene count from bakta JSON results
    gene_counts = {}
    bakta_results_dir = 'scripts/bakta_results'

    # Handle both relative path (from root directory) and absolute path
    if not os.path.exists(bakta_results_dir):
        bakta_results_dir = 'bakta_results'  # If running from scripts/ directory

    if os.path.exists(bakta_results_dir):
        json_files = glob.glob(os.path.join(bakta_results_dir, '*_bakta_results.json'))
        total_files = len(json_files)
        print(f"Processing {total_files} bakta result files...")

        for i, json_file in enumerate(json_files, 1):
            # Extract contig name from filename (remove _bakta_results.json)
            basename = os.path.basename(json_file)
            contig = basename.replace('_bakta_results.json', '')

            count = 0
            try:
                with open(json_file, 'r') as jf:
                    data = json.load(jf)
                    # Count CDS features in the bakta results
                    for feature in data.get('features', []):
                        if feature.get('type') == 'cds':
                            count += 1
            except (json.JSONDecodeError, IOError, KeyError):
                continue
            gene_counts[contig] = count

            # Progress reporting
            if i % 100 == 0 or i == total_files:
                print(f"  Processed {i}/{total_files} bakta files...")

        print(f"Found gene counts for {len(gene_counts)} contigs")
    else:
        print(f"Warning: Bakta results directory not found at '{bakta_results_dir}'")
        print("No gene counts will be added.")

    # Annotate each metadata CSV
    metadata_files = glob.glob('metadata*.csv')
    print(f"\nUpdating {len(metadata_files)} metadata files...")

    for meta_file in metadata_files:
        print(f"  Processing {meta_file}...")
        # Read existing data
        with open(meta_file, newline='') as mf:
            reader = csv.DictReader(mf)
            rows = list(reader)
            fieldnames = reader.fieldnames.copy()

        # Ensure gene_count column exists
        if 'gene_count' not in fieldnames:
            fieldnames.append('gene_count')

        # Update rows
        updated_count = 0
        for row in rows:
            contig = row.get('Contig', '')
            gene_count = gene_counts.get(contig, 0)
            row['gene_count'] = str(gene_count)
            if gene_count > 0:
                updated_count += 1

        # Write back in place
        with open(meta_file, 'w', newline='') as mf:
            writer = csv.DictWriter(mf, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
            writer.writeheader()
            writer.writerows(rows)

        print(f"    Updated {updated_count}/{len(rows)} rows with gene counts")

    print(f"\nCompleted! Updated gene counts in {len(metadata_files)} metadata files.")

if __name__ == '__main__':
    main()