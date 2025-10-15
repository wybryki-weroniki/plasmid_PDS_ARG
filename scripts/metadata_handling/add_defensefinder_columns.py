#!/usr/bin/env python3
"""
Add defense finder summary columns to the metadata CSV.
 - systems_count: number of detected defense systems
 - defensefinder: comma-separated list of system subtypes
 - systems_gene_count: sum of genes_count across systems
"""
import csv
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', default='metadata_updated.csv',
                        help='Input metadata CSV file')
    parser.add_argument('-o', '--output', default='metadata_updated.csv',
                        help='Output CSV file to write updates')
    args = parser.parse_args()

    with open(args.input, newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames.copy()
    # Ensure new columns exist
    new_cols = ['systems_count', 'defensefinder', 'systems_gene_count']
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    for row in rows:
        contig = row.get('Contig', '')
        fname = f"{contig}_defense_finder_systems.tsv"
        if not os.path.exists(fname):
            row['systems_count'] = '0'
            row['defensefinder'] = ''
            row['systems_gene_count'] = '0'
            continue
        with open(fname, newline='') as f2:
            tsv = csv.DictReader(f2, delimiter='\t')
            subtypes = []
            gene_counts = []
            for rec in tsv:
                st = rec.get('subtype', '').strip()
                if st:
                    subtypes.append(st)
                try:
                    gene_counts.append(int(rec.get('genes_count', '0')))
                except ValueError:
                    gene_counts.append(0)
            row['systems_count'] = str(len(subtypes))
            row['defensefinder'] = ','.join(subtypes)
            row['systems_gene_count'] = str(sum(gene_counts))

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)

if __name__ == '__main__':
    main()