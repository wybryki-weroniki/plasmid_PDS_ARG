#!/usr/bin/env python3
"""
Split E. coli metadata into plasmid and chromosome subsets.

Given a metadata CSV filtered to E. coli (e.g., metadata_updated_ecoli.csv),
this script writes two files:
  - metadata_updated_ecoli_plasmid.csv  (rows where Type == 'Plasmid')
  - metadata_updated_ecoli_chromosome.csv  (rows where Type == 'Chromosome')
"""
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', default='metadata_updated_ecoli.csv',
                        help='Input E. coli metadata CSV file')
    parser.add_argument('-p', '--plasmid_out', default='metadata_updated_ecoli_plasmid.csv',
                        help='Output CSV for Plasmid rows')
    parser.add_argument('-c', '--chromosome_out', default='metadata_updated_ecoli_chromosome.csv',
                        help='Output CSV for Chromosome rows')
    args = parser.parse_args()

    with open(args.input, newline='') as fin:
        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames
        plasmid_rows = []
        chrom_rows = []
        for row in reader:
            t = row.get('Type', '')
            if t == 'Plasmid':
                plasmid_rows.append(row)
            elif t == 'Chromosome':
                chrom_rows.append(row)

    with open(args.plasmid_out, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(plasmid_rows)

    with open(args.chromosome_out, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(chrom_rows)

if __name__ == '__main__':
    main()