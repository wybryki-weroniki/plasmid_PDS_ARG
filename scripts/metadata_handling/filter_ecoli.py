#!/usr/bin/env python3
"""
Filter the metadata CSV to only rows where mlst-PubMLST == 'ecoli'.
"""
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', default='metadata_updated.csv',
                        help='Input metadata CSV file')
    parser.add_argument('-o', '--output', default='metadata_updated_ecoli.csv',
                        help='Output CSV file for filtered records')
    args = parser.parse_args()

    with open(args.input, newline='') as fin:
        reader = csv.DictReader(fin)
        rows = [r for r in reader if r.get('mlst-PubMLST') == 'ecoli']
        fieldnames = reader.fieldnames

    with open(args.output, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)

if __name__ == '__main__':
    main()