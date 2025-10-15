#!/usr/bin/env python3
"""
Add a 'system_length' column to the detailed metadata CSV by summing the amino acid
lengths of all proteins listed in 'protein_in_syst' for each system.
Proteins are read from Prodigal .prt FASTA output files in the repository.
"""
import csv
import os
import sys
import argparse


def parse_prt_lengths(path, wanted, lengths):
    """Parse a Prodigal .prt FASTA file and record lengths for wanted IDs."""
    with open(path) as fh:
        seq_id = None
        seq_len = 0
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if seq_id and seq_id in wanted:
                    lengths[seq_id] = seq_len
                seq_id = line[1:].split()[0]
                seq_len = 0
            else:
                if seq_id in wanted:
                    seq_len += len(line.strip().replace('*', ''))
        if seq_id and seq_id in wanted:
            lengths[seq_id] = seq_len


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', default='metadata_updated_ecoli_systems_detail.csv',
                        help='Input detailed metadata CSV file')
    parser.add_argument('-o', '--output', default='metadata_updated_ecoli_systems_detail.csv',
                        help='Output CSV file with system_length added')
    args = parser.parse_args()

    with open(args.input, newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames.copy()

    col_name = 'system_length'
    if col_name not in fieldnames:
        fieldnames.append(col_name)

    wanted = set()
    for row in rows:
        prot_str = row.get('protein_in_syst', '').strip()
        if not prot_str or prot_str.upper() == 'NA':
            continue
        for pid in prot_str.split(','):
            if pid:
                wanted.add(pid)

    lengths = {}
    for root, dirs, files in os.walk('.'):
        for fname in files:
            if fname.endswith('.prt'):
                parse_prt_lengths(os.path.join(root, fname), wanted, lengths)

    missing = set()
    for row in rows:
        prot_str = row.get('protein_in_syst', '').strip()
        if not prot_str or prot_str.upper() == 'NA':
            row[col_name] = 'NA'
        else:
            total = 0
            for pid in prot_str.split(','):
                if pid in lengths:
                    total += lengths[pid]
                else:
                    missing.add(pid)
            row[col_name] = str(total)

    if missing:
        sys.stderr.write(f"Warning: {len(missing)} protein IDs not found in .prt files:\n")
        for pid in sorted(missing):
            sys.stderr.write(f"  {pid}\n")

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()