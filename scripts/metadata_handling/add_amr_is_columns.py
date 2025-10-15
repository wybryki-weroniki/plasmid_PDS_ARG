#!/usr/bin/env python3
"""
Add AMR and IS summary columns to the metadata CSV.
 - AMR_count: count of AMR genes from AMRFinderPlus
 - IS_count: count of insertion sequences from Abricate-ISfinder
 - IS_density: IS_count per 10kb of genome length
 - AMR_binary_presence: 1 if AMR_count > 0 else 0
 - AMR_presence: 'present' or 'absent'
 - systems_binary_presence: 1 if systems_count > 0 else 0
 - systems_presence: 'present' or 'absent'
"""
import csv
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

    # New columns to add
    new_cols = [
        'AMR_count', 'IS_count', 'IS_density',
        'AMR_binary_presence', 'AMR_presence',
        'systems_binary_presence', 'systems_presence'
    ]
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    for row in rows:
        # AMR_count
        amr = row.get('AMRFinderPlus', '').strip()
        if not amr or amr.upper() == 'NA':
            amr_count = 0
        else:
            amr_count = sum(1 for x in amr.split(',') if x)
        row['AMR_count'] = str(amr_count)
        # AMR binary and presence
        row['AMR_binary_presence'] = '1' if amr_count > 0 else '0'
        row['AMR_presence'] = 'present' if amr_count > 0 else 'absent'

        # IS_count
        isf = row.get('Abricate-ISfinder', '').strip()
        if not isf or isf.upper() == 'NA':
            is_count = 0
        else:
            is_count = sum(1 for x in isf.split(',') if x)
        row['IS_count'] = str(is_count)

        # IS_density per 10kb
        try:
            length = float(row.get('Length-(bp)', 0))
            density = is_count * 10000.0 / length if length > 0 else 0.0
        except ValueError:
            density = 0.0
        # format with six decimal places
        row['IS_density'] = f"{density:.6f}"

        # systems presence re-check
        try:
            sys_count = int(row.get('systems_count', '0'))
        except ValueError:
            sys_count = 0
        row['systems_binary_presence'] = '1' if sys_count > 0 else '0'
        row['systems_presence'] = 'present' if sys_count > 0 else 'absent'

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)

if __name__ == '__main__':
    main()