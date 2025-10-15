#!/usr/bin/env python3
"""
Generate a detailed metadata CSV by merging metadata_updated_ecoli.csv with defense finder systems results.
"""
import os
import pandas as pd


def main():
    meta = pd.read_csv("metadata_updated_ecoli.csv", dtype=str, keep_default_na=False)
    tsv_files = {}
    for root, dirs, files in os.walk('.'):
        for fname in files:
            if fname.endswith('_defense_finder_systems.tsv'):
                contig = fname[: -len('_defense_finder_systems.tsv')]
                tsv_files[contig] = os.path.join(root, fname)
    new_cols = [
        'type', 'subtype', 'activity', 'sys_beg', 'sys_end',
        'protein_in_syst', 'genes_count_per_system', 'name_of_profiles_in_sys'
    ]
    combined_records = []

    for _, meta_row in meta.iterrows():
        contig = meta_row['Contig']
        systems_path = tsv_files.get(contig)
        if systems_path:
            df_sys = pd.read_csv(systems_path, sep='\t', dtype=str, keep_default_na=False)
            if 'genes_count' in df_sys.columns:
                df_sys = df_sys.rename(columns={'genes_count': 'genes_count_per_system'})
            df_sys = df_sys[new_cols]
            if not df_sys.empty:
                for _, sys_row in df_sys.iterrows():
                    rec = meta_row.to_dict()
                    rec.update(sys_row.to_dict())
                    combined_records.append(rec)
                continue
        rec = meta_row.to_dict()
        for col in new_cols:
            rec[col] = 'NA'
        combined_records.append(rec)
    final_df = pd.DataFrame(combined_records, columns=list(meta.columns) + new_cols)
    final_df.to_csv("metadata_updated_ecoli_systems_detail.csv", index=False)


if __name__ == '__main__':
    main()