#!/usr/bin/env python3
"""
Run AMRFinderPlus directly on FASTA files to get detailed Element Type and Subtype information.

This script:
1. Finds all FASTA files in the dataset
2. Runs AMRFinderPlus directly on each FASTA file
3. Collects all results into a comprehensive CSV with detailed AMR information
4. Provides Element Type and Subtype columns for each AMR finding
"""
import subprocess
import pandas as pd
import argparse
from pathlib import Path
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import os

def run_amrfinder_on_fasta(fasta_file, output_dir, organism="Escherichia"):
    """Run AMRFinderPlus on a single FASTA file."""
    fasta_path = Path(fasta_file)
    contig_name = fasta_path.stem
    output_file = output_dir / f"{contig_name}_amrfinder.tsv"

    # Adaptive timeout based on file size and type
    file_size_mb = fasta_path.stat().st_size / (1024 * 1024)

    if 'chromo' in contig_name:
        # Chromosome files need more time
        timeout = max(900, int(file_size_mb * 60))  # 15 min minimum, +1 min per MB
    else:
        # Plasmid files are usually smaller
        timeout = max(600, int(file_size_mb * 30))  # 10 min minimum, +30s per MB

    timeout = min(timeout, 1800)  # Cap at 30 minutes

    try:
        # Activate conda environment and run AMRFinderPlus
        cmd = [
            'bash', '-c',
            f'source ~/miniconda3/etc/profile.d/conda.sh && conda activate amrfinder && '
            f'amrfinder --nucleotide "{fasta_file}" --organism {organism} --output "{output_file}"'
        ]

        result = subprocess.run(cmd,
                              capture_output=True,
                              text=True,
                              timeout=timeout)

        if result.returncode == 0:
            # Check if output file has results (more than just header)
            if output_file.exists():
                with open(output_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:  # More than just header
                        return {'contig': contig_name, 'status': 'success', 'results_file': output_file, 'amr_count': len(lines) - 1}
                    else:
                        return {'contig': contig_name, 'status': 'no_amr', 'results_file': None, 'amr_count': 0}
            else:
                return {'contig': contig_name, 'status': 'no_output', 'results_file': None, 'amr_count': 0}
        else:
            print(f"Error processing {contig_name}: {result.stderr}")
            return {'contig': contig_name, 'status': 'error', 'results_file': None, 'amr_count': 0}

    except subprocess.TimeoutExpired:
        print(f"Timeout processing {contig_name} (timeout: {timeout//60}min {timeout%60}s)")
        return {'contig': contig_name, 'status': 'timeout', 'results_file': None, 'amr_count': 0, 'timeout': timeout}
    except Exception as e:
        print(f"Exception processing {contig_name}: {e}")
        return {'contig': contig_name, 'status': 'exception', 'results_file': None, 'amr_count': 0}

def collect_fasta_files(assemblies_dir, filter_metadata_file=None):
    """Collect FASTA files from the assemblies directory, optionally filtered by metadata."""
    assemblies_path = Path(assemblies_dir)
    all_fasta_files = list(assemblies_path.glob('*.fasta'))
    print(f"Found {len(all_fasta_files)} total FASTA files")

    if filter_metadata_file:
        print(f"Filtering FASTA files based on metadata: {filter_metadata_file}")

        # Read the metadata file to get the Contig column
        try:
            metadata_df = pd.read_csv(filter_metadata_file)
            if 'Contig' not in metadata_df.columns:
                print(f"Warning: 'Contig' column not found in {filter_metadata_file}")
                print(f"Available columns: {list(metadata_df.columns)}")
                return all_fasta_files

            # Get the set of contig names from metadata
            contig_names = set(metadata_df['Contig'].astype(str))
            print(f"Found {len(contig_names)} contigs in metadata file")

            # Filter FASTA files to only include those with matching names
            filtered_fasta_files = []
            for fasta_file in all_fasta_files:
                contig_name = fasta_file.stem  # filename without extension
                if contig_name in contig_names:
                    filtered_fasta_files.append(fasta_file)

            print(f"Filtered to {len(filtered_fasta_files)} FASTA files matching metadata")
            print(f"Excluded: {len(all_fasta_files) - len(filtered_fasta_files)} files not in metadata")

            return filtered_fasta_files

        except Exception as e:
            print(f"Error reading metadata file {filter_metadata_file}: {e}")
            print(f"Using all FASTA files instead")
            return all_fasta_files

    return all_fasta_files

def combine_amr_results(output_dir, combined_output_file):
    """Combine all individual AMR results into a single CSV."""
    all_results = []

    # AMRFinderPlus output columns
    expected_columns = [
        'Protein id', 'Contig id', 'Start', 'Stop', 'Strand',
        'Element symbol', 'Element name', 'Scope', 'Type', 'Subtype',
        'Class', 'Subclass', 'Method', 'Target length',
        'Reference sequence length', '% Coverage of reference',
        '% Identity to reference', 'Alignment length',
        'Closest reference accession', 'Closest reference name',
        'HMM accession', 'HMM description'
    ]

    result_files = list(Path(output_dir).glob('*_amrfinder.tsv'))

    for result_file in result_files:
        contig_name = result_file.stem.replace('_amrfinder', '')

        try:
            # Read the TSV file
            df = pd.read_csv(result_file, sep='\t')

            if len(df) > 0:
                # Add contig name as a column for reference
                df['Source_Contig'] = contig_name
                all_results.append(df)

        except Exception as e:
            print(f"Error reading {result_file}: {e}")

    if all_results:
        # Combine all results
        combined_df = pd.concat(all_results, ignore_index=True)

        # Reorder columns to put Source_Contig first
        cols = ['Source_Contig'] + [col for col in combined_df.columns if col != 'Source_Contig']
        combined_df = combined_df[cols]

        # Save to CSV
        combined_df.to_csv(combined_output_file, index=False, quoting=csv.QUOTE_ALL)

        print(f"\n COMBINED RESULTS SUMMARY:")
        print(f"  • Total AMR findings: {len(combined_df)}")
        print(f"  • Contigs with AMR: {combined_df['Source_Contig'].nunique()}")
        print(f"  • Results saved to: {combined_output_file}")

        # Show some statistics
        if len(combined_df) > 0:
            print(f"\n AMR STATISTICS:")
            print(f"  • Element Types found: {combined_df['Type'].value_counts().to_dict()}")
            print(f"  • Classes found: {combined_df['Class'].value_counts().head().to_dict()}")

        return combined_df
    else:
        print("No AMR results found in any files")
        return None

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assemblies-dir',
                       default='../assemblies',
                       help='Directory containing FASTA files (default: ../assemblies)')
    parser.add_argument('--output-dir',
                       default='amrfinder_results',
                       help='Directory for individual AMRFinderPlus result files (default: amrfinder_results)')
    parser.add_argument('--combined-output',
                       default='amrfinder_combined_results.csv',
                       help='Combined output CSV file (default: amrfinder_combined_results.csv)')
    parser.add_argument('--max-workers', type=int, default=4,
                       help='Maximum number of parallel workers')
    parser.add_argument('--organism', default='Escherichia',
                       help='Organism for AMRFinderPlus (default: Escherichia)')
    parser.add_argument('--test-run', action='store_true',
                       help='Test run on first 10 files only')
    parser.add_argument('--filter-metadata',
                       help='Filter FASTA files to only process those listed in the Contig column of this metadata CSV file')

    args = parser.parse_args()

    print("AMRFinderPlus Direct Runner")
    print("=" * 50)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Collect FASTA files
    fasta_files = collect_fasta_files(args.assemblies_dir, args.filter_metadata)

    if args.test_run:
        fasta_files = fasta_files[:10]
        print(f"Test run mode: processing first {len(fasta_files)} files")

    if not fasta_files:
        print("No FASTA files found!")
        return

    print(f"Output directory: {output_dir}")
    print(f"Combined output: {args.combined_output}")
    print(f"Organism: {args.organism}")
    print(f"Max workers: {args.max_workers}")

    # Process files in parallel
    results = []
    successful = 0
    failed = 0
    no_amr = 0

    print(f"\n Processing {len(fasta_files)} FASTA files...")

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        # Submit all jobs
        future_to_file = {
            executor.submit(run_amrfinder_on_fasta, fasta_file, output_dir, args.organism): fasta_file
            for fasta_file in fasta_files
        }

        # Collect results as they complete
        for i, future in enumerate(as_completed(future_to_file), 1):
            result = future.result()
            results.append(result)

            if result['status'] == 'success':
                successful += 1
                print(f"  {result['contig']} ({result['amr_count']} AMR genes)")
            elif result['status'] == 'no_amr':
                no_amr += 1
                if i % 50 == 0:  # Less verbose for no AMR cases
                    print(f"  ... processed {i}/{len(fasta_files)} files")
            else:
                failed += 1
                print(f"  {result['contig']} ({result['status']})")

            # Progress update every 100 files
            if i % 100 == 0:
                print(f"  Progress: {i}/{len(fasta_files)} ({i/len(fasta_files)*100:.1f}%)")

    print(f"\n PROCESSING SUMMARY:")
    print(f"  • Total files: {len(fasta_files)}")
    print(f"  • Successful with AMR: {successful}")
    print(f"  • No AMR found: {no_amr}")
    print(f"  • Failed/Error: {failed}")

    # Combine all results
    print(f"\n Combining results...")
    combined_df = combine_amr_results(output_dir, args.combined_output)

    print(f"\n Complete! Direct AMRFinderPlus analysis finished.")
    print(f"Individual results in: {output_dir}")
    print(f"Combined results in: {args.combined_output}")

if __name__ == "__main__":
    main()