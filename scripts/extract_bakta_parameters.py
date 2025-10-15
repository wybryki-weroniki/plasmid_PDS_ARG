#!/usr/bin/env python3
"""
Script to extract useful parameters from FASTA files and metadata for Bakta API annotation.

This script extracts:
1. Circular/linear topology from FASTA headers
2. Taxonomic information from metadata
3. Replicon types (chromosome/plasmid)
4. Sequence lengths
5. Creates replicon table for Bakta
"""

import os
import pandas as pd
import re
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('extract_bakta_parameters.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def parse_fasta_header(header_line: str) -> Dict[str, str]:
    """
    Parse FASTA header to extract sequence metadata.

    Example header: ">RHB01-C04_1 1 length=4843636 depth=1.00x circular=true"

    Args:
        header_line: FASTA header line starting with '>'

    Returns:
        Dictionary with parsed metadata
    """
    metadata = {}

    # Remove '>' and split by spaces
    parts = header_line.strip('>').split()

    if len(parts) > 0:
        metadata['sequence_id'] = parts[0]

    # Extract key=value pairs
    for part in parts:
        if '=' in part:
            key, value = part.split('=', 1)
            metadata[key] = value

    return metadata


def extract_fasta_info(fasta_path: str) -> Dict[str, str]:
    """
    Extract information from a single FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary with extracted information
    """
    info = {
        'file_path': fasta_path,
        'filename': os.path.basename(fasta_path),
        'sequence_id': None,
        'length': None,
        'circular': False,
        'depth': None
    }

    try:
        with open(fasta_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('>'):
                metadata = parse_fasta_header(first_line)

                info['sequence_id'] = metadata.get('sequence_id', os.path.basename(fasta_path).replace('.fasta', ''))
                info['length'] = metadata.get('length')
                info['circular'] = metadata.get('circular', 'false').lower() == 'true'
                info['depth'] = metadata.get('depth')

        logger.info(f"Processed FASTA: {info['filename']} - Circular: {info['circular']}, Length: {info['length']}")

    except Exception as e:
        logger.error(f"Error processing FASTA {fasta_path}: {str(e)}")

    return info


def load_metadata(metadata_path: str) -> pd.DataFrame:
    """
    Load and process metadata CSV file.

    Args:
        metadata_path: Path to metadata CSV file

    Returns:
        Processed metadata DataFrame
    """
    try:
        metadata = pd.read_csv(metadata_path)
        logger.info(f"Loaded metadata with {len(metadata)} rows and columns: {list(metadata.columns)}")
        return metadata
    except Exception as e:
        logger.error(f"Error loading metadata {metadata_path}: {str(e)}")
        return pd.DataFrame()


def create_replicon_table(fasta_info_list: List[Dict], metadata_df: pd.DataFrame, output_path: str) -> None:
    """
    Create replicon table for Bakta API.

    Bakta replicon table format:
    locus,new_locus,type,topology,name

    Args:
        fasta_info_list: List of FASTA information dictionaries
        metadata_df: Metadata DataFrame
        output_path: Path to save replicon table
    """
    replicon_data = []

    for info in fasta_info_list:
        sequence_id = info['sequence_id']
        filename_base = info['filename'].replace('.fasta', '')

        # Find matching metadata row
        metadata_row = None
        if not metadata_df.empty:
            # Try matching by Contig column
            matches = metadata_df[metadata_df['Contig'] == filename_base]
            if not matches.empty:
                metadata_row = matches.iloc[0]
            else:
                # Try partial matching
                matches = metadata_df[metadata_df['Contig'].str.contains(filename_base.split('_')[0], na=False)]
                if not matches.empty:
                    metadata_row = matches.iloc[0]

        # Determine replicon type
        replicon_type = "contig"  # default
        if metadata_row is not None:
            type_col = metadata_row.get('Type', '').lower()
            if 'chromosome' in type_col:
                replicon_type = "chromosome"
            elif 'plasmid' in type_col:
                replicon_type = "plasmid"

        # Determine topology
        topology = "circular" if info['circular'] else "linear"

        # Create name
        name = sequence_id
        if metadata_row is not None:
            taxon = metadata_row.get('mlst.PubMLST', '')
            if taxon and taxon != 'nan':
                name = f"{taxon}_{sequence_id}"

        replicon_data.append({
            'locus': sequence_id,
            'new_locus': sequence_id,
            'type': replicon_type,
            'topology': topology,
            'name': name
        })

    # Save replicon table
    replicon_df = pd.DataFrame(replicon_data)
    replicon_df.to_csv(output_path, index=False)
    logger.info(f"Created replicon table with {len(replicon_data)} entries: {output_path}")

    return replicon_df


def create_bakta_parameters_json(fasta_info_list: List[Dict], metadata_df: pd.DataFrame, output_path: str) -> None:
    """
    Create JSON file with extracted parameters for Bakta API calls.

    Args:
        fasta_info_list: List of FASTA information dictionaries
        metadata_df: Metadata DataFrame
        output_path: Path to save JSON file
    """
    bakta_parameters = {}

    for info in fasta_info_list:
        filename_base = info['filename'].replace('.fasta', '')

        # Find matching metadata
        metadata_row = None
        if not metadata_df.empty:
            matches = metadata_df[metadata_df['Contig'] == filename_base]
            if not matches.empty:
                metadata_row = matches.iloc[0]

        # Extract taxon information
        genus = None
        species = None
        if metadata_row is not None:
            taxon = metadata_row.get('mlst.PubMLST', '')
            if taxon and taxon != 'nan':
                if taxon == 'ecoli':
                    genus = "Escherichia"
                    species = "coli"
                elif taxon == 'koxytoca':
                    genus = "Klebsiella"
                    species = "oxytoca"
                # Add more mappings as needed

        # Create parameter set
        params = {
            'fasta_file': info['file_path'],
            'sequence_id': info['sequence_id'],
            'length': info['length'],
            'circular': info['circular'],
            'replicon_type': metadata_row.get('Type', 'Unknown') if metadata_row is not None else 'Unknown',
            'genus': genus,
            'species': species,
            'taxon': metadata_row.get('mlst.PubMLST', '') if metadata_row is not None else '',
        }

        bakta_parameters[filename_base] = params

    # Save parameters
    with open(output_path, 'w') as f:
        json.dump(bakta_parameters, f, indent=2)

    logger.info(f"Created Bakta parameters JSON with {len(bakta_parameters)} entries: {output_path}")
    return bakta_parameters


def main():
    """Main function to extract all Bakta parameters."""
    # Set up paths
    assemblies_dir = Path(__file__).parent.parent
    metadata_path = assemblies_dir / "metadata_updated.csv"

    # Find all FASTA files
    fasta_files = list(assemblies_dir.glob("*.fasta"))
    logger.info(f"Found {len(fasta_files)} FASTA files")

    if len(fasta_files) == 0:
        logger.error("No FASTA files found in assemblies directory")
        return

    # Load metadata
    metadata_df = load_metadata(str(metadata_path))

    # Process all FASTA files
    fasta_info_list = []
    for fasta_file in fasta_files:
        info = extract_fasta_info(str(fasta_file))
        fasta_info_list.append(info)

    # Create output directory
    output_dir = Path(__file__).parent / "bakta_params"
    output_dir.mkdir(exist_ok=True)

    # Create replicon table
    replicon_table_path = output_dir / "replicons.csv"
    replicon_df = create_replicon_table(fasta_info_list, metadata_df, str(replicon_table_path))

    # Create parameters JSON
    params_json_path = output_dir / "bakta_parameters.json"
    bakta_params = create_bakta_parameters_json(fasta_info_list, metadata_df, str(params_json_path))

    # Create summary report
    summary_path = output_dir / "extraction_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("BAKTA PARAMETER EXTRACTION SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total FASTA files processed: {len(fasta_files)}\n")
        f.write(f"Metadata entries: {len(metadata_df)}\n\n")

        # Count by type
        type_counts = replicon_df['type'].value_counts()
        f.write("Replicon types:\n")
        for replicon_type, count in type_counts.items():
            f.write(f"  {replicon_type}: {count}\n")

        # Count by topology
        topo_counts = replicon_df['topology'].value_counts()
        f.write("\nTopology distribution:\n")
        for topology, count in topo_counts.items():
            f.write(f"  {topology}: {count}\n")

        # Taxon distribution
        taxon_counts = {}
        for params in bakta_params.values():
            taxon = params.get('taxon', 'Unknown')
            taxon_counts[taxon] = taxon_counts.get(taxon, 0) + 1

        f.write("\nTaxon distribution:\n")
        for taxon, count in sorted(taxon_counts.items()):
            f.write(f"  {taxon}: {count}\n")

    logger.info(f"Parameter extraction complete! Check {output_dir} for results.")
    print(f"\nSUMMARY:")
    print(f"Processed {len(fasta_files)} FASTA files")
    print(f"Created replicon table: {replicon_table_path}")
    print(f"Created parameters JSON: {params_json_path}")
    print(f"Created summary report: {summary_path}")


if __name__ == "__main__":
    main()