#!/usr/bin/env python3

"""
Script to add GTEx expression data as attributes to GENCODE GFF3 annotations.
The script matches genes based on Ensembl IDs (stripped of version numbers) and adds
expression data for selected tissues as GFF attributes in the form gtex_{tissue}: {value}.
"""

import argparse
import gzip
import csv
import re
from typing import Dict, Tuple, List, Set

def normalize_tissue_name(tissue: str) -> str:
    """
    Normalize tissue name by:
    1. Converting to lowercase
    2. Replacing non-letter chars with underscores
    3. Collapsing multiple sequential underscores
    """
    # Convert to lowercase
    tissue = tissue.lower()
    # Replace non-letter chars with underscore
    tissue = re.sub(r'[^a-z]+', '_', tissue)
    # Collapse multiple sequential underscores
    tissue = re.sub(r'_+', '_', tissue)
    # Remove leading/trailing underscores
    tissue = tissue.strip('_')
    return tissue

def strip_version(ensembl_id: str) -> str:
    """Remove version number from Ensembl ID."""
    return ensembl_id.split('.')[0]

def parse_gff_info(info_field: str) -> Dict[str, str]:
    """Parse the GFF3 info field into a dictionary."""
    info_dict = {}
    for item in info_field.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def format_gff_info(info_dict: Dict[str, str]) -> str:
    """Format dictionary back into GFF3 info field string."""
    return ';'.join(f"{key}={value}" for key, value in info_dict.items())

def process_gtex(gtex_file: str, selected_tissues: Set[str]) -> Tuple[Dict[str, int], Dict[str, Dict[str, float]]]:
    """
    Process GTEx file and extract expression data for selected tissues.
    Returns:
    - Dict mapping normalized tissue name to its column index
    - Dict mapping stripped Ensembl ID to dict of {normalized_tissue: expression_value}
    """
    tissue_indices = {}
    expression_data = {}
    
    with open(gtex_file, 'r') as f:
        # Skip version and dimensions lines
        f.readline()
        f.readline()
        
        # Read header and find indices of selected tissues
        header = f.readline().strip().split('\t')
        for i, tissue in enumerate(header[2:], 2):  # Start from 2 to account for Name and Description
            norm_tissue = normalize_tissue_name(tissue)
            if norm_tissue in selected_tissues:
                tissue_indices[norm_tissue] = i
        
        # Read expression data
        for line in f:
            fields = line.strip().split('\t')
            ensembl_id = strip_version(fields[0])
            
            # Create dict of tissue:expression for this gene
            gene_expression = {}
            for norm_tissue, idx in tissue_indices.items():
                try:
                    value = float(fields[idx]) if fields[idx] else 0.0
                    gene_expression[norm_tissue] = value
                except (IndexError, ValueError):
                    gene_expression[norm_tissue] = 0.0
            
            expression_data[ensembl_id] = gene_expression
    
    return tissue_indices, expression_data

def main():
    parser = argparse.ArgumentParser(description='Add GTEx expression data as attributes to GENCODE GFF3')
    parser.add_argument('gff_file', help='Input GENCODE GFF3 file (gzipped)')
    parser.add_argument('gtex_file', help='Input GTEx expression file')
    parser.add_argument('output_file', help='Output GFF3 file (will be gzipped)')
    parser.add_argument('--tissues', nargs='+', required=True,
                      help='One or more tissue names to include from GTEx data')
    args = parser.parse_args()

    # Normalize selected tissue names
    selected_tissues = {normalize_tissue_name(t) for t in args.tissues}
    
    # Process GTEx file to get expression data for selected tissues
    tissue_indices, expression_data = process_gtex(args.gtex_file, selected_tissues)
    
    # Process GFF file and add expression data
    with gzip.open(args.gff_file, 'rt') as fin, gzip.open(args.output_file, 'wt') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) != 9:
                fout.write(line)
                continue
            
            # Parse info field
            info_dict = parse_gff_info(fields[8])
            ensembl_id = info_dict.get('ID', '')
            if not ensembl_id:
                fout.write(line)
                continue
            
            # Get stripped version for matching
            stripped_id = strip_version(ensembl_id)
            
            # If this gene has expression data and it's a gene feature, add it
            if stripped_id in expression_data and fields[2] == 'gene':
                gene_expression = expression_data[stripped_id]
                for tissue, value in gene_expression.items():
                    info_dict[f'gtex_{tissue}'] = f'{value:.6f}'
                
                # Update info field with new attributes
                fields[8] = format_gff_info(info_dict)
                fout.write('\t'.join(fields) + '\n')
            else:
                fout.write(line)

if __name__ == '__main__':
    main() 