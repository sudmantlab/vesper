#!/usr/bin/env python3

"""
Script to combine GENCODE GFF3 annotations with GTEx expression data into a BED format file.
The script matches genes based on Ensembl IDs (stripped of version numbers) and combines
positional information from GFF3 with expression data from GTEx.
"""

import argparse
import gzip
import csv
import re
from typing import Dict, Tuple, List

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

def process_gff(gff_file: str) -> Dict[str, Tuple[str, int, int, str]]:
    """
    Process GFF3 file and extract gene information.
    Returns dict mapping stripped Ensembl ID to (chromosome, start, end, gene_name).
    Only processes gene features.
    """
    gene_info = {}
    with gzip.open(gff_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) != 9 or fields[2] != 'gene':
                continue
                
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            info = parse_gff_info(fields[8])
            
            # Extract Ensembl ID and gene name
            ensembl_id = info.get('ID', '')
            if not ensembl_id:
                continue
                
            stripped_id = strip_version(ensembl_id)
            gene_name = info.get('gene_name', stripped_id)
            
            gene_info[stripped_id] = (chrom, start, end, gene_name)
    
    return gene_info

def process_gtex(gtex_file: str) -> Tuple[List[str], Dict[str, List[float]]]:
    """
    Process GTEx file and extract expression data.
    Returns tuple of (tissue_names, expression_data_dict).
    Expression data dict maps stripped Ensembl ID to list of expression values.
    """
    expression_data = {}
    tissue_names = []
    
    with open(gtex_file, 'r') as f:
        # Skip version and dimensions lines
        f.readline()
        f.readline()
        
        # Read header
        header = f.readline().strip().split('\t')
        tissue_names = header[2:]  # Skip Name and Description columns
        
        # Read expression data
        for line in f:
            fields = line.strip().split('\t')
            ensembl_id = strip_version(fields[0])
            expression_values = [float(x) if x != '' else 0.0 for x in fields[2:]]
            expression_data[ensembl_id] = expression_values
    
    return tissue_names, expression_data

def main():
    parser = argparse.ArgumentParser(description='Combine GENCODE GFF3 and GTEx data into BED format')
    parser.add_argument('gff_file', help='Input GENCODE GFF3 file (gzipped)')
    parser.add_argument('gtex_file', help='Input GTEx expression file')
    parser.add_argument('output_file', help='Output BED file')
    args = parser.parse_args()

    # Process input files
    gene_info = process_gff(args.gff_file)
    tissue_names, expression_data = process_gtex(args.gtex_file)
    
    # Write output BED file
    with open(args.output_file, 'w') as f:
        # Write header
        header = ['#chrom', 'start', 'end', 'ensembl_id', 'gene_name'] + tissue_names
        f.write('\t'.join(header) + '\n')
        
        # Write data
        for ensembl_id in sorted(gene_info.keys()):
            if ensembl_id in expression_data:
                chrom, start, end, gene_name = gene_info[ensembl_id]
                expression_values = expression_data[ensembl_id]
                
                output_fields = [
                    chrom,
                    str(start),
                    str(end),
                    ensembl_id,
                    gene_name
                ] + [str(x) for x in expression_values]
                
                f.write('\t'.join(output_fields) + '\n')

if __name__ == '__main__':
    main() 