#!/usr/bin/env python3
import os
from pathlib import Path
import logging
import argparse
import sys
from rich_argparse import RawDescriptionRichHelpFormatter

from vesper.commands import run_refine, run_annotate
from vesper.utils import setup_file_logging

def display_ascii():
        print("""
                    『 ° °*.ﾟ°*. ﾟ｡｡*･｡° °
                     .ﾟ°*. ﾟ \033[1mvesper\033[0m ﾟ｡｡*･｡
                      ﾟ° ﾟ.ﾟ° ﾟ｡｡*･｡° ﾟ*° 』

        \033[3ma haplotype-aware de novo structural variant caller\033[0m
        """)

class VesperArgumentParser(argparse.ArgumentParser):
    """Custom argument parser that shows program-specific help on error."""
    
    def error(self, message):
        """Upon error, prints help message and error."""
        display_ascii()
        self.print_help()
        self.exit(2, f'\n\033[31mERROR\033[0m: {message}\n')

def parse_args():
    """Parse command line arguments."""
    parser = VesperArgumentParser(
        formatter_class=RawDescriptionRichHelpFormatter,
        epilog="""
    Recommended usage:
    vesper call (optional + not yet available) → vesper annotate → vesper refine

    - vesper call: provide a FASTQ file of PacBio HiFi reads for alignment and variant calling.
    - vesper annotate: provide a VCF file of candidate variants to be annotated.
    - vesper refine: provide a VCF file of annotated candidate variants and the origin BAM file.

    View inputs & arguments for each command with vesper {command} --help.
        """
    )
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # vesper call
    call_parser = subparsers.add_parser('call',
        help='Align reads* and call de novo structural variant candidates',
        description='Align reads* and call de novo structural variant candidates.',
        formatter_class=parser.formatter_class,
        epilog = """
Example:
  vesper call --fastq input.fastq --output-dir output/
        """
        )
    call_parser.add_argument("--fastq", "-f", required=True,
                            help="fastq/fastq.gz file for alignment (required)")
    call_parser.add_argument("--output-dir", "-o", required=True,
                            help="Output directory (required)")
    call_parser.add_argument("--threads", "-t", 
                            default="8",
                            help="Number of threads to use (integer or 'auto', default: 8)")
    call_parser.add_argument("--logging",
                            help="Log directory (default: output/logs)")
    call_parser.add_argument("--debug", action="store_true",
                            help="Enable debug logging")
    call_parser.add_argument("--console-output", action="store_true",
                            help="Enable logging to console (default: False)")

    # vesper annotate
    annotate_parser = subparsers.add_parser('annotate',
        help='Annotate variants with genomic features',
        description='Annotate structural variants with overlapping genomic features.',
        formatter_class=parser.formatter_class,
        epilog = """
Usage notes:
  - Provide annotation files with --files (file format is auto-detected by extension)
  - Use --names to provide shorthand names for each file (optional, defaults to filenames)
  - Names must be unique across all annotation files

Examples:
  # Single annotation file
  vesper annotate --vcf input.vcf --files annotations.bed --names repeats --output-dir output/
  vesper annotate --vcf input.vcf --files annotations.gff3.gz --names genes --output-dir output/
  
  # Multiple annotation files
  vesper annotate --vcf input.vcf --files annotations1.bed annotations2.bed --names repeats centromeres --output-dir output/
  vesper annotate --vcf input.vcf --files annotations1.gff3.gz annotations2.tsv --names genes features --output-dir output/
        """
        )
    annotate_parser.add_argument("--vcf", "-v", required=True,
                                help="Input VCF file (required)")
    annotate_parser.add_argument("--output-dir", "-o", required=True,
                                help="Output directory (required)")
    
    annotate_files = annotate_parser.add_argument_group('annotation files')
    annotate_files.add_argument("--files", "-f", nargs='+', default=[],
                                help="Annotation files (BED, GFF/GTF, TSV). File format is auto-detected by extension.")
    annotate_files.add_argument("--names", "-n", nargs='+', default=[],
                                help="Shorthand names for annotation files. Must match number of files.")
    annotate_parser.add_argument("--proximal-span", type=int, default=100,
                                help="Distance (+/-) in base pairs to search for proximal features (default: 100)")
    annotate_parser.add_argument("--repeatmasker-n", type=int, default=0,
                                help="Number of top-scoring repeat annotations to return per insertion (default: return all annotations)")
    annotate_parser.add_argument("--threads", "-t", type=str, default="8",
                                help="Number of threads to use (integer or 'auto' (defaults to 2x CPU count, up to 32), default: 8)")
    annotate_parser.add_argument("--test-mode", type=int, default=None,
                                help="Run in test mode with limited variants. Specify the number of variants to process (default: disabled)")
    annotate_parser.add_argument("--logging",
                                help="Log directory (default: output/logs)")
    annotate_parser.add_argument("--debug", action="store_true",
                                help="Enable debug logging")
    annotate_parser.add_argument("--console-output", action="store_true",
                                help="Enable logging to console (default: False)")
    annotate_parser.add_argument("--rebuild", action="store_true",
                                help="Force rebuild of all annotation databases")
    
    # vesper refine
    refine_parser = subparsers.add_parser('refine', 
        help='Refine variants using read evidence',
        description='Refine structural variants using supporting read evidence.',
        formatter_class=parser.formatter_class,
        epilog="""
Example:
  vesper refine --vcf input.vcf --bam input.bam --output-dir output/
  vesper refine --vcf input.vcf --bam input.bam --output-dir output/ --test-mode 50
        """)
    refine_parser.add_argument("--vcf", "-v", required=True,
                              help="Input VCF file (required)")
    refine_parser.add_argument("--bam", "-b", required=True,
                              help="Input BAM file (required)")
    refine_parser.add_argument("--output-dir", "-o", required=True,
                              help="Output directory (required)")
    refine_parser.add_argument("--min-support", type=int, default=1,
                              help="Minimum supporting reads (default: 1)")
    refine_parser.add_argument("--max-af", type=float, default=0.1,
                              help="Maximum allele frequency (default: 0.1)")
    refine_parser.add_argument("--test-mode", type=int, default=None,
                              help="Run in test mode with limited variants. Specify the number of variants to process (default: disabled)")
    refine_parser.add_argument("--threads", type=str, default="8",
                              help="Number of threads to use (integer or 'auto', default: 8)")
    refine_parser.add_argument("--auto-load-registry", choices=['True', 'False'], default='True',
                              help="Whether to automatically load existing registry if found (default: True)")
    refine_parser.add_argument("--force-new-registry", action="store_true",
                              help="Forces rebuilding of registry even if one exists.")
    refine_parser.add_argument("--logging",
                              help="Log directory (default: output/logs)")
    refine_parser.add_argument("--debug", action="store_true",
                              help="Enable debug logging")
    refine_parser.add_argument("--console-output", action="store_true",
                              help="Enable logging to console (default: False)")

    args, _ = parser.parse_known_args()
    if args.command is None:
        display_ascii()
        parser.print_help()
        sys.exit(0)

    if args.command == 'annotate':
        subparser = annotate_parser
    elif args.command == 'refine':
        subparser = refine_parser

    
    args = parser.parse_args()
    
    if args.command == 'annotate':
        if not args.files:
            subparser.error("ERROR: At least one annotation file must be provided with --files")
        
        # Check that number of names matches number of files when names are provided
        if args.names and len(args.names) != len(args.files):
            subparser.error(f"ERROR: Number of names ({len(args.names)}) must match number of files ({len(args.files)})")
        
        # Check for duplicate names
        if args.names and len(set(args.names)) != len(args.names):
            subparser.error("ERROR: Duplicate names found in --names. Names must be unique.")
        
        # Validate that all files exist
        for file_path in args.files:
            if not os.path.exists(file_path):
                subparser.error(f"ERROR: File does not exist: {file_path}")
    
    return args

def main():
    args = parse_args()
    
    if args.command == 'refine':
        log_dir = Path(args.logging) if args.logging else Path(args.output_dir) / 'logs'
        logger = setup_file_logging(log_dir, 'refine', args.debug, args.console_output)
        run_refine(args, logger)
        
    elif args.command == 'annotate':
        log_dir = Path(args.logging) if args.logging else Path(args.output_dir) / 'logs'
        logger = setup_file_logging(log_dir, 'annotate', args.debug, args.console_output)
        run_annotate(args, logger)

if __name__ == "__main__":
    main()
