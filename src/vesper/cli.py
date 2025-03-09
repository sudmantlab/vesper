#!/usr/bin/env python3
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
        self.print_help()
        self.exit(2, f'\n\033[31mERROR\033[0m: {message}\n')
        
    def print_help(self, file=None):
        super().print_help(file)
        
    def format_help(self):
        """Format help message with ASCII art header."""
        display_ascii()
        return super().format_help()

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
  - At least one of --bed or --gff must be provided.
  - Use --bed-names and/or --gff-names to provide shorthand names for each file.
  - Names must be unique across all annotation files.

Examples:
  # Single annotation file
  vesper annotate --vcf input.vcf --bed annotations.bed --bed-names repeats --output-dir output/
  vesper annotate --vcf input.vcf --gff annotations.gff3.gz --gff-names genes --output-dir output/
  
  # Multiple annotation files
  vesper annotate --vcf input.vcf --bed annotations1.bed annotations2.bed --bed-names repeats centromeres --output-dir output/
  vesper annotate --vcf input.vcf --gff annotations1.gff3.gz annotations2.gff3.gz --gff-names genes regulatory --output-dir output/
        """
        )
    annotate_parser.add_argument("--vcf", "-v", required=True,
                                help="Input VCF file (required)")
    annotate_parser.add_argument("--output-dir", "-o", required=True,
                                help="Output directory (required)")
    
    annotate_files = annotate_parser.add_argument_group('annotation files')
    annotate_files.add_argument("--bed", "-b", nargs='+', default=[],
                                help="BED file(s) for annotations. Multiple files can be provided.")
    annotate_files.add_argument("--bed-names", "-bn", nargs='+', default=[],
                                help="Shorthand names for BED files (required). Must match number of BED files.")
    annotate_files.add_argument("--gff", "-g", nargs='+', default=[],
                                help="GFF/GTF file(s) for annotations. Multiple files can be provided.")
    annotate_files.add_argument("--gff-names", "-gn", nargs='+', default=[],
                                help="Shorthand names for GFF files (required). Must match number of GFF files.")
    
    annotate_parser.add_argument("--proximal-span", type=int, default=100,
                                help="Distance (+/-) in base pairs to search for proximal features (default: 100)")
    annotate_parser.add_argument("--test-mode", type=int, default=None,
                                help="Run in test mode with limited variants. Specify the number of variants to process (default: disabled)")
    annotate_parser.add_argument("--logging",
                                help="Log directory (default: output/logs)")
    annotate_parser.add_argument("--debug", action="store_true",
                                help="Enable debug logging")
    annotate_parser.add_argument("--console-output", action="store_true",
                                help="Enable logging to console (default: False)")
    
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
        parser.print_help()
        sys.exit(0)

    if args.command == 'annotate':
        subparser = annotate_parser
    elif args.command == 'refine':
        subparser = refine_parser

    
    args = parser.parse_args()
    
    if args.command == 'annotate':
        if not (args.bed or args.gff):
            subparser.error("ERROR: At least one of --bed or --gff must be provided")
        
        if (args.bed and args.gff) and not (args.bed_names and args.gff_names):
            subparser.error("ERROR: Use --bed-names (-bn) and --gff-names (-gn) to provide shorthand names for each file.")
        if args.bed and not args.bed_names:
            subparser.error("ERROR: Add shorthand name(s) for BED files with --bed-names (-bn)")
        if args.gff and not args.gff_names:
            subparser.error("ERROR: Add shorthand name(s) for GFF files with --gff-names (-gn)")
            
        # Check that number of names matches number of files
        if args.bed_names and len(args.bed_names) != len(args.bed):
            subparser.error(f"ERROR: Number of BED names ({len(args.bed_names)}) must match number of BED files ({len(args.bed)})")
        if args.gff_names and len(args.gff_names) != len(args.gff):
            subparser.error(f"ERROR: Number of GFF names ({len(args.gff_names)}) must match number of GFF files ({len(args.gff)})")

        # Check for duplicate names
        if args.bed_names and len(set(args.bed_names)) != len(args.bed_names):
            subparser.error("ERROR: Duplicate names found in --bed-names. Names must be unique.")
        if args.gff_names and len(set(args.gff_names)) != len(args.gff_names):
            subparser.error("ERROR: Duplicate names found in --gff-names. Names must be unique.")
        if args.bed_names and args.gff_names:
            all_names = args.bed_names + args.gff_names
            if len(set(all_names)) != len(all_names):
                subparser.error("ERROR: Names must be unique across both BED and GFF files.")
    
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
