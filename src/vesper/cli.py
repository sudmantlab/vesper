#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path
import argparse
import sys
from rich_argparse import RawDescriptionRichHelpFormatter

from vesper.commands import run_refine, run_annotate, run_summarize
from vesper.utils import setup_file_logging, SummarizeConfig

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
    - vesper annotate: provide a VCF file to annotate insertions with RepeatMasker (external annotations optional).
    - vesper refine: provide a VCF file and BAM to calculate read-based confidence scores.
    - vesper summarize: provide one or more refined VCFs to generate a TSV summary.

    View inputs & arguments for each command with vesper {command} --help.
        """
    )
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # vesper annotate
    annotate_parser = subparsers.add_parser('annotate',
        help='Annotate insertion variants with RepeatMasker and optionally with genomic features',
        description='Annotate insertion sequences with RepeatMasker and optionally overlay genomic features.',
        formatter_class=parser.formatter_class,
        epilog = """
Examples:
vesper annotate --vcf input.vcf --output-dir output/
vesper annotate --vcf input.vcf --files annotations/segdups.gff annotations/repeats.gff --names segdups repeats --output-dir output/
        """
        )
    annotate_parser.add_argument("--vcf", "-v", required=True,
                                help="Input VCF file (required)")
    annotate_parser.add_argument("--output-dir", "-o", required=True,
                                help="Output directory (required)")
    
    annotate_files = annotate_parser.add_argument_group('optional external annotations')
    annotate_files.add_argument("--files", "-f", nargs='+', default=[],
                                help="GFF/GTF annotation file paths")
    annotate_files.add_argument("--names", "-n", nargs='+', default=[],
                                help="Space-delimited shorthand names for annotation files. Required if --files is used, must match number of files.")
    annotate_parser.add_argument("--proximal-span", type=int, default=100,
                                help="Distance (+/-) in base pairs to search for proximal features (default: 100)")
    annotate_parser.add_argument("--repeatmasker-n", type=int, default=0,
                                help="Number of top-scoring repeat annotations to return per insertion (default: all)")
    annotate_parser.add_argument("--threads", "-t", type=int, default=4,
                                help="Number of threads (default: 4)")
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
    refine_parser.add_argument("--threads", type=int, default=4,
                              help="Number of threads (default: 4)")
    refine_parser.add_argument("--logging",
                              help="Log directory (default: output/logs)")
    refine_parser.add_argument("--debug", action="store_true",
                              help="Enable debug logging")
    refine_parser.add_argument("--console-output", action="store_true",
                              help="Enable logging to console (default: False)")

    # vesper summarize
    summarize_parser = subparsers.add_parser('summarize',
        help='Summarize refined VCFs into a TSV',
        description='Summarize refined VCFs into a TSV with confidence and annotation fields.',
        formatter_class=parser.formatter_class,
        epilog="""
Examples:
  vesper summarize --input sample1.vcf.gz sample2.vcf.gz --output summary.tsv
  vesper summarize --input sample1.vcf.gz --output summary.tsv --sample-names SAMPLE1
        """
        )
    summarize_parser.add_argument("--input", "-i", required=True, nargs='+',
                                  help="Input VCF file(s) (bgzipped and indexed). Provide one or more paths.")
    summarize_parser.add_argument("--output", "-o", required=True,
                                  help="Output TSV file.")
    summarize_parser.add_argument("--sample-names", "-s", nargs='+',
                                  help="Optional sample names corresponding to each input. Must match number of inputs.")
    summarize_parser.add_argument("--logging",
                                  help="Log directory (default: <output-dir>/logs)")
    summarize_parser.add_argument("--debug", action="store_true",
                                  help="Enable debug logging")
    summarize_parser.add_argument("--console-output", action="store_true",
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
    elif args.command == 'summarize':
        subparser = summarize_parser

    
    args = parser.parse_args()
    
    if args.command == 'annotate':
        # Only validate if files are provided
        if args.files:
            # Check that names are provided when files are provided
            if not args.names:
                subparser.error("ERROR: --names is required when --files is provided")
            
            # Check that number of names matches number of files
            if len(args.names) != len(args.files):
                subparser.error(f"ERROR: Number of names ({len(args.names)}) must match number of files ({len(args.files)})")
            
            # Check for duplicate names
            if len(set(args.names)) != len(args.names):
                subparser.error("ERROR: Duplicate names found in --names. Names must be unique.")
            
            # Validate that all files exist
            for file_path in args.files:
                if not os.path.exists(file_path):
                    subparser.error(f"ERROR: File does not exist: {file_path}")
    elif args.command == 'summarize':
        if args.sample_names and len(args.sample_names) != len(args.input):
            subparser.error(f"ERROR: Number of sample names ({len(args.sample_names)}) must match number of inputs ({len(args.input)})")

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
    elif args.command == 'summarize':
        config = SummarizeConfig.from_args(args)
        logger = setup_file_logging(config.log_dir, 'summarize', config.debug, config.console_output)
        run_summarize(config, logger)

if __name__ == "__main__":
    main()
