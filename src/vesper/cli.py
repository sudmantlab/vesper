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


def add_common_args(parser):
    """Add arguments common to all commands."""
    parser.add_argument("--vcf", "-v", required=True, help="Input VCF file (required)")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory (required)")
    parser.add_argument("--logging", help="Log directory (default: output/logs)")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--console-output", action="store_true", help="Enable logging to console (default: False)")

class VesperArgumentParser(argparse.ArgumentParser):
    """Custom argument parser that shows program-specific help on error."""
    
    def error(self, message):
        """Upon error, prints help message and error."""
        display_ascii()
        self.print_help()
        self.exit(2, f'\n\033[31mERROR\033[0m: {message}\n')

def parse_args():
    """Parse command line arguments."""
    # Main parser
    parser = VesperArgumentParser(
        formatter_class=RawDescriptionRichHelpFormatter,
        epilog="View inputs & arguments for each command with vesper {command} --help"
    )
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # vesper annotate
    annotate_parser = subparsers.add_parser('annotate',
        help='Annotate variants with genomic features',
        description='Annotate structural variants with overlapping genomic features.',
        formatter_class=parser.formatter_class,
        epilog = """
Example:
  vesper annotate --vcf input.vcf --bed annotations.bed --output-dir output/
  vesper annotate --vcf input.vcf --output-dir output/ --test-mode 50
        """
        )
    add_common_args(annotate_parser)
    annotate_parser.add_argument("--bed", "-b", default="annotations/hg38/GRCH38_repeatmasker.bed",
                                help="BED file for annotations (default: annotations/hg38/GRCH38_repeatmasker.bed)")
    annotate_parser.add_argument("--test-mode", type=int, default=None,
                              help="Run in test mode with limited variants. Specify the number of variants to process (default: disabled)")
    
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
    add_common_args(refine_parser)
    refine_parser.add_argument("--bam", "-b", required=True, help="Input BAM file (required)")
    refine_parser.add_argument("--min-support", type=int, default=1,
                              help="Minimum supporting reads (default: 1)")
    refine_parser.add_argument("--max-af", type=float, default=0.1,
                              help="Maximum allele frequency (default: 0.1)")
    refine_parser.add_argument("--test-mode", type=int, default=None,
                              help="Run in test mode with limited variants. Specify the number of variants to process (default: disabled)")
    refine_parser.add_argument("--auto-load-registry", choices=['True', 'False'], default='True',
                              help="Whether to automatically load existing registry if found (default: True)")
    refine_parser.add_argument("--force-new-registry", choices=['True', 'False'], default='False',
                              help="Forces creation of new registry even if one exists (default: False)")

    # First parse args without handling errors to see what subcommand was chosen
    args, _ = parser.parse_known_args()
    
    # Show main help if no command specified
    if args.command is None:
        display_ascii()
        parser.print_help()
        sys.exit(0)
    
    # Get the appropriate subparser based on the command
    if args.command == 'annotate':
        subparser = annotate_parser
    elif args.command == 'refine':
        subparser = refine_parser
    
    # Make subparser also use our custom error handling
    subparser.__class__ = VesperArgumentParser
    
    # Parse all arguments (will show appropriate help on error)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    if args.command == 'refine':
        # Get the log directory from args or use default
        log_dir = Path(args.logging) if args.logging else Path(args.output_dir) / 'logs'
        
        # Set up file logging with the command name
        logger = setup_file_logging(log_dir, 'refine', args.debug, args.console_output)
        
        # Run the command
        run_refine(args, logger)
        
    elif args.command == 'annotate':
        # Get the log directory from args or use default
        log_dir = Path(args.logging) if args.logging else Path(args.output_dir) / 'logs'
        
        # Set up file logging with the command name
        logger = setup_file_logging(log_dir, 'annotate', args.debug, args.console_output)
        
        # Run the command
        run_annotate(args, logger)

if __name__ == "__main__":
    main()
