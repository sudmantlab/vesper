import argparse
from pathlib import Path
from typing import Optional
import logging
import sys
import os

from .utils.config import Config

def setup_logging(log_path: Optional[Path] = None, debug: bool = False) -> None:
    """Configure logging for the application.
    
    Args:
        log_path: Optional path to write log file. If None, only console logging is enabled.
        debug: If True, set log level to DEBUG
    """
    log_level = logging.DEBUG if debug else logging.INFO
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Console handler
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(log_level)
    console.setFormatter(logging.Formatter(log_format))
    root_logger.addHandler(console)
    
    # File handler if log path specified
    if log_path:
        # Create logs directory if it doesn't exist
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(str(log_path))
        file_handler.setLevel(log_level)
        file_handler.setFormatter(logging.Formatter(log_format))
        root_logger.addHandler(file_handler)

def parse_args() -> Config:
    """Parse command line arguments and return Config instance."""
    parser = argparse.ArgumentParser(
        description='Vesper: Variant Effect Predictor',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('-i', '--input', type=Path, required=True,
                       help='Input VCF file (must end with .vcf.gz)')
    parser.add_argument('-b', '--bam', type=Path, required=True,
                       help='Input BAM file (must end with .bam)')
    parser.add_argument('-o', '--output', type=Path, required=True,
                       help='Output VCF file (must end with .vcf.gz)')

    # Optional arguments
    parser.add_argument('--min-support', type=int, default=1,
                       help='Minimum support threshold (Default: 1)')
    parser.add_argument('--max-af', type=float, default=0.1,
                       help='Maximum allele frequency threshold (Default: 0.1)')
    parser.add_argument('--logging', type=Path, default=None, metavar='LOG_DIR',
                       help='(Directory) Saves logs to file in specified directory')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug-level logging (Default: False)')

    args = parser.parse_args()
    
    # Setup logging if enabled
    if args.logging:
        log_path = args.logging / f'vesper_{os.getpid()}.log'
        setup_logging(log_path, args.debug)
    else:
        setup_logging(debug=args.debug)
    
    return Config.from_args(args)

def main():
    """Main entry point for the CLI."""
    config = parse_args()
    # TODO: Add main program logic here
    pass
