from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Union
import os

def get_threads(threads_arg):
    """Convert threads argument to an integer.
    
    Args:
        threads_arg: String or integer specifying threads ('auto' or int)
        
    Returns:
        int: Number of threads to use
    """
    if threads_arg == 'auto':
        return min(32, max(8, os.cpu_count() * 2))
    else:
        try:
            return int(threads_arg)
        except ValueError:
            try:
                parsed = int(float(threads_arg))
                print(f"WARNING: Only integer values are supported for threads. Converting {threads_arg} to {parsed}.")
                return parsed
            except ValueError:
                print(f"WARNING: Invalid threads argument {threads_arg} could not be parsed, falling back to default (8)")
                return 8

@dataclass
class CallConfig:
    """Configuration for the call command."""
    # Required arguments
    fastq_file: Path
    output_dir: Path
    
    # Optional arguments
    log_dir: Optional[Path] = None  # output_dir/logs if not specified
    debug: bool = False
    threads: int = 8

    @classmethod
    def from_args(cls, args):
        """Create CallConfig instance from parsed command line arguments."""
        return cls(
            fastq_file=Path(args.fastq),
            output_dir=Path(args.output_dir),
            log_dir=Path(args.logging) if args.logging else Path(args.output_dir) / 'logs',
            debug=args.debug,
            threads=get_threads(args.threads)
        )

@dataclass
class RefineConfig:
    """Configuration for the refine command."""
    # Required arguments
    vcf_input: Path
    output_dir: Path
    bam_file: Path
    
    # Optional arguments
    log_dir: Optional[Path] = None  # output_dir/logs if not specified
    debug: bool = False
    min_support: int = 1
    max_af: float = 0.1
    auto_load_registry: bool = True
    force_new_registry: bool = False
    test_mode: Optional[int] = None
    threads: int = 8

    @classmethod
    def from_args(cls, args):
        """Create RefineConfig instance from parsed command line arguments."""
        return cls(
            vcf_input=Path(args.vcf),
            output_dir=Path(args.output_dir),
            bam_file=Path(args.bam),
            log_dir=Path(args.logging) if args.logging else Path(args.output_dir) / 'logs',
            min_support=args.min_support,
            max_af=args.max_af,
            auto_load_registry=args.auto_load_registry == 'True',
            force_new_registry=args.force_new_registry,
            debug=args.debug,
            test_mode=args.test_mode,
            threads=get_threads(args.threads)
        )

@dataclass
class AnnotateConfig:
    """Configuration for the annotate command."""
    vcf_input: Path
    output_dir: Path
    
    # Optional additional annotation
    gff_files: List[Path] = None
    gff_names: List[str] = None

    # RepeatMasker
    proximal_span: int = 100
    repeatmasker_n: Optional[int] = None # number of top-scoring repeat annotations to return (None = all)
    
    # TODO: Trace logging inheritance 
    log_dir: Optional[Path] = None  # output_dir/logs if not specified
    debug: bool = False
    test_mode: Optional[int] = None
    threads: int = 8
    rebuild: bool = False
    
    @classmethod
    def from_args(cls, args):
        """Create AnnotateConfig instance from parsed command line arguments."""
        config = cls(
            vcf_input=Path(args.vcf),
            output_dir=Path(args.output_dir),
            log_dir=Path(args.logging) if args.logging else Path(args.output_dir) / 'logs',
            debug=args.debug,
            test_mode=args.test_mode,
            proximal_span=args.proximal_span,
            repeatmasker_n=args.repeatmasker_n if args.repeatmasker_n > 0 else None, # Interpret 0 from CLI as None (keep all)
            threads=get_threads(args.threads),
            rebuild=args.rebuild if hasattr(args, 'rebuild') else False
        )
        
        # Process GFF files and names
        config.gff_files = [Path(file) for file in args.files]
        
        if args.names:
            if len(args.names) != len(args.files):
                raise ValueError("Number of names must match number of files")
            config.gff_names = args.names
        else:
            # Default to using the file stem (filename without extension)
            config.gff_names = [Path(file).stem for file in args.files]
            
        return config