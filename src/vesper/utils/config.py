from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

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
    threads: int = 4

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
            threads=args.threads
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
    threads: int = 4
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
            threads=args.threads,
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

@dataclass
class SummarizeConfig:
    """Configuration for the summarize command."""
    vcf_inputs: List[Path]
    output_path: Path
    sample_names: Optional[List[str]] = None
    log_dir: Optional[Path] = None
    debug: bool = False
    console_output: bool = False

    @classmethod
    def from_args(cls, args):
        """Create SummarizeConfig instance from parsed command line arguments."""
        input_paths = [Path(path) for path in args.input]
        sample_names = args.sample_names if getattr(args, 'sample_names', None) else None

        if sample_names and len(sample_names) != len(input_paths):
            raise ValueError("Number of sample names must match number of input files")

        output_path = Path(args.output)
        log_dir = Path(args.logging) if getattr(args, 'logging', None) else output_path.parent / 'logs'

        return cls(
            vcf_inputs=input_paths,
            output_path=output_path,
            sample_names=sample_names,
            log_dir=log_dir,
            debug=getattr(args, 'debug', False),
            console_output=getattr(args, 'console_output', False)
        )
