from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class Config:
    """Configuration class for parameters."""
    # Input/output files
    vcf_input: Path
    bam_file: Path
    output_dir: Path
    
    # Registry configuration
    auto_load_registry: bool = True  # Whether to automatically load existing registry if found
    force_new_registry: bool = False  # Force creation of new registry even if one exists
    
    # Parameters
    min_support: int = 1
    max_af: float = 0.1
    
    # Logging configuration
    log_dir: Optional[Path] = None
    debug: bool = False
    
    @classmethod
    def from_args(cls, args):
        """Create Config instance from parsed command line arguments."""
        return cls(
            vcf_input=Path(args.vcf),
            bam_file=Path(args.bam),
            output_dir=Path(args.output_dir),
            auto_load_registry=args.auto_load_registry == 'True',  # Convert string to bool
            force_new_registry=args.force_new_registry,
            min_support=args.min_support,
            max_af=args.max_af,
            log_dir=args.logging,
            debug=args.debug
        )
