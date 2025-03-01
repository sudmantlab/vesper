from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class Config:
    """Configuration class for parameters."""
    # Input/output files
    vcf_input: Path
    bam_file: Path
    vcf_output: Path
    
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
            vcf_input=Path(args.input),
            bam_file=Path(args.bam),
            vcf_output=Path(args.output),
            min_support=args.min_support,
            max_af=args.max_af,
            log_dir=args.logging,
            debug=args.debug
        )
