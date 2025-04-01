"""Logging utilities for Vesper."""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


def setup_file_logging(log_dir: Path, command_name: str, debug: bool = False, console_output: bool = False) -> logging.Logger:
    """Configure logging to write to a timestamped log file and optionally to console.
    
    Args:
        log_dir: Directory to write log files to
        command_name: Name of the command being run (e.g., 'annotate', 'refine')
        debug: Whether to enable debug logging
        console_output: Whether to also log to console (default: False)
    
    Returns:
        Configured logger instance
    """
    # Create log directory if it doesn't exist
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Create timestamped log filename
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    log_file = log_dir / f"{timestamp}.{command_name}.log"
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if debug else logging.INFO)
    
    # Clear any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create console handler if requested
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG if debug else logging.INFO)
        console_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_format)
        root_logger.addHandler(console_handler)
    
    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_format)
    root_logger.addHandler(file_handler)
    
    # Create and return the vesper logger
    logger = logging.getLogger('vesper')
    print(f"{timestamp} - Logging to file: {log_file}")
    
    return logger 