"""Module for handling RepeatMasker processing for insertion sequences."""

import os
import re
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import shutil
import tempfile
import uuid
import json
from glob import glob
from dataclasses import dataclass, field

from ..models.variants import VariantAnalysis
from ..models.repeatmasker import RepeatMaskerResult

class RepeatMaskerProcessor:
    """Process insertion sequences using RepeatMasker.
    
    This processor extracts insertion sequences from VariantAnalysis objects and
    runs RepeatMasker to identify repeat elements. Uses temporary directories
    for batch processing to avoid collisions during multithreading.
    """
    
    def __init__(self, output_dir: Path):
        """Initialize the RepeatMasker processor.
        
        Args:
            output_dir: Base output directory path
        """
        self.logger = logging.getLogger(__name__)
        self.output_dir = output_dir
        self.repeatmasker_dir = output_dir / "repeatmasker"
        self.current_temp_dir = None
        
        self.repeatmasker_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Created RepeatMasker output directory: {self.repeatmasker_dir}")
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._cleanup_temp_dir()
        
    def _create_temp_dir(self) -> Path:
        """Create a temporary directory for batch processing.
        
        Returns:
            Path to temporary directory
        """
        if self.current_temp_dir and self.current_temp_dir.exists():
            self._cleanup_temp_dir() # remove any existing temp dir
        
        self.current_temp_dir = self.repeatmasker_dir / f"batch_{uuid.uuid4().hex[:8]}"
        self.current_temp_dir.mkdir(parents=True)
        return self.current_temp_dir
    
    def _cleanup_temp_dir(self) -> None:
        """Clean up the current temporary directory."""
        if self.current_temp_dir and self.current_temp_dir.exists():
            shutil.rmtree(self.current_temp_dir)
            self.logger.debug(f"Cleaned up temporary directory: {self.current_temp_dir}")
            self.current_temp_dir = None
    
    def _write_insertion_fasta(self, variants: List[VariantAnalysis], temp_dir: Path) -> Path:
        """Extract insertion sequences to FASTA file. 
        Each insertion sequence is written as a separate record in the FASTA file,
        using the variant ID as the record header.
        
        Args:
            variants: List of VariantAnalysis objects to process
            temp_dir: Temporary directory for the batch processing
        
        Returns:
            Path to the created FASTA file
        """
        insertion_count = 0
        
        fasta_path = temp_dir / "insertions.fa"
        with open(fasta_path, 'w') as fasta_file:
            for v in variants:
                alt = v.variant.alt
                if alt not in ["N", ".", "<INS>"]:
                    fasta_file.write(f">{v.variant.ID}\n")
                    fasta_file.write(f"{alt}\n")
                    insertion_count += 1
        
        self.logger.info(f"Wrote {insertion_count} insertion sequences to {fasta_path}")
        return fasta_path
    
    def _run_repeatmasker(self, fasta_path: Path, temp_dir: Path) -> None:
        """Run RepeatMasker on a batch of insertion sequences.
        
        Args:
            fasta_path: Path to the batch insertion sequences FASTA file
            temp_dir: Temporary directory for the batch processing
        """
        if not fasta_path.exists():
            raise FileNotFoundError(f"Insertion FASTA file not found: {fasta_path}")
        
        command = [
            "RepeatMasker",
            "-engine", "rmblast",
            "-nocut",
            "-gff",
            "-species", "human",
            "-dir", str(temp_dir),
            str(fasta_path)
        ]
        
        self.logger.info(f"Running RepeatMasker with command: {' '.join(command)}")
        
        try:
            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            self.logger.debug(f"RepeatMasker output: {result.stdout}")
            self.logger.info("RepeatMasker completed successfully")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"RepeatMasker failed: {e.stderr}")
            raise RuntimeError(f"RepeatMasker failed: {e.stderr}")
        except FileNotFoundError:
            self.logger.error("RepeatMasker executable not found. Make sure it's installed and in your PATH.")
            raise RuntimeError("RepeatMasker executable not found. Make sure it's installed and in your PATH.")
        
        expected_files = [
            fasta_path.with_suffix(".fa.out"),
            fasta_path.with_suffix(".fa.tbl"),
            fasta_path.with_suffix(".fa.cat.gz"),
            fasta_path.with_suffix(".fa.masked"),
            fasta_path.with_suffix(".fa.out.gff")
        ]
        
        missing_files = [f for f in expected_files if not f.exists()]
        
        if missing_files:
            self.logger.warning(f"Missing expected RepeatMasker output files: {', '.join(str(f) for f in missing_files)}")

    def _parse_output(self, temp_dir: Path) -> Dict[str, List[RepeatMaskerResult]]:
        """Parse a RepeatMasker .out file and return all annotations.
        
        Args:
            temp_dir: Temporary directory for the batch processing
                
        Returns:
            Dictionary mapping sequence IDs to repeat annotations
        """
        outfile = temp_dir / "insertions.fa.out"
        results = {}
        if not outfile.exists():
            self.logger.debug(f"No RepeatMasker output file found: {outfile}")
            return results
        
        with open(outfile, 'r') as f:
            for line in f: # skip header
                if "score" in line.lower():
                    break
            next(f, None)
            
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                try:
                    # note and remove trailing asterisks (doc info below)
                    # asterisk (*) in the final column = annotation partially overlaps a higher-scoring match (<80%)
                    # results are sorted in the end anyway
                    if line.endswith('*'):
                        line = line[:-1].strip()
                    fields = re.split(r'\s+', line.strip())
                    if not fields or len(fields) < 14:
                        self.logger.warning(f"Skipping invalid line (not enough fields): {line}")
                        continue
                    
                    sw_score = int(fields[0])
                    perc_div = float(fields[1])
                    perc_del = float(fields[2])
                    perc_ins = float(fields[3])
                    query_name = fields[4]
                    query_begin = int(fields[5])
                    query_end = int(fields[6])
                    query_left_raw = fields[7]
                    query_left = int(query_left_raw.strip('()')) if query_left_raw.strip('()').isdigit() else 0 # ex. (100) -> 100
                    strand = fields[8]
                    
                    # begin annoying repeat class/family parsing
                    repeat_name = fields[9]
                    repeat_class = None
                    position_shift = 0
                    
                    # look for class/family field in the next 3 fields
                    # easiest detection: repeat class/family indicated by slash for many categories
                    for i in range(10, min(13, len(fields))):
                        if '/' in fields[i]: # ex. LTR/ERVL, Retroposon/SVA, LINE/L1...
                            repeat_class = fields[i]
                            position_shift = i - 10 + 1 
                            break
                    
                    # if not found, likely Simple_repeat, Low_complexity, etc...
                    # we then have an extra field that's not part of position-in-repeat coordinates
                    if repeat_class is None and len(fields) > 14 and not all(f.strip('()').isdigit() for f in fields[10:13]):
                        repeat_class = fields[10]
                        position_shift = 1
                    
                    # adjust field index for parsing remaining fields based on position_shift
                    repeat_start_idx = 10 + position_shift
                    
                    if repeat_start_idx + 2 >= len(fields):
                        self.logger.warning(f"Skipping invalid line (can't locate position-in-repeat fields): {line}")
                        continue
                    
                    # finally, parse the position-in-repeat fields (next three after the adjusted index)
                    repeat_start_raw = fields[repeat_start_idx]
                    repeat_start = int(repeat_start_raw.strip('()')) if repeat_start_raw.strip('()').isdigit() else 0
                    
                    repeat_end = int(fields[repeat_start_idx+1])
                    
                    repeat_left_raw = fields[repeat_start_idx+2]
                    repeat_left = int(repeat_left_raw.strip('()')) if repeat_left_raw.strip('()').isdigit() else 0
                    
                    result = RepeatMaskerResult(
                        repeat_name=repeat_name,
                        repeat_class=repeat_class,
                        sw_score=sw_score,
                        divergence=perc_div,
                        deletion=perc_del,
                        insertion=perc_ins,
                        query_start=query_begin,
                        query_end=query_end,
                        query_left=query_left,
                        strand=strand,
                        repeat_start=repeat_start,
                        repeat_end=repeat_end,
                        repeat_left=repeat_left,
                        match_length=round(query_end - query_begin + 1, 2),
                        match_coverage=(query_end - query_begin + 1) / (query_end + query_left),
                    )
                    
                    if query_name not in results:
                        results[query_name] = []
                    results[query_name].append(result)
                    
                except Exception as e:
                    self.logger.warning(f"Error parsing line: {line}. Error: {e}")
                    continue
        
        return results
    
    def _parse_and_sort(self, temp_dir: Path) -> Dict[str, List[RepeatMaskerResult]]:
        """For each query sequence, sort and return annotations.
        
        Args:
            temp_dir: Path to the directory containing RepeatMasker output files
            
        Returns:
            List of annotations (dictionaries) for each query sequence.
        """
        try:
            all_results = self._parse_output(temp_dir)
        except Exception as e:
            self.logger.warning(f"Error parsing RepeatMasker output: {e}. Continuing...")
            return {}
        
        sorted_results = {}
        
        for query_name, annotations in all_results.items():
            if not annotations: # insertion seq not IDed w/ repetitive motif
                continue
                
            sorted_annotations = sorted(annotations, key=lambda x: (x.match_length), reverse=True)
            sorted_results[query_name] = sorted_annotations
        
        return sorted_results
    
    def assign_results(self, variant_id: str, results: Dict[str, List[RepeatMaskerResult]], n: Optional[int] = None) -> List[RepeatMaskerResult]:
        """Query RepeatMasker results for the given variant ID and optionally return the top n results.
        
        Args:
            variant_id: ID of the variant to query
            results: Dictionary mapping sequence IDs to repeat annotations
            n: Optional number of top results to return. If None, returns all results.
            
        Returns:
            List of RepeatMaskerResult objects corresponding to the variant ID
        """
        if variant_id in results:
            if n is not None:
                return results[variant_id][:n]
            return results[variant_id]
        return []
    
    def batch_analysis(self, variants: List[VariantAnalysis], chunk_idx: int, n: Optional[int] = None) -> None:
        """Analyze a batch of variants with RepeatMasker.
        
        Args:
            variants: List of variants to analyze
            n: Optional number of top results to return. If None, returns all results.
        """
        if not variants:
            return
            
        temp_dir = self._create_temp_dir()
        try:
            fasta_path = self._write_insertion_fasta(variants, temp_dir)
            if fasta_path.stat().st_size == 0:
                self.logger.info("No insertions found in batch - skipping RepeatMasker analysis")
                return
            
            self._run_repeatmasker(fasta_path, temp_dir)
            results = self._parse_and_sort(temp_dir)
            
            for variant in variants:
                variant.repeatmasker_results = self.assign_results(variant.variant.ID, results, n)

            # TODO: can this be done more cleanly?
            temp_json = self.repeatmasker_dir / f"rm.chunk_{chunk_idx}.json"
            results_dict = {
                var_id: [result.to_dict() for result in result_list] 
                for var_id, result_list in results.items()
            }
            with open(temp_json, 'w') as f:
                json.dump(results_dict, f)
                    
            self.logger.info(f"Processed {len(variants)} variants through RepeatMasker, {len(results)} results found")
        except Exception as e:
            self.logger.error(f"Error processing variants with RepeatMasker: {e}")
            raise e
        finally:
            self._cleanup_temp_dir()

    @staticmethod
    def merge_temp_jsons(directory: Path, final_json_path: Path) -> None:
        """Merge all temporary JSON files into a single final JSON file.
        
        Args:
            directory: Directory containing temporary JSON files
            final_json_path: Path for final merged JSON file
        """
        merged_results = {}
        temp_jsons = glob(str(directory / "*.chunk_*.json"))
        
        for temp_json in temp_jsons:
            with open(temp_json, 'r') as f:
                chunk_results = json.load(f)
                merged_results.update(chunk_results)
            os.remove(temp_json)
            
        with open(final_json_path, 'w') as f:
            json.dump(merged_results, f, indent=4)