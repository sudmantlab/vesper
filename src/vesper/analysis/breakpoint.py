from dataclasses import dataclass
from typing import Optional, Tuple, List
import re
import pysam
from collections import Counter
import edlib

from ..models.variants import Variant, VariantAnalysis, SVType
from ..models.reads import AlignedRead

@dataclass
class BreakpointContext:
    """Represents the referencesequence context around a variant breakpoint."""
    left_context: str  # Sequence upstream of breakpoint
    right_context: str  # Sequence downstream of breakpoint
    breakpoint_pos: int  # Position of the breakpoint
    microhomology: Optional[str] = None  # Microhomology sequence if present
    tsd: Optional[str] = None  # Target site duplication if present
    poly_a_tail: Optional[str] = None  # Poly-A tail if present

# pysam CIGAR operation codes: BAM_CMATCH 0, BAM_CINS 1, BAM_CDEL 2, BAM_CREF_SKIP 3,
# BAM_CSOFT_CLIP 4, BAM_CHARD_CLIP 5, BAM_CPAD 6, BAM_CEQUAL 7, BAM_CDIFF 8
CIGAR_CONSUMES_QUERY = {0, 1, 4, 7, 8} # M, I, S, =, X


class BreakpointAnalyzer:
    """Analyzes breakpoint contexts for structural variants."""
    
    def __init__(self, reference_fasta: str, context_size: int = 25):
        """
        Args:
            reference_fasta: Path to reference genome FASTA file (must be indexed, e.g., with samtools faidx)
            context_size: Number of bases to extract on each side of breakpoint
        """
        self.context_size = context_size
        self.reference_fasta = pysam.FastaFile(reference_fasta)
        
    def _get_reference_sequence(self, chrom: str, start: int, end: int) -> str:
        """Get reference sequence for a genomic region using pysam."""
        try:
            return self.reference_fasta.fetch(chrom, start, end)
        except ValueError as e:
            # pysam raises ValueError for fetch on invalid region (e.g., start < 0)
            # or if chromosome not found. Re-raise or handle appropriately.
            # For now, return empty string, but logging/error might be better.
            print(f"Error fetching sequence for {chrom}:{start}-{end}: {e}") # Consider proper logging
            return ""
        except KeyError as e:
            print(f"Error: Chromosome '{chrom}' not found in FASTA file {self.reference_fasta.filename}. Ensure FASTA and index are correct.")
            return ""
        
    def _find_insertion_in_read(self,read: AlignedRead, expected_len: int) -> Optional[Tuple[int, int, int]]:
        """
        Parses CIGAR to find the main insertion operation and its position in the read sequence.
        Returns (insertion_start_in_read, insertion_end_in_read, actual_ins_len) or None.
        """
        if not read.cigartuples:
            return None

        # Simple approach: Find the largest 'I' operation close to expected_len
        # More complex logic might be needed for fragmented insertions or complex CIGARs
        best_ins_op = None # (op_len, query_start)
        current_query_idx = 0

        for op, length in read.cigartuples:
            if op == 1: # Insertion
                # Check if length is reasonably close to expected_len (e.g., within 10% or 10bp)
                if abs(length - expected_len) <= max(10, expected_len * 0.1):
                    if best_ins_op is None or length > best_ins_op[0]:
                        best_ins_op = (length, current_query_idx)

            if op in CIGAR_CONSUMES_QUERY:
                current_query_idx += length

        if best_ins_op:
            ins_len, ins_start_in_read = best_ins_op
            ins_end_in_read = ins_start_in_read + ins_len
            return ins_start_in_read, ins_end_in_read + 1, ins_len # 0-based, exclusive end
        else:
            # TODO: Could add logic here to check for large soft clips as alternative evidence
            pass

        return None

    def _find_insertion_by_sequence(self,read: AlignedRead, ins_seq: str, max_edit_distance_ratio: float = 0.10) -> Optional[Tuple[int, int]]:
        """
        Finds the insertion sequence within the read sequence using edlib.
        Uses infix mode ('HW') to find the best match for the entire ins_seq within the read.
        Returns (start_pos_in_read, end_pos_in_read) (exclusive end) if found within threshold, else None.
        """
        if not ins_seq or not read.sequence:
            return None

        # Remove base 'N' from insertion sequence if present, as edlib might handle it differently
        # Depending on alignment needs, more sophisticated handling might be required.
        query_seq = ins_seq.replace('N', '').replace('n','')
        if not query_seq: # Handle case where insertion is all Ns
            return None

        read_seq = read.sequence

        try:
            result = edlib.align(query_seq, read_seq, mode="HW", task="locations")
        except ValueError:
            # Edlib can raise ValueError if sequences are empty after N removal, etc.
            return None

        edit_distance = result.get('editDistance', -1)
        locations = result.get('locations')
        alignment_length = result.get('alignmentLength')

        if edit_distance == -1 or not locations:
            # No alignment found
            return None

        # Check if edit distance is within the allowed threshold
        max_dist = int(len(query_seq) * max_edit_distance_ratio)
        if edit_distance > max_dist:
            return None

        # Edlib HW mode returns the end position (inclusive) of the query match in the target
        # We only expect one location for HW mode
        end_pos_inclusive = locations[0][1]
        # Approximate start position based on query length (adjusts for internal indels)
        # Note: This start position might be slightly off if there are many indels, but good for flanking.
        start_pos = end_pos_inclusive - len(query_seq) + 1 

        # Return 0-based, exclusive end coordinates
        return start_pos, end_pos_inclusive + 1, alignment_length
    
    def extract_context(self, variant: Variant) -> BreakpointContext:
        """Extract sequence context around variant breakpoint."""
        # For insertions, the breakpoint is at the insertion point
        # For deletions, it's at the start and end of the deletion
        # For other SV types, we need to handle them appropriately
        # TODO: actually handle other SV types
        if variant.sv_type == SVType.INS:
            start = max(0, variant.position - self.context_size)
            end = variant.position + self.context_size
            ref_seq = self._get_reference_sequence(variant.chrom, start, end)
            left_context = ref_seq[:self.context_size]
            right_context = ref_seq[-self.context_size:]
            breakpoint_pos = variant.position
        else:
            raise NotImplementedError(f"SV type {variant.sv_type} not supported (yet).") # TODO: idk this shouldn't be a problem but hold for now
            
        return BreakpointContext(
            left_context=left_context,
            right_context=right_context,
            breakpoint_pos=breakpoint_pos
        )
    
    def _longest_common_substring(self, left: str, right: str) -> Optional[str]:
        """Find longest common substring between two strings.
        
        Args:
            left: First sequence (left context sequence or window)
            right: Second sequence (right context sequence or window)

        Returns:
            A tuple of (longest_common_substring, length, left_start, right_start)
            Note: left_start and right_start are the positions of the longest common substring in the provided left and right sequences, respectively.
            These are *relative* positions and should be explicitly converted to absolute positions separately if/when needed.
        """
        max_substring = ""
        left_start, right_start = 0, 0
        max_len = 0
        for i in range(len(left)):
            for j in range(len(right)):
                k = 0
                while (i + k < len(left) and 
                    j + k < len(right) and 
                    left[i + k] == right[j + k]):
                    k += 1
                if k > max_len:
                    max_len = k
                    max_substring = left[i:i + k]
                    left_start, right_start = i, j

        return (max_substring, len(max_substring), left_start, right_start)
    
    def detect_microhomology(self, context: BreakpointContext, min_len: int = 4, max_len: int = 50) -> Optional[str]:
        """
        Detect microhomology at breakpoint by finding longest common substring.

        Args:
            context: BreakpointContext object containing left and right context sequences
            min_len: Minimum length of microhomology to detect
            max_len: Maximum length of microhomology to detect

        Returns:
            A tuple of (True, microhomology_sequence, microhomology_length, (left_microhomology_start, right_microhomology_start)) if microhomology is found.
            Else (False, 0, 0, (0, 0)).
        """
        left = context.left_context
        right = context.right_context
        
        LCS, LCS_len, left_start, right_start = self._longest_common_substring(left, right)

        if LCS and max_len >= LCS_len >= min_len:
            # previously returned absolute coordinates, switch to relative
            # return (True, LCS, LCS_len, (context.breakpoint_pos - self.context_size + left_start, context.breakpoint_pos + right_start))
            return (True, LCS, LCS_len, (-self.context_size + left_start, right_start))
        else:
            return (False, 0, 0, (0, 0))
    
    def detect_tsd(self, variant_analysis: VariantAnalysis, min_len: int = 5, max_len: int = 15) -> Optional[Tuple[str, int]]:
        """
        Detect target site duplication (TSD) by checking two possibilities:
        1. TSD is included at the start of the called insertion sequence.
        2. TSD flanks the called insertion sequence externally.
        Returns the first (longest) TSD found in any supporting read.

        Args:
            variant_analysis: The VariantAnalysis object containing the variant and reads.

        Returns:
            A tuple of (True, tsd_sequence, tsd_len, (left_tsd_start, right_tsd_start)) if a TSD is found.
            left_tsd_start and right_tsd_start are the positions of the TSD relative to the breakpoint.
            Else (False, 0, 0, (0, 0)).
        """
        variant = variant_analysis.variant
        if variant.sv_type != SVType.INS or not variant_analysis.support_reads or not variant_analysis.support_reads.reads:
            return None

        for read in variant_analysis.support_reads.reads:
            coords = None
            cigar_insertion = self._find_insertion_in_read(read, variant.sv_length)
            if cigar_insertion:
                coords = (cigar_insertion[0], cigar_insertion[1]) # [start, end)
            else: # Fallback to sequence alignment
                align_insertion = self._find_insertion_by_sequence(read, variant.alt)
                if align_insertion:
                    coords = (align_insertion[0], align_insertion[1]) # [start, end)

            if coords is None:
                continue # Cannot determine insertion boundaries in this read

            ins_start_in_read, ins_end_in_read = coords
            read_seq = read.sequence
            read_len = len(read_seq)

            lw_start, lw_end = max(0, ins_start_in_read - self.context_size), ins_start_in_read + self.context_size
            left_window = read_seq[lw_start:lw_end]
            rw_start, rw_end = ins_end_in_read - self.context_size, min(ins_end_in_read + self.context_size, read_len)
            right_window = read_seq[rw_start:rw_end]

            tsd_seq, tsd_len, tsd_left_start, tsd_right_start = self._longest_common_substring(left_window, right_window)
            if tsd_len >= min_len and tsd_len <= max_len:
                # distance relative to breakpoint
                left_tsd_start = tsd_left_start - self.context_size
                right_tsd_start = tsd_right_start - self.context_size

                return (True, tsd_seq, tsd_len, (left_tsd_start, right_tsd_start))
        return (False, 0, 0, (0, 0))
    
    def detect_poly_a_tail(self, variant: Variant) -> Optional[str]:
        """Detect poly-A tail in insertion sequence."""
        if variant.sv_type != SVType.INS:
            print('test')
            return None

        alt_seq = variant.alt
        seed_pattern = re.compile(r'A{10}')
        match = seed_pattern.search(alt_seq)
        
        if match:
            start = match.start()
            end = match.end()
            
            while end < len(alt_seq):
                if alt_seq[end:end+2] == 'AA':
                    end += 2
                elif alt_seq[end] == 'A':
                    end += 1
                    break
                else:
                    break
            
            return (True, len(alt_seq[start:end]), (start, end))
        else:
            # Try searching for poly-T stretches
            seed_pattern_t = re.compile(r'T{10}')
            match_t = seed_pattern_t.search(alt_seq)
            
            if match_t:
                start = match_t.start()
                end = match_t.end()
                
                while end < len(alt_seq):
                    if alt_seq[end:end+2] == 'TT':
                        end += 2
                    elif alt_seq[end] == 'T':
                        end += 1
                        break
                    else:
                        break
                
                return (True, len(alt_seq[start:end]), (start, end))
            
            return (False, 0, (0, 0))
    
    def analyze_breakpoint(self, variant_analysis: VariantAnalysis) -> BreakpointContext:
        """Perform comprehensive breakpoint analysis."""
        context = self.extract_context(variant_analysis.variant)
        
        context.microhomology = self.detect_microhomology(context)
        context.poly_a_tail = self.detect_poly_a_tail(variant_analysis.variant)
        context.tsd = self.detect_tsd(variant_analysis)
            
        return context 