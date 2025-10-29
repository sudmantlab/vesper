from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, List, TYPE_CHECKING
import re
import edlib

from ..models.variants import SVType

if TYPE_CHECKING:
    from ..models.variants import VariantAnalysis
    from ..models.reads import AlignedRead

@dataclass
class BreakpointContext:
    """Represents the sequence context around a variant breakpoint."""
    read_name: Optional[str] = None  # Name of supporting read
    read_seq: Optional[str] = None  # Sequence of supporting read
    read_ins_coords: Optional[Tuple[int, int]] = None  # Insertion coordinates in supporting read
    tsd: Optional[Tuple[bool, str, str, int, Tuple[int, int]]] = None  # Target site duplication if present
    poly_a_tail: Optional[Tuple[bool, str, Tuple[int, int], int]] = None  # Poly-A tail if present

    support_reads: Optional[List[AlignedRead]] = None  # (Debug only) List of supporting reads

# pysam CIGAR operation codes: BAM_CMATCH 0, BAM_CINS 1, BAM_CDEL 2, BAM_CREF_SKIP 3,
# BAM_CSOFT_CLIP 4, BAM_CHARD_CLIP 5, BAM_CPAD 6, BAM_CEQUAL 7, BAM_CDIFF 8
CIGAR_CONSUMES_QUERY = {0, 1, 4, 7, 8} # M, I, S, =, X


class BreakpointAnalyzer:
    """Analyzes breakpoint contexts for structural variants."""
    
    def __init__(self, context_size: int = 25, debug: bool = False):
        """
        Args:
            context_size: Number of bases to extract on each side of breakpoint
        """
        self.context_size = context_size
        self.debug = debug # returns read sequences in BreakpointContext
        
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
        
    def _is_homopolymer(self, sequence, threshold=0.7):
        """
        Check if any single base makes up more than threshold% of the sequence.
        
        Args:
            sequence (str): DNA sequence to check
            threshold (float): Fraction threshold for homopolymer detection (default 0.7)
        Returns:
            bool: True if any base exceeds threshold fraction, False otherwise
        """
        if not sequence:
            return False
        
        seq_len = len(sequence)
        base_counts = {
            'A': sequence.upper().count('A') / seq_len,
            'C': sequence.upper().count('C') / seq_len, 
            'G': sequence.upper().count('G') / seq_len,
            'T': sequence.upper().count('T') / seq_len
        }

        for base, frac in base_counts.items():
            if frac >= threshold:
                stretch = base * int(frac * seq_len)
                if stretch in sequence:
                    return True
        return False
    
    def _get_supporting_read(self, context: BreakpointContext, variant_analysis: VariantAnalysis) -> Tuple[int, int]:
        """Find valid supporting read for analysis by searching for valid insertion coordinates in read."""
        variant, coords = variant_analysis.variant, None
        if variant.sv_type != SVType.INS:
            return coords
        
        for read in variant_analysis.support_reads.reads:
            read_name = read.name
            read_seq = read.sequence
            # Prioritize CIGAR-based insertion finding
            cigar_insertion = self._find_insertion_in_read(read, variant.sv_length)
            if cigar_insertion:
                coords = (cigar_insertion[0] - 1, cigar_insertion[1] - 1) # [start, end)
            else: # Fallback to sequence alignment if CIGAR doesn't give a clear insertion
                align_insertion = self._find_insertion_by_sequence(read, variant.alt)
                if align_insertion:
                    coords = (align_insertion[0] - 1, align_insertion[1] - 1) # [start, end)
        
        context.read_name = read_name
        context.read_seq = read_seq
        context.read_ins_coords = coords
        return context
    
    def detect_tsd(self, context: BreakpointContext, variant_analysis: VariantAnalysis, min_len: int = 7, max_len: int = 20) -> Tuple[bool, bool, str, int, Tuple[int, int]]:
        """
        Detect target site duplication (TSD) allowing for an edit distance of 1.
        Checks for potential TSDs flanking the insertion site in supporting reads.
        Returns the first (longest) TSD found in any supporting read within the edit distance threshold.

        Args:
            variant_analysis: The VariantAnalysis object containing the variant and reads.
            min_len: Minimum length of TSD to consider.
            max_len: Maximum length of TSD to consider.

        Returns:
            1) A tuple of (True, left_tsd_sequence, right_tsd_sequence, tsd_len, (left_tsd_start_rel, right_tsd_start_rel)) if a TSD is found.
            left_tsd_start_rel and right_tsd_start_rel are the positions of the TSD relative to the breakpoint context start.
            Else 3) (False, None, None, None, 0, (0, 0)) if no TSD is found.
        """
        variant = variant_analysis.variant
        if variant.sv_type != SVType.INS or not variant_analysis.support_reads or not variant_analysis.support_reads.reads:
            return (False, None, "", 0, (0, 0)) # Return type adjusted

        coords = context.read_ins_coords
        if coords is None:
            return (False, None, "", 0, (0, 0)) # Return type adjusted

        read_seq = context.read_seq
        read_len = len(read_seq)

        ins_start_in_read, ins_end_in_read = coords
        query_start, query_end = variant_analysis.repeatmasker_results[0].query_start, variant_analysis.repeatmasker_results[0].query_end

        lw_start = max(0, ins_start_in_read - self.context_size) # start upstream of ins start - in case query is downstream of ins (ex. minus strand insertion)
        lw_end = ins_start_in_read + query_start # end at query start (match does not include TSD)
        left_window = read_seq[lw_start:lw_end]
        
        # use query end as window start
        rw_start = ins_start_in_read + query_end
        rw_end = min(read_len, max(ins_end_in_read + self.context_size, rw_start + self.context_size)) # end downstream of query end
        right_window = read_seq[rw_start:rw_end]

        if not left_window or not right_window: # Skip if windows are empty
            return (False, None, "", 0, (0, 0)) # Return type adjusted

        best_exact_match = None

        for tsd_len in range(max_len, min_len - 1, -1):
            if best_exact_match:
                break  # Already found the longest possible exact match
            # Iterate starting positions in right_window
            for j in range(len(right_window) - tsd_len + 1):
                substring_right = right_window[j : j + tsd_len]
                try:
                    # Align substring_right (query) against left_window (target)
                    result = edlib.align(substring_right, left_window, mode="HW", task="locations", k=0)
                except ValueError:
                    return (False, None, "", 0, (0, 0)) # Return type adjusted

                edit_distance = result.get('editDistance', -1)
                if edit_distance == 0: # Found exact match
                    locations = result.get('locations')
                    if locations:
                        # Match position is in left_window
                        left_tsd_start_in_window = locations[0][0] 
                        # Query position is in right_window
                        right_tsd_start_in_window = j
                        
                        left_tsd_start_rel = (lw_start + left_tsd_start_in_window) - ins_start_in_read
                        right_tsd_start_rel = (rw_start + right_tsd_start_in_window) - ins_end_in_read
                        
                        # Extract the matched sequence from the left window
                        target_start, target_end = locations[0]
                        left_tsd_sequence = left_window[target_start : target_end + 1]

                        # Store if neither sequence is a homopolymer
                        if not self._is_homopolymer(left_tsd_sequence, threshold = 1.0) and not self._is_homopolymer(substring_right, threshold = 0.8):
                            best_exact_match = (True, left_tsd_sequence, substring_right, tsd_len, (left_tsd_start_rel, right_tsd_start_rel))
                            break # Found longest non-homopolymer exact match
        
        # If an exact match was found, return it immediately for this read
        if best_exact_match:
            return best_exact_match

        return (False, None, "", 0, (0, 0)) # Return type adjusted

    def _gap_extend(self, seed_char: str, seed_len: int, sequence: str, min_total_length: int, max_impurity: float, left_bound: int, right_bound: int, max_local_impurity: int = 4) -> Optional[Tuple[bool, str, Tuple[int, int], int]]:
        """Gap extend a seed pattern in a sequence."""
        seed_pattern = re.compile(seed_char * seed_len)
        match = seed_pattern.search(sequence)

         # TODO: Fix poly-A tail detection to use full read seq rather than insertion seq to account for imperfect breakpoint calling

        if match:
            start = match.start()
            end = match.end()

            # Extend backwards
            pos = max(start, left_bound + 1)
            impure_count = 0
            impure_ratio = 0
            local_impure = 0
            while pos > left_bound and (impure_count/(end - pos + 1)) <= max_impurity and local_impure <= max_local_impurity:
                if sequence[pos] == seed_char:
                    start = pos
                    local_impure = 0
                else:
                    local_impure += 1
                    impure_count += 1
                    impure_ratio = impure_count/(end - pos + 1)
                    if impure_ratio > max_impurity:
                        break
                    start = pos  # Include the gap
                pos -= 1

            # Extend forwards
            pos = end
            impure_count = 0
            impure_ratio = 0
            local_impure = 0
            while pos < right_bound and (impure_count/(pos - start + 1)) <= max_impurity and local_impure <= max_local_impurity:
                if sequence[pos] == seed_char:
                    end = pos
                    local_impure = 0
                else:
                    local_impure += 1
                    impure_count += 1
                    impure_ratio = impure_count/(pos - start + 1)
                    if impure_ratio > max_impurity:
                        break
                    end = pos
                pos += 1
        
            if end - start >= min_total_length:
                # adjust end so full tail included with end-exclusive range
                end += 1
                matched = sequence[start:end]

                # Shrink ends to exclude non-matching (impure) bases
                while matched[0] != seed_char:
                    matched = matched[1:]
                    start += 1
                while matched[-1] != seed_char:
                    matched = matched[:-1]
                    end -= 1
                
                # Calculate final impurity ratio
                impure_count = len(matched) - matched.count(seed_char)
                impure_ratio = impure_count/len(matched)

                # # adjust end so full tail included with end-exclusive range
                # end += 1
                matched = sequence[start:end]

                return (True, matched, (start, end), len(matched), impure_ratio)
        else:
            return (False, None, (0, 0), 0, 0)

    def detect_poly_a_tail(self, context: BreakpointContext, variant_analysis: VariantAnalysis, min_total_length: int = 10, max_impurity: float = 0.20) -> Optional[Tuple[bool, str, Tuple[int, int], int]]:
        """Detect poly-A tail in insertion sequence."""
        variant = variant_analysis.variant
        if variant.sv_type != SVType.INS:
            return (False, None, (0, 0), 0)
        
        alt_seq = variant.alt
        tsd_len, left_tsd_start = context.tsd[3], context.tsd[4][0]
        left_tsd_end = left_tsd_start + tsd_len
        original_strand = variant_analysis.repeatmasker_results[0].strand

        # take into account TSD positions
        if original_strand == 'C': # prevent poly-A tail from extending into left TSD by setting bounds
            left_bound, right_bound = left_tsd_end, len(alt_seq)
            tail = self._gap_extend('T', 10, alt_seq, min_total_length, max_impurity, left_bound, right_bound)  
        else:
            left_bound, right_bound = 0, len(alt_seq)
            tail = self._gap_extend('A', 10, alt_seq, min_total_length, max_impurity, left_bound, right_bound)
          
        return tail
    
    def analyze_breakpoint(self, variant_analysis: VariantAnalysis) -> BreakpointContext:
        """Perform comprehensive breakpoint analysis."""
        context = BreakpointContext()
        context = self._get_supporting_read(context, variant_analysis)
        context.tsd = self.detect_tsd(context, variant_analysis)
        context.poly_a_tail = self.detect_poly_a_tail(context, variant_analysis)
        if self.debug:
            context.support_reads = variant_analysis.support_reads

        return context 
