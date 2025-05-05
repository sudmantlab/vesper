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
    
    def __init__(self, reference_fasta: str, context_size: int = 25, max_edit: int = 1):
        """
        Args:
            reference_fasta: Path to reference genome FASTA file (must be indexed, e.g., with samtools faidx)
            context_size: Number of bases to extract on each side of breakpoint
        """
        self.context_size = context_size
        self.max_edit = max_edit
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
    
    def detect_tsd(self, context: BreakpointContext, variant_analysis: VariantAnalysis, min_len: int = 7, max_len: int = 20, max_edit: int = 1) -> Tuple[bool, bool, str, int, Tuple[int, int]]:
        """
        Detect target site duplication (TSD) allowing for an edit distance of 1.
        Checks for potential TSDs flanking the insertion site in supporting reads.
        Returns the first (longest) TSD found in any supporting read within the edit distance threshold.

        Args:
            variant_analysis: The VariantAnalysis object containing the variant and reads.
            min_len: Minimum length of TSD to consider.
            max_len: Maximum length of TSD to consider.

        Returns:
            1) A tuple of (True, edit_distance, left_tsd_sequence, right_tsd_sequence, tsd_len, (left_tsd_start_rel, right_tsd_start_rel)) if a TSD is found.
            left_tsd_start_rel and right_tsd_start_rel are the positions of the TSD relative to the breakpoint context start.
            Else 3) (False, None, None, None, 0, (0, 0)) if no TSD is found.
        """
        variant = variant_analysis.variant
        if variant.sv_type != SVType.INS or not variant_analysis.support_reads or not variant_analysis.support_reads.reads:
            return (False, False, "", 0, (0, 0)) # Return type adjusted

        for read in variant_analysis.support_reads.reads:
            coords = None
            # Prioritize CIGAR-based insertion finding
            cigar_insertion = self._find_insertion_in_read(read, variant.sv_length)
            if cigar_insertion:
                coords = (cigar_insertion[0], cigar_insertion[1]) # [start, end)
            else: # Fallback to sequence alignment if CIGAR doesn't give a clear insertion
                align_insertion = self._find_insertion_by_sequence(read, variant.alt)
                if align_insertion:
                    coords = (align_insertion[0], align_insertion[1]) # [start, end)

            if coords is None:
                continue # Cannot determine insertion boundaries in this read

            ins_start_in_read, ins_end_in_read = coords
            ins_len = ins_end_in_read - ins_start_in_read # avoid bugs where SVLEN != read ins seq
            query_start, query_end = variant_analysis.repeatmasker_results[0].query_start, variant_analysis.repeatmasker_results[0].query_end
            read_seq = read.sequence
            read_len = len(read_seq)

            lw_start = max(0, ins_start_in_read + query_start - self.context_size) # start upstream of query start
            lw_end = lw_start + self.context_size 
            left_window = read_seq[lw_start:lw_end]

            if context.poly_a_tail[0]: # use poly-A tail as window start
                poly_a_start, poly_a_end = context.poly_a_tail[2][0], context.poly_a_tail[2][1]
                if variant_analysis.repeatmasker_results[0].strand == 'C': # orient poly-A relative to ins
                    poly_a_start, poly_a_end = ins_len - poly_a_end, ins_len - poly_a_start
                rw_start = ins_start_in_read + poly_a_end # start at poly-A tail end
                rw_end = min(read_len, max(ins_end_in_read + self.context_size, rw_start + self.context_size)) # end downstream of poly-A tail end
                right_window = read_seq[rw_start:rw_end]
            else: # use query end as window start
                rw_start = ins_start_in_read + query_end # start at query end
                rw_end =min(read_len, max(ins_end_in_read + self.context_size, rw_start + self.context_size)) # end downstream of query end
                right_window = read_seq[rw_start:rw_end]

            if not left_window or not right_window: # Skip if windows are empty
                continue

            best_exact_match = None
            best_inexact_match = None

            # --- Pass 1: Search for the longest EXACT match (k=0) ---
            for tsd_len in range(max_len, min_len - 1, -1):
                if best_exact_match: break # Already found the longest possible exact match
                # Iterate starting positions in right_window
                for j in range(len(right_window) - tsd_len + 1):
                    substring_right = right_window[j : j + tsd_len]
                    try:
                        # Align substring_right (query) against left_window (target)
                        result = edlib.align(substring_right, left_window, mode="HW", task="locations", k=0)
                    except ValueError:
                        continue

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
                                best_exact_match = (True, 0, left_tsd_sequence, substring_right, tsd_len, (left_tsd_start_rel, right_tsd_start_rel))
                                break # Found longest non-homopolymer exact match
            
            # If an exact match was found, return it immediately for this read
            if best_exact_match:
                return best_exact_match
            
            # --- Pass 2: Search for the longest INEXACT match (k=max_edit) --- 
            # This part only runs if no exact match was found in Pass 1
            for tsd_len in range(max_len, min_len - 1, -1):
                if best_inexact_match: break # Already found the longest possible inexact match
                for j in range(len(right_window) - tsd_len + 1):
                    substring_right = right_window[j : j + tsd_len]
                    try:
                        # Align substring_right (query) against left_window (target)
                        result = edlib.align(substring_right, left_window, mode="HW", task="locations", k=max_edit)
                    except ValueError:
                        continue 

                    edit_distance = result.get('editDistance', -1)
                    # Check if a valid *inexact* alignment was found (must be > 0 and <= max_edit)
                    if 0 < edit_distance <= max_edit:
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
                            
                            # # Store if neither sequence is a homopolymer
                            if not self._is_homopolymer(left_tsd_sequence, threshold = 1.0) and not self._is_homopolymer(substring_right, threshold = 0.8):
                                best_inexact_match = (True, edit_distance, left_tsd_sequence, substring_right, tsd_len, (left_tsd_start_rel, right_tsd_start_rel))
                                break # Found longest non-homopolymer inexact match

            # If an inexact match was found (and no exact one was), return it
            if best_inexact_match:
                return best_inexact_match

            # If loop finishes for THIS read without finding any suitable TSD, continue to next read

        # If loop finishes for ALL reads without finding a TSD:
        return (False, None, None, None, 0, (0, 0))

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
            
            return (False, None, (0, 0))
    
    def analyze_breakpoint(self, variant_analysis: VariantAnalysis) -> BreakpointContext:
        """Perform comprehensive breakpoint analysis."""
        context = self.extract_context(variant_analysis.variant)
        
        context.microhomology = self.detect_microhomology(context)
        context.poly_a_tail = self.detect_poly_a_tail(variant_analysis.variant)
        context.tsd = self.detect_tsd(context, variant_analysis, max_edit = self.max_edit)
        return context 