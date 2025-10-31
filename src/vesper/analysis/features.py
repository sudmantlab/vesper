from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Protocol, Sequence, Tuple

import edlib
import re

from ..models.contexts import BreakpointReadContext, ExtractedSupportRead, VariantContext
from ..models.variants import SVType

LEGACY_CONTEXT_DEFAULT = 30
LEGACY_MIN_TSD_LEN = 7
LEGACY_MAX_TSD_LEN = 20
LEGACY_SEED_LEN = 10
LEGACY_MAX_LOCAL_IMPURITY = 4
CIGAR_CONSUMES_QUERY = {0, 1, 4, 7, 8}


def _get_support_read(variant_ctx: VariantContext, support_rank: int) -> Optional[ExtractedSupportRead]:
    try:
        return variant_ctx.support_reads[support_rank]
    except (AttributeError, IndexError):
        return None


def _legacy_locate_insertion(read: ExtractedSupportRead, expected_len: int) -> Optional[Tuple[int, int]]:
    if not read or not read.cigartuples:
        return None

    tolerance = max(10, int(expected_len * 0.1)) if expected_len else 10
    best: Optional[Tuple[int, int]] = None
    query_idx = 0

    for op, length in read.cigartuples:
        if op == 1:
            within_tolerance = (
                not expected_len or abs(length - expected_len) <= tolerance
            )
            if within_tolerance:
                if best is None or length > best[0]:
                    best = (length, query_idx)
        if op in CIGAR_CONSUMES_QUERY:
            query_idx += length

    if best is None:
        return None

    ins_len, ins_start = best
    ins_end = ins_start + ins_len
    return max(0, ins_start - 1), ins_end


def _legacy_extract_windows(
    sequence: str,
    ins_start: int,
    ins_end: int,
    repeat_start: Optional[int],
    repeat_end: Optional[int],
    context_size: int,
) -> Tuple[str, str, int, int]:
    left_start = max(0, ins_start - context_size)
    if isinstance(repeat_start, int) and repeat_start >= 0:
        left_end = min(len(sequence), ins_start + repeat_start)
    else:
        left_end = ins_start

    if isinstance(repeat_end, int) and repeat_end >= 0:
        right_start = min(len(sequence), ins_start + repeat_end)
    else:
        right_start = ins_end

    right_end = min(len(sequence), max(ins_end + context_size, right_start + context_size))

    left_window = sequence[left_start:left_end]
    right_window = sequence[right_start:right_end]
    return left_window, right_window, left_start, right_start


def _legacy_detect_tsd(
    variant_ctx: VariantContext,
    context: BreakpointReadContext,
    min_len: int,
    max_len: int,
) -> Optional[Dict[str, Any]]:
    support_read = _get_support_read(variant_ctx, context.support_rank)
    if not support_read or support_read.sequence is None:
        return None

    coords = _legacy_locate_insertion(support_read, variant_ctx.variant.sv_length)
    if not coords:
        return None

    ins_start, ins_end = coords

    context_size = context.metadata.get("context_size", LEGACY_CONTEXT_DEFAULT)
    repeat_start = variant_ctx.metadata.get("repeat_query_start")
    repeat_end = variant_ctx.metadata.get("repeat_query_end")

    (
        left_window,
        right_window,
        left_window_start,
        right_window_start,
    ) = _legacy_extract_windows(
        support_read.sequence,
        ins_start,
        ins_end,
        repeat_start,
        repeat_end,
        context_size,
    )

    if not left_window or not right_window:
        return None

    max_candidate_len = min(max_len, len(left_window), len(right_window))
    if max_candidate_len < min_len:
        return None

    best = None
    for tsd_len in range(max_candidate_len, min_len - 1, -1):
        for j in range(len(right_window) - tsd_len + 1):
            substring_right = right_window[j : j + tsd_len]
            try:
                result = edlib.align(
                    substring_right,
                    left_window,
                    mode="HW",
                    task="locations",
                    k=0,
                )
            except ValueError:
                continue

            if result.get("editDistance", -1) != 0:
                continue

            locations = result.get("locations")
            if not locations:
                continue

            left_start = locations[0][0]
            left_seq = left_window[left_start : left_start + tsd_len]

            if _is_homopolymer(left_seq, threshold=1.0) or _is_homopolymer(substring_right, threshold=0.8):
                continue

            best = (
                left_seq,
                tsd_len,
                left_start,
                j,
                left_window,
                right_window,
                ins_start,
                ins_end,
                repeat_start,
                repeat_end,
            )
            break
        if best:
            break

    if not best:
        return None

    (
        sequence,
        tsd_len,
        left_start,
        right_start,
        left_window,
        right_window,
        ins_start,
        ins_end,
        repeat_start,
        repeat_end,
    ) = best

    left_tsd_start_rel = (left_window_start + left_start) - ins_start
    right_tsd_start_rel = (right_window_start + right_start) - ins_end

    return {
        "tsd_present": True,
        "tsd_sequence": sequence,
        "tsd_length": tsd_len,
        "tsd_left_offset": left_tsd_start_rel,
        "tsd_right_offset": right_tsd_start_rel,
        "tsd_edit_distance": 0,
    }


def _legacy_gap_extend(
    seed_char: str,
    seed_len: int,
    sequence: str,
    min_total_length: int,
    max_impurity: float,
    left_bound: int,
    right_bound: int,
) -> Optional[Tuple[str, float]]:
    seed_pattern = re.compile(seed_char * seed_len)
    match = seed_pattern.search(sequence)
    if not match:
        return None

    start = match.start()
    end = match.end()

    pos = max(start, left_bound + 1)
    impure_count = 0
    local_impure = 0
    while pos > left_bound and (impure_count / (end - pos + 1)) <= max_impurity and local_impure <= LEGACY_MAX_LOCAL_IMPURITY:
        if sequence[pos] == seed_char:
            start = pos
            local_impure = 0
        else:
            local_impure += 1
            impure_count += 1
            impurity_ratio = impure_count / (end - pos + 1)
            if impurity_ratio > max_impurity:
                break
            start = pos
        pos -= 1

    pos = end
    impure_count = 0
    local_impure = 0
    while pos < right_bound and (impure_count / (pos - start + 1)) <= max_impurity and local_impure <= LEGACY_MAX_LOCAL_IMPURITY:
        if sequence[pos] == seed_char:
            end = pos
            local_impure = 0
        else:
            local_impure += 1
            impure_count += 1
            impurity_ratio = impure_count / (pos - start + 1)
            if impurity_ratio > max_impurity:
                break
            end = pos
        pos += 1

    if end - start < min_total_length:
        return None

    end += 1
    matched = sequence[start:end]

    while matched and matched[0] != seed_char:
        matched = matched[1:]
        start += 1
    while matched and matched[-1] != seed_char:
        matched = matched[:-1]
        end -= 1

    if not matched or len(matched) < min_total_length:
        return None

    impure_count = len(matched) - matched.count(seed_char)
    impurity_ratio = impure_count / len(matched) if matched else 0.0

    return matched, impurity_ratio


class FeatureDetector(Protocol):
    name: str
    output_columns: Sequence[str]

    def analyze(self, context: BreakpointReadContext, variant_ctx: VariantContext) -> Mapping[str, Any]:
        ...


BASE_CONTEXT_COLUMNS = [
    "sample_id",
    "variant_id",
    "read_id",
    "chrom",
    "position",
    "sv_type",
    "sv_length",
    "strand",
    "support_rank",
    "ins_start",
    "ins_end",
    "left_seq",
    "insert_seq",
    "right_seq",
    "mapq",
    "edit_distance",
    "is_primary",
]


@dataclass
class FeatureSet:
    name: str
    detectors: Sequence[FeatureDetector]

    def columns(self, context_columns: Sequence[str] | None = None) -> List[str]:
        base = list(context_columns or BASE_CONTEXT_COLUMNS)
        seen = set(base)
        for detector in self.detectors:
            for column in detector.output_columns:
                if column not in seen:
                    base.append(column)
                    seen.add(column)
        return base

    def analyze(self, context: BreakpointReadContext, variant_ctx: VariantContext) -> Dict[str, Any]:
        row: MutableMapping[str, Any] = context.to_dict()
        for detector in self.detectors:
            row.update(detector.analyze(context, variant_ctx))
        return dict(row)


def _is_homopolymer(sequence: str, threshold: float = 0.8) -> bool:
    if not sequence:
        return False
    seq_len = len(sequence)
    for base in "ACGT":
        if sequence.count(base) / seq_len >= threshold:
            return True
    return False


class TsdDetector:
    name = "tsd"
    output_columns = (
        "tsd_present",
        "tsd_sequence",
        "tsd_length",
        "tsd_left_offset",
        "tsd_right_offset",
        "tsd_edit_distance",
    )

    def __init__(
        self,
        min_len: int = 7,
        max_len: int = 20,
        max_edit_distance: int = 0,
    ) -> None:
        self.min_len = min_len
        self.max_len = max_len
        self.max_edit_distance = max_edit_distance

    def analyze(self, context: BreakpointReadContext, variant_ctx: VariantContext) -> Mapping[str, Any]:
        if variant_ctx.variant.sv_type != SVType.INS:
            return self._empty()
        support_read = _get_support_read(variant_ctx, context.support_rank)
        if not support_read or support_read.sequence is None:
            return self._empty()

        tsd_result = _legacy_detect_tsd(
            variant_ctx,
            context,
            min_len=self.min_len,
            max_len=self.max_len,
        )
        if not tsd_result:
            return self._empty()

        return {
            "tsd_present": tsd_result["tsd_present"],
            "tsd_sequence": tsd_result["tsd_sequence"],
            "tsd_length": tsd_result["tsd_length"],
            "tsd_left_offset": tsd_result["tsd_left_offset"],
            "tsd_right_offset": tsd_result["tsd_right_offset"],
            "tsd_edit_distance": tsd_result["tsd_edit_distance"],
        }

    def _empty(self) -> Mapping[str, Any]:
        return {
            "tsd_present": False,
            "tsd_sequence": ".",
            "tsd_length": 0,
            "tsd_left_offset": 0,
            "tsd_right_offset": 0,
            "tsd_edit_distance": 0,
        }


class PolyATailDetector:
    name = "poly_a"
    output_columns = (
        "poly_a_present",
        "poly_a_sequence",
        "poly_a_length",
        "poly_a_impurity",
    )

    def __init__(self, min_total_length: int = 10, max_impurity: float = 0.2) -> None:
        self.min_total_length = min_total_length
        self.max_impurity = max_impurity

    def analyze(self, context: BreakpointReadContext, variant_ctx: VariantContext) -> Mapping[str, Any]:
        if variant_ctx.variant.sv_type != SVType.INS:
            return self._empty()
        alt_seq = variant_ctx.variant.alt
        if not alt_seq or alt_seq == ".":
            return self._empty()

        strand = context.metadata.get("repeat_strand") or variant_ctx.metadata.get("repeat_strand") or "+"
        seed_char = "T" if strand in {"-", "C"} else "A"

        tsd_result = _legacy_detect_tsd(
            variant_ctx,
            context,
            min_len=LEGACY_MIN_TSD_LEN,
            max_len=LEGACY_MAX_TSD_LEN,
        )

        if tsd_result and tsd_result["tsd_present"]:
            left_tsd_start = tsd_result["tsd_left_offset"]
            tsd_len = tsd_result["tsd_length"]
        else:
            left_tsd_start = 0
            tsd_len = 0

        if strand in {"-", "C"}:
            left_bound = left_tsd_start + tsd_len
        else:
            left_bound = 0
        right_bound = len(alt_seq)

        match = _legacy_gap_extend(
            seed_char,
            LEGACY_SEED_LEN,
            alt_seq,
            self.min_total_length,
            self.max_impurity,
            left_bound,
            right_bound,
        )
        if not match:
            return self._empty()

        matched, impurity = match
        return {
            "poly_a_present": True,
            "poly_a_sequence": matched,
            "poly_a_length": len(matched),
            "poly_a_impurity": impurity,
        }

    def _find_poly_tail(self, sequence: str, seed_char: str) -> Optional[Tuple[str, float]]:
        best_seq = None
        best_impurity = 1.0
        length = len(sequence)

        for start in range(length):
            impurity = 0
            count = 0
            for end in range(start, length):
                if sequence[end] != seed_char:
                    impurity += 1
                count += 1
                impurity_ratio = impurity / count
                if impurity_ratio > self.max_impurity:
                    break
                if count >= self.min_total_length:
                    candidate = sequence[start : end + 1]
                    candidate_impurity = impurity_ratio
                    if best_seq is None or (len(candidate), -candidate_impurity) > (
                        len(best_seq), -best_impurity
                    ):
                        best_seq = candidate
                        best_impurity = candidate_impurity

        if best_seq is None:
            return None
        return best_seq, best_impurity

    def _empty(self) -> Mapping[str, Any]:
        return {
            "poly_a_present": False,
            "poly_a_sequence": ".",
            "poly_a_length": 0,
            "poly_a_impurity": 0.0,
        }


BASE_FEATURE_SET = FeatureSet(name="base", detectors=())
ALU_FEATURE_SET = FeatureSet(name="alu", detectors=(TsdDetector(), PolyATailDetector()))


FEATURE_SETS = {
    feature_set.name: feature_set
    for feature_set in [BASE_FEATURE_SET, ALU_FEATURE_SET]
}


def get_feature_set(name: str) -> FeatureSet:
    try:
        return FEATURE_SETS[name]
    except KeyError as exc:
        raise KeyError(f"Unknown feature set '{name}'") from exc


def available_feature_sets() -> Iterable[str]:
    return FEATURE_SETS.keys()
