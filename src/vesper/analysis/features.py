from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Protocol, Sequence, Tuple

import edlib

from ..models.contexts import BreakpointReadContext, ExtractedSupportRead, VariantContext
from ..models.variants import SVType


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
        if context.left_seq == "." or context.right_seq == "." or context.ins_start is None:
            return self._empty()

        left_window = context.left_seq
        right_window = context.right_seq
        best = None

        max_len = min(self.max_len, len(left_window), len(right_window))
        if max_len < self.min_len:
            return self._empty()

        for tsd_len in range(max_len, self.min_len - 1, -1):
            candidate = self._find_match(left_window, right_window, tsd_len)
            if candidate:
                best = candidate
                break

        if not best:
            return self._empty()

        sequence, left_start, right_start, edit_distance = best
        if _is_homopolymer(sequence):
            return self._empty()

        left_offset = left_start - len(left_window)
        right_offset = right_start

        return {
            "tsd_present": True,
            "tsd_sequence": sequence,
            "tsd_length": len(sequence),
            "tsd_left_offset": left_offset,
            "tsd_right_offset": right_offset,
            "tsd_edit_distance": edit_distance,
        }

    def _find_match(
        self,
        left_window: str,
        right_window: str,
        length: int,
    ) -> Optional[Tuple[str, int, int, int]]:
        if len(right_window) < length or len(left_window) < length:
            return None

        max_start = len(right_window) - length
        for right_start in range(max_start + 1):
            substring_right = right_window[right_start : right_start + length]
            try:
                result = edlib.align(
                    substring_right,
                    left_window,
                    mode="HW",
                    task="locations",
                    k=self.max_edit_distance if self.max_edit_distance >= 0 else -1,
                )
            except ValueError:
                continue

            edit_distance = result.get("editDistance", -1)
            if edit_distance == -1:
                continue

            if not result.get("locations"):
                continue

            left_start = result["locations"][0][0]
            if edit_distance > self.max_edit_distance:
                continue

            return substring_right, left_start, right_start, edit_distance

        return None

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
        if context.insert_seq == ".":
            return self._empty()

        strand = context.metadata.get("repeat_strand") or variant_ctx.metadata.get("repeat_strand") or "+"
        seed_char = "T" if strand in {"-", "C"} else "A"
        sequence = context.insert_seq

        match = self._find_poly_tail(sequence, seed_char)
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
            end = start
            while end < length:
                if sequence[end] != seed_char:
                    impurity += 1
                count += 1
                impurity_ratio = impurity / count
                if impurity_ratio > self.max_impurity:
                    break
                end += 1
            if count >= self.min_total_length and impurity / count <= self.max_impurity:
                candidate = sequence[start:end]
                candidate_impurity = impurity / count
                if best_seq is None or (len(candidate), -candidate_impurity) > (len(best_seq), -best_impurity):
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
