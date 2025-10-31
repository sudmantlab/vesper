from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Optional, Tuple

from ..models.contexts import BreakpointReadContext, ExtractedSupportRead, VariantContext
from ..models.variants import SVType

CIGAR_CONSUMES_QUERY = {0, 1, 4, 7, 8}


class ContextExtractionError(Exception):
    pass


class MissingSequenceError(ContextExtractionError):
    pass


class InsertionNotFound(ContextExtractionError):
    pass


@dataclass
class ContextExtractor:
    context_size: int = 25
    length_tolerance: float = 0.1
    strict: bool = False

    def iter_read_contexts(self, variant_ctx: VariantContext) -> Iterator[BreakpointReadContext]:
        if variant_ctx.variant.sv_type not in {SVType.INS, SVType.DEL}:  # scaffolding for future buildout on DEL types
            return iter(())
        return self._iter_breakpoints(variant_ctx)

    def _iter_breakpoints(self, variant_ctx: VariantContext) -> Iterator[BreakpointReadContext]:
        for idx, read in enumerate(variant_ctx.support_reads):
            try:
                yield self._analyze_read(variant_ctx, read, idx)
            except ContextExtractionError as exc:
                if isinstance(exc, MissingSequenceError):
                    raise
                if self.strict:
                    raise
                yield self._placeholder_context(variant_ctx, read, idx, reason="error")

    def _analyze_read(
        self,
        variant_ctx: VariantContext,
        read: ExtractedSupportRead,
        support_rank: int,
    ) -> BreakpointReadContext:
        if read.sequence is None:
            raise MissingSequenceError(f"Sequence missing for read {read.read_id}")

        coords = self._locate_event(read, variant_ctx.variant.sv_type, variant_ctx.variant.sv_length)
        if coords is None:
            if self.strict:
                raise InsertionNotFound(f"Event not found in read {read.read_id}")
            return self._placeholder_context(variant_ctx, read, support_rank, reason="missing_event")

        ins_start, ins_end, event_length = coords
        legacy_ins_start = max(0, ins_start - 1)

        repeat_start_rel = variant_ctx.metadata.get("repeat_query_start")
        repeat_end_rel = variant_ctx.metadata.get("repeat_query_end")

        if variant_ctx.variant.sv_type == SVType.INS:  # may also work for DUP, sniffles2 offers conversion of DUP to DEL
            left_start = max(0, legacy_ins_start - self.context_size)
            if isinstance(repeat_start_rel, int) and repeat_start_rel >= 0:
                left_end = min(len(read.sequence), legacy_ins_start + repeat_start_rel)
            else:
                left_end = legacy_ins_start

            if isinstance(repeat_end_rel, int) and repeat_end_rel >= 0:
                right_start_candidate = min(len(read.sequence), legacy_ins_start + repeat_end_rel)
                if right_start_candidate > ins_end:
                    right_start = ins_end
                else:
                    right_start = right_start_candidate
            else:
                right_start = ins_end

            right_extension_base = ins_end + self.context_size
            right_end = min(len(read.sequence), max(right_extension_base, right_start + self.context_size))

            left_seq = read.sequence[left_start:left_end]
            insert_seq = read.sequence[ins_start:ins_end]
            right_seq = read.sequence[right_start:right_end]
        else:  # scaffolding for future buildout on DEL types
            left_start = max(0, legacy_ins_start - self.context_size)
            left_end = legacy_ins_start
            right_start = legacy_ins_start
            right_end = min(len(read.sequence), legacy_ins_start + self.context_size)
            left_seq = read.sequence[left_start:left_end]
            insert_seq = "."
            right_seq = read.sequence[right_start:right_end]

        metadata = dict(variant_ctx.metadata)
        metadata["context_status"] = "ok"
        metadata["event_length"] = event_length
        metadata.setdefault("context_size", self.context_size)

        return BreakpointReadContext(
            sample_id=variant_ctx.sample_id,
            variant_id=variant_ctx.variant.ID,
            read_id=read.read_id,
            chrom=variant_ctx.variant.chrom,
            position=variant_ctx.variant.position,
            sv_type=variant_ctx.variant.sv_type,
            sv_length=variant_ctx.variant.sv_length,
            strand=read.strand,
            support_rank=support_rank,
            ins_start=legacy_ins_start,
            ins_end=ins_end,
            left_seq=left_seq if left_seq else ".",
            insert_seq=insert_seq if insert_seq else ".",
            right_seq=right_seq if right_seq else ".",
            mapq=read.mapq,
            edit_distance=read.edit_distance,
            is_primary=not (read.is_secondary or read.is_supplementary),
            metadata=metadata,
        )

    def _placeholder_context(
        self,
        variant_ctx: VariantContext,
        read: ExtractedSupportRead,
        support_rank: int,
        reason: str,
        ) -> BreakpointReadContext:
        metadata = dict(variant_ctx.metadata)
        metadata["context_status"] = reason
        metadata.setdefault("event_length", variant_ctx.variant.sv_length)
        metadata.setdefault("context_size", self.context_size)
        return BreakpointReadContext(
            sample_id=variant_ctx.sample_id,
            variant_id=variant_ctx.variant.ID,
            read_id=read.read_id,
            chrom=variant_ctx.variant.chrom,
            position=variant_ctx.variant.position,
            sv_type=variant_ctx.variant.sv_type,
            sv_length=variant_ctx.variant.sv_length,
            strand=read.strand,
            support_rank=support_rank,
            ins_start=None,
            ins_end=None,
            left_seq=".",
            insert_seq=".",
            right_seq=".",
            mapq=read.mapq,
            edit_distance=read.edit_distance,
            is_primary=not (read.is_secondary or read.is_supplementary),
            metadata=metadata,
        )

    def _locate_event(
        self,
        read: ExtractedSupportRead,
        sv_type: SVType,
        expected_length: int,
    ) -> Optional[Tuple[int, int, int]]:
        return self._find_event_in_read(read, sv_type, expected_length)

    def _find_event_in_read(
        self,
        read: ExtractedSupportRead,
        sv_type: SVType,
        expected_len: int,
    ) -> Optional[Tuple[int, int, int]]:
        if not read.cigartuples:
            return None

        if sv_type == SVType.INS:
            target_op = 1
        elif sv_type == SVType.DEL:
            target_op = 2
        else:
            return None

        best: Optional[Tuple[int, int, int]] = None
        query_idx = 0
        tolerance = max(10, int(expected_len * self.length_tolerance)) if expected_len else 10

        for op, length in read.cigartuples:
            if op == target_op:
                skip = False
                if expected_len and abs(length - expected_len) > tolerance and best is not None:
                    skip = True

                if not skip:
                    if sv_type == SVType.INS:
                        candidate_start = query_idx
                        candidate_end = query_idx + length
                    else:
                        candidate_start = query_idx
                        candidate_end = query_idx

                    candidate = (length, candidate_start, candidate_end)

                    if best is None:
                        best = candidate
                    else:
                        if expected_len:
                            if abs(length - expected_len) < abs(best[0] - expected_len):
                                best = candidate
                        elif length > best[0]:
                            best = candidate
            if op in CIGAR_CONSUMES_QUERY:
                query_idx += length

        if best:
            length, start, end = best
            return start, end, length

        return None
