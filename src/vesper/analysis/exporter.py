from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, MutableMapping, Optional, Sequence

import pandas as pd

from .context_extractor import ContextExtractor
from .features import FeatureSet, BASE_CONTEXT_COLUMNS
from ..models.contexts import VariantContext


def iter_feature_rows(
    variant_contexts: Iterable[VariantContext],
    extractor: ContextExtractor,
    feature_set: FeatureSet,
):
    for variant_ctx in variant_contexts:
        for read_ctx in extractor.iter_read_contexts(variant_ctx):
            yield feature_set.analyze(read_ctx, variant_ctx)


def export_breakpoint_tsv(
    rows: Iterable[Mapping[str, object]],
    columns: Sequence[str],
    output_path: Path,
    drop_columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    row_list = [dict(row) for row in rows]

    if not row_list:
        exported = pd.DataFrame(columns=list(columns))
    else:
        normalized = []
        for row in row_list:
            normalized_row: MutableMapping[str, object] = dict(row)
            for column in columns:
                normalized_row.setdefault(column, None)
            normalized.append(dict(normalized_row))

        exported = pd.DataFrame(normalized)
        ordered = [c for c in columns if c in exported.columns]
        extras = [c for c in exported.columns if c not in ordered]
        exported = exported.loc[:, ordered + extras]

    if drop_columns:
        exported = exported.drop(columns=list(drop_columns))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    exported.to_csv(output_path, sep="\t", index=False)
    return exported


def default_columns(feature_set: FeatureSet) -> Sequence[str]:
    return feature_set.columns(BASE_CONTEXT_COLUMNS)
