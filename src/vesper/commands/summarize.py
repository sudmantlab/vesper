from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Iterable, List

from vesper.utils.config import SummarizeConfig

HEADER = (
    "SAMPLE\tCHROM\tPOS\tID\tQUAL\tFILTER\tCONFIDENCE\tCONFIDENCE_FLAGS\t"
    "RSFTCLIP\tNSFTCLIP\tRMAPQ\tNSMAPQ\tSVLEN\tREPEATMASKER_RESULTS\tOVERLAPPING\n"
)

QUERY_FORMAT = (
    r"%CHROM\t%POS\t%ID\t%QUAL\t%FILTER\t%INFO/CONFIDENCE\t%INFO/CONFIDENCE_FLAGS\t"
    r"%INFO/RSFTCLIP\t%INFO/NSFTCLIP\t%INFO/RMAPQ\t%INFO/NSMAPQ\t%INFO/SVLEN\t"
    r"%INFO/REPEATMASKER_RESULTS\t%INFO/OVERLAPPING\n"
)


def _infer_sample_name(vcf_path: Path) -> str:
    name = vcf_path.name
    for suffix in (".vcf.gz", ".vcf"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if "." in name:
        name = name.split(".")[0]
    return name


def _gather_rows(vcf_path: Path, sample_name: str, output_handle, logger: logging.Logger) -> None:
    command = ["bcftools", "query", "-f", QUERY_FORMAT, str(vcf_path)]
    logger.debug(f"Running command: {' '.join(command)}")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        output_handle.write(f"{sample_name}\t{line}\n")

    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()
    if return_code != 0:
        message = stderr.strip() or "bcftools query failed"
        raise RuntimeError(f"{vcf_path}: {message}")


def run_summarize(config: SummarizeConfig, logger: logging.Logger) -> None:
    if not config.vcf_inputs:
        raise ValueError("No input VCF files provided.")

    for path in config.vcf_inputs:
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")

    config.output_path.parent.mkdir(parents=True, exist_ok=True)

    if config.sample_names:
        sample_names: Iterable[str] = config.sample_names
    else:
        sample_names = [_infer_sample_name(path) for path in config.vcf_inputs]

    logger.info(f"Writing summary to {config.output_path}")
    with config.output_path.open("w", encoding="utf-8") as handle:
        handle.write(HEADER)
        for vcf_path, sample_name in zip(config.vcf_inputs, sample_names):
            logger.info(f"Summarizing {vcf_path} (sample '{sample_name}')")
            _gather_rows(vcf_path, sample_name, handle, logger)

    logger.info("Summarization completed successfully.")
