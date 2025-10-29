import io
import logging
from types import SimpleNamespace

import pytest

from vesper.commands.summarize import run_summarize, HEADER
from vesper.utils import SummarizeConfig


class DummyProcess:
    def __init__(self, stdout_text: str, stderr_text: str = "", return_code: int = 0):
        self.stdout = io.StringIO(stdout_text)
        self.stderr = io.StringIO(stderr_text)
        self._return_code = return_code

    def wait(self) -> int:
        return self._return_code


def test_run_summarize_infers_samples(tmp_path, monkeypatch):
    vcf1 = tmp_path / "sample1.annotated.refined.vcf.gz"
    vcf2 = tmp_path / "sample2.annotated.refined.vcf.gz"
    vcf1.write_text("dummy")
    vcf2.write_text("dummy")
    output_path = tmp_path / "summary.tsv"

    outputs = iter([
        "chr1\t1\tsample1_var\t10\tPASS\t1\tNONE\t0\t0\t10\t20\t100\t.\t.\n",
        "chr2\t2\tsample2_var\t20\tPASS\t0.8\tMAPQ_DIFF\t5\t2\t25\t12\t200\t.\t.\n"
    ])

    monkeypatch.setattr(
        "vesper.commands.summarize.subprocess.Popen",
        lambda *args, **kwargs: DummyProcess(next(outputs))
    )

    args = SimpleNamespace(
        input=[str(vcf1), str(vcf2)],
        output=str(output_path),
        sample_names=None,
        logging=None,
        debug=False,
        console_output=False
    )
    config = SummarizeConfig.from_args(args)

    logger = logging.getLogger("test_run_summarize_infers_samples")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())

    run_summarize(config, logger)

    contents = output_path.read_text().strip().splitlines()
    assert contents[0] == HEADER.strip()
    assert contents[1].startswith("sample1\tchr1\t1")
    assert contents[2].startswith("sample2\tchr2\t2")


def test_summarize_config_sample_name_mismatch():
    args = SimpleNamespace(
        input=["a", "b"],
        output="summary.tsv",
        sample_names=["sample_a"],
        logging=None,
        debug=False,
        console_output=False
    )
    with pytest.raises(ValueError):
        SummarizeConfig.from_args(args)


def test_run_summarize_missing_bcftools(tmp_path, monkeypatch):
    vcf_path = tmp_path / "sample.annotated.refined.vcf.gz"
    vcf_path.write_text("dummy")
    config = SummarizeConfig(
        vcf_inputs=[vcf_path],
        output_path=tmp_path / "summary.tsv",
        sample_names=None,
        log_dir=tmp_path / "logs",
        debug=False,
        console_output=False
    )

    def raising_popen(*args, **kwargs):
        raise FileNotFoundError("bcftools not found")

    monkeypatch.setattr("vesper.commands.summarize.subprocess.Popen", raising_popen)

    logger = logging.getLogger("test_run_summarize_missing_bcftools")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())

    with pytest.raises(RuntimeError):
        run_summarize(config, logger)
