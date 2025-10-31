# vesper
`vesper` is a toolkit for annotating and refining *de novo* transposable elements (TEs) from PacBio HiFi data.

```
vesper/
├── src/vesper/
│   ├── cli.py
│   ├── analysis/
│   ├── commands/
│   ├── misc/
│   ├── models/
│   ├── processors/
│   └── utils/
│       └── config.py
├── notebooks/
├── tests/           
├── environment.yml  
└── pyproject.toml   
```

## Usage

`vesper` performs two main functions: `annotate`, and `refine`.

1. `vesper annotate` - Provide a VCF file to annotate insertion sequences with RepeatMasker. External BED/GFF/TSV annotations are optional.
    - Yields `{filename}.annotated.vcf.gz`: variants with RepeatMasker results and optional external annotations in the INFO field
2. `vesper refine` - Provide a VCF file and the supporting BAM file to calculate read-based confidence scores.
    - Yields `{filename}.refined.vcf.gz`: variants with confidence scores based on read quality metrics

## Installation

We strongly suggest using `conda` (or equiv. `mamba`) to install `vesper`


1. Directly create and install the environment from `environment.yml`:
```bash
conda env create -f environment.yml
```

*Note*: This environment should build without any issues on Linux and osx-64. However, not all packages are compiled for osx-arm64 (Apple Silicon): if you are using an M-series Mac, you can use the `--subdir` option to emulate osx-64:

```bash
conda env create -f environment.yml --subdir osx-64
```

2. Activate the environment:
```bash
conda activate vesper
```

3. Install `vesper` via `pip` or `uv` (`bioconda` upload soon!)
```bash
# pip
pip install -e .

# uv
uv pip install -e .
```

4. Verify that everything is installed:

```bash
vesper --help
```

## Workflow

A typical `vesper` run takes an input VCF from `sniffles2` using `--qc-all` mode, which allows `vesper` to screen for ultra-low frequency variants. An example workflow is provided in `Snakefile` and can be run with `snakemake`.

## Documentation

Extended documentation for `vesper annotate` and `vesper refine` are available in the `docs/` folder.

## Troubleshooting

`environment.yml` requests `repeatmasker==4.1.2.p1` which may cause some errors due to `python` shebang line at the head of certain files (`'share/RepeatMasker/famdb.py`, `share/RepeatMasker/util/RM2Bed.py`)

```
#!/usr/bin/env python3
```

You can remedy this by simply modifying it to refer to `python` instead of `python3`. You may be able to install a newer version of RepeatMasker as well: this behavior has been tricky to replicate, so please file an issue if you find a way to consistently resolve it.