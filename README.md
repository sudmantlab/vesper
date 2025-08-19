# vesper
`vesper` is a haplotype-aware *de novo* structural variant caller for PacBio HiFi data.

```
vesper/
├── src/vesper/
│   ├── cli.py
│   ├── commands/
│   ├── models/
│   ├── processors/
│   └── utils/
│       └── config.py
├── notebooks/
├── tests/           
├── environment.yml  
├── requirements.txt 
└── pyproject.toml   
```

## Usage

`vesper` performs three main functions: `call` (not implemented yet), `annotate`, and `refine`. You have two options for running `vesper`:

1. `vesper annotate` - Provide a VCF file to annotate insertion sequences with RepeatMasker. External BED/GFF/TSV annotations are optional.
    - Yields `{filename}.annotated.vcf.gz`: variants with RepeatMasker results and optional external annotations in the INFO field
2. `vesper refine` - Provide a VCF file and the supporting BAM file to calculate read-based confidence scores.
    - Yields `{filename}.refined.vcf.gz`: variants with confidence scores based on read quality metrics

### Utilities (WIP)

- `vesper construct` (not implemented yet) – perform local reassembly of complex mosaic variants.
- `vesper merge` (not implemented yet) – merge variants across multiple VCF files.

## Installation

Use `venv` (Python's built-in virtual environment) or `conda` to install. (TODO: build Docker image before release)

### Option 1: Using venv

TODO: add instructions for using `conda`, proper `requirements.txt` and `environment.yml` files, `pip` install via pypi, etc...

1. Create and activate a virtual environment:
```bash
# Create virtual environment
python -m venv venv

# Activate it (choose based on your shell)
# On Unix/macOS:
source venv/bin/activate
# On Windows:
.\venv\Scripts\activate
```

2. Install dependencies:
```bash
# pip (dev install)
pip install -e ".[dev]"

# uv (slightly faster)
pip install uv
uv pip install -e ".[dev]"
```

`requirements.txt` is generated via `pip freeze > requirements.txt` for noting specific versions in the most recent functional build and can also be used to install dependencies.

### Option 2: Using conda

#TODO: Fix overlapping dependency specifications in `environment.yml` and `pyproject.toml`.

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

3. Install `vesper` via `pip` or `uv`:
```bash
# pip
pip install -e ".[dev]"

# uv
uv pip install -e ".[dev]"
```

`environment-list.txt` is generated via `conda env export > environment-list.txt` and can be used to re-create the environment with hard pins.

## Testing

After setting up either environment, you can run tests with:
```bash
pytest
```

This will:
- Run all tests in the `tests/` directory
- Generate coverage reports
- Show test results in verbose mode