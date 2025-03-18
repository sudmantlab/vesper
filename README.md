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

1. `vesper annotate` - Provide a VCF file of candidate variants to be annotated and one or more BED/GFF files of annotations.
    - Yields `{filename}.annotated.vcf.gz`: variants with annotations in the INFO field
2. `vesper refine` - Provide a VCF file of (annotated) candidate variants to be refined and the supporting BAM file.
    - Yields `{filename}.refined.vcf.gz`: variants with both annotations and confidence scores

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

## Testing

After setting up either environment, you can run tests with:
```bash
pytest
```

This will:
- Run all tests in the `tests/` directory
- Generate coverage reports
- Show test results in verbose mode