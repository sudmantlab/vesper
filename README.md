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

`vesper` performs three main functions: `call`, `annotate`, and `refine`. You have two options for running `vesper`:

1. Provide a FASTQ file of PacBio HiFi reads. `vesper call` will align your reads and perform first-pass structural variant calling using one or more of the following approaches.
    - Reference genome alignment
    - Haplotype-aware de novo assembly with reference-guided scaffolding
    - Pangenome graph alignment and surjection to a reference genome
2. Provide a BAM file, reference genome + annotations, and a VCF file of candidate variants to be refined. 
    - `vesper annotate` will annotate the variants.
    - `vesper refine` will score variants and filter for low-confidence regions in the reference and/or assembly.

> **Note:** At the moment, only retrotransposon discovery is supported 😔

The output of all three functions is a VCF file with varying levels of confidence:
- WIP: vesper call` - `{filename}.candidates.vcf.gz`, unfiltered candidate variants
- `vesper annotate` - `{filename}.annotated.vcf.gz`, variants with reference annotations
- `vesper refine` - `{filename}.refined.vcf.gz`, variants with both annotations and confidence scores

### Utilities (WIP)

`vesper` provides several utilities for examining and manipulating variants.
- `vesper construct` can be used to perform local reassembly of complex mosaic variants.
- To do...

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

*Note*: The `requirements.txt` file is generated via `pip freeze > requirements.txt` for noting specific versions in the most recent functional build, but isn't used for repopulating the environment.

## Testing

After setting up either environment, you can run tests with:
```bash
pytest
```

This will:
- Run all tests in the `tests/` directory
- Generate coverage reports
- Show test results in verbose mode