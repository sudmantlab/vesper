# vesper
                                                                                
                                                                   ..           
                                                                   .-           
                                                                   --           
                                                                  .--           
                                                                  ----          
                                                                ..----..        
                                                            .--------------.    
                                 ...........                  ..--------..      
                                              ...........         ----          
                                                         ........ .--.          
     -##+    ##. -#-+#+   .#--+#+  ##+#####.    +#-##-  -##-###....--.          
      +##    #. -#-  +##  ##   .+  +#+   +##.  ##. .##-  ###  +   .--....       
       ##-  +- .##.-+#-.  ####.    +#+    ##+ -## -#+-   ##-       .-......     
       -## .#  .##.        .####-  +#+    +## ###        ##-       .........    
        ##+#.  .##+     . .   -##. +#+    +#- -##.     . ##-       . .......    
        .##+    -##+  .+  #.   ##  +##   .#+   ###-  .+  ##-        .......     
         -#.     .+###-   -+###-   +#####+.     .####.  -##+    ........        
                                   +#+                     ........             
                                   ###              .........                   
                                             ........                           
                                 .............                                  
                          ........                                              
                       .......                                                  
                       .......                                                  
                            ...........                                          

`vesper` is a haplotype-aware *de novo* structural variant caller  for PacBio HiFi data.

```
vesper/
├── src/vesper/
│   ├── models/      # Data models and types
│   ├── processors/  # File processing and analysis
│   └── utils/       # Utility functions and config
├── notebooks/       # Jupyter notebooks for examples and analysis
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
- `vesper call` - `{filename}.candidates.vcf.gz`, unfiltered candidate variants
- `vesper annotate` - `{filename}.annotated.vcf.gz`, variants with reference annotations
- `vesper refine` - `{filename}.refined.vcf.gz`, variants with both annotations and confidence scores

### Utilities

`vesper` provides several utilities for examining and manipulating variants.
- `vesper construct` can be used to perform local reassembly of complex mosaic variants.
- To do...

## Installation

Use `venv` (Python's built-in virtual environment) or `conda` to install. (TODO: build Docker image before release)

### Option 1: Using venv + pip/uv

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

2. Install dependencies (use `pip` or `uv`):
```bash
# pip (dev install)
pip install -e ".[dev]"

# pip
pip install .

# uv
uv pip install .
```

### Option 2: Using conda

Create and activate the conda environment (use `conda` or `mamba`):
```bash
# Create environment from specification
conda env create -f environment.yml

# Activate the environment
conda activate vesper
```


## Testing

After setting up either environment, you can run tests with:
```bash
pytest
```

This will:
- Run all tests in the `tests/` directory
- Generate coverage reports
- Show test results in verbose mode

## Environment Management Tips

### Updating Dependencies

If you add new dependencies:
1. Add them to `pyproject.toml`
2. Update `requirements.txt`:
```bash
pip freeze > requirements.txt
```
3. Update `environment.yml` manually with the new package

### Recreating Environments

For venv:
```bash
# if in venv, deactivate it
deactivate
# Remove old environment
rm -rf venv/
# Create new one following installation steps
```

For conda:
```bash
# Remove old environment
conda deactivate
conda env remove -n vesper
# Create new one
conda env create -f environment.yml
``` 