# AutoREACTER v0.2.0 - README

**AutoREACTER** is a Python-based toolkit for managing and automating reaction modeling in LAMMPS, developed as part of the **Multiscale Polymer Toolkit (MuPT)**.  
This repository is in **beta** and under active development — APIs and functionality may change without notice.

## Current Reaction Support (Beta)

For now, the package supports **only the reaction types below** and their relevant functional groups.

### Polyesterification (Polycondensation)
1. Hydroxy–carboxylic acid polycondensation *(including amino acid self-/co-polycondensation)*   
2. Hydroxy acid halide polycondensation *(including self-condensation and mixed halide cases)*  
3. Diol + di-acid halide polycondensation  
4. Diol + di-carboxylic acid polycondensation  

### Polyamidation (Polycondensation)
5. Amino acid polycondensation *(including amino acid self-/co-polycondensation)*  
6. Diamine + di-carboxylic acid polycondensation  
7. Diamine + di-carboxylic acid halide polycondensation  

## Quick Start

### Installation

```bash
git clone https://github.com/NanoCIPHER-Lab/AutoREACTER
```
### How to Use

AutoREACTER can be run either directly from the command line for automated workflows or via Jupyter Notebooks for interactive, step-by-step execution with visual feedback. You will need to provide a JSON input file to run AutoREACTER. 

#### 1. Jupyter Notebook Mode (Interactive & Visual)
See: **[examples/README.md](https://github.com/NanoCIPHER-Lab/AutoREACTER/blob/main/examples/README.md)** for usage guidelines and examples.\
For a step-by-step workflow with detailed visualizations of monomers, functional groups, and templates, use the provided Jupyter notebooks. This mode is highly recommended for visualizing your reactions interactively before generating LAMMPS files.

#### 2. Command-Line Interface (CLI) Mode
Use the main `AutoREACTER.py` executable script to run the full workflow. You will need to provide a JSON input file detailing your system 
### How to Use

AutoREACTER can be run either directly from the command line for automated workflows or via Jupyter Notebooks for interactive, step-by-step execution with visual feedback.

First, set up and activate your Conda environment to ensure all necessary dependencies are installed:

```bash
conda create -n autoRX -y -c conda-forge python=3.13 numpy pandas rdkit ipykernel networkx
conda activate autoRX
```
Once activated, use the main AutoREACTER.py
```bash
# Run the automated workflow with your configuration file
python AutoREACTER.py -i path/to/your/input.json or python AutoREACTER.py --input path/to/your/input.jsonf

# View all available commands and flags
python AutoREACTER.py --help or python AutoREACTER.py -h
``

# Run the interactive cleanup utility to manage old cache/run directories
python AutoREACTER.py --cleanup N or python AutoREACTER.py -c N          → delete runs older than N days e.g., 7, 30,
python AutoREACTER.py --cleanup all or python AutoREACTER.py -c all      → delete all cached runs