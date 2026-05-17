# Welcome to AutoREACTER

**AutoREACTER v0.2-beta.0**

AutoREACTER is a Python-based toolkit for managing and automating reaction modeling workflows in LAMMPS. It is developed as part of the Multiscale Polymer Toolkit (MuPT).

**Note: AutoREACTER is currently in v0.2-beta.0. It is under active development, and APIs or functionality may change as we continue to expand the reaction library and force field support.**

AutoREACTER is designed to automate the setup and modeling of complex chemical reactions in LAMMPS. It bridges the gap between raw chemical structures and REACTER-ready LAMMPS input files for atomistic simulations.

For **[REACTER](https://www.reacter.org/)**, users typically need to manually prepare molecule files, pre-reaction templates, post-reaction templates, and map files. Then, atom types must be assigned to these files, and input scripts must be written for LAMMPS.

AutoREACTER reduces this manual process by generating reaction templates, preparing MAP files, building topology files, and writing LAMMPS input files using simplified configuration files, as described in the *Input Configuration Documentation*.

---

## Advantages of AutoREACTER

1. **Dual Workflows:**  
   Run automated pipelines through the CLI, or use interactive Jupyter Notebooks to visually validate intermediate steps, functional groups, and template mappings before writing a single LAMMPS file.

2. **Automatic Reaction Detection:**  
   Molecules are first filtered based on functional group availability. Based on these functional groups, the necessary reaction templates for REACTER are automatically generated.

3. **Automated Force Field Integration:**  
   Fully integrated with the LUNAR API to automatically assign partial charges for all reactant and product templates. Available force field options are selected through the input configuration workflow.

4. **Smart File Generation:**  
   Automatically generates fully configured LAMMPS input scripts, including system densification, pre-equilibration, and multi-stage reaction files based on temperature, density, and monomer/molecule ratios or counts.

---

## How AutoREACTER Works

AutoREACTER is built on a clean, modular, class-based architecture. When you submit your simulation setup, the toolkit processes it through five core stages:

### 1. Input Processing

AutoREACTER reads the `input.json` file and validates the simulation setup, including temperatures, densities, monomer SMILES strings, stoichiometric ratios, monomer or molecule counts, and force field selections.

### 2. Functional Group and Reaction Detection

AutoREACTER analyzes the input molecules to identify available functional groups. These functional groups are cross-referenced with the internal reaction library to determine all possible reactive pathways. Non-reactive molecules are identified separately.

### 3. Reaction Template Generation

AutoREACTER uses RDKit and mapped SMARTS strings to generate pre-reaction and post-reaction molecules and identify initiator atoms.

A 1:1 atom mapping is performed between pre-reaction and post-reaction structures. The `walker.py` module traverses the molecular graph from initiator atoms to extract templates and edges. RDKit is also used to identify by-products throughout the reaction process.

### 4. Force Field Assignment

AutoREACTER communicates with the LUNAR API to assign force field parameters, atom types, and partial charges.

LUNAR repository: [https://github.com/CMMRLab/LUNAR](https://github.com/CMMRLab/LUNAR)

### 5. LAMMPS File Generation

Based on the provided inputs and generated chemical data, AutoREACTER writes the required LAMMPS input scripts for the full simulation workflow.

This includes files for system densification, pre-reaction equilibration, reaction stages, and post-reaction equilibration.

---

## Terminology Note

In the documentation and within AutoREACTER, molecules are often referred to as *monomers*, as this package is primarily designed for polymer systems.