![AutoREACTER Logo](_static/logo.png)
# Welcome to AutoREACTER

**AutoREACTER v0.2.0**

AutoREACTER is a Python-based toolkit for managing and automating reaction modeling in LAMMPS, developed as part of the Multiscale Polymer Toolkit (MuPT).

**Note: AutoREACTER is currently in v0.2-beta. It is under active development, and APIs or functionality may change as we continue to expand the reaction library and force field support.**

AutoREACTER is a powerful Python-based toolkit designed to automate the setup and modeling of complex chemical reactions in LAMMPS. AutoREACTER bridges the gap between raw chemical structures and REACTER-ready LAMMPS input files for atomistic simulations.

For **[REACTER](https://www.reacter.org/)**, you need to manually prepare molecule files, pre-reaction templates, post-reaction templates, and map files. Then, you must assign atom types to these files and write input scripts for LAMMPS.

AutoREACTER removes this manual process of defining reaction templates, preparing MAP files, and building topology files—using simplified input files as described in *Input Configuration Documentation*.

## Advantages of AutoREACTER

1. **Dual Workflows:** Run automated pipelines via the CLI, or use interactive Jupyter Notebooks to visually validate intermediate steps, functional groups, and template mappings before writing a single LAMMPS file.

2. **Automatic Reaction Detection:** Molecules are first filtered based on functional group availability. Based on these functional groups, the necessary reaction templates for REACTER are automatically generated.

3. **Automated Force Field Integration:** Fully integrated with the LUNAR API to automatically assign partial charges, and for all reactant and product templates. The available force fields are listed in AutoREACTER-Supported Force Fields documentation.

4. **Smart File Generation:** Automatically generates fully configured LAMMPS input scripts, including system densification, pre-equilibration, and multi-stage reaction files based on temperature, density, and monomer/molecule ratios or counts.

## How AutoREACTER Works

AutoREACTER is built on a clean, modular, class-based architecture. When you submit your simulation setup, the toolkit processes it through five core stages:

1. **Input Processing:**  
   Reads your `input.json` file and validates system replicas, temperatures, densities, monomer SMILES strings, stoichiometric ratios, and force field selections.

2. **Functional Group & Reaction Detection:**  
   Analyzes monomers to identify available functional groups. These are cross-referenced with the internal reaction library to determine all possible reactive pathways. Non-reactive molecules are identified separately.

3. **Reaction Template Generation:**  
   Uses RDKit and mapped SMARTS strings to generate pre- and post-reaction molecules and identify initiator atoms.  
   A 1:1 atom mapping is performed between pre- and post-reaction structures. The `walker.py` module traverses the molecular graph from initiator atoms to extract templates and edges. RDKit is also used to identify by-products throughout the reaction process.

4. **Force Field Assignment:**  
   Communicates with the LUNAR API to assign force field parameters, atom types, and partial charges.  
   LUNAR repository: [https://github.com/CMMRLab/LUNAR](https://github.com/CMMRLab/LUNAR)

5. **LAMMPS File Generation:**  
   Based on the provided inputs and generated chemical data, AutoREACTER writes the required LAMMPS input scripts.

---

**Note:** In the documentation and within AutoREACTER, molecules are often referred to as *monomers*, as this package is primarily designed for polymer systems.
