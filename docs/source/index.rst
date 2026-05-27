.. AutoREACTER documentation master file, created by
   sphinx-quickstart on Tue Apr 14 10:57:05 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


AutoREACTER documentation
=========================

Welcome to the AutoREACTER documentation.

AutoREACTER is an open-source Python toolkit for automating the setup of input files for reactive molecular dynamics (MD) in LAMMPS. AutoREACTER fully automates the creation of all input files necessary to run REACTER, requiring as user input only the SMILES strings of the reacting molecules. Reactions are suggested from a customizable database currently populated with polymerization reactions.

GitHub repository: https://github.com/NanoCIPHER-Lab/AutoREACTER

**Note: AutoREACTER is currently in v0.2-beta and under active development. APIs, configuration schemas, and core functionality may change or break without notice as we expand reaction library and force field support**

REACTER is a method for modeling chemical reactions using classical fixed-bond force fields. The required inputs files are pre- and post-reaction templates and a map file providing the atom-to-atom mapping for the reaction.

AutoREACTER bridges the gap between raw chemical structures and REACTER-ready LAMMPS input files for atomistic simulations. AutoREACTER automatically generates atom-typed reaction templates, map files, and the initial monomer melt using simplified input files as described in Input Configuration Documentation.

.. toctree::
   :maxdepth: 1
   :caption: User Guide:

   overview.md
   getting-started.md
   input-configuration.md
   supported-reactions.md
   supported-force-fields.md
   clean_up.md
   change_log.md
   contact.md

.. toctree::
   :maxdepth: 2
   :caption: Developer API:

   modules