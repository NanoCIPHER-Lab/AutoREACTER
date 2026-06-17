#!/usr/bin/env python3
"""
AutoREACTER Workflow: Molecular Reaction Screening & Analysis
==============================================================

This script demonstrates a typical AutoREACTER (Automated REACTion ExploreR)
workflow for screening functional groups, selecting reaction partners, and
processing reaction templates in a combinatorial or high-throughput manner.

The workflow proceeds through the following stages:
    1. Import and version check
    2. Load a configuration / input file (JSON format)
    3. Visualize molecules, functional groups, and reactions
    4. Select reactive species (reactions and non-reactants)
    5. Prepare and display reaction templates
    6. Execute the final processing step

Usage:
    python auto_reacter_workflow.py

Requirements:
    - AutoREACTER package installed (imported as ``arx``)
    - A valid JSON input file (e.g., ``example_1_inputs_count_mode.json``)
"""

# ---------------------------------------------------------------------------
# 1. Import AutoREACTER
# ---------------------------------------------------------------------------
try:
    import AutoREACTER as arx
    print("AutoREACTER imported successfully.")
except ImportError:
    print("Failed to import AutoREACTER. "
          "Ensure the package is installed in your environment.")
    # Re-raise to halt execution if the core dependency is missing
    raise

# ---------------------------------------------------------------------------
# 2. Version check & configuration load
# ---------------------------------------------------------------------------
# Display the installed version for reproducibility / debugging
print(arx.__version__)

# Load the JSON input file that defines molecules, functional groups,
# reaction rules, and processing parameters for the current run.
arx.load("example_1_inputs_count_mode.json")

# ---------------------------------------------------------------------------
# 3. Inspection: visualise the loaded data
# ---------------------------------------------------------------------------

# Display all molecules parsed from the input file.
# Useful for verifying that SMILES strings, names, and counts are correct.
arx.show_molecules()

# Display functional groups detected / defined on each molecule.
# Helps confirm that the correct reactive sites have been identified.
arx.show_functional_groups()

# Display the full list of possible reactions (before any filtering).
# Each reaction is typically a tuple of (reactants, products, template).
arx.show_reactions()

# ---------------------------------------------------------------------------
# 4. Selection: narrow down the reaction space
# ---------------------------------------------------------------------------

# Interactively or automatically select which reactions to keep.
# This may filter by energetics, user-defined rules, or duplicate removal.
arx.select_reactions()

# Specify molecules that should NOT participate as reactants
# (e.g., solvents, catalysts, or spectator species).
arx.select_non_reactants()

# ---------------------------------------------------------------------------
# 5. Template preparation
# ---------------------------------------------------------------------------

# Prepare reactions for template extraction: canonicalise SMILES,
# map atom indices, and group equivalent transformations.
arx.prepare_reactions()

# Display reaction templates of a given type.
# The "delete" argument likely shows templates that involve bond-breaking
# or atom removal patterns.
arx.show_reaction_templates("delete")

# ---------------------------------------------------------------------------
# 6. Final processing
# ---------------------------------------------------------------------------

# Execute the full pipeline: enumerate products, score candidates,
# and write output files (CSV / SDF / HTML reports).
arx.process()
