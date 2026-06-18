#!/usr/bin/env python3
"""
AutoREACTER Workflow: Required Pipeline Example
===============================================

This script demonstrates the required AutoREACTER workflow only.
Visualization images are generated internally and saved in the session image directory.

Usage:
    python auto_reacter_workflow.py -i example_1_inputs_count_mode.json
"""

import argparse
import sys


# ---------------------------------------------------------------------------
# 1. Command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Run the required AutoREACTER workflow from an input JSON file."
)

parser.add_argument(
    "-i",
    "--input",
    "-in",
    required=True,
    help="Path to the AutoREACTER input JSON file."
)

args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Import AutoREACTER
# ---------------------------------------------------------------------------
try:
    import AutoREACTER as arx
    print("AutoREACTER imported successfully.")
except ImportError:
    print(
        "Failed to import AutoREACTER. "
        "Ensure the package is installed in your environment."
    )
    raise


# ---------------------------------------------------------------------------
# 3. Version check & input load
# ---------------------------------------------------------------------------
# Print the installed version for reproducibility and debugging.
print(arx.__version__)

# Load the JSON input file and initialize the AutoREACTER session.
arx.run(args.input)
session = arx.session()


# ---------------------------------------------------------------------------
# 4. Select reactions and non-reactants
# ---------------------------------------------------------------------------
# Select which detected reactions should be processed.
arx.select_reactions()

# Select non-reactant molecules, if any are detected.
arx.select_non_reactants()


# ---------------------------------------------------------------------------
# 5. Reaction template preparation
# ---------------------------------------------------------------------------
# Prepare reaction templates for downstream REACTER file generation.
arx.prepare_reactions()


# ---------------------------------------------------------------------------
# 6. Review checkpoint
# ---------------------------------------------------------------------------
# AutoREACTER saves visualization images internally in the session image directory.
print(
    "\n[INFO] Reaction preparation completed."
    "\n[INFO] Please review the generated images in the session image directory."
)

# Try to print the image directory if it is exposed by AutoREACTER.
session = arx.session()
img_dir = getattr(session, "img_dir", None) or getattr(session, "images_dir", None)

if img_dir is not None:
    print(f"[INFO] Image directory: {img_dir}")

ok_pass = input("\nType 'ok' to continue with final processing: ").strip().lower()

if ok_pass != "ok":
    print("[EXIT] Workflow stopped by user.")
    sys.exit(0)


# ---------------------------------------------------------------------------
# 7. Final processing
# ---------------------------------------------------------------------------
# Run 3D geometry setup, force-field generation, REACTER file building,
# and LAMMPS simulation setup.
arx.process()

print("\n[INFO] AutoREACTER workflow completed successfully.")
