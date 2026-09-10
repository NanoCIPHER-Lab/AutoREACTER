#!/usr/bin/env python3
"""
AutoREACTER Workflow Example
============================

This script demonstrates the standard AutoREACTER workflow using a JSON
input file.

It is intended for users running AutoREACTER from the source repository.

Example
-------
From the AutoREACTER repository root:

    python examples/example_1.py \
        -i examples/polyamide_count_mode_basic.json

AutoREACTER can also be installed from PyPI and used directly from a
user-created Python script:

    python -m pip install AutoREACTER
"""

import argparse
import sys


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run the standard AutoREACTER workflow "
            "from an input JSON file."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        "-in",
        required=True,
        help="Path to the AutoREACTER input JSON file.",
    )

    return parser.parse_args()


def main():
    """Run the standard AutoREACTER workflow."""

    args = parse_arguments()

    # -------------------------------------------------------------------------
    # 1. Import AutoREACTER
    # -------------------------------------------------------------------------

    try:
        import AutoREACTER as arx
    except ImportError:
        print(
            "[ERROR] Failed to import AutoREACTER.\n"
            "Install it with:\n\n"
            "    python -m pip install AutoREACTER\n"
        )
        raise

    print("[OK] AutoREACTER imported successfully.")
    print(f"[INFO] AutoREACTER version: {arx.__version__}")

    # -------------------------------------------------------------------------
    # 2. Initialize workflow
    # -------------------------------------------------------------------------

    arx.run(args.input)

    # -------------------------------------------------------------------------
    # 3. Select reactions
    # -------------------------------------------------------------------------

    arx.select_reactions()

    # -------------------------------------------------------------------------
    # 4. Select non-reactants
    # -------------------------------------------------------------------------

    arx.select_non_reactants()

    # -------------------------------------------------------------------------
    # 5. Prepare reaction templates
    # -------------------------------------------------------------------------

    arx.prepare_reactions()

    # -------------------------------------------------------------------------
    # 6. Review checkpoint
    # -------------------------------------------------------------------------

    session = arx.session()

    print(
        "\n[INFO] Reaction preparation completed."
        f"\n[INFO] Review generated images in: {session.images_dir}"
    )

    confirmation = input(
        "\nType 'ok' to continue with final processing: "
    ).strip().lower()

    if confirmation != "ok":
        print("[EXIT] Workflow stopped by user.")
        return 0

    # -------------------------------------------------------------------------
    # 7. Final processing
    # -------------------------------------------------------------------------

    arx.process()

    print(
        "\n[INFO] AutoREACTER workflow completed successfully."
        f"\n[INFO] Output directory: {session.output_dir}"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())