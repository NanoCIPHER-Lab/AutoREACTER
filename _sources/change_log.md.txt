![AutoREACTER Logo](_static/logo.png)
# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
<!--
DRAFT - USE THIS FOR FUTURE RELEASES

Use this section to track changes you are currently working on but haven't "officially" labeled as a version yet.

## [X.Y.Z] - YYYY-MM- DD keep data format ISO standard
Added : for new features.
Changed:  for changes in existing functionality.
Deprecated: for soon-to-be removed features.
Removed: for now removed features.
Fixed: for any bug fixes.
Security: in case of vulnerabilities.


You can keep Unreleased changes here until you are ready to label a new version. Once you label a version, move the unreleased changes to a new section with the version number and date.
This serves two purposes:
    People can see what changes they might expect in upcoming releases
    At release time, you can move the Unreleased section changes into a new release version section.
 -->
---

## v0.2.0 - 2026-04-17

### BREAKING CHANGES
1. **Architectural Shift: Migrated from a script-based workflow to a Class-Based Architecture. Not backward compatible**

2. Input Schema (input.json): Replicas now require explicit per-instance density and temperature definitions. Added force field specification on input.json.

3. Atom Mapping: Switched from RDKit map numbers to isotope numbers. 0.1.0v relying on map numbers for atom tracking will require migration.

4. Directory Structure: Outputs are no longer saved to the root directory; they are now stored in dated run folders.

5. Library-Driven Workflow: Included functional_groups_library.py and reactions_library.py" to store structured dictionary definitions and SMIRKS transformations for automated reaction discovery.

### Added
1. Command-Line Interface (CLI): Introduced AutoREACTER.py as the primary entry point, supporting automated workflows with --help and --cleanup flags.

2. Interactive Visualizations: Added Jupyter Notebook integration for step-by-step visual validation of monomers, functional groups, and reactant/product templates.

3. Modern Atom Typing: Native support for string-based LAMMPS atom types (e.g., PCFF) through internal integration.

4. Force Field Selection: Users can now explicitly specify supported force fields directly within the input.json configuration.

5. Retention Utility: Added RetentionCleanup to manage disk space by removing aged simulation runs.

6. Full Logic Separation: Each stage now has its own dedicated detector and library (e.g., functional_groups, reactions), allowing for isolated debugging and cleaner updates.

7. Non-Monomer Control: Added a specialized detector for non-reactive species, giving users the choice to keep or discard them in the final simulation.

### Changed
1. Core Logic: Refactored workflow into modular classes: InputParser, Detector, ReactionSelection, TemplateBuilder, and LUNARClient.

2. Mapping Robustness: Atom tracking now utilizes isotopes (will be deleted after completion of atom mapping) to prevent the loss of metadata during chemical reactions.

3. Input Optimization: Refactored the preparation logic to generate only necessary LAMMPS input files based on explicit replica definitions, reducing file bloat.

### Fixed
1. Path Validation: Implemented strict directory validation for LUNAR, the system now rejects invalid paths during the input phase.

2. Index Corruption: Resolved an issue where reaction indices were lost when dropping reactions. By utilizing Python dataclasses and setting an `active_status = False`, the system now preserves the original data structure and prevents data leaks or data losses.

3. Duplicate Handling: Fixed bugs causing duplicate reaction instances during the discovery phase. The system now uses dataclass-based validation to flag duplicates as inactive `(active_status = False)` rather than removing them, preventing data crashes and ensuring a clean run.

---

## v0.1.0 - 2026-02-17

### Added

1. Initial Release: Formal launch of AutoREACTER as a script-based automation workflow for polymer simulation preparation.

2. Jupyter-Exclusive Execution: Initial workflow designed specifically for interactive execution within Jupyter Notebook environments.

3. SMARTS-Based Detection: Implementation of a core Functional Group Detector using SMARTS-based pattern matching to classify monomers (only supported for step growth monomers only).

4. Library-Driven Workflow: Included functional_groups_library.py and reactions_library.py to store structured dictionary definitions and SMIRKS transformations for automated reaction discovery.

5. LUNAR Integration: Automated passing of monomer and non-monomer MOL files to LUNAR for 3D geometry and initial property preparation.

6. Deterministic Reindexing: Implementation of NumPy-based reindexing to generate consistent reaction template molecules and mapping files for REACTER.

7. LAMMPS Support: Generation of ready-to-run REACTER input files for LAMMPS, utilizing legacy integer-based atom typing.

8. Unified Input Parsing: Introduced input_parser.py to validate a single structured inputs.json file defining temperatures, densities, and SMILES-based monomer compositions.