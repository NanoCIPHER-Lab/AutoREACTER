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

## [0.2.0] - 2026-04-17

### Breaking
- Migrated to a class-based architecture (not backward compatible).
- Updated input schema: explicit density and temperature per replica required.
- Replaced RDKit atom map numbers with isotope-based tracking.
- Changed output structure to dated run directories.
- Introduced library-driven workflow for reactions and functional groups.

### Added
- CLI interface (AutoREACTER.py) with --help and --cleanup support.
- Jupyter-based interactive visualization tools.
- Support for string-based LAMMPS atom types (e.g., PCFF).
- Force field selection via input.json.
- Retention cleanup utility for managing old runs.
- Modular detectors and libraries for workflow stages.
- Control over non-reactive species retention.

### Changed
- Refactored core workflow into modular classes.
- Improved atom mapping robustness using isotopes.
- Optimized input generation to reduce unnecessary LAMMPS files.

### Fixed
- Enforced strict validation of external tool paths.

---

## [0.1.0] - 2026-02-17

### Added

- Initial release of AutoREACTER as a script-based workflow for polymer simulation preparation.
- Support for Jupyter Notebook–based interactive execution.
- Functional group detection using SMARTS-based pattern matching (step-growth monomers only).
- Library-driven workflow with structured definitions for functional groups and reactions.
- Deterministic reindexing for consistent reaction templates and mapping files.
- Generation of REACTER-compatible LAMMPS input files using integer-based atom typing.
- Unified input parsing via a single structured input file for simulation configuration.