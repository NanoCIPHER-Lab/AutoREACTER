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

# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.3.0] - 2026-08-21

### Added

- **Reaction progression workflow:** Added iterative reaction-product generation so AutoREACTER can generate templates from products formed in earlier reaction steps. This improves support for multistage reactions, copolymerizations, and small-molecule systems where the initial monomers alone are not sufficient to generate all required polymerization templates.
- **Reaction iteration control:** Added `reaction_iteration_depth` to control how many reaction-product iterations are attempted. The default value is `5`, and the loop can be disabled by setting `reaction_iteration_depth` to `false`.
- **Index-based functional-group detection:** Added index-based functional-group detection to track reactive atom positions more accurately during reaction progression.
- **Index-based reaction detection:** Added index-based reaction detection so products generated in earlier reaction steps can be checked for additional valid reactions.
- **RDKit/NetworkX reaction deduplication:** Added graph-based reaction deduplication to reduce redundant reaction pathways during reaction progression.
- **LAMMPS template deduplication:** Added NetworkX-based LAMMPS template deduplication using generated pre-reaction templates, post-reaction templates, and map files.
- **Wildcard template support:** Added support for LAMMPS wildcard map generation to reduce duplicate edge-template cases when using supported LAMMPS versions.
- **Expanded reaction libraries:** Added and reorganized polymer reaction libraries for epoxy-amine, polyamide, polyester, polycarbonate, polysiloxane, polyurea, polyurethane, vinyl, and related polymerization chemistries.
- **Expanded functional-group libraries:** Added modular functional-group libraries for nitrogen, oxygen, carboxyl/carbonyl, vinyl/alkene, sulfur, silicon, ring, active-center, and mixed AB-type groups.
- **TFE/vinyl support:** Added support for tetrafluoroethylene and additional vinyl polymerization workflows.
- **Public API session access:** Added `arx.session()` to expose the active AutoREACTER session for workflow inspection.
- **Unit tests:** Added broader unit-test coverage for input parsing, ARX CLI behavior, reaction preparation, LUNAR client utilities, force-field wrapper components, walkers, and LAMMPS writers.
- **Documentation pages:** Added new documentation pages for advanced options and template deduplication.

### Changed

- **Input parser:** Refactored input parsing and validation, including cleaner handling of workflow options, simulation setup fields, force-field aliases, and input schema checks.
- **Reaction-library organization:** Refactored reaction libraries into dedicated modules and registries instead of relying on one monolithic reaction-library file.
- **Functional-group organization:** Refactored functional-group definitions into modular registry-based libraries.
- **Reaction preparation:** Refactored reaction preparation to support reaction progression, deduplication, inactive-template filtering, and clearer error handling.
- **REACTER file handling:** Simplified REACTER file metadata storage by storing LAMMPS molecule paths directly on monomer entries and template/map paths directly on reaction metadata.
- **LAMMPS writers:** Updated LAMMPS input writers to use the refactored REACTER metadata and filter inactive reaction templates.
- **3D molecule preparation:** Improved 3D embedding and repair handling for congested or difficult polymer structures.
- **Examples:** Replaced older test-style example JSON files with cleaner v0.3 example inputs.
- **Documentation structure:** Reorganized documentation into clearer user-facing and developer-facing pages.

### Fixed

- Fixed reaction progression issues where products from earlier steps were not correctly reused for later reaction detection.
- Fixed index-alignment issues during functional-group and reaction detection.
- Fixed radical handling for vinyl polymerization products and reaction deduplication.
- Fixed duplicate-template detection behavior for both RDKit-level reaction metadata and LAMMPS-level template files.
- Fixed handling of inactive or duplicate reaction templates so they are skipped in later workflow stages.
- Fixed empty-reaction cases so AutoREACTER raises clearer errors when no valid reaction instances are found.
- Fixed LAMMPS molecule/template file path handling after the REACTER metadata refactor.
- Fixed cache staging behavior during output directory preparation.
- Fixed documentation heading and toctree issues for cleaner Sphinx builds.

### Removed

- Removed legacy compatibility shims.
- Removed unused placeholder detector, fragment-comparison, and legacy library files.
- Removed older cluttered example JSON files in favor of focused v0.3 examples.

## [0.2.3] - [2026-06-24]

### Added
* Added a new `Session` class to centralize session-based input/output handling across the entire package.
* Added a new `ARXCLI` module to provide a streamlined Command-Line Interface for running the AutoREACTER workflow with a single command.
  - The CLI supports input file parsing, session management, and orchestrating the end-to-end workflow.


### Changed
* Centralized state management across the entire package. The `Session` object is now the single source of truth. All modules and functions have been refactored to use the `Session` for accessing inputs, outputs, and configuration.
* Centralized the versioning system to dynamically pull from the package's `__version__` attribute.
* Updated all documentation, examples, and CI configuration to reflect the new API and workflow changes.

## [0.2.2] - 2026-06-09

### Added

* Added a new PyPI-ready version of **AutoREACTER**.
* Added new reactions to the **reaction library**.
* Added the `ARXSession` module for centralized session-based input/output handling.
* Added unit tests for the LUNAR API wrapper.
* Added `MANIFEST.in` to support packaging of required non-Python files.
* Added outputs save direction same as the input file directory. 

### Changed

* Refactored cache and staging logic to use temporary/session-based directories instead of manually managed `base_dir` paths.
* Updated `.frc` file handling to use package-included force field files.
* Moved `frc_files` inside the AutoREACTER package for more reliable access after installation.
* Improved PCFF `.frc` directory discovery to be more robust against repository layout changes.
* Updated `prepare_reactions.py`, `molecule_3d_preparation.py`, `REACTER_files_builder.py`, and `lunar_api_wrapper.py` to use session-based inputs.
* Updated documentation, examples, CI configuration, and getting-started instructions.

### Deprecated

* Deprecated direct `base_dir`-style path handling in favor of `ARXSession`-based workflow management.

### Removed

* Removed older manual cache directory resolution logic.
* Removed temporary/editor backup files from the repository.
* Removed outdated path-resolution behavior that depended on repository-relative assumptions.
* Removed the old dated folder output method.

### Fixed

* Fixed `.frc` file path resolution issues when running AutoREACTER from different locations.
* Improved error messages for missing or unresolved `.frc` files.

### Security

* No security-related changes in this release.

---
## [0.2.1] - 2026-05-04

### Fixed
- Fixed issue where `cache/00_cache` was not cleared when initiating a new AutoREACTER run.
- Improved cache handling to prevent bugs caused by deleting cache files in the middle of a run.
- Added a warning message when the `--cleanup` command is used to delete all cached runs, to prevent accidental data loss.

## [0.2.0] - 2026-05-01

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