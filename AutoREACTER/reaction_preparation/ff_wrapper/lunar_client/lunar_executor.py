"""
LUNAR Executor Module

Handles the direct execution of LUNAR's underlying Python scripts 
(atom_typing.py, all2lmp.py, and bond_react_merge.py) via subprocesses.
"""

import os
import sys
import time
import subprocess
from pathlib import Path
from dataclasses import dataclass

from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata
from AutoREACTER.reaction_preparation.ff_wrapper.ff_validator import FFValidator

# Import central/shared data structures
from AutoREACTER.reaction_preparation.ff_wrapper.ff_wrapper import FFFiles, MoleculeFile, TemplateFile, DataFiles

# Local utility imports
from AutoREACTER.reaction_preparation.ff_wrapper.lunar_client.lunar_utils import move_merge_outputs, get_ending_integer

# =============================================================================
# LUNAR-Specific Data Structures
# =============================================================================

@dataclass(slots=True)
class AtomTypingResult:
    """Container for atom_typing.py output file paths."""
    id: str                  # Identifier: "data1", "pre1", "post1", etc.
    molecule: bool           # True for data molecules, False for reaction templates
    typed_data_file: Path    # Path to generated *_typed.data file
    nta_file: Path           # Path to generated *.nta file


@dataclass(slots=True)
class All2LMPResult:
    """Container for all2lmp.py conversion results."""
    id: str                  # Identifier matching the input
    molecule: bool           # True for molecules, False for templates
    all2lmp_data_file: Path  # Path to generated *_IFF.data file

class LunarExecutor:
    """
    Executes LUNAR scripts and manages their immediate input/output files.
    """
    def __init__(self, lunar_location: Path, cache_dir: Path):
        self.lunar_location = Path(lunar_location)
        self.atom_typing_py = self.lunar_location / "atom_typing.py"
        self.all2lmp_py = self.lunar_location / "all2lmp.py"
        self.bond_react_merge_py = self.lunar_location / "bond_react_merge.py"

        # Setup Stage Caches
        self.cache_atom_typing = Path(cache_dir) / "atom_typing"
        self.cache_all2lmp = Path(cache_dir) / "all2lmp"
        self.cache_bond_react_merge = Path(cache_dir) / "bond_react_merge"

        for p in [self.cache_atom_typing, self.cache_all2lmp, self.cache_bond_react_merge]:
            p.mkdir(parents=True, exist_ok=True)

    def run_atom_typing(
        self, 
        updated_inputs: SimulationSetup,
        prepared_reactions: list[ReactionMetadata],
        force_field: str
    ) -> list[AtomTypingResult]:
        """Runs LUNAR's atom_typing.py on all molecules and templates."""
        atom_typing_result = []

        def run_command(mol_file: Path, output_dir: Path):
            subprocess.run([
                sys.executable,
                str(self.atom_typing_py),
                '-topo', str(mol_file),
                '-dir', str(output_dir),
                '-ff', force_field,
                '-del-method', 'mass',
                '-del-crit', '0'
            ], check=True)

        # 1. Process Data Molecules (Monomers)
        for entry in updated_inputs.monomers:
            if getattr(entry, "status", True) is False:
                continue

            if not entry.molecule_3Dmol_path:
                raise ValueError(
                    f"Active monomer '{entry.name}' is missing molecule_3Dmol_path. "
                    "Ensure 3D geometry preparation completed successfully."
                )

            run_command(entry.molecule_3Dmol_path, self.cache_atom_typing)
            
            atom_typing_result.append(
                AtomTypingResult(
                    id=entry.name,
                    molecule=True,
                    typed_data_file=self.cache_atom_typing / f"{entry.name}_typed.data",
                    nta_file=self.cache_atom_typing / f"{entry.name}_typed.nta"
                )
            )
            time.sleep(0.1)

        # 2. Process Reaction Templates
        for reaction in prepared_reactions:
            for mol_path, label in [
                (reaction.reactant_combined_3Dmol_path, "pre"),
                (reaction.product_combined_3Dmol_path, "post"),
            ]:
                output_dir = self.cache_atom_typing / f"{label}{reaction.reaction_id}"
                output_dir.mkdir(parents=True, exist_ok=True)

                run_command(mol_path, output_dir)

                atom_typing_result.append(
                    AtomTypingResult(
                        id=f"{label}{reaction.reaction_id}",
                        molecule=False,
                        typed_data_file=output_dir / f"{mol_path.stem}_typed.data",
                        nta_file=output_dir / f"{mol_path.stem}_typed.nta"
                    )
                )
                time.sleep(0.1)

        # Verify outputs
        for files in atom_typing_result:
            if not files.typed_data_file.is_file() or not files.nta_file.is_file():
                raise FileNotFoundError(f"Expected LUNAR output not found for {files.id}")
            print(f"[LUNAR atom_typing] Generated files for {files.id}")

        return atom_typing_result

    def run_all2lmp(
        self, 
        atom_typing_results: list[AtomTypingResult], 
        frc_file: Path
    ) -> list[All2LMPResult]:
        """Runs LUNAR's all2lmp.py to convert typed data into LAMMPS format."""
        all2lmp_results = []

        for entry in atom_typing_results:
            subprocess.run([
                sys.executable,
                str(self.all2lmp_py),
                '-topo', str(entry.typed_data_file),
                '-nta', str(entry.nta_file),
                '-frc', str(frc_file),
                '-asm', 'T',
                '-dir', str(self.cache_all2lmp)
            ], check=True)
            
            all2lmp_results.append(All2LMPResult(
                id=entry.id,
                molecule=entry.molecule,
                all2lmp_data_file=Path(f"{entry.id}_typed_IFF.data")
            ))

        time.sleep(0.1)
        
        # Verify outputs
        for result in all2lmp_results:
            path = self.cache_all2lmp / result.all2lmp_data_file
            if not path.is_file():
                raise FileNotFoundError(f"Expected output file not found at {path}")
            print(f"[LUNAR all2lmp] Generated file for {result.id}")

        return all2lmp_results

    def run_bond_react_merge(
        self,
        merge_input_file_path: Path,
        all2lmp_results: list[All2LMPResult]
    ) -> FFFiles:
        """Executes bond_react_merge.py to create the final unified simulation setup."""
        env = os.environ.copy()
        env["QT_QPA_PLATFORM"] = "offscreen"
        subprocess.run(
            [
                sys.executable,
                str(self.bond_react_merge_py),
                "-files", f"infile:{merge_input_file_path.name}", 
                "-atomstyle", "full",
                "-tl", "T",
                "-wrd", "T",
            ],
            cwd=str(self.cache_bond_react_merge),
            env=env,
            check=True
        )

        move_merge_outputs(self.cache_all2lmp, self.cache_bond_react_merge)
        output_dir = self.cache_bond_react_merge

        molecule_files = []
        template_files = []

        # Parse outputs into unified FFFiles structure
        for result in all2lmp_results:
            name = result.id
            if result.molecule:
                molecule_files.append(MoleculeFile(
                    id=name,
                    molecule_files=DataFiles(
                        data_file=output_dir / f"{name}_typed_IFF_merged.data",
                        lmp_molecule_file=output_dir / f"{name}_typed_IFF_merged.lmpmol"
                    )
                ))
            else:
                if name.startswith("pre"):
                    rid = get_ending_integer(name)
                    template_files.append(TemplateFile(
                        reaction_id=rid,
                        pre_reaction_file=DataFiles(
                            data_file=output_dir / f"{name}_typed_IFF_merged.data",
                            lmp_molecule_file=output_dir / f"{name}_typed_IFF_merged.lmpmol"
                        ),
                        post_reaction_file=DataFiles(
                            data_file=output_dir / f"post{rid}_typed_IFF_merged.data",
                            lmp_molecule_file=output_dir / f"post{rid}_typed_IFF_merged.lmpmol"
                        )
                    ))

        final_files = FFFiles(
            force_field_data=output_dir / "force_field.data",
            in_file=self.cache_all2lmp / "in.create_atoms.script",
            molecule_files=molecule_files,
            template_files=template_files
        )

        FFValidator(final_files)

        return final_files