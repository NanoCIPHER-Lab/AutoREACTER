"""Shared parameterization workflow for AutoREACTER.

Coordinates optimized RDKit molecules, reaction templates, backend
parameterization, and LAMMPS output merging.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Iterable

from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import (
    ReactionMetadata,
)
from AutoREACTER.reaction_preparation.ff_wrapper.parameterization_client.lammps_merge_builder import (
    LammpsMergeBuilder,
)

if TYPE_CHECKING:
    from AutoREACTER.session import Session

AVAILABLE_FORCE_FIELDS = {"OPLSAA", "GAFF", "OPENFF"}

@dataclass(slots=True)
class ParameterizationResult:
    """One parameterized structure and its original molecule/template ID."""

    id: str
    molecule: bool
    typed_data_file: Path


class ParameterizationAPIWrapper:
    """Shared parameterization workflow independent of the backend."""

    def __init__(
        self,
        ARX: "Session",
        force_field_name=None,
    ) -> None:
        self.session = ARX
        self.inputs = ARX.inputs
        self.reaction_metadata = ARX.reaction_metadata
        self.prepared_reactions_with_3d_mols = ARX.reaction_metadata
        # Default to Foyer until additional backends are registered.
        if force_field_name in AVAILABLE_FORCE_FIELDS and force_field_name == "OPLSAA":
            from AutoREACTER.reaction_preparation.ff_wrapper.parameterization_client.foyer.foyer_backend import (
                FoyerBackend,
            )
            self.backend = FoyerBackend(self.inputs.force_field)

        elif force_field_name in AVAILABLE_FORCE_FIELDS and force_field_name == "OPENFF":
            from AutoREACTER.reaction_preparation.ff_wrapper.parameterization_client.openff.openff_backend import (
                OpenFFBackend,
            )
            self.backend = OpenFFBackend(self.inputs.force_field)

        # Preserve the original Foyer staging directory.
        self.cache_dir = Path(ARX.staging_dir) / self.backend.name
        self.cache_atom_typing = self.cache_dir / "atom_typing"
        self.cache_molecules = self.cache_atom_typing / "molecules"
        self.cache_templates = self.cache_atom_typing / "templates"
        self.cache_merge_data = self.cache_dir / "merge_data"

        for directory in (
            self.cache_dir,
            self.cache_atom_typing,
            self.cache_molecules,
            self.cache_templates,
            self.cache_merge_data,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        self.atom_typing_results = (
            self._run_atom_type_monomers(self.inputs.monomers)
            + self._atom_type_reactions(self.prepared_reactions_with_3d_mols)
        )

        self.final_parameterization_files = LammpsMergeBuilder(
            self.cache_merge_data
        ).build(self.atom_typing_results)

        # Temporary compatibility with the previous Foyer interface.
        self.final_foyer_files = self.final_parameterization_files

    def _run_atom_type_monomers(
        self, monomers: Iterable[object]
    ) -> list[ParameterizationResult]:
        """Type every active input monomer, using its optimized RDKit Mol."""
        results: list[ParameterizationResult] = []
        for entry in monomers:
            if getattr(entry, "status", True) is False:
                continue
            molecule_id = str(entry.name)
            mol = entry.rdkit_mol
            if mol is None:
                raise ValueError(f"Missing RDKit molecule for '{molecule_id}'.")

            path = self.backend.parameterize(
                mol, molecule_id=molecule_id, output_dir=self.cache_molecules
            )
            results.append(
                ParameterizationResult(
                    id=molecule_id, molecule=True, typed_data_file=path
                )
            )
        return results

    def _atom_type_reactions(
        self, prepared_reactions: list[ReactionMetadata]
    ) -> list[ParameterizationResult]:
        """Type both sides of each active reaction under preN/postN IDs."""
        results: list[ParameterizationResult] = []
        for reaction in prepared_reactions:
            if getattr(reaction, "activity_stats", True) is False:
                continue
            for prefix, mol in (
                ("pre", reaction.reactant_combined_RDmol),
                ("post", reaction.product_combined_RDmol),
            ):
                molecule_id = f"{prefix}{reaction.reaction_id}"
                if mol is None:
                    raise ValueError(
                        f"Missing RDKit molecule for template '{molecule_id}'."
                    )
                path = self.backend.parameterize(
                    mol, molecule_id=molecule_id, output_dir=self.cache_templates
                )
                results.append(
                    ParameterizationResult(
                        id=molecule_id, molecule=False, typed_data_file=path
                    )
                )
        return results

