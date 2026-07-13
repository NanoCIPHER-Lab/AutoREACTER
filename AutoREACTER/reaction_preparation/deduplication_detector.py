from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx
from rdkit import Chem

if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.reaction_processor.reaction_progression import (
        ReactionMetadata,
    )


class DeduplicationDetector:
    """
    Detect duplicate pre/post reaction graph pairs.

    Coordinates, atom IDs, and bond IDs are ignored during graph
    isomorphism comparison.

    A coupled graph is created for every reaction pair so that the same
    atom mapping must be valid for both the reactant and product graphs.
    """

    NODE_ATTRIBUTE = "atom_label"
    EDGE_ATTRIBUTE = "bond_label"

    LAMMPS_COMPARISON_GROUP = "lammps"
    RDKIT_COMPARISON_GROUP = "rdkit"

    def __init__(self) -> None:
        """
        Initialize independent deduplication caches.
        """
        self.seen_reactions: dict[str, list[nx.Graph]] = {
            self.LAMMPS_COMPARISON_GROUP: [],
            self.RDKIT_COMPARISON_GROUP: [],
        }

    def is_duplicate(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
        comparison_group: str,
    ) -> bool:
        """
        Check whether an equivalent pre/post reaction pair has been seen.

        Args:
            pre_template_graph:
                Graph representing the reactant state.

            post_template_graph:
                Graph representing the product state.

            comparison_group:
                Cache group used for the comparison, such as ``lammps``
                or ``rdkit``.

        Returns:
            True if an equivalent reaction pair has already been seen.
            Otherwise, stores the reaction pair and returns False.
        """
        coupled_graph = self._couple_graphs(
            pre_template_graph=pre_template_graph,
            post_template_graph=post_template_graph,
        )

        node_match = nx.algorithms.isomorphism.categorical_node_match(
            ["phase", self.NODE_ATTRIBUTE],
            [None, None],
        )

        edge_match = nx.algorithms.isomorphism.categorical_edge_match(
            ["relationship", self.EDGE_ATTRIBUTE],
            [None, None],
        )

        seen_graphs = self.seen_reactions.setdefault(
            comparison_group,
            [],
        )

        for seen_graph in seen_graphs:
            # Cheap checks before running graph isomorphism.
            if coupled_graph.number_of_nodes() != seen_graph.number_of_nodes():
                continue

            if coupled_graph.number_of_edges() != seen_graph.number_of_edges():
                continue

            if nx.is_isomorphic(
                coupled_graph,
                seen_graph,
                node_match=node_match,
                edge_match=edge_match,
            ):
                return True

        seen_graphs.append(coupled_graph.copy())

        return False

    def _couple_graphs(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
    ) -> nx.Graph:
        """
        Combine the reactant and product graphs into one graph.

        Every atom is represented twice:

            ("pre", atom_id)
            ("post", atom_id)

        A correspondence edge connects the same atom ID in the pre- and
        post-reaction graphs. This requires one consistent atom mapping
        across both reaction states.
        """
        pre_atom_ids = set(pre_template_graph.nodes)
        post_atom_ids = set(post_template_graph.nodes)

        if pre_atom_ids != post_atom_ids:
            missing_from_post = sorted(
                pre_atom_ids - post_atom_ids
            )
            missing_from_pre = sorted(
                post_atom_ids - pre_atom_ids
            )

            raise ValueError(
                "Pre- and post-reaction graphs must contain matching "
                "atom IDs. "
                f"Missing from post graph: {missing_from_post}. "
                f"Missing from pre graph: {missing_from_pre}."
            )

        coupled_graph = nx.Graph()

        self._add_phase_to_coupled_graph(
            source_graph=pre_template_graph,
            coupled_graph=coupled_graph,
            phase="pre",
        )

        self._add_phase_to_coupled_graph(
            source_graph=post_template_graph,
            coupled_graph=coupled_graph,
            phase="post",
        )

        for atom_id in pre_atom_ids:
            coupled_graph.add_edge(
                ("pre", atom_id),
                ("post", atom_id),
                relationship="atom_correspondence",
                **{self.EDGE_ATTRIBUTE: None},
            )

        return coupled_graph

    def _add_phase_to_coupled_graph(
        self,
        source_graph: nx.Graph,
        coupled_graph: nx.Graph,
        phase: str,
    ) -> None:
        """
        Add one reaction phase to a coupled pre/post graph.
        """
        for atom_id, attributes in source_graph.nodes(data=True):
            atom_label = attributes.get(self.NODE_ATTRIBUTE)

            if atom_label is None:
                raise ValueError(
                    f"Node {atom_id} is missing the required "
                    f"{self.NODE_ATTRIBUTE!r} attribute."
                )

            coupled_graph.add_node(
                (phase, atom_id),
                phase=phase,
                **{self.NODE_ATTRIBUTE: atom_label},
            )

        for atom1_id, atom2_id, attributes in source_graph.edges(
            data=True
        ):
            bond_label = attributes.get(self.EDGE_ATTRIBUTE)

            if bond_label is None:
                raise ValueError(
                    f"Edge {atom1_id}-{atom2_id} is missing the required "
                    f"{self.EDGE_ATTRIBUTE!r} attribute."
                )

            coupled_graph.add_edge(
                (phase, atom1_id),
                (phase, atom2_id),
                relationship="bond",
                **{self.EDGE_ATTRIBUTE: bond_label},
            )

    def rdkit_mol_to_networkx(
        self,
        molecule: Chem.Mol,
    ) -> nx.Graph:
        """
        Convert an in-memory RDKit molecule into a NetworkX graph.

        Coordinates are not read or stored.

        Node attribute:
            atom_label:
                Chemical element symbol.

        Edge attribute:
            bond_label:
                RDKit bond type.
        """
        if molecule is None:
            raise ValueError(
                "Cannot create a graph from a None RDKit molecule."
            )

        graph = nx.Graph()

        for atom in molecule.GetAtoms():
            atom_id = atom.GetIdx()

            graph.add_node(
                atom_id,
                **{
                    self.NODE_ATTRIBUTE: atom.GetSymbol(),
                },
            )

        for bond in molecule.GetBonds():
            atom1_id = bond.GetBeginAtomIdx()
            atom2_id = bond.GetEndAtomIdx()

            graph.add_edge(
                atom1_id,
                atom2_id,
                **{
                    self.EDGE_ATTRIBUTE: str(
                        bond.GetBondType()
                    ),
                },
            )

        return graph

    def lammps_molecule_to_networkx(
        self,
        file_path: str | Path,
    ) -> nx.Graph:
        """
        Convert a LAMMPS molecule-template file into a NetworkX graph.

        Only the following sections are used:

            Types
            Bonds

        Coordinates, charges, angles, dihedrals, impropers, atom IDs,
        and bond IDs are not used in graph isomorphism comparison.
        """
        file_path = Path(file_path)

        if not file_path.is_file():
            raise FileNotFoundError(
                f"LAMMPS molecule file does not exist: {file_path}"
            )

        sections = self._read_lammps_sections(file_path)

        if "Types" not in sections:
            raise ValueError(
                f"Types section was not found in {file_path}."
            )

        graph = nx.Graph()

        self._add_lammps_atoms(
            graph=graph,
            type_lines=sections["Types"],
            file_path=file_path,
        )

        self._add_lammps_bonds(
            graph=graph,
            bond_lines=sections.get("Bonds", []),
            file_path=file_path,
        )

        print(
            f"Graph created from {file_path.name}: "
            f"{graph.number_of_nodes()} atoms and "
            f"{graph.number_of_edges()} bonds."
        )

        return graph

    def _add_lammps_atoms(
        self,
        graph: nx.Graph,
        type_lines: list[str],
        file_path: Path,
    ) -> None:
        """
        Add atoms from a LAMMPS Types section to a graph.
        """
        for line in type_lines:
            parts = line.split()

            if len(parts) < 2:
                raise ValueError(
                    f"Invalid Types line in {file_path}: {line!r}"
                )

            try:
                atom_id = int(parts[0])
            except ValueError as error:
                raise ValueError(
                    f"Invalid atom ID in {file_path}: {line!r}"
                ) from error

            atom_type = parts[1]

            if atom_id in graph:
                raise ValueError(
                    f"Duplicate atom ID {atom_id} in {file_path}."
                )

            graph.add_node(
                atom_id,
                **{
                    self.NODE_ATTRIBUTE: atom_type,
                },
            )

    def _add_lammps_bonds(
        self,
        graph: nx.Graph,
        bond_lines: list[str],
        file_path: Path,
    ) -> None:
        """
        Add bonds from a LAMMPS Bonds section to a graph.
        """
        for line in bond_lines:
            parts = line.split()

            if len(parts) < 4:
                raise ValueError(
                    f"Invalid Bonds line in {file_path}: {line!r}"
                )

            try:
                bond_id = int(parts[0])
                atom1_id = int(parts[2])
                atom2_id = int(parts[3])
            except ValueError as error:
                raise ValueError(
                    f"Invalid Bonds line in {file_path}: {line!r}"
                ) from error

            bond_type = parts[1]

            self._validate_bond_atoms(
                graph=graph,
                bond_id=bond_id,
                atom1_id=atom1_id,
                atom2_id=atom2_id,
                source=file_path,
            )

            graph.add_edge(
                atom1_id,
                atom2_id,
                **{
                    self.EDGE_ATTRIBUTE: bond_type,
                },
            )

    def _read_lammps_sections(
        self,
        file_path: Path,
    ) -> dict[str, list[str]]:
        """
        Read relevant sections from a LAMMPS molecule-template file.
        """
        relevant_sections = {
            "Types",
            "Bonds",
        }

        all_section_headers = {
            "Coords",
            "Types",
            "Charges",
            "Molecules",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
            "Special Bond Counts",
            "Special Bonds",
        }

        sections: dict[str, list[str]] = {}
        current_section: str | None = None

        with file_path.open(
            "r",
            encoding="utf-8",
        ) as file:
            for raw_line in file:
                # Remove comments while preserving the actual data.
                line = raw_line.split(
                    "#",
                    maxsplit=1,
                )[0].strip()

                if not line:
                    continue

                if line in all_section_headers:
                    if line in relevant_sections:
                        current_section = line
                        sections.setdefault(
                            current_section,
                            [],
                        )
                    else:
                        current_section = None

                    continue

                if current_section is not None:
                    sections[current_section].append(line)

        return sections

    def _validate_bond_atoms(
        self,
        graph: nx.Graph,
        bond_id: int,
        atom1_id: int,
        atom2_id: int,
        source: Path | str,
    ) -> None:
        """
        Ensure that both atoms referenced by a bond exist.
        """
        undefined_atoms = [
            atom_id
            for atom_id in (atom1_id, atom2_id)
            if atom_id not in graph
        ]

        if undefined_atoms:
            raise ValueError(
                f"Bond {bond_id} references undefined atom IDs "
                f"{undefined_atoms} in {source}."
            )

    def compare_graphs(
        self,
        molecule_file_paths: list[str | Path],
    ) -> dict[str, bool]:
        """
        Compare LAMMPS pre/post molecule-template pairs.

        A pre-template filename must contain ``pre``. Its corresponding
        post-template filename is found by replacing the first occurrence
        of ``pre`` with ``post``.

        Returns:
            Mapping from each pre-template path to its duplicate status.
        """
        results: dict[str, bool] = {}

        for file_path_value in molecule_file_paths:
            pre_file_path = Path(file_path_value)

            if "pre" not in pre_file_path.name:
                continue

            post_file_path = pre_file_path.with_name(
                pre_file_path.name.replace(
                    "pre",
                    "post",
                    1,
                )
            )

            if not post_file_path.is_file():
                print(
                    "Skipping reaction because its post-template file "
                    f"does not exist: {post_file_path}"
                )
                continue

            pre_graph = self.lammps_molecule_to_networkx(
                pre_file_path
            )

            post_graph = self.lammps_molecule_to_networkx(
                post_file_path
            )

            duplicate = self.is_duplicate(
                pre_template_graph=pre_graph,
                post_template_graph=post_graph,
                comparison_group=self.LAMMPS_COMPARISON_GROUP,
            )

            results[str(pre_file_path)] = duplicate

            status = (
                "Duplicate"
                if duplicate
                else "Unique"
            )

            print(
                f"{status} reaction: "
                f"{pre_file_path.name} -> "
                f"{post_file_path.name}"
            )

        return results

    def compare_graphs_mol(
        self,
        reaction_metadata_items: list[ReactionMetadata],
    ) -> list[ReactionMetadata]:
        """
        Detect duplicate reactions using in-memory RDKit molecules.

        Reactions whose ``activity_stats`` value is already False are
        skipped.

        Duplicate reactions are disabled by setting:

            reaction_metadata.activity_stats = False

        Args:
            reaction_metadata_items:
                Prepared reaction metadata objects.

        Returns:
            The original metadata list with duplicate reactions disabled.
        """
        for reaction_index, reaction_metadata in enumerate(
            reaction_metadata_items,
            start=1,
        ):
            if reaction_metadata.activity_stats is False:
                continue

            reactant_mol = (
                reaction_metadata.reactant_combined_RDmol
            )
            product_mol = (
                reaction_metadata.product_combined_RDmol
            )
            template_reactant_to_product_mapping = (
                reaction_metadata.template_reactant_to_product_mapping
            )
            if reactant_mol is None:
                raise ValueError(
                    f"Reaction {reaction_index} does not contain a "
                    "combined reactant RDKit molecule."
                )

            if product_mol is None:
                raise ValueError(
                    f"Reaction {reaction_index} does not contain a "
                    "combined product RDKit molecule."
                )

            pre_graph = self.rdkit_mol_to_networkx(
                reactant_mol
            )

            post_graph = self.rdkit_mol_to_networkx(
                product_mol
            )

            duplicate = self.is_duplicate(
                pre_template_graph=pre_graph,
                post_template_graph=post_graph,
                comparison_group=self.RDKIT_COMPARISON_GROUP,
            )

            if duplicate:
                reaction_metadata.activity_stats = False

                print(
                    f"Reaction {reaction_index}: "
                    "duplicate reaction detected and disabled."
                )
            else:
                print(
                    f"Reaction {reaction_index}: "
                    "unique reaction retained."
                )

        return reaction_metadata_items

    def clear_cache(
        self,
        comparison_group: str | None = None,
    ) -> None:
        """
        Clear stored reaction graphs.

        Args:
            comparison_group:
                Clear only one cache group. When None, all cache groups
                are cleared.
        """
        if comparison_group is None:
            for seen_graphs in self.seen_reactions.values():
                seen_graphs.clear()

            return

        self.seen_reactions.setdefault(
            comparison_group,
            [],
        ).clear()


if __name__ == "__main__":
    folder_path = Path(
        "/mnt/c/Users/janit/Documents/GitHub/AutoREACTER/"
        "examples/AutoREACTER_outputs/"
        "Epoxy_Test_Primary_Diamine_Diepoxy/"
        "LAMMPS_input_files/"
        "Epoxy_Test_Primary_Diamine_Diepoxy_epoxy_test"
    )

    if not folder_path.is_dir():
        raise NotADirectoryError(
            f"Invalid folder path: {folder_path}"
        )

    pre_template_files = sorted(
        file_path
        for file_path in folder_path.glob("*.molecule")
        if "pre" in file_path.name
    )

    print(
        f"Found {len(pre_template_files)} "
        "pre-reaction molecule files."
    )

    detector = DeduplicationDetector()

    results = detector.compare_graphs(
        pre_template_files
    )

    print("\nDeduplication results:")

    for file_path, duplicate in results.items():
        status = (
            "duplicate"
            if duplicate
            else "unique"
        )

        print(
            f"{Path(file_path).name}: {status}"
        )