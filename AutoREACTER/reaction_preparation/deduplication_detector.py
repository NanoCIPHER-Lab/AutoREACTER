from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx
from rdkit import Chem

if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.reaction_processor.reaction_progression import (
        ReactionMetadata,
    )


"""
In here we do not add any coordinations, atom IDs, or bond IDs to the graphs. This is because these
There tempalate can generate from different reactanats but can have same template
How ever the is a possible herdle in here that element can be same but atom types in 
LAMMPS can be different. In future release we can add one more attribute to 
the graph which is atom type in LAMMPS and bond type in LAMMPS. 
This will make sure that we are comparing the same template and not different template with same element.
"""

class DeduplicationDetector:
    """
    Detect duplicate pre/post-reaction graph pairs.

    Coordinates, atom IDs, and bond IDs are ignored during graph
    isomorphism comparison.

    Each pre/post pair is combined into a coupled graph. Correspondence
    edges connect the same atom across the reactant and product states,
    ensuring that one consistent atom mapping must satisfy both graphs.
    """

    NODE_ATTRIBUTE = "atom_label"
    EDGE_ATTRIBUTE = "bond_label"

    LAMMPS_COMPARISON_GROUP = "lammps"
    RDKIT_COMPARISON_GROUP = "rdkit"

    _PRE_PHASE = "pre"
    _POST_PHASE = "post"

    _BOND_RELATIONSHIP = "bond"
    _ATOM_CORRESPONDENCE_RELATIONSHIP = "atom_correspondence"

    _LAMMPS_RELEVANT_SECTIONS = {
        "Types",
        "Bonds",
    }

    _LAMMPS_SECTION_HEADERS = {
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

    def __init__(self) -> None:
        """Initialize independent deduplication caches."""
        self.seen_reactions: dict[str, list[nx.Graph]] = {
            self.LAMMPS_COMPARISON_GROUP: [],
            self.RDKIT_COMPARISON_GROUP: [],
        }

    # ------------------------------------------------------------------
    # Public duplicate-detection API
    # ------------------------------------------------------------------

    def is_duplicate(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
        comparison_group: str,
    ) -> bool:
        """
        Check whether an equivalent pre/post-reaction pair was seen.

        Args:
            pre_template_graph:
                Graph representing the reactant state.
            post_template_graph:
                Graph representing the product state.
            comparison_group:
                Cache group used for the comparison, such as ``lammps``
                or ``rdkit``.

        Returns:
            True when an equivalent reaction pair is already cached.
            Otherwise, caches the reaction pair and returns False.
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

    def compare_graphs(
        self,
        molecule_file_paths: list[str | Path],
    ) -> dict[str, bool]:
        """
        Compare LAMMPS pre/post molecule-template pairs.

        A pre-template filename must contain ``pre``. The corresponding
        post-template filename is determined by replacing the first
        occurrence of ``pre`` with ``post``.

        Args:
            molecule_file_paths:
                Candidate LAMMPS molecule-template file paths.

        Returns:
            A mapping from each processed pre-template path to its
            duplicate status.
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

            status = "Duplicate" if duplicate else "Unique"

            print(
                f"{status} reaction: "
                f"{pre_file_path.name} -> {post_file_path.name}"
            )

        return results

    def compare_graphs_mol(
        self,
        reaction_metadata_items: list[ReactionMetadata],
    ) -> list[ReactionMetadata]:
        """
        Detect duplicate reactions using in-memory RDKit molecules.

        Reactions whose ``activity_stats`` value is already False are
        skipped. Detected duplicates are disabled by setting
        ``activity_stats`` to False.

        Only atoms included in
        ``template_reactant_to_product_mapping`` participate in duplicate
        detection. Product atom indices are relabeled into reactant-index
        space before the pre/post graphs are coupled.

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

            reactant_mol = reaction_metadata.reactant_combined_RDmol
            product_mol = reaction_metadata.product_combined_RDmol

            reactant_to_product_mapping = (
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

            if not reactant_to_product_mapping:
                raise ValueError(
                    f"Reaction {reaction_index} does not contain a "
                    "template_reactant_to_product_mapping; duplicate "
                    "comparison cannot be restricted to template atoms."
                )

            reactant_template_indices = set(
                reactant_to_product_mapping
            )

            product_to_reactant_mapping = {
                product_index: reactant_index
                for reactant_index, product_index
                in reactant_to_product_mapping.items()
            }

            product_template_indices = set(
                product_to_reactant_mapping
            )

            pre_graph = self.rdkit_mol_to_networkx(
                molecule=reactant_mol,
                atom_idxs=reactant_template_indices,
            )

            post_graph = self.rdkit_mol_to_networkx(
                molecule=product_mol,
                atom_idxs=product_template_indices,
                idx_relabel=product_to_reactant_mapping,
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
                Cache group to clear. When None, all cache groups are
                cleared.
        """
        if comparison_group is None:
            for seen_graphs in self.seen_reactions.values():
                seen_graphs.clear()

            return

        self.seen_reactions.setdefault(
            comparison_group,
            [],
        ).clear()

    # ------------------------------------------------------------------
    # Public graph-conversion API
    # ------------------------------------------------------------------

    def rdkit_mol_to_networkx(
        self,
        molecule: Chem.Mol,
        atom_idxs: set[int] | None = None,
        idx_relabel: dict[int, int] | None = None,
    ) -> nx.Graph:
        """
        Convert an in-memory RDKit molecule into a NetworkX graph.

        Coordinates are not read or stored.

        Node attributes:
            atom_label:
                Chemical element symbol.

        Edge attributes:
            bond_label:
                RDKit bond type represented as a string.

        Args:
            molecule:
                RDKit molecule to convert.
            atom_idxs:
                Optional atom-index set defining the induced subgraph.
                Atoms outside this set and bonds touching excluded atoms
                are omitted.
            idx_relabel:
                Optional mapping from RDKit atom indices to graph node
                IDs. This can be used to express a product graph in the
                corresponding reactant-index space.

        Returns:
            NetworkX representation of the selected molecule region.
        """
        if molecule is None:
            raise ValueError(
                "Cannot create a graph from a None RDKit molecule."
            )

        included_atom_indices = (
            atom_idxs
            if atom_idxs is not None
            else {
                atom.GetIdx()
                for atom in molecule.GetAtoms()
            }
        )

        if idx_relabel is not None:
            missing_relabels = sorted(
                included_atom_indices - idx_relabel.keys()
            )

            if missing_relabels:
                raise ValueError(
                    "The atom-index relabel mapping does not contain "
                    f"entries for atom indices {missing_relabels}."
                )

        graph = nx.Graph()

        for atom in molecule.GetAtoms():
            atom_index = atom.GetIdx()

            if atom_index not in included_atom_indices:
                continue

            node_id = self._resolve_node_id(
                atom_index=atom_index,
                idx_relabel=idx_relabel,
            )

            graph.add_node(
                node_id,
                **{
                    self.NODE_ATTRIBUTE: atom.GetSymbol(),
                },
            )

        for bond in molecule.GetBonds():
            atom1_index = bond.GetBeginAtomIdx()
            atom2_index = bond.GetEndAtomIdx()

            if (
                atom1_index not in included_atom_indices
                or atom2_index not in included_atom_indices
            ):
                continue

            node1_id = self._resolve_node_id(
                atom_index=atom1_index,
                idx_relabel=idx_relabel,
            )
            node2_id = self._resolve_node_id(
                atom_index=atom2_index,
                idx_relabel=idx_relabel,
            )

            graph.add_edge(
                node1_id,
                node2_id,
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

        Only the ``Types`` and ``Bonds`` sections are used. Coordinates,
        charges, angles, dihedrals, impropers, atom IDs, and bond IDs do
        not participate in graph-isomorphism comparison.

        Args:
            file_path:
                Path to the LAMMPS molecule-template file.

        Returns:
            NetworkX representation of the molecule template.
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

    # ------------------------------------------------------------------
    # Private coupled-graph helpers
    # ------------------------------------------------------------------

    def _couple_graphs(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
    ) -> nx.Graph:
        """
        Combine reactant and product graphs into one coupled graph.

        Every atom is represented twice:

            ("pre", atom_id)
            ("post", atom_id)

        A correspondence edge connects matching atom IDs across the two
        phases, requiring one consistent mapping for both states.
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
            phase=self._PRE_PHASE,
        )

        self._add_phase_to_coupled_graph(
            source_graph=post_template_graph,
            coupled_graph=coupled_graph,
            phase=self._POST_PHASE,
        )

        for atom_id in pre_atom_ids:
            coupled_graph.add_edge(
                (self._PRE_PHASE, atom_id),
                (self._POST_PHASE, atom_id),
                relationship=(
                    self._ATOM_CORRESPONDENCE_RELATIONSHIP
                ),
                **{
                    self.EDGE_ATTRIBUTE: None,
                },
            )

        return coupled_graph

    def _add_phase_to_coupled_graph(
        self,
        source_graph: nx.Graph,
        coupled_graph: nx.Graph,
        phase: str,
    ) -> None:
        """Add one reaction phase to a coupled pre/post graph."""
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
                **{
                    self.NODE_ATTRIBUTE: atom_label,
                },
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
                relationship=self._BOND_RELATIONSHIP,
                **{
                    self.EDGE_ATTRIBUTE: bond_label,
                },
            )

    @staticmethod
    def _resolve_node_id(
        atom_index: int,
        idx_relabel: dict[int, int] | None,
    ) -> int:
        """Resolve an RDKit atom index to its output graph node ID."""
        if idx_relabel is None:
            return atom_index

        return idx_relabel[atom_index]

    # ------------------------------------------------------------------
    # Private LAMMPS parsing helpers
    # ------------------------------------------------------------------

    def _read_lammps_sections(
        self,
        file_path: Path,
    ) -> dict[str, list[str]]:
        """Read relevant sections from a LAMMPS molecule-template file."""
        sections: dict[str, list[str]] = {}
        current_section: str | None = None

        with file_path.open(
            "r",
            encoding="utf-8",
        ) as file:
            for raw_line in file:
                line = raw_line.split(
                    "#",
                    maxsplit=1,
                )[0].strip()

                if not line:
                    continue

                if line in self._LAMMPS_SECTION_HEADERS:
                    if line in self._LAMMPS_RELEVANT_SECTIONS:
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

    def _add_lammps_atoms(
        self,
        graph: nx.Graph,
        type_lines: list[str],
        file_path: Path,
    ) -> None:
        """Add atoms from a LAMMPS ``Types`` section."""
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
        """Add bonds from a LAMMPS ``Bonds`` section."""
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

    @staticmethod
    def _validate_bond_atoms(
        graph: nx.Graph,
        bond_id: int,
        atom1_id: int,
        atom2_id: int,
        source: Path | str,
    ) -> None:
        """Ensure that both atoms referenced by a bond exist."""
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


def _main() -> None:
    """Run the standalone LAMMPS deduplication example."""
    folder_path = Path(
        "/mnt/c/Users/janit/Documents/GitHub/AutoREACTER/"
        "examples/AutoREACTER_outputs/"
        "Epoxy_Test_Primary_Diamine_Diepoxy/"
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
    results = detector.compare_graphs(pre_template_files)

    print("\nDeduplication results:")

    for file_path, duplicate in results.items():
        status = "duplicate" if duplicate else "unique"

        print(
            f"{Path(file_path).name}: {status}"
        )


if __name__ == "__main__":
    _main()