"""
Graph-based reaction deduplication for RDKit molecules and LAMMPS
molecule templates.

The comparison intentionally ignores coordinates, atom IDs, and bond IDs.
RDKit comparisons use chemical element, radical state, and bond type.
LAMMPS comparisons use atom type and bond type.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx
from rdkit import Chem

if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import (
        ReactionMetadata,
    )


class DeduplicationDetector:
    """Detect duplicate pre/post-reaction graph pairs."""

    DEEP_CHECK = True

    NODE_ATTRIBUTE = "atom_label"
    EDGE_ATTRIBUTE = "bond_label"

    RADICAL_COUNT_ATTRIBUTE = "radical_count"
    RADICAL_PRESENT_ATTRIBUTE = "contains_radical"
    RADICAL_SIGNATURE_ATTRIBUTE = "radical_signature"

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

    _LAMMPS_MAP_SECTION_HEADERS = {
        "InitiatorIDs",
        "EdgeIDs",
        "Equivalences",
        "DeleteIDs",
        "Wildcards",
    }

    def __init__(self) -> None:
        """Initialize independent graph-comparison caches."""
        self.seen_reactions: dict[str, list[nx.Graph]] = {
            self.LAMMPS_COMPARISON_GROUP: [],
            self.RDKIT_COMPARISON_GROUP: [],
        }

        self.seen_reaction_pairs: dict[
            str,
            list[tuple[nx.Graph, nx.Graph]],
        ] = {
            self.LAMMPS_COMPARISON_GROUP: [],
            self.RDKIT_COMPARISON_GROUP: [],
        }

    # ------------------------------------------------------------------
    # Duplicate-detection API
    # ------------------------------------------------------------------

    def is_duplicate(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
        comparison_group: str,
    ) -> bool:
        """
        Check whether an equivalent coupled pre/post pair was previously
        cached.

        The pre- and post-reaction graphs are coupled using atom
        correspondence edges. This requires one consistent atom mapping
        to satisfy both reaction phases.

        This path is mainly used for RDKit molecules where pre and post
        graphs have already been relabeled into the same atom-index space.
        """
        coupled_graph = self._couple_graphs(
            pre_template_graph=pre_template_graph,
            post_template_graph=post_template_graph,
        )

        return self.is_duplicate_coupled_graph(
            coupled_graph=coupled_graph,
            comparison_group=comparison_group,
        )

    def is_duplicate_coupled_graph(
        self,
        coupled_graph: nx.Graph,
        comparison_group: str,
    ) -> bool:
        """
        Check whether an already-coupled pre/post graph was previously cached.

        This is used for LAMMPS templates where atom correspondence comes from
        the RXN_*.map Equivalences section instead of matching atom IDs.
        """
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
            if (
                coupled_graph.graph.get(
                    self.RADICAL_SIGNATURE_ATTRIBUTE
                )
                != seen_graph.graph.get(
                    self.RADICAL_SIGNATURE_ATTRIBUTE
                )
            ):
                continue

            if (
                coupled_graph.number_of_nodes()
                != seen_graph.number_of_nodes()
            ):
                continue

            if (
                coupled_graph.number_of_edges()
                != seen_graph.number_of_edges()
            ):
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

    def is_duplicate_lammps_template_pair(
        self,
        pre_file_path: str | Path,
        post_file_path: str | Path,
        comparison_group: str | None = None,
        map_file_path: str | Path | None = None,
        wildcards: bool = False,
    ) -> bool:
        """
        Check LAMMPS template duplication using:
            - template_pre_*.molecule
            - template_post_*.molecule
            - RXN_*.map Equivalences

        This ignores coordinates, atom IDs, and bond IDs, but preserves
        the pre-to-post atom correspondence from the LAMMPS map file.

        When ``wildcards`` is True, the EdgeIDs from the map file are treated
        as wildcard atoms. Those pre-template atoms and their equivalent
        post-template atoms are removed before graph comparison.
        """
        pre_file_path = Path(pre_file_path)
        post_file_path = Path(post_file_path)

        if comparison_group is None:
            comparison_group = self.LAMMPS_COMPARISON_GROUP

        pre_graph = self.lammps_molecule_to_networkx(pre_file_path)
        post_graph = self.lammps_molecule_to_networkx(post_file_path)

        if map_file_path is None:
            map_file_path = self._lammps_map_path_from_template_path(
                pre_file_path
            )
        else:
            map_file_path = Path(map_file_path)

        if not map_file_path.is_file():
            return self.is_duplicate_pair(
                pre_graph=pre_graph,
                post_graph=post_graph,
                comparison_group=comparison_group,
            )

        pre_to_post_mapping = self._read_lammps_equivalences(
            map_file_path
        )

        if not pre_to_post_mapping:
            raise ValueError(
                f"No Equivalences mapping found in {map_file_path}."
            )

        if wildcards:
            pre_graph, post_graph, pre_to_post_mapping = (
                self._apply_lammps_wildcards(
                    pre_graph=pre_graph,
                    post_graph=post_graph,
                    pre_to_post_mapping=pre_to_post_mapping,
                    map_file_path=map_file_path,
                )
            )

        coupled_graph = self._couple_lammps_graphs(
            pre_graph=pre_graph,
            post_graph=post_graph,
            pre_to_post_mapping=pre_to_post_mapping,
            source=map_file_path,
        )

        return self.is_duplicate_coupled_graph(
            coupled_graph=coupled_graph,
            comparison_group=comparison_group,
        )

    @staticmethod
    def _lammps_reaction_id_from_template_path(
        file_path: str | Path,
    ) -> str:
        """
        Extract reaction ID from names like:
            template_pre_22.molecule
            template_post_22.molecule
            template_pre_1_homo2.molecule
        """
        file_path = Path(file_path)

        match = re.match(
            r"template_(?:pre|post)_(.+)\.molecule$",
            file_path.name,
        )

        if match is None:
            raise ValueError(
                f"Could not infer reaction ID from template file name: "
                f"{file_path.name}"
            )

        return match.group(1)

    def _lammps_map_path_from_template_path(
        self,
        file_path: str | Path,
    ) -> Path:
        """Return the RXN_*.map path matching a template molecule path."""
        file_path = Path(file_path)
        reaction_id = self._lammps_reaction_id_from_template_path(
            file_path
        )
        return file_path.with_name(f"RXN_{reaction_id}.map")

    @classmethod
    def _is_lammps_map_section_header(cls, line: str) -> bool:
        """Return True if a stripped line starts a map-file section."""
        if not line:
            return False

        first_token = line.split()[0]
        return (
            line in cls._LAMMPS_MAP_SECTION_HEADERS
            or first_token in cls._LAMMPS_MAP_SECTION_HEADERS
        )

    @staticmethod
    def _is_lammps_count_line(line: str) -> bool:
        """Return True for count lines such as '11 edgeIDs'."""
        return re.match(
            r"^\s*\d+\s+(?:edgeIDs|equivalences|wildcards)\s*$",
            line,
            flags=re.IGNORECASE,
        ) is not None

    @classmethod
    def _read_lammps_map_integer_section(
        cls,
        map_file_path: str | Path,
        section_name: str,
    ) -> list[list[int]]:
        """Read integer rows from a LAMMPS bond/react map-file section."""
        map_file_path = Path(map_file_path)
        current_section: str | None = None
        rows: list[list[int]] = []

        with map_file_path.open("r", encoding="utf-8") as file:
            for raw_line in file:
                line = raw_line.split("#", maxsplit=1)[0].strip()

                if not line:
                    continue

                first_token = line.split()[0]

                if cls._is_lammps_map_section_header(line):
                    current_section = first_token
                    continue

                if current_section != section_name:
                    continue

                try:
                    rows.append([int(part) for part in line.split()])
                except ValueError:
                    continue

        return rows

    @classmethod
    def _read_lammps_equivalences(
        cls,
        map_file_path: str | Path,
    ) -> dict[int, int]:
        """
        Read pre-to-post atom equivalences from a LAMMPS bond/react map file.

        Returns:
            {pre_atom_id: post_atom_id}
        """
        mapping: dict[int, int] = {}

        for row in cls._read_lammps_map_integer_section(
            map_file_path,
            "Equivalences",
        ):
            if len(row) < 2:
                continue

            pre_atom_id = row[0]
            post_atom_id = row[1]

            if (
                pre_atom_id in mapping
                and mapping[pre_atom_id] != post_atom_id
            ):
                raise ValueError(
                    f"Conflicting Equivalences entry in "
                    f"{map_file_path}: pre atom {pre_atom_id} maps to "
                    f"both {mapping[pre_atom_id]} and {post_atom_id}."
                )

            mapping[pre_atom_id] = post_atom_id

        post_to_pre: dict[int, int] = {}

        for pre_atom_id, post_atom_id in mapping.items():
            if (
                post_atom_id in post_to_pre
                and post_to_pre[post_atom_id] != pre_atom_id
            ):
                raise ValueError(
                    f"Non-bijective Equivalences section in "
                    f"{map_file_path}: post atom {post_atom_id} is mapped "
                    f"from both {post_to_pre[post_atom_id]} and "
                    f"{pre_atom_id}."
                )

            post_to_pre[post_atom_id] = pre_atom_id

        return mapping

    @classmethod
    def _read_lammps_edge_ids(
        cls,
        map_file_path: str | Path,
    ) -> list[int]:
        """Read pre-template EdgeIDs from a LAMMPS bond/react map file."""
        edge_ids: list[int] = []

        for row in cls._read_lammps_map_integer_section(
            map_file_path,
            "EdgeIDs",
        ):
            if not row:
                continue

            edge_ids.append(row[0])

        return edge_ids

    @staticmethod
    def _remove_lammps_map_section(
        lines: list[str],
        section_name: str,
    ) -> list[str]:
        """Remove a named map-file section and its body."""
        section_headers = {
            "InitiatorIDs",
            "EdgeIDs",
            "Equivalences",
            "DeleteIDs",
            "Wildcards",
        }

        output: list[str] = []
        inside_target = False

        for line in lines:
            stripped = line.strip()
            first_token = stripped.split()[0] if stripped else ""

            is_header = (
                stripped in section_headers
                or first_token in section_headers
            )

            if is_header and first_token == section_name:
                inside_target = True
                continue

            if inside_target and is_header:
                inside_target = False

            if not inside_target:
                output.append(line)

        return output

    def _replace_lammps_count_lines(
        self,
        lines: list[str],
        edge_count: int,
        equivalence_count: int,
        wildcard_count: int,
    ) -> list[str]:
        """
        Replace the map-file count block.

        The final order is:
            N edgeIDs
            N equivalences
            N wildcards
        """
        filtered_lines = [
            line
            for line in lines
            if not self._is_lammps_count_line(line.strip())
        ]

        insert_idx = 0
        for idx, line in enumerate(filtered_lines):
            if line.strip().startswith("#") or not line.strip():
                insert_idx = idx + 1
                continue
            break

        count_lines = [
            f"{edge_count} edgeIDs\n",
            f"{equivalence_count} equivalences\n",
            f"{wildcard_count} wildcards\n",
        ]

        return (
            filtered_lines[:insert_idx]
            + count_lines
            + filtered_lines[insert_idx:]
        )

    @staticmethod
    def _read_lammps_single_column_section(
        map_file_path: str | Path,
        section_name: str,
    ) -> list[int]:
        """
        Read a single-column integer section from a LAMMPS bond/react map file.
        """
        map_file_path = Path(map_file_path)

        section_headers = {
            "InitiatorIDs",
            "EdgeIDs",
            "Equivalences",
            "DeleteIDs",
            "Wildcards",
        }

        current_section: str | None = None
        values: list[int] = []

        with map_file_path.open("r", encoding="utf-8") as file:
            for raw_line in file:
                line = raw_line.split("#", maxsplit=1)[0].strip()

                if not line:
                    continue

                first_token = line.split()[0]

                if line in section_headers:
                    current_section = line
                    continue

                if first_token in section_headers:
                    current_section = first_token
                    continue

                if current_section != section_name:
                    continue

                try:
                    values.append(int(first_token))
                except ValueError:
                    continue

        return values

    def _write_lammps_wildcard_map_file(
        self,
        map_file_path: str | Path,
        wildcard_ids: list[int],
    ) -> None:
        """
        Replace the original RXN_*.map with a wildcard-enabled map.

        The map path is unchanged. The Wildcards section is written from the
        pre-template EdgeIDs.
        """
        map_file_path = Path(map_file_path)

        # Preserve input order while removing repeated IDs.
        seen_wildcards: set[int] = set()
        wildcard_ids = [
            atom_id
            for atom_id in wildcard_ids
            if not (
                atom_id in seen_wildcards
                or seen_wildcards.add(atom_id)
            )
        ]

        original_lines = map_file_path.read_text(
            encoding="utf-8",
        ).splitlines()

        header_lines: list[str] = []

        for line in original_lines:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("#"):
                header_lines.append(stripped)
                continue

            break

        initiator_ids = self._read_lammps_single_column_section(
            map_file_path,
            "InitiatorIDs",
        )
        edge_ids = self._read_lammps_edge_ids(map_file_path)
        equivalences = self._read_lammps_equivalences(map_file_path)
        delete_ids = self._read_lammps_single_column_section(
            map_file_path,
            "DeleteIDs",
        )

        lines: list[str] = []

        for header_line in header_lines:
            lines.append(f"{header_line}\n")

        lines.append("\n")
        lines.append(f"{len(edge_ids)} edgeIDs\n")
        lines.append(f"{len(equivalences)} equivalences\n")

        if delete_ids:
            lines.append(f"{len(delete_ids)} deleteIDs\n")

        lines.append(f"{len(wildcard_ids)} wildcards\n")
        lines.append("\n")

        lines.append("InitiatorIDs\n")
        lines.append("\n")
        for atom_id in initiator_ids:
            lines.append(f"{atom_id}\n")
        lines.append("\n")

        lines.append("EdgeIDs\n")
        lines.append("\n")
        for atom_id in edge_ids:
            lines.append(f"{atom_id}\n")
        lines.append("\n")

        lines.append("Equivalences\n")
        lines.append("\n")
        for pre_atom_id, post_atom_id in equivalences.items():
            lines.append(f"{pre_atom_id:<5} {post_atom_id}\n")
        lines.append("\n")

        if delete_ids:
            lines.append("DeleteIDs\n")
            lines.append("\n")
            for atom_id in delete_ids:
                lines.append(f"{atom_id}\n")
            lines.append("\n")

        lines.append("Wildcards\n")
        lines.append("\n")
        for atom_id in wildcard_ids:
            lines.append(f"{atom_id}\n")

        map_file_path.write_text(
            "".join(lines),
            encoding="utf-8",
        )

    @staticmethod
    def _remove_nodes_if_present(
        graph: nx.Graph,
        node_ids: list[int],
    ) -> nx.Graph:
        """Return a graph copy with selected nodes removed."""
        graph_copy = graph.copy()

        graph_copy.remove_nodes_from(
            node_id
            for node_id in node_ids
            if node_id in graph_copy
        )

        return graph_copy

    def _apply_lammps_wildcards(
        self,
        pre_graph: nx.Graph,
        post_graph: nx.Graph,
        pre_to_post_mapping: dict[int, int],
        map_file_path: Path,
    ) -> tuple[nx.Graph, nx.Graph, dict[int, int]]:
        """
        Apply wildcard-style template comparison.

        EdgeIDs are defined in pre-template atom ID space. Their equivalent
        post-template atom IDs are obtained from the Equivalences section.
        These atoms are removed before graph comparison.
        """
        pre_edge_ids = self._read_lammps_edge_ids(map_file_path)

        post_edge_ids = [
            pre_to_post_mapping[pre_atom_id]
            for pre_atom_id in pre_edge_ids
            if pre_atom_id in pre_to_post_mapping
        ]

        pre_graph = self._remove_nodes_if_present(
            pre_graph,
            pre_edge_ids,
        )
        post_graph = self._remove_nodes_if_present(
            post_graph,
            post_edge_ids,
        )

        pre_edge_id_set = set(pre_edge_ids)
        post_edge_id_set = set(post_edge_ids)

        pre_to_post_mapping = {
            pre_atom_id: post_atom_id
            for pre_atom_id, post_atom_id in pre_to_post_mapping.items()
            if pre_atom_id not in pre_edge_id_set
            and post_atom_id not in post_edge_id_set
        }

        return pre_graph, post_graph, pre_to_post_mapping

    def _couple_lammps_graphs(
        self,
        pre_graph: nx.Graph,
        post_graph: nx.Graph,
        pre_to_post_mapping: dict[int, int],
        source: Path | str,
    ) -> nx.Graph:
        """
        Couple LAMMPS pre/post template graphs using RXN_*.map equivalences.

        Unlike the RDKit coupling path, this does not require matching atom IDs
        in pre and post files. The map file defines correspondence.
        """
        coupled_graph = nx.Graph()

        coupled_graph.graph[self.RADICAL_SIGNATURE_ATTRIBUTE] = (
            pre_graph.graph.get(self.RADICAL_COUNT_ATTRIBUTE, 0),
            post_graph.graph.get(self.RADICAL_COUNT_ATTRIBUTE, 0),
            pre_graph.graph.get(self.RADICAL_PRESENT_ATTRIBUTE, False),
            post_graph.graph.get(self.RADICAL_PRESENT_ATTRIBUTE, False),
        )

        self._add_phase_to_coupled_graph(
            source_graph=pre_graph,
            coupled_graph=coupled_graph,
            phase=self._PRE_PHASE,
        )

        self._add_phase_to_coupled_graph(
            source_graph=post_graph,
            coupled_graph=coupled_graph,
            phase=self._POST_PHASE,
        )

        for pre_atom_id, post_atom_id in pre_to_post_mapping.items():
            if pre_atom_id not in pre_graph:
                raise ValueError(
                    f"Map file {source} references pre atom "
                    f"{pre_atom_id}, but that atom is not in the "
                    "pre-template graph."
                )

            if post_atom_id not in post_graph:
                raise ValueError(
                    f"Map file {source} references post atom "
                    f"{post_atom_id}, but that atom is not in the "
                    "post-template graph."
                )

            coupled_graph.add_edge(
                (self._PRE_PHASE, pre_atom_id),
                (self._POST_PHASE, post_atom_id),
                relationship=self._ATOM_CORRESPONDENCE_RELATIONSHIP,
                **{
                    self.EDGE_ATTRIBUTE: None,
                },
            )

        return coupled_graph

    def is_duplicate_pair(
        self,
        pre_graph: nx.Graph,
        post_graph: nx.Graph,
        comparison_group: str,
    ) -> bool:
        """
        Check whether an equivalent uncoupled pre/post graph pair was
        previously cached.

        A reaction is considered a duplicate only when both its reactant
        graph and product graph match the same cached reaction entry.
        """
        node_match = nx.algorithms.isomorphism.categorical_node_match(
            self.NODE_ATTRIBUTE,
            None,
        )

        edge_match = nx.algorithms.isomorphism.categorical_edge_match(
            self.EDGE_ATTRIBUTE,
            None,
        )

        seen_pairs = self.seen_reaction_pairs.setdefault(
            comparison_group,
            [],
        )

        for seen_pre_graph, seen_post_graph in seen_pairs:
            if (
                pre_graph.number_of_nodes()
                != seen_pre_graph.number_of_nodes()
            ):
                continue

            if (
                pre_graph.number_of_edges()
                != seen_pre_graph.number_of_edges()
            ):
                continue

            if (
                post_graph.number_of_nodes()
                != seen_post_graph.number_of_nodes()
            ):
                continue

            if (
                post_graph.number_of_edges()
                != seen_post_graph.number_of_edges()
            ):
                continue

            pre_matches = nx.is_isomorphic(
                pre_graph,
                seen_pre_graph,
                node_match=node_match,
                edge_match=edge_match,
            )

            if not pre_matches:
                continue

            post_matches = nx.is_isomorphic(
                post_graph,
                seen_post_graph,
                node_match=node_match,
                edge_match=edge_match,
            )

            if post_matches:
                return True

        seen_pairs.append(
            (
                pre_graph.copy(),
                post_graph.copy(),
            )
        )

        return False

    def compare_graphs(
        self,
        molecule_file_paths: list[str | Path],
        wildcards: bool = False,
    ) -> dict[str, bool]:
        """
        Compare LAMMPS pre/post molecule-template pairs.

        A pre-template filename must contain ``pre``. Its post-template
        path is determined by replacing the first occurrence of ``pre``
        with ``post``.
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

            duplicate = self.is_duplicate_lammps_template_pair(
                pre_file_path=pre_file_path,
                post_file_path=post_file_path,
                comparison_group=self.LAMMPS_COMPARISON_GROUP,
                wildcards=wildcards,
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
        reaction_metadata_items: list["ReactionMetadata"],
        index_source: str = "template",
        deep_check: bool = True,
    ) -> list["ReactionMetadata"]:
        
        
        """
        Detect duplicate reactions using in-memory RDKit molecules.

        Each call performs one independent deduplication pass over the
        supplied accumulated reaction pool. The RDKit coupled-graph cache
        is therefore cleared before the comparison starts.

        Reactions that are already inactive are ignored. Repeated references
        to the exact same ReactionMetadata object are removed from the
        returned list without disabling the retained object. This matters
        because setting ``activity_stats`` to False on one repeated reference
        would otherwise disable every occurrence of that same object.

        For distinct ReactionMetadata objects, the first unique reaction is
        retained. Later equivalent reactions are disabled by setting
        ``activity_stats`` to False and are excluded from the returned pool.
        """
        self.clear_cache(self.RDKIT_COMPARISON_GROUP)

        unique_reactions: list["ReactionMetadata"] = []
        retained_object_ids: set[int] = set()

        for reaction_index, reaction_metadata in enumerate(
            reaction_metadata_items,
            start=1,
        ):
            if not reaction_metadata.activity_stats:
                continue

            reaction_object_id = id(reaction_metadata)

            if reaction_object_id in retained_object_ids:
                continue

            reactant_mol = reaction_metadata.reactant_combined_RDmol
            product_mol = reaction_metadata.product_combined_RDmol

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

            reactant_to_product_mapping = (
                self._select_reactant_to_product_mapping(
                    reaction_metadata=reaction_metadata,
                    reaction_index=reaction_index,
                    index_source=index_source,
                )
            )

            reactant_indices = set(reactant_to_product_mapping)
            product_indices = set(
                reactant_to_product_mapping.values()
            )

            product_to_reactant_mapping = {
                product_idx: reactant_idx
                for reactant_idx, product_idx
                in reactant_to_product_mapping.items()
            }

            if len(product_to_reactant_mapping) != len(
                reactant_to_product_mapping
            ):
                raise ValueError(
                    f"Reaction {reaction_index} contains a non-bijective "
                    "reactant-to-product mapping."
                )

            pre_graph = self.rdkit_mol_to_networkx(
                molecule=reactant_mol,
                atom_idxs=reactant_indices,
                deep_check=deep_check,
            )

            post_graph = self.rdkit_mol_to_networkx(
                molecule=product_mol,
                atom_idxs=product_indices,
                idx_relabel=product_to_reactant_mapping,
                deep_check=deep_check
            )

            reactant_radical_count = self._count_radical_atoms(
                reactant_mol
            )
            product_radical_count = self._count_radical_atoms(
                product_mol
            )

            pre_graph.graph[self.RADICAL_COUNT_ATTRIBUTE] = (
                reactant_radical_count
            )
            post_graph.graph[self.RADICAL_COUNT_ATTRIBUTE] = (
                product_radical_count
            )

            pre_graph.graph[self.RADICAL_PRESENT_ATTRIBUTE] = (
                reactant_radical_count > 0
            )
            post_graph.graph[self.RADICAL_PRESENT_ATTRIBUTE] = (
                product_radical_count > 0
                or bool(getattr(reaction_metadata, "is_radical", False))
            )

            duplicate = self.is_duplicate(
                pre_template_graph=pre_graph,
                post_template_graph=post_graph,
                comparison_group=self.RDKIT_COMPARISON_GROUP,
            )

            if duplicate:
                reaction_metadata.activity_stats = False
                continue

            retained_object_ids.add(reaction_object_id)
            unique_reactions.append(reaction_metadata)

        return unique_reactions

    @classmethod
    def _one_neighbor_edge_environment_signature(
        cls,
        atom: Chem.Atom,
        included_atom_indices: set[int],
    ) -> tuple[tuple[str, int, bool, str, int, bool, str], ...]:
        """
        Return a one-bond external chemical-environment signature for
        boundary atoms.

        Only neighbors outside the restricted comparison graph are used.
        This avoids walking into the next molecule or building a larger
        shell. Internal atoms return an empty signature.
        """
        external_neighbor_signatures = []

        try:
            atom.GetOwningMol().UpdatePropertyCache(strict=False)
        except RuntimeError:
            pass

        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)

            if neighbor.GetIdx() in included_atom_indices:
                continue

            external_neighbor_signatures.append(
                (
                    neighbor.GetSymbol(),
                    neighbor.GetFormalCharge(),
                    neighbor.GetIsAromatic(),
                    str(neighbor.GetHybridization()),
                    cls._safe_total_hydrogen_count(neighbor),
                    cls._is_radical_atom(neighbor),
                    str(bond.GetBondType()),
                )
            )

        return tuple(sorted(external_neighbor_signatures))

    @staticmethod
    def _safe_total_hydrogen_count(atom: Chem.Atom) -> int:
        """Return total hydrogen count without failing on unsanitized mols."""
        try:
            return atom.GetTotalNumHs()
        except RuntimeError:
            explicit_h_neighbors = sum(
                neighbor.GetAtomicNum() == 1
                for neighbor in atom.GetNeighbors()
            )
            return explicit_h_neighbors + atom.GetNumExplicitHs()

    def clear_cache(
        self,
        comparison_group: str | None = None,
    ) -> None:
        """
        Clear graph-comparison caches.

        Args:
            comparison_group:
                Specific cache group to clear. All groups are cleared
                when omitted.
        """
        if comparison_group is None:
            for seen_graphs in self.seen_reactions.values():
                seen_graphs.clear()

            for seen_pairs in self.seen_reaction_pairs.values():
                seen_pairs.clear()

            return

        self.seen_reactions.setdefault(
            comparison_group,
            [],
        ).clear()

        self.seen_reaction_pairs.setdefault(
            comparison_group,
            [],
        ).clear()

    # ------------------------------------------------------------------
    # RDKit graph conversion
    # ------------------------------------------------------------------

    def rdkit_mol_to_networkx(
        self,
        molecule: Chem.Mol,
        atom_idxs: set[int] | None = None,
        idx_relabel: dict[int, int] | None = None,
        deep_check: bool = True,
    ) -> nx.Graph:
        """
        Convert an RDKit molecule into a NetworkX graph.

        Coordinates are not read or stored.
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

            is_radical = self._is_radical_atom(atom)

            
            if not deep_check:
                atom_label = (
                    atom.GetSymbol(),
                    is_radical,
                )
            else:
                atom_label = (
                    atom.GetSymbol(),
                    is_radical,
                    self._one_neighbor_edge_environment_signature(
                        atom=atom,
                        included_atom_indices=included_atom_indices,
                    ),
                )

            graph.add_node(
                node_id,
                **{
                    self.NODE_ATTRIBUTE: atom_label,
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

    # ------------------------------------------------------------------
    # LAMMPS graph conversion
    # ------------------------------------------------------------------

    def lammps_molecule_to_networkx(
        self,
        file_path: str | Path,
    ) -> nx.Graph:
        """
        Convert a LAMMPS molecule-template file into a NetworkX graph.

        Only the ``Types`` and ``Bonds`` sections are included.
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

        return graph

    def write_wildcard_maps(
        self,
        template_files: list["ReactionMetadata"],
    ) -> list["ReactionMetadata"]:
        """
        Write Wildcards sections for active LAMMPS templates without
        performing duplicate filtering.
        """
        active_templates: list["ReactionMetadata"] = []

        for template in template_files:
            if not template.activity_stats:
                continue

            if template.map_file is None:
                template.activity_stats = False
                continue

            map_file_path = Path(template.map_file)

            if not map_file_path.is_file():
                template.activity_stats = False
                continue

            wildcard_ids = self._read_lammps_edge_ids(map_file_path)

            self._write_lammps_wildcard_map_file(
                map_file_path=map_file_path,
                wildcard_ids=wildcard_ids,
            )

            template.map_file = map_file_path
            active_templates.append(template)

        return active_templates

    def compare_lammps_templates(
        self,
        template_files: list["ReactionMetadata"],
        wildcards: bool = False,
    ) -> list["ReactionMetadata"]:
        """
        Compare LAMMPS pre/post molecule-template pairs.

        Filters a list of ReactionMetadata objects by generating LAMMPS graphs
        for their pre and post reaction files, detecting duplicates, setting
        duplicate reactions inactive, and returning only active unique templates.
        """
        self.clear_cache(self.LAMMPS_COMPARISON_GROUP)

        unique_templates: list["ReactionMetadata"] = []

        for template in template_files:
            if not template.activity_stats:
                continue

            if (
                template.pre_reaction_file is None
                or template.post_reaction_file is None
            ):
                print(
                    f"Skipping template ID {template.reaction_id}: "
                    "Missing pre or post reaction file definitions."
                )
                template.activity_stats = False
                continue

            pre_file_path = Path(template.pre_reaction_file)
            post_file_path = Path(template.post_reaction_file)

            if not pre_file_path.is_file():
                print(
                    f"Skipping template ID {template.reaction_id}: "
                    f"Pre-template file does not exist: {pre_file_path}"
                )
                template.activity_stats = False
                continue

            if not post_file_path.is_file():
                print(
                    f"Skipping template ID {template.reaction_id}: "
                    f"Post-template file does not exist: {post_file_path}"
                )
                template.activity_stats = False
                continue

            map_file_value = getattr(template, "map_file", None)
            map_file_path = (
                Path(map_file_value)
                if map_file_value is not None
                else self._lammps_map_path_from_template_path(
                    pre_file_path
                )
            )

            duplicate = self.is_duplicate_lammps_template_pair(
                pre_file_path=pre_file_path,
                post_file_path=post_file_path,
                comparison_group=self.LAMMPS_COMPARISON_GROUP,
                map_file_path=map_file_path,
                wildcards=wildcards,
            )

            if duplicate:
                template.activity_stats = False
                print(
                    f"Duplicate template disabled: "
                    f"RXN_{template.reaction_id}"
                )
                continue

            if wildcards and map_file_path.is_file():
                wildcard_ids = self._read_lammps_edge_ids(map_file_path)

                self._write_lammps_wildcard_map_file(
                    map_file_path=map_file_path,
                    wildcard_ids=wildcard_ids,
                )

                template.map_file = map_file_path

            template.activity_stats = True
            unique_templates.append(template)

        return unique_templates

    # ------------------------------------------------------------------
    # Reaction-index selection
    # ------------------------------------------------------------------

    def _select_reactant_to_product_mapping(
        self,
        reaction_metadata: "ReactionMetadata",
        reaction_index: int,
        index_source: str,
    ) -> dict[int, int]:
        """Select the atom-index mapping used for graph restriction."""
        if index_source == "template":
            mapping = (
                reaction_metadata.template_reactant_to_product_mapping
            )

            if not mapping:
                raise ValueError(
                    f"Reaction {reaction_index} does not contain a "
                    "template_reactant_to_product_mapping."
                )

            return mapping

        if index_source == "first_shell":
            first_shell_indices = reaction_metadata.first_shell
            full_mapping = (
                reaction_metadata.reactant_to_product_mapping
            )

            if not first_shell_indices:
                raise ValueError(
                    f"Reaction {reaction_index} does not contain "
                    "first_shell indices."
                )

            if not full_mapping:
                raise ValueError(
                    f"Reaction {reaction_index} does not contain a "
                    "reactant_to_product_mapping."
                )

            return {
                reactant_index: full_mapping[reactant_index]
                for reactant_index in first_shell_indices
                if reactant_index in full_mapping
            }

        raise ValueError(
            f"Unsupported index_source {index_source!r}. "
            "Expected 'template' or 'first_shell'."
        )

    # ------------------------------------------------------------------
    # Coupled-graph helpers
    # ------------------------------------------------------------------

    def _couple_graphs(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
    ) -> nx.Graph:
        """
        Combine reactant and product graphs using atom-correspondence
        edges.

        This RDKit path expects pre and post graph atom IDs to already match.
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

        coupled_graph.graph[self.RADICAL_SIGNATURE_ATTRIBUTE] = (
            pre_template_graph.graph.get(
                self.RADICAL_COUNT_ATTRIBUTE,
                0,
            ),
            post_template_graph.graph.get(
                self.RADICAL_COUNT_ATTRIBUTE,
                0,
            ),
            pre_template_graph.graph.get(
                self.RADICAL_PRESENT_ATTRIBUTE,
                False,
            ),
            post_template_graph.graph.get(
                self.RADICAL_PRESENT_ATTRIBUTE,
                False,
            ),
        )

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
                relationship=self._ATOM_CORRESPONDENCE_RELATIONSHIP,
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
        """Add one reaction phase to a coupled graph."""
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
        """Resolve an RDKit atom index to its graph node ID."""
        if idx_relabel is None:
            return atom_index

        return idx_relabel[atom_index]

    # ------------------------------------------------------------------
    # LAMMPS parsing helpers
    # ------------------------------------------------------------------

    def _read_lammps_sections(
        self,
        file_path: Path,
    ) -> dict[str, list[str]]:
        """Read relevant sections from a LAMMPS molecule file."""
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
        """Ensure both atoms referenced by a bond exist."""
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

    @classmethod
    def _count_radical_atoms(
        cls,
        molecule: Chem.Mol,
    ) -> int:
        """Count radical atoms in a complete RDKit molecule."""
        return sum(
            cls._is_radical_atom(atom)
            for atom in molecule.GetAtoms()
        )

    @staticmethod
    def _is_radical_atom(atom: Chem.Atom) -> bool:
        """Return True for an explicit or structurally under-valent radical atom.

        Progression normally sets ``NumRadicalElectrons`` on the cleaned
        product. The structural fallback is needed because the original
        ``RunReactants`` product stored in ``ReactionMetadata`` can retain
        incomplete valence bookkeeping before later sanitization.

        The fallback is intentionally limited to neutral, non-aromatic carbon
        atoms used by the current vinyl-radical implementation.
        """
        if atom.GetNumRadicalElectrons() > 0:
            return True

        if atom.GetAtomicNum() != 6:
            return False

        if atom.GetFormalCharge() != 0:
            return False

        if atom.GetIsAromatic():
            return False

        try:
            atom.GetOwningMol().UpdatePropertyCache(strict=False)
            if atom.GetNumImplicitHs() > 0:
                return False
        except RuntimeError:
            pass

        graph_bond_valence = sum(
            bond.GetBondTypeAsDouble()
            for bond in atom.GetBonds()
        )

        explicit_hydrogen_neighbors = sum(
            neighbor.GetAtomicNum() == 1
            for neighbor in atom.GetNeighbors()
        )

        effective_valence = graph_bond_valence

        if explicit_hydrogen_neighbors == 0:
            effective_valence += atom.GetNumExplicitHs()

        return abs(effective_valence - 3.0) < 1.0e-6


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

        print(f"{Path(file_path).name}: {status}")


if __name__ == "__main__":
    _main()
