import os
from pathlib import Path

import networkx as nx


class DeduplicationDetector:
    def __init__(self):
        # Graphs are not hashable, so store reaction graph pairs in a list.
        self.seen_reactions: list[tuple[nx.Graph, nx.Graph]] = []

    def is_duplicate(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
    ) -> bool:
        """
        Check whether a pre/post reaction graph pair has already been seen.

        Atom types and bond types are included in the graph comparison.
        Coordinates are intentionally ignored because the same molecular
        topology can have different coordinates.
        """
        node_match = nx.algorithms.isomorphism.categorical_node_match(
            "atom_type",
            None,
        )
        edge_match = nx.algorithms.isomorphism.categorical_edge_match(
            "bond_type",
            None,
        )

        for seen_pre_graph, seen_post_graph in self.seen_reactions:
            pre_matches = nx.is_isomorphic(
                pre_template_graph,
                seen_pre_graph,
                node_match=node_match,
                edge_match=edge_match,
            )

            if not pre_matches:
                continue

            post_matches = nx.is_isomorphic(
                post_template_graph,
                seen_post_graph,
                node_match=node_match,
                edge_match=edge_match,
            )

            if post_matches:
                return True

        self.seen_reactions.append(
            (
                pre_template_graph.copy(),
                post_template_graph.copy(),
            )
        )

        return False

    def lammps_molecule_to_networkx(
        self,
        file_path: str | Path,
    ) -> nx.Graph:
        """
        Convert a LAMMPS molecule-template file into a NetworkX graph.

        Supported sections:
            Types
            Coords
            Bonds
        """
        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(
                f"Molecule file does not exist: {file_path}"
            )

        sections = self._read_sections(file_path)

        atom_mapping: dict[int, str] = {}
        coord_mapping: dict[int, tuple[float, float, float]] = {}
        bond_mapping: dict[int, tuple[str, int, int]] = {}

        for line in sections.get("Types", []):
            atom_data = line.split()

            if len(atom_data) < 2:
                raise ValueError(
                    f"Invalid Types line in {file_path}: {line!r}"
                )

            atom_id = int(atom_data[0])
            atom_type = atom_data[1]
            atom_mapping[atom_id] = atom_type

        for line in sections.get("Coords", []):
            coord_data = line.split()

            if len(coord_data) < 4:
                raise ValueError(
                    f"Invalid Coords line in {file_path}: {line!r}"
                )

            atom_id = int(coord_data[0])
            x = float(coord_data[1])
            y = float(coord_data[2])
            z = float(coord_data[3])

            coord_mapping[atom_id] = (x, y, z)

        for line in sections.get("Bonds", []):
            bond_data = line.split()

            if len(bond_data) < 4:
                raise ValueError(
                    f"Invalid Bonds line in {file_path}: {line!r}"
                )

            bond_id = int(bond_data[0])
            bond_type = bond_data[1]
            atom1_id = int(bond_data[2])
            atom2_id = int(bond_data[3])

            bond_mapping[bond_id] = (
                bond_type,
                atom1_id,
                atom2_id,
            )

        graph = nx.Graph()

        for atom_id, atom_type in atom_mapping.items():
            if atom_id not in coord_mapping:
                raise ValueError(
                    f"Atom {atom_id} has a type but no coordinates "
                    f"in {file_path}."
                )

            graph.add_node(
                atom_id,
                atom_type=atom_type,
                coords=coord_mapping[atom_id],
            )

        for bond_id, bond_data in bond_mapping.items():
            bond_type, atom1_id, atom2_id = bond_data

            if atom1_id not in graph or atom2_id not in graph:
                raise ValueError(
                    f"Bond {bond_id} references an undefined atom "
                    f"in {file_path}."
                )

            graph.add_edge(
                atom1_id,
                atom2_id,
                bond_id=bond_id,
                bond_type=bond_type,
            )

        print(
            f"Graph created from {file_path} with "
            f"{graph.number_of_nodes()} nodes and "
            f"{graph.number_of_edges()} edges."
        )

        return graph

    def _read_sections(
        self,
        file_path: Path,
    ) -> dict[str, list[str]]:
        """
        Read relevant sections from a LAMMPS molecule-template file.
        """
        supported_sections = {
            "Types",
            "Coords",
            "Bonds",
            "Charges",
            "Molecules",
            "Angles",
            "Dihedrals",
            "Impropers",
            "Special Bond Counts",
            "Special Bonds",
        }

        sections: dict[str, list[str]] = {}
        current_section: str | None = None

        with file_path.open("r", encoding="utf-8") as file:
            for raw_line in file:
                # Remove comments but preserve the actual data.
                line = raw_line.split("#", maxsplit=1)[0].strip()

                if not line:
                    continue

                if line in supported_sections:
                    current_section = line
                    sections.setdefault(current_section, [])
                    continue

                if current_section is not None:
                    # Stop collecting when another unsupported header/count
                    # line is encountered only through recognized sections.
                    sections[current_section].append(line)

        return sections

    def _couple_graphs(
        self,
        pre_template_graph: nx.Graph,
        post_template_graph: nx.Graph,
    ) -> dict[str, nx.Graph]:
        """
        Store the pre- and post-template graphs together.
        """
        return {
            "pre_template_graph": pre_template_graph,
            "post_template_graph": post_template_graph,
        }

    def _compare_graphs(
        self,
        mol_file_paths: list[str | Path],
    ) -> dict[str, bool]:
        """
        Process pre-template files and report whether each reaction is
        a duplicate of a previously processed pre/post pair.

        Returns:
            Mapping from the pre-template file path to duplicate status.
        """
        results: dict[str, bool] = {}

        for file_path_value in mol_file_paths:
            file_path = Path(file_path_value)

            if "pre" not in file_path.name:
                continue

            post_file_name = file_path.name.replace("pre", "post", 1)
            post_file_path = file_path.with_name(post_file_name)

            if not post_file_path.exists():
                print(
                    f"Post-template file does not exist for "
                    f"{file_path}."
                )
                continue

            pre_template_graph = self.lammps_molecule_to_networkx(
                file_path
            )
            post_template_graph = self.lammps_molecule_to_networkx(
                post_file_path
            )

            is_duplicate = self.is_duplicate(
                pre_template_graph,
                post_template_graph,
            )

            results[str(file_path)] = is_duplicate

            if is_duplicate:
                print(
                    f"Duplicate reaction: {file_path.name} and "
                    f"{post_file_path.name}"
                )
            else:
                print(
                    f"Unique reaction: {file_path.name} and "
                    f"{post_file_path.name}"
                )

        return results


if __name__ == "__main__":
    deduplication_detector = DeduplicationDetector()

    folder_path = Path(
        "/mnt/c/Users/Janitha/Documents/AutoREACTER/examples/"
        "AutoREACTER_outputs/"
        "Epoxy_Test_Primary_Diamine_Diepoxy"
    )

    pre_template_files = sorted(
        file_path
        for file_path in folder_path.glob("*.molecule")
        if "pre" in file_path.name
    )

    results = deduplication_detector._compare_graphs(
        pre_template_files
    )

    print("\nDeduplication results:")

    for file_path, is_duplicate in results.items():
        status = "duplicate" if is_duplicate else "unique"
        print(f"{Path(file_path).name}: {status}")