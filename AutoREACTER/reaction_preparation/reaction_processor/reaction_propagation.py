from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Iterable, Set, Tuple, TYPE_CHECKING
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from pathlib import Path
import pandas as pd

if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata

from AutoREACTER.detectors.functional_groups_detector import (
    FunctionalGroupsDetector,
    MonomerRole,
)
from AutoREACTER.detectors.reaction_detector import (
    ReactionDetector,
    ReactionInstance,
)
from AutoREACTER.reaction_preparation.reaction_processor.utils import (
    add_dict_as_new_columns, add_column_safe, prepare_paths
)
from AutoREACTER.reaction_preparation.reaction_processor.walker import reaction_atom_walker
from AutoREACTER.reaction_preparation.reaction_processor.reaction_processing_support import (
    assign_first_shell_and_initiators,
    assign_atom_map_numbers_and_set_isotopes,
    build_atom_index_mapping,
    build_reaction,
    build_reactants,
    build_reaction_tuple,
    clear_isotopes,
    detect_byproducts,
    detect_duplicates,
    reassign_atom_map_numbers_by_isotope,
    reveal_template_map_numbers,
    validate_mapping,
)

MAX_ITERATIONS = 6  # Maximum iterations for reaction pooling to prevent infinite loops


@dataclass(slots=True)
class TemplateIndexedMolecule:
    mol: Chem.Mol
    indexes: List[int]
    name: Optional[str] = None


@dataclass(slots=True)
class ReactantPool:
    """
    A pool of reactants for possible reactions for multiple reaction loops.
    """
    loop_no: int
    pool: List[MonomerRole]


@dataclass(slots=True, frozen=True)
class CombinedReactionSMILES:
    reactant_key: str
    product_key: str


class InfiniteReactionLoopError(Exception):
    """Raised when the reaction pooling process exceeds a reasonable number of iterations, indicating a potential infinite loop."""


class ReactionPropagation:
    def __init__(self, cache: Path):
        """
        Initialize ReactionPropagation with a cache directory.

        Args:
            cache: Path to the cache directory for storing intermediate CSV files.
        """
        self.cache = Path(cache)
        self.csv_cache = prepare_paths(self.cache, "csv_cache")
        self.fg_detector = FunctionalGroupsDetector()
        self.rxn_detector = ReactionDetector()


    #  PUBLIC                                                              

    def run_propagation_loop(
        self,
        reactions_metadata: list[ReactionMetadata],
        monomer_roles: list[MonomerRole],
    ) -> list[ReactionMetadata]:
        """
        Propagation loop that iteratively discovers new reactions from the
        products of previous reactions.

        Each iteration:
          1. Converts prior reaction products into TemplateIndexedMolecules.
          2. Detects functional groups restricted to the mapped product indices.
          3. Detects which new polymerization reactions are now possible.
          4. Filters reactions already seen (by canonical SMILES + reaction name).
          5. Runs prepare_reactions on the new batch and appends results.
          6. Stops when no new reactions are found or MAX_ITERATIONS is hit.

        Args:
            reactions_metadata: Seed metadata from the first reaction pass.
            monomer_roles: Original monomer roles (used for guard / future use).

        Returns:
            The fully accumulated list of ReactionMetadata (seed + all propagation rounds).

        Raises:
            ValueError: If monomer_roles is None.
            InfiniteReactionLoopError: If MAX_ITERATIONS is exceeded without convergence.
        """
        if monomer_roles is None:
            raise ValueError(
                "Monomer roles must be provided for the propagation loop."
            )
        # seen_instance_keys prevents re-running identical (mol-pair, rxn-name) combos
        seen_instance_keys: Set[Tuple] = self._build_seen_instance_keys(reactions_metadata)

        for iteration in range(MAX_ITERATIONS + 1):
            if iteration >= MAX_ITERATIONS:
                raise InfiniteReactionLoopError(
                    f"Exceeded maximum iterations ({MAX_ITERATIONS}) in reaction "
                    "propagation loop. Potential infinite loop detected. "
                    "Try increasing MAX_ITERATIONS if this is a false positive. "
                    "If the issue persists, please raise an issue at: "
                    "https://github.com/NanoCIPHER-Lab/AutoREACTER/issues"
                )

            # --- Step 1: build indexed product molecules from current metadata ---
            template_indexed_molecules = self._prepare_for_second_loop(reactions_metadata)

            if not template_indexed_molecules:
                # Nothing left to propagate from
                return reactions_metadata

            
            # --- Step 2: index-restricted functional-group detection ---
            detected_roles = self.fg_detector.index_based_functional_groups_detector(
                template_indexed_molecules
            )

            if not detected_roles:
                # No reactive sites remaining — propagation converged
                return reactions_metadata
            import json
            from dataclasses import dataclass, asdict
            
            monomer_roles.extend(detected_roles)  # Add original monomer roles back in for reaction detection
            
            roles_dict_list = [asdict(role) for role in monomer_roles]
            json_data = json.dumps(roles_dict_list, indent=4)
            
            print(f"Iteration {iteration}: Detected {len(monomer_roles)} functional group roles for reaction detection.")
            print(json_data)  # Print the JSON data for debugging purposes
            # --- Step 3: detect candidate reactions from the new roles ---
            candidate_instances = self.rxn_detector.reaction_detector(
                monomer_roles
            )
            # to begin with the loop we will need to  self.fg_detector.index_based_functional_groups_detector( give bothe template indexes molecules 
            # and original monomer roles so in that script we can flag weather a molecule a original monomer or a product of a reaction so that the 
            # reaction_detector can detect reactions from both original monomers and new products. 
            # This is because the reaction_detector is currently designed to take monomer roles as 
            # input and produce reaction instances, but in the propagation loop we want to be able to detect 
            # reactions that can occur between newly formed products as well as between products and original monomers. 
            # Therefore, the functional group detection step needs to be aware of which molecules are original monomers 
            # and which are new products, so it can assign roles accordingly for the reaction detection step. Additionally, 
            # the reaction detection logic may need to be enhanced to consider reactions between products, which may not have 
            # been necessary in the initial pass where only original monomers were considered. This is a critical next step to 
            # enable the propagation loop to discover new reactions iteratively from the evolving set of molecules.
            # we need reaction_detector.index_based_reaction_detector to return the reaction with indexes weather monomer
            # or a product to begin with also we can not reconstruct the molecules which may lead to errors reconstructing 
            # the same molecule multiple times 
            raise NotImplementedError("Reaction detection from detected roles is not yet implemented. "
                "This is a critical next step to enable the propagation loop. "
                "Please implement the logic to convert detected functional group roles "
                "into ReactionInstance objects that can be processed in the next steps."
            )
            # --- Step 4: filter out already-seen (mol-pair, reaction-name) combos ---
            new_instances = self._filter_seen_instances(
                candidate_instances, seen_instance_keys
            )

            if not new_instances:
                # All candidates were already processed — propagation converged
                return reactions_metadata

            # --- Step 5: prepare new reactions and append to master list ---
            new_metadata = self.prepare_reactions(new_instances)
            reactions_metadata.extend(new_metadata)

            # Register newly processed instances so they are not repeated
            for instance in new_instances:
                key = self._build_instance_key(instance)
                seen_instance_keys.add(key)

        # Unreachable — loop exits via return or raises above
        return reactions_metadata

    def prepare_reactions(
        self, reaction_instances: list[ReactionInstance]
    ) -> list[ReactionMetadata]:
        """
        Main pipeline: processes reaction instances, detects duplicates, and
        enriches metadata with template mappings and edge atoms.

        Args:
            reaction_instances: List of detected reaction instances to process.

        Returns:
            List of processed ReactionMetadata objects with template mappings
            and edge atoms.
        """
        reactions_metadata = self._process_reaction_instances(reaction_instances)
        reaction_metadata = detect_duplicates(reactions_metadata)

        for reaction in reaction_metadata:
            if not reaction.activity_stats:
                continue  # Skip reactions marked as duplicates

            combined_reactant_molecule_object = reaction.reactant_combined_RDmol
            reaction_dataframe = reaction.reaction_dataframe
            csv_save_path = reaction.csv_path

            # Build full atom mapping dictionary from dataframe
            fully_mapped_dict = (
                reaction_dataframe.set_index("reactant_idx")["product_idx"].to_dict()
            )
            first_shell = reaction_dataframe["first_shell"].dropna().tolist()

            # Generate template mapping by walking reaction graph
            template_mapped_dict, edge_atoms = reaction_atom_walker(
                combined_reactant_molecule_object,
                first_shell,
                fully_mapped_dict,
            )

            # Add template mapping and edge atoms to dataframe
            reaction_dataframe = add_dict_as_new_columns(
                reaction_dataframe,
                template_mapped_dict,
                titles=["template_reactant_idx", "template_product_idx"],
            )
            reaction_dataframe = add_column_safe(
                reaction_dataframe, edge_atoms, "edge_atoms"
            )

            # Persist updated dataframe and update metadata fields
            reaction.reaction_dataframe = reaction_dataframe.copy()
            reaction_dataframe.to_csv(csv_save_path, index=False)
            reaction.edge_atoms = edge_atoms
            reaction.template_reactant_to_product_mapping = template_mapped_dict

        return reaction_metadata

    #  PRIVATE — reaction processing

    def _process_reaction_instances(
        self, detected_reactions: list[ReactionInstance]
    ) -> list[ReactionMetadata]:
        """
        Converts ReactionInstance objects into ReactionMetadata by building
        molecules and running reactions.

        Args:
            detected_reactions: List of detected reaction instances.

        Returns:
            List of ReactionMetadata objects with atom mappings.
        """
        reaction_metadata: list[ReactionMetadata] = []

        for reaction in detected_reactions:
            rxn_smarts = reaction.reaction_smarts
            reactant_smiles_1 = reaction.monomer_1.smiles
            same_reactants = reaction.same_reactants

            reactant_smiles_2 = (
                reactant_smiles_1 if same_reactants else reaction.monomer_2.smiles
            )
            delete_atoms = reaction.delete_atom

            rxn = build_reaction(rxn_smarts)
            mol_reactant_1, mol_reactant_2 = build_reactants(
                reactant_smiles_1, reactant_smiles_2
            )
            reaction_tuple = build_reaction_tuple(
                same_reactants, mol_reactant_1, mol_reactant_2
            )

            reaction_metadata = self._process_reaction_products(
                rxn,
                self.csv_cache,
                reaction_tuple,
                delete_atoms,
                reaction_metadata,
            )

        return reaction_metadata

    def _process_reaction_products(
        self,
        rxn: ChemicalReaction,
        csv_cache: Path,
        reaction_tuple: list,
        delete_atoms: bool = True,
        reaction_metadata: Optional[list[ReactionMetadata]] = None,
    ) -> list[ReactionMetadata]:
        """
        Runs reactions on reactant pairs and builds metadata for each product set.

        For full (original) molecules the standard RunReactants path is used.
        For index-restricted molecules the forced-reaction path is used, keeping
        only product sets where the available_sites actually participated in the
        SMARTS match.

        Args:
            rxn: RDKit ChemicalReaction object.
            csv_cache: Path to cache directory for saving CSVs.
            reaction_tuple: List of reactant pairs to process.
                Each element is either:
                  - [mol1, mol2]               — standard full-molecule reaction
                  - [mol1, mol2, available_sites] — forced / index-restricted reaction
            delete_atoms: Whether to detect and track byproducts.
            reaction_metadata: Accumulator list for metadata objects.

        Returns:
            Updated list of ReactionMetadata objects.
        """
        if reaction_metadata is None:
            reaction_metadata = []

        for pair in reaction_tuple:
            r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])

            # Determine whether this is a forced (index-restricted) reaction
            available_sites: Optional[List[int]] = pair[2] if len(pair) > 2 else None

            assign_atom_map_numbers_and_set_isotopes(r1, r2)

            if available_sites is not None:
                # --- Forced / index-restricted path ---
                raw_product_sets = self._forced_reaction(
                    rxn, (r1, r2), available_sites
                )
                # raw_product_sets: list of (reactant_tuple, product_set)
                product_sets = [ps for _, ps in raw_product_sets]
            else:
                # --- Standard full-molecule path ---
                products = rxn.RunReactants((r1, r2))
                product_sets = list(products) if products else []

            if not product_sets:
                continue

            for product_set in product_sets:
                df = pd.DataFrame(columns=["reactant_idx", "product_idx"])

                reactant_combined = Chem.CombineMols(r1, r2)
                if len(product_set) == 1:
                    product_combined = product_set[0]
                else:
                    product_combined = Chem.CombineMols(product_set[0], product_set[1])

                reassign_atom_map_numbers_by_isotope(product_combined)

                mapping_dict, df = build_atom_index_mapping(
                    reactant_combined, product_combined
                )
                reverse_mapping = {v: k for k, v in mapping_dict.items()}

                reveal_template_map_numbers(product_combined)
                validate_mapping(df, reactant_combined, product_combined)

                first_shell, initiator_idxs = assign_first_shell_and_initiators(
                    reactant_combined, product_combined, reverse_mapping
                )

                byproduct_reactant_idxs = detect_byproducts(
                    product_combined, reverse_mapping, delete_atoms
                )

                df_combined = pd.concat(
                    [
                        df,
                        pd.Series(first_shell, name="first_shell"),
                        pd.Series(initiator_idxs, name="initiators"),
                        pd.Series(byproduct_reactant_idxs, name="byproduct_idx"),
                    ],
                    axis=1,
                ).astype(pd.Int64Dtype())

                total_products = len(reaction_metadata) + 1

                clear_isotopes(reactant_combined, product_combined)

                df_combined.to_csv(
                    csv_cache / f"reaction_{total_products}.csv", index=False
                )

                from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import ReactionMetadata

                reaction_metadata.append(
                    ReactionMetadata(
                        reaction_id=total_products,
                        reactant_combined_RDmol=reactant_combined,
                        product_combined_RDmol=product_combined,
                        reactant_to_product_mapping=mapping_dict,
                        product_to_reactant_mapping=reverse_mapping,
                        first_shell=first_shell,
                        initiators=initiator_idxs,
                        byproduct_indices=byproduct_reactant_idxs,
                        csv_path=csv_cache / f"reaction_{total_products}.csv",
                        reaction_dataframe=df_combined,
                        delete_atom=delete_atoms,
                        delete_atom_idx=(
                            byproduct_reactant_idxs[0]
                            if byproduct_reactant_idxs
                            else None
                        ),
                        activity_stats=True,
                    )
                )

        return reaction_metadata

    def _forced_reaction(
        self,
        rxn: ChemicalReaction,
        reactant_tuple: tuple[Chem.Mol, Chem.Mol],
        available_sites: List[int],
    ) -> list[tuple[tuple[Chem.Mol, Chem.Mol], list[Chem.Mol]]]:
        """
        Run a reaction and keep only product sets where the reacting atoms
        overlap with the caller-supplied available_sites indices.

        The SMARTS match atom indices returned by RunReactants tell us exactly
        which atoms in the reactant molecules were mapped by the reaction
        template.  We discard any product set whose matched atoms have NO
        overlap with available_sites, ensuring that only reactions occurring at
        the designated propagation sites are accepted.

        Args:
            rxn: The RDKit ChemicalReaction object.
            reactant_tuple: (mol1, mol2) to react.
            available_sites: Atom indices (in the combined / individual reactant
                space) that are permitted reaction sites for this propagation
                step.

        Returns:
            List of (reactant_tuple, product_set) for every accepted product set.
        """
        r1, r2 = reactant_tuple

        # RunReactants returns product sets; GetReactingAtoms returns, per
        # reactant, the atom indices that participated in the reaction.
        raw_products = rxn.RunReactants((r1, r2))

        if not raw_products:
            return []

        available_set = set(available_sites)
        accepted: list[tuple[tuple[Chem.Mol, Chem.Mol], list[Chem.Mol]]] = []

        for product_set in raw_products:
            # GetReactingAtoms(mapped=False) gives per-reactant lists of atom
            # indices that were actually transformed by the reaction SMARTS.
            reacting_per_reactant: tuple[tuple[int, ...], ...] = (
                rxn.GetReactingAtoms(mappedAtomsOnly=False)
            )

            # Flatten all reacting atom indices across both reactants.
            # r2 atoms are offset by r1.GetNumAtoms() in the combined molecule,
            # but here we compare against the per-molecule indices directly, so
            # we check r1 sites and r2 sites separately against available_set.
            r1_reacting = set(reacting_per_reactant[0]) if len(reacting_per_reactant) > 0 else set()
            r2_reacting = set(reacting_per_reactant[1]) if len(reacting_per_reactant) > 1 else set()

            # Offset r2 indices into the combined-molecule index space
            r2_offset = r1.GetNumAtoms()
            r2_reacting_combined = {idx + r2_offset for idx in r2_reacting}

            all_reacting_combined = r1_reacting | r2_reacting_combined

            # Accept the product set only if at least one reacting atom is an
            # available (propagation) site
            if all_reacting_combined & available_set:
                accepted.append((reactant_tuple, list(product_set)))

        return accepted

    #  PRIVATE — second-loop preparation                                  #

    def _prepare_for_second_loop(
        self,
        reaction_metadata_list: list[ReactionMetadata],
    ) -> list[TemplateIndexedMolecule]:
        """
        Converts active reaction products into TemplateIndexedMolecule objects
        so that the index-restricted functional-group detector can operate on
        the correct atom positions in the next propagation iteration.

        The indexes stored on each TemplateIndexedMolecule correspond to the
        *product* atom indices that descended from the template reaction sites
        (i.e. the keys of template_reactant_to_product_mapping inverted to
        product space).  These are the only atoms that should be examined for
        new reactive functionality in the propagation loop.

        Args:
            reaction_metadata_list: Current full list of reaction metadata.

        Returns:
            List of TemplateIndexedMolecule ready for index-based FG detection.
        """
        template_indexed_molecules: list[TemplateIndexedMolecule] = []

        for reaction in reaction_metadata_list:
            if not reaction.activity_stats:
                continue  # Skip duplicates / inactive reactions

            if (
                reaction.product_combined_RDmol is not None
                and reaction.template_reactant_to_product_mapping is not None
            ):
                # Invert mapping: template reactant idx -> template product idx
                # Keys of the inverted map are the relevant product-space indices
                template_product_to_reactant_to_product_mapping = {
                    v: k
                    for k, v in reaction.template_reactant_to_product_mapping.items()
                }

                template_indexed_molecules.append(
                    TemplateIndexedMolecule(
                        mol=reaction.product_combined_RDmol,
                        indexes=list(
                            template_product_to_reactant_to_product_mapping.keys()
                        ),
                        name=f"reaction_{reaction.reaction_id}_product",
                    )
                )

        return template_indexed_molecules

    #  PRIVATE — seen-instance deduplication                              #

    def _build_instance_key(self, instance: ReactionInstance) -> Tuple:
        """
        Build a canonical deduplication key for a ReactionInstance.

        The key is (reaction_name, canonical_smiles_1, canonical_smiles_2)
        where the two SMILES are sorted so that A+B and B+A produce the same key.

        Args:
            instance: The reaction instance to key.

        Returns:
            A hashable tuple uniquely representing this (mol-pair, reaction) combo.
        """
        smiles_1 = self._canonical_smiles(instance.monomer_1.smiles)
        smiles_2 = self._canonical_smiles(
            instance.monomer_2.smiles if instance.monomer_2 else instance.monomer_1.smiles
        )

        # Sort so order of reactants doesn't produce a duplicate key
        ordered = tuple(sorted([smiles_1, smiles_2]))
        return (instance.reaction_name, ordered)

    def _build_seen_instance_keys(
        self,
        reactions_metadata: Iterable[ReactionMetadata],
    ) -> Set[Tuple]:
        """
        Pre-populate the seen-instance set from the seed metadata produced by
        the first reaction pass.

        We reconstruct keys from the combined-reactant SMILES because
        ReactionMetadata does not store the originating ReactionInstance.

        Args:
            reactions_metadata: Seed reaction metadata.

        Returns:
            Set of seen keys (reaction_name is unavailable in metadata, so the
            key here is the frozenset of canonical reactant fragment SMILES,
            which is sufficient to prevent re-running the same mol-pair).
        """
        seen: Set[Tuple] = set()

        for reaction in reactions_metadata:
            if not reaction.activity_stats:
                continue

            key = self._build_metadata_key(reaction)
            seen.add(key)

        return seen

    def _build_metadata_key(self, reaction: ReactionMetadata) -> Tuple:
        """
        Build a deduplication key directly from a ReactionMetadata object.

        Since ReactionMetadata does not store the reaction name, the key uses
        the canonical SMILES of both the reactant and product combined molecules,
        which uniquely identifies the transformation.

        Args:
            reaction: ReactionMetadata to key.

        Returns:
            Hashable tuple (reactant_canonical_smiles, product_canonical_smiles).
        """
        r_smiles = self._canonical_smiles_from_mol(reaction.reactant_combined_RDmol)
        p_smiles = self._canonical_smiles_from_mol(reaction.product_combined_RDmol)
        return (r_smiles, p_smiles)

    def _filter_seen_instances(
        self,
        candidate_instances: list[ReactionInstance],
        seen_keys: Set[Tuple],
    ) -> list[ReactionInstance]:
        """
        Remove candidate ReactionInstances whose (mol-pair, reaction-name) key
        has already been processed.

        Args:
            candidate_instances: Newly detected reaction instances.
            seen_keys: Set of already-processed keys.

        Returns:
            Filtered list containing only genuinely new instances.
        """
        new_instances: list[ReactionInstance] = []

        for instance in candidate_instances:
            key = self._build_instance_key(instance)
            if key not in seen_keys:
                new_instances.append(instance)

        return new_instances

    #  PRIVATE — SMILES helpers

    def _canonical_smiles(self, smiles: str) -> str:
        """
        Convert a SMILES string to its canonical RDKit form.
        Falls back to the original string if parsing fails.

        Args:
            smiles: Input SMILES string.

        Returns:
            Canonical SMILES string.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    def _canonical_smiles_from_mol(self, mol: Chem.Mol) -> str:
        """
        Build a canonical SMILES string from an RDKit Mol object.

        Args:
            mol: RDKit molecule.

        Returns:
            Canonical SMILES string, or empty string if mol is None.
        """
        if mol is None:
            return ""
        return Chem.MolToSmiles(
            Chem.Mol(mol), canonical=True, isomericSmiles=True
        )

    def _build_combined_molecule_SMILES(self, mol: Chem.Mol) -> str:
        """
        Build a stable key for a combined molecule.
        """
        if mol is None:
            raise ValueError("Combined molecule cannot be None.")

        return Chem.MolToSmiles(
            Chem.Mol(mol),
            canonical=True,
            isomericSmiles=True,
        )

    def _build_reaction_state(
        self,
        reactant_combined_mol: Chem.Mol,
        product_combined_mol: Chem.Mol,
    ) -> CombinedReactionSMILES:
        """
        Build a comparable reactant/product state from combined molecules.
        """
        return CombinedReactionSMILES(
            reactant_key=self._build_combined_molecule_SMILES(reactant_combined_mol),
            product_key=self._build_combined_molecule_SMILES(product_combined_mol),
        )

    def _build_seen_reaction_combinations(
        self,
        reactions_metadata: Iterable[ReactionMetadata],
    ) -> set[CombinedReactionSMILES]:
        """
        Build the initial seen-state set before entering the propagation loop.
        """
        seen_combinations: set[CombinedReactionSMILES] = set()

        for reaction in reactions_metadata:
            if not reaction.activity_stats:
                continue

            state = self._build_reaction_state(
                reaction.reactant_combined_RDmol,
                reaction.product_combined_RDmol,
            )
            seen_combinations.add(state)

        return seen_combinations

    def _is_seen_reaction(
        self,
        reactant_combined_mol: Chem.Mol,
        product_combined_mol: Chem.Mol,
        seen_reactions: set[CombinedReactionSMILES],
    ) -> bool:
        """
        Check whether the reactant/product combined-molecule pair has already been seen.
        """
        reaction = self._build_reaction_state(
            reactant_combined_mol,
            product_combined_mol,
        )
        return reaction in seen_reactions

    def _register_reaction_state(
        self,
        reactant_combined_mol: Chem.Mol,
        product_combined_mol: Chem.Mol,
        seen_reactions: set[CombinedReactionSMILES],
    ) -> None:
        """
        Register a newly accepted reactant/product state.
        """
        reaction = self._build_reaction_state(
            reactant_combined_mol,
            product_combined_mol,
        )
        seen_reactions.add(reaction)

    def _compare_combined_molecules_placeholder_networkx(
        self,
        reactant_combined_mol: Chem.Mol,
        product_combined_mol: Chem.Mol,
    ) -> bool:
        """
        PLACEHOLDER:
        Compare two combined molecules using a graph-based approach (e.g., NetworkX).
        This is a more computationally intensive method that can be used as a fallback
        when SMILES-based comparison is inconclusive.
        """
        return False