"""
* Monomer Functionality Detection Module
--------------------------------------
Purpose
This module takes a set of monomer SMILES strings and determines which monomers
are chemically “eligible” to participate in polymerization reactions. It does
this by scanning each monomer for known reactive functional groups using RDKit
SMARTS substructure matching. The output is a structured dictionary describing
which functional group(s) were found in each monomer and how many reactive sites
of each type exist.
* What “functionality” means here
- A “functional group” is a chemical motif that can react (e.g., -OH, -COOH,
  epoxide rings, vinyl double bonds, etc.).
- A “reactive site count” is the number of times that motif appears in the
  monomer. For example, ethylene glycol has two -OH sites; adipic acid has two
  -COOH sites.
- A “functionality type” is a rule for deciding whether the monomer is allowed
  to proceed to the polymerization engine. This module is intentionally a
  FILTER + ANNOTATOR:
  (1) Filter: remove monomers that would cause chain termination or invalid
      step-growth logic.
  (2) Annotator: retain counts so later steps can handle multi-site reactions.
* Core approach
1) Convert SMILES → RDKit Mol
2) For each defined functional group class (stored in monomer_types):
   - Build RDKit pattern(s) from SMARTS
   - Find ALL matches using substructure matching
   - Count how many matches exist
   - Decide whether the monomer “qualifies” under the rules below
   - If qualified, record:
       • functionality_type
       • functional_group_name
       • SMARTS pattern(s)
       • site count(s)
* Functionality Elaboration Logic (rules enforced)
1) Vinyl and Mono (presence-based eligibility)
   - Goal: identify monomers that can undergo chain-growth / ring-opening style
     Reactions where a single reactive motif is sufficient to participate.
   - Pass condition:
       • at least one match to the SMARTS pattern (count >= 1)
   - Notes:
       • If a monomer has multiple vinyl sites (e.g., divinyl crosslinker),
         it still “passes” here because it is reactive. The multiplicity is
         tracked in the site count and is intended to be exploited later by
         the reaction engine (crosslinking / branching logic).
       • This stage does NOT decide how many reactions will happen; it only
         verifies that the reaction is possible and records the capacity.
* 2) Di_Identical (step-growth safety rule: prevent termination)
   - Goal: enforce correct step-growth behavior for monomers expected to link
     through identical end-groups (A-A type monomers such as diols, diacids,
     diamines, diisocyanates, diepoxides).
   - Pass condition:
       • at least two occurrences of the SAME functional group (count >= 2)
   - Fail behavior:
       • if only one match is found (count == 1), the monomer is treated as
         non-functional for this category (i.e., rejected for step-growth).
   - Why this matters:
       • A monomer with only one reactive end group behaves like a chain
         stopper in step-growth polymerization, which can silently break
         propagation and produce wrong molecular weight distributions.
       • This module intentionally blocks those “A1” species from entering
         The step-growth reaction engine under a “di_identical” assumption.
* 3) Di_Different (heterobifunctional eligibility)
   - Goal: identify monomers that carry two distinct reactive groups (A-B type)
     such as amino acids (NH2 + COOH), hydroxy acids (OH + COOH), etc.
   - Pass condition:
       • at least one match to SMARTS_1 (count1 >= 1)
       • AND at least one match to SMARTS_2 (count2 >= 1)
   - Notes:
       • The module records BOTH counts. This is important because a monomer
         may have uneven stoichiometry (e.g., OH=1 but COOH=2), which can enable
         branching or different pairing possibilities later.
       • This stage still does not enumerate pairings; it only guarantees both
         chemical capabilities exist and preserve the information required for
         Later enumeration.
* What this module guarantees downstream
- Every monomer returned by this detector has at least one valid reactive motif
  according to the category rules above.
- Step-growth “di_identical” monomers have at least two identical reactive sites,
  reducing the risk of chain termination artifacts.
- “di_different” monomers contain both required group types.
- Reactive site multiplicities are explicitly captured (counts), enabling later
  logic to generate many possible reactions for monomers with multiple sites.
* What this module does NOT attempt to do (by design)
    - It does not enumerate all possible site-to-site reaction pairings.
    - It does not resolve conflicts between overlapping SMARTS matches.
* Output concept (high-level)
For each monomer index, the detector returns:
- "smiles": the original monomer SMILES
- One or more detected functional group entries (indexed 1..N), each including:
    • functionality_type (vinyl / mono / di_identical / di_different)
    • functional_group_name (e.g., diol, di_carboxylic_acid, epoxide, vinyl)
    • functional_group_smarts_1 (+ smarts_2 if applicable)
    • functional_count_1 (+ count_2 if applicable)
* This output is the “functional inventory” that the polymerization engine uses
later to determine how many reactions are possible and how to build reaction
templates safely.
"""

from __future__ import annotations # Enable postponed evaluation of annotations for forward references.
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, rdchem
import logging
import json
from typing import Dict, List, Tuple, Optional
import pathlib, os
from PIL.Image import Image
from AutoREACTER.input_parser import MonomerEntry
from AutoREACTER.detectors.functional_groups_library import FunctionalGroupsLibrary
from AutoREACTER.reaction_preparation.reaction_processor.prepare_reactions import TemplateIndexedMolecule

logger = logging.getLogger(__name__)  # Module-level logger for future diagnostics.

# Dictionary defining various types of monomers would typically be loaded here, but is accessed via FunctionalGroupsLibrary.
# Includes functionality types (e.g., 'vinyl', 'mono', 'di_different', 'di_identical'), SMARTS patterns,
# and group names. Expandable based on literature (e.g., "J. Chem. Inf. Model. 2023, 63, 5539−5548").
# TODO: Address monomers with mixed groups like COCl and COOH.
# TODO: Add more functional groups as needed from literature/user requirements.


@dataclass(slots=True)
class FunctionalGroupInfo:
    """
    Immutable dataclass storing information about a detected functional group.

    Attributes:
        functionality_type (str): Type of functionality (e.g., 'vinyl', 'mono', 'di_identical', 'di_different').
        fg_name (str): Name of the functional group (e.g., 'acrylate').
        fg_smarts_1 (str): Primary SMARTS pattern for matching.
        fg_count_1 (int): Number of matches for fg_smarts_1.
        fg_smarts_2 (Optional[str]): Secondary SMARTS pattern (for 'di_different' types).
        fg_count_2 (Optional[int]): Number of matches for fg_smarts_2.
    """
    functionality_type: str
    fg_name: str
    fg_smarts_1: str
    fg_count_1: int
    fg_smarts_2: Optional[str] = None
    fg_count_2: Optional[int] = None


@dataclass(slots=True, frozen=True)
class MonomerRole:
    """
    Immutable dataclass representing a monomer with its detected functional groups.

    Attributes:
        smiles (str): SMILES string of the monomer.
        name (str): Human-readable name of the monomer.
        functionalities (Tuple[FunctionalGroupInfo, ...]): Tuple of all detected functional groups.
    """
    smiles: str
    name: str
    indexes: Tuple[int, ...]
    functionalities: Tuple[FunctionalGroupInfo, ...]  # Tuple of detected functionalities for the monomer

class ReactantPool:
    """
    A pool of reactants for possible reactions for multiple reaction loops.
    """
    loop_no : int
    pool : List[MonomerRole] 

@dataclass(slots=True)
class FunctionalGroupVisualization:
    """
    Dataclass for storing visualization information about detected functional groups.

    Attributes:
        monomer: rdchem.Mol - RDKit molecule object for the monomer.
        name: str - Name of the monomer.
        indexes_to_highlight: Tuple[int, ...] - Atom indices to highlight in visualization.
    """
    monomer: rdchem.Mol
    name: str
    indexes_to_highlight: Tuple[int, ...]

class FunctionalGroupsDetector:
    """
    Detector class for identifying functional groups in monomers using SMARTS-based substructure matching.

    Uses a predefined FunctionalGroupsLibrary to match monomers against known patterns and classify them
    into roles like mono-functional, di-functional (identical or different), or vinyl.

    Attributes:
        functional_group_library (FunctionalGroupsLibrary): Library of predefined monomer types and SMARTS.
        monomer_types (dict): Dictionary of monomer types from the library.
        selected_monomers (dict): Internal storage for selected monomers (currently unused in return).
    """

    def __init__(self):
        """Initialize the detector with the functional groups library."""
        self.functional_group_library = FunctionalGroupsLibrary()
        self.monomer_types = self.functional_group_library.monomer_types

    def _detect_monomer_functionality(
        self, mol: Chem.Mol, functionality_type: str, smarts_1: str, smarts_2: Optional[str] = None
        ) -> Tuple[int, Optional[int], Optional[int], Optional[Tuple[Tuple[int, ...], ...]]]:
        """
        Detect the functionality type of a monomer based on its SMILES string using RDKit substructure matching
        with SMARTS patterns.

        Performs substructure matching and counts occurrences to determine if the monomer qualifies as
        mono-functional (1 match), di-functional identical (>=2 matches of same pattern), di-functional different
        (>=1 match each of two patterns), or no match.

        Args:
            mol (Chem.Mol): RDKit molecule object.
            functionality_type (str): Expected type (e.g., 'mono', 'di_identical', 'di_different', 'vinyl').
            smarts_1 (str): Primary SMARTS pattern for matching.
            smarts_2 (str, optional): Secondary SMARTS pattern for 'di_different' types. Defaults to None.

        Returns:
            Tuple[int, Optional[int], Optional[int], Optional[Tuple[Tuple[int, ...], ...]]]:
                - functionality_count: 0 (no match), 1 (single match), or 2 (dual matches).
                - count_1: Number of matches for smarts_1 (or 0 if invalid).
                - count_2: Number of matches for smarts_2 (None if not applicable).
                - functional_matches: Tuple of tuples containing atom indices of matched functional groups.

        Notes:
            - Invalid SMILES or SMARTS patterns return (0, None, None).
            - For 'di_identical': Requires >=2 matches of smarts_1.
            - Duplicate logic noted in original; this implementation handles all cases.
            - No validation on SMARTS complexity or performance for large molecules.
        """
        # Convert SMILES to RDKit molecule object for substructure matching.
        
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return 0, None, None, None

        # Create pattern from primary SMARTS and validate.
        patt1 = Chem.MolFromSmarts(smarts_1)
        if patt1 is None:
            logger.warning(f"Invalid primary SMARTS: {smarts_1}")
            return 0, None, None, None

        if smarts_2:
            # Handle 'di_different' types: Match two distinct patterns.
            patt2 = Chem.MolFromSmarts(smarts_2)
            if patt2 is None:
                logger.warning(f"Invalid secondary SMARTS: {smarts_2}")
                return 0, None, None, None

            # Count matches once (avoid duplicate calls)
            matches1 = mol.GetSubstructMatches(patt1)
            matches2 = mol.GetSubstructMatches(patt2)

            functional_count_1 = len(matches1)
            functional_count_2 = len(matches2)

            # Combine safely (always tuple of tuples)
            functional_matches = matches1 + matches2

            # Debug
            # print(
            #     f"SMILES: {smiles}, "
            #     f"Pattern 1 Matches: {functional_count_1}, "
            #     f"Pattern 2 Matches: {functional_count_2}"
            # )

            if functional_count_1 >= 1 and functional_count_2 >= 1:
                return 2, functional_count_1, functional_count_2, functional_matches

            return 0, functional_count_1, functional_count_2, functional_matches
        else:
            # Single pattern matching for 'vinyl', 'mono', or 'di_identical'.
            functional_matches = mol.GetSubstructMatches(patt1)
            functional_count_1 = len(functional_matches)

            if functionality_type == "di_identical":
                # Require at least two identical functional groups.
                if functional_count_1 >= 2:
                    return 2, functional_count_1, None, functional_matches
                return 0, functional_count_1, None, functional_matches
            
            else:
                # For 'vinyl' or 'mono': Single match suffices.
                if functional_count_1 >= 1:
                    return 1, functional_count_1, None, functional_matches  
                
        return 0, 0, None, ()  # No matches found
    
    def _index_based_detect_monomer_functionality(self, mol: Chem.Mol, indexes: list[int], functionality_type: str, smarts_1: str, smarts_2: Optional[str] = None
        ) -> Tuple[int, Optional[int], Optional[int], Optional[Tuple[Tuple[int, ...], ...]]]:
        """
        Detect functional groups in a post-reaction template by restricting the search to specific atom indices.

        This function is intended for a second-pass functionality check after reaction-template generation.
        It focuses matching on mapped atom indices from a pre-reacted molecule so that functional groups are
        identified only at the relevant reacted sites. This is especially useful when the generated
        post-reaction template contains multiple functional groups, and only the secondary or site-specific
        functionalities need to be detected. It is also useful when atom types change during the reaction
        and the functional groups must be re-identified in the post-reaction template.

        Args:
            mol (Chem.Mol):
                RDKit molecule object representing the post-reaction template molecule.
            indexes (list[int]):
                List of atom indices used to restrict matching to specific substructures or mapped sites.
            functionality_type (str):
                Expected functionality type. Examples include 'mono', 'di_identical',
                'di_different', and 'vinyl'.
            smarts_1 (str):
                Primary SMARTS pattern used for functional group matching.
            smarts_2 (str, optional):
                Secondary SMARTS pattern used only when `functionality_type` is
                'di_different'. Defaults to None.

        Returns:
            Tuple[int, Optional[int], Optional[int], Optional[Tuple[Tuple[int, ...], ...]]]:
                A tuple containing:
                    - functionality_count:
                        0 if no valid match is found,
                        1 if a single functional group is found,
                        2 if two functional groups are found.
                    - count_1:
                        Number of matches for `smarts_1`, or 0 if `smarts_1` is invalid.
                    - count_2:
                        Number of matches for `smarts_2`, or None if not applicable.
                    - functional_matches:
                        Tuple of tuples containing the atom indices of the matched functional groups.

        Notes:
            - If the molecule or SMARTS patterns are invalid, the function returns
            (0, None, None, None).
            - For 'di_identical', at least two matches of `smarts_1` are required.
            - For 'di_different', both SMARTS patterns may be required depending on the logic.
            - This function is designed for targeted substructure detection and does not validate
            SMARTS complexity or optimize performance for large molecules.
        """
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return 0, None, None, None

        # Create pattern from primary SMARTS and validate.
        patt1 = Chem.MolFromSmarts(smarts_1)
        if patt1 is None:
            logger.warning(f"Invalid primary SMARTS: {smarts_1}")
            return 0, None, None, None

        if smarts_2:
            patt2 = Chem.MolFromSmarts(smarts_2)
            if patt2 is None:
                logger.warning(f"Invalid secondary SMARTS: {smarts_2}")
                return 0, None, None, None
            matches2 = mol.GetSubstructMatches(patt2)
        else:
            matches2 = ()
            # Count matches once (avoid duplicate calls)
        matches1 = mol.GetSubstructMatches(patt1)
            

        functional_count_1 = len(matches1)
        functional_count_2 = len(matches2)

        # Combine safely (always tuple of tuples)
        functional_matches = matches1 + matches2

        # Debug
        # print(
        #     f"SMILES: {smiles}, "
        #     f"Pattern 1 Matches: {functional_count_1}, "
        #     f"Pattern 2 Matches: {functional_count_2}"
        # )

        if functional_count_1 >= 1 or functional_count_2 >= 1:
            matches1_valid = self._find_matching_tuple(indexes, matches1)
            matches2_valid = self._find_matching_tuple(indexes, matches2)
            valid_matches = tuple(
                m for m in (matches1_valid, matches2_valid) if m is not None
            )
            if matches1_valid or matches2_valid:
                return 2, functional_count_1, functional_count_2, valid_matches
                
        return 0, 0, None, ()  # No matches found

    def functional_groups_detector(self, monomers: list[MonomerEntry]) -> Tuple[list[MonomerRole], list[FunctionalGroupVisualization]]:
        """
        Detect functional groups across a list of monomers and categorize them into roles.

        Iterates over predefined monomer_types, matches each against the monomer's SMILES,
        and collects valid detections. Prints matches for debugging/user feedback.

        Args:
            monomers (list[MonomerEntry]): List of MonomerEntry objects (each with 'smiles' and 'name').

        Returns:
            list[MonomerRole]: List of monomers with detected functionalities.
                Each MonomerRole contains smiles, name, and tuple of FunctionalGroupInfo.
            list[FunctionalGroupVisualization]: List of visualization data for detected functional groups.

        Notes:
            - Matches criteria: 'vinyl'/'mono' (>=1 primary), 'di_identical' (>=2 primary),
              'di_different' (>=1 each pattern).
        """

        # NOTE: Original code notes "from here code is still messed up; need to fix the loops".
        # Current implementation iterates monomers outer, types inner; consider inverting for multi-match optimization.

        # Collect roles for all monomers with detections.
        monomer_roles = []
        monomer_roles_visualization = []  # For potential future visualization output.
        

        for monomer in monomers:
            smiles = monomer.smiles
            detected_functionalities = []
            all_matches = []
            # Check against each predefined functional group type.
            for functional_group in self.monomer_types.values():
                ftype = functional_group["functionality_type"]
                smarts_1 = functional_group["smarts_1"]
                smarts_2 = functional_group.get("smarts_2")  # May be None for non-di_different types.
                functionality_count, count_1, count_2, functional_matches = self._detect_monomer_functionality(
                    monomer.rdkit_mol, ftype, smarts_1, smarts_2
                )

                # Determine if this functionality matches criteria.
                if functionality_count > 0:
                    # Add matches to the overall list for visualization.
                    all_matches.extend(functional_matches)
                    # Log match for user feedback.
                    print(f"{smiles} has functionality: {functional_group['group_name']}")
                    if functional_group.get("comments"):
                        print(f"Note: {smiles} - {functional_group['comments']}")
                    
                    # Store detection details.
                    detected_functionalities.append(
                        FunctionalGroupInfo(
                            functionality_type=ftype,
                            fg_name=functional_group['group_name'],
                            fg_smarts_1=smarts_1,
                            fg_count_1=count_1,
                            fg_smarts_2=smarts_2,
                            fg_count_2=count_2
                        )
                    )
                    # Debug print for match details.
                    # print(
                    #     f"Monomer {monomer.name} (SMILES: {smiles}) matches {ftype} "
                    #     f"with {functional_group['group_name']} (Count 1: {count_1}, Count 2: {count_2})"
                    # )

            # Add to roles if any functionalities detected.
            if detected_functionalities:
                monomer_roles.append(
                    MonomerRole(
                        smiles=smiles,
                        name=monomer.name,
                        indexes=(),  # Placeholder: no index-based detection in the first loop. 
                        functionalities=tuple(detected_functionalities)
                    )
                )
                monomer_roles_visualization.append(
                    FunctionalGroupVisualization(
                        monomer=monomer.rdkit_mol,
                        name=monomer.name,
                        indexes_to_highlight=tuple(all_matches)
                    ))
        return monomer_roles, monomer_roles_visualization
    
    def index_based_functional_groups_detector(
        self,
        dimers: list[TemplateIndexedMolecule]
    ) -> Tuple[list[MonomerRole], list[FunctionalGroupVisualization]]:
        """
        Detect functional groups in template/product molecules by restricting
        matching to user-provided atom indices.
        copy of functional_groups_detector but with index-based detection for post-reaction templates.

        Args:
            dimers (list[TemplateIndexedMolecule]):
                List of objects containing:
                    - mol: RDKit molecule
                    - indexes: atom indices to focus detection on

        Returns:
            Tuple[list[MonomerRole], list[FunctionalGroupVisualization]]:
                - Detected functional roles
                - Visualization metadata with matched atom indices
        """
        monomer_roles = []

        for monomer in dimers:
            mol = monomer.mol
            indexes = monomer.indexes

            detected_functionalities = []
            all_matches = []

            for functional_group in self.monomer_types.values():
                ftype = functional_group["functionality_type"]
                smarts_1 = functional_group["smarts_1"]
                smarts_2 = functional_group.get("smarts_2")

                functionality_count, count_1, count_2, functional_matches = (
                    self._index_based_detect_monomer_functionality(
                        mol=mol,
                        indexes=indexes,
                        functionality_type=ftype,
                        smarts_1=smarts_1,
                        smarts_2=smarts_2,
                    )
                )

                if functionality_count > 0:
                    all_matches.extend(functional_matches)

                    detected_functionalities.append(
                        FunctionalGroupInfo(
                            functionality_type=ftype,
                            fg_name=functional_group["group_name"],
                            fg_smarts_1=smarts_1,
                            fg_count_1=count_1,
                            fg_smarts_2=smarts_2,
                            fg_count_2=count_2,
                        )
                    )

            if detected_functionalities:
                monomer_roles.append(
                    MonomerRole(
                        smiles="",   # placeholder: user does not want MolToSmiles()
                        name=monomer.name,     # placeholder: user does not want to use name even though it's available in the object; can be adjusted as needed
                        indexes=tuple(indexes),
                        functionalities=tuple(detected_functionalities),
                    )
                )

        return monomer_roles
    
    def _find_matching_tuple(self, indexes, tuple_list):
        """
        Find the first tuple in tuple_list that contains at least one index from indexes.

        Args:
            indexes (list[int]): Atom indices to match against.
            tuple_list (list[tuple[int, ...]]): Candidate tuples of atom indices.

        Returns:
            tuple[int, ...] | None: First matching tuple, or None if no match.
        """
        index_set = set(indexes)
        for tup in tuple_list:
            if index_set.intersection(tup):
                return tup
        return None
    


    def functional_group_highlighted_molecules_image_grid(self, monomer_roles_visualization: list[FunctionalGroupVisualization]) -> Image:
            """Convert monomer roles with detected functionalities into visualizations.
            Args:
                monomer_roles_visualization (list[FunctionalGroupVisualization]): List of monomer roles with detected functionalities for visualization.
            Returns:
                PIL.Image.Image: A single grid image of the monomers with highlighted functional groups.
            """
            molecules = []
            indextohighlight = []
            names = []
            for viz in monomer_roles_visualization:  # <-- use actual objects, not dict version
                mol = viz.monomer
                name = viz.name
                matches = viz.indexes_to_highlight or ()

                # Flatten tuple-of-tuples into unique atom indices
                highlight_atoms = sorted(
                    {atom for match in matches for atom in match}
                )

                molecules.append(mol)
                indextohighlight.append(highlight_atoms)
                names.append(name)

            img = Draw.MolsToGridImage(
                molecules,
                molsPerRow=3,
                legends=names,
                subImgSize=(500, 500),
                highlightAtomLists=indextohighlight
            )

            return img

if __name__ == "__main__":

    monomer_dictionary = {
        1:  "C=CCO",                 # vinyl monomer # veryfied
        2:  "CCCCCCC=C",           
        3:  " CC=C1CC2CC1C=C2",           # cyclic olefin monomer (norbornene) # veryfied
        4: "C1=CC2CCC1C2",           # cyclic olefin monomer (norbornadiene) # veryfied
        5:  "C1CCC(=O)OCC1",           # lactone monomer (ε-caprolactone) # veryfied
        6:  "OCCC(O)=O",           # hydroxy carboxylic acid (lactic acid) # veryfied
        7:  "C1=CC=C(C(=C1)C(=O)O)O", # hydroxy carboxylic acid (p-hydroxybenzoic acid) # veryfied
        8:  "O=C(O)CCCCC(=O)O",      # di-carboxylic acid (adipic acid) # veryfied
        9:  "O=C(Cl)c1ccc(cc1)C(=O)Cl", # di-carboxylic acid halide (terephthaloyl chloride) # veryfied
        10: "OCCO",                  # diol (ethylene glycol) # veryfied
        11: "O=C1OC(=O)C=C1",        # cyclic anhydride (maleic anhydride) # veryfied
        12: "C1CO1",                 # epoxide monomer (ethylene oxide) # veryfied
        13: "O=C1CCCCN1",            # lactam monomer (caprolactam) # veryfied
        14: "NCC(=O)O",              # amino acid (glycine) # veryfied
        15: "NCCCC(=O)O",            # amino acid (aminocaproic acid) # veryfied
        16: "NCCN",                  # di-amine (ethylenediamine)
        17: "NCCCCCCN",              # primary di-amine (hexamethylenediamine)
        18: "O=C1OC(=O)c2cc(C3OC(=O)C=CC3=O)ccc2C1=O",  # di-cyclic anhydride (pyromellitic dianhydride)
        19: "O=C=N-C-N=C=O",                # di-isocyanate (simplified; represents TDI core)
        20: "C1=CC(=CC=C1C(C2=CC=CC=C2)(O)CC3CO3)CC4CO4" # di-epoxide (bisphenol A diglycidyl ether)
    }
    monomer_dictionary = {1: "O=C(Cl)c1cc(C(=O)Cl)cc(C(=O)Cl)c1",
                          2: "Nc1cccc(N)c1"}

    monomers = []

    for i, smiles in monomer_dictionary.items():
        smiles = smiles.strip()
        mol = Chem.MolFromSmiles(smiles)

        monomers.append(
            MonomerEntry(
                id=i,
                data_id=f"TEST_{i}",
                name=f"M{i}",
                smiles=smiles,
                count=None,
                ratio=1.0,
                atom_count=mol.GetNumAtoms() if mol else 0,
                molar_mass=Descriptors.MolWt(mol) if mol else 0.0,
                rdkit_mol=mol
            )
        )

    detector = FunctionalGroupsDetector()
    results = detector.functional_groups_detector(monomers)

    from dataclasses import asdict
    
    results_dict = [asdict(role) for role in results[0]]
    visualization_dict = [asdict(viz) for viz in results[1]]

    print(json.dumps(results_dict, indent=4))
    print(visualization_dict)

    from rdkit.Chem import Draw
    autoreacter_dir = pathlib.Path(__file__).parent.parent.parent.resolve()
    save_dir = autoreacter_dir / "cache" / "_cache_test" / "functional_groups_visualization"
    os.makedirs(save_dir, exist_ok=True)
    for viz in results[1]:  # <-- use actual objects, not dict version
        mol = viz.monomer
        name = viz.name
        matches = viz.indexes_to_highlight or ()

        # Flatten tuple-of-tuples into unique atom indices
        highlight_atoms = sorted(
            {atom for match in matches for atom in match}
        )

        file_path =  save_dir / f"{name}_functional_groups.png"

        Draw.MolToFile(
            mol,
            str(file_path),
            size=(500, 500),
            highlightAtoms=highlight_atoms
        )

        print(f"Saved: {file_path}")


        