import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AutoREACTER.detectors.reactions_library import registry
from AutoREACTER.detectors.reactions_library.polyamides import (
    REACTIONS,
)


AMINO_HOMO = (
    "Amino Acid Polycondensation (Polyamidation)"
)

AMINO_CO = (
    "Amino Acid and Amino Acid "
    "Polycondensation (Polyamidation)"
)

DIAMINE_ACID = (
    "Di-Amine and Di-Carboxylic Acid "
    "Polycondensation (Polyamidation)"
)

DIAMINE_HALIDE = (
    "Di-Amine and Di-Carboxylic Acid Halide "
    "Polycondensation (Polyamidation)"
)

CAPROLACTAM_HYDROLYSIS = (
    "Hydrolytic Initiation of Caprolactam"
)


EXPECTED_REACTIONS = {
    AMINO_HOMO,
    AMINO_CO,
    DIAMINE_ACID,
    DIAMINE_HALIDE,
    CAPROLACTAM_HYDROLYSIS,
}


def run_reaction(reaction_name, smiles_1, smiles_2):
    reaction_info = REACTIONS[reaction_name]

    rxn = rdChemReactions.ReactionFromSmarts(
        reaction_info["reaction"]
    )

    assert rxn is not None

    mol_1 = Chem.AddHs(
        Chem.MolFromSmiles(smiles_1)
    )

    mol_2 = Chem.AddHs(
        Chem.MolFromSmiles(smiles_2)
    )

    return rxn.RunReactants(
        (
            mol_1,
            mol_2,
        )
    )


# =============================================================================
# Structure
# =============================================================================

def test_expected_active_reactions():
    assert set(REACTIONS) == EXPECTED_REACTIONS


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_required_fields_exist(reaction_name):
    reaction = REACTIONS[reaction_name]

    required = {
        "same_reactants",
        "reactant_1",
        "product",
        "delete_atom",
        "reaction",
        "reference",
        "comments",
    }

    assert required.issubset(
        reaction.keys()
    )


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_reaction_smarts_parses(reaction_name):
    rxn = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    assert rxn is not None


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_reaction_passes_registry_validation(
    reaction_name,
):
    errors = registry._validate_reaction_smarts(
        reaction_name,
        REACTIONS[reaction_name],
    )

    assert errors == []


# =============================================================================
# Metadata
# =============================================================================

def test_amino_acid_homopolymerization_metadata():
    reaction = REACTIONS[AMINO_HOMO]

    assert reaction["same_reactants"] is True
    assert reaction["reactant_1"] == "amino_acid"

    assert reaction.get(
        "reactant_2"
    ) is None

    assert reaction["product"] == "polyamide_chain"
    assert reaction["delete_atom"] is True


def test_amino_acid_copolymerization_metadata():
    reaction = REACTIONS[AMINO_CO]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "amino_acid"
    assert reaction["reactant_2"] == "amino_acid"
    assert reaction["product"] == "polyamide_chain"
    assert reaction["delete_atom"] is True


def test_diamine_diacid_metadata():
    reaction = REACTIONS[DIAMINE_ACID]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "di_amine"

    assert (
        reaction["reactant_2"]
        == "di_carboxylic_acid"
    )

    assert reaction["product"] == "polyamide_chain"
    assert reaction["delete_atom"] is True


def test_diamine_diacid_halide_metadata():
    reaction = REACTIONS[DIAMINE_HALIDE]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "di_amine"

    assert (
        reaction["reactant_2"]
        == "di_carboxylic_acid_halide"
    )

    assert reaction["product"] == "polyamide_chain"
    assert reaction["delete_atom"] is True


def test_caprolactam_hydrolysis_metadata():
    reaction = REACTIONS[
        CAPROLACTAM_HYDROLYSIS
    ]

    assert reaction["same_reactants"] is False
    assert reaction["reactant_1"] == "water"
    assert reaction["reactant_2"] == "lactam"
    assert reaction["product"] == "polyamide_chain"
    assert reaction["delete_atom"] is False


# =============================================================================
# Initiator maps
# =============================================================================

@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_reactions_contain_initiator_maps_1_and_2(
    reaction_name,
):
    reaction = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    reactants = [
        Chem.Mol(
            reaction.GetReactantTemplate(i)
        )
        for i in range(
            reaction.GetNumReactantTemplates()
        )
    ]

    products = [
        Chem.Mol(
            reaction.GetProductTemplate(i)
        )
        for i in range(
            reaction.GetNumProductTemplates()
        )
    ]

    reactant_maps = (
        registry._atom_maps_in_templates(
            reactants
        )
    )

    product_maps = (
        registry._atom_maps_in_templates(
            products
        )
    )

    assert {
        1,
        2,
    }.issubset(
        reactant_maps
    )

    assert {
        1,
        2,
    }.issubset(
        product_maps
    )


@pytest.mark.parametrize(
    "reaction_name",
    sorted(EXPECTED_REACTIONS),
)
def test_products_contain_initiator_bond(
    reaction_name,
):
    reaction = rdChemReactions.ReactionFromSmarts(
        REACTIONS[reaction_name]["reaction"]
    )

    products = [
        Chem.Mol(
            reaction.GetProductTemplate(i)
        )
        for i in range(
            reaction.GetNumProductTemplates()
        )
    ]

    assert registry._has_bond_between_atom_maps(
        products,
        1,
        2,
    ) is True


# =============================================================================
# Amino-acid polycondensation
# =============================================================================

def test_glycine_homopolycondensation_executes():
    products = run_reaction(
        AMINO_HOMO,
        "NCC(=O)O",
        "NCC(=O)O",
    )

    assert len(products) > 0


def test_glycine_condensation_produces_two_product_fragments():
    """
    The SMARTS generates the amide-containing product plus the
    eliminated water fragment.
    """
    products = run_reaction(
        AMINO_HOMO,
        "NCC(=O)O",
        "NCC(=O)O",
    )

    assert len(
        products[0]
    ) == 2


def test_glycine_condensation_products_are_sanitizable():
    products = run_reaction(
        AMINO_HOMO,
        "NCC(=O)O",
        "NCC(=O)O",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


# =============================================================================
# Diamine + dicarboxylic acid
# =============================================================================

def test_diamine_diacid_polycondensation_executes():
    products = run_reaction(
        DIAMINE_ACID,
        "NCCN",
        "O=C(O)CCC(=O)O",
    )

    assert len(products) > 0


def test_diamine_diacid_products_are_sanitizable():
    products = run_reaction(
        DIAMINE_ACID,
        "NCCN",
        "O=C(O)CCC(=O)O",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_diamine_diacid_reaction_produces_byproduct():
    products = run_reaction(
        DIAMINE_ACID,
        "NCCN",
        "O=C(O)CCC(=O)O",
    )

    assert len(
        products[0]
    ) == 2


# =============================================================================
# Diamine + acid halide
# =============================================================================

def test_diamine_acid_halide_polycondensation_executes():
    products = run_reaction(
        DIAMINE_HALIDE,
        "NCCN",
        "O=C(Cl)CCC(=O)Cl",
    )

    assert len(products) > 0


def test_diamine_acid_halide_products_are_sanitizable():
    products = run_reaction(
        DIAMINE_HALIDE,
        "NCCN",
        "O=C(Cl)CCC(=O)Cl",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_diamine_acid_halide_reaction_produces_byproduct():
    products = run_reaction(
        DIAMINE_HALIDE,
        "NCCN",
        "O=C(Cl)CCC(=O)Cl",
    )

    assert len(
        products[0]
    ) == 2


# =============================================================================
# Hydrolytic initiation of caprolactam
# =============================================================================

def test_caprolactam_hydrolysis_executes():
    products = run_reaction(
        CAPROLACTAM_HYDROLYSIS,
        "O",
        "O=C1CCCCCN1",
    )

    assert len(products) > 0


def test_caprolactam_hydrolysis_product_is_sanitizable():
    products = run_reaction(
        CAPROLACTAM_HYDROLYSIS,
        "O",
        "O=C1CCCCCN1",
    )

    for product in products[0]:
        Chem.SanitizeMol(product)


def test_caprolactam_hydrolysis_produces_single_open_chain_product():
    products = run_reaction(
        CAPROLACTAM_HYDROLYSIS,
        "O",
        "O=C1CCCCCN1",
    )

    assert len(
        products[0]
    ) == 1


# =============================================================================
# Negative chemistry tests
# =============================================================================

def test_diamine_diacid_reaction_does_not_accept_simple_alcohol():
    products = run_reaction(
        DIAMINE_ACID,
        "CCO",
        "O=C(O)CCC(=O)O",
    )

    assert len(products) == 0


def test_diamine_acid_halide_reaction_does_not_accept_dicarboxylic_acid():
    products = run_reaction(
        DIAMINE_HALIDE,
        "NCCN",
        "O=C(O)CCC(=O)O",
    )

    assert len(products) == 0