from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, Draw
import pandas as pd
from pathlib import Path
import os


def is_number_in_set(set_of_tuples, reactant):
    print(set_of_tuples)
    for atom in reactant.GetAtoms():
        if atom.GetAtomMapNum() == 101 or atom.GetAtomMapNum() == 102:
            x = atom.GetIdx()
            for t in set_of_tuples:
                if x in t:
                    print(f"Found matching atom index {x} in tuple {t}")
                    return t
    raise ValueError("No matching atom found in the provided set of tuples.")


def smart_mapping(reactant, smarts_template, match_tuple):
    if not match_tuple:
        return
    SMARTS_index_to_Map_Number = {}
    for smarts_atom in smarts_template.GetAtoms():
        map_num = smarts_atom.GetAtomMapNum()
        if map_num != 0:
            SMARTS_index_to_Map_Number[smarts_atom.GetIdx()] = map_num
    for smarts_pos, map_num in SMARTS_index_to_Map_Number.items():
        if smarts_pos < len(match_tuple):
            atom_index_in_mol = match_tuple[smarts_pos]
            atom = reactant.GetAtomWithIdx(atom_index_in_mol)
            atom.SetAtomMapNum(map_num)


def is_consecutive(num_list):
    if not num_list:
        return False
    return len(set(num_list)) == len(num_list) and max(num_list) - min(num_list) + 1 == len(num_list)


def prepare_paths(cache):
    csv_cache = Path(cache) / "reactant_product_sets" / "csv"
    os.makedirs(csv_cache, exist_ok=True)
    return csv_cache


def build_reaction(rxn_smarts):
    return AllChem.ReactionFromSmarts(rxn_smarts)


def build_reactants(reactant_smiles_1, reactant_smiles_2):
    mol_reactant_1 = Chem.MolFromSmiles(reactant_smiles_1)
    mol_reactant_1 = Chem.AddHs(mol_reactant_1)
    mol_reactant_2 = Chem.MolFromSmiles(reactant_smiles_2)
    mol_reactant_2 = Chem.AddHs(mol_reactant_2)
    return mol_reactant_1, mol_reactant_2


def process_reactions(rxn, csv_cache, reaction_tuple, key=None, molecule_and_csv_path_dict=None, delete_atoms=False):
    if molecule_and_csv_path_dict is None:
        molecule_and_csv_path_dict = {}
    total_products = 0
    mols = []
    output = ""

    for j, pair in enumerate(reaction_tuple):
        r1, r2 = Chem.Mol(pair[0]), Chem.Mol(pair[1])
        products = rxn.RunReactants((r1, r2))
        matches_1 = r1.GetSubstructMatches(rxn.GetReactants()[0])
        matches_2 = r2.GetSubstructMatches(rxn.GetReactants()[1])

        for product_set in products:
            try:
                del df
            except:
                pass
            df = pd.DataFrame(columns=["reactant_index", "product_idx"])
            first_shell = []
            initatiator_idxs = []
            mapping_dict = {}

            for i, product in enumerate(product_set):
                num_total_atoms = 0
                if i == 1:
                    num_total_atoms = product_set[0].GetNumAtoms()
                for atom in product.GetAtoms():
                    if atom.HasProp("react_idx"):
                        r_idx = atom.GetIntProp("react_idx")
                        a_idx = atom.GetIntProp("react_atom_idx")
                        r = (r1, r2)[r_idx]
                        r_atom = r.GetAtomWithIdx(a_idx)
                        map_num = atom.GetIdx() + 1 + num_total_atoms + 100
                        r_atom.SetAtomMapNum(map_num)
                        atom.SetAtomMapNum(map_num)

            reactant_combined = Chem.CombineMols(r1, r2)
            product_combined = Chem.CombineMols(*product_set)

            for r_atom in reactant_combined.GetAtoms():
                for p_atom in product_combined.GetAtoms():
                    if r_atom.GetAtomMapNum() == p_atom.GetAtomMapNum():
                        new_row = pd.DataFrame([{"reactant_index": r_atom.GetIdx(), "product_idx": p_atom.GetIdx()}])
                        df = pd.concat([df, new_row], ignore_index=True)
                        break

            num_reactant_atoms = reactant_combined.GetNumAtoms()
            num_product_atoms = product_combined.GetNumAtoms()
            reactant_mapped = df["reactant_index"].notna().sum()
            product_mapped = df["product_idx"].notna().sum()

            if reactant_mapped != product_mapped:
                raise ValueError(
                    "Mismatch in mapped atom counts between columns: "
                    f"reactant_index mapped={reactant_mapped}, product_idx mapped={product_mapped}"
                )

            if reactant_mapped != num_reactant_atoms:
                df.to_csv(csv_cache / f"debug_mapping_{key}_{total_products}.csv", index=False)
                raise ValueError(
                    "Mapping does not cover all reactant atoms: "
                    f"reactant atoms={num_reactant_atoms}, mapped={reactant_mapped}"
                )

            if product_mapped != num_product_atoms:
                raise ValueError(
                    "Mapping does not cover all product atoms: "
                    f"product atoms={num_product_atoms}, mapped={product_mapped}"
                )

            if is_consecutive(df["reactant_index"].tolist()) is False or is_consecutive(df["product_idx"].tolist()) is False:
                raise ValueError(
                    "Mapping indices are not consecutive: "
                    f"reactant indices={df['reactant_index'].tolist()}, product indices={df['product_idx'].tolist()}"
                )

            sub_set1 = is_number_in_set(matches_1, r1)
            smart_mapping(reactant=r1, smarts_template=rxn.GetReactants()[0], match_tuple=sub_set1 if sub_set1 else ())

            sub_set2 = is_number_in_set(matches_2, r2)
            smart_mapping(reactant=r2, smarts_template=rxn.GetReactants()[1], match_tuple=sub_set2 if sub_set2 else ())

            reactant_combined = Chem.CombineMols(r1, r2)
            for atom in reactant_combined.GetAtoms():
                if atom.GetAtomMapNum() <= 99:
                    first_shell.append(atom.GetIdx())
                if atom.GetAtomMapNum() in [1, 2]:
                    initatiator_idxs.append(atom.GetIdx())

            byproduct_indexs = []
            if delete_atoms:
                frags = rdmolops.GetMolFrags(product_combined, asMols=True)
                smallest_mol = min(frags, key=lambda m: m.GetNumAtoms())
                for atom in smallest_mol.GetAtoms():
                    byproduct_map_number = atom.GetAtomMapNum()
                    for p_atom in product_combined.GetAtoms():
                        if p_atom.GetAtomMapNum() == byproduct_map_number:
                            byproduct_indexs.append(p_atom.GetIdx())
                            break

            total_products += 1
            mols.append(reactant_combined)
            mols.append(product_combined)
            if key not in molecule_and_csv_path_dict:
                molecule_and_csv_path_dict[key] = {}
            sub_dict = molecule_and_csv_path_dict[key][total_products] = {}
            sub_dict["reactant"] = reactant_combined
            sub_dict["product"] = product_combined
            sub_dict["csv_path"] = csv_cache / f"reaction_{key}_{total_products}.csv"
            sub_dict["mapping_dataframe"] = df
            sub_dict["delete_atoms"] = delete_atoms

            for r_idx, p_idx in mapping_dict.items():
                new_row = pd.DataFrame([{"reactant_index": r_idx, "product_idx": p_idx}])
                df = pd.concat([df, new_row], ignore_index=True)

            first_shell_column = pd.Series(first_shell, name="first_shell")
            initatiator_idxs_column = pd.Series(initatiator_idxs, name="initiators")
            if len(initatiator_idxs) != 2:
                raise ValueError(f"Expected 2 initiators, got {len(initatiator_idxs)}: {initatiator_idxs}")
            byproduct_indexs_column = pd.Series(byproduct_indexs, name="byproduct_indices")
            df_combined = pd.concat([df, first_shell_column, initatiator_idxs_column, byproduct_indexs_column], axis=1)
            df_combined.to_csv(csv_cache / f"reaction_{key}_{total_products}.csv", index=False)
            print(f"Saved reaction {key}_{total_products} to CSV")

    return mols, output, molecule_and_csv_path_dict


def save_output(output, out_path):
    with open(out_path, "w") as f:
        f.write(output)


def save_grid_image(mols, csv_cache, key=None):
    if mols:
        img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(900, 900))
        img_path = Path(csv_cache) / (f"reaction_grid_{key}.png" if key is not None else "reaction_grid.png")
        img.save(str(img_path))
        return img_path
    return None


def reaction_tuples(same_reactants, mol_reactant_1, mol_reactant_2):
    if not mol_reactant_2:
        mol_reactant_2 = mol_reactant_1
    if same_reactants:
        return [[mol_reactant_1, mol_reactant_2], [mol_reactant_2, mol_reactant_1]]
    else:
        return [[mol_reactant_1, mol_reactant_2]]


def processing_dict(detected_reactions, cache):
    molecule_and_csv_path_dict = {}
    csv_cache = prepare_paths(cache)
    for key in detected_reactions:
        reaction_dict = detected_reactions[key]
        rxn_smarts = reaction_dict["reaction"]
        reactant_smiles_1 = reaction_dict["monomer_1"]["smiles"]
        same_reactants = reaction_dict["same_reactants"]
        if same_reactants:
            reactant_smiles_2 = reactant_smiles_1
        else:
            reactant_smiles_2 = reaction_dict["monomer_2"]["smiles"]
        delete_atoms = reaction_dict.get("delete_atom", False)
        rxn = build_reaction(rxn_smarts)
        mol_reactant_1, mol_reactant_2 = build_reactants(reactant_smiles_1, reactant_smiles_2)
        reaction_tuple = reaction_tuples(same_reactants, mol_reactant_1, mol_reactant_2)

        mols, output, molecule_and_csv_path_dict = process_reactions(rxn, csv_cache, reaction_tuple, key, molecule_and_csv_path_dict, delete_atoms=delete_atoms)
        save_output(output, f"reaction_mapping_output_{key}.txt")
        save_grid_image(mols, csv_cache, key)
    return molecule_and_csv_path_dict


def run_all(cache, rxn_smarts, reactant_smiles_1, reactant_smiles_2):
    molecule_and_csv_path_dict = {}
    csv_cache = prepare_paths(cache)
    rxn = build_reaction(rxn_smarts)
    mol_reactant_1, mol_reactant_2 = build_reactants(reactant_smiles_1, reactant_smiles_2)
    reaction_tuple = [[mol_reactant_1, mol_reactant_2]]
    delete_atoms = True
    mols, output, molecule_and_csv_path_dict = process_reactions(rxn, csv_cache, reaction_tuple, None, molecule_and_csv_path_dict, delete_atoms)
    save_output(output, "reaction_mapping_output.txt")
    save_grid_image(mols, csv_cache, None)
    return molecule_and_csv_path_dict


if __name__ == "__main__":
    detected_reactions = {
        1: {
            "reaction_name": "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)",
            "same_reactants": True,
            "reactant_1": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "reference": {
                "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                "reaction_and_mechanism": [
                    "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                    "https://pubs.acs.org/doi/10.1021/ed073pA312",
                ],
            },
            "monomer_1": {
                "smiles": "O=C(O)c1cc(O)c(Cl)c(C(=O)O)c1",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 2,
            },
        },
        2: {
            "reaction_name": "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)",
            "same_reactants": True,
            "reactant_1": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "reference": {
                "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                "reaction_and_mechanism": [
                    "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                    "https://pubs.acs.org/doi/10.1021/ed073pA312",
                ],
            },
            "monomer_1": {
                "smiles": "O=C(O)CCCC(O)CCCO",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 2,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1,
            },
        },
        3: {
            "reaction_name": "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)",
            "same_reactants": False,
            "reactant_1": "hydroxy_carboxylic_acid",
            "reactant_2": "hydroxy_carboxylic_acid",
            "product": "polyester_chain",
            "delete_atom": True,
            "reaction": "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
            "reference": {
                "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
                "reaction_and_mechanism": [
                    "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                    "https://pubs.acs.org/doi/10.1021/ed073pA312",
                ],
            },
            "monomer_1": {
                "smiles": "O=C(O)c1cc(O)c(Cl)c(C(=O)O)c1",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 1,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 2,
            },
            "monomer_2": {
                "smiles": "O=C(O)CCCC(O)CCCO",
                "functionality_type": "di_different",
                "functional_group_name": "hydroxy_carboxylic_acid",
                "functional_group_smarts_1": "[OX2H1;!$(OC=*):1]",
                "functional_count_1": 2,
                "functional_group_smarts_2": "[CX3:2](=[O])[OX2H1]",
                "functional_count_2": 1,
            },
        },
    }

    cache = "C:\\Users\\Janitha\\Documents\\GitHub\\reaction_lammps_mupt\\cache\\00_cache"
    rxn_smarts = "[O;!$(OC=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]"
    reactant_smiles_1 = "O=C(O)c1cc(O)c(Cl)c(C(=O)O)c1"
    reactant_smiles_2 = "O=C(O)CCCC(O)CCCO"

    molecule_dict_csv_path_dict = run_all(cache, rxn_smarts, reactant_smiles_1, reactant_smiles_2)
    molecule_dict_csv_path_dict = processing_dict(detected_reactions, cache)
    import pprint
    pprint.pprint(molecule_dict_csv_path_dict)
    print(molecule_dict_csv_path_dict)
