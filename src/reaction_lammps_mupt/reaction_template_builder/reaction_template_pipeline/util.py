from rdkit import Chem

def compare_products(dict, _prod2):
    prod2 = Chem.Mol(_prod2)
    for key, dict_2 in dict.items():
        prod = dict_2.get("product")
        if prod is not None:
            prod1 = Chem.Mol(prod)

            for atom in prod1.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in prod2.GetAtoms():
                atom.SetAtomMapNum(0)

            smi1 = Chem.MolToSmiles(prod1, canonical=True)
            smi2 = Chem.MolToSmiles(prod2, canonical=True)

            if smi1 == smi2:
                print("DUPLICATE PRODUCT FOUND")
                return False
    return True

def format_detected_reactions_dict(detected_reactions):
    reactions_names = ""
    smart_references = ""
    mechanism_references = ""
    formtted_dict = {}

    for key, reaction in detected_reactions.items():
        reactions_name = reaction.get("reaction_name")
        if reactions_name in reactions_names:
            continue
        if reactions_names == "":
            reactions_names += reactions_name
        else:
            reactions_names += f", {reactions_name}"

        reference = reaction.get("reference")
        if reference is None:
            continue

        smarts_ref = reference.get("smarts")
        mech_refs = reference.get("reaction_and_mechanism")

        if smarts_ref:
            if smarts_ref not in smart_references:
                if smart_references == "":
                    smart_references += smarts_ref
                else:
                    smart_references += f", {smarts_ref}"

        if mech_refs:
            mech_block = ", ".join(mech_refs)
            if mech_block not in mechanism_references:
                if mechanism_references == "":
                    mechanism_references += mech_block
                else:
                    mechanism_references += f", {mech_block}"

    formtted_dict = {
        "reactions_names": reactions_names,
        "smart_references": smart_references,
        "mechanism_references": mechanism_references,
    }

    return formtted_dict


