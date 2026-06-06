"""
Merge Builder Module

Generates the text-based merge_input.txt configuration file required 
by LUNAR's bond_react_merge.py script.
"""

from pathlib import Path
from typing import TYPE_CHECKING, Any

# Local utility imports
from AutoREACTER.reaction_preparation.lunar_client.lunar_utils import normalize_path, get_ending_integer

if TYPE_CHECKING:
    # Forward declaration for type hinting to avoid circular imports during refactor
    from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import All2LMPResult

def write_bond_react_merge_input(
    cache_bond_react_merge: Path,
    cache_all2lmp: Path,
    all2lmp_results: list['All2LMPResult']
) -> Path:
    """
    Generate the merge_input.txt file for bond_react_merge.py.
    
    This file tells LUNAR which data files are molecules vs. reaction
    templates and how they should be combined. Molecules get 'dataN' 
    tags while reactions get 'preN'/'postN' pairs.
    
    Args:
        cache_bond_react_merge: Directory to save the merge_input.txt
        cache_all2lmp: Directory containing the all2lmp output data files
        all2lmp_results: Results from the all2lmp conversion stage
        
    Returns:
        Path to the generated merge_input.txt file
    """
    cache_dir = Path(cache_bond_react_merge)
    cache_dir.mkdir(parents=True, exist_ok=True)

    merge_file = cache_dir / "merge_input.txt"

    # Header with column descriptions
    merge_files = f"""# anything following the "#" character will be ignored

{"# file-tag":<10}{"filename":<150}{"comment (required)"}
"""

    # Tag data molecules as data1, data2, etc.
    data_counter = 1
    for r in all2lmp_results:
        if not r.molecule:
            continue

        full_path = Path(cache_all2lmp) / r.all2lmp_data_file
        full_path = normalize_path(full_path)

        tag = f"data{data_counter}"
        comment = "# This datafile will have all coeffs in it"

        merge_files += f"{tag:<10}{full_path:<150}{comment}\n"
        data_counter += 1

    # Group reaction templates by ID and tag as preN/postN pairs
    reaction_pairs = {}

    for r in all2lmp_results:
        if r.molecule:
            continue

        rid = get_ending_integer(r.id)
        if rid is None:
            continue

        if rid not in reaction_pairs:
            reaction_pairs[rid] = {}

        if r.id.startswith("pre"):
            reaction_pairs[rid]["pre"] = r
        elif r.id.startswith("post"):
            reaction_pairs[rid]["post"] = r

    # Write pre/post pairs in order
    rxn_counter = 1
    for rid in sorted(reaction_pairs):
        pair = reaction_pairs[rid]

        pre = pair.get("pre")
        post = pair.get("post")

        if not pre or not post:
            raise ValueError(f"Incomplete reaction pair for reaction ID {rid}")

        pre_path = normalize_path(Path(cache_all2lmp) / pre.all2lmp_data_file)
        post_path = normalize_path(Path(cache_all2lmp) / post.all2lmp_data_file)

        merge_files += f"{f'pre{rxn_counter}':<10}{pre_path:<150}# for rxn{rxn_counter}\n"
        merge_files += f"{f'post{rxn_counter}':<10}{post_path:<150}# for rxn{rxn_counter}\n"

        rxn_counter += 1

    merge_files += f"\n# Specify the parent_directory of where to write results (optional)\n"

    with open(merge_file, "w") as f:
        f.write(merge_files)

    return merge_file