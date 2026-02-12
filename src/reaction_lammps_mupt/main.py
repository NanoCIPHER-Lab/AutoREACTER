import json
from pathlib import Path
from initialization import initialize
from input_parser import InputParser
from cache import get_default_cache_dir, retention_cleanup, get_current_cache_dir
from cache import delete_default_cache_files as delete_cache_dir
from cache import copy_to_date_folder as Copy

initialize()
cache_dir = get_default_cache_dir()
final_cache_dir = get_current_cache_dir()


if __name__ == "__main__":
    current_dir = Path(__file__).resolve().parents[2]  # same as parent.parent.parent, but clearer
    input_file = current_dir / "inputs.json"

    with input_file.open("r", encoding="utf-8") as f:
        inputs = json.load(f)

    # modify this to a cmd line argument later and to support this way as well
    parser = InputParser()
    validated_inputs = parser.validate_inputs(inputs)
    print(validated_inputs)



# then will go to the reaction detector part 
# if we failed to generate reactions we say "No reactions detected with given monomers"
# we generate scripts without automatic reactions we request for templates from user and the map file -> we say this to lunar and The INPUT Script Generator

# then will come back to main.py file 

# Next part will be to go to the LUNAR for parameterization

# then will come back to main.py file 

# Then with all the data will go to INPUT Script Generator

# and finally will result in OUTPUT files as a ready Folder

# MAIN.PY name will changed in the future to another name 