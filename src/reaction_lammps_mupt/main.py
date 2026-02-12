import json
from pathlib import Path
from initialization import initialize
from input_parser import InputParser
from cache import GetCacheDir, RunDirectoryManager, RetentionCleanup
from detectors.detector import Detector

"""
Main script for initializing the reaction LAMMPS MUPT pipeline.
"""
initialize()
get_cache_dir = GetCacheDir()
cache_dir = get_cache_dir.staging_dir
print(f"Base cache directory: {cache_dir}")
dated_cache_dir = RunDirectoryManager.make_dated_run_dir(get_cache_dir.cache_base_dir, chdir_to="none")
print(f"Staging directory created at: {dated_cache_dir}")

# optional use
"""
TODO: Use this in CMD line interface and in the future when we want to support both ways
    of running the code (with a json file or with cmd line arguments)"""

# RunDirectoryManager.copy_into_run(cache_dir, dated_cache_dir)



if __name__ == "__main__":
    current_dir = Path(__file__).resolve().parents[2]  # same as parent.parent.parent, but clearer
    input_file = current_dir / "inputs.json"

    with input_file.open("r", encoding="utf-8") as f:
        inputs = json.load(f)

    # modify this to a cmd line argument later and to support this way as well
    parser = InputParser()
    validated_inputs = parser.validate_inputs(inputs)
    print(validated_inputs)
    # RetentionCleanup.run(GetCacheDir().cache_base_dir)  # Optional: Clean up old cache files based on retention policy.
    detector = Detector(inputs)
    detected_reactions = detector.reactions_dict
    non_monomers = detector.non_reactants_list
    RunDirectoryManager.copy_into_run(cache_dir, dated_cache_dir)

