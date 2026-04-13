import os
from pathlib import Path

from AutoREACTER.sim_setup.writers.lammps_settings import LammpsInitialSettings, LammpsSettings
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles
from AutoREACTER.sim_setup.writers.densification_writer import DensificationWriter

class Writer:
    def __init__(self, reacter_files: REACTERFiles):
        self.reacter_files = reacter_files
        self.lammps_initial_setup = LammpsInitialSettings(reacter_files)
        self.settings = self.lammps_initial_setup.get_LUNAR_lammps_settings()
        self.densification_writer = DensificationWriter(self.settings, reacter_files)

    def write_all_files(self, run_dir: Path, replica):
        lammps_dir = run_dir / "LAMMPS_files"
        lammps_dir.mkdir(parents=True, exist_ok=True)
        self.densification_writer.write_files(lammps_dir, replica)
        self.densification_writer.write_densification_script(lammps_dir, replica)