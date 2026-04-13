import os
from pathlib import Path

from AutoREACTER.sim_setup.writers.lammps_settings import LammpsInitialSettings
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles
from AutoREACTER.sim_setup.writers.densification_writer import DensificationWriter
from AutoREACTER.input_parser import SimulationSetup

class Writer:

    def __init__(self, reacter_files: REACTERFiles):
        self.reacter_files = reacter_files
        self.lammps_initial_setup = LammpsInitialSettings(reacter_files)
        self.settings = self.lammps_initial_setup.get_LUNAR_lammps_settings()
        self.densification_writer = DensificationWriter(self.settings, reacter_files)

    def write_all_files(self, run_dir: Path, simulation_setup: SimulationSetup) -> None:
        lammps_dir = run_dir / "LAMMPS_files"
        lammps_dir.mkdir(parents=True, exist_ok=True)

        replicas = simulation_setup.replicas
        for replica in replicas:
            self.densification_writer.write_lammps_densification_files(lammps_dir, replica)
        