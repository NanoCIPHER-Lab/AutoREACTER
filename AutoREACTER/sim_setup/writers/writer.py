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

    def write_all_files(self, run_dir: Path, simulation_setup: SimulationSetup) -> None:
        sim_name = simulation_setup.simulation_name
        print(f"\n[INFO] Writing all files for simulation: {sim_name}")
        lammps_dir = run_dir / "LAMMPS_input_files"
        lammps_dir.mkdir(parents=True, exist_ok=True)

        replicas = simulation_setup.replicas
        for replica in replicas:
            sub_dir = lammps_dir / f"{sim_name}_{replica.tag}"
            sub_dir.mkdir(parents=True, exist_ok=True)
            DensificationWriter(
                out_dir=sub_dir, 
                settings=self.settings, 
                reacter_files=self.reacter_files, 
                replica=replica, 
                sim_name=sim_name
            )
        