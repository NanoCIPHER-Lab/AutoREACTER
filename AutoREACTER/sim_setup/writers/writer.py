from AutoREACTER.sim_setup.writers.lammps_settings import LammpsInitialSettings, LammpsSettings
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles


class Writer:
    def __init__(self, reacter_files: REACTERFiles):
        self.reacter_files = reacter_files
        self.lammps_initial_setup = LammpsInitialSettings(reacter_files)
        self.settings: LammpsSettings = self.lammps_initial_setup.get_LUNAR_lammps_settings()