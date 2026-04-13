from pathlib import Path
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.sim_setup.system_property_calculations import SystemPropertyCalculations
from AutoREACTER.sim_setup.writers.writer import Writer
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles

class SimulationSetupManager:


    def setup_and_write_simulation(self, setup: SimulationSetup, reacter_files: REACTERFiles, run_dir: Path) -> SimulationSetup:
        calculator = SystemPropertyCalculations(setup)
        updated_setup = calculator.process_all()
        # Debugging output to verify the updated setup
        # for replica in updated_setup.replicas:
        #     print(f"\n[INFO] Replica: {replica.tag}")
        #     print(f"  - Calculated Monomer Counts: {replica.monomer_counts}")
        #     print(f"  - Initial Box Length: {replica.initial_box_length:.2f} Å")
        #     print(f"  - Target Density: {replica.density} g/cm³")

        self.generate_input_files(setup=updated_setup, reacter_files=reacter_files, run_dir=run_dir)
        print(f"\n[SUCCESS] All 5 simulation stages written to: {run_dir / 'LAMMPS_input_files'}")
        return updated_setup

    def generate_input_files(self, setup: SimulationSetup, reacter_files: REACTERFiles, run_dir: Path) -> None:
        writer = Writer(reacter_files=reacter_files)
        writer.write_all_files(run_dir, setup)