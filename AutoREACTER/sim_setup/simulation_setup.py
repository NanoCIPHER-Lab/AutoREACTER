from pathlib import Path
from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.sim_setup.system_property_calculations import SystemPropertyCalculations
from AutoREACTER.sim_setup.writers.writer import Writer
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles

class SimulationSetupManager:
    def populate_physical_parameters(self, setup: SimulationSetup, reacter_files: REACTERFiles, run_dir: Path) -> SimulationSetup:
        calculator = SystemPropertyCalculations(setup)
        updated_setup = calculator.process_all()

        for replica in updated_setup.replicas:
            print(f"\n[INFO] Replica: {replica.tag}")
            print(f"  - Calculated Counts: {replica.monomer_counts}")
            print(f"  - Initial Box Volume: {replica.initial_box_volume:.2f} Å³ (Length: {replica.initial_box_length:.2f} Å)")

        writer = Writer(reacter_files)
        writer.write_all_files(run_dir, updated_setup)
        
        print(f"\n[SUCCESS] LAMMPS densification packages written to: {run_dir}")
        
        return updated_setup
    
    def write_densification_files(self, setup: SimulationSetup, reacter_files: REACTERFiles, run_dir: Path) -> None:
        writer = Writer(reacter_files)
        writer.write_all_files(run_dir, setup)
        print(f"\n[SUCCESS] LAMMPS densification packages written to: {run_dir}")