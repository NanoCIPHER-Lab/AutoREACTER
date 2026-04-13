from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.sim_setup.system_property_calculations import SystemPropertyCalculations

class SimulationSetupManager:
    def populate_physical_parameters(self, setup: SimulationSetup) -> SimulationSetup:
        
        calculator = SystemPropertyCalculations(setup)
        updated_setup = calculator.process_all()

        for replica in updated_setup.replicas:
            # FIXED: Changed replica.name to replica.tag 
            print(f"\n[INFO] Replica: {replica.tag}")
            
            # Print target atoms if ratio mode, otherwise note count mode
            target_atoms = replica.total_atoms if replica.total_atoms else "N/A (Count Mode)"
            print(f"  - Target Total Atoms: {target_atoms}")
            print(f"  - Calculated Counts: {replica.monomer_counts}")
            print(f"  - Target Box Volume: {replica.box_volume:.2f} Å³ (Length: {replica.box_length:.2f} Å)")
            print(f"  - Sparse Box Volume: {replica.initial_box_volume:.2f} Å³ (Length: {replica.initial_box_length:.2f} Å)")
        
        return updated_setup