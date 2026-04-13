from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.sim_setup.system_property_calculations import SystemPropertyCalculations

class SimulationSetupManager:
    def populate_physical_parameters(self, setup: SimulationSetup) -> SimulationSetup:
        
        calculator = SystemPropertyCalculations(setup)
        updated_setup = calculator.process_all()

        for replica in updated_setup.replicas:
            print(f"\n[INFO] Replica: {replica.tag}")
            print(f"  - Calculated Counts: {replica.monomer_counts}")
            print(f"  - Initial Box Volume: {replica.initial_box_volume:.2f} Å³ (Length: {replica.initial_box_length:.2f} Å)")
        
        return updated_setup