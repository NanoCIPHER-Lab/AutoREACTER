import math
from rdkit.Chem import Descriptors
from AutoREACTER.input_parser import SimulationSetup

N_A = 6.02214076e23 # Avogadro's constant
CM_2_A3 = 1e24  # Conversion factor from cm^3 to Angstroms^3

class NoneMonomerError(Exception):
    """Custom exception for cases where a monomer is expected but not found."""
    pass

class SystemPropertyCalculations:

    def __init__(self, simulation_setup: SimulationSetup):
        self.simulation_setup = simulation_setup

    def process_all(self) -> SimulationSetup:
        self._populate_monomer_properties()
        self._calculate_replica_properties()
        
        # Fixed typo and properly returning the setup
        return self.simulation_setup
        
    def _populate_monomer_properties(self) -> None:
        print("\n[INFO] Populating monomer properties based on RDKit Mol objects...")
        for monomer in self.simulation_setup.monomers:
            if not monomer.status:
                continue
            
            if monomer.rdkit_mol is None:
                raise NoneMonomerError(f"Monomer with ID {monomer.id} has no RDKit Mol object.")
            
            monomer.num_atoms = monomer.rdkit_mol.GetNumAtoms()
            monomer.molecular_weight = Descriptors.MolWt(monomer.rdkit_mol)

    def _calculate_replica_properties(self) -> None:
        active_monomers = {}

        for monomer in self.simulation_setup.monomers:
            if monomer.status:
                # FIXED: Key this by name instead of id so it matches ratio lookups
                active_monomers[monomer.name] = monomer
        
        for replica in self.simulation_setup.replicas:
            # 1. Get monomer counts if ratio mode is used
            if self.simulation_setup.composition_method == "ratio":
                self._get_monomer_counts(replica, active_monomers)
            # 2. Calculate box dimensions based on monomer counts and properties
            self._calculate_box_dimensions(replica, active_monomers)

    def _get_monomer_counts(self, replica, active_monomers: dict) -> None:
        if replica.total_atoms is None:
            raise ValueError(f"Replica '{replica.tag}' is missing 'total_atoms'.")
        
        atoms_per_unit = 0.0
        for m_name, ratio in replica.monomer_ratios.items():
            if m_name in active_monomers:
                atoms_per_unit += ratio * active_monomers[m_name].num_atoms
                
        if atoms_per_unit == 0:
            raise ValueError(f"Replica '{replica.tag}' has no active monomers to calculate ratios.")
            
        multiplier = replica.total_atoms / atoms_per_unit
        
        if replica.monomer_counts is None:
            replica.monomer_counts = {}
            
        for m_name, ratio in replica.monomer_ratios.items():
            if m_name in active_monomers:
                calculated_count = max(1, math.ceil(ratio * multiplier))
                replica.monomer_counts[m_name] = calculated_count
                
                if active_monomers[m_name].count is None:
                    active_monomers[m_name].count = {}
                active_monomers[m_name].count[replica.tag] = calculated_count

    def _calculate_box_dimensions(self, replica, active_monomers: dict) -> None:
        if replica.monomer_counts is None:
            raise ValueError(f"Replica '{replica.tag}' lacks monomer counts.")

        total_mass_g = 0.0
        for m_name, count in replica.monomer_counts.items():
            if m_name in active_monomers:
                # Mass of this species = (Count * MW) / Avogadro
                mass_g = (count * active_monomers[m_name].molecular_weight) / N_A
                total_mass_g += mass_g
        
        target_volume_cm3 = total_mass_g / replica.density
        target_volume_A3 = target_volume_cm3 * CM_2_A3

        # Initial (Sparse) Dimensions
        initial_density = replica.density / 4.0
        initial_volume_cm3 = total_mass_g / initial_density
        initial_volume_A3 = initial_volume_cm3 * CM_2_A3
        initial_box_length = round(math.pow(initial_volume_A3, 1.0/3.0), 2)
        replica.initial_box_volume = initial_volume_A3
        replica.initial_box_length = initial_box_length