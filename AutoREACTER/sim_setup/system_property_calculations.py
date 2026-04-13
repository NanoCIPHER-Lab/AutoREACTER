import math
from rdkit.Chem import Descriptors

from AutoREACTER.input_parser import SimulationSetup

# Physical constants
N_A = 6.02214076e23      # Avogadro's constant (mol⁻¹)
CM_2_A3 = 1e24           # Conversion factor: 1 cm³ = 10²⁴ Å³


class NoneMonomerError(Exception):
    """Custom exception raised when a monomer is expected to have an RDKit molecule but does not."""
    pass


class SystemPropertyCalculations:
    """
    Computes essential physical and stoichiometric properties for each simulation replica.
    
    This class is responsible for two main tasks:
      1. Populating monomer-level properties (atom count, molecular weight) from RDKit objects.
      2. Calculating replica-level properties including monomer counts (when using ratio mode),
         total system mass, target volume, and initial (low-density) box dimensions.
    
    All calculations are performed in-place on the provided SimulationSetup object.
    The class follows a clear separation between monomer-level and replica-level logic.
    """

    def __init__(self, simulation_setup: SimulationSetup) -> None:
        """
        Initialize with a SimulationSetup object that will be mutated with calculated values.
        
        Args:
            simulation_setup: Fully parsed SimulationSetup containing monomers,
                              replicas, composition method, and target densities.
        """
        self.simulation_setup = simulation_setup

    def process_all(self) -> SimulationSetup:
        """
        Run the complete property calculation pipeline.
        
        This is the main public method called by SimulationSetupManager.
        It populates monomer properties first, then calculates all replica-specific
        values (counts, box dimensions, etc.).
        
        Returns:
            SimulationSetup: The same object reference now populated with all calculated
                             properties (monomer counts, initial_box_length, etc.).
        """
        self._populate_monomer_properties()
        self._calculate_replica_properties()
        return self.simulation_setup

    def _populate_monomer_properties(self) -> None:
        """
        Populate each active monomer with basic structural properties derived from its RDKit molecule.
        
        For every monomer whose status is True:
          - num_atoms is set to the heavy atom count.
          - molecular_weight is set using RDKit's Descriptors.MolWt (g/mol).
        
        Raises:
            NoneMonomerError: If a monomer is marked active but has no RDKit Mol object.
        """
        print("\n[INFO] Populating monomer properties based on RDKit Mol objects...")

        for monomer in self.simulation_setup.monomers:
            if not monomer.status:
                continue

            if monomer.rdkit_mol is None:
                raise NoneMonomerError(
                    f"Monomer with ID {monomer.id} has no RDKit Mol object."
                )

            monomer.num_atoms = monomer.rdkit_mol.GetNumAtoms()
            monomer.molecular_weight = Descriptors.MolWt(monomer.rdkit_mol)

    def _calculate_replica_properties(self) -> None:
        """
        Compute all properties required for each replica (monomer counts and box dimensions).
        
        This method first builds a lookup dictionary of active monomers (keyed by name
        to match ratio definitions), then iterates over every replica and calls the
        appropriate calculation helpers based on the composition_method.
        """
        # Build lookup of only active monomers for fast access during ratio calculations
        active_monomers = {
            monomer.name: monomer
            for monomer in self.simulation_setup.monomers
            if monomer.status
        }

        for replica in self.simulation_setup.replicas:
            # 1. Determine integer monomer counts when user supplied stoichiometric ratios
            if self.simulation_setup.composition_method == "ratio":
                self._get_monomer_counts(replica, active_monomers)

            # 2. Calculate initial (sparse) and target box dimensions from mass and density
            self._calculate_box_dimensions(replica, active_monomers)

    def _get_monomer_counts(self, replica, active_monomers: dict) -> None:
        """
        Convert monomer ratios into integer counts that sum close to the requested total_atoms.
        
        The algorithm calculates a multiplier so the weighted atom count matches the
        target total_atoms, then ceilings each contribution (minimum of 1 monomer).
        The resulting counts are stored both on the replica and back on the monomer
        objects for later use in file writing.
        
        Args:
            replica: The Replica object being processed.
            active_monomers: Dictionary mapping monomer names to Monomer objects.
        
        Raises:
            ValueError: If total_atoms is not set or if no valid monomers are present.
        """
        if replica.total_atoms is None:
            raise ValueError(f"Replica '{replica.tag}' is missing 'total_atoms'.")

        # Calculate total atoms per "unit" of the stoichiometric ratio
        atoms_per_unit = sum(
            ratio * active_monomers[m_name].num_atoms
            for m_name, ratio in replica.monomer_ratios.items()
            if m_name in active_monomers
        )

        if atoms_per_unit == 0:
            raise ValueError(f"Replica '{replica.tag}' has no active monomers to calculate ratios.")

        multiplier = replica.total_atoms / atoms_per_unit

        if replica.monomer_counts is None:
            replica.monomer_counts = {}

        for m_name, ratio in replica.monomer_ratios.items():
            if m_name in active_monomers:
                calculated_count = max(1, math.ceil(ratio * multiplier))
                replica.monomer_counts[m_name] = calculated_count

                # Also record this count on the monomer for easier lookup later
                if active_monomers[m_name].count is None:
                    active_monomers[m_name].count = {}
                active_monomers[m_name].count[replica.tag] = calculated_count

    def _calculate_box_dimensions(self, replica, active_monomers: dict) -> None:
        """
        Compute the initial (low-density) simulation box size based on total mass and target density.
        
        The method first calculates the total mass of the system using monomer counts
        and molecular weights. It then derives:
          - The target volume at the user-specified density.
          - An initial sparse box at 25% of the target density (to avoid atomic overlaps
            during the early minimization/equilibration stages).
        
        All results are stored directly on the replica object.
        
        Args:
            replica: The Replica object being processed.
            active_monomers: Dictionary of active monomers for quick property lookup.
        
        Raises:
            ValueError: If monomer_counts have not been calculated yet.
        """
        if replica.monomer_counts is None:
            raise ValueError(f"Replica '{replica.tag}' lacks monomer counts.")

        # Calculate total system mass in grams
        total_mass_g = 0.0
        for m_name, count in replica.monomer_counts.items():
            if m_name in active_monomers:
                mass_g = (count * active_monomers[m_name].molecular_weight) / N_A
                total_mass_g += mass_g

        # Target volume at user-defined density
        target_volume_cm3 = total_mass_g / replica.density
        target_volume_A3 = target_volume_cm3 * CM_2_A3

        # Initial sparse box (25% of target density) to reduce starting overlaps
        initial_density = replica.density / 4.0
        initial_volume_cm3 = total_mass_g / initial_density
        initial_volume_A3 = initial_volume_cm3 * CM_2_A3

        initial_box_length = round(math.pow(initial_volume_A3, 1.0 / 3.0), 2)

        replica.initial_box_volume = initial_volume_A3
        replica.initial_box_length = initial_box_length
