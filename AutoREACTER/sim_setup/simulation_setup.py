from pathlib import Path

from AutoREACTER.input_parser import SimulationSetup
from AutoREACTER.sim_setup.system_property_calculations import SystemPropertyCalculations
from AutoREACTER.sim_setup.writers.writer import Writer
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles


class SimulationSetupManager:
    """
    Orchestrates the final stage of the AutoREACTER workflow: system property
    calculation and LAMMPS input file generation.
    
    This manager coordinates two primary responsibilities:
      1. Computing critical simulation parameters (monomer counts, box dimensions,
         densities, etc.) using SystemPropertyCalculations.
      2. Generating all required LAMMPS input files for the five simulation stages
         (minimization, equilibration, production, etc.) via the Writer component.
    
    The class acts as a thin, high-level coordinator that ensures the updated
    SimulationSetup object (with calculated physical properties) is passed
    correctly to the file-writing layer.
    """

    def setup_and_write_simulation(
        self,
        setup: SimulationSetup,
        reacter_files: REACTERFiles,
        run_dir: Path
    ) -> SimulationSetup:
        """
        Execute the complete simulation setup pipeline and write all LAMMPS files.
        
        This is the main entry point called by the top-level AutoREACTER workflow.
        It performs system property calculations on all replicas and then delegates
        file generation to generate_input_files().
        
        Args:
            setup: SimulationSetup object containing user-defined parameters,
                   monomer definitions, and reaction data.
            reacter_files: REACTERFiles object containing all processed molecular
                           templates, reaction files, and data files generated
                           by the Lunar API stage.
            run_dir: Path to the final output directory where the 'LAMMPS_input_files'
                     subdirectory will be created.
        
        Returns:
            SimulationSetup: The same object with all calculated properties populated
                            (monomer counts, box lengths, densities, etc.). This
                            updated object can be used for logging or further processing.
        """
        # Calculate physical properties (monomer counts, box size, density, etc.)
        calculator = SystemPropertyCalculations(setup)
        updated_setup = calculator.process_all()

        # Optional debugging block (uncomment during development to inspect calculations)
        # for replica in updated_setup.replicas:
        #     print(f"\n[INFO] Replica: {replica.tag}")
        #     print(f"  - Calculated Monomer Counts: {replica.monomer_counts}")
        #     print(f"  - Initial Box Length: {replica.initial_box_length:.2f} Å")
        #     print(f"  - Target Density: {replica.density} g/cm³")

        # Generate all LAMMPS input files for the five simulation stages
        self.generate_input_files(
            setup=updated_setup,
            reacter_files=reacter_files,
            run_dir=run_dir
        )

        print(f"\n[SUCCESS] All 5 simulation stages written to: {run_dir / 'LAMMPS_input_files'}")
        return updated_setup

    def generate_input_files(
        self,
        setup: SimulationSetup,
        reacter_files: REACTERFiles,
        run_dir: Path
    ) -> None:
        """
        Instantiate the Writer and delegate creation of all LAMMPS input files.
        
        This method follows the single-responsibility principle by keeping file I/O
        logic inside the dedicated Writer class.
        
        Args:
            setup: SimulationSetup containing all calculated physical properties.
            reacter_files: REACTERFiles object with molecule templates and reaction data.
            run_dir: Base output directory. The Writer will create a 'LAMMPS_input_files'
                     subdirectory inside it.
        """
        writer = Writer(reacter_files=reacter_files)
        writer.write_all_files(run_dir, setup)
