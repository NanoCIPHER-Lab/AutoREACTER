"""
Module for generating LAMMPS input scripts for second-stage molecular dynamics reaction simulations.

This module handles the creation of LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator)
input files for the second reaction stage of AutoREACTER workflows, including reaction template setup,
simulation parameters, and required file management.
"""

import random
import shutil
from pathlib import Path
from datetime import datetime
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles
from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica

# Generate current date string for file headers and timestamping
now = datetime.now().strftime("%Y-%m-%d")


class RxnSecondStageWriter:
    """
    Generates LAMMPS input scripts for second-stage reaction simulations.
    
    This class creates a complete LAMMPS input script for the second reaction stage,
    configures inter-molecular reaction parameters, and manages the copying of required
    reaction template files to the output directory.
    
    Attributes:
        settings (LammpsSettings): LAMMPS simulation settings and parameters
        out_dir (Path): Output directory for generated files
        sim_name (str): Base name for the simulation
        reacter_files (REACTERFiles): Container for reaction template files
        second_stage_file_name (str): Name of the generated second-stage input file
    """

    def __init__(self, out_dir: Path, settings: LammpsSettings, reacter_files: REACTERFiles, replica: Replica, sim_name: str):
        """
        Initialize the RxnSecondStageWriter and generate second-stage reaction files.
        
        Args:
            out_dir (Path): Root output directory for simulation files
            settings (LammpsSettings): LAMMPS simulation settings and force field parameters
            reacter_files (REACTERFiles): Container with reaction template and molecule files
            replica (Replica): Replica configuration with temperature and tag information
            sim_name (str): Base simulation name used for file naming
        """
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.reacter_files = reacter_files
        # Generate and store the second-stage input file name
        self.second_stage_file_name = self.write_second_stage_reaction_files(replica=replica)

    def write_second_stage_reaction_files(self, replica: Replica) -> str:
        """
        Generate LAMMPS input script for second-stage inter-molecular reactions.
        
        Creates a complete LAMMPS input file with reaction definitions, simulation
        parameters, and thermal dynamics settings. This stage simulates reactions
        between molecules at medium intermolecular distances (1M-3.5 Angstroms).
        
        Args:
            replica (Replica): Replica configuration containing temperature and unique tag
            
        Returns:
            str: Filename of the generated second-stage input script
            
        Raises:
            FileNotFoundError: If required reaction template files are missing
        """
        tag = f"{self.sim_name}_{replica.tag}"
        rxn_dir = self.out_dir / "4_reaction_second_stage"
        rxn_dir.mkdir(parents=True, exist_ok=True)
        
        # Reference commonly used objects for cleaner code
        s = self.settings
        rf = self.reacter_files

        # Define input and output data file names
        input_data = f"{tag}_reacted_0M-1M_3.5A.data"
        output_base = f"{tag}_reacted_1M-3.5_5.0A"
        
        # Build LAMMPS script header and core simulation parameters
        lines = [
            f"# {tag} Second Reaction Stage Script - Generated {now} by AutoREACTER\n",
            "#------------Initialization------------",
            f"{'units':<16} {s.units}",
            f"{'dimension':<16} {s.dimension}",
            f"{'boundary':<16} {s.boundary}",
            f"{'atom_style':<16} {s.atom_style}",
            "",
            # Force field style definitions
            "#------------Force Field Styles------------",
            f"{'angle_style':<16} {s.angle_style}",
            f"{'bond_style':<16} {s.bond_style}",
            f"{'dihedral_style':<16} {s.dihedral_style}",
            f"{'improper_style':<16} {s.improper_style}",
            "",
            # Pair and electrostatic interaction styles
            f"{'pair_style':<16} {s.pair_style}",
            f"{'kspace_style':<16} {s.kspace_style}",
            f"{'pair_modify':<16} {s.pair_modify}"
        ]
        
        # Add optional neighbor list settings if specified
        if s.neighbor:
            lines.append(f"{'neighbor':<16} {s.neighbor}")
        if s.neigh_modify:
            lines.append(f"{'neigh_modify':<16} {s.neigh_modify}")

        # Define data file reading with extra bonding capacity for reactions
        lines.extend([
            "",
            "#------------Read First Stage reacted Box------------",
            f"{'read_data':<16} \"{input_data}\" &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            # Initial energy minimization and simulation setup
            "#------------Minimization and Velocity Initialization------------",
            f"{'minimize':<16} 1.0e-4 1.0e-6 1000 10000",
            f"{'velocity':<16} all create {replica.temperature} {random.randint(10000, 9999999)} loop geom",
            f"{'timestep':<16} 1.0",
            f"{'thermo':<16} 100",
            f"{'reset_timestep':<16} 1000000",
            ""
        ])

        # Define reaction templates and build reaction fix commands
        lines.append("#------------Define Reaction Templates------------")
        rxn_commands = []
        for i, template in enumerate(rf.template_files, 1):
            # Create unique identifiers for pre- and post-reaction molecule templates
            pre_id = f"mol_pre_{i}"
            post_id = f"mol_post_{i}"
            
            # Extract molecule and mapping file names from template objects
            pre_file = template.pre_reaction_file.lmp_molecule_file.name
            post_file = template.post_reaction_file.lmp_molecule_file.name
            map_file = template.map_file.name
            
            # Register molecule templates in LAMMPS
            lines.append(f"{'molecule':<16} {pre_id} {pre_file}")
            lines.append(f"{'molecule':<16} {post_id} {post_file}")
            
            # Build reaction command with cutoff at 5.0 Angstroms for inter-molecular reactions
            # stabilize_steps: equilibration steps after reaction to prevent explosion
            rxn_str = f"react rxn_stp_{i} all 1 0.0 5.0 {pre_id} {post_id} {map_file} stabilize_steps 200 rescale_charges yes"
            rxn_commands.append(rxn_str)

        # Combine all reaction commands with line continuation
        all_reactions = " & \n                ".join(rxn_commands)
        
        # Configure bond/react fix with stabilization groups and define fixes for simulation
        lines.extend([
            "",
            f"{'fix':<16} rxns all bond/react stabilization yes statted_grp 0.03 &",
            f"{'':<16} {all_reactions}",
            "",
            "# Note: If atoms are being deleted during the reaction, ensure you use the correct Map file",
            "#       (e.g., RXN_i_with_delete_ids.map). ", 
            "#       NPT is recommended for deletion to account for density changes.\n",
            # NVT thermostat for stabilized group (uncomment NPT if density changes expected)
            f"{'fix':<16} 1 statted_grp_REACT nvt temp {replica.temperature} {replica.temperature} 100.0\n",   
            f"#{'fix':<16} 1 statted_grp_REACT npt temp {replica.temperature} {replica.temperature} 100.0 iso 0.0 0.0 1000.0",
            "",
            # Output configuration
            f"{'thermo_style':<16} custom step time temp f_rxns[*] press density vol pe ke etotal",
            f"{'dump':<16} traj all xyz 1000 {output_base}.xyz",
            f"{'dump_modify':<16} traj types labels",
            "",
            # Run simulation for 2.5 million timesteps and save periodic backups
            f"{'run':<16} 2500000",
            f"{'restart':<16} 100 {output_base}_backup1.restart {output_base}_backup2.restart",
            f"{'write_restart':<16} {output_base}.restart",
            f"{'write_data':<16} {output_base}.data nofix"
        ])

        # Write assembled LAMMPS input script to file
        in_file_path = rxn_dir / f"in.{tag}_reaction_stage_2"
        with open(in_file_path, 'w') as f:
            f.write("\n".join(lines))
        
        # Copy all required reaction template files to output directory
        self._copy_required_files(dest_dir=rxn_dir)

        return in_file_path.name

    def _copy_required_files(self, dest_dir: Path) -> None:
        """
        Copy reaction template and molecule definition files to the output directory.
        
        Copies pre-reaction molecules, post-reaction molecules, and atom mapping files
        required by LAMMPS to the simulation directory for execution.
        
        Args:
            dest_dir (Path): Destination directory for copied files
            
        Raises:
            FileNotFoundError: If any required reaction file is missing or cannot be accessed
        """
        rf = self.reacter_files
        for template in rf.template_files:
            # Collect all files associated with this reaction template
            files = [
                template.map_file,
                template.pre_reaction_file.lmp_molecule_file,
                template.post_reaction_file.lmp_molecule_file
            ]
            
            # Copy each file to the destination, preserving metadata
            for file in files:
                if file is None or not file.exists():
                    raise FileNotFoundError(f"Required reaction file not found: {file}")
                shutil.copy2(file, dest_dir / file.name)
