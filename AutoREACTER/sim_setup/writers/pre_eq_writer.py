from pathlib import Path
from datetime import datetime
from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica


class PreEqWriter:
    """
    Generates a LAMMPS input script for the pre-equilibration stage of a polymer
    simulation. This stage follows densification and consists of:
    
    1. Energy minimization
    2. NVT temperature ramp from 298.15 K to target temperature
    3. NPT equilibration at target temperature and 0 atm pressure
    4. Final NVT equilibration at target temperature
    
    The script reads the shrinked (densified) box produced by the previous stage.
    """

    def __init__(
        self,
        out_dir: Path,
        settings: LammpsSettings,
        replica: Replica,
        sim_name: str,
    ):
        """
        Initialize the pre-equilibration writer and immediately generate the LAMMPS script.

        Args:
            out_dir: Base output directory for simulation files
            settings: LAMMPS settings configuration (units, styles, pair potentials, etc.)
            replica: Contains replica-specific parameters such as target temperature and tag
            sim_name: Name prefix used for all output files and the input script
        """
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.in_pre_eq_file_name = self.write_pre_eq_file(replica=replica)

    def write_pre_eq_file(self, replica: Replica) -> str:
        """
        Creates the complete LAMMPS input script for pre-equilibration and writes
        it to the '2_pre_equilibration' subdirectory.

        The generated script performs a multi-stage equilibration protocol:
        - Minimization of the densified configuration
        - NVT ramp from room temperature (298.15 K) to the target temperature
        - NPT equilibration at target temperature and zero pressure (isotropic)
        - Extended NVT equilibration at the final temperature

        Args:
            replica: Replica object containing temperature and simulation tag

        Returns:
            str: The filename of the generated LAMMPS input script
        """
        tag = f"{self.sim_name}_{replica.tag}"
        pre_eq_dir = self.out_dir / "2_pre_equilibration"
        pre_eq_dir.mkdir(parents=True, exist_ok=True)

        s = self.settings
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Filenames used throughout the LAMMPS script
        input_data = f"{tag}_shrinked_box.data"
        output_xyz = f"{tag}_pre_equilibration.xyz"
        output_data = f"{tag}_pre_equilibrated.data"
        restart_1 = f"{tag}_pre_equilibration_backup1.restart"
        restart_2 = f"{tag}_pre_equilibration_backup2.restart"

        lines = [
            f"# {tag} Pre-Equilibration Script - Generated {now} by AutoREACTER",
            "",
            "#------------Initialization------------",
            f"{'units':<16} {s.units}",
            f"{'dimension':<16} {s.dimension}",
            f"{'boundary':<16} {s.boundary}",
            f"{'atom_style':<16} {s.atom_style}",
            "",
            "#------------Force Field Styles------------",
            f"{'bond_style':<16} {s.bond_style}",
            f"{'angle_style':<16} {s.angle_style}",
            f"{'dihedral_style':<16} {s.dihedral_style}",
            f"{'improper_style':<16} {s.improper_style}",
            "",
            f"{'pair_style':<16} {s.pair_style}",
            f"{'kspace_style':<16} {s.kspace_style}",
            f"{'pair_modify':<16} {s.pair_modify}",
        ]

        # Add optional neighbor list settings if they were configured
        if s.neighbor:
            lines.append(f"{'neighbor':<16} {s.neighbor}")
        if s.neigh_modify:
            lines.append(f"{'neigh_modify':<16} {s.neigh_modify}")

        lines.extend([
            "",
            "#------------Read Densified Configuration------------",
            "# Reads the final configuration from the densification stage.",
            f"{'read_data':<16} {input_data} &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            "#------------Simulation Setup------------",
            "# Perform an initial energy minimization to remove bad contacts",
            f"{'minimize':<16} 0.0 1.0e-8 1000 100000",
            f"{'reset_timestep':<16} 0",
            "",
            "# Use a 1 fs timestep for all stages of equilibration",
            f"{'timestep':<16} 1",
            f"{'thermo':<16} 1000",
            f"{'thermo_style':<16} custom step time temp press density vol pe ke etotal",
            "",
            "# Trajectory output (XYZ format with atom type labels)",
            f"{'dump':<16} dump_2 all xyz 1000 {output_xyz}",
            f"{'dump_modify':<16} dump_2 types labels",
            "",
            "# Automatic restart files every 100 steps",
            f"{'restart':<16} 100 {restart_1} {restart_2}",
            "",
            "#------------Stage 1: NVT Temperature Ramp------------",
            "# (25,000 steps × 1 fs timestep)",
            f"{'fix':<16} nvt_1 all nvt temp 298.15 {replica.temperature} 100.0",
            f"{'run':<16} 25000",
            f"{'unfix':<16} nvt_1",
            "",
            "#------------Stage 2: NPT Equilibration------------",
            "# Isotropic pressure control at 0 atm while maintaining target temperature.",
            "# This allows the box volume to adjust to the correct density at the",
            "# desired temperature and pressure.",
            f"{'fix':<16} npt_2 all npt temp {replica.temperature} {replica.temperature} 100.0 iso 0.0 0.0 1000.0",
            f"{'run':<16} 25000",
            f"{'unfix':<16} npt_2",
            "",
            "#------------Stage 3: Final NVT Equilibration------------",
            "# Extended constant-volume equilibration at the final temperature.",
            "# This stabilizes the system after the density adjustment in the NPT stage.",
            f"{'fix':<16} nvt_3 all nvt temp {replica.temperature} {replica.temperature} 100.0",
            f"{'run':<16} 50000",
            f"{'unfix':<16} nvt_3",
            "",
            "#------------Cleanup and Output------------",
            f"{'undump':<16} dump_2",
            f"{'write_data':<16} {output_data}"
        ])

        in_dest_location = pre_eq_dir / f"in.{tag}_pre_equilibration"
        with open(in_dest_location, 'w') as f:
            f.write("\n".join(lines))

        return in_dest_location.name
