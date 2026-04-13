"""
Generate the LAMMPS input script for the post-equilibration stage.

This writer creates the script used after the reaction stage has completed.
It starts from the reacted data file, performs a short sequence of
temperature/pressure relaxation steps, and writes the final equilibrated
structure and restart files to the post-equilibration output directory.
"""

from pathlib import Path
from datetime import datetime

from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica


# Timestamp used in generated file headers so output scripts are traceable.
now = datetime.now().strftime("%Y-%m-%d")


class PostEqWriter:
    """
    Build the LAMMPS input file for the post-equilibration stage.

    The constructor writes the script immediately, so instantiating this class
    is enough to generate the file in the correct output directory.

    Attributes:
        settings (LammpsSettings): Shared LAMMPS simulation settings.
        out_dir (Path): Root output directory for the run.
        sim_name (str): Base simulation name used in file naming.
    """

    def __init__(self, out_dir: Path, settings: LammpsSettings, replica: Replica, sim_name: str):
        """
        Initialize the writer and generate the post-equilibration input script.

        Args:
            out_dir (Path): Root directory where stage-specific folders are created.
            settings (LammpsSettings): LAMMPS force-field and run parameters.
            replica (Replica): Replica-specific settings such as tag and temperature.
            sim_name (str): Simulation name prefix used to build output filenames.
        """
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name

        # Write the post-equilibration file immediately during construction.
        self.write_post_eq_file(replica=replica)

    def write_post_eq_file(self, replica: Replica) -> str:
        """
        Create the LAMMPS input script for the post-equilibration stage.

        The script:
        1. Reads the reacted data from the previous stage.
        2. Minimizes the structure.
        3. Runs NVT to move the system toward 300 K.
        4. Runs NPT to relax pressure and density.
        5. Finishes with a final NVT stabilization.
        6. Writes the final data and restart files.

        Args:
            replica (Replica): Replica configuration, including the run tag and temperature.

        Returns:
            str: The filename of the generated LAMMPS input script.
        """
        tag = f"{self.sim_name}_{replica.tag}"

        # Store all post-equilibration artifacts in a dedicated stage folder.
        post_dir = self.out_dir / "5_post_equilibration"
        post_dir.mkdir(parents=True, exist_ok=True)

        s = self.settings

        # Input is the reacted structure produced by the second reaction stage.
        input_data = f"{tag}_reacted_1M-3.5_5.0A.data"

        # Output files written after equilibration completes.
        output_xyz = f"{tag}_post_equilibration.xyz"
        output_data = f"{tag}_post_equilibrated.data"

        # Restart files are written in rotation so the simulation can be resumed safely.
        restart_1 = f"{tag}_post_equilibration_backup1.restart"
        restart_2 = f"{tag}_post_equilibration_backup2.restart"

        # Base LAMMPS setup: units, topology style, and force-field definitions.
        lines = [
            f"# Script for {tag}:",
            "# Post-reaction equilibration",
            f"# Generated {now}\n",
            f"{'units':<16} {s.units}",
            f"{'dimension':<16} {s.dimension}",
            f"{'boundary':<16} {s.boundary}",
            f"{'atom_style':<16} {s.atom_style}",
            "",
            # Bonded interaction styles must match the force field used upstream.
            f"{'angle_style':<16} {s.angle_style}",
            f"{'bond_style':<16} {s.bond_style}",
            f"{'dihedral_style':<16} {s.dihedral_style}",
            f"{'improper_style':<16} {s.improper_style}",
            "",
            # Pair and long-range interaction settings are inherited from the simulation config.
            f"{'pair_style':<16} {s.pair_style}",
            f"{'kspace_style':<16} {s.kspace_style}",
            f"{'pair_modify':<16} {s.pair_modify}",
        ]

        # Optional neighbor settings are only included when configured.
        if s.neighbor:
            lines.append(f"{'neighbor':<16} {s.neighbor}")
        if s.neigh_modify:
            lines.append(f"{'neigh_modify':<16} {s.neigh_modify}")

        lines.extend([
            "",
            # Read the reacted data file and reserve extra bonded capacity for future reactions.
            f"{'read_data':<16} {input_data} &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            # Light minimization removes any bad contacts before dynamics begins.
            f"{'minimize':<16} 0.0 1.0e-8 1000 100000",
            f"{'reset_timestep':<16} 0",
            f"{'timestep':<16} 1",
            f"{'thermo':<16} 1000",
            f"{'thermo_style':<16} custom step time temp press density vol pe ke etotal",
            f"{'dump':<16} dump_1 all xyz 100 {output_xyz}",
            f"{'dump_modify':<16} dump_1 types labels",
            # Create rolling restart files so the run can be resumed if needed.
            f"{'restart':<16} 100 {restart_1} {restart_2}",
            "",
            # Stage 1: bring the system from the replica temperature to 300 K under NVT.
            f"{'fix':<16} nvt_1 all nvt temp {replica.temperature} 300.0 100.0",
            f"{'run':<16} 100000",
            f"{'unfix':<16} nvt_1",
            "",
            # Stage 2: relax pressure and density at 300 K using NPT.
            f"{'fix':<16} npt_2 all npt temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0",
            f"{'run':<16} 100000",
            f"{'unfix':<16} npt_2",
            "",
            # Stage 3: final thermal stabilization at 300 K under NVT.
            f"{'fix':<16} nvt_3 all nvt temp 300.0 300.0 100.0",
            f"{'run':<16} 100000",
            f"{'unfix':<16} nvt_3",
            "",
            # Stop collecting trajectory frames before saving the final structure.
            f"{'undump':<16} dump_1",
            f"{'write_data':<16} {output_data}",
        ])

        # Write the generated LAMMPS input script to the post-equilibration folder.
        in_file_path = post_dir / f"in.{tag}_post_equilibration"
        with open(in_file_path, "w") as f:
            f.write("\n".join(lines))

        return in_file_path.name
