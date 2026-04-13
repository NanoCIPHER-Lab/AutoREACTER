from pathlib import Path
from datetime import datetime
from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica

now = datetime.now().strftime("%Y-%m-%d")


class PostEqWriter:
    def __init__(self, out_dir: Path, settings: LammpsSettings, replica: Replica, sim_name: str):
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.write_post_eq_file(replica=replica)

    def write_post_eq_file(self, replica: Replica) -> str:
        tag = f"{self.sim_name}_{replica.tag}"
        post_dir = self.out_dir / "5_post_equilibration"
        post_dir.mkdir(parents=True, exist_ok=True)
        
        s = self.settings

        input_data = f"{tag}_reacted_1M-3.5_5.0A.data"
        output_xyz = f"{tag}_post_equilibration.xyz"
        output_data = f"{tag}_post_equilibrated.data"
        restart_1 = f"{tag}_post_equilibration_backup1.restart"
        restart_2 = f"{tag}_post_equilibration_backup2.restart"

        lines = [
            f"# Script for {tag}:",
            "# Post-reaction equilibration",
            f"# Generated {now}\n",
            f"{'units':<16} {s.units}",
            f"{'dimension':<16} {s.dimension}",
            f"{'boundary':<16} {s.boundary}",
            f"{'atom_style':<16} {s.atom_style}",
            "",
            f"{'angle_style':<16} {s.angle_style}",
            f"{'bond_style':<16} {s.bond_style}",
            f"{'dihedral_style':<16} {s.dihedral_style}",
            f"{'improper_style':<16} {s.improper_style}",
            "",
            f"{'pair_style':<16} {s.pair_style}",
            f"{'kspace_style':<16} {s.kspace_style}",
            f"{'pair_modify':<16} {s.pair_modify}"
        ]

        if s.neighbor:
            lines.append(f"{'neighbor':<16} {s.neighbor}")
        if s.neigh_modify:
            lines.append(f"{'neigh_modify':<16} {s.neigh_modify}")

        lines.extend([
            "",
            f"{'read_data':<16} {input_data} &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            f"{'minimize':<16} 0.0 1.0e-8 1000 100000",
            f"{'reset_timestep':<16} 0",
            f"{'timestep':<16} 1",
            f"{'thermo':<16} 1000",
            f"{'thermo_style':<16} custom step time temp press density vol pe ke etotal",
            f"{'dump':<16} dump_1 all xyz 100 {output_xyz}",
            f"{'dump_modify':<16} dump_1 types labels",
            f"{'restart':<16} 100 {restart_1} {restart_2}",
            "",
            f"{'fix':<16} nvt_1 all nvt temp {replica.temperature} 300.0 100.0",
            f"{'run':<16} 100000",
            f"{'unfix':<16} nvt_1",
            "",
            f"{'fix':<16} npt_2 all npt temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0",
            f"{'run':<16} 100000",
            f"{'unfix':<16} npt_2",
            "",
            f"{'fix':<16} nvt_3 all nvt temp 300.0 300.0 100.0",
            f"{'run':<16} 100000",
            f"{'unfix':<16} nvt_3",
            "",
            f"{'undump':<16} dump_1",
            f"{'write_data':<16} {output_data}"
        ])

        in_file_path = post_dir / f"in.{tag}_post_equilibration"
        with open(in_file_path, 'w') as f:
            f.write("\n".join(lines))
        
        return in_file_path.name