from pathlib import Path
from datetime import datetime
from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica

now = datetime.now().strftime("%Y-%m-%d")

class PreEqWriter:
    def __init__(self, out_dir: Path, settings: LammpsSettings, replica: Replica, sim_name: str):
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.in_pre_eq_file_name = self.write_pre_eq_file(replica=replica)

    def write_pre_eq_file(self, replica: Replica) -> str:
        tag = f"{self.sim_name}_{replica.tag}"
        dens_dir = self.out_dir / "2_pre_equilibration"
        dens_dir.mkdir(parents=True, exist_ok=True)
        
        s = self.settings
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        input_data = f"{tag}_shrinked_box.data"
        output_xyz = f"{tag}_pre_equilibration.xyz"
        output_data = f"{tag}_pre_equilibrated.data"
        restart_1 = f"{tag}_pre_equilibration_backup1.restart"
        restart_2 = f"{tag}_pre_equilibration_backup2.restart"

        lines = [
            f"# {tag} Pre-Equilibration Script - Generated {now} by AutoREACTER\n",
            f"#------------Initialization------------",
            f"{'units':<16} {s.units}",
            f"{'dimension':<16} {s.dimension}",
            f"{'boundary':<16} {s.boundary}",
            f"{'atom_style':<16} {s.atom_style}",
            "",
            "# ------------Force Field Styles------------",
            f"{'angle_style':<16} {s.angle_style}",
            f"{'bond_style':<16} {s.bond_style}",
            f"{'dihedral_style':<16} {s.dihedral_style}",
            f"{'improper_style':<16} {s.improper_style}",
            "",
            f"{'pair_style':<16} {s.pair_style}",
            f"{'kspace_style':<16} {s.kspace_style}",
            f"{'pair_modify':<16} {s.pair_modify}",
        ]

        if s.neighbor:
            lines.append(f"{'neighbor':<16} {s.neighbor}")
        if s.neigh_modify:
            lines.append(f"{'neigh_modify':<16} {s.neigh_modify}")

        lines.extend([
            "",
            "#------------Read densified box------------",
            f"{'read_data':<16} {input_data} &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            "#------------Simulation Setup------------",
            f"{'minimize':<16} 0.0 1.0e-8 1000 100000",
            f"{'reset_timestep':<16} 0",
            f"{'timestep':<16} 1",
            f"{'thermo':<16} 1000",
            f"{'thermo_style':<16} custom step time temp press density vol pe ke etotal",
            f"{'dump':<16} dump_2 all xyz 1000 {output_xyz}",
            f"{'dump_modify':<16} dump_2 types labels",
            f"{'restart':<16} 100 {restart_1} {restart_2}",
            "",
            "# First NVT ramp/hold",
            f"{'fix':<16} nvt_1 all nvt temp 298.15 {replica.temperature} 100.0",
            f"{'run':<16} 25000",
            f"{'unfix':<16} nvt_1",
            "",
            "# NPT equilibration",
            f"{'fix':<16} npt_2 all npt temp {replica.temperature} {replica.temperature} 100.0 iso 0.0 0.0 1000.0",
            f"{'run':<16} 25000",
            f"{'unfix':<16} npt_2",
            "",
            "# Final NVT equilibration",
            f"{'fix':<16} nvt_3 all nvt temp {replica.temperature} {replica.temperature} 100.0",
            f"{'run':<16} 50000",
            f"{'unfix':<16} nvt_3",
            "",
            f"{'undump':<16} dump_2",
            f"{'write_data':<16} {output_data}"
        ])

        in_dest_location = dens_dir / f"in.{tag}_pre_equilibration"
        with open(in_dest_location, 'w') as f:
            f.write("\n".join(lines))
    
        return in_dest_location.name

