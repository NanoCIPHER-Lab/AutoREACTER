import random
import shutil
from pathlib import Path
from datetime import datetime
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles
from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica

now = datetime.now().strftime("%Y-%m-%d")

class RxnSecondStageWriter:
    def __init__(self, out_dir: Path, settings: LammpsSettings, reacter_files: REACTERFiles, replica: Replica, sim_name: str):
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.reacter_files = reacter_files
        self.second_stage_file_name = self.write_second_stage_reaction_files(replica=replica)

    def write_second_stage_reaction_files(self, replica: Replica) -> str:
        tag = f"{self.sim_name}_{replica.tag}"
        rxn_dir = self.out_dir / "4_reaction_second_stage"
        rxn_dir.mkdir(parents=True, exist_ok=True)
        
        s = self.settings
        rf = self.reacter_files
        

        input_data = f"{tag}_reacted_0M-1M_3.5A.data"
        output_base = f"{tag}_reacted_1M-3.5_5.0A"
        
        lines = [
            f"# {tag} Second Reaction Stage Script - Generated {now} by AutoREACTER\n",
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
            f"{'pair_modify':<16} tail yes mix sixthpower",
            "",
            f"{'read_data':<16} \"{input_data}\" &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            f"{'minimize':<16} 1.0e-4 1.0e-6 1000 10000",
            f"{'velocity':<16} all create {replica.temperature} {random.randint(10000, 9999999)} loop geom",
            f"{'timestep':<16} 1.0",
            f"{'thermo':<16} 100",
            f"{'reset_timestep':<16} 1000000",
            ""
        ]

        lines.append("#------------Define Reaction Templates------------")
        rxn_commands = []
        for i, template in enumerate(rf.template_files, 1):
            pre_id = f"mol_pre_{i}"
            post_id = f"mol_post_{i}"
            
            pre_file = template.pre_reaction_file.lmp_molecule_file.name
            post_file = template.post_reaction_file.lmp_molecule_file.name
            map_file = template.map_file.name
            
            lines.append(f"{'molecule':<16} {pre_id} {pre_file}")
            lines.append(f"{'molecule':<16} {post_id} {post_file}")
            
            # Distance updated to 5.0 and stabilize_steps to 200 as per your template
            rxn_str = f"react rxn_stp_{i} all 1 0.0 5.0 {pre_id} {post_id} {map_file} stabilize_steps 200 rescale_charges yes"
            rxn_commands.append(rxn_str)

        all_reactions = " & \n                ".join(rxn_commands)
        
        lines.extend([
            "",
            f"{'fix':<16} rxns all bond/react stabilization yes statted_grp 0.03 &",
            f"{'':<16} {all_reactions}",
            "",
            "# Note: If atoms are being deleted during the reaction, ensure you use the correct Map file",
            "#       (e.g., RXN_i_with_delete_ids.map). ", 
            "#       NPT is recommended for deletion to account for density changes.\n",
            f"{'fix':<16} 1 statted_grp_REACT nvt temp {replica.temperature} {replica.temperature} 100.0\n",   
            f"#{'fix':<16} 1 statted_grp_REACT npt temp {replica.temperature} {replica.temperature} 100.0 iso 0.0 0.0 1000.0",
            "",
            f"{'thermo_style':<16} custom step time temp f_rxns[*] press density vol pe ke etotal",
            f"{'dump':<16} traj all xyz 1000 {output_base}.xyz",
            f"{'dump_modify':<16} traj types labels",
            "",
            f"{'run':<16} 2500000",
            f"{'restart':<16} 100 {output_base}_backup1.restart {output_base}_backup2.restart",
            f"{'write_restart':<16} {output_base}.restart",
            f"{'write_data':<16} {output_base}.data nofix"
        ])

        in_file_path = rxn_dir / f"in.{tag}_reaction_stage_2"
        with open(in_file_path, 'w') as f:
            f.write("\n".join(lines))
        
        self._copy_required_files(dest_dir=rxn_dir)

        return in_file_path.name

    def _copy_required_files(self, dest_dir: Path) -> None:
        rf = self.reacter_files
        for template in rf.template_files:
            files = [
                template.map_file,
                template.pre_reaction_file.lmp_molecule_file,
                template.post_reaction_file.lmp_molecule_file
            ]
            
            for file in files:
                if file is None or not file.exists():
                    raise FileNotFoundError(f"Required reaction file not found: {file}")
                shutil.copy2(file, dest_dir / file.name)