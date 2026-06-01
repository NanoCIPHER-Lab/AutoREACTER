import random
import shutil
from pathlib import Path
from datetime import datetime
from AutoREACTER.reaction_preparation.lunar_client.REACTER_files_builder import REACTERFiles
from AutoREACTER.sim_setup.writers.lammps_settings import LammpsSettings
from AutoREACTER.input_parser import Replica

now = datetime.now().strftime("%Y-%m-%d")


class RxnFirstStageWriter:
    """Generates LAMMPS input scripts and copies required files for the first reaction stage.

    This class constructs a complete LAMMPS input file (.in) that sets up a bond/react
    simulation, including initialization, force-field styles, molecule template reads,
    reaction fix definitions, thermostat/barostat settings, output dumps, and run/restart
    commands. All referenced template, map, and molecule files are also copied into the
    reaction output directory.

    Parameters
    ----------
    out_dir : Path
        Root output directory under which ``3_reaction_first_stage`` will be created.
    settings : LammpsSettings
        LAMMPS simulation settings (units, atom_style, pair_style, kspace_style, etc.).
    reacter_files : REACTERFiles
        Container for reaction template files (pre/post molecule files, map files).
    replica : Replica
        Replica metadata including temperature and a human-readable tag.
    sim_name : str
        Base simulation name used to prefix output file names.
    """

    def __init__(
        self,
        out_dir: Path,
        settings: LammpsSettings,
        reacter_files: REACTERFiles,
        replica: Replica,
        sim_name: str,
    ):
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.reacter_files = reacter_files
        self.first_stage_file_name = self.write_first_stage_reaction_files(replica=replica)


    #  Public API
    def write_first_stage_reaction_files(self, replica: Replica) -> str:
        """Build and write the first-stage LAMMPS reaction input file.

        The generated input file follows this general structure:

        1. **Initialization** – units, dimension, boundary, atom_style.
        2. **Force-field styles** – angle/bond/dihedral/improper/pair/kspace styles.
        3. **Read equilibrated box** – ``read_data`` on the pre-equilibrated
           structure with extra bond/angle/dihedral/improper/special slots.
        4. **Minimization & velocity** – energy minimization followed by
           Gaussian velocity initialization at *replica.temperature*.
        5. **Reaction templates** – ``molecule`` commands for every pre/post
           template pair plus the corresponding ``fix bond/react`` commands
           joined with LAMMPS line-continuation (``&``).
        6. **Thermostat / barostat** – NVT (or commented-out NPT) on the
           ``statted_grp_REACT`` group.
        7. **Output & run** – thermo style, XYZ trajectory dump, a 1 000 000-step
           run, restart files, and a final ``write_data`` call.

        All auxiliary files (map files, pre/post molecule files) are copied
        into the reaction directory via :meth:`_copy_required_files`.

        Parameters
        ----------
        replica : Replica
            Replica carrying the temperature and tag for this stage.

        Returns
        -------
        str
            The **file name** (not path) of the generated ``in.*`` script.
        """
        tag = f"{self.sim_name}_{replica.tag}"
        rxn_dir = self.out_dir / "3_reaction_first_stage"
        rxn_dir.mkdir(parents=True, exist_ok=True)

        s = self.settings
        rf = self.reacter_files
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        input_data = f"{tag}_pre_equilibrated.data"
        output_base = f"{tag}_reacted_0M-1M_3.5A"

        # ----- Header / Initialization ---------------------------------
        lines = [
            f"# {tag} First Reaction Stage Script - Generated {now} by AutoREACTER\n",
            "#------------Initialization------------",
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

        # ----- Read equilibrated structure ------------------------------
        lines.extend([
            "",
            "#------------Read Equilibrated Box------------",
            f"{'read_data':<16} {input_data} &",
            f"{'':<16} extra/bond/per/atom 50 &",
            f"{'':<16} extra/angle/per/atom 50 &",
            f"{'':<16} extra/dihedral/per/atom 50 &",
            f"{'':<16} extra/improper/per/atom 50 &",
            f"{'':<16} extra/special/per/atom 50",
            "",
            "#------------Minimization and Velocity------------",
            f"{'minimize':<16} 1.0e-4 1.0e-6 1000 10000",
            f"{'velocity':<16} all create {replica.temperature} {random.randint(10_000, 9_999_999)} dist gaussian",
            f"{'timestep':<16} 1.0",
            f"{'thermo':<16} 100",
            f"{'reset_timestep':<16} 0",
            "",
        ])

        # ----- Reaction templates & fix bond/react ---------------------
        lines.append("#------------Define Reaction Templates------------")
        rxn_commands: list[str] = []

        for i, template in enumerate(rf.template_files, 1):
            pre_id = f"mol_pre_{i}"
            post_id = f"mol_post_{i}"

            # Extract filenames from the dataclass fields
            pre_file = template.pre_reaction_file.lmp_molecule_file.name
            post_file = template.post_reaction_file.lmp_molecule_file.name
            map_file = template.map_file.name

            lines.append(f"{'molecule':<16} {pre_id} {pre_file}")
            lines.append(f"{'molecule':<16} {post_id} {post_file}\n")

            rxn_str = (
                f"react rxn_stp_{i} all 1 0.0 3.5 {pre_id} {post_id} {map_file} "
                f"stabilize_steps 60 rescale_charges yes"
            )
            rxn_commands.append(rxn_str)

        all_reactions = " & \n                ".join(rxn_commands)

        lines.extend([
            "",
            f"{'fix':<16} rxns all bond/react stabilization yes statted_grp 0.03 &",
            f"{'':<16} {all_reactions}",
            "",
            "",
            "# Note: If atoms are being deleted during the reaction, ensure you use the correct Map file",
            "#       (e.g., RXN_i_with_delete_ids.map).",
            "        NPT is recommended for deletion to account for density changes.\n",
            f"{'fix':<16} 1 statted_grp_REACT nvt temp {replica.temperature} {replica.temperature} 100.0\n",
            f"#{'fix':<16} 1 statted_grp_REACT npt temp {replica.temperature} {replica.temperature} 100.0 iso 0.0 0.0 1000.0",
            "",
            f"{'thermo_style':<16} custom step time temp f_rxns[*] press density vol pe ke etotal",
            f"{'dump':<16} traj all xyz 1000 {output_base}.xyz",
            f"{'dump_modify':<16} traj types labels",
            "",
            f"{'run':<16} 1000000",
            f"{'restart':<16} 100 {output_base}_backup1.restart {output_base}_backup2.restart",
            f"{'write_restart':<16} {output_base}.restart",
            f"{'write_data':<16} {output_base}.data nofix",
        ])

        # ----- Write the .in file --------------------------------------
        in_file_path = rxn_dir / f"in.{tag}_reaction"
        with open(in_file_path, "w") as f:
            f.write("\n".join(lines))

        # ----- Copy auxiliary files into the reaction directory ---------
        self._copy_required_files(dest_dir=rxn_dir)

        return in_file_path.name

    #  Internal helpers
    def _copy_required_files(self, dest_dir: Path) -> None:
        """Copy every map and molecule file referenced by the reaction templates.

        Parameters
        ----------
        dest_dir : Path
            Target directory (``3_reaction_first_stage``).

        Raises
        ------
        FileNotFoundError
            If any referenced file does not exist on disk.
        """
        rf = self.reacter_files

        for template in rf.template_files:
            files: list[Path] = [
                template.map_file,
                template.pre_reaction_file.lmp_molecule_file,
                template.post_reaction_file.lmp_molecule_file,
            ]

            for file in files:
                if file is None or not file.exists():
                    raise FileNotFoundError(f"Required reaction file not found: {file}")
                shutil.copy2(file, dest_dir / file.name)
