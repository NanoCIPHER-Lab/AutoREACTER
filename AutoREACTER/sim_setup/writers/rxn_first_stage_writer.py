from __future__ import annotations

import random
import shutil

from pathlib import Path
from datetime import datetime
from typing import TYPE_CHECKING

from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import (
    REACTERFiles,
)
from AutoREACTER.sim_setup.writers.lammps_settings import (
    LammpsSettings,
)
from AutoREACTER.input_parser import Simulation


if TYPE_CHECKING:
    from AutoREACTER.kinetics.estimator import KineticsEstimator


now = datetime.now().strftime("%Y-%m-%d")


class RxnFirstStageWriter:
    """Generates LAMMPS input scripts and copies required files for the first reaction stage.

    This class constructs a complete LAMMPS input file (.in) that sets up a bond/react
    simulation, including initialization, force-field styles, molecule template reads,
    reaction fix definitions, thermostat/barostat settings, output dumps, and run/restart
    commands. All referenced template, map, and molecule files are also copied into the
    reaction output directory.

    When kinetics estimation is enabled, RMG kinetics are estimated for each active
    reaction at the current simulation conditions and written directly into the copied
    LAMMPS reaction map file. No reaction-by-simulation kinetics matrix is stored.

    Parameters
    ----------
    out_dir : Path
        Root output directory under which ``3_reaction_first_stage`` will be created.
    settings : LammpsSettings
        LAMMPS simulation settings (units, atom_style, pair_style, kspace_style, etc.).
    reacter_files : REACTERFiles
        Container for reaction template files (pre/post molecule files, map files).
    simulation : Simulation
        Simulation metadata including temperature and a human-readable tag.
    sim_name : str
        Base simulation name used to prefix output file names.
    kinetics_estimator : KineticsEstimator | None
        Optional initialized RMG kinetics estimator. When provided, kinetics are
        estimated and written directly into each copied active reaction map file.
    """

    def __init__(
        self,
        out_dir: Path,
        settings: LammpsSettings,
        reacter_files: REACTERFiles,
        simulation: Simulation,
        sim_name: str,
        kinetics_estimator: "KineticsEstimator | None" = None,
    ):
        self.settings = settings
        self.out_dir = out_dir
        self.sim_name = sim_name
        self.reacter_files = reacter_files
        self.kinetics_estimator = kinetics_estimator

        self.first_stage_file_name = (
            self.write_first_stage_reaction_files(
                simulation=simulation
            )
        )

    #  Public API
    def write_first_stage_reaction_files(
        self,
        simulation: Simulation,
    ) -> str:
        """Build and write the first-stage LAMMPS reaction input file.

        The generated input file follows this general structure:

        1. **Initialization** – units, dimension, boundary, atom_style.
        2. **Force-field styles** – angle/bond/dihedral/improper/pair/kspace styles.
        3. **Read equilibrated box** – ``read_data`` on the pre-equilibrated
           structure with extra bond/angle/dihedral/improper/special slots.
        4. **Minimization & velocity** – energy minimization followed by
           Gaussian velocity initialization at *simulation.temperature*.
        5. **Reaction templates** – ``molecule`` commands for every pre/post
           template pair plus the corresponding ``fix bond/react`` commands
           joined with LAMMPS line-continuation (``&``).
        6. **Thermostat / barostat** – NVT (or commented-out NPT) on the
           ``statted_grp_REACT`` group.
        7. **Output & run** – thermo style, XYZ trajectory dump, a 1 000 000-step
           run, restart files, and a final ``write_data`` call.

        All auxiliary files (map files, pre/post molecule files) are copied
        into the reaction directory via :meth:`_copy_required_files`.

        If a kinetics estimator is provided, kinetics are estimated after the
        reaction files are copied, and the Arrhenius constraint is written only
        to the copied active map file.

        Parameters
        ----------
        simulation : Simulation
            Simulation object containing the temperature and tag for this stage.

        Returns
        -------
        str
            The **file name** (not path) of the generated ``in.*`` script.
        """
        tag = (
            f"{self.sim_name}_{simulation.tag}"
        )

        rxn_dir = (
            self.out_dir
            / "3_reaction_first_stage"
        )

        rxn_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        s = self.settings
        rf = self.reacter_files

        now = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

        input_data = (
            f"{tag}_pre_equilibrated.data"
        )

        output_base = (
            f"{tag}_reacted_0M-1M_3.5A"
        )

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
            lines.append(
                f"{'neighbor':<16} {s.neighbor}"
            )

        if s.neigh_modify:
            lines.append(
                f"{'neigh_modify':<16} {s.neigh_modify}"
            )

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
            (
                f"{'velocity':<16} "
                f"all create {simulation.temperature} "
                f"{random.randint(10_000, 9_999_999)} "
                f"dist gaussian"
            ),
            f"{'timestep':<16} 1.0",
            f"{'thermo':<16} 100",
            f"{'reset_timestep':<16} 0",
            "",
        ])

        # ----- Reaction templates & fix bond/react ---------------------
        lines.append(
            "#------------Define Reaction Templates------------"
        )

        rxn_commands: list[str] = []

        for template in [
            t for t in rf.template_files
            if getattr(t, "activity_stats", True)
        ]:

            # Extract filenames from the dataclass fields
            pre_file = (
                template.pre_reaction_file.name
            )

            post_file = (
                template.post_reaction_file.name
            )

            # Use the active map generated by AutoREACTER. For reactions with
            # deleted atoms, this map includes the required DeleteIDs section.
            map_file = (
                template.map_file.name
            )

            reaction_id = (
                template.reaction_id
            )

            pre_id = (
                f"mol_pre_{reaction_id}"
            )

            post_id = (
                f"mol_post_{reaction_id}"
            )

            lines.append(
                f"{'molecule':<16} "
                f"{pre_id:<16} "
                f"{pre_file}"
            )

            lines.append(
                f"{'molecule':<16} "
                f"{post_id:<16} "
                f"{post_file}\n"
            )

            rxn_stp = (
                f"rxn_stp_{reaction_id}"
            )

            rxn_str = (
                f"react "
                f"{rxn_stp:<15} "
                f"all 1 0.0 3.5 "
                f"{pre_id:<14} "
                f"{post_id:<15} "
                f"{map_file:<15} "
                f"stabilize_steps 60 "
                f"rescale_charges yes"
            )

            rxn_commands.append(
                rxn_str
            )

        if not rxn_commands:
            raise ValueError(
                "No active reaction templates were found."
            )

        all_reactions = (
            " & \n                "
        ).join(rxn_commands)

        lines.extend([
            "",
            (
                f"{'fix':<16}"
                f"rxns all bond/react "
                f"stabilization yes statted_grp 0.03 &"
            ),
            f"{'':<16}{all_reactions}",
            "",
            "",
            "# Note: If atoms are being deleted during the reaction, ensure you use the correct Map file",
            "#       (e.g., RXN_i_with_delete_ids.map).",
            "#       NPT is recommended for deletion to account for density changes.",
            "#       If kinetics estimation is enabled, Arrhenius constraints are written",
            "#       directly into the copied active reaction map file.\n",
            (
                f"{'fix':<16} "
                f"1 statted_grp_REACT nvt "
                f"temp {simulation.temperature} "
                f"{simulation.temperature} "
                f"100.0\n"
            ),
            (
                f"#{'fix':<16} "
                f"1 statted_grp_REACT npt "
                f"temp {simulation.temperature} "
                f"{simulation.temperature} "
                f"100.0 iso 0.0 0.0 1000.0"
            ),
            "",
            (
                f"{'thermo_style':<16} "
                f"custom step time temp "
                f"f_rxns[*] press density "
                f"vol pe ke etotal"
            ),
            (
                f"{'dump':<16} "
                f"traj all xyz 1000 "
                f"{output_base}.xyz"
            ),
            f"{'dump_modify':<16} traj types labels",
            "",
            f"{'run':<16} 1000000",
            (
                f"{'restart':<16} "
                f"100 "
                f"{output_base}_backup1.restart "
                f"{output_base}_backup2.restart"
            ),
            (
                f"{'write_restart':<16} "
                f"{output_base}.restart"
            ),
            (
                f"{'write_data':<16} "
                f"{output_base}.data nofix"
            ),
        ])

        # ----- Write the .in file --------------------------------------
        in_file_path = (
            rxn_dir
            / f"in.{tag}_reaction"
        )

        with open(
            in_file_path,
            "w",
        ) as f:
            f.write(
                "\n".join(lines)
            )

        # ----- Copy auxiliary files into the reaction directory ---------
        self._copy_required_files(
            dest_dir=rxn_dir,
            simulation=simulation,
        )

        return in_file_path.name

    #  Internal helpers
    def _copy_required_files(
        self,
        dest_dir: Path,
        simulation: Simulation,
    ) -> None:
        """Copy every map and molecule file referenced by the reaction templates.

        The standard RXN_N.map file is always copied and is the map used by
        AutoREACTER's generated LAMMPS script.

        When an optional RXN_N_with_delete_ids.map file is available, it is
        also copied into the reaction directory for the user. AutoREACTER does
        not automatically use the supplementary DeleteIDs map.

        When kinetics estimation is enabled, RMG kinetics are estimated for
        the current reaction and simulation after the files are copied. The
        resulting Arrhenius constraint is written directly into the copied
        active map file. The original AutoREACTER map file is never modified.

        Parameters
        ----------
        dest_dir : Path
            Target directory (``3_reaction_first_stage``).
        simulation : Simulation
            Simulation object containing the temperature and other conditions
            required for kinetics estimation.

        Raises
        ------
        FileNotFoundError
            If any required standard reaction file does not exist on disk.
        """
        rf = self.reacter_files

        for template in [
            t for t in rf.template_files
            if getattr(t, "activity_stats", True)
        ]:
            files: list[Path] = [
                template.map_file,
                template.pre_reaction_file,
                template.post_reaction_file,
            ]

            # Supplementary DeleteIDs map is optional.
            map_file_with_delete_ids = getattr(
                template,
                "map_file_with_delete_ids",
                None,
            )

            if map_file_with_delete_ids is not None:
                files.append(
                    map_file_with_delete_ids
                )

            # ----- Copy all required reaction files ---------------------
            for file in files:
                if (
                    file is None
                    or not file.exists()
                ):
                    raise FileNotFoundError(
                        f"Required reaction file not found: {file}"
                    )

                shutil.copy2(
                    file,
                    dest_dir / file.name,
                )

            # ----- Estimate kinetics and modify the copied map file -----
            #
            # Kinetics are not stored on ReactionMetadata. Each reaction is
            # evaluated for the current simulation conditions, written directly
            # into this simulation's copied map file, and then discarded.
            #
            # If kinetics estimation was not enabled in Writer, the estimator
            # is None and the normal map file is left unchanged.
            if self.kinetics_estimator is None:
                continue

            copied_map_file = (
                dest_dir
                / template.map_file.name
            )

            self.kinetics_estimator.estimate_and_write_map(
                reaction=template,
                simulation=simulation,
                map_file=copied_map_file,
                lammps_units=self.settings.units,
            )