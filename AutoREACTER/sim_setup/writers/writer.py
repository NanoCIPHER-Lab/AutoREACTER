from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from AutoREACTER.sim_setup.writers.lammps_settings import (
    LammpsInitialSettings,
)
from AutoREACTER.reaction_preparation.ff_wrapper.REACTER_files_builder import (
    REACTERFiles,
)
from AutoREACTER.sim_setup.writers.densification_writer import (
    DensificationWriter,
)
from AutoREACTER.sim_setup.writers.pre_eq_writer import (
    PreEqWriter,
)
from AutoREACTER.sim_setup.writers.rxn_first_stage_writer import (
    RxnFirstStageWriter,
)
from AutoREACTER.sim_setup.writers.rxn_second_stage_writer import (
    RxnSecondStageWriter,
)
from AutoREACTER.sim_setup.writers.post_eq_writer import (
    PostEqWriter,
)
from AutoREACTER.input_parser import SimulationSetup


if TYPE_CHECKING:
    from AutoREACTER.kinetics.estimator import KineticsEstimator
    from AutoREACTER.session import Session


class Writer:

    def __init__(
        self,
        reacter_files: REACTERFiles,
        session: "Session",
    ):
        self.reacter_files = reacter_files

        self.lammps_initial_setup = LammpsInitialSettings(
            reacter_files
        )

        self.settings = (
            self.lammps_initial_setup.get_LUNAR_lammps_settings()
        )

        # ----- Optional kinetics estimation -----------------------------
        #
        # RMG is an optional dependency. It is imported and initialized only
        # when kinetics estimation is explicitly enabled for the session.
        #
        # The estimator is created once here and reused for every simulation
        # and reaction. This prevents the RMG kinetics database from being
        # repeatedly loaded for every reaction map.
        self.kinetics_estimator: "KineticsEstimator | None" = None

        if session.kinetics:
            from AutoREACTER.kinetics.estimator import (
                KineticsEstimator,
            )

            self.kinetics_estimator = KineticsEstimator()

    def write_all_files(
        self,
        run_dir: Path,
        simulation_setup: SimulationSetup,
    ) -> None:
        sim_name = simulation_setup.simulation_name

        lammps_dir = (
            run_dir
            / "LAMMPS_input_files"
        )

        lammps_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        simulations = simulation_setup.simulations

        for simulation in simulations:
            sub_dir = (
                lammps_dir
                / f"{sim_name}_{simulation.tag}"
            )

            sub_dir.mkdir(
                parents=True,
                exist_ok=True,
            )

            # ----- Densification ----------------------------------------
            DensificationWriter(
                out_dir=sub_dir,
                settings=self.settings,
                reacter_files=self.reacter_files,
                simulation=simulation,
                sim_name=sim_name,
            )

            # ----- Pre-equilibration -----------------------------------
            PreEqWriter(
                out_dir=sub_dir,
                settings=self.settings,
                simulation=simulation,
                sim_name=sim_name,
            )

            # ----- First reaction stage --------------------------------
            #
            # The optional kinetics estimator is passed to the first-stage
            # writer. If kinetics estimation is disabled, this value is None
            # and the reaction maps are copied without modification.
            #
            # If enabled, the already-loaded estimator calculates the kinetics
            # for each active reaction at the current simulation conditions and
            # writes the Arrhenius constraint directly into the copied map file.
            RxnFirstStageWriter(
                out_dir=sub_dir,
                settings=self.settings,
                reacter_files=self.reacter_files,
                simulation=simulation,
                sim_name=sim_name,
                kinetics_estimator=self.kinetics_estimator,
            )

            # ----- Optional second reaction stage -----------------------
            if simulation_setup.write_second_reaction_stage:
                RxnSecondStageWriter(
                    out_dir=sub_dir,
                    settings=self.settings,
                    reacter_files=self.reacter_files,
                    simulation=simulation,
                    sim_name=sim_name,
                )

            # ----- Post-equilibration ----------------------------------
            PostEqWriter(
                out_dir=sub_dir,
                settings=self.settings,
                simulation=simulation,
                sim_name=sim_name,
                write_second_reaction_stage=(
                    simulation_setup.write_second_reaction_stage
                ),
            )