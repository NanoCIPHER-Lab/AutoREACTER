from __future__ import annotations

import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem


# ============================================================================
# Exceptions
# ============================================================================


class KineticsEstimationError(Exception):
    """Raised when complete and reliable kinetics cannot be resolved."""


# ============================================================================
# Normalized kinetics container
# ============================================================================


@dataclass(frozen=True)
class ResolvedKinetics:
    """Normalized kinetics used internally by AutoREACTER.

    Parameters
    ----------
    A
        Arrhenius pre-exponential factor.

    n
        Temperature exponent.

    Ea_j_mol
        Activation energy in J/mol.

    T0_K
        Reference temperature used by modified Arrhenius expressions.

    source
        Human-readable provenance.

    source_kind
        Usually ``"curated"`` or ``"rmg"``.

    mode
        ``"lammps_probability"``
            Parameters already represent an Arrhenius probability model
            appropriate for the LAMMPS bond/react map constraint.

        ``"physical_rate"``
            Parameters represent a physical kinetic rate constant and require
            an explicit probability conversion/calibration factor before they
            can be used by LAMMPS.

    A_units
        Original units of A when known.

    probability_scale
        Conversion/calibration factor used to map physical kinetics to the
        dimensionless LAMMPS reaction probability.

    family
        RMG family when applicable.

    estimated
        True when RMG obtained the kinetics through a rate-rule estimate rather
        than direct kinetics.
    """

    A: float
    n: float
    Ea_j_mol: float
    T0_K: float

    source: str
    source_kind: str
    mode: str

    A_units: str | None = None
    probability_scale: float | None = None

    family: str | None = None
    estimated: bool = False
    comment: str | None = None


# ============================================================================
# Main estimator
# ============================================================================


class KineticsEstimator:
    """Resolve reaction kinetics and write LAMMPS Arrhenius constraints.

    Resolution order
    ----------------
    1. User / curated AutoREACTER kinetics.
    2. Exact RMG reaction match.
    3. RMG reactant-only generation followed by strict product isomorphism.
    4. Failure.

    AutoREACTER never accepts:
    - fingerprint-only similarity,
    - formula-only similarity,
    - nearest chemical reaction,
    - incomplete kinetics,
    - ambiguous competing kinetics.

    RMG rate-rule interpolation is allowed only after RMG has identified the
    correct chemical transformation.

    Notes
    -----
    RMG returns physical rate constants. LAMMPS ``fix bond/react`` uses its
    Arrhenius constraint as a reaction acceptance probability.

    Therefore physical RMG kinetics require an explicit ``probability_scale``
    before they are written to a map file.
    """

    AMBIENT_PRESSURE = 1.0e5  # Pa

    # Boltzmann constant expressed in kcal/mol/K for LAMMPS "real" units.
    KB_REAL = 0.00198720425864083

    CURATED_REACTION_ATTRIBUTES = (
        "curated_kinetics",
        "kinetics_data",
        "user_kinetics",
    )

    CURATED_SESSION_ATTRIBUTES = (
        "curated_kinetics",
        "kinetics_registry",
        "kinetics_data",
    )

    def __init__(
        self,
        session: Any | None = None,
    ):
        """Initialize the kinetics resolver and RMG database."""

        self.session = session

        self.database: Any | None = None

        self._rmg_molecule_class: Any | None = None
        self._rmg_species_class: Any | None = None
        self._from_rdkit_mol: Any | None = None

        self._info(
            "Initializing kinetics estimator."
        )

        # ------------------------------------------------------------------
        # RMG is optional from AutoREACTER's perspective, but if kinetics
        # resolution reaches the RMG path it must be installed correctly.
        # ------------------------------------------------------------------
        try:
            from rmgpy import settings
            from rmgpy.data.kinetics.database import KineticsDatabase
            from rmgpy.molecule.converter import from_rdkit_mol
            from rmgpy.molecule.molecule import Molecule
            from rmgpy.species import Species

        except ImportError as exc:
            self._fail(
                "RMG is required for automatic kinetics estimation. "
                "Install RMG-Py and configure RMG-database, or provide "
                "curated kinetics for every reaction.",
                cause=exc,
            )

        self._rmg_molecule_class = Molecule
        self._rmg_species_class = Species
        self._from_rdkit_mol = from_rdkit_mol

        # ------------------------------------------------------------------
        # Locate RMG kinetics database.
        # ------------------------------------------------------------------
        kinetics_path = os.path.join(
            settings["database.directory"],
            "kinetics",
        )

        self._info(
            f"RMG kinetics database path: {kinetics_path}"
        )

        if not os.path.isdir(
            kinetics_path
        ):
            self._fail(
                "RMG kinetics database directory does not exist: "
                f"{kinetics_path}"
            )

        # ------------------------------------------------------------------
        # Load family rate rules.
        #
        # We currently do not load every kinetics library because family
        # matching / rate-rule estimation is the main automatic path.
        # ------------------------------------------------------------------
        self._info(
            "Loading all RMG reaction families. "
            "Kinetics libraries are disabled for this integration."
        )

        self.database = KineticsDatabase()

        try:
            self.database.load(
                kinetics_path,
                families="all",
                libraries=[],
            )

        except Exception as exc:
            self._fail(
                "RMG kinetics database could not be loaded.",
                cause=exc,
            )

        family_count = len(
            getattr(
                self.database,
                "families",
                {},
            )
        )

        self._info(
            "RMG kinetics database loaded successfully with "
            f"{family_count} reaction families."
        )

    # ======================================================================
    # Logging
    # ======================================================================

    @staticmethod
    def _info(
        message: str,
    ) -> None:
        """Print an informational message."""

        print(
            f"[INFO] {message}"
        )

    @staticmethod
    def _warning(
        message: str,
    ) -> None:
        """Print a warning message."""

        print(
            f"[WARNING] {message}"
        )

    @staticmethod
    def _fail(
        message: str,
        cause: Exception | None = None,
    ) -> None:
        """Print a warning and raise a kinetics error."""

        print(
            f"[WARNING] {message}"
        )

        if cause is None:
            raise KineticsEstimationError(
                message
            )

        raise KineticsEstimationError(
            message
        ) from cause

    # ======================================================================
    # Public API
    # ======================================================================

    def estimate_and_write_map(
        self,
        reaction: Any,
        simulation: Any,
        map_file: Path,
        lammps_units: str,
    ) -> bool:
        """Resolve kinetics and write the Arrhenius map constraint.

        The map file is modified only after:
        - chemistry resolution succeeds,
        - kinetics are complete,
        - unit conversion succeeds,
        - LAMMPS probability conversion succeeds,
        - the resulting probability at the simulation temperature is valid.
        """

        reaction_id = getattr(
            reaction,
            "reaction_id",
            "unknown",
        )

        try:
            temperature = float(
                simulation.temperature
            )

        except (
            TypeError,
            ValueError,
            AttributeError,
        ) as exc:
            self._fail(
                f"Reaction {reaction_id} does not have a valid "
                "simulation temperature.",
                cause=exc,
            )

        if (
            not math.isfinite(
                temperature
            )
            or temperature <= 0.0
        ):
            self._fail(
                f"Reaction {reaction_id} has invalid temperature "
                f"{temperature} K."
            )

        self._info(
            f"Estimating kinetics for AutoREACTER reaction "
            f"{reaction_id} at {temperature:.2f} K."
        )

        # ------------------------------------------------------------------
        # Resolve kinetics from curated data first, then RMG.
        # ------------------------------------------------------------------
        kinetics = self._resolve_kinetics(
            reaction=reaction,
        )

        self._validate_resolved_kinetics(
            kinetics=kinetics,
            reaction_id=reaction_id,
        )

        # ------------------------------------------------------------------
        # Convert normalized kinetics to LAMMPS map parameters.
        # ------------------------------------------------------------------
        A_lammps, n_lammps, Ea_lammps = (
            self._convert_to_lammps_probability_parameters(
                kinetics=kinetics,
                reaction=reaction,
                lammps_units=lammps_units,
                reaction_id=reaction_id,
            )
        )

        # ------------------------------------------------------------------
        # Validate actual probability at the simulation temperature.
        # ------------------------------------------------------------------
        probability = self._evaluate_lammps_probability(
            A=A_lammps,
            n=n_lammps,
            Ea=Ea_lammps,
            temperature=temperature,
            lammps_units=lammps_units,
            reaction_id=reaction_id,
        )

        self._info(
            f"Reaction {reaction_id}: Arrhenius probability at "
            f"{temperature:.2f} K = {probability:.6e}"
        )

        seed = random.randint(
            10_000,
            2_147_483_646,
        )

        # ------------------------------------------------------------------
        # Nothing has touched the map until this point.
        # ------------------------------------------------------------------
        self._write_arrhenius_constraint(
            map_file=map_file,
            A=A_lammps,
            n=n_lammps,
            Ea=Ea_lammps,
            seed=seed,
        )

        self._info(
            f"Reaction {reaction_id}: kinetics successfully written."
        )

        self._info(
            f"Kinetics source: {kinetics.source}"
        )

        if kinetics.family:
            self._info(
                f"RMG reaction family: {kinetics.family}"
            )

        self._info(
            "LAMMPS Arrhenius parameters: "
            f"A={A_lammps:.8e}, "
            f"n={n_lammps:.8g}, "
            f"Ea={Ea_lammps:.8g}"
        )

        self._info(
            f"Updated reaction map: {map_file.name}"
        )

        return True

    # ======================================================================
    # Main resolver
    # ======================================================================

    def _resolve_kinetics(
        self,
        reaction: Any,
    ) -> ResolvedKinetics:
        """Resolve the best scientifically supported kinetics source."""

        reaction_id = getattr(
            reaction,
            "reaction_id",
            "unknown",
        )

        # ==================================================================
        # 1. Curated / user-provided kinetics
        # ==================================================================
        curated_payload = self._find_curated_kinetics(
            reaction=reaction,
        )

        if curated_payload is not None:
            self._info(
                f"Reaction {reaction_id}: curated kinetics found."
            )

            resolved = self._normalize_curated_kinetics(
                payload=curated_payload,
                reaction_id=reaction_id,
            )

            self._info(
                f"Reaction {reaction_id}: using curated kinetics "
                f"from {resolved.source}."
            )

            return resolved

        # ==================================================================
        # 2. RMG
        # ==================================================================
        self._info(
            f"Reaction {reaction_id}: no curated kinetics found."
        )

        self._info(
            f"Reaction {reaction_id}: attempting RMG kinetics resolution."
        )

        try:
            return self._estimate_rmg_kinetics(
                reaction=reaction,
            )

        except KineticsEstimationError as exc:
            self._warning(
                f"Reaction {reaction_id}: RMG could not provide "
                f"reliable kinetics."
            )

            self._warning(
                str(
                    exc
                )
            )

        # ==================================================================
        # 3. Nothing reliable
        # ==================================================================
        self._fail(
            f"No reliable kinetics source is available for reaction "
            f"{reaction_id}. Provide curated/user kinetics for this "
            "reaction or add support for the appropriate reaction family."
        )

    # ======================================================================
    # Curated kinetics
    # ======================================================================

    def _find_curated_kinetics(
        self,
        reaction: Any,
    ) -> dict[str, Any] | None:
        """Look for curated kinetics on the reaction or Session.

        Supported reaction attributes
        -----------------------------
        reaction.curated_kinetics
        reaction.kinetics_data
        reaction.user_kinetics

        Supported Session registries
        ----------------------------
        session.curated_kinetics
        session.kinetics_registry
        session.kinetics_data

        Session registries may be dictionaries keyed by reaction ID.
        """

        reaction_id = getattr(
            reaction,
            "reaction_id",
            None,
        )

        # ------------------------------------------------------------------
        # Direct reaction metadata.
        # ------------------------------------------------------------------
        for attribute in self.CURATED_REACTION_ATTRIBUTES:
            payload = getattr(
                reaction,
                attribute,
                None,
            )

            if isinstance(
                payload,
                dict,
            ):
                return payload

        # ------------------------------------------------------------------
        # Optional Session-level registry.
        # ------------------------------------------------------------------
        if self.session is None:
            return None

        for attribute in self.CURATED_SESSION_ATTRIBUTES:
            registry = getattr(
                self.session,
                attribute,
                None,
            )

            if not isinstance(
                registry,
                dict,
            ):
                continue

            # --------------------------------------------------------------
            # Registry may itself be one kinetics payload.
            # --------------------------------------------------------------
            if self._looks_like_kinetics_payload(
                registry
            ):
                return registry

            # --------------------------------------------------------------
            # Or it may be keyed by reaction ID.
            # --------------------------------------------------------------
            keys_to_try = [
                reaction_id,
                str(
                    reaction_id
                ),
            ]

            for key in keys_to_try:
                if key in registry:
                    payload = registry[
                        key
                    ]

                    if isinstance(
                        payload,
                        dict,
                    ):
                        return payload

        return None

    @staticmethod
    def _looks_like_kinetics_payload(
        payload: dict[str, Any],
    ) -> bool:
        """Return True when a dictionary appears to contain kinetics."""

        return all(
            key in payload
            for key in (
                "A",
                "n",
                "Ea",
            )
        )

    def _normalize_curated_kinetics(
        self,
        payload: dict[str, Any],
        reaction_id: Any,
    ) -> ResolvedKinetics:
        """Validate and normalize user/curated kinetics."""

        required = (
            "A",
            "n",
            "Ea",
            "Ea_units",
            "source",
            "mode",
        )

        missing = [
            key
            for key in required
            if key not in payload
        ]

        if missing:
            self._fail(
                f"Curated kinetics for reaction {reaction_id} are "
                f"incomplete. Missing: {missing}"
            )

        try:
            A = float(
                payload["A"]
            )

            n = float(
                payload["n"]
            )

            Ea_input = float(
                payload["Ea"]
            )

            T0 = float(
                payload.get(
                    "T0",
                    1.0,
                )
            )

        except (
            TypeError,
            ValueError,
        ) as exc:
            self._fail(
                f"Curated kinetics for reaction {reaction_id} contain "
                "non-numeric Arrhenius parameters.",
                cause=exc,
            )

        source = str(
            payload["source"]
        ).strip()

        mode = str(
            payload["mode"]
        ).strip().lower()

        Ea_units = str(
            payload["Ea_units"]
        ).strip()

        A_units_raw = payload.get(
            "A_units"
        )

        A_units = (
            str(
                A_units_raw
            ).strip()
            if A_units_raw is not None
            else None
        )

        if not source:
            self._fail(
                f"Curated kinetics for reaction {reaction_id} must "
                "include a non-empty source."
            )

        if mode not in (
            "lammps_probability",
            "physical_rate",
        ):
            self._fail(
                f"Curated kinetics for reaction {reaction_id} use "
                f"unsupported mode {mode!r}. Expected "
                "'lammps_probability' or 'physical_rate'."
            )

        if mode == "physical_rate" and not A_units:
            self._fail(
                f"Physical-rate kinetics for reaction {reaction_id} "
                "must include A_units."
            )

        Ea_j_mol = self._activation_energy_to_j_per_mol(
            value=Ea_input,
            units=Ea_units,
            reaction_id=reaction_id,
        )

        probability_scale = payload.get(
            "probability_scale"
        )

        if probability_scale is not None:
            try:
                probability_scale = float(
                    probability_scale
                )

            except (
                TypeError,
                ValueError,
            ) as exc:
                self._fail(
                    f"probability_scale for reaction {reaction_id} "
                    "must be numeric.",
                    cause=exc,
                )

        estimated = bool(
            payload.get(
                "estimated",
                False,
            )
        )

        comment = payload.get(
            "comment"
        )

        if comment is not None:
            comment = str(
                comment
            )

        resolved = ResolvedKinetics(
            A=A,
            n=n,
            Ea_j_mol=Ea_j_mol,
            T0_K=T0,
            source=source,
            source_kind="curated",
            mode=mode,
            A_units=A_units,
            probability_scale=probability_scale,
            family=None,
            estimated=estimated,
            comment=comment,
        )

        self._validate_resolved_kinetics(
            kinetics=resolved,
            reaction_id=reaction_id,
        )

        if estimated:
            self._warning(
                f"Reaction {reaction_id}: curated kinetics are marked "
                "as estimated."
            )

        return resolved

    # ======================================================================
    # RMG kinetics
    # ======================================================================

    def _estimate_rmg_kinetics(
        self,
        reaction: Any,
    ) -> ResolvedKinetics:
        """Find the exact AutoREACTER reaction in RMG."""

        if self.database is None:
            self._fail(
                "RMG kinetics database is not initialized."
            )

        reaction_id = getattr(
            reaction,
            "reaction_id",
            "unknown",
        )

        reactant_rdkit = getattr(
            reaction,
            "reactant_combined_RDmol",
            None,
        )

        product_rdkit = getattr(
            reaction,
            "product_combined_RDmol",
            None,
        )

        if reactant_rdkit is None:
            self._fail(
                f"Reaction {reaction_id} is missing "
                "reactant_combined_RDmol."
            )

        if product_rdkit is None:
            self._fail(
                f"Reaction {reaction_id} is missing "
                "product_combined_RDmol."
            )

        reactants = self._combined_rdkit_to_rmg_molecules(
            reactant_rdkit
        )

        products = self._combined_rdkit_to_rmg_molecules(
            product_rdkit
        )

        self._info(
            f"Reaction {reaction_id}: converted "
            f"{len(reactants)} reactant molecule(s) and "
            f"{len(products)} product molecule(s) to RMG."
        )

        # ------------------------------------------------------------------
        # Optional known family.
        # ------------------------------------------------------------------
        family_hint = getattr(
            reaction,
            "rmg_family",
            None,
        )

        only_families = None

        if family_hint:
            family_hint = str(
                family_hint
            )

            if family_hint not in self.database.families:
                self._fail(
                    f"Reaction {reaction_id} specifies RMG family "
                    f"{family_hint!r}, but that family is not loaded."
                )

            only_families = [
                family_hint
            ]

            self._info(
                f"Reaction {reaction_id}: restricting RMG search to "
                f"family {family_hint}."
            )

        # ==================================================================
        # Stage 1 — exact product-constrained RMG search
        # ==================================================================
        self._info(
            f"Reaction {reaction_id}: attempting normal RMG "
            "product-constrained matching."
        )

        try:
            kwargs = {
                "reactants": reactants,
                "products": products,
            }

            if only_families is not None:
                kwargs["only_families"] = only_families

            exact_candidates = self.database.generate_reactions(
                **kwargs
            )

        except Exception as exc:
            self._warning(
                f"Reaction {reaction_id}: direct RMG matching raised "
                f"an exception: {exc}"
            )

            exact_candidates = []

        if exact_candidates:
            self._info(
                f"Reaction {reaction_id}: RMG returned "
                f"{len(exact_candidates)} direct candidate(s)."
            )

            return self._select_rmg_candidate(
                candidates=exact_candidates,
                reaction_id=reaction_id,
                source_description="direct RMG product match",
            )

        # ==================================================================
        # Stage 2 — reactants only, then strict product graph matching
        # ==================================================================
        self._warning(
            f"Reaction {reaction_id}: no direct product-constrained "
            "RMG match was found."
        )

        self._warning(
            "Falling back to RMG reactant-only generation followed by "
            "strict Molecule-level product isomorphism matching "
            "(including resonance structures)."
        )

        self._warning(
            "This fallback does NOT accept unrelated reactions, "
            "formula-only matches, or fingerprint similarity."
        )

        try:
            kwargs = {
                "reactants": reactants,
            }

            if only_families is not None:
                kwargs["only_families"] = only_families

            generated_candidates = self.database.generate_reactions(
                **kwargs
            )

        except Exception as exc:
            self._fail(
                f"RMG could not generate reaction candidates for "
                f"reaction {reaction_id}.",
                cause=exc,
            )

        if not generated_candidates:
            self._fail(
                f"RMG generated no reaction candidates for "
                f"reaction {reaction_id}."
            )

        self._info(
            f"Reaction {reaction_id}: RMG generated "
            f"{len(generated_candidates)} candidate reaction(s) "
            "from the reactants."
        )

        matching_candidates = []

        for candidate in generated_candidates:
            if self._candidate_products_match(
                candidate=candidate,
                target_products=products,
            ):
                matching_candidates.append(
                    candidate
                )

        if not matching_candidates:
            self._warning(
                f"Reaction {reaction_id}: RMG generated "
                f"{len(generated_candidates)} candidate reaction(s), "
                "but none reproduced the intended AutoREACTER "
                "product graphs."
            )

            self._warning(
                "This indicates that the intended chemistry is not "
                "represented by the currently loaded RMG reaction-family "
                "database, or that its mechanism differs from the "
                "AutoREACTER transformation."
            )

            self._fail(
                f"No reliable RMG kinetics are available for reaction "
                f"{reaction_id}. Curated kinetics are required."
            )

        self._info(
            f"Reaction {reaction_id}: "
            f"{len(matching_candidates)} candidate(s) passed strict "
            "product-isomorphism validation."
        )

        return self._select_rmg_candidate(
            candidates=matching_candidates,
            reaction_id=reaction_id,
            source_description="strict RMG isomorphism fallback",
        )

    # ======================================================================
    # Strict product matching
    # ======================================================================

    def _candidate_products_match(
        self,
        candidate: Any,
        target_products: list[Any],
    ) -> bool:
        """Require exact graph-isomorphic product sets.

        Resonance structures are considered, but molecular formulas or
        fingerprint similarity are not sufficient.
        """

        candidate_products = getattr(
            candidate,
            "products",
            None,
        )

        if candidate_products is None:
            return False

        if len(
            candidate_products
        ) != len(
            target_products
        ):
            return False

        candidate_option_sets = [
            self._rmg_entity_molecule_options(
                entity
            )
            for entity in candidate_products
        ]

        target_option_sets = [
            self._rmg_molecule_resonance_options(
                molecule
            )
            for molecule in target_products
        ]

        return self._match_product_sets(
            candidate_sets=candidate_option_sets,
            target_sets=target_option_sets,
        )

    def _match_product_sets(
        self,
        candidate_sets: list[list[Any]],
        target_sets: list[list[Any]],
    ) -> bool:
        """Find a one-to-one isomorphic assignment between product sets."""

        if not candidate_sets:
            return not target_sets

        first_candidate = candidate_sets[0]

        for index, target_options in enumerate(
            target_sets
        ):
            matched = False

            for candidate_molecule in first_candidate:
                for target_molecule in target_options:
                    try:
                        if candidate_molecule.is_isomorphic(
                            target_molecule
                        ):
                            matched = True
                            break

                    except Exception:
                        continue

                if matched:
                    break

            if not matched:
                continue

            remaining_targets = (
                target_sets[:index]
                + target_sets[
                    index + 1:
                ]
            )

            if self._match_product_sets(
                candidate_sets=candidate_sets[1:],
                target_sets=remaining_targets,
            ):
                return True

        return False

    def _rmg_entity_molecule_options(
        self,
        entity: Any,
    ) -> list[Any]:
        """Extract resonance molecule options from Species or Molecule."""

        species_molecules = getattr(
            entity,
            "molecule",
            None,
        )

        if species_molecules:
            return list(
                species_molecules
            )

        if hasattr(
            entity,
            "is_isomorphic",
        ):
            return self._rmg_molecule_resonance_options(
                entity
            )

        return []

    def _rmg_molecule_resonance_options(
        self,
        molecule: Any,
    ) -> list[Any]:
        """Generate resonance structures for one RMG Molecule."""

        if self._rmg_species_class is None:
            return [
                molecule
            ]

        try:
            species = self._rmg_species_class(
                molecule=[
                    molecule.copy(
                        deep=True
                    )
                ]
            )

            species.generate_resonance_structures(
                keep_isomorphic=True,
                filter_structures=True,
                save_order=False,
            )

            if species.molecule:
                return list(
                    species.molecule
                )

        except Exception:
            pass

        return [
            molecule
        ]

    # ======================================================================
    # RMG candidate selection
    # ======================================================================

    def _select_rmg_candidate(
        self,
        candidates: list[Any],
        reaction_id: Any,
        source_description: str,
    ) -> ResolvedKinetics:
        """Select one complete and non-conflicting RMG result."""

        valid: list[
            tuple[Any, ResolvedKinetics]
        ] = []

        for index, candidate in enumerate(
            candidates,
            start=1,
        ):
            kinetics = getattr(
                candidate,
                "kinetics",
                None,
            )

            family = self._reaction_family_label(
                candidate
            )

            if kinetics is None:
                self._warning(
                    f"Reaction {reaction_id}: RMG candidate {index} "
                    f"from family {family} contains no kinetics."
                )

                continue

            try:
                normalized = self._normalize_rmg_kinetics(
                    candidate=candidate,
                    kinetics=kinetics,
                    reaction_id=reaction_id,
                )

            except KineticsEstimationError as exc:
                self._warning(
                    f"Reaction {reaction_id}: RMG candidate {index} "
                    f"from family {family} is incomplete: {exc}"
                )

                continue

            valid.append(
                (
                    candidate,
                    normalized,
                )
            )

        if not valid:
            self._fail(
                f"RMG matched reaction {reaction_id}, but none of the "
                "matching reactions supplied complete Arrhenius "
                "A, n, and Ea data."
            )

        # ------------------------------------------------------------------
        # One valid result.
        # ------------------------------------------------------------------
        if len(valid) == 1:
            result = valid[0][1]

            self._report_rmg_result(
                kinetics=result,
                reaction_id=reaction_id,
                source_description=source_description,
            )

            return result

        # ------------------------------------------------------------------
        # Multiple RMG representations can be resonance / degeneracy
        # duplicates. Collapse them only if the normalized kinetics agree.
        # ------------------------------------------------------------------
        signatures = {
            self._resolved_signature(
                kinetics
            )
            for _, kinetics in valid
        }

        if len(signatures) == 1:
            self._warning(
                f"Reaction {reaction_id}: RMG returned "
                f"{len(valid)} equivalent complete kinetics candidates."
            )

            self._warning(
                "Their normalized Arrhenius data are equivalent; "
                "using the first representation."
            )

            result = valid[0][1]

            self._report_rmg_result(
                kinetics=result,
                reaction_id=reaction_id,
                source_description=source_description,
            )

            return result

        families = sorted(
            {
                kinetics.family or "unknown"
                for _, kinetics in valid
            }
        )

        self._fail(
            f"Reaction {reaction_id} has {len(valid)} conflicting "
            f"complete RMG kinetics candidates from families "
            f"{families}. AutoREACTER will not guess between them."
        )

    def _normalize_rmg_kinetics(
        self,
        candidate: Any,
        kinetics: Any,
        reaction_id: Any,
    ) -> ResolvedKinetics:
        """Convert one RMG kinetics object into normalized data."""

        model_name = type(
            kinetics
        ).__name__

        # ------------------------------------------------------------------
        # AutoREACTER requires explicit A/n/Ea.
        # ------------------------------------------------------------------
        for parameter in (
            "A",
            "n",
            "Ea",
        ):
            if getattr(
                kinetics,
                parameter,
                None,
            ) is None:
                raise KineticsEstimationError(
                    f"RMG kinetics model {model_name} does not "
                    f"provide required parameter {parameter}."
                )

        A = self._rmg_quantity_value_si(
            kinetics.A,
            "A",
        )

        n = self._rmg_quantity_value_si(
            kinetics.n,
            "n",
        )

        Ea = self._rmg_quantity_value_si(
            kinetics.Ea,
            "Ea",
        )

        T0_quantity = getattr(
            kinetics,
            "T0",
            None,
        )

        if T0_quantity is None:
            T0 = 1.0

        else:
            T0 = self._rmg_quantity_value_si(
                T0_quantity,
                "T0",
            )

        A_units = getattr(
            kinetics.A,
            "units",
            None,
        )

        if A_units is not None:
            A_units = str(
                A_units
            )

        family = self._reaction_family_label(
            candidate
        )

        comment = str(
            getattr(
                kinetics,
                "comment",
                "",
            )
            or ""
        ).strip()

        lower_comment = comment.lower()

        estimation_terms = (
            "estimated",
            "average",
            "averaged",
            "rate rule",
            "rate rules",
            "template",
        )

        estimated = any(
            term in lower_comment
            for term in estimation_terms
        )

        probability_scale = self._find_probability_scale(
            reaction_id=reaction_id,
        )

        source = (
            f"RMG family {family}"
        )

        resolved = ResolvedKinetics(
            A=A,
            n=n,
            Ea_j_mol=Ea,
            T0_K=T0,
            source=source,
            source_kind="rmg",
            mode="physical_rate",
            A_units=A_units,
            probability_scale=probability_scale,
            family=family,
            estimated=estimated,
            comment=comment or None,
        )

        self._validate_resolved_kinetics(
            kinetics=resolved,
            reaction_id=reaction_id,
        )

        return resolved

    def _report_rmg_result(
        self,
        kinetics: ResolvedKinetics,
        reaction_id: Any,
        source_description: str,
    ) -> None:
        """Print provenance for an accepted RMG result."""

        self._info(
            f"Reaction {reaction_id}: selected RMG kinetics using "
            f"{source_description}."
        )

        self._info(
            f"Reaction {reaction_id}: RMG family = "
            f"{kinetics.family}"
        )

        if kinetics.A_units:
            self._info(
                f"Reaction {reaction_id}: RMG A units = "
                f"{kinetics.A_units}"
            )

        if kinetics.estimated:
            self._warning(
                f"Reaction {reaction_id}: RMG kinetics were estimated "
                "from its rate-rule/template hierarchy."
            )

            self._warning(
                "The chemical transformation is matched, but the "
                "numerical kinetics are estimated rather than direct "
                "experimental data."
            )

        else:
            self._info(
                f"Reaction {reaction_id}: no RMG rate-rule estimation "
                "warning was detected."
            )

        if kinetics.comment:
            comment = " ".join(
                kinetics.comment.split()
            )

            if len(
                comment
            ) > 400:
                comment = (
                    comment[:397]
                    + "..."
                )

            self._info(
                f"RMG kinetics source: {comment}"
            )

    # ======================================================================
    # RMG helpers
    # ======================================================================

    @staticmethod
    def _rmg_quantity_value_si(
        quantity: Any,
        name: str,
    ) -> float:
        """Return an RMG quantity in SI units."""

        if quantity is None:
            raise KineticsEstimationError(
                f"RMG parameter {name} is missing."
            )

        value = getattr(
            quantity,
            "value_si",
            quantity,
        )

        try:
            value = float(
                value
            )

        except (
            TypeError,
            ValueError,
        ) as exc:
            raise KineticsEstimationError(
                f"RMG parameter {name} is not numeric: "
                f"{value!r}"
            ) from exc

        if not math.isfinite(
            value
        ):
            raise KineticsEstimationError(
                f"RMG parameter {name} is not finite."
            )

        return value

    @staticmethod
    def _reaction_family_label(
        reaction: Any,
    ) -> str:
        """Return a readable RMG reaction family."""

        family = getattr(
            reaction,
            "family",
            None,
        )

        if family is None:
            return "unknown"

        if isinstance(
            family,
            str,
        ):
            return family

        label = getattr(
            family,
            "label",
            None,
        )

        if label:
            return str(
                label
            )

        return str(
            family
        )

    @staticmethod
    def _resolved_signature(
        kinetics: ResolvedKinetics,
    ) -> tuple[str, ...]:
        """Create a stable comparison signature."""

        return (
            kinetics.family or "",
            kinetics.mode,
            f"{kinetics.A:.12g}",
            f"{kinetics.n:.12g}",
            f"{kinetics.Ea_j_mol:.12g}",
            f"{kinetics.T0_K:.12g}",
            kinetics.A_units or "",
        )

    # ======================================================================
    # Probability-scale lookup
    # ======================================================================

    def _find_probability_scale(
        self,
        reaction_id: Any,
    ) -> float | None:
        """Find an optional physical-rate -> LAMMPS probability scale.

        Supported Session forms include:

            session.kinetics_probability_scale = 1.0e-12

        or

            session.kinetics_probability_scale = {
                1: 1.0e-12,
                2: 5.0e-13,
            }

        or equivalent string reaction IDs.
        """

        if self.session is None:
            return None

        raw = getattr(
            self.session,
            "kinetics_probability_scale",
            None,
        )

        if raw is None:
            return None

        if isinstance(
            raw,
            dict,
        ):
            value = raw.get(
                reaction_id,
                raw.get(
                    str(
                        reaction_id
                    )
                ),
            )

        else:
            value = raw

        if value is None:
            return None

        try:
            value = float(
                value
            )

        except (
            TypeError,
            ValueError,
        ):
            return None

        return value

    # ======================================================================
    # Normalized kinetics validation
    # ======================================================================

    def _validate_resolved_kinetics(
        self,
        kinetics: ResolvedKinetics,
        reaction_id: Any,
    ) -> None:
        """Require complete, finite kinetics."""

        numeric_values = {
            "A": kinetics.A,
            "n": kinetics.n,
            "Ea": kinetics.Ea_j_mol,
            "T0": kinetics.T0_K,
        }

        for name, value in numeric_values.items():
            if not math.isfinite(
                value
            ):
                self._fail(
                    f"Reaction {reaction_id}: kinetics parameter "
                    f"{name} is not finite."
                )

        if kinetics.A <= 0.0:
            self._fail(
                f"Reaction {reaction_id}: Arrhenius A must be "
                f"positive. Received {kinetics.A}."
            )

        if kinetics.T0_K <= 0.0:
            self._fail(
                f"Reaction {reaction_id}: Arrhenius T0 must be "
                f"positive. Received {kinetics.T0_K}."
            )

        if not kinetics.source:
            self._fail(
                f"Reaction {reaction_id}: kinetics have no source "
                "provenance."
            )

        if kinetics.mode not in (
            "lammps_probability",
            "physical_rate",
        ):
            self._fail(
                f"Reaction {reaction_id}: unsupported kinetics mode "
                f"{kinetics.mode!r}."
            )

    # ======================================================================
    # Activation-energy conversion
    # ======================================================================

    def _activation_energy_to_j_per_mol(
        self,
        value: float,
        units: str,
        reaction_id: Any,
    ) -> float:
        """Convert common activation-energy units to J/mol."""

        normalized = (
            units.strip()
            .lower()
            .replace(
                " ",
                "",
            )
        )

        conversions = {
            "j/mol": 1.0,
            "jmol^-1": 1.0,
            "jmol-1": 1.0,

            "kj/mol": 1000.0,
            "kjmol^-1": 1000.0,
            "kjmol-1": 1000.0,

            "cal/mol": 4.184,
            "calmol^-1": 4.184,
            "calmol-1": 4.184,

            "kcal/mol": 4184.0,
            "kcalmol^-1": 4184.0,
            "kcalmol-1": 4184.0,
        }

        if normalized not in conversions:
            self._fail(
                f"Reaction {reaction_id}: unsupported Ea units "
                f"{units!r}. Supported units are J/mol, kJ/mol, "
                "cal/mol, and kcal/mol."
            )

        result = (
            value
            * conversions[
                normalized
            ]
        )

        if not math.isfinite(
            result
        ):
            self._fail(
                f"Reaction {reaction_id}: converted activation "
                "energy is not finite."
            )

        return result

    # ======================================================================
    # Physical kinetics -> LAMMPS probability parameters
    # ======================================================================

    def _convert_to_lammps_probability_parameters(
        self,
        kinetics: ResolvedKinetics,
        reaction: Any,
        lammps_units: str,
        reaction_id: Any,
    ) -> tuple[float, float, float]:
        """Convert normalized kinetics to LAMMPS Arrhenius parameters."""

        if not isinstance(
            lammps_units,
            str,
        ):
            self._fail(
                f"Reaction {reaction_id}: invalid LAMMPS units "
                f"{lammps_units!r}."
            )

        if lammps_units.lower() != "real":
            self._fail(
                "Automatic kinetics conversion currently supports "
                "only LAMMPS 'real' units. "
                f"Reaction {reaction_id} uses {lammps_units!r}."
            )

        # ------------------------------------------------------------------
        # RMG / modified Arrhenius form:
        #
        #     A * (T / T0)^n * exp(...)
        #
        # LAMMPS:
        #
        #     A_lammps * T^n * exp(...)
        #
        # therefore:
        #
        #     A_lammps = A / T0^n
        # ------------------------------------------------------------------
        try:
            normalized_A = (
                kinetics.A
                / math.pow(
                    kinetics.T0_K,
                    kinetics.n,
                )
            )

        except (
            ValueError,
            OverflowError,
            ZeroDivisionError,
        ) as exc:
            self._fail(
                f"Reaction {reaction_id}: unable to normalize "
                "Arrhenius A/T0.",
                cause=exc,
            )

        # ------------------------------------------------------------------
        # Direct LAMMPS probability data needs no physical-rate scaling.
        # ------------------------------------------------------------------
        if kinetics.mode == "lammps_probability":
            A_lammps = normalized_A

            self._info(
                f"Reaction {reaction_id}: kinetics are already "
                "defined as a LAMMPS probability model."
            )

        # ------------------------------------------------------------------
        # Physical rates require explicit calibration.
        # ------------------------------------------------------------------
        else:
            scale = kinetics.probability_scale

            # Allow reaction metadata to provide the scale too.
            if scale is None:
                scale = getattr(
                    reaction,
                    "kinetics_probability_scale",
                    None,
                )

                if scale is not None:
                    try:
                        scale = float(
                            scale
                        )

                    except (
                        TypeError,
                        ValueError,
                    ):
                        self._fail(
                            f"Reaction {reaction_id}: "
                            "kinetics_probability_scale is not numeric."
                        )

            if scale is None:
                self._warning(
                    f"Reaction {reaction_id}: RMG supplied complete "
                    "physical Arrhenius kinetics."
                )

                if kinetics.A_units:
                    self._warning(
                        f"RMG A units are {kinetics.A_units}."
                    )

                self._warning(
                    "LAMMPS bond/react uses its Arrhenius expression "
                    "as an acceptance probability, not directly as a "
                    "physical rate coefficient."
                )

                self._fail(
                    f"Reaction {reaction_id}: a calibrated "
                    "kinetics_probability_scale is required before "
                    "physical RMG kinetics can be written to the "
                    "LAMMPS map."
                )

            if (
                not math.isfinite(
                    scale
                )
                or scale <= 0.0
            ):
                self._fail(
                    f"Reaction {reaction_id}: probability scale must "
                    f"be positive and finite. Received {scale}."
                )

            self._warning(
                f"Reaction {reaction_id}: converting physical "
                "kinetics to LAMMPS probability using calibrated "
                f"scale {scale:.8e}."
            )

            A_lammps = (
                normalized_A
                * scale
            )

        if (
            not math.isfinite(
                A_lammps
            )
            or A_lammps <= 0.0
        ):
            self._fail(
                f"Reaction {reaction_id}: resulting LAMMPS "
                f"pre-exponential factor is invalid: {A_lammps}"
            )

        # ------------------------------------------------------------------
        # J/mol -> kcal/mol for units real.
        # ------------------------------------------------------------------
        Ea_lammps = (
            kinetics.Ea_j_mol
            / 4184.0
        )

        if not math.isfinite(
            Ea_lammps
        ):
            self._fail(
                f"Reaction {reaction_id}: resulting LAMMPS "
                "activation energy is invalid."
            )

        return (
            A_lammps,
            kinetics.n,
            Ea_lammps,
        )

    # ======================================================================
    # LAMMPS probability validation
    # ======================================================================

    def _evaluate_lammps_probability(
        self,
        A: float,
        n: float,
        Ea: float,
        temperature: float,
        lammps_units: str,
        reaction_id: Any,
    ) -> float:
        """Evaluate the map Arrhenius probability at simulation temperature."""

        if lammps_units.lower() != "real":
            self._fail(
                "Probability validation currently supports only "
                "LAMMPS 'real' units."
            )

        try:
            probability = (
                A
                * math.pow(
                    temperature,
                    n,
                )
                * math.exp(
                    -Ea
                    / (
                        self.KB_REAL
                        * temperature
                    )
                )
            )

        except (
            OverflowError,
            ValueError,
            ZeroDivisionError,
        ) as exc:
            self._fail(
                f"Reaction {reaction_id}: could not evaluate "
                "LAMMPS Arrhenius probability.",
                cause=exc,
            )

        if not math.isfinite(
            probability
        ):
            self._fail(
                f"Reaction {reaction_id}: Arrhenius probability "
                "is not finite."
            )

        if probability <= 0.0:
            self._fail(
                f"Reaction {reaction_id}: Arrhenius probability "
                f"is non-positive ({probability})."
            )

        if probability > 1.0:
            self._fail(
                f"Reaction {reaction_id}: Arrhenius expression "
                f"evaluates to {probability:.6e} at "
                f"{temperature:.2f} K, which exceeds 1.0. "
                "The probability calibration is invalid."
            )

        return probability

    # ======================================================================
    # RDKit -> RMG
    # ======================================================================

    def _combined_rdkit_to_rmg_molecules(
        self,
        combined_mol: Chem.Mol,
    ) -> list[Any]:
        """Split combined RDKit structures into separate RMG molecules."""

        if combined_mol is None:
            self._fail(
                "Combined RDKit molecule is missing."
            )

        if self._rmg_molecule_class is None:
            self._fail(
                "RMG Molecule class is unavailable."
            )

        if self._from_rdkit_mol is None:
            self._fail(
                "RMG RDKit converter is unavailable."
            )

        try:
            fragments = Chem.GetMolFrags(
                combined_mol,
                asMols=True,
                sanitizeFrags=True,
            )

        except Exception as exc:
            self._fail(
                "RDKit could not split the combined reaction molecule.",
                cause=exc,
            )

        if not fragments:
            self._fail(
                "Combined RDKit molecule contains no fragments."
            )

        molecules = []

        for fragment_index, fragment in enumerate(
            fragments,
            start=1,
        ):
            molecule = self._rmg_molecule_class()

            try:
                self._from_rdkit_mol(
                    molecule,
                    fragment,
                )

            except Exception as exc:
                self._fail(
                    f"RMG could not convert RDKit fragment "
                    f"{fragment_index}.",
                    cause=exc,
                )

            vertices = getattr(
                molecule,
                "vertices",
                None,
            )

            if not vertices:
                self._fail(
                    f"RMG conversion produced an empty molecule for "
                    f"fragment {fragment_index}."
                )

            molecules.append(
                molecule
            )

        return molecules

    # ======================================================================
    # LAMMPS map writer
    # ======================================================================

    def _write_arrhenius_constraint(
        self,
        map_file: Path,
        A: float,
        n: float,
        Ea: float,
        seed: int,
    ) -> None:
        """Insert or replace one Arrhenius constraint in a map file."""

        map_file = Path(
            map_file
        )

        if not map_file.exists():
            self._fail(
                f"Reaction map file does not exist: {map_file}"
            )

        if not map_file.is_file():
            self._fail(
                f"Reaction map path is not a file: {map_file}"
            )

        values = {
            "A": A,
            "n": n,
            "Ea": Ea,
        }

        for name, value in values.items():
            if not math.isfinite(
                float(
                    value
                )
            ):
                self._fail(
                    f"Cannot write Arrhenius constraint because "
                    f"{name} is invalid: {value}"
                )

        if A <= 0.0:
            self._fail(
                f"Cannot write Arrhenius constraint because "
                f"A={A} is not positive."
            )

        if seed <= 0:
            self._fail(
                f"LAMMPS Arrhenius random seed is invalid: {seed}"
            )

        try:
            lines = (
                map_file.read_text()
                .splitlines()
            )

        except Exception as exc:
            self._fail(
                f"Could not read reaction map file {map_file}.",
                cause=exc,
            )

        if not lines:
            self._fail(
                f"Reaction map file is empty: {map_file}"
            )

        arrhenius_line = (
            f"arrhenius "
            f"{A:.8e} "
            f"{n:.8g} "
            f"{Ea:.8g} "
            f"{seed}"
        )

        # ==================================================================
        # Replace existing Arrhenius constraint
        # ==================================================================
        existing_arrhenius = [
            index
            for index, line in enumerate(
                lines
            )
            if line.strip()
            .lower()
            .startswith(
                "arrhenius "
            )
        ]

        if len(
            existing_arrhenius
        ) > 1:
            self._fail(
                f"Reaction map {map_file.name} already contains "
                "multiple Arrhenius constraints."
            )

        if existing_arrhenius:
            self._warning(
                f"{map_file.name} already contains an Arrhenius "
                "constraint. Replacing it with validated kinetics."
            )

            lines[
                existing_arrhenius[0]
            ] = arrhenius_line

            self._write_and_validate_map(
                map_file=map_file,
                lines=lines,
            )

            return

        # ==================================================================
        # Locate constraints header
        # ==================================================================
        constraints_header_index = None
        constraints_count = None

        for index, line in enumerate(
            lines
        ):
            parts = line.split()

            if (
                len(parts) == 2
                and parts[0].isdigit()
                and parts[1].lower()
                == "constraints"
            ):
                constraints_header_index = index
                constraints_count = int(
                    parts[0]
                )

                break

        # ==================================================================
        # Locate Constraints section
        # ==================================================================
        constraints_section_index = None

        for index, line in enumerate(
            lines
        ):
            if line.strip() == "Constraints":
                constraints_section_index = index
                break

        # ------------------------------------------------------------------
        # Existing header claims constraints but section is absent.
        # ------------------------------------------------------------------
        if (
            constraints_count is not None
            and constraints_count > 0
            and constraints_section_index is None
        ):
            self._fail(
                f"Reaction map {map_file.name} declares "
                f"{constraints_count} constraint(s) but contains "
                "no Constraints section."
            )

        # ==================================================================
        # Increment / create constraints count
        # ==================================================================
        if constraints_header_index is not None:
            lines[
                constraints_header_index
            ] = (
                f"{constraints_count + 1} constraints"
            )

        else:
            try:
                initiator_index = next(
                    index
                    for index, line in enumerate(
                        lines
                    )
                    if line.strip()
                    == "InitiatorIDs"
                )

            except StopIteration as exc:
                self._fail(
                    f"Invalid reaction map {map_file.name}: "
                    "InitiatorIDs section was not found.",
                    cause=exc,
                )

            insert_index = initiator_index

            while (
                insert_index > 0
                and not lines[
                    insert_index - 1
                ].strip()
            ):
                insert_index -= 1

            lines.insert(
                insert_index,
                "1 constraints",
            )

            # --------------------------------------------------------------
            # Header insertion shifts body indices.
            # --------------------------------------------------------------
            constraints_section_index = None

            for index, line in enumerate(
                lines
            ):
                if line.strip() == "Constraints":
                    constraints_section_index = index
                    break

        # ==================================================================
        # Add Arrhenius line
        # ==================================================================
        if constraints_section_index is None:
            lines.extend(
                [
                    "",
                    "Constraints",
                    "",
                    arrhenius_line,
                ]
            )

        else:
            section_names = {
                "InitiatorIDs",
                "Equivalences",
                "EdgeIDs",
                "Wildcards",
                "DeleteIDs",
                "CreateIDs",
                "ChiralIDs",
                "Constraints",
            }

            next_section_index = None

            for index in range(
                constraints_section_index + 1,
                len(
                    lines
                ),
            ):
                if (
                    lines[index].strip()
                    in section_names
                ):
                    next_section_index = index
                    break

            if next_section_index is None:
                if (
                    lines
                    and lines[-1].strip()
                ):
                    lines.append(
                        ""
                    )

                lines.append(
                    arrhenius_line
                )

            else:
                insert_index = next_section_index

                while (
                    insert_index
                    > constraints_section_index + 1
                    and not lines[
                        insert_index - 1
                    ].strip()
                ):
                    insert_index -= 1

                lines.insert(
                    insert_index,
                    arrhenius_line,
                )

        self._write_and_validate_map(
            map_file=map_file,
            lines=lines,
        )

    # ======================================================================
    # Map validation
    # ======================================================================

    def _write_and_validate_map(
        self,
        map_file: Path,
        lines: list[str],
    ) -> None:
        """Write a map and verify its Arrhenius constraint structure."""

        try:
            map_file.write_text(
                "\n".join(
                    lines
                )
                + "\n"
            )

        except Exception as exc:
            self._fail(
                f"Could not write reaction map {map_file}.",
                cause=exc,
            )

        try:
            written_lines = (
                map_file.read_text()
                .splitlines()
            )

        except Exception as exc:
            self._fail(
                f"Could not re-read modified map {map_file}.",
                cause=exc,
            )

        arrhenius_lines = [
            line
            for line in written_lines
            if line.strip()
            .lower()
            .startswith(
                "arrhenius "
            )
        ]

        if len(
            arrhenius_lines
        ) != 1:
            self._fail(
                f"Reaction map validation failed for "
                f"{map_file.name}: expected exactly one "
                "Arrhenius constraint but found "
                f"{len(arrhenius_lines)}."
            )

        if not any(
            line.strip()
            == "Constraints"
            for line in written_lines
        ):
            self._fail(
                f"Reaction map validation failed for "
                f"{map_file.name}: Constraints section is missing."
            )

        constraint_headers = []

        for line in written_lines:
            parts = line.split()

            if (
                len(parts) == 2
                and parts[0].isdigit()
                and parts[1].lower()
                == "constraints"
            ):
                constraint_headers.append(
                    int(
                        parts[0]
                    )
                )

        if len(
            constraint_headers
        ) != 1:
            self._fail(
                f"Reaction map validation failed for "
                f"{map_file.name}: expected exactly one "
                "'<N> constraints' header."
            )

        if constraint_headers[0] < 1:
            self._fail(
                f"Reaction map validation failed for "
                f"{map_file.name}: constraint count is zero."
            )

        self._info(
            f"Validated modified reaction map: {map_file.name}"
        )