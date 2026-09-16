from __future__ import annotations

import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem


class KineticsEstimationError(Exception):
    """Raised when usable kinetics cannot be obtained."""


@dataclass(frozen=True)
class ResolvedKinetics:
    A: float
    n: float
    Ea_j_mol: float
    T0_K: float
    source: str
    family: str | None = None
    A_units: str | None = None
    estimated: bool = False
    mode: str = "physical_rate"
    probability_scale: float | None = None
    comment: str | None = None


class KineticsEstimator:
    """Resolve RMG kinetics and write LAMMPS bond/react Arrhenius constraints.

    Strict chemistry matching is attempted first. If that fails, the optional
    heuristic path treats monomer/small-molecule chemistry as a surrogate for
    the corresponding local polymer reaction and calibrates the resulting RMG
    physical kinetics to a chosen LAMMPS acceptance probability.
    """

    HEURISTIC_MODE = True
    HEURISTIC_ALLOW_UNHINTED_FALLBACK = True
    HEURISTIC_REFERENCE_PROBABILITY = 0.10

    KB_REAL = 0.00198720425864083       # kcal/(mol K)
    R_SI = 8.31446261815324             # J/(mol K)

    def __init__(self, session: Any | None = None):
        self.session = session
        self.database = None
        self._Molecule = None
        self._Reaction = None
        self._from_rdkit_mol = None

        target_p = getattr(
            session,
            "kinetics_target_probability",
            self.HEURISTIC_REFERENCE_PROBABILITY,
        ) if session is not None else self.HEURISTIC_REFERENCE_PROBABILITY

        try:
            target_p = float(target_p)
        except (TypeError, ValueError):
            target_p = self.HEURISTIC_REFERENCE_PROBABILITY

        if not 0.0 < target_p <= 1.0:
            target_p = self.HEURISTIC_REFERENCE_PROBABILITY

        self.target_probability = target_p

        self._info("Initializing kinetics estimator.")

        try:
            from rmgpy import settings
            from rmgpy.data.kinetics.database import KineticsDatabase
            from rmgpy.molecule.converter import from_rdkit_mol
            from rmgpy.molecule.molecule import Molecule
            from rmgpy.reaction import Reaction
        except ImportError as exc:
            self._fail(
                "RMG-Py is required when kinetics=True.",
                cause=exc,
            )

        self._Molecule = Molecule
        self._Reaction = Reaction
        self._from_rdkit_mol = from_rdkit_mol

        kinetics_path = os.path.join(settings["database.directory"], "kinetics")
        self._info(f"RMG kinetics database path: {kinetics_path}")

        if not os.path.isdir(kinetics_path):
            self._fail(f"RMG kinetics directory does not exist: {kinetics_path}")

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
            self._fail("RMG kinetics database could not be loaded.", cause=exc)

        self._info(
            f"RMG kinetics database loaded successfully with "
            f"{len(self.database.families)} reaction families."
        )

        if self.HEURISTIC_MODE:
            self._warning(
                "HEURISTIC KINETICS FALLBACK IS ENABLED. "
                "Exact chemistry is still attempted first."
            )

    # ------------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------------

    @staticmethod
    def _info(message: str) -> None:
        print(f"[INFO] {message}")

    @staticmethod
    def _warning(message: str) -> None:
        print(f"[WARNING] {message}")

    @staticmethod
    def _fail(message: str, cause: Exception | None = None) -> None:
        print(f"[WARNING] {message}")
        if cause is None:
            raise KineticsEstimationError(message)
        raise KineticsEstimationError(message) from cause

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def estimate_and_write_map(
        self,
        reaction: Any,
        simulation: Any,
        map_file: Path,
        lammps_units: str,
    ) -> bool:
        reaction_id = getattr(reaction, "reaction_id", "unknown")

        try:
            temperature = float(simulation.temperature)
        except (TypeError, ValueError, AttributeError) as exc:
            self._fail(
                f"Reaction {reaction_id}: invalid simulation temperature.",
                cause=exc,
            )

        if not math.isfinite(temperature) or temperature <= 0.0:
            self._fail(
                f"Reaction {reaction_id}: invalid temperature {temperature!r} K."
            )

        self._info(
            f"Estimating kinetics for AutoREACTER reaction {reaction_id} "
            f"at {temperature:.2f} K."
        )

        kinetics = self._resolve_kinetics(reaction)
        self._validate_resolved_kinetics(kinetics, reaction_id)

        A_lmp, n_lmp, Ea_lmp = self._to_lammps_probability(
            kinetics=kinetics,
            reaction=reaction,
            temperature=temperature,
            lammps_units=lammps_units,
            reaction_id=reaction_id,
        )

        probability = self._lammps_probability(
            A_lmp,
            n_lmp,
            Ea_lmp,
            temperature,
            reaction_id,
        )

        self._info(f"Kinetics source: {kinetics.source}")
        if kinetics.family:
            self._info(f"RMG family: {kinetics.family}")
        self._info(f"LAMMPS A  = {A_lmp:.8e}")
        self._info(f"LAMMPS n  = {n_lmp:.8g}")
        self._info(f"LAMMPS Ea = {Ea_lmp:.8g} kcal/mol")
        self._info(
            f"Arrhenius probability at {temperature:.2f} K = {probability:.8e}"
        )

        seed = random.randint(10_000, 2_147_483_646)
        self._write_arrhenius_constraint(
            Path(map_file),
            A_lmp,
            n_lmp,
            Ea_lmp,
            seed,
        )

        self._info(f"Updated reaction map: {Path(map_file).name}")
        return True

    # ------------------------------------------------------------------
    # Resolution order
    # ------------------------------------------------------------------

    def _resolve_kinetics(self, reaction: Any) -> ResolvedKinetics:
        reaction_id = getattr(reaction, "reaction_id", "unknown")

        curated = self._find_curated_kinetics(reaction)
        if curated is not None:
            self._info(f"Reaction {reaction_id}: using curated kinetics.")
            return self._normalize_curated(curated, reaction_id)

        self._info(f"Reaction {reaction_id}: no curated kinetics found.")
        self._info(f"Reaction {reaction_id}: attempting RMG kinetics resolution.")
        return self._estimate_rmg_kinetics(reaction)

    def _find_curated_kinetics(self, reaction: Any) -> dict[str, Any] | None:
        for attr in ("curated_kinetics", "kinetics_data", "user_kinetics"):
            value = getattr(reaction, attr, None)
            if isinstance(value, dict):
                return value

        if self.session is None:
            return None

        rid = getattr(reaction, "reaction_id", None)
        for attr in ("curated_kinetics", "kinetics_registry", "kinetics_data"):
            registry = getattr(self.session, attr, None)
            if not isinstance(registry, dict):
                continue
            if all(k in registry for k in ("A", "n", "Ea")):
                return registry
            for key in (rid, str(rid)):
                if key in registry and isinstance(registry[key], dict):
                    return registry[key]

        return None

    def _normalize_curated(
        self,
        data: dict[str, Any],
        reaction_id: Any,
    ) -> ResolvedKinetics:
        for key in ("A", "n", "Ea", "Ea_units", "source"):
            if key not in data:
                self._fail(
                    f"Reaction {reaction_id}: curated kinetics missing {key!r}."
                )

        try:
            A = float(data["A"])
            n = float(data["n"])
            Ea = float(data["Ea"])
            T0 = float(data.get("T0", 1.0))
        except (TypeError, ValueError) as exc:
            self._fail(
                f"Reaction {reaction_id}: non-numeric curated kinetics.",
                cause=exc,
            )

        mode = str(data.get("mode", "physical_rate")).lower()
        if mode not in {"physical_rate", "lammps_probability"}:
            self._fail(
                f"Reaction {reaction_id}: unsupported kinetics mode {mode!r}."
            )

        probability_scale = data.get("probability_scale")
        if probability_scale is not None:
            probability_scale = float(probability_scale)

        return ResolvedKinetics(
            A=A,
            n=n,
            Ea_j_mol=self._ea_to_jmol(Ea, str(data["Ea_units"]), reaction_id),
            T0_K=T0,
            source=str(data["source"]),
            family=(str(data["family"]) if data.get("family") else None),
            A_units=(str(data["A_units"]) if data.get("A_units") else None),
            estimated=bool(data.get("estimated", False)),
            mode=mode,
            probability_scale=probability_scale,
            comment=(str(data["comment"]) if data.get("comment") else None),
        )

    # ------------------------------------------------------------------
    # RMG search
    # ------------------------------------------------------------------

    def _estimate_rmg_kinetics(self, reaction: Any) -> ResolvedKinetics:
        reaction_id = getattr(reaction, "reaction_id", "unknown")
        reactant_mol = getattr(reaction, "reactant_combined_RDmol", None)
        product_mol = getattr(reaction, "product_combined_RDmol", None)

        if reactant_mol is None or product_mol is None:
            self._fail(
                f"Reaction {reaction_id}: missing combined RDKit reactant/product."
            )

        reactants = self._rdkit_combined_to_rmg(reactant_mol)
        products = self._rdkit_combined_to_rmg(product_mol)

        self._info(
            f"Reaction {reaction_id}: converted {len(reactants)} reactant "
            f"molecule(s) and {len(products)} product molecule(s) to RMG."
        )

        family_hint = getattr(reaction, "rmg_family", None)
        if family_hint is not None:
            family_hint = str(family_hint)
            if family_hint not in self.database.families:
                self._warning(
                    f"Reaction {reaction_id}: RMG family hint {family_hint!r} "
                    "is not loaded; ignoring it."
                )
                family_hint = None

        if family_hint:
            self._info(
                f"Reaction {reaction_id}: restricting RMG search to family "
                f"{family_hint}."
            )
        else:
            self._warning(
                f"Reaction {reaction_id}: no rmg_family hint reached the estimator."
            )

        only_families = [family_hint] if family_hint else None
        target = self._Reaction(reactants=reactants, products=products)

        # Exact product-constrained RMG search.
        self._info(
            f"Reaction {reaction_id}: attempting normal RMG "
            "product-constrained matching."
        )

        try:
            kwargs = {"reactants": reactants, "products": products}
            if only_families:
                kwargs["only_families"] = only_families
            direct = self.database.generate_reactions(**kwargs)
        except Exception as exc:
            self._warning(
                f"Reaction {reaction_id}: product-constrained RMG search failed: {exc}"
            )
            direct = []

        if direct:
            self._info(
                f"Reaction {reaction_id}: RMG returned {len(direct)} direct candidate(s)."
            )
            return self._select_candidate(
                direct,
                target_products=products,
                reaction_id=reaction_id,
                heuristic=False,
                label="direct RMG product match",
            )

        # Reactant-only search.
        self._warning(
            f"Reaction {reaction_id}: no direct product-constrained RMG match."
        )

        try:
            kwargs = {"reactants": reactants}
            if only_families:
                kwargs["only_families"] = only_families
            generated = self.database.generate_reactions(**kwargs)
        except Exception as exc:
            self._fail(
                f"Reaction {reaction_id}: RMG reactant-only generation failed.",
                cause=exc,
            )

        if not generated:
            self._fail(
                f"Reaction {reaction_id}: RMG generated no candidate reactions."
            )

        self._info(
            f"Reaction {reaction_id}: RMG generated {len(generated)} "
            "candidate reaction(s) from the reactants."
        )

        exact = []
        for candidate in generated:
            try:
                if candidate.is_isomorphic(
                    target,
                    either_direction=False,
                    check_identical=False,
                    check_template_rxn_products=True,
                    strict=False,
                ):
                    exact.append(candidate)
            except Exception:
                continue

        if exact:
            self._info(
                f"Reaction {reaction_id}: {len(exact)} reactant-generated "
                "candidate(s) passed RMG reaction isomorphism."
            )
            return self._select_candidate(
                exact,
                target_products=products,
                reaction_id=reaction_id,
                heuristic=False,
                label="RMG reaction-isomorphism fallback",
            )

        if not self.HEURISTIC_MODE:
            self._fail(
                f"Reaction {reaction_id}: no exact RMG transformation found."
            )

        if not family_hint and not self.HEURISTIC_ALLOW_UNHINTED_FALLBACK:
            self._fail(
                f"Reaction {reaction_id}: heuristic fallback requires rmg_family."
            )

        self._warning("=" * 68)
        self._warning("HEURISTIC MONOMER / SMALL-MOLECULE SURROGATE MODE")
        self._warning("=" * 68)
        self._warning(
            "Assumption: monomer-monomer intrinsic kinetics approximate the "
            "same local reaction when one reactant is part of a polymer chain."
        )
        self._warning(
            "Full product-graph equality is being relaxed for this fallback."
        )
        if not family_hint:
            self._warning(
                "No family hint is available, so broad RMG candidates will be "
                "ranked by product formula/count and kinetics availability."
            )

        return self._select_candidate(
            generated,
            target_products=products,
            reaction_id=reaction_id,
            heuristic=True,
            label=(
                f"heuristic surrogate in family {family_hint}"
                if family_hint
                else "heuristic surrogate from broad RMG search"
            ),
        )

    # ------------------------------------------------------------------
    # Candidate kinetics
    # ------------------------------------------------------------------

    def _select_candidate(
        self,
        candidates: list[Any],
        target_products: list[Any],
        reaction_id: Any,
        heuristic: bool,
        label: str,
    ) -> ResolvedKinetics:
        usable: list[tuple[float, int, ResolvedKinetics]] = []

        for index, candidate in enumerate(candidates, start=1):
            family = self._family_label(candidate)
            try:
                kinetics_obj, provenance = self._get_candidate_kinetics(candidate)
                resolved = self._normalize_rmg_kinetics(
                    candidate,
                    kinetics_obj,
                    provenance,
                    reaction_id,
                )
            except KineticsEstimationError as exc:
                self._warning(
                    f"Reaction {reaction_id}: candidate {index} ({family}) "
                    f"has no usable kinetics: {exc}"
                )
                continue

            score = (
                self._heuristic_score(candidate, target_products, resolved)
                if heuristic
                else 0.0
            )
            usable.append((score, index, resolved))

        if not usable:
            self._fail(
                f"Reaction {reaction_id}: RMG candidates were found but none "
                "yielded usable Arrhenius kinetics."
            )

        if heuristic:
            usable.sort(key=lambda item: (item[0], -item[1]), reverse=True)
            score, index, chosen = usable[0]
            self._warning(
                f"Reaction {reaction_id}: heuristic candidate #{index} selected "
                f"with score {score:.1f}; family={chosen.family}."
            )
            if len(usable) > 1:
                self._warning(
                    f"Runner-up: score={usable[1][0]:.1f}, "
                    f"family={usable[1][2].family}."
                )
        else:
            chosen = usable[0][2]
            if len(usable) > 1:
                self._warning(
                    f"Reaction {reaction_id}: multiple exact candidates have "
                    "usable kinetics; using the first one."
                )

        self._info(f"Reaction {reaction_id}: selected kinetics using {label}.")
        self._report_rmg(chosen, reaction_id)
        return chosen

    def _get_candidate_kinetics(self, candidate: Any) -> tuple[Any, str]:
        kinetics = getattr(candidate, "kinetics", None)
        if kinetics is not None:
            return kinetics, "kinetics attached to generated RMG reaction"

        family_label = self._family_label(candidate)
        family = self.database.families.get(family_label)
        if family is None:
            raise KineticsEstimationError(
                f"family {family_label!r} is unavailable"
            )

        template = getattr(candidate, "template", None)
        degeneracy = getattr(candidate, "degeneracy", 1) or 1
        if template is None:
            raise KineticsEstimationError("candidate has no RMG template")

        # Normal RMG family resolver: training/depository/rules as available.
        try:
            result = family.get_kinetics(
                candidate,
                template_labels=template,
                degeneracy=degeneracy,
                estimator="",
                return_all_kinetics=False,
            )
            if result:
                kinetics_obj = result[0] if isinstance(result, tuple) else result
                source = (
                    str(result[1])
                    if isinstance(result, tuple) and len(result) > 1
                    else "RMG family kinetics"
                )
                if kinetics_obj is not None:
                    candidate.kinetics = kinetics_obj
                    return kinetics_obj, source
        except Exception as exc:
            self._warning(
                f"RMG family.get_kinetics failed for {family_label}: {exc}"
            )

        # Final RMG-native fallback: hierarchy/rate rules.
        try:
            result = family.estimate_kinetics_using_rate_rules(
                template,
                degeneracy=degeneracy,
            )
            kinetics_obj = result[0] if isinstance(result, tuple) else result
            if kinetics_obj is not None:
                candidate.kinetics = kinetics_obj
                return kinetics_obj, "RMG hierarchical rate-rule estimate"
        except Exception as exc:
            raise KineticsEstimationError(
                f"rate-rule estimation failed: {exc}"
            ) from exc

        raise KineticsEstimationError("RMG returned no kinetics")

    def _normalize_rmg_kinetics(
        self,
        candidate: Any,
        kinetics: Any,
        provenance: str,
        reaction_id: Any,
    ) -> ResolvedKinetics:
        arrhenius = getattr(kinetics, "arrhenius", None)
        if isinstance(arrhenius, (list, tuple)) and arrhenius:
            self._warning(
                f"Reaction {reaction_id}: MultiArrhenius-like kinetics found; "
                "using the first expression."
            )
            kinetics = arrhenius[0]

        A_q = getattr(kinetics, "A", None)
        n_q = getattr(kinetics, "n", None)
        if A_q is None or n_q is None:
            raise KineticsEstimationError(
                f"{type(kinetics).__name__} does not expose A and n"
            )

        A = self._quantity_si(A_q, "A")
        n = self._quantity_si(n_q, "n")

        Ea_q = getattr(kinetics, "Ea", None)
        estimated = False
        if Ea_q is not None:
            Ea = self._quantity_si(Ea_q, "Ea")
        else:
            E0_q = getattr(kinetics, "E0", None)
            if E0_q is None or not self.HEURISTIC_MODE:
                raise KineticsEstimationError(
                    f"{type(kinetics).__name__} does not expose Ea"
                )
            Ea = self._quantity_si(E0_q, "E0")
            estimated = True
            self._warning(
                f"Reaction {reaction_id}: {type(kinetics).__name__} has E0 but "
                "not Ea; using E0 as the heuristic activation-energy surrogate."
            )

        T0_q = getattr(kinetics, "T0", None)
        T0 = self._quantity_si(T0_q, "T0") if T0_q is not None else 1.0
        units = getattr(A_q, "units", None)
        comment = str(getattr(kinetics, "comment", "") or "").strip()

        lower = (comment + " " + provenance).lower()
        if any(x in lower for x in ("estimate", "rate rule", "template", "average")):
            estimated = True

        return ResolvedKinetics(
            A=A,
            n=n,
            Ea_j_mol=Ea,
            T0_K=T0,
            source=f"RMG {self._family_label(candidate)}: {provenance}",
            family=self._family_label(candidate),
            A_units=(str(units) if units is not None else None),
            estimated=estimated,
            mode="physical_rate",
            comment=(comment or None),
        )

    def _report_rmg(self, kinetics: ResolvedKinetics, reaction_id: Any) -> None:
        self._info(f"Reaction {reaction_id}: RMG family = {kinetics.family}")
        self._info(f"Reaction {reaction_id}: physical A = {kinetics.A:.8e}")
        if kinetics.A_units:
            self._info(f"Reaction {reaction_id}: A units = {kinetics.A_units}")
        self._info(f"Reaction {reaction_id}: n = {kinetics.n:.8g}")
        self._info(
            f"Reaction {reaction_id}: physical Ea = "
            f"{kinetics.Ea_j_mol / 1000.0:.6f} kJ/mol"
        )
        if kinetics.estimated:
            self._warning(
                f"Reaction {reaction_id}: RMG kinetics include an estimated "
                "or hierarchical component."
            )

    # ------------------------------------------------------------------
    # Heuristic candidate ranking
    # ------------------------------------------------------------------

    def _heuristic_score(
        self,
        candidate: Any,
        target_products: list[Any],
        kinetics: ResolvedKinetics,
    ) -> float:
        score = 0.0
        candidate_products = list(getattr(candidate, "products", []) or [])

        if len(candidate_products) == len(target_products):
            score += 25.0

        target_formulas = sorted(self._formula(x) for x in target_products)
        candidate_formulas = sorted(self._entity_formula(x) for x in candidate_products)

        if candidate_formulas == target_formulas and candidate_formulas:
            score += 100.0
        elif self._formula_signature(candidate_formulas) == self._formula_signature(target_formulas):
            score += 35.0

        if kinetics.family == "Diels_alder_addition":
            score += 5.0
        if not kinetics.estimated:
            score += 3.0
        score += 2.0
        return score

    @staticmethod
    def _formula_signature(formulas: list[str]) -> str:
        return "+".join(sorted(formulas))

    def _entity_formula(self, entity: Any) -> str:
        molecules = getattr(entity, "molecule", None)
        if molecules:
            return self._formula(molecules[0])
        return self._formula(entity)

    @staticmethod
    def _formula(molecule: Any) -> str:
        try:
            return str(molecule.get_formula())
        except Exception:
            return ""

    # ------------------------------------------------------------------
    # Physical kinetics -> LAMMPS probability
    # ------------------------------------------------------------------

    def _to_lammps_probability(
        self,
        kinetics: ResolvedKinetics,
        reaction: Any,
        temperature: float,
        lammps_units: str,
        reaction_id: Any,
    ) -> tuple[float, float, float]:
        if not isinstance(lammps_units, str) or lammps_units.lower() != "real":
            self._fail("Automatic kinetics conversion currently supports units real only.")

        Ea_lmp = kinetics.Ea_j_mol / 4184.0
        n_lmp = kinetics.n

        if kinetics.mode == "lammps_probability":
            A_lmp = kinetics.A / math.pow(kinetics.T0_K, kinetics.n)
            return A_lmp, n_lmp, Ea_lmp

        physical_rate = self._physical_rate(kinetics, temperature)
        self._info(
            f"Reaction {reaction_id}: physical surrogate rate at "
            f"{temperature:.2f} K = {physical_rate:.8e}"
        )

        scale = kinetics.probability_scale
        if scale is None:
            scale = getattr(reaction, "kinetics_probability_scale", None)
            if scale is not None:
                scale = float(scale)

        if scale is not None:
            if not math.isfinite(scale) or scale <= 0.0:
                self._fail(
                    f"Reaction {reaction_id}: invalid kinetics_probability_scale."
                )
            A_lmp = kinetics.A / math.pow(kinetics.T0_K, kinetics.n) * scale
            return A_lmp, n_lmp, Ea_lmp

        if not self.HEURISTIC_MODE:
            self._fail(
                f"Reaction {reaction_id}: physical kinetics require a "
                "kinetics_probability_scale."
            )

        # Preserve n and Ea, choose A so P(T_reference) = target probability.
        thermal = math.pow(temperature, n_lmp) * math.exp(
            -Ea_lmp / (self.KB_REAL * temperature)
        )
        if not math.isfinite(thermal) or thermal <= 0.0:
            self._fail(
                f"Reaction {reaction_id}: invalid Arrhenius thermal factor."
            )

        A_lmp = self.target_probability / thermal
        implied_scale = self.target_probability / physical_rate

        self._warning("=" * 68)
        self._warning("HEURISTIC PHYSICAL-RATE -> MD-PROBABILITY CALIBRATION")
        self._warning("=" * 68)
        self._warning(
            f"Reaction {reaction_id}: preserving RMG n={n_lmp:.8g} and "
            f"Ea={Ea_lmp:.8g} kcal/mol."
        )
        self._warning(
            f"Reaction {reaction_id}: choosing A so P={self.target_probability:.3f} "
            f"at {temperature:.2f} K."
        )
        self._warning(
            f"Reaction {reaction_id}: implied scale = {implied_scale:.8e}."
        )
        self._warning(
            "This is a modeling shortcut, not an experimentally validated "
            "condensed-phase rate-to-probability conversion."
        )

        return A_lmp, n_lmp, Ea_lmp

    def _physical_rate(self, k: ResolvedKinetics, temperature: float) -> float:
        rate = (
            k.A
            * math.pow(temperature / k.T0_K, k.n)
            * math.exp(-k.Ea_j_mol / (self.R_SI * temperature))
        )
        if not math.isfinite(rate) or rate <= 0.0:
            raise KineticsEstimationError(
                f"Invalid physical Arrhenius rate at {temperature} K: {rate!r}"
            )
        return rate

    def _lammps_probability(
        self,
        A: float,
        n: float,
        Ea: float,
        temperature: float,
        reaction_id: Any,
    ) -> float:
        probability = (
            A
            * math.pow(temperature, n)
            * math.exp(-Ea / (self.KB_REAL * temperature))
        )
        if not math.isfinite(probability) or probability <= 0.0:
            self._fail(
                f"Reaction {reaction_id}: invalid LAMMPS probability {probability!r}."
            )
        if probability > 1.0 + 1.0e-12:
            self._fail(
                f"Reaction {reaction_id}: LAMMPS probability {probability:.8e} > 1."
            )
        return min(probability, 1.0)

    # ------------------------------------------------------------------
    # RDKit -> RMG
    # ------------------------------------------------------------------

    def _rdkit_combined_to_rmg(self, combined: Chem.Mol) -> list[Any]:
        try:
            fragments = Chem.GetMolFrags(
                combined,
                asMols=True,
                sanitizeFrags=True,
            )
        except Exception as exc:
            self._fail("RDKit could not split combined molecule.", cause=exc)

        output = []
        for i, fragment in enumerate(fragments, start=1):
            molecule = self._Molecule()
            try:
                self._from_rdkit_mol(molecule, fragment)
            except Exception as exc:
                self._fail(f"RMG conversion failed for fragment {i}.", cause=exc)
            if not getattr(molecule, "vertices", None):
                self._fail(f"RMG conversion produced empty fragment {i}.")
            output.append(molecule)

        return output

    # ------------------------------------------------------------------
    # Validation / utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _quantity_si(quantity: Any, name: str) -> float:
        value = getattr(quantity, "value_si", quantity)
        try:
            value = float(value)
        except (TypeError, ValueError) as exc:
            raise KineticsEstimationError(
                f"RMG parameter {name} is not numeric: {value!r}"
            ) from exc
        if not math.isfinite(value):
            raise KineticsEstimationError(f"RMG parameter {name} is not finite")
        return value

    @staticmethod
    def _family_label(reaction: Any) -> str:
        family = getattr(reaction, "family", None)
        if family is None:
            return "unknown"
        if isinstance(family, str):
            return family
        return str(getattr(family, "label", family))

    def _validate_resolved_kinetics(
        self,
        kinetics: ResolvedKinetics,
        reaction_id: Any,
    ) -> None:
        for name, value in {
            "A": kinetics.A,
            "n": kinetics.n,
            "Ea": kinetics.Ea_j_mol,
            "T0": kinetics.T0_K,
        }.items():
            if not math.isfinite(value):
                self._fail(
                    f"Reaction {reaction_id}: kinetics parameter {name} is not finite."
                )
        if kinetics.A <= 0.0:
            self._fail(f"Reaction {reaction_id}: A must be positive.")
        if kinetics.T0_K <= 0.0:
            self._fail(f"Reaction {reaction_id}: T0 must be positive.")

    def _ea_to_jmol(self, value: float, units: str, reaction_id: Any) -> float:
        key = units.strip().lower().replace(" ", "")
        factors = {
            "j/mol": 1.0,
            "kj/mol": 1000.0,
            "cal/mol": 4.184,
            "kcal/mol": 4184.0,
            "jmol^-1": 1.0,
            "kjmol^-1": 1000.0,
            "calmol^-1": 4.184,
            "kcalmol^-1": 4184.0,
        }
        if key not in factors:
            self._fail(
                f"Reaction {reaction_id}: unsupported Ea units {units!r}."
            )
        return value * factors[key]

    # ------------------------------------------------------------------
    # Map patching
    # ------------------------------------------------------------------

    def _write_arrhenius_constraint(
        self,
        map_file: Path,
        A: float,
        n: float,
        Ea: float,
        seed: int,
    ) -> None:
        if not map_file.is_file():
            self._fail(f"Reaction map file does not exist: {map_file}")

        lines = map_file.read_text().splitlines()
        if not lines:
            self._fail(f"Reaction map is empty: {map_file}")

        arr_line = f"arrhenius {A:.8e} {n:.8g} {Ea:.8g} {seed}"
        existing = [
            i for i, line in enumerate(lines)
            if line.strip().lower().startswith("arrhenius ")
        ]

        if len(existing) > 1:
            self._fail(
                f"{map_file.name} contains multiple Arrhenius constraints."
            )

        if existing:
            lines[existing[0]] = arr_line
        else:
            header_idx, count = self._constraints_header(lines)
            section_idx = self._section(lines, "Constraints")

            if count is not None and count > 0 and section_idx is None:
                self._fail(
                    f"{map_file.name} declares constraints but has no Constraints section."
                )

            if header_idx is None:
                initiator_idx = self._section(lines, "InitiatorIDs")
                if initiator_idx is None:
                    self._fail(
                        f"{map_file.name}: InitiatorIDs section not found."
                    )
                insert_at = initiator_idx
                while insert_at > 0 and not lines[insert_at - 1].strip():
                    insert_at -= 1
                lines.insert(insert_at, "1 constraints")
            else:
                lines[header_idx] = f"{count + 1} constraints"

            section_idx = self._section(lines, "Constraints")
            if section_idx is None:
                if lines and lines[-1].strip():
                    lines.append("")
                lines.extend(["Constraints", "", arr_line])
            else:
                lines.insert(self._constraints_end(lines, section_idx), arr_line)

        self._validate_map(lines, map_file.name)

        tmp = map_file.with_name(map_file.name + ".tmp")
        try:
            tmp.write_text("\n".join(lines) + "\n")
            tmp.replace(map_file)
        except Exception as exc:
            try:
                if tmp.exists():
                    tmp.unlink()
            except OSError:
                pass
            self._fail(f"Could not update {map_file} atomically.", cause=exc)

        self._info(f"Validated modified reaction map: {map_file.name}")

    @staticmethod
    def _constraints_header(lines: list[str]) -> tuple[int | None, int | None]:
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit() and parts[1].lower() == "constraints":
                return i, int(parts[0])
        return None, None

    @staticmethod
    def _section(lines: list[str], name: str) -> int | None:
        for i, line in enumerate(lines):
            if line.strip() == name:
                return i
        return None

    @staticmethod
    def _constraints_end(lines: list[str], start: int) -> int:
        sections = {
            "InitiatorIDs",
            "Equivalences",
            "EdgeIDs",
            "Wildcards",
            "DeleteIDs",
            "CreateIDs",
            "ChiralIDs",
            "Constraints",
        }
        for i in range(start + 1, len(lines)):
            if lines[i].strip() in sections:
                return i
        if lines and lines[-1].strip():
            lines.append("")
        return len(lines)

    def _validate_map(self, lines: list[str], name: str) -> None:
        arr = [
            line for line in lines
            if line.strip().lower().startswith("arrhenius ")
        ]
        if len(arr) != 1:
            self._fail(
                f"Map validation failed for {name}: found {len(arr)} Arrhenius lines."
            )
        if self._section(lines, "Constraints") is None:
            self._fail(
                f"Map validation failed for {name}: Constraints section missing."
            )
        _, count = self._constraints_header(lines)
        if count is None or count < 1:
            self._fail(
                f"Map validation failed for {name}: invalid constraints count."
            )
