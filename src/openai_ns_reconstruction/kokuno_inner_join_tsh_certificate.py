"""Numerical PA.10 T_sh certificate for the selected Appendix-B join.

Pinned provenance is KokunoYumeto/yang-mills-interacting-workbench at commit
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.
The source first bounds ``||ell_i||_{C^0_eta} <= B_0`` and then requires

    T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf),
    f=(1+eta^2)^(-1),

with ``X_sep=110 exp(T_sh)`` and ``x_sep=X_sep/X_R < exp(-8)`` before the
axial-restoration interval.  For the source flat step used in this repository,
``||sigma'||_inf=8`` and ``||log f||_inf=log(2)``.

This module makes that inequality executable for the *selected autonomous*
Appendix-B realization.  It estimates a conservative C0 envelope for ell_i on
Chebyshev-Lobatto construction nodes and checks it on disjoint midpoint nodes.
That finite check is deliberately labelled numerical rather than a proof of the
source C0 bound.  Consequently it can safely choose a source-form T_sh and
route the selected data into PA.16 when the radial geometry permits, but it does
not promote ``source_T_sh_lower_bound_verified`` or the global leading field.

The separate source axis-pressure datum exposed by #438 also remains a hard
blocker.  No hidden Kokuno/OpenAI parameter is inferred.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_appendix_b_boundary import X_I
from .kokuno_appendix_b_incoming_moments import KokunoAppendixBIncomingMoments
from .kokuno_inner_join_exit import KokunoInnerJoinExit, InnerJoinSolveResult

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-inner-join-tsh-certificate-v1"
SIGMA_PRIME_SUP = 8.0
LOG_F_SUP = math.log(2.0)
LOG_X_RESTORE_START = -8.0

_SOURCE_FORMULAS = {
    "B0_bound": "||ell_i||_{C^0_eta} <= B_0",
    "PA10_T_sh": "T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf)",
    "source_f": "f(eta)=(1+eta^2)^(-1); ||log f||_inf=log 2 on [-1,1]",
    "source_step": "source flat sigma; ||sigma'||_inf=8",
    "separation": "X_sep=110*exp(T_sh); x_sep=X_sep/X_R<exp(-8)",
    "outer_scale": "X_R=110*(C*P_*)^10",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "PA10_T_sh_formula_executable": True,
    "selected_B0_C0_envelope_numerically_checked": True,
    "selected_T_sh_lower_bound_numerically_instantiated": True,
    "selected_inner_join_route_executable_when_geometry_allows": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_outer_pressure_datum_bound": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _chebyshev_lobatto_nodes(count: int) -> np.ndarray:
    k = np.arange(count, dtype=float)
    return np.cos(np.pi * k / float(count - 1))[::-1]


@dataclass(frozen=True)
class KokunoInnerJoinTshCertificate:
    """Selected-realization numerical instantiation of the public PA.10 bound."""

    incoming: KokunoAppendixBIncomingMoments = field(
        default_factory=KokunoAppendixBIncomingMoments
    )
    eta_nodes: int = 17
    envelope_relative_padding: float = 0.02
    envelope_absolute_padding: float = 1.0e-6

    def __post_init__(self) -> None:
        if not isinstance(self.incoming, KokunoAppendixBIncomingMoments):
            raise TypeError("incoming must be a KokunoAppendixBIncomingMoments")
        if isinstance(self.eta_nodes, bool) or not isinstance(
            self.eta_nodes, (int, np.integer)
        ):
            raise TypeError("eta_nodes must be an integer")
        n = int(self.eta_nodes)
        if not 9 <= n <= 65 or n % 2 == 0:
            raise ValueError("eta_nodes must be an odd integer in [9,65]")
        rel = float(self.envelope_relative_padding)
        absolute = float(self.envelope_absolute_padding)
        if not math.isfinite(rel) or not 0.0 < rel <= 0.2:
            raise ValueError("envelope_relative_padding must lie in (0,0.2]")
        if not math.isfinite(absolute) or not 0.0 < absolute <= 1.0e-2:
            raise ValueError("envelope_absolute_padding must lie in (0,1e-2]")
        object.__setattr__(self, "eta_nodes", n)
        object.__setattr__(self, "envelope_relative_padding", rel)
        object.__setattr__(self, "envelope_absolute_padding", absolute)

    @property
    def construction_eta(self) -> np.ndarray:
        return _chebyshev_lobatto_nodes(self.eta_nodes)

    @property
    def holdout_eta(self) -> np.ndarray:
        nodes = self.construction_eta
        return 0.5 * (nodes[:-1] + nodes[1:])

    def _ell_values(self, eta: np.ndarray) -> np.ndarray:
        values = np.asarray(self.incoming.boundary.boundary_values(eta)["ell_i"], dtype=float)
        if values.shape != eta.shape or np.any(~np.isfinite(values)):
            raise RuntimeError("Appendix-B ell_i evaluation returned invalid values")
        return values

    @cached_property
    def _envelope_cache(self) -> dict[str, Any]:
        # Appendix-B boundary evaluation integrates a nontrivial ODE for every
        # eta.  A certificate instance is immutable, so evaluating this finite
        # construction/holdout set once is both exact replay and a substantial
        # CI/runtime improvement; it does not change the numerical contract.
        construction = self._ell_values(self.construction_eta)
        holdout = self._ell_values(self.holdout_eta)
        construction_max = float(np.max(np.abs(construction)))
        holdout_max = float(np.max(np.abs(holdout)))
        B0 = construction_max * (1.0 + self.envelope_relative_padding) + self.envelope_absolute_padding
        passed = bool(holdout_max <= B0)
        return {
            "eta_construction_nodes": self.eta_nodes,
            "eta_holdout_nodes": int(self.holdout_eta.size),
            "construction_max_abs_ell_i": construction_max,
            "holdout_max_abs_ell_i": holdout_max,
            "selected_numerical_B0_envelope": B0,
            "envelope_relative_padding": self.envelope_relative_padding,
            "envelope_absolute_padding": self.envelope_absolute_padding,
            "disjoint_holdout_within_envelope": passed,
            "analytic_source_B0_bound_proved": False,
        }

    def envelope_report(self) -> dict[str, Any]:
        return dict(self._envelope_cache)

    @property
    def selected_B0(self) -> float:
        report = self._envelope_cache
        if not report["disjoint_holdout_within_envelope"]:
            raise RuntimeError("selected numerical B0 envelope failed its disjoint eta holdout")
        return float(report["selected_numerical_B0_envelope"])

    @property
    def T_sh_lower_bound(self) -> float:
        return 20.0 * SIGMA_PRIME_SUP * (self.selected_B0 + LOG_F_SUP)

    @property
    def selected_T_sh(self) -> float:
        # The source permits equality.  nextafter makes the floating comparison
        # robust while changing no mathematical scale.
        return math.nextafter(self.T_sh_lower_bound, math.inf)

    @property
    def max_T_sh_for_separation(self) -> float:
        # log x_sep = log(X_i/X_R)+T_sh must be strictly below -8.
        return float(self.incoming.outer_schedule.log_X_R) - math.log(X_I) + LOG_X_RESTORE_START

    @property
    def separation_geometry_feasible(self) -> bool:
        return self.selected_T_sh < self.max_T_sh_for_separation

    @property
    def required_delta_log_X_R(self) -> float:
        needed = self.selected_T_sh - self.max_T_sh_for_separation
        return max(0.0, needed + (1.0e-9 if needed >= 0.0 else 0.0))

    @property
    def required_log_C_multiplier_if_only_X_R_is_enlarged(self) -> float:
        # X_R=110(CP_*)^10, with P_* fixed in this diagnostic.
        return self.required_delta_log_X_R / 10.0

    @property
    def required_log10_C_multiplier_if_only_X_R_is_enlarged(self) -> float:
        return self.required_log_C_multiplier_if_only_X_R_is_enlarged / math.log(10.0)

    @property
    def required_C_multiplier_float64_saturated(self) -> bool:
        return self.required_log_C_multiplier_if_only_X_R_is_enlarged > math.log(
            np.finfo(float).max
        )

    @property
    def required_C_multiplier_if_only_X_R_is_enlarged(self) -> float:
        """Finite float64 compatibility view; exact size is reported in log space.

        The selected numerical PA.10 scale can require a C multiplier far beyond
        float64.  Saturating this compatibility property avoids an overflow while
        ``required_log_C_multiplier_if_only_X_R_is_enlarged`` and its log10 form
        retain the actual finite diagnostic scale without pretending a clipped
        float is the exact multiplier.
        """
        log_multiplier = self.required_log_C_multiplier_if_only_X_R_is_enlarged
        if self.required_C_multiplier_float64_saturated:
            return float(np.finfo(float).max)
        return math.exp(log_multiplier)

    def geometry_report(self) -> dict[str, Any]:
        log_x_i = math.log(X_I) - float(self.incoming.outer_schedule.log_X_R)
        return {
            "sigma_prime_sup": SIGMA_PRIME_SUP,
            "log_f_sup": LOG_F_SUP,
            "T_sh_lower_bound_from_selected_numerical_B0": self.T_sh_lower_bound,
            "selected_T_sh": self.selected_T_sh,
            "log_x_i": log_x_i,
            "selected_log_x_sep": log_x_i + self.selected_T_sh,
            "required_log_x_sep_upper_bound": LOG_X_RESTORE_START,
            "max_T_sh_for_current_outer_schedule": self.max_T_sh_for_separation,
            "separation_geometry_feasible": self.separation_geometry_feasible,
            "required_delta_log_X_R": self.required_delta_log_X_R,
            "required_log_C_multiplier_if_only_X_R_is_enlarged": self.required_log_C_multiplier_if_only_X_R_is_enlarged,
            "required_log10_C_multiplier_if_only_X_R_is_enlarged": self.required_log10_C_multiplier_if_only_X_R_is_enlarged,
            "required_C_multiplier_if_only_X_R_is_enlarged_float64": self.required_C_multiplier_if_only_X_R_is_enlarged,
            "required_C_multiplier_float64_saturated": self.required_C_multiplier_float64_saturated,
        }

    def build_selected_inner_join(self, *, quadrature_points: int = 96) -> KokunoInnerJoinExit:
        if not self.separation_geometry_feasible:
            raise ValueError(
                "selected PA.10 T_sh does not fit before log x=-8 under the current outer schedule"
            )
        return KokunoInnerJoinExit(
            T_sh=self.selected_T_sh,
            outer_schedule=self.incoming.outer_schedule,
            quadrature_points=quadrature_points,
        )

    def solve_selected_at_eta(
        self,
        eta: float,
        *,
        quadrature_points: int = 96,
        absolute_tolerance: float = 2.0e-11,
    ) -> InnerJoinSolveResult:
        join = self.build_selected_inner_join(quadrature_points=quadrature_points)
        return self.incoming.solve_inner_join_at_eta(
            join,
            eta=float(eta),
            absolute_tolerance=float(absolute_tolerance),
        )

    def report(self) -> dict[str, Any]:
        pressure = self.incoming.pressure_datum_report()
        envelope = self.envelope_report()
        return {
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "selected_B0": envelope,
            "geometry": self.geometry_report(),
            "pressure_datum": pressure,
            "selected_route_ready": bool(
                envelope["disjoint_holdout_within_envelope"]
                and self.separation_geometry_feasible
            ),
            "source_route_ready": False,
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "parameters": {
                "eta_nodes": self.eta_nodes,
                "envelope_relative_padding": self.envelope_relative_padding,
                "envelope_absolute_padding": self.envelope_absolute_padding,
            },
            "dependencies": {"incoming": self.incoming.to_payload()},
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self._unsigned_payload()).encode()).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoInnerJoinTshCertificate":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected inner-join T_sh certificate schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("inner-join T_sh source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("inner-join T_sh truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode()).hexdigest()
        if claimed != expected:
            raise ValueError("inner-join T_sh payload SHA mismatch")
        params = payload.get("parameters", {})
        dep = payload.get("dependencies", {})
        obj = cls(
            incoming=KokunoAppendixBIncomingMoments.from_payload(dep.get("incoming")),
            eta_nodes=params.get("eta_nodes"),
            envelope_relative_padding=params.get("envelope_relative_padding"),
            envelope_absolute_padding=params.get("envelope_absolute_padding"),
        )
        if obj.to_payload() != payload:
            raise ValueError("inner-join T_sh payload does not replay exactly")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoInnerJoinTshCertificate":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta-nodes", type=int, default=17)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    obj = KokunoInnerJoinTshCertificate(eta_nodes=args.eta_nodes)
    report = obj.report()
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
