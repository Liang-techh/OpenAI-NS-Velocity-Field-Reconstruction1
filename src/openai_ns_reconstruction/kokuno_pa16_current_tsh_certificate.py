"""Current-lineage numerical PA.10 ``T_sh`` certificate for the PA.16 join.

Pinned public provenance is the corrected 2026-09-09 KokunoYumeto
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.
The public Appendix-B / PA.10 continuation first has an analytic bound

    ||ell_i||_{C^0_eta} <= B_0

and then requires

    T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf),
    f(eta)=(1+eta^2)^(-1).

For the source flat step already used by this repository,
``||sigma'||_inf=8``.  On the source eta interval ``[-1,1]``,
``||log f||_inf=log(2)``.  The angular transition must also finish before the
axial-restoration interval,

    log x_sep = log(X_i/X_R) + T_sh < -8,
    X_i=110,  X_R=110(C P_*)^10.

This module instantiates those formulas on the *current* Agent-1 lineage from
``KokunoPA10ActualXiPrefixMoments``.  It evaluates the current Xi boundary
``ell_i(eta)`` on Chebyshev-Lobatto construction nodes and on disjoint midpoint
nodes, then builds a padded numerical envelope.  That finite sampling is useful
candidate-side evidence but is intentionally not promoted to the source's
analytic ``B_0`` theorem; consequently ``source_T_sh_lower_bound_verified``
stays false.

The class can build the existing source-form PA.16 join only when the selected
numerical lower bound fits the public separation geometry.  It does not claim
that PA.16 has solved, that the global inner/outer join is complete, or that a
matched pressure/forcing/full NS residual exists.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from functools import cached_property
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_inner_join_exit import KokunoInnerJoinExit, InnerJoinSolveResult
from .kokuno_pa10_actual_xi_prefix_moments import KokunoPA10ActualXiPrefixMoments
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)


SCHEMA = "kokuno-pa16-current-tsh-certificate-v1"
SOURCE_ETA_INTERVAL = (-1.0, 1.0)
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

_NUMERICAL_REALIZATION = {
    "B0_envelope": (
        "finite Chebyshev-Lobatto construction nodes on [-1,1] plus disjoint "
        "midpoint holdout nodes, with fixed relative/absolute padding"
    ),
    "upstream": (
        "current candidate-side #947 Xi ell_i/G_i/five-moment data; autonomous "
        "pressure, kappa_0 and outer scale remain explicitly candidate-side"
    ),
    "T_sh_selection": (
        "nextafter of the public lower-bound formula evaluated on the selected "
        "numerical B0 envelope; not a recovered hidden source value"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_lineage_Xi_data_bound": True,
    "PA10_T_sh_formula_executable": True,
    "current_source_eta_interval_covered": True,
    "selected_B0_C0_envelope_numerically_checked": True,
    "selected_T_sh_lower_bound_numerically_instantiated": True,
    "selected_join_route_executable_when_geometry_allows": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_prepared_appendixA_pressure_stress_materialized": False,
    "source_admitted_kappa0_materialized": False,
    "source_prepared_upstream_five_moment_discrepancy_materialized": False,
    "five_moment_repair_applied": False,
    "inner_to_outer_join_completed": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _chebyshev_lobatto_nodes(count: int) -> np.ndarray:
    k = np.arange(count, dtype=float)
    return np.cos(np.pi * k / float(count - 1))[::-1]


@dataclass(frozen=True)
class KokunoPA16CurrentTshCertificate:
    """Numerically instantiate the public PA.10 ``T_sh`` bound on #947 data."""

    moments: KokunoPA10ActualXiPrefixMoments = field(
        default_factory=KokunoPA10ActualXiPrefixMoments,
        repr=False,
        compare=False,
    )
    eta_nodes: int = 17
    envelope_relative_padding: float = 0.02
    envelope_absolute_padding: float = 1.0e-6

    def __post_init__(self) -> None:
        if not isinstance(self.moments, KokunoPA10ActualXiPrefixMoments):
            raise TypeError("moments must be KokunoPA10ActualXiPrefixMoments")
        if isinstance(self.eta_nodes, bool) or not isinstance(
            self.eta_nodes, (int, np.integer)
        ):
            raise TypeError("eta_nodes must be an integer")
        n = int(self.eta_nodes)
        if n < 9 or n > 65 or n % 2 == 0:
            raise ValueError("eta_nodes must be an odd integer in [9,65]")
        rel = float(self.envelope_relative_padding)
        absolute = float(self.envelope_absolute_padding)
        if not math.isfinite(rel) or not 0.0 < rel <= 0.2:
            raise ValueError("envelope_relative_padding must lie in (0,0.2]")
        if not math.isfinite(absolute) or not 0.0 < absolute <= 1.0e-2:
            raise ValueError("envelope_absolute_padding must lie in (0,1e-2]")
        lo, hi = self.moments.eta_interval
        if lo > SOURCE_ETA_INTERVAL[0] or hi < SOURCE_ETA_INTERVAL[1]:
            raise ValueError(
                "current Xi evaluator must cover the full source eta interval [-1,1]"
            )
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
        eta = np.asarray(eta, dtype=float)
        values = np.asarray(self.moments.bridge.handoff_at_Xi(eta)["ell_i"], dtype=float)
        if values.shape != eta.shape or np.any(~np.isfinite(values)):
            raise RuntimeError("current Xi ell_i evaluator returned invalid values")
        return values

    @cached_property
    def _envelope_cache(self) -> dict[str, Any]:
        construction = self._ell_values(self.construction_eta)
        holdout = self._ell_values(self.holdout_eta)
        construction_max = float(np.max(np.abs(construction)))
        holdout_max = float(np.max(np.abs(holdout)))
        B0 = (
            construction_max * (1.0 + self.envelope_relative_padding)
            + self.envelope_absolute_padding
        )
        return {
            "eta_interval": list(SOURCE_ETA_INTERVAL),
            "eta_construction_nodes": self.eta_nodes,
            "eta_holdout_nodes": int(self.holdout_eta.size),
            "construction_max_abs_ell_i": construction_max,
            "holdout_max_abs_ell_i": holdout_max,
            "selected_numerical_B0_envelope": B0,
            "envelope_relative_padding": self.envelope_relative_padding,
            "envelope_absolute_padding": self.envelope_absolute_padding,
            "disjoint_holdout_within_envelope": bool(holdout_max <= B0),
            "analytic_source_B0_bound_proved": False,
        }

    def envelope_report(self) -> dict[str, Any]:
        return dict(self._envelope_cache)

    @property
    def selected_B0(self) -> float:
        report = self._envelope_cache
        if not report["disjoint_holdout_within_envelope"]:
            raise RuntimeError("selected numerical B0 envelope failed disjoint eta holdout")
        return float(report["selected_numerical_B0_envelope"])

    @property
    def T_sh_lower_bound(self) -> float:
        return 20.0 * SIGMA_PRIME_SUP * (self.selected_B0 + LOG_F_SUP)

    @property
    def selected_T_sh(self) -> float:
        return math.nextafter(self.T_sh_lower_bound, math.inf)

    @property
    def log_x_i(self) -> float:
        return math.log(self.moments.X_i) - float(self.moments.outer_schedule.log_X_R)

    @property
    def max_T_sh_for_separation(self) -> float:
        # log(X_i/X_R)+T_sh < -8.
        return LOG_X_RESTORE_START - self.log_x_i

    @property
    def separation_geometry_feasible(self) -> bool:
        return self.selected_T_sh < self.max_T_sh_for_separation

    @property
    def geometry_margin(self) -> float:
        return self.max_T_sh_for_separation - self.selected_T_sh

    @property
    def required_delta_log_X_R(self) -> float:
        deficit = self.selected_T_sh - self.max_T_sh_for_separation
        return max(0.0, deficit + (1.0e-9 if deficit >= 0.0 else 0.0))

    @property
    def required_log_C_multiplier_if_only_X_R_is_enlarged(self) -> float:
        # X_R = 110 (C P_*)^10, holding P_* fixed.
        return self.required_delta_log_X_R / 10.0

    @property
    def required_log10_C_multiplier_if_only_X_R_is_enlarged(self) -> float:
        return self.required_log_C_multiplier_if_only_X_R_is_enlarged / math.log(10.0)

    def geometry_report(self) -> dict[str, Any]:
        return {
            "X_i": self.moments.X_i,
            "C": self.moments.C,
            "log_P_star": float(self.moments.outer_schedule.log_P_star),
            "log_X_R": float(self.moments.outer_schedule.log_X_R),
            "log_x_i": self.log_x_i,
            "sigma_prime_sup": SIGMA_PRIME_SUP,
            "log_f_sup": LOG_F_SUP,
            "T_sh_lower_bound_from_selected_numerical_B0": self.T_sh_lower_bound,
            "selected_T_sh": self.selected_T_sh,
            "selected_log_x_sep": self.log_x_i + self.selected_T_sh,
            "required_log_x_sep_upper_bound": LOG_X_RESTORE_START,
            "max_T_sh_for_current_outer_schedule": self.max_T_sh_for_separation,
            "geometry_margin": self.geometry_margin,
            "separation_geometry_feasible": self.separation_geometry_feasible,
            "required_delta_log_X_R": self.required_delta_log_X_R,
            "required_log_C_multiplier_if_only_X_R_is_enlarged": (
                self.required_log_C_multiplier_if_only_X_R_is_enlarged
            ),
            "required_log10_C_multiplier_if_only_X_R_is_enlarged": (
                self.required_log10_C_multiplier_if_only_X_R_is_enlarged
            ),
        }

    def build_selected_inner_join(self, *, quadrature_points: int = 96) -> KokunoInnerJoinExit:
        if not self.separation_geometry_feasible:
            raise ValueError(
                "selected current-lineage PA.10 T_sh does not fit before log x=-8 "
                "under the current outer schedule"
            )
        return KokunoInnerJoinExit(
            T_sh=self.selected_T_sh,
            outer_schedule=self.moments.outer_schedule,
            quadrature_points=quadrature_points,
        )

    def solve_selected_at_eta(
        self,
        eta: float,
        *,
        prefix_order: int | None = None,
        quadrature_points: int = 96,
        absolute_tolerance: float = 2.0e-11,
    ) -> InnerJoinSolveResult:
        """Run the existing PA.16 repair only after this numerical geometry gate.

        This method is candidate-side and does not promote the source analytic
        ``B_0``/``T_sh`` theorem.  The returned solve is likewise not a global
        Navier-Stokes validation.
        """
        join = self.build_selected_inner_join(quadrature_points=quadrature_points)
        data = self.moments.pa16_input_at_eta(float(eta), order=prefix_order)
        return join.solve_at_eta(
            eta=float(eta),
            ell_i=float(data["ell_i"]),
            G_i=float(data["G_i"]),
            incoming_scaled_discrepancy=np.asarray(
                data["incoming_scaled_discrepancy"], dtype=float
            ),
            absolute_tolerance=float(absolute_tolerance),
        )

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
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
            "numerical_realization": dict(_NUMERICAL_REALIZATION),
            "current_parent_semantic_sha256": self.moments.semantic_sha256,
            "selected_B0": envelope,
            "geometry": self.geometry_report(),
            "selected_route_ready": bool(
                envelope["disjoint_holdout_within_envelope"]
                and self.separation_geometry_feasible
            ),
            "source_route_ready": False,
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "moments": self.moments.configuration(),
            "eta_nodes": self.eta_nodes,
            "envelope_relative_padding": self.envelope_relative_padding,
            "envelope_absolute_padding": self.envelope_absolute_padding,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentTshCertificate":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current T_sh certificate schema")
        return cls(
            moments=KokunoPA10ActualXiPrefixMoments.from_configuration(
                payload.get("moments")
            ),
            eta_nodes=payload.get("eta_nodes"),
            envelope_relative_padding=payload.get("envelope_relative_padding"),
            envelope_absolute_padding=payload.get("envelope_absolute_padding"),
        )

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "truth_boundary": _TRUTH_BOUNDARY,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()

    def save_configuration(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.configuration(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA16CurrentTshCertificate":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta-nodes", type=int, default=17)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    certificate = KokunoPA16CurrentTshCertificate(eta_nodes=args.eta_nodes)
    text = json.dumps(certificate.report(), indent=2, sort_keys=True, allow_nan=False)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
