"""Shared-C PA.10 ``T_sh`` screen for the selected rescaled leading path.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected 2026-09-09 reader
(Zenodo 22678406).

The reconstruction requires

    ||ell_i||_{C^0_eta} <= B_0,
    T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf),
    f(eta)=(1+eta^2)^(-1),

and then ``X_sep=110 exp(T_sh)`` must satisfy
``x_sep=X_sep/X_R < exp(-8)``.  For the source flat step used by the
repository, ``||sigma'||_inf=8`` and ``||log f||_inf=log(2)`` on
``eta in [-1,1]``.

Earlier finite-C code made this inequality executable on the old diagnostic
outer scale.  PRs #481/#491 establish that the selected source-rescaled path
must instead use one shared ``C`` in both ``F=phi/C`` and
``X_R=110(C P_*)^10``.  This module performs the same local numerical B_0
screen against that *same* rescaled scale.  It deliberately does not propose a
post-hoc C multiplier: ``ell_i=log(C E_i)`` uses the same C, so changing C also
changes the quantity entering the PA.10 lower bound.  Any future C adjustment
must therefore recompute the coupled upstream path rather than holding B_0
fixed.

The finite eta screen below is numerical evidence for the selected autonomous
realization, not an analytic proof of the source B_0 bound.  Consequently a
PASS would only authorize the selected PA.16 handoff; it would not set
``source_T_sh_lower_bound_verified`` or establish a global leading profile.
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

from .kokuno_rescaled_appendix_b_boundary import X_I
from .kokuno_rescaled_appendix_b_incoming_moments import (
    KokunoRescaledAppendixBIncomingMoments,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-rescaled-inner-join-tsh-certificate-v1"
SIGMA_PRIME_SUP = 8.0
LOG_F_SUP = math.log(2.0)
LOG_X_RESTORE_START = -8.0

_SOURCE_FORMULAS = {
    "B0_bound": "||ell_i||_{C^0_eta} <= B_0",
    "PA10_T_sh": "T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf)",
    "source_f": "f(eta)=(1+eta^2)^(-1); ||log f||_inf=log 2 on [-1,1]",
    "source_step": "source flat sigma; ||sigma'||_inf=8",
    "boundary_log": "ell_i(eta)=log(C E(X_i,eta)); X_i=110",
    "shared_C": "the same fixed C in F=phi/C is used in X_R=110(C P_*)^10",
    "separation": "X_sep=110 exp(T_sh); x_sep=X_sep/X_R<exp(-8)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "PA10_T_sh_formula_executable": True,
    "selected_rescaled_B0_C0_envelope_numerically_checked": True,
    "selected_rescaled_T_sh_lower_bound_numerically_instantiated": True,
    "shared_C_outer_scale_used": True,
    "coupled_C_repair_not_inferred_from_fixed_B0": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
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
class KokunoRescaledInnerJoinTshCertificate:
    """Numerically instantiate PA.10 on the #491 shared-C rescaled path."""

    incoming: KokunoRescaledAppendixBIncomingMoments = field(
        default_factory=KokunoRescaledAppendixBIncomingMoments
    )
    eta_nodes: int = 5
    envelope_relative_padding: float = 0.02
    envelope_absolute_padding: float = 1.0e-6

    def __post_init__(self) -> None:
        if not isinstance(self.incoming, KokunoRescaledAppendixBIncomingMoments):
            raise TypeError("incoming must be a KokunoRescaledAppendixBIncomingMoments")
        if isinstance(self.eta_nodes, bool) or not isinstance(
            self.eta_nodes, (int, np.integer)
        ):
            raise TypeError("eta_nodes must be an integer")
        n = int(self.eta_nodes)
        if not 5 <= n <= 65 or n % 2 == 0:
            raise ValueError("eta_nodes must be an odd integer in [5,65]")
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
        values = np.asarray(
            self.incoming.boundary.boundary_values(eta)["ell_i"], dtype=float
        )
        if values.shape != eta.shape or np.any(~np.isfinite(values)):
            raise RuntimeError("rescaled Appendix-B ell_i evaluation returned invalid values")
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
    def max_T_sh_for_separation(self) -> float:
        return (
            float(self.incoming.binding.log_X_R)
            - math.log(X_I)
            + LOG_X_RESTORE_START
        )

    @property
    def separation_geometry_margin(self) -> float:
        return self.max_T_sh_for_separation - self.selected_T_sh

    @property
    def separation_geometry_feasible(self) -> bool:
        return bool(self.separation_geometry_margin > 0.0)

    @property
    def selected_log_x_sep(self) -> float:
        return float(self.incoming.binding.log_x_sep(self.selected_T_sh))

    def geometry_report(self) -> dict[str, Any]:
        return {
            "selected_log_C": float(self.incoming.binding.selected_log_C),
            "selected_log_X_R": float(self.incoming.binding.log_X_R),
            "selected_log_x_i": float(self.incoming.binding.log_x_i),
            "sigma_prime_sup": SIGMA_PRIME_SUP,
            "log_f_sup": LOG_F_SUP,
            "T_sh_lower_bound_from_selected_numerical_B0": self.T_sh_lower_bound,
            "selected_T_sh": self.selected_T_sh,
            "max_T_sh_for_current_shared_C_scale": self.max_T_sh_for_separation,
            "selected_log_x_sep": self.selected_log_x_sep,
            "required_log_x_sep_upper_bound": LOG_X_RESTORE_START,
            "separation_geometry_margin": self.separation_geometry_margin,
            "separation_geometry_feasible": self.separation_geometry_feasible,
            "fixed_B0_C_multiplier_repair_reported": False,
            "reason_no_fixed_B0_C_repair": (
                "ell_i=log(C E_i) uses the same C; changing C requires recomputing "
                "the coupled Appendix-B boundary/B0 rather than holding B0 fixed"
            ),
        }

    @property
    def selected_pa16_handoff_allowed(self) -> bool:
        return bool(
            self._envelope_cache["disjoint_holdout_within_envelope"]
            and self.separation_geometry_feasible
        )

    def require_selected_pa16_handoff(self) -> None:
        if not self.selected_pa16_handoff_allowed:
            raise ValueError(
                "selected shared-C PA.10 geometry is not admissible; PA.16 handoff remains closed"
            )

    def pa16_inputs_at_eta(self, eta: float) -> dict[str, Any]:
        self.require_selected_pa16_handoff()
        payload = dict(self.incoming.pa16_inputs_at_eta(float(eta)))
        payload["selected_T_sh"] = self.selected_T_sh
        payload["selected_rescaled_T_sh_certificate_sha256"] = self.sha256
        payload["source_T_sh_lower_bound_verified"] = False
        return payload

    def report(self) -> dict[str, Any]:
        envelope = self.envelope_report()
        return {
            "schema": SCHEMA,
            "certificate_sha256": self.sha256,
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
            "selected_pa16_handoff_allowed": self.selected_pa16_handoff_allowed,
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
        return hashlib.sha256(
            _canonical_json(self._unsigned_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(
        cls, payload: dict[str, Any]
    ) -> "KokunoRescaledInnerJoinTshCertificate":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected rescaled inner-join T_sh certificate schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("rescaled T_sh source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("rescaled T_sh truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("rescaled T_sh payload SHA-256 mismatch")
        params = payload.get("parameters")
        deps = payload.get("dependencies")
        if not isinstance(params, dict) or not isinstance(deps, dict):
            raise ValueError("missing rescaled T_sh parameters/dependencies")
        obj = cls(
            incoming=KokunoRescaledAppendixBIncomingMoments.from_payload(
                deps.get("incoming")
            ),
            eta_nodes=int(params["eta_nodes"]),
            envelope_relative_padding=float(params["envelope_relative_padding"]),
            envelope_absolute_padding=float(params["envelope_absolute_padding"]),
        )
        if obj.to_payload() != payload:
            raise ValueError("rescaled T_sh replay changed payload")
        return obj

    def save_json(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoRescaledInnerJoinTshCertificate":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit the selected shared-C PA.10 T_sh numerical certificate"
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--eta-nodes", type=int, default=5)
    parser.add_argument("--axis-quadrature-points", type=int, default=64)
    args = parser.parse_args(argv)
    certificate = KokunoRescaledInnerJoinTshCertificate(
        incoming=KokunoRescaledAppendixBIncomingMoments(
            axis_quadrature_points=args.axis_quadrature_points
        ),
        eta_nodes=args.eta_nodes,
    )
    text = json.dumps(certificate.report(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
