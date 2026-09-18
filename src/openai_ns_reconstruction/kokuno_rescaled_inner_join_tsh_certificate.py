"""Shared-C PA.10 obstruction screen for the selected rescaled leading path.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected 2026-09-09 reader
(Zenodo 22678406).

The public reconstruction requires

    ||ell_i||_{C^0_eta} <= B_0,
    T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf),
    f(eta)=(1+eta^2)^(-1),

and then ``X_sep=110 exp(T_sh)`` must satisfy
``x_sep=X_sep/X_R < exp(-8)``.  For the source flat step used by the
repository, ``||sigma'||_inf=8`` and ``||log f||_inf=log(2)`` on
``eta in [-1,1]``.

PRs #481/#491 establish that the selected source-rescaled path must use one
shared ``C`` in both ``F=phi/C`` and ``X_R=110(C P_*)^10``.  A first version
of this module sampled an eta grid in an attempt to build a numerical upper
envelope for ``B_0``.  That was unnecessarily expensive and, more
importantly, stronger than is needed to reject the current selected scale.

For any single eta,

    B_0 >= ||ell_i||_inf >= |ell_i(eta)|.

Therefore one point gives a rigorous *necessary* lower bound for every
admissible ``T_sh``.  We evaluate that witness at the stationary real-axis
phase maximum used by the selected normalization.  If even

    20 ||sigma'||_inf (|ell_i(eta_w)| + ||log f||_inf)

is not strictly below the geometric ceiling
``log(X_R)-log(110)-8``, then no unknown larger source ``B_0`` can rescue the
current selected path.  This one-point obstruction is both cheaper and more
logically conservative than pretending a finite eta sample proves the source
``C^0`` bound.

This module still does not prove the source analytic ``B_0`` estimate, recover
hidden parameters, complete PA.16, match global pressure, validate the PDE, or
identify the OpenAI field.  It also deliberately does not propose a post-hoc
``C`` multiplier: ``ell_i=log(C E_i)`` belongs to the same coupled upstream
path, so any new ``C`` must be propagated from Appendix B again.
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
SCHEMA = "kokuno-rescaled-inner-join-tsh-certificate-v2"
SIGMA_PRIME_SUP = 8.0
LOG_F_SUP = math.log(2.0)
LOG_X_RESTORE_START = -8.0
WITNESS_POLICY = "selected-real-axis-phase-stationary-eta"

_SOURCE_FORMULAS = {
    "B0_bound": "||ell_i||_{C^0_eta} <= B_0",
    "pointwise_B0_lower_bound": "B_0 >= ||ell_i||_inf >= |ell_i(eta_w)|",
    "PA10_T_sh": "T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf)",
    "pointwise_T_sh_lower_bound": (
        "T_sh >= 20 ||sigma'||_inf (|ell_i(eta_w)| + ||log f||_inf)"
    ),
    "source_f": "f(eta)=(1+eta^2)^(-1); ||log f||_inf=log 2 on [-1,1]",
    "source_step": "source flat sigma; ||sigma'||_inf=8",
    "boundary_log": "ell_i(eta)=log(C E(X_i,eta)); X_i=110",
    "shared_C": "the same fixed C in F=phi/C is used in X_R=110(C P_*)^10",
    "separation": "X_sep=110 exp(T_sh); x_sep=X_sep/X_R<exp(-8)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "PA10_T_sh_formula_executable": True,
    "selected_rescaled_pointwise_B0_lower_bound_executable": True,
    "selected_rescaled_PA10_separation_obstruction_executable": True,
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


@dataclass(frozen=True)
class KokunoRescaledInnerJoinTshCertificate:
    """Fail closed when one PA.10 pointwise witness already breaks separation."""

    incoming: KokunoRescaledAppendixBIncomingMoments = field(
        default_factory=KokunoRescaledAppendixBIncomingMoments
    )

    def __post_init__(self) -> None:
        if not isinstance(self.incoming, KokunoRescaledAppendixBIncomingMoments):
            raise TypeError("incoming must be a KokunoRescaledAppendixBIncomingMoments")

    @property
    def witness_eta(self) -> float:
        value = float(
            self.incoming.boundary.reference.seed.phase_stationary_eta
        )
        if not math.isfinite(value) or not -1.0 <= value <= 1.0:
            raise RuntimeError("selected phase-stationary eta is invalid")
        return value

    @cached_property
    def _witness_cache(self) -> dict[str, float]:
        eta = self.witness_eta
        values = self.incoming.boundary.boundary_values(np.asarray(eta))
        ell_i = float(np.asarray(values["ell_i"]))
        if not math.isfinite(ell_i):
            raise RuntimeError("rescaled Appendix-B ell_i witness became non-finite")
        return {
            "eta": eta,
            "ell_i": ell_i,
            "abs_ell_i": abs(ell_i),
        }

    @property
    def pointwise_B0_lower_bound(self) -> float:
        return float(self._witness_cache["abs_ell_i"])

    @property
    def pointwise_T_sh_lower_bound(self) -> float:
        return 20.0 * SIGMA_PRIME_SUP * (
            self.pointwise_B0_lower_bound + LOG_F_SUP
        )

    @property
    def max_T_sh_for_separation(self) -> float:
        return (
            float(self.incoming.binding.log_X_R)
            - math.log(X_I)
            + LOG_X_RESTORE_START
        )

    @property
    def separation_margin_upper_bound(self) -> float:
        """Upper bound on possible PA.10 geometry margin from the point witness."""
        return self.max_T_sh_for_separation - self.pointwise_T_sh_lower_bound

    @property
    def selected_path_pointwise_obstructed(self) -> bool:
        # Separation requires strict T_sh < max_T_sh_for_separation.
        return bool(self.pointwise_T_sh_lower_bound >= self.max_T_sh_for_separation)

    @property
    def minimum_log_x_sep_from_witness(self) -> float:
        return (
            math.log(X_I)
            + self.pointwise_T_sh_lower_bound
            - float(self.incoming.binding.log_X_R)
        )

    def pointwise_report(self) -> dict[str, Any]:
        return {
            "witness_policy": WITNESS_POLICY,
            "witness_eta": float(self._witness_cache["eta"]),
            "ell_i_at_witness": float(self._witness_cache["ell_i"]),
            "abs_ell_i_at_witness": self.pointwise_B0_lower_bound,
            "pointwise_B0_lower_bound": self.pointwise_B0_lower_bound,
            "analytic_source_B0_bound_proved": False,
        }

    def geometry_report(self) -> dict[str, Any]:
        return {
            "selected_log_C": float(self.incoming.binding.selected_log_C),
            "selected_log_X_R": float(self.incoming.binding.log_X_R),
            "selected_log_x_i": float(self.incoming.binding.log_x_i),
            "sigma_prime_sup": SIGMA_PRIME_SUP,
            "log_f_sup": LOG_F_SUP,
            "pointwise_T_sh_lower_bound": self.pointwise_T_sh_lower_bound,
            "max_T_sh_for_current_shared_C_scale": self.max_T_sh_for_separation,
            "minimum_log_x_sep_from_witness": self.minimum_log_x_sep_from_witness,
            "required_log_x_sep_upper_bound": LOG_X_RESTORE_START,
            "separation_margin_upper_bound": self.separation_margin_upper_bound,
            "selected_path_pointwise_obstructed": self.selected_path_pointwise_obstructed,
            "fixed_B0_C_multiplier_repair_reported": False,
            "reason_no_fixed_B0_C_repair": (
                "ell_i=log(C E_i) belongs to the same coupled Appendix-B path; "
                "changing C requires recomputing upstream data rather than holding B0 fixed"
            ),
        }

    @property
    def selected_pa16_handoff_allowed(self) -> bool:
        # This screening object never promotes a handoff.  The current selected
        # path is rejected when the point witness obstructs it; if a future path
        # is not obstructed, the full source B0 certificate is still required.
        return False

    def require_selected_pa16_handoff(self) -> None:
        if self.selected_path_pointwise_obstructed:
            raise ValueError(
                "selected shared-C PA.10 path is obstructed by a pointwise B0 lower bound; "
                "PA.16 handoff remains closed"
            )
        raise ValueError(
            "pointwise PA.10 screen is not an analytic B0 certificate; "
            "PA.16 handoff remains closed"
        )

    def pa16_inputs_at_eta(self, eta: float) -> dict[str, Any]:
        self.require_selected_pa16_handoff()
        raise AssertionError("unreachable")

    def report(self) -> dict[str, Any]:
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
            "pointwise_B0_lower_bound": self.pointwise_report(),
            "geometry": self.geometry_report(),
            "selected_pa16_handoff_allowed": False,
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
            "parameters": {"witness_policy": WITNESS_POLICY},
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
        params = payload.get("parameters")
        if params != {"witness_policy": WITNESS_POLICY}:
            raise ValueError("rescaled T_sh witness policy changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("rescaled T_sh payload SHA-256 mismatch")
        deps = payload.get("dependencies")
        if not isinstance(deps, dict):
            raise ValueError("missing rescaled T_sh dependencies")
        obj = cls(
            incoming=KokunoRescaledAppendixBIncomingMoments.from_payload(
                deps.get("incoming")
            )
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
        description="Emit the selected shared-C PA.10 pointwise obstruction receipt"
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--axis-quadrature-points", type=int, default=64)
    args = parser.parse_args(argv)
    certificate = KokunoRescaledInnerJoinTshCertificate(
        incoming=KokunoRescaledAppendixBIncomingMoments(
            axis_quadrature_points=args.axis_quadrature_points
        )
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
