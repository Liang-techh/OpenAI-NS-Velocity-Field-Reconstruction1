"""Pressure/forcing-free transport-viscous block for the PA.10 center velocity.

Pinned public provenance:
KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

Agent-1 #839 exposes the same inner PA.10 contraction-center velocity together
with velocity_dt, analytic self-advection, and a frozen controlled Cartesian
Laplacian.  This module assembles the repository target operator

    T0 = partial_t u + (u . grad)u - nu Delta u,   nu = 0.01.

The viscosity is the preregistered repository problem value, not a hidden
Kokuno/OpenAI parameter.  No pressure or forcing is inserted here, so T0 is not
a complete Navier--Stokes residual.  The field remains inner-only and fails
closed whenever the inherited Laplacian stencil leaves the source inner domain.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_cartesian_center_laplacian import (
    KokunoPA10CartesianCenterLaplacian,
)
from .kokuno_pa10_cartesian_center_velocity import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)


SCHEMA = "kokuno-pa10-cartesian-center-transport-viscous-v1"
VISCOSITY = 1.0e-2

_SOURCE_FORMULAS = {
    "source_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "repository_transport_operator": "T0=partial_t u+(u.grad)u-nu Delta u",
    "repository_viscosity": "nu=0.01 from the frozen constrained problem",
    "pressure_forcing_status": "grad p and restricted forcing are not included",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "inner_cartesian_center_velocity_inherited": True,
    "inner_cartesian_center_velocity_dt_inherited": True,
    "inner_cartesian_center_self_advection_inherited": True,
    "inner_cartesian_center_velocity_laplacian_inherited": True,
    "inner_cartesian_center_transport_viscous_executable": True,
    "transport_viscosity_is_repository_choice": True,
    "transport_viscosity_caller_tunable": False,
    "matched_pressure_gradient_included": False,
    "restricted_forcing_included": False,
    "complete_momentum_residual": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "complete_candidate_api_ready": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _assemble_transport(
    velocity_dt: Any, self_advection: Any, velocity_laplacian: Any
) -> np.ndarray:
    """Assemble the frozen pressure/forcing-free transport-viscous vector."""
    dt, adv, lap = np.broadcast_arrays(
        np.asarray(velocity_dt, dtype=float),
        np.asarray(self_advection, dtype=float),
        np.asarray(velocity_laplacian, dtype=float),
    )
    if dt.shape[-1:] != (3,):
        raise ValueError("transport inputs must have trailing Cartesian component axis of length 3")
    if np.any(~np.isfinite(dt)) or np.any(~np.isfinite(adv)) or np.any(~np.isfinite(lap)):
        raise ValueError("transport inputs must contain only finite values")
    return dt + adv - VISCOSITY * lap


@dataclass(frozen=True)
class KokunoPA10CartesianCenterTransportViscous:
    """Expose T0 = u_t + (u.grad)u - 0.01 Delta u for the current inner center."""

    laplacian: KokunoPA10CartesianCenterLaplacian = field(
        default_factory=KokunoPA10CartesianCenterLaplacian
    )

    def __post_init__(self) -> None:
        if not isinstance(self.laplacian, KokunoPA10CartesianCenterLaplacian):
            raise TypeError("laplacian must be KokunoPA10CartesianCenterLaplacian")

    @property
    def field(self):
        return self.laplacian.field

    @property
    def spatial(self):
        return self.laplacian.spatial

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.laplacian.velocity(x, y, z, t)

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.laplacian.velocity_dt(x, y, z, t)

    def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.laplacian.self_advection(x, y, z, t)

    def velocity_laplacian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.laplacian.velocity_laplacian(x, y, z, t)

    def transport_viscous(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        out = _assemble_transport(
            self.velocity_dt(x, y, z, t),
            self.self_advection(x, y, z, t),
            self.velocity_laplacian(x, y, z, t),
        )
        if np.any(~np.isfinite(out)):
            raise RuntimeError("PA.10 center transport-viscous block became nonfinite")
        return out

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "laplacian_configuration": self.laplacian.configuration(),
            "viscosity": VISCOSITY,
            "operator": "velocity_dt+self_advection-viscosity*velocity_laplacian",
        }

    @property
    def field_sha256(self) -> str:
        return self.laplacian.field_sha256

    @property
    def temporal_derivative_sha256(self) -> str:
        return self.laplacian.temporal_derivative_sha256

    @property
    def spatial_derivative_sha256(self) -> str:
        return self.laplacian.spatial_derivative_sha256

    @property
    def self_advection_sha256(self) -> str:
        return self.laplacian.self_advection_sha256

    @property
    def laplacian_sha256(self) -> str:
        return self.laplacian.laplacian_sha256

    @property
    def transport_viscous_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10CartesianCenterTransportViscous":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        if payload.get("viscosity") != VISCOSITY:
            raise ValueError("viscosity drifted from the frozen repository value")
        if payload.get("operator") != "velocity_dt+self_advection-viscosity*velocity_laplacian":
            raise ValueError("unsupported transport-viscous operator")
        return cls(
            laplacian=KokunoPA10CartesianCenterLaplacian.from_configuration(
                payload["laplacian_configuration"]
            )
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10CartesianCenterTransportViscous":
        return cls.from_configuration(json.loads(Path(path).read_text()))


def build_report() -> dict[str, Any]:
    """Build a deterministic scoped engineering receipt, not a PDE receipt."""
    candidate = KokunoPA10CartesianCenterTransportViscous()
    x0, x1 = candidate.spatial.source_X_interval
    span = x1 - x0
    if not span > 0.0:
        raise RuntimeError("source inner X interval is degenerate")

    X = x0 + span * np.array([0.17, 0.21, 0.25, 0.29, 0.19, 0.27])
    eta = np.array([-0.29, -0.17, -0.04, 0.08, 0.20, 0.29])
    t = np.array([0.465, 0.485, 0.505, 0.525, 0.545, 0.495])
    theta = np.array([0.31, 0.83, 1.41, 2.03, 2.67, 3.29])
    physical = candidate.field.cartesian_from_similarity(X, eta, t, theta)

    args = (physical["x"], physical["y"], physical["z"], physical["t"])
    dt = candidate.velocity_dt(*args)
    adv = candidate.self_advection(*args)
    lap = candidate.velocity_laplacian(*args)
    transport = candidate.transport_viscous(*args)
    reconstructed = dt + adv - VISCOSITY * lap

    boundary = candidate.field.cartesian_from_similarity(x1, 0.0, 0.5, 0.0)
    boundary_fail_closed = False
    try:
        candidate.transport_viscous(
            boundary["x"], boundary["y"], boundary["z"], boundary["t"]
        )
    except ValueError:
        boundary_fail_closed = True

    checks = {
        "viscosity": VISCOSITY,
        "assembly_closure_max_abs": float(np.max(np.abs(transport - reconstructed))),
        "transport_nontrivial_rms": float(np.sqrt(np.mean(transport * transport))),
        "time_term_rms": float(np.sqrt(np.mean(dt * dt))),
        "self_advection_rms": float(np.sqrt(np.mean(adv * adv))),
        "viscous_term_rms": float(np.sqrt(np.mean((VISCOSITY * lap) ** 2))),
        "inner_boundary_stencil_fails_closed": boundary_fail_closed,
        "all_probe_values_finite": bool(
            np.all(np.isfinite(dt))
            and np.all(np.isfinite(adv))
            and np.all(np.isfinite(lap))
            and np.all(np.isfinite(transport))
        ),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
        },
        "source_formulas": _SOURCE_FORMULAS,
        "configuration": candidate.configuration(),
        "field_sha256": candidate.field_sha256,
        "temporal_derivative_sha256": candidate.temporal_derivative_sha256,
        "spatial_derivative_sha256": candidate.spatial_derivative_sha256,
        "self_advection_sha256": candidate.self_advection_sha256,
        "laplacian_sha256": candidate.laplacian_sha256,
        "transport_viscous_sha256": candidate.transport_viscous_sha256,
        "machine_checks": checks,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
        "scientific_gates": {
            "momentum_max_l2": 1.0e-3,
            "divergence_max_l2": 1.0e-5,
            "free_residual_defined_forcing_forbidden": True,
        },
    }
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path, required=True)
    args = parser.parse_args()

    candidate = KokunoPA10CartesianCenterTransportViscous()
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    candidate.save_configuration(args.config_output)
    print(json.dumps(report["machine_checks"], indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
