"""Composite-ready public bundle for the frozen Kokuno Agent-2 oscillatory field.

This module introduces no new velocity, coefficient, carrier, phase, support,
complete-curl formula, pressure, forcing, or correction.  It packages the four
already materialized public Agent-2 interfaces under one immutable provenance
identity:

    A_osc, u_osc, partial_t A_osc, partial_t u_osc.

The admitted physical velocity remains PR #561 and its independently audited
public time derivative remains PR #579.  The public vector-potential interfaces
from PRs #616/#625 remain Agent-2 self-check interfaces pending a separate
Agent-4 public black-box audit.  Bundling them here is an integration/provenance
increment only; it does not promote those scientific truth states.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_oscillatory_vector_potential import evaluate_vector_potential_osc
from .kokuno_public_oscillatory_vector_potential_time_derivative import (
    evaluate_vector_potential_osc_dt,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-PUBLIC-OSCILLATORY-BUNDLE-044"
SCHEMA = "kokuno-a2-public-oscillatory-bundle-v1"
PARENT_AGENT2_PR = 625
PARENT_AGENT2_HEAD = "e12bfb0a65031e44b1556443d88d11f2e9e92904"

UPSTREAM = {
    "velocity": {
        "pr": 561,
        "head": "732800ce4990464b49c8aa32d0dff4580f6684d4",
        "provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
        "status": "independently_admitted_public_oscillatory_velocity",
    },
    "velocity_dt": {
        "pr": 579,
        "head": "6c8e71c800a17c0e9df1802042f9feebf9dbce04",
        "provider": (
            "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:"
            "velocity_osc_dt"
        ),
        "status": "independently_admitted_public_time_derivative",
    },
    "vector_potential": {
        "pr": 616,
        "head": "aa79cb5eb601ad179c79fb8a4f3e21a093d94c31",
        "provider": (
            "openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential:"
            "vector_potential_osc"
        ),
        "status": "agent2_self_check_only_pending_independent_agent4_audit",
    },
    "vector_potential_dt": {
        "pr": 625,
        "head": PARENT_AGENT2_HEAD,
        "provider": (
            "openai_ns_reconstruction."
            "kokuno_public_oscillatory_vector_potential_time_derivative:"
            "vector_potential_osc_dt"
        ),
        "status": "agent2_self_check_only_pending_independent_agent4_audit",
    },
}


@dataclass(frozen=True)
class KokunoPublicOscillatoryBundleIdentity:
    """Immutable identity of the frozen four-interface Agent-2 handoff."""

    schema: str = SCHEMA
    parent_agent2_pr: int = PARENT_AGENT2_PR
    parent_agent2_head: str = PARENT_AGENT2_HEAD

    def payload(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "parent_agent2_pr": self.parent_agent2_pr,
            "parent_agent2_head": self.parent_agent2_head,
            "upstream": UPSTREAM,
            "coordinate_contract": (
                "public R=hypot(x,y), theta=atan2(y,x), z; per-beta source "
                "Z_beta=epsilon_beta*z so public partial_z equals source D_z"
            ),
            "representation_contract": (
                "localized m=+/-1 real pair assembled at vector-potential level, "
                "then complete-curl velocity; no free Cartesian component fitting"
            ),
            "provenance_boundary": {
                "kokuno_corrected_source_formula_reused": True,
                "repository_public_coordinate_gauge_realization": True,
                "repository_autonomous_time_modulation": True,
                "paper_exact": False,
                "openai_field_identified": False,
            },
        }

    def sha256(self) -> str:
        raw = json.dumps(
            self.payload(), sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()


IDENTITY = KokunoPublicOscillatoryBundleIdentity()


def evaluate_oscillatory_bundle(x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
    """Evaluate all frozen Agent-2 public oscillatory interfaces in one call.

    Inputs follow the existing public providers and are NumPy-broadcastable.
    The returned arrays all have trailing Cartesian vector dimension three.
    No numerical spatial or temporal differentiation is performed in this
    adapter; it delegates to the already frozen public interfaces.
    """
    potential = evaluate_vector_potential_osc(x, y, z, t)
    potential_dt = evaluate_vector_potential_osc_dt(x, y, z, t)
    u = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    u_t = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)
    A = np.asarray(potential["vector_potential_cartesian_total"], dtype=float)
    A_t = np.asarray(
        potential_dt["vector_potential_dt_cartesian_total"], dtype=float
    )
    support = np.asarray(potential["support_mask"], dtype=bool)

    shapes = {tuple(value.shape) for value in (A, u, A_t, u_t)}
    if len(shapes) != 1:
        raise RuntimeError("frozen Agent-2 public interfaces returned inconsistent shapes")
    shape = next(iter(shapes))
    if not shape or shape[-1] != 3:
        raise RuntimeError("frozen Agent-2 public interfaces lost Cartesian vector shape")
    if support.shape != shape[:-1]:
        raise RuntimeError("frozen Agent-2 support mask shape changed")
    if not all(np.all(np.isfinite(value)) for value in (A, u, A_t, u_t)):
        raise RuntimeError("frozen Agent-2 public bundle produced non-finite values")

    return {
        "vector_potential": A,
        "velocity": u,
        "vector_potential_dt": A_t,
        "velocity_dt": u_t,
        "support_mask": support,
        "bundle_identity_sha256": IDENTITY.sha256(),
        "bundle_identity": IDENTITY.payload(),
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "vector_potential_candidate_changed": False,
            "source_formula_changed": False,
            "integration_adapter_only": True,
            "independent_agent4_vector_potential_audit_required": True,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = np.asarray(
        [
            (0.51, 0.13, -0.47, 0.36),
            (-0.31, 0.66, -0.47, 0.36),
            (0.82, -0.36, 0.43, 0.36),
            (-0.74, -0.55, 0.43, 0.36),
            (0.55, -0.18, -0.47, 0.50),
            (-0.29, -0.70, -0.47, 0.50),
            (0.85, 0.31, 0.43, 0.50),
            (-0.78, 0.49, 0.43, 0.50),
            (0.48, 0.24, -0.47, 0.64),
            (-0.38, 0.62, -0.47, 0.64),
            (0.79, -0.42, 0.43, 0.64),
            (-0.70, -0.61, 0.43, 0.64),
        ],
        dtype=float,
    )
    return rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]


def _max_abs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def verification_receipt() -> dict[str, Any]:
    """Verify exact API replay and frozen support without new scientific claims."""
    x, y, z, t = _verification_cloud()
    bundle = evaluate_oscillatory_bundle(x, y, z, t)
    direct_potential = np.asarray(
        evaluate_vector_potential_osc(x, y, z, t)["vector_potential_cartesian_total"],
        dtype=float,
    )
    direct_potential_dt = np.asarray(
        evaluate_vector_potential_osc_dt(x, y, z, t)[
            "vector_potential_dt_cartesian_total"
        ],
        dtype=float,
    )
    direct_velocity = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    direct_velocity_dt = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)

    replay_errors = {
        "vector_potential": _max_abs(bundle["vector_potential"], direct_potential),
        "velocity": _max_abs(bundle["velocity"], direct_velocity),
        "vector_potential_dt": _max_abs(
            bundle["vector_potential_dt"], direct_potential_dt
        ),
        "velocity_dt": _max_abs(bundle["velocity_dt"], direct_velocity_dt),
    }

    exterior_x = np.asarray((0.0, 1.85, 0.8), dtype=float)
    exterior_y = np.asarray((0.0, 0.0, 0.0), dtype=float)
    exterior_z = np.asarray((0.0, 0.0, 2.15), dtype=float)
    exterior_t = np.asarray((0.5, 0.5, 0.5), dtype=float)
    exterior = evaluate_oscillatory_bundle(
        exterior_x, exterior_y, exterior_z, exterior_t
    )
    exterior_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in ("vector_potential", "velocity", "vector_potential_dt", "velocity_dt")
    )

    failed_guards: list[str] = []
    if any(value != 0.0 for value in replay_errors.values()):
        failed_guards.append("exact_public_api_replay")
    if exterior_max != 0.0:
        failed_guards.append("registered_support_exterior_zero")
    if not np.any(bundle["support_mask"]):
        failed_guards.append("nontrivial_support_intersection")
    if bundle["bundle_identity_sha256"] != IDENTITY.sha256():
        failed_guards.append("bundle_identity_replay")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "bundle_identity_sha256": IDENTITY.sha256(),
        "sample_count": int(t.size),
        "exact_public_api_replay_max_abs": replay_errors,
        "support_interior_count": int(np.count_nonzero(bundle["support_mask"])),
        "support_exterior_absolute_max": exterior_max,
        "failed_guards": failed_guards,
        "upstream": UPSTREAM,
        "truth_boundary": bundle["truth_boundary"],
        "scientific_scope": (
            "integration/provenance replay only; no new curl/divergence/NS residual "
            "admission and no replacement for independent Agent-4 vector-potential audit"
        ),
    }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = verification_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
