"""CR002 audit for the PA.10 inner velocity_dt domain-neighborhood scope.

The analytic ``velocity_dt`` introduced on Agent-1 PR #811 is a pointwise
fixed-Cartesian derivative on the current *inner* PA.10 representation.
Agent-4 PR #814 independently checks it with finite differences, but its random
protocol is deliberately well inside the source-X interval.

This audit machine-locks the distinction between those statements.  In
particular, a point can be valid for the analytic derivative while a particular
finite-difference stencil around that same fixed Cartesian point leaves the
currently materialized inner-X representation.  That is a representation-domain
fact, not a physical singularity and not evidence against a future global join.
"""

from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import kokuno_a4_pa10_cartesian_center_velocity_dt_independent_audit as agent4
from .kokuno_pa10_cartesian_center_velocity_dt import (
    KokunoPA10CartesianCenterVelocityTimeDerivative,
)

SCHEMA = "cr002-kokuno-velocity-dt-domain-neighborhood-scope-v1"
TASK = "CR002-KOKUNO-VELOCITY-DT-DOMAIN-NEIGHBORHOOD-073"
ROOT = Path(__file__).resolve().parents[2]
SCOPE_PATH = ROOT / "configs" / "kokuno_velocity_dt_domain_neighborhood_scope.json"
CR001_PATH = ROOT / "configs" / "constraints.json"

EXPECTED_AGENT4_HEAD = "6f2c9789921885a2c21bb69c642d7f820fdeed41"
EXPECTED_AGENT1_HEAD = "6f16e837c3a95a6e54af5c2cb7725d094f667026"
EXPECTED_AGENT4_SOURCE_BLOB = "210f445053f027034f3fe20dc7765cd32be1e74d"
EXPECTED_AGENT1_SOURCE_BLOB = "c1d576820416945443bbc24cc931cb64e9e032d9"


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_scope() -> dict[str, Any]:
    return _read_json(SCOPE_PATH)


def _validate_cr001() -> None:
    data = _read_json(CR001_PATH)
    if data["nu"] != 0.01:
        raise ValueError("CR001 nu drift")
    domain = data["domain"]
    if domain["physical"] != "R^3":
        raise ValueError("CR001 physical-domain drift")
    if domain["evaluation_box"] != [[-2, 2], [-2, 2], [-2, 2]]:
        raise ValueError("CR001 evaluation-box drift")
    if domain["support"] != "r < 2 and abs(z) < 2":
        raise ValueError("CR001 support drift")
    if domain["time_interval"] != [0.25, 0.75]:
        raise ValueError("CR001 time-window drift")

    forcing = data["forcing"]
    if forcing["mode"] != "restricted_two_parameter_family":
        raise ValueError("CR001 forcing-family drift")
    if forcing["parameters"] != {"a": [0.0, 10.0], "c": [0.0, 10.0]}:
        raise ValueError("CR001 forcing-bound drift")
    restriction = forcing["restriction"]
    if "No residual-dependent basis or pointwise free force" not in restriction:
        raise ValueError("CR001 free-force prohibition drift")

    nontriviality = data["nontriviality"]
    if nontriviality["reference_energy"] != 1.0:
        raise ValueError("CR001 reference-energy drift")
    if nontriviality["reference_energy_abs_tolerance"] != 0.001:
        raise ValueError("CR001 energy-tolerance drift")
    if "reject collapsed candidates" not in nontriviality["enforcement"]:
        raise ValueError("CR001 amplitude-collapse prohibition drift")

    validation = data["validation"]
    if validation["seed"] != 914027 or validation["held_out_points"] != 4096:
        raise ValueError("CR001 validation-sample drift")
    if validation["derivative_steps"] != [0.02, 0.01, 0.005]:
        raise ValueError("CR001 derivative-ladder drift")
    if validation["quadrature_orders_per_axis"] != [24, 48, 96]:
        raise ValueError("CR001 quadrature-ladder drift")
    thresholds = validation["thresholds"]
    expected = {
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
        "pde_residual_max": 1.0e-3,
        "pde_residual_L2": 1.0e-3,
    }
    for key, value in expected.items():
        if thresholds[key] != value:
            raise ValueError(f"CR001 threshold drift: {key}")
    if validation["failure_policy"] != (
        "retain failed results; changing thresholds requires a new experiment version"
    ):
        raise ValueError("CR001 failure-policy drift")
    if validation["time_derivatives"] != (
        "one-sided calibrated stencils at time endpoints; centered stencils only "
        "strictly within the declared window"
    ):
        raise ValueError("CR001 time-derivative policy drift")


def _validate_agent4_protocol_source() -> None:
    if agent4.SEED != 9173471:
        raise ValueError("Agent-4 seed drift")
    if agent4.SAMPLE_COUNT != 1024:
        raise ValueError("Agent-4 sample-count drift")
    if tuple(agent4.TIME_STEPS) != (4.0e-4, 2.0e-4, 1.0e-4):
        raise ValueError("Agent-4 FD4 step drift")

    source = inspect.getsource(agent4._sample_fixed_cartesian_points)
    required_fragments = (
        "rng.uniform(0.06 * xmax, 0.30 * xmax, size=SAMPLE_COUNT)",
        "rng.uniform(-0.42, 0.42, size=SAMPLE_COUNT)",
        "rng.uniform(0.34, 0.66, size=SAMPLE_COUNT)",
    )
    for fragment in required_fragments:
        if fragment not in source:
            raise ValueError(f"Agent-4 interior sampling source drift: {fragment}")


def validate_scope(payload: Mapping[str, Any]) -> None:
    if not isinstance(payload, Mapping):
        raise TypeError("scope payload must be a mapping")
    expected = load_scope()
    if dict(payload) != expected:
        raise ValueError("domain-neighborhood scope contract drift")

    if payload.get("schema") != SCHEMA or payload.get("task") != TASK:
        raise ValueError("scope identity drift")
    stack = payload["stack"]
    if stack != {
        "agent1_velocity_dt_head": EXPECTED_AGENT1_HEAD,
        "agent1_velocity_dt_pr": 811,
        "agent1_velocity_dt_source_blob": EXPECTED_AGENT1_SOURCE_BLOB,
        "agent4_audit_source_blob": EXPECTED_AGENT4_SOURCE_BLOB,
        "audited_head": EXPECTED_AGENT4_HEAD,
        "audited_pr": 814,
    }:
        raise ValueError("upstream binding drift")

    classes = payload["classification"]
    if set(classes) != {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }:
        raise ValueError("four-way provenance classification drift")

    domain = payload["domain_scope"]
    required_true = {
        "inner_source_X_only",
        "outside_inner_X_raises",
        "analytic_velocity_dt_is_fixed_cartesian",
        "pointwise_velocity_dt_evaluable_does_not_imply_frozen_fd_stencil_fits",
        "agent4_fd4_evidence_is_interior_subdomain_evidence",
    }
    required_false = {
        "agent4_fd4_evidence_establishes_entire_inner_X_domain",
        "agent4_fd4_evidence_establishes_inner_boundary_behavior",
        "agent4_fd4_evidence_establishes_global_join_behavior",
        "full_inner_domain_pde_derivative_validated",
        "global_pde_derivative_validated",
    }
    for key in required_true:
        if domain.get(key) is not True:
            raise ValueError(f"required domain-scope truth lost: {key}")
    for key in required_false:
        if domain.get(key) is not False:
            raise ValueError(f"premature domain-scope promotion: {key}")

    states = payload["independent_states"]
    if states["inner_velocity_dt_interface_available"] is not True:
        raise ValueError("inner velocity_dt availability lost")
    if states["inner_velocity_dt_agent4_interior_audit_registered"] is not True:
        raise ValueError("Agent-4 audit registration lost")
    for key in (
        "velocity_export_ready_promoted_by_this_contract",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if states[key] is not False:
            raise ValueError(f"truth-state promotion forbidden: {key}")

    _validate_cr001()
    _validate_agent4_protocol_source()


def boundary_witness() -> dict[str, Any]:
    """Return a deterministic witness separating pointwise and stencil domains."""
    scope = load_scope()
    spec = scope["boundary_witness"]
    provider = KokunoPA10CartesianCenterVelocityTimeDerivative()
    x_min, x_max = provider.source_X_interval
    if not (x_max > x_min >= 0.0):
        raise RuntimeError("unexpected source-X interval")

    fraction = float(spec["source_X_fraction_of_upper_bound"])
    X0 = fraction * x_max
    eta0 = float(spec["eta"])
    t0 = float(spec["t"])
    h = float(spec["fd4_step"])
    point = provider.field.cartesian_from_similarity(X0, eta0, t0, 0.0)

    x = point["x"]
    y = point["y"]
    z = point["z"]
    dt_value = np.asarray(provider.velocity_dt(x, y, z, t0), dtype=float)
    if dt_value.shape != (3,) or not np.all(np.isfinite(dt_value)):
        raise RuntimeError("pointwise inner velocity_dt witness is not finite")

    plus2_t = t0 + 2.0 * h
    shifted = provider.field.similarity_coordinates(x, y, z, plus2_t)
    shifted_X = float(np.asarray(shifted["X"]))
    shifted_fraction = shifted_X / x_max
    if not shifted_fraction > 1.0:
        raise RuntimeError("configured boundary witness did not leave inner X")

    outside_error = False
    try:
        provider.velocity(x, y, z, plus2_t)
    except ValueError:
        outside_error = True
    if not outside_error:
        raise RuntimeError("current inner velocity unexpectedly accepted outer-X stencil point")

    q0 = 1.0 - t0
    crossing_delta = q0 * (1.0 - fraction)
    expected_delta = float(spec["expected_crossing_delta_at_eta_zero"])
    if not np.isclose(crossing_delta, expected_delta, rtol=0.0, atol=5.0e-16):
        raise RuntimeError("boundary-witness crossing delta drift")
    if not (crossing_delta < 2.0 * h):
        raise RuntimeError("boundary witness no longer crossed by fine FD4 outer point")

    return {
        "source_X_interval": [float(x_min), float(x_max)],
        "source_X_fraction_at_center_time": fraction,
        "pointwise_velocity_dt_finite": True,
        "pointwise_velocity_dt_norm": float(np.linalg.norm(dt_value)),
        "fd4_step": h,
        "positive_outer_offset": 2.0 * h,
        "crossing_delta_at_eta_zero": crossing_delta,
        "shifted_source_X_fraction": shifted_fraction,
        "velocity_at_t_plus_2h_rejected_by_inner_domain": outside_error,
        "interpretation": (
            "current representation-domain witness only; not a physical singularity "
            "and not a statement about a future global join"
        ),
    }


def audit() -> dict[str, Any]:
    scope = load_scope()
    validate_scope(scope)
    witness = boundary_witness()
    return {
        "schema": SCHEMA,
        "task": TASK,
        "scope_validated": True,
        "boundary_witness": witness,
        "truth_boundary": copy.deepcopy(scope["independent_states"]),
    }


def main() -> int:
    print(json.dumps(audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
