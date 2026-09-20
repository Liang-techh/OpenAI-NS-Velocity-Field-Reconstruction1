"""CR002 audit for the ST052-M endpoint temporal-curvature derivative scope.

Agent-7 PR #818 preregisters an autonomous response direction

    g2(t) = 16 (t - 0.25) (0.75 - t)

whose *velocity-value* contribution vanishes at the registered time endpoints.
That statement must not be silently upgraded to preservation of ``velocity_dt``
or of the Navier-Stokes momentum residual: ``g2'`` is +8 at t=.25 and -8
at t=.75.  This module locks that distinction without fitting a coefficient or
changing any candidate field.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-st052m-endpoint-temporal-derivative-scope-v1"
TASK = "CR002-ST052M-ENDPOINT-TEMPORAL-DERIVATIVE-SCOPE-117"
ROOT = Path(__file__).resolve().parents[2]
SCOPE_PATH = ROOT / "configs" / "st052m_endpoint_temporal_derivative_scope.json"
CR001_PATH = ROOT / "configs" / "constraints.json"
UPSTREAM_PATH = (
    ROOT
    / "experiments"
    / "root_st052"
    / "agent7_st052m_endpoint_temporal_curvature_preflight.py"
)
EXPECTED_UPSTREAM_BLOB = "1d5e5353e3792ceb0a3f99d16749cbb63e3cc301"
EXPECTED_CR001_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _git_blob_sha1(path: Path) -> str:
    raw = path.read_bytes()
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def load_scope() -> dict[str, Any]:
    return _read_json(SCOPE_PATH)


def _g2(time: float) -> float:
    t = float(time)
    return 16.0 * (t - 0.25) * (0.75 - t)


def _g2_prime(time: float) -> float:
    t = float(time)
    return 16.0 * (1.0 - 2.0 * t)


def derivative_witness() -> dict[str, float]:
    values = {
        "g2_start": _g2(0.25),
        "g2_mid": _g2(0.50),
        "g2_end": _g2(0.75),
        "g2_prime_start": _g2_prime(0.25),
        "g2_prime_mid": _g2_prime(0.50),
        "g2_prime_end": _g2_prime(0.75),
    }
    expected = {
        "g2_start": 0.0,
        "g2_mid": 1.0,
        "g2_end": 0.0,
        "g2_prime_start": 8.0,
        "g2_prime_mid": 0.0,
        "g2_prime_end": -8.0,
    }
    for key, target in expected.items():
        if not math.isclose(values[key], target, rel_tol=0.0, abs_tol=1.0e-15):
            raise ValueError(f"temporal derivative witness drift: {key}")
    return values


def _validate_upstream_source() -> None:
    if _git_blob_sha1(UPSTREAM_PATH) != EXPECTED_UPSTREAM_BLOB:
        raise ValueError("audited Agent-7 #818 source blob drift")
    source = UPSTREAM_PATH.read_text(encoding="utf-8")
    required_fragments = (
        "TIME_START = 0.25",
        "TIME_END = 0.75",
        "TIME_MID = 0.50",
        "return float(16.0 * (t - TIME_START) * (TIME_END - t))",
        '"temporal_coefficient_selected": False',
        '"held_out_pde_residual_evaluated": False',
        '"public_image_numeric_target_used": False',
        '"pde_validated": False',
        '"paper_exact": False',
        '"openai_field_identified": False',
    )
    for fragment in required_fragments:
        if fragment not in source:
            raise ValueError(f"Agent-7 #818 source semantics drift: {fragment}")


def _validate_cr001() -> None:
    if _git_blob_sha1(CR001_PATH) != EXPECTED_CR001_BLOB:
        raise ValueError("CR001 canonical contract blob drift")
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
        raise ValueError("CR001 forcing bounds drift")
    if "No residual-dependent basis or pointwise free force" not in forcing["restriction"]:
        raise ValueError("CR001 free-force prohibition drift")

    nontriviality = data["nontriviality"]
    if nontriviality["reference_energy"] != 1.0:
        raise ValueError("CR001 reference-energy drift")
    if nontriviality["reference_energy_abs_tolerance"] != 0.001:
        raise ValueError("CR001 reference-energy tolerance drift")
    if "reject collapsed candidates" not in nontriviality["enforcement"]:
        raise ValueError("CR001 amplitude-collapse prohibition drift")

    validation = data["validation"]
    if validation["seed"] != 914027 or validation["held_out_points"] != 4096:
        raise ValueError("CR001 held-out sample drift")
    if validation["times"] != [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]:
        raise ValueError("CR001 validation-time drift")
    if validation["derivative_steps"] != [0.02, 0.01, 0.005]:
        raise ValueError("CR001 derivative-ladder drift")
    if validation["quadrature_orders_per_axis"] != [24, 48, 96]:
        raise ValueError("CR001 quadrature-ladder drift")
    if validation["time_derivatives"] != (
        "one-sided calibrated stencils at time endpoints; centered stencils only "
        "strictly within the declared window"
    ):
        raise ValueError("CR001 endpoint derivative policy drift")
    thresholds = validation["thresholds"]
    expected_thresholds = {
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
        "pde_residual_max": 1.0e-3,
        "pde_residual_L2": 1.0e-3,
    }
    for key, value in expected_thresholds.items():
        if thresholds[key] != value:
            raise ValueError(f"CR001 threshold drift: {key}")
    if validation["failure_policy"] != (
        "retain failed results; changing thresholds requires a new experiment version"
    ):
        raise ValueError("CR001 failure-policy drift")


def validate_scope(payload: Mapping[str, Any]) -> None:
    if not isinstance(payload, Mapping):
        raise TypeError("scope payload must be a mapping")
    expected = load_scope()
    if dict(payload) != expected:
        raise ValueError("endpoint temporal derivative scope contract drift")
    if payload.get("schema") != SCHEMA or payload.get("task") != TASK:
        raise ValueError("scope identity drift")

    stack = payload["stack"]
    if stack != {
        "audited_pr": 818,
        "audited_head": "c5acc7f8f65b64885ec5660ed9c991696668e300",
        "audited_source_path": (
            "experiments/root_st052/agent7_st052m_endpoint_temporal_curvature_preflight.py"
        ),
        "audited_source_blob": EXPECTED_UPSTREAM_BLOB,
        "source_parent_pr": 810,
        "frozen_child_pr": 775,
    }:
        raise ValueError("upstream stack binding drift")

    classes = payload["classification"]
    if set(classes) != {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }:
        raise ValueError("four-way provenance classification drift")
    if classes["public_source_fact"] != []:
        raise ValueError("autonomous g2 data laundered into public-source facts")

    basis = payload["time_basis"]
    if basis["coefficient_selected"] is not False:
        raise ValueError("temporal coefficient prematurely selected")
    if basis["coefficient_may_not_be_selected_from_cr001_held_out_pde_samples"] is not True:
        raise ValueError("held-out tuning firewall lost")
    if basis["public_openai_numeric_time_law_claimed"] is not False:
        raise ValueError("autonomous g2 promoted to OpenAI numeric time law")
    witness = derivative_witness()
    if basis["witness"] != witness:
        raise ValueError("stored temporal derivative witness drift")

    scope = payload["scope"]
    required_true = {
        "future_child_assumes_static_spatial_channel_for_this_identity",
        "g2_zero_endpoints_preserve_additive_velocity_values_at_endpoints",
        "g2_zero_at_t025_preserves_initial_kinetic_energy_value_if_parent_and_domain_are_unchanged",
        "nonzero_g2_prime_endpoints_require_explicit_velocity_dt_accounting_if_beta_times_C_is_nonzero",
        "initial_energy_value_preservation_is_not_pde_validation",
        "initial_energy_value_preservation_is_not_new_nontriviality_evidence",
        "endpoint_velocity_value_preservation_is_not_visual_correspondence",
        "identifiability_pass_would_only_establish_local_expression_capacity",
        "identifiability_pass_would_not_select_beta",
    }
    required_false = {
        "preflight_changes_candidate_velocity",
        "g2_zero_endpoints_imply_velocity_dt_endpoint_preservation",
        "g2_zero_endpoints_imply_momentum_residual_endpoint_neutrality",
        "held_out_pde_residual_evaluated_by_this_increment",
    }
    for key in required_true:
        if scope.get(key) is not True:
            raise ValueError(f"required temporal-scope truth lost: {key}")
    for key in required_false:
        if scope.get(key) is not False:
            raise ValueError(f"forbidden temporal-scope promotion: {key}")

    states = payload["independent_states"]
    if states["existing_callable_velocity_delivery_invalidated"] is not False:
        raise ValueError("PDE/temporal audit cannot invalidate callable delivery")
    for key in (
        "velocity_export_ready_promoted_by_this_contract",
        "temporal_child_materialized",
        "temporal_coefficient_selected",
        "visualization_ready_promoted_by_this_contract",
        "visual_correspondence_verified",
        "source_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if states[key] is not False:
            raise ValueError(f"truth-state promotion forbidden: {key}")

    _validate_upstream_source()
    _validate_cr001()


def audit() -> dict[str, Any]:
    scope = load_scope()
    validate_scope(scope)
    return {
        "schema": SCHEMA,
        "task": TASK,
        "scope_validated": True,
        "derivative_witness": derivative_witness(),
        "interpretation": (
            "g2 preserves the additive velocity value at t=.25/.75 but has nonzero "
            "endpoint time derivative; velocity_dt and momentum-residual consequences "
            "therefore require separate evaluation after any coefficient is frozen"
        ),
        "truth_boundary": copy.deepcopy(scope["independent_states"]),
    }


def main() -> int:
    print(json.dumps(audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
