"""Fail-closed Agent-3 admission from an independent oscillatory audit.

Agent 3 must not turn a source-compatible oscillatory provider into a same-cycle
mean correction merely because it exposes ``velocity_osc(x,y,z,t)``.  The
public provider must first pass the independently frozen Agent-4 black-box
preflight for divergence, phase-mean covariance rank, support, nontriviality,
parameter responsiveness, and mutation sensitivity.

This module consumes the *audit receipt*, recomputes the frozen guard verdict,
and only opens the next materialization stage when every guard passes.  A PASS
here still does not materialize ``requestedStress``, a finite-head mean debt,
a velocity correction, or a Navier--Stokes residual reduction.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

TASK = "KOKUNO-A3-INDEPENDENT-OSCILLATORY-ADMISSION-035"
SCHEMA = "kokuno-agent3-independent-oscillatory-admission-v1"
EXPECTED_AUDIT_SCHEMA = "kokuno-agent4-public-oscillatory-independent-audit-v1"
EXPECTED_AGENT2_HEAD = "46c1699c7fc7ab2dc103d5cf430ead3a39bbbc1e"
EXPECTED_SEED = 9173241
EXPECTED_FD4_STEPS = (0.02, 0.01, 0.005)
EXPECTED_PHASE_RESOLUTIONS = (32, 64, 128)

# Copied exactly from the independently preregistered Agent-4 audit.  The
# caller cannot weaken these values through its receipt.
FROZEN_GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 3.0,
    "minimum_covariance_rank_ratio": 2.0e-2,
    "maximum_covariance_resolution_drift": 1.0e-6,
    "axis_near_absolute_max": 1.0e-14,
    "radial_exterior_absolute_max": 1.0e-14,
    "project_axial_exterior_absolute_max": 1.0e-12,
    "minimum_nontrivial_velocity_rms": 1.0e-8,
    "minimum_parameter_perturbation_relative_change": 1.0e-4,
    "minimum_divergence_mutation_relative_rms": 5.0e-2,
    "maximum_duplicated_covariance_rank_ratio": 1.0e-8,
}


def _mapping(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _sequence(value: object, name: str, length: int) -> Sequence[Any]:
    if not isinstance(value, (list, tuple)) or len(value) != length:
        raise ValueError(f"{name} must have length {length}")
    return value


def _finite(value: object, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _strict_bool(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be boolean")
    return value


def _same_float(actual: object, expected: float, name: str) -> None:
    value = _finite(actual, name)
    if value != float(expected):
        raise ValueError(f"{name} changed from the frozen independent-audit value")


def _check_receipt_identity(audit: Mapping[str, Any]) -> None:
    if audit.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("unexpected independent oscillatory audit schema")
    if audit.get("parent_agent2_head") != EXPECTED_AGENT2_HEAD:
        raise ValueError("audit is not bound to the frozen Agent-2 public provider head")
    if int(audit.get("seed", -1)) != EXPECTED_SEED:
        raise ValueError("independent audit seed changed")

    protocol = _mapping(audit.get("protocol"), "protocol")
    if tuple(protocol.get("fd4_steps", ())) != EXPECTED_FD4_STEPS:
        raise ValueError("independent FD4 ladder changed")
    if tuple(protocol.get("phase_resolutions", ())) != EXPECTED_PHASE_RESOLUTIONS:
        raise ValueError("independent phase-resolution ladder changed")
    if protocol.get("primary_interface") != "module-level velocity_osc(x,y,z,t) only":
        raise ValueError("audit did not use the public oscillatory interface")
    for key in (
        "pressure_or_forcing_fit",
        "candidate_internal_derivative_oracle_used",
        "training_tensor_or_loss_used",
    ):
        if protocol.get(key) is not False:
            raise ValueError(f"audit protocol violated fail-closed rule: {key}")

    guards = _mapping(audit.get("guards_frozen_before_actions"), "guards")
    if set(guards) != set(FROZEN_GUARDS):
        raise ValueError("audit guard set differs from the frozen Agent-4 protocol")
    for key, expected in FROZEN_GUARDS.items():
        _same_float(guards[key], expected, f"guards.{key}")

    truth = _mapping(audit.get("truth_boundary"), "truth_boundary")
    if truth.get("oscillatory_component_only") is not True:
        raise ValueError("audit must remain oscillatory-component-only")
    for key in (
        "global_leading_velocity_pressure_available",
        "materialized_correction_available",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"independent audit truth boundary was improperly promoted: {key}")
    _same_float(truth.get("formal_momentum_normalized_max_l2_gate"), 1.0e-3, "momentum gate")
    _same_float(truth.get("formal_divergence_max_l2_gate"), 1.0e-5, "divergence gate")


def evaluate_independent_oscillatory_admission(audit: Mapping[str, Any]) -> dict[str, Any]:
    """Recompute the independent preflight and gate Agent-3 mean-debt work."""

    audit = _mapping(audit, "audit")
    _check_receipt_identity(audit)
    failures: list[str] = []

    divergence = _mapping(audit.get("divergence"), "divergence")
    div_rows = _sequence(divergence.get("rows"), "divergence.rows", 3)
    for index, (row_raw, step) in enumerate(zip(div_rows, EXPECTED_FD4_STEPS, strict=True)):
        row = _mapping(row_raw, f"divergence.rows[{index}]")
        _same_float(row.get("step"), step, f"divergence.rows[{index}].step")
    finest = _mapping(div_rows[-1], "divergence.rows[-1]")
    finest_rms = _finite(finest.get("relative_rms"), "finest relative divergence rms")
    finest_max = _finite(finest.get("relative_max"), "finest relative divergence max")
    refinement = tuple(
        _finite(v, "divergence refinement ratio")
        for v in _sequence(
            divergence.get("relative_rms_refinement_ratios"),
            "divergence.relative_rms_refinement_ratios",
            2,
        )
    )
    if finest_rms > FROZEN_GUARDS["finest_relative_divergence_rms"]:
        failures.append("finest_relative_divergence_rms")
    if finest_max > FROZEN_GUARDS["finest_relative_divergence_max"]:
        failures.append("finest_relative_divergence_max")
    if min(refinement) < FROZEN_GUARDS["minimum_divergence_refinement_ratio"]:
        failures.append("minimum_divergence_refinement_ratio")
    divergence_passed = not any(name.startswith("finest_") or name == "minimum_divergence_refinement_ratio" for name in failures)
    if _strict_bool(divergence.get("local_divergence_passed"), "local_divergence_passed") != divergence_passed:
        raise ValueError("reported local divergence verdict disagrees with frozen guard recomputation")

    covariance = _mapping(audit.get("phase_mean_covariance"), "phase_mean_covariance")
    cov_rows = _sequence(covariance.get("rows"), "phase_mean_covariance.rows", 3)
    rank_ratios: list[float] = []
    for index, (row_raw, resolution) in enumerate(zip(cov_rows, EXPECTED_PHASE_RESOLUTIONS, strict=True)):
        row = _mapping(row_raw, f"phase_mean_covariance.rows[{index}]")
        if int(row.get("phase_resolution", -1)) != resolution:
            raise ValueError("phase covariance resolution changed")
        rank_ratios.append(_finite(row.get("rank_ratio"), "covariance rank ratio"))
    drift = tuple(
        _finite(v, "covariance resolution drift")
        for v in _sequence(
            covariance.get("relative_drift_to_finest"),
            "phase_mean_covariance.relative_drift_to_finest",
            2,
        )
    )
    duplicate_ratio = _finite(
        covariance.get("duplicated_first_column_rank_ratio"),
        "duplicated covariance rank ratio",
    )
    covariance_passed = True
    if min(rank_ratios) < FROZEN_GUARDS["minimum_covariance_rank_ratio"]:
        failures.append("minimum_covariance_rank_ratio")
        covariance_passed = False
    if max(drift) > FROZEN_GUARDS["maximum_covariance_resolution_drift"]:
        failures.append("maximum_covariance_resolution_drift")
        covariance_passed = False
    if duplicate_ratio > FROZEN_GUARDS["maximum_duplicated_covariance_rank_ratio"]:
        failures.append("maximum_duplicated_covariance_rank_ratio")
        covariance_passed = False
    if _strict_bool(covariance.get("covariance_rank_preflight_passed"), "covariance verdict") != covariance_passed:
        raise ValueError("reported covariance verdict disagrees with frozen guard recomputation")

    support = _mapping(audit.get("support"), "support")
    axis_max = _finite(support.get("axis_near_absolute_max"), "axis-near absolute max")
    radial_max = _finite(support.get("radial_exterior_absolute_max"), "radial exterior absolute max")
    axial_max = _finite(support.get("project_axial_exterior_absolute_max"), "axial exterior absolute max")
    radial_axis_support_passed = True
    if axis_max > FROZEN_GUARDS["axis_near_absolute_max"]:
        failures.append("axis_near_absolute_max")
        radial_axis_support_passed = False
    if radial_max > FROZEN_GUARDS["radial_exterior_absolute_max"]:
        failures.append("radial_exterior_absolute_max")
        radial_axis_support_passed = False
    axial_support_passed = axial_max <= FROZEN_GUARDS["project_axial_exterior_absolute_max"]
    if not axial_support_passed:
        failures.append("project_axial_exterior_absolute_max")
    if _strict_bool(support.get("radial_axis_support_passed"), "radial support verdict") != radial_axis_support_passed:
        raise ValueError("reported radial/axis support verdict disagrees with recomputation")
    if _strict_bool(support.get("project_axial_support_passed"), "axial support verdict") != axial_support_passed:
        raise ValueError("reported axial support verdict disagrees with recomputation")

    nontriviality = _mapping(audit.get("nontriviality"), "nontriviality")
    velocity_rms = _finite(nontriviality.get("heldout_velocity_rms"), "heldout velocity rms")
    nontriviality_passed = velocity_rms >= FROZEN_GUARDS["minimum_nontrivial_velocity_rms"]
    if not nontriviality_passed:
        failures.append("minimum_nontrivial_velocity_rms")
    if _strict_bool(nontriviality.get("passed"), "nontriviality verdict") != nontriviality_passed:
        raise ValueError("reported nontriviality verdict disagrees with recomputation")

    perturbation = _mapping(audit.get("parameter_perturbation"), "parameter_perturbation")
    perturbation_values = (
        _finite(perturbation.get("h_minus_10pct_relative_change"), "h- perturbation"),
        _finite(perturbation.get("h_plus_10pct_relative_change"), "h+ perturbation"),
    )
    perturbation_passed = min(perturbation_values) >= FROZEN_GUARDS["minimum_parameter_perturbation_relative_change"]
    if not perturbation_passed:
        failures.append("minimum_parameter_perturbation_relative_change")
    if _strict_bool(perturbation.get("passed"), "parameter perturbation verdict") != perturbation_passed:
        raise ValueError("reported parameter-perturbation verdict disagrees with recomputation")

    mutation = _mapping(audit.get("divergence_mutation"), "divergence_mutation")
    mutation_rms = _finite(mutation.get("relative_rms"), "divergence mutation relative rms")
    mutation_passed = mutation_rms >= FROZEN_GUARDS["minimum_divergence_mutation_relative_rms"]
    if not mutation_passed:
        failures.append("minimum_divergence_mutation_relative_rms")
    if _strict_bool(mutation.get("caught"), "divergence mutation verdict") != mutation_passed:
        raise ValueError("reported divergence-mutation verdict disagrees with recomputation")

    local_passed = bool(
        divergence_passed
        and covariance_passed
        and radial_axis_support_passed
        and nontriviality_passed
        and perturbation_passed
        and mutation_passed
    )
    project_support_passed = bool(radial_axis_support_passed and axial_support_passed)
    public_passed = bool(local_passed and project_support_passed)
    if _strict_bool(audit.get("local_curl_covariance_preflight_passed"), "local preflight verdict") != local_passed:
        raise ValueError("reported local oscillatory verdict disagrees with recomputation")
    if _strict_bool(audit.get("project_support_preflight_passed"), "project support verdict") != project_support_passed:
        raise ValueError("reported project support verdict disagrees with recomputation")
    if _strict_bool(audit.get("public_oscillatory_preflight_passed"), "public preflight verdict") != public_passed:
        raise ValueError("reported public oscillatory verdict disagrees with recomputation")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "consumed_audit_schema": EXPECTED_AUDIT_SCHEMA,
        "consumed_agent2_head": EXPECTED_AGENT2_HEAD,
        "frozen_guards": dict(FROZEN_GUARDS),
        "failed_guards": tuple(failures),
        "independent_public_oscillatory_preflight_passed": public_passed,
        "correction_ingest_allowed": public_passed,
        "same_cycle_requested_stress_materialization_allowed": public_passed,
        "candidate_finite_head_mean_debt_materialization_allowed": public_passed,
        "requested_stress_actual_state_values_materialized": False,
        "finite_head_mean_debt_materialized": False,
        "real_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_rerun_allowed": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "surrogate_defect_used": False,
        "scope": (
            "independent oscillatory-handoff admission only; a PASS merely permits Agent 3 "
            "to materialize same-cycle requestedStress/mean debt next"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    audit = json.loads(args.audit.read_text(encoding="utf-8"))
    report = evaluate_independent_oscillatory_admission(audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
