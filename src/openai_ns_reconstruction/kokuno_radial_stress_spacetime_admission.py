"""Agent-3 admission for Agent-4's spacetime radial-stress generalization.

Agent 4 PR #600 independently rebuilt the admitted oscillatory self-defect with
Cartesian FD8 at two new off-grid spacetime states, reconstructed the compact
moment-complement radial stress, and checked the radial operator with a separate
FD8 derivative.  This module pins that exact execution contract and recomputes
all scientific guards from the returned report instead of trusting a PASS bit.

The narrow promotion is correction-side only: Agent 3 may reuse the already
admitted compact/radial operator without retuning it merely because the cycle
state moves away from the original #595 anchor.  A future cycle still must feed
an actual same-cycle full composite defect satisfying the same support/moment
contract.  This module does not create Agent-1 leading/cross data, pressure,
forcing, a signed inverse, a correction velocity, or a Navier--Stokes PASS.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .kokuno_agent4_radial_stress_spacetime_generalization import (
    audit_spacetime_generalization,
)

TASK = "KOKUNO-A3-RADIAL-STRESS-SPACETIME-ADMISSION-043"
SCHEMA = "kokuno-a3-radial-stress-spacetime-admission-v1"
EXPECTED_AUDIT_SCHEMA = "kokuno-a4-radial-stress-spacetime-generalization-v2"

PARENT_AGENT3_PR = 598
PARENT_AGENT3_HEAD = "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233"
ACTUAL_AUDITOR_AGENT4_PR = 600
ACTUAL_AUDITOR_AGENT4_HEAD = "3c9934ba390e28597ad32fcaf0e320af8da990ad"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
PRIOR_AGENT4_PR = 595
PRIOR_AGENT4_HEAD = "19296acad4f84ba05d3b095bb8c01bf5a2c95891"
AUDIT_WORKFLOW_RUN_ID = 35434164786
AUDIT_STANDARD_WORKFLOW_RUN_ID = 35434164785
AUDIT_ARTIFACT_ID = 10581683366
AUDIT_ARTIFACT_ZIP_DIGEST = (
    "sha256:c6a7617d45abf893e07708ea3fb13a963d7c3e3d3621f14300a85d007685cebc"
)

RADIAL_COUNTS = (25, 49, 97)
ANGULAR_COUNT = 16
NEW_SPACETIME_STATES = ((0.37, -0.31), (0.63, 0.31))
PRIOR_ANCHOR_STATE = (0.50, 0.08)
CARTESIAN_FD8_STEP = 0.0025
NU = 0.01
FROZEN_GUARDS = {
    "finest_relative_rms_max": 2.0e-2,
    "finest_relative_max_max": 5.0e-2,
    "minimum_refinement_ratio": 2.0,
    "moment_complement_relative_max": 1.0e-10,
    "edge_relative_max": 1.0e-8,
    "sign_flip_mutation_min": 5.0e-1,
    "nontrivial_self_defect_rms_min": 1.0e-8,
}
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


def _approx_equal(
    left: float,
    right: float,
    *,
    rel: float = 1.0e-12,
    abs_: float = 1.0e-15,
) -> bool:
    return math.isclose(float(left), float(right), rel_tol=rel, abs_tol=abs_)


def _require_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _require_exact_protocol(report: Mapping[str, Any]) -> None:
    if report.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("Agent-4 spacetime radial-stress schema changed")

    provenance = _require_mapping(report.get("provenance"), "provenance")
    expected_provenance = {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "prior_agent4_pr": PRIOR_AGENT4_PR,
        "prior_agent4_head": PRIOR_AGENT4_HEAD,
    }
    for key, expected in expected_provenance.items():
        if provenance.get(key) != expected:
            raise ValueError(f"Agent-4 provenance changed: {key}")
    anchor = _require_mapping(provenance.get("prior_anchor_state"), "prior_anchor_state")
    if not _approx_equal(float(anchor.get("t", math.nan)), PRIOR_ANCHOR_STATE[0]):
        raise ValueError("prior anchor time changed")
    if not _approx_equal(float(anchor.get("z", math.nan)), PRIOR_ANCHOR_STATE[1]):
        raise ValueError("prior anchor z changed")

    repair = _require_mapping(report.get("protocol_repair"), "protocol_repair")
    if int(repair.get("first_attempt_workflow", -1)) != 35433365648:
        raise ValueError("pre-result stencil-feasibility workflow identity changed")
    if repair.get("first_attempt_scientific_values_emitted") is not False:
        raise ValueError("pre-result repair can only be admitted when no scientific values were emitted")

    protocol = _require_mapping(report.get("protocol"), "protocol")
    if tuple(int(value) for value in protocol.get("radial_counts", ())) != RADIAL_COUNTS:
        raise ValueError("radial resolution ladder changed")
    if int(protocol.get("angular_count", -1)) != ANGULAR_COUNT:
        raise ValueError("angular count changed")
    if not _approx_equal(float(protocol.get("cartesian_fd8_step", math.nan)), CARTESIAN_FD8_STEP):
        raise ValueError("Cartesian FD8 step changed")
    if not _approx_equal(float(protocol.get("nu", math.nan)), NU):
        raise ValueError("viscosity changed")
    if protocol.get("cartesian_derivative") != "centered FD8, independently implemented":
        raise ValueError("Cartesian derivative implementation label changed")
    if protocol.get("radial_derivative") != "separate centered FD8":
        raise ValueError("radial derivative implementation label changed")

    states = protocol.get("new_spacetime_states")
    if not isinstance(states, list) or len(states) != len(NEW_SPACETIME_STATES):
        raise ValueError("new spacetime-state set changed")
    for row, expected in zip(states, NEW_SPACETIME_STATES):
        mapping = _require_mapping(row, "new_spacetime_state")
        if not _approx_equal(float(mapping.get("t", math.nan)), expected[0]):
            raise ValueError("new spacetime time changed")
        if not _approx_equal(float(mapping.get("z", math.nan)), expected[1]):
            raise ValueError("new spacetime z changed")

    protocol_anchor = _require_mapping(protocol.get("prior_anchor_state"), "protocol.prior_anchor_state")
    if not _approx_equal(float(protocol_anchor.get("t", math.nan)), PRIOR_ANCHOR_STATE[0]):
        raise ValueError("protocol anchor time changed")
    if not _approx_equal(float(protocol_anchor.get("z", math.nan)), PRIOR_ANCHOR_STATE[1]):
        raise ValueError("protocol anchor z changed")

    guards = _require_mapping(protocol.get("guards"), "protocol.guards")
    if set(guards) != set(FROZEN_GUARDS):
        raise ValueError("frozen scientific guard set changed")
    for key, expected in FROZEN_GUARDS.items():
        if not _approx_equal(float(guards[key]), expected):
            raise ValueError(f"frozen scientific guard changed: {key}")


def _recompute_channel(
    levels: list[Mapping[str, Any]],
    convergence: Mapping[str, Any],
    channel: str,
    exponent: int,
) -> dict[str, Any]:
    errors: list[float] = []
    rows: list[Mapping[str, Any]] = []
    for level in levels:
        row = _require_mapping(level.get(channel), f"{channel} level")
        if int(row.get("exponent", -1)) != exponent:
            raise ValueError(f"{channel} exponent changed")
        errors.append(float(row["operator_relative_rms"]))
        rows.append(row)
    ratios = [errors[index] / max(errors[index + 1], 1.0e-300) for index in range(len(errors) - 1)]
    finest = rows[-1]

    declared = _require_mapping(convergence.get(channel), f"convergence.{channel}")
    declared_errors = [float(value) for value in declared.get("relative_rms_by_radial_count", ())]
    declared_ratios = [float(value) for value in declared.get("refinement_ratios", ())]
    if len(declared_errors) != len(errors) or any(
        not _approx_equal(a, b) for a, b in zip(declared_errors, errors)
    ):
        raise ValueError(f"{channel} convergence RMS ladder disagrees with levels")
    if len(declared_ratios) != len(ratios) or any(
        not _approx_equal(a, b) for a, b in zip(declared_ratios, ratios)
    ):
        raise ValueError(f"{channel} refinement ratios disagree with levels")

    derived = {
        "relative_rms_by_radial_count": tuple(errors),
        "refinement_ratios": tuple(ratios),
        "finest_relative_rms": errors[-1],
        "finest_relative_max": float(finest["operator_relative_max"]),
        "finest_moment_complement_relative": float(finest["moment_complement_relative"]),
        "finest_edge_relative": float(finest["edge_relative"]),
        "finest_sign_flip_mutation_relative_rms": float(finest["sign_flip_mutation_relative_rms"]),
        "finest_weighted_moment": float(finest["weighted_moment"]),
    }
    for key in (
        "finest_relative_max",
        "finest_moment_complement_relative",
        "finest_edge_relative",
        "finest_sign_flip_mutation_relative_rms",
        "finest_weighted_moment",
    ):
        if not _approx_equal(float(declared[key]), float(derived[key])):
            raise ValueError(f"{channel} declared {key} disagrees with finest level")

    failed: list[str] = []
    if derived["finest_relative_rms"] > FROZEN_GUARDS["finest_relative_rms_max"]:
        failed.append("finest_relative_rms")
    if derived["finest_relative_max"] > FROZEN_GUARDS["finest_relative_max_max"]:
        failed.append("finest_relative_max")
    if min(ratios) < FROZEN_GUARDS["minimum_refinement_ratio"]:
        failed.append("minimum_refinement_ratio")
    if derived["finest_moment_complement_relative"] > FROZEN_GUARDS[
        "moment_complement_relative_max"
    ]:
        failed.append("moment_complement")
    if derived["finest_edge_relative"] > FROZEN_GUARDS["edge_relative_max"]:
        failed.append("edge_relative")
    if derived["finest_sign_flip_mutation_relative_rms"] < FROZEN_GUARDS[
        "sign_flip_mutation_min"
    ]:
        failed.append("sign_flip_mutation")
    derived["failed_guards"] = tuple(failed)
    derived["passed"] = not failed
    return derived


def _recompute_state(state: Mapping[str, Any], expected_index: int) -> dict[str, Any]:
    if int(state.get("state_index", -1)) != expected_index:
        raise ValueError("state index changed")
    expected_t, expected_z = NEW_SPACETIME_STATES[expected_index]
    if not _approx_equal(float(state.get("t", math.nan)), expected_t):
        raise ValueError("state time changed")
    if not _approx_equal(float(state.get("z", math.nan)), expected_z):
        raise ValueError("state z changed")

    raw_levels = state.get("levels")
    if not isinstance(raw_levels, list) or len(raw_levels) != len(RADIAL_COUNTS):
        raise ValueError("state radial levels changed")
    levels: list[Mapping[str, Any]] = []
    for raw_level, expected_count in zip(raw_levels, RADIAL_COUNTS):
        level = _require_mapping(raw_level, "state level")
        if int(level.get("radial_count", -1)) != expected_count:
            raise ValueError("state radial count changed")
        if float(level.get("boundary_margin", 0.0)) <= float(level.get("stencil_reach", math.inf)):
            raise ValueError("FD8 Cartesian stencil no longer fits radial support margin")
        levels.append(level)

    convergence = _require_mapping(state.get("convergence"), "state.convergence")
    theta = _recompute_channel(levels, convergence, "theta_e2", 2)
    axial = _recompute_channel(levels, convergence, "axial_e1", 1)
    finest_self_defect_rms = float(levels[-1]["self_defect_rms"])

    failed: list[str] = []
    if finest_self_defect_rms < FROZEN_GUARDS["nontrivial_self_defect_rms_min"]:
        failed.append("nontrivial_self_defect_rms")
    failed.extend(f"theta_e2:{name}" for name in theta["failed_guards"])
    failed.extend(f"axial_e1:{name}" for name in axial["failed_guards"])

    declared_failed = tuple(str(value) for value in state.get("failed_guards", ()))
    if declared_failed != tuple(failed):
        raise ValueError("state declared failed_guards disagrees with recomputed guards")
    if bool(state.get("state_passed")) is not (not failed):
        raise ValueError("state declared PASS disagrees with recomputed guards")

    return {
        "state_index": expected_index,
        "t": expected_t,
        "z": expected_z,
        "finest_self_defect_rms": finest_self_defect_rms,
        "theta_e2": theta,
        "axial_e1": axial,
        "failed_guards": tuple(failed),
        "passed": not failed,
    }


def admit_spacetime_radial_stress_report(report: Mapping[str, Any]) -> dict[str, Any]:
    """Fail closed unless the exact #600 scientific contract still passes."""
    _require_exact_protocol(report)
    raw_states = report.get("states")
    if not isinstance(raw_states, list) or len(raw_states) != len(NEW_SPACETIME_STATES):
        raise ValueError("Agent-4 state report set changed")
    states = [_recompute_state(_require_mapping(row, "state"), index) for index, row in enumerate(raw_states)]
    failed = tuple(
        f"state{state['state_index']}:{name}"
        for state in states
        for name in state["failed_guards"]
    )
    declared_failed = tuple(str(value) for value in report.get("failed_guards", ()))
    if declared_failed != failed:
        raise ValueError("Agent-4 declared global failed_guards disagrees with recomputed guards")
    if bool(report.get("spacetime_radial_stress_generalization_passed")) is not (not failed):
        raise ValueError("Agent-4 declared spacetime PASS disagrees with recomputed guards")
    if failed:
        raise ValueError("Agent-4 #600 spacetime radial-stress audit does not pass frozen guards")

    truth = _require_mapping(report.get("truth_boundary"), "truth_boundary")
    required_true = (
        "public_black_box_velocity_consumed",
        "public_black_box_time_derivative_consumed",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"independent spacetime audit lost required evidence: {key}")
    required_false = (
        "agent3_self_defect_operator_used",
        "agent3_compact_stress_constructor_used",
        "agent3_radial_derivative_used",
        "pressure_fitted",
        "forcing_fitted",
        "full_same_cycle_composite_requested_stress_materialized",
        "public_velocity_correction_materialized",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    )
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"spacetime operator audit cannot promote {key}")
    if not _approx_equal(float(truth.get("final_normalized_momentum_gate", math.nan)), FINAL_NORMALIZED_MOMENTUM_GATE):
        raise ValueError("final normalized momentum gate changed")
    if not _approx_equal(float(truth.get("final_normalized_divergence_gate", math.nan)), FINAL_NORMALIZED_DIVERGENCE_GATE):
        raise ValueError("final normalized divergence gate changed")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "source_agent4_pr": ACTUAL_AUDITOR_AGENT4_PR,
        "source_agent4_head": ACTUAL_AUDITOR_AGENT4_HEAD,
        "source_agent4_dedicated_workflow": AUDIT_WORKFLOW_RUN_ID,
        "source_agent4_standard_workflow": AUDIT_STANDARD_WORKFLOW_RUN_ID,
        "source_agent4_artifact_id": AUDIT_ARTIFACT_ID,
        "source_agent4_artifact_zip_digest": AUDIT_ARTIFACT_ZIP_DIGEST,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "states": states,
        "failed_guards": failed,
        "independent_spacetime_radial_stress_generalization_admitted": True,
        "compact_radial_stress_operator_independently_admitted": True,
        "full_composite_radial_stress_execution_allowed_when_actual_defect_available": True,
        "finite_correction_cycle_radial_operator_reuse_allowed_without_retuning": True,
        "radial_operator_reuse_scope": (
            "same admitted compact/moment operator; every future cycle must still provide an actual "
            "same-cycle defect and satisfy support/moment/regularity and gain guards"
        ),
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "real_full_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def evaluate_spacetime_radial_stress_admission() -> dict[str, Any]:
    """Re-run Agent 4's public-only audit and apply the Agent-3 admission gate."""
    return admit_spacetime_radial_stress_report(audit_spacetime_generalization())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate_spacetime_radial_stress_admission()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
