"""Agent-5 routing checkpoint after Agent-4's bounded-coordinate audit.

This is a narrow integration overlay on the v22 unit-routing checkpoint.  It
binds Agent 4 PR #402's independent black-box PASS for Agent 2 PR #400's
repository-autonomous bounded coordinate interface.  It does not reinterpret
the still-rank-deficient one-band Agent-3 calibration as a correction PASS.

The formal full-domain gates remain unchanged: normalized momentum max/L2
<=1e-3 and divergence max/L2 <=1e-5.  The Agent-4 result bound here is only a
local coordinate/value/tangent implementation audit.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
from typing import Any

from .kokuno_bounded_coordinate_routing_checkpoint import (
    checkpoint_sha256,
    build_checkpoint as build_v22_checkpoint,
)

SCHEMA = "kokuno-agent5-bounded-coordinate-independent-routing-checkpoint-v23"
TASK_ID = "KOKUNO-A5-BOUNDED-COORDINATE-INDEPENDENT-ROUTING-023"

AGENT4_COORDINATE_AUDIT_RECEIPT = {
    "task_id": "KOKUNO-A4-BOUNDED-FAMILY-COORDINATE-INDEPENDENT-AUDIT-022",
    "source_pr": 402,
    "source_head": "a41358f1c11bb9af35505747a4fb2e524e63e1d9",
    "dedicated_run": 35346814135,
    "standard_run": 35346814165,
    "artifact_id": 10546873045,
    "artifact_digest": "sha256:aa31aed8f10c39773b062be45185664849e9214b99ccedec54ac0f6e4b01ea02",
    "seed": 9173131,
    "held_out_radii": 7,
    "held_out_angles": 24,
    "physical_labels": 6,
    "slow_bands": [5, 6, 7],
    "forward_steps": [0.04, 0.02, 0.01],
    "common_finest_covariance_response_relative_rms": 0.005000000000002316,
    "common_refinement_ratios": [1.9999999999999036, 1.9999999999991926],
    "common_quadratic_identity_relative_error": 1.3405568623477447e-15,
    "band_finest_covariance_response_relative_rms": 0.003408850955119712,
    "band_refinement_ratios": [1.9999999999992453, 1.999999999991075],
    "simultaneous_label_velocity_permutation_relative_rms": 1.0613488584753089e-16,
    "minimum_label_multiplier_at_l1_0125": 0.875,
    "misaligned_outer_band_response_relative_change": 1.9999999999999993,
    "local_guards": {
        "finest_covariance_response_relative_rms": 1.0e-2,
        "minimum_refinement_ratio": 1.8,
        "common_quadratic_identity_relative_error": 5.0e-13,
        "permutation_invariance_relative_rms": 5.0e-13,
        "minimum_label_multiplier": 0.875 - 5.0e-13,
        "misalignment_mutation_relative_change_min": 0.1,
    },
    "repository_bounded_coordinate_independent_preflight_passed": True,
    "actual_source_multiband_family_audited": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def build_checkpoint() -> dict[str, Any]:
    payload = deepcopy(build_v22_checkpoint())
    payload.pop("checkpoint_sha256", None)
    payload["schema"] = SCHEMA
    payload["task_id"] = TASK_ID
    payload["agent4_bounded_coordinate_audit_receipt"] = AGENT4_COORDINATE_AUDIT_RECEIPT

    payload["upstream"]["agent4"] = {
        "latest_pr": 402,
        "head_sha": AGENT4_COORDINATE_AUDIT_RECEIPT["source_head"],
        "dedicated_run": AGENT4_COORDINATE_AUDIT_RECEIPT["dedicated_run"],
        "standard_run": AGENT4_COORDINATE_AUDIT_RECEIPT["standard_run"],
        "artifact_id": AGENT4_COORDINATE_AUDIT_RECEIPT["artifact_id"],
        "artifact_digest": AGENT4_COORDINATE_AUDIT_RECEIPT["artifact_digest"],
        "source_compatible_partition_independent_preflight_passed": True,
        "bounded_coordinate_repository_independent_preflight_passed": True,
        "actual_multiband_physical_family_independent_audit_run": False,
        "consumed_in_executable_ancestry": False,
    }

    payload["states"]["repository_bounded_coordinate_independently_audited"] = True
    payload["states"]["actual_multiband_physical_family_independently_audited"] = False

    for stage in payload["pipeline_frontier"]:
        if stage["stage"] == "oscillatory_augmentation":
            stage.setdefault("evidence", []).append(
                "Agent-4 #402 independently clears the repository-autonomous delta_common/delta_band value/tangent handoff on a fresh three-band black-box family"
            )
            stage["status"] = "coordinate_contract_independently_cleared_waiting_for_actual_multiband_physical_family"
        if stage["stage"] == "mean_radial_corrections":
            stage.setdefault("evidence", []).append(
                "Agent-4 #402 clears the generic bounded-coordinate interface but does not audit an actual source/background multi-band family"
            )

    payload["routing"] = {
        "new_fact": (
            "Agent-4 #402 independently clears Agent-2 #400's repository-autonomous bounded-coordinate value/tangent interface under frozen local guards.  Together with Agent-3 #401, coordinate semantics and units are no longer the correction-lane blocker."
        ),
        "interpretation": (
            "The remaining blocker is physical, not bookkeeping: the one-band real calibration still has an inactive delta_band direction and 0/25 rank-two coverage.  Agent 2 must now instantiate a real multi-band positive-order/background family before Agent 3 can retest rank and Agent 4 can audit the actual-family receipt."
        ),
        "do_not_do": [
            "do not rerun the generic bounded-coordinate audit after #402 unless the interface changes",
            "do not call #402 an audit of the missing actual-source multi-band physical family",
            "do not reuse the physical-amplitude budget numerically in fractional coordinates without the #401 conversion",
            "do not call delta_band a second correction direction when only one slow band is active",
            "do not fabricate per-label free coefficients to manufacture rank",
            "do not rerun the finite correction cycle while required-node rank-two coverage is 0/25",
            "do not instantiate/export the final Kokuno candidate before leading, oscillatory and correction stages coexist",
            "do not weaken the fixed 1e-3 momentum or 1e-5 divergence gates",
            "do not use residual-defined free forcing",
        ],
        "shortest_next_closure": [
            "Agent 2: bind positive-order/background + mode data on the independently cleared source-compatible partition and produce a public Q-scaled by-beta/total xyz,t family with at least two active slow bands; retain source/autonomous provenance.",
            "Agent 3: consume that exact multi-band family through #401 and #393; require rank-two coverage, algebraic fit, aggregate budget, spacetime and radial d_z sigma_1 guards before correction materialization.",
            "Agent 4: independently audit that first actual-family tangent/rank receipt, then later own the formal held-out full-domain PDE gate on the frozen complete candidate.",
            "Agent 1: finish global leading assembly from #399, including remaining overlays, matched pressure and only the fixed/restricted forcing contract.",
            "Agent 5: instantiate the deterministic candidate artifact plus Python/MATLAB smoke only after leading + oscillatory + correction coexist in one executable ancestry.",
        ],
    }

    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = build_checkpoint()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent5/bounded_coordinate_independent_routing_checkpoint_v23/checkpoint.json",
    )
    args = parser.parse_args(argv)
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
