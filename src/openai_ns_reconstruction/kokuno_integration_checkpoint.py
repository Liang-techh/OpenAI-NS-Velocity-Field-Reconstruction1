"""Agent-5 integration checkpoint for the executable Kokuno reconstruction path.

This module is glue only.  It consumes the public interfaces delivered by the
Kokuno lanes and writes a reproducible checkpoint bundle without reimplementing
their mathematics:

* Agent 1: leading-axis profile plus native q/X/eta coordinates;
* Agent 2: complete-curl oscillatory correction (through Agent 3/4 consumers);
* Agent 3: measured phase-mean defect and compact radial-stress inverse;
* Agent 4: independent PCHIP/FD4 validation of that radial inverse.

The current base velocity used by Agents 3/4 remains the serialized capped
bipolar candidate and is explicitly an engineering bridge, not the Kokuno
leading field.  A complete 3-D Kokuno leading velocity and an applied finite
mean-velocity correction cycle do not yet exist in this ancestry, so the fixed
held-out normalized Navier--Stokes gate of 1e-3 is recorded but unassessed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_independent_radial_stress import run_real_candidate_preflight
from .kokuno_leading_axis_profiles import KokunoLeadingAxisProfile
from .kokuno_radial_stress import generate_current_candidate_radial_stress_report
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SCHEMA = "kokuno-agent5-integration-checkpoint-v2"
REGISTERED_PDE_THRESHOLD = 1.0e-3
BRIDGE_CANDIDATE_SHA256 = "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609"
DEFAULT_BRIDGE_CANDIDATE = "artifacts/bipolar_joint_capped/candidate.json"

UPSTREAM_HEADS = {
    "agent1_native_coordinates": "eabba21ca1f1375abbbf0b33537d66b82c5530e8",
    "agent2_complete_curl_core": "7e3a21e8f3b1950e28490bc35124571f285abe8a",
    "agent3_radial_stress_inverse": "cac3f4de3de2ca47364423e9a1b70407dd819526",
    "agent4_independent_radial_validator": "92755841265df56ff81b23ac820d2de4d2a6080f",
}

OPEN_UNCONSUMED_SIBLINGS = {
    "agent2_analytic_complete_curl_derivatives": {
        "pr": 222,
        "head": "930e84e333a829bee0e84d0f60f14e4e7ad26907",
        "reason": (
            "Useful for the later full-momentum composite, but not required to validate the "
            "current Agent-3 radial primitive; this checkpoint does not pretend it is in ancestry."
        ),
    }
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    """Hash a checkpoint payload while excluding its own digest field."""
    canonical = dict(payload)
    canonical.pop("receipt_sha256", None)
    return hashlib.sha256(_canonical_json(canonical).encode("utf-8")).hexdigest()


def _validate_checkpoint(payload: dict[str, Any]) -> None:
    if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
        raise ValueError("unsupported Kokuno Agent-5 integration checkpoint schema")
    if payload.get("upstream_heads") != UPSTREAM_HEADS:
        raise ValueError("Kokuno upstream-head provenance changed")
    if payload.get("open_unconsumed_siblings") != OPEN_UNCONSUMED_SIBLINGS:
        raise ValueError("Kokuno sibling-consumption provenance changed")

    bridge = payload.get("engineering_bridge", {})
    if bridge.get("candidate_sha256") != BRIDGE_CANDIDATE_SHA256:
        raise ValueError("engineering bridge candidate identity changed")
    if bridge.get("role") != "temporary_engineering_bridge_not_kokuno_leading":
        raise ValueError("engineering bridge was relabeled as Kokuno leading data")

    state = payload.get("stage_state", {})
    required_true = (
        "leading_axis_profile_ready",
        "native_coordinates_ready",
        "oscillatory_ready",
        "mean_defect_ready",
        "radial_inverse_ready",
        "radial_inverse_independent_validation_ready",
    )
    if not all(state.get(key) is True for key in required_true):
        raise ValueError("a prerequisite Kokuno checkpoint stage is not ready")
    required_false = (
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "pde_validated",
    )
    if not all(state.get(key) is False for key in required_false):
        raise ValueError("Kokuno scientific stage was promoted beyond current evidence")

    gate = payload.get("pde_gate", {})
    if gate.get("metric") != "held_out_normalized_full_momentum_residual":
        raise ValueError("Kokuno PDE gate metric changed")
    if gate.get("threshold") != REGISTERED_PDE_THRESHOLD:
        raise ValueError("Kokuno PDE threshold changed")
    if gate.get("assessed") is not False or gate.get("passed") is not False:
        raise ValueError("unassessed Kokuno PDE gate was promoted")

    truth = payload.get("truth_boundary", {})
    for key in ("pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary {key} was promoted")

    expected_sha = checkpoint_sha256(payload)
    if payload.get("receipt_sha256") != expected_sha:
        raise ValueError("Kokuno integration checkpoint SHA mismatch")


def load_integration_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    _validate_checkpoint(payload)
    return payload


def build_integration_checkpoint(
    *,
    candidate_path: str | Path = DEFAULT_BRIDGE_CANDIDATE,
    output_dir: str | Path = "artifacts/kokuno_agent5/integration_checkpoint_v2",
    time: float = 0.5,
    radial_count: int = 65,
    angular_count: int = 16,
) -> dict[str, Any]:
    """Run the currently available Agent-1--4 chain and write one checkpoint bundle."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = Path(candidate_path)

    leading = KokunoLeadingAxisProfile()
    coordinates = KokunoNativeSimilarityCoordinates(h=leading.h)
    leading_path = leading.save_json(output_dir / "leading_axis_profile.json")
    coordinates_path = coordinates.save_json(output_dir / "native_similarity_coordinates.json")

    # Exercise the native source relation on off-axis/off-time points so the
    # checkpoint records an executable coordinate bridge rather than metadata only.
    z_probe = np.array([-0.35, -0.12, 0.0, 0.18, 0.41], dtype=float)
    t_probe = np.array([0.31, 0.43, 0.55, 0.67, 0.79], dtype=float)
    native_relation_max_abs = float(np.max(np.abs(coordinates.relation_residual(z_probe, t_probe))))

    agent3_path = output_dir / "agent3_radial_stress_report.json"
    agent3 = generate_current_candidate_radial_stress_report(
        candidate_path=candidate_path,
        output=agent3_path,
        time=float(time),
        radial_count=int(radial_count),
        angular_count=int(angular_count),
    )
    if agent3.get("candidate_sha256") != BRIDGE_CANDIDATE_SHA256:
        raise ValueError("Agent-3 report is not tied to the frozen engineering bridge candidate")

    agent4_path = output_dir / "agent4_independent_radial_stress_report.json"
    agent4 = run_real_candidate_preflight(
        candidate_path=candidate_path,
        output=agent4_path,
        radial_counts=(33, 49, 65),
        angular_count=int(angular_count),
        time=float(time),
    )
    if agent4.get("candidate_sha256") != BRIDGE_CANDIDATE_SHA256:
        raise ValueError("Agent-4 report is not tied to the frozen engineering bridge candidate")

    projection = agent3["projection"]
    theta = agent3["theta_e2"]
    finest_independent_theta = agent4["levels"][-1]["theta_e2"]
    radial_independent_ready = bool(agent4["structural_preflight_passed"])

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream_heads": dict(UPSTREAM_HEADS),
        "open_unconsumed_siblings": OPEN_UNCONSUMED_SIBLINGS,
        "source_bundle": {
            "leading_axis_profile_path": str(leading_path),
            "leading_axis_profile_sha256": leading.sha256,
            "native_coordinates_path": str(coordinates_path),
            "native_coordinates_sha256": coordinates.sha256,
            "native_relation_probe_max_abs": native_relation_max_abs,
            "agent3_report_path": str(agent3_path),
            "agent4_report_path": str(agent4_path),
        },
        "engineering_bridge": {
            "candidate_path": str(candidate_path),
            "candidate_sha256": BRIDGE_CANDIDATE_SHA256,
            "role": "temporary_engineering_bridge_not_kokuno_leading",
            "used_for": "real phase-mean defect measurement until Agent 1 exposes a complete 3-D Kokuno leading velocity",
        },
        "stage_metrics": {
            "full_sampled_mean_defect_increment_rms": float(projection["full_mean_defect_increment_rms"]),
            "ring_theta_raw_rms": float(projection["ring_theta_raw_rms"]),
            "theta_gate_capture_rms_ratio": float(projection["theta_gate_capture_rms_ratio"]),
            "theta_radial_stress_rms": float(theta["stress_rms"]),
            "theta_radial_identity_agent3_rms": float(theta["radial_identity_fd_rms"]),
            "theta_radial_independent_stress_disagreement_normalized_rms": float(
                finest_independent_theta["stress_disagreement_normalized_rms"]
            ),
            "theta_radial_independent_identity_rms": float(finest_independent_theta["identity_rms"]),
            "theta_radial_independent_mutation_ratio": float(
                finest_independent_theta["mutation_to_baseline_identity_ratio"]
            ),
            "full_stage_residual_comparison_available": False,
        },
        "stage_state": {
            "leading_axis_profile_ready": True,
            "native_coordinates_ready": True,
            "leading_ready": False,
            "oscillatory_ready": True,
            "mean_defect_ready": True,
            "radial_inverse_ready": True,
            "radial_inverse_independent_validation_ready": radial_independent_ready,
            "correction_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
        "pde_gate": {
            "metric": "held_out_normalized_full_momentum_residual",
            "threshold": REGISTERED_PDE_THRESHOLD,
            "assessed": False,
            "passed": False,
            "reason": (
                "A complete 3-D Kokuno leading velocity, an applied finite mean/radial velocity correction cycle, "
                "and the compatible pressure/fixed-or-restricted forcing composite do not yet exist in this ancestry."
            ),
        },
        "blockers": [
            "Agent 1 still lacks the full radial leading profile and complete 3-D leading velocity.",
            "Agent 3 radial stress is a validated primitive but has not been converted into a public mean-velocity correction.",
            "The finite correction cycle has not been run on a complete Kokuno composite candidate.",
            "No compatible pressure plus fixed-or-preregistered restricted forcing has been attached to that absent composite.",
            "Only after those blockers close may Agent 4 assess the unchanged 1e-3 held-out normalized full-momentum gate.",
        ],
        "truth_boundary": {
            "kokuno_replay_is_independent_pde_validation": False,
            "engineering_bridge_is_kokuno_leading": False,
            "free_forcing_used": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    payload["receipt_sha256"] = checkpoint_sha256(payload)
    _validate_checkpoint(payload)

    receipt_path = output_dir / "integration_checkpoint.json"
    receipt_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the current Kokuno Agent-1--4 integration checkpoint")
    parser.add_argument("--candidate", default=DEFAULT_BRIDGE_CANDIDATE)
    parser.add_argument("--output-dir", default="artifacts/kokuno_agent5/integration_checkpoint_v2")
    parser.add_argument("--time", type=float, default=0.5)
    parser.add_argument("--radial-count", type=int, default=65)
    parser.add_argument("--angular-count", type=int, default=16)
    args = parser.parse_args()
    payload = build_integration_checkpoint(
        candidate_path=args.candidate,
        output_dir=args.output_dir,
        time=args.time,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
    )
    print(
        json.dumps(
            {
                "receipt_sha256": payload["receipt_sha256"],
                "stage_state": payload["stage_state"],
                "stage_metrics": payload["stage_metrics"],
                "pde_gate": payload["pde_gate"],
                "blockers": payload["blockers"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
