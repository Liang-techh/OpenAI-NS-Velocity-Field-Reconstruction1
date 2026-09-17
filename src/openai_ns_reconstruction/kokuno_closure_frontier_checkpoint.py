"""Agent-5 checkpoint for the current Kokuno closure frontier.

This module deliberately integrates *interfaces and routing decisions*, not new
mathematics.  It keeps the last replayable reference-stage velocity from Agent 5,
binds Agent 1's executable five-moment repair primitive, Agent 2's source
phase/covector contract, Agent 3's now-materialized-but-rejected signed-curl
correction cycle, and Agent 4's independent reference-stage seam audit.

No global leading field is invented here.  The signed Agent-3 correction is not
added to the routed velocity because its disjoint held-out finite-cycle screen
worsened the measured mean defect.  The source phase contract is likewise not
composed into the complete curl until Agent 1 -> source V/G mapping exists.  The
registered held-out full-momentum gate remains exactly 1e-3 and unassessed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_five_moment_repair import KokunoFiveMomentRepair
from .kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract
from .kokuno_reference_composite_checkpoint import KokunoReferenceCompositeCandidate
from .kokuno_signed_amplitude_curl_cycle import KokunoSignedAmplitudeCurlCorrection


SCHEMA = "kokuno-agent5-closure-frontier-checkpoint-v6"
TASK_ID = "KOKUNO-A5-CLOSURE-FRONTIER-CHECKPOINT-006"
REGISTERED_PDE_THRESHOLD = 1.0e-3
REGISTERED_DIVERGENCE_THRESHOLD = 1.0e-5
REFERENCE_COMPOSITE_SHA256 = (
    "f6e04bcd4ffbdb175b4336556ecf420905cc8edb374cbea4633f8feae34da2dd"
)

UPSTREAM = {
    "agent1_five_moment_repair": {
        "pr": 261,
        "head_sha": "fb53f6d984577386ac8beb02499690d099ded55d",
        "standard_run": 35282087830,
        "constrained_tests_passed": 229,
        "source_outer_target_discrepancy_supplied": False,
        "global_leading_profile_reconstructed": False,
    },
    "agent2_source_phase_contract": {
        "pr": 262,
        "head_sha": "71078654aadadcab1549189c06c2737e4595544f",
        "dedicated_run": 35283120571,
        "standard_run": 35283120773,
        "agent1_profile_to_source_V_G_mapping_completed": False,
        "complete_curl_velocity_changed": False,
    },
    "agent3_signed_amplitude_curl_cycle": {
        "pr": 263,
        "head_sha": "04f1335dd73eecffb95d18c523f1afb81ee89195",
        "dedicated_run": 35283445927,
        "standard_run": 35283445968,
        "artifact_id": 10524121068,
        "artifact_digest": (
            "sha256:5c01745ea7962f02be5cf319924c5929e186e8a57f894eb4afecaf234efcce8a"
        ),
        "held_in_mean_defect_ratio": 1.15246146,
        "held_out_mean_defect_ratio": 1.15068886,
        "held_out_theta_ratio": 1.08769840,
        "pressure_inclusive_raw_operator_ratio": 0.99992353,
        "correction_divergence_max": 1.950e-9,
        "public_velocity_correction_materialized": True,
        "accepted_for_next_cycle": False,
    },
    "agent4_independent_reference_audit": {
        "pr": 254,
        "head_sha": "7fcbb55b493b498b707045333122c5dd15a79aba",
        "dedicated_run": 35278817731,
        "standard_run": 35278817695,
        "artifact_id": 10521976161,
        "finest_step": 0.001,
        "finest_equal_stratum_time_sample_l2": 53.37532045280016,
        "finest_worst_residual_max": 141.70334110256638,
        "transition_worst_divergence_max": 0.18240010126668693,
        "post_core_worst_divergence_max": 0.35519874322564604,
        "formal_full_domain_gate_assessed": False,
    },
}

STATES = {
    "reference_stage_velocity_export_ready": True,
    "five_moment_repair_primitive_ready": True,
    "real_five_moment_target_ready": False,
    "source_phase_contract_ready": True,
    "agent1_to_source_V_G_mapping_ready": False,
    "source_phase_composed_into_complete_curl": False,
    "signed_curl_correction_materialized": True,
    "signed_curl_correction_accepted": False,
    "independent_reference_validation_ready": True,
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _five_moment_probe(repair: KokunoFiveMomentRepair) -> dict[str, Any]:
    zero = np.zeros(5, dtype=float)
    jacobian = repair.transformed_jacobian(
        zero, eta=0.23, P_star=1.15, f_eta=0.85
    )
    singular = np.linalg.svd(jacobian, compute_uv=False)
    rank = int(np.linalg.matrix_rank(jacobian))
    condition = float(singular[0] / singular[-1])
    return {
        "primitive_sha256": repair.sha256,
        "probe": {
            "eta": 0.23,
            "P_star": 1.15,
            "f_eta": 0.85,
            "linearized_rank": rank,
            "singular_values": [float(v) for v in singular],
            "condition_number": condition,
            "real_target_used": False,
        },
    }


def _phase_probe(contract: KokunoOscillatoryPhaseContract) -> dict[str, Any]:
    values = contract.evaluate_from_V(
        R=np.asarray([1.2]),
        theta=np.asarray([0.3]),
        Z=np.asarray([0.2]),
        pulse_v=np.asarray([0.1]),
        V=np.asarray([0.4]),
        V_R=np.asarray([0.2]),
        V_Z=np.asarray([0.05]),
        G=np.asarray([0.1]),
        G_R=np.asarray([0.02]),
        G_Z=np.asarray([0.03]),
    )
    n = np.asarray(values["n_Phi"], dtype=float)
    return {
        "contract_sha256": contract.sha256,
        "manufactured_interface_probe": {
            "Phi": float(np.asarray(values["Phi"]).reshape(-1)[0]),
            "n_Phi": [float(v) for v in n.reshape(-1, 3)[0]],
            "k": int(values["k"]),
            "k_m": int(values["k_m"]),
            "used_for_parameter_selection": False,
            "used_for_pde_acceptance": False,
        },
    }


def build_checkpoint() -> dict[str, Any]:
    candidate = KokunoReferenceCompositeCandidate()
    if candidate.sha256 != REFERENCE_COMPOSITE_SHA256:
        raise RuntimeError("reference composite identity drifted from Agent-5 v5 checkpoint")

    repair = KokunoFiveMomentRepair()
    phase = KokunoOscillatoryPhaseContract()
    # Importing the public Agent-3 type is part of the integration contract.  It
    # is intentionally not instantiated/routed because its exact-head cycle was
    # rejected on disjoint held-out mean-defect evidence.
    if not callable(getattr(KokunoSignedAmplitudeCurlCorrection, "velocity", None)):
        raise RuntimeError("Agent-3 signed-curl public velocity interface is unavailable")

    point = np.asarray([[0.1, 0.0, 0.05]], dtype=float)
    replay = np.asarray(candidate.at_points(point, 0.5), dtype=float)
    if replay.shape != (1, 3) or not np.all(np.isfinite(replay)):
        raise RuntimeError("reference-stage [u,v,w] replay failed")

    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "upstream": UPSTREAM,
        "fixed_gates": {
            "held_out_normalized_full_momentum_residual": REGISTERED_PDE_THRESHOLD,
            "divergence_max": REGISTERED_DIVERGENCE_THRESHOLD,
            "changed": False,
        },
        "reference_stage": {
            "candidate_sha256": candidate.sha256,
            "velocity_probe": {
                "point_xyz": point[0].tolist(),
                "time": 0.5,
                "uvw": replay[0].tolist(),
                "used_for_pde_acceptance": False,
            },
        },
        "agent1_five_moment_interface": _five_moment_probe(repair),
        "agent2_source_phase_interface": _phase_probe(phase),
        "routing": {
            "reference_stage_remains_replayable_not_promoted_global": True,
            "apply_agent3_signed_curl_correction": False,
            "reason_agent3_not_routed": (
                "materialized correction worsened held-in and disjoint held-out mean defect"
            ),
            "compose_source_phase_into_agent2_complete_curl": False,
            "reason_source_phase_not_composed": (
                "Agent-1 profile to source V/G mapping is not implemented"
            ),
            "run_formal_full_domain_gate": False,
            "reason_full_gate_not_run": (
                "global matched leading field and accepted finite correction cycle do not exist"
            ),
        },
        "states": STATES,
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "free_residual_defined_forcing_used": False,
            "surrogate_promoted_to_final_candidate": False,
            "kokuno_replay_called_independent_pde_validation": False,
            "threshold_relaxed": False,
        },
        "next_shortest_closure": [
            "Agent 1: supply the real outer/cone mismatch and complete matched core-to-heat leading field",
            "Agent 2: map Agent-1/base symbols to source V/G and compose the source phase/covector into the native-X complete curl",
            "Agent 3: add the missing independent covariance/stress column or resolve source sign/frame before another finite correction cycle",
            "Agent 4: assess the unchanged 1e-3 held-out full-domain gate only after a global composite exists",
        ],
    }
    payload["checkpoint_sha256"] = _sha(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise ValueError("checkpoint must be a JSON object")
    supplied = dict(payload)
    claimed = supplied.pop("checkpoint_sha256", None)
    if claimed != _sha(supplied):
        raise ValueError("checkpoint SHA mismatch")
    if supplied.get("schema") != SCHEMA or supplied.get("task_id") != TASK_ID:
        raise ValueError("checkpoint schema/task mismatch")
    if supplied.get("upstream") != UPSTREAM:
        raise ValueError("upstream evidence/provenance changed")
    if supplied.get("fixed_gates") != {
        "held_out_normalized_full_momentum_residual": REGISTERED_PDE_THRESHOLD,
        "divergence_max": REGISTERED_DIVERGENCE_THRESHOLD,
        "changed": False,
    }:
        raise ValueError("registered PDE/divergence gates changed")
    if supplied.get("states") != STATES:
        raise ValueError("frontier scientific states changed")
    routing = supplied.get("routing", {})
    if routing.get("apply_agent3_signed_curl_correction") is not False:
        raise ValueError("rejected Agent-3 correction cannot be routed")
    if routing.get("compose_source_phase_into_agent2_complete_curl") is not False:
        raise ValueError("unmapped source phase cannot be composed")
    if routing.get("run_formal_full_domain_gate") is not False:
        raise ValueError("formal full-domain gate cannot be marked run yet")
    truth = supplied.get("truth_boundary", {})
    required_truth = {
        "source_version_and_provenance_preserved": True,
        "free_residual_defined_forcing_used": False,
        "surrogate_promoted_to_final_candidate": False,
        "kokuno_replay_called_independent_pde_validation": False,
        "threshold_relaxed": False,
    }
    if truth != required_truth:
        raise ValueError("truth boundary changed")
    return payload


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_checkpoint(payload)


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    target = Path(output_dir)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError("output directory must be absent or empty")
    target.mkdir(parents=True, exist_ok=True)

    candidate = KokunoReferenceCompositeCandidate()
    repair = KokunoFiveMomentRepair()
    phase = KokunoOscillatoryPhaseContract()
    candidate_path = candidate.save_json(target / "reference_composite_candidate.json")
    repair_path = repair.save_json(target / "five_moment_repair_primitive.json")
    phase_path = phase.save_json(target / "source_phase_contract.json")

    checkpoint = build_checkpoint()
    checkpoint_path = target / "closure_frontier_checkpoint.json"
    checkpoint_path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # Fail closed on the same public loaders used downstream.
    if KokunoReferenceCompositeCandidate.load_json(candidate_path).sha256 != candidate.sha256:
        raise RuntimeError("reference candidate save/load identity changed")
    if KokunoFiveMomentRepair.load_json(repair_path).sha256 != repair.sha256:
        raise RuntimeError("five-moment primitive save/load identity changed")
    if KokunoOscillatoryPhaseContract.load_json(phase_path).sha256 != phase.sha256:
        raise RuntimeError("source phase contract save/load identity changed")
    load_checkpoint(checkpoint_path)
    return checkpoint


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent5/closure_frontier_checkpoint_v6",
    )
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
