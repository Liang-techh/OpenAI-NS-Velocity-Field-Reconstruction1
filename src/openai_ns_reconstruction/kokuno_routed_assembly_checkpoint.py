"""Kokuno Agent-5 routed assembly checkpoint.

This glue layer freezes the currently accepted routing decisions without inventing
the still-missing global core-to-heat splice:

* Agent 1 leading core is the callable inner field.
* Agent 2 oscillatory correction is retained at amplitude 0.125 and phase 0.0;
  the later training-selected phase pi/4 did not generalize on held-out data.
* Agent 3's +0.005 stress-shaped mean-swirl lift on the actual leading core was
  rejected and is not applied.
* Agent 1's source heat-exterior component is separately executable and
  serializable, but its c_inf=1 value is still an autonomous placeholder.
* Agent 4's independent three-stage core audit is bound as evidence only.

No global velocity is exposed because the source outer profile and compensated
core-to-heat matching are not complete.  The fixed full-domain normalized
Navier--Stokes gate remains 1e-3 and unassessed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_core_composite_checkpoint import KokunoCoreCompositeCandidate
from .kokuno_heat_exterior_profile import KokunoHeatExteriorProfile


SCHEMA = "kokuno-agent5-routed-assembly-checkpoint-v1"
TASK = "KOKUNO-A5-ROUTED-ASSEMBLY-CHECKPOINT-004"
PDE_THRESHOLD = 1.0e-3
DIVERGENCE_THRESHOLD = 1.0e-5
OSCILLATORY_AMPLITUDE = 0.125
OSCILLATORY_PHASE = 0.0

UPSTREAM = {
    "agent1_core": {
        "pr": 229,
        "head": "6f9cf3c4874087c0218a729df3c84a1fb45fc3a4",
        "role": "callable source-native finite leading core",
    },
    "agent1_heat_exterior": {
        "pr": 238,
        "head": "5c7ffdb9851e03d86394ed53c13043278473bb6b",
        "role": "source heat-factor and pure-swirl exterior component",
        "harvest": "exact upstream Git blobs; no Agent-5 mathematical rewrite",
        "c_inf_status": "autonomous_placeholder_pending_source_outer_assembly",
    },
    "agent2_phase_screen": {
        "pr": 239,
        "head": "b049601cd368d56a47db32e62ee4dfe28fae1f13",
        "role": "bounded phase-offset screen on actual leading core",
        "training_selected_phase": 0.7853981633974483,
        "retained_phase": OSCILLATORY_PHASE,
        "amplitude": OSCILLATORY_AMPLITUDE,
        "selected_generalizes_vs_phase0": False,
        "accepted_for_next_cycle": False,
        "heldout_finest_mean_rms": {
            "leading_only": 8.168462285885303,
            "phase0": 8.016485041340767,
            "training_selected_pi_over_4": 8.043024902731021,
        },
    },
    "agent3_core_mean_cycle": {
        "pr": 241,
        "head": "23b4860f2d1ea514aed4a1cce71e9c06b7b596a3",
        "role": "mean/radial correction screen on actual leading core",
        "selected_amplitude": 0.005,
        "accepted_for_next_cycle": False,
        "held_in_mean_defect_ratio": 2.18984,
        "held_out_total_mean_defect_ratio": 2.43062,
        "held_out_theta_ratio": 5.58634,
        "pressure_inclusive_raw_phase_mean_operator_ratio": 0.999477,
        "applied_to_routed_core": False,
        "source_exact_stress_to_velocity_lift": False,
    },
    "agent4_core_independent": {
        "pr": 242,
        "head": "9056f79cb83aed06df508b1d0787d71a23af59a2",
        "run": 35273292374,
        "artifact_id": 10519059689,
        "artifact_digest": "sha256:b77201c133ab4f6e1fc92ea9caff3df9a0571e6514e700b08a9d558d6bc25d1b",
        "role": "fresh held-out independent FD4 three-stage core audit",
        "validation_seed": 9172941,
        "finest_step": 0.005,
        "stage_ratios": {
            "leading_plus_oscillatory_over_leading": 0.9857381469915921,
            "after_rejected_mean_correction_over_leading_plus_oscillatory": 0.9994501078005282,
            "after_rejected_mean_correction_over_leading": 0.9851960972738397,
        },
        "local_worst_normalized_residual_max": 19.86354779674492,
        "local_worst_divergence_max": 6.244116819864587e-05,
        "local_reference_threshold_met": False,
        "local_divergence_reference_met": False,
        "formal_full_domain_gate_assessed": False,
        "pde_validated": False,
    },
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _receipt_sha(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("receipt_sha256", None)
    return hashlib.sha256(_canonical_json(copy).encode("utf-8")).hexdigest()


def routed_core_candidate() -> KokunoCoreCompositeCandidate:
    """Return the currently retained core-only velocity route.

    The numerical core is leading + oscillatory(amplitude=.125, phase=0).
    Agent-3's rejected mean correction is intentionally absent.
    """
    return KokunoCoreCompositeCandidate(
        oscillatory=KokunoCompleteCurlCorrection(
            amplitude=OSCILLATORY_AMPLITUDE,
            phase=OSCILLATORY_PHASE,
        )
    )


def global_velocity(*args: Any, **kwargs: Any) -> np.ndarray:
    """Fail closed rather than naively adding the unmatched core and heat pieces."""
    del args, kwargs
    raise RuntimeError(
        "global Kokuno velocity is unavailable: source outer profile and compensated "
        "core-to-heat matching are not complete"
    )


def write_routed_assembly_checkpoint(
    output_dir: str | Path = "artifacts/kokuno_agent5/routed_assembly_checkpoint_v4",
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    core = routed_core_candidate()
    heat = KokunoHeatExteriorProfile()

    core_path = core.save_json(output_dir / "core_candidate.json")
    heat_path = heat.save_json(output_dir / "heat_exterior_component.json")

    core_replay = KokunoCoreCompositeCandidate.load_json(core_path)
    heat_replay = KokunoHeatExteriorProfile.load_json(heat_path)
    if core_replay.sha256 != core.sha256:
        raise RuntimeError("core candidate save/load changed identity")
    if heat_replay.sha256 != heat.sha256:
        raise RuntimeError("heat-exterior save/load changed identity")

    core_probe = np.asarray(core_replay.velocity(0.1, 0.05, 0.0, 0.5), dtype=float)
    heat_probe = np.asarray(heat_replay.velocity(2.0, 0.75, 0.0, 0.4), dtype=float)
    if core_probe.shape != (3,) or heat_probe.shape != (3,):
        raise RuntimeError("component velocity API must return Cartesian [u,v,w]")
    if not np.all(np.isfinite(core_probe)) or not np.all(np.isfinite(heat_probe)):
        raise RuntimeError("component smoke probe became non-finite")
    if float(np.linalg.norm(core_probe)) <= 0.0 or float(np.linalg.norm(heat_probe)) <= 0.0:
        raise RuntimeError("component smoke probe unexpectedly collapsed to zero")

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "upstream": UPSTREAM,
        "routing": {
            "core_velocity": "Agent1 leading core + Agent2 complete curl",
            "oscillatory_amplitude": OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": OSCILLATORY_PHASE,
            "phase0_retained_after_agent2_holdout": True,
            "agent3_mean_correction_applied": False,
            "agent3_rejection_respected": True,
            "heat_exterior_combined_with_core": False,
            "reason_not_combined": (
                "source outer profile and compensated core-to-heat matching are incomplete"
            ),
        },
        "artifacts": {
            "core_candidate": {
                "path": str(core_path),
                "sha256": core.sha256,
                "file_sha256": _file_sha256(core_path),
                "api": "velocity(x,y,z,t)->Cartesian [u,v,w]",
            },
            "heat_exterior_component": {
                "path": str(heat_path),
                "sha256": heat.sha256,
                "file_sha256": _file_sha256(heat_path),
                "api": "velocity(x,y,z,t)->Cartesian [u,v,w] for X>0 exterior component",
            },
        },
        "smoke": {
            "core_probe_velocity": core_probe.tolist(),
            "heat_probe_velocity": heat_probe.tolist(),
            "visualization_or_delivery_only": True,
            "used_for_pde_acceptance": False,
        },
        "registered_gate": {
            "held_out_normalized_ns_residual_threshold": PDE_THRESHOLD,
            "divergence_reference": DIVERGENCE_THRESHOLD,
            "assessed": False,
            "passed": False,
            "why_unassessed": [
                "global source outer profile is not reconstructed",
                "compensated core-to-heat matching is not implemented",
                "Agent-3 current stress-to-velocity mean correction failed its own held-out gate",
                "compatible final fixed/restricted forcing contract is not attached",
            ],
        },
        "stage_state": {
            "leading_core_ready": True,
            "heat_exterior_component_ready": True,
            "core_to_heat_matching_complete": False,
            "leading_ready": False,
            "oscillatory_ready": True,
            "correction_ready": False,
            "core_velocity_export_ready": True,
            "velocity_export_ready": False,
            "independent_core_validation_ready": True,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "truth_boundary": {
            "paper_exact_reconstruction_claimed": False,
            "heat_component_is_global_matched_field": False,
            "rejected_mean_correction_promoted": False,
            "visual_success_implies_pde_success": False,
            "free_residual_defined_forcing_used": False,
        },
    }
    receipt["receipt_sha256"] = _receipt_sha(receipt)
    receipt_path = output_dir / "assembly_checkpoint.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def load_routed_assembly_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != SCHEMA:
        raise ValueError("unsupported Kokuno Agent-5 routed assembly schema")
    if payload.get("receipt_sha256") != _receipt_sha(payload):
        raise ValueError("routed assembly receipt SHA mismatch")
    if payload.get("upstream") != UPSTREAM:
        raise ValueError("routed assembly upstream provenance changed")
    if payload["registered_gate"][
        "held_out_normalized_ns_residual_threshold"
    ] != PDE_THRESHOLD:
        raise ValueError("registered 1e-3 PDE threshold changed")
    if payload["registered_gate"]["assessed"] or payload["registered_gate"]["passed"]:
        raise ValueError("global PDE gate was promoted without a global candidate")
    states = payload["stage_state"]
    if states["leading_ready"] or states["correction_ready"] or states["velocity_export_ready"]:
        raise ValueError("incomplete Kokuno assembly was promoted")
    if states["pde_validated"]:
        raise ValueError("PDE validation was promoted")
    if payload["routing"]["agent3_mean_correction_applied"]:
        raise ValueError("rejected Agent-3 mean correction was silently applied")
    if payload["routing"]["heat_exterior_combined_with_core"]:
        raise ValueError("unmatched heat exterior was silently combined with the core")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write the routed Kokuno core + heat-component integration checkpoint"
    )
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent5/routed_assembly_checkpoint_v4",
    )
    args = parser.parse_args()
    receipt = write_routed_assembly_checkpoint(args.output_dir)
    print(
        json.dumps(
            {
                "receipt_sha256": receipt["receipt_sha256"],
                "routing": receipt["routing"],
                "stage_state": receipt["stage_state"],
                "registered_gate": receipt["registered_gate"],
                "agent4_evidence": receipt["upstream"]["agent4_core_independent"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
