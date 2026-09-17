"""Agent-5 integration checkpoint for the executable Kokuno reference stage.

This module performs one deliberately narrow integration step.  It combines
Agent 1's source-derived *reference continuation* with the already retained
Agent-2 complete-curl correction under one public Cartesian ``[u,v,w]`` API,
then binds the exact-head evidence from Agents 2--4 into a deterministic
checkpoint receipt.

The reference continuation is not the final Kokuno leading field.  Independent
Agent-4 validation finds a large transition/post-core divergence and momentum
obstruction, while Agent 3 has only established rank-one signed-covariance
*realizability* and has not materialized that correction as a public velocity.
Accordingly this artifact is exportable as an intermediate research stage but
fails closed for ``leading_ready``, full ``velocity_export_ready`` and
``pde_validated``.  The registered 1e-3 PDE gate is not changed or assessed.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate


CANDIDATE_SCHEMA = "kokuno-reference-composite-candidate-v1"
CHECKPOINT_SCHEMA = "kokuno-agent5-reference-integration-checkpoint-v5"
TASK_ID = "KOKUNO-A5-REFERENCE-INTEGRATION-CHECKPOINT-005"
ROUTED_AMPLITUDE = 0.125
ROUTED_PHASE = 0.0
REGISTERED_PDE_THRESHOLD = 1.0e-3
REGISTERED_DIVERGENCE_THRESHOLD = 1.0e-5

# Exact-head evidence is copied from successful upstream workflow artifacts, not
# recomputed or selected by this integration lane.  Every value is tied to an
# immutable head SHA plus workflow artifact digest so later stages can tell
# evidence binding from independent validation.
UPSTREAM_EVIDENCE: dict[str, Any] = {
    "agent1_reference_continuation": {
        "pr": 251,
        "head_sha": "035ef85a2d48d00e907afcbace20932542ee0b6a",
        "standard_test_run": 35276765796,
        "candidate_sha256": "5042bbb390ec2d2e26ddb8133d887587b8321a71880d0a0b478cc9151843c2d7",
        "source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
        "corrected_release": "zenodo:22678406",
        "stage_only": True,
    },
    "agent2_reference_oscillatory_screen": {
        "pr": 252,
        "head_sha": "defb572e954e2e225b14588403823a1757458fcf",
        "standard_test_run": 35277801428,
        "dedicated_run": 35277801462,
        "artifact_id": 10521706314,
        "artifact_digest": "sha256:a69b644342c65d65daaa9eaeb2aaa7895980a5169ee944c5d6d823c64f83b79b",
        "validation_seed": 9172981,
        "finest_step": 0.005,
        "routed_amplitude": ROUTED_AMPLITUDE,
        "routed_phase": ROUTED_PHASE,
        "inner_mean_residual_rms_ratio": 0.9698914432979021,
        "post_core_mean_residual_rms_ratio": 0.9998713322004918,
        "all_strata_mean_residual_rms_ratio": 0.9971060818942786,
        "post_core_worst_composite_divergence_max": 0.612132220610281,
        "post_core_worst_composite_residual_max": 145.4876121326331,
        "parameter_selection_performed": False,
        "formal_full_domain_gate_assessed": False,
    },
    "agent3_signed_covariance_inverse": {
        "pr": 253,
        "head_sha": "17f15998218f885b072e2e483de4aa26f2f80268",
        "standard_test_run": 35278036364,
        "dedicated_run": 35278036467,
        "artifact_id": 10521700003,
        "artifact_digest": "sha256:b634cb07afcf07f3698fd8cda0a67573e89090af618016b093aabdf2e551b853",
        "input_leading_core_sha256": "936768138f959e969caa5a8537647022f17818e65304ba7f60200d067c9bff67",
        "active_target_nodes": 26,
        "unresolved_active_nodes": 0,
        "rank1_realizable_on_active_nodes": True,
        "delta_amplitude_max_abs": 0.006635614730692422,
        "within_existing_agent2_absolute_amplitude_bound": True,
        "public_velocity_correction_materialized": False,
        "residual_reduction_claimed": False,
        "formal_full_domain_gate_assessed": False,
    },
    "agent4_independent_reference_audit": {
        "pr": 254,
        "head_sha": "7fcbb55b493b498b707045333122c5dd15a79aba",
        "standard_test_run": 35278817695,
        "dedicated_run": 35278817731,
        "artifact_id": 10521976161,
        "artifact_digest": "sha256:331ed0cb2524d3783ecaaf8db5aaaaec449f0e0622397a97d48bca422f5d56f2",
        "validation_seed": 9172991,
        "finest_step": 0.001,
        "candidate_sha256": "5042bbb390ec2d2e26ddb8133d887587b8321a71880d0a0b478cc9151843c2d7",
        "finest_equal_stratum_time_sample_l2": 53.37532045280016,
        "finest_worst_residual_max": 141.70334110256638,
        "transition_worst_divergence_max": 0.18240010126668693,
        "post_core_worst_divergence_max": 0.35519874322564604,
        "post_core_radial_component_rms": 70.1736988132618,
        "formal_full_domain_gate_assessed": False,
    },
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_axis(values: Any, name: str) -> np.ndarray:
    axis = np.asarray(values, dtype=float)
    if axis.ndim != 1 or axis.size == 0 or not np.all(np.isfinite(axis)):
        raise ValueError(f"{name} must be a nonempty finite one-dimensional axis")
    return axis


@dataclass(frozen=True)
class KokunoReferenceCompositeCandidate:
    """Callable Agent-1 reference stage plus the frozen routed Agent-2 curl.

    This class is intentionally fixed to the routed amplitude/phase.  It is an
    intermediate stage artifact, not the final global Kokuno candidate.
    """

    amplitude: float = ROUTED_AMPLITUDE
    phase: float = ROUTED_PHASE
    _leading: KokunoReferenceContinuationCandidate = field(
        init=False, repr=False, compare=False
    )
    _correction: KokunoCompleteCurlCorrection = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if float(self.amplitude) != ROUTED_AMPLITUDE or float(self.phase) != ROUTED_PHASE:
            raise ValueError("reference checkpoint is frozen to routed amplitude=.125, phase=0")
        object.__setattr__(self, "amplitude", ROUTED_AMPLITUDE)
        object.__setattr__(self, "phase", ROUTED_PHASE)
        object.__setattr__(self, "_leading", KokunoReferenceContinuationCandidate())
        object.__setattr__(
            self,
            "_correction",
            KokunoCompleteCurlCorrection(amplitude=ROUTED_AMPLITUDE, phase=ROUTED_PHASE),
        )
        if self._leading.sha256 != UPSTREAM_EVIDENCE["agent1_reference_continuation"]["candidate_sha256"]:
            raise RuntimeError("reference-continuation identity drifted from bound Agent-1 evidence")

    @property
    def leading(self) -> KokunoReferenceContinuationCandidate:
        return self._leading

    @property
    def correction(self) -> KokunoCompleteCurlCorrection:
        return self._correction

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        leading = np.asarray(self._leading.velocity(x, y, z, t), dtype=float)
        correction = np.asarray(self._correction.velocity(x, y, z, t), dtype=float)
        if leading.shape != correction.shape or leading.shape[-1] != 3:
            raise RuntimeError("leading and oscillatory velocity layouts do not match")
        result = leading + correction
        if not np.all(np.isfinite(result)):
            raise RuntimeError("reference composite velocity became non-finite")
        return result

    __call__ = velocity

    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        # No autonomous pressure fit is introduced by the routed curl screen.
        return self._leading.pressure(x, y, z, t)

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = np.asarray(points_xyz, dtype=float)
        if points.ndim < 2 or points.shape[-1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points_xyz must be finite with shape (...,3)")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def grid(self, x: Any, y: Any, z: Any, times: Any) -> np.ndarray:
        x_axis = _finite_axis(x, "x")
        y_axis = _finite_axis(y, "y")
        z_axis = _finite_axis(z, "z")
        time_axis = _finite_axis(times, "times")
        Xg, Yg, Zg = np.meshgrid(x_axis, y_axis, z_axis, indexing="ij")
        output = np.empty(
            (time_axis.size, x_axis.size, y_axis.size, z_axis.size, 3), dtype=float
        )
        for index, time in enumerate(time_axis):
            output[index] = self.velocity(Xg, Yg, Zg, float(time))
        return output

    def to_payload(self) -> dict[str, Any]:
        correction_metadata = self._correction.metadata()
        return {
            "schema": CANDIDATE_SCHEMA,
            "composition": "Agent-1 reference continuation + routed Agent-2 complete curl",
            "leading": {
                "schema": self._leading.to_payload()["schema"],
                "sha256": self._leading.sha256,
                "source": self._leading.to_payload()["source"],
            },
            "oscillatory": {
                "family": correction_metadata["family"],
                "parameters": correction_metadata["parameters"],
                "source": correction_metadata["source"],
                "autonomous_choices": correction_metadata["autonomous_choices"],
            },
            "routing": {
                "amplitude": self.amplitude,
                "phase": self.phase,
                "agent2_head_sha": UPSTREAM_EVIDENCE["agent2_reference_oscillatory_screen"]["head_sha"],
                "parameter_selection_performed_here": False,
            },
            "api": {
                "velocity_layout": "Cartesian [u,v,w] on final component axis",
                "grid_layout": "(time,x,y,z,component)",
                "pressure": "Agent-1 reference-continuation pressure unchanged",
            },
            "truth_boundary": {
                "reference_stage_only": True,
                "reference_composite_velocity_export_ready": True,
                "cone_modulation_completed": False,
                "five_moment_repair_completed": False,
                "heat_compensation_completed": False,
                "core_to_heat_matching_completed": False,
                "mean_radial_correction_materialized": False,
                "global_leading_profile_reconstructed": False,
                "formal_full_domain_pde_gate_assessed": False,
                "pde_validated": False,
                "paper_exact": False,
                "openai_field_identified": False,
                "blowup_proved": False,
            },
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoReferenceCompositeCandidate":
        if not isinstance(payload, dict):
            raise ValueError("candidate payload must be a JSON object")
        candidate = cls()
        canonical = json.loads(_canonical_json(candidate.to_payload()))
        supplied = dict(payload)
        supplied_sha = supplied.pop("sha256", None)
        if supplied != canonical:
            raise ValueError("reference-composite payload metadata changed")
        if supplied_sha is not None and supplied_sha != candidate.sha256:
            raise ValueError("reference-composite SHA mismatch")
        return candidate

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoReferenceCompositeCandidate":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def build_checkpoint(candidate: KokunoReferenceCompositeCandidate | None = None) -> dict[str, Any]:
    candidate = KokunoReferenceCompositeCandidate() if candidate is None else candidate
    if candidate.leading.sha256 != UPSTREAM_EVIDENCE["agent4_independent_reference_audit"]["candidate_sha256"]:
        raise RuntimeError("Agent-4 evidence is not bound to this reference-continuation identity")

    a2 = UPSTREAM_EVIDENCE["agent2_reference_oscillatory_screen"]
    a3 = UPSTREAM_EVIDENCE["agent3_signed_covariance_inverse"]
    a4 = UPSTREAM_EVIDENCE["agent4_independent_reference_audit"]
    large_sampled_divergence = bool(
        a4["post_core_worst_divergence_max"] > REGISTERED_DIVERGENCE_THRESHOLD
        and a4["transition_worst_divergence_max"] > REGISTERED_DIVERGENCE_THRESHOLD
    )
    large_sampled_momentum = bool(
        a4["finest_worst_residual_max"] > REGISTERED_PDE_THRESHOLD
    )

    return {
        "schema": CHECKPOINT_SCHEMA,
        "task_id": TASK_ID,
        "candidate_sha256": candidate.sha256,
        "upstream_evidence": UPSTREAM_EVIDENCE,
        "fixed_gates": {
            "held_out_normalized_full_momentum_residual": REGISTERED_PDE_THRESHOLD,
            "divergence_max": REGISTERED_DIVERGENCE_THRESHOLD,
            "changed": False,
        },
        "integration_findings": {
            "reference_continuation_identity_consistent_across_agents_1_2_4": True,
            "retained_oscillation_improves_inner_mean_rms": bool(
                a2["inner_mean_residual_rms_ratio"] < 1.0
            ),
            "retained_oscillation_materially_repairs_post_core_seam": False,
            "agent2_post_core_ratio": a2["post_core_mean_residual_rms_ratio"],
            "agent4_large_sampled_transition_and_post_core_divergence": large_sampled_divergence,
            "agent4_large_sampled_momentum_obstruction": large_sampled_momentum,
            "agent3_rank1_signed_covariance_realizability_ready": bool(
                a3["rank1_realizable_on_active_nodes"]
                and a3["within_existing_agent2_absolute_amplitude_bound"]
            ),
            "agent3_public_velocity_correction_materialized": bool(
                a3["public_velocity_correction_materialized"]
            ),
            "reference_stage_eligible_for_global_promotion": False,
        },
        "routing": {
            "retain_reference_stage_as_replayable_source_checkpoint": True,
            "retain_oscillatory_parameters_without_reselection": {
                "amplitude": ROUTED_AMPLITUDE,
                "phase": ROUTED_PHASE,
            },
            "apply_agent3_correction_to_reference_composite": False,
            "reason_agent3_not_applied": (
                "Agent-3 established rank-one covariance realizability on the finite leading core, "
                "but no vector-potential-level public velocity correction or held-out residual "
                "validation exists yet"
            ),
            "next_shortest_blockers": [
                "Agent 1: continue source cone modulation, five-moment repair, heat compensation and matched core-to-heat assembly while addressing the sampled seam divergence",
                "Agent 3: materialize the signed differential amplitude at vector-potential level and rerun disjoint held-in/held-out residual checks",
                "Agent 4: assess the unchanged 1e-3 full-domain gate only after a globally matched composite and compatible fixed/restricted forcing contract exist",
            ],
        },
        "states": {
            "reference_continuation_ready": True,
            "oscillatory_ready": True,
            "signed_covariance_realizability_ready": True,
            "independent_reference_validation_ready": True,
            "reference_composite_velocity_export_ready": True,
            "leading_ready": False,
            "correction_ready": False,
            "velocity_export_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def write_checkpoint(output_dir: str | Path) -> dict[str, Any]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    candidate = KokunoReferenceCompositeCandidate()
    candidate_path = candidate.save_json(output / "reference_composite_candidate.json")
    replay = KokunoReferenceCompositeCandidate.load_json(candidate_path)
    if replay.sha256 != candidate.sha256:
        raise RuntimeError("reference-composite save/load changed candidate identity")

    probes = np.array(
        [[0.11, 0.03, 0.02], [0.31, -0.17, -0.03], [0.64, 0.12, 0.02]],
        dtype=float,
    )
    replay_values = replay.at_points(probes, 0.5)
    direct_values = candidate.at_points(probes, 0.5)
    if not np.allclose(replay_values, direct_values, rtol=0.0, atol=0.0):
        raise RuntimeError("reference-composite replay changed public [u,v,w]")

    checkpoint = build_checkpoint(candidate)
    checkpoint["delivery_smoke"] = {
        "probe_points": probes.tolist(),
        "time": 0.5,
        "velocity": direct_values.tolist(),
        "finite": bool(np.all(np.isfinite(direct_values))),
        "nontrivial": bool(np.max(np.linalg.norm(direct_values, axis=1)) > 0.0),
        "save_load_same_sha": True,
        "used_for_pde_acceptance": False,
    }
    checkpoint_text = _canonical_json(checkpoint)
    checkpoint["checkpoint_sha256"] = hashlib.sha256(checkpoint_text.encode("utf-8")).hexdigest()
    (output / "integration_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return checkpoint


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/kokuno_agent5/reference_integration_checkpoint_v5"),
    )
    args = parser.parse_args(argv)
    checkpoint = write_checkpoint(args.output_dir)
    print(
        json.dumps(
            {
                "task_id": checkpoint["task_id"],
                "candidate_sha256": checkpoint["candidate_sha256"],
                "checkpoint_sha256": checkpoint["checkpoint_sha256"],
                "reference_composite_velocity_export_ready": checkpoint["states"][
                    "reference_composite_velocity_export_ready"
                ],
                "velocity_export_ready": checkpoint["states"]["velocity_export_ready"],
                "pde_validated": checkpoint["states"]["pde_validated"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
