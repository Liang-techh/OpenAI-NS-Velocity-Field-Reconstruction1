"""Agent-5 integration checkpoint for the first callable Kokuno-derived core composite.

This module harvests, without reimplementing their mathematics:
- Agent 1: source-native finite leading-core velocity/pressure;
- Agent 2: exact-curl oscillatory correction and analytic differential operators;
- Agent 3: the first bounded mean-correction-cycle receipt (provenance only here);
- Agent 4: independent fourth-order black-box momentum/divergence evaluator.

The callable artifact produced here is deliberately *core only*: Agent-1 has not
completed the global heat/exterior matching, and Agent-3's selected mean-swirl
correction was measured against the capped-bipolar engineering bridge rather than
against the Agent-1 leading core.  Consequently that correction is recorded but
not transplanted into this composite, and the registered full-domain 1e-3 PDE
gate remains unassessed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_independent_leading_core import evaluate_fd4
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate

SCHEMA = "kokuno-agent5-core-composite-checkpoint-v1"
PDE_THRESHOLD = 1.0e-3
DIVERGENCE_THRESHOLD = 1.0e-5
NU = 0.01
VALIDATION_SEED = 9_172_861
DEFAULT_TIMES = (0.375, 0.5, 0.625)
DEFAULT_STEPS = (0.02, 0.01, 0.005)
DIAGNOSTIC_OSCILLATORY_AMPLITUDE = 0.125

UPSTREAM = {
    "agent1": {
        "pr": 229,
        "head": "6f9cf3c4874087c0218a729df3c84a1fb45fc3a4",
        "role": "source-native finite leading-core velocity/pressure",
    },
    "agent2": {
        "pr": 230,
        "head": "5eec6b7f81831dcb8a8c03f89e99bd84cfb597d4",
        "analytic_operator_head": "930e84e333a829bee0e84d0f60f14e4e7ad26907",
        "role": "complete-curl oscillatory correction plus analytic operators",
    },
    "agent3": {
        "pr": 232,
        "head": "3a7c7598ad8f0b9323b9b9ca025eb1d720106007",
        "role": "one-step stress-shaped mean-correction screen",
        "bridge_candidate_sha256": "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609",
        "selected_amplitude": 0.005,
        "held_out_mean_defect_rms_ratio": 0.96010719,
        "held_out_raw_operator_rms_ratio": 0.99920872,
        "independent_correction_divergence_max": 7.0863e-11,
        "accepted_for_next_cycle": True,
        "stress_to_velocity_lift_source_exact": False,
    },
    "agent4": {
        "pr": 231,
        "head": "46705613bf41486e5d52da8c73a6fb25eb16eeca",
        "role": "independent black-box FD4 leading-core validator",
    },
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _result_metrics(result: dict[str, np.ndarray]) -> dict[str, float]:
    divergence = np.asarray(result["divergence"], dtype=float)
    momentum = np.asarray(result["momentum"], dtype=float)
    norms = np.linalg.norm(momentum, axis=-1)
    return {
        "divergence_max": float(np.max(np.abs(divergence))),
        "divergence_rms": float(np.sqrt(np.mean(divergence * divergence))),
        "momentum_max": float(np.max(norms)),
        "momentum_rms": float(np.sqrt(np.mean(norms * norms))),
    }


@dataclass(frozen=True)
class KokunoCoreCompositeCandidate:
    """Replayable core-only ``leading + complete-curl oscillatory`` candidate."""

    leading: KokunoLeadingCoreSeriesCandidate = field(
        default_factory=KokunoLeadingCoreSeriesCandidate
    )
    oscillatory: KokunoCompleteCurlCorrection = field(
        default_factory=lambda: KokunoCompleteCurlCorrection(
            amplitude=DIAGNOSTIC_OSCILLATORY_AMPLITUDE
        )
    )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        base = np.asarray(self.leading.velocity(x, y, z, t), dtype=float)
        delta = np.asarray(self.oscillatory.velocity(x, y, z, t), dtype=float)
        if base.shape != delta.shape or base.shape[-1] != 3:
            raise ValueError("leading and oscillatory velocities must share (...,3) layout")
        result = base + delta
        if not np.all(np.isfinite(result)):
            raise RuntimeError("core composite velocity became non-finite")
        return result

    __call__ = velocity

    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Use Agent-1 core pressure unchanged; no pressure correction is fitted."""
        return self.leading.pressure(x, y, z, t)

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = np.asarray(points_xyz, dtype=float)
        if points.ndim < 2 or points.shape[-1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points_xyz must be a finite array with final dimension 3")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def grid(self, x: Any, y: Any, z: Any, times: Any) -> np.ndarray:
        x_axis = np.asarray(x, dtype=float)
        y_axis = np.asarray(y, dtype=float)
        z_axis = np.asarray(z, dtype=float)
        t_axis = np.asarray(times, dtype=float)
        for name, axis in (("x", x_axis), ("y", y_axis), ("z", z_axis), ("times", t_axis)):
            if axis.ndim != 1 or axis.size == 0 or not np.all(np.isfinite(axis)):
                raise ValueError(f"{name} must be a nonempty finite one-dimensional axis")
        X, Y, Z = np.meshgrid(x_axis, y_axis, z_axis, indexing="ij")
        out = np.empty((len(t_axis), len(x_axis), len(y_axis), len(z_axis), 3), dtype=float)
        for index, time in enumerate(t_axis):
            out[index] = self.velocity(X, Y, Z, float(time))
        return out

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "upstream": UPSTREAM,
            "components": {
                "leading": {
                    "sha256": self.leading.sha256,
                    "payload": self.leading.to_payload(),
                },
                "oscillatory": {
                    "metadata": self.oscillatory.metadata(),
                    "diagnostic_amplitude_not_promoted": True,
                },
                "agent3_mean_correction": {
                    "upstream_selected_amplitude": UPSTREAM["agent3"]["selected_amplitude"],
                    "measured_against": "capped_bipolar_engineering_bridge",
                    "applied_to_agent1_core": False,
                    "reason": "must be remeasured/reselected on Agent-1 leading core before composition",
                },
            },
            "api": {
                "velocity": "velocity(x,y,z,t)->[...,3] Cartesian [u,v,w]",
                "pressure": "Agent-1 core pressure unchanged",
                "grid_layout": "(time,x,y,z,component)",
            },
            "truth_boundary": {
                "leading_core_ready": True,
                "leading_ready": False,
                "oscillatory_ready": True,
                "correction_one_step_engineering_receipt_ready": True,
                "correction_ready": False,
                "core_velocity_export_ready": True,
                "velocity_export_ready": False,
                "global_outer_matching_complete": False,
                "compatible_final_forcing_attached": False,
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
    def load_json(cls, path: str | Path) -> "KokunoCoreCompositeCandidate":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        expected_sha = payload.pop("sha256", None)
        if payload.get("schema") != SCHEMA or expected_sha is None:
            raise ValueError("invalid Kokuno Agent-5 core composite artifact")
        leading_payload = payload["components"]["leading"]["payload"]
        leading = KokunoLeadingCoreSeriesCandidate.from_payload(leading_payload)
        osc_parameters = payload["components"]["oscillatory"]["metadata"]["parameters"]
        oscillatory = KokunoCompleteCurlCorrection(**osc_parameters)
        candidate = cls(leading=leading, oscillatory=oscillatory)
        if candidate.sha256 != expected_sha:
            raise ValueError("core composite artifact SHA mismatch")
        if _canonical_json(payload) != _canonical_json(candidate.to_payload()):
            raise ValueError("core composite artifact metadata changed")
        return candidate


def _held_out_core_points(seed: int, count: int) -> np.ndarray:
    if count <= 0:
        raise ValueError("count must be positive")
    rng = np.random.default_rng(seed)
    # Chosen so all FD4 +/-2h probes for h<=.02 remain safely inside Lambda*X<=4.1.
    xy = rng.uniform(-0.12, 0.12, size=(count, 2))
    z = rng.uniform(-0.16, 0.16, size=(count, 1))
    return np.concatenate((xy, z), axis=1)


def generate_core_composite_checkpoint(
    *,
    output_dir: str | Path = "artifacts/kokuno_agent5/core_composite_checkpoint_v3",
    seed: int = VALIDATION_SEED,
    point_count: int = 24,
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    leading = KokunoLeadingCoreSeriesCandidate()
    composite = KokunoCoreCompositeCandidate(leading=leading)
    candidate_path = composite.save_json(output_dir / "core_composite_candidate.json")
    replay = KokunoCoreCompositeCandidate.load_json(candidate_path)
    if replay.sha256 != composite.sha256:
        raise RuntimeError("core composite save/load replay changed candidate identity")

    points = _held_out_core_points(seed, point_count)
    rows: list[dict[str, Any]] = []
    for time in DEFAULT_TIMES:
        for step in DEFAULT_STEPS:
            baseline = evaluate_fd4(leading.velocity, leading.pressure, points, time, step, nu=NU)
            augmented = evaluate_fd4(composite.velocity, composite.pressure, points, time, step, nu=NU)
            baseline_metrics = _result_metrics(baseline)
            augmented_metrics = _result_metrics(augmented)
            rows.append({
                "time": float(time),
                "step": float(step),
                "leading_only": baseline_metrics,
                "leading_plus_oscillatory": augmented_metrics,
                "momentum_rms_ratio": augmented_metrics["momentum_rms"] / max(baseline_metrics["momentum_rms"], np.finfo(float).tiny),
            })

    finest = min(DEFAULT_STEPS)
    finest_rows = [row for row in rows if row["step"] == finest]
    baseline_finest = float(np.mean([row["leading_only"]["momentum_rms"] for row in finest_rows]))
    augmented_finest = float(np.mean([row["leading_plus_oscillatory"]["momentum_rms"] for row in finest_rows]))
    max_divergence_finest = float(max(row["leading_plus_oscillatory"]["divergence_max"] for row in finest_rows))

    receipt = {
        "task": "KOKUNO-A5-CORE-COMPOSITE-CHECKPOINT-003",
        "candidate_path": str(candidate_path),
        "candidate_sha256": composite.sha256,
        "upstream": UPSTREAM,
        "sampling": {
            "seed": int(seed),
            "point_count": int(point_count),
            "times": list(DEFAULT_TIMES),
            "steps": list(DEFAULT_STEPS),
            "nu": NU,
            "operator": "Agent-4 independent black-box fourth-order finite differences",
            "pressure": "Agent-1 leading-core pressure, unchanged after oscillatory augmentation",
            "forcing": "zero diagnostic only; no force fitted",
        },
        "baseline_vs_core_composite": {
            "rows": rows,
            "finest_step": finest,
            "mean_leading_only_momentum_rms": baseline_finest,
            "mean_leading_plus_oscillatory_momentum_rms": augmented_finest,
            "mean_ratio": augmented_finest / max(baseline_finest, np.finfo(float).tiny),
            "max_composite_divergence": max_divergence_finest,
        },
        "agent3_bridge_cycle_evidence": {
            **UPSTREAM["agent3"],
            "consumed_as_provenance_only": True,
            "mean_correction_applied_to_core_composite": False,
        },
        "registered_gate": {
            "held_out_normalized_ns_residual_threshold": PDE_THRESHOLD,
            "divergence_reference": DIVERGENCE_THRESHOLD,
            "assessed": False,
            "passed": False,
            "why_unassessed": [
                "Agent-1 candidate remains core-only; global outer/heat matching is missing",
                "Agent-3 selected mean correction was measured against a different engineering bridge and is not transplanted",
                "compatible final fixed/restricted forcing contract is not attached",
            ],
        },
        "stage_state": composite.to_payload()["truth_boundary"],
    }
    receipt_path = output_dir / "integration_checkpoint.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Agent-5 Kokuno core-composite checkpoint")
    parser.add_argument("--output-dir", default="artifacts/kokuno_agent5/core_composite_checkpoint_v3")
    parser.add_argument("--seed", type=int, default=VALIDATION_SEED)
    parser.add_argument("--points", type=int, default=24)
    args = parser.parse_args()
    receipt = generate_core_composite_checkpoint(output_dir=args.output_dir, seed=args.seed, point_count=args.points)
    print(json.dumps({
        "candidate_sha256": receipt["candidate_sha256"],
        "baseline_vs_core_composite": receipt["baseline_vs_core_composite"],
        "registered_gate": receipt["registered_gate"],
        "stage_state": receipt["stage_state"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
