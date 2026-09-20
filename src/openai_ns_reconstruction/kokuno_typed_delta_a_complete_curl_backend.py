"""Agent-2 backend for Agent-3 typed signed-amplitude correction handoff.

This module is intentionally narrow. Agent 3 owns the no-surrogate chain that
produces a typed ``delta_a_sigma(R)`` evaluation. Agent 2 owns the localized
vector-potential / complete-curl realization. The backend below translates only
that already-typed amplitude differential into the Agent-2 correction kernel.

It does not import Agent-3 modules, recompute a defect, fit velocity components,
or accept pressure/forcing/residual/gain/target inputs. This keeps the ownership
boundary usable even while the Agent-3 handoff PR is stacked on a different
branch lineage.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from .kokuno_public_signed_amplitude_complete_curl import (
    correction_vector_potential,
    correction_velocity,
    correction_velocity_dt,
    profile_from_delta_a,
)
from .kokuno_public_z_pullback_velocity import default_field

TASK = "KOKUNO-A2-TYPED-DELTA-A-COMPLETE-CURL-BACKEND-050"
SCHEMA = "kokuno-a2-typed-delta-a-complete-curl-backend-v1"
PARENT_AGENT2_PR = 725
PARENT_AGENT2_HEAD = "1930af781dbfff63a0ebac7b58eb102eb473fd73"
AGENT3_HANDOFF_PR = 717
AGENT3_HANDOFF_HEAD = "dd7edcc13fc16e59e4d53d5e4e36646c0c827526"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"


def _finite_1d(values: Any, *, label: str) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.ndim != 1 or out.size == 0 or not np.all(np.isfinite(out)):
        raise ValueError(f"{label} must be a finite nonempty one-dimensional profile")
    return out


@dataclass(frozen=True)
class TypedDeltaACompleteCurlBackend:
    """A2-owned backend matching Agent-3's external evaluator seam.

    ``velocity_evaluator`` and ``velocity_dt_evaluator`` have the exact callable
    shape expected by Agent-3 #717:

        evaluator(typed_amplitude_evaluation, x, y, z, t) -> 3-vector

    The backend deliberately keeps Agent-3 certification hard-false. Independent
    Agent-4 audit is still required before a real full-candidate adapter may be
    labelled source-certified.
    """

    radii: tuple[float, ...]
    reference_time: float
    provenance: str
    producer_kind: str = "agent2-complete-curl-typed-delta-a-backend"

    def __post_init__(self) -> None:
        field = default_field()
        radii = _finite_1d(self.radii, label="radii")
        if radii.size < 8:
            raise ValueError("radii must contain at least 8 nodes")
        if np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        if not (field.radial_inner < radii[0] < radii[-1] < field.radial_outer):
            raise ValueError("typed delta-a radial grid must stay strictly inside A2 support")
        reference_time = float(self.reference_time)
        if not math.isfinite(reference_time) or not (
            field.time_min <= reference_time <= field.time_max
        ):
            raise ValueError("reference_time must lie in the frozen A2 time interval")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if not isinstance(self.producer_kind, str) or not self.producer_kind.strip():
            raise ValueError("producer_kind must be nonempty")
        object.__setattr__(self, "radii", tuple(float(v) for v in radii))
        object.__setattr__(self, "reference_time", reference_time)

    def _profile(self, amplitude_evaluation: Any):
        if not hasattr(amplitude_evaluation, "delta_a"):
            raise TypeError("typed amplitude evaluation must expose delta_a")
        if not hasattr(amplitude_evaluation, "base_amplitudes"):
            raise TypeError("typed amplitude evaluation must expose base_amplitudes")
        base = amplitude_evaluation.base_amplitudes
        if not hasattr(base, "radii"):
            raise TypeError("typed amplitude base witness must expose radii")

        incoming_radii = _finite_1d(base.radii, label="typed amplitude radii")
        expected_radii = np.asarray(self.radii, dtype=float)
        if incoming_radii.shape != expected_radii.shape or not np.array_equal(
            incoming_radii, expected_radii
        ):
            raise ValueError("typed delta-a radial grid must match the A2 backend grid exactly")

        delta_a = np.asarray(amplitude_evaluation.delta_a, dtype=float)
        if delta_a.shape != (expected_radii.size, 2) or not np.all(np.isfinite(delta_a)):
            raise ValueError("typed delta_a must be finite with shape (len(radii),2)")

        return profile_from_delta_a(
            expected_radii,
            delta_a,
            reference_time=self.reference_time,
            producer_kind=self.producer_kind,
            provenance=self.provenance,
            # Hard-false until an independent A4 audit certifies the actual
            # full-candidate reference inputs and this A2 correction backend.
            source_mean_amplitude_differential_certified=False,
        )

    def vector_potential_evaluator(
        self, amplitude_evaluation: Any, x: float, y: float, z: float, t: float
    ) -> np.ndarray:
        return np.asarray(
            correction_vector_potential(self._profile(amplitude_evaluation), x, y, z, t),
            dtype=float,
        )

    def velocity_evaluator(
        self, amplitude_evaluation: Any, x: float, y: float, z: float, t: float
    ) -> np.ndarray:
        return np.asarray(
            correction_velocity(self._profile(amplitude_evaluation), x, y, z, t),
            dtype=float,
        )

    def velocity_dt_evaluator(
        self, amplitude_evaluation: Any, x: float, y: float, z: float, t: float
    ) -> np.ndarray:
        return np.asarray(
            correction_velocity_dt(self._profile(amplitude_evaluation), x, y, z, t),
            dtype=float,
        )

    def agent3_adapter_kwargs(self, identity: Any) -> dict[str, Any]:
        """Return kwargs consumable by Agent-3 ``CompleteCurlCorrectionAdapter``.

        No Agent-3 import is needed on this branch. The certification bit remains
        hard-false; a later independent audit may promote it in a separate PR.
        """

        return {
            "identity": identity,
            "radii": self.radii,
            "producer_kind": self.producer_kind,
            "provenance": self.provenance,
            "velocity_evaluator": self.velocity_evaluator,
            "velocity_dt_evaluator": self.velocity_dt_evaluator,
            "source_agent2_complete_curl_certified": False,
        }


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(TypedDeltaACompleteCurlBackend)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
        "delta_a",
    }
    return {
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "agent3_handoff_pr": AGENT3_HANDOFF_PR,
        "agent3_handoff_head": AGENT3_HANDOFF_HEAD,
        "typed_delta_a_consumed_only_after_agent3_materialization": True,
        "agent3_mean_radial_chain_reimplemented": False,
        "vector_potential_first_complete_curl_reused": True,
        "public_z_pullback_inherited": True,
        "forbidden_constructor_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "source_mean_amplitude_differential_certified_for_full_candidate": False,
        "independent_agent4_correction_audit_required": True,
        "real_full_candidate_bound": False,
        "real_finite_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def _regression_amplitude(radii: np.ndarray, delta_a: np.ndarray) -> Any:
    return SimpleNamespace(
        delta_a=np.asarray(delta_a, dtype=float),
        base_amplitudes=SimpleNamespace(radii=tuple(float(v) for v in radii)),
    )


def verification_receipt() -> dict[str, Any]:
    radii = np.linspace(0.24, 1.26, 17)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 6
    delta_a = np.stack(
        (
            0.010 * bump * (1.0 + 0.12 * np.cos(2.0 * math.pi * s)),
            -0.008 * bump * (1.0 - 0.09 * np.sin(2.0 * math.pi * s)),
        ),
        axis=-1,
    )
    delta_a[[0, -1], :] = 0.0
    amplitude = _regression_amplitude(radii, delta_a)
    provenance = "typed adapter regression; no Agent-3 defect or residual input"
    backend = TypedDeltaACompleteCurlBackend(
        radii=tuple(float(v) for v in radii),
        reference_time=0.50,
        provenance=provenance,
    )
    direct = profile_from_delta_a(
        radii,
        delta_a,
        reference_time=0.50,
        producer_kind=backend.producer_kind,
        provenance=provenance,
        source_mean_amplitude_differential_certified=False,
    )

    rows = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for idx, (radius, z, t) in enumerate(
        zip(
            (0.41, 0.58, 0.76, 0.94) * 3,
            (-0.55, -0.17, 0.23, 0.59) * 3,
            (0.37,) * 4 + (0.50,) * 4 + (0.63,) * 4,
        )
    ):
        theta = 0.213 + (idx + 1) * golden
        rows.append((radius * math.cos(theta), radius * math.sin(theta), z, t))
    points = np.asarray(rows, dtype=float)
    x, y, z, t = points.T

    backend_velocity = np.stack(
        [backend.velocity_evaluator(amplitude, *row) for row in points], axis=0
    )
    backend_dt = np.stack(
        [backend.velocity_dt_evaluator(amplitude, *row) for row in points], axis=0
    )
    backend_A = np.stack(
        [backend.vector_potential_evaluator(amplitude, *row) for row in points], axis=0
    )
    direct_velocity = correction_velocity(direct, x, y, z, t)
    direct_dt = correction_velocity_dt(direct, x, y, z, t)
    direct_A = correction_vector_potential(direct, x, y, z, t)

    payload = backend.agent3_adapter_kwargs(identity="typed-regression-cycle")
    failed = []
    if not np.array_equal(backend_velocity, direct_velocity):
        failed.append("velocity_exact_replay")
    if not np.array_equal(backend_dt, direct_dt):
        failed.append("velocity_dt_exact_replay")
    if not np.array_equal(backend_A, direct_A):
        failed.append("vector_potential_exact_replay")
    if payload["source_agent2_complete_curl_certified"] is not False:
        failed.append("premature_agent2_certification")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "sample_count": int(points.shape[0]),
        "exact_replay": {
            "velocity_max_abs_difference": float(np.max(np.abs(backend_velocity - direct_velocity))),
            "velocity_dt_max_abs_difference": float(np.max(np.abs(backend_dt - direct_dt))),
            "vector_potential_max_abs_difference": float(np.max(np.abs(backend_A - direct_A))),
        },
        "agent3_adapter_payload": {
            "producer_kind": payload["producer_kind"],
            "radii_count": len(payload["radii"]),
            "source_agent2_complete_curl_certified": payload[
                "source_agent2_complete_curl_certified"
            ],
            "velocity_callable": callable(payload["velocity_evaluator"]),
            "velocity_dt_callable": callable(payload["velocity_dt_evaluator"]),
        },
        "failed_guards": failed,
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(verification_receipt(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
