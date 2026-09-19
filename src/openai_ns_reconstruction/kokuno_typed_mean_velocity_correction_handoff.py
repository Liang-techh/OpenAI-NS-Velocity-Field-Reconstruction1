"""Typed handoff from Agent-3 signed amplitude differentials to a correction field.

This is intentionally an adapter seam, not another curl implementation.  Agent 3
owns the mean/radial correction chain through the typed signed amplitude
differential ``delta_a``.  Agent 2 owns the public oscillatory complete-curl
family.  This module keeps that ownership boundary explicit:

    raw same-cycle fields
      -> actual NS defect
      -> cylindrical mean channels
      -> compact radial stress
      -> finite-head debt
      -> Delta C / epsilon
      -> H_ref^{-1}
      -> delta_y
      -> delta_a
      -> EXTERNAL complete-curl adapter
      -> CorrectionFieldProvider(delta_u, delta_u_t).

The external adapter receives the *typed* amplitude-differential evaluation, not
caller-supplied residual/defect/stress/target arrays.  It must provide both the
velocity correction and its time derivative because the finite correction-cycle
contract audits both identities before recomputing the true NS defect.

No Agent-2 curl formula is copied here.  The deterministic receipt uses a plainly
labelled analytic regression adapter only to test the typed plumbing; it is not a
source-certified Kokuno correction and is not evidence of residual reduction.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from .kokuno_finite_correction_cycle_contract import CorrectionFieldProvider
from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from .kokuno_typed_mean_amplitude_differential import (
    MeanReferenceBaseAmplitudes,
    TypedMeanAmplitudeDifferentialEvaluation,
    materialize_typed_mean_amplitude_differential,
)
from .kokuno_typed_mean_inverse_rhs import MeanInverseReferenceScale
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry
from .kokuno_typed_mean_reference_inverse import MeanReferenceInverseOperator

TASK = "KOKUNO-A3-TYPED-VELOCITY-CORRECTION-HANDOFF-057"
SCHEMA = "kokuno-a3-typed-velocity-correction-handoff-v1"
PARENT_AGENT3_PR = 709
PARENT_AGENT3_HEAD = "bfaebe3339c984b78e1a797c634b588d61a15cad"
AGENT2_PUBLIC_COMPLETE_CURL_PR = 616
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

AmplitudeAwareVectorEvaluator = Callable[
    [TypedMeanAmplitudeDifferentialEvaluation, float, float, float, float],
    Sequence[float],
]


@dataclass(frozen=True)
class CompleteCurlCorrectionAdapter:
    """External complete-curl provider owned outside Agent 3.

    ``source_agent2_complete_curl_certified`` may be true only for an adapter that
    actually wraps the admitted Agent-2 public complete-curl family (including the
    matching time-derivative semantics).  The flag is deliberately false in the
    deterministic regression in this module.
    """

    identity: CycleIdentity
    radii: tuple[float, ...]
    producer_kind: str
    provenance: str
    velocity_evaluator: AmplitudeAwareVectorEvaluator
    velocity_dt_evaluator: AmplitudeAwareVectorEvaluator
    source_agent2_complete_curl_certified: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.identity, CycleIdentity):
            raise TypeError("identity must be a CycleIdentity")
        radii = np.asarray(self.radii, dtype=float)
        if radii.ndim != 1 or radii.size == 0 or not np.all(np.isfinite(radii)):
            raise ValueError("radii must be a finite nonempty one-dimensional profile")
        if np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        if not isinstance(self.producer_kind, str) or not self.producer_kind.strip():
            raise ValueError("producer_kind must be nonempty")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if not callable(self.velocity_evaluator) or not callable(self.velocity_dt_evaluator):
            raise TypeError("complete-curl velocity and time-derivative evaluators must be callable")
        if not isinstance(self.source_agent2_complete_curl_certified, bool):
            raise TypeError("source_agent2_complete_curl_certified must be bool")
        if self.source_agent2_complete_curl_certified:
            normalized = self.producer_kind.lower().replace("_", "-")
            if "agent2" not in normalized or "complete-curl" not in normalized:
                raise ValueError(
                    "source-certified complete-curl adapters must identify Agent 2 complete-curl provenance"
                )
        object.__setattr__(self, "radii", tuple(float(x) for x in radii))


@dataclass(frozen=True)
class TypedMeanVelocityCorrectionHandoff:
    """Typed amplitude differential plus a finite-cycle correction provider."""

    identity: CycleIdentity
    to_identity: CycleIdentity
    amplitude_differential: TypedMeanAmplitudeDifferentialEvaluation
    adapter: CompleteCurlCorrectionAdapter
    correction_provider: CorrectionFieldProvider

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "to_identity": asdict(self.to_identity),
            "amplitude_differential": self.amplitude_differential.to_receipt(),
            "adapter": {
                "identity": asdict(self.adapter.identity),
                "radii": list(self.adapter.radii),
                "producer_kind": self.adapter.producer_kind,
                "provenance": self.adapter.provenance,
                "source_agent2_complete_curl_certified": (
                    self.adapter.source_agent2_complete_curl_certified
                ),
            },
            "correction_provider": {
                "from_identity": asdict(self.correction_provider.from_identity),
                "to_identity": asdict(self.correction_provider.to_identity),
                "source_ref": self.correction_provider.source_ref,
            },
        }


def _finite_vector(value: Sequence[float], *, label: str) -> tuple[float, float, float]:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must return a finite 3-vector")
    return tuple(float(x) for x in arr)


def materialize_typed_mean_velocity_correction_handoff(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    geometry: RadialMeanStressGeometry,
    reference_scale: MeanInverseReferenceScale,
    reference_operator: MeanReferenceInverseOperator,
    base_amplitudes: MeanReferenceBaseAmplitudes,
    complete_curl_adapter: CompleteCurlCorrectionAdapter,
    to_identity: CycleIdentity,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanVelocityCorrectionHandoff:
    """Materialize a finite-cycle correction provider without implementing curl.

    The no-surrogate Agent-3 chain is recomputed internally through ``delta_a``.
    Only then is the typed amplitude differential passed to the external
    complete-curl adapter.  The resulting ``CorrectionFieldProvider`` is directly
    consumable by the existing finite correction-cycle audit contract.
    """

    amplitude = materialize_typed_mean_amplitude_differential(
        velocity,
        velocity_dt,
        pressure,
        restricted_forcing,
        geometry,
        reference_scale,
        reference_operator,
        base_amplitudes,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    if not isinstance(complete_curl_adapter, CompleteCurlCorrectionAdapter):
        raise TypeError("complete_curl_adapter must be a CompleteCurlCorrectionAdapter")
    if complete_curl_adapter.identity != amplitude.identity:
        raise ValueError("complete-curl adapter identity must match the correction-cycle identity")

    geometry_radii = np.asarray(geometry.radii, dtype=float)
    adapter_radii = np.asarray(complete_curl_adapter.radii, dtype=float)
    if geometry_radii.shape != adapter_radii.shape or not np.array_equal(
        geometry_radii, adapter_radii
    ):
        raise ValueError("complete-curl adapter radial grid must match the delta-a radial grid exactly")

    if not isinstance(to_identity, CycleIdentity):
        raise TypeError("to_identity must be a CycleIdentity")
    if to_identity.cycle_id != amplitude.identity.cycle_id:
        raise ValueError("velocity correction must remain within one cycle_id")
    if to_identity.cycle_index != amplitude.identity.cycle_index + 1:
        raise ValueError("velocity correction must advance exactly one cycle index")
    if to_identity.state_token == amplitude.identity.state_token:
        raise ValueError("velocity correction must change the cycle state token")

    def correction_velocity(x: float, y: float, z: float, t: float) -> tuple[float, float, float]:
        return _finite_vector(
            complete_curl_adapter.velocity_evaluator(amplitude, x, y, z, t),
            label="complete-curl correction velocity",
        )

    def correction_velocity_dt(
        x: float, y: float, z: float, t: float
    ) -> tuple[float, float, float]:
        return _finite_vector(
            complete_curl_adapter.velocity_dt_evaluator(amplitude, x, y, z, t),
            label="complete-curl correction velocity_dt",
        )

    provider = CorrectionFieldProvider(
        from_identity=amplitude.identity,
        to_identity=to_identity,
        source_ref=(
            f"typed-mean-delta-a->{complete_curl_adapter.producer_kind}:"
            f"{complete_curl_adapter.provenance}"
        ),
        velocity_evaluator=correction_velocity,
        velocity_dt_evaluator=correction_velocity_dt,
    )
    return TypedMeanVelocityCorrectionHandoff(
        identity=amplitude.identity,
        to_identity=to_identity,
        amplitude_differential=amplitude,
        adapter=complete_curl_adapter,
        correction_provider=provider,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_typed_mean_velocity_correction_handoff)
    forbidden = {
        "residual",
        "defect",
        "mean_source",
        "source",
        "stress",
        "requested_stress",
        "debt",
        "delta_c",
        "inverse_rhs",
        "rhs",
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
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "agent2_public_complete_curl_pr": AGENT2_PUBLIC_COMPLETE_CURL_PR,
        "raw_same_cycle_providers_required": True,
        "typed_delta_a_recomputed_internally": True,
        "external_complete_curl_adapter_required": True,
        "velocity_and_time_derivative_required_from_adapter": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_delta_y_allowed": False,
        "caller_supplied_delta_a_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "agent3_reimplements_agent2_complete_curl": False,
        "typed_velocity_correction_handoff_executable": True,
        "source_epsilon_materialized_for_full_candidate": False,
        "source_h_ref_materialized_for_full_candidate": False,
        "source_base_amplitudes_materialized_for_full_candidate": False,
        "source_agent2_complete_curl_adapter_bound_for_full_candidate": False,
        "real_candidate_delta_amplitudes_materialized": False,
        "public_velocity_correction_materialized_for_real_candidate": False,
        "real_full_candidate_defect_consumed": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
    }


def _compact_profile(radius: float, *, center: float, halfwidth: float) -> float:
    s = (radius - center) / halfwidth
    if abs(s) >= 1.0:
        return 0.0
    return math.cos(0.5 * math.pi * s) ** 8


def _analytic_regression_handoff() -> TypedMeanVelocityCorrectionHandoff:
    identity = CycleIdentity("analytic-typed-velocity-correction-handoff-v1", 0, "t0")
    to_identity = CycleIdentity(identity.cycle_id, 1, "t1")
    center = 0.50
    halfwidth = 0.30
    radii = tuple(float(value) for value in np.linspace(0.20, 0.80, 49))

    def profile(x: float, y: float) -> float:
        return _compact_profile(math.hypot(x, y), center=center, halfwidth=halfwidth)

    velocity = VectorFieldProvider(
        identity,
        "analytic:u=t*psi(r)*(-y,x,1)",
        lambda x, y, z, t: (
            -t * y * profile(x, y),
            t * x * profile(x, y),
            t * profile(x, y),
        ),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        "analytic:u_t=psi(r)*(-y,x,1)",
        lambda x, y, z, t: (
            -y * profile(x, y),
            x * profile(x, y),
            profile(x, y),
        ),
    )
    pressure = ScalarFieldProvider(identity, "analytic:p=0", lambda x, y, z, t: 0.0)
    forcing = RestrictedForcingProvider(
        identity,
        "analytic:f=0",
        "regression-zero-forcing; fixed independently of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    geometry = RadialMeanStressGeometry(
        time=0.0,
        axial_z=0.13,
        radii=radii,
        bump_center=center,
        bump_halfwidth=halfwidth,
        angular_count=32,
        phase=0.017,
    )
    reference_scale = MeanInverseReferenceScale(
        identity=identity,
        radii=radii,
        epsilon=tuple(0.25 + 0.05 * radius for radius in radii),
        producer_kind="analytic-regression",
        provenance="fixed positive regression scale independent of computed defect/debt",
        source_reference_scale_certified=False,
    )
    reference_operator = MeanReferenceInverseOperator(
        identity=identity,
        radii=radii,
        A_c=tuple(2.0 for _ in radii),
        u_star=tuple(3.0 for _ in radii),
        h_plus=tuple(1.0 for _ in radii),
        h_minus=tuple(1.0 for _ in radii),
        producer_kind="analytic-regression",
        provenance="displayed H_ref regression with A_c=2,u_*=3,h_+=h_-=1",
        source_reference_operator_certified=False,
    )
    base_amplitudes = MeanReferenceBaseAmplitudes(
        identity=identity,
        radii=radii,
        a_plus=tuple(2.0 for _ in radii),
        a_minus=tuple(3.0 for _ in radii),
        producer_kind="analytic-regression",
        provenance="fixed positive base amplitudes independent of computed defect/debt",
        source_reference_base_amplitudes_certified=False,
    )

    def regression_velocity(
        amplitude: TypedMeanAmplitudeDifferentialEvaluation,
        x: float,
        y: float,
        z: float,
        t: float,
    ) -> tuple[float, float, float]:
        # Analytic divergence-free plumbing probe only; NOT the Agent-2 curl.
        scale = float(np.sqrt(np.mean(amplitude.delta_a * amplitude.delta_a)))
        return (-scale * y, scale * x, 0.0)

    def regression_velocity_dt(
        amplitude: TypedMeanAmplitudeDifferentialEvaluation,
        x: float,
        y: float,
        z: float,
        t: float,
    ) -> tuple[float, float, float]:
        return (0.0, 0.0, 0.0)

    adapter = CompleteCurlCorrectionAdapter(
        identity=identity,
        radii=radii,
        producer_kind="analytic-regression-not-agent2-complete-curl",
        provenance="divergence-free typed-plumbing probe; no Agent-2 curl formula reproduced",
        velocity_evaluator=regression_velocity,
        velocity_dt_evaluator=regression_velocity_dt,
        source_agent2_complete_curl_certified=False,
    )
    return materialize_typed_mean_velocity_correction_handoff(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        reference_scale,
        reference_operator,
        base_amplitudes,
        adapter,
        to_identity,
        viscosity=0.01,
        spatial_step=0.005,
    )


def deterministic_receipt() -> dict[str, object]:
    handoff = _analytic_regression_handoff()
    provider = handoff.correction_provider
    points = np.asarray(
        [
            (0.20, 0.11, -0.07, 0.0),
            (-0.31, 0.19, 0.08, 0.0),
            (0.14, -0.27, 0.03, 0.0),
        ],
        dtype=float,
    )
    values = np.asarray([provider.velocity_evaluator(*p) for p in points], dtype=float)
    dt_values = np.asarray([provider.velocity_dt_evaluator(*p) for p in points], dtype=float)
    norms = np.linalg.norm(values, axis=1)
    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "agent2_public_complete_curl_pr": AGENT2_PUBLIC_COMPLETE_CURL_PR,
            "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "analytic_regression": {
            "forcing_is_zero_and_fixed_independently_of_residual": True,
            "adapter_is_source_agent2_certified": False,
            "adapter_is_only_typed_plumbing_probe": True,
            "correction_vector_rms": float(np.sqrt(np.mean(norms * norms))),
            "correction_vector_max": float(np.max(norms)),
            "correction_dt_vector_max": float(np.max(np.linalg.norm(dt_values, axis=1))),
            "correction_nontrivial": bool(np.max(norms) > 0.0),
            "from_identity": asdict(provider.from_identity),
            "to_identity": asdict(provider.to_identity),
        },
        "handoff": handoff.to_receipt(),
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
