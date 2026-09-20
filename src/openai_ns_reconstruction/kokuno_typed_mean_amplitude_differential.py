"""Convert the typed signed covariance differential into amplitude differentials.

This is the narrow Kokuno Agent-3 seam after
:mod:`kokuno_typed_mean_reference_inverse`. The upstream chain recomputes the
same-cycle momentum defect from raw fields, projects the actual cylindrical
mean, reconstructs compact radial stress, forms the frozen finite-head debt,
normalizes it by a provenance-bearing reference scale, and solves

    H_ref delta_y = Delta C / epsilon.

For the corrected Kokuno signed-pair coordinates ``y_sigma=a_sigma**2``, the
linear differential relation is

    delta_y_sigma = 2 a_sigma delta_a_sigma,

hence

    delta_a_sigma = delta_y_sigma / (2 a_sigma).

This module implements only that differential coordinate conversion. It does
not rebuild Agent-2 complete curls, does not materialize a velocity correction,
and does not claim that the regression amplitudes are source-certified values
for the real full candidate.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from .kokuno_typed_mean_inverse_rhs import MeanInverseReferenceScale
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry
from .kokuno_typed_mean_reference_inverse import (
    MeanReferenceInverseOperator,
    TypedMeanReferenceInverseEvaluation,
    materialize_typed_mean_reference_inverse,
)

TASK = "KOKUNO-A3-TYPED-MEAN-AMPLITUDE-DIFFERENTIAL-056"
SCHEMA = "kokuno-a3-typed-mean-amplitude-differential-v1"
PARENT_AGENT3_PR = 703
PARENT_AGENT3_HEAD = "56957cd51c7e6638608f28f8fa21df590194dc84"
SOURCE_SIGNED_PAIR_PR = 472
SOURCE_DIFFERENTIAL_FORMULA_PR = 253
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
AMPLITUDE_DIFFERENTIAL_FORMULA = "delta_a_sigma = delta_y_sigma / (2 * a_sigma)"


def _positive_profile(
    values: tuple[float, ...], radii: np.ndarray, name: str
) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.shape != radii.shape or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite profile on the exact radial grid")
    if np.any(out <= 0.0):
        raise ValueError(f"{name} must be strictly positive pointwise")
    return out


@dataclass(frozen=True)
class MeanReferenceBaseAmplitudes:
    """Provenance-bearing positive base amplitudes for the signed pair."""

    identity: CycleIdentity
    radii: tuple[float, ...]
    a_plus: tuple[float, ...]
    a_minus: tuple[float, ...]
    producer_kind: str
    provenance: str
    source_reference_base_amplitudes_certified: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.identity, CycleIdentity):
            raise TypeError("identity must be a CycleIdentity")
        radii = np.asarray(self.radii, dtype=float)
        if radii.ndim != 1 or radii.size == 0 or not np.all(np.isfinite(radii)):
            raise ValueError("radii must be a finite nonempty one-dimensional profile")
        if np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        a_plus = _positive_profile(self.a_plus, radii, "a_plus")
        a_minus = _positive_profile(self.a_minus, radii, "a_minus")
        if not isinstance(self.producer_kind, str) or not self.producer_kind.strip():
            raise ValueError("producer_kind must be nonempty")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if not isinstance(self.source_reference_base_amplitudes_certified, bool):
            raise TypeError("source_reference_base_amplitudes_certified must be bool")
        object.__setattr__(self, "radii", tuple(float(x) for x in radii))
        object.__setattr__(self, "a_plus", tuple(float(x) for x in a_plus))
        object.__setattr__(self, "a_minus", tuple(float(x) for x in a_minus))

    def stacked(self) -> np.ndarray:
        return np.stack(
            (np.asarray(self.a_plus, dtype=float), np.asarray(self.a_minus, dtype=float)),
            axis=-1,
        )


@dataclass(frozen=True)
class TypedMeanAmplitudeDifferentialEvaluation:
    """Typed reference inverse plus ``delta_y -> delta_a`` conversion."""

    identity: CycleIdentity
    reference_inverse: TypedMeanReferenceInverseEvaluation
    base_amplitudes: MeanReferenceBaseAmplitudes
    delta_y: np.ndarray
    delta_a: np.ndarray
    multiplication_back_closure_max_abs: float
    maximum_relative_differential: float
    minimum_full_step_amplitude: float
    full_step_positive_everywhere: bool

    def to_receipt(self) -> dict[str, object]:
        base = self.base_amplitudes
        return {
            "identity": asdict(self.identity),
            "reference_inverse": self.reference_inverse.to_receipt(),
            "base_amplitudes": {
                "identity": asdict(base.identity),
                "radii": list(base.radii),
                "a_plus": list(base.a_plus),
                "a_minus": list(base.a_minus),
                "producer_kind": base.producer_kind,
                "provenance": base.provenance,
                "source_reference_base_amplitudes_certified": (
                    base.source_reference_base_amplitudes_certified
                ),
            },
            "amplitude_differential": {
                "formula": AMPLITUDE_DIFFERENTIAL_FORMULA,
                "ordering": ["sigma_plus", "sigma_minus"],
                "delta_y": self.delta_y.tolist(),
                "delta_a": self.delta_a.tolist(),
                "delta_a_vector_rms": float(
                    np.sqrt(np.mean(np.sum(self.delta_a * self.delta_a, axis=-1)))
                ),
                "multiplication_back_closure_max_abs": (
                    self.multiplication_back_closure_max_abs
                ),
                "maximum_relative_differential": self.maximum_relative_differential,
                "minimum_full_step_amplitude": self.minimum_full_step_amplitude,
                "full_step_positive_everywhere": self.full_step_positive_everywhere,
                "full_step_positivity_is_diagnostic_not_admission_gate": True,
            },
        }


def materialize_typed_mean_amplitude_differential(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    geometry: RadialMeanStressGeometry,
    reference_scale: MeanInverseReferenceScale,
    reference_operator: MeanReferenceInverseOperator,
    base_amplitudes: MeanReferenceBaseAmplitudes,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanAmplitudeDifferentialEvaluation:
    """Recompute the no-surrogate chain and materialize signed ``delta_a``."""

    inverse = materialize_typed_mean_reference_inverse(
        velocity,
        velocity_dt,
        pressure,
        restricted_forcing,
        geometry,
        reference_scale,
        reference_operator,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    if not isinstance(base_amplitudes, MeanReferenceBaseAmplitudes):
        raise TypeError("base_amplitudes must be a MeanReferenceBaseAmplitudes witness")
    if base_amplitudes.identity != inverse.identity:
        raise ValueError("base-amplitude identity must match the correction-cycle identity")

    geometry_radii = np.asarray(geometry.radii, dtype=float)
    base_radii = np.asarray(base_amplitudes.radii, dtype=float)
    if geometry_radii.shape != base_radii.shape or not np.array_equal(
        geometry_radii, base_radii
    ):
        raise ValueError("base-amplitude radial grid must match the debt radial grid exactly")

    base = base_amplitudes.stacked()
    delta_y = np.asarray(inverse.delta_y, dtype=float)
    if base.shape != delta_y.shape:
        raise RuntimeError("base-amplitude shape drifted from the signed reference inverse")

    delta_a = delta_y / (2.0 * base)
    if not np.all(np.isfinite(delta_a)):
        raise ArithmeticError("signed amplitude differential produced non-finite coordinates")

    reconstructed = 2.0 * base * delta_a
    closure = float(np.max(np.abs(reconstructed - delta_y)))
    scale = max(1.0, float(np.max(np.abs(delta_y))))
    if closure > 1024.0 * math.ulp(scale):
        raise ArithmeticError("signed amplitude differential closure exceeded guard")

    relative = np.abs(delta_a) / base
    maximum_relative = float(np.max(relative))
    full_step = base + delta_a
    minimum_full_step = float(np.min(full_step))
    full_step_positive = bool(np.all(full_step > 0.0))

    return TypedMeanAmplitudeDifferentialEvaluation(
        identity=inverse.identity,
        reference_inverse=inverse,
        base_amplitudes=base_amplitudes,
        delta_y=delta_y,
        delta_a=delta_a,
        multiplication_back_closure_max_abs=closure,
        maximum_relative_differential=maximum_relative,
        minimum_full_step_amplitude=minimum_full_step,
        full_step_positive_everywhere=full_step_positive,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_typed_mean_amplitude_differential)
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
        "correction",
    }
    return {
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "source_formula_provenance_prs": [
            SOURCE_SIGNED_PAIR_PR,
            SOURCE_DIFFERENTIAL_FORMULA_PR,
        ],
        "amplitude_differential_formula": AMPLITUDE_DIFFERENTIAL_FORMULA,
        "typed_reference_inverse_provider": (
            "kokuno_typed_mean_reference_inverse:materialize_typed_mean_reference_inverse"
        ),
        "raw_same_cycle_providers_required": True,
        "caller_supplied_delta_y_allowed": False,
        "caller_supplied_delta_a_allowed": False,
        "provenance_bearing_positive_base_amplitudes_required": True,
        "source_base_amplitude_certification_required_for_scientific_use": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "typed_amplitude_differential_executable": True,
        "source_epsilon_materialized_for_full_candidate": False,
        "source_h_ref_materialized_for_full_candidate": False,
        "source_base_amplitudes_materialized_for_full_candidate": False,
        "real_candidate_delta_amplitudes_materialized": False,
        "public_velocity_correction_materialized": False,
        "agent2_complete_curl_reimplemented": False,
        "signed_mean_inverse_input_ready_for_real_candidate": False,
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


def deterministic_receipt() -> dict[str, object]:
    """Non-``f=R`` mechanics regression; base amplitudes are not source-certified."""

    identity = CycleIdentity("analytic-typed-mean-amplitude-differential-v1", 0, "t0")
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
    evaluation = materialize_typed_mean_amplitude_differential(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        reference_scale,
        reference_operator,
        base_amplitudes,
        viscosity=0.01,
        spatial_step=0.005,
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "source_signed_pair_pr": SOURCE_SIGNED_PAIR_PR,
            "source_differential_formula_pr": SOURCE_DIFFERENTIAL_FORMULA_PR,
            "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "analytic_regression": {
            "forcing_is_zero_and_fixed_independently_of_residual": True,
            "source_reference_scale_certified": False,
            "source_reference_operator_certified": False,
            "source_reference_base_amplitudes_certified": False,
            "evaluation": evaluation.to_receipt(),
        },
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    payload = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
