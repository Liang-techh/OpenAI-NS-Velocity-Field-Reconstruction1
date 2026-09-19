"""Prepare the profile-level ``Delta C / epsilon`` mean-inverse right-hand side.

This is the narrow Agent-3 seam after :mod:`kokuno_typed_mean_finite_head_debt`.
The finite-head debt is recomputed from raw same-cycle fields; the caller cannot
supply a residual, defect, stress, debt, target, gain, or success threshold.
The only additional input is a provenance-bearing positive reference-scale
profile ``epsilon`` bound to the same cycle identity and radial grid.

The historical signed-mean correction gate (#493) fixes the physical-to-
reference conversion as

    H_ref y = Delta C / epsilon

profile-wise.  This module materializes only that RHS.  It does *not* claim that
the repository-autonomous finite-head factor equals Kokuno's formal
``missingWeight``; it does not materialize ``H_ref`` or solve the signed inverse;
and it does not authorize a velocity correction or residual claim.
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
from .kokuno_typed_mean_finite_head_debt import (
    TypedMeanFiniteHeadDebtEvaluation,
    materialize_typed_mean_finite_head_debt,
)
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry

TASK = "KOKUNO-A3-TYPED-MEAN-INVERSE-RHS-054"
SCHEMA = "kokuno-a3-typed-mean-inverse-rhs-v1"
PARENT_AGENT3_PR = 686
PARENT_AGENT3_HEAD = "df05669aaee94cd8de2463a35c99a73c476118c1"
HISTORICAL_INVERSE_READINESS_PR = 493
PHYSICAL_TO_REFERENCE_FORMULA = "H_ref y = Delta C / epsilon"


@dataclass(frozen=True)
class MeanInverseReferenceScale:
    """Positive reference-scale profile with explicit provenance and identity.

    A naked epsilon array is intentionally not accepted by the public execution
    API.  This witness is only a typed/provenance boundary; constructing one is
    not, by itself, theorem evidence or proof that the scale is the paper-exact
    Kokuno epsilon for the full candidate.
    """

    identity: CycleIdentity
    radii: tuple[float, ...]
    epsilon: tuple[float, ...]
    producer_kind: str
    provenance: str
    source_reference_scale_certified: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.identity, CycleIdentity):
            raise TypeError("identity must be a CycleIdentity")
        radii = np.asarray(self.radii, dtype=float)
        epsilon = np.asarray(self.epsilon, dtype=float)
        if radii.ndim != 1 or radii.size == 0 or not np.all(np.isfinite(radii)):
            raise ValueError("radii must be a finite nonempty one-dimensional profile")
        if epsilon.shape != radii.shape or not np.all(np.isfinite(epsilon)):
            raise ValueError("epsilon must be a finite profile with the same shape as radii")
        if np.any(epsilon <= 0.0):
            raise ValueError("epsilon must be strictly positive pointwise")
        if np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        if not isinstance(self.producer_kind, str) or not self.producer_kind.strip():
            raise ValueError("producer_kind must be nonempty")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if not isinstance(self.source_reference_scale_certified, bool):
            raise TypeError("source_reference_scale_certified must be bool")
        object.__setattr__(self, "radii", tuple(float(x) for x in radii))
        object.__setattr__(self, "epsilon", tuple(float(x) for x in epsilon))


@dataclass(frozen=True)
class TypedMeanInverseRHSEvaluation:
    """Actual typed debt plus its profile-wise reference-scale normalization."""

    identity: CycleIdentity
    debt: TypedMeanFiniteHeadDebtEvaluation
    reference_scale: MeanInverseReferenceScale
    theta_rhs: np.ndarray
    axial_rhs: np.ndarray
    rhs_rms: float
    rhs_max_abs: float
    reconstruction_closure_max_abs: float

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "reference_scale": {
                "identity": asdict(self.reference_scale.identity),
                "radii": list(self.reference_scale.radii),
                "epsilon": list(self.reference_scale.epsilon),
                "producer_kind": self.reference_scale.producer_kind,
                "provenance": self.reference_scale.provenance,
                "source_reference_scale_certified": self.reference_scale.source_reference_scale_certified,
            },
            "finite_head_debt": self.debt.to_receipt(),
            "mean_inverse_rhs": {
                "formula": PHYSICAL_TO_REFERENCE_FORMULA,
                "ordering": ["theta_e2", "axial_e1"],
                "theta_e2_values": self.theta_rhs.tolist(),
                "axial_e1_values": self.axial_rhs.tolist(),
                "rhs_rms": self.rhs_rms,
                "rhs_max_abs": self.rhs_max_abs,
                "reconstruction_closure_max_abs": self.reconstruction_closure_max_abs,
            },
        }


def materialize_typed_mean_inverse_rhs(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    geometry: RadialMeanStressGeometry,
    reference_scale: MeanInverseReferenceScale,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanInverseRHSEvaluation:
    """Recompute actual typed debt, then form ``Delta C / epsilon`` pointwise."""

    debt = materialize_typed_mean_finite_head_debt(
        velocity,
        velocity_dt,
        pressure,
        restricted_forcing,
        geometry,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    if not isinstance(reference_scale, MeanInverseReferenceScale):
        raise TypeError("reference_scale must be a MeanInverseReferenceScale")
    if reference_scale.identity != debt.identity:
        raise ValueError("reference-scale identity must match the same correction-cycle identity")
    geometry_radii = np.asarray(geometry.radii, dtype=float)
    reference_radii = np.asarray(reference_scale.radii, dtype=float)
    if geometry_radii.shape != reference_radii.shape or not np.array_equal(
        geometry_radii, reference_radii
    ):
        raise ValueError("reference-scale radial grid must match the debt radial grid exactly")

    epsilon = np.asarray(reference_scale.epsilon, dtype=float)
    theta_debt = np.asarray(debt.theta_debt, dtype=float)
    axial_debt = np.asarray(debt.axial_debt, dtype=float)
    if theta_debt.shape != epsilon.shape or axial_debt.shape != epsilon.shape:
        raise RuntimeError("finite-head debt profile shape drifted from the reference scale")

    theta_rhs = theta_debt / epsilon
    axial_rhs = axial_debt / epsilon
    if not np.all(np.isfinite(theta_rhs)) or not np.all(np.isfinite(axial_rhs)):
        raise ArithmeticError("Delta C / epsilon produced non-finite inverse RHS values")
    rhs = np.stack((theta_rhs, axial_rhs), axis=-1)
    reconstructed = np.stack((theta_rhs * epsilon, axial_rhs * epsilon), axis=-1)
    debt_matrix = np.stack((theta_debt, axial_debt), axis=-1)
    closure = float(np.max(np.abs(reconstructed - debt_matrix)))
    rhs_rms = float(np.sqrt(np.mean(np.sum(rhs * rhs, axis=-1))))
    rhs_max = float(np.max(np.abs(rhs)))
    scale = max(1.0, float(np.max(np.abs(debt_matrix))))
    if closure > 256.0 * math.ulp(scale):
        raise ArithmeticError("Delta C / epsilon reconstruction closure exceeded floating-point guard")

    return TypedMeanInverseRHSEvaluation(
        identity=debt.identity,
        debt=debt,
        reference_scale=reference_scale,
        theta_rhs=theta_rhs,
        axial_rhs=axial_rhs,
        rhs_rms=rhs_rms,
        rhs_max_abs=rhs_max,
        reconstruction_closure_max_abs=closure,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_typed_mean_inverse_rhs)
    forbidden = {
        "residual",
        "defect",
        "mean_source",
        "source",
        "stress",
        "requested_stress",
        "debt",
        "delta_c",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "h_ref",
        "inverse_matrix",
        "correction",
    }
    return {
        "physical_to_reference_formula": PHYSICAL_TO_REFERENCE_FORMULA,
        "historical_formula_provenance_pr": HISTORICAL_INVERSE_READINESS_PR,
        "typed_actual_debt_provider": "kokuno_typed_mean_finite_head_debt:materialize_typed_mean_finite_head_debt",
        "raw_same_cycle_providers_required": True,
        "naked_epsilon_public_parameter_allowed": False,
        "provenance_bearing_reference_scale_required": True,
        "caller_supplied_precomputed_defect_allowed": False,
        "caller_supplied_precomputed_debt_allowed": False,
        "caller_supplied_inverse_rhs_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "profile_level_delta_c_over_epsilon_materialized": True,
        "reference_scale_source_certification_required_for_scientific_use": True,
        "formal_missing_weight_materialized": False,
        "formal_missing_weight_equality_claimed": False,
        "h_ref_materialized": False,
        "signed_mean_inverse_solved": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
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
    """Analytic non-f=R mechanics regression for the normalization seam only."""

    identity = CycleIdentity("analytic-typed-mean-inverse-rhs-v1", 0, "t0")
    center = 0.50
    halfwidth = 0.30
    radii = tuple(float(value) for value in np.linspace(0.20, 0.80, 49))

    def profile(x: float, y: float) -> float:
        return _compact_profile(math.hypot(x, y), center=center, halfwidth=halfwidth)

    velocity = VectorFieldProvider(
        identity,
        "analytic:u=t*psi(r)*(-y,x,1)",
        lambda x, y, z, t: (-t * y * profile(x, y), t * x * profile(x, y), t * profile(x, y)),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        "analytic:u_t=psi(r)*(-y,x,1)",
        lambda x, y, z, t: (-y * profile(x, y), x * profile(x, y), profile(x, y)),
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
    reference = MeanInverseReferenceScale(
        identity=identity,
        radii=radii,
        epsilon=tuple(0.25 + 0.05 * radius for radius in radii),
        producer_kind="analytic-regression",
        provenance="fixed positive regression scale independent of computed defect/debt",
        source_reference_scale_certified=False,
    )
    evaluation = materialize_typed_mean_inverse_rhs(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        reference,
        viscosity=0.01,
        spatial_step=0.005,
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "structural_source": "Kokuno corrected reconstruction signed mean-correction architecture",
            "physical_to_reference_formula_provenance_pr": HISTORICAL_INVERSE_READINESS_PR,
            "repository_role": "typed actual finite-head debt -> profile-level Delta C / epsilon inverse RHS",
        },
        "analytic_regression": {
            "field": "u=t*psi(r)*(-y,x,1), u_t=psi(r)*(-y,x,1), p=0, f=0 at t=0",
            "reference_scale_kind": "analytic-regression; not source-certified",
            "evaluation": evaluation.to_receipt(),
        },
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
