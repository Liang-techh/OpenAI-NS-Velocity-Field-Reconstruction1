"""Solve the typed mean-correction reference 2x2 inverse without surrogate RHS input.

This is the narrow Agent-3 seam after :mod:`kokuno_typed_mean_inverse_rhs`.
The upstream module recomputes the actual same-cycle NS defect, angular mean
channels, compact radial stress, finite-head debt, and profile-level
``Delta C / epsilon`` right-hand side.  This module adds only the displayed
two-sign reference covariance inverse from the corrected Kokuno reconstruction:

    H_ref =
      [[-A_c h_+, -A_c h_-],
       [-u_* h_+, +u_* h_-]]

and solves, profile-wise,

    H_ref delta_y = Delta C / epsilon.

The solve is intentionally expressed in signed differential coordinates
``delta_y``.  It does not identify them with a new positive base amplitude,
does not divide by ``2 a_sigma`` to form ``delta a_sigma``, and does not
materialize a correction velocity.  A source-certified operator/reference
scale and the real full candidate remain separate prerequisites.
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
from .kokuno_typed_mean_inverse_rhs import (
    MeanInverseReferenceScale,
    TypedMeanInverseRHSEvaluation,
    materialize_typed_mean_inverse_rhs,
)
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry

TASK = "KOKUNO-A3-TYPED-MEAN-REFERENCE-INVERSE-055"
SCHEMA = "kokuno-a3-typed-mean-reference-inverse-v1"
PARENT_AGENT3_PR = 696
PARENT_AGENT3_HEAD = "a8b87637710df0469333e1c3fec03f3713fe4e0e"
SOURCE_SIGNED_PAIR_PR = 472
SOURCE_UNIT_BRIDGE_PR = 474
UPSTREAM_AGENT3_REFERENCE_INVERSE_PR = 703
UPSTREAM_AGENT3_REFERENCE_INVERSE_HEAD = "56957cd51c7e6638608f28f8fa21df590194dc84"
A5_NUMPY2_COMPATIBILITY_NOTE = (
    "Agent-5 integration preserves the Agent-3 matrix/RHS algebra and only makes "
    "the stacked vector RHS axis explicit for NumPy >=2 batch solve semantics."
)
REFERENCE_MATRIX_FORMULA = (
    "H_ref=[[-A_c*h_plus,-A_c*h_minus],"
    "[-u_star*h_plus,+u_star*h_minus]]"
)
SIGNED_REFERENCE_SOLVE_FORMULA = "H_ref * delta_y = Delta C / epsilon"


def _profile(values: tuple[float, ...], radii: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.shape != radii.shape or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite profile on the exact radial grid")
    return out


@dataclass(frozen=True)
class MeanReferenceInverseOperator:
    """Provenance-bearing profile of the displayed two-sign reference map.

    The operator witness is a typed boundary, not a theorem certificate by
    itself.  Scientific use requires ``source_reference_operator_certified`` to
    come from an independently admitted source binding.
    """

    identity: CycleIdentity
    radii: tuple[float, ...]
    A_c: tuple[float, ...]
    u_star: tuple[float, ...]
    h_plus: tuple[float, ...]
    h_minus: tuple[float, ...]
    producer_kind: str
    provenance: str
    source_reference_operator_certified: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.identity, CycleIdentity):
            raise TypeError("identity must be a CycleIdentity")
        radii = np.asarray(self.radii, dtype=float)
        if radii.ndim != 1 or radii.size == 0 or not np.all(np.isfinite(radii)):
            raise ValueError("radii must be a finite nonempty one-dimensional profile")
        if np.any(np.diff(radii) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        profiles: dict[str, np.ndarray] = {}
        for name in ("A_c", "u_star", "h_plus", "h_minus"):
            profiles[name] = _profile(tuple(getattr(self, name)), radii, name)
            if np.any(profiles[name] <= 0.0):
                raise ValueError(f"{name} must be strictly positive pointwise")
        if not isinstance(self.producer_kind, str) or not self.producer_kind.strip():
            raise ValueError("producer_kind must be nonempty")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if not isinstance(self.source_reference_operator_certified, bool):
            raise TypeError("source_reference_operator_certified must be bool")
        object.__setattr__(self, "radii", tuple(float(x) for x in radii))
        for name, values in profiles.items():
            object.__setattr__(self, name, tuple(float(x) for x in values))

    def matrices(self) -> np.ndarray:
        """Return the profile of source-displayed 2x2 reference matrices."""
        A_c = np.asarray(self.A_c, dtype=float)
        u_star = np.asarray(self.u_star, dtype=float)
        h_plus = np.asarray(self.h_plus, dtype=float)
        h_minus = np.asarray(self.h_minus, dtype=float)
        matrices = np.empty((A_c.size, 2, 2), dtype=float)
        matrices[:, 0, 0] = -A_c * h_plus
        matrices[:, 0, 1] = -A_c * h_minus
        matrices[:, 1, 0] = -u_star * h_plus
        matrices[:, 1, 1] = u_star * h_minus
        return matrices


@dataclass(frozen=True)
class TypedMeanReferenceInverseEvaluation:
    """Typed actual RHS plus the signed differential reference-coordinate solve."""

    identity: CycleIdentity
    rhs: TypedMeanInverseRHSEvaluation
    operator: MeanReferenceInverseOperator
    theta_axial_rhs: np.ndarray
    delta_y: np.ndarray
    determinants: np.ndarray
    condition_numbers: np.ndarray
    reconstruction_closure_max_abs: float

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "inverse_rhs": self.rhs.to_receipt(),
            "reference_operator": {
                "identity": asdict(self.operator.identity),
                "radii": list(self.operator.radii),
                "A_c": list(self.operator.A_c),
                "u_star": list(self.operator.u_star),
                "h_plus": list(self.operator.h_plus),
                "h_minus": list(self.operator.h_minus),
                "producer_kind": self.operator.producer_kind,
                "provenance": self.operator.provenance,
                "source_reference_operator_certified": (
                    self.operator.source_reference_operator_certified
                ),
                "matrix_formula": REFERENCE_MATRIX_FORMULA,
            },
            "signed_reference_inverse": {
                "formula": SIGNED_REFERENCE_SOLVE_FORMULA,
                "channel_ordering": ["theta_e2", "axial_e1"],
                "rhs": self.theta_axial_rhs.tolist(),
                "delta_y_ordering": ["sigma_plus", "sigma_minus"],
                "delta_y": self.delta_y.tolist(),
                "determinants": self.determinants.tolist(),
                "condition_numbers": self.condition_numbers.tolist(),
                "minimum_abs_determinant": float(np.min(np.abs(self.determinants))),
                "maximum_condition_number": float(np.max(self.condition_numbers)),
                "delta_y_vector_rms": float(
                    np.sqrt(np.mean(np.sum(self.delta_y * self.delta_y, axis=-1)))
                ),
                "reconstruction_closure_max_abs": self.reconstruction_closure_max_abs,
            },
        }


def materialize_typed_mean_reference_inverse(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    geometry: RadialMeanStressGeometry,
    reference_scale: MeanInverseReferenceScale,
    reference_operator: MeanReferenceInverseOperator,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanReferenceInverseEvaluation:
    """Recompute the typed RHS and solve the displayed 2x2 reference map."""

    rhs = materialize_typed_mean_inverse_rhs(
        velocity,
        velocity_dt,
        pressure,
        restricted_forcing,
        geometry,
        reference_scale,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    if not isinstance(reference_operator, MeanReferenceInverseOperator):
        raise TypeError("reference_operator must be a MeanReferenceInverseOperator")
    if reference_operator.identity != rhs.identity:
        raise ValueError("reference-operator identity must match the correction-cycle identity")

    geometry_radii = np.asarray(geometry.radii, dtype=float)
    operator_radii = np.asarray(reference_operator.radii, dtype=float)
    if geometry_radii.shape != operator_radii.shape or not np.array_equal(
        geometry_radii, operator_radii
    ):
        raise ValueError("reference-operator radial grid must match the debt radial grid exactly")

    matrices = reference_operator.matrices()
    target = np.stack(
        (np.asarray(rhs.theta_rhs, dtype=float), np.asarray(rhs.axial_rhs, dtype=float)),
        axis=-1,
    )
    if matrices.shape != (target.shape[0], 2, 2):
        raise RuntimeError("reference-operator shape drifted from the typed inverse RHS")

    determinants = np.linalg.det(matrices)
    if not np.all(np.isfinite(determinants)) or np.any(determinants == 0.0):
        raise ArithmeticError("reference operator is singular or non-finite")
    singular_values = np.linalg.svd(matrices, compute_uv=False)
    if not np.all(np.isfinite(singular_values)) or np.any(singular_values[:, 1] <= 0.0):
        raise ArithmeticError("reference operator singular values are invalid")
    condition_numbers = singular_values[:, 0] / singular_values[:, 1]

    # NumPy >=2 treats a two-dimensional RHS as a stack of (M, K) matrices.
    # Make the one-vector-per-radius solve explicit, then remove only that axis.
    # This is a representation/API compatibility repair; the H_ref and RHS
    # mathematics are unchanged from Agent 3 PR #703 exact head above.
    delta_y = np.linalg.solve(matrices, target[..., np.newaxis])[..., 0]
    if not np.all(np.isfinite(delta_y)):
        raise ArithmeticError("signed reference inverse produced non-finite coordinates")
    reconstructed = np.einsum("nij,nj->ni", matrices, delta_y)
    closure = float(np.max(np.abs(reconstructed - target)))
    scale = max(1.0, float(np.max(np.abs(target))))
    if closure > 1024.0 * math.ulp(scale):
        raise ArithmeticError("signed reference inverse reconstruction closure exceeded guard")

    return TypedMeanReferenceInverseEvaluation(
        identity=rhs.identity,
        rhs=rhs,
        operator=reference_operator,
        theta_axial_rhs=target,
        delta_y=delta_y,
        determinants=determinants,
        condition_numbers=condition_numbers,
        reconstruction_closure_max_abs=closure,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_typed_mean_reference_inverse)
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
        "reference_matrix_formula": REFERENCE_MATRIX_FORMULA,
        "signed_reference_solve_formula": SIGNED_REFERENCE_SOLVE_FORMULA,
        "source_formula_provenance_prs": [SOURCE_SIGNED_PAIR_PR, SOURCE_UNIT_BRIDGE_PR],
        "upstream_agent3_reference_inverse_pr": UPSTREAM_AGENT3_REFERENCE_INVERSE_PR,
        "upstream_agent3_reference_inverse_head": UPSTREAM_AGENT3_REFERENCE_INVERSE_HEAD,
        "a5_numpy2_compatibility_note": A5_NUMPY2_COMPATIBILITY_NOTE,
        "typed_actual_inverse_rhs_provider": (
            "kokuno_typed_mean_inverse_rhs:materialize_typed_mean_inverse_rhs"
        ),
        "raw_same_cycle_providers_required": True,
        "caller_supplied_inverse_rhs_allowed": False,
        "caller_supplied_signed_coordinates_allowed": False,
        "provenance_bearing_reference_operator_required": True,
        "source_reference_operator_certification_required_for_scientific_use": True,
        "source_reference_scale_certification_required_for_scientific_use": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "typed_reference_inverse_algebra_executable": True,
        "source_h_ref_materialized_for_full_candidate": False,
        "source_epsilon_materialized_for_full_candidate": False,
        "base_positive_amplitudes_materialized": False,
        "delta_amplitudes_materialized": False,
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
    """Non-f=R mechanics regression for the typed reference inverse only."""

    identity = CycleIdentity("analytic-typed-mean-reference-inverse-v1", 0, "t0")
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
    operator = MeanReferenceInverseOperator(
        identity=identity,
        radii=radii,
        A_c=tuple(2.0 for _ in radii),
        u_star=tuple(3.0 for _ in radii),
        h_plus=tuple(1.0 for _ in radii),
        h_minus=tuple(1.0 for _ in radii),
        producer_kind="analytic-regression",
        provenance="displayed H_ref formula with fixed A_c=2,u_*=3,h_+=h_-=1",
        source_reference_operator_certified=False,
    )
    evaluation = materialize_typed_mean_reference_inverse(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        reference_scale,
        operator,
        viscosity=0.01,
        spatial_step=0.005,
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "upstream_agent3_reference_inverse_pr": UPSTREAM_AGENT3_REFERENCE_INVERSE_PR,
            "upstream_agent3_reference_inverse_head": UPSTREAM_AGENT3_REFERENCE_INVERSE_HEAD,
            "a5_numpy2_compatibility_note": A5_NUMPY2_COMPATIBILITY_NOTE,
            "source_signed_pair_pr": SOURCE_SIGNED_PAIR_PR,
            "source_unit_bridge_pr": SOURCE_UNIT_BRIDGE_PR,
            "structural_source": "Kokuno corrected reconstruction signed two-column mean correction",
            "repository_role": "typed actual DeltaC/epsilon -> signed H_ref reference coordinates",
        },
        "analytic_regression": {
            "field": "u=t*psi(r)*(-y,x,1), u_t=psi(r)*(-y,x,1), p=0, f=0 at t=0",
            "reference_scale_kind": "analytic-regression; not source-certified",
            "reference_operator_kind": "analytic-regression; not source-certified",
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
