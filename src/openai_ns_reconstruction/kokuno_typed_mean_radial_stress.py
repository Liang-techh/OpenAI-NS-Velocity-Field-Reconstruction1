"""Connect the typed actual mean defect directly to compact radial stress.

This is the narrow Agent-3 seam after ``kokuno_typed_mean_defect_projection``.
Callers provide raw same-cycle field providers and physical radial geometry; no
precomputed residual/defect/stress/target/gain enters the public API.  At each
radius the module recomputes

    R = u_t + (u . grad)u + grad(p) - nu Delta(u) - f,

projects its angular mean channels, and sends only the resulting actual
``bar R_theta`` and ``bar R_z`` profiles into the already-admitted compact
moment-complement radial inverse

    M_e F = int r^e F dr,
    P_e F = F - b_e M_e F,
    sigma_e(r) = -r^(-e) int_0^r s^e P_e(s) ds,

with e=2 for theta and e=1 for axial.  The existing radial constructor is reused
rather than reimplemented or retuned.

This adapter is repository engineering machinery.  It does not independently
validate restricted-forcing semantics, create a correction velocity, run a
finite correction cycle, assess the final normalized NS residual, or prove a
blow-up theorem.
"""

from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .kokuno_compact_radial_stress_adapter import _compact_radial_stress
from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from .kokuno_typed_mean_defect_projection import (
    CylindricalRing,
    MeanDefectRingEvaluation,
    evaluate_actual_mean_defect_ring,
)

TASK = "KOKUNO-A3-TYPED-MEAN-RADIAL-STRESS-052"
SCHEMA = "kokuno-a3-typed-mean-radial-stress-v1"
PARENT_AGENT3_PR = 670
PARENT_AGENT3_HEAD = "56833156dc0db1ac03ceb42bf49fa2c55935c85e"


@dataclass(frozen=True)
class RadialMeanStressGeometry:
    """One fixed-(t,z) radial profile and its compact correction bump."""

    time: float
    axial_z: float
    radii: tuple[float, ...]
    bump_center: float
    bump_halfwidth: float
    angular_count: int = 32
    phase: float = 0.0

    def __post_init__(self) -> None:
        radii = tuple(float(value) for value in self.radii)
        object.__setattr__(self, "radii", radii)
        scalars = (self.time, self.axial_z, self.bump_center, self.bump_halfwidth, self.phase)
        if not all(math.isfinite(float(value)) for value in scalars):
            raise ValueError("radial-stress geometry must be finite")
        if len(radii) < 9:
            raise ValueError("radial-stress geometry requires at least nine radii")
        if not all(math.isfinite(value) and value > 0.0 for value in radii):
            raise ValueError("radii must be finite and strictly positive")
        if not all(right > left for left, right in zip(radii, radii[1:])):
            raise ValueError("radii must be strictly increasing")
        if self.bump_halfwidth <= 0.0:
            raise ValueError("bump_halfwidth must be strictly positive")
        if self.angular_count < 8:
            raise ValueError("angular_count must be at least 8")
        support_min = self.bump_center - self.bump_halfwidth
        support_max = self.bump_center + self.bump_halfwidth
        tolerance = 64.0 * math.ulp(max(1.0, abs(self.bump_center), abs(self.bump_halfwidth)))
        if support_min < radii[0] - tolerance or support_max > radii[-1] + tolerance:
            raise ValueError("compact bump support must lie inside the sampled radial interval")


@dataclass(frozen=True)
class TypedMeanRadialStressEvaluation:
    """Actual mean-defect profiles and their compact radial stresses."""

    identity: CycleIdentity
    geometry: RadialMeanStressGeometry
    viscosity: float
    spatial_step: float
    rings: tuple[MeanDefectRingEvaluation, ...]
    radial_mean_profile: np.ndarray
    theta_mean_profile: np.ndarray
    axial_mean_profile: np.ndarray
    theta_stress: dict[str, object]
    axial_stress: dict[str, object]

    def to_receipt(self) -> dict[str, object]:
        def clean_stress(report: dict[str, object]) -> dict[str, object]:
            cleaned: dict[str, object] = {}
            for key, value in report.items():
                if isinstance(value, np.ndarray):
                    cleaned[key] = value.tolist()
                elif isinstance(value, np.floating):
                    cleaned[key] = float(value)
                else:
                    cleaned[key] = value
            return cleaned

        return {
            "identity": asdict(self.identity),
            "geometry": asdict(self.geometry),
            "viscosity": self.viscosity,
            "spatial_step": self.spatial_step,
            "rings": [ring.to_receipt() for ring in self.rings],
            "radial_mean_profile": self.radial_mean_profile.tolist(),
            "theta_mean_profile": self.theta_mean_profile.tolist(),
            "axial_mean_profile": self.axial_mean_profile.tolist(),
            "theta_e2": clean_stress(self.theta_stress),
            "axial_e1": clean_stress(self.axial_stress),
        }


def materialize_actual_mean_radial_stress(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    geometry: RadialMeanStressGeometry,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanRadialStressEvaluation:
    """Recompute actual mean defect and immediately apply the admitted inverse.

    The source profiles are intentionally not public arguments.  They are
    materialized internally from raw providers through #670's typed mean
    projection, then passed unchanged to the existing compact stress operator.
    """

    rings = tuple(
        evaluate_actual_mean_defect_ring(
            velocity,
            velocity_dt,
            pressure,
            restricted_forcing,
            CylindricalRing(
                time=geometry.time,
                radius=radius,
                axial_z=geometry.axial_z,
                angular_count=geometry.angular_count,
                phase=geometry.phase,
            ),
            viscosity=viscosity,
            spatial_step=spatial_step,
        )
        for radius in geometry.radii
    )
    identities = {item.identity for item in rings}
    if len(identities) != 1:
        raise RuntimeError("typed radial profile produced mixed cycle identities")

    radial = np.asarray([item.radial_mean for item in rings], dtype=float)
    theta = np.asarray([item.theta_mean for item in rings], dtype=float)
    axial = np.asarray([item.axial_mean for item in rings], dtype=float)
    radii = np.asarray(geometry.radii, dtype=float)

    theta_stress = _compact_radial_stress(
        radii,
        theta,
        exponent=2,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    axial_stress = _compact_radial_stress(
        radii,
        axial,
        exponent=1,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    return TypedMeanRadialStressEvaluation(
        identity=next(iter(identities)),
        geometry=geometry,
        viscosity=float(viscosity),
        spatial_step=float(spatial_step),
        rings=rings,
        radial_mean_profile=radial,
        theta_mean_profile=theta,
        axial_mean_profile=axial,
        theta_stress=theta_stress,
        axial_stress=axial_stress,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_actual_mean_radial_stress)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "source",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
    }
    return {
        "source_defect_formula": "u_t + (u dot grad)u + grad(p) - nu*Delta(u) - f",
        "actual_mean_projection_provider": "kokuno_typed_mean_defect_projection:evaluate_actual_mean_defect_ring",
        "radial_stress_constructor": "kokuno_compact_radial_stress_adapter:_compact_radial_stress (exact A3 operator extraction)",
        "theta_exponent": 2,
        "axial_exponent": 1,
        "radial_inverse_formula": "sigma_e=-r^(-e) integral r^e(F-b_e*M_e(F)) dr",
        "raw_same_cycle_providers_required": True,
        "caller_supplied_precomputed_defect_allowed": False,
        "caller_supplied_mean_source_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "actual_mean_channels_fed_directly_to_existing_radial_operator": True,
        "compact_support_geometry_checked": True,
        "restricted_forcing_semantics_independently_validated_here": False,
        "radial_operator_revalidated_here": False,
        "kokuno_theorem_machine_replayed_here": False,
        "real_full_candidate_defect_consumed": False,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
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
    """Nonzero typed-seam regression with fixed zero forcing, never ``f=R``."""

    identity = CycleIdentity("analytic-typed-radial-stress-v1", 0, "t0")
    center = 0.50
    halfwidth = 0.30

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
        "regression-zero-forcing; independent of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    radii = tuple(float(value) for value in np.linspace(0.20, 0.80, 49))
    geometry = RadialMeanStressGeometry(
        time=0.0,
        axial_z=0.13,
        radii=radii,
        bump_center=center,
        bump_halfwidth=halfwidth,
        angular_count=32,
        phase=0.017,
    )
    evaluation = materialize_actual_mean_radial_stress(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        viscosity=0.01,
        spatial_step=0.005,
    )
    expected_psi = np.asarray(
        [_compact_profile(radius, center=center, halfwidth=halfwidth) for radius in radii],
        dtype=float,
    )
    expected_theta = np.asarray(radii, dtype=float) * expected_psi
    theta_error = float(np.max(np.abs(evaluation.theta_mean_profile - expected_theta)))
    axial_error = float(np.max(np.abs(evaluation.axial_mean_profile - expected_psi)))
    radial_error = float(np.max(np.abs(evaluation.radial_mean_profile)))
    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "structural_source": "Kokuno corrected reconstruction mean correction / compact radial stress architecture",
            "repository_role": "actual typed mean defect -> existing compact radial-stress operator",
        },
        "analytic_regression": {
            "field": "u=t*psi(r)*(-y,x,1), u_t=psi(r)*(-y,x,1), p=0, f=0 at t=0",
            "expected_actual_defect": "R=psi(r)*(-y,x,1)",
            "expected_mean_channels": "bar_R_theta=r*psi(r), bar_R_z=psi(r), bar_R_r=0",
            "maximum_theta_mean_error": theta_error,
            "maximum_axial_mean_error": axial_error,
            "maximum_radial_mean_error": radial_error,
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
