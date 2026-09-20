"""Route the strict-inner transport mean into the admitted compact radial stress.

Kokuno Agent 3 owns the mean-correction / compact-radial-stress lane.  This
module consumes the exact typed strict-inner transport mean produced by Agent-3
PR #850 from Agent-2 #840 + Agent-1 #839 and immediately applies the already
admitted compact moment-complement radial inverse

    M_e(F) = integral r^e F(r) dr,
    P_e(F) = F - b_e M_e(F),
    sigma_e(r) = -r^(-e) integral_0^r s^e P_e(F)(s) ds,

with e=2 for the tangential channel and e=1 for the axial channel.

The routed source is still only

    T = u_t + (u . grad)u - nu Delta u,  nu=0.01,

for the current strict-inner ``u_inner + u_osc`` field.  Matched pressure,
restricted forcing, the corrected/global leading join, and Agent-3 correction
transport are absent.  Consequently these stresses are scoped transport
diagnostics and are not authorized as a complete NS-defect correction target.
"""

from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_actual_oscillatory_mean_stress import _compact_radial_stress
from .kokuno_strict_inner_transport_mean import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    REPOSITORY_VISCOSITY,
    ExactAgent1InnerTransportBackend,
    ExactAgent2StrictInnerTransportBackend,
    StrictInnerTransportMeanWitness,
    _materialize_from_transport_provider,
    materialize_strict_inner_transport_mean,
)

TASK = "KOKUNO-A3-STRICT-INNER-TRANSPORT-RADIAL-STRESS-073"
SCHEMA = "kokuno-a3-strict-inner-transport-radial-stress-v1"
PARENT_AGENT3_PR = 850
PARENT_AGENT3_HEAD = "7c3e5c23df9125722f055df4696f593ea17dd0e8"
RADIAL_OPERATOR_AGENT3_PR = 677
RADIAL_OPERATOR_AGENT3_HEAD = "56024553b981833a283f98648c8599011d6f38ab"
RADIAL_OPERATOR_INDEPENDENT_ADMISSION_AGENT3_PR = 607
RADIAL_OPERATOR_INDEPENDENT_AUDIT_AGENT4_PR = 600

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"


@dataclass(frozen=True)
class StrictInnerTransportRadialGeometry:
    """Fixed-(t,z) radial profile used by the compact stress constructor."""

    time: float
    axial_z: float
    radii: tuple[float, ...]
    bump_center: float
    bump_halfwidth: float

    def __post_init__(self) -> None:
        radii = tuple(float(value) for value in self.radii)
        object.__setattr__(self, "radii", radii)
        scalars = (
            float(self.time),
            float(self.axial_z),
            float(self.bump_center),
            float(self.bump_halfwidth),
        )
        if not all(math.isfinite(value) for value in scalars):
            raise ValueError("radial geometry values must be finite")
        if len(radii) < 9:
            raise ValueError("radial geometry requires at least nine radii")
        if not all(math.isfinite(value) and value > 0.0 for value in radii):
            raise ValueError("radii must be finite and strictly positive")
        if not all(right > left for left, right in zip(radii, radii[1:])):
            raise ValueError("radii must be strictly increasing")
        if self.bump_halfwidth <= 0.0:
            raise ValueError("bump_halfwidth must be strictly positive")
        support_min = self.bump_center - self.bump_halfwidth
        support_max = self.bump_center + self.bump_halfwidth
        tolerance = 64.0 * math.ulp(
            max(1.0, abs(self.bump_center), abs(self.bump_halfwidth))
        )
        if support_min < radii[0] - tolerance or support_max > radii[-1] + tolerance:
            raise ValueError(
                "compact bump support must lie inside the sampled radial interval"
            )


def _clean_stress(report: dict[str, object]) -> dict[str, object]:
    cleaned: dict[str, object] = {}
    for key, value in report.items():
        if isinstance(value, np.ndarray):
            cleaned[key] = np.asarray(value, dtype=float).tolist()
        elif isinstance(value, np.floating):
            cleaned[key] = float(value)
        else:
            cleaned[key] = value
    return cleaned


@dataclass(frozen=True)
class StrictInnerTransportRadialStressWitness:
    """Scoped strict-inner transport mean and its compact theta/z stresses."""

    geometry: StrictInnerTransportRadialGeometry
    mean_witness: StrictInnerTransportMeanWitness
    radial_mean_profile: np.ndarray
    theta_mean_profile: np.ndarray
    axial_mean_profile: np.ndarray
    theta_stress: dict[str, object]
    axial_stress: dict[str, object]
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "geometry": asdict(self.geometry),
            "backend_kind": self.backend_kind,
            "mean_witness": self.mean_witness.to_receipt(),
            "radial_mean_profile": np.asarray(
                self.radial_mean_profile, dtype=float
            ).tolist(),
            "theta_mean_profile": np.asarray(
                self.theta_mean_profile, dtype=float
            ).tolist(),
            "axial_mean_profile": np.asarray(
                self.axial_mean_profile, dtype=float
            ).tolist(),
            "theta_e2": _clean_stress(self.theta_stress),
            "axial_e1": _clean_stress(self.axial_stress),
            "truth_boundary": truth_boundary(),
        }


def _materialize_from_mean_witness(
    witness: StrictInnerTransportMeanWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    backend_kind: str,
) -> StrictInnerTransportRadialStressWitness:
    """Apply the admitted radial operator to one already-typed A3 mean witness.

    This helper is intentionally internal.  The public scientific materializer
    below recomputes the mean witness from exact pinned upstream backends, so a
    caller cannot inject a precomputed mean profile as success evidence.
    """

    if not isinstance(witness, StrictInnerTransportMeanWitness):
        raise TypeError("witness must be StrictInnerTransportMeanWitness")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    radii = np.asarray(geometry.radii, dtype=float)
    expected_shape = (radii.size, 3)
    mean_transport = np.asarray(witness.mean_transport_cylindrical, dtype=float)
    if mean_transport.shape != expected_shape:
        raise ValueError(
            f"transport mean shape {mean_transport.shape} does not match {expected_shape}"
        )
    if not np.all(np.isfinite(mean_transport)):
        raise ValueError("transport mean contains non-finite values")

    witness_radii = np.asarray(witness.radius, dtype=float)
    witness_z = np.asarray(witness.z, dtype=float)
    witness_t = np.asarray(witness.t, dtype=float)
    if witness_radii.shape != radii.shape or not np.array_equal(witness_radii, radii):
        raise ValueError("mean witness radial grid does not match radial geometry")
    if witness_z.shape != radii.shape or not np.array_equal(
        witness_z, np.full_like(radii, float(geometry.axial_z))
    ):
        raise ValueError("mean witness z values do not match fixed radial geometry")
    if witness_t.shape != radii.shape or not np.array_equal(
        witness_t, np.full_like(radii, float(geometry.time))
    ):
        raise ValueError("mean witness time values do not match fixed radial geometry")
    if witness.viscosity != REPOSITORY_VISCOSITY:
        raise ValueError("mean witness viscosity drifted from repository value 0.01")
    if tuple(witness.angular_orders) != tuple(ANGULAR_ORDERS):
        raise ValueError("mean witness angular ladder drifted")

    radial = np.array(mean_transport[:, 0], copy=True)
    theta = np.array(mean_transport[:, 1], copy=True)
    axial = np.array(mean_transport[:, 2], copy=True)

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

    return StrictInnerTransportRadialStressWitness(
        geometry=geometry,
        mean_witness=witness,
        radial_mean_profile=radial,
        theta_mean_profile=theta,
        axial_mean_profile=axial,
        theta_stress=theta_stress,
        axial_stress=axial_stress,
        backend_kind=backend_kind,
    )


def materialize_strict_inner_transport_radial_stress(
    transport_backend: ExactAgent2StrictInnerTransportBackend,
    inner_backend: ExactAgent1InnerTransportBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> StrictInnerTransportRadialStressWitness:
    """Recompute the exact #850 scoped mean, then apply compact radial stress."""

    if not isinstance(transport_backend, ExactAgent2StrictInnerTransportBackend):
        raise TypeError(
            "transport_backend must be ExactAgent2StrictInnerTransportBackend"
        )
    if not isinstance(inner_backend, ExactAgent1InnerTransportBackend):
        raise TypeError("inner_backend must be ExactAgent1InnerTransportBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    radii = np.asarray(geometry.radii, dtype=float)
    witness = materialize_strict_inner_transport_mean(
        transport_backend,
        inner_backend,
        radii,
        np.full_like(radii, float(geometry.axial_z)),
        np.full_like(radii, float(geometry.time)),
    )
    return _materialize_from_mean_witness(
        witness,
        geometry,
        backend_kind="exact-a3-850-mean-to-admitted-radial-stress",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_strict_inner_transport_radial_stress)
    forbidden = {
        "residual",
        "defect",
        "mean",
        "source",
        "stress",
        "inverse",
        "pressure",
        "forcing",
        "target",
        "gain",
        "alpha",
        "damping",
        "angular_order",
        "spatial_step",
        "time_step",
        "viscosity",
        "nu",
        "delta_y",
        "delta_a",
        "normalized_score",
        "scientific_threshold",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
        "radial_operator_agent3_head": RADIAL_OPERATOR_AGENT3_HEAD,
        "radial_operator_independent_admission_agent3_pr": (
            RADIAL_OPERATOR_INDEPENDENT_ADMISSION_AGENT3_PR
        ),
        "radial_operator_independent_audit_agent4_pr": (
            RADIAL_OPERATOR_INDEPENDENT_AUDIT_AGENT4_PR
        ),
        "scoped_transport_formula": "u_t + (u dot grad)u - nu*Delta(u)",
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "angular_orders": list(ANGULAR_ORDERS),
        "radial_inverse_formula": (
            "M_e F=int r^e F; P_e F=F-b_e M_e F; "
            "sigma_e=-r^(-e) int_0^r s^e P_e F ds"
        ),
        "theta_exponent": 2,
        "axial_exponent": 1,
        "radial_component_recorded_but_not_inverted": True,
        "mean_recomputed_from_exact_upstream_in_public_api": True,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "compact_support_geometry_checked": True,
        "moment_complement_recorded": True,
        "strict_inner_transport_radial_stress_materialized": True,
        "complete_ns_defect": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "global_corrected_leading_join_materialized": False,
        "correction_transport_included": False,
        "scoped_transport_stress_authorized_as_correction_target": False,
        "finite_head_debt_authorized_from_scoped_transport": False,
        "real_agent3_delta_a_bound": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def _mechanics_transport_provider():
    """Manufactured rotating-frame transport with nonzero theta/z mean profiles."""

    def provider(x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray):
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        radius = np.sqrt(xb * xb + yb * yb)
        theta = np.arctan2(yb, xb)
        c = np.cos(theta)
        s = np.sin(theta)

        mr = 0.05 * radius
        mt = 0.20 * (1.0 - ((radius - 0.50) / 0.30) ** 2) ** 2
        mz = 0.12 + 0.03 * radius

        mode = np.cos(2.0 * theta)
        mr = mr + 0.04 * mode
        mt = mt - 0.03 * np.sin(2.0 * theta)
        mz = mz + 0.02 * mode

        transport = np.stack(
            (mr * c - mt * s, mr * s + mt * c, mz), axis=-1
        )
        zero = np.zeros_like(transport)
        return transport, zero, zero, transport, REPOSITORY_VISCOSITY

    return provider


def build_mechanics_report() -> dict[str, object]:
    """Manufactured interface/radial-operator regression, never candidate evidence."""

    radii = np.linspace(0.20, 0.80, 49)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(value) for value in radii),
        bump_center=0.50,
        bump_halfwidth=0.30,
    )
    mean_witness = _materialize_from_transport_provider(
        _mechanics_transport_provider(),
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        agent1_backend=None,
        agent2_backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    stress_witness = _materialize_from_mean_witness(
        mean_witness,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    expected_theta = 0.20 * (1.0 - ((radii - 0.50) / 0.30) ** 2) ** 2
    expected_axial = 0.12 + 0.03 * radii
    expected_radial = 0.05 * radii
    return {
        "schema": SCHEMA,
        "task": TASK,
        "mechanics_only": True,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
            "structural_source": (
                "Kokuno corrected reconstruction mean correction / "
                "compact radial stress architecture"
            ),
        },
        "expected_radial_mean_profile": expected_radial.tolist(),
        "expected_theta_mean_profile": expected_theta.tolist(),
        "expected_axial_mean_profile": expected_axial.tolist(),
        "evaluation": stress_witness.to_receipt(),
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_mechanics_report()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
