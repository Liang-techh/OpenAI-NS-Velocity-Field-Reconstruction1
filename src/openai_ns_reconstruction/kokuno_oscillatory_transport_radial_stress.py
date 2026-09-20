"""Route the typed oscillatory transport mean through the compact radial inverse.

Agent 2 PR #901 exposes the strict-inner Cartesian oscillatory transport increment

    delta T_osc = d_t u_osc
                  + (u_inner . grad) u_osc
                  + (u_osc . grad) u_inner
                  + (u_osc . grad) u_osc
                  - 0.01 Delta u_osc.

Agent 3 PR #902 projects that typed quantity, together with T_inner and
T_(inner+osc), onto the rotating cylindrical m=0 frame.  This module performs
one next Agent-3-owned step: on a fixed-(z,t) radial profile, apply the already
admitted polynomial-first-cell compact radial inverse from PR #875 to the
cylindrical tangential (e=2) and axial (e=1) means.

The identical linear radial operator is also applied to the #902 inner and
inner+oscillatory transport means, allowing the attribution identity

    sigma(delta T_osc) = sigma(T_(inner+osc)) - sigma(T_inner)

to be checked at stress level on the same radial grid and compact moment
functional.  The radial cylindrical component is recorded but is not inverted.

This remains strict-inner pressure/forcing-free transport.  It is not a
complete NS defect and is not authorized as a real correction target.
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

from .kokuno_oscillatory_transport_mean import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    REPOSITORY_VISCOSITY,
    ExactAgent2OscillatoryTransportHandoff,
    OscillatoryTransportMeanWitness,
    _materialize_from_provider,
    materialize_oscillatory_transport_mean,
)
from .kokuno_strict_inner_transport_first_cell_polynomial_radial_stress import (
    _polynomial_first_cell_compact_radial_stress,
)
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)

TASK = "KOKUNO-A3-OSCILLATORY-TRANSPORT-RADIAL-STRESS-085"
SCHEMA = "kokuno-a3-oscillatory-transport-radial-stress-v1"
PARENT_AGENT3_PR = 902
PARENT_AGENT3_HEAD = "baa36f77039423ec5c6b86a8d0a60da5d2fa9811"
PARENT_AGENT3_SOURCE_BLOB = "08933745fa25d0f1a8a136843f855cee7675bf75"
RADIAL_OPERATOR_AGENT3_PR = 875
RADIAL_OPERATOR_AGENT3_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

MEAN_ATTRIBUTION_ABSOLUTE_GATE = 1.0e-12
STRESS_ATTRIBUTION_RELATIVE_GATE = 2.0e-11


def _clean_report(report: dict[str, object]) -> dict[str, object]:
    cleaned: dict[str, object] = {}
    for key, value in report.items():
        if isinstance(value, np.ndarray):
            cleaned[key] = np.asarray(value, dtype=float).tolist()
        elif isinstance(value, np.floating):
            cleaned[key] = float(value)
        else:
            cleaned[key] = value
    return cleaned


def _relative_max_closure(lhs: np.ndarray, rhs: np.ndarray) -> float:
    a = np.asarray(lhs, dtype=float)
    b = np.asarray(rhs, dtype=float)
    scale = max(1.0, float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    return float(np.max(np.abs(a - b)) / scale)


def _stress_array(report: dict[str, object], label: str) -> np.ndarray:
    value = np.asarray(report.get("stress"), dtype=float)
    if value.ndim != 1 or value.size < 9 or not np.all(np.isfinite(value)):
        raise ValueError(f"{label} stress is malformed or non-finite")
    return value


@dataclass(frozen=True)
class OscillatoryTransportRadialStressWitness:
    """Scoped oscillatory transport mean plus theta/z compact radial stresses."""

    geometry: StrictInnerTransportRadialGeometry
    mean_witness: OscillatoryTransportMeanWitness
    radial_delta_mean_profile: np.ndarray
    theta_delta_mean_profile: np.ndarray
    axial_delta_mean_profile: np.ndarray
    theta_delta_stress: dict[str, object]
    theta_inner_stress: dict[str, object]
    theta_combined_stress: dict[str, object]
    axial_delta_stress: dict[str, object]
    axial_inner_stress: dict[str, object]
    axial_combined_stress: dict[str, object]
    mean_attribution_closure_absolute_max: float
    theta_stress_attribution_relative_max: float
    axial_stress_attribution_relative_max: float
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "geometry": asdict(self.geometry),
            "backend_kind": self.backend_kind,
            "mean_witness": self.mean_witness.to_receipt(),
            "radial_delta_mean_profile": np.asarray(
                self.radial_delta_mean_profile, dtype=float
            ).tolist(),
            "theta_delta_mean_profile": np.asarray(
                self.theta_delta_mean_profile, dtype=float
            ).tolist(),
            "axial_delta_mean_profile": np.asarray(
                self.axial_delta_mean_profile, dtype=float
            ).tolist(),
            "theta_e2": {
                "delta": _clean_report(self.theta_delta_stress),
                "inner": _clean_report(self.theta_inner_stress),
                "inner_plus_oscillatory": _clean_report(self.theta_combined_stress),
                "stress_attribution_relative_max": self.theta_stress_attribution_relative_max,
            },
            "axial_e1": {
                "delta": _clean_report(self.axial_delta_stress),
                "inner": _clean_report(self.axial_inner_stress),
                "inner_plus_oscillatory": _clean_report(self.axial_combined_stress),
                "stress_attribution_relative_max": self.axial_stress_attribution_relative_max,
            },
            "mean_attribution_closure_absolute_max": self.mean_attribution_closure_absolute_max,
            "truth_boundary": truth_boundary(),
        }


def _apply_channel(
    radii: np.ndarray,
    delta: np.ndarray,
    inner: np.ndarray,
    combined: np.ndarray,
    *,
    exponent: int,
    geometry: StrictInnerTransportRadialGeometry,
) -> tuple[dict[str, object], dict[str, object], dict[str, object], float]:
    kwargs = {
        "exponent": exponent,
        "bump_center": geometry.bump_center,
        "bump_halfwidth": geometry.bump_halfwidth,
    }
    delta_report = _polynomial_first_cell_compact_radial_stress(
        radii, delta, **kwargs
    )
    inner_report = _polynomial_first_cell_compact_radial_stress(
        radii, inner, **kwargs
    )
    combined_report = _polynomial_first_cell_compact_radial_stress(
        radii, combined, **kwargs
    )
    sigma_delta = _stress_array(delta_report, "delta")
    sigma_inner = _stress_array(inner_report, "inner")
    sigma_combined = _stress_array(combined_report, "combined")
    closure = _relative_max_closure(sigma_delta, sigma_combined - sigma_inner)
    return delta_report, inner_report, combined_report, closure


def _materialize_from_mean_witness(
    witness: OscillatoryTransportMeanWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    backend_kind: str,
) -> OscillatoryTransportRadialStressWitness:
    if not isinstance(witness, OscillatoryTransportMeanWitness):
        raise TypeError("witness must be OscillatoryTransportMeanWitness")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    radii = np.asarray(geometry.radii, dtype=float)
    expected = (radii.size, 3)
    arrays = {
        "delta": np.asarray(witness.mean_transport_increment_cylindrical, dtype=float),
        "inner": np.asarray(witness.mean_inner_transport_cylindrical, dtype=float),
        "combined": np.asarray(
            witness.mean_inner_plus_oscillatory_transport_cylindrical, dtype=float
        ),
    }
    for label, value in arrays.items():
        if value.shape != expected:
            raise ValueError(f"{label} mean shape {value.shape} does not match {expected}")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{label} mean contains non-finite values")

    witness_r = np.asarray(witness.radius, dtype=float)
    witness_z = np.asarray(witness.z, dtype=float)
    witness_t = np.asarray(witness.t, dtype=float)
    if witness_r.shape != radii.shape or not np.array_equal(witness_r, radii):
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

    mean_closure = float(np.max(np.abs(arrays["delta"] - (arrays["combined"] - arrays["inner"]))))
    if mean_closure > MEAN_ATTRIBUTION_ABSOLUTE_GATE:
        raise ValueError("oscillatory mean before/after attribution closure failed")

    theta_delta, theta_inner, theta_combined, theta_closure = _apply_channel(
        radii,
        arrays["delta"][:, 1],
        arrays["inner"][:, 1],
        arrays["combined"][:, 1],
        exponent=2,
        geometry=geometry,
    )
    axial_delta, axial_inner, axial_combined, axial_closure = _apply_channel(
        radii,
        arrays["delta"][:, 2],
        arrays["inner"][:, 2],
        arrays["combined"][:, 2],
        exponent=1,
        geometry=geometry,
    )

    return OscillatoryTransportRadialStressWitness(
        geometry=geometry,
        mean_witness=witness,
        radial_delta_mean_profile=np.array(arrays["delta"][:, 0], copy=True),
        theta_delta_mean_profile=np.array(arrays["delta"][:, 1], copy=True),
        axial_delta_mean_profile=np.array(arrays["delta"][:, 2], copy=True),
        theta_delta_stress=theta_delta,
        theta_inner_stress=theta_inner,
        theta_combined_stress=theta_combined,
        axial_delta_stress=axial_delta,
        axial_inner_stress=axial_inner,
        axial_combined_stress=axial_combined,
        mean_attribution_closure_absolute_max=mean_closure,
        theta_stress_attribution_relative_max=theta_closure,
        axial_stress_attribution_relative_max=axial_closure,
        backend_kind=backend_kind,
    )


def materialize_oscillatory_transport_radial_stress(
    handoff_backend: ExactAgent2OscillatoryTransportHandoff,
    geometry: StrictInnerTransportRadialGeometry,
) -> OscillatoryTransportRadialStressWitness:
    """Recompute the exact #902 mean and route it through the #875 radial inverse."""
    if not isinstance(handoff_backend, ExactAgent2OscillatoryTransportHandoff):
        raise TypeError("handoff_backend must be ExactAgent2OscillatoryTransportHandoff")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    radii = np.asarray(geometry.radii, dtype=float)
    witness = materialize_oscillatory_transport_mean(
        handoff_backend,
        radii,
        np.full_like(radii, float(geometry.axial_z)),
        np.full_like(radii, float(geometry.time)),
    )
    return _materialize_from_mean_witness(
        witness,
        geometry,
        backend_kind="exact-a3-902-oscillatory-mean-to-a3-875-radial-stress",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_oscillatory_transport_radial_stress)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse",
        "pressure", "forcing", "target", "gain", "alpha", "damping",
        "angular_order", "spatial_step", "time_step", "viscosity", "nu",
        "delta_y", "delta_a", "normalized_score", "scientific_threshold",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
        "radial_operator_agent3_head": RADIAL_OPERATOR_AGENT3_HEAD,
        "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "angular_orders": list(ANGULAR_ORDERS),
        "strict_inner_only": True,
        "oscillatory_transport_increment_materialized": True,
        "cylindrical_m0_projection_materialized": True,
        "scoped_oscillatory_transport_radial_inverse_performed": True,
        "radial_inverse_rule": "polynomial_first_cell_compact_moment_complement",
        "theta_exponent": 2,
        "axial_exponent": 1,
        "radial_component_recorded_but_not_inverted": True,
        "stress_level_before_after_attribution_checked": True,
        "mean_attribution_absolute_gate": MEAN_ATTRIBUTION_ABSOLUTE_GATE,
        "stress_attribution_relative_gate": STRESS_ATTRIBUTION_RELATIVE_GATE,
        "mean_recomputed_from_exact_upstream_in_public_api": True,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "global_corrected_leading_join_materialized": False,
        "complete_ns_defect": False,
        "scoped_transport_stress_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
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


def _cartesian_from_cylindrical(
    radial: np.ndarray,
    tangential: np.ndarray,
    axial: np.ndarray,
    theta: np.ndarray,
) -> np.ndarray:
    c = np.cos(theta)
    s = np.sin(theta)
    return np.stack((radial * c - tangential * s, radial * s + tangential * c, axial), axis=-1)


def _mechanics_provider():
    """Manufactured linear-attribution fixture; never candidate evidence."""
    def provider(x: np.ndarray, y: np.ndarray, z: np.ndarray, t: np.ndarray):
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        r = np.sqrt(xb * xb + yb * yb)
        theta = np.arctan2(yb, xb)
        mode_c = np.cos(2.0 * theta)
        mode_s = np.sin(2.0 * theta)

        inner = _cartesian_from_cylindrical(
            0.03 * r + 0.01 * mode_c,
            0.11 + 0.04 * r + 0.02 * mode_s,
            -0.06 + 0.03 * r + 0.01 * mode_c,
            theta,
        )
        delta = _cartesian_from_cylindrical(
            -0.02 * r + 0.015 * mode_s,
            0.08 * (1.0 + r) - 0.01 * mode_c,
            0.05 - 0.02 * r + 0.008 * mode_s,
            theta,
        )
        combined = inner + delta
        time_part = 0.25 * delta
        nonlinear = 0.50 * delta
        viscous = 0.25 * delta
        return time_part, nonlinear, viscous, delta, inner, combined, REPOSITORY_VISCOSITY
    return provider


def build_mechanics_report() -> dict[str, object]:
    radii = np.linspace(0.02, 0.80, 79)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(v) for v in radii),
        bump_center=0.50,
        bump_halfwidth=0.28,
    )
    mean = _materialize_from_provider(
        _mechanics_provider(),
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        backend=None,
        backend_kind="manufactured-mechanics-only",
    )
    result = _materialize_from_mean_witness(
        mean, geometry, backend_kind="manufactured-mechanics-only"
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "witness": result.to_receipt(),
        "truth_boundary": truth_boundary(),
    }


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_mechanics_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
