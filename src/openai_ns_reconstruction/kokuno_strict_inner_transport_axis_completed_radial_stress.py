"""Axis-complete the scoped Kokuno radial-stress lower-limit discretization.

Agent-3 PR #858 routes the exact typed strict-inner transport mean through the
compact moment-complement radial inverse

    M_e(F) = int r^e F(r) dr,
    P_e(F) = F - b_e M_e(F),
    sigma_e(r) = -r^(-e) int_0^r s^e P_e(F)(s) ds,

with e=2 for the tangential channel and e=1 for the axial channel.  The reused
discrete helper in #858 receives only strictly positive radii and initializes
its cumulative trapezoid at the first supplied radius.  Consequently its first
reported stress is zero by construction and does not include the formal
``0 -> r_min`` primitive.

This increment keeps the same typed source and the same compact bump but uses an
explicit *axis-completed discrete weighted primitive*: for e>=1 and a bounded
source, the weighted integrand r^e P_e(F) has the formal endpoint value zero at
r=0.  We prepend that weighted endpoint to the quadrature before computing the
moment, moment complement, and cumulative primitive.  The source itself is not
invented at the axis, and no caller may supply a defect/mean/stress surrogate.

The result fixes only the discrete lower-limit realization on the current
strict-inner sampled radial interval.  It does not certify the unsampled global
radial domain, axis regularity of the future global candidate, matched pressure
or forcing, a complete NS defect, or a finite correction-cycle improvement.
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

from .kokuno_actual_oscillatory_mean_stress import (
    _compact_cos8_bump,
    _compact_radial_stress,
    _scalar_rms,
)
from .kokuno_strict_inner_transport_mean import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    REPOSITORY_VISCOSITY,
    ExactAgent1InnerTransportBackend,
    ExactAgent2StrictInnerTransportBackend,
    StrictInnerTransportMeanWitness,
    _materialize_from_transport_provider,
)
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
    StrictInnerTransportRadialStressWitness,
    _materialize_from_mean_witness as _materialize_legacy_from_mean_witness,
    _mechanics_transport_provider,
    materialize_strict_inner_transport_radial_stress,
)

TASK = "KOKUNO-A3-AXIS-COMPLETED-RADIAL-LOWER-LIMIT-074"
SCHEMA = "kokuno-a3-axis-completed-radial-lower-limit-v1"
PARENT_AGENT3_PR = 858
PARENT_AGENT3_HEAD = "9140b0bab1bfc34d496460675e92e63fc6dd721b"
INDEPENDENT_SCOPE_AUDIT_PR = 863
INDEPENDENT_SCOPE_AUDIT_HEAD = "7dddb17c9a853f0aa328c251050901104628ba9c"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"


def _axis_completed_cumulative_trapezoid(
    weighted_values: np.ndarray, radii: np.ndarray
) -> np.ndarray:
    """Integrate a weighted radial integrand from the formal lower limit zero.

    ``weighted_values`` means the already-weighted quantity ``r^e G(r)`` on a
    strictly-positive radial grid.  For the admitted e=1,2 channels, bounded G
    implies the formal weighted endpoint is zero.  This helper prepends exactly
    that endpoint and returns the cumulative primitive on the positive nodes.

    This is a numerical quadrature convention, not an independent certificate
    that the future global source is bounded/regular all the way to the axis.
    """

    values = np.asarray(weighted_values, dtype=float)
    r = np.asarray(radii, dtype=float)
    if values.ndim != 1 or r.ndim != 1 or values.shape != r.shape:
        raise ValueError("axis-completed primitive requires matching 1D arrays")
    if r.size < 3 or np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("axis-completed radii must be positive and strictly increasing")
    if not (np.all(np.isfinite(values)) and np.all(np.isfinite(r))):
        raise ValueError("axis-completed primitive inputs must be finite")

    augmented_r = np.concatenate((np.array([0.0]), r))
    augmented_values = np.concatenate((np.array([0.0]), values))
    increments = 0.5 * (augmented_values[1:] + augmented_values[:-1]) * np.diff(
        augmented_r
    )
    return np.cumsum(increments)


def _axis_completed_compact_radial_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    bump_center: float,
    bump_halfwidth: float,
) -> dict[str, Any]:
    """Apply the moment-complement inverse with an explicit discrete r=0 endpoint.

    The compact bump remains supported inside the supplied positive-radius
    interval, so no bump value is needed on ``0 <= r < r_min``.  The formal
    weighted endpoint is inserted as zero.  For transparency the legacy #858
    positive-radius-truncated result is computed from the same source and the
    delta is returned as a diagnostic.
    """

    r = np.asarray(radii, dtype=float)
    f = np.asarray(source, dtype=float)
    if r.ndim != 1 or f.shape != r.shape or r.size < 9:
        raise ValueError(
            "axis-completed radial stress requires matching 1D arrays with >=9 nodes"
        )
    if exponent not in (1, 2):
        raise ValueError("only the axial e=1 and tangential e=2 channels are supported")
    if np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("radii must be positive and strictly increasing")
    if not (np.all(np.isfinite(r)) and np.all(np.isfinite(f))):
        raise ValueError("radii and source must be finite")
    center = float(bump_center)
    halfwidth = float(bump_halfwidth)
    if not (math.isfinite(center) and math.isfinite(halfwidth) and halfwidth > 0.0):
        raise ValueError("compact bump geometry must be finite with positive halfwidth")
    support_min = center - halfwidth
    support_max = center + halfwidth
    tolerance = 64.0 * math.ulp(max(1.0, abs(center), abs(halfwidth)))
    if support_min < r[0] - tolerance or support_max > r[-1] + tolerance:
        raise ValueError("compact bump support must remain inside positive sampled radii")

    weight = r**exponent
    source_weighted = weight * f
    source_primitive = _axis_completed_cumulative_trapezoid(source_weighted, r)
    moment = float(source_primitive[-1])

    bump_raw = _compact_cos8_bump(r, center, halfwidth)
    bump_weighted_raw = weight * bump_raw
    bump_weighted_integral = float(
        _axis_completed_cumulative_trapezoid(bump_weighted_raw, r)[-1]
    )
    if not math.isfinite(bump_weighted_integral) or bump_weighted_integral <= 0.0:
        raise RuntimeError("axis-completed compact bump lost positive normalization")
    bump = bump_raw / bump_weighted_integral

    complement = f - bump * moment
    complement_weighted = weight * complement
    primitive = _axis_completed_cumulative_trapezoid(complement_weighted, r)
    complement_moment = float(primitive[-1])
    stress = -primitive / weight
    reconstructed_force = -f + bump * moment

    legacy = _compact_radial_stress(
        r,
        f,
        exponent=exponent,
        bump_center=center,
        bump_halfwidth=halfwidth,
    )
    legacy_stress = np.asarray(legacy["stress"], dtype=float)
    stress_delta = stress - legacy_stress

    return {
        "exponent": int(exponent),
        "quadrature_lower_limit": 0.0,
        "first_positive_radius": float(r[0]),
        "virtual_axis_weighted_integrand": 0.0,
        "weighted_moment": moment,
        "bump_weighted_integral": float(
            _axis_completed_cumulative_trapezoid(weight * bump, r)[-1]
        ),
        "moment_complement_weighted_moment": complement_moment,
        "axis_to_first_radius_weighted_primitive": float(primitive[0]),
        "stress": stress,
        "stress_rms": _scalar_rms(stress),
        "stress_max_abs": float(np.max(np.abs(stress))),
        "stress_inner_positive_edge": float(stress[0]),
        "stress_outer_edge": float(stress[-1]),
        "reconstructed_force_rms": _scalar_rms(reconstructed_force),
        "legacy_truncated_stress_inner_edge": float(legacy_stress[0]),
        "legacy_truncated_stress_outer_edge": float(legacy_stress[-1]),
        "stress_delta_vs_legacy_truncated": stress_delta,
        "stress_delta_vs_legacy_truncated_rms": _scalar_rms(stress_delta),
        "stress_delta_vs_legacy_truncated_max_abs": float(
            np.max(np.abs(stress_delta))
        ),
        "legacy_truncated_weighted_moment": float(legacy["weighted_moment"]),
        "legacy_truncated_moment_complement_weighted_moment": float(
            legacy["moment_complement_weighted_moment"]
        ),
    }


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


@dataclass(frozen=True)
class AxisCompletedStrictInnerTransportRadialStressWitness:
    """#858 scoped stress plus an explicit discrete 0->r lower-limit realization."""

    geometry: StrictInnerTransportRadialGeometry
    legacy_witness: StrictInnerTransportRadialStressWitness
    theta_axis_completed_stress: dict[str, object]
    axial_axis_completed_stress: dict[str, object]
    backend_kind: str

    @property
    def mean_witness(self) -> StrictInnerTransportMeanWitness:
        return self.legacy_witness.mean_witness

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "geometry": asdict(self.geometry),
            "backend_kind": self.backend_kind,
            "mean_witness": self.mean_witness.to_receipt(),
            "radial_mean_profile": np.asarray(
                self.legacy_witness.radial_mean_profile, dtype=float
            ).tolist(),
            "theta_mean_profile": np.asarray(
                self.legacy_witness.theta_mean_profile, dtype=float
            ).tolist(),
            "axial_mean_profile": np.asarray(
                self.legacy_witness.axial_mean_profile, dtype=float
            ).tolist(),
            "theta_e2_axis_completed": _clean_report(
                self.theta_axis_completed_stress
            ),
            "axial_e1_axis_completed": _clean_report(self.axial_axis_completed_stress),
            "legacy_positive_radius_truncated": {
                "theta_e2": _clean_report(self.legacy_witness.theta_stress),
                "axial_e1": _clean_report(self.legacy_witness.axial_stress),
            },
            "truth_boundary": truth_boundary(),
        }


def _materialize_axis_completed_from_mean_witness(
    witness: StrictInnerTransportMeanWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    backend_kind: str,
) -> AxisCompletedStrictInnerTransportRadialStressWitness:
    """Internal mechanics helper; public scientific path recomputes exact upstream."""

    legacy = _materialize_legacy_from_mean_witness(
        witness,
        geometry,
        backend_kind=backend_kind + ":legacy-858-check",
    )
    radii = np.asarray(geometry.radii, dtype=float)
    theta = np.asarray(legacy.theta_mean_profile, dtype=float)
    axial = np.asarray(legacy.axial_mean_profile, dtype=float)
    theta_stress = _axis_completed_compact_radial_stress(
        radii,
        theta,
        exponent=2,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    axial_stress = _axis_completed_compact_radial_stress(
        radii,
        axial,
        exponent=1,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    return AxisCompletedStrictInnerTransportRadialStressWitness(
        geometry=geometry,
        legacy_witness=legacy,
        theta_axis_completed_stress=theta_stress,
        axial_axis_completed_stress=axial_stress,
        backend_kind=backend_kind,
    )


def materialize_axis_completed_strict_inner_transport_radial_stress(
    transport_backend: ExactAgent2StrictInnerTransportBackend,
    inner_backend: ExactAgent1InnerTransportBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> AxisCompletedStrictInnerTransportRadialStressWitness:
    """Recompute exact typed #858 source, then axis-complete the discrete primitive."""

    if not isinstance(transport_backend, ExactAgent2StrictInnerTransportBackend):
        raise TypeError(
            "transport_backend must be ExactAgent2StrictInnerTransportBackend"
        )
    if not isinstance(inner_backend, ExactAgent1InnerTransportBackend):
        raise TypeError("inner_backend must be ExactAgent1InnerTransportBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    legacy = materialize_strict_inner_transport_radial_stress(
        transport_backend,
        inner_backend,
        geometry,
    )
    radii = np.asarray(geometry.radii, dtype=float)
    theta_stress = _axis_completed_compact_radial_stress(
        radii,
        np.asarray(legacy.theta_mean_profile, dtype=float),
        exponent=2,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    axial_stress = _axis_completed_compact_radial_stress(
        radii,
        np.asarray(legacy.axial_mean_profile, dtype=float),
        exponent=1,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    return AxisCompletedStrictInnerTransportRadialStressWitness(
        geometry=geometry,
        legacy_witness=legacy,
        theta_axis_completed_stress=theta_stress,
        axial_axis_completed_stress=axial_stress,
        backend_kind="exact-a3-858-typed-source-axis-completed-discrete-lower-limit",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_axis_completed_strict_inner_transport_radial_stress
    )
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
        "independent_scope_audit_pr": INDEPENDENT_SCOPE_AUDIT_PR,
        "independent_scope_audit_head": INDEPENDENT_SCOPE_AUDIT_HEAD,
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "angular_orders": list(ANGULAR_ORDERS),
        "theta_exponent": 2,
        "axial_exponent": 1,
        "axis_completed_discrete_lower_limit_node_included": True,
        "virtual_axis_weighted_integrand_zero_for_e_ge_1": True,
        "axis_to_first_positive_radius_primitive_recorded": True,
        "legacy_positive_radius_truncation_difference_recorded": True,
        "scoped_discrete_lower_limit_realization_improved": True,
        "formal_axis_based_inverse_verified": False,
        "formal_full_domain_moment_verified": False,
        "axis_regularity_verified": False,
        "bounded_source_to_axis_independently_certified": False,
        "global_radial_domain_materialized": False,
        "mean_recomputed_from_exact_upstream_in_public_api": True,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(
            signature.parameters
        ),
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


def build_mechanics_report() -> dict[str, object]:
    """Manufactured lower-limit regression; never candidate residual evidence."""

    radii = np.linspace(0.02, 0.80, 79)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(value) for value in radii),
        bump_center=0.50,
        bump_halfwidth=0.28,
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
    witness = _materialize_axis_completed_from_mean_witness(
        mean_witness,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    receipt = witness.to_receipt()
    receipt["mechanics_only"] = True
    receipt["candidate_residual_evidence"] = False
    receipt["formal_axis_inverse_certificate"] = False
    return receipt


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
