"""Polynomially reconstruct the Kokuno radial inverse first cell.

PR #868 repaired the missing ``0 -> r_min`` interval by prepending a virtual
weighted endpoint and applying a trapezoid.  Independent CR002 PR #871 then
showed that this is still only a discrete convention: for a constant unweighted
source and the e=2 channel, the first-cell primitive is 50% too large.

This Agent-3 increment keeps the same exact typed strict-inner transport mean,
compact moment complement, and positive-radius grid.  It changes only the first
cell of the weighted radial quadrature.  The unweighted source ``G(r)`` is
reconstructed from the first three positive nodes by a quadratic Lagrange
polynomial, and ``r^e G(r)`` is integrated analytically from 0 to ``r_min``.
The same rule is used for the source moment, compact-bump normalization, and
moment-complement primitive, so outer-edge closure is preserved by construction.
A two-node affine reconstruction is recorded independently as a local
sensitivity diagnostic.

The rule is exact for any quadratic unweighted source on the first cell.  That
mechanics fact directly removes the constant-source e=2 defect identified by
CR002 #871.  It is still an extrapolatory repository numerical realization: it
does not certify the real typed source's formal first-cell error/convergence,
full radial-domain moment semantics, global axis regularity, a complete NS
defect, or a finite correction-cycle residual reduction.
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
from numpy.polynomial import Polynomial

from .kokuno_actual_oscillatory_mean_stress import _compact_cos8_bump, _scalar_rms
from .kokuno_strict_inner_transport_axis_completed_radial_stress import (
    AxisCompletedStrictInnerTransportRadialStressWitness,
    _axis_completed_compact_radial_stress,
    _materialize_axis_completed_from_mean_witness,
    materialize_axis_completed_strict_inner_transport_radial_stress,
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
    _mechanics_transport_provider,
)

TASK = "KOKUNO-A3-FIRST-CELL-POLYNOMIAL-QUADRATURE-075"
SCHEMA = "kokuno-a3-first-cell-polynomial-radial-stress-v1"
PARENT_AGENT3_PR = 868
PARENT_AGENT3_HEAD = "c8268100a36942de9fa966d76f22ad6148ceaa02"
PARENT_AGENT3_SOURCE_BLOB = "041d4298b5350b25a5699492ddf39c35fae4aaab"
INDEPENDENT_SCOPE_AUDIT_PR = 871
INDEPENDENT_SCOPE_AUDIT_HEAD = "594af2273c2c27a7a64a9a4340da43472343b83a"
INDEPENDENT_SCOPE_AUDIT_BLOB = "dc92c3537c286af8fd2f801584d43a6e6c58f782"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"


def _validate_radial_inputs(
    values: np.ndarray, radii: np.ndarray, *, minimum_nodes: int = 3
) -> tuple[np.ndarray, np.ndarray]:
    f = np.asarray(values, dtype=float)
    r = np.asarray(radii, dtype=float)
    if f.ndim != 1 or r.ndim != 1 or f.shape != r.shape:
        raise ValueError("first-cell quadrature requires matching 1D arrays")
    if r.size < minimum_nodes:
        raise ValueError(f"first-cell quadrature requires at least {minimum_nodes} nodes")
    if np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("radii must be positive and strictly increasing")
    if not (np.all(np.isfinite(f)) and np.all(np.isfinite(r))):
        raise ValueError("first-cell quadrature inputs must be finite")
    return f, r


def _first_cell_lagrange_weights(
    radii: np.ndarray, *, exponent: int, degree: int
) -> np.ndarray:
    """Return weights integrating r^e times the local Lagrange extrapolant.

    The interpolation nodes are the first ``degree+1`` strictly-positive radii,
    while the integration interval is ``[0, r[0]]``.  Degree 1 and 2 are used in
    this module.  The operation is deliberately recorded as extrapolation, not
    as a source theorem or an axis-regularity proof.
    """

    r = np.asarray(radii, dtype=float)
    if exponent not in (1, 2):
        raise ValueError("only e=1 and e=2 channels are supported")
    if degree not in (1, 2):
        raise ValueError("first-cell polynomial degree must be 1 or 2")
    if r.ndim != 1 or r.size < degree + 1:
        raise ValueError("insufficient positive-radius nodes for interpolation")
    nodes = r[: degree + 1]
    if np.any(nodes <= 0.0) or np.any(np.diff(nodes) <= 0.0):
        raise ValueError("interpolation nodes must be positive and increasing")

    monomial = Polynomial([0.0] * exponent + [1.0])
    weights: list[float] = []
    upper = float(nodes[0])
    for j, node in enumerate(nodes):
        basis = Polynomial([1.0])
        denominator = 1.0
        for k, other in enumerate(nodes):
            if k == j:
                continue
            basis = basis * Polynomial([-float(other), 1.0])
            denominator *= float(node - other)
        basis = basis / denominator
        primitive = (basis * monomial).integ()
        weights.append(float(primitive(upper) - primitive(0.0)))
    out = np.asarray(weights, dtype=float)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("non-finite first-cell interpolation weights")
    return out


def _first_cell_polynomial_integral(
    values: np.ndarray,
    radii: np.ndarray,
    *,
    exponent: int,
    degree: int,
) -> float:
    f, r = _validate_radial_inputs(values, radii, minimum_nodes=degree + 1)
    weights = _first_cell_lagrange_weights(r, exponent=exponent, degree=degree)
    return float(np.dot(weights, f[: degree + 1]))


def _cumulative_with_polynomial_first_cell(
    values: np.ndarray,
    radii: np.ndarray,
    *,
    exponent: int,
    degree: int = 2,
) -> np.ndarray:
    """Integrate r^e G(r), using exact polynomial integration on first cell.

    Only the ``0 -> r_min`` cell changes relative to #868.  Every positive-node
    interval retains the inherited weighted trapezoid, making the increment
    minimal and isolating the CR002 #871 defect.
    """

    f, r = _validate_radial_inputs(values, radii, minimum_nodes=max(3, degree + 1))
    first = _first_cell_polynomial_integral(
        f, r, exponent=exponent, degree=degree
    )
    weighted = (r**exponent) * f
    increments = 0.5 * (weighted[1:] + weighted[:-1]) * np.diff(r)
    cumulative = np.empty_like(r)
    cumulative[0] = first
    if r.size > 1:
        cumulative[1:] = first + np.cumsum(increments)
    return cumulative


def _relative_difference(a: float, b: float) -> float:
    scale = max(abs(float(a)), abs(float(b)), np.finfo(float).tiny)
    return abs(float(a) - float(b)) / scale


def _polynomial_first_cell_compact_radial_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    bump_center: float,
    bump_halfwidth: float,
) -> dict[str, Any]:
    """Apply the compact inverse with a quadratic-extrapolated first cell."""

    f, r = _validate_radial_inputs(source, radii, minimum_nodes=9)
    if exponent not in (1, 2):
        raise ValueError("only the axial e=1 and tangential e=2 channels are supported")
    center = float(bump_center)
    halfwidth = float(bump_halfwidth)
    if not (math.isfinite(center) and math.isfinite(halfwidth) and halfwidth > 0.0):
        raise ValueError("compact bump geometry must be finite with positive halfwidth")
    support_min = center - halfwidth
    support_max = center + halfwidth
    tolerance = 64.0 * math.ulp(max(1.0, abs(center), abs(halfwidth)))
    if support_min < r[0] - tolerance or support_max > r[-1] + tolerance:
        raise ValueError("compact bump support must remain inside positive sampled radii")

    source_quadratic_primitive = _cumulative_with_polynomial_first_cell(
        f, r, exponent=exponent, degree=2
    )
    source_linear_first = _first_cell_polynomial_integral(
        f, r, exponent=exponent, degree=1
    )
    source_quadratic_first = float(source_quadratic_primitive[0])
    moment = float(source_quadratic_primitive[-1])

    bump_raw = _compact_cos8_bump(r, center, halfwidth)
    bump_quadratic_primitive = _cumulative_with_polynomial_first_cell(
        bump_raw, r, exponent=exponent, degree=2
    )
    bump_integral = float(bump_quadratic_primitive[-1])
    if not math.isfinite(bump_integral) or bump_integral <= 0.0:
        raise RuntimeError("polynomial-first-cell compact bump lost positive normalization")
    bump = bump_raw / bump_integral

    complement = f - bump * moment
    complement_primitive = _cumulative_with_polynomial_first_cell(
        complement, r, exponent=exponent, degree=2
    )
    complement_linear_first = _first_cell_polynomial_integral(
        complement, r, exponent=exponent, degree=1
    )
    complement_quadratic_first = float(complement_primitive[0])
    complement_moment = float(complement_primitive[-1])
    weight = r**exponent
    stress = -complement_primitive / weight
    reconstructed_force = -f + bump * moment

    parent = _axis_completed_compact_radial_stress(
        r,
        f,
        exponent=exponent,
        bump_center=center,
        bump_halfwidth=halfwidth,
    )
    parent_stress = np.asarray(parent["stress"], dtype=float)
    stress_delta = stress - parent_stress

    return {
        "exponent": int(exponent),
        "quadrature_lower_limit": 0.0,
        "first_positive_radius": float(r[0]),
        "first_cell_rule": "quadratic_lagrange_extrapolation_of_unweighted_source",
        "first_cell_polynomial_degree": 2,
        "first_cell_quadratic_weights": _first_cell_lagrange_weights(
            r, exponent=exponent, degree=2
        ),
        "first_cell_linear_weights": _first_cell_lagrange_weights(
            r, exponent=exponent, degree=1
        ),
        "source_first_cell_linear_estimate": source_linear_first,
        "source_first_cell_quadratic_estimate": source_quadratic_first,
        "source_first_cell_linear_quadratic_relative_difference": _relative_difference(
            source_linear_first, source_quadratic_first
        ),
        "weighted_moment": moment,
        "bump_weighted_integral": float(
            _cumulative_with_polynomial_first_cell(
                bump, r, exponent=exponent, degree=2
            )[-1]
        ),
        "moment_complement_weighted_moment": complement_moment,
        "axis_to_first_radius_weighted_primitive": complement_quadratic_first,
        "complement_first_cell_linear_estimate": complement_linear_first,
        "complement_first_cell_quadratic_estimate": complement_quadratic_first,
        "complement_first_cell_linear_quadratic_relative_difference": _relative_difference(
            complement_linear_first, complement_quadratic_first
        ),
        "stress": stress,
        "stress_rms": _scalar_rms(stress),
        "stress_max_abs": float(np.max(np.abs(stress))),
        "stress_inner_positive_edge": float(stress[0]),
        "stress_outer_edge": float(stress[-1]),
        "reconstructed_force_rms": _scalar_rms(reconstructed_force),
        "parent_868_stress": parent_stress,
        "stress_delta_vs_parent_868": stress_delta,
        "stress_delta_vs_parent_868_rms": _scalar_rms(stress_delta),
        "stress_delta_vs_parent_868_max_abs": float(np.max(np.abs(stress_delta))),
        "parent_868_first_cell_primitive": float(
            parent["axis_to_first_radius_weighted_primitive"]
        ),
        "parent_868_stress_inner_positive_edge": float(
            parent["stress_inner_positive_edge"]
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
class PolynomialFirstCellStrictInnerTransportRadialStressWitness:
    geometry: StrictInnerTransportRadialGeometry
    parent_868_witness: AxisCompletedStrictInnerTransportRadialStressWitness
    theta_polynomial_first_cell_stress: dict[str, object]
    axial_polynomial_first_cell_stress: dict[str, object]
    backend_kind: str

    @property
    def mean_witness(self) -> StrictInnerTransportMeanWitness:
        return self.parent_868_witness.mean_witness

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "geometry": asdict(self.geometry),
            "backend_kind": self.backend_kind,
            "mean_witness": self.mean_witness.to_receipt(),
            "radial_mean_profile": np.asarray(
                self.parent_868_witness.legacy_witness.radial_mean_profile, dtype=float
            ).tolist(),
            "theta_mean_profile": np.asarray(
                self.parent_868_witness.legacy_witness.theta_mean_profile, dtype=float
            ).tolist(),
            "axial_mean_profile": np.asarray(
                self.parent_868_witness.legacy_witness.axial_mean_profile, dtype=float
            ).tolist(),
            "theta_e2_polynomial_first_cell": _clean_report(
                self.theta_polynomial_first_cell_stress
            ),
            "axial_e1_polynomial_first_cell": _clean_report(
                self.axial_polynomial_first_cell_stress
            ),
            "parent_868_axis_completed": {
                "theta_e2": _clean_report(
                    self.parent_868_witness.theta_axis_completed_stress
                ),
                "axial_e1": _clean_report(
                    self.parent_868_witness.axial_axis_completed_stress
                ),
            },
            "truth_boundary": truth_boundary(),
        }


def _materialize_polynomial_first_cell_from_mean_witness(
    witness: StrictInnerTransportMeanWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    backend_kind: str,
) -> PolynomialFirstCellStrictInnerTransportRadialStressWitness:
    parent = _materialize_axis_completed_from_mean_witness(
        witness,
        geometry,
        backend_kind=backend_kind + ":parent-868-axis-completed",
    )
    radii = np.asarray(geometry.radii, dtype=float)
    theta = np.asarray(parent.legacy_witness.theta_mean_profile, dtype=float)
    axial = np.asarray(parent.legacy_witness.axial_mean_profile, dtype=float)
    theta_stress = _polynomial_first_cell_compact_radial_stress(
        radii,
        theta,
        exponent=2,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    axial_stress = _polynomial_first_cell_compact_radial_stress(
        radii,
        axial,
        exponent=1,
        bump_center=geometry.bump_center,
        bump_halfwidth=geometry.bump_halfwidth,
    )
    return PolynomialFirstCellStrictInnerTransportRadialStressWitness(
        geometry=geometry,
        parent_868_witness=parent,
        theta_polynomial_first_cell_stress=theta_stress,
        axial_polynomial_first_cell_stress=axial_stress,
        backend_kind=backend_kind,
    )


def materialize_polynomial_first_cell_strict_inner_transport_radial_stress(
    transport_backend: ExactAgent2StrictInnerTransportBackend,
    inner_backend: ExactAgent1InnerTransportBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> PolynomialFirstCellStrictInnerTransportRadialStressWitness:
    """Recompute exact typed source through #868, then replace only first cell."""

    if not isinstance(transport_backend, ExactAgent2StrictInnerTransportBackend):
        raise TypeError("transport_backend must be ExactAgent2StrictInnerTransportBackend")
    if not isinstance(inner_backend, ExactAgent1InnerTransportBackend):
        raise TypeError("inner_backend must be ExactAgent1InnerTransportBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    parent = materialize_axis_completed_strict_inner_transport_radial_stress(
        transport_backend,
        inner_backend,
        geometry,
    )
    radii = np.asarray(geometry.radii, dtype=float)
    theta = np.asarray(parent.legacy_witness.theta_mean_profile, dtype=float)
    axial = np.asarray(parent.legacy_witness.axial_mean_profile, dtype=float)
    return PolynomialFirstCellStrictInnerTransportRadialStressWitness(
        geometry=geometry,
        parent_868_witness=parent,
        theta_polynomial_first_cell_stress=_polynomial_first_cell_compact_radial_stress(
            radii,
            theta,
            exponent=2,
            bump_center=geometry.bump_center,
            bump_halfwidth=geometry.bump_halfwidth,
        ),
        axial_polynomial_first_cell_stress=_polynomial_first_cell_compact_radial_stress(
            radii,
            axial,
            exponent=1,
            bump_center=geometry.bump_center,
            bump_halfwidth=geometry.bump_halfwidth,
        ),
        backend_kind="exact-a3-868-typed-source-polynomial-first-cell",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_polynomial_first_cell_strict_inner_transport_radial_stress
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
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "independent_scope_audit_pr": INDEPENDENT_SCOPE_AUDIT_PR,
        "independent_scope_audit_head": INDEPENDENT_SCOPE_AUDIT_HEAD,
        "independent_scope_audit_blob": INDEPENDENT_SCOPE_AUDIT_BLOB,
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "angular_orders": list(ANGULAR_ORDERS),
        "first_cell_rule": "quadratic_lagrange_extrapolation_of_unweighted_source",
        "first_cell_degree": 2,
        "first_cell_affine_diagnostic_degree": 1,
        "constant_source_cr002_e2_mechanics_defect_removed": True,
        "quadratic_unweighted_first_cell_mechanics_exact": True,
        "real_source_first_cell_linear_quadratic_disagreement_recorded": True,
        "real_source_first_cell_convergence_verified": False,
        "formal_first_cell_accuracy_verified": False,
        "formal_axis_based_inverse_verified": False,
        "formal_full_domain_moment_verified": False,
        "axis_regularity_verified": False,
        "bounded_source_to_axis_independently_certified": False,
        "global_radial_domain_materialized": False,
        "mean_recomputed_from_exact_upstream_in_public_api": True,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
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


def _formal_polynomial_first_cell(
    coefficients: tuple[float, ...], r_min: float, exponent: int
) -> float:
    return float(
        sum(
            coefficient * r_min ** (exponent + power + 1) / (exponent + power + 1)
            for power, coefficient in enumerate(coefficients)
        )
    )


def build_mechanics_report() -> dict[str, object]:
    """Mechanics-only exactness and typed-routing receipt."""

    radii = np.linspace(0.02, 0.80, 79)
    coefficients = (1.2, -0.7, 0.4)
    source = coefficients[0] + coefficients[1] * radii + coefficients[2] * radii**2
    polynomial_checks: dict[str, object] = {}
    for exponent in (1, 2):
        quadratic = _first_cell_polynomial_integral(
            source, radii, exponent=exponent, degree=2
        )
        linear = _first_cell_polynomial_integral(
            source, radii, exponent=exponent, degree=1
        )
        formal = _formal_polynomial_first_cell(coefficients, float(radii[0]), exponent)
        constant = np.ones_like(radii)
        constant_estimate = _first_cell_polynomial_integral(
            constant, radii, exponent=exponent, degree=2
        )
        constant_formal = float(radii[0] ** (exponent + 1) / (exponent + 1))
        polynomial_checks[f"e{exponent}"] = {
            "formal_quadratic_source_first_cell": formal,
            "quadratic_rule_first_cell": quadratic,
            "quadratic_rule_absolute_error": abs(quadratic - formal),
            "linear_rule_first_cell": linear,
            "linear_quadratic_relative_difference": _relative_difference(linear, quadratic),
            "constant_source_formal_first_cell": constant_formal,
            "constant_source_quadratic_rule_first_cell": constant_estimate,
            "constant_source_absolute_error": abs(constant_estimate - constant_formal),
        }

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
    witness = _materialize_polynomial_first_cell_from_mean_witness(
        mean_witness,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "formal_axis_inverse_certificate": False,
        "polynomial_first_cell_checks": polynomial_checks,
        "typed_routing_receipt": witness.to_receipt(),
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
