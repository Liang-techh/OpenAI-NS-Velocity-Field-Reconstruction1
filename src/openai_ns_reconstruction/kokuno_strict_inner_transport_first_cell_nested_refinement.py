"""Nested positive-radius refinement for the Kokuno first-cell radial inverse.

Agent-3 PR #875 replaces the first radial cell ``[0,r_min]`` by analytic
integration of a quadratic extrapolant built from the first three positive
radii.  Independent CR002 PR #879 then fixes the observability boundary: close
agreement of finite extrapolants that see only positive-radius samples cannot,
by itself, certify the real source on the unsampled interval ``[0,r_min]``.

This increment therefore does one narrower thing.  Starting from one already
registered strict-inner radial geometry, it evaluates the *same typed upstream
transport source* on a frozen nested positive-radius ladder

    base:   r_min, ...
    medium: r_min/2, r_min, ...
    fine:   r_min/4, r_min/2, r_min, ...

and applies the unchanged #875 polynomial-first-cell radial reconstruction at
each level.  Stresses are compared only on the common base radii.  The resulting
finite-grid refinement ratio is an engineering stability diagnostic, not a
formal first-cell error bound: every finite ladder still leaves an unsampled
interval next to the axis, exactly as CR002 #879 warns.

No caller may inject a residual, defect, source, mean, stress, pressure,
forcing, gain, damping, viscosity, derivative step, or scientific threshold.
The public path recomputes the exact typed strict-inner source at every level.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_strict_inner_transport_first_cell_polynomial_radial_stress import (
    PolynomialFirstCellStrictInnerTransportRadialStressWitness,
    _materialize_polynomial_first_cell_from_mean_witness,
    materialize_polynomial_first_cell_strict_inner_transport_radial_stress,
)
from .kokuno_strict_inner_transport_mean import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    REPOSITORY_VISCOSITY,
    ExactAgent1InnerTransportBackend,
    ExactAgent2StrictInnerTransportBackend,
    _materialize_from_transport_provider,
)
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)

TASK = "KOKUNO-A3-FIRST-CELL-NESTED-REFINEMENT-076"
SCHEMA = "kokuno-a3-first-cell-nested-refinement-v1"
PARENT_AGENT3_PR = 875
PARENT_AGENT3_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
PARENT_AGENT3_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"
INDEPENDENT_OBSERVABILITY_AUDIT_PR = 879
INDEPENDENT_OBSERVABILITY_AUDIT_HEAD = "1763282ab49c72a8998628805b53ed69f5e3890b"
INDEPENDENT_OBSERVABILITY_AUDIT_BLOB = "34e5eef092cc7b6f52341d25a9f8d92fdfc3de62"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

# Frozen before exact-head numerical output.  This is an engineering refinement
# diagnostic only, not the repository momentum/divergence acceptance gate.
REFINEMENT_RATIO_GATE = 1.5
ABSOLUTE_STABILITY_FLOOR = 1.0e-10


def _scalar_rms(values: np.ndarray) -> float:
    data = np.asarray(values, dtype=float)
    if data.size == 0 or not np.all(np.isfinite(data)):
        raise ValueError("RMS input must be nonempty and finite")
    return float(np.sqrt(np.mean(data * data)))


def _nested_geometries(
    geometry: StrictInnerTransportRadialGeometry,
) -> tuple[
    StrictInnerTransportRadialGeometry,
    StrictInnerTransportRadialGeometry,
    StrictInnerTransportRadialGeometry,
]:
    """Return frozen base / half-first-radius / quarter-first-radius grids."""

    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")
    radii = np.asarray(geometry.radii, dtype=float)
    if radii.ndim != 1 or radii.size < 9:
        raise ValueError("nested refinement requires at least nine base radii")
    r0 = float(radii[0])
    medium_radii = np.concatenate(([0.5 * r0], radii))
    fine_radii = np.concatenate(([0.25 * r0, 0.5 * r0], radii))

    def make(values: np.ndarray) -> StrictInnerTransportRadialGeometry:
        return StrictInnerTransportRadialGeometry(
            time=float(geometry.time),
            axial_z=float(geometry.axial_z),
            radii=tuple(float(value) for value in values),
            bump_center=float(geometry.bump_center),
            bump_halfwidth=float(geometry.bump_halfwidth),
        )

    return make(radii), make(medium_radii), make(fine_radii)


def _stress_refinement_metrics(
    coarse: dict[str, object],
    medium: dict[str, object],
    fine: dict[str, object],
) -> dict[str, object]:
    """Compare three #875 stresses on the common coarse positive radii."""

    c = np.asarray(coarse["stress"], dtype=float)
    m_all = np.asarray(medium["stress"], dtype=float)
    f_all = np.asarray(fine["stress"], dtype=float)
    if c.ndim != 1 or m_all.ndim != 1 or f_all.ndim != 1:
        raise ValueError("stress arrays must be one-dimensional")
    if m_all.size != c.size + 1 or f_all.size != c.size + 2:
        raise ValueError("nested stress arrays do not match the frozen ladder")
    if not (np.all(np.isfinite(c)) and np.all(np.isfinite(m_all)) and np.all(np.isfinite(f_all))):
        raise ValueError("nested stress arrays must be finite")

    m = m_all[1:]
    f = f_all[2:]
    coarse_medium = m - c
    medium_fine = f - m
    cm_rms = _scalar_rms(coarse_medium)
    mf_rms = _scalar_rms(medium_fine)
    cm_max = float(np.max(np.abs(coarse_medium)))
    mf_max = float(np.max(np.abs(medium_fine)))
    fine_rms = _scalar_rms(f)
    floor_reached = mf_rms <= ABSOLUTE_STABILITY_FLOOR
    raw_ratio = None if mf_rms == 0.0 else float(cm_rms / mf_rms)
    ratio_gate_passed = bool(
        floor_reached
        or (
            mf_rms <= cm_rms
            and raw_ratio is not None
            and raw_ratio >= REFINEMENT_RATIO_GATE
        )
    )

    return {
        "coarse_medium_common_radii_rms": cm_rms,
        "medium_fine_common_radii_rms": mf_rms,
        "coarse_medium_common_radii_max_abs": cm_max,
        "medium_fine_common_radii_max_abs": mf_max,
        "raw_refinement_ratio": raw_ratio,
        "fine_common_stress_rms": fine_rms,
        "medium_fine_relative_to_fine_rms": mf_rms / max(fine_rms, np.finfo(float).tiny),
        "absolute_stability_floor": ABSOLUTE_STABILITY_FLOOR,
        "absolute_stability_floor_reached": floor_reached,
        "refinement_ratio_gate": REFINEMENT_RATIO_GATE,
        "engineering_stability_gate_passed": ratio_gate_passed,
        "coarse_common_first_edge_stress": float(c[0]),
        "medium_common_first_edge_stress": float(m[0]),
        "fine_common_first_edge_stress": float(f[0]),
        "medium_new_subcell_edge_stress": float(m_all[0]),
        "fine_new_subcell_edge_stress": float(f_all[0]),
        "finite_positive_radius_refinement_is_formal_error_bound": False,
    }


def _common_mean_replay_max_abs(
    coarse: PolynomialFirstCellStrictInnerTransportRadialStressWitness,
    medium: PolynomialFirstCellStrictInnerTransportRadialStressWitness,
    fine: PolynomialFirstCellStrictInnerTransportRadialStressWitness,
) -> float:
    c = np.asarray(coarse.mean_witness.mean_transport_cylindrical, dtype=float)
    m = np.asarray(medium.mean_witness.mean_transport_cylindrical, dtype=float)
    f = np.asarray(fine.mean_witness.mean_transport_cylindrical, dtype=float)
    if c.ndim != 2 or c.shape[1] != 3:
        raise ValueError("coarse mean transport must have shape (N,3)")
    if m.shape != (c.shape[0] + 1, 3) or f.shape != (c.shape[0] + 2, 3):
        raise ValueError("nested mean witnesses do not match frozen radial ladder")
    return float(
        max(
            np.max(np.abs(c - m[1:])),
            np.max(np.abs(c - f[2:])),
            np.max(np.abs(m[1:] - f[2:])),
        )
    )


@dataclass(frozen=True)
class NestedFirstCellRefinementWitness:
    base_geometry: StrictInnerTransportRadialGeometry
    coarse: PolynomialFirstCellStrictInnerTransportRadialStressWitness
    medium: PolynomialFirstCellStrictInnerTransportRadialStressWitness
    fine: PolynomialFirstCellStrictInnerTransportRadialStressWitness
    theta_e2_refinement: dict[str, object]
    axial_e1_refinement: dict[str, object]
    common_mean_replay_max_abs: float
    backend_kind: str

    @property
    def engineering_stability_passed(self) -> bool:
        return bool(
            self.theta_e2_refinement["engineering_stability_gate_passed"]
            and self.axial_e1_refinement["engineering_stability_gate_passed"]
        )

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "base_geometry": asdict(self.base_geometry),
            "backend_kind": self.backend_kind,
            "nested_first_radii": {
                "coarse": float(self.coarse.geometry.radii[0]),
                "medium": float(self.medium.geometry.radii[0]),
                "fine": float(self.fine.geometry.radii[0]),
            },
            "common_mean_replay_max_abs": float(self.common_mean_replay_max_abs),
            "theta_e2_refinement": dict(self.theta_e2_refinement),
            "axial_e1_refinement": dict(self.axial_e1_refinement),
            "engineering_stability_passed": self.engineering_stability_passed,
            "coarse": self.coarse.to_receipt(),
            "medium": self.medium.to_receipt(),
            "fine": self.fine.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _assemble_nested_witness(
    geometry: StrictInnerTransportRadialGeometry,
    coarse: PolynomialFirstCellStrictInnerTransportRadialStressWitness,
    medium: PolynomialFirstCellStrictInnerTransportRadialStressWitness,
    fine: PolynomialFirstCellStrictInnerTransportRadialStressWitness,
    *,
    backend_kind: str,
) -> NestedFirstCellRefinementWitness:
    theta = _stress_refinement_metrics(
        coarse.theta_polynomial_first_cell_stress,
        medium.theta_polynomial_first_cell_stress,
        fine.theta_polynomial_first_cell_stress,
    )
    axial = _stress_refinement_metrics(
        coarse.axial_polynomial_first_cell_stress,
        medium.axial_polynomial_first_cell_stress,
        fine.axial_polynomial_first_cell_stress,
    )
    replay = _common_mean_replay_max_abs(coarse, medium, fine)
    return NestedFirstCellRefinementWitness(
        base_geometry=geometry,
        coarse=coarse,
        medium=medium,
        fine=fine,
        theta_e2_refinement=theta,
        axial_e1_refinement=axial,
        common_mean_replay_max_abs=replay,
        backend_kind=backend_kind,
    )


def materialize_nested_first_cell_refinement(
    transport_backend: ExactAgent2StrictInnerTransportBackend,
    inner_backend: ExactAgent1InnerTransportBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> NestedFirstCellRefinementWitness:
    """Recompute exact typed source on the frozen base/half/quarter radial ladder."""

    if not isinstance(transport_backend, ExactAgent2StrictInnerTransportBackend):
        raise TypeError("transport_backend must be ExactAgent2StrictInnerTransportBackend")
    if not isinstance(inner_backend, ExactAgent1InnerTransportBackend):
        raise TypeError("inner_backend must be ExactAgent1InnerTransportBackend")
    coarse_geometry, medium_geometry, fine_geometry = _nested_geometries(geometry)
    coarse = materialize_polynomial_first_cell_strict_inner_transport_radial_stress(
        transport_backend, inner_backend, coarse_geometry
    )
    medium = materialize_polynomial_first_cell_strict_inner_transport_radial_stress(
        transport_backend, inner_backend, medium_geometry
    )
    fine = materialize_polynomial_first_cell_strict_inner_transport_radial_stress(
        transport_backend, inner_backend, fine_geometry
    )
    return _assemble_nested_witness(
        geometry,
        coarse,
        medium,
        fine,
        backend_kind="exact-typed-source-positive-radius-nested-refinement",
    )


def _materialize_nested_from_transport_provider(
    provider: Callable[..., tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]],
    geometry: StrictInnerTransportRadialGeometry,
    *,
    backend_kind: str,
) -> NestedFirstCellRefinementWitness:
    """Internal mechanics/CI path; public scientific API remains backend-typed."""

    levels: list[PolynomialFirstCellStrictInnerTransportRadialStressWitness] = []
    for level_name, level_geometry in zip(
        ("coarse", "medium", "fine"), _nested_geometries(geometry)
    ):
        radii = np.asarray(level_geometry.radii, dtype=float)
        mean = _materialize_from_transport_provider(
            provider,
            radii,
            np.full_like(radii, float(level_geometry.axial_z)),
            np.full_like(radii, float(level_geometry.time)),
            agent1_backend=None,
            agent2_backend=None,
            backend_kind=f"{backend_kind}:{level_name}",
        )
        levels.append(
            _materialize_polynomial_first_cell_from_mean_witness(
                mean,
                level_geometry,
                backend_kind=f"{backend_kind}:{level_name}:polynomial-first-cell",
            )
        )
    return _assemble_nested_witness(
        geometry,
        levels[0],
        levels[1],
        levels[2],
        backend_kind=backend_kind,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_nested_first_cell_refinement)
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
        "refinement_ratio_gate",
        "stability_floor",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "independent_observability_audit_pr": INDEPENDENT_OBSERVABILITY_AUDIT_PR,
        "independent_observability_audit_head": INDEPENDENT_OBSERVABILITY_AUDIT_HEAD,
        "independent_observability_audit_blob": INDEPENDENT_OBSERVABILITY_AUDIT_BLOB,
        "nested_first_radius_factors": [1.0, 0.5, 0.25],
        "refinement_ratio_gate": REFINEMENT_RATIO_GATE,
        "absolute_stability_floor": ABSOLUTE_STABILITY_FLOOR,
        "engineering_gate_frozen_before_exact_head_output": True,
        "angular_orders": list(ANGULAR_ORDERS),
        "repository_viscosity": REPOSITORY_VISCOSITY,
        "public_api_recomputes_exact_typed_source": True,
        "caller_supplied_source_allowed": False,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "real_source_positive_radius_nested_refinement_assessed": True,
        "finite_positive_radius_refinement_stability_is_formal_error_bound": False,
        "finite_nested_refinement_closes_hidden_cell_observability": False,
        "real_source_first_cell_convergence_verified": False,
        "formal_first_cell_accuracy_verified": False,
        "formal_axis_based_inverse_verified": False,
        "formal_full_domain_moment_verified": False,
        "axis_regularity_verified": False,
        "bounded_source_to_axis_independently_certified": False,
        "global_radial_domain_materialized": False,
        "complete_ns_defect": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "global_corrected_leading_join_materialized": False,
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


def _mechanics_provider():
    """Smooth cubic radial mean plus m=2 modes; never candidate evidence."""

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
        mr = 0.03 * radius
        mt = 0.20 + 0.40 * radius - 0.30 * radius**2 + 0.15 * radius**3
        mz = 0.30 + 0.20 * radius + 0.08 * radius**3
        mr = mr + 0.02 * np.cos(2.0 * theta)
        mt = mt - 0.015 * np.sin(2.0 * theta)
        mz = mz + 0.01 * np.cos(2.0 * theta)
        transport = np.stack((mr * c - mt * s, mr * s + mt * c, mz), axis=-1)
        zero = np.zeros_like(transport)
        return transport, zero, zero, transport, REPOSITORY_VISCOSITY

    return provider


def build_mechanics_report() -> dict[str, object]:
    radii = np.linspace(0.20, 1.00, 17)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(value) for value in radii),
        bump_center=0.65,
        bump_halfwidth=0.18,
    )
    witness = _materialize_nested_from_transport_provider(
        _mechanics_provider(),
        geometry,
        backend_kind="manufactured-smooth-cubic-mechanics-only",
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "formal_first_cell_certificate": False,
        "witness": witness.to_receipt(),
        "truth_boundary": truth_boundary(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = build_mechanics_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
