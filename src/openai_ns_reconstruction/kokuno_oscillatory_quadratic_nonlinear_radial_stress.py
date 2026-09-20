"""Route the current pure oscillatory quadratic mean through compact radial stress.

Agent 2 PR #908 exposes the current checksum-bound strict-inner nonlinear pieces
without performing cylindrical averaging.  Agent 3 PR #915 projects those pieces
and isolates

    Q_bar = <(u_osc . grad) u_osc>_theta,
    M_bar = <(u_inner . grad) u_osc + (u_osc . grad) u_inner>_theta,
    A_bar = M_bar + Q_bar.

This module performs one next Agent-3-owned step: apply the already-audited
polynomial-first-cell compact radial inverse from PR #875 to the tangential and
axial channels of Q_bar, M_bar and A_bar.  Because the same radial operator is
linear, it also checks the stress-level attribution

    sigma(A_bar) = sigma(M_bar) + sigma(Q_bar).

The radial cylindrical component is recorded but deliberately not inverted by
this two-channel operator.  This remains strict-inner nonlinear transport only:
there is no pressure, forcing, corrected/global leading join, complete NS defect,
correction velocity, finite correction-cycle run, or held-out residual claim.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import inspect
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_oscillatory_nonlinear_mean_attribution import (
    ANGULAR_ORDERS,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    ExactAgent2OscillatoryNonlinearHandoff,
    OscillatoryNonlinearMeanAttributionWitness,
    _analytic_provider,
    _materialize_from_provider,
    materialize_oscillatory_nonlinear_mean_attribution,
)
from .kokuno_oscillatory_transport_radial_stress import (
    _clean_report,
    _relative_max_closure,
    _stress_array,
)
from .kokuno_strict_inner_transport_first_cell_polynomial_radial_stress import (
    _polynomial_first_cell_compact_radial_stress,
)
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)

TASK = "KOKUNO-A3-OSCILLATORY-QUADRATIC-NONLINEAR-RADIAL-STRESS-087"
SCHEMA = "kokuno-a3-oscillatory-quadratic-nonlinear-radial-stress-v1"
PARENT_AGENT3_PR = 915
PARENT_AGENT3_HEAD = "278b98dbda2f8de8b21a86f7e7e4740939387dd2"
PARENT_AGENT3_SOURCE_BLOB = "ff3bb9d29ecab1c8469dc2599ce67dcc9cc13805"
RADIAL_OPERATOR_AGENT3_PR = 875
RADIAL_OPERATOR_AGENT3_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"

MEAN_PIECE_CLOSURE_ABSOLUTE_GATE = 1.0e-12
STRESS_PIECE_CLOSURE_RELATIVE_GATE = 2.0e-11


@dataclass(frozen=True)
class OscillatoryQuadraticNonlinearRadialStressWitness:
    """Pure-quadratic/mixed/aggregate nonlinear means and compact stresses."""

    geometry: StrictInnerTransportRadialGeometry
    mean_witness: OscillatoryNonlinearMeanAttributionWitness
    radial_quadratic_mean_profile: np.ndarray
    radial_mixed_mean_profile: np.ndarray
    radial_aggregate_mean_profile: np.ndarray
    theta_quadratic_stress: dict[str, object]
    theta_mixed_stress: dict[str, object]
    theta_aggregate_stress: dict[str, object]
    axial_quadratic_stress: dict[str, object]
    axial_mixed_stress: dict[str, object]
    axial_aggregate_stress: dict[str, object]
    mean_piece_closure_absolute_max: float
    theta_stress_piece_closure_relative_max: float
    axial_stress_piece_closure_relative_max: float
    backend_kind: str

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "geometry": asdict(self.geometry),
            "backend_kind": self.backend_kind,
            "mean_witness": self.mean_witness.to_receipt(),
            "radial_quadratic_mean_profile": np.asarray(
                self.radial_quadratic_mean_profile, dtype=float
            ).tolist(),
            "radial_mixed_mean_profile": np.asarray(
                self.radial_mixed_mean_profile, dtype=float
            ).tolist(),
            "radial_aggregate_mean_profile": np.asarray(
                self.radial_aggregate_mean_profile, dtype=float
            ).tolist(),
            "theta_e2": {
                "quadratic": _clean_report(self.theta_quadratic_stress),
                "mixed": _clean_report(self.theta_mixed_stress),
                "aggregate": _clean_report(self.theta_aggregate_stress),
                "stress_piece_closure_relative_max": self.theta_stress_piece_closure_relative_max,
            },
            "axial_e1": {
                "quadratic": _clean_report(self.axial_quadratic_stress),
                "mixed": _clean_report(self.axial_mixed_stress),
                "aggregate": _clean_report(self.axial_aggregate_stress),
                "stress_piece_closure_relative_max": self.axial_stress_piece_closure_relative_max,
            },
            "mean_piece_closure_absolute_max": self.mean_piece_closure_absolute_max,
            "truth_boundary": truth_boundary(),
        }


def _apply_piece_channel(
    radii: np.ndarray,
    quadratic: np.ndarray,
    mixed: np.ndarray,
    aggregate: np.ndarray,
    *,
    exponent: int,
    geometry: StrictInnerTransportRadialGeometry,
) -> tuple[dict[str, object], dict[str, object], dict[str, object], float]:
    kwargs = {
        "exponent": exponent,
        "bump_center": geometry.bump_center,
        "bump_halfwidth": geometry.bump_halfwidth,
    }
    quadratic_report = _polynomial_first_cell_compact_radial_stress(
        radii, quadratic, **kwargs
    )
    mixed_report = _polynomial_first_cell_compact_radial_stress(
        radii, mixed, **kwargs
    )
    aggregate_report = _polynomial_first_cell_compact_radial_stress(
        radii, aggregate, **kwargs
    )

    sigma_quadratic = _stress_array(quadratic_report, "quadratic")
    sigma_mixed = _stress_array(mixed_report, "mixed")
    sigma_aggregate = _stress_array(aggregate_report, "aggregate")
    closure = _relative_max_closure(
        sigma_aggregate, sigma_mixed + sigma_quadratic
    )
    return quadratic_report, mixed_report, aggregate_report, closure


def _materialize_from_mean_witness(
    witness: OscillatoryNonlinearMeanAttributionWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    backend_kind: str,
) -> OscillatoryQuadraticNonlinearRadialStressWitness:
    if not isinstance(witness, OscillatoryNonlinearMeanAttributionWitness):
        raise TypeError("witness must be OscillatoryNonlinearMeanAttributionWitness")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    radii = np.asarray(geometry.radii, dtype=float)
    expected = (radii.size, 3)
    arrays = {
        "quadratic": np.asarray(
            witness.mean_oscillatory_self_advection_cylindrical, dtype=float
        ),
        "mixed": np.asarray(witness.mean_mixed_cross_cylindrical, dtype=float),
        "aggregate": np.asarray(
            witness.mean_aggregate_nonlinear_cylindrical, dtype=float
        ),
    }
    for label, value in arrays.items():
        if value.shape != expected:
            raise ValueError(
                f"{label} nonlinear mean shape {value.shape} does not match {expected}"
            )
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{label} nonlinear mean contains non-finite values")

    witness_r = np.asarray(witness.radius, dtype=float)
    witness_z = np.asarray(witness.z, dtype=float)
    witness_t = np.asarray(witness.t, dtype=float)
    if witness_r.shape != radii.shape or not np.array_equal(witness_r, radii):
        raise ValueError("nonlinear mean witness radial grid does not match radial geometry")
    expected_z = np.full_like(radii, float(geometry.axial_z))
    expected_t = np.full_like(radii, float(geometry.time))
    if witness_z.shape != radii.shape or not np.array_equal(witness_z, expected_z):
        raise ValueError("nonlinear mean witness z values do not match radial geometry")
    if witness_t.shape != radii.shape or not np.array_equal(witness_t, expected_t):
        raise ValueError("nonlinear mean witness time values do not match radial geometry")
    if tuple(witness.angular_orders) != tuple(ANGULAR_ORDERS):
        raise ValueError("nonlinear mean witness angular ladder drifted")

    mean_closure = float(
        np.max(np.abs(arrays["aggregate"] - (arrays["mixed"] + arrays["quadratic"])))
    )
    if mean_closure > MEAN_PIECE_CLOSURE_ABSOLUTE_GATE:
        raise ValueError("quadratic/mixed nonlinear mean attribution closure failed")

    theta_quadratic, theta_mixed, theta_aggregate, theta_closure = _apply_piece_channel(
        radii,
        arrays["quadratic"][:, 1],
        arrays["mixed"][:, 1],
        arrays["aggregate"][:, 1],
        exponent=2,
        geometry=geometry,
    )
    axial_quadratic, axial_mixed, axial_aggregate, axial_closure = _apply_piece_channel(
        radii,
        arrays["quadratic"][:, 2],
        arrays["mixed"][:, 2],
        arrays["aggregate"][:, 2],
        exponent=1,
        geometry=geometry,
    )

    if theta_closure > STRESS_PIECE_CLOSURE_RELATIVE_GATE:
        raise ValueError("theta radial-stress nonlinear piece attribution closure failed")
    if axial_closure > STRESS_PIECE_CLOSURE_RELATIVE_GATE:
        raise ValueError("axial radial-stress nonlinear piece attribution closure failed")

    return OscillatoryQuadraticNonlinearRadialStressWitness(
        geometry=geometry,
        mean_witness=witness,
        radial_quadratic_mean_profile=np.array(arrays["quadratic"][:, 0], copy=True),
        radial_mixed_mean_profile=np.array(arrays["mixed"][:, 0], copy=True),
        radial_aggregate_mean_profile=np.array(arrays["aggregate"][:, 0], copy=True),
        theta_quadratic_stress=theta_quadratic,
        theta_mixed_stress=theta_mixed,
        theta_aggregate_stress=theta_aggregate,
        axial_quadratic_stress=axial_quadratic,
        axial_mixed_stress=axial_mixed,
        axial_aggregate_stress=axial_aggregate,
        mean_piece_closure_absolute_max=mean_closure,
        theta_stress_piece_closure_relative_max=theta_closure,
        axial_stress_piece_closure_relative_max=axial_closure,
        backend_kind=backend_kind,
    )


def materialize_oscillatory_quadratic_nonlinear_radial_stress(
    handoff_backend: ExactAgent2OscillatoryNonlinearHandoff,
    geometry: StrictInnerTransportRadialGeometry,
) -> OscillatoryQuadraticNonlinearRadialStressWitness:
    """Recompute #915 typed means and route quadratic/mixed/aggregate pieces."""
    if not isinstance(handoff_backend, ExactAgent2OscillatoryNonlinearHandoff):
        raise TypeError("handoff_backend must be ExactAgent2OscillatoryNonlinearHandoff")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    radii = np.asarray(geometry.radii, dtype=float)
    witness = materialize_oscillatory_nonlinear_mean_attribution(
        handoff_backend,
        radii,
        np.full_like(radii, float(geometry.axial_z)),
        np.full_like(radii, float(geometry.time)),
    )
    return _materialize_from_mean_witness(
        witness,
        geometry,
        backend_kind="exact-a3-915-quadratic-mean-to-a3-875-radial-stress",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_oscillatory_quadratic_nonlinear_radial_stress)
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
        "angular_orders": list(ANGULAR_ORDERS),
        "strict_inner_only": True,
        "current_typed_nonlinear_decomposition_consumed": True,
        "oscillatory_quadratic_self_mean_materialized": True,
        "mixed_cross_mean_materialized": True,
        "aggregate_nonlinear_mean_materialized": True,
        "scoped_quadratic_nonlinear_radial_inverse_performed": True,
        "radial_inverse_rule": "polynomial_first_cell_compact_moment_complement",
        "theta_exponent": 2,
        "axial_exponent": 1,
        "radial_component_recorded_but_not_inverted": True,
        "stress_level_quadratic_mixed_attribution_checked": True,
        "mean_piece_closure_absolute_gate": MEAN_PIECE_CLOSURE_ABSOLUTE_GATE,
        "stress_piece_closure_relative_gate": STRESS_PIECE_CLOSURE_RELATIVE_GATE,
        "mean_recomputed_from_exact_upstream_in_public_api": True,
        "agent2_curl_or_advection_reimplemented_by_agent3": False,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "global_corrected_leading_join_materialized": False,
        "complete_ns_defect": False,
        "scoped_quadratic_radial_stress_authorized_as_correction_target": False,
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


def build_mechanics_report() -> dict[str, object]:
    """Deterministic mechanics-only replay; never candidate residual evidence."""
    radii = np.linspace(0.02, 0.80, 79)
    geometry = StrictInnerTransportRadialGeometry(
        time=0.50,
        axial_z=0.13,
        radii=tuple(float(v) for v in radii),
        bump_center=0.50,
        bump_halfwidth=0.28,
    )
    mean = _materialize_from_provider(
        _analytic_provider,
        radii,
        np.full_like(radii, geometry.axial_z),
        np.full_like(radii, geometry.time),
        backend=None,
        backend_kind="analytic-projection-regression-only",
    )
    witness = _materialize_from_mean_witness(
        mean,
        geometry,
        backend_kind="manufactured-mechanics-only",
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
        "radial_operator_agent3_head": RADIAL_OPERATOR_AGENT3_HEAD,
        "registered_quadratic_mean": [0.17, -0.03, 0.09],
        "registered_mixed_mean": [0.13, -0.04, 0.25],
        "registered_aggregate_mean": [0.30, -0.07, 0.34],
        "witness": witness.to_receipt(),
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "truth_boundary": truth_boundary(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = build_mechanics_report()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
