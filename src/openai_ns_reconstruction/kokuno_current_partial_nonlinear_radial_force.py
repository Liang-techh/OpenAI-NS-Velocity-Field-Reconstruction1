"""Materialize the current partial radial force implied by compact axial stress.

The corrected 2026-09-09 Kokuno reconstruction uses a symmetric compact stress.
For the two scalar compact profiles already reconstructed by Agent 3, its
cylindrical divergence contains the radial component

    (div T)_r = partial_z sigma_1,

where ``sigma_1`` is the axial/e=1 compact radial stress.  Therefore the radial
force is not a third independently chosen moment-complement inverse.  It must be
retained by differentiating the same axial stress in the physical axial
coordinate.  Historical Agent-3 PR #357 exposed this source seam on an older
candidate; this module applies it to the current #976 lineage without reusing
those historical candidate numbers.

This remains a partial-domain nonlinear-transport diagnostic through the current
leading-field scope.  It is not a complete Navier--Stokes defect, not an
authorized correction target, and not a correction velocity.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import inspect
from pathlib import Path
from typing import Callable, Mapping

import numpy as np

from .kokuno_current_partial_nonlinear_mean_attribution import (
    ExactCurrentPartialNonlinearBackend,
    _git_blob_sha1,
)
from .kokuno_current_partial_nonlinear_radial_stress import (
    CurrentPartialNonlinearRadialStressWitness,
    materialize_current_partial_nonlinear_radial_stress,
)
from .kokuno_oscillatory_nonlinear_mean_attribution import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
)
from .kokuno_oscillatory_transport_radial_stress import _stress_array
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)

TASK = "KOKUNO-A3-CURRENT-PARTIAL-NONLINEAR-RADIAL-FORCE-104"
SCHEMA = "kokuno-a3-current-partial-nonlinear-radial-force-v1"

PARENT_AGENT3_PR = 976
PARENT_AGENT3_HEAD = "0a777668050ff9e73d6c8db1f72e5fabfa9b04a6"
PARENT_AGENT3_SOURCE_BLOB = "9afa50c573186b86cd86fa7ffe44a5986a808b6d"

HISTORICAL_RADIAL_FORCE_AGENT3_PR = 357
HISTORICAL_RADIAL_FORCE_AGENT3_HEAD = "169558a6994a255119d6af97430e717a221c034c"
CR002_CHANNEL_COMPLETENESS_PR = 979

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_RADIAL_FORCE_FORMULA = "(div T)_r = partial_z sigma_1"
SOURCE_FORMULA_LOCATION = "R33-R34; retained under the R41 chart map"

# Reuse the historical A3 radial-force audit ladder.  These are fixed numerical
# differentiation controls, not Navier--Stokes acceptance thresholds.
Z_DERIVATIVE_STEP_LADDER = (0.02, 0.01, 0.005)
FINE_PAIR_RELATIVE_STABILITY_GATE = 5.0e-2
RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE = 5.0e-10


def _rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not np.all(np.isfinite(array)):
        raise ValueError("RMS input must be a finite one-dimensional array")
    return float(np.sqrt(np.mean(array * array)))


def _relative_max_closure(lhs: np.ndarray, rhs: np.ndarray) -> float:
    a = np.asarray(lhs, dtype=float)
    b = np.asarray(rhs, dtype=float)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("closure arrays must be matching one-dimensional arrays")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("closure arrays must be finite")
    scale = max(1.0, float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    return float(np.max(np.abs(a - b)) / scale)


def centered_axial_stress_derivative(
    sigma_plus: np.ndarray,
    sigma_minus: np.ndarray,
    step: float,
) -> np.ndarray:
    """Centered physical-z derivative of one axial compact-stress profile."""
    plus = np.asarray(sigma_plus, dtype=float)
    minus = np.asarray(sigma_minus, dtype=float)
    step = float(step)
    if plus.shape != minus.shape or plus.ndim != 1 or plus.size < 9:
        raise ValueError("stress profiles must be matching 1-D arrays with >=9 nodes")
    if not np.all(np.isfinite(plus)) or not np.all(np.isfinite(minus)):
        raise ValueError("stress profiles must be finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("z derivative step must be positive and finite")
    return (plus - minus) / (2.0 * step)


def _validate_step_ladder(steps: tuple[float, ...]) -> tuple[float, ...]:
    checked = tuple(float(value) for value in steps)
    if len(checked) < 2:
        raise ValueError("z derivative ladder requires at least two levels")
    if not all(np.isfinite(value) and value > 0.0 for value in checked):
        raise ValueError("z derivative steps must be positive and finite")
    if not all(right < left for left, right in zip(checked, checked[1:])):
        raise ValueError("z derivative steps must be strictly decreasing")
    return checked


def _differentiate_axial_stress_ladder(
    geometry: StrictInnerTransportRadialGeometry,
    sampler: Callable[[float], Mapping[str, np.ndarray]],
    *,
    steps: tuple[float, ...] = Z_DERIVATIVE_STEP_LADDER,
) -> dict[str, object]:
    """Differentiate quadratic/mixed/aggregate sigma_1 on a fixed z ladder.

    ``sampler`` is internal/test-only.  The public materializer below always
    supplies it by recomputing exact #976 witnesses from the typed current
    candidate backend, so callers cannot inject a radial force or stress.
    """
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")
    checked_steps = _validate_step_ladder(steps)
    expected_shape = (len(geometry.radii),)
    pieces = ("quadratic", "mixed", "aggregate")
    levels: list[dict[str, object]] = []

    for step in checked_steps:
        plus = sampler(float(geometry.axial_z + step))
        minus = sampler(float(geometry.axial_z - step))
        derivatives: dict[str, np.ndarray] = {}
        for piece in pieces:
            if piece not in plus or piece not in minus:
                raise ValueError(f"axial-stress sampler omitted {piece}")
            plus_array = np.asarray(plus[piece], dtype=float)
            minus_array = np.asarray(minus[piece], dtype=float)
            if plus_array.shape != expected_shape or minus_array.shape != expected_shape:
                raise ValueError("axial-stress sampler returned the wrong radial shape")
            derivatives[piece] = centered_axial_stress_derivative(
                plus_array, minus_array, step
            )

        closure = _relative_max_closure(
            derivatives["aggregate"],
            derivatives["mixed"] + derivatives["quadratic"],
        )
        if closure > RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE:
            raise ValueError("radial-force nonlinear piece attribution closure failed")
        levels.append(
            {
                "z_step": step,
                "z_minus": float(geometry.axial_z - step),
                "z_plus": float(geometry.axial_z + step),
                "quadratic": derivatives["quadratic"],
                "mixed": derivatives["mixed"],
                "aggregate": derivatives["aggregate"],
                "piece_closure_relative_max": closure,
            }
        )

    stability: dict[str, tuple[float, ...]] = {}
    for piece in pieces:
        comparisons: list[float] = []
        for coarse, fine in zip(levels, levels[1:]):
            coarse_array = np.asarray(coarse[piece], dtype=float)
            fine_array = np.asarray(fine[piece], dtype=float)
            scale = max(_rms(fine_array), 1.0e-14)
            comparisons.append(_rms(coarse_array - fine_array) / scale)
        stability[piece] = tuple(comparisons)

    aggregate_finest_pair = float(stability["aggregate"][-1])
    return {
        "steps": checked_steps,
        "levels": tuple(levels),
        "stability": stability,
        "aggregate_finest_pair_relative_rms_difference": aggregate_finest_pair,
        "aggregate_derivative_stability_preflight_passed": bool(
            aggregate_finest_pair <= FINE_PAIR_RELATIVE_STABILITY_GATE
        ),
    }


def _authenticate_parent_source() -> str:
    source = inspect.getsourcefile(materialize_current_partial_nonlinear_radial_stress)
    if source is None:
        raise RuntimeError("current partial radial-stress source is unavailable")
    blob = _git_blob_sha1(Path(source))
    if blob != PARENT_AGENT3_SOURCE_BLOB:
        raise RuntimeError("A3 #976 current partial radial-stress source blob drifted")
    return blob


def _axial_stress_pieces(
    witness: CurrentPartialNonlinearRadialStressWitness,
) -> dict[str, np.ndarray]:
    if not isinstance(witness, CurrentPartialNonlinearRadialStressWitness):
        raise TypeError("witness must be CurrentPartialNonlinearRadialStressWitness")
    inherited = witness.inherited
    return {
        "quadratic": _stress_array(inherited.axial_quadratic_stress, "quadratic axial"),
        "mixed": _stress_array(inherited.axial_mixed_stress, "mixed axial"),
        "aggregate": _stress_array(inherited.axial_aggregate_stress, "aggregate axial"),
    }


@dataclass(frozen=True)
class CurrentPartialNonlinearRadialForceWitness:
    """Current partial compact stress plus source-required radial force."""

    geometry: StrictInnerTransportRadialGeometry
    center_stress: CurrentPartialNonlinearRadialStressWitness
    z_steps: tuple[float, ...]
    quadratic_radial_force_levels: tuple[np.ndarray, ...]
    mixed_radial_force_levels: tuple[np.ndarray, ...]
    aggregate_radial_force_levels: tuple[np.ndarray, ...]
    piece_closure_relative_max_levels: tuple[float, ...]
    quadratic_stability_relative_rms: tuple[float, ...]
    mixed_stability_relative_rms: tuple[float, ...]
    aggregate_stability_relative_rms: tuple[float, ...]
    aggregate_derivative_stability_preflight_passed: bool
    parent_source_blob: str

    @property
    def finest_quadratic_radial_force(self) -> np.ndarray:
        return np.asarray(self.quadratic_radial_force_levels[-1], dtype=float)

    @property
    def finest_mixed_radial_force(self) -> np.ndarray:
        return np.asarray(self.mixed_radial_force_levels[-1], dtype=float)

    @property
    def finest_aggregate_radial_force(self) -> np.ndarray:
        return np.asarray(self.aggregate_radial_force_levels[-1], dtype=float)

    def _piece_receipt(self, piece: str, levels: tuple[np.ndarray, ...]) -> dict[str, object]:
        arrays = tuple(np.asarray(value, dtype=float) for value in levels)
        finest = arrays[-1]
        center_means = {
            "quadratic": self.center_stress.inherited.radial_quadratic_mean_profile,
            "mixed": self.center_stress.inherited.radial_mixed_mean_profile,
            "aggregate": self.center_stress.inherited.radial_aggregate_mean_profile,
        }
        radial_mean = np.asarray(center_means[piece], dtype=float)
        mean_scale = max(_rms(radial_mean), 1.0e-14)
        stability = {
            "quadratic": self.quadratic_stability_relative_rms,
            "mixed": self.mixed_stability_relative_rms,
            "aggregate": self.aggregate_stability_relative_rms,
        }[piece]
        return {
            "levels": [
                {
                    "z_step": float(step),
                    "radial_force": array.tolist(),
                    "radial_force_rms": _rms(array),
                    "radial_force_max_abs": float(np.max(np.abs(array))),
                }
                for step, array in zip(self.z_steps, arrays)
            ],
            "coarse_to_fine_relative_rms_differences": [float(v) for v in stability],
            "finest_radial_force_rms": _rms(finest),
            "finest_radial_force_max_abs": float(np.max(np.abs(finest))),
            "center_recorded_radial_mean_rms": _rms(radial_mean),
            "force_to_recorded_radial_mean_rms_ratio": _rms(finest) / mean_scale,
        }

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "geometry": {
                "time": float(self.geometry.time),
                "axial_z": float(self.geometry.axial_z),
                "radii": [float(v) for v in self.geometry.radii],
                "bump_center": float(self.geometry.bump_center),
                "bump_halfwidth": float(self.geometry.bump_halfwidth),
            },
            "source_radial_force_formula": SOURCE_RADIAL_FORCE_FORMULA,
            "source_formula_location": SOURCE_FORMULA_LOCATION,
            "z_derivative_step_ladder": [float(v) for v in self.z_steps],
            "quadratic": self._piece_receipt(
                "quadratic", self.quadratic_radial_force_levels
            ),
            "mixed": self._piece_receipt("mixed", self.mixed_radial_force_levels),
            "aggregate": self._piece_receipt(
                "aggregate", self.aggregate_radial_force_levels
            ),
            "piece_closure_relative_max_levels": [
                float(v) for v in self.piece_closure_relative_max_levels
            ],
            "aggregate_derivative_stability_preflight_passed": bool(
                self.aggregate_derivative_stability_preflight_passed
            ),
            "parent_center_radial_stress_receipt": self.center_stress.to_receipt(),
            "parent_agent3_976_source_blob": self.parent_source_blob,
            "provenance": {
                "parent_agent3_pr": PARENT_AGENT3_PR,
                "parent_agent3_head": PARENT_AGENT3_HEAD,
                "historical_radial_force_agent3_pr": HISTORICAL_RADIAL_FORCE_AGENT3_PR,
                "historical_radial_force_agent3_head": HISTORICAL_RADIAL_FORCE_AGENT3_HEAD,
                "cr002_channel_completeness_pr": CR002_CHANNEL_COMPLETENESS_PR,
                "source_reader_repository": SOURCE_READER_REPO,
                "source_reader_head": SOURCE_READER_HEAD,
                "source_reader_path": SOURCE_READER_PATH,
                "source_reader_date": SOURCE_READER_DATE,
            },
            "truth_boundary": truth_boundary(),
        }


def materialize_current_partial_nonlinear_radial_force(
    backend: ExactCurrentPartialNonlinearBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> CurrentPartialNonlinearRadialForceWitness:
    """Recompute #976 stresses at neighboring z and form ``partial_z sigma_1``."""
    if not isinstance(backend, ExactCurrentPartialNonlinearBackend):
        raise TypeError("backend must be ExactCurrentPartialNonlinearBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    parent_blob = _authenticate_parent_source()
    center = materialize_current_partial_nonlinear_radial_stress(backend, geometry)

    def sampler(z_value: float) -> Mapping[str, np.ndarray]:
        shifted = replace(geometry, axial_z=float(z_value))
        witness = materialize_current_partial_nonlinear_radial_stress(backend, shifted)
        return _axial_stress_pieces(witness)

    differentiated = _differentiate_axial_stress_ladder(geometry, sampler)
    levels = differentiated["levels"]
    return CurrentPartialNonlinearRadialForceWitness(
        geometry=geometry,
        center_stress=center,
        z_steps=tuple(float(v) for v in differentiated["steps"]),
        quadratic_radial_force_levels=tuple(
            np.asarray(level["quadratic"], dtype=float) for level in levels
        ),
        mixed_radial_force_levels=tuple(
            np.asarray(level["mixed"], dtype=float) for level in levels
        ),
        aggregate_radial_force_levels=tuple(
            np.asarray(level["aggregate"], dtype=float) for level in levels
        ),
        piece_closure_relative_max_levels=tuple(
            float(level["piece_closure_relative_max"]) for level in levels
        ),
        quadratic_stability_relative_rms=tuple(
            float(v) for v in differentiated["stability"]["quadratic"]
        ),
        mixed_stability_relative_rms=tuple(
            float(v) for v in differentiated["stability"]["mixed"]
        ),
        aggregate_stability_relative_rms=tuple(
            float(v) for v in differentiated["stability"]["aggregate"]
        ),
        aggregate_derivative_stability_preflight_passed=bool(
            differentiated["aggregate_derivative_stability_preflight_passed"]
        ),
        parent_source_blob=parent_blob,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_partial_nonlinear_radial_force)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "radial_force",
        "pressure", "forcing", "target", "gain", "alpha", "damping",
        "spatial_step", "derivative_step", "z_step", "angular_order",
        "viscosity", "nu", "correction", "stage_budget", "normalized_score",
        "scientific_threshold",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_path": SOURCE_READER_PATH,
        "source_reader_date": SOURCE_READER_DATE,
        "source_radial_force_formula": SOURCE_RADIAL_FORCE_FORMULA,
        "source_formula_location": SOURCE_FORMULA_LOCATION,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "historical_radial_force_agent3_pr": HISTORICAL_RADIAL_FORCE_AGENT3_PR,
        "cr002_channel_completeness_pr": CR002_CHANNEL_COMPLETENESS_PR,
        "z_derivative_step_ladder": list(Z_DERIVATIVE_STEP_LADDER),
        "fine_pair_relative_stability_gate": FINE_PAIR_RELATIVE_STABILITY_GATE,
        "radial_force_piece_closure_relative_gate": RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE,
        "current_partial_axial_e1_stress_consumed": True,
        "current_partial_radial_force_from_dz_sigma1_materialized": True,
        "quadratic_mixed_aggregate_radial_force_attribution_materialized": True,
        "independent_third_radial_moment_inverse_introduced": False,
        "recorded_radial_mean_equated_to_radial_force": False,
        "agent2_curl_or_jacobian_reimplemented_by_agent3": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_radial_force_allowed": False,
        "caller_supplied_defect_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "partial_domain_through_current_Xh_only": True,
        "velocity_beyond_Xh_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_defect": False,
        "full_three_component_ns_correction_object_materialized": False,
        "scoped_current_partial_stress_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
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
