"""Materialize current-I4 radial force from exact compact axial stress.

The corrected 2026-09-09 Kokuno symmetric compact-stress representation uses

    (div T)_r = partial_z sigma_1,

where ``sigma_1`` is the axial/e=1 compact stress. Agent-3 PR #1109 already
materializes exact current-I4 mixed/quadratic/aggregate compact stresses on the
identity-bound Agent-2 #1080 / Agent-1 #1079 leading+oscillatory candidate.
This module performs only that next divergence step by recomputing the same
exact stress at neighboring physical-z slices.

No independent third radial moment inverse is introduced. This remains scoped
correction-side machinery: it is not the source RF30--RF49 five-row mean
correction, not a complete Navier--Stokes defect, not an authorized correction
target, not a Cartesian correction velocity, and not a finite correction cycle.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import inspect
from pathlib import Path
from typing import Mapping

import numpy as np

from .kokuno_current_partial_nonlinear_radial_force import (
    FINE_PAIR_RELATIVE_STABILITY_GATE,
    RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE,
    Z_DERIVATIVE_STEP_LADDER,
    _differentiate_axial_stress_ladder,
    _rms,
)
from .kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_PR,
    AGENT1_LEADING_SOURCE_BLOB,
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_PR,
    AGENT2_COMPOSITE_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
    SOURCE_BLOB_SHA1,
    SOURCE_MEAN_DEFECT_FORMULAS,
    SOURCE_MEAN_INTERVAL_FORMULA,
    SOURCE_PATH,
    SOURCE_READER_DATE,
    SOURCE_READER_HEAD,
    SOURCE_REPOSITORY,
    _git_blob_sha1,
)
from .kokuno_current_i4_nonlinear_radial_stress import (
    CurrentI4NonlinearRadialStressWitness,
    materialize_current_i4_nonlinear_radial_stress,
)
from .kokuno_oscillatory_nonlinear_mean_attribution import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
)
from .kokuno_oscillatory_transport_radial_stress import _stress_array
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)

TASK = "KOKUNO-A3-CURRENT-I4-NONLINEAR-RADIAL-FORCE-121"
SCHEMA = "kokuno-a3-current-i4-nonlinear-radial-force-v1"

PARENT_AGENT3_PR = 1109
PARENT_AGENT3_HEAD = "49590ff311fef1dec4bd850b3013fd989485c87a"
PARENT_AGENT3_SOURCE_BLOB = "b8d9928934d4df74b52e8cd0d2c27c468257cdc9"

HISTORICAL_CURRENT_I2_RADIAL_FORCE_AGENT3_PR = 1090
HISTORICAL_CURRENT_I2_RADIAL_FORCE_AGENT3_HEAD = (
    "5e8512392df683375fabd132f35577e1b94773ff"
)

SOURCE_RADIAL_FORCE_FORMULA = "(div T)_r = partial_z sigma_1"
SOURCE_FORMULA_LOCATION = "R33-R34; retained under the R41 chart map"


def _authenticate_parent_source() -> str:
    source = inspect.getsourcefile(materialize_current_i4_nonlinear_radial_stress)
    if source is None:
        raise RuntimeError("current-I4 nonlinear radial-stress source is unavailable")
    blob = _git_blob_sha1(Path(source))
    if blob != PARENT_AGENT3_SOURCE_BLOB:
        raise RuntimeError("A3 #1109 current-I4 radial-stress source blob drifted")
    return blob


def _axial_stress_pieces(
    witness: CurrentI4NonlinearRadialStressWitness,
) -> dict[str, np.ndarray]:
    if not isinstance(witness, CurrentI4NonlinearRadialStressWitness):
        raise TypeError("witness must be CurrentI4NonlinearRadialStressWitness")
    inherited = witness.inherited
    return {
        "quadratic": _stress_array(
            inherited.axial_quadratic_stress, "quadratic axial"
        ),
        "mixed": _stress_array(inherited.axial_mixed_stress, "mixed axial"),
        "aggregate": _stress_array(
            inherited.axial_aggregate_stress, "aggregate axial"
        ),
    }


@dataclass(frozen=True)
class CurrentI4NonlinearRadialForceWitness:
    """Exact current-I4 compact stress plus source-required radial force."""

    geometry: StrictInnerTransportRadialGeometry
    center_stress: CurrentI4NonlinearRadialStressWitness
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

    def _piece_receipt(
        self, piece: str, levels: tuple[np.ndarray, ...]
    ) -> dict[str, object]:
        arrays = tuple(np.asarray(value, dtype=float) for value in levels)
        finest = arrays[-1]
        inherited = self.center_stress.inherited
        center_means = {
            "quadratic": inherited.radial_quadratic_mean_profile,
            "mixed": inherited.radial_mixed_mean_profile,
            "aggregate": inherited.radial_aggregate_mean_profile,
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
            "coarse_to_fine_relative_rms_differences": [
                float(value) for value in stability
            ],
            "finest_radial_force_rms": _rms(finest),
            "finest_radial_force_max_abs": float(np.max(np.abs(finest))),
            "center_recorded_radial_mean_rms": _rms(radial_mean),
            "force_to_recorded_radial_mean_rms_ratio": _rms(finest) / mean_scale,
        }

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "backend_kind": "exact-current-a3-pr1109-i4-axial-stress-to-radial-force",
            "geometry": {
                "time": float(self.geometry.time),
                "axial_z": float(self.geometry.axial_z),
                "radii": [float(value) for value in self.geometry.radii],
                "bump_center": float(self.geometry.bump_center),
                "bump_halfwidth": float(self.geometry.bump_halfwidth),
            },
            "source_radial_force_formula": SOURCE_RADIAL_FORCE_FORMULA,
            "source_formula_location": SOURCE_FORMULA_LOCATION,
            "z_derivative_step_ladder": [float(value) for value in self.z_steps],
            "quadratic": self._piece_receipt(
                "quadratic", self.quadratic_radial_force_levels
            ),
            "mixed": self._piece_receipt("mixed", self.mixed_radial_force_levels),
            "aggregate": self._piece_receipt(
                "aggregate", self.aggregate_radial_force_levels
            ),
            "piece_closure_relative_max_levels": [
                float(value) for value in self.piece_closure_relative_max_levels
            ],
            "aggregate_derivative_stability_preflight_passed": bool(
                self.aggregate_derivative_stability_preflight_passed
            ),
            "parent_center_radial_stress_receipt": self.center_stress.to_receipt(),
            "parent_agent3_1109_source_blob": self.parent_source_blob,
            "provenance": {
                "parent_agent3_pr": PARENT_AGENT3_PR,
                "parent_agent3_head": PARENT_AGENT3_HEAD,
                "historical_current_i2_radial_force_agent3_pr": (
                    HISTORICAL_CURRENT_I2_RADIAL_FORCE_AGENT3_PR
                ),
                "historical_current_i2_radial_force_agent3_head": (
                    HISTORICAL_CURRENT_I2_RADIAL_FORCE_AGENT3_HEAD
                ),
                "agent2_composite_pr": AGENT2_COMPOSITE_PR,
                "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
                "agent2_composite_source_blob": AGENT2_COMPOSITE_SOURCE_BLOB,
                "agent1_leading_pr": AGENT1_LEADING_PR,
                "agent1_leading_head": AGENT1_LEADING_HEAD,
                "agent1_leading_source_blob": AGENT1_LEADING_SOURCE_BLOB,
                "source_reader_repository": SOURCE_REPOSITORY,
                "source_reader_head": SOURCE_READER_HEAD,
                "source_reader_path": SOURCE_PATH,
                "source_reader_blob": SOURCE_BLOB_SHA1,
                "source_reader_date": SOURCE_READER_DATE,
                "source_mean_interval_formula": SOURCE_MEAN_INTERVAL_FORMULA,
                "source_mean_defect_formulas": SOURCE_MEAN_DEFECT_FORMULAS,
            },
            "truth_boundary": truth_boundary(),
        }


def materialize_current_i4_nonlinear_radial_force(
    backend: ExactCurrentI4NonlinearBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> CurrentI4NonlinearRadialForceWitness:
    """Recompute exact #1109 stresses at neighboring z and form ``partial_z sigma_1``."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    parent_blob = _authenticate_parent_source()
    center = materialize_current_i4_nonlinear_radial_stress(backend, geometry)

    def sampler(z_value: float) -> Mapping[str, np.ndarray]:
        shifted = replace(geometry, axial_z=float(z_value))
        witness = materialize_current_i4_nonlinear_radial_stress(backend, shifted)
        return _axial_stress_pieces(witness)

    differentiated = _differentiate_axial_stress_ladder(geometry, sampler)
    levels = differentiated["levels"]
    return CurrentI4NonlinearRadialForceWitness(
        geometry=geometry,
        center_stress=center,
        z_steps=tuple(float(value) for value in differentiated["steps"]),
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
            float(value) for value in differentiated["stability"]["quadratic"]
        ),
        mixed_stability_relative_rms=tuple(
            float(value) for value in differentiated["stability"]["mixed"]
        ),
        aggregate_stability_relative_rms=tuple(
            float(value) for value in differentiated["stability"]["aggregate"]
        ),
        aggregate_derivative_stability_preflight_passed=bool(
            differentiated["aggregate_derivative_stability_preflight_passed"]
        ),
        parent_source_blob=parent_blob,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_nonlinear_radial_force)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "radial_force",
        "pressure", "forcing", "target", "gain", "alpha", "damping",
        "spatial_step", "derivative_step", "z_step", "angular_order",
        "viscosity", "nu", "correction", "stage_budget", "normalized_score",
        "scientific_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "source_reader_repository": SOURCE_REPOSITORY,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_path": SOURCE_PATH,
        "source_reader_blob": SOURCE_BLOB_SHA1,
        "source_reader_date": SOURCE_READER_DATE,
        "source_reader_classification": "structural_math_provenance_only",
        "source_radial_force_formula": SOURCE_RADIAL_FORCE_FORMULA,
        "source_formula_location": SOURCE_FORMULA_LOCATION,
        "source_mean_interval_formula": SOURCE_MEAN_INTERVAL_FORMULA,
        "source_mean_defect_formulas": SOURCE_MEAN_DEFECT_FORMULAS,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "z_derivative_step_ladder": list(Z_DERIVATIVE_STEP_LADDER),
        "fine_pair_relative_stability_gate": FINE_PAIR_RELATIVE_STABILITY_GATE,
        "radial_force_piece_closure_relative_gate": (
            RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE
        ),
        "current_I4_leading_plus_oscillatory_identity_consumed": True,
        "current_I4_axial_e1_stress_consumed": True,
        "current_I4_radial_force_from_dz_sigma1_materialized": True,
        "quadratic_mixed_aggregate_radial_force_attribution_materialized": True,
        "independent_third_radial_moment_inverse_introduced": False,
        "recorded_radial_mean_equated_to_radial_force": False,
        "agent2_curl_or_jacobian_reimplemented_by_agent3": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_radial_force_allowed": False,
        "caller_supplied_defect_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "public_parameters": tuple(signature.parameters),
        "source_I4_five_row_mean_correction_materialized": False,
        "public_float64_I2_delta_verified": False,
        "post_I4_or_pulse_leading_identity_consumed": False,
        "outer_global_velocity_consumed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_defect": False,
        "scoped_current_I4_radial_force_authorized_as_complete_ns_defect": False,
        "scoped_current_I4_radial_force_authorized_as_correction_target": False,
        "current_I4_correction_velocity_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proof_claimed": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
