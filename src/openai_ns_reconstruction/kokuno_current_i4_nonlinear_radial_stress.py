"""Route exact current-I4 nonlinear m=0 means through compact radial stress.

This Kokuno Agent-3 increment consumes the exact current-I4 nonlinear mean
attribution from PR #1101,

    M_bar = <(u_lead . grad)u_osc + (u_osc . grad)u_lead>_theta,
    Q_bar = <(u_osc . grad)u_osc>_theta,
    N_bar = M_bar + Q_bar,

for the identity-bound Agent-2 PR #1080 / Agent-1 PR #1079 candidate.  It
reuses the already-reviewed Agent-3 compact moment-complement radial-stress
operator from PR #920 / PR #875.  No Agent-2 curl/Jacobian machinery is
reimplemented.

The inherited source-specific symmetric representation reconstructs tangential
``e=2`` and axial ``e=1`` scalar stresses.  The radial force belongs to the
later ``partial_z sigma_1`` divergence step; this module deliberately does not
combine that derivative into the same increment and does not invent an
independent third radial moment inverse.

RF40a records I4 as the reserved mean-correction interval, while RF30--RF49
specify the later source five-row correction/recomputation cycle.  This module
only reconstructs compact radial stresses from the actual current-I4 nonlinear
mean.  It does not claim that the RF30--RF49 five-row correction has been
executed, and it remains scoped correction-side machinery rather than a
complete Navier--Stokes defect or correction target.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path

import numpy as np

from .kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_PR,
    AGENT1_LEADING_SOURCE_BLOB,
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_PR,
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_HEAD,
    AGENT2_DIFFERENTIAL_PR,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
    SOURCE_BLOB_SHA1,
    SOURCE_MEAN_DEFECT_FORMULAS,
    SOURCE_MEAN_INTERVAL_FORMULA,
    SOURCE_PATH,
    SOURCE_READER_DATE,
    SOURCE_READER_HEAD,
    SOURCE_REPOSITORY,
    _git_blob_sha1,
    materialize_current_i4_nonlinear_mean_attribution,
)
from .kokuno_oscillatory_nonlinear_mean_attribution import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    OscillatoryNonlinearMeanAttributionWitness,
)
from .kokuno_oscillatory_quadratic_nonlinear_radial_stress import (
    MEAN_PIECE_CLOSURE_ABSOLUTE_GATE,
    STRESS_PIECE_CLOSURE_RELATIVE_GATE,
    OscillatoryQuadraticNonlinearRadialStressWitness,
    _materialize_from_mean_witness,
)
from .kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)

TASK = "KOKUNO-A3-CURRENT-I4-NONLINEAR-RADIAL-STRESS-120"
SCHEMA = "kokuno-a3-current-i4-nonlinear-radial-stress-v1"

PARENT_AGENT3_PR = 1101
PARENT_AGENT3_HEAD = "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b"
PARENT_AGENT3_SOURCE_BLOB = "10d067570ddba834db3f19c9859c2efdad16b3e8"

RADIAL_ROUTING_AGENT3_PR = 920
RADIAL_ROUTING_AGENT3_HEAD = "3ca100009e770f2aecb5859c721bdd99ffb8552c"
RADIAL_OPERATOR_AGENT3_PR = 875
RADIAL_OPERATOR_AGENT3_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"


def _authenticate_parent_source() -> str:
    source = inspect.getsourcefile(materialize_current_i4_nonlinear_mean_attribution)
    if source is None:
        raise RuntimeError("current-I4 nonlinear mean source is unavailable")
    blob = _git_blob_sha1(Path(source))
    if blob != PARENT_AGENT3_SOURCE_BLOB:
        raise RuntimeError("A3 #1101 current-I4 nonlinear mean source blob drifted")
    return blob


@dataclass(frozen=True)
class CurrentI4NonlinearRadialStressWitness:
    """Exact current-I4 mean attribution plus inherited compact radial stresses."""

    inherited: OscillatoryQuadraticNonlinearRadialStressWitness
    parent_source_blob: str

    @property
    def mean_witness(self) -> OscillatoryNonlinearMeanAttributionWitness:
        return self.inherited.mean_witness

    @property
    def mean_piece_closure_absolute_max(self) -> float:
        return float(self.inherited.mean_piece_closure_absolute_max)

    @property
    def theta_stress_piece_closure_relative_max(self) -> float:
        return float(self.inherited.theta_stress_piece_closure_relative_max)

    @property
    def axial_stress_piece_closure_relative_max(self) -> float:
        return float(self.inherited.axial_stress_piece_closure_relative_max)

    def to_receipt(self) -> dict[str, object]:
        base = self.inherited.to_receipt()
        return {
            "schema": SCHEMA,
            "task": TASK,
            "backend_kind": "exact-current-a3-pr1101-mean-to-inherited-compact-radial-stress",
            "geometry": base["geometry"],
            "current_i4_mean_witness": base["mean_witness"],
            "radial_quadratic_mean_profile": base["radial_quadratic_mean_profile"],
            "radial_mixed_mean_profile": base["radial_mixed_mean_profile"],
            "radial_aggregate_mean_profile": base["radial_aggregate_mean_profile"],
            "theta_e2": base["theta_e2"],
            "axial_e1": base["axial_e1"],
            "mean_piece_closure_absolute_max": self.mean_piece_closure_absolute_max,
            "parent_agent3_1101_source_blob": self.parent_source_blob,
            "provenance": {
                "parent_agent3_pr": PARENT_AGENT3_PR,
                "parent_agent3_head": PARENT_AGENT3_HEAD,
                "radial_routing_agent3_pr": RADIAL_ROUTING_AGENT3_PR,
                "radial_routing_agent3_head": RADIAL_ROUTING_AGENT3_HEAD,
                "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
                "radial_operator_agent3_head": RADIAL_OPERATOR_AGENT3_HEAD,
                "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
                "agent2_composite_pr": AGENT2_COMPOSITE_PR,
                "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
                "agent2_composite_source_blob": AGENT2_COMPOSITE_SOURCE_BLOB,
                "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
                "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
                "agent2_differential_source_blob": AGENT2_DIFFERENTIAL_SOURCE_BLOB,
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


def _materialize_current_mean_witness(
    witness: OscillatoryNonlinearMeanAttributionWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    parent_source_blob: str,
) -> CurrentI4NonlinearRadialStressWitness:
    """Apply only the inherited compact radial operator to an exact #1101 mean."""
    if not isinstance(witness, OscillatoryNonlinearMeanAttributionWitness):
        raise TypeError("witness must be OscillatoryNonlinearMeanAttributionWitness")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")
    if parent_source_blob != PARENT_AGENT3_SOURCE_BLOB:
        raise ValueError("parent source blob does not match exact A3 #1101")

    inherited = _materialize_from_mean_witness(
        witness,
        geometry,
        backend_kind="current-a3-pr1101-mean-reusing-a3-pr920-a3-pr875-radial-operator",
    )
    if inherited.mean_witness is not witness:
        raise RuntimeError("inherited radial routing detached the current #1101 mean witness")
    return CurrentI4NonlinearRadialStressWitness(
        inherited=inherited,
        parent_source_blob=parent_source_blob,
    )


def materialize_current_i4_nonlinear_radial_stress(
    backend: ExactCurrentI4NonlinearBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> CurrentI4NonlinearRadialStressWitness:
    """Recompute exact #1101 means and route them through compact radial stress."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    parent_blob = _authenticate_parent_source()
    radii = np.asarray(geometry.radii, dtype=float)
    z = np.full_like(radii, float(geometry.axial_z))
    t = np.full_like(radii, float(geometry.time))
    witness = materialize_current_i4_nonlinear_mean_attribution(backend, radii, z, t)
    if witness.backend_kind != "exact-current-a2-pr1080-through-i4-nonlinear-attribution":
        raise RuntimeError("current-I4 nonlinear mean backend identity drifted")
    return _materialize_current_mean_witness(
        witness,
        geometry,
        parent_source_blob=parent_blob,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_nonlinear_radial_stress)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "spatial_step",
        "derivative_step", "angular_order", "time_step", "viscosity", "nu",
        "correction", "stage_budget", "normalized_score", "scientific_threshold",
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
        "source_mean_interval_formula": SOURCE_MEAN_INTERVAL_FORMULA,
        "source_mean_defect_formulas": SOURCE_MEAN_DEFECT_FORMULAS,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "radial_routing_agent3_pr": RADIAL_ROUTING_AGENT3_PR,
        "radial_routing_agent3_head": RADIAL_ROUTING_AGENT3_HEAD,
        "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
        "radial_operator_agent3_head": RADIAL_OPERATOR_AGENT3_HEAD,
        "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_composite_source_blob": AGENT2_COMPOSITE_SOURCE_BLOB,
        "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "agent2_differential_source_blob": AGENT2_DIFFERENTIAL_SOURCE_BLOB,
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "agent1_leading_source_blob": AGENT1_LEADING_SOURCE_BLOB,
        "current_I4_leading_plus_oscillatory_identity_consumed": True,
        "current_I4_mixed_nonlinear_mean_consumed": True,
        "current_I4_quadratic_oscillatory_mean_consumed": True,
        "current_I4_aggregate_nonlinear_mean_consumed": True,
        "source_I4_reserved_mean_correction_interval_recorded": True,
        "current_I4_nonlinear_radial_inverse_performed": True,
        "current_I4_quadratic_radial_stress_materialized": True,
        "current_I4_mixed_radial_stress_materialized": True,
        "current_I4_aggregate_radial_stress_materialized": True,
        "radial_inverse_rule": "polynomial_first_cell_compact_moment_complement",
        "theta_exponent": 2,
        "axial_exponent": 1,
        "independent_third_radial_moment_inverse_required": False,
        "radial_force_requires_later_partial_z_sigma_1_step": True,
        "moment_complement_preserved_by_inherited_operator": True,
        "polynomial_first_cell_axis_extrapolation_inherited": True,
        "outer_edge_stress_closure_inherited": True,
        "mean_piece_closure_absolute_gate": MEAN_PIECE_CLOSURE_ABSOLUTE_GATE,
        "stress_piece_closure_relative_gate": STRESS_PIECE_CLOSURE_RELATIVE_GATE,
        "agent2_curl_or_jacobian_reimplemented_by_agent3": False,
        "caller_supplied_mean_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_defect_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "public_parameters": tuple(signature.parameters),
        "source_I4_five_row_mean_correction_materialized": False,
        "current_I4_correction_velocity_materialized": False,
        "post_I4_or_pulse_leading_identity_consumed": False,
        "outer_global_velocity_consumed": False,
        "radial_force_materialized_in_this_increment": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_defect": False,
        "scoped_current_I4_radial_stress_authorized_as_complete_ns_defect": False,
        "scoped_current_I4_radial_stress_authorized_as_correction_target": False,
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
