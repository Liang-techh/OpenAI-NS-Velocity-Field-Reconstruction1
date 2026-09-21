"""Route current RF40 axial-shutdown nonlinear m=0 means through compact radial stress.

This Agent-3 increment consumes the exact current-lineage mean attribution from
PR #1011,

    M_bar = <(u_lead . grad)u_osc + (u_osc . grad)u_lead>_theta,
    Q_bar = <(u_osc . grad)u_osc>_theta,
    N_bar = M_bar + Q_bar,

for the identity-bound Agent-2 PR #999 composite through the RF40 axial-shutdown
endpoint ``X_2``.  It reuses the already-reviewed Agent-3 compact
moment-complement radial-stress operator from PR #920 / PR #875.  No Agent-2
curl/Jacobian machinery and no Agent-1 RF40 leading construction is reproduced.

The inherited source-specific symmetric representation reconstructs tangential
``e=2`` and axial ``e=1`` scalar stresses.  The radial force belongs to the later
``partial_z sigma_1`` divergence step; this module deliberately does not invent
an independent third radial moment inverse.

This remains scoped correction-side machinery.  It is not a complete
Navier--Stokes defect, not an authorized correction target, not a Cartesian
correction velocity, and not held-out residual evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path

import numpy as np

from .kokuno_current_rf40_axial_shutdown_nonlinear_mean_attribution import (
    ExactCurrentRF40AxialShutdownNonlinearBackend,
    _git_blob_sha1,
    materialize_current_rf40_axial_shutdown_nonlinear_mean_attribution,
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

TASK = "KOKUNO-A3-CURRENT-RF40-AXIAL-SHUTDOWN-NONLINEAR-RADIAL-STRESS-110"
SCHEMA = "kokuno-a3-current-rf40-axial-shutdown-nonlinear-radial-stress-v1"

PARENT_AGENT3_PR = 1011
PARENT_AGENT3_HEAD = "6aabf7b6dd28ae683a3776d518e9a977507911e9"
PARENT_AGENT3_SOURCE_BLOB = "1ccef44607487d7c3cbdfa6dbce0067fdd302ea5"

RADIAL_ROUTING_AGENT3_PR = 920
RADIAL_ROUTING_AGENT3_HEAD = "3ca100009e770f2aecb5859c721bdd99ffb8552c"
RADIAL_OPERATOR_AGENT3_PR = 875
RADIAL_OPERATOR_AGENT3_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"

AGENT2_COMPOSITE_PR = 999
AGENT2_COMPOSITE_HEAD = "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5"
AGENT2_COMPOSITE_SOURCE_BLOB = "255e65ce0c37652858509c33e7c9fad40b73ca97"
AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"
AGENT1_LEADING_PR = 993
AGENT1_LEADING_HEAD = "2ac6460b483efb1c07f2fa65e7fed781a32f2718"
AGENT1_LEADING_SOURCE_BLOB = "057a0514c8a941c3e59922b8158f480434b4441e"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_READER_DATE = "2026-09-09"


def _authenticate_parent_source() -> str:
    source = inspect.getsourcefile(
        materialize_current_rf40_axial_shutdown_nonlinear_mean_attribution
    )
    if source is None:
        raise RuntimeError("RF40 axial-shutdown nonlinear mean source is unavailable")
    blob = _git_blob_sha1(Path(source))
    if blob != PARENT_AGENT3_SOURCE_BLOB:
        raise RuntimeError("A3 #1011 RF40 axial-shutdown nonlinear mean source blob drifted")
    return blob


@dataclass(frozen=True)
class CurrentRF40AxialShutdownNonlinearRadialStressWitness:
    """Current through-X2 mean attribution plus inherited compact radial stresses."""

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
            "backend_kind": (
                "exact-current-a3-pr1011-mean-to-inherited-compact-radial-stress"
            ),
            "geometry": base["geometry"],
            "current_rf40_axial_shutdown_mean_witness": base["mean_witness"],
            "radial_quadratic_mean_profile": base["radial_quadratic_mean_profile"],
            "radial_mixed_mean_profile": base["radial_mixed_mean_profile"],
            "radial_aggregate_mean_profile": base["radial_aggregate_mean_profile"],
            "theta_e2": base["theta_e2"],
            "axial_e1": base["axial_e1"],
            "mean_piece_closure_absolute_max": self.mean_piece_closure_absolute_max,
            "parent_agent3_1011_source_blob": self.parent_source_blob,
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
                "source_reader_repository": SOURCE_READER_REPO,
                "source_reader_head": SOURCE_READER_HEAD,
                "source_reader_path": SOURCE_READER_PATH,
                "source_reader_blob": SOURCE_READER_BLOB,
                "source_reader_date": SOURCE_READER_DATE,
            },
            "truth_boundary": truth_boundary(),
        }


def _materialize_current_mean_witness(
    witness: OscillatoryNonlinearMeanAttributionWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    parent_source_blob: str,
) -> CurrentRF40AxialShutdownNonlinearRadialStressWitness:
    """Apply only the inherited compact radial operator to an exact #1011 mean witness."""
    if not isinstance(witness, OscillatoryNonlinearMeanAttributionWitness):
        raise TypeError("witness must be OscillatoryNonlinearMeanAttributionWitness")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")
    if parent_source_blob != PARENT_AGENT3_SOURCE_BLOB:
        raise ValueError("parent source blob does not match exact A3 #1011")

    inherited = _materialize_from_mean_witness(
        witness,
        geometry,
        backend_kind=(
            "current-a3-pr1011-mean-reusing-a3-pr920-a3-pr875-radial-operator"
        ),
    )
    if inherited.mean_witness is not witness:
        raise RuntimeError("inherited radial routing detached the current #1011 mean witness")
    return CurrentRF40AxialShutdownNonlinearRadialStressWitness(
        inherited=inherited,
        parent_source_blob=parent_source_blob,
    )


def materialize_current_rf40_axial_shutdown_nonlinear_radial_stress(
    backend: ExactCurrentRF40AxialShutdownNonlinearBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> CurrentRF40AxialShutdownNonlinearRadialStressWitness:
    """Recompute #1011 means and route them through the compact radial inverse."""
    if not isinstance(backend, ExactCurrentRF40AxialShutdownNonlinearBackend):
        raise TypeError("backend must be ExactCurrentRF40AxialShutdownNonlinearBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    parent_blob = _authenticate_parent_source()
    radii = np.asarray(geometry.radii, dtype=float)
    z = np.full_like(radii, float(geometry.axial_z))
    t = np.full_like(radii, float(geometry.time))
    witness = materialize_current_rf40_axial_shutdown_nonlinear_mean_attribution(
        backend, radii, z, t
    )
    if witness.backend_kind != (
        "exact-current-a2-pr999-through-rf40-axial-shutdown-nonlinear-attribution"
    ):
        raise RuntimeError("RF40 axial-shutdown nonlinear mean backend identity drifted")
    return _materialize_current_mean_witness(
        witness,
        geometry,
        parent_source_blob=parent_blob,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_current_rf40_axial_shutdown_nonlinear_radial_stress
    )
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "spatial_step",
        "derivative_step", "angular_order", "time_step", "viscosity", "nu",
        "correction", "stage_budget", "normalized_score", "scientific_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_path": SOURCE_READER_PATH,
        "source_reader_blob": SOURCE_READER_BLOB,
        "source_reader_date": SOURCE_READER_DATE,
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
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "current_leading_plus_oscillation_through_RF40_axial_shutdown_consumed": True,
        "current_RF40_axial_shutdown_mixed_nonlinear_mean_consumed": True,
        "current_RF40_axial_shutdown_quadratic_nonlinear_mean_consumed": True,
        "current_RF40_axial_shutdown_aggregate_nonlinear_mean_consumed": True,
        "current_RF40_axial_shutdown_nonlinear_radial_inverse_performed": True,
        "current_RF40_axial_shutdown_quadratic_radial_stress_materialized": True,
        "current_RF40_axial_shutdown_mixed_radial_stress_materialized": True,
        "current_RF40_axial_shutdown_aggregate_radial_stress_materialized": True,
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
        "partial_domain_through_RF40_axial_shutdown_only": True,
        "velocity_after_RF40_axial_shutdown_materialized": False,
        "RF40_lambda_turn_current_lineage_consumed": False,
        "RF40_power_law_current_lineage_consumed": False,
        "full_post_XR_RF40_current_lineage_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "pressure_gradient_included": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_included": False,
        "complete_ns_defect": False,
        "scoped_current_RF40_axial_shutdown_radial_stress_authorized_as_correction_target": False,
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
