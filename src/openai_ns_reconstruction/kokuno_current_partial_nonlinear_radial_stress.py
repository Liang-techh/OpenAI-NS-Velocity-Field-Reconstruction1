"""Route the current partial nonlinear m=0 means through compact radial stress.

This Agent-3 increment consumes the exact current-lineage mean attribution from
PR #971,

    A_bar = <(u_lead . grad) u_osc + (u_osc . grad) u_lead>_theta,
    Q_bar = <(u_osc . grad) u_osc>_theta,
    N_bar = A_bar + Q_bar,

and applies the already-existing Agent-3 polynomial-first-cell compact radial
moment-complement operator used by PR #920 / PR #875.  No Agent-2 curl or
Jacobian implementation is reproduced here.

Only the tangential and axial cylindrical channels are inverted, using the
already-registered exponents e=2 and e=1.  The radial channel is retained as a
diagnostic because the inherited two-channel operator does not reconstruct it.
This remains a current *partial-domain* leading+oscillatory transport
reconstruction through the present X_h.  It is not a complete Navier--Stokes
defect, an authorized correction target, or a correction velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path

import numpy as np

from .kokuno_current_partial_nonlinear_mean_attribution import (
    ExactCurrentPartialNonlinearBackend,
    _git_blob_sha1,
    materialize_current_partial_nonlinear_mean_attribution,
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

TASK = "KOKUNO-A3-CURRENT-PARTIAL-NONLINEAR-RADIAL-STRESS-103"
SCHEMA = "kokuno-a3-current-partial-nonlinear-radial-stress-v1"

PARENT_AGENT3_PR = 971
PARENT_AGENT3_HEAD = "bdb02ad487a175ab62748687c3172127d09aeb9c"
PARENT_AGENT3_SOURCE_BLOB = "9c62ae3c3b4b3ca2d9082f822d57feef700ea6cf"

RADIAL_ROUTING_AGENT3_PR = 920
RADIAL_ROUTING_AGENT3_HEAD = "3ca100009e770f2aecb5859c721bdd99ffb8552c"
RADIAL_OPERATOR_AGENT3_PR = 875
RADIAL_OPERATOR_AGENT3_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"

AGENT2_COMPOSITE_PR = 970
AGENT2_COMPOSITE_HEAD = "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT1_LEADING_PR = 965
AGENT1_LEADING_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"


def _authenticate_parent_source() -> str:
    source = inspect.getsourcefile(materialize_current_partial_nonlinear_mean_attribution)
    if source is None:
        raise RuntimeError("current partial nonlinear mean source is unavailable")
    blob = _git_blob_sha1(Path(source))
    if blob != PARENT_AGENT3_SOURCE_BLOB:
        raise RuntimeError("A3 #971 current partial nonlinear mean source blob drifted")
    return blob


@dataclass(frozen=True)
class CurrentPartialNonlinearRadialStressWitness:
    """Current-lineage mean attribution plus inherited compact radial stresses."""

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
        # Reuse only the inherited numerical radial operator.  Replace its
        # historical strict-inner truth boundary with this current partial one.
        return {
            "schema": SCHEMA,
            "task": TASK,
            "backend_kind": "exact-current-a3-pr971-mean-to-inherited-compact-radial-stress",
            "geometry": base["geometry"],
            "current_partial_mean_witness": base["mean_witness"],
            "radial_quadratic_mean_profile": base["radial_quadratic_mean_profile"],
            "radial_mixed_mean_profile": base["radial_mixed_mean_profile"],
            "radial_aggregate_mean_profile": base["radial_aggregate_mean_profile"],
            "theta_e2": base["theta_e2"],
            "axial_e1": base["axial_e1"],
            "mean_piece_closure_absolute_max": self.mean_piece_closure_absolute_max,
            "parent_agent3_971_source_blob": self.parent_source_blob,
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
                "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
                "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
                "agent1_leading_pr": AGENT1_LEADING_PR,
                "agent1_leading_head": AGENT1_LEADING_HEAD,
                "source_reader_repository": SOURCE_READER_REPO,
                "source_reader_head": SOURCE_READER_HEAD,
                "source_reader_date": SOURCE_READER_DATE,
            },
            "truth_boundary": truth_boundary(),
        }


def _materialize_current_mean_witness(
    witness: OscillatoryNonlinearMeanAttributionWitness,
    geometry: StrictInnerTransportRadialGeometry,
    *,
    parent_source_blob: str,
) -> CurrentPartialNonlinearRadialStressWitness:
    """Apply only the inherited compact radial operator to a current mean witness."""
    if not isinstance(witness, OscillatoryNonlinearMeanAttributionWitness):
        raise TypeError("witness must be OscillatoryNonlinearMeanAttributionWitness")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")
    if parent_source_blob != PARENT_AGENT3_SOURCE_BLOB:
        raise ValueError("parent source blob does not match exact A3 #971")

    inherited = _materialize_from_mean_witness(
        witness,
        geometry,
        backend_kind="current-a3-pr971-mean-reusing-a3-pr920-a3-pr875-radial-operator",
    )
    if inherited.mean_witness is not witness:
        raise RuntimeError("inherited radial routing detached the current mean witness")
    return CurrentPartialNonlinearRadialStressWitness(
        inherited=inherited,
        parent_source_blob=parent_source_blob,
    )


def materialize_current_partial_nonlinear_radial_stress(
    backend: ExactCurrentPartialNonlinearBackend,
    geometry: StrictInnerTransportRadialGeometry,
) -> CurrentPartialNonlinearRadialStressWitness:
    """Recompute #971 current means and route them through the compact radial inverse."""
    if not isinstance(backend, ExactCurrentPartialNonlinearBackend):
        raise TypeError("backend must be ExactCurrentPartialNonlinearBackend")
    if not isinstance(geometry, StrictInnerTransportRadialGeometry):
        raise TypeError("geometry must be StrictInnerTransportRadialGeometry")

    parent_blob = _authenticate_parent_source()
    radii = np.asarray(geometry.radii, dtype=float)
    z = np.full_like(radii, float(geometry.axial_z))
    t = np.full_like(radii, float(geometry.time))
    witness = materialize_current_partial_nonlinear_mean_attribution(
        backend, radii, z, t
    )
    if witness.backend_kind != "exact-current-a2-pr970-partial-nonlinear-attribution":
        raise RuntimeError("current nonlinear mean backend identity drifted")
    return _materialize_current_mean_witness(
        witness,
        geometry,
        parent_source_blob=parent_blob,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_partial_nonlinear_radial_stress)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "spatial_step",
        "derivative_step", "angular_order", "time_step", "viscosity", "nu",
        "correction", "stage_budget", "normalized_score", "scientific_threshold",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "radial_routing_agent3_pr": RADIAL_ROUTING_AGENT3_PR,
        "radial_routing_agent3_head": RADIAL_ROUTING_AGENT3_HEAD,
        "radial_operator_agent3_pr": RADIAL_OPERATOR_AGENT3_PR,
        "radial_operator_agent3_head": RADIAL_OPERATOR_AGENT3_HEAD,
        "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
        "current_partial_leading_plus_oscillation_consumed": True,
        "current_partial_mixed_nonlinear_mean_consumed": True,
        "current_partial_quadratic_nonlinear_mean_consumed": True,
        "current_partial_aggregate_nonlinear_mean_consumed": True,
        "current_partial_nonlinear_radial_inverse_performed": True,
        "current_partial_quadratic_radial_stress_materialized": True,
        "current_partial_mixed_radial_stress_materialized": True,
        "current_partial_aggregate_radial_stress_materialized": True,
        "radial_inverse_rule": "polynomial_first_cell_compact_moment_complement",
        "theta_exponent": 2,
        "axial_exponent": 1,
        "radial_component_recorded_but_not_inverted": True,
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
        "partial_domain_through_current_Xh_only": True,
        "strict_inner_only": False,
        "velocity_beyond_Xh_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "pressure_gradient_included": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_included": False,
        "complete_ns_defect": False,
        "scoped_current_partial_radial_stress_authorized_as_correction_target": False,
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
