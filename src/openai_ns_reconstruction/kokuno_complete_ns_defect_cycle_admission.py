"""Fail-closed admission gate before any real Kokuno finite correction cycle.

This module deliberately does *not* accept residual/defect arrays, pressure,
forcing, correction fields, gains, sample sets, or thresholds from callers.
Instead it authenticates the exact current A3 #1045 correction-side source and
turns its truth boundary into an executable blocker receipt.  A later adapter
may only replace this blocker when the repository has one identity-bound global
candidate carrying a matched pressure, preregistered restricted forcing, a
complete NS defect, disjoint held-in/held-out evidence, an independent second
covariance column with bounded-inverse preflight, and an accepted correction-
cycle gain receipt.

The finite-stage engineering iteration is
    u^(k+1) = u^k + delta_u_k.
The corrected 2026-09-09 Kokuno reader is structural provenance only; this gate
is not independent PDE validation and does not run an iteration.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path

from .kokuno_current_rf40_power_law_nonlinear_mean_attribution import _git_blob_sha1
from .kokuno_current_rf40_power_law_nonlinear_radial_force import (
    materialize_current_rf40_power_law_nonlinear_radial_force,
    truth_boundary as parent_truth_boundary,
)
from .kokuno_oscillatory_nonlinear_mean_attribution import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
)

TASK = "KOKUNO-A3-COMPLETE-NS-DEFECT-CYCLE-ADMISSION-115"
SCHEMA = "kokuno-a3-complete-ns-defect-cycle-admission-v1"
PARENT_AGENT3_PR = 1045
PARENT_AGENT3_HEAD = "3ec20a77b551be819a71e3308e9ecab2a2226406"
PARENT_AGENT3_SOURCE_BLOB = "39c98d3f4df05846848b897c240b3b56fb6fbe2c"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_READER_DATE = "2026-09-09"
FINITE_CYCLE_FORMULA = "u^(k+1)=u^k+delta_u_k"
GAIN_GUARD_PRIOR_A3_PR = 293
SECOND_COLUMN_COVERAGE_PRIOR_A3_PR = 282
SECOND_COLUMN_BOUNDED_INVERSE_PRIOR_A3_PR = 302

REQUIRED_EVIDENCE = (
    "one_global_candidate_semantic_identity",
    "outer_global_velocity_materialized",
    "matched_global_pressure_materialized",
    "preregistered_restricted_forcing_materialized",
    "forcing_proven_not_residual_defined",
    "complete_identity_bound_ns_defect_materialized",
    "held_in_held_out_sample_sets_disjoint",
    "independent_second_covariance_column_passed",
    "second_column_bounded_inverse_preflight_passed",
    "cartesian_correction_velocity_materialized",
    "correction_cycle_gain_guard_accepted",
)


def _authenticate_parent_source() -> str:
    source = inspect.getsourcefile(materialize_current_rf40_power_law_nonlinear_radial_force)
    if source is None:
        raise RuntimeError("A3 #1045 radial-force source is unavailable")
    blob = _git_blob_sha1(Path(source))
    if blob != PARENT_AGENT3_SOURCE_BLOB:
        raise RuntimeError("A3 #1045 radial-force source blob drifted")
    return blob


@dataclass(frozen=True)
class CompleteNSDefectCycleAdmissionWitness:
    parent_source_blob: str
    parent_truth: dict[str, object]
    missing_required_evidence: tuple[str, ...]
    finite_cycle_admitted: bool

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "finite_cycle_formula": FINITE_CYCLE_FORMULA,
            "finite_cycle_admitted": self.finite_cycle_admitted,
            "missing_required_evidence": list(self.missing_required_evidence),
            "required_evidence": list(REQUIRED_EVIDENCE),
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "parent_agent3_source_blob": self.parent_source_blob,
            "prior_a3_gain_guard_pr": GAIN_GUARD_PRIOR_A3_PR,
            "prior_a3_second_column_coverage_pr": SECOND_COLUMN_COVERAGE_PRIOR_A3_PR,
            "prior_a3_second_column_bounded_inverse_pr": SECOND_COLUMN_BOUNDED_INVERSE_PRIOR_A3_PR,
            "source_reader": {
                "repository": SOURCE_READER_REPO,
                "head": SOURCE_READER_HEAD,
                "path": SOURCE_READER_PATH,
                "blob": SOURCE_READER_BLOB,
                "date": SOURCE_READER_DATE,
                "role": "structural-source-not-independent-validation",
            },
            "fixed_scientific_gates": {
                "normalized_momentum": FINAL_NORMALIZED_MOMENTUM_GATE,
                "normalized_divergence": FINAL_NORMALIZED_DIVERGENCE_GATE,
            },
            "parent_truth_boundary": dict(self.parent_truth),
            "truth_boundary": truth_boundary(),
        }


def materialize_current_complete_ns_defect_cycle_admission() -> CompleteNSDefectCycleAdmissionWitness:
    """Authenticate #1045 and fail closed on every missing cycle prerequisite."""
    parent_blob = _authenticate_parent_source()
    parent = parent_truth_boundary()

    # The exact #1045 lineage is correction-side only.  It has no global
    # velocity/pressure/forcing/complete defect or real Cartesian correction.
    # Second-column and gain-guard receipts are also deliberately not inherited
    # merely because historical A3 experiments exist.
    missing = list(REQUIRED_EVIDENCE)
    invariant_false = (
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_included",
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
    )
    for key in invariant_false:
        if parent.get(key) is not False:
            raise RuntimeError(f"parent truth boundary unexpectedly promoted {key}")

    return CompleteNSDefectCycleAdmissionWitness(
        parent_source_blob=parent_blob,
        parent_truth=dict(parent),
        missing_required_evidence=tuple(missing),
        finite_cycle_admitted=False,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_complete_ns_defect_cycle_admission)
    forbidden = {
        "residual", "defect", "pressure", "forcing", "velocity", "correction",
        "gain", "alpha", "damping", "held_in", "held_out", "samples",
        "scientific_threshold", "momentum_threshold", "divergence_threshold",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "caller_supplied_surrogate_defect_allowed": False,
        "caller_supplied_residual_arrays_allowed": False,
        "caller_supplied_forcing_allowed": False,
        "residual_defined_forcing_allowed": False,
        "same_samples_allowed_for_held_in_and_held_out": False,
        "scoped_mean_stress_force_counts_as_complete_defect": False,
        "historical_gain_guard_counts_as_current_acceptance": False,
        "historical_second_column_preflight_counts_as_current_evidence": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "complete_ns_defect": False,
        "complete_defect_identity_bound_to_global_candidate": False,
        "matched_global_pressure_materialized": False,
        "preregistered_restricted_forcing_materialized": False,
        "forcing_proven_not_residual_defined": False,
        "held_in_held_out_sample_sets_disjoint": False,
        "independent_second_covariance_column_passed": False,
        "second_column_bounded_inverse_preflight_passed": False,
        "current_real_ns_correction_velocity_materialized": False,
        "correction_cycle_gain_guard_accepted": False,
        "finite_correction_cycle_admitted": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
