from __future__ import annotations

from copy import deepcopy
import math

import numpy as np

from openai_ns_reconstruction.kokuno_a4_current_tsh_independent_audit import (
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    GRID_SIZES,
    MAX_ENVELOPE_STABILITY,
    OFFGRID_COUNT,
    SEED,
    _assess_report_from_observation,
    independent_ell_observation,
    materialize_current_tsh_independent_audit,
)


def _manufactured_report(*, B0: float, X_i: float, log_X_R: float) -> dict:
    lower = 160.0 * (float(B0) + math.log(2.0))
    selected_T = math.nextafter(lower, math.inf)
    log_x_i = math.log(float(X_i)) - float(log_X_R)
    max_T = -8.0 - log_x_i
    margin = max_T - selected_T
    feasible = selected_T < max_T
    return {
        "selected_B0": {
            "selected_numerical_B0_envelope": float(B0),
            "disjoint_holdout_within_envelope": True,
            "analytic_source_B0_bound_proved": False,
        },
        "geometry": {
            "T_sh_lower_bound_from_selected_numerical_B0": lower,
            "selected_T_sh": selected_T,
            "log_X_R": float(log_X_R),
            "log_x_i": log_x_i,
            "max_T_sh_for_current_outer_schedule": max_T,
            "selected_log_x_sep": log_x_i + selected_T,
            "geometry_margin": margin,
            "separation_geometry_feasible": feasible,
        },
        "selected_route_ready": feasible,
    }


def test_frozen_protocol_and_repository_gates_are_not_relaxed() -> None:
    assert SEED == 9173641
    assert GRID_SIZES == (129, 257, 513)
    assert OFFGRID_COUNT == 257
    assert MAX_ENVELOPE_STABILITY == 2.0e-3
    assert FINAL_MOMENTUM_GATE == 1.0e-3
    assert FINAL_DIVERGENCE_GATE == 1.0e-5


def test_uniform_sampler_is_nested_and_locks_a_manufactured_endpoint_maximum() -> None:
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assert observation["nested_uniform_grids_exact"] is True
    assert observation["nested_max_monotone"] is True
    assert observation["medium_to_fine_scale_normalized_max_change"] == 0.0
    assert [level["eta_points"] for level in observation["levels"]] == [129, 257, 513]
    assert all(level["max_abs_ell_i"] == 1.25 for level in observation["levels"])
    assert observation["combined_observed_max_abs_ell_i"] == 1.25
    assert abs(observation["combined_worst_eta"]) == 1.0


def test_independent_public_algebra_accepts_a_consistent_manufactured_report() -> None:
    X_i = 110.0
    log_X_R = 400.0
    report = _manufactured_report(B0=1.30, X_i=X_i, log_X_R=log_X_R)
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assessed = _assess_report_from_observation(
        report, observation, X_i=X_i, log_X_R=log_X_R
    )
    assert assessed["audit_pass"] is True
    assert all(assessed["checks"].values())
    assert assessed["independent"]["observed_max_abs_ell_i"] == 1.25


def test_understated_B0_is_rejected_by_independent_public_samples() -> None:
    X_i = 110.0
    log_X_R = 400.0
    report = _manufactured_report(B0=1.20, X_i=X_i, log_X_R=log_X_R)
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assessed = _assess_report_from_observation(
        report, observation, X_i=X_i, log_X_R=log_X_R
    )
    assert assessed["checks"]["selected_B0_covers_independent_probes"] is False
    assert assessed["audit_pass"] is False


def test_T_sh_formula_drift_is_rejected() -> None:
    X_i = 110.0
    log_X_R = 400.0
    report = _manufactured_report(B0=1.30, X_i=X_i, log_X_R=log_X_R)
    report["geometry"]["T_sh_lower_bound_from_selected_numerical_B0"] += 1.0e-3
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assessed = _assess_report_from_observation(
        report, observation, X_i=X_i, log_X_R=log_X_R
    )
    assert assessed["checks"]["T_sh_lower_bound_formula_matches"] is False
    assert assessed["audit_pass"] is False


def test_outer_scale_identity_drift_is_rejected() -> None:
    X_i = 110.0
    log_X_R = 400.0
    report = _manufactured_report(B0=1.30, X_i=X_i, log_X_R=log_X_R)
    report["geometry"]["log_X_R"] += 1.0e-3
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assessed = _assess_report_from_observation(
        report, observation, X_i=X_i, log_X_R=log_X_R
    )
    assert assessed["checks"]["log_X_R_identity_matches"] is False
    assert assessed["audit_pass"] is False


def test_feasibility_boolean_laundering_is_rejected_even_if_algebra_is_unchanged() -> None:
    X_i = 110.0
    log_X_R = 400.0
    report = _manufactured_report(B0=1.30, X_i=X_i, log_X_R=log_X_R)
    report["geometry"]["separation_geometry_feasible"] = not report["geometry"][
        "separation_geometry_feasible"
    ]
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assessed = _assess_report_from_observation(
        report, observation, X_i=X_i, log_X_R=log_X_R
    )
    assert assessed["checks"]["feasibility_boolean_matches"] is False
    assert assessed["audit_pass"] is False


def test_negative_route_geometry_can_be_a_valid_consistency_audit() -> None:
    # A deliberately small outer scale makes the candidate-side route infeasible.
    # Correctly reporting that negative state is not an audit failure.
    X_i = 110.0
    log_X_R = math.log(X_i) + 20.0
    report = _manufactured_report(B0=1.30, X_i=X_i, log_X_R=log_X_R)
    assert report["geometry"]["separation_geometry_feasible"] is False
    observation = independent_ell_observation(lambda eta: 1.0 + 0.25 * eta * eta)
    assessed = _assess_report_from_observation(
        report, observation, X_i=X_i, log_X_R=log_X_R
    )
    assert assessed["audit_pass"] is True
    assert assessed["candidate_side_route_feasible"] is False


def test_real_artifact_receipt_keeps_full_ns_promotions_fail_closed() -> None:
    receipt = materialize_current_tsh_independent_audit()
    assert receipt["save_reload_semantic_identity_exact"] is True
    assert isinstance(receipt["audit_pass"], bool)
    assert isinstance(receipt["candidate_side_route_feasible"], bool)
    assert receipt["source_B0_analytic_bound_proved"] is False
    assert receipt["source_T_sh_lower_bound_verified"] is False
    assert receipt["leading_only_ns_residual_assessed"] is False
    assert receipt["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert receipt["after_correction_ns_residual_assessed"] is False
    assert receipt["same_protocol_comparable_to_st006"] is False
    assert receipt["pde_validated"] is False
    assert set(receipt["mutation_firewall"]) == {
        "understated_B0_rejected",
        "T_sh_formula_drift_rejected",
        "outer_scale_log_X_R_drift_rejected",
        "feasibility_boolean_flip_rejected",
        "public_ell_amplitude_mutation_rejected",
    }
    assert np.isfinite(
        receipt["independent_ell_observation"]["combined_observed_max_abs_ell_i"]
    )
