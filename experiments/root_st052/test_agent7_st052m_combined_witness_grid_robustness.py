from __future__ import annotations

import math

import agent7_st052m_combined_witness_grid_robustness as audit


def _record(tip: float, aspect: float, central: float = 0.0) -> dict:
    return {
        "increments_child_vs_linear": {
            "tip_thinning_increment": tip,
            "aspect_gain_increment": aspect,
            "central_radial_rms_increment": central,
        }
    }


def test_frozen_source_identity_and_no_retuning_contract():
    assert audit.PREREG_ISSUE == 658
    assert audit.SOURCE_WITNESS_PR == 652
    assert audit.SOURCE_WITNESS_HEAD == "b0d98d0d5fdb3ca6a438fdf6ac0dbaaef7781faf"
    assert audit.GRID_RESOLUTIONS == (33, 41, 49)
    assert audit.MID_TIME == 0.50
    assert audit.PARETO_TOL == 5.0e-5
    assert audit.MARGIN_RESOLUTION_FRACTION == 0.25
    assert audit.witness.SWIRL_A == 0.065899695471146
    assert math.isclose(audit.witness.SHOULDER_LAMBDA, -0.02, rel_tol=0.0, abs_tol=1e-15)
    assert audit.TRUTH["witness_retuned"] is False
    assert audit.TRUTH["new_spatial_basis_added"] is False
    assert audit.TRUTH["new_temporal_basis_added"] is False
    assert audit.TRUTH["optimization_performed"] is False
    assert audit.TRUTH["parameter_grid_scan_performed"] is False


def test_resolution_claim_requires_drift_small_relative_to_existing_margin():
    # 41-grid miss is 1e-6 beyond the frozen boundary; 49-grid drift is only
    # 2e-7, within the preregistered 25% margin-resolution budget.
    records = {
        33: _record(-5.08e-5, 2.0e-5),
        41: _record(-5.10e-5, 2.1e-5),
        49: _record(-5.12e-5, 2.2e-5),
    }
    d = audit.robustness_decision(records)
    assert d["both_41_and_49_below_boundary"] is True
    assert d["drift_small_enough_to_resolve_41_margin"] is True
    assert d["aspect_increment_same_sign_all_grids"] is True
    assert d["pareto_rejection_grid_resolved"] is True


def test_large_grid_drift_keeps_tolerance_edge_rejection_unresolved():
    records = {
        33: _record(-4.7e-5, 2.0e-5),
        41: _record(-5.006e-5, 2.1e-5),
        49: _record(-4.95e-5, 2.2e-5),
    }
    d = audit.robustness_decision(records)
    assert d["pareto_rejection_grid_resolved"] is False
    assert d["drift_small_enough_to_resolve_41_margin"] is False


def test_nonfinite_or_aspect_sign_flip_cannot_be_called_resolved():
    sign_flip = {
        33: _record(-5.2e-5, 2.0e-5),
        41: _record(-5.2e-5, 2.0e-5),
        49: _record(-5.2e-5, -1.0e-6),
    }
    assert audit.robustness_decision(sign_flip)["pareto_rejection_grid_resolved"] is False

    nonfinite = {
        33: _record(-5.2e-5, 2.0e-5),
        41: _record(-5.2e-5, 2.0e-5),
        49: _record(float("nan"), 2.0e-5),
    }
    assert audit.robustness_decision(nonfinite)["pareto_rejection_grid_resolved"] is False
