from __future__ import annotations

import numpy as np

from agent7_morphology_transfer import CHILD_IDS, PARENT_ID, _weighted_quantile, audit_morphology_transfer


def test_weighted_quantile_is_weight_sensitive():
    values = np.array([0.0, 1.0, 2.0])
    weights = np.array([1.0, 1.0, 8.0])
    assert _weighted_quantile(values, weights, 0.50) == 2.0
    assert _weighted_quantile(values, weights, 0.90) == 2.0


def test_real_frozen_transfer_smoke_is_fail_closed():
    report = audit_morphology_transfer(times=(0.50,), grid_sizes=(17, 19))
    assert report["parent_id"] == PARENT_ID
    assert tuple(report["child_ids"]) == CHILD_IDS
    assert report["contract"]["fit_or_parameter_selection"] is False
    assert report["contract"]["clean_transfer_rule_reused_from_pr_378"] is True
    assert report["upstream_residual_receipt"]["classification"].startswith("frozen upstream")
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["held_out_pde_residual_recomputed"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    for ident in (PARENT_ID,) + CHILD_IDS:
        assert report["structure"][ident]["outside_support_max_abs_velocity"] == 0.0
        assert report["structure"][ident]["representative_inward_swirl_bipolar_signs_all_pass"]
