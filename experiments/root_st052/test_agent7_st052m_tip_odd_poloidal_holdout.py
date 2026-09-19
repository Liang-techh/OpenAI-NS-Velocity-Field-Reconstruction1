import numpy as np

import agent7_st052m_tip_odd_poloidal_holdout as audit


def test_fresh_seed_protocol_is_disjoint_from_development():
    seeds, metadata = audit._fresh_seed_table()
    assert seeds.shape == (64, 3)
    assert len(metadata) == 64
    assert audit._development_seed_overlap_count(seeds) == 0
    assert sorted({m["radius"] for m in metadata}) == list(audit.HOLDOUT_RADII)
    assert sorted({m["z_abs"] for m in metadata}) == [0.75, 1.25]
    assert {m["band"] for m in metadata} == {"shoulder", "tip"}


def test_fresh_numerical_protocol_does_not_reuse_selection_protocol():
    base = audit.screen.loc.parent.base
    assert audit.OUTPUT_SAMPLES == 37
    assert base.OUTPUT_SAMPLES == 33
    assert audit.SOLVER_MAX_STEP == 0.008
    assert base.SOLVER_MAX_STEP == 0.01
    assert audit.SOLVER_RTOL == 5.0e-10
    assert base.SOLVER_RTOL == 1.0e-9
    assert audit.CHECKPOINT_INDICES == (0, 9, 18, 27, 36)
    times = np.linspace(*audit.TIME_INTERVAL, audit.OUTPUT_SAMPLES)
    for index, expected in zip(audit.CHECKPOINT_INDICES, audit.CHECKPOINT_TIMES):
        assert abs(float(times[index]) - expected) < 1.0e-15


def _record(tip_dr, tip_ur, tip_count, shoulder_dr, shoulder_count, turns):
    return {
        "all": {
            "path_count": 64,
            "inward_path_count": 32,
            "mean_delta_r": -0.01,
            "mean_absolute_turns": float(turns),
            "maximum_absolute_turns": float(turns),
        },
        "shoulder": {
            "path_count": 32,
            "inward_path_count": int(shoulder_count),
            "mean_delta_r": float(shoulder_dr),
            "min_delta_r": float(shoulder_dr),
            "max_delta_r": float(shoulder_dr),
            "mean_start_radial_velocity": -0.1,
            "mean_absolute_turns": float(turns),
            "maximum_absolute_turns": float(turns),
        },
        "tip": {
            "path_count": 32,
            "inward_path_count": int(tip_count),
            "mean_delta_r": float(tip_dr),
            "min_delta_r": float(tip_dr),
            "max_delta_r": float(tip_dr),
            "mean_start_radial_velocity": float(tip_ur),
            "mean_absolute_turns": float(turns),
            "maximum_absolute_turns": float(turns),
        },
    }


def _records(row):
    return {"segments": {f"s{i}": row for i in range(4)}, "whole_interval": row}


def test_directional_gate_requires_actual_inward_tip_mean():
    control = _records(_record(-0.01, -0.1, 32, -0.01, 32, 1.0))
    linear = _records(_record(0.02, 0.2, 0, -0.02, 32, 0.95))
    child = _records(_record(0.01, 0.1, 0, -0.02, 32, 0.96))

    good = _records(_record(-0.002, -0.02, 8, -0.02, 32, 0.98))
    decision = audit.directional_decision(control, linear, child, good)
    assert decision["fresh_tip_path_directional_preferred"] is True

    still_outward = _records(_record(0.001, -0.02, 8, -0.02, 32, 0.98))
    rejected = audit.directional_decision(control, linear, child, still_outward)
    assert rejected["fresh_tip_path_directional_preferred"] is False
    assert all(
        row["candidate_tip_mean_inward"] is False
        for row in rejected["segments"].values()
    )


def test_truth_boundary_freezes_no_retuning_or_promotion():
    assert audit.PREREG_ISSUE == 713
    assert audit.SOURCE_PARENT_PR == 701
    assert audit.SOURCE_PARENT_HEAD == "e3ecc8654f2c9a2324265bbd01df596b4e503404"
    assert audit.TRUTH["source_701_candidate_retuned"] is False
    assert audit.TRUTH["heldout_evaluation_feedback_into_alpha"] is False
    assert audit.TRUTH["new_spatial_basis_added"] is False
    assert audit.TRUTH["new_temporal_basis_added"] is False
    assert audit.TRUTH["parameter_grid_scan_performed"] is False
    assert audit.TRUTH["optimization_performed"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["openai_field_identified"] is False
