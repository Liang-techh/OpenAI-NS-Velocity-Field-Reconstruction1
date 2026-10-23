from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction import constrained_st050rc_energy_neutral_redist_material_path as m


def _measurement():
    per_path = []
    pair_rows = []
    idx = 0
    for radius in m.SEED_RADII:
        for angle in range(8):
            for z in (-0.3, 0.3):
                per_path.append(
                    {
                        "path_index": idx,
                        "seed": {"radius": radius, "angle_index": angle, "z": z},
                        "absolute_turns": radius + 0.01 * angle,
                        "radius_change": -0.1 * radius,
                    }
                )
                idx += 1
            pair_rows.append(
                {
                    "radius": radius,
                    "angle_index": angle,
                    "axial_separation_change": 0.02 * radius,
                }
            )
    return {"per_path": per_path, "pair_rows": pair_rows}


def test_source_and_truth_boundary_are_fail_closed():
    assert m.TASK_ID == "CR-A9-051"
    assert m.BASE_MAIN_SHA == "f0193d66c9d92948b4820ebcb70263673995b324"
    assert m.REDISTRIBUTION_GAIN == pytest.approx(0.025)
    assert m.SCIPY_SOURCE["license"] == "BSD-3-Clause"
    assert m.SCIPY_SOURCE["classification"].startswith("direct migration")
    for key in (
        "canonical_velocity_changed",
        "production_candidate_selected",
        "production_redistribution_gain_selected",
        "pressure_or_force_changed",
        "held_out_pde_residual_evaluated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert m.TRUTH_BOUNDARY[key] is False
    assert m.TRUTH_BOUNDARY["comparison_is_descriptive_not_acceptance"] is True


def test_seed_radius_summary_requires_the_frozen_population():
    out = m._summarize_by_seed_radius(_measurement())
    assert set(out) == {"0.6", "0.9", "1.2"}
    for radius in m.SEED_RADII:
        row = out[f"{radius:.1f}"]
        assert row["path_count"] == 16
        assert row["pair_count"] == 8
        assert row["mean_radius_change"] == pytest.approx(-0.1 * radius)
        assert row["mean_pair_axial_separation_change"] == pytest.approx(0.02 * radius)

    broken = _measurement()
    broken["per_path"].pop()
    with pytest.raises(ValueError, match="population drift"):
        m._summarize_by_seed_radius(broken)


def test_comparison_is_directional_and_does_not_define_acceptance():
    parent = {
        "mean_radius_change": -0.08,
        "mean_absolute_turns": 0.02,
        "maximum_absolute_turns": 0.03,
        "mean_pair_axial_separation_change": 0.07,
        "inward_path_count": 48,
        "pair_axial_separation_growth_count": 16,
        "pair_axial_separation_shrink_count": 8,
    }
    child = dict(parent)
    child.update(
        mean_radius_change=-0.079,
        mean_absolute_turns=0.021,
        maximum_absolute_turns=0.031,
        mean_pair_axial_separation_change=0.069,
    )
    got = m._compare(child, parent)
    assert got["mean_absolute_turns_relative_change"] == pytest.approx(0.05)
    assert got["maximum_absolute_turns_relative_change"] == pytest.approx(1 / 30)
    assert got["radial_contraction_magnitude_relative_change"] < 0
    assert got["mean_pair_separation_change_relative_change"] < 0
    assert got["inward_path_count_delta"] == 0


def test_truth_audit_rejects_promotion():
    report = {
        "task_id": m.TASK_ID,
        "base_main_sha": m.BASE_MAIN_SHA,
        "frozen_material_path_contract": {"path_count": 48},
        "parent_material_paths": {"path_count": 48, "paired_material_line_count": 24},
        "child_material_paths": {"path_count": 48, "paired_material_line_count": 24},
        "truth_boundary": copy.deepcopy(m.TRUTH_BOUNDARY),
    }
    m._audit_truth(report)
    report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary drift"):
        m._audit_truth(report)


def test_relative_change_rejects_zero_reference():
    assert m._relative_change(1.1, 1.0) == pytest.approx(0.1)
    with pytest.raises(ValueError, match="zero comparison denominator"):
        m._relative_change(1.0, 0.0)
