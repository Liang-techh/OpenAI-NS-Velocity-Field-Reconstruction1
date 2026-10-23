from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import constrained_st051b_temporal_material_path as m


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


def test_source_identity_and_truth_boundary_are_fail_closed():
    assert m.TASK_ID == "CR-A9-056"
    assert m.BASE_MAIN_SHA == "f0193d66c9d92948b4820ebcb70263673995b324"
    assert m.ST051_PARENT_HEAD == "4b784f1b8457af2ead49295631d834d4e882000b"
    assert m.AGENT7_SOURCE_HEAD == "f2828ae259eebc867162ef5842a99f0860053229"
    assert m.AGENT9_PATH_SOURCE_HEAD == "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
    assert m.TEMPORAL_GAMMA == pytest.approx(0.025)
    assert m.TEMPORAL_GAINS == pytest.approx((0.025, 0.0375, 0.05))
    assert m.FROZEN_SOURCE_ALPHA == pytest.approx(2.520520814687742)
    assert m.SCIPY_SOURCE["license"] == "BSD-3-Clause"
    assert m.SCIPY_SOURCE["classification"].startswith("direct migration")
    for key in (
        "canonical_velocity_changed",
        "production_candidate_selected",
        "production_static_gain_selected",
        "production_temporal_coefficient_selected",
        "redistribution_rebalanced_on_st051b",
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


def test_seed_radius_summary_requires_frozen_population():
    out = m._summarize_by_seed_radius(_measurement())
    assert set(out) == {"0.6", "0.9", "1.2"}
    for radius in m.SEED_RADII:
        row = out[f"{radius:.1f}"]
        assert row["path_count"] == 16
        assert row["pair_count"] == 8
        assert row["mean_radius_change"] == pytest.approx(-0.1 * radius)
        assert row["mean_pair_axial_separation_change"] == pytest.approx(0.02 * radius)

    broken = _measurement()
    broken["pair_rows"].pop()
    with pytest.raises(ValueError, match="population drift"):
        m._summarize_by_seed_radius(broken)


def test_comparison_is_descriptive_and_directional():
    low = {
        "mean_radius_change": -0.08,
        "mean_absolute_turns": 0.02,
        "maximum_absolute_turns": 0.03,
        "mean_pair_axial_separation_change": 0.07,
        "inward_path_count": 48,
        "pair_axial_separation_growth_count": 16,
        "pair_axial_separation_shrink_count": 8,
    }
    temporal = dict(low)
    temporal.update(
        mean_radius_change=-0.0804,
        mean_absolute_turns=0.0202,
        maximum_absolute_turns=0.0303,
        mean_pair_axial_separation_change=0.0702,
    )
    got = m._compare(temporal, low)
    assert got["mean_absolute_turns_relative_change"] == pytest.approx(0.01)
    assert got["maximum_absolute_turns_relative_change"] == pytest.approx(0.01)
    assert got["radial_contraction_magnitude_relative_change"] > 0
    assert got["mean_pair_separation_change_relative_change"] > 0
    assert got["inward_path_count_delta"] == 0


def test_by_radius_comparison_preserves_registered_bands():
    old = m._summarize_by_seed_radius(_measurement())
    new = copy.deepcopy(old)
    for row in new.values():
        row["mean_absolute_turns"] *= 1.02
        row["maximum_absolute_turns"] *= 1.03
        row["mean_radius_change"] *= 1.01
        row["mean_pair_axial_separation_change"] *= 0.99
    got = m._compare_by_radius(new, old)
    assert set(got) == {"0.6", "0.9", "1.2"}
    for row in got.values():
        assert row["mean_absolute_turns_relative_change"] == pytest.approx(0.02)
        assert row["maximum_absolute_turns_relative_change"] == pytest.approx(0.03)
        assert row["radial_contraction_magnitude_relative_change"] == pytest.approx(0.01)
        assert row["mean_pair_separation_change_relative_change"] == pytest.approx(-0.01)


def test_truth_audit_rejects_promotion_and_schedule_drift():
    path = {"path_count": 48, "paired_material_line_count": 24}
    report = {
        "task_id": m.TASK_ID,
        "base_main_sha": m.BASE_MAIN_SHA,
        "static_low_material_paths": dict(path),
        "static_high_material_paths": dict(path),
        "temporal_material_paths": dict(path),
        "redistribution": {
            "static_low_gain": m.STATIC_LOW_GAIN,
            "static_high_gain": m.STATIC_HIGH_GAIN,
            "temporal_gamma": m.TEMPORAL_GAMMA,
            "temporal_gains_at_registered_times": list(m.TEMPORAL_GAINS),
            "gamma_values_outside_preregistered_grid_evaluated": False,
            "transform_rebalanced_on_st051b": False,
        },
        "truth_boundary": copy.deepcopy(m.TRUTH_BOUNDARY),
    }
    m._audit_truth(report)
    report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary drift"):
        m._audit_truth(report)

    report["truth_boundary"]["pde_validated"] = False
    report["redistribution"]["temporal_gains_at_registered_times"][-1] = 0.06
    with pytest.raises(ValueError, match="temporal schedule identity drift"):
        m._audit_truth(report)


def test_static_reference_controls_are_fixed():
    assert m.STATIC_CONTROL_REFERENCE["0.025"]["mean_absolute_turns"] == pytest.approx(
        0.02229297191872165
    )
    assert m.STATIC_CONTROL_REFERENCE["0.05"]["mean_absolute_turns"] == pytest.approx(
        0.02249374673546196
    )
    assert m._relative_change(1.1, 1.0) == pytest.approx(0.1)
    with pytest.raises(ValueError, match="zero comparison denominator"):
        m._relative_change(1.0, 0.0)
