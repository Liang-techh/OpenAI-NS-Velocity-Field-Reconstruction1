from __future__ import annotations

import math

import pytest

from openai_ns_reconstruction import constrained_st052_material_path as m


def _measurement(turn_scale: float = 1.0, radius_scale: float = 1.0) -> dict:
    rows = []
    for radius in m.SEED_RADII:
        for z in (-0.3, 0.3):
            for angle in range(8):
                rows.append(
                    {
                        "seed": {"radius": radius, "z": z, "angle_index": angle},
                        "absolute_turns": turn_scale * radius,
                        "radius_change": -radius_scale * 0.1 * radius,
                        "abs_z_change": 0.01 * radius,
                    }
                )
    return {
        "per_path": rows,
        "mean_absolute_turns": turn_scale * 0.9,
        "maximum_absolute_turns": turn_scale * 1.2,
        "mean_radius_change": -radius_scale * 0.09,
        "mean_pair_axial_separation_change": turn_scale * 0.07,
        "inward_path_count": 48,
        "pair_axial_separation_growth_count": 16,
        "pair_axial_separation_shrink_count": 8,
    }


def test_frozen_source_and_truth_contract() -> None:
    assert m.TASK_ID == "CR-A9-058"
    assert m.BASE_MAIN_SHA == "f0193d66c9d92948b4820ebcb70263673995b324"
    assert m.ST052_PR == 508
    assert m.ST052_HEAD == "b3b8bfdbe1077f9ec967d158602951997d81e17d"
    assert m.AGENT9_PATH_PR == 435
    assert m.AGENT9_PATH_HEAD == "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
    assert m.SEED_RADII == (0.6, 0.9, 1.2)
    assert m.SCIPY_SOURCE["classification"].startswith("direct migration")
    assert m.SCIPY_SOURCE["copied_upstream_implementation"] is False
    assert all(
        m.TRUTH_BOUNDARY[key] is False
        for key in (
            "canonical_velocity_changed",
            "production_candidate_selected",
            "held_out_pde_residual_recomputed_in_this_increment",
            "parent_pde_receipt_transferred",
            "visualization_ready",
            "visual_correspondence_verified",
            "pde_validated",
            "source_correspondence_verified",
            "paper_exact",
            "openai_field_identified",
            "blowup_proved",
        )
    )
    assert m.TRUTH_BOUNDARY["comparison_is_descriptive_not_acceptance"] is True


def test_by_seed_radius_uses_exact_frozen_population() -> None:
    out = m.by_seed_radius(_measurement()["per_path"] and _measurement())
    assert set(out) == {"0.6", "0.9", "1.2"}
    assert all(row["path_count"] == 16 for row in out.values())
    assert out["0.6"]["mean_absolute_turns"] == pytest.approx(0.6)
    assert out["1.2"]["mean_radius_change"] == pytest.approx(-0.12)


def test_by_seed_radius_fails_closed_on_seed_drift() -> None:
    bad = _measurement()
    bad["per_path"] = bad["per_path"][:-1]
    with pytest.raises(ValueError, match="unexpected frozen seed population"):
        m.by_seed_radius(bad)


def test_comparison_reports_signed_changes_without_acceptance_rule() -> None:
    parent = _measurement(turn_scale=1.0, radius_scale=1.0)
    child = _measurement(turn_scale=1.1, radius_scale=0.95)
    out = m.comparison(child, parent)
    assert out["mean_absolute_turns_relative"] == pytest.approx(0.1)
    assert out["maximum_absolute_turns_relative"] == pytest.approx(0.1)
    assert out["radial_contraction_magnitude_relative"] == pytest.approx(-0.05)
    assert out["mean_pair_axial_separation_change_relative"] == pytest.approx(0.1)
    assert out["inward_path_count_delta"] == 0
    assert out["pair_growth_count_delta"] == 0
    assert out["pair_shrink_count_delta"] == 0
    assert math.isfinite(out["mean_absolute_turns_relative"])
