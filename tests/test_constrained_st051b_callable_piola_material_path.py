from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.constrained_st051b_callable_piola_material_path import (
    AGENT7_HEAD,
    AGENT7_TASK_ID,
    AGENT9_PATH_HEAD,
    BASE_MAIN_SHA,
    PIOLA_BETA,
    REDISTRIBUTION_GAIN,
    TASK_ID,
    TRUTH_BOUNDARY,
    by_seed_radius,
    comparison,
)


def _fake_measurement(scale: float = 1.0):
    per_path = []
    for radius in (0.6, 0.9, 1.2):
        for angle in range(8):
            for z in (-0.3, 0.3):
                per_path.append(
                    {
                        "seed": {"radius": radius, "angle_index": angle, "z": z},
                        "absolute_turns": scale * (0.01 + 0.002 * radius),
                        "radius_change": -scale * (0.05 + 0.01 * radius),
                        "abs_z_change": scale * 0.02,
                    }
                )
    return {
        "mean_radius_change": -0.08 * scale,
        "mean_absolute_turns": 0.02 * scale,
        "maximum_absolute_turns": 0.03 * scale,
        "mean_pair_axial_separation_change": 0.06 * scale,
        "inward_path_count": 48,
        "pair_axial_separation_growth_count": 16,
        "pair_axial_separation_shrink_count": 8,
        "per_path": per_path,
    }


def test_frozen_identities_and_truth_boundary():
    assert TASK_ID == "CR-A9-057"
    assert BASE_MAIN_SHA == "f0193d66c9d92948b4820ebcb70263673995b324"
    assert AGENT7_HEAD == "c87ffa3f1a798df0429f16f5ee12dfb41d6c8ded"
    assert AGENT7_TASK_ID == "CR003-ST051B-CALLABLE-PIOLA-TRANSFER-076"
    assert AGENT9_PATH_HEAD == "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
    assert REDISTRIBUTION_GAIN == 0.025
    assert PIOLA_BETA == 0.075
    for key in (
        "canonical_velocity_changed",
        "production_beta_selected",
        "production_candidate_selected",
        "held_out_pde_residual_evaluated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert TRUTH_BOUNDARY[key] is False


def test_seed_radius_summary_preserves_frozen_population():
    rows = by_seed_radius(_fake_measurement())
    assert set(rows) == {"0.6", "0.9", "1.2"}
    assert all(row["path_count"] == 16 for row in rows.values())
    assert rows["0.6"]["mean_absolute_turns"] == pytest.approx(0.0112)
    assert rows["1.2"]["mean_radius_change"] == pytest.approx(-0.062)


def test_seed_radius_summary_fails_closed_on_population_drift():
    measurement = _fake_measurement()
    measurement["per_path"].pop()
    with pytest.raises(ValueError, match="population"):
        by_seed_radius(measurement)


def test_comparison_uses_contraction_magnitude_not_signed_radius_ratio():
    base = _fake_measurement(1.0)
    child = _fake_measurement(1.1)
    comp = comparison(child, base)
    assert comp["mean_absolute_turns_relative"] == pytest.approx(0.1)
    assert comp["maximum_absolute_turns_relative"] == pytest.approx(0.1)
    assert comp["radial_contraction_magnitude_relative"] == pytest.approx(0.1)
    assert comp["mean_pair_axial_separation_change_relative"] == pytest.approx(0.1)
    assert comp["inward_path_count_delta"] == 0
    assert comp["pair_growth_count_delta"] == 0
    assert comp["pair_shrink_count_delta"] == 0


def test_truth_boundary_is_json_serializable():
    encoded = json.dumps(TRUTH_BOUNDARY, sort_keys=True)
    assert "pde_validated" in encoded
