import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("agent7_st052m_linear_temporal_activation.py")
SPEC = importlib.util.spec_from_file_location("agent7_linear_temporal", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_preregistered_identity_and_truth_boundary():
    assert mod.PREREG_ISSUE == 585
    assert mod.SOURCE_AGENT7_HEAD == "39b106ad8cb8df2064cabead3a12682089575e74"
    assert mod.SOURCE_MORPH_HEAD == "0ddf5f6b321ba78e27aa3ca3e073ef593e55869a"
    assert mod.SOURCE_PATH_PR == 574
    assert mod.SOURCE_RENDER_PR == 583
    assert mod.TAPER_TAU == 0.05
    assert mod.EXPECTED_BETA == 0.08837490297155456
    assert mod.GRID_RESOLUTION == 41
    assert mod.TRUTH["canonical_velocity_changed"] is False
    assert mod.TRUTH["saved_velocity_changed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False


def test_linear_activation_has_frozen_endpoints():
    assert mod.activation(0.25) == 0.0
    assert mod.activation(0.50) == 0.5
    assert mod.activation(0.75) == 1.0
    with np.testing.assert_raises(ValueError):
        mod.activation(0.20)
    with np.testing.assert_raises(ValueError):
        mod.activation(0.80)


def test_offmidplane_seed_table_matches_frozen_agent9_protocol():
    seeds, metadata = mod._seed_table()
    assert seeds.shape == (48, 3)
    assert len(metadata) == 48
    assert set(row["band"] for row in metadata) == {"shoulder", "tip"}
    assert sum(row["band"] == "shoulder" for row in metadata) == 24
    assert sum(row["band"] == "tip" for row in metadata) == 24
    assert set(abs(float(row["z_abs"])) for row in metadata) == {0.85, 1.15}
    assert set(int(row["sign"]) for row in metadata) == {-1, 1}


def test_path_comparison_uses_same_signed_relative_convention():
    control = {
        "mean_absolute_turns": 2.0,
        "contraction_magnitude": 4.0,
        "mean_pair_axial_separation_change": 5.0,
        "inward_path_count": 40,
    }
    candidate = {
        "mean_absolute_turns": 1.98,
        "contraction_magnitude": 4.4,
        "mean_pair_axial_separation_change": 4.75,
        "inward_path_count": 41,
    }
    out = mod.compare_paths(candidate, control)
    assert np.isclose(out["mean_absolute_turns_relative"], -0.01)
    assert np.isclose(out["contraction_magnitude_relative"], 0.10)
    assert np.isclose(out["mean_pair_axial_separation_change_relative"], -0.05)
    assert out["inward_path_count_delta"] == 1


def test_clean_rule_thresholds_are_frozen_to_issue_585():
    c = mod.CRITERIA
    assert c["reference_energy_relative_error_max"] == 1e-10
    assert c["late_static_identity_abs_max"] == 1e-11
    assert c["aggregate_mean_turns_relative_min"] == -0.015
    assert c["aggregate_pair_axial_change_relative_min"] == -0.05
    assert c["tip_pair_axial_change_relative_min"] == -0.25
    assert c["late_effect_retention_min"] == 0.999
    assert c["mid_effect_retention_min"] == 0.30
    assert c["support_max_abs"] == 1e-12
    assert c["divergence_fd_max"] == 1e-5
    assert c["response_rank_min"] == 2
    assert c["response_condition_max"] == 8.0
