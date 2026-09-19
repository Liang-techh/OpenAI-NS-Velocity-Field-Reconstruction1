import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("agent7_st052m_linear_ramp_convergence.py")
SPEC = importlib.util.spec_from_file_location("agent7_linear_ramp_convergence", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_preregistered_identity_and_truth_boundary():
    assert mod.PREREG_ISSUE == 593
    assert mod.SOURCE_TEMPORAL_PR == 587
    assert mod.SOURCE_TEMPORAL_HEAD == "0b93819095f6c8576a7571bdc2d2fbef4154944d"
    assert mod.SOURCE_STATIC_PR == 559
    assert mod.SOURCE_STATIC_HEAD == "39b106ad8cb8df2064cabead3a12682089575e74"
    assert mod.EXPECTED_BETA == 0.08837490297155456
    assert mod.GRID_LEVELS == (25, 33, 41)
    assert mod.REFERENCE_TIMES == (0.25, 0.375, 0.50, 0.625, 0.75)
    assert mod.TRUTH["basis_changed"] is False
    assert mod.TRUTH["canonical_velocity_changed"] is False
    assert mod.TRUTH["saved_velocity_changed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False


def test_clean_rule_thresholds_are_frozen_to_issue_593():
    c = mod.CRITERIA
    assert c["identity_relative_abs_max"] == 1e-12
    assert c["active_aspect_change_min"] == 0.0
    assert c["active_tip_radial_change_max"] == 0.0
    assert c["central_radial_abs_change_max"] == 0.005
    assert c["tip_enstrophy_fraction_retention_min"] == 0.95
    assert c["fine_effect_drift_abs_max"] == 0.0075
    assert c["finest_adjacent_effect_backslide_max"] == 0.001
    assert c["late_static_identity_abs_max"] == 1e-11


def test_active_point_rule_accepts_only_desired_direction_and_locality():
    good = {
        "full_aspect_ratio": 0.01,
        "smooth_tip_radial_rms": -0.01,
        "central_radial_rms": 0.001,
        "tip_enstrophy_fraction_retention": 1.05,
    }
    assert mod._active_point_pass(good)
    bad_aspect = dict(good, full_aspect_ratio=-1e-4)
    bad_tip = dict(good, smooth_tip_radial_rms=1e-4)
    bad_center = dict(good, central_radial_rms=0.006)
    bad_retention = dict(good, tip_enstrophy_fraction_retention=0.94)
    assert not mod._active_point_pass(bad_aspect)
    assert not mod._active_point_pass(bad_tip)
    assert not mod._active_point_pass(bad_center)
    assert not mod._active_point_pass(bad_retention)
