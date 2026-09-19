import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("agent7_st052m_threshold_free_morphology.py")
SPEC = importlib.util.spec_from_file_location("agent7_threshold_free", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_frozen_contract_and_truth_boundary():
    assert mod.PREREG_ISSUE == 575
    assert mod.SOURCE_PRIOR_HEAD == "85f2f0940d80757313651d2ab50c6ef89dc22914"
    assert mod.SOURCE_COMPENSATED_HEAD == "39b106ad8cb8df2064cabead3a12682089575e74"
    assert mod.SOURCE_PATH_HEAD == "196bcc5445039a2b3533128daa5c20cbb151cd58"
    assert mod.EXPECTED_BETA == 0.08837490297155456
    assert mod.GRID_LEVELS == (25, 33, 41)
    assert mod.REFERENCE_TIMES == (0.25, 0.50, 0.75)
    assert mod.TIP_WEIGHT_WINDOW == (0.50, 0.80)
    assert mod.TRUTH["basis_changed"] is False
    assert mod.TRUTH["velocity_coefficient_changed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False


def test_tensor_trapezoid_weights_have_expected_mass():
    for n in (5, 9):
        w = mod.tensor_trapezoid_weights(n)
        assert w.shape == (n, n, n)
        assert np.isclose(np.sum(w), float((n - 1) ** 3))
        assert w[0, 0, 0] == 0.125
        assert w[1, 1, 1] == 1.0


def test_smooth_tip_weight_is_fixed_compact_bump():
    s = np.array([0.0, 0.50, 0.575, 0.65, 0.725, 0.80, 1.0])
    w = mod.smooth_tip_weight(s)
    assert w[0] == 0.0
    assert w[1] == 0.0
    assert w[5] == 0.0
    assert w[6] == 0.0
    assert np.isclose(w[3], 1.0)
    assert 0.0 < w[2] < 1.0
    assert np.isclose(w[2], w[4])


def test_enstrophy_moments_are_amplitude_scale_invariant():
    axis = np.linspace(-2.0, 2.0, 25)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    magnitude = np.exp(-0.7 * (xx**2 + yy**2) - 0.25 * zz**2) + 0.05
    m1 = mod.enstrophy_moment_metrics(magnitude, axis)
    m2 = mod.enstrophy_moment_metrics(3.25 * magnitude, axis)
    for key in (
        "full_axial_rms",
        "full_radial_rms",
        "full_aspect_ratio",
        "smooth_tip_radial_rms",
        "smooth_tip_axial_rms",
        "smooth_tip_enstrophy_fraction",
        "central_radial_rms",
    ):
        assert np.isclose(m1[key], m2[key], rtol=2e-14, atol=2e-14)
    assert np.isclose(m2["full_enstrophy_trapezoid"] / m1["full_enstrophy_trapezoid"], 3.25**2)


def test_compare_metrics_and_point_rule_use_preregistered_directions():
    control = {
        "full_axial_rms": 1.0,
        "full_radial_rms": 1.0,
        "full_aspect_ratio": 1.0,
        "smooth_tip_radial_rms": 1.0,
        "smooth_tip_axial_rms": 1.0,
        "smooth_tip_enstrophy_fraction": 0.20,
        "central_radial_rms": 1.0,
        "full_enstrophy_trapezoid": 2.0,
    }
    child = dict(control)
    child.update(
        full_axial_rms=1.01,
        full_radial_rms=0.995,
        full_aspect_ratio=1.015,
        smooth_tip_radial_rms=0.99,
        smooth_tip_axial_rms=1.002,
        smooth_tip_enstrophy_fraction=0.198,
        central_radial_rms=1.001,
        full_enstrophy_trapezoid=2.04,
    )
    rel = mod.compare_metrics(control, child)
    assert rel["smooth_tip_radial_rms"] < 0.0
    assert rel["full_aspect_ratio"] > 0.0
    assert rel["tip_enstrophy_fraction_retention"] == 0.99
    assert mod._point_rule(rel)
    assert not mod._point_rule(dict(rel, full_aspect_ratio=0.0))
