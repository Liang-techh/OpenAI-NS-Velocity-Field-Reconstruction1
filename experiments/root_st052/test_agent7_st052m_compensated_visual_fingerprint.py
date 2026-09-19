import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("agent7_st052m_compensated_visual_fingerprint.py")
SPEC = importlib.util.spec_from_file_location("agent7_visual_fp", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_frozen_contract_and_truth_boundary():
    assert mod.PREREG_ISSUE == 567
    assert mod.SOURCE_COMPENSATED_HEAD == "39b106ad8cb8df2064cabead3a12682089575e74"
    assert mod.SOURCE_RENDER_HEAD == "619e3d2ff8c0231c9aced4ce20216c03d78cdef8"
    assert mod.SOURCE_PATH_HEAD == "196bcc5445039a2b3533128daa5c20cbb151cd58"
    assert mod.EXPECTED_BETA == 0.08837490297155456
    assert mod.GRID_LEVELS == (25, 33)
    assert mod.REFERENCE_TIMES == (0.25, 0.50, 0.75)
    assert mod.VORTICITY_QUANTILE == 0.985
    assert mod.TRUTH["canonical_velocity_changed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False


def test_vorticity_recovers_solid_rotation():
    axis = np.linspace(-1.0, 1.0, 17)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    field = np.stack((-yy, xx, np.zeros_like(zz)), axis=-1)
    omega, magnitude = mod.vorticity(field, axis[1] - axis[0])
    assert np.max(np.abs(omega[..., 0])) < 1e-12
    assert np.max(np.abs(omega[..., 1])) < 1e-12
    assert np.max(np.abs(omega[..., 2] - 2.0)) < 1e-12
    assert np.max(np.abs(magnitude - 2.0)) < 1e-12


def test_compare_metrics_reports_expected_signs():
    control = {
        "selected_axial_rms": 1.0,
        "selected_radial_rms": 1.0,
        "tip_enstrophy_weighted_radial_rms": 1.0,
        "central_enstrophy_weighted_radial_rms": 1.0,
        "full_vorticity_rms": 1.0,
        "selected_max_abs_z": 1.0,
    }
    child = dict(control)
    child.update(
        selected_axial_rms=1.08,
        selected_radial_rms=0.96,
        tip_enstrophy_weighted_radial_rms=0.98,
        central_enstrophy_weighted_radial_rms=1.001,
        full_vorticity_rms=1.03,
        selected_max_abs_z=1.2,
    )
    rel = mod.compare_metrics(control, child)
    assert rel["top_axial_rms"] > 0
    assert rel["top_radial_rms"] < 0
    assert rel["tip_radial_rms"] < 0
    assert abs(rel["central_radial_rms"]) < 0.005
    assert rel["top_max_abs_z_delta"] > 0


def test_sign_agreement_requires_desired_direction_on_both_grids():
    good = {"top_axial_rms": 0.05, "top_radial_rms": -0.02, "tip_radial_rms": -0.01}
    bad = dict(good, top_radial_rms=0.01)
    assert mod._same_required_sign(good, good)
    assert not mod._same_required_sign(good, bad)
