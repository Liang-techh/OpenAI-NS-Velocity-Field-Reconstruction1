import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("agent7_st052m_vector_potential_poloidal_screen.py")
SPEC = importlib.util.spec_from_file_location("agent7_vector_potential_poloidal_screen", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_preregistered_identity_and_truth_boundary():
    assert mod.PREREG_ISSUE == 612
    assert mod.SOURCE_TEMPORAL_PR == 587
    assert mod.SOURCE_TEMPORAL_HEAD == "0b93819095f6c8576a7571bdc2d2fbef4154944d"
    assert mod.GRID_RESOLUTION == 41
    assert mod.PERTURBATION_RMS_FRACTION == 0.005
    assert mod.R_SUPPORT == 1.55
    assert mod.Z_SUPPORT_INNER == 0.70
    assert mod.Z_SUPPORT_OUTER == 1.45
    assert mod.TRUTH["candidate_basis_changed"] is False
    assert mod.TRUTH["parameter_scan_performed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False


def test_endpoint_zero_envelope_is_frozen():
    assert mod.endpoint_zero_envelope(0.25) == 0.0
    assert mod.endpoint_zero_envelope(0.50) == 1.0
    assert mod.endpoint_zero_envelope(0.75) == 0.0


def test_vector_potential_correction_is_compact_and_center_exact():
    pts = np.array([
        [0.8, 0.0, 0.0],
        [0.8, 0.0, 0.65],
        [0.8, 0.0, 1.46],
        [1.56, 0.0, 1.0],
        [0.8, 0.0, 1.0],
    ])
    corr = mod.vector_potential_correction(pts)
    assert np.max(np.abs(corr[:4])) == 0.0
    assert np.linalg.norm(corr[4]) > 0.0


def test_vector_potential_correction_has_small_cartesian_fd_divergence():
    rng = np.random.default_rng(6127)
    pts = rng.uniform(-1.25, 1.25, size=(24, 3))
    h = 1.0e-5
    div = np.zeros(len(pts))
    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = h
        up = mod.vector_potential_correction(pts + shift)
        um = mod.vector_potential_correction(pts - shift)
        div += (up[:, axis] - um[:, axis]) / (2.0 * h)
    assert np.max(np.abs(div)) < 1.0e-6


def test_pareto_rule_requires_one_common_sign():
    good = {
        "plus": {
            "aspect_gain_increment": 1e-3,
            "tip_thinning_increment": 2e-3,
            "turns_fidelity_increment": 1e-3,
            "axial_pair_fidelity_increment": 1e-3,
        },
        "minus": {
            "aspect_gain_increment": -1e-3,
            "tip_thinning_increment": -2e-3,
            "turns_fidelity_increment": -1e-3,
            "axial_pair_fidelity_increment": -1e-3,
        },
    }
    r = mod._pareto(good, True)
    assert r["target_free_pareto_direction"] == "+epsilon"

    mixed = {
        "plus": dict(good["plus"], turns_fidelity_increment=-1e-3),
        "minus": dict(good["minus"], turns_fidelity_increment=1e-3),
    }
    r = mod._pareto(mixed, True)
    assert r["target_free_pareto_direction"] is None
