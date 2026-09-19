import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("agent7_st052m_temporal_curvature_sensitivity.py")
SPEC = importlib.util.spec_from_file_location("agent7_temporal_curvature_sensitivity", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_preregistered_identity_and_truth_boundary():
    assert mod.PREREG_ISSUE == 603
    assert mod.SOURCE_TEMPORAL_PR == 587
    assert mod.SOURCE_TEMPORAL_HEAD == "0b93819095f6c8576a7571bdc2d2fbef4154944d"
    assert mod.SOURCE_STATIC_PR == 559
    assert mod.SOURCE_STATIC_HEAD == "39b106ad8cb8df2064cabead3a12682089575e74"
    assert mod.EPSILON == 0.02
    assert mod.ZERO_TOL == 1e-6
    assert mod.GRID_RESOLUTION == 41
    assert mod.MID_TIME == 0.50
    assert mod.TRUTH["candidate_basis_changed"] is False
    assert mod.TRUTH["parameter_scan_performed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False


def test_curvature_basis_preserves_endpoints_and_midpoint_normalization():
    assert mod.curvature_basis(0.25) == 0.0
    assert mod.curvature_basis(0.50) == 1.0
    assert mod.curvature_basis(0.75) == 0.0
    assert mod.curved_activation(0.25, mod.EPSILON) == 0.0
    assert mod.curved_activation(0.25, -mod.EPSILON) == 0.0
    assert mod.curved_activation(0.75, mod.EPSILON) == 1.0
    assert mod.curved_activation(0.75, -mod.EPSILON) == 1.0
    assert mod.curved_activation(0.50, mod.EPSILON) == 0.52
    assert mod.curved_activation(0.50, -mod.EPSILON) == 0.48


def test_local_pareto_rule_requires_one_common_improving_sign():
    positive = {
        "aspect_desirability": 1.0,
        "tip_thinning_desirability": 2.0,
        "path_turns_desirability": 0.5,
        "path_pair_axial_desirability": 0.25,
    }
    r = mod.local_pareto(positive)
    assert r["target_free_local_pareto_direction"] == "+lambda"

    negative = {k: -v for k, v in positive.items()}
    r = mod.local_pareto(negative)
    assert r["target_free_local_pareto_direction"] == "-lambda"

    mixed = dict(positive, path_turns_desirability=-0.5)
    r = mod.local_pareto(mixed)
    assert r["target_free_local_pareto_direction"] is None


def test_curved_activation_rejects_unregistered_amplitude():
    try:
        mod.curved_activation(0.5, 0.021)
    except ValueError:
        pass
    else:
        raise AssertionError("unregistered lambda must be rejected")
