import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_axial_taper_delivery_capsule import (
    _validate_recipe,
    build_candidate,
    delivery_capsule,
    load_recipe,
    write_delivery_bundle,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_axial_taper_mode import (
    BASE_AXIAL_PLATEAU_Q,
    Eq45SupportedPhi10BlendAxialTaperCandidate,
)


def test_recipe_keeps_selection_and_scientific_states_separate():
    recipe = load_recipe()
    assert recipe["status"]["velocity_export_ready"] is True
    assert recipe["status"]["caller_declared_blend_weight"] is True
    assert recipe["status"]["caller_declared_axial_taper"] is True
    assert recipe["status"]["blend_weight_selected"] is False
    assert recipe["status"]["axial_taper_value_selected"] is False
    assert recipe["status"]["visualization_ready"] is False
    assert recipe["status"]["pde_validated"] is False
    assert recipe["sibling_evidence"]["restricted_force_status"] == "open_unconsumed"
    assert recipe["sibling_evidence"]["divergence_status"] == "open_unconsumed"
    assert recipe["sibling_evidence"]["morphology_status"] == (
        "open_unconsumed_calibration_pending"
    )


def test_base_taper_endpoint_replays_declared_blend_exactly():
    candidate = build_candidate(0.5, BASE_AXIAL_PLATEAU_Q)
    points = np.asarray(
        [
            [0.31, 0.17, 0.22],
            [0.55, -0.21, 1.67],
            [1.71, 0.12, 0.0],
        ],
        dtype=float,
    )
    for time in (0.3125, 0.6875):
        actual = candidate.at_points(points, time)
        expected = candidate.base.at_points(points, time)
        np.testing.assert_array_equal(actual, expected)


def test_delivery_bundle_roundtrips_exact_candidate_and_public_velocity(tmp_path):
    paths = write_delivery_bundle(tmp_path, blend_weight=0.5, axial_plateau_q=0.7225)
    candidate = Eq45SupportedPhi10BlendAxialTaperCandidate.load_json(paths["candidate"])
    capsule = json.loads(paths["capsule"].read_text(encoding="utf-8"))

    assert capsule["candidate_sha256"] == candidate.sha256
    assert capsule["base_blend_sha256"] == candidate.base_sha256
    assert capsule["base_supported_sha256"] == candidate.supported_base_sha256
    assert capsule["declared_blend_weight"] == 0.5
    assert capsule["declared_axial_plateau_q"] == 0.7225
    assert capsule["blend_weight_selected_by_capsule"] is False
    assert capsule["axial_taper_selected_by_capsule"] is False
    assert capsule["status"]["visualization_ready"] is False
    assert capsule["status"]["pde_validated"] is False

    rebuilt = build_candidate(0.5, 0.7225)
    points = np.asarray(
        [[0.29, -0.13, 1.74], [0.83, 0.26, 0.41], [1.74, 0.0, 1.73]],
        dtype=float,
    )
    for time in (0.3125, 0.5, 0.6875):
        np.testing.assert_array_equal(
            candidate.at_points(points, time), rebuilt.at_points(points, time)
        )


def test_capsule_public_grid_layout_and_outer_support_face():
    candidate = build_candidate(0.25, 0.81)
    grid = candidate.grid(
        np.asarray([-0.4, 0.4]),
        np.asarray([-0.2, 0.2]),
        np.asarray([-2.0, 0.0, 2.0]),
        np.asarray([0.3125, 0.6875]),
    )
    assert grid.shape == (2, 2, 2, 3, 3)
    assert np.all(np.isfinite(grid))
    np.testing.assert_array_equal(grid[:, :, :, 0, :], 0.0)
    np.testing.assert_array_equal(grid[:, :, :, -1, :], 0.0)
    assert np.any(np.abs(grid[:, :, :, 1, :]) > 0.0)

    capsule = delivery_capsule(0.25, 0.81)
    assert capsule["delivery"]["grid_layout"] == "time,x,y,z,component"
    assert capsule["delivery"]["component_order"] == ["u", "v", "w"]
    assert capsule["axial_identity_half_height"] == pytest.approx(1.8)


@pytest.mark.parametrize(
    ("blend_weight", "axial_plateau_q"),
    [(-0.01, 0.64), (1.01, 0.64), (0.5, 0.63), (0.5, 0.82), (np.nan, 0.64)],
)
def test_caller_controls_remain_bounded(blend_weight, axial_plateau_q):
    with pytest.raises(ValueError):
        build_candidate(blend_weight, axial_plateau_q)


def test_recipe_mutations_fail_closed():
    recipe = load_recipe()

    promoted = copy.deepcopy(recipe)
    promoted["status"]["visualization_ready"] = True
    with pytest.raises(ValueError):
        _validate_recipe(promoted)

    promoted = copy.deepcopy(recipe)
    promoted["status"]["pde_validated"] = True
    with pytest.raises(ValueError):
        _validate_recipe(promoted)

    selected = copy.deepcopy(recipe)
    selected["representation"]["axial_plateau_q_selected"] = True
    with pytest.raises(ValueError):
        _validate_recipe(selected)

    consumed = copy.deepcopy(recipe)
    consumed["sibling_evidence"]["divergence_status"] = "consumed"
    with pytest.raises(ValueError):
        _validate_recipe(consumed)
