import copy

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_axial_taper_mode import (
    BASE_AXIAL_PLATEAU_Q,
    MAX_AXIAL_PLATEAU_Q,
    Eq45SupportedPhi10BlendAxialTaperCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)


def _blend(weight=0.5):
    return Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=governed_supported_seed(),
        blend_weight=weight,
    )


def test_base_axial_taper_endpoint_exactly_replays_existing_blend():
    base = _blend()
    candidate = Eq45SupportedPhi10BlendAxialTaperCandidate(
        base=base,
        axial_plateau_q=BASE_AXIAL_PLATEAU_Q,
    )
    points = np.asarray(
        [
            [0.35, 0.0, 1.68],
            [0.80, 0.0, -1.82],
            [1.72, 0.0, 0.0],
            [0.45, 0.25, 0.8],
        ]
    )
    for time in (0.25, 0.3125, 0.5, 0.75):
        np.testing.assert_array_equal(
            candidate.at_points(points, time),
            base.at_points(points, time),
        )

    assert candidate.axial_identity_half_height == pytest.approx(1.6)
    assert candidate.taper.radial_plateau_q == base.base.taper.radial_plateau_q
    assert candidate.taper.radial_support == base.base.taper.radial_support
    assert candidate.taper.axial_half_height == base.base.taper.axial_half_height


def test_extended_axial_plateau_changes_only_axial_collar_and_keeps_support_faces_zero():
    base = _blend()
    baseline = Eq45SupportedPhi10BlendAxialTaperCandidate(base=base, axial_plateau_q=0.64)
    extended = Eq45SupportedPhi10BlendAxialTaperCandidate(base=base, axial_plateau_q=0.81)
    time = 0.3125

    unchanged_points = np.asarray(
        [
            [0.4, 0.2, 0.7],
            [1.55, 0.0, 0.0],
            [1.72, 0.0, 0.0],
            [1.84, 0.0, 0.0],
        ]
    )
    np.testing.assert_array_equal(
        extended.at_points(unchanged_points, time),
        baseline.at_points(unchanged_points, time),
    )

    axial_collar = np.asarray(
        [
            [0.35, 0.0, 1.68],
            [0.80, 0.0, 1.75],
            [1.25, 0.0, 1.82],
            [0.35, 0.0, -1.68],
            [0.80, 0.0, -1.75],
            [1.25, 0.0, -1.82],
        ]
    )
    delta = extended.at_points(axial_collar, time) - baseline.at_points(axial_collar, time)
    assert float(np.linalg.norm(delta)) > 1.0e-10
    assert extended.axial_identity_half_height == pytest.approx(1.8)

    support_faces = np.asarray(
        [
            [0.5, 0.0, 2.0],
            [0.5, 0.0, -2.0],
            [2.0, 0.0, 0.5],
            [-2.0, 0.0, -0.5],
            [2.1, 0.0, 0.0],
            [0.0, 0.0, 2.1],
        ]
    )
    np.testing.assert_array_equal(
        extended.at_points(support_faces, time),
        np.zeros((support_faces.shape[0], 3)),
    )


def test_axial_taper_candidate_mixed_time_grid_and_json_roundtrip(tmp_path):
    candidate = Eq45SupportedPhi10BlendAxialTaperCandidate(
        base=_blend(0.5), axial_plateau_q=0.75
    )
    points = np.asarray([[0.4, 0.0, 1.7], [0.8, 0.2, 1.72], [1.0, 0.0, 0.0]])
    times = np.asarray([0.25, 0.3125, 0.75])
    values = candidate.at_points(points, times)
    assert values.shape == (3, 3)
    assert np.all(np.isfinite(values))
    assert float(np.linalg.norm(values)) > 0.0

    grid = candidate.grid(
        np.linspace(-0.5, 0.5, 3),
        np.linspace(-0.25, 0.25, 2),
        np.linspace(-1.8, 1.8, 4),
        np.asarray([0.25, 0.5]),
    )
    assert grid.shape == (2, 3, 2, 4, 3)
    assert np.all(np.isfinite(grid))

    path = tmp_path / "axial_taper_candidate.json"
    candidate.save_json(path)
    loaded = Eq45SupportedPhi10BlendAxialTaperCandidate.load_json(path)
    assert loaded.sha256 == candidate.sha256
    assert loaded.to_dict() == candidate.to_dict()
    np.testing.assert_array_equal(loaded.at_points(points, times), values)


def test_axial_taper_bounds_and_truth_metadata_fail_closed():
    base = _blend()
    for value in (BASE_AXIAL_PLATEAU_Q - 1.0e-6, MAX_AXIAL_PLATEAU_Q + 1.0e-6, np.nan):
        with pytest.raises(ValueError, match="autonomous trial bound"):
            Eq45SupportedPhi10BlendAxialTaperCandidate(base=base, axial_plateau_q=value)

    candidate = Eq45SupportedPhi10BlendAxialTaperCandidate(base=base, axial_plateau_q=0.75)
    truth = candidate.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["visualization_candidate_only"] is True
    for key in (
        "axial_taper_value_selected_by_this_module",
        "canonical_velocity_changed",
        "spatial_profile_basis_grown",
        "temporal_family_grown",
        "radial_taper_changed",
        "physical_support_outer_faces_changed",
        "force_or_pressure_fitted",
        "pde_objective_used_to_choose_axial_taper",
        "public_image_fitted",
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False

    payload = candidate.to_dict()
    modified = copy.deepcopy(payload)
    modified["axial_taper_control"]["bounds"] = [0.5, 0.9]
    with pytest.raises(ValueError, match="bounds were modified"):
        Eq45SupportedPhi10BlendAxialTaperCandidate.from_dict(modified)

    modified = copy.deepcopy(payload)
    modified["truth_boundary"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="truth-boundary metadata"):
        Eq45SupportedPhi10BlendAxialTaperCandidate.from_dict(modified)
