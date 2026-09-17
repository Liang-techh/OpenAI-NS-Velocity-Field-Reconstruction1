import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_c2_compact_temporal_capacity import (
    compact_phi10_snapshot,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


PROBES = np.array(
    [
        [0.22, 0.11, -0.37],
        [0.71, -0.24, 0.18],
        [1.23, 0.17, -0.53],
        [1.72, -0.08, 0.41],
    ],
    dtype=float,
)


def test_blend_candidate_materializes_screened_endpoint_families_and_midpoint():
    base = governed_supported_seed()
    quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)
    blend0 = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base, blend_weight=0.0
    )
    blend1 = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base, blend_weight=1.0
    )
    blend50 = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base, blend_weight=0.5
    )

    for time in (0.25, 0.3125, 0.5, 0.625, 0.75):
        np.testing.assert_allclose(
            blend0.at_points(PROBES, time), quartic.at_points(PROBES, time), rtol=0.0, atol=0.0
        )
        compact = compact_phi10_snapshot(base, time)
        np.testing.assert_allclose(
            blend1.at_points(PROBES, time), compact.at_points(PROBES, time), rtol=0.0, atol=0.0
        )

    # Every endpoint family shares the selected early snapshot and static anchors.
    for candidate in (blend0, blend50, blend1):
        np.testing.assert_allclose(
            candidate.at_points(PROBES, 0.25), quartic.at_points(PROBES, 0.25), rtol=0.0, atol=0.0
        )
        for time in (0.5, 0.625, 0.75):
            np.testing.assert_allclose(
                candidate.at_points(PROBES, time), base.at_points(PROBES, time), rtol=0.0, atol=0.0
            )

    # At the screened off-keyframe the one-mode public velocity is affine in lambda.
    time = 0.3125
    expected = 0.5 * (
        blend0.at_points(PROBES, time) + blend1.at_points(PROBES, time)
    )
    np.testing.assert_allclose(
        blend50.at_points(PROBES, time), expected, rtol=0.0, atol=2.0e-15
    )
    assert blend50.coefficient_at(time) == pytest.approx(-1.1380789910456575)


def test_blend_candidate_exposes_mixed_time_and_regular_grid_velocity():
    candidate = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=governed_supported_seed(), blend_weight=0.35
    )
    times = np.array([0.25, 0.3125, 0.4375, 0.6875])
    velocity = candidate.at_points(PROBES, times)
    assert velocity.shape == PROBES.shape
    assert np.all(np.isfinite(velocity))

    grid = candidate.grid(
        np.array([-0.4, 0.2]),
        np.array([0.0]),
        np.array([-0.3, 0.5]),
        np.array([0.25, 0.5, 0.75]),
    )
    assert grid.shape == (3, 2, 1, 2, 3)
    assert np.all(np.isfinite(grid))


def test_blend_candidate_json_sha_roundtrip_and_fail_closed_metadata(tmp_path):
    candidate = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=governed_supported_seed(), blend_weight=0.4
    )
    path = tmp_path / "blend.json"
    candidate.save_json(path)
    restored = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.load_json(path)
    assert restored.sha256 == candidate.sha256
    assert restored.to_dict() == candidate.to_dict()
    np.testing.assert_allclose(
        restored.at_points(PROBES, 0.333), candidate.at_points(PROBES, 0.333), rtol=0.0, atol=0.0
    )

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["temporal_mode"]["blend_weight_selected"] = True
    with pytest.raises(ValueError, match="selection provenance"):
        Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.from_dict(payload)

    payload = candidate.to_dict()
    payload["truth_boundary"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.from_dict(payload)

    payload = candidate.to_dict()
    payload["temporal_mode"]["compact_return_tau"] += 1.0e-3
    with pytest.raises(ValueError, match="compact return"):
        Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.from_dict(payload)


def test_blend_candidate_rejects_invalid_weight_time_and_degenerate_amplitude():
    base = governed_supported_seed()
    for value in (-0.01, 1.01, np.nan, np.inf):
        with pytest.raises(ValueError, match=r"\[0,1\]"):
            Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
                base=base, blend_weight=value
            )
    with pytest.raises(ValueError, match="nonzero"):
        Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
            base=base, blend_weight=0.5, early_delta=0.0
        )

    candidate = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base, blend_weight=0.5
    )
    with pytest.raises(ValueError, match="declared interval"):
        candidate.at_points(PROBES, 0.8)
    with pytest.raises(ValueError, match="scalar"):
        candidate.snapshot(np.array([0.25, 0.5]))

    truth = candidate.to_dict()["truth_boundary"]
    assert truth["screened_temporal_blend_materialized"] is True
    assert truth["velocity_export_ready"] is True
    assert truth["visualization_candidate_only"] is True
    for key in (
        "canonical_velocity_changed",
        "blend_weight_selected_by_this_module",
        "physical_support_validated",
        "spatial_basis_grown",
        "force_or_pressure_fitted",
        "pde_objective_used_to_choose_blend",
        "public_image_fitted",
        "production_temporal_shape_promoted",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
