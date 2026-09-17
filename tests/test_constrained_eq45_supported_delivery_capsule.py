from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_delivery_capsule import (
    EXPECTED_CHILD_SHA256,
    EXPECTED_PARENT_SHA256,
    build_eq45_supported_delivery_capsule,
    load_checked_eq45_supported_delivery_capsule,
    public_velocity_probe,
    validate_eq45_supported_delivery_capsule,
    write_eq45_supported_delivery_capsule,
)
from openai_ns_reconstruction.eq45_supported_delivery import (
    Eq45SupportedDeliveryField,
    default_field,
)


ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / "artifacts" / "constrained" / "eq45_supported_delivery_capsule.json"


def test_checked_capsule_matches_current_supported_delivery(tmp_path: Path) -> None:
    checked = load_checked_eq45_supported_delivery_capsule(CAPSULE, ROOT)
    rebuilt = build_eq45_supported_delivery_capsule(ROOT)
    assert checked == rebuilt
    assert checked["candidate"]["child_sha256"] == EXPECTED_CHILD_SHA256
    assert checked["candidate"]["parent_sha256"] == EXPECTED_PARENT_SHA256
    assert checked["states"]["velocity_export_ready"] is True
    assert checked["states"]["visualization_ready"] is False
    assert checked["states"]["pde_validated"] is False

    regenerated = tmp_path / "capsule.json"
    write_eq45_supported_delivery_capsule(regenerated, ROOT)
    assert json.loads(regenerated.read_text(encoding="utf-8")) == checked


def test_capsule_public_field_replays_velocity_grid_and_save_load(tmp_path: Path) -> None:
    field = default_field()
    probe = public_velocity_probe()
    assert probe.shape == (3,)
    assert np.all(np.isfinite(probe))

    support_face = field.velocity(2.0, 0.0, 0.0, 0.5)
    np.testing.assert_array_equal(support_face, np.zeros(3))

    axes = np.array([-1.0, 0.0, 1.0])
    times = np.array([0.25, 0.5, 0.75])
    grid = field.grid(axes, axes, axes, times)
    assert grid.shape == (3, 3, 3, 3, 3)
    assert np.all(np.isfinite(grid))

    path = tmp_path / "supported_candidate.json"
    field.save_candidate(path)
    reloaded = Eq45SupportedDeliveryField.load_candidate(path)
    assert reloaded.sha256 == field.sha256 == EXPECTED_CHILD_SHA256

    points = np.array([[0.37, -0.29, 0.41], [1.72, 0.11, -1.66]])
    np.testing.assert_allclose(
        reloaded.at_points(points, 0.5), field.at_points(points, 0.5), rtol=0.0, atol=0.0
    )


def test_capsule_fails_closed_on_identity_claim_or_asset_promotion() -> None:
    capsule = build_eq45_supported_delivery_capsule(ROOT)

    bad_identity = deepcopy(capsule)
    bad_identity["candidate"]["child_sha256"] = EXPECTED_PARENT_SHA256
    with pytest.raises(ValueError, match="does not match governed state"):
        validate_eq45_supported_delivery_capsule(bad_identity, ROOT)

    bad_claim = deepcopy(capsule)
    bad_claim["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="does not match governed state"):
        validate_eq45_supported_delivery_capsule(bad_claim, ROOT)

    bad_asset = deepcopy(capsule)
    bad_asset["integrated_assets"]["matlab_mat_exporter"] = True
    with pytest.raises(ValueError, match="does not match governed state"):
        validate_eq45_supported_delivery_capsule(bad_asset, ROOT)
