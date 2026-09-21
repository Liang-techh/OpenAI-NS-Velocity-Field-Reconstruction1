from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_current_partial_leading_oscillatory_velocity as mod


class _FakeLeading:
    def __init__(self, x_limit: float = 10.0) -> None:
        self.x_limit = float(x_limit)

    def velocity(self, x, y, z, t):
        xx, yy, zz, tt = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if np.any(np.abs(xx) > self.x_limit):
            raise ValueError("outside fake leading domain")
        return np.stack((xx + 2.0 * yy, zz - tt, xx - yy + 0.5 * zz), axis=-1)


def test_private_composition_replays_exact_addition_and_shape() -> None:
    points = np.asarray(
        [
            [[0.32, 0.11, -0.20], [0.41, -0.17, -0.08]],
            [[-0.36, 0.24, 0.05], [-0.52, -0.16, 0.16]],
        ],
        dtype=float,
    )
    times = np.asarray([[0.31, 0.39], [0.47, 0.55]])
    out = mod._evaluate_with_backend(
        _FakeLeading(),
        points[..., 0],
        points[..., 1],
        points[..., 2],
        times,
    )
    assert out.leading.shape == (2, 2, 3)
    assert out.oscillatory.shape == (2, 2, 3)
    assert out.velocity.shape == (2, 2, 3)
    np.testing.assert_allclose(out.velocity, out.leading + out.oscillatory, rtol=0.0, atol=0.0)
    assert np.sqrt(np.mean(np.sum(out.oscillatory**2, axis=-1))) > mod.OSCILLATORY_SIGNAL_FLOOR


def test_axis_mask_reduces_composition_to_leading() -> None:
    x = np.zeros(3)
    y = np.zeros(3)
    z = np.asarray((-0.1, 0.0, 0.12))
    t = np.asarray((0.37, 0.51, 0.67))
    out = mod._evaluate_with_backend(_FakeLeading(), x, y, z, t)
    np.testing.assert_array_equal(out.oscillatory, np.zeros((3, 3)))
    np.testing.assert_array_equal(out.velocity, out.leading)


def test_leading_domain_failure_is_not_masked_by_oscillatory_support() -> None:
    with pytest.raises(ValueError, match="outside fake leading domain"):
        mod._evaluate_with_backend(_FakeLeading(x_limit=1.0), 2.0, 0.0, 0.0, 0.5)


def test_parent_time_contract_is_preserved_even_on_axis() -> None:
    with pytest.raises(ValueError, match="time is outside"):
        mod._evaluate_with_backend(_FakeLeading(), 0.0, 0.0, 0.0, 0.20)


@pytest.mark.parametrize(
    "args",
    [
        (np.nan, 0.0, 0.0, 0.5),
        (0.0, np.inf, 0.0, 0.5),
        (np.ones(2), np.ones(3), 0.0, 0.5),
    ],
)
def test_bad_coordinates_fail_closed(args) -> None:
    with pytest.raises(ValueError):
        mod._evaluate_with_backend(_FakeLeading(), *args)


def test_public_velocity_has_no_scientific_tuning_inputs() -> None:
    assert list(inspect.signature(mod.velocity).parameters) == ["x", "y", "z", "t"]
    contract = mod.public_contract()
    assert contract["forbidden_inputs_present"] == []
    assert contract["new_oscillatory_parameters"] is False
    assert contract["leading_profile_reimplemented"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["pressure_or_forcing_added"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_exact_agent1_runtime_is_lazy_on_a2_parent_lineage() -> None:
    assert mod.AGENT1_PR == 965
    assert mod.AGENT1_HEAD == "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
    assert mod.AGENT1_SOURCE_BLOB_SHA1 == "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"


def test_git_blob_hash_helper_uses_git_blob_object_convention(tmp_path: Path) -> None:
    path = tmp_path / "sample.py"
    path.write_bytes(b"abc\n")
    import hashlib

    expected = hashlib.sha1(b"blob 4\0abc\n").hexdigest()
    assert mod._git_blob_sha1(path) == expected


def test_required_truth_boundary_keeps_partial_scope_fail_closed() -> None:
    required = mod._required_agent1_truth()
    assert required["current_joined_cartesian_spacetime_leading_velocity_through_Xh_materialized"] is True
    for key in (
        "outer_global_leading_velocity_materialized",
        "velocity_beyond_Xh_materialized",
        "unified_global_cartesian_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
    ):
        assert required[key] is False
