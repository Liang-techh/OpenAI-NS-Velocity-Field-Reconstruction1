from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_partial_composite_batch_differentials import (
    LEADING_FD6_SPATIAL_STEP,
    _curl_from_jacobian,
    _evaluate_with_interfaces,
    _fd6_jacobian,
    public_contract,
    semantic_payload,
)


class _AffineLeading:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack(
            (
                2.0 * x - y + 0.5 * z + 0.1 * t,
                x + 3.0 * y - 2.0 * z - 0.2 * t,
                -x + 0.25 * y + 4.0 * z + 0.05 * t,
            ),
            axis=-1,
        )


_LEAD_JAC = np.asarray(
    [[2.0, -1.0, 0.5], [1.0, 3.0, -2.0], [-1.0, 0.25, 4.0]], dtype=float
)
_OSC_JAC = np.asarray(
    [[0.0, 0.2, 0.0], [-0.2, 0.0, 0.0], [0.0, 0.0, 0.1]], dtype=float
)


def _osc_velocity(points):
    p = np.asarray(points, dtype=float)
    return np.stack((0.2 * p[..., 1], -0.2 * p[..., 0], 0.1 * p[..., 2]), axis=-1)


def _fake_oscillatory_differentials(points, time):
    p = np.asarray(points, dtype=float)
    shape = p.shape[:-1]
    np.broadcast_to(np.asarray(time, dtype=float), shape)
    jac = np.broadcast_to(_OSC_JAC, shape + (3, 3)).copy()
    return SimpleNamespace(
        velocity=_osc_velocity(p),
        velocity_jacobian=jac,
        divergence=np.trace(jac, axis1=-2, axis2=-1),
        vorticity=_curl_from_jacobian(jac),
        interior_mask=np.ones(shape, dtype=bool),
    )


class _FakeComposite:
    semantic_sha256 = "a" * 64

    def __init__(self):
        self.leading_backend = _AffineLeading()

    def velocity(self, x, y, z, t):
        points = np.stack(np.broadcast_arrays(x, y, z), axis=-1)
        return self.leading_backend.velocity(x, y, z, t) + _osc_velocity(points)


def test_fd6_jacobian_reconstructs_affine_leading_field():
    leading = _AffineLeading()
    x = np.asarray([0.21, -0.37, 0.52])
    y = np.asarray([-0.11, 0.26, 0.07])
    z = np.asarray([0.19, -0.08, 0.31])
    t = np.asarray([0.33, 0.47, 0.61])
    jac = _fd6_jacobian(leading.velocity, x, y, z, t)
    expected = np.broadcast_to(_LEAD_JAC, jac.shape)
    assert LEADING_FD6_SPATIAL_STEP == pytest.approx(1.0e-3)
    assert np.max(np.abs(jac - expected)) <= 2.0e-11


def test_staged_partial_composite_differentials_close_exactly():
    points = np.asarray(
        [[0.31, 0.14, -0.22], [-0.43, 0.27, 0.16], [0.62, -0.19, 0.08]],
        dtype=float,
    )
    time = np.asarray([0.35, 0.49, 0.63])
    out = _evaluate_with_interfaces(
        _FakeComposite(), _fake_oscillatory_differentials, points, time
    )
    expected_lead = np.broadcast_to(_LEAD_JAC, out.leading_velocity_jacobian.shape)
    expected_osc = np.broadcast_to(_OSC_JAC, out.oscillatory_velocity_jacobian.shape)
    expected_total = expected_lead + expected_osc

    assert np.max(np.abs(out.leading_velocity_jacobian - expected_lead)) <= 2.0e-11
    assert np.max(np.abs(out.oscillatory_velocity_jacobian - expected_osc)) == 0.0
    assert np.max(np.abs(out.velocity_jacobian - expected_total)) <= 2.0e-11
    assert np.max(
        np.abs(out.velocity - (out.leading_velocity + out.oscillatory_velocity))
    ) <= 1.0e-14
    assert np.max(
        np.abs(out.divergence - np.trace(out.velocity_jacobian, axis1=-2, axis2=-1))
    ) == 0.0
    assert np.max(np.abs(out.vorticity - _curl_from_jacobian(out.velocity_jacobian))) == 0.0
    assert np.all(out.oscillatory_interior_mask)


def test_shape_preserving_broadcast_and_empty_batches():
    points = np.asarray(
        [
            [[0.2, 0.1, 0.0], [0.3, -0.2, 0.1]],
            [[-0.4, 0.15, -0.1], [0.5, 0.25, 0.2]],
        ],
        dtype=float,
    )
    out = _evaluate_with_interfaces(
        _FakeComposite(), _fake_oscillatory_differentials, points, 0.47
    )
    assert out.velocity.shape == (2, 2, 3)
    assert out.velocity_jacobian.shape == (2, 2, 3, 3)
    assert out.divergence.shape == (2, 2)
    assert out.vorticity.shape == (2, 2, 3)
    assert out.oscillatory_interior_mask.shape == (2, 2)

    empty = np.empty((0, 3), dtype=float)
    empty_out = _evaluate_with_interfaces(
        _FakeComposite(), _fake_oscillatory_differentials, empty, np.empty((0,))
    )
    assert empty_out.velocity.shape == (0, 3)
    assert empty_out.velocity_jacobian.shape == (0, 3, 3)
    assert empty_out.divergence.shape == (0,)


def test_bad_points_time_and_backend_fail_closed():
    field = _FakeComposite()
    with pytest.raises(ValueError):
        _evaluate_with_interfaces(field, _fake_oscillatory_differentials, [1.0, 2.0], 0.5)
    with pytest.raises(ValueError):
        _evaluate_with_interfaces(
            field,
            _fake_oscillatory_differentials,
            [[np.nan, 0.0, 0.0]],
            0.5,
        )
    with pytest.raises(ValueError):
        _evaluate_with_interfaces(
            field,
            _fake_oscillatory_differentials,
            np.zeros((2, 3)),
            np.zeros((3,)),
        )
    with pytest.raises(TypeError):
        _evaluate_with_interfaces(
            object(), _fake_oscillatory_differentials, np.zeros((1, 3)), 0.5
        )


def test_composition_drift_is_rejected():
    class _BadComposite(_FakeComposite):
        def velocity(self, x, y, z, t):
            out = super().velocity(x, y, z, t)
            return out + np.asarray([1.0e-5, 0.0, 0.0])

    with pytest.raises(RuntimeError, match="leading \\+ oscillatory"):
        _evaluate_with_interfaces(
            _BadComposite(),
            _fake_oscillatory_differentials,
            np.asarray([[0.3, 0.2, 0.1]]),
            0.5,
        )


def test_public_contract_has_no_scientific_tuning_escape_hatch():
    contract = public_contract()
    assert contract["public_inputs"] == ["points", "time"]
    assert contract["forbidden_inputs_present"] == []
    assert contract["caller_tunable_derivative_step"] is False
    assert contract["leading_profile_reimplemented"] is False
    assert contract["complete_curl_reimplemented"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["radial_inverse_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["pressure_or_forcing_added"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_semantic_truth_boundary_stays_partial_and_non_pde():
    payload = semantic_payload(_FakeComposite())
    truth = payload["truth_boundary"]
    assert truth["current_partial_leading_plus_oscillatory_velocity_materialized"] is True
    assert truth["current_partial_composite_spatial_differentials_materialized"] is True
    for key in (
        "velocity_beyond_Xh_materialized",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_composite_save_load_available",
        "mean_projection_performed",
        "radial_inverse_performed",
        "correction_velocity_constructed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_velocity_pressure_forcing_api",
        "complete_ns_residual",
        "same_protocol_comparable_to_st006",
        "residual_reduction_claimed",
        "velocity_export_ready",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    ):
        assert truth[key] is False
