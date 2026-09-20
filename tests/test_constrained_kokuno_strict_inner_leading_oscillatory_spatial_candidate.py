from __future__ import annotations

import copy
import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_diagnostic import (
    evaluate_vorticity_osc_fd6,
)
from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_spatial_candidate import (
    OSCILLATORY_SPATIAL_STEP,
    KokunoStrictInnerLeadingOscillatorySpatialCandidate,
    public_contract,
)


class _AffineBackend:
    field_sha256 = "1" * 64
    temporal_derivative_sha256 = "2" * 64
    spatial_derivative_sha256 = "3" * 64

    _J = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [-1.0, 1.0, -0.25],
            [0.7, -0.4, 2.0],
        ],
        dtype=float,
    )
    _DT = np.asarray([0.5, -0.3, 0.1], dtype=float)

    @staticmethod
    def _broadcast(x, y, z, t):
        return np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )

    def velocity(self, x, y, z, t):
        x, y, z, t = self._broadcast(x, y, z, t)
        return np.stack(
            (
                x + 2.0 * y + 3.0 * z + 0.5 * t,
                -x + y - 0.25 * z - 0.3 * t,
                0.7 * x - 0.4 * y + 2.0 * z + 0.1 * t,
            ),
            axis=-1,
        )

    def velocity_dt(self, x, y, z, t):
        x, _, _, _ = self._broadcast(x, y, z, t)
        return np.broadcast_to(self._DT, x.shape + (3,)).copy()

    def velocity_jacobian(self, x, y, z, t):
        x, _, _, _ = self._broadcast(x, y, z, t)
        return np.broadcast_to(self._J, x.shape + (3, 3)).copy()


def _points():
    x = np.asarray([0.43, 0.58, 0.71, 0.84], dtype=float)
    y = np.asarray([0.17, -0.26, 0.31, -0.39], dtype=float)
    z = np.asarray([-0.46, -0.11, 0.27, 0.63], dtype=float)
    t = np.asarray([0.44, 0.48, 0.52, 0.56], dtype=float)
    return x, y, z, t


def test_spatial_candidate_adds_analytic_inner_and_fixed_fd6_oscillatory_jacobians():
    backend = _AffineBackend()
    candidate = KokunoStrictInnerLeadingOscillatorySpatialCandidate(backend)
    x, y, z, t = _points()

    evaluated = candidate.evaluate(x, y, z, t)
    oscillatory = evaluate_vorticity_osc_fd6(
        x, y, z, t, spatial_step=OSCILLATORY_SPATIAL_STEP
    )
    expected_osc = np.asarray(oscillatory["velocity_gradient_fd6"], dtype=float)
    expected_inner = backend.velocity_jacobian(x, y, z, t)

    np.testing.assert_allclose(evaluated.inner_leading_jacobian, expected_inner, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(evaluated.oscillatory_jacobian, expected_osc, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        evaluated.velocity_jacobian,
        expected_inner + expected_osc,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        candidate.velocity_jacobian(x, y, z, t),
        evaluated.velocity_jacobian,
        rtol=0.0,
        atol=0.0,
    )
    assert evaluated.jacobian_additive_closure_abs_max <= 1.0e-13
    np.testing.assert_allclose(
        evaluated.divergence,
        np.trace(evaluated.velocity_jacobian, axis1=-2, axis2=-1),
        rtol=0.0,
        atol=0.0,
    )
    expected_curl = np.stack(
        (
            evaluated.velocity_jacobian[..., 2, 1] - evaluated.velocity_jacobian[..., 1, 2],
            evaluated.velocity_jacobian[..., 0, 2] - evaluated.velocity_jacobian[..., 2, 0],
            evaluated.velocity_jacobian[..., 1, 0] - evaluated.velocity_jacobian[..., 0, 1],
        ),
        axis=-1,
    )
    np.testing.assert_allclose(evaluated.vorticity, expected_curl, rtol=0.0, atol=0.0)


def test_velocity_and_time_derivative_semantic_identity_are_preserved():
    backend = _AffineBackend()
    candidate = KokunoStrictInnerLeadingOscillatorySpatialCandidate(backend)
    base = candidate.differentiable_candidate
    x, y, z, t = _points()

    np.testing.assert_allclose(candidate.velocity(x, y, z, t), base.velocity(x, y, z, t), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(candidate.velocity_dt(x, y, z, t), base.velocity_dt(x, y, z, t), rtol=0.0, atol=0.0)
    assert candidate.velocity_candidate_sha256 == base.velocity_candidate_sha256
    assert candidate.differentiable_sha256 == base.differentiable_sha256


def test_manifest_roundtrip_and_identity_firewall(tmp_path):
    backend = _AffineBackend()
    candidate = KokunoStrictInnerLeadingOscillatorySpatialCandidate(backend)
    path = tmp_path / "spatial_candidate.json"
    manifest = candidate.save_manifest(path)
    rebound = KokunoStrictInnerLeadingOscillatorySpatialCandidate.load_manifest(backend, path)
    x, y, z, t = _points()

    assert rebound.spatial_sha256 == candidate.spatial_sha256
    np.testing.assert_allclose(
        rebound.velocity_jacobian(x, y, z, t),
        candidate.velocity_jacobian(x, y, z, t),
        rtol=0.0,
        atol=0.0,
    )

    mutated = copy.deepcopy(manifest)
    mutated["spatial_configuration"]["oscillatory_spatial_derivative"]["fixed_spatial_step"] = 2.0e-3
    stripped = dict(mutated)
    stripped.pop("manifest_sha256")
    import hashlib

    raw = json.dumps(stripped, sort_keys=True, separators=(",", ":"), allow_nan=False)
    mutated["manifest_sha256"] = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    with pytest.raises(ValueError, match="configuration/provenance changed"):
        KokunoStrictInnerLeadingOscillatorySpatialCandidate.from_manifest(backend, mutated)


def test_public_jacobian_surface_has_no_scientific_tuning_knobs():
    parameters = set(
        inspect.signature(
            KokunoStrictInnerLeadingOscillatorySpatialCandidate.velocity_jacobian
        ).parameters
    )
    assert parameters == {"self", "x", "y", "z", "t"}
    contract = public_contract()
    assert contract["forbidden_velocity_jacobian_inputs_present"] is False
    assert contract["oscillatory_spatial_step"] == 1.0e-3
    truth = contract["truth_boundary"]
    assert truth["strict_inner_composite_velocity_jacobian_executable"] is True
    assert truth["inner_leading_velocity_jacobian_analytic"] is True
    assert truth["oscillatory_velocity_jacobian_uses_fixed_fd6"] is True
    assert truth["global_leading_velocity_materialized"] is False
    assert truth["agent3_correction_velocity_included"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False


def test_backend_contract_fails_closed_when_spatial_surface_is_missing():
    class MissingSpatial:
        field_sha256 = "1" * 64
        temporal_derivative_sha256 = "2" * 64
        spatial_derivative_sha256 = "3" * 64

        def velocity(self, x, y, z, t):
            x = np.asarray(x, dtype=float)
            return np.zeros(x.shape + (3,), dtype=float)

        def velocity_dt(self, x, y, z, t):
            x = np.asarray(x, dtype=float)
            return np.zeros(x.shape + (3,), dtype=float)

    with pytest.raises(TypeError, match="velocity_jacobian"):
        KokunoStrictInnerLeadingOscillatorySpatialCandidate(MissingSpatial())


def test_backend_jacobian_shape_is_checked():
    class BadShape(_AffineBackend):
        def velocity_jacobian(self, x, y, z, t):
            x, _, _, _ = self._broadcast(x, y, z, t)
            return np.zeros(x.shape + (3,), dtype=float)

    x, y, z, t = _points()
    candidate = KokunoStrictInnerLeadingOscillatorySpatialCandidate(BadShape())
    with pytest.raises(ValueError, match="velocity_jacobian must have shape"):
        candidate.velocity_jacobian(x, y, z, t)
