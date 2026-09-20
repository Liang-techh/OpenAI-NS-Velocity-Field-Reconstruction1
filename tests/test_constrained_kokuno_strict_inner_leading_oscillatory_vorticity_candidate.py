from __future__ import annotations

import copy
import hashlib
import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_vorticity_candidate import (
    KokunoStrictInnerLeadingOscillatoryVorticityCandidate,
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
    return (
        np.asarray([0.43, 0.58, 0.71, 0.84], dtype=float),
        np.asarray([0.17, -0.26, 0.31, -0.39], dtype=float),
        np.asarray([-0.46, -0.11, 0.27, 0.63], dtype=float),
        np.asarray([0.44, 0.48, 0.52, 0.56], dtype=float),
    )


def _curl_from_jacobian(jacobian):
    j = np.asarray(jacobian, dtype=float)
    return np.stack(
        (
            j[..., 2, 1] - j[..., 1, 2],
            j[..., 0, 2] - j[..., 2, 0],
            j[..., 1, 0] - j[..., 0, 1],
        ),
        axis=-1,
    )


def test_vorticity_is_exact_curl_of_existing_spatial_candidate():
    candidate = KokunoStrictInnerLeadingOscillatoryVorticityCandidate(_AffineBackend())
    x, y, z, t = _points()
    evaluated = candidate.evaluate(x, y, z, t)
    expected = _curl_from_jacobian(evaluated.velocity_jacobian)

    np.testing.assert_allclose(evaluated.vorticity, expected, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(candidate.vorticity(x, y, z, t), expected, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        evaluated.divergence,
        np.trace(evaluated.velocity_jacobian, axis1=-2, axis2=-1),
        rtol=0.0,
        atol=0.0,
    )


def test_parent_velocity_time_and_spatial_identities_are_preserved():
    candidate = KokunoStrictInnerLeadingOscillatoryVorticityCandidate(_AffineBackend())
    parent = candidate.spatial_candidate
    x, y, z, t = _points()

    np.testing.assert_allclose(candidate.velocity(x, y, z, t), parent.velocity(x, y, z, t), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(candidate.velocity_dt(x, y, z, t), parent.velocity_dt(x, y, z, t), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(candidate.velocity_jacobian(x, y, z, t), parent.velocity_jacobian(x, y, z, t), rtol=0.0, atol=0.0)
    assert candidate.velocity_candidate_sha256 == parent.velocity_candidate_sha256
    assert candidate.differentiable_sha256 == parent.differentiable_sha256
    assert candidate.spatial_sha256 == parent.spatial_sha256


def test_batch_and_scalar_vorticity_match():
    candidate = KokunoStrictInnerLeadingOscillatoryVorticityCandidate(_AffineBackend())
    x, y, z, t = _points()
    batch = candidate.vorticity(x, y, z, t)
    scalar = np.stack(
        [candidate.vorticity(x[i], y[i], z[i], t[i]) for i in range(x.size)], axis=0
    )
    np.testing.assert_allclose(batch, scalar, rtol=0.0, atol=1.0e-12)


def test_manifest_roundtrip_and_vorticity_identity_firewall(tmp_path):
    backend = _AffineBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryVorticityCandidate(backend)
    path = tmp_path / "vorticity_candidate.json"
    manifest = candidate.save_manifest(path)
    rebound = KokunoStrictInnerLeadingOscillatoryVorticityCandidate.load_manifest(backend, path)
    x, y, z, t = _points()

    assert rebound.vorticity_sha256 == candidate.vorticity_sha256
    np.testing.assert_allclose(rebound.vorticity(x, y, z, t), candidate.vorticity(x, y, z, t), rtol=0.0, atol=0.0)

    mutated = copy.deepcopy(manifest)
    mutated["vorticity_configuration"]["composition"] = "mutated"
    stripped = dict(mutated)
    stripped.pop("manifest_sha256")
    raw = json.dumps(stripped, sort_keys=True, separators=(",", ":"), allow_nan=False)
    mutated["manifest_sha256"] = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    with pytest.raises(ValueError, match="vorticity_configuration changed"):
        KokunoStrictInnerLeadingOscillatoryVorticityCandidate.from_manifest(backend, mutated)


def test_public_vorticity_surface_has_no_scientific_tuning_knobs():
    parameters = set(inspect.signature(KokunoStrictInnerLeadingOscillatoryVorticityCandidate.vorticity).parameters)
    assert parameters == {"self", "x", "y", "z", "t"}
    contract = public_contract()
    assert contract["forbidden_vorticity_inputs_present"] is False
    truth = contract["truth_boundary"]
    assert truth["strict_inner_composite_vorticity_executable"] is True
    assert truth["vorticity_derived_from_existing_spatial_candidate"] is True
    assert truth["new_derivative_realization_introduced"] is False
    assert truth["three_resolution_morphology_verified"] is False
    assert truth["global_leading_velocity_materialized"] is False
    assert truth["agent3_correction_velocity_included"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False


def test_nonfinite_backend_jacobian_fails_closed():
    class NonFiniteBackend(_AffineBackend):
        def velocity_jacobian(self, x, y, z, t):
            out = super().velocity_jacobian(x, y, z, t)
            out[..., 0, 0] = np.nan
            return out

    candidate = KokunoStrictInnerLeadingOscillatoryVorticityCandidate(NonFiniteBackend())
    x, y, z, t = _points()
    with pytest.raises(ValueError, match="velocity_jacobian must be finite"):
        candidate.vorticity(x, y, z, t)
