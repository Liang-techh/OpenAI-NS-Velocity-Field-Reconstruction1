from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_nonlinear_handoff import (
    KokunoOscillatoryNonlinearHandoff,
    public_contract,
)


class PolynomialInnerBackend:
    field_sha256 = "1" * 64

    @staticmethod
    def _broadcast(x, y, z, t):
        return np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )

    def velocity(self, x, y, z, t):
        x, y, z, _ = self._broadcast(x, y, z, t)
        return np.stack((0.2 * x + 0.1, -0.15 * y + 0.05, 0.1 * z - 0.02), axis=-1)

    def velocity_dt(self, x, y, z, t):
        x, _, _, _ = self._broadcast(x, y, z, t)
        return np.zeros(x.shape + (3,), dtype=float)

    def self_advection(self, x, y, z, t):
        u = self.velocity(x, y, z, t)
        return np.stack((0.2 * u[..., 0], -0.15 * u[..., 1], 0.1 * u[..., 2]), axis=-1)

    def velocity_laplacian(self, x, y, z, t):
        return self.velocity_dt(x, y, z, t)


def _probes():
    return (
        np.asarray([0.42, 0.58, 0.73, 0.91]),
        np.asarray([0.17, -0.21, 0.26, -0.18]),
        np.asarray([-0.45, -0.12, 0.23, 0.51]),
        np.asarray([0.47, 0.49, 0.51, 0.53]),
    )


def test_nonlinear_handoff_decomposes_parent_aggregate():
    handoff = KokunoOscillatoryNonlinearHandoff(PolynomialInnerBackend())
    result = handoff.evaluate(*_probes())
    reconstructed = (
        result.inner_advects_oscillation
        + result.oscillation_advects_inner
        + result.oscillatory_self_advection
    )
    assert result.oscillatory_nonlinear_increment.shape == (4, 3)
    assert np.all(np.isfinite(reconstructed))
    assert np.max(np.abs(reconstructed - result.oscillatory_nonlinear_increment)) <= 5.0e-13
    assert np.max(
        np.abs(
            result.mixed_cross_advection
            - result.inner_advects_oscillation
            - result.oscillation_advects_inner
        )
    ) <= 5.0e-13
    assert result.decomposition_closure_abs_max <= 5.0e-13
    assert np.sqrt(np.mean(np.sum(result.oscillatory_self_advection**2, axis=-1))) > 1.0e-12


def test_contract_stays_raw_cartesian_and_parameter_closed():
    contract = public_contract()
    assert contract["forbidden_evaluate_inputs_present"] == []
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False
    assert contract["decomposition"] == [
        "inner_advects_oscillation",
        "oscillation_advects_inner",
        "oscillatory_self_advection",
    ]


def test_manifest_rebind_and_semantic_tamper_fail_closed(tmp_path):
    backend = PolynomialInnerBackend()
    handoff = KokunoOscillatoryNonlinearHandoff(backend)
    path = tmp_path / "nonlinear_handoff.json"
    manifest = handoff.save_manifest(path)
    rebound = KokunoOscillatoryNonlinearHandoff.load_manifest(backend, path)
    assert rebound.handoff_sha256 == handoff.handoff_sha256
    assert manifest["payload"]["scientific_boundary"]["paper_exact"] is False

    raw = json.loads(path.read_text())
    raw["payload"]["scientific_boundary"]["paper_exact"] = True
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="checksum mismatch"):
        KokunoOscillatoryNonlinearHandoff.load_manifest(backend, path)
