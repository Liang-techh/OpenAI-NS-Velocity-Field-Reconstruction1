from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_nonlinear_attribution import (
    KokunoOscillatoryNonlinearAttribution,
    public_contract,
)


class PolynomialInnerBackend:
    field_sha256 = "2" * 64

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
        return np.stack((0.21 * x + 0.09, -0.13 * y + 0.04, 0.11 * z - 0.03), axis=-1)

    def velocity_dt(self, x, y, z, t):
        x, _, _, _ = self._broadcast(x, y, z, t)
        return np.zeros(x.shape + (3,), dtype=float)

    def self_advection(self, x, y, z, t):
        u = self.velocity(x, y, z, t)
        return np.stack((0.21 * u[..., 0], -0.13 * u[..., 1], 0.11 * u[..., 2]), axis=-1)

    def velocity_laplacian(self, x, y, z, t):
        return self.velocity_dt(x, y, z, t)


def _probes():
    return (
        np.asarray([0.41, 0.57, 0.72, 0.88, 1.03]),
        np.asarray([0.16, -0.20, 0.25, -0.17, 0.12]),
        np.asarray([-0.43, -0.11, 0.21, 0.49, 0.66]),
        np.asarray([0.472, 0.486, 0.500, 0.514, 0.528]),
    )


def _vrms(value):
    value = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(value * value, axis=-1))))


def _canonical_sha(payload):
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def test_attribution_replays_parent_magnitudes_and_normalized_shares():
    attribution = KokunoOscillatoryNonlinearAttribution(PolynomialInnerBackend())
    parent = attribution._parent.evaluate(*_probes())
    result = attribution.evaluate(*_probes())

    expected = [
        _vrms(parent.inner_advects_oscillation),
        _vrms(parent.oscillation_advects_inner),
        _vrms(parent.oscillatory_self_advection),
    ]
    denominator = sum(expected)
    assert denominator > 0.0
    assert result.inner_advects_oscillation_rms == pytest.approx(expected[0])
    assert result.oscillation_advects_inner_rms == pytest.approx(expected[1])
    assert result.oscillatory_self_advection_rms == pytest.approx(expected[2])
    assert result.mixed_cross_advection_rms == pytest.approx(_vrms(parent.mixed_cross_advection))
    assert result.aggregate_nonlinear_increment_rms == pytest.approx(
        _vrms(parent.oscillatory_nonlinear_increment)
    )
    shares = np.asarray([
        result.inner_advects_oscillation_share,
        result.oscillation_advects_inner_share,
        result.oscillatory_self_advection_share,
    ])
    assert np.all(np.isfinite(shares))
    assert np.all((shares >= 0.0) & (shares <= 1.0))
    assert float(np.sum(shares)) == pytest.approx(1.0, abs=1.0e-15)
    assert result.mixed_magnitude_share == pytest.approx(shares[0] + shares[1], abs=1.0e-15)
    assert result.decomposition_closure_abs_max <= 5.0e-13
    assert result.mixed_closure_abs_max <= 5.0e-13


def test_contract_is_diagnostic_not_mean_projection_or_residual_gate():
    contract = public_contract()
    assert contract["forbidden_evaluate_inputs_present"] == []
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["shares_are_acceptance_gates"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    payload = KokunoOscillatoryNonlinearAttribution(PolynomialInnerBackend()).semantic_payload()
    metrics = payload["metric_contract"]
    assert metrics["shares_are_non_cancelling_magnitude_diagnostics"] is True
    assert metrics["shares_are_ns_residual_fractions"] is False
    assert metrics["shares_are_m0_mean_defect_fractions"] is False
    assert metrics["shares_are_acceptance_gates"] is False
    boundary = payload["scientific_boundary"]
    assert boundary["mean_projection_performed"] is False
    assert boundary["same_protocol_st006_comparison_valid"] is False
    assert boundary["residual_reduction_claimed"] is False


def test_manifest_rebind_and_checksum_valid_semantic_tamper_fail_closed(tmp_path):
    backend = PolynomialInnerBackend()
    attribution = KokunoOscillatoryNonlinearAttribution(backend)
    path = tmp_path / "nonlinear_attribution.json"
    attribution.save_manifest(path)
    rebound = KokunoOscillatoryNonlinearAttribution.load_manifest(backend, path)
    assert rebound.attribution_sha256 == attribution.attribution_sha256

    raw = json.loads(path.read_text())
    raw["payload"]["scientific_boundary"]["mean_projection_performed"] = True
    raw["attribution_sha256"] = _canonical_sha(raw["payload"])
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="semantics mismatch"):
        KokunoOscillatoryNonlinearAttribution.load_manifest(backend, path)
