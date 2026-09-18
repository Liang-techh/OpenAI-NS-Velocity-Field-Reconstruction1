import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy import (
    KokunoHeatReplacementDiscrepancy,
)


def _tail_normalized(model: KokunoHeatReplacementDiscrepancy, eta: float) -> np.ndarray:
    scales = model.reference_tail_scales()
    e_k = scales["e_K"]
    X_k = scales["X_K"]
    return model.discrepancy(eta) / np.asarray(
        [e_k**2, X_k * e_k**2, X_k**1.5 * e_k]
    )


def test_source_tail_geometry_and_terminal_guard_are_executable():
    model = KokunoHeatReplacementDiscrepancy()
    np.testing.assert_allclose(model.chi_K(np.asarray([0.1, 0.2])), 0.0, atol=0.0)
    np.testing.assert_allclose(model.chi_K(np.asarray([0.5, 0.7])), 1.0, atol=0.0)
    assert model.X_K == pytest.approx(np.exp(0.2) * model.X_tail)
    assert model.X_b == pytest.approx(np.exp(3.0) * model.X_tail)
    assert 0.0 <= model.max_terminal_log_slope() < model.h / 4.0
    with pytest.raises(ValueError, match="terminal bound"):
        KokunoHeatReplacementDiscrepancy(rho_o=0.001)


def test_actual_heat_replacement_discrepancy_is_even_finite_and_endpoint_zero():
    model = KokunoHeatReplacementDiscrepancy()
    eta = np.asarray([-1.0, -0.2, 0.0, 0.2, 1.0])
    discrepancy = model.discrepancy(eta)
    assert discrepancy.shape == (5, 3)
    assert np.all(np.isfinite(discrepancy))
    np.testing.assert_allclose(discrepancy[0], 0.0, atol=0.0)
    np.testing.assert_allclose(discrepancy[-1], 0.0, atol=0.0)
    np.testing.assert_allclose(discrepancy[1], discrepancy[3], rtol=2.0e-11, atol=2.0e-12)
    np.testing.assert_allclose(
        discrepancy[3],
        np.asarray([-2.20937473e-09, 6.25423006e-06, -2.63159207]),
        rtol=2.0e-6,
        atol=2.0e-10,
    )
    assert discrepancy[3, 0] < 0.0
    assert discrepancy[3, 1] > 0.0
    assert discrepancy[3, 2] < 0.0


def test_source_large_scale_normalized_discrepancy_has_expected_inverse_scale():
    first = KokunoHeatReplacementDiscrepancy(X_tail=1000.0)
    second = KokunoHeatReplacementDiscrepancy(X_tail=2000.0)
    d1 = _tail_normalized(first, 0.2)
    d2 = _tail_normalized(second, 0.2)
    ratio = np.abs(d2 / d1)
    assert np.all((ratio > 0.49) & (ratio < 0.51))


def test_normalized_target_exactly_cancels_physical_discrepancy_scaling():
    model = KokunoHeatReplacementDiscrepancy()
    eta = np.asarray([-0.4, 0.2])
    X_star = 17.0
    e_star = 0.3
    target = model.normalized_repair_target(eta, X_star=X_star, e_star=e_star)
    physical_scales = np.asarray(
        [e_star**2, X_star * e_star**2, X_star**1.5 * e_star]
    )
    np.testing.assert_allclose(
        target * physical_scales + model.discrepancy(eta),
        0.0,
        rtol=0.0,
        atol=2.0e-14,
    )


def test_quadrature_stability_and_serialization_truth_are_fail_closed(tmp_path):
    base = KokunoHeatReplacementDiscrepancy()
    refined = KokunoHeatReplacementDiscrepancy(
        heat_quadrature_order=128,
        transition_quadrature_order=256,
        tail_quadrature_order=128,
    )
    np.testing.assert_allclose(
        base.discrepancy(0.37), refined.discrepancy(0.37), rtol=4.0e-7, atol=2.0e-10
    )

    destination = base.save_json(tmp_path / "heat_discrepancy.json")
    loaded = KokunoHeatReplacementDiscrepancy.load_json(destination)
    assert loaded.sha256 == base.sha256
    assert loaded.to_payload() == base.to_payload()

    payload = json.loads(destination.read_text(encoding="utf-8"))
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoHeatReplacementDiscrepancy.from_payload(payload)
