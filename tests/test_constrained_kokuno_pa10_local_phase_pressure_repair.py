import numpy as np

from openai_ns_reconstruction.kokuno_pa10_selected_pressure_primitive import (
    KokunoPA10SelectedPressurePrimitive,
)


def _fd4_eta(primitive, field, y, eta, h):
    def values(offset):
        return np.asarray(primitive.evaluate(y, eta + offset * h)[field], dtype=float)

    return (-values(2.0) + 8.0 * values(1.0) - 8.0 * values(-1.0) + values(-2.0)) / (
        12.0 * h
    )


def test_stationary_local_phase_resolves_lambda_minus_half_layer():
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=64)
    root = float(primitive.core.phase_stationary_eta)
    width = float(primitive.core.rescaling_lambda ** -0.5)
    offsets = np.asarray([-0.60, -0.30, 0.0, 0.30, 0.60], dtype=float)
    eta = root + offsets * width

    log_g = np.asarray(primitive.evaluate(np.ones_like(eta), eta)["log_g"], dtype=float)

    assert log_g[2] == 0.0
    assert log_g[0] < -1.0e-8
    assert log_g[1] < -1.0e-8
    assert log_g[3] < -1.0e-8
    assert log_g[4] < -1.0e-8
    assert np.ptp(log_g) > 1.0e-6


def test_stationary_anchored_zeta_is_zero_at_selected_root():
    primitive = KokunoPA10SelectedPressurePrimitive()
    root = float(primitive.core.phase_stationary_eta)
    got = primitive.evaluate(1.7, root)
    assert float(got["zeta_star_local"]) == 0.0
    assert float(got["log_g"]) == 0.0


def test_pressure_eta_matches_value_only_fd4_on_natural_layer():
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=64)
    root = float(primitive.core.phase_stationary_eta)
    width = float(primitive.core.rescaling_lambda ** -0.5)
    y = np.asarray([0.73, 1.61, 2.84, 3.70], dtype=float)
    eta = root + np.asarray([-0.45, -0.20, 0.25, 0.50], dtype=float) * width

    declared = np.asarray(primitive.evaluate(y, eta)["p_eta"], dtype=float)
    coarse = _fd4_eta(primitive, "p", y, eta, 0.10 * width)
    fine = _fd4_eta(primitive, "p", y, eta, 0.05 * width)

    def relative_rms(estimate):
        error = np.asarray(estimate) - declared
        return float(np.sqrt(np.mean(error * error)) / np.sqrt(np.mean(declared * declared)))

    coarse_error = relative_rms(coarse)
    fine_error = relative_rms(fine)
    assert fine_error < 5.0e-3
    assert fine_error < coarse_error


def test_payload_records_precision_safe_value_path(tmp_path):
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=32)
    payload = primitive.to_payload()
    assert payload["schema"] == "kokuno-pa10-selected-pressure-primitive-v2"
    numerics = payload["numerics"]
    assert numerics["stationary_local_phase_quadrature_order"] == 64
    assert numerics["global_phase_subtraction_used_for_public_pressure_values"] is False
    assert primitive.truth_boundary["selected_stationary_binary64_root_is_source_exact_hidden_value"] is False

    path = tmp_path / "selected-pressure-v2.json"
    primitive.save(path)
    restored = KokunoPA10SelectedPressurePrimitive.load(path)
    assert restored.sha256 == primitive.sha256
