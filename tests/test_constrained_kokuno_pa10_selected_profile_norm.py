from __future__ import annotations

import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_selected_profile_norm import (
    KokunoPA10SelectedProfileNormEnvelope,
)


def _fast() -> KokunoPA10SelectedProfileNormEnvelope:
    return KokunoPA10SelectedProfileNormEnvelope(
        coarse_eta_points=17,
        fine_eta_points=33,
        safety_fraction=1.0e-8,
        max_refinement_drift_fraction=1.0e-5,
    )


def test_selected_combined_profile_is_vectorized_finite_and_nontrivial() -> None:
    model = _fast()
    eta = np.asarray([-1.0, -0.25, 0.0, model.seed.phase_stationary_eta, 0.5, 1.0])
    values = model.combined_log_profile(eta)
    derivative = model.combined_log_profile_eta(eta)
    assert values.shape == eta.shape
    assert derivative.shape == eta.shape
    assert np.all(np.isfinite(values))
    assert np.all(np.isfinite(derivative))
    assert np.max(values) - np.min(values) > 1.0e10


def test_stable_formula_avoids_logC_logF_cancellation_at_X0() -> None:
    model = _fast()
    eta = np.asarray([-1.0, -0.5, 0.0, 0.5, 1.0])
    center = model.seed.profile_values_scaled(np.full(eta.shape, 4.0), eta)
    # For the selected real-axis normalization, log(phi)=log(C)+log(F), but
    # that representation catastrophically cancels the O(1) center value at
    # eta=0 because log(C) and log(F) are both O(1e27).  The direct source
    # rescaling used by combined_log_profile must retain that information.
    replay = model.seed.log_C_real_axis + center["log_F"]
    stable = model.combined_log_profile(eta)

    large = np.abs(stable) > 1.0e20
    scale = np.abs(stable[large])
    assert np.all(np.abs(replay[large] - stable[large]) / scale < 5.0e-15)

    zero_index = 2
    assert replay[zero_index] == 0.0
    assert stable[zero_index] != 0.0
    expected_zero = np.log(model.seed.f0(4.0 * model.seed.chi(np.asarray(0.0))))
    assert stable[zero_index] == pytest.approx(float(expected_zero), rel=0.0, abs=1.0e-15)


def test_analytic_eta_derivative_agrees_with_centered_difference() -> None:
    model = _fast()
    eta = 0.23
    h = 2.0e-6
    fd = float(
        (model.combined_log_profile(eta + h) - model.combined_log_profile(eta - h))
        / (2.0 * h)
    )
    exact = float(model.combined_log_profile_eta(np.asarray(eta)))
    assert np.isfinite(fd)
    assert abs(fd - exact) / max(1.0, abs(exact)) < 2.0e-7


def test_envelope_is_padded_stable_and_fail_closed() -> None:
    model = _fast()
    envelope = model.envelope
    fine = envelope["fine"]
    assert envelope["combined_profile_C0_candidate_upper_envelope"] > fine[
        "sample_max_abs_combined_log_profile"
    ]
    assert envelope["nested_grid_relative_drift"] <= envelope[
        "max_allowed_refinement_drift_fraction"
    ]
    assert envelope["numerical_envelope_guard_passed"] is True
    assert envelope["continuum_supremum_proved"] is False
    report = model.report(include_selected_B0_diagnostic=False)
    truth = report["truth_boundary"]
    assert truth["selected_center_combined_profile_executable"] is True
    assert truth["selected_center_combined_profile_eta_derivative_executable"] is True
    assert truth["selected_center_is_source_fixed_point"] is False
    assert truth["selected_profile_norm_envelope_is_continuum_interval_proof"] is False
    assert truth["combined_profile_C0_norm_machine_bound"] is False
    assert truth["source_B0_dependencies_machine_bound"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_candidate_B0_diagnostic_consumes_explicit_engineering_inputs_only() -> None:
    model = _fast()
    value = model.candidate_B0_diagnostic(M0=4.0, L0=66.0)
    profile = model.envelope["combined_profile_C0_candidate_upper_envelope"]
    # The additive O(1e2) correction is below one binary64 ulp at the selected
    # O(1e27) profile scale, so equality is an expected floating-point outcome.
    assert value >= profile
    assert np.isfinite(value)


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path) -> None:
    model = _fast()
    payload = model.to_payload()
    replay = KokunoPA10SelectedProfileNormEnvelope.from_payload(payload)
    assert replay.sha256 == model.sha256
    assert replay.to_payload() == payload

    path = tmp_path / "profile_norm.json"
    model.save_json(path)
    loaded = KokunoPA10SelectedProfileNormEnvelope.load_json(path)
    assert loaded.sha256 == model.sha256
    assert json.loads(path.read_text())["sha256"] == model.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_B0_dependencies_machine_bound"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoPA10SelectedProfileNormEnvelope.from_payload(tampered)


def test_domain_rejects_invalid_eta() -> None:
    model = _fast()
    with pytest.raises(ValueError, match="eta must lie"):
        model.combined_log_profile(np.asarray([1.001]))
    with pytest.raises(ValueError, match="eta must lie"):
        model.combined_log_profile_eta(np.asarray([-1.001]))
