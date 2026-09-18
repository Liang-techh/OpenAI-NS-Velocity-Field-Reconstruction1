from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_coupled_normalization_pa10_screen import (
    KokunoCoupledNormalizationPA10Screen,
)
from openai_ns_reconstruction.kokuno_rescaled_inner_join_tsh_certificate import (
    KokunoRescaledInnerJoinTshCertificate,
)


@pytest.fixture(scope="module")
def selected_certificate() -> KokunoRescaledInnerJoinTshCertificate:
    return KokunoRescaledInnerJoinTshCertificate()


@pytest.fixture(scope="module")
def selected_screen(selected_certificate) -> KokunoCoupledNormalizationPA10Screen:
    return KokunoCoupledNormalizationPA10Screen.from_selected_certificate(
        selected_certificate
    )


def test_selected_obstruction_becomes_explicit_upstream_cancellation_band(
    selected_certificate, selected_screen
):
    geometry = selected_screen.geometry_report()
    selected = selected_screen.selected_realization_report(selected_certificate)

    assert geometry["t_sh_strict_upper"] == pytest.approx(
        selected_certificate.max_T_sh_for_separation, rel=2e-16
    )
    assert geometry["pointwise_feasible_band_nonempty"] is True
    assert geometry["ell_abs_strict_cap"] > 0.0

    lo, hi = geometry["required_log_E_i_open_interval"]
    assert lo < -selected_screen.log_C < hi
    assert selected["selected_log_E_i"] > hi
    assert selected["selected_realization_excluded"] is True
    assert selected["evaluation"]["strict_margin"] < 0.0
    assert selected["selected_abs_ell_over_strict_cap"] == pytest.approx(16.0)
    assert selected["source_T_sh_lower_bound_verified"] is False


def test_perfect_cancellation_is_only_not_excluded_not_certified():
    screen = KokunoCoupledNormalizationPA10Screen(
        log_C=100.0,
        log_P_star=5.0,
    )
    result = screen.evaluation_report(-100.0)
    assert result["ell_i"] == pytest.approx(0.0, abs=0.0)
    assert result["strict_margin"] > 0.0
    assert result["not_excluded_by_necessary_screen"] is True
    assert screen.geometry_report()["source_T_sh_lower_bound_verified"] is False


def test_vectorized_screen_rejects_values_outside_open_band():
    screen = KokunoCoupledNormalizationPA10Screen(
        log_C=100.0,
        log_P_star=5.0,
    )
    cap = screen.ell_abs_strict_cap
    values = np.array(
        [
            -screen.log_C,
            -screen.log_C + 0.5 * cap,
            -screen.log_C + 1.1 * cap,
            -screen.log_C - 1.1 * cap,
        ]
    )
    result = screen.evaluate_log_E_i(values)
    np.testing.assert_array_equal(
        result["not_excluded_by_necessary_screen"],
        np.array([True, True, False, False]),
    )
    assert np.all(result["strict_margin"][:2] > 0.0)
    assert np.all(result["strict_margin"][2:] < 0.0)


def test_empty_band_and_serialization_truth_fail_closed(selected_screen):
    impossible = KokunoCoupledNormalizationPA10Screen(
        log_C=1.0,
        log_P_star=-1.0,
    )
    assert impossible.pointwise_feasible_band_nonempty is False
    with pytest.raises(ValueError, match="no pointwise PA.10-feasible"):
        impossible.required_log_E_i_open_interval()

    payload = selected_screen.to_payload()
    replay = KokunoCoupledNormalizationPA10Screen.from_payload(payload)
    assert replay.to_payload() == payload

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_T_sh_lower_bound_verified"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoCoupledNormalizationPA10Screen.from_payload(tampered)

    assert len(selected_screen.sha256) == 64
