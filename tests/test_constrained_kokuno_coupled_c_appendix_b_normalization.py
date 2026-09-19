from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_coupled_c_appendix_b_normalization import (
    KokunoCoupledCAppendixBNormalization,
    KokunoCoupledCRescaledCoreSeed,
)


def test_log_c_oversize_is_applied_before_profile_evaluation():
    seed = KokunoCoupledCRescaledCoreSeed(log_C_factor=4.0)
    eta = seed.phase_stationary_eta
    minimal = seed.minimal_real_axis_log_C
    assert seed.log_C_real_axis == pytest.approx(4.0 * minimal, rel=2e-16)
    assert seed.log_C_oversize == pytest.approx(3.0 * minimal, rel=2e-16)
    assert float(np.asarray(seed.log_g(eta))) == pytest.approx(
        -seed.log_C_oversize, rel=2e-15
    )


@pytest.fixture(scope="module")
def coupled_candidate() -> KokunoCoupledCAppendixBNormalization:
    return KokunoCoupledCAppendixBNormalization(log_C_factor=32.0)


def test_dyadic_oversize_repropagates_to_a_nonexcluded_pointwise_pa10_witness(
    coupled_candidate,
):
    report = coupled_candidate.pointwise_report()
    assert report["log_C_factor"] == pytest.approx(32.0)
    assert report["selected_log_C"] > report["minimal_real_axis_log_C"]
    assert abs(report["ell_i"]) < report["ell_abs_strict_cap"]
    assert report["strict_margin"] > 0.0
    assert report["not_excluded_by_necessary_screen"] is True
    assert abs(report["G_i"]) > 1.0e-3
    assert report["source_B0_analytic_bound_proved"] is False
    assert report["source_T_sh_lower_bound_verified"] is False
    assert report["selected_pa16_handoff_allowed"] is False


def test_velocity_stays_nontrivial_and_payload_truth_boundary_fails_closed(
    coupled_candidate,
):
    value = np.asarray(coupled_candidate.velocity(0.0, 0.0, 0.0, 0.0), dtype=float)
    assert value.shape == (3,)
    assert np.all(np.isfinite(value))
    assert np.linalg.norm(value) > 1.0e-6

    payload = coupled_candidate.to_payload()
    replay = KokunoCoupledCAppendixBNormalization.from_payload(payload)
    assert replay.to_payload() == payload
    assert len(coupled_candidate.sha256) == 64

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_T_sh_lower_bound_verified"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoCoupledCAppendixBNormalization.from_payload(tampered)
