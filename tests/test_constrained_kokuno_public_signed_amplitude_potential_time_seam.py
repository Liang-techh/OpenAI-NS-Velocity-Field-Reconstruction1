from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    evaluate_signed_amplitude_correction,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_potential_time_seam import (
    correction_vector_potential_dt,
    truth_boundary,
    verification_receipt,
    _verification_profile,
)


def test_public_vector_potential_dt_replays_parent_materialization_exactly() -> None:
    correction = _verification_profile()
    x = np.asarray((0.43, -0.61, 0.82), dtype=float)
    y = np.asarray((0.17, 0.26, -0.31), dtype=float)
    z = np.asarray((-0.41, 0.19, 0.55), dtype=float)
    t = np.asarray((0.38, 0.50, 0.62), dtype=float)
    public = correction_vector_potential_dt(correction, x, y, z, t)
    parent = evaluate_signed_amplitude_correction(correction, x, y, z, t)[
        "vector_potential_dt_cartesian_total"
    ]
    assert np.array_equal(public, parent)


def test_potential_time_commutation_receipt_passes_preregistered_guards() -> None:
    receipt = verification_receipt()
    assert receipt["schema"] == "kokuno-a2-signed-amplitude-potential-time-seam-v1"
    assert receipt["sample_count"] == 18
    assert receipt["failed_guards"] == []
    assert min(receipt["potential_time_refinement_ratios"]) >= 20.0
    assert min(receipt["curl_rms_refinement_ratios"]) >= 8.0
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12


def test_truth_boundary_keeps_scientific_promotions_closed() -> None:
    boundary = truth_boundary()
    assert boundary["public_vector_potential_time_derivative_executable"] is True
    assert boundary["potential_time_commutation_auditable"] is True
    assert boundary["curl_time_commutation_auditable"] is True
    assert boundary["vector_potential_first_complete_curl_reused"] is True
    assert boundary["candidate_amplitude_phase_support_retuned"] is False
    assert boundary["agent3_mean_radial_chain_reimplemented"] is False
    assert boundary["source_agent2_complete_curl_certified_for_full_candidate"] is False
    assert boundary["independent_agent4_correction_vector_potential_audit_required"] is True
    assert boundary["real_agent3_delta_a_bound"] is False
    assert boundary["heldout_ns_momentum_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["paper_exact"] is False


def test_public_signature_has_no_residual_or_forcing_escape_hatch() -> None:
    parameters = set(inspect.signature(correction_vector_potential_dt).parameters)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
    }
    assert parameters.isdisjoint(forbidden)


def test_parent_time_domain_guard_is_inherited() -> None:
    correction = _verification_profile()
    with pytest.raises(ValueError):
        correction_vector_potential_dt(correction, 0.5, 0.0, 0.0, 0.20)
