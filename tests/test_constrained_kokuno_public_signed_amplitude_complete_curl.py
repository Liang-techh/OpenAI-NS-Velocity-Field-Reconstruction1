from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    SignedRadialAmplitudeCorrection,
    correction_velocity,
    correction_velocity_dt,
    correction_vector_potential,
    profile_from_delta_a,
    truth_boundary,
    verification_receipt,
)


def _profile(scale: float = 1.0) -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.24, 1.26, 13)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(np.pi * s) ** 6
    delta = np.stack(
        (scale * 0.01 * bump, -scale * 0.007 * bump),
        axis=-1,
    )
    delta[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta,
        reference_time=0.5,
        producer_kind="focused-regression",
        provenance="test-only typed signed radial profile",
        source_mean_amplitude_differential_certified=False,
    )


def test_public_profile_builder_is_typed_and_fail_closed() -> None:
    signature = inspect.signature(profile_from_delta_a)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "forcing",
        "pressure",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(signature.parameters)

    radii = np.linspace(0.24, 1.26, 13)
    bad = np.zeros((13, 2))
    bad[0, 0] = 1.0e-3
    with pytest.raises(ValueError, match="endpoints"):
        profile_from_delta_a(
            radii,
            bad,
            reference_time=0.5,
            producer_kind="bad",
            provenance="bad noncompact profile",
        )

    with pytest.raises(ValueError):
        profile_from_delta_a(
            radii[::-1],
            np.zeros((13, 2)),
            reference_time=0.5,
            producer_kind="bad",
            provenance="bad radii",
        )


def test_zero_profile_is_exact_zero_for_potential_velocity_and_dt() -> None:
    profile = _profile(scale=0.0)
    x = np.asarray((0.42, 0.71, 0.98))
    y = np.asarray((0.13, -0.19, 0.21))
    z = np.asarray((-0.4, 0.2, 0.7))
    t = np.asarray((0.37, 0.5, 0.63))

    for evaluator in (
        correction_vector_potential,
        correction_velocity,
        correction_velocity_dt,
    ):
        values = evaluator(profile, x, y, z, t)
        assert values.shape == (3, 3)
        assert np.array_equal(values, np.zeros_like(values))


def test_nonzero_profile_is_finite_nontrivial_and_compact() -> None:
    profile = _profile()
    x = np.asarray((0.45, 0.72, 0.99))
    y = np.asarray((0.11, -0.17, 0.09))
    z = np.asarray((-0.5, 0.0, 0.6))
    t = np.asarray((0.39, 0.5, 0.61))

    velocity = correction_velocity(profile, x, y, z, t)
    velocity_dt = correction_velocity_dt(profile, x, y, z, t)
    potential = correction_vector_potential(profile, x, y, z, t)
    assert np.all(np.isfinite(velocity))
    assert np.all(np.isfinite(velocity_dt))
    assert np.all(np.isfinite(potential))
    assert np.max(np.linalg.norm(velocity, axis=-1)) > 1.0e-8
    assert np.max(np.linalg.norm(potential, axis=-1)) > 1.0e-8

    exterior = correction_velocity(
        profile,
        np.asarray((0.0, 0.18, 1.31, 0.75)),
        np.asarray((0.0, 0.0, 0.0, 0.0)),
        np.asarray((0.0, 0.0, 0.0, 2.05)),
        np.full(4, 0.5),
    )
    assert np.array_equal(exterior, np.zeros_like(exterior))


def test_truth_boundary_keeps_source_and_repository_claims_separate() -> None:
    truth = truth_boundary()
    assert truth["source_complete_curl_product_rule_reused"] is True
    assert truth["vector_potential_first_contract"] is True
    assert truth["agent3_mean_radial_chain_reimplemented"] is False
    assert truth["shared_across_beta_profile_is_repository_realization"] is True
    assert truth["quintic_profile_lift_is_repository_realization"] is True
    assert truth["autonomous_time_lift_is_repository_realization"] is True
    assert truth["source_agent2_complete_curl_certified_for_full_candidate"] is False
    assert truth["independent_agent4_correction_curl_audit_required"] is True
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_preregistered_complete_curl_receipt() -> None:
    receipt = verification_receipt()
    assert receipt["sample_count"] == 18
    assert receipt["correction"]["source_mean_amplitude_differential_certified"] is False
    assert receipt["failed_guards"] == []
    assert min(receipt["curl_rms_refinement_ratios"]) >= 8.0
    assert min(receipt["time_rms_refinement_ratios"]) >= 20.0
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12
    assert receipt["compact_c2_boundary_max_abs"] <= 1.0e-12
    assert receipt["truth_boundary"]["pde_validated"] is False
