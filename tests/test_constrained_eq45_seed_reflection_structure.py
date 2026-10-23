from dataclasses import replace

import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_seed_reflection_structure import (
    FROZEN_CANDIDATE_SHA256,
    PUBLIC_CROSSCHECK_TOLERANCE,
    audit_frozen_eq45_seed_reflection,
    certify_even_eta_seed,
    public_velocity_reflection_crosscheck,
)


def test_frozen_eq45_seed_has_strict_reflection_certificate():
    report = audit_frozen_eq45_seed_reflection()

    assert report["schema"] == "eq45_seed_reflection_structure_v1"
    assert report["status"] == "strict_seed_reflection_structure_verified"
    assert report["candidate_sha256"] == FROZEN_CANDIDATE_SHA256

    strict = report["strict_evidence"]
    assert strict["odd_eta_coefficients_exactly_zero"] is True
    assert strict["odd_eta_modes"] == [
        {"mode": [0, 1], "phi": 0.0, "F": 0.0},
        {"mode": [1, 1], "phi": 0.0, "F": 0.0},
    ]
    assert strict["derived_profile_parity"]["v0"] == "odd_in_eta"
    assert strict["derived_profile_parity"]["U"] == "even_in_eta"
    assert strict["derived_physical_channel_parity"]["u_r"] == "odd_in_z"
    assert strict["derived_physical_channel_parity"]["u_theta"] == "even_in_z"
    assert strict["derived_physical_channel_parity"]["w"] == "even_in_z"
    assert (
        strict["derived_physical_channel_parity"]["midplane_radial_velocity"]
        == "exactly_zero_at_z_equals_0"
    )

    crosscheck = report["public_velocity_crosscheck"]
    assert crosscheck["tolerance"] == PUBLIC_CROSSCHECK_TOLERANCE
    for name in (
        "radial_odd_max_abs",
        "swirl_even_max_abs",
        "axial_even_max_abs",
        "midplane_radial_max_abs",
    ):
        assert crosscheck[name] <= PUBLIC_CROSSCHECK_TOLERANCE

    states = report["states"]
    assert states["velocity_export_ready"] is True
    for name in (
        "visualization_ready",
        "visual_correspondence_verified",
        "physical_support_validated",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert states[name] is False


def test_odd_eta_mutation_invalidates_strict_certificate():
    candidate = Eq45VelocityCandidate.load_json(
        "artifacts/constrained/eq45_velocity_candidate_seed.json"
    )
    basis = candidate.profile_basis
    coefficients = list(basis.phi_coefficients)
    odd_index = basis.mode_indices.index((0, 1))
    coefficients[odd_index] = 0.125
    mutant_basis = replace(basis, phi_coefficients=tuple(coefficients))
    mutant = replace(candidate, profile_basis=mutant_basis)

    with pytest.raises(ValueError, match="odd-eta profile coefficients"):
        certify_even_eta_seed(mutant)


def test_public_reflection_crosscheck_fails_closed_on_invalid_tolerance():
    candidate = Eq45VelocityCandidate.load_json(
        "artifacts/constrained/eq45_velocity_candidate_seed.json"
    )
    with pytest.raises(ValueError, match="tolerance"):
        public_velocity_reflection_crosscheck(candidate, tolerance=0.0)
