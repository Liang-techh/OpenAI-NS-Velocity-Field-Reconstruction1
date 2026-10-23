from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_swirl_capacity import (
    audit_supported_swirl_capacity,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
PARENT_SHA = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
CHILD_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"


def _child() -> Eq45SupportedVelocityCandidate:
    return Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.load_json(SEED))


def test_supported_swirl_capacity_is_independent_stable_and_swirl_selective():
    report = audit_supported_swirl_capacity(_child())

    assert report["parent_sha256"] == PARENT_SHA
    assert report["child_sha256"] == CHILD_SHA
    assert report["numerical_rank"] == 2
    assert 4.0 < report["condition_number"] < 6.0
    assert report["response_refinement_relative_change"] < 5e-4
    assert 0.50 < report["f12_novelty_outside_f10_span"] < 0.60

    singular = report["singular_values"]
    assert 1.30 < singular[0] < 1.40
    assert 0.27 < singular[1] < 0.30

    norms = report["response_norms"]
    assert 1.20 < norms["(1, 0)"] < 1.35
    assert 0.50 < norms["(1, 2)"] < 0.60

    for mode in ("(1, 0)", "(1, 2)"):
        cross = report["swirl_poloidal_cross_talk"][mode]
        assert cross["swirl_response_fraction"] > 1.0 - 1e-12
        assert cross["poloidal_cross_talk_fraction"] < 1e-12
        assert report["plateau_supported_vs_parent_response_relative_difference"][mode] < 1e-12

        collar = report["collar_supported_vs_parent_response_relative_difference"][mode]
        assert 0.40 < collar["radial_collar"] < 0.55
        assert 0.25 < collar["axial_collar"] < 0.40

    shares = report["region_response_share"]
    assert shares["(1, 0)"]["radial_collar"] > 0.90
    assert shares["(1, 2)"]["radial_collar"] > 0.90
    assert shares["(1, 2)"]["axial_collar"] > shares["(1, 0)"]["axial_collar"]


def test_supported_swirl_capacity_rejects_bad_steps_and_times():
    child = _child()
    with pytest.raises(ValueError):
        audit_supported_swirl_capacity(child, coefficient_steps=(0.02, 0.04))
    with pytest.raises(ValueError):
        audit_supported_swirl_capacity(child, coefficient_steps=(0.02, 0.02))
    with pytest.raises(ValueError):
        audit_supported_swirl_capacity(child, times=(0.20, 0.50))


def test_supported_swirl_capacity_does_not_mutate_or_promote_candidate():
    child = _child()
    before = child.sha256
    report = audit_supported_swirl_capacity(child)

    assert child.sha256 == before
    assert report["parameter_count_added"] == 0
    assert report["existing_parameters_audited"] == 2
    assert report["interpretation"]["new_basis_required_by_this_audit"] is False
    assert report["interpretation"]["public_visual_target_fitted"] is False

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["new_basis_added"] is False
    assert truth["forcing_or_pressure_refit"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False
