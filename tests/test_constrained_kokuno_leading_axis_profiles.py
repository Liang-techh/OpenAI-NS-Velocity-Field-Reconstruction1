import json
from dataclasses import fields

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_leading_axis_profiles import (
    CORRECTED_RELEASE,
    KokunoLeadingAxisProfile,
    SOURCE_COMMIT,
)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"h": 0.0},
        {"h": -1e-3},
        {"h": 0.01},
        {"h": 0.0100001},
        {"h": np.nan},
        {"h": np.inf},
        {"j0": 0.0},
        {"j0": -1e-3},
        {"j0": 0.050001},
        {"j0": np.nan},
        {"j0": np.inf},
    ],
)
def test_source_parameter_bounds_fail_closed(kwargs):
    with pytest.raises(ValueError):
        KokunoLeadingAxisProfile(**kwargs)


def test_vectorized_formulas_and_simplified_identity():
    profile = KokunoLeadingAxisProfile(h=0.006, j0=0.021)
    eta = np.array([-1.0, -0.3, 0.0, 0.4, 1.0])
    out = profile.evaluate(eta)

    expected_u = 4.0 * eta + profile.j0
    d = 1.0 - eta**2
    expected_h = profile.D * eta + d * expected_u
    expected_w = 1.0 - 4.0 * d - 2.0 * profile.D * eta * expected_u

    np.testing.assert_allclose(out["U_star"], expected_u, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(out["H_star"], expected_h, rtol=0.0, atol=2e-16)
    np.testing.assert_allclose(out["W_star"], expected_w, rtol=0.0, atol=3e-16)
    np.testing.assert_allclose(profile.w_star(eta), profile.w_star_simplified(eta), rtol=0.0, atol=8e-16)
    np.testing.assert_allclose(out["dU_star_deta"], 4.0, rtol=0.0, atol=0.0)


def test_scalar_and_vector_evaluation_agree():
    profile = KokunoLeadingAxisProfile()
    eta = 0.37
    scalar = profile.evaluate(eta)
    vector = profile.evaluate(np.array([eta]))
    for key in scalar:
        np.testing.assert_allclose(np.asarray(scalar[key]), np.asarray(vector[key])[0], rtol=0.0, atol=0.0)


def test_analytic_derivatives_match_centered_differences():
    profile = KokunoLeadingAxisProfile(h=0.004, j0=0.031)
    eta = np.array([-0.8, -0.2, 0.1, 0.65])
    eps = 1e-6

    fd_h = (profile.h_star(eta + eps) - profile.h_star(eta - eps)) / (2.0 * eps)
    fd_w = (profile.w_star(eta + eps) - profile.w_star(eta - eps)) / (2.0 * eps)

    np.testing.assert_allclose(profile.dh_star_deta(eta), fd_h, rtol=2e-9, atol=2e-9)
    np.testing.assert_allclose(profile.dw_star_deta(eta), fd_w, rtol=2e-8, atol=2e-9)


def test_positive_j0_preserves_source_asymmetry_and_nontriviality():
    profile = KokunoLeadingAxisProfile()
    assert profile.u_star(0.0) == pytest.approx(profile.j0)
    assert profile.u_star(0.0) > 0.0
    assert profile.dw_star_deta(0.0) < 0.0
    assert not np.allclose(profile.u_star(np.array([-0.5, 0.0, 0.5])), 0.0)


def test_no_amplitude_collapse_parameter_is_exposed():
    assert [field.name for field in fields(KokunoLeadingAxisProfile)] == ["h", "j0"]


def test_payload_pins_source_and_coordinate_discrepancy():
    payload = KokunoLeadingAxisProfile().to_payload()
    assert payload["source"]["commit"] == SOURCE_COMMIT
    assert payload["source"]["corrected_release"] == CORRECTED_RELEASE
    assert payload["source"]["pdf_independently_parsed"] is False
    assert payload["coordinates"]["native_time_relation"] == "tau=1-t=q(1-eta^2)=q-z^2*q^(2h)"
    assert payload["coordinates"]["current_repo_eq45_user_relation"] == "q-z^2*q^(-2h)=1-t"
    assert payload["coordinates"]["mapping_status"] == "pending_explicit_convention_reconciliation"
    assert payload["truth_boundary"]["full_3d_velocity_candidate"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False


def test_serialization_roundtrip_preserves_sha_and_values(tmp_path):
    profile = KokunoLeadingAxisProfile(h=0.007, j0=0.019)
    path = profile.save_json(tmp_path / "profile.json")
    reloaded = KokunoLeadingAxisProfile.load_json(path)
    assert reloaded == profile
    assert reloaded.sha256 == profile.sha256
    eta = np.linspace(-1.0, 1.0, 17)
    for key, value in profile.evaluate(eta).items():
        np.testing.assert_array_equal(reloaded.evaluate(eta)[key], value)


@pytest.mark.parametrize("section", ["source", "formulas", "coordinates", "truth_boundary"])
def test_serialization_rejects_truth_or_provenance_mutation(tmp_path, section):
    profile = KokunoLeadingAxisProfile()
    path = profile.save_json(tmp_path / "profile.json")
    payload = json.loads(path.read_text())
    first_key = next(iter(payload[section]))
    payload[section][first_key] = "tampered"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError):
        KokunoLeadingAxisProfile.load_json(path)


def test_serialization_rejects_sha_mutation(tmp_path):
    profile = KokunoLeadingAxisProfile()
    path = profile.save_json(tmp_path / "profile.json")
    payload = json.loads(path.read_text())
    payload["sha256"] = "0" * 64
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="SHA mismatch"):
        KokunoLeadingAxisProfile.load_json(path)
