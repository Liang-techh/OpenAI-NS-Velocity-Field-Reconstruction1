import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_reference_plateau_x100 import (
    SCHEMA,
    KokunoPA10ReferencePlateauToX100,
)


def test_reference_plateau_extends_parent_exactly_and_stays_nontrivial():
    profile = KokunoPA10ReferencePlateauToX100()
    eta = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])

    inherited_x = np.full_like(eta, 0.5 * (profile.X_0 + profile.X_2))
    got = profile.values(inherited_x, eta)
    parent = profile.collar.values(inherited_x, eta)
    np.testing.assert_array_equal(got["F_reference"], parent["F_reference"])
    np.testing.assert_array_equal(got["U_reference"], parent["U_reference"])
    np.testing.assert_array_equal(got["E_reference"], parent["E_reference"])

    endpoint = profile.collar.values(np.full_like(eta, profile.X_2), eta)
    for X in (2.0 * profile.X_2, 1.0, 10.0, 50.0, profile.X_b_ref):
        values = profile.values(np.full_like(eta, X), eta)
        np.testing.assert_array_equal(values["F_reference"], endpoint["F_reference"])
        np.testing.assert_array_equal(values["U_reference"], endpoint["U_reference"])
        np.testing.assert_allclose(
            values["E_reference"],
            np.sqrt(2.0 * X) * endpoint["F_reference"],
            rtol=0.0,
            atol=0.0,
        )

    handoff = profile.handoff_at_X100(eta)
    assert np.all(np.isfinite(handoff["F_reference"]))
    assert np.all(np.isfinite(handoff["U_reference"]))
    assert np.all(np.isfinite(handoff["E_reference"]))
    assert np.linalg.norm(handoff["F_reference"]) > 1.0e-12
    assert np.linalg.norm(handoff["U_reference"]) > 1.0e-12


def test_plateau_radial_jets_and_public_log_slope_are_exact():
    profile = KokunoPA10ReferencePlateauToX100()
    eta = np.asarray([-0.6, -0.2, 0.0, 0.2, 0.6])
    X = np.asarray([0.5, 1.0, 5.0, 20.0, 99.0])
    assert np.all(X > profile.X_2)

    derivatives = profile.radial_derivatives(X, eta)
    slopes = profile.radial_log_slopes(X, eta)
    values = profile.values(X, eta)

    np.testing.assert_array_equal(derivatives["F_reference_X"], np.zeros_like(X))
    np.testing.assert_array_equal(derivatives["U_reference_X"], np.zeros_like(X))
    np.testing.assert_array_equal(slopes["D_X_log_F_reference"], np.zeros_like(X))
    np.testing.assert_array_equal(slopes["D_X_U_reference"], np.zeros_like(X))
    np.testing.assert_array_equal(slopes["l_reference"], np.ones_like(X))
    np.testing.assert_allclose(
        derivatives["E_reference_X"],
        values["F_reference"] / np.sqrt(2.0 * X),
        rtol=5.0e-15,
        atol=0.0,
    )


def test_plateau_derivatives_replay_with_independent_centered_difference():
    profile = KokunoPA10ReferencePlateauToX100()
    X = np.asarray([0.5, 2.0, 10.0, 70.0])
    eta = np.asarray([-0.5, -0.1, 0.25, 0.65])
    h = 2.0e-5
    plus = profile.values(X + h, eta)
    minus = profile.values(X - h, eta)
    finite_F = (plus["F_reference"] - minus["F_reference"]) / (2.0 * h)
    finite_U = (plus["U_reference"] - minus["U_reference"]) / (2.0 * h)
    finite_E = (plus["E_reference"] - minus["E_reference"]) / (2.0 * h)
    analytic = profile.radial_derivatives(X, eta)

    np.testing.assert_allclose(finite_F, analytic["F_reference_X"], atol=2.0e-11, rtol=0.0)
    np.testing.assert_allclose(finite_U, analytic["U_reference_X"], atol=2.0e-11, rtol=0.0)
    np.testing.assert_allclose(finite_E, analytic["E_reference_X"], atol=2.0e-9, rtol=2.0e-9)


def test_X2_handoff_has_zero_reference_radial_jet_and_X100_matches_it():
    profile = KokunoPA10ReferencePlateauToX100()
    eta = np.linspace(-0.8, 0.8, 9)
    x2 = profile.collar.handoff_at_collar_end(eta)
    x100 = profile.handoff_at_X100(eta)

    np.testing.assert_allclose(x2["F_reference_X"], 0.0, atol=1.0e-13, rtol=0.0)
    np.testing.assert_allclose(x2["U_reference_X"], 0.0, atol=1.0e-13, rtol=0.0)
    np.testing.assert_array_equal(x100["F_reference"], x2["F_reference"])
    np.testing.assert_array_equal(x100["U_reference"], x2["U_reference"])
    np.testing.assert_array_equal(x100["F_reference_X"], np.zeros_like(eta))
    np.testing.assert_array_equal(x100["U_reference_X"], np.zeros_like(eta))
    np.testing.assert_array_equal(x100["l_reference"], np.ones_like(eta))


def test_vectorized_shape_domain_and_source_scale_guards():
    profile = KokunoPA10ReferencePlateauToX100()
    X = np.asarray([[profile.X_2, 1.0], [20.0, 100.0]])
    eta = np.asarray([[-0.2], [0.3]])
    values = profile.values(X, eta)
    assert values["F_reference"].shape == (2, 2)
    assert values["U_reference"].shape == (2, 2)
    assert values["E_reference"].shape == (2, 2)
    assert profile.X_2 < profile.X_b_ref == 100.0 < profile.X_i == 110.0

    with pytest.raises(ValueError):
        profile.values(100.0001, 0.0)
    with pytest.raises(ValueError):
        profile.values(-1.0e-12, 0.0)
    with pytest.raises(ValueError):
        profile.values(1.0, profile.eta_interval[1] + 0.1)


def test_configuration_roundtrip_semantic_identity_and_truth_boundary(tmp_path):
    profile = KokunoPA10ReferencePlateauToX100()
    path = tmp_path / "reference_plateau_x100.json"
    payload = profile.save_configuration(path)
    loaded = KokunoPA10ReferencePlateauToX100.load_configuration(path)

    assert json.loads(path.read_text()) == payload
    assert payload["schema"] == SCHEMA
    assert loaded.configuration() == profile.configuration()
    assert loaded.semantic_sha256 == profile.semantic_sha256
    assert profile.truth_boundary["public_reference_plateau_to_X100_executable"] is True
    assert profile.truth_boundary["reference_stress_primitives_p1r_nsr_materialized"] is False
    assert profile.truth_boundary["kappa0_continuation_to_X100_materialized"] is False
    assert profile.truth_boundary["outer_global_leading_velocity_materialized"] is False
    assert profile.truth_boundary["pde_validated"] is False
    assert profile.truth_boundary["paper_exact"] is False
    assert profile.truth_boundary["openai_field_identified"] is False

    bad = dict(payload)
    bad["source_X_b_ref"] = 99.0
    with pytest.raises(ValueError):
        KokunoPA10ReferencePlateauToX100.from_configuration(bad)


def test_report_records_zero_plateau_drift_without_promoting_science():
    profile = KokunoPA10ReferencePlateauToX100()
    report = profile.report()
    handoff = report["X100_handoff_probe"]
    assert report["schema"] == SCHEMA
    assert report["geometry"]["X_b_ref"] == 100.0
    assert handoff["max_abs_F_plateau_value_drift"] == 0.0
    assert handoff["max_abs_U_plateau_value_drift"] == 0.0
    assert max(abs(v) for v in handoff["F_X"]) == 0.0
    assert max(abs(v) for v in handoff["U_X"]) == 0.0
    assert all(math.isclose(v, 1.0, rel_tol=0.0, abs_tol=0.0) for v in handoff["l_reference"])
    assert report["truth_boundary"]["actual_inner_to_outer_bridge_materialized"] is False
    assert report["truth_boundary"]["heldout_ns_residual_assessed"] is False
