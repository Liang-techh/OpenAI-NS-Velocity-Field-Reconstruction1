import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_reference_p1_x100 import (
    SCHEMA,
    ETA_FD_STEP,
    KokunoPA10ReferenceP1ToX100,
)


def test_p1_X0_is_exact_stress_free_initial_datum():
    primitive = KokunoPA10ReferenceP1ToX100()
    eta = np.asarray([-0.6, 0.0, 0.6])
    X = np.full_like(eta, primitive.X_0)
    got = primitive.values(X, eta)["p1_reference"]
    expected = primitive.initial_p1(eta)
    np.testing.assert_allclose(got, expected, rtol=2.0e-13, atol=2.0e-13)
    assert np.all(np.isfinite(got))
    assert np.linalg.norm(got) > 1.0e-12


def test_reference_M_reuses_natural_primitive_at_X0():
    primitive = KokunoPA10ReferenceP1ToX100()
    eta = np.asarray([-0.5, 0.1, 0.55])
    X = np.full_like(eta, primitive.X_0)
    got = primitive.M_reference(X, eta)
    natural = primitive.reference.collar.source_normalized.values(X, eta)["M_0"]
    np.testing.assert_allclose(got, natural, rtol=0.0, atol=2.0e-14)


def test_p1_radial_ode_replays_with_independent_centered_difference():
    primitive = KokunoPA10ReferenceP1ToX100()
    X = 4.0
    eta = 0.25
    step = 2.0e-3
    plus = primitive.values(X + step, eta)["p1_reference"]
    minus = primitive.values(X - step, eta)["p1_reference"]
    fd = (plus - minus) / (2.0 * step)
    ode = primitive.radial_derivative(X, eta)
    scale = max(1.0, abs(float(ode)))
    assert abs(float(fd - ode)) / scale < 2.0e-3


def test_source_terms_and_X100_handoff_are_finite_nontrivial():
    primitive = KokunoPA10ReferenceP1ToX100()
    eta = np.asarray([-0.5, 0.0, 0.5])
    X = np.asarray([primitive.X_0, 10.0, 100.0])
    terms = primitive.source_terms(X, eta)
    for key in (
        "M_reference",
        "M_reference_eta",
        "log_E_reference_eta",
        "W_reference",
        "H_c_reference",
        "H_reference",
        "ell_reference",
        "L",
        "S_q_reference",
    ):
        assert np.all(np.isfinite(terms[key]))
    handoff = primitive.handoff_at_X100(eta)
    assert np.all(np.isfinite(handoff["p1_reference"]))
    assert np.all(np.isfinite(handoff["p1_reference_X"]))
    assert np.linalg.norm(handoff["p1_reference"]) > 1.0e-12


def test_vectorization_domain_serialization_and_truth_boundary(tmp_path):
    primitive = KokunoPA10ReferenceP1ToX100()
    X = np.asarray([[primitive.X_0, 2.0], [10.0, 100.0]])
    eta = np.asarray([[-0.2], [0.3]])
    got = primitive.values(X, eta)["p1_reference"]
    assert got.shape == (2, 2)

    lo, hi = primitive.eta_interval
    with pytest.raises(ValueError):
        primitive.values(primitive.X_0 - 1.0e-9, 0.0)
    with pytest.raises(ValueError):
        primitive.values(100.001, 0.0)
    with pytest.raises(ValueError):
        primitive.values(1.0, lo + ETA_FD_STEP)
    with pytest.raises(ValueError):
        primitive.values(1.0, hi - ETA_FD_STEP)

    path = tmp_path / "reference_p1_x100.json"
    payload = primitive.save_configuration(path)
    loaded = KokunoPA10ReferenceP1ToX100.load_configuration(path)
    assert json.loads(path.read_text()) == payload
    assert payload["schema"] == SCHEMA
    assert loaded.configuration() == primitive.configuration()
    assert loaded.semantic_sha256 == primitive.semantic_sha256

    truth = primitive.truth_boundary
    assert truth["public_reference_p1_ode_executable_to_X100"] is True
    assert truth["reference_eta_jets_analytic"] is False
    assert truth["reference_nsr_materialized"] is False
    assert truth["kappa0_continuation_to_X100_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    bad = dict(payload)
    bad["eta_fd_step"] = 1.0e-4
    with pytest.raises(ValueError):
        KokunoPA10ReferenceP1ToX100.from_configuration(bad)
