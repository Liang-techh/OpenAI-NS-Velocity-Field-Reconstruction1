import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_ideal_exterior_qs_schedule import (
    KokunoPublicIdealExteriorQsSchedule,
)


def _centered_derivative(fun, x: float, eps: float = 2.0e-6) -> float:
    return float((np.asarray(fun(x + eps)) - np.asarray(fun(x - eps))) / (2.0 * eps))


def test_source_ideal_entry_and_stage_endpoints_are_deterministic():
    q = KokunoPublicIdealExteriorQsSchedule()
    assert q.lambda_value == pytest.approx(0.05, rel=0.0, abs=0.0)
    assert q.h_value == pytest.approx(0.005, rel=0.0, abs=0.0)
    assert q.q_exterior_entry_source_ideal == pytest.approx(
        (0.05 - 0.005) / 0.95, rel=0.0, abs=2e-16
    )
    assert q.q_release1_end_source_ideal == pytest.approx(
        0.5435933032548154, rel=2e-13, abs=2e-15
    )
    assert q.hold_length == pytest.approx(
        4.0 * math.log(200.0), rel=0.0, abs=2e-15
    )
    assert q.q_hold_end_source_ideal == pytest.approx(
        21.630896422116002, rel=2e-13, abs=2e-14
    )
    assert q.q_in_source_ideal == pytest.approx(
        13.455991081964195, rel=2e-13, abs=2e-14
    )


def test_release1_and_release2_satisfy_public_ode_away_from_flat_endpoints():
    q = KokunoPublicIdealExteriorQsSchedule()
    for s in (0.12, 0.37, 0.63, 0.88):
        fd1 = _centered_derivative(q.q_release1, s)
        rhs1 = float(np.asarray(q.release1_ode_rhs(s)))
        assert fd1 == pytest.approx(rhs1, rel=2e-8, abs=2e-9)

        fd2 = _centered_derivative(q.q_release2, s)
        rhs2 = float(np.asarray(q.release2_ode_rhs(s)))
        assert fd2 == pytest.approx(rhs2, rel=2e-8, abs=2e-9)


def test_lminus1_hold_is_exactly_linear_source_transport():
    q = KokunoPublicIdealExteriorQsSchedule()
    pts = np.array([0.0, 0.1, 0.5, 0.9]) * q.hold_length
    vals = q.q_hold(pts)
    expected = q.q_release1_end_source_ideal + (1.0 - q.h_value) * pts
    assert np.array_equal(vals, expected)


def test_ideal_matching_bridge_hits_exact_parent_qp_target():
    q = KokunoPublicIdealExteriorQsSchedule()
    assert q.q_in_source_ideal > q.q_p > 0.0
    expected = math.log(q.q_in_source_ideal / q.q_p) / (1.0 - q.h_value)
    assert q.ideal_matching_length == pytest.approx(expected, rel=0.0, abs=0.0)
    assert q.ideal_matching_length == pytest.approx(10.1116603621, rel=2e-11)
    end = float(np.asarray(q.q_matching_bridge_source_ideal(q.ideal_matching_length)))
    assert end == pytest.approx(q.q_p, rel=3e-15, abs=1e-18)


def test_matching_bridge_is_exponential_and_vectorized():
    q = KokunoPublicIdealExteriorQsSchedule()
    y = np.array([0.0, 0.25, 0.5, 1.0]) * q.ideal_matching_length
    values = q.q_matching_bridge_source_ideal(y)
    expected = q.q_in_source_ideal * np.exp(-(1.0 - q.h_value) * y)
    assert np.array_equal(values, expected)
    with pytest.raises(ValueError):
        q.q_matching_bridge_source_ideal(q.ideal_matching_length + 0.01)


def test_report_and_truth_boundary_do_not_launder_ideal_transport_into_current_candidate():
    q = KokunoPublicIdealExteriorQsSchedule()
    report = q.schedule_report()
    truth = report["truth_boundary"]
    assert report["source_ideal_assumptions"] == [
        "I=XH/(1-lambda)",
        "U=0",
        "M=0",
        "J=0",
    ]
    assert report["current_candidate_q_s_claimed"] is False
    assert report["current_cartesian_bridge_composed"] is False
    assert report["matching_end_relative_error_to_q_p"] < 5e-15
    assert truth["public_source_ideal_exterior_qs_schedule_materialized"] is True
    assert truth["public_source_ideal_qs_release2_endpoint_materialized"] is True
    assert truth["public_source_ideal_l_minus_h_matching_length_materialized"] is True
    assert truth["current_q_s_release2_endpoint_materialized"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["current_cartesian_terminal_multiplier_composed"] is False
    assert truth["source_exterior_heat_replacement_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_configuration_roundtrip_and_mutation_fail_closed(tmp_path):
    q = KokunoPublicIdealExteriorQsSchedule()
    path = tmp_path / "qs.json"
    q.save_configuration(path)
    loaded = KokunoPublicIdealExteriorQsSchedule.load_configuration(path)
    assert loaded.semantic_sha256 == q.semantic_sha256
    assert loaded.configuration() == q.configuration()

    payload = json.loads(path.read_text())
    payload["truth_boundary"]["current_q_s_release2_endpoint_materialized"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPublicIdealExteriorQsSchedule.from_configuration(payload)


def test_no_scientific_tuning_or_cartesian_velocity_surface_is_exposed():
    q = KokunoPublicIdealExteriorQsSchedule()
    public = {
        name
        for name, value in inspect.getmembers(type(q))
        if not name.startswith("_") and callable(value)
    }
    assert "velocity" not in public
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "viscosity",
        "optimizer",
        "gain",
        "threshold",
        "q_in",
        "q_p_input",
        "matching_length_input",
    }
    assert not (public & forbidden)
