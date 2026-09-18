import copy
import hashlib
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pressure_datum_binding import (
    KokunoPressureDatumBinding,
)


@pytest.fixture(scope="module")
def binding():
    return KokunoPressureDatumBinding()


def test_public_axis_pressure_bound_is_bound_exactly_without_hidden_recovery(binding):
    expected_required = math.sqrt(2.5) * math.exp(binding.outer_schedule.log_P_star)
    assert binding.required_pressure_scale == pytest.approx(expected_required, rel=2e-15)
    assert binding.selected_pressure_scale >= binding.required_pressure_scale
    assert binding.pressure_ratio_to_lower_bound >= 1.0

    eta = np.linspace(-1.0, 1.0, 17)
    selected = binding.selected_axis_pressure(eta)
    upper = binding.source_axis_pressure_upper_bound(eta)
    assert np.all(selected <= upper)
    assert np.all(binding.pressure_bound_satisfied(eta))

    report = binding.datum_report()
    assert report["all_sample_bounds_satisfied"] is True
    assert report["algebraic_bound_satisfied_by_construction"] is True
    assert report["required_pressure_scale"] > 1e13
    assert report["selected_to_previous_baseline_ratio"] > 1e13
    assert report["hidden_source_value_recovered"] is False


def test_selected_B13_reference_profile_and_velocity_are_executable(binding):
    eta = np.asarray([-0.5, 0.0, 0.5])
    values = binding.profile_values(np.zeros_like(eta), eta)
    assert set(values) == {"F", "U", "Pi"}
    assert np.all(np.isfinite(values["F"]))
    assert np.all(np.isfinite(values["U"]))
    assert np.all(np.isfinite(values["Pi"]))
    assert np.allclose(values["U"], 4.0 * eta + binding.j0, rtol=0.0, atol=1e-15)
    assert np.allclose(values["Pi"], binding.selected_axis_pressure(eta), rtol=2e-15)

    velocity = binding.velocity(
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.0]),
        np.asarray([0.0, 0.1]),
        np.asarray([0.5, 0.5]),
    )
    assert velocity.shape == (2, 3)
    assert np.all(np.isfinite(velocity))
    assert np.max(np.linalg.norm(velocity, axis=-1)) > 0.0


def test_binary64_core_probe_is_an_explicit_receipt_not_a_source_claim(binding):
    # Degree zero is guaranteed to exercise the selected pressure datum without
    # entering the pressure-driven nonlinear recurrence.
    finite = binding.series_probe(maxdegree=0, eta_nodes=33)
    assert finite["status"] == "finite"
    assert finite["maxdegree"] == 0
    assert finite["eta_nodes"] == 33
    assert finite["candidate_sha256"]

    # The production-degree probe may succeed or expose a binary64 dynamic-range
    # barrier.  Either result is a valid executable receipt; neither promotes
    # the scientific truth boundary.
    production = binding.series_probe(maxdegree=14, eta_nodes=33)
    assert production["status"] in {"finite", "nonfinite_or_failed"}
    json.dumps(production, sort_keys=True, allow_nan=False)

    truth = binding.report()["truth_boundary"]
    assert truth["selected_source_compatible_pressure_datum_executable"] is True
    assert truth["source_hidden_pressure_scale_recovered"] is False
    assert truth["actual_source_appendix_B_trajectory_reconstructed"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path, binding):
    path = binding.save_json(tmp_path / "pressure_datum.json")
    replay = KokunoPressureDatumBinding.load_json(path)
    assert replay.to_payload() == binding.to_payload()
    assert replay.sha256 == binding.sha256

    payload = copy.deepcopy(binding.to_payload())
    payload["truth_boundary"]["global_pressure_matched"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoPressureDatumBinding.from_payload(payload)


def test_margin_and_domain_guards_fail_closed():
    with pytest.raises(ValueError, match="pressure_margin"):
        KokunoPressureDatumBinding(pressure_margin=0.999)
    obj = KokunoPressureDatumBinding(pressure_margin=1.1)
    assert obj.pressure_square_ratio_to_source_bound > 1.0
    with pytest.raises(ValueError, match="eta"):
        obj.selected_axis_pressure(1.01)
