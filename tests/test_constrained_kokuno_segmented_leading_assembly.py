import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_segmented_leading_assembly import (
    KokunoSegmentedLeadingAssembly,
    UnreconstructedKokunoStageError,
)


def _physical_point_for_X(X: float, *, t: float = 0.5) -> tuple[float, float, float, float]:
    # z=0 gives eta=0 and q=1-t exactly in the native Kokuno coordinates.
    q = 1.0 - t
    return math.sqrt(2.0 * q * X), 0.0, 0.0, t


def test_reference_stage_delegates_exact_velocity_and_pressure():
    assembly = KokunoSegmentedLeadingAssembly()
    point = (0.2, -0.1, 0.05, 0.5)
    coordinates = assembly._coordinates.evaluate(*point)
    assert float(coordinates["X"]) < assembly.reference.X2
    assert assembly.stage_for_similarity(coordinates["X"]).item() == "reference_continuation"

    expected_velocity = assembly.reference.velocity(*point)
    actual_velocity = assembly.velocity(*point)
    np.testing.assert_array_equal(actual_velocity, expected_velocity)

    expected_pressure = assembly.reference.pressure(*point)
    actual_pressure = assembly.pressure(*point)
    np.testing.assert_array_equal(actual_pressure, expected_pressure)


def test_repaired_i2_stage_is_reached_and_delegates_velocity():
    assembly = KokunoSegmentedLeadingAssembly()
    X_star = math.exp(assembly.repaired_i2.outer_schedule.log_X_star)
    point = _physical_point_for_X(X_star)
    coordinates = assembly._coordinates.evaluate(*point)
    log_X = math.log(float(coordinates["X"]))
    low, high = assembly.i2_log_interval
    assert low < log_X < high
    assert assembly.stage_for_similarity(coordinates["X"]).item() == "audited_repaired_I2"

    expected = assembly.repaired_i2.velocity(*point)
    actual = assembly.velocity(*point)
    np.testing.assert_array_equal(actual, expected)
    assert np.all(np.isfinite(actual))
    assert np.linalg.norm(actual) > 0.0


def test_unreconstructed_gap_fails_closed_instead_of_bridging():
    assembly = KokunoSegmentedLeadingAssembly()
    # X=1 is immediately beyond the short reference continuation and vastly before I2.
    point = _physical_point_for_X(1.0)
    assert assembly.stage_for_similarity(np.asarray(1.0)).item() == "unreconstructed"
    with pytest.raises(UnreconstructedKokunoStageError, match="unreconstructed Kokuno radial stage"):
        assembly.velocity(*point)


def test_pressure_refuses_i2_until_global_matching_exists():
    assembly = KokunoSegmentedLeadingAssembly()
    X_star = math.exp(assembly.repaired_i2.outer_schedule.log_X_star)
    point = _physical_point_for_X(X_star)
    with pytest.raises(UnreconstructedKokunoStageError, match="global matched pressure"):
        assembly.pressure(*point)


def test_coverage_report_exposes_both_large_missing_gaps_and_heat_scale():
    assembly = KokunoSegmentedLeadingAssembly()
    report = assembly.coverage_report()
    low, high = assembly.i2_log_interval

    assert report["global_coverage_complete"] is False
    assert report["reference"]["executable_velocity"] is True
    assert report["repaired_I2"]["executable_velocity"] is True
    assert report["repaired_I2"]["independent_audit_pr"] == 339
    assert report["unreconstructed_log_X_gaps"][0] == [
        math.log(assembly.reference.X2),
        low,
    ]
    assert report["unreconstructed_log_X_gaps"][1] == [high, assembly.heat_log_X_K]
    assert assembly.heat_log_X_K > math.log(np.finfo(float).max)
    assert report["terminal_heat"]["matched_velocity_routed_here"] is False


def test_mixed_supported_vector_routes_each_executable_segment():
    assembly = KokunoSegmentedLeadingAssembly()
    core_point = _physical_point_for_X(0.2)
    X_star = math.exp(assembly.repaired_i2.outer_schedule.log_X_star)
    i2_point = _physical_point_for_X(X_star)

    x = np.asarray([core_point[0], i2_point[0]])
    y = np.asarray([0.0, 0.0])
    z = np.asarray([0.0, 0.0])
    t = np.asarray([0.5, 0.5])
    actual = assembly.velocity(x, y, z, t)

    expected_core = np.asarray(assembly.reference.velocity(*core_point), dtype=float)
    expected_i2 = np.asarray(assembly.repaired_i2.velocity(*i2_point), dtype=float)
    np.testing.assert_array_equal(actual[0], expected_core)
    np.testing.assert_array_equal(actual[1], expected_i2)


def test_payload_roundtrip_and_truth_boundary_fail_closed(tmp_path):
    assembly = KokunoSegmentedLeadingAssembly()
    path = tmp_path / "segmented.json"
    assembly.save_json(path)
    loaded = KokunoSegmentedLeadingAssembly.load_json(path)
    assert loaded.to_payload() == assembly.to_payload()
    assert loaded.sha256 == assembly.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["truth_boundary"]["global_leading_profile_reconstructed"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False
    payload["truth_boundary"]["global_leading_profile_reconstructed"] = True
    with pytest.raises(ValueError, match="payload hash or content mismatch"):
        KokunoSegmentedLeadingAssembly.from_payload(payload)
