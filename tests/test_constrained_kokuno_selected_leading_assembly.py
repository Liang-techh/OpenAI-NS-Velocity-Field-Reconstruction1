import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_selected_leading_assembly import (
    KokunoSelectedLeadingAssembly,
    UnreconstructedSelectedKokunoStageError,
)


@pytest.fixture(scope="module")
def assembly():
    return KokunoSelectedLeadingAssembly()


def _physical_point(log_X, eta, h):
    eta = float(eta)
    q = 1.0 / (1.0 - eta * eta)
    D = 0.5 - float(h)
    z = q**D * eta
    radius = math.sqrt(2.0 * q) * math.exp(0.5 * float(log_X))
    return radius, 0.0, z, 0.0


def test_selected_leading_stage_coverage_is_fail_closed(assembly):
    intervals = assembly._intervals()
    mod_left, mod_right = intervals["modulation"]
    i1_low, i1_high = intervals["I1"]
    i2_low, i2_high = intervals["I2"]
    i3_low, i3_high = intervals["I3"]

    probes_log = np.asarray(
        [
            0.5 * (assembly.outer_base.log_X_R + mod_left),
            0.5 * (mod_left + mod_right),
            0.5 * (mod_right + i1_low),
            0.5 * (i1_low + i1_high),
            0.5 * (i2_low + i2_high),
            0.5 * (i2_high + i3_low),
            0.5 * (i3_low + i3_high),
        ]
    )
    stages = assembly.stage_for_similarity(np.exp(probes_log)).tolist()
    assert stages == [
        "rf40_outer_base",
        "selected_autonomous_modulation",
        "post_modulation_prefix_memory",
        "selected_eta_smooth_I1_repair",
        "audited_repaired_I2",
        "post_I2_base_before_I3",
        "unreconstructed",
    ]

    gap_log = 0.5 * (math.log(assembly.reference.X2) + assembly.outer_base.log_X_R)
    assert assembly.stage_for_similarity(math.exp(gap_log)).item() == "unreconstructed"
    report = assembly.coverage_report()
    assert report["global_coverage_complete"] is False
    assert report["global_pressure_matched"] is False
    assert report["selected_outer"]["modulation_prefix_M_memory_carried"] is True


def test_reference_axis_and_rf40_base_delegate_exactly(assembly):
    reference_value = np.asarray(assembly.reference.velocity(0.0, 0.0, 0.0, 0.5))
    routed_reference = np.asarray(assembly.velocity(0.0, 0.0, 0.0, 0.5))
    assert np.all(np.isfinite(routed_reference))
    assert np.array_equal(routed_reference, reference_value)

    mod_left, _ = assembly._intervals()["modulation"]
    log_X = 0.5 * (assembly.outer_base.log_X_R + mod_left)
    point = _physical_point(log_X, 0.2, assembly.h)
    expected = np.asarray(assembly.outer_base.velocity(*point))
    routed = np.asarray(assembly.velocity(*point))
    assert np.array_equal(routed, expected)


def test_modulation_prefix_memory_survives_after_local_support(assembly):
    mod_left, mod_right = assembly._intervals()["modulation"]
    i1_low, _ = assembly._intervals()["I1"]
    eta = 0.2

    inside = _physical_point(0.5 * (mod_left + mod_right), eta, assembly.h)
    inside_value = np.asarray(assembly.velocity(*inside))
    inside_base = np.asarray(assembly.outer_base.velocity(*inside))
    assert np.all(np.isfinite(inside_value))
    assert not np.array_equal(inside_value, inside_base)

    # Outside the compact A/B support, E and U have returned to the RF40 base,
    # but incompressibility still sees the nonzero prefix M until I1 repairs it.
    after_log = 0.5 * (mod_right + i1_low)
    after = _physical_point(after_log, eta, assembly.h)
    after_value = np.asarray(assembly.velocity(*after))
    after_base = np.asarray(assembly.outer_base.velocity(*after))
    assert np.all(np.isfinite(after_value))
    assert not np.array_equal(after_value, after_base)

    memory = assembly.moment_memory_report(eta)
    assert memory["normalized_modulation_M"] != 0.0


def test_I1_residual_memory_matches_eta_smooth_family_closure(assembly):
    eta = 0.2
    memory = assembly.moment_memory_report(eta)
    closure = assembly.i1_family.closure_report(eta)
    residual_M = float(memory["normalized_post_I1_M_residual"])
    closure_M = float(closure["residual_normalized"][0])
    assert residual_M == pytest.approx(closure_M, rel=1.0e-10, abs=1.0e-15)

    i1_low, i1_high = assembly._intervals()["I1"]
    point = _physical_point(0.5 * (i1_low + i1_high), eta, assembly.h)
    value = np.asarray(assembly.velocity(*point))
    assert np.all(np.isfinite(value))


def test_I2_places_audited_heat_delta_on_full_rf40_base(assembly):
    i2_low, i2_high = assembly._intervals()["I2"]
    log_X = 0.5 * (i2_low + i2_high)
    point = _physical_point(log_X, 0.0, assembly.h)

    routed = np.asarray(assembly.velocity(*point))
    full_base = np.asarray(assembly.outer_base.velocity(*point))
    local_patch = np.asarray(assembly.repaired_i2.velocity(*point))
    assert np.all(np.isfinite(routed))

    # At eta=0 and y=0, the selected M-memory radial term vanishes.  The local
    # I2 helper is pure swirl and has zero x component, whereas RF40 retains its
    # incompressibility-derived historical v0.  The assembly must keep RF40 v0.
    assert full_base[0] != 0.0
    assert local_patch[0] == 0.0
    assert routed[0] == pytest.approx(full_base[0], rel=2.0e-13, abs=0.0)


def test_missing_join_I3_and_outer_pressure_fail_closed(assembly):
    join_log = 0.5 * (math.log(assembly.reference.X2) + assembly.outer_base.log_X_R)
    join_point = _physical_point(join_log, 0.0, assembly.h)
    with pytest.raises(UnreconstructedSelectedKokunoStageError):
        assembly.velocity(*join_point)

    i3_low, i3_high = assembly._intervals()["I3"]
    i3_point = _physical_point(0.5 * (i3_low + i3_high), 0.0, assembly.h)
    with pytest.raises(UnreconstructedSelectedKokunoStageError):
        assembly.velocity(*i3_point)

    base_point = _physical_point(assembly.outer_base.log_X_R + 0.5, 0.0, assembly.h)
    with pytest.raises(UnreconstructedSelectedKokunoStageError):
        assembly.pressure(*base_point)


def test_selected_leading_payload_round_trip_and_truth_tamper(assembly):
    payload = assembly.to_payload()
    restored = KokunoSelectedLeadingAssembly.from_payload(copy.deepcopy(payload))
    assert restored.sha256 == assembly.sha256
    assert restored.coverage_report() == assembly.coverage_report()

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["global_leading_profile_reconstructed"] = True
    with pytest.raises(ValueError):
        KokunoSelectedLeadingAssembly.from_payload(tampered)
