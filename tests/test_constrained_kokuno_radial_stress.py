import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_mean_defect import PhaseMeanDefectContract
from openai_ns_reconstruction.kokuno_radial_stress import (
    CompactRadialStressInverse,
    RadialMeanDefectProfile,
    _compact_c4_annulus_window,
    sample_real_phase_mean_radial_profiles,
)


def _manufactured_profile(exponent: int, component: str, count: int = 65) -> RadialMeanDefectProfile:
    r0, r1 = 0.2, 0.8
    radii = np.linspace(r0, r1, count)
    window = _compact_c4_annulus_window(radii, r0, r1)
    raw = 0.7 + 0.4 * radii - 0.15 * radii**2
    return RadialMeanDefectProfile(
        radii=radii,
        raw_values=raw,
        values=raw * window,
        exponent=exponent,
        component=component,
        time=0.5,
        z=0.0,
        angular_count=16,
        annulus=(r0, r1),
    )


@pytest.mark.parametrize("exponent,component", [(2, "theta"), (1, "z")])
def test_moment_complement_closes_support_and_radial_identity(exponent, component):
    inverse = CompactRadialStressInverse(_manufactured_profile(exponent, component))
    metrics = inverse.metrics(query_count=401, fd_step=1.0e-5)

    assert abs(inverse.moment) > 1.0e-4
    assert metrics["normalized_bump_weighted_moment"] == 1.0
    assert metrics["moment_complement_weighted_moment_abs"] < 1.0e-12
    assert metrics["inner_edge_stress_abs"] == 0.0
    assert metrics["outer_edge_stress_abs"] == 0.0
    assert metrics["outside_stress_max_abs"] == 0.0
    assert metrics["raw_inverse_outer_tail_without_moment_subtraction"] > 1.0e-4
    assert metrics["radial_identity_fd_rms"] < 2.0e-8
    assert metrics["radial_identity_fd_max_abs"] < 1.0e-7


def test_equal_moment_subtraction_is_not_optional():
    profile = _manufactured_profile(2, "theta")
    inverse = CompactRadialStressInverse(profile)
    raw_outer_tail = abs(inverse.moment) / inverse.r_outer**profile.exponent
    assert raw_outer_tail > 1.0e-4
    assert inverse.stress(np.array([inverse.r_outer, inverse.r_outer + 0.1])).tolist() == [0.0, 0.0]


def test_profile_contract_fails_closed_on_wrong_source_pairing_and_edge_support():
    profile = _manufactured_profile(2, "theta")
    with pytest.raises(ValueError, match="source pairing"):
        RadialMeanDefectProfile(
            radii=profile.radii,
            raw_values=profile.raw_values,
            values=profile.values,
            exponent=1,
            component="theta",
            time=0.5,
            z=0.0,
            angular_count=16,
            annulus=profile.annulus,
        )

    bad_values = profile.values.copy()
    bad_values[0] = 0.1
    with pytest.raises(ValueError, match="vanish"):
        RadialMeanDefectProfile(
            radii=profile.radii,
            raw_values=profile.raw_values,
            values=bad_values,
            exponent=2,
            component="theta",
            time=0.5,
            z=0.0,
            angular_count=16,
            annulus=profile.annulus,
        )


def test_real_phase_mean_defect_feeds_radial_profiles_without_surrogate_tensor():
    correction = KokunoCompleteCurlCorrection()
    contract = PhaseMeanDefectContract(phase_count=8, step=0.01)

    def zero_base(points, time):
        del time
        return np.zeros_like(points, dtype=float)

    theta, axial, receipt = sample_real_phase_mean_radial_profiles(
        zero_base,
        correction,
        contract,
        time=0.5,
        z=0.0,
        r_inner=0.12,
        r_outer=0.72,
        radial_count=17,
        angular_count=8,
    )

    assert theta.source_kind == "real_phase_mean_defect_ring_average"
    assert axial.source_kind == "real_phase_mean_defect_ring_average"
    assert theta.exponent == 2 and theta.component == "theta"
    assert axial.exponent == 1 and axial.component == "z"
    assert max(np.max(np.abs(theta.raw_values)), np.max(np.abs(axial.raw_values))) > 1.0e-8
    assert receipt["full_mean_defect_increment_rms"] > 1.0e-8
    assert 0.0 <= receipt["theta_gate_capture_rms_ratio"] <= 1.0
    assert 0.0 <= receipt["z_gate_capture_rms_ratio"] <= 1.0


def test_annulus_must_stay_inside_agent2_transverse_box():
    correction = KokunoCompleteCurlCorrection()
    contract = PhaseMeanDefectContract(phase_count=8, step=0.01)

    def zero_base(points, time):
        del time
        return np.zeros_like(points, dtype=float)

    with pytest.raises(ValueError, match="transverse box"):
        sample_real_phase_mean_radial_profiles(
            zero_base,
            correction,
            contract,
            r_inner=0.12,
            r_outer=0.85,
            radial_count=17,
            angular_count=8,
        )
