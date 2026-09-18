import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rescaled_appendix_b_selected import (
    make_source_scale_aware_appendix_b_boundary,
)
from openai_ns_reconstruction.kokuno_rescaled_outer_scale_binding import (
    KokunoRescaledOuterScaleBinding,
)


@pytest.fixture(scope="module")
def binding():
    boundary = make_source_scale_aware_appendix_b_boundary(max_step=0.5)
    return KokunoRescaledOuterScaleBinding(boundary=boundary)


def test_selected_rescaled_C_is_not_the_old_finite_outer_C(binding):
    assert binding.selected_log_C > 1.0e10
    assert binding.template_log_C == pytest.approx(math.log(2.0), abs=2.0e-15)
    assert binding.selected_C_materializable_binary64 is False
    assert binding.finite_template_C_matches_selected_C is False
    with pytest.raises(ValueError, match="does not match the finite outer schedule C"):
        binding.require_finite_template_C_match()


def test_public_XR_formula_is_applied_with_the_same_selected_log_C(binding):
    schedule = binding.template_outer_schedule
    expected = math.log(110.0) + 10.0 * (
        binding.selected_log_C + schedule.log_P_star
    )
    assert binding.log_X_R == pytest.approx(expected, rel=2.0e-15)
    expected_shift = 10.0 * (binding.selected_log_C - math.log(schedule.C))
    assert binding.delta_log_X_R_from_finite_template == pytest.approx(
        expected_shift, rel=2.0e-15
    )
    assert binding.log_X_R - binding.template_log_X_R == pytest.approx(
        expected_shift, rel=2.0e-15, abs=2.0e-3
    )
    assert binding.log_X_R > binding.template_log_X_R + 1.0e10


def test_reserved_intervals_shift_with_XR_but_local_patch_amplitude_is_stable(binding):
    shift = binding.delta_log_X_R_from_finite_template
    selected = binding.reserved_log_intervals()
    template = binding.template_outer_schedule.reserved_log_intervals()
    assert selected.keys() == template.keys()
    for key in selected:
        for new_value, old_value in zip(selected[key], template[key], strict=True):
            assert new_value - old_value == pytest.approx(
                shift, rel=2.0e-15, abs=2.0e-3
            )
    assert binding.log_e_star == pytest.approx(
        binding.template_outer_schedule.log_e_star, rel=2.0e-13, abs=2.0e-13
    )


def test_pa15_log_scaling_and_separation_are_available_without_materializing_XR(binding):
    scales = binding.pa15_log_inverse_scales()
    assert scales == {
        "M": -binding.log_X_R,
        "I": -1.5 * binding.log_X_R,
        "J": -1.5 * binding.log_X_R,
        "S": -binding.log_X_R,
        "C_p": 0.0,
    }
    assert all(math.isfinite(value) for value in scales.values())
    assert binding.log_x_i == pytest.approx(
        math.log(110.0) - binding.log_X_R, rel=2.0e-15
    )
    assert binding.separation_geometry_feasible(1.0e6) is True
    with pytest.raises(ValueError, match="nonnegative"):
        binding.log_x_sep(-1.0)


def test_report_records_unproved_source_C_bounds(binding):
    report = binding.scale_report()
    assert report["finite_template_C_matches_selected_C"] is False
    assert report["selected_C_materializable_binary64"] is False
    assert report["source_complex_C_bound_verified"] is False
    assert report["source_all_PA11_C_bounds_verified"] is False
    assert np.isfinite(report["selected_log_X_R"])
    assert np.isfinite(report["selected_log_c_patch"])


def test_serialization_and_truth_boundary_fail_closed(binding, tmp_path):
    payload = binding.to_payload()
    replay = KokunoRescaledOuterScaleBinding.from_payload(payload)
    assert replay.sha256 == binding.sha256
    truth = payload["truth_boundary"]
    assert truth["shared_C_dependency_recorded"] is True
    assert truth["source_complex_C_bound_verified"] is False
    assert truth["actual_source_incoming_five_moment_discrepancy_bound"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    path = tmp_path / "rescaled-outer-scale-binding.json"
    binding.save_json(path)
    assert KokunoRescaledOuterScaleBinding.load_json(path).sha256 == binding.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_complex_C_bound_verified"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoRescaledOuterScaleBinding.from_payload(tampered)
