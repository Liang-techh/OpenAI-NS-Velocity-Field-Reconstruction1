import copy

import pytest

import openai_ns_reconstruction.kokuno_radial_force_spacetime_envelope as mod


def _cell(time, z, *, rms=1.0e-4, stable=True, sha="candidate-sha"):
    return {
        "source": {
            "repository": "KokunoYumeto/yang-mills-interacting-workbench",
            "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
            "path": "navier-stokes/navier_stokes_workbench.tex",
            "formulas": "R33-R34 and R41",
        },
        "inputs": {
            "time": float(time),
            "z": float(z),
            "leading_candidate_sha256": sha,
            "surrogate_defect_used": False,
            "z_derivative_steps": list(mod.Z_DERIVATIVE_STEP_LADDER),
        },
        "radial_force_derivative_audit": {
            "fine_pair_relative_stability_tolerance": mod.FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
            "finest_pair_relative_rms_difference": 1.0e-3,
            "radial_force_derivative_stability_preflight_passed": bool(stable),
        },
        "full_tensor_force_completeness": {
            "finest_radial_force_rms": float(rms),
            "finest_radial_force_max_abs": 2.0 * float(rms),
            "radial_force_over_tangential_force_rms": 0.4,
            "radial_force_over_raw_ring_radial_defect_rms": 0.01,
            "radial_force_over_gated_ring_radial_defect_rms": 0.02,
            "radial_force_machine_nonzero": True,
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "surrogate_defect_used": False,
        },
    }


def _grid(times=(0.375, 0.5, 0.625), z_values=(0.06, 0.08, 0.10)):
    return [_cell(time, z, rms=1.0e-4 + 1.0e-5 * i) for i, (z, time) in enumerate(
        (z, time) for z in z_values for time in times
    )]


def test_default_spacetime_grid_reuses_agent3_frozen_screen():
    assert tuple(mod.HELD_OUT_CYCLE_TIMES) == (0.375, 0.5, 0.625)
    assert tuple(mod.SPATIAL_SCREEN_Z) == (0.06, 0.08, 0.10)
    assert mod.Z_DERIVATIVE_STEP_LADDER == (0.02, 0.01, 0.005)
    assert mod.FINE_PAIR_RELATIVE_STABILITY_TOLERANCE == 0.05


def test_summary_builds_conservative_envelope_over_complete_grid():
    reports = _grid()
    summary = mod.summarize_radial_force_spacetime(reports)
    assert summary["cell_count"] == 9
    assert summary["all_cells_derivative_stable"] is True
    assert summary["all_cells_radial_force_machine_nonzero"] is True
    assert summary["radial_force_rms_min"] == pytest.approx(1.0e-4)
    assert summary["radial_force_rms_max"] == pytest.approx(1.8e-4)
    assert summary["radial_force_max_abs_across_cells"] == pytest.approx(3.6e-4)
    assert summary["strongest_radial_force_cell"]["time"] == pytest.approx(0.625)
    assert summary["strongest_radial_force_cell"]["z"] == pytest.approx(0.10)


def test_summary_retains_scientific_derivative_rejection_without_crashing():
    reports = _grid()
    reports[4]["radial_force_derivative_audit"][
        "radial_force_derivative_stability_preflight_passed"
    ] = False
    reports[4]["radial_force_derivative_audit"][
        "finest_pair_relative_rms_difference"
    ] = 0.08
    summary = mod.summarize_radial_force_spacetime(reports)
    assert summary["all_cells_derivative_stable"] is False
    assert summary["finest_pair_relative_rms_difference_max"] == pytest.approx(0.08)


def test_summary_fails_closed_on_missing_or_duplicate_cells():
    reports = _grid()
    with pytest.raises(ValueError, match="missing frozen spacetime cells"):
        mod.summarize_radial_force_spacetime(reports[:-1])
    duplicate = reports + [copy.deepcopy(reports[0])]
    with pytest.raises(ValueError, match="unique"):
        mod.summarize_radial_force_spacetime(duplicate)


def test_summary_rejects_surrogate_changed_ladder_or_provenance():
    reports = _grid()
    bad = copy.deepcopy(reports)
    bad[0]["inputs"]["surrogate_defect_used"] = True
    with pytest.raises(ValueError, match="real candidate defect"):
        mod.summarize_radial_force_spacetime(bad)

    bad = copy.deepcopy(reports)
    bad[0]["inputs"]["z_derivative_steps"] = [0.03, 0.01, 0.005]
    with pytest.raises(ValueError, match="z-derivative ladder"):
        mod.summarize_radial_force_spacetime(bad)

    bad = copy.deepcopy(reports)
    bad[-1]["source"]["commit"] = "different"
    with pytest.raises(ValueError, match="provenance"):
        mod.summarize_radial_force_spacetime(bad)


def test_summary_rejects_candidate_identity_drift():
    reports = _grid()
    reports[-1]["inputs"]["leading_candidate_sha256"] = "other-sha"
    with pytest.raises(ValueError, match="leading candidate changed"):
        mod.summarize_radial_force_spacetime(reports)


def test_generate_report_keeps_cycle_fail_closed(monkeypatch, tmp_path):
    times = (0.375, 0.625)
    z_values = (0.06, 0.10)

    def fake_cell(*, output, time, z, z_steps):
        assert tuple(z_steps) == mod.Z_DERIVATIVE_STEP_LADDER
        return _cell(time, z, rms=1.0e-4 + 1.0e-5 * (time + z))

    monkeypatch.setattr(mod, "generate_radial_force_cell_report", fake_cell)
    report = mod.generate_actual_core_report(
        output_dir=tmp_path / "report",
        times=times,
        z_values=z_values,
    )
    assert report["radial_force_spacetime_envelope"]["cell_count"] == 4
    assert report["routing"]["finite_correction_cycle_rerun_allowed"] is False
    assert report["routing"]["second_public_covariance_column_available"] is False
    assert report["truth_boundary"]["surrogate_defect_used"] is False
    assert report["truth_boundary"]["residual_reduction_claimed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert (tmp_path / "report" / "summary.json").is_file()
