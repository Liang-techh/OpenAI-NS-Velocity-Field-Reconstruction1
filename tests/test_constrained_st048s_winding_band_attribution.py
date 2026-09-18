import pytest

from openai_ns_reconstruction import constrained_st048s_winding_band_attribution as mod


def _paths(*, st006: bool, turns_by_radius=(0.04, 0.025, 0.005), radius_changes=(-0.08, -0.10, -0.06)):
    rows = []
    index = 0
    for radius, turns, dr in zip(mod.SEED_RADII, turns_by_radius, radius_changes):
        for angle_index in range(mod.SEED_ANGLES):
            for z in mod.SEED_Z:
                row = {
                    "path_index": index,
                    "seed": {
                        "radius": radius,
                        "angle_index": angle_index,
                        "angle_radians": 2.0 * 3.141592653589793 * angle_index / mod.SEED_ANGLES,
                        "z": z,
                    },
                    "radius_change": dr,
                    "abs_z_change": 0.02 if radius < 1.2 else -0.01,
                }
                if st006:
                    row["signed_turns"] = turns
                else:
                    row["absolute_turns"] = turns
                rows.append(row)
                index += 1
    return rows


def _st006():
    return {
        "task_id": mod.ST006_TASK_ID,
        "schema": mod.ST006_SCHEMA,
        "report_sha256": mod.ST006_REPORT_SHA256,
        "candidate": "ST006",
        "candidate_sha256": mod.ST006_CANDIDATE_SHA256,
        "registered_contract": {
            "seed_radii": list(mod.SEED_RADII),
            "seed_z": list(mod.SEED_Z),
            "seed_angles": mod.SEED_ANGLES,
            "seed_count": 48,
            "time_interval": list(mod.TIME_INTERVAL),
            "output_samples": mod.OUTPUT_SAMPLES,
            "solver": {
                "method": mod.SOLVER_METHOD,
                "rtol": mod.SOLVER_RTOL,
                "atol": mod.SOLVER_ATOL,
                "max_step": mod.SOLVER_MAX_STEP,
            },
        },
        "interpretation_boundary": {
            "all_sampled_points_inside_registered_box": True,
            "pde_validated": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "openai_field_identified": False,
            "paper_exact": False,
            "blowup_proved": False,
            "visual_acceptance_threshold_defined": False,
            "image_fit_used": False,
            "hidden_time_alignment_used": False,
            "camera_registration_used": False,
            "openai_seed_locations_used": False,
        },
        "per_path": _paths(st006=True),
    }


def _swirl():
    return {
        "task_id": mod.SWIRL_TASK_ID,
        "schema": mod.SWIRL_SCHEMA,
        "report_sha256": mod.SWIRL_REPORT_SHA256,
        "candidate": {
            "candidate_id": "ST048-S",
            "original_raw_candidate_sha256": mod.ST048S_RAW_SHA256,
            "pde_validated": False,
            "source_correspondence_verified": False,
        },
        "frozen_contract": {
            "seed_radii": list(mod.SEED_RADII),
            "seed_z": list(mod.SEED_Z),
            "seed_angles": mod.SEED_ANGLES,
            "path_count": 48,
            "time_interval": list(mod.TIME_INTERVAL),
            "output_samples": mod.OUTPUT_SAMPLES,
            "solver_method": mod.SOLVER_METHOD,
            "solver_rtol": mod.SOLVER_RTOL,
            "solver_atol": mod.SOLVER_ATOL,
            "solver_max_step": mod.SOLVER_MAX_STEP,
        },
        "truth_boundary": {
            "canonical_velocity_changed": False,
            "production_candidate_selected": False,
            "production_kappa_selected": False,
            "held_out_pde_residual_evaluated": False,
            "visual_acceptance_threshold_defined": False,
            "hidden_openai_time_camera_seed_or_velocity_used": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "comparison_is_descriptive_not_acceptance": True,
        },
        "measurements": {
            "temporal_piola_kappa_0": {
                "per_path": _paths(st006=False, turns_by_radius=(0.029, 0.021, 0.006)),
            },
            "temporal_piola_kappa_005": {
                "per_path": _paths(st006=False, turns_by_radius=(0.030, 0.022, 0.0062), radius_changes=(-0.078, -0.098, -0.058)),
            },
        },
    }


def test_stratification_routes_largest_deficit_to_inner_radius():
    report = mod.build_report(_st006(), _swirl())
    mod.audit_truth_boundary(report)
    assert report["routing"]["largest_candidate_side_winding_deficit_radius"] == pytest.approx(0.6)
    assert report["routing"]["outer_r12_is_winding_surplus_vs_st006"] is True
    assert report["truth_boundary"]["production_radial_profile_selected"] is False
    assert report["analysis_contract"]["new_integration"] is False


def test_rejects_contract_drift():
    bad = _swirl()
    bad["frozen_contract"]["seed_radii"] = [0.55, 0.9, 1.2]
    with pytest.raises(ValueError, match="seed-radii"):
        mod.build_report(_st006(), bad)


def test_rejects_truth_promotion():
    bad = _swirl()
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        mod.build_report(_st006(), bad)


def test_rejects_report_identity_drift():
    bad = _st006()
    bad["report_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="report SHA"):
        mod.build_report(bad, _swirl())


def test_rejects_missing_path_seed():
    bad = _st006()
    bad["per_path"] = bad["per_path"][:-1]
    with pytest.raises(ValueError, match="exactly 48 paths"):
        mod.build_report(bad, _swirl())
