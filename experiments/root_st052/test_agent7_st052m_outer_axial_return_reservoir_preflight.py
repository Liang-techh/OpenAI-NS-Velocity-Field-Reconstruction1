from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "agent7_st052m_outer_axial_return_reservoir_preflight.py"
SPEC = importlib.util.spec_from_file_location("agent7_outer_axial_return_reservoir", MODULE_PATH)
assert SPEC and SPEC.loader
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_frozen_lineage_and_truth_boundary() -> None:
    report = audit.evaluate()
    assert report["task_id"] == "CR003-ST052M-OUTER-AXIAL-RETURN-RESERVOIR-PREFLIGHT-108"
    assert report["issue"] == 766
    assert report["source_parent"] == {
        "pr": 759,
        "head": "40c7c1cd0ed611788b6766a0e5b3a5782a0418d7",
    }
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["coefficient_fit_performed"] is False
    assert report["fresh_714_path_data_used"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["pressure_or_force_changed"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_visible_tip_is_inward_and_outer_reservoir_is_outward() -> None:
    report = audit.evaluate()
    directional = report["directional_spot_checks"]
    assert all(v < 0.0 for v in directional["visible_inward_radial_components"].values())
    assert all(v > 0.0 for v in directional["reservoir_outward_radial_components"].values())
    assert directional["positive_negative_z_max_radial_mismatch"] == 0.0
    assert directional["p9_radial_at_abs_z_1p75"] > 0.0
    assert directional["reservoir_radial_at_abs_z_1p75"] < 0.0
    assert directional["p9_outward_reservoir_inward_at_1p75"] is True


def test_support_and_endpoint_semantics_are_frozen() -> None:
    report = audit.evaluate()
    assert report["exact_endpoint_join_checks"]["potential"] == [0.0, 1.0, 0.0, -0.0, -1.0, -0.0]
    assert max(abs(v) for v in report["exact_endpoint_join_checks"]["derivative"]) == 0.0
    support = report["support_checks"]
    assert support["shoulder_max_abs"] == 0.0
    assert support["axial_outside_max_abs"] == 0.0
    assert support["radial_outside_max_abs"] == 0.0


def test_signed_integral_debt_is_relocated_not_removed() -> None:
    report = audit.evaluate()
    gate = report["signed_integral_refinement"]["gate"]
    assert gate["passes"] is True
    assert gate["max_abs_finest_full_integral"] <= audit.INTEGRAL_TOL
    assert gate["max_abs_finest_visible_plus_reservoir"] <= audit.INTEGRAL_TOL
    assert gate["max_abs_finest_endpoint_error"] <= audit.INTEGRAL_TOL
    assert gate["max_medium_to_fine_change"] <= audit.REFINEMENT_TOL
    fine = report["signed_integral_refinement"]["levels"][str(audit.GRID_SIZES[-1])]
    for row in fine.values():
        assert row["visible_signed_integral"] < 0.0
        assert row["reservoir_signed_integral"] > 0.0
    scope = report["structural_scope"]
    assert scope["return_flow_removed"] is False
    assert scope["return_flow_relocated_outside_visible_1_to_1p8_lobe"] is True


def test_reservoir_shape_adds_axial_response_capacity_without_dimension_growth() -> None:
    report = audit.evaluate()
    geometry = report["response_geometry_vs_p9"]
    assert geometry["passes"] is True
    assert geometry["rank"] == 2
    assert abs(geometry["column_cosine"]) < audit.COSINE_MAX_ABS
    assert geometry["condition_number"] <= audit.CONDITION_MAX
    decision = report["decision"]
    assert decision["outer_axial_return_reservoir_is_structural_escape_from_1p8_endpoint_obstruction"] is True
    assert decision["one_later_frozen_nonlinear_child_replay_justified"] is True
    assert decision["add_second_basis_dimension_now"] is False
    assert decision["closer_visualization_delivery_this_increment"] is False
