from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "agent7_st052m_fixed_radius_return_balance_audit.py"
SPEC = importlib.util.spec_from_file_location("agent7_fixed_radius_return_balance", MODULE_PATH)
assert SPEC and SPEC.loader
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_frozen_lineage_and_truth_boundary() -> None:
    report = audit.evaluate()
    assert report["task_id"] == "CR003-ST052M-FIXED-RADIUS-RETURN-BALANCE-AUDIT-106"
    assert report["issue"] == 758
    assert report["source_parent"] == {
        "pr": 750,
        "head": "376223a3bb2dae0d53756ea48da26937a95b6e58",
    }
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["second_poloidal_basis_added_to_candidate"] is False
    assert report["fresh_714_path_data_used"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_endpoint_identity_is_exact_at_every_frozen_radius() -> None:
    report = audit.evaluate()
    checks = report["exact_fixed_radius_identity"]["endpoint_checks"]
    for radius_record in checks.values():
        for channel_record in radius_record.values():
            assert channel_record["F_inner"] == 0.0
            assert channel_record["F_outer"] == 0.0
            assert channel_record["exact_integral_from_endpoints"] == 0.0
    assert report["exact_fixed_radius_identity"]["applies_to_arbitrary_finite_linear_combination"] is True


def test_refined_signed_integrals_obey_preregistered_gate() -> None:
    report = audit.evaluate()
    gate = report["numerical_gate"]
    assert gate["passes"] is True
    assert gate["max_abs_finest_signed_integral"] <= audit.SIGNED_INTEGRAL_TOL
    assert gate["max_medium_to_fine_change"] <= audit.REFINEMENT_CHANGE_TOL


def test_two_channel_witness_is_independent_and_radially_distinct() -> None:
    report = audit.evaluate()
    geometry = report["response_geometry"]
    assert geometry["rank"] == 2
    assert abs(geometry["column_cosine"]) < audit.RANK_COSINE_MAX_ABS
    witness = report["radial_localization_witness"]
    assert witness["annular_zero_at_inner_probe_radii"] is True
    assert witness["annular_active_at_outer_probe_radii"] is True


def test_radial_only_growth_is_not_promoted_as_structural_escape() -> None:
    report = audit.evaluate()
    decision = report["decision"]
    assert decision["independent_two_channel_witness"] is True
    assert decision["fixed_radius_signed_balance_survives_basis_growth"] is True
    assert decision["radial_only_second_channel_can_move_return_debt_between_radii"] is False
    assert decision["radial_localization_only_is_structural_escape"] is False
    assert decision["add_second_radial_localized_channel_now"] is False
    assert decision["closer_visualization_delivery_this_increment"] is False
    assert "1.8<|z|<2.0" in decision["next_if_740_fails"]
