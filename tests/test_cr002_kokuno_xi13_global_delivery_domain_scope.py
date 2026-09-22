from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_xi13_global_delivery_domain_scope import (
    EXPECTED_CONSTRAINTS_BLOB,
    EXPECTED_PARENT_BLOB,
    assert_contract,
    audit_contract,
    prefix_nonuniqueness_witness,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "kokuno_xi13_global_delivery_domain_scope.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"
CANDIDATE = ROOT / "src" / "openai_ns_reconstruction" / "kokuno_pa16_current_cartesian_pulse_end_compensated.py"
PROJECT_STATUS = ROOT / "project_status.json"


def _write_json(tmp_path: Path, payload: dict, name: str = "mutated.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def _contract() -> dict:
    return json.loads(CONTRACT.read_text())


def test_prefix_nonuniqueness_witness_is_nonvacuous() -> None:
    receipt = prefix_nonuniqueness_witness(13.0, 13.5)
    assert receipt["same_prefix_value"] is True
    assert receipt["same_endpoint_value"] is True
    assert receipt["different_outer_value"] is True
    assert receipt["outer_delta"] > 0.0


def test_prefix_nonuniqueness_witness_rejects_invalid_probe() -> None:
    with pytest.raises(ValueError, match="strictly beyond"):
        prefix_nonuniqueness_witness(13.0, 13.0)


def test_repository_contract_is_fail_closed() -> None:
    assert_contract(CONTRACT, CONSTRAINTS, CANDIDATE, PROJECT_STATUS)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("stage_delivery_scope", "xi13_is_canonical_physical_support_endpoint"), True),
        (("stage_delivery_scope", "xi13_is_canonical_evaluation_box_boundary"), True),
        (("stage_delivery_scope", "terminal_exterior_tail_materialized"), True),
        (("stage_delivery_scope", "outer_global_leading_velocity_materialized"), True),
        (("stage_delivery_scope", "global_project_domain_coverage_verified"), True),
        (("stage_delivery_scope", "unified_global_kokuno_velocity_export_ready"), True),
        (("stage_delivery_scope", "matched_global_pressure_materialized"), True),
        (("stage_delivery_scope", "restricted_forcing_materialized"), True),
        (("stage_delivery_scope", "heldout_complete_ns_residual_assessed"), True),
        (("independent_delivery_states", "kokuno_global_velocity_export_ready"), True),
        (("independent_delivery_states", "kokuno_visual_correspondence_verified"), True),
        (("independent_delivery_states", "kokuno_pde_validated"), True),
        (("independent_delivery_states", "kokuno_paper_exact"), True),
        (("independent_delivery_states", "kokuno_openai_field_identified"), True),
        (("future_extension_promotion_rule", "reuse_xi13_prefix_receipts_as_outer_domain_evidence"), True),
        (("future_extension_promotion_rule", "new_identity_required_for_terminal_or_exterior_extension"), False),
        (("future_extension_promotion_rule", "fresh_save_load_audit_required"), False),
        (("future_extension_promotion_rule", "fresh_finite_cartesian_domain_audit_required"), False),
        (("future_extension_promotion_rule", "fresh_downstream_numerical_audit_required"), False),
        (("cr001_invariants", "forcing_family_change_this_increment"), True),
        (("cr001_invariants", "threshold_relaxation_this_increment"), True),
    ],
)
def test_truth_and_scope_promotions_fail_closed(tmp_path: Path, path: tuple[str, str], value: object) -> None:
    payload = _contract()
    payload[path[0]][path[1]] = value
    errors = audit_contract(_write_json(tmp_path, payload), CONSTRAINTS, CANDIDATE, PROJECT_STATUS)
    assert errors


def test_public_pulse_endpoint_cannot_be_relabelled(tmp_path: Path) -> None:
    payload = _contract()
    payload["stage_delivery_scope"]["public_source_pulse_endpoint_xi"] = 14.0
    errors = audit_contract(_write_json(tmp_path, payload), CONSTRAINTS, CANDIDATE, PROJECT_STATUS)
    assert any("public pulse endpoint" in error for error in errors)


def test_four_way_provenance_is_mandatory(tmp_path: Path) -> None:
    payload = _contract()
    payload["four_way_provenance"].pop("pending_unknown")
    errors = audit_contract(_write_json(tmp_path, payload), CONSTRAINTS, CANDIDATE, PROJECT_STATUS)
    assert any("four_way_provenance" in error for error in errors)


def test_forbidden_implication_set_is_mandatory(tmp_path: Path) -> None:
    payload = _contract()
    payload["evidence_transfer_firewall"]["forbidden_implications"] = [
        "source-coordinate pulse endpoint xi=13 -> terminal/exterior/global leading completion"
    ]
    errors = audit_contract(_write_json(tmp_path, payload), CONSTRAINTS, CANDIDATE, PROJECT_STATUS)
    assert any("forbidden_implications" in error for error in errors)


def test_cr001_threshold_mutation_is_rejected(tmp_path: Path) -> None:
    constraints = json.loads(CONSTRAINTS.read_text())
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    mutated = _write_json(tmp_path, constraints, "constraints.json")
    errors = audit_contract(CONTRACT, mutated, CANDIDATE, PROJECT_STATUS)
    assert any(
        "canonical configs/constraints.json blob changed" in error
        or "live momentum max gate" in error
        for error in errors
    )


def test_canonical_eq45_delivery_cannot_be_downgraded(tmp_path: Path) -> None:
    status = json.loads(PROJECT_STATUS.read_text())
    status["states"]["velocity_export_ready"] = False
    mutated = _write_json(tmp_path, status, "project_status.json")
    errors = audit_contract(CONTRACT, CONSTRAINTS, CANDIDATE, mutated)
    assert any("velocity_export_ready" in error for error in errors)


def test_parent_blobs_are_pinned() -> None:
    payload = _contract()
    assert payload["exact_parent"]["candidate_module_blob_sha"] == EXPECTED_PARENT_BLOB
    assert payload["cr001_invariants"]["constraints_blob_sha"] == EXPECTED_CONSTRAINTS_BLOB
