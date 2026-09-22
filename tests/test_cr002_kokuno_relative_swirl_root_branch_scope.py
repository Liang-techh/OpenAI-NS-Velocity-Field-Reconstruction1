from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.audit_kokuno_relative_swirl_root_branch_scope import (
    assert_contract,
    audit_contract,
    multi_root_witness,
)
from openai_ns_reconstruction.kokuno_public_relative_swirl_compensator import (
    KokunoPublicRelativeSwirlCompensator,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/kokuno_relative_swirl_root_branch_scope.json"
CONSTRAINTS = ROOT / "configs/constraints.json"
PARENT = ROOT / "src/openai_ns_reconstruction/kokuno_public_relative_swirl_compensator.py"
STATUS = ROOT / "project_status.json"


def test_cr002_relative_swirl_root_branch_contract_passes_exact_parent() -> None:
    assert_contract(CONTRACT, CONSTRAINTS, PARENT, STATUS)


def test_autonomous_mechanics_witness_has_two_exact_roots() -> None:
    payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    witness = payload["mechanics_witness"]
    receipt = multi_root_witness(
        witness["angular_row"],
        witness["pressure_linear_row"],
        witness["pressure_quadratic_matrix"],
        witness["angular_target"],
        witness["root_a"],
        witness["root_b"],
    )
    assert receipt["two_distinct_exact_roots"] is True
    assert receipt["root_distance"] == 2.0


def test_exact_1154_autonomous_rows_have_zero_and_far_zero_target_roots() -> None:
    comp = KokunoPublicRelativeSwirlCompensator()
    a = comp.angular_row_scaled
    b = comp.pressure_linear_row_scaled
    qmat = comp.pressure_quadratic_matrix_scaled
    direction = np.array([-a[1] / a[0], 1.0], dtype=float)
    qa = float(direction @ qmat @ direction)
    qb = float(2.0 * (b @ direction))
    assert abs(qa) > 1e-12
    far = (-qb / qa) * direction
    assert float(np.linalg.norm(far)) > 1e-6
    assert abs(float(a @ far)) < 1e-10
    assert abs(float(2.0 * (b @ far) + far @ qmat @ far)) < 1e-10
    selected = comp.solve(0.0)
    chosen = np.array([float(selected.c1), float(selected.c2)])
    assert np.allclose(chosen, np.zeros(2), rtol=0.0, atol=1e-14)
    assert float(np.sum(np.abs(far))) >= 1.0


def test_mutating_branch_provenance_fails_closed(tmp_path: Path) -> None:
    payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    payload["root_branch_scope"]["current_root_selection_is_repository_autonomous"] = False
    mutated = tmp_path / "contract.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")
    errors = audit_contract(mutated, CONSTRAINTS, PARENT, STATUS)
    assert any("current_root_selection_is_repository_autonomous must remain true" in error for error in errors)


def test_mutating_evidence_transfer_rule_fails_closed(tmp_path: Path) -> None:
    payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    payload["branch_rule_change_promotion"]["reuse_old_coefficient_receipts_after_branch_rule_change"] = True
    mutated = tmp_path / "contract.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")
    errors = audit_contract(mutated, CONSTRAINTS, PARENT, STATUS)
    assert any("reuse_old_coefficient_receipts_after_branch_rule_change must remain false" in error for error in errors)
