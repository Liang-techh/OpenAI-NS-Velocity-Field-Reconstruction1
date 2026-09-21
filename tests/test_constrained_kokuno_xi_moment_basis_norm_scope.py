from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_xi_moment_basis_norm_scope import (
    MomentBasisNormScopeError,
    _smallness,
    audit_scope,
)

SCOPE_PATH = Path("configs/kokuno_xi_moment_basis_norm_scope.json")
CONSTRAINTS_PATH = Path("configs/constraints.json")


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(tmp_path: Path, name: str, payload: dict) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def test_registered_scope_passes_without_promoting_science() -> None:
    receipt = audit_scope()
    assert receipt["status"] == "PASS"
    assert receipt["current_xi_gain_authorized"] is False
    assert receipt["pde_validated"] is False


def test_mechanics_witness_distinguishes_consistent_vs_d_only_scaling() -> None:
    original = _smallness(1.0, 1.0, 0.2)
    consistent = _smallness(0.1, 0.1, 0.02)
    inconsistent = _smallness(1.0, 1.0, 0.02)
    assert original == pytest.approx(1.6)
    assert consistent == pytest.approx(original)
    assert inconsistent == pytest.approx(0.16)
    assert inconsistent < 1.0 < original


@pytest.mark.parametrize(
    ("mutator", "match"),
    [
        (
            lambda s: s["moment_coordinate_contract"].__setitem__(
                "current_xi_discrepancy_authorized_for_gain_stage_now", True
            ),
            "prematurely authorized|unsafe promotion",
        ),
        (
            lambda s: s["pinned_inputs"]["current_xi_row_transform"].__setitem__(
                "discrepancy_from_complete_ns_defect", True
            ),
            "complete-NS defect",
        ),
        (
            lambda s: s["pinned_inputs"]["existing_gain_firewall"].__setitem__(
                "basis_or_row_normalization_identity_field_present", True
            ),
            "basis-binding gap",
        ),
        (
            lambda s: s["pinned_inputs"]["current_xi_discrepancy"].__setitem__(
                "source_normalization_C", 1.0
            ),
            "outer-normalization identity drift",
        ),
        (
            lambda s: s["delivery_state"].__setitem__("pde_validated", True),
            "unsafe scientific promotion",
        ),
        (
            lambda s: s["delivery_state"].__setitem__(
                "canonical_eq45_velocity_export_ready", False
            ),
            "must not cancel Eq45 delivery",
        ),
    ],
)
def test_scope_mutations_fail_closed(tmp_path: Path, mutator, match: str) -> None:
    scope = _load(SCOPE_PATH)
    mutator(scope)
    mutated = _write(tmp_path, "scope.json", scope)
    with pytest.raises(MomentBasisNormScopeError, match=match):
        audit_scope(mutated, CONSTRAINTS_PATH, require_parent_blob=False)


def test_d_only_rescaling_cannot_be_declared_same_system(tmp_path: Path) -> None:
    scope = _load(SCOPE_PATH)
    promotions = scope["moment_coordinate_contract"]["forbidden_promotions"]
    promotions.remove("rescale discrepancy_d alone => same five-moment correction system")
    mutated = _write(tmp_path, "scope.json", scope)
    with pytest.raises(MomentBasisNormScopeError, match="forbidden-promotion firewall drift"):
        audit_scope(mutated, CONSTRAINTS_PATH, require_parent_blob=False)


def test_future_binding_must_cover_b_q_and_d(tmp_path: Path) -> None:
    scope = _load(SCOPE_PATH)
    requirements = scope["moment_coordinate_contract"]["future_gain_authorization_requires"]
    requirements.remove(
        "one checksum-bound moment_basis_id shared by discrepancy_d, matrix_B and bilinear_Q"
    )
    mutated = _write(tmp_path, "scope.json", scope)
    with pytest.raises(MomentBasisNormScopeError, match="future basis-binding requirements drift"):
        audit_scope(mutated, CONSTRAINTS_PATH, require_parent_blob=False)


def test_autonomous_outer_normalization_cannot_be_laundered_as_public_source(tmp_path: Path) -> None:
    scope = _load(SCOPE_PATH)
    sentence = (
        "current outer normalization C=1000 is a repository realization, "
        "not a recovered hidden source value"
    )
    scope["provenance_classes"]["autonomous_design"].remove(sentence)
    scope["provenance_classes"]["public_source_fact"].append(sentence)
    mutated = _write(tmp_path, "scope.json", scope)
    with pytest.raises(MomentBasisNormScopeError, match="autonomous C=1000 provenance lost"):
        audit_scope(mutated, CONSTRAINTS_PATH, require_parent_blob=False)


def test_cr001_threshold_relaxation_is_rejected(tmp_path: Path) -> None:
    constraints = _load(CONSTRAINTS_PATH)
    constraints["validation"]["thresholds"]["pde_residual_max"] = 2e-3
    mutated_constraints = _write(tmp_path, "constraints.json", constraints)
    with pytest.raises(MomentBasisNormScopeError, match="momentum max gate drift"):
        audit_scope(SCOPE_PATH, mutated_constraints, require_parent_blob=False)


def test_cr001_free_force_shortcut_is_rejected(tmp_path: Path) -> None:
    constraints = _load(CONSTRAINTS_PATH)
    constraints["forcing"]["restriction"] = (
        "Only a,c may be fitted. Residual-dependent pointwise free force allowed."
    )
    mutated_constraints = _write(tmp_path, "constraints.json", constraints)
    with pytest.raises(
        MomentBasisNormScopeError,
        match="residual-defined free-force prohibition missing",
    ):
        audit_scope(SCOPE_PATH, mutated_constraints, require_parent_blob=False)


def test_pending_basis_identity_cannot_be_silently_removed(tmp_path: Path) -> None:
    scope = _load(SCOPE_PATH)
    scope["provenance_classes"]["pending_unknown"].remove(
        "the exact moment-coordinate/basis/row-normalization identity tying d, B and Q together is not yet bound"
    )
    mutated = _write(tmp_path, "scope.json", scope)
    with pytest.raises(MomentBasisNormScopeError, match="pending basis identity silently promoted"):
        audit_scope(mutated, CONSTRAINTS_PATH, require_parent_blob=False)
