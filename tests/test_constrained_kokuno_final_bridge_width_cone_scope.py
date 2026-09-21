from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_final_bridge_width_cone_scope import validate_scope


ROOT = Path(__file__).resolve().parents[1]
SCOPE_PATH = ROOT / "configs" / "kokuno_final_bridge_width_cone_scope.json"
PARENT_PATH = ROOT / "src" / "openai_ns_reconstruction" / "kokuno_a5_actual_final_bridge_xi110_ingest_contract.py"
CONSTRAINTS_PATH = ROOT / "configs" / "constraints.json"


def _fixture() -> tuple[dict, str, dict]:
    scope = json.loads(SCOPE_PATH.read_text(encoding="utf-8"))
    parent = PARENT_PATH.read_text(encoding="utf-8")
    constraints = json.loads(CONSTRAINTS_PATH.read_text(encoding="utf-8"))
    return scope, parent, constraints


def _must_fail(scope: dict, parent: str, constraints: dict) -> None:
    with pytest.raises(ValueError):
        validate_scope(scope, parent, constraints)


def test_registered_scope_passes_without_scientific_promotion() -> None:
    scope, parent, constraints = _fixture()
    result = validate_scope(scope, parent, constraints)
    assert result["status"] == "pass"
    assert result["scientific_promotion"] is False
    assert scope["registered_realization"]["widths_are_repository_autonomous"] is True
    assert scope["evidence_scope"]["final_bridge_cone_admissibility_independently_certified"] is False
    assert scope["state_boundary"]["canonical_eq45_velocity_export_ready"] is True
    assert scope["state_boundary"]["kokuno_velocity_export_ready"] is False


def test_mechanics_witness_distinguishes_existence_from_selected_width_admission() -> None:
    scope, _, _ = _fixture()
    witness = scope["mechanics_witness"]
    lo, hi = witness["toy_allowed_interval"]
    bound = witness["toy_admissible_upper_bound"]
    existing = witness["toy_existing_admissible_width"]
    selected = witness["toy_autonomous_selected_width"]
    assert lo < existing < bound
    assert lo < selected < hi
    assert selected > bound
    # This is logic-only.  No toy number is a Kokuno/OpenAI source bound.
    assert "not a Kokuno/OpenAI bound" in witness["scope"]


@pytest.mark.parametrize(
    "mutator",
    [
        lambda s: s["provenance"]["public_source_fact"].append(
            "the public source fixes both final-bridge logarithmic widths to 0.02"
        ),
        lambda s: s["registered_realization"].__setitem__("widths_recovered_from_public_source", True),
        lambda s: s["evidence_scope"].__setitem__("geometric_width_fit_implies_source_cone_admission", True),
        lambda s: s["evidence_scope"].__setitem__("profile_derivative_consistency_implies_source_cone_admission", True),
        lambda s: s["evidence_scope"].__setitem__("endpoint_handoff_consistency_implies_source_cone_admission", True),
        lambda s: s["evidence_scope"].__setitem__("autonomous_widths_identify_source_hidden_widths", True),
        lambda s: s["evidence_scope"].__setitem__("final_bridge_cone_admissibility_independently_certified", True),
        lambda s: s["evidence_scope"].__setitem__("source_admitted_global_kappa0_smallness", True),
        lambda s: s["evidence_scope"].__setitem__("agent4_exact_head_scientific_admission_observed", True),
        lambda s: s["state_boundary"].__setitem__("kokuno_velocity_export_ready", True),
        lambda s: s["state_boundary"].__setitem__("pde_validated", True),
        lambda s: s["state_boundary"].__setitem__("paper_exact", True),
        lambda s: s["state_boundary"].__setitem__("openai_field_identified", True),
        lambda s: s["state_boundary"].__setitem__("pending_kokuno_cone_or_pde_work_blocks_canonical_callable_velocity_delivery", True),
        lambda s: s["cr001_lock"].__setitem__("momentum_max_gate", 1.0e-2),
        lambda s: s["cr001_lock"].__setitem__("divergence_volume_l2_gate", 1.0e-4),
        lambda s: s["cr001_lock"].__setitem__("residual_defined_free_force_forbidden", False),
    ],
)
def test_truth_boundary_mutations_fail_closed(mutator) -> None:
    scope, parent, constraints = _fixture()
    mutated = copy.deepcopy(scope)
    mutator(mutated)
    _must_fail(mutated, parent, constraints)


def test_missing_pending_source_width_identity_fails_closed() -> None:
    scope, parent, constraints = _fixture()
    mutated = copy.deepcopy(scope)
    mutated["provenance"]["pending_unknown"] = [
        item
        for item in mutated["provenance"]["pending_unknown"]
        if "source-hidden numerical values" not in item
    ]
    _must_fail(mutated, parent, constraints)


def test_parent_truth_boundary_promotion_is_detected() -> None:
    scope, parent, constraints = _fixture()
    promoted_parent = parent.replace(
        '"final_bridge_cone_admissibility_independently_certified": False',
        '"final_bridge_cone_admissibility_independently_certified": True',
        1,
    )
    _must_fail(scope, promoted_parent, constraints)


def test_canonical_threshold_relaxation_is_detected() -> None:
    scope, parent, constraints = _fixture()
    relaxed = copy.deepcopy(constraints)
    relaxed["validation"]["thresholds"]["pde_residual_L2"] = 1.0e-2
    _must_fail(scope, parent, relaxed)
