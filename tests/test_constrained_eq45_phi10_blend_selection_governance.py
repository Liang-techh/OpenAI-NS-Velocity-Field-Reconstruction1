from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_phi10_blend_selection_governance import (
    audit_phi10_blend_selection_evidence_scope,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/eq45_phi10_blend_selection_evidence_scope.json"
CONSTRAINTS = ROOT / "configs/constraints.json"
RECIPE = ROOT / "artifacts/constrained/eq45_supported_phi10_compact_quartic_blend_recipe.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(tmp_path: Path, name: str, value: dict) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def test_blend_selection_evidence_scope_passes() -> None:
    result = audit_phi10_blend_selection_evidence_scope()
    assert result["velocity_export_ready"] is True
    assert result["blend_weight_selected"] is False
    assert result["candidate_selection_resolved"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["unconsumed_sibling_prs"] == [184, 188, 189]


def test_rejects_source_laundering(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["source_classification"]["compact_quartic_phi10_temporal_interpolation"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification drifted"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))


def test_rejects_threshold_relaxation(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["cr001_binding"]["pde_residual_L2_threshold"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_L2 threshold drifted"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))


def test_rejects_free_force_drift(tmp_path: Path) -> None:
    constraints = copy.deepcopy(_load(CONSTRAINTS))
    constraints["forcing"]["restriction"] = "Fit an arbitrary pointwise force to the residual."
    with pytest.raises(ValueError, match="no-free-force guard drifted"):
        audit_phi10_blend_selection_evidence_scope(constraints_path=_write(tmp_path, "constraints.json", constraints))


def test_rejects_delivery_recipe_selecting_lambda(tmp_path: Path) -> None:
    recipe = copy.deepcopy(_load(RECIPE))
    recipe["representation"]["blend_weight_selected"] = True
    with pytest.raises(ValueError, match="selected a blend weight"):
        audit_phi10_blend_selection_evidence_scope(recipe_path=_write(tmp_path, "recipe.json", recipe))


def test_rejects_single_lane_selection_authority(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["evidence_ledger"]["target_free_morphology"]["may_select_blend_weight"] = True
    with pytest.raises(ValueError, match="may not select blend weight"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))


def test_rejects_public_comparison_bypass(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["selection_contract"]["public_observable_comparison_required_for_visual_correspondence_selection"] = False
    with pytest.raises(ValueError, match="must remain true"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))


def test_rejects_visual_state_promotion(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="truth-state boundary drifted"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))


def test_rejects_pde_state_promotion(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-state boundary drifted"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))


def test_pde_pending_does_not_become_export_blocker(tmp_path: Path) -> None:
    contract = copy.deepcopy(_load(CONTRACT))
    contract["selection_contract"]["formal_pde_pass_required_for_candidate_local_save_load_export"] = True
    with pytest.raises(ValueError, match="must remain false"):
        audit_phi10_blend_selection_evidence_scope(contract_path=_write(tmp_path, "contract.json", contract))
