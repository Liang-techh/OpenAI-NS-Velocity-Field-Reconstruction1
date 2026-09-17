import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_shared_frame_visual_governance import (
    audit_shared_frame_visual_evidence_scope,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/shared_frame_visual_evidence_scope.json"


def _mutated_contract(tmp_path, mutate):
    data = json.loads(CONTRACT.read_text(encoding="utf-8"))
    mutate(data)
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_live_shared_frame_visual_evidence_scope_passes():
    result = audit_shared_frame_visual_evidence_scope()
    assert result["velocity_export_ready"] is True
    assert result["candidate_selection_resolved"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["per_panel_autoscaling"] is False
    assert result["projected_streamlines_are_true_3d_streamlines"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        lambda data: data["comparison_semantics"].__setitem__("per_panel_autoscaling", True),
        lambda data: data["comparison_semantics"].__setitem__(
            "projected_streamlines_are_true_3d_streamlines", True
        ),
        lambda data: data["activity_guard"].__setitem__("scientific_acceptance_threshold", True),
        lambda data: data["activity_guard"].__setitem__("may_replace_cr001_energy_nontriviality", True),
        lambda data: data["promotion_contract"].__setitem__("candidate_selection_resolved", True),
        lambda data: data["promotion_contract"].__setitem__("visual_correspondence_verified", True),
        lambda data: data["promotion_contract"].__setitem__(
            "formal_pde_pass_required_for_candidate_local_save_load_export", True
        ),
        lambda data: data["source_classification"].__setitem__(
            "single_global_speed_normalization", "public_source_fact"
        ),
        lambda data: data["cr001_binding"].__setitem__("pde_residual_L2_threshold", 0.01),
    ],
)
def test_governance_mutations_fail_closed(tmp_path, mutate):
    path = _mutated_contract(tmp_path, mutate)
    with pytest.raises(ValueError):
        audit_shared_frame_visual_evidence_scope(contract_path=path)
