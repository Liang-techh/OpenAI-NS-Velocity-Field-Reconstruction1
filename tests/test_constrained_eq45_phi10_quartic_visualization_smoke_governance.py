from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_phi10_quartic_visualization_smoke_governance import (
    audit_quartic_visualization_smoke_scope,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/eq45_phi10_quartic_visualization_smoke_scope.json"
CONSTRAINTS = ROOT / "configs/constraints.json"
DELIVERY = ROOT / "configs/delivery_state_contract.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(tmp_path: Path, name: str, value: dict) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    return path


def _audit_mutated(tmp_path: Path, contract: dict, *, constraints: dict | None = None):
    return audit_quartic_visualization_smoke_scope(
        contract_path=_write(tmp_path, "contract.json", contract),
        constraints_path=_write(
            tmp_path,
            "constraints.json",
            constraints if constraints is not None else _load(CONSTRAINTS),
        ),
        delivery_state_path=DELIVERY,
    )


def test_current_quartic_visualization_smoke_scope_passes():
    report = audit_quartic_visualization_smoke_scope()
    assert report["velocity_export_ready"] is True
    assert report["visualization_smoke_generated"] is True
    assert report["saved_slice_contains_full_cartesian_velocity_at_each_sample"] is True
    assert report["sampled_y0_max_abs_v"] > 0.0
    assert report["line_overlay_is_true_3d_streamline"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_rejects_promoting_two_component_overlay_to_true_3d_streamline(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["representation_semantics"]["line_overlay_is_true_3d_streamline"] = True
    with pytest.raises(ValueError, match="representation semantics drifted"):
        _audit_mutated(tmp_path, contract)


def test_rejects_promoting_two_component_overlay_to_full_trajectory_projection(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["representation_semantics"]["line_overlay_is_full_3d_trajectory_projection"] = True
    with pytest.raises(ValueError, match="representation semantics drifted"):
        _audit_mutated(tmp_path, contract)


def test_rejects_render_smoke_as_visualization_or_correspondence_acceptance(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="truth-state boundary"):
        _audit_mutated(tmp_path, contract)

    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="truth-state boundary"):
        _audit_mutated(tmp_path, contract)


def test_rejects_render_smoke_as_pde_acceptance_or_export_blocker(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["governed_conclusions"]["smoke_may_establish_pde_validated"] = True
    with pytest.raises(ValueError, match="must remain false"):
        _audit_mutated(tmp_path, contract)

    contract = copy.deepcopy(_load(CONTRACT))
    contract["truth_states"]["velocity_export_ready"] = False
    with pytest.raises(ValueError, match="truth-state boundary"):
        _audit_mutated(tmp_path, contract)


def test_rejects_reclassifying_renderer_choice_as_public_source_fact(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["source_classification"]["two_component_u_w_line_overlay"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification drifted"):
        _audit_mutated(tmp_path, contract)


def test_rejects_candidate_or_overlay_component_drift(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    contract["smoke_binding"]["candidate_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="candidate identity drifted"):
        _audit_mutated(tmp_path, contract)

    contract = copy.deepcopy(_load(CONTRACT))
    contract["smoke_binding"]["line_overlay_components"] = ["u", "v", "w"]
    with pytest.raises(ValueError, match="line overlay components drifted"):
        _audit_mutated(tmp_path, contract)


def test_rejects_posthoc_cr001_threshold_or_energy_tolerance_relaxation(tmp_path: Path):
    contract = copy.deepcopy(_load(CONTRACT))
    constraints = copy.deepcopy(_load(CONSTRAINTS))
    constraints["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_L2 threshold drifted"):
        _audit_mutated(tmp_path, contract, constraints=constraints)

    constraints = copy.deepcopy(_load(CONSTRAINTS))
    constraints["nontriviality"]["reference_energy_abs_tolerance"] = 0.1
    with pytest.raises(ValueError, match="reference energy tolerance drifted"):
        _audit_mutated(tmp_path, contract, constraints=constraints)
