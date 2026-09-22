from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_source_ideal_qs_evidence_audit import audit

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = Path("configs/kokuno_source_ideal_qs_evidence_scope.json")
SOURCE = Path("src/openai_ns_reconstruction/kokuno_public_ideal_exterior_qs_schedule.py")
CONSTRAINTS = Path("configs/constraints.json")
STATUS = Path("project_status.json")


def _fixture(tmp_path: Path) -> Path:
    for rel in (CONTRACT, SOURCE, CONSTRAINTS, STATUS):
        target = tmp_path / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, target)
    return tmp_path


def _mutate_json(root: Path, rel: Path, mutate) -> None:
    path = root / rel
    payload = json.loads(path.read_text(encoding="utf-8"))
    mutate(payload)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_exact_repository_scope_passes() -> None:
    report = audit(ROOT)
    assert report["dependency_pr"] == 1217
    assert report["source_ideal_reference_materialized"] is True
    assert report["current_candidate_qs_materialized"] is False
    assert report["source_ideal_to_current_transfer_allowed"] is False
    assert report["canonical_eq45_velocity_export_ready"] is True
    assert report["kokuno_unified_global_velocity_export_ready"] is False


def test_source_ideal_tmatch_cannot_be_promoted_to_current_bridge(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONTRACT,
        lambda p: p["evidence_transfer"].__setitem__(
            "source_ideal_t_match_may_be_used_as_current_cartesian_bridge_length", True
        ),
    )
    with pytest.raises(AssertionError, match="forbidden evidence transfer"):
        audit(root)


def test_current_candidate_promotion_without_materialization_fails(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONTRACT,
        lambda p: p["truth_state"].__setitem__(
            "current_candidate_qs_release2_endpoint_materialized", True
        ),
    )
    with pytest.raises(AssertionError, match="current_candidate_qs_release2_endpoint_materialized"):
        audit(root)


def test_source_module_cannot_silently_promote_current_global_truth(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    source_path = root / SOURCE
    text = source_path.read_text(encoding="utf-8")
    text = text.replace(
        '"current_q_s_release2_endpoint_materialized": False',
        '"current_q_s_release2_endpoint_materialized": True',
        1,
    )
    source_path.write_text(text, encoding="utf-8")
    with pytest.raises(AssertionError, match="source-ideal module blob|current_q_s_release2_endpoint_materialized"):
        audit(root)


def test_cr001_threshold_relaxation_is_not_a_way_around_the_bridge(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONSTRAINTS,
        lambda p: p["validation"]["thresholds"].__setitem__("pde_residual_max", 0.01),
    )
    with pytest.raises(AssertionError, match="canonical CR001 constraints blob|momentum max threshold"):
        audit(root)


def test_kokuno_incompleteness_does_not_revoke_canonical_eq45_delivery(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        STATUS,
        lambda p: p["states"].__setitem__("velocity_export_ready", False),
    )
    with pytest.raises(AssertionError, match="canonical Eq45 velocity_export_ready"):
        audit(root)


def test_same_number_is_not_same_state_identity() -> None:
    # Autonomous mechanics-only witness: identical scalar output can arise from
    # states with different carried M history, so scalar equality cannot prove
    # source-ideal/current-candidate identity.
    ideal = {"q_s": 0.25, "M": 0.0}
    carried = {"q_s": 0.25, "M": 0.125}
    assert ideal["q_s"] == carried["q_s"]
    assert ideal != carried
