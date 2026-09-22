from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_current_qs_eta_bridge_scope import (
    _mechanics_eta_bridge_witness,
    audit,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = Path("configs/kokuno_current_qs_eta_bridge_scope.json")
SOURCE = Path(
    "src/openai_ns_reconstruction/kokuno_current_exterior_qs_release2_state.py"
)
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
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_exact_repository_scope_passes() -> None:
    report = audit(ROOT)
    assert report["dependency_pr"] == 1225
    assert report["current_q_s_release2_endpoint_materialized"] is True
    assert report["current_q_s_eta_independence_proven"] is False
    assert report["current_scalar_matching_length_materialized"] is False
    assert report["current_cartesian_matching_bridge_materialized"] is False
    assert report["mechanics_only_matching_time_spread"] > 0.0
    assert (
        report[
            "mechanics_only_max_abs_terminal_mismatch_after_scalarization"
        ]
        > 0.0
    )
    assert report["canonical_eq45_velocity_export_ready"] is True
    assert report["kokuno_unified_global_velocity_export_ready"] is False


def test_eta_resolved_endpoint_cannot_be_scalarized_by_declaration(
    tmp_path: Path,
) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONTRACT,
        lambda p: p["evidence_transfer"].__setitem__(
            "eta_resolved_qs_endpoint_implies_one_scalar_matching_length", True
        ),
    )
    with pytest.raises(AssertionError, match="forbidden evidence transfer"):
        audit(root)


def test_finite_probe_agreement_cannot_be_upgraded_to_exact_eta_independence(
    tmp_path: Path,
) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONTRACT,
        lambda p: p["truth_state"].__setitem__(
            "current_q_s_proven_eta_independent", True
        ),
    )
    with pytest.raises(AssertionError, match="current_q_s_proven_eta_independent"):
        audit(root)


def test_source_module_bridge_promotion_fails_closed(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    source_path = root / SOURCE
    text = source_path.read_text(encoding="utf-8")
    text = text.replace(
        '"current_l_minus_h_matching_bridge_materialized": False',
        '"current_l_minus_h_matching_bridge_materialized": True',
        1,
    )
    source_path.write_text(text, encoding="utf-8")
    with pytest.raises(
        AssertionError,
        match="current-Qs module blob|current_l_minus_h_matching_bridge_materialized",
    ):
        audit(root)


def test_scalar_summary_cannot_define_current_bridge_identity(
    tmp_path: Path,
) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONTRACT,
        lambda p: p["evidence_transfer"].__setitem__(
            "min_max_mean_or_rms_implied_time_may_define_current_bridge_identity",
            True,
        ),
    )
    with pytest.raises(AssertionError, match="forbidden evidence transfer"):
        audit(root)


def test_cr001_threshold_relaxation_is_not_a_bridge_escape(tmp_path: Path) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        CONSTRAINTS,
        lambda p: p["validation"]["thresholds"].__setitem__(
            "pde_residual_max", 0.01
        ),
    )
    with pytest.raises(
        AssertionError,
        match="canonical CR001 constraints blob|momentum max threshold",
    ):
        audit(root)


def test_kokuno_bridge_incompleteness_does_not_revoke_eq45_delivery(
    tmp_path: Path,
) -> None:
    root = _fixture(tmp_path)
    _mutate_json(
        root,
        STATUS,
        lambda p: p["states"].__setitem__("velocity_export_ready", False),
    )
    with pytest.raises(
        AssertionError, match="project status blob|canonical Eq45 velocity_export_ready"
    ):
        audit(root)


def test_mechanics_only_eta_bridge_scalarization_witness() -> None:
    payload = json.loads((ROOT / CONTRACT).read_text(encoding="utf-8"))
    report = _mechanics_eta_bridge_witness(payload)
    times = report["pointwise_matching_time"]
    assert max(times) - min(times) > 0.0
    assert report["representative_eta"] == 0.0
    assert report["max_abs_terminal_mismatch"] > 0.0

    # Logic-only witness: choosing the eta=0 time makes that point hit q_p,
    # while the eta +/- 0.5 endpoints do not both hit q_p.
    q_p = payload["mechanics_witness"]["q_p"]
    terminal = report["terminal_after_representative_scalar_time"]
    assert terminal[1] == pytest.approx(q_p, rel=0.0, abs=1e-15)
    assert terminal[0] != pytest.approx(q_p, rel=0.0, abs=1e-8)
    assert terminal[2] != pytest.approx(q_p, rel=0.0, abs=1e-8)
