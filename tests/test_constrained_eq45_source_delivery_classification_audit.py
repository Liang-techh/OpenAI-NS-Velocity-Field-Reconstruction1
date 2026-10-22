from __future__ import annotations

import json
from pathlib import Path
import shutil

import pytest

from openai_ns_reconstruction.constrained_eq45_source_delivery_classification_audit import (
    audit_eq45_source_delivery_classification,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILES = (
    "eq45_source_delivery_classification_bridge.json",
    "eq45_source_contract.json",
    "velocity_delivery_contract.json",
    "constraints.json",
)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    configs = tmp_path / "configs"
    configs.mkdir(parents=True)
    for name in CONFIG_FILES:
        shutil.copy2(REPO_ROOT / "configs" / name, configs / name)
    shutil.copy2(REPO_ROOT / "project_status.json", tmp_path / "project_status.json")
    return tmp_path


def _classification_row(payload: dict, prefix: str) -> dict:
    matches = [
        row
        for row in payload["source_classification"]
        if row["item"].startswith(prefix)
    ]
    assert len(matches) == 1
    return matches[0]


def test_live_eq45_source_delivery_classification_bridge_passes() -> None:
    receipt = audit_eq45_source_delivery_classification(REPO_ROOT)
    assert receipt["legacy_source_labels_mapped_to_canonical_four_classes"] is True
    assert receipt["public_source_and_autonomous_numeric_instantiation_kept_distinct"] is True
    assert receipt["callable_delivery_kept_independent_from_pde_and_exactness"] is True
    assert receipt["cr001_checked_unchanged"] is True


def test_bridge_uses_exact_canonical_four_classes_and_forbids_zero_profile_shortcut() -> None:
    contract = _read_json(
        REPO_ROOT / "configs" / "eq45_source_delivery_classification_bridge.json"
    )
    assert set(contract["canonical_classification_vocabulary"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert contract["legacy_eq45_source_label_map"] == {
        "user_required": "user_requirement",
        "public_source": "public_source_fact",
        "autonomous": "autonomous_design",
        "pending": "pending_unknown",
    }
    forbidden = "\n".join(contract["forbidden_inferences"])
    assert "Do not set unresolved public-source profile data to zero" in forbidden
    assert "does not block" not in forbidden.lower()  # prose is not used as a hidden promotion gate
    assert contract["source_to_delivery_bridge"][
        "callable_delivery_requires_public_source_profile_identification"
    ] is False
    assert contract["source_to_delivery_bridge"]["callable_delivery_requires_pde_validation"] is False


def test_rejects_legacy_source_label_laundering(tmp_path: Path) -> None:
    root = _minimal_repo(tmp_path)
    path = root / "configs" / "eq45_source_contract.json"
    payload = _read_json(path)
    row = _classification_row(
        payload,
        "finite numerical choices such as h, sigma, radial degree, eta nodes, fit weights and visualization tolerances",
    )
    row["classification"] = "public_source"
    _write_json(path, payload)
    with pytest.raises(AssertionError, match="legacy classification"):
        audit_eq45_source_delivery_classification(root)


def test_rejects_public_source_profile_recovery_promotion(tmp_path: Path) -> None:
    root = _minimal_repo(tmp_path)
    path = root / "configs" / "eq45_source_contract.json"
    payload = _read_json(path)
    payload["claim_status"]["final_profiles_numerically_identified"] = True
    _write_json(path, payload)
    with pytest.raises(AssertionError, match="public-source claim promoted"):
        audit_eq45_source_delivery_classification(root)


def test_rejects_autonomous_candidate_relabelled_as_public_source(tmp_path: Path) -> None:
    root = _minimal_repo(tmp_path)
    path = root / "configs" / "velocity_delivery_contract.json"
    payload = _read_json(path)
    row = _classification_row(
        payload,
        "Eq45 support-connected candidate identity, API, export format, and delivery wrapper",
    )
    row["classification"] = "public_source_fact"
    _write_json(path, payload)
    with pytest.raises(AssertionError, match="canonical candidate provenance class"):
        audit_eq45_source_delivery_classification(root)


def test_rejects_delivery_to_exact_openai_promotion(tmp_path: Path) -> None:
    root = _minimal_repo(tmp_path)
    path = root / "configs" / "eq45_source_delivery_classification_bridge.json"
    payload = _read_json(path)
    payload["source_to_delivery_bridge"]["canonical_candidate_may_be_called_exact_openai_field"] = True
    _write_json(path, payload)
    with pytest.raises(AssertionError, match="exact OpenAI field"):
        audit_eq45_source_delivery_classification(root)


def test_rejects_project_truth_state_laundering(tmp_path: Path) -> None:
    root = _minimal_repo(tmp_path)
    path = root / "project_status.json"
    payload = _read_json(path)
    payload["states"]["visual_correspondence_verified"] = True
    _write_json(path, payload)
    with pytest.raises(AssertionError, match="project truth state visual_correspondence_verified"):
        audit_eq45_source_delivery_classification(root)


def test_rejects_cr001_threshold_drift(tmp_path: Path) -> None:
    root = _minimal_repo(tmp_path)
    path = root / "configs" / "constraints.json"
    payload = _read_json(path)
    payload["validation"]["thresholds"]["pde_residual_max"] = 0.002
    _write_json(path, payload)
    with pytest.raises(AssertionError, match="pde_residual_max"):
        audit_eq45_source_delivery_classification(root)
