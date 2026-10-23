from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_acceptance_ledger import build_acceptance_ledger

ROOT = Path(__file__).resolve().parents[1]
TASK_STATES = {"CR005": "IN_PROGRESS", "CR009": "IN_PROGRESS", "CR011": "TODO", "CR012": "TODO"}


def _load(path: str):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _current_ledger():
    return build_acceptance_ledger(
        _load("configs/constraints.json"),
        _load("artifacts/constrained/coupled_joint/validation.json"),
        _load("artifacts/constrained/structure_identities.json"),
        candidate_path="artifacts/constrained/coupled_joint/candidate.json",
        task_states=TASK_STATES,
    )


def test_current_coupled_joint_snapshot_is_fail_closed():
    ledger = _current_ledger()
    status = {item["id"]: item for item in ledger["evidence"]}
    assert status["pde_residual_max"]["observed"] == pytest.approx(1.0022183043210762)
    assert status["pde_residual_L2"]["observed"] == pytest.approx(2.2215324099175593)
    assert status["divergence_max"]["observed"] == pytest.approx(0.0036918882656411522)
    assert status["divergence_L2"]["observed"] == pytest.approx(0.0006086624401033801)
    assert status["pde_residual_max"]["status"] == "fail"
    assert status["divergence_max"]["status"] == "fail"
    assert status["core_drift"]["status"] == "pass"
    assert status["structure_identities"]["status"] == "pass"
    assert status["boundary_support_current_candidate"]["status"] == "pending"
    assert status["cr009_generalization"]["status"] == "pending"
    assert ledger["dependency_gate_open"] is False
    assert ledger["acceptance_ready"] is False
    assert ledger["truth_boundary"] == {
        "paper_exact": False,
        "full_blowup_proof": False,
        "small_residual_implies_singularity": False,
    }


def test_checked_snapshot_matches_regeneration():
    assert _current_ledger() == _load("artifacts/constrained/coupled_joint/acceptance_ledger.json")


def test_false_pass_label_is_rejected():
    config = _load("configs/constraints.json")
    validation = _load("artifacts/constrained/coupled_joint/validation.json")
    structure = _load("artifacts/constrained/structure_identities.json")
    mutated = copy.deepcopy(validation)
    mutated["status"] = "passed_validation"
    with pytest.raises(ValueError, match="claims pass"):
        build_acceptance_ledger(
            config,
            mutated,
            structure,
            candidate_path="artifacts/constrained/coupled_joint/candidate.json",
            task_states={"CR005": "DONE", "CR009": "DONE", "CR011": "DONE", "CR012": "DONE"},
            boundary_support_status="pass",
            generalization_status="pass",
            clean_replay_status="pass",
        )
