from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_acceptance import (
    build_supported_acceptance_state,
)
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)

ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
SNAPSHOT = ROOT / "artifacts" / "constrained" / "eq45_supported_acceptance_state.json"


def _load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _inputs():
    child = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.load_json(SEED))
    return (
        child,
        _load("configs/constraints.json"),
        _load("artifacts/constrained/eq45_supported_energy_audit.json"),
        _load("configs/delivery_state_contract.json"),
    )


def _state() -> dict:
    return build_supported_acceptance_state(*_inputs())


def test_supported_acceptance_snapshot_matches_checked_evidence() -> None:
    state = _state()
    checked = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    assert state == checked

    evidence = {item["id"]: item for item in state["evidence"]}
    assert state["candidate_sha256"] == (
        "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
    )
    assert evidence["velocity_export_ready"]["status"] == "pass"
    assert evidence["cr001_reference_energy"]["status"] == "fail"
    assert evidence["energy_quadrature_convergence"]["status"] == "pass"
    assert evidence["validation_time_energy_range"]["status"] == "pass"
    assert evidence["supported_child_divergence_validation"]["status"] == "pending"
    assert evidence["supported_child_pde_validation"]["status"] == "pending"
    assert evidence["visualization_ready"]["status"] == "pending"
    assert state["states"]["velocity_export_ready"] is True
    assert state["states"]["visualization_ready"] is False
    assert state["states"]["pde_validated"] is False
    assert state["acceptance_ready"] is False
    assert state["interpretation"]["energy_failure_blocks_velocity_export"] is False


def test_false_reference_energy_pass_label_is_rejected() -> None:
    child, constraints, energy, delivery = _inputs()
    mutated = deepcopy(energy)
    mutated["reference_energy_gate_pass"] = True
    with pytest.raises(ValueError, match="inconsistent"):
        build_supported_acceptance_state(child, constraints, mutated, delivery)


def test_reference_energy_threshold_drift_is_rejected() -> None:
    child, constraints, energy, delivery = _inputs()
    mutated = deepcopy(constraints)
    mutated["nontriviality"]["reference_energy_abs_tolerance"] = 0.1
    with pytest.raises(ValueError, match="reference tolerance"):
        build_supported_acceptance_state(child, mutated, energy, delivery)


def test_delivery_contract_cannot_merge_pde_failure_into_export_blocking() -> None:
    child, constraints, energy, delivery = _inputs()
    mutated = deepcopy(delivery)
    mutated["forbidden_inferences"] = [
        item
        for item in mutated["forbidden_inferences"]
        if item.get("from") != ["pde_validation_failed"]
    ]
    with pytest.raises(ValueError, match="PDE-failure/export separation"):
        build_supported_acceptance_state(child, constraints, energy, mutated)


def test_missing_independent_delivery_state_is_rejected() -> None:
    child, constraints, energy, delivery = _inputs()
    mutated = deepcopy(delivery)
    del mutated["states"]["visualization_ready"]
    with pytest.raises(ValueError, match="visualization_ready"):
        build_supported_acceptance_state(child, constraints, energy, mutated)
