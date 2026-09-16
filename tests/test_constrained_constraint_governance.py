from __future__ import annotations

import copy
import json
from pathlib import Path

from openai_ns_reconstruction.constrained_constraint_governance import (
    audit_constraint_directory,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _minimal_config(experiment_id: str, train_seed: int = 11, validation_seed: int = 29) -> dict:
    return {
        "schema_version": 1,
        "experiment_id": experiment_id,
        "nu": 0.01,
        "domain": {"time_interval": [0.25, 0.75]},
        "optimization": {"seed": train_seed},
        "validation": {"seed": validation_seed, "thresholds": {"pde_residual_max": 1e-3}},
    }


def _minimal_governance() -> dict:
    return {
        "schema_version": 1,
        "baseline": "constraints.json",
        "protected_paths": [
            "schema_version",
            "nu",
            "domain.time_interval",
            "validation.thresholds",
        ],
        "seed_separation": {
            "training_path": "optimization.seed",
            "validation_path": "validation.seed",
            "must_differ": True,
        },
        "known_findings": [],
    }


def test_repository_contract_has_no_silent_threshold_or_problem_drift():
    audit = audit_constraint_directory(REPO_ROOT / "configs")
    assert audit.config_count == 10
    assert audit.acceptance_contract_pass
    assert audit.seed_separation_pass
    assert audit.protected_drift == ()

    # This is a deliberately surfaced metadata blocker, not a silent pass.
    assert not audit.metadata_identity_pass
    assert audit.known_findings_match
    assert audit.duplicate_experiment_ids == {
        "tensor_correction_stage1": (
            "constraints_axial_swirl.json",
            "constraints_inner_swirl.json",
            "constraints_localized_swirl.json",
            "constraints_quintic_swirl.json",
            "constraints_temporal_swirl.json",
            "constraints_tensor.json",
        )
    }


def test_threshold_relaxation_is_detected(tmp_path: Path):
    baseline = _minimal_config("base")
    variant = copy.deepcopy(baseline)
    variant["experiment_id"] = "variant"
    variant["validation"]["thresholds"]["pde_residual_max"] = 1.0
    _write_json(tmp_path / "constraints.json", baseline)
    _write_json(tmp_path / "constraints_variant.json", variant)
    _write_json(tmp_path / "constraint_governance.json", _minimal_governance())

    audit = audit_constraint_directory(tmp_path)
    assert not audit.acceptance_contract_pass
    assert [(item.config, item.path) for item in audit.protected_drift] == [
        ("constraints_variant.json", "validation.thresholds")
    ]


def test_training_validation_seed_reuse_is_detected(tmp_path: Path):
    baseline = _minimal_config("base", train_seed=17, validation_seed=17)
    variant = copy.deepcopy(baseline)
    variant["experiment_id"] = "variant"
    _write_json(tmp_path / "constraints.json", baseline)
    _write_json(tmp_path / "constraints_variant.json", variant)
    _write_json(tmp_path / "constraint_governance.json", _minimal_governance())

    audit = audit_constraint_directory(tmp_path)
    assert not audit.seed_separation_pass
    assert not audit.acceptance_contract_pass


def test_unregistered_identity_collision_is_not_silently_accepted(tmp_path: Path):
    baseline = _minimal_config("duplicate")
    variant = copy.deepcopy(baseline)
    _write_json(tmp_path / "constraints.json", baseline)
    _write_json(tmp_path / "constraints_variant.json", variant)
    _write_json(tmp_path / "constraint_governance.json", _minimal_governance())

    audit = audit_constraint_directory(tmp_path)
    assert not audit.metadata_identity_pass
    assert not audit.known_findings_match
    assert audit.duplicate_experiment_ids == {
        "duplicate": ("constraints.json", "constraints_variant.json")
    }
