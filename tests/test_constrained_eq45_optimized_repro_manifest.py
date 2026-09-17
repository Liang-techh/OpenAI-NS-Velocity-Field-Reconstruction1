import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_optimized_repro_manifest import (
    MANIFEST_SCHEMA,
    audit_optimized_repro_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "artifacts" / "constrained" / "eq45_profile_force_optimized_repro_manifest.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"
OPTIMIZATION = ROOT / "artifacts" / "constrained" / "eq45_profile_force_optimization.json"
CANDIDATE = ROOT / "artifacts" / "constrained" / "eq45_profile_force_optimized_candidate.json"
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _payload(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_optimized_repro_manifest_binds_candidate_optimizer_and_project_contract():
    report = audit_optimized_repro_manifest(MANIFEST, CONSTRAINTS, OPTIMIZATION, CANDIDATE, SEED)

    assert report["schema"] == MANIFEST_SCHEMA
    assert report["candidate_sha256"] == "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"
    assert report["velocity_probe_rms"] > 0.0
    assert report["fit_seed"] == 20260917
    assert report["diagnostic_holdout_seed"] == 914117
    assert report["project_validation_seed"] == 914027
    assert len({report["fit_seed"], report["diagnostic_holdout_seed"], report["project_validation_seed"]}) == 3
    assert report["velocity_export_ready"] is True
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False


def test_manifest_mutations_fail_closed(tmp_path):
    base = _payload(MANIFEST)
    mutations = []

    wrong_sha = copy.deepcopy(base)
    wrong_sha["candidate"]["sha256"] = "0" * 64
    mutations.append(wrong_sha)

    wrong_fit_seed = copy.deepcopy(base)
    wrong_fit_seed["optimization_receipt"]["fit_seed"] = 20260918
    mutations.append(wrong_fit_seed)

    conflated_seed = copy.deepcopy(base)
    conflated_seed["project_contract"]["project_validation_seed"] = 914117
    mutations.append(conflated_seed)

    relaxed_threshold = copy.deepcopy(base)
    relaxed_threshold["project_contract"]["validation_thresholds"]["pde_residual_L2"] = 0.01
    mutations.append(relaxed_threshold)

    promoted_pde = copy.deepcopy(base)
    promoted_pde["delivery_state"]["pde_validated"] = True
    mutations.append(promoted_pde)

    promoted_visual = copy.deepcopy(base)
    promoted_visual["delivery_state"]["visualization_ready"] = True
    mutations.append(promoted_visual)

    claimed_frame_mapping = copy.deepcopy(base)
    claimed_frame_mapping["reference_visualization"]["public_openai_frame_time_mapping"] = "matched"
    mutations.append(claimed_frame_mapping)

    for index, payload in enumerate(mutations):
        path = tmp_path / f"manifest_mutation_{index}.json"
        _write(path, payload)
        with pytest.raises(ValueError):
            audit_optimized_repro_manifest(path, CONSTRAINTS, OPTIMIZATION, CANDIDATE, SEED)


def test_manifest_detects_underlying_receipt_or_constraint_drift(tmp_path):
    optimization = _payload(OPTIMIZATION)
    optimization["contract"]["holdout_seed"] = 914118
    drifted_optimization = tmp_path / "optimization.json"
    _write(drifted_optimization, optimization)
    with pytest.raises(ValueError, match="diagnostic_holdout_seed"):
        audit_optimized_repro_manifest(MANIFEST, CONSTRAINTS, drifted_optimization, CANDIDATE, SEED)

    constraints = _payload(CONSTRAINTS)
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    drifted_constraints = tmp_path / "constraints.json"
    _write(drifted_constraints, constraints)
    with pytest.raises(ValueError, match="project validation thresholds"):
        audit_optimized_repro_manifest(MANIFEST, drifted_constraints, OPTIMIZATION, CANDIDATE, SEED)
