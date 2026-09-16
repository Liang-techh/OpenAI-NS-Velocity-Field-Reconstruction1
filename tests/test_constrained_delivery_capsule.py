import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_delivery_capsule import (
    build_delivery_capsule,
    validate_delivery_capsule,
)


def _write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _fixture_repo(root: Path, *, validation_status="failed_validation") -> Path:
    candidate = {"family": "coupled_velocity_v1", "status": "candidate", "coefficients": [0.2]}
    _write_json(root / "artifacts/constrained/coupled_joint/candidate.json", candidate)
    _write_json(root / "src/openai_ns_reconstruction/data/velocity_candidate.json", candidate)
    (root / "src/openai_ns_reconstruction/velocity_components.py").parent.mkdir(parents=True, exist_ok=True)
    (root / "src/openai_ns_reconstruction/velocity_components.py").write_text("def velocity(x,y,z,t): return x,y,z\n")
    (root / "docs").mkdir(parents=True, exist_ok=True)
    (root / "docs/VELOCITY_API.md").write_text("velocity API\n")
    (root / "docs/STRUCTURE_IDENTITIES.md").write_text("structure\n")
    _write_json(
        root / "configs/constraints_coupled.json",
        {
            "units": "dimensionless",
            "domain": {
                "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
                "time_interval": [0.25, 0.75],
            },
        },
    )
    _write_json(
        root / "artifacts/constrained/coupled_joint/training.json",
        {"seed": 20260916, "calls": 29, "max_nfev": 80},
    )
    _write_json(
        root / "artifacts/constrained/coupled_joint/validation.json",
        {
            "status": validation_status,
            "candidate": "artifacts/constrained/coupled_joint/candidate.json",
            "seed": 914027,
            "points": 4096,
        },
    )
    _write_json(
        root / "artifacts/constrained/structure_identities.json",
        {"status": "five_symbolic_identities_verified"},
    )
    grid = root / "artifacts/visual/velocity_api/grid.npz"
    grid.parent.mkdir(parents=True, exist_ok=True)
    grid.write_bytes(b"reference-grid")
    _write_json(grid.parent / "metadata.json", {"scope": "reference grid"})
    _write_json(root / "project_status.json", {"visual_correspondence": "not_yet_verified"})
    return root


def test_capsule_records_reproducibility_without_promoting_validation(tmp_path):
    root = _fixture_repo(tmp_path)
    capsule = build_delivery_capsule(root)
    assert capsule["candidate"]["family"] == "coupled_velocity_v1"
    assert capsule["reproducibility"]["training_seed"] == 20260916
    assert capsule["reproducibility"]["validation_seed"] == 914027
    assert capsule["recommended_evaluation"]["reference_times"] == [0.25, 0.5, 0.75]
    assert capsule["claim_states"]["velocity_export_ready"] is True
    assert capsule["claim_states"]["visualization_ready"] is False
    assert capsule["claim_states"]["pde_validated"] is False
    assert capsule["claim_states"]["paper_exact"] is False
    assert capsule["final_artifact_complete"] is False
    assert capsule["inventory"]["matlab_export_entry"]["status"] == "pending"


def test_capsule_rejects_candidate_drift_and_seed_reuse(tmp_path):
    root = _fixture_repo(tmp_path)
    _write_json(
        root / "src/openai_ns_reconstruction/data/velocity_candidate.json",
        {"family": "coupled_velocity_v1", "status": "candidate", "coefficients": [0.3]},
    )
    with pytest.raises(ValueError, match="candidate bytes differ"):
        build_delivery_capsule(root)

    root = _fixture_repo(tmp_path / "second")
    _write_json(
        root / "artifacts/constrained/coupled_joint/validation.json",
        {
            "status": "failed_validation",
            "candidate": "artifacts/constrained/coupled_joint/candidate.json",
            "seed": 20260916,
            "points": 4096,
        },
    )
    with pytest.raises(ValueError, match="seeds must be distinct"):
        build_delivery_capsule(root)


def test_capsule_rejects_stale_validation_candidate_path(tmp_path):
    root = _fixture_repo(tmp_path)
    _write_json(
        root / "artifacts/constrained/coupled_joint/validation.json",
        {
            "status": "failed_validation",
            "candidate": "artifacts/constrained/older/candidate.json",
            "seed": 914027,
            "points": 4096,
        },
    )
    with pytest.raises(ValueError, match="does not point"):
        build_delivery_capsule(root)


def test_validator_cannot_promote_pde_or_paper_exact(tmp_path):
    capsule = build_delivery_capsule(_fixture_repo(tmp_path))
    capsule["claim_states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        validate_delivery_capsule(capsule)

    capsule = build_delivery_capsule(_fixture_repo(tmp_path / "second"))
    capsule["claim_states"]["paper_exact"] = True
    with pytest.raises(ValueError, match="paper_exact"):
        validate_delivery_capsule(capsule)


def test_current_repository_delivery_capsule_is_truth_bounded():
    root = Path(__file__).resolve().parents[1]
    capsule = build_delivery_capsule(root)
    assert capsule["candidate"]["family"] == "coupled_velocity_v1"
    assert capsule["reproducibility"]["training_seed"] != capsule["reproducibility"]["validation_seed"]
    assert capsule["claim_states"]["velocity_export_ready"] is True
    assert capsule["claim_states"]["pde_validated"] is False
    assert capsule["claim_states"]["paper_exact"] is False
    assert capsule["current_validation_label"] == "failed_validation"
