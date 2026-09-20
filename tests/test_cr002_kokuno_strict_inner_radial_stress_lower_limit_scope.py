from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_strict_inner_radial_stress_lower_limit_scope import (
    AuditError,
    audit,
)

CONTRACT_REL = Path("configs/kokuno_strict_inner_radial_stress_lower_limit_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
TARGET_REL = Path(
    "src/openai_ns_reconstruction/kokuno_strict_inner_transport_radial_stress.py"
)
HELPER_REL = Path(
    "src/openai_ns_reconstruction/kokuno_actual_oscillatory_mean_stress.py"
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _sandbox(tmp_path: Path) -> Path:
    root = _repo_root()
    for rel in (CONTRACT_REL, CONSTRAINTS_REL, TARGET_REL, HELPER_REL):
        destination = tmp_path / rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / rel, destination)
    return tmp_path


def _mutate_json(path: Path, mutator) -> None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    mutator(payload)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_cr002_radial_stress_lower_limit_scope_happy_path() -> None:
    report = audit(_repo_root())
    assert report["ok"] is True
    assert report["target_pr"] == 858
    assert report["formal_lower_limit"] == 0.0
    assert (
        report["realized_lower_limit"]
        == "r_min = first supplied strictly-positive radial sample"
    )
    assert report["formal_axis_based_inverse_verified"] is False
    assert report["pde_validated"] is False
    assert report["pde_pending_blocks_callable_velocity_delivery"] is False


def test_manufactured_lower_limit_witness_is_deterministic() -> None:
    report = audit(_repo_root())
    witness = report["manufactured_formal_sigma_at_rmin"]
    assert witness["e1"] == pytest.approx(-0.1, abs=1e-15)
    assert witness["e2"] == pytest.approx(-1.0 / 15.0, abs=1e-15)


@pytest.mark.parametrize(
    ("name", "mutator"),
    [
        (
            "realized_lower_limit_laundered_to_axis",
            lambda data: data["radial_inverse_scope"].__setitem__(
                "realized_discrete_lower_limit", "0"
            ),
        ),
        (
            "inner_edge_zero_promoted_to_axis_regularity",
            lambda data: data["radial_inverse_scope"].__setitem__(
                "stress_inner_edge_zero_is_axis_regularity_evidence", True
            ),
        ),
        (
            "sampled_moment_promoted_to_full_domain",
            lambda data: data["radial_inverse_scope"].__setitem__(
                "sampled_weighted_moment_is_formal_full_domain_moment", True
            ),
        ),
        (
            "formal_inverse_promoted_without_certificate",
            lambda data: data["truth_boundary"].__setitem__(
                "formal_axis_based_inverse_verified", True
            ),
        ),
        (
            "axis_regularity_promoted",
            lambda data: data["truth_boundary"].__setitem__(
                "axis_regularity_verified", True
            ),
        ),
        (
            "correction_target_promoted",
            lambda data: data["truth_boundary"].__setitem__(
                "correction_target_admitted", True
            ),
        ),
        (
            "pde_promoted",
            lambda data: data["truth_boundary"].__setitem__(
                "pde_validated", True
            ),
        ),
        (
            "pde_pending_made_delivery_blocker",
            lambda data: data["truth_boundary"].__setitem__(
                "pde_pending_blocks_callable_velocity_delivery", True
            ),
        ),
        (
            "autonomous_sampling_laundered_as_public",
            lambda data: (
                data["classification"]["autonomous_design"].remove(
                    "strict_inner_positive_radius_sampling"
                ),
                data["classification"]["public_source_fact"].append(
                    "strict_inner_positive_radius_sampling"
                ),
            ),
        ),
        (
            "provenance_overlap",
            lambda data: data["classification"]["public_source_fact"].append(
                "manufactured_lower_limit_witness"
            ),
        ),
        (
            "momentum_threshold_relaxed",
            lambda data: data["cr001_snapshot"].__setitem__(
                "pde_residual_max", 0.002
            ),
        ),
        (
            "free_residual_force_enabled",
            lambda data: data["cr001_snapshot"].__setitem__(
                "free_residual_force_allowed", True
            ),
        ),
        (
            "amplitude_collapse_enabled",
            lambda data: data["cr001_snapshot"].__setitem__(
                "amplitude_collapse_allowed", True
            ),
        ),
    ],
)
def test_contract_mutations_fail_closed(
    tmp_path: Path, name: str, mutator
) -> None:
    sandbox = _sandbox(tmp_path)
    _mutate_json(sandbox / CONTRACT_REL, mutator)
    with pytest.raises(AuditError):
        audit(sandbox)


def test_live_cr001_threshold_relaxation_fails_closed(tmp_path: Path) -> None:
    sandbox = _sandbox(tmp_path)

    def relax(data):
        data["validation"]["thresholds"]["pde_residual_max"] = 0.002

    _mutate_json(sandbox / CONSTRAINTS_REL, relax)
    with pytest.raises(AuditError, match="pde_residual_max"):
        audit(sandbox)


def test_live_cr001_free_force_laundering_fails_closed(tmp_path: Path) -> None:
    sandbox = _sandbox(tmp_path)

    def permit(data):
        data["forcing"]["restriction"] = "pointwise free force allowed"

    _mutate_json(sandbox / CONSTRAINTS_REL, permit)
    with pytest.raises(AuditError, match="free-force"):
        audit(sandbox)


def test_helper_lower_limit_mechanics_drift_fails_closed(tmp_path: Path) -> None:
    sandbox = _sandbox(tmp_path)
    path = sandbox / HELPER_REL
    text = path.read_text(encoding="utf-8")
    text = text.replace(
        "out = np.zeros_like(values)",
        "out = np.ones_like(values)",
        1,
    )
    path.write_text(text, encoding="utf-8")
    with pytest.raises(AuditError, match="helper"):
        audit(sandbox)
