import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.eq45_supported_delivery import (
    Eq45SupportedDeliveryField,
    default_field,
    velocity,
)


ROOT = Path(__file__).resolve().parents[1]


def _json(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def test_current_delivery_identity_agrees_across_governed_surfaces():
    contract = _json("configs/velocity_delivery_contract.json")
    status = _json("project_status.json")
    capsule = _json("artifacts/constrained/eq45_supported_delivery_capsule.json")
    primary = contract["primary_deliverable"]

    assert contract["schema_version"] == 2
    assert primary["canonical"] is True
    assert primary["candidate_family"] == status["candidate_family"] == capsule["candidate"]["family"]
    assert primary["candidate_sha256"] == status["candidate_sha256"] == capsule["candidate"]["child_sha256"]
    assert primary["api"] == status["velocity_api"] == capsule["public_interface"]["velocity"]
    assert primary["grid_api"] == status["grid_api"] == capsule["public_interface"]["grid"]
    assert primary["load_api"] == status["candidate_save_load_api"] == capsule["public_interface"]["load_candidate"]


def test_readme_points_to_active_eq45_delivery_and_scopes_legacy_surface():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    contract = _json("configs/velocity_delivery_contract.json")
    legacy = contract["legacy_compatibility_delivery"]

    assert "from openai_ns_reconstruction.eq45_supported_delivery import velocity" in readme
    assert "python -m openai_ns_reconstruction.eq45_export_bundle" in readme
    assert "### Legacy compatibility surface" in readme
    assert legacy["canonical"] is False
    assert legacy["candidate_family"] == "coupled_velocity_v1"
    assert legacy["api"] in readme
    assert legacy["cli"] in readme
    assert "not** the current canonical Eq45 candidate" in readme


def test_velocity_export_readiness_cannot_promote_independent_science_states():
    contract = _json("configs/velocity_delivery_contract.json")
    status = _json("project_status.json")
    capsule = _json("artifacts/constrained/eq45_supported_delivery_capsule.json")

    assert contract["claim_status"]["velocity_export_ready"] is True
    assert status["states"]["velocity_export_ready"] is True
    assert capsule["states"]["velocity_export_ready"] is True

    independent_false_states = (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for state in independent_false_states:
        assert contract["claim_status"][state] is False
        assert status["states"][state] is False
        assert capsule["states"][state] is False


def test_canonical_candidate_save_load_replays_same_velocity(tmp_path):
    contract = _json("configs/velocity_delivery_contract.json")
    expected_sha = contract["primary_deliverable"]["candidate_sha256"]

    field = default_field()
    assert field.sha256 == expected_sha

    point = (0.1, 0.0, 0.1, 0.5)
    direct = velocity(*point)
    wrapped = field.velocity(*point)
    assert direct.shape == (3,)
    assert np.isfinite(direct).all()
    np.testing.assert_array_equal(direct, wrapped)

    path = tmp_path / "candidate.json"
    field.save_candidate(path)
    reloaded = Eq45SupportedDeliveryField.load_candidate(path)
    assert reloaded.sha256 == expected_sha
    np.testing.assert_array_equal(reloaded.velocity(*point), direct)


def test_api_governance_does_not_mutate_cr001_contract():
    contract = _json("configs/velocity_delivery_contract.json")
    constraints = _json("configs/constraints.json")
    guard = contract["cr001_nonmutation"]
    thresholds = constraints["validation"]["thresholds"]

    assert guard == {
        "constraints_path": "configs/constraints.json",
        "scientific_parameters_changed": False,
        "thresholds_changed": False,
        "forcing_contract_changed": False,
        "validation_sample_changed": False,
    }
    assert constraints["nu"] == 0.01
    assert constraints["validation"]["seed"] == 914027
    assert constraints["validation"]["held_out_points"] == 4096
    assert thresholds["pde_residual_max"] == 1e-3
    assert thresholds["pde_residual_L2"] == 1e-3
