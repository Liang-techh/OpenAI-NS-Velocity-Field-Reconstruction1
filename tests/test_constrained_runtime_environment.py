from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_runtime_environment import (
    PROJECT_DISTRIBUTION,
    build_runtime_environment_receipt,
    validate_runtime_environment_receipt,
)


def test_checked_in_velocity_has_resolved_runtime_receipt():
    report = build_runtime_environment_receipt()

    assert report["candidate"]["family"] == "coupled_velocity_v1"
    assert len(report["candidate"]["sha256"]) == 64
    assert report["candidate"]["public_evaluator"].endswith("VelocityField.at_points")
    assert report["runtime"]["python_version"]
    required = report["dependencies"]["required_resolved_versions"]
    assert required[PROJECT_DISTRIBUTION] == "0.2.0"
    assert required["numpy"]
    assert required["scipy"]
    assert report["delivery_state"] == {
        "velocity_export_ready": "not_assessed_here",
        "visualization_ready": "not_assessed_here",
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    assert len(report["receipt_sha256"]) == 64


def test_receipt_is_deterministic_within_one_runtime():
    first = build_runtime_environment_receipt()
    second = build_runtime_environment_receipt()
    assert first == second


def test_missing_required_distribution_fails_closed():
    with pytest.raises(ValueError, match="required distribution is not installed"):
        build_runtime_environment_receipt(
            required_distributions=(
                PROJECT_DISTRIBUTION,
                "numpy",
                "scipy",
                "definitely-not-an-installed-distribution-cr011",
            )
        )


def test_claim_promotion_and_version_tampering_are_rejected():
    report = build_runtime_environment_receipt()

    promoted = deepcopy(report)
    promoted["delivery_state"]["pde_validated"] = True
    with pytest.raises(ValueError, match="cannot promote"):
        validate_runtime_environment_receipt(promoted)

    tampered = deepcopy(report)
    tampered["dependencies"]["required_resolved_versions"]["numpy"] = "0.0.0-tampered"
    with pytest.raises(ValueError, match="digest mismatch"):
        validate_runtime_environment_receipt(tampered)
