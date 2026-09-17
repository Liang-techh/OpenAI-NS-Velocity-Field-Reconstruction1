from copy import deepcopy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_bipolar_capped_delivery import (
    EXPECTED_CANDIDATE_SHA256,
    load_field,
    load_recipe,
    validate_recipe,
    write_delivery_bundle,
)
from openai_ns_reconstruction.eq45_supported_delivery import Eq45SupportedDeliveryField


def test_capped_delivery_recipe_binds_candidate_and_failure_boundary():
    recipe = load_recipe()
    field = load_field()

    assert field.sha256 == EXPECTED_CANDIDATE_SHA256
    assert recipe["status"]["velocity_export_ready"] is True
    assert recipe["status"]["visualization_ready"] is False
    assert recipe["status"]["pde_validated"] is False
    receipt = recipe["fresh_validation_receipt"]
    assert receipt["fresh_sample_divergence_thresholds_passed"] is True
    assert receipt["formal_registered_divergence_gate_assessed"] is False
    assert receipt["necessary_axisymmetric_pressure_PDE_condition_passed"] is False

    points = np.array(
        [[0.1, 0.0, 0.1], [0.4, -0.3, 0.2], [-0.7, 0.2, -0.5]], dtype=float
    )
    values = field.at_points(points, 0.5)
    assert values.shape == (3, 3)
    assert np.all(np.isfinite(values))
    assert np.linalg.norm(values) > 0.0


def test_capped_delivery_bundle_roundtrips_public_velocity_and_grid(tmp_path):
    paths = write_delivery_bundle(tmp_path)
    restored = Eq45SupportedDeliveryField.load_candidate(paths["candidate"])
    capsule = json.loads(paths["capsule"].read_text())

    assert restored.sha256 == EXPECTED_CANDIDATE_SHA256
    assert capsule["candidate_sha256"] == EXPECTED_CANDIDATE_SHA256
    assert capsule["status"]["velocity_export_ready"] is True
    assert capsule["status"]["visualization_ready"] is False
    assert capsule["status"]["pde_validated"] is False

    field = load_field()
    points = np.array(
        [[0.2, 0.1, 0.3], [-0.5, 0.25, -0.4], [0.9, -0.1, 0.6]], dtype=float
    )
    for time in (0.25, 0.5, 0.75):
        np.testing.assert_array_equal(
            restored.at_points(points, time), field.at_points(points, time)
        )


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("optimizer", "energy_defect_cap_role"), "acceptance_threshold"),
        (("optimizer", "acceptance_threshold_changed"), True),
        (("fresh_validation_receipt", "formal_registered_divergence_gate_assessed"), True),
        (
            ("fresh_validation_receipt", "necessary_axisymmetric_pressure_PDE_condition_passed"),
            True,
        ),
        (("status", "pde_validated"), True),
        (("status", "visualization_ready"), True),
        (("sibling_evidence", "spatial_redistribution_status"), "consumed"),
    ],
)
def test_capped_delivery_recipe_fails_closed_on_claim_laundering(path, value):
    recipe = deepcopy(load_recipe())
    cursor = recipe
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(ValueError):
        validate_recipe(recipe)
