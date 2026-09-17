import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_velocity_parameter_perturbation import (
    audit_velocity_parameter_perturbations,
)
from openai_ns_reconstruction.velocity_components import DEFAULT_CANDIDATE


def test_packaged_candidate_perturbations_use_public_velocity_artifacts():
    report = audit_velocity_parameter_perturbations(
        DEFAULT_CANDIDATE,
        delta=1.0e-3,
        coefficient_indices=(0, 58, 116),
        seed=1907,
        held_out_points=32,
        times=(0.35, 0.65),
    )
    assert report["candidate_sha256"]
    assert report["held_out_points"] == 32
    assert len(report["variants"]) == 3
    for row in report["variants"]:
        assert row["base_velocity_rms"] > 0
        assert row["vector_error_rms"] > 0
        assert row["vector_error_max"] >= row["vector_error_rms"]
        assert np.isfinite(row["relative_rms_to_base"])
        assert row["relative_rms_to_base"] > 0
        assert abs(row["signed_delta"]) == pytest.approx(1.0e-3)


def test_same_seed_is_bitwise_reproducible():
    kwargs = dict(
        delta=2.0e-3,
        coefficient_indices=(0,),
        seed=991,
        held_out_points=16,
        times=(0.4, 0.6),
    )
    first = audit_velocity_parameter_perturbations(DEFAULT_CANDIDATE, **kwargs)
    second = audit_velocity_parameter_perturbations(DEFAULT_CANDIDATE, **kwargs)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


@pytest.mark.parametrize(
    "invalid",
    [
        {"delta": 0.0},
        {"coefficient_indices": ()},
        {"coefficient_indices": (0, 0)},
        {"coefficient_indices": (117,)},
        {"held_out_points": 7},
        {"times": (0.2,)},
        {"box": ((0, 0), (-1, 1), (-1, 1))},
    ],
)
def test_invalid_requests_fail_closed(invalid):
    kwargs = dict(
        delta=1.0e-3,
        coefficient_indices=(0,),
        seed=1,
        held_out_points=8,
        times=(0.5,),
    )
    kwargs.update(invalid)
    with pytest.raises((TypeError, ValueError)):
        audit_velocity_parameter_perturbations(DEFAULT_CANDIDATE, **kwargs)
