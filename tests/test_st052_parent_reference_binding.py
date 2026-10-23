from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.st052_parent_reference_binding import (
    DEFAULT_PARENT_REFERENCE,
    TRUTH_BOUNDARY,
    St052ParentBindingError,
    St052ParentReferenceSpec,
    St052ReferenceBoundParent,
    behavioral_binding_id,
    verify_parent_callable,
)


def _reference_callable(points, time):
    points = np.asarray(points, dtype=float)
    refs = np.asarray(DEFAULT_PARENT_REFERENCE.points, dtype=float)
    times = np.asarray(DEFAULT_PARENT_REFERENCE.times, dtype=float)
    velocity = np.asarray(DEFAULT_PARENT_REFERENCE.velocity, dtype=float)
    out = np.empty_like(points)
    for i, point in enumerate(points):
        match = np.where(
            (times == float(time)) & np.all(np.isclose(refs, point, atol=0.0, rtol=0.0), axis=1)
        )[0]
        if len(match) == 1:
            out[i] = velocity[match[0]]
        else:
            out[i] = np.array((0.1, -0.2, 0.3))
    return out


def test_reference_spec_round_trip_and_stable_binding_id(tmp_path):
    path = tmp_path / "parent-reference.json"
    DEFAULT_PARENT_REFERENCE.save(path)
    loaded = St052ParentReferenceSpec.load(path)
    assert loaded == DEFAULT_PARENT_REFERENCE
    assert loaded.sha256() == DEFAULT_PARENT_REFERENCE.sha256()
    assert len(behavioral_binding_id(loaded)) == 64

    obj = json.loads(path.read_text())
    obj["spec"]["velocity"][0][0] += 1.0e-4
    path.write_text(json.dumps(obj))
    with pytest.raises(St052ParentBindingError, match="checksum mismatch"):
        St052ParentReferenceSpec.load(path)


def test_exact_reference_callable_binds_and_receipt_is_fail_closed(tmp_path):
    receipt = verify_parent_callable(_reference_callable)
    assert receipt["passed"] is True
    assert receipt["probe_count"] == 10
    assert receipt["max_abs"] == 0.0
    assert receipt["full_parent_function_equality_proved"] is False
    assert receipt["full_parent_artifact_materialized"] is False

    bound = St052ReferenceBoundParent(_reference_callable)
    value = bound(np.array([[0.123, -0.2, 0.4]]), 0.5)
    np.testing.assert_allclose(value, [[0.1, -0.2, 0.3]], atol=0.0, rtol=0.0)
    path = tmp_path / "binding-receipt.json"
    bound.save_receipt(path)
    saved = json.loads(path.read_text())
    assert saved["behavioral_binding_id"] == bound.binding_id
    assert saved["truth_boundary"]["pde_validated"] is False
    assert saved["truth_boundary"]["full_st052_parent_materialized_here"] is False


def test_perturbed_parent_fails_runtime_behavior_binding():
    def perturbed(points, time):
        out = _reference_callable(points, time).copy()
        out[:, 0] += 2.0e-6
        return out

    with pytest.raises(St052ParentBindingError, match="behavioral reference mismatch"):
        verify_parent_callable(perturbed)
    with pytest.raises(St052ParentBindingError):
        St052ReferenceBoundParent(perturbed)


def test_shape_nonfinite_and_time_guards():
    def wrong_shape(points, time):
        return np.zeros((len(points), 2))

    with pytest.raises(St052ParentBindingError, match="shape"):
        verify_parent_callable(wrong_shape)

    bound = St052ReferenceBoundParent(_reference_callable)
    with pytest.raises(ValueError, match="shape"):
        bound(np.zeros((3, 2)), 0.5)
    with pytest.raises(ValueError, match="time outside"):
        bound(np.zeros((1, 3)), 0.8)


def test_truth_boundary_does_not_promote_reference_binding():
    assert TRUTH_BOUNDARY["runtime_parent_behavior_bound_on_frozen_reference"] is True
    assert TRUTH_BOUNDARY["full_st052_parent_materialized_here"] is False
    assert TRUTH_BOUNDARY["parent_artifact_sha256_assigned"] is False
    assert TRUTH_BOUNDARY["whole_domain_parent_equivalence_proved"] is False
    assert TRUTH_BOUNDARY["complete_child_candidate_materialized"] is False
    assert TRUTH_BOUNDARY["complete_candidate_save_load_ready"] is False
    assert TRUTH_BOUNDARY["visualization_ready"] is False
    assert TRUTH_BOUNDARY["pde_validated"] is False
