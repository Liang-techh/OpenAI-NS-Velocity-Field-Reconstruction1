import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_swirl_capacity import (
    audit_supported_swirl_capacity,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _child() -> Eq45SupportedVelocityCandidate:
    return Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.load_json(SEED))


def test_supported_swirl_capacity_calibration():
    report = audit_supported_swirl_capacity(_child())
    # One-shot exact-CI calibration.  This assertion is replaced by fail-closed
    # numerical regressions after harvesting the report from the Actions log.
    assert False, json.dumps(report, sort_keys=True)


def test_supported_swirl_capacity_rejects_bad_steps_and_times():
    child = _child()
    with pytest.raises(ValueError):
        audit_supported_swirl_capacity(child, coefficient_steps=(0.02, 0.04))
    with pytest.raises(ValueError):
        audit_supported_swirl_capacity(child, coefficient_steps=(0.02, 0.02))
    with pytest.raises(ValueError):
        audit_supported_swirl_capacity(child, times=(0.20, 0.50))


def test_supported_swirl_capacity_does_not_mutate_candidate():
    child = _child()
    before = child.sha256
    audit = audit_supported_swirl_capacity(child)
    assert child.sha256 == before
    assert audit["parent_sha256"] == child.parent.sha256
    assert audit["child_sha256"] == child.sha256
    assert audit["parameter_count_added"] == 0
    assert audit["truth_boundary"]["velocity_changed"] is False
    assert audit["truth_boundary"]["new_basis_added"] is False
