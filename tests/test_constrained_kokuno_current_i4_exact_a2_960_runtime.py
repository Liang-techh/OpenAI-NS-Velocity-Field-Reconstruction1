from __future__ import annotations

import inspect
import re

import numpy as np

from openai_ns_reconstruction import kokuno_oscillatory_batch_axis_safety as axis
from openai_ns_reconstruction import kokuno_oscillatory_batch_differentials as diff
from openai_ns_reconstruction import kokuno_public_candidate_velocity as public_candidate
from openai_ns_reconstruction import kokuno_public_carrier_resolved_velocity as carrier
from openai_ns_reconstruction import kokuno_public_oscillatory_vorticity_diagnostic as diagnostic
from openai_ns_reconstruction import kokuno_public_z_pullback_velocity as z_pullback
from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT2_DIFFERENTIAL_FUNCTION,
    AGENT2_DIFFERENTIAL_MODULE,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    _git_blob_sha1,
)


EXPECTED_BLOBS = {
    "differential": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "axis_safety": "599baa190ec742d0e532e5517078534124e68d6b",
    "vorticity_diagnostic": "4ae525bad1e4b83c9dcabc1a97d3931562d099a5",
    "z_pullback": "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653",
    "public_candidate": "a0e1b07981ba23a6bd547df14e8cdc89ad0d89ca",
    "carrier_resolved": "77edc85e4ac7aced997c90c2d93e7552e7d0c39b",
}


def _blob(module: object) -> str:
    source = inspect.getsourcefile(module)
    assert source is not None
    return _git_blob_sha1(source)


def test_exact_a2_960_runtime_closure_blobs_are_unchanged() -> None:
    assert _blob(diff) == EXPECTED_BLOBS["differential"]
    assert _blob(axis) == EXPECTED_BLOBS["axis_safety"]
    assert _blob(diagnostic) == EXPECTED_BLOBS["vorticity_diagnostic"]
    assert _blob(z_pullback) == EXPECTED_BLOBS["z_pullback"]
    assert _blob(public_candidate) == EXPECTED_BLOBS["public_candidate"]
    assert _blob(carrier) == EXPECTED_BLOBS["carrier_resolved"]


def test_current_a3_firewall_accepts_the_exact_a2_960_function_identity() -> None:
    fn = diff.evaluate_oscillatory_batch_differentials
    assert fn.__module__ == AGENT2_DIFFERENTIAL_MODULE
    assert fn.__name__ == AGENT2_DIFFERENTIAL_FUNCTION
    source = inspect.getsourcefile(fn)
    assert source is not None
    assert _git_blob_sha1(source) == AGENT2_DIFFERENTIAL_SOURCE_BLOB
    assert AGENT2_DIFFERENTIAL_SOURCE_BLOB == EXPECTED_BLOBS["differential"]
    semantic = diff.batch_differentials_sha256()
    assert re.fullmatch(r"[0-9a-f]{64}", semantic)
    payload = diff.semantic_payload()
    assert payload["truth_boundary"]["oscillatory_velocity_changed"] is False
    assert payload["truth_boundary"]["correction_velocity_constructed"] is False
    assert payload["truth_boundary"]["complete_ns_residual"] is False
    assert diff.public_contract()["forbidden_inputs_present"] == []


def test_axis_and_support_exterior_remain_exact_zero() -> None:
    t0, t1 = axis._time_interval()
    time = 0.5 * (t0 + t1)
    _, r1, _, z1 = axis._support()
    points = np.asarray(
        [
            (0.0, 0.0, 0.0),
            (2.0 * r1, 0.0, 0.0),
            (0.0, 0.0, 2.0 * max(abs(z1), 1.0)),
        ],
        dtype=float,
    )
    out = diff.evaluate_oscillatory_batch_differentials(points, time)
    assert not np.any(out.interior_mask)
    assert np.array_equal(out.velocity, np.zeros((3, 3)))
    assert np.array_equal(out.velocity_jacobian, np.zeros((3, 3, 3)))
    assert np.array_equal(out.divergence, np.zeros(3))
    assert np.array_equal(out.vorticity, np.zeros((3, 3)))
