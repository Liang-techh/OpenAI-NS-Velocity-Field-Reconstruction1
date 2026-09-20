from __future__ import annotations

import importlib
import math
import os
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_agent2_complete_curl_binding import (
    AGENT2_BACKEND_HEAD,
    bind_agent2_complete_curl_backend,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import CycleIdentity


class _StubBackend:
    def __init__(self, *, extra: bool = False, mutate_identity: bool = False) -> None:
        self.extra = extra
        self.mutate_identity = mutate_identity

    def agent3_adapter_kwargs(self, identity: CycleIdentity):
        chosen = (
            CycleIdentity(identity.cycle_id, identity.cycle_index, identity.state_token + "-bad")
            if self.mutate_identity
            else identity
        )
        payload = {
            "identity": chosen,
            "radii": (0.2, 0.3),
            "producer_kind": "agent2-complete-curl-test-backend",
            "provenance": "schema-only unit stub",
            "velocity_evaluator": lambda amplitude, x, y, z, t: (0.0, 0.0, 0.0),
            "velocity_dt_evaluator": lambda amplitude, x, y, z, t: (0.0, 0.0, 0.0),
            "source_agent2_complete_curl_certified": False,
        }
        if self.extra:
            payload["residual"] = 0.0
        return payload


def test_binding_is_fail_closed_on_payload_schema_and_identity() -> None:
    identity = CycleIdentity("binding-schema-test", 0, "s0")
    adapter = bind_agent2_complete_curl_backend(_StubBackend(), identity)
    assert adapter.identity == identity
    assert adapter.source_agent2_complete_curl_certified is False

    with pytest.raises(ValueError, match="schema mismatch"):
        bind_agent2_complete_curl_backend(_StubBackend(extra=True), identity)
    with pytest.raises(ValueError, match="changed the requested cycle identity"):
        bind_agent2_complete_curl_backend(_StubBackend(mutate_identity=True), identity)
    with pytest.raises(TypeError, match="agent3_adapter_kwargs"):
        bind_agent2_complete_curl_backend(object(), identity)


def _load_exact_agent2_backend():
    package_path = os.environ.get("KOKUNO_A2_EXACT_PACKAGE_PATH")
    if not package_path:
        pytest.skip("exact Agent-2 worktree is supplied only by dedicated cross-lineage CI")

    # Agent-3 modules above are already imported from this branch.  Prepend only
    # now so the sibling A2-only backend and its private dependencies resolve from
    # the pinned external worktree without copying any curl implementation here.
    import openai_ns_reconstruction as package

    if package_path not in package.__path__:
        package.__path__.insert(0, package_path)
    module = importlib.import_module(
        "openai_ns_reconstruction.kokuno_typed_delta_a_complete_curl_backend"
    )
    assert module.AGENT3_HANDOFF_PR == 717
    return module.TypedDeltaACompleteCurlBackend


def test_exact_agent2_734_backend_binds_and_replays_without_agent3_curl_copy() -> None:
    Backend = _load_exact_agent2_backend()
    radii = np.linspace(0.24, 1.26, 17)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 6
    delta_a = np.stack(
        (
            0.010 * bump * (1.0 + 0.12 * np.cos(2.0 * math.pi * s)),
            -0.008 * bump * (1.0 - 0.09 * np.sin(2.0 * math.pi * s)),
        ),
        axis=-1,
    )
    delta_a[[0, -1], :] = 0.0
    amplitude = SimpleNamespace(
        delta_a=delta_a,
        base_amplitudes=SimpleNamespace(radii=tuple(float(v) for v in radii)),
    )
    identity = CycleIdentity("a2-a3-exact-backend-binding-v1", 0, "s0")
    backend = Backend(
        radii=tuple(float(v) for v in radii),
        reference_time=0.50,
        provenance=f"exact Agent-2 backend {AGENT2_BACKEND_HEAD}",
    )
    adapter = bind_agent2_complete_curl_backend(backend, identity)

    points = (
        (0.41, 0.07, -0.31, 0.37),
        (0.58, -0.13, 0.17, 0.50),
        (0.76, 0.11, 0.29, 0.63),
    )
    saw_nonzero = False
    for point in points:
        direct_u = np.asarray(backend.velocity_evaluator(amplitude, *point), dtype=float)
        bound_u = np.asarray(adapter.velocity_evaluator(amplitude, *point), dtype=float)
        direct_ut = np.asarray(backend.velocity_dt_evaluator(amplitude, *point), dtype=float)
        bound_ut = np.asarray(adapter.velocity_dt_evaluator(amplitude, *point), dtype=float)
        assert np.array_equal(bound_u, direct_u)
        assert np.array_equal(bound_ut, direct_ut)
        assert np.all(np.isfinite(bound_u))
        assert np.all(np.isfinite(bound_ut))
        saw_nonzero = saw_nonzero or bool(np.linalg.norm(bound_u) > 0.0)
    assert saw_nonzero
    assert adapter.producer_kind == "agent2-complete-curl-typed-delta-a-backend"
    assert adapter.source_agent2_complete_curl_certified is False


def test_binding_truth_boundary_does_not_promote_scientific_state() -> None:
    boundary = truth_boundary()
    assert boundary["agent2_backend_head"] == AGENT2_BACKEND_HEAD
    assert boundary["agent2_complete_curl_reimplemented_by_agent3"] is False
    assert boundary["payload_schema_fail_closed"] is True
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["caller_supplied_defect_allowed"] is False
    assert boundary["caller_supplied_stress_allowed"] is False
    assert boundary["caller_supplied_target_allowed"] is False
    assert boundary["caller_supplied_gain_allowed"] is False
    assert boundary["caller_supplied_scientific_threshold_allowed"] is False
    assert boundary["source_agent2_complete_curl_certified_by_this_binding"] is False
    assert boundary["real_full_candidate_bound"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
