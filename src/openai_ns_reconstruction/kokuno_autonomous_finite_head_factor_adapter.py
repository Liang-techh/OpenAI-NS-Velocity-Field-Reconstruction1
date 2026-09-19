"""Pure Agent-5 adapter for the frozen autonomous finite-head factor.

This module extracts only the repository-autonomous finite-head partition replay
used by Agent-3 PR #586 and later PR #686.  The original #586 module also
imports the oscillatory-candidate-specific Agent-2 graph, which is intentionally
not part of the current Agent-5 integration ancestry.  Keeping the scalar factor
replay here avoids importing that unrelated graph while preserving the exact
profile, frozen inputs, arithmetic, and provenance.

This is repository engineering.  It is not Kokuno's theorem-machine
``missingWeight`` and does not claim paper-exact identity.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any

AUTONOMOUS_FACTOR_ORIGIN_PR = 514
AUTONOMOUS_FACTOR_ORIGIN_HEAD = "50d1b8ca965b35d5e0a7210c22e2302cd30c20e8"
AUTONOMOUS_FACTOR_REPLAY_PR = 586
AUTONOMOUS_FACTOR_REPLAY_HEAD = "ca4da20c90ccd7a5046f79fb2e7965ec05b6086f"
AUTONOMOUS_PROFILE_VERSION = "agent3-cinf-logdyadic-l2-v1"
AUTONOMOUS_PROFILE_FORMULA = (
    "b(s)=exp(1-1/(1-s^2)) for |s|<1 else 0; "
    "raw_m(q)=b(log2(q)+m); mask_m=raw_m/sqrt(sum_j raw_j^2)"
)
AUTONOMOUS_PROFILE_SUPPORT = "|log2(q)+m|<1"
PREPARED_N = 5
BAND = 5
COORDINATE_Q = 1.3
EXPECTED_AUTONOMOUS_MISSING_WEIGHT = 0.28410624176179694


def _autonomous_profile(s: float) -> float:
    value = float(s)
    if not math.isfinite(value):
        raise ValueError("profile coordinate must be finite")
    if abs(value) >= 1.0:
        return 0.0
    return math.exp(1.0 - 1.0 / (1.0 - value * value))


def replay_autonomous_finite_head_factor() -> dict[str, Any]:
    """Replay the exact frozen PR #514/#586 autonomous factor."""
    q_phys = math.ldexp(COORDINATE_Q, -BAND)
    log_q = math.log2(q_phys)
    center = -log_q
    lo = math.floor(center) - 2
    hi = math.ceil(center) + 2
    raw = {m: _autonomous_profile(log_q + m) for m in range(lo, hi + 1)}
    raw = {m: value for m, value in raw.items() if value > 0.0}
    norm_sq = math.fsum(value * value for value in raw.values())
    if not math.isfinite(norm_sq) or norm_sq <= 0.0:
        raise ArithmeticError("autonomous partition failed to cover physical q")
    norm = math.sqrt(norm_sq)
    masks = {m: value / norm for m, value in raw.items()}
    closure = math.fsum(value * value for value in masks.values())
    if abs(closure - 1.0) > 128.0 * math.ulp(1.0):
        raise ArithmeticError("autonomous squared partition lost normalization")
    finite_head = set(range(-1, PREPARED_N))
    weight = math.fsum(value * value for m, value in masks.items() if m in finite_head)
    if abs(weight - EXPECTED_AUTONOMOUS_MISSING_WEIGHT) > 5.0e-15:
        raise RuntimeError("autonomous PR #514/#586 factor replay drifted")
    profile_payload = json.dumps(
        {
            "version": AUTONOMOUS_PROFILE_VERSION,
            "formula": AUTONOMOUS_PROFILE_FORMULA,
            "support": AUTONOMOUS_PROFILE_SUPPORT,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        "origin_pr": AUTONOMOUS_FACTOR_ORIGIN_PR,
        "origin_head": AUTONOMOUS_FACTOR_ORIGIN_HEAD,
        "replay_pr": AUTONOMOUS_FACTOR_REPLAY_PR,
        "replay_head": AUTONOMOUS_FACTOR_REPLAY_HEAD,
        "profile_version": AUTONOMOUS_PROFILE_VERSION,
        "profile_sha256": hashlib.sha256(profile_payload).hexdigest(),
        "profile_formula": AUTONOMOUS_PROFILE_FORMULA,
        "profile_support": AUTONOMOUS_PROFILE_SUPPORT,
        "prepared_N": PREPARED_N,
        "band": BAND,
        "coordinate_q": COORDINATE_Q,
        "physical_q": q_phys,
        "finite_head_indices": tuple(range(-1, PREPARED_N)),
        "active_mask_indices": tuple(sorted(masks)),
        "active_mask_values": {str(m): masks[m] for m in sorted(masks)},
        "squared_partition_sum": closure,
        "autonomous_missing_weight": weight,
        "profile_is_repository_autonomous": True,
        "formal_theorem_machine_bump_identity_claimed": False,
        "formal_missing_weight_equality_claimed": False,
        "theorem_missing_weight_replaced": False,
        "theorem_missing_weight_materialized": False,
        "integration_adapter_only": True,
    }
