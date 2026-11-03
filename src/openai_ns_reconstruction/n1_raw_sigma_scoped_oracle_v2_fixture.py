"""Synthetic checks for :mod:`n1_raw_sigma_scoped_oracle_v1`.

The fixture receives the candidate module from its caller.  It imports no
graph producer, old interval leaf, or evaluator and does not touch a saved
receipt.  The sidecar and lexical environments below are synthetic copies of
the narrow contract only.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SCHEMA = "n1_raw_sigma_scoped_oracle_v2_fixture"


class _Env:
    def __init__(self, value: Any, *, shadow_value: Any = None, missing: bool = False) -> None:
        self.value = value
        self.shadow_value = shadow_value
        self.missing = missing
        self.calls: list[Any] = []

    def lookup(self, node_id: int) -> Any:
        self.calls.append(node_id)
        if self.missing:
            return None
        if node_id == 9000:
            return self.shadow_value
        if node_id == 35848:
            return self.value
        return None


class _Backend:
    def __init__(self, module: Any) -> None:
        self.module = module

    def validate_scalar(self, value: Any) -> bool:
        return isinstance(value, self.module.ScopedSigmaInterval)


def _copy(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _copy(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy(item) for item in value]
    return value


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _make_inputs(module: Any) -> tuple[Any, Any, Any, Any, Any, Any]:
    F = Fraction
    interval_type = module.FractionInterval
    safe = interval_type(F(0), F(1, 2**128))
    bounds = module.AuthenticatedSigmaBounds(
        graph_sha256=module.EXPECTED_GRAPH_SHA256,
        argument_interval=safe,
        sigma_interval=interval_type(F(0), F(1, 2**1900)),
        sigma_s_interval=interval_type(F(0), F(1, 2**1600)),
        authenticated=True,
        source_result_sha256="a" * 64,
        receipt_sha256="b" * 64,
        authentication_basis=module.AUTHENTICATION_BASIS,
        coordinate=module.RAW_S_COORDINATE,
        normalized=False,
        source_identity_verified=False,
        owner_binding=False,
    )
    rows = {node_id: _copy(dict(row)) for node_id, row in module.EXPECTED_GRAPH_ROWS.items()}
    sidecar = {
        "authenticated": True,
        "graph_sha256": module.EXPECTED_GRAPH_SHA256,
        "rows": rows,
    }
    saved_result = {
        "authenticated": True,
        "authentication_basis": module.AUTHENTICATION_BASIS,
        "graph_sha256": module.EXPECTED_GRAPH_SHA256,
        "source_result_sha256": "a" * 64,
        "receipt_sha256": "b" * 64,
        "direct_lower_collar": {
            "authenticated": True,
            "evaluated_argument_interval": {"lower": F(0), "upper": F(1, 2**128)},
            "interval_leaf_input": (F(0), F(1, 2**128)),
            "sigma": {"lower": F(0), "upper": F(1, 2**1900)},
            "d_sigma_ds": (F(0), F(1, 2**1600)),
            "coordinate": module.RAW_S_COORDINATE,
            "normalized": False,
        },
    }
    oracle = module.build_oracle(saved_result, sidecar)
    env = _Env(safe, shadow_value=interval_type(F(0), F(0)))
    backend = _Backend(module)
    return oracle, env, backend, rows, safe, saved_result


def run_fixture(module: Any) -> dict[str, Any]:
    F = Fraction
    names = (
        "api",
        "sigma_leaf",
        "sigma_prime_leaf",
        "real_argument_aliases",
        "interval_alias_conflict_rejected",
        "backend_accepts_frozen_interval",
        "node_swap_rejected",
        "unauthorized_node_rejected",
        "interval_below_zero_rejected",
        "interval_above_safe_rejected",
        "missing_binder_rejected",
        "lexical_shadowing_uses_node_id",
        "malicious_record_rejected",
        "normalized_unit_rejected",
        "metadata_conflict_rejected",
        "immutable_result",
        "saved_result_requires_authentication",
        "graph_hash_rejected",
        "no_normalized_claim",
    )
    checks = {name: False for name in names}
    diagnostics: dict[str, Any] = {}
    required = (
        "FractionInterval",
        "AuthenticatedSigmaBounds",
        "ScopedSigmaInterval",
        "build_oracle",
        "EXPECTED_GRAPH_ROWS",
    )
    checks["api"] = all(hasattr(module, name) for name in required)
    if not checks["api"]:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "diagnostics": {"missing_api": [name for name in required if not hasattr(module, name)]},
            "synthetic_only": True,
            "used_old_fixture": False,
            "used_old_module": False,
        }

    try:
        oracle, env, backend, rows, safe, saved_result = _make_inputs(module)
        sigma = oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", env)
        prime = oracle(module.SIGMA_PRIME_NODE_ID, rows[module.SIGMA_PRIME_NODE_ID], "unknown", env)
        checks["sigma_leaf"] = (
            sigma.lower == Fraction(0)
            and sigma.upper == Fraction(1, 2**1900)
            and sigma.quantity == "sigma"
            and sigma.input_interval == safe
            and sigma.source_identity_verified is False
            and sigma.owner_binding is False
            and sigma.normalized is False
            and sigma.certified is False
            and sigma.conditional is True
        )
        direct = saved_result["direct_lower_collar"]
        checks["real_argument_aliases"] = (
            isinstance(direct["evaluated_argument_interval"], dict)
            and isinstance(direct["interval_leaf_input"], tuple)
            and direct["evaluated_argument_interval"]["lower"] == direct["interval_leaf_input"][0]
            and direct["evaluated_argument_interval"]["upper"] == direct["interval_leaf_input"][1]
        )
        alias_conflict = _copy(saved_result)
        alias_conflict["direct_lower_collar"]["interval_leaf_input"] = (F(0), F(1, 2**127))
        checks["interval_alias_conflict_rejected"] = _rejected(
            lambda: module.build_oracle(
                alias_conflict,
                {"authenticated": True, "graph_sha256": module.EXPECTED_GRAPH_SHA256, "rows": rows},
            )
        )
        checks["sigma_prime_leaf"] = (
            prime.lower == Fraction(0)
            and prime.upper == Fraction(1, 2**1600)
            and prime.quantity == "sigma_s"
            and prime.coordinate == module.RAW_S_COORDINATE
        )
        checks["backend_accepts_frozen_interval"] = backend.validate_scalar(sigma) is True
        checks["node_swap_rejected"] = (
            _rejected(lambda: oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_PRIME_NODE_ID], "unknown", env))
            and _rejected(lambda: oracle(module.SIGMA_PRIME_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", env))
        )
        checks["unauthorized_node_rejected"] = _rejected(
            lambda: oracle(35950, {"operation": "product", "arguments": [1, 2]}, "unknown", env)
        )
        below = _Env(module.FractionInterval(F(-1), Fraction(0)))
        above = _Env(module.FractionInterval(F(0), F(1, 2**128) + F(1)))
        checks["interval_below_zero_rejected"] = _rejected(
            lambda: oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", below)
        )
        checks["interval_above_safe_rejected"] = _rejected(
            lambda: oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", above)
        )
        checks["missing_binder_rejected"] = _rejected(
            lambda: oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", _Env(None, missing=True))
        )
        shadow = _Env(safe, shadow_value=module.FractionInterval(F(1), F(1)))
        shadow_result = oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", shadow)
        checks["lexical_shadowing_uses_node_id"] = shadow_result.input_interval == safe and shadow.calls == [module.ARGUMENT_NODE_ID]
        forged = _copy(rows[module.SIGMA_NODE_ID])
        forged["argument"] = 9000
        checks["malicious_record_rejected"] = _rejected(
            lambda: oracle(module.SIGMA_NODE_ID, forged, "unknown", env)
        )
        normalized_env = _Env({"lower": Fraction(0), "upper": F(1), "unit": "normalized_q"})
        checks["normalized_unit_rejected"] = _rejected(
            lambda: oracle(module.SIGMA_NODE_ID, rows[module.SIGMA_NODE_ID], "unknown", normalized_env)
        )
        conflicting = _copy(rows[module.SIGMA_NODE_ID])
        conflicting["op"] = "product"
        checks["metadata_conflict_rejected"] = _rejected(
            lambda: oracle(module.SIGMA_NODE_ID, conflicting, "unknown", env)
        )
        original_lower = sigma.lower
        try:
            sigma.lower = Fraction(1)  # type: ignore[misc]
            immutable = False
        except Exception:
            immutable = sigma.lower == original_lower
        checks["immutable_result"] = immutable and sigma.interval.lower == Fraction(0)
        checks["saved_result_requires_authentication"] = _rejected(
            lambda: module.build_oracle(
                {
                    "graph_sha256": module.EXPECTED_GRAPH_SHA256,
                    "argument_interval": safe,
                    "sigma_interval": module.FractionInterval(F(0), F(1)),
                    "d_sigma_ds": module.FractionInterval(F(0), F(1)),
                },
                {"authenticated": True, "graph_sha256": module.EXPECTED_GRAPH_SHA256, "rows": rows},
            )
        )
        checks["graph_hash_rejected"] = _rejected(
            lambda: module.build_oracle(
                module.AuthenticatedSigmaBounds(
                    graph_sha256=module.EXPECTED_GRAPH_SHA256,
                    argument_interval=safe,
                    sigma_interval=module.FractionInterval(F(0), F(1)),
                    sigma_s_interval=module.FractionInterval(F(0), F(1)),
                    authenticated=True,
                ),
                {"authenticated": True, "graph_sha256": "c" * 64, "rows": rows},
            )
        )
        checks["no_normalized_claim"] = (
            sigma.full_g5 is False
            and sigma.physical_m15 is False
            and sigma.source_form == module.RAW_S_FORM
        )
        diagnostics = {
            "sigma": {"lower": sigma.lower, "upper": sigma.upper},
            "sigma_prime": {"lower": prime.lower, "upper": prime.upper},
            "env_calls": env.calls,
        }
    except Exception as exc:
        diagnostics = {"exception": type(exc).__name__, "message": str(exc)}

    return {
        "schema": SCHEMA,
        "component": "raw_sigma_scoped_interval_oracle_fixture",
        "passed": all(checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "synthetic_only": True,
        "used_old_fixture": False,
        "used_old_module": False,
        "source_identity_verified": False,
        "owner_binding": False,
        "full_g5": False,
        "physical_m15": False,
    }


run = run_fixture


__all__ = ["SCHEMA", "run", "run_fixture"]
