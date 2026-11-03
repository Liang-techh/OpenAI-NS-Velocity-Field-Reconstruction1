"""Synthetic checks for :mod:`n1_collar_sigma_actual_domain_v1`.

The fixture supplies a tiny in-memory saved mapping and an injected interval
stub.  It never imports the previous interval leaf, opens the real JSON file,
reads source bytes, computes a live graph, or constructs a value from the
real enormous exponent.  The positive case records exactly two evaluator
calls, first on ``[0, 2**(-128)]`` for direct ``sigma(s)`` and then on
``[0, 1]`` for the active normalized argument; negative tuple cases must be
rejected before the stub is entered.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SCHEMA = "n1_collar_sigma_actual_domain_v1_fixture"

_SYNTHETIC_SC = (0, 15, -200, 4)
_SYNTHETIC_SAFE_UPPER = Fraction(1, 2**128)


def _saved_mapping(
    lower: tuple[int, int, int, int] = _SYNTHETIC_SC,
    upper: tuple[int, int, int, int] | None = None,
) -> dict[str, Any]:
    if upper is None:
        upper = lower
    return {
        "explicit_current_inner_exit_strict_collar": {
            "selected_first_phase_endpoint": {
                "lower_exact_mpf_tuple": list(lower),
                "upper_exact_mpf_tuple": list(upper),
            }
        }
    }


class _SyntheticIntervalEvaluator:
    """A callable stub with a visible call ledger and rational intervals."""

    def __init__(
        self,
        *,
        ordinary_x_derivative: bool = True,
        z_derivative: bool = False,
        source_identity_verified: bool = False,
        owner_binding: bool = False,
        derivative_key: str = "d_sigma_dx",
    ) -> None:
        self.calls: list[tuple[Fraction, Fraction]] = []
        self.ordinary_x_derivative = ordinary_x_derivative
        self.z_derivative = z_derivative
        self.source_identity_verified = source_identity_verified
        self.owner_binding = owner_binding
        self.derivative_key = derivative_key

    def __call__(self, lower: Fraction, upper: Fraction) -> dict[str, Any]:
        self.calls.append((lower, upper))
        return {
            "input_interval": {"lower": lower, "upper": upper},
            "sigma": {"lower": Fraction(0), "upper": Fraction(1, 2)},
            self.derivative_key: {
                "lower": Fraction(0),
                "upper": Fraction(1),
            },
            "ordinary_x_derivative": self.ordinary_x_derivative,
            "z_derivative": self.z_derivative,
            "source_identity_verified": self.source_identity_verified,
            "owner_binding": self.owner_binding,
        }


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def run_fixture(module: Any) -> dict[str, Any]:
    """Run tuple, call-count, cap, and fixed-s ordinary-Z checks."""

    checks: dict[str, bool] = {
        "consumer_api": False,
        "tuple_api": False,
        "tuple_proof": False,
        "exact_lower_upper_match": False,
        "safe_input_once": False,
        "normalized_active_input_once": False,
        "two_distinct_call_objects": False,
        "cap_intersection_is_bound": False,
        "normalized_piecewise_windows": False,
        "normalized_chain_rule_bound": False,
        "fixed_s_sigma_Z_zero": False,
        "fixed_s_over_sc_sigma_Z_zero": False,
        "fixed_s_without_ratio_not_normalized_zero": False,
        "unverified_s_Z_not_asserted": False,
        "lower_upper_mismatch_rejected": False,
        "sign_error_rejected": False,
        "bitcount_error_rejected": False,
        "large_positive_exponent_rejected": False,
        "wrong_derivative_semantics_rejected": False,
        "wrong_derivative_key_rejected": False,
        "source_identity_true_rejected": False,
        "owner_binding_true_rejected": False,
        "source_identity_not_claimed": False,
        "no_physical_claims": False,
        "direct_collar_scope_only": False,
    }
    consume = getattr(module, "evaluate_actual_collar_sigma", None)
    validate = getattr(module, "validate_sc_tuple", None)
    prove = getattr(module, "prove_half_sc_domain", None)
    checks["consumer_api"] = callable(consume)
    checks["tuple_api"] = callable(validate)
    checks["tuple_proof"] = callable(prove)
    selected: dict[str, Any] = {}

    if checks["consumer_api"] and checks["tuple_api"] and checks["tuple_proof"]:
        try:
            evaluator = _SyntheticIntervalEvaluator()
            positive = consume(
                _saved_mapping(),
                evaluator,
                fixed_s_in_z=True,
                fixed_s_over_sc_in_z=True,
            )
            selected = {
                "schema": positive.get("schema"),
                "component": positive.get("component"),
                "sc_tuple": positive.get("sc_tuple"),
                "sigma": positive.get("sigma"),
                "d_sigma_ds": positive.get("d_sigma_ds"),
                "normalized_inner_cutoff": positive.get("normalized_inner_cutoff"),
                "derivative_cap": positive.get("derivative_cap"),
                "interval_leaf_call_count": positive.get("interval_leaf_call_count"),
            }
            checks["exact_lower_upper_match"] = positive.get("sc_tuple") == _SYNTHETIC_SC
            checks["safe_input_once"] = (
                len(evaluator.calls) == 2
                and evaluator.calls[0] == (Fraction(0), _SYNTHETIC_SAFE_UPPER)
            )
            checks["normalized_active_input_once"] = evaluator.calls[1:] == [
                (Fraction(0), Fraction(1))
            ]
            checks["two_distinct_call_objects"] = (
                positive.get("direct_lower_collar") is not None
                and positive.get("normalized_inner_cutoff") is not None
                and positive.get("direct_lower_collar")
                is not positive.get("normalized_inner_cutoff")
                and positive.get("interval_leaf_call_count") == 2
                and positive.get("expected_interval_leaf_call_count") == 2
            )
            cap = getattr(module, "ACTUAL_DERIVATIVE_CAP", None)
            checks["cap_intersection_is_bound"] = (
                positive.get("d_sigma_ds", {}).get("upper") == cap
                and positive.get("derivative_cap", {}).get("is_selected_value") is False
                and positive.get("sigma_value") is None
                and positive.get("d_sigma_ds_value") is None
            )
            normalized = positive.get("normalized_inner_cutoff", {})
            flat = normalized.get("flat_windows", {})
            checks["normalized_piecewise_windows"] = (
                normalized.get("q_domain") == {
                    "lower": Fraction(0),
                    "upper": Fraction(1),
                }
                and normalized.get("active_window", {}).get("q_domain") == {
                    "lower": Fraction(1, 2),
                    "upper": Fraction(3, 4),
                }
                and flat.get("q_le_half", {}).get("sigma")
                == {"lower": Fraction(0), "upper": Fraction(0)}
                and flat.get("q_ge_three_quarters", {}).get("sigma")
                == {"lower": Fraction(1), "upper": Fraction(1)}
            )
            checks["normalized_chain_rule_bound"] = (
                normalized.get("d_sigma_dq", {}).get("lower") == Fraction(0)
                and normalized.get("d_sigma_dq", {}).get("upper") == 4
                and normalized.get("sigma_x_not_Z1") is True
                and normalized.get("argument_derivative_is_not_Z1") is True
            )
            checks["fixed_s_sigma_Z_zero"] = (
                positive.get("fixed_s_in_z") is True
                and positive.get("sigma_Z") == Fraction(0)
                and positive.get("sigma_z") == Fraction(0)
            )
            checks["fixed_s_over_sc_sigma_Z_zero"] = (
                positive.get("fixed_s_over_sc_in_z") is True
                and positive.get("normalized_sigma_Z") == Fraction(0)
                and normalized.get("sigma_Z") == Fraction(0)
            )

            direct_only_evaluator = _SyntheticIntervalEvaluator()
            direct_only = consume(
                _saved_mapping(), direct_only_evaluator, fixed_s_in_z=True
            )
            checks["fixed_s_without_ratio_not_normalized_zero"] = (
                direct_only.get("fixed_s_in_z") is True
                and direct_only.get("sigma_Z") == Fraction(0)
                and direct_only.get("fixed_s_over_sc_in_z") is False
                and direct_only.get("normalized_sigma_Z") is None
                and direct_only.get("normalized_sigma_z") is None
                and direct_only.get("normalized_inner_cutoff", {}).get("sigma_Z") is None
                and direct_only_evaluator.calls == [
                    (Fraction(0), _SYNTHETIC_SAFE_UPPER),
                    (Fraction(0), Fraction(1)),
                ]
            )
            checks["source_identity_not_claimed"] = (
                positive.get("source_identity_verified") is False
                and positive.get("owner_binding") is False
                and positive.get("source_files_authenticated_here") is False
                and positive.get("source_files_and_sha_authentication") == "runner_only"
            )
            checks["no_physical_claims"] = (
                positive.get("physical_m15") is False
                and positive.get("five_moment_closure") is False
                and positive.get("graph_binding") is False
            )
            checks["direct_collar_scope_only"] = (
                positive.get("cutoff_scope") == "direct lower-collar chi sigma(s) only"
                and positive.get("normalized_inner_cutoff_claimed") is True
                and positive.get("normalized_inner_cutoff_live_binding_claimed") is False
                and positive.get("normalized_inner_cutoff_node_binding") is False
            )

            no_z_evaluator = _SyntheticIntervalEvaluator()
            no_z = consume(_saved_mapping(), no_z_evaluator)
            checks["unverified_s_Z_not_asserted"] = (
                no_z.get("sigma_Z") is None
                and no_z.get("sigma_z") is None
                and no_z.get("fixed_s_in_z") is False
                and no_z.get("normalized_sigma_Z") is None
                and no_z.get("normalized_sigma_z") is None
                and no_z_evaluator.calls == [
                    (Fraction(0), _SYNTHETIC_SAFE_UPPER),
                    (Fraction(0), Fraction(1)),
                ]
            )

            mismatch_evaluator = _SyntheticIntervalEvaluator()
            mismatch = _saved_mapping(upper=(0, 14, -200, 4))
            checks["lower_upper_mismatch_rejected"] = _rejected(
                lambda: consume(mismatch, mismatch_evaluator)
            ) and not mismatch_evaluator.calls

            sign_evaluator = _SyntheticIntervalEvaluator()
            checks["sign_error_rejected"] = _rejected(
                lambda: consume(_saved_mapping(lower=(1, 15, -200, 4)), sign_evaluator)
            ) and not sign_evaluator.calls

            bitcount_evaluator = _SyntheticIntervalEvaluator()
            checks["bitcount_error_rejected"] = _rejected(
                lambda: consume(_saved_mapping(lower=(0, 15, -200, 5)), bitcount_evaluator)
            ) and not bitcount_evaluator.calls

            huge_exponent_evaluator = _SyntheticIntervalEvaluator()
            checks["large_positive_exponent_rejected"] = _rejected(
                lambda: consume(
                    _saved_mapping(lower=(0, 15, 10**300, 4)),
                    huge_exponent_evaluator,
                )
            ) and not huge_exponent_evaluator.calls

            wrong_semantics_evaluator = _SyntheticIntervalEvaluator(
                ordinary_x_derivative=False
            )
            checks["wrong_derivative_semantics_rejected"] = _rejected(
                lambda: consume(_saved_mapping(), wrong_semantics_evaluator)
            ) and wrong_semantics_evaluator.calls == [
                (Fraction(0), _SYNTHETIC_SAFE_UPPER)
            ]

            wrong_key_evaluator = _SyntheticIntervalEvaluator(
                derivative_key="d_sigma_ds"
            )
            checks["wrong_derivative_key_rejected"] = _rejected(
                lambda: consume(_saved_mapping(), wrong_key_evaluator)
            ) and wrong_key_evaluator.calls == [
                (Fraction(0), _SYNTHETIC_SAFE_UPPER)
            ]

            source_identity_evaluator = _SyntheticIntervalEvaluator(
                source_identity_verified=True
            )
            checks["source_identity_true_rejected"] = _rejected(
                lambda: consume(_saved_mapping(), source_identity_evaluator)
            ) and source_identity_evaluator.calls == [
                (Fraction(0), _SYNTHETIC_SAFE_UPPER)
            ]

            owner_binding_evaluator = _SyntheticIntervalEvaluator(
                owner_binding=True
            )
            checks["owner_binding_true_rejected"] = _rejected(
                lambda: consume(_saved_mapping(), owner_binding_evaluator)
            ) and owner_binding_evaluator.calls == [
                (Fraction(0), _SYNTHETIC_SAFE_UPPER)
            ]
        except Exception:
            for name in checks:
                if name not in {"consumer_api", "tuple_api", "tuple_proof"}:
                    checks[name] = False

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "synthetic_only": True,
        "used_saved_real_mapping": False,
        "used_previous_interval_leaf": False,
        "source_identity_verified": False,
        "selected": selected,
    }


run = run_fixture


__all__ = ["SCHEMA", "run", "run_fixture"]
