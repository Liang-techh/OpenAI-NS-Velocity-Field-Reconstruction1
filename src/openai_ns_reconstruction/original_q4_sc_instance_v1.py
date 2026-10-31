"""Bind one new, explicit Q4 ``sc`` dyadic instance.

The caller authenticates the archive receipt before this module sees it.  This
adapter reads only the named receipt fields and the selected inner-exit
endpoint.  It does not replay a producer, walk a graph, evaluate a field, or
infer a historical assignment.

The endpoint is a finite-precision exported dyadic.  A singleton exported
interval is exact for this *new selection*, not a claim about an analytic
value, an interval midpoint, or an old source owner.
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCHEMA = "original_q4_sc_instance_v1"


class Missing(ValueError):
    """A required field was not supplied."""


class Mismatch(ValueError):
    """A supplied field failed a strict admission check."""


def _field(value: Any, key: str, label: str) -> Any:
    if not isinstance(value, Mapping) or key not in value:
        raise Missing(label)
    return value[key]


def _strict_int(value: Any) -> bool:
    """A built-in int only; bool is deliberately excluded."""

    return type(value) is int


def _source_family_scope(scope: Mapping[str, Any]) -> dict[str, Any]:
    """Retain only an abstract receipt-scope family marker.

    The real archived sc receipt does not expose a family value in ``scope``.
    Presence of the authenticated scope is retained as an abstraction; no
    family is reconstructed from graph contents or from an unrelated field.
    """

    if not isinstance(scope, Mapping):
        raise Missing("result.archive_result.scope")
    exposed = "source_family" in scope or "family" in scope
    value = scope.get("source_family", scope.get("family"))
    if exposed and (value is None or value is False or value == ""):
        raise Mismatch("receipt scope source_family is empty")
    return {
        "receipt_scope_present": True,
        "source_family_value_exposed": exposed,
        "source_family_value": value if exposed else None,
        "reconstructed": False,
    }


def _graph_id(result: Mapping[str, Any]) -> dict[str, str]:
    """Use the authenticated Q4/LOOP raw input hashes as graph identity."""

    input_hashes = _field(result, "input_hashes", "receipt.result.input_hashes")
    if not isinstance(input_hashes, Mapping):
        raise Missing("receipt.result.input_hashes")
    q4_keys = []
    loop_keys = []
    for key in input_hashes:
        if not isinstance(key, str):
            continue
        lower = key.lower()
        # Receipts may also pin source/fixture/runner files.  Only the two
        # archived raw JSON gzip inputs are eligible graph-identity material.
        if not lower.endswith(".json.gz"):
            continue
        is_q4 = "q4" in lower
        is_loop = "loop" in lower
        if is_q4 and is_loop:
            raise Mismatch("one input hash cannot be both Q4 and LOOP")
        if is_q4:
            q4_keys.append(key)
        if is_loop:
            loop_keys.append(key)
    if len(q4_keys) != 1 or len(loop_keys) != 1:
        raise Missing("receipt.result.input_hashes must expose one Q4 and one LOOP hash")
    q4_hash = input_hashes[q4_keys[0]]
    loop_hash = input_hashes[loop_keys[0]]
    if not isinstance(q4_hash, str) or not q4_hash:
        raise Mismatch("Q4 input hash must be a non-empty string")
    if not isinstance(loop_hash, str) or not loop_hash:
        raise Mismatch("LOOP input hash must be a non-empty string")
    return {
        "q4": q4_hash,
        "loop": loop_hash,
        "q4_input_hash": q4_hash,
        "loop_input_hash": loop_hash,
    }


def _sha_field(evidence: Mapping[str, Any]) -> Any:
    raw = _field(
        evidence,
        "selected_shared_node_sha256",
        "result.archive_result.evidence.selected_shared_node_sha256",
    )
    if isinstance(raw, Mapping):
        raw = _field(
            raw,
            "sc",
            "result.archive_result.evidence.selected_shared_node_sha256.sc",
        )
    if not isinstance(raw, str) or not raw:
        raise Mismatch("selected sc node sha256 must be a non-empty string")
    return raw


def _canonical_positive_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    """Validate a positive canonical ``[sign, mantissa, exponent, bitcount]``.

    The representation invariant is checked without constructing any power of
    two.  Canonical nonzero mpmath mantissas are positive, odd, and have
    ``mantissa.bit_length() == bitcount``.
    """

    if not isinstance(value, (list, tuple)) or len(value) != 4:
        raise Missing(label + " must be a four-element exact_mpf_tuple")
    sign, mantissa, exponent, bitcount = value
    if not all(_strict_int(item) for item in value):
        raise Mismatch(label + " contains a non-int or bool")
    if sign not in (0, 1):
        raise Mismatch(label + ".sign must be 0 or 1")
    if mantissa < 0:
        raise Mismatch(label + ".mantissa must be non-negative")
    if bitcount < 0:
        raise Mismatch(label + ".bitcount must be non-negative")
    if mantissa == 0:
        if (sign, mantissa, exponent, bitcount) != (0, 0, 0, 0):
            raise Mismatch(label + " has a non-canonical zero")
        raise Mismatch(label + " must be strictly positive")
    if sign != 0:
        raise Mismatch(label + " is not strictly positive")
    if bitcount != mantissa.bit_length():
        raise Mismatch(label + ".bitcount is not canonical")
    if mantissa % 2 != 1:
        raise Mismatch(label + ".mantissa is not canonical (not odd)")
    if bitcount <= 0:
        raise Mismatch(label + ".bitcount must be positive")
    return sign, mantissa, exponent, bitcount


def _scope(source_family: Any = None, graph_id: Any = None,
           *, binding_materialized: bool = False) -> dict[str, Any]:
    """The deliberately narrow claim scope used by every returned result."""

    return {
        "newly_chosen_binding": True,
        "existing_archive_assignment_proved": False,
        "corrected_owner_identity_proved": False,
        "old_errors_removed": False,
        "binding_materialized": bool(binding_materialized),
        "source_family": source_family,
        "graph_id": graph_id,
        # This is a policy boundary, not a claim that bind consumed a flat
        # receipt; ``prior_q4_flat_conclusion_used`` remains false below.
        "q4_flat_conclusion_used_only_conditionally": True,
        "prior_q4_flat_result_requires_external_authentication": True,
        "prior_q4_flat_conclusion_used": False,
        "flat_conclusion_proved_here": False,
        "numeric_values_computed": False,
        "numeric_field_evaluation": False,
        "analytic_endpoint_value_claimed": False,
        "analytic_value_equal_interval_midpoint_claimed": False,
        "old_same_source_claimed": False,
        "source_graph_walked": False,
    }


def _dyadic_record(name: str, bits: tuple[int, int, int, int],
                   expression: str) -> dict[str, Any]:
    sign, mantissa, exponent, bitcount = bits
    return {
        "name": name,
        "exact_mpf_tuple": [sign, mantissa, exponent, bitcount],
        "sign": sign,
        "mantissa": mantissa,
        "integer_mantissa": mantissa,
        "binary_exponent": exponent,
        "exponent": exponent,
        "bitcount": bitcount,
        "symbolic_expression": expression,
        "power_of_two_materialized": False,
    }


def _derived_dyadics(sc: tuple[int, int, int, int]) -> dict[str, Any]:
    """Build sc/2 and 3sc/4 using integer mantissa/exponent updates only."""

    _sign, mantissa, exponent, bitcount = sc
    half = (0, mantissa, exponent - 1, bitcount)
    three_mantissa = 3 * mantissa
    three_quarters = (
        0, three_mantissa, exponent - 2, three_mantissa.bit_length()
    )
    _canonical_positive_mpf_tuple(half, "derived sc_over_two")
    _canonical_positive_mpf_tuple(three_quarters, "derived three_sc_over_four")
    half_record = _dyadic_record(
        "sc_over_two", half, "canonical_mantissa * 2**binary_exponent"
    )
    three_quarters_record = _dyadic_record(
        "three_sc_over_four", three_quarters,
        "canonical_mantissa * 2**binary_exponent",
    )
    comparison = {
        "sc_over_two_lt_three_sc_over_four": True,
        "three_sc_over_four_lt_sc": True,
        "proof": (
            "at common exponent e-2, 2*m < 3*m < 4*m for positive m; "
            "only integer mantissas and exponent shifts are used"
        ),
        "common_exponent": exponent - 2,
        "left_common_integer_mantissa": 2 * mantissa,
        "middle_integer_mantissa": three_mantissa,
        "right_common_integer_mantissa": 4 * mantissa,
        "power_of_two_materialized": False,
    }
    return {
        "sc_over_two": half_record,
        "sc_half": half_record,
        "three_sc_over_four": three_quarters_record,
        "sc_three_quarters": three_quarters_record,
        "comparison": comparison,
        "sc_less_than_one_eighth_proof": {
            "bitcount_plus_exponent": bitcount + exponent,
            "criterion": "bitcount + exponent <= -3",
            "passed": bitcount + exponent <= -3,
            "power_of_two_materialized": False,
        },
    }


def _failure(status: str, checks: dict[str, bool], errors: list[str],
             source_family: Any, graph_id: Any) -> dict[str, Any]:
    claim_scope = _scope(source_family, graph_id, binding_materialized=False)
    return {
        "schema": SCHEMA,
        "status": status,
        "passed": False,
        "checks": checks,
        "errors": errors,
        "binding": None,
        "selected_source_parameter": None,
        "derived_dyadics": None,
        "scope": claim_scope,
        "proof_scope": dict(claim_scope),
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
    }


def bind(sc_receipt: Any, endpoint: Any, expected_sc_id: Any,
         expected_sc_hash: Any) -> dict[str, Any]:
    """Return a new explicit Q4 ``sc`` binding, or a fail-closed result.

    ``sc_receipt`` is already authenticated by its caller.  This function
    checks its verified status, source-family scope, graph id, selected node
    id and selected-node digest.  ``endpoint`` is the direct
    ``selected_first_phase_endpoint`` object; a wrapper of that name is also
    accepted when the direct tuple keys are absent.
    """

    checks: dict[str, bool] = {}
    errors: list[str] = []
    source_family: Any = None
    graph_id: Any = None

    def check(label: str, condition: bool) -> None:
        checks[label] = bool(condition)
        if not condition:
            raise Mismatch(label)

    try:
        check("expected_sc_id_is_strict_int", _strict_int(expected_sc_id))
        if not isinstance(expected_sc_hash, str) or not expected_sc_hash:
            raise Mismatch("expected_sc_hash must be a non-empty string")
        checks["expected_sc_hash_is_nonempty_string"] = True

        result = _field(sc_receipt, "result", "receipt.result")
        archive = _field(result, "archive_result", "receipt.result.archive_result")
        check("archive_status_verified",
              _field(archive, "status", "archive_result.status") == "verified")
        # Archived receipts in this repository put ``passed`` on the outer
        # result, not always on archive_result.  A present field is checked;
        # authenticated status=verified is sufficient when it is absent.
        if "passed" in archive:
            check("archive_passed_true", archive["passed"] is True)
        else:
            checks["archive_passed_field_absent_allowed_by_verified_status"] = True
        evidence = _field(archive, "evidence", "archive_result.evidence")
        scope = _field(archive, "scope", "archive_result.scope")
        if not isinstance(evidence, Mapping):
            raise Missing("archive_result.evidence")
        if not isinstance(scope, Mapping):
            raise Missing("archive_result.scope")
        source_family = _source_family_scope(scope)
        checks["receipt_scope_present_without_family_reconstruction"] = True
        archive_checks = archive.get("checks")
        if isinstance(archive_checks, Mapping) and \
                "archive_family_metadata_matches" in archive_checks:
            check("archive_family_metadata_matches",
                  archive_checks["archive_family_metadata_matches"] is True)
        else:
            checks["archive_family_metadata_check_absent_no_family_claim"] = True
        graph_id = _graph_id(result)
        checks["graph_identifier_explicit"] = True

        receipt_sc_id = _field(evidence, "sc", "archive_result.evidence.sc")
        check("receipt_sc_id_is_strict_int", _strict_int(receipt_sc_id))
        check("receipt_sc_id_matches_expected", receipt_sc_id == expected_sc_id)
        receipt_sc_hash = _sha_field(evidence)
        check("selected_sc_hash_matches_expected",
              receipt_sc_hash == expected_sc_hash)

        if not isinstance(endpoint, Mapping):
            raise Missing("selected_first_phase_endpoint")
        endpoint_obj = endpoint
        if ("lower_exact_mpf_tuple" not in endpoint_obj
                and "upper_exact_mpf_tuple" not in endpoint_obj
                and "selected_first_phase_endpoint" in endpoint_obj):
            endpoint_obj = _field(endpoint_obj, "selected_first_phase_endpoint",
                                  "selected_first_phase_endpoint")
        lower_raw = _field(
            endpoint_obj, "lower_exact_mpf_tuple",
            "selected_first_phase_endpoint.lower_exact_mpf_tuple",
        )
        upper_raw = _field(
            endpoint_obj, "upper_exact_mpf_tuple",
            "selected_first_phase_endpoint.upper_exact_mpf_tuple",
        )
        lower = _canonical_positive_mpf_tuple(
            lower_raw, "selected_first_phase_endpoint.lower_exact_mpf_tuple"
        )
        checks["lower_exact_mpf_tuple_canonical_positive"] = True
        upper = _canonical_positive_mpf_tuple(
            upper_raw, "selected_first_phase_endpoint.upper_exact_mpf_tuple"
        )
        checks["upper_exact_mpf_tuple_canonical_positive"] = True
        check("endpoint_is_exact_singleton_dyadic", lower == upper)
        sc = lower
        _sign, mantissa, exponent, bitcount = sc
        # For normalized positive m*2**e with m<2**b, b+e <= -3 proves
        # m*2**e < 2**(-3), with no construction of 2**abs(e).
        check("sc_strictly_below_one_eighth_without_expansion",
              bitcount + exponent <= -3)
        checks["sc_positive_proved_from_canonical_tuple"] = True

        derived = _derived_dyadics(sc)
        checks["derived_sc_over_two_is_canonical"] = True
        checks["derived_three_sc_over_four_is_canonical"] = True
        checks["derived_order_comparison_proved"] = (
            derived["comparison"]["sc_over_two_lt_three_sc_over_four"]
            and derived["comparison"]["three_sc_over_four_lt_sc"]
        )
        checks["derived_power_of_two_not_materialized"] = (
            derived["comparison"]["power_of_two_materialized"] is False
        )

        claim_scope = _scope(source_family, graph_id, binding_materialized=True)
        selected_source = {
            "name": "sc",
            "graph_id": graph_id,
            "node_id": expected_sc_id,
            "sc_node_id": expected_sc_id,
            "exact_mpf_tuple": list(sc),
            "sc_exact_mpf_tuple": list(sc),
            "finite_precision_export_dyadic": True,
            "selected_from_singleton_endpoint": True,
            "analytic_value_claimed": False,
            "interval_midpoint_equality_claimed": False,
            "old_same_source_claimed": False,
        }
        binding = {
            "binding_kind": "new_explicit_q4_sc_parameter_instance",
            "graph_id": graph_id,
            "node_id": expected_sc_id,
            "sc_node_id": expected_sc_id,
            "source_family": source_family,
            "selected_source_parameter": selected_source,
            "sc_exact_mpf_tuple": list(sc),
            "endpoint": {
                "lower_exact_mpf_tuple": list(lower),
                "upper_exact_mpf_tuple": list(upper),
                "is_exact_singleton": True,
                "finite_precision_export_dyadic": True,
                "analytic_value_equal_interval_midpoint_claimed": False,
            },
            "derived_dyadics": derived,
            "scope": dict(claim_scope),
        }
        flat_scope = {
            "used": False,
            "conditional_only": True,
            "requires_external_authenticated_receipt": True,
            "statement": (
                "prior Q4 flat-at-sc/2 conclusion is not authenticated or "
                "used by bind; an outer runner must supply that receipt"
            ),
            "sc_positive_for_new_selection": True,
            "numeric_field_evaluation": False,
            "source_family_rebound_or_old_identity": False,
            "scope": dict(claim_scope),
        }
        return {
            "schema": SCHEMA,
            "status": "verified",
            "passed": True,
            "checks": checks,
            "errors": [],
            "binding": binding,
            "selected_source_parameter": selected_source,
            "derived_dyadics": derived,
            "q4_flat_conclusion": flat_scope,
            "scope": claim_scope,
            "proof_scope": dict(claim_scope),
            "source_code_executed": False,
            "source_graph_walked": False,
            "source_values_evaluated": False,
        }
    except Missing as exc:
        errors.append(str(exc))
        return _failure("unresolved", checks, errors, source_family, graph_id)
    except Mismatch as exc:
        errors.append(str(exc))
        return _failure("mismatch", checks, errors, source_family, graph_id)
    except Exception as exc:  # defensive fail-closed boundary for JSON-like input
        errors.append(type(exc).__name__ + ": " + str(exc))
        return _failure("unresolved", checks, errors, source_family, graph_id)


__all__ = ["SCHEMA", "Missing", "Mismatch", "bind"]
