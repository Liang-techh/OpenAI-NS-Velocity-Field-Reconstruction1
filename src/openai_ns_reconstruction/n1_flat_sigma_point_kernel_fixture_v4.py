"""Synthetic-only fixture for :mod:`n1_flat_sigma_point_kernel_v4`.

The fixture builds a fictional v4 receipt in memory and receives the leaf
module from its caller.  It does not read an archive, import an old source
module, consume a saved numerical result, or recompute an AST.  Its numeric
points, precisions, and synthetic source metadata are new to this revision;
the receipt schema and graph node IDs are the fixed real v1 input contract.
"""

from __future__ import annotations

from copy import deepcopy
from decimal import Decimal, localcontext
from typing import Any


__all__ = ["SCHEMA", "run_fixture", "run"]


SCHEMA = "n1_flat_sigma_point_kernel_fixture_v4"

_SIGMA_MODULE = "synthetic_v4_outer_receipt.py"
_SIGMA_FUNCTION = "stable_sigma_v4"
_SIGMA_MODULE_SHA = "e" * 64
_SIGMA_AST_SHA = "f" * 64
_PRIME_MODULE = "synthetic_v4_derivative_receipt.py"
_PRIME_FUNCTION = "sigma_jets_v4"
_PRIME_MODULE_SHA = "g" * 64
_PRIME_AST_SHA = "h" * 64
_SIGMA_DEFINITION = (
    "0 for x<=0; 1 for x>=1; "
    "exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside"
)
_PRIME_ORDINARY_DERIVATIVE = "same original sigma Taylor coefficient[1]"
_ARGUMENT_VARIABLE_NAME = "original_lower_collar_s"

# The receipt's exact graph IDs remain fixed across execution revisions.
_ARGUMENT_NODE_ID = 35848
_SIGMA_NODE_ID = 35949
_SIGMA_PRIME_NODE_ID = 35953
_UNRELATED_NODE_ID = 48
_UNRELATED_OPERATION = "exact_rational"
_PRIOR_SCHEMA = "n1_raw_g5_dag_program_v1"


def _synthetic_program(
    module: Any,
    *,
    sigma_argument: int = _ARGUMENT_NODE_ID,
    sigma_definition: str = _SIGMA_DEFINITION,
    prime_ordinary_derivative: str = _PRIME_ORDINARY_DERIVATIVE,
    argument_name: str = _ARGUMENT_VARIABLE_NAME,
    sigma_status: str = "opaque",
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build a fictional v4 receipt and independent caller digests."""

    # Keep the synthetic receipt independent of leaf constants.  These are
    # literal IDs/schema from the real prior receipt and intentionally remain
    # fixed even if a future leaf accidentally changes its declarations.
    sigma_id = _SIGMA_NODE_ID
    prime_id = _SIGMA_PRIME_NODE_ID
    argument_id = _ARGUMENT_NODE_ID
    prior_schema = _PRIOR_SCHEMA
    sigma_source = {
        "module": _SIGMA_MODULE,
        "function": _SIGMA_FUNCTION,
        "module_sha256": _SIGMA_MODULE_SHA,
        "AST_sha256": _SIGMA_AST_SHA,
    }
    prime_source = {
        "module": _PRIME_MODULE,
        "function": _PRIME_FUNCTION,
        "module_sha256": _PRIME_MODULE_SHA,
        "AST_sha256": _PRIME_AST_SHA,
    }
    records = {
        argument_id: {
            "operation": "bound_variable",
            "name": argument_name,
        },
        sigma_id: {
            "operation": "exact_original_source_flat_sigma",
            "argument": sigma_argument,
            "definition": sigma_definition,
            "exact_function_not_its_interval_cap": True,
            "source_binding": sigma_source,
        },
        prime_id: {
            "operation": "analytic_unary",
            "function": "original_flat_sigma_prime",
            "argument": sigma_argument,
            "ordinary_derivative": prime_ordinary_derivative,
            "exact_flat_source_derivative_not_a_bound": True,
            "derivative_source": prime_source,
        },
        _UNRELATED_NODE_ID: {
            "operation": _UNRELATED_OPERATION,
            "argument": argument_id,
        },
    }
    # Place the unrelated supported wrapper first so the target-only filter is
    # exercised before either target opaque check.
    ordered = (_UNRELATED_NODE_ID, argument_id, sigma_id, prime_id)
    instructions = [
        {
            "node_id": node,
            "record": dict(records[node]),
            "original_record": dict(records[node]),
            "status": (
                sigma_status
                if node == sigma_id
                else "supported"
                if node == _UNRELATED_NODE_ID
                else "opaque"
            ),
            "op": "supported" if node == _UNRELATED_NODE_ID else "opaque",
            "opaque_operation": records[node]["operation"],
        }
        for node in ordered
    ]
    program = {
        "schema": prior_schema,
        "component": "n1_raw_g5_dag_program_v1",
        "status": "compiled",
        "compile_status": "compiled",
        "receipt_version": 4,
        "graph_hash": "synthetic-v4-receipt-m8q5",
        "conditional": True,
        "certified": False,
        "instructions": instructions,
        "nodes": {node: dict(row) for node, row in zip(ordered, instructions)},
        "records": {node: dict(records[node]) for node in ordered},
    }
    digests = {
        "sigma": dict(sigma_source),
        "sigma_prime": dict(prime_source),
    }
    return program, digests


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _independent_formula(x: Decimal, precision: int) -> tuple[Decimal, Decimal]:
    """Direct Decimal formula at a fresh moderate point for a fixture check."""

    with localcontext() as context:
        context.prec = precision
        one = Decimal(1)
        a = (-one / (x * x)).exp()
        complement = one - x
        b = (-one / (complement * complement)).exp()
        sigma = a / (a + b)
        derivative = sigma * (one - sigma) * (
            Decimal(2) / (x * x * x)
            + Decimal(2) / (complement * complement * complement)
        )
        return +sigma, +derivative


def _copy_program(program: dict[str, Any]) -> dict[str, Any]:
    """Copy wrapper layers so one negative case cannot mutate another."""
    return deepcopy(program)


def _ulp(value: Decimal, precision: int) -> Decimal:
    if value.is_zero():
        return Decimal(1).scaleb(1 - precision)
    return Decimal(1).scaleb(value.copy_abs().adjusted() - precision + 1)


def _close_decimal(
    actual: Decimal,
    expected: Decimal,
    *,
    precision: int,
    absolute_tolerance: Decimal,
    ulps: int = 4,
) -> bool:
    with localcontext() as context:
        context.prec = precision + 24
        difference = abs(actual - expected)
        ulp_bound = Decimal(ulps) * max(
            _ulp(actual, precision), _ulp(expected, precision)
        )
        return difference <= max(absolute_tolerance, ulp_bound)


def _all_false(checks: dict[str, bool], names: tuple[str, ...]) -> None:
    for name in names:
        checks[name] = False


def run_fixture(module: Any) -> dict[str, Any]:
    """Exercise new v4 point and binding corners using synthetic data only."""

    checks: dict[str, bool] = {}
    bind_entry = getattr(module, "bind_flat_sigma_point_kernel", None)
    evaluator = getattr(module, "evaluate_point", None)
    checks["binding_api"] = callable(bind_entry)
    checks["point_api"] = callable(evaluator)
    selected: dict[str, Any] = {}
    names = (
        "bound",
        "schema",
        "component",
        "prior_schema",
        "same_argument",
        "source_digest_provenance_only",
        "semantic_source_records",
        "old_program_not_claimed_executable",
        "point_only_metadata",
        "endpoints",
        "endpoint_metadata",
        "reflection_symmetry",
        "interior_signs",
        "ordinary_derivative_formula",
        "bound_callables_execute",
        "extreme_precision_fails_closed",
        "nonfinite_rejected",
        "non_decimal_rejected",
        "wrong_source_digest_rejected",
        "wrong_source_module_rejected",
        "prior_missing_status_rejected",
        "metadata_alias_conflict_rejected",
        "record_metadata_alias_conflict_rejected",
        "wrong_argument_rejected",
        "wrong_definition_rejected",
        "wrong_prime_semantics_rejected",
        "wrong_argument_name_rejected",
        "supported_operation_rejected",
        "single_wrapper_supported_rejected",
        "record_instruction_conflict_rejected",
        "unrelated_supported_allowed",
        "target_supported_operation_rejected",
        "record_int_string_conflict_rejected",
        "wrapper_op_marker_rejected",
        "prime_source_alias_conflict_rejected",
    )

    if checks["binding_api"] and checks["point_api"]:
        program, digests = _synthetic_program(module)
        try:
            bound = bind_entry(program, digests)
            selected = {
                "schema": bound.get("schema"),
                "component": bound.get("component"),
                "status": bound.get("status"),
                "bound_argument_node_id": bound.get("bound_argument_node_id"),
                "source_identity_verified": bound.get("source_identity_verified"),
                "caller_supplied_source_records_match": bound.get(
                    "caller_supplied_source_records_match"
                ),
                "ast_assertion_status": bound.get("ast_assertion_status"),
                "prior_program_executable": bound.get("prior_program_executable"),
            }
            # These checks intentionally use the independent fixed receipt
            # literals, not declarations imported from the leaf under test.
            sigma_id = _SIGMA_NODE_ID
            prime_id = _SIGMA_PRIME_NODE_ID
            argument_id = _ARGUMENT_NODE_ID
            checks["bound"] = bound.get("status") == "bound"
            unrelated_rows = [
                row
                for row in program["instructions"]
                if row.get("node_id") == _UNRELATED_NODE_ID
            ]
            checks["unrelated_supported_allowed"] = (
                checks["bound"]
                and len(unrelated_rows) == 1
                and unrelated_rows[0].get("status") == "supported"
                and unrelated_rows[0].get("op") == "supported"
                and unrelated_rows[0].get("opaque_operation")
                == _UNRELATED_OPERATION
            )
            checks["schema"] = bound.get("schema") == getattr(
                module, "SCHEMA", "n1_flat_sigma_point_kernel_v4"
            )
            checks["component"] = bound.get("component") == getattr(
                module, "COMPONENT", "n1_flat_sigma_point_kernel_v4"
            )
            checks["prior_schema"] = bound.get("prior_schema") == _PRIOR_SCHEMA
            checks["same_argument"] = (
                bound["bound_nodes"]["sigma"]["argument"] == argument_id
                and bound["bound_nodes"]["sigma_prime"]["argument"] == argument_id
                and bound.get("bound_argument_node_id") == argument_id
                and bound["bound_nodes"]["sigma"]["node_id"] == sigma_id
                and bound["bound_nodes"]["sigma_prime"]["node_id"] == prime_id
            )
            checks["source_digest_provenance_only"] = (
                bound.get("source_identity_verified") is False
                and bound.get("caller_supplied_source_records_match") is True
                and bound.get("ast_assertion_status") == "saved_assertion_not_recomputed"
                and bound.get("claims", {}).get("ast_recomputed") is False
                and bound["source_identities"]["sigma"]["module"] == _SIGMA_MODULE
                and bound["source_identities"]["sigma_prime"]["module"] == _PRIME_MODULE
            )
            sigma_record = bound["bound_nodes"]["sigma"]["original_record"]
            prime_record = bound["bound_nodes"]["sigma_prime"]["original_record"]
            checks["semantic_source_records"] = (
                sigma_record["definition"] == _SIGMA_DEFINITION
                and sigma_record["exact_function_not_its_interval_cap"] is True
                and prime_record["ordinary_derivative"] == _PRIME_ORDINARY_DERIVATIVE
                and prime_record["exact_flat_source_derivative_not_a_bound"] is True
                and sigma_record["argument"] == argument_id
                and prime_record["argument"] == argument_id
                and bound["bound_nodes"]["sigma"]["source_binding"]["AST_sha256"]
                == _SIGMA_AST_SHA
            )
            checks["old_program_not_claimed_executable"] = (
                bound.get("prior_program_executable") is False
                and bound.get("prior_program_execution_supported") is False
                and bound.get("claims", {}).get("prior_program_executable") is False
            )
            checks["point_only_metadata"] = (
                bound.get("certified") is False
                and bound.get("point_evaluation_only") is True
                and bound.get("interval") is False
                and bound.get("full_Navier_Stokes") is False
            )

            endpoint_cases = (
                (Decimal("-4.625"), Decimal(0), Decimal(0)),
                (Decimal("4.625"), Decimal(1), Decimal(0)),
            )
            endpoint_results = [evaluator(x, 97) for x, _, _ in endpoint_cases]
            checks["endpoints"] = all(
                result["value"] == expected
                and result["derivative"] == expected_prime
                and isinstance(result["value"], Decimal)
                and isinstance(result["derivative"], Decimal)
                for result, (_, expected, expected_prime) in zip(
                    endpoint_results, endpoint_cases
                )
            )
            checks["endpoint_metadata"] = all(
                result.get("precision") == 97
                and result.get("certified") is False
                and result.get("point_evaluation_only") is True
                and result.get("interval") is False
                for result in endpoint_results
            )

            reflection_points = (
                Decimal("0.1836"),
                Decimal("0.3467"),
                Decimal("0.8275"),
            )
            reflection_precision = 163
            reflection_tolerance = Decimal("1e-145")
            reflection_checks = []
            for point in reflection_points:
                with localcontext() as context:
                    context.prec = reflection_precision + 8
                    reflected = +(Decimal(1) - point)
                left = evaluator(point, reflection_precision)
                right = evaluator(reflected, reflection_precision)
                with localcontext() as context:
                    context.prec = reflection_precision
                    expected_value = +(Decimal(1) - left["value"])
                reflection_checks.append(
                    _close_decimal(
                        right["value"],
                        expected_value,
                        precision=reflection_precision,
                        absolute_tolerance=reflection_tolerance,
                    )
                    and _close_decimal(
                        right["derivative"],
                        left["derivative"],
                        precision=reflection_precision,
                        absolute_tolerance=reflection_tolerance,
                    )
                    and left.get("certified") is False
                    and right.get("certified") is False
                    and left.get("point_evaluation_only") is True
                    and right.get("point_evaluation_only") is True
                )
            checks["reflection_symmetry"] = all(reflection_checks)

            sign_points = (
                Decimal("0.1836"),
                Decimal("0.3467"),
                Decimal("0.5478"),
                Decimal("0.8275"),
            )
            signs = [evaluator(point, 157) for point in sign_points]
            checks["interior_signs"] = all(
                Decimal(0) < item["value"] < Decimal(1)
                and item["derivative"] > Decimal(0)
                for item in signs
            )

            ordinary_point = Decimal("0.4126")
            ordinary = evaluator(ordinary_point, 167)
            expected_sigma, expected_prime = _independent_formula(ordinary_point, 167)
            with localcontext() as context:
                context.prec = 152
                tolerance = Decimal("1e-140")
                checks["ordinary_derivative_formula"] = (
                    abs(ordinary["value"] - expected_sigma) < tolerance
                    and abs(ordinary["derivative"] - expected_prime) < tolerance
                )

            checks["bound_callables_execute"] = (
                bound["sigma"](ordinary_point, 167) == ordinary["value"]
                and bound["sigma_prime"](ordinary_point, 167) == ordinary["derivative"]
                and bound(ordinary_point, 167)["value"] == ordinary["value"]
            )
            checks["extreme_precision_fails_closed"] = _rejected(
                lambda: evaluator(Decimal("1e-69"), 61)
            )
            checks["nonfinite_rejected"] = _rejected(
                lambda: evaluator(Decimal("NaN"), 97)
            ) and _rejected(lambda: evaluator(Decimal("Infinity"), 97))
            checks["non_decimal_rejected"] = _rejected(lambda: evaluator(3, 97))

            wrong_digest = {
                "sigma": dict(digests["sigma"]),
                "sigma_prime": dict(digests["sigma_prime"]),
            }
            wrong_digest["sigma"]["AST_sha256"] = "i" * 64
            checks["wrong_source_digest_rejected"] = _rejected(
                lambda: bind_entry(program, wrong_digest)
            )
            wrong_module = {
                "sigma": dict(digests["sigma"]),
                "sigma_prime": dict(digests["sigma_prime"]),
            }
            wrong_module["sigma"]["module_sha256"] = "j" * 64
            checks["wrong_source_module_rejected"] = _rejected(
                lambda: bind_entry(program, wrong_module)
            )

            missing_status_program = dict(program)
            missing_status_program.pop("status", None)
            missing_status_program.pop("compile_status", None)
            checks["prior_missing_status_rejected"] = _rejected(
                lambda: bind_entry(missing_status_program, digests)
            )

            alias_conflict_sigma = dict(digests["sigma"])
            alias_conflict_sigma["source_module"] = "synthetic_v4_conflicting_module.py"
            alias_conflict_digests = {
                "sigma": alias_conflict_sigma,
                "sigma_prime": dict(digests["sigma_prime"]),
            }
            checks["metadata_alias_conflict_rejected"] = _rejected(
                lambda: bind_entry(program, alias_conflict_digests)
            )

            # Exercise the same fail-closed rule against a preserved receipt
            # record, using a records-only shape so no earlier copy masks the
            # conflicting alias.
            record_alias_program = {
                "schema": _PRIOR_SCHEMA,
                "status": "compiled",
                "compile_status": "compiled",
                "records": {
                    node: dict(row) for node, row in program["records"].items()
                },
            }
            conflict_sigma_record = dict(
                record_alias_program["records"][_SIGMA_NODE_ID]
            )
            conflict_sigma_source = dict(conflict_sigma_record["source_binding"])
            conflict_sigma_source["source_module"] = (
                "synthetic_v4_conflicting_record_module.py"
            )
            conflict_sigma_record["source_binding"] = conflict_sigma_source
            record_alias_program["records"][_SIGMA_NODE_ID] = conflict_sigma_record
            checks["record_metadata_alias_conflict_rejected"] = _rejected(
                lambda: bind_entry(record_alias_program, digests)
            )

            bad_program, bad_digests = _synthetic_program(
                module, sigma_argument=argument_id - 11
            )
            checks["wrong_argument_rejected"] = _rejected(
                lambda: bind_entry(bad_program, bad_digests)
            )
            bad_definition_program, bad_definition_digests = _synthetic_program(
                module, sigma_definition="synthetic v4 interval cap"
            )
            checks["wrong_definition_rejected"] = _rejected(
                lambda: bind_entry(bad_definition_program, bad_definition_digests)
            )
            bad_prime_program, bad_prime_digests = _synthetic_program(
                module, prime_ordinary_derivative="synthetic v4 derivative bound"
            )
            checks["wrong_prime_semantics_rejected"] = _rejected(
                lambda: bind_entry(bad_prime_program, bad_prime_digests)
            )
            bad_variable_program, bad_variable_digests = _synthetic_program(
                module, argument_name="synthetic_v4_s"
            )
            checks["wrong_argument_name_rejected"] = _rejected(
                lambda: bind_entry(bad_variable_program, bad_variable_digests)
            )
            supported_program, supported_digests = _synthetic_program(
                module, sigma_status="supported"
            )
            checks["supported_operation_rejected"] = _rejected(
                lambda: bind_entry(supported_program, supported_digests)
            )
            checks["target_supported_operation_rejected"] = (
                checks["supported_operation_rejected"]
            )

            # Only one instruction wrapper is marked supported; nodes and the
            # bare records remain opaque/unchanged.  This must still fail.
            single_supported_program = _copy_program(program)
            for row in single_supported_program["instructions"]:
                if row.get("node_id") == _SIGMA_NODE_ID:
                    row["status"] = "supported"
                    break
            checks["single_wrapper_supported_rejected"] = _rejected(
                lambda: bind_entry(single_supported_program, digests)
            )

            # A bare record that disagrees with an instruction record must not
            # be hidden by wrapper selection.
            record_instruction_conflict = _copy_program(program)
            record_instruction_conflict["records"][_SIGMA_NODE_ID]["definition"] = (
                "synthetic v4 raw-record/instruction conflict"
            )
            checks["record_instruction_conflict_rejected"] = _rejected(
                lambda: bind_entry(record_instruction_conflict, digests)
            )

            # Both integer and string aliases for one node must be checked;
            # the string alias is deliberately inconsistent with the integer.
            record_int_string_conflict = _copy_program(program)
            string_sigma_record = dict(
                record_int_string_conflict["records"][_SIGMA_NODE_ID]
            )
            string_sigma_record["definition"] = (
                "synthetic v4 integer-string record alias conflict"
            )
            record_int_string_conflict["records"][str(_SIGMA_NODE_ID)] = (
                string_sigma_record
            )
            checks["record_int_string_conflict_rejected"] = _rejected(
                lambda: bind_entry(record_int_string_conflict, digests)
            )

            # An opaque wrapper cannot carry an op marker claiming support.
            wrapper_op_conflict = _copy_program(program)
            for row in wrapper_op_conflict["instructions"]:
                if row.get("node_id") == _SIGMA_NODE_ID:
                    row["status"] = "opaque"
                    row["op"] = "supported"
                    break
            checks["wrapper_op_marker_rejected"] = _rejected(
                lambda: bind_entry(wrapper_op_conflict, digests)
            )

            # sigma_prime cannot silently prefer derivative_source when a
            # conflicting source_binding alias is also present.
            prime_source_alias_conflict = _copy_program(program)
            for row in prime_source_alias_conflict["instructions"]:
                if row.get("node_id") == _SIGMA_PRIME_NODE_ID:
                    prime_record = row["record"]
                    alias = dict(prime_record["derivative_source"])
                    alias["module"] = "synthetic_v4_prime_alias_conflict.py"
                    prime_record["source_binding"] = alias
                    row["original_record"]["source_binding"] = dict(alias)
                    break
            for row in prime_source_alias_conflict["nodes"].values():
                if row.get("node_id") == _SIGMA_PRIME_NODE_ID:
                    prime_record = row["record"]
                    alias = dict(prime_record["derivative_source"])
                    alias["module"] = "synthetic_v4_prime_alias_conflict.py"
                    prime_record["source_binding"] = alias
                    row["original_record"]["source_binding"] = dict(alias)
                    break
            prime_record = prime_source_alias_conflict["records"][_SIGMA_PRIME_NODE_ID]
            alias = dict(prime_record["derivative_source"])
            alias["module"] = "synthetic_v4_prime_alias_conflict.py"
            prime_record["source_binding"] = alias
            checks["prime_source_alias_conflict_rejected"] = _rejected(
                lambda: bind_entry(prime_source_alias_conflict, digests)
            )
        except Exception:
            _all_false(checks, names)
    else:
        _all_false(checks, names)

    return {
        "schema": SCHEMA,
        "component": "n1_flat_sigma_point_kernel_v4",
        "passed": all(checks.values()),
        "checks": checks,
        "synthetic_only": True,
        "used_old_results": False,
        "used_old_module": False,
        "ast_recomputed": False,
        "numeric_equivalence_certified": False,
        "selected": selected,
    }


run = run_fixture





