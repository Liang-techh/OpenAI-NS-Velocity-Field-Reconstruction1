"""Synthetic-only checks for :mod:`n1_raw_g5_dag_program_v1`.

The fixture creates a fictional indexed graph in memory.  It does not read an
archive, import a producer, or use any saved result.  ``run_fixture(module)``
receives the new compiler/executor module explicitly and is the only module
implementation it calls.
"""

from __future__ import annotations

from fractions import Fraction
import copy
import json
from typing import Any


__all__ = ["SCHEMA", "run_fixture", "run"]


SCHEMA = "n1_raw_g5_dag_program_fixture_v1"


class _ExactBackend:
    """Small exact backend used only by this fictional fixture."""

    def rational(self, numerator: int, denominator: int) -> Fraction:
        return Fraction(numerator, denominator)

    def add(self, left: Fraction, right: Fraction) -> Fraction:
        return left + right

    def multiply(self, left: Fraction, right: Fraction) -> Fraction:
        return left * right

    def negative(self, value: Fraction) -> Fraction:
        return -value

    def divide(self, numerator: Fraction, denominator: Fraction) -> Fraction:
        return numerator / denominator

    def is_positive(self, value: Fraction) -> bool:
        return value > 0

    def unary(self, function: str, value: Fraction) -> Fraction:
        if function == "abs":
            return abs(value)
        raise ValueError(f"fixture backend does not provide {function}")


class _RejectingScalarBackend(_ExactBackend):
    """Backend negative case: advertised scalar validation must be honoured."""

    def validate_scalar(self, _value: Any) -> bool:
        return False


class _PositiveDivWithoutPositiveCheckBackend(_ExactBackend):
    """positive_div must not bypass the compiler's strict-positive guard."""

    def positive_div(self, numerator: Fraction, denominator: Fraction) -> Fraction:
        return numerator / denominator

    def is_positive(self, _value: Fraction) -> bool:
        return False


class _ArbitraryScalar:
    """Non-Real object used to check the explicit scalar validator gate."""


class _AcceptingArbitraryBackend(_ExactBackend):
    def validate_scalar(self, value: Any) -> bool:
        return isinstance(value, _ArbitraryScalar)


def _synthetic_graph() -> dict[str, Any]:
    """Return a fresh graph with shared nodes, multiple roots, and one Z1."""

    return {
        "exact_graph_nodes": [
            {"operation": "exact_rational", "numerator": 1, "denominator": 2},  # 0
            {"operation": "exact_rational", "numerator": 1, "denominator": 3},  # 1
            {"operation": "sum", "arguments": [0, 1]},  # 2 = 5/6
            {"operation": "product", "arguments": [2, 0]},  # 3 = 5/12
            {"operation": "negative", "argument": 3},  # 4 = -5/12
            {"operation": "positive_quotient", "numerator": 4, "denominator": 0},  # 5
            {"operation": "analytic_unary", "function": "abs", "argument": 5},  # 6
            {
                "operation": "fictional_future_source",
                "ordinary_Z": "opaque",
                # This nested payload must not be interpreted as a graph edge.
                "payload": {"argument": 999, "binder": {"children": [998]}},
            },  # 7
            {"operation": "sum", "arguments": [6, 7]},  # 8
            {
                "operation": "exact_rational",
                "numerator": 7,
                "denominator": 11,
                "ordinary_Z": "Z1",
            },  # 9
        ],
        "source_family": "fictional-g5-fixture",
    }


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def run_fixture(module: Any) -> dict[str, Any]:
    """Run synthetic acceptance and fail-closed checks against ``module``."""

    checks: dict[str, bool] = {}
    selected: dict[str, Any] = {}
    compile_program = getattr(module, "compile_program", None)
    execute = getattr(module, "execute", None)
    checks["compiler_api"] = callable(compile_program)
    checks["executor_api"] = callable(execute)

    graph = _synthetic_graph()
    roots = [3, 8, 9]
    backend = _ExactBackend()
    oracle_calls: list[tuple[int, str, Any]] = []

    def oracle(node_id: int, original_record: Any, operation: str) -> Fraction:
        oracle_calls.append((node_id, operation, copy.deepcopy(original_record)))
        if node_id != 7 or operation != "unknown":
            raise ValueError(f"unexpected opaque node {node_id}/{operation}")
        return Fraction(1, 12)

    if checks["compiler_api"] and checks["executor_api"]:
        try:
            program = compile_program(graph, roots, "fixture-g5-graph-hash")
            result = execute(program, backend, oracle)
            root_values = result.get("root_values", [])
            selected = {
                "program": {
                    "schema": program.get("schema"),
                    "roots": list(program.get("roots", [])),
                    "topological_order": list(program.get("topological_order", [])),
                    "missing_oracle_nodes": list(program.get("missing_oracle_nodes", [])),
                },
                "result": {
                    "schema": result.get("schema"),
                    "status": result.get("status"),
                    "root_values": [str(value) for value in root_values],
                    "opaque_nodes_used": list(result.get("opaque_nodes_used", [])),
                },
            }
            checks["compiled"] = program.get("status") == "compiled"
            checks["graph_hash_preserved"] = program.get("graph_hash") == "fixture-g5-graph-hash"
            checks["multiple_roots"] = program.get("roots") == roots
            checks["shared_node_once"] = (
                list(program.get("topological_order", [])).count(0) == 1
                and len(program.get("instructions", []))
                == len(set(program.get("topological_order", [])))
            )
            checks["exact_fraction_arithmetic"] = root_values[:1] == [Fraction(5, 12)]
            checks["ordinary_Z1_root_retained"] = (
                program.get("records", {}).get(9, {}).get("ordinary_Z") == "Z1"
                or program.get("records", {}).get("9", {}).get("ordinary_Z") == "Z1"
            )
            checks["unknown_only_through_oracle"] = (
                bool(oracle_calls)
                and [node_id for node_id, _, _ in oracle_calls] == [7]
                and [operation for _, operation, _ in oracle_calls] == ["unknown"]
                and program.get("missing_oracle_nodes") == [7]
            )
            checks["conditional_not_certified"] = (
                result.get("status") == "conditional"
                and result.get("conditional") is True
                and result.get("certified") is False
            )
            claims = result.get("claims", {})
            checks["no_physical_claims"] = all(
                claims.get(name) is False
                for name in ("function_equality", "physical_m15", "Kchi", "Navier_Stokes")
            )
        except Exception:
            for name in (
                "compiled",
                "graph_hash_preserved",
                "multiple_roots",
                "shared_node_once",
                "exact_fraction_arithmetic",
                "ordinary_Z1_root_retained",
                "unknown_only_through_oracle",
                "conditional_not_certified",
                "no_physical_claims",
            ):
                checks[name] = False

        # Missing oracle is an execution error, never a zero/default value.
        try:
            program = compile_program(graph, roots, "fixture-g5-graph-hash")
            checks["missing_oracle_fails_closed"] = _rejected(
                lambda: execute(program, backend, None)
            )
        except Exception:
            checks["missing_oracle_fails_closed"] = False

        cycle_graph = {
            "exact_graph_nodes": [
                {"operation": "sum", "arguments": [1]},
                {"operation": "negative", "argument": 0},
            ]
        }
        checks["cycle_rejected"] = _rejected(
            lambda: compile_program(cycle_graph, [0], "fixture-cycle")
        )

        invalid_graph = {
            "exact_graph_nodes": [
                {"operation": "sum", "arguments": [99]},
            ]
        }
        checks["invalid_id_rejected"] = _rejected(
            lambda: compile_program(invalid_graph, [0], "fixture-invalid-id")
        )
        checks["boolean_root_rejected"] = _rejected(
            lambda: compile_program(graph, [True], "fixture-bool-id")
        )
        rational_conflict_graph = {
            "exact_graph_nodes": [
                {
                    "operation": "exact_rational",
                    "numerator": 1,
                    "value": 2,
                    "denominator": 3,
                }
            ]
        }
        checks["rational_alias_conflict_rejected"] = _rejected(
            lambda: compile_program(rational_conflict_graph, [0], "fixture-rational-conflict")
        )

        # Oracle outputs are scalar values, never source records/cap metadata,
        # non-finite floats, or an index placeholder.
        opaque_program = compile_program(graph, [8], "fixture-opaque-negative")
        checks["oracle_mapping_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda *_args: {"cap": 1})
        )
        checks["oracle_list_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda *_args: [1, 2])
        )
        checks["oracle_tuple_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda *_args: (1, 2))
        )
        checks["oracle_set_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda *_args: {1, 2})
        )
        checks["oracle_nan_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda *_args: float("nan"))
        )
        checks["oracle_inf_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda *_args: float("inf"))
        )
        checks["oracle_node_id_rejected"] = _rejected(
            lambda: execute(opaque_program, backend, lambda node_id, *_args: node_id)
        )
        checks["backend_scalar_validator_enforced"] = _rejected(
            lambda: execute(opaque_program, _RejectingScalarBackend(), oracle)
        )
        checks["oracle_complex_rejected"] = _rejected(
            lambda: execute(
                compile_program(graph, [7], "fixture-complex"),
                backend,
                lambda *_args: complex(1, 0),
            )
        )
        arbitrary_program = compile_program(graph, [7], "fixture-arbitrary")
        checks["oracle_arbitrary_object_requires_validator"] = _rejected(
            lambda: execute(arbitrary_program, backend, lambda *_args: _ArbitraryScalar())
        )
        try:
            arbitrary_result = execute(
                arbitrary_program,
                _AcceptingArbitraryBackend(),
                lambda *_args: _ArbitraryScalar(),
            )
            checks["oracle_arbitrary_object_validator_gate"] = (
                len(arbitrary_result.get("root_values", [])) == 1
                and isinstance(arbitrary_result["root_values"][0], _ArbitraryScalar)
            )
        except Exception:
            checks["oracle_arbitrary_object_validator_gate"] = False

        positive_div_program = compile_program(graph, [5], "fixture-positive-div")
        checks["positive_div_still_checks_positive"] = _rejected(
            lambda: execute(positive_div_program, _PositiveDivWithoutPositiveCheckBackend(), None)
        )

        integral_graph = {
            "exact_graph_nodes": [
                {
                    "operation": "definite_integral",
                    "binder": {"kind": "fresh_binder", "token": "fixture-t"},
                    "lower": {"operation": "exact_rational", "numerator": 0},
                    "upper": {"operation": "exact_rational", "numerator": 1},
                    "integrand": {"operation": "bound_variable", "name": "fixture-t"},
                }
            ]
        }
        integral_program = compile_program(integral_graph, [0], "fixture-integral")
        checks["integral_runtime_rejected"] = _rejected(
            lambda: execute(
                integral_program,
                backend,
                {"integral": lambda *_args: Fraction(1, 1)},
            )
        )
        substitution_graph = {
            "exact_graph_nodes": [
                {
                    "operation": "function_substitution",
                    "expression": 1001,
                    "variable": 1002,
                    "value": 1003,
                }
            ]
        }
        substitution_program = compile_program(substitution_graph, [0], "fixture-substitution")
        checks["substitution_runtime_rejected"] = _rejected(
            lambda: execute(
                substitution_program,
                backend,
                {"substitution": lambda *_args: Fraction(1, 1)},
            )
        )

        regular_component_graph = {
            "exact_graph_nodes": [
                {"operation": "current_original_unique_regular_n1_component"}
            ]
        }
        regular_component_program = compile_program(
            regular_component_graph, [0], "fixture-regular-component"
        )
        regular_calls: list[str] = []

        def regular_component_oracle(
            _node_id: int, _record: Any, operation: str
        ) -> Fraction:
            regular_calls.append(operation)
            return Fraction(2, 1)

        regular_result = execute(regular_component_program, backend, regular_component_oracle)
        checks["regular_component_category_explicit"] = (
            regular_calls == ["regular_component"]
            and regular_result.get("root_values") == [Fraction(2, 1)]
        )

        variable_graph = {
            "exact_graph_nodes": [{"operation": "bound_variable", "name": "fixture-x"}]
        }
        variable_program = compile_program(variable_graph, [0], "fixture-variable")
        checks["wrong_category_rejected"] = _rejected(
            lambda: execute(variable_program, backend, {"unknown": lambda *_args: Fraction(1, 1)})
        )
        checks["implicit_resolve_rejected"] = _rejected(
            lambda: execute(
                variable_program,
                backend,
                {"resolve": lambda *_args: Fraction(1, 1)},
            )
        )
    else:
        for name in (
            "compiled",
            "graph_hash_preserved",
            "multiple_roots",
            "shared_node_once",
            "exact_fraction_arithmetic",
            "ordinary_Z1_root_retained",
            "unknown_only_through_oracle",
            "conditional_not_certified",
            "no_physical_claims",
            "missing_oracle_fails_closed",
            "cycle_rejected",
            "invalid_id_rejected",
            "boolean_root_rejected",
            "rational_alias_conflict_rejected",
            "oracle_mapping_rejected",
            "oracle_list_rejected",
            "oracle_tuple_rejected",
            "oracle_set_rejected",
            "oracle_nan_rejected",
            "oracle_inf_rejected",
            "oracle_node_id_rejected",
            "backend_scalar_validator_enforced",
            "oracle_complex_rejected",
            "oracle_arbitrary_object_requires_validator",
            "oracle_arbitrary_object_validator_gate",
            "positive_div_still_checks_positive",
            "integral_runtime_rejected",
            "substitution_runtime_rejected",
            "regular_component_category_explicit",
            "wrong_category_rejected",
            "implicit_resolve_rejected",
        ):
            checks[name] = False

    output: dict[str, Any] = {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "synthetic_only": True,
        "used_old_results": False,
        "selected": selected,
    }
    try:
        json.dumps(output, sort_keys=True)
    except (TypeError, ValueError):
        output["json_safe"] = False
        output["passed"] = False
    else:
        output["json_safe"] = True
    return output


run = run_fixture

