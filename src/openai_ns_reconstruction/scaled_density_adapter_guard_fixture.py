"""Direct adapter-boundary fixture for the scaled density input adapter.

This fixture exercises only adapter selection and row/coordinate validation.
It deliberately does not call a density, endpoint, correlation, graph, source
provider, or arithmetic consumer.
"""

from __future__ import annotations

from typing import Any, Callable

from openai_ns_reconstruction.original_finite_N_density_kernel import (
    PYTHON_OPERATIONS,
    ScalarOperations,
    ordinary_jet,
)
from openai_ns_reconstruction.original_scaled_signed_density import (
    ScaledDensityJetError,
    _as_jet,
    _constant_jet,
    _coordinate_jet,
    _row,
    _select_operations,
)


def _no_arithmetic(*_args: Any, **_kwargs: Any) -> Any:
    raise AssertionError("scaled adapter guard fixture must not evaluate arithmetic")


def _adapter(label: str, constant_log: list[Any] | None = None) -> ScalarOperations:
    def constant(value: Any) -> tuple[str, Any]:
        if constant_log is not None:
            constant_log.append(value)
        return ("constant", value)

    return ScalarOperations(
        _no_arithmetic,
        _no_arithmetic,
        _no_arithmetic,
        _no_arithmetic,
        _no_arithmetic,
        constant,
    )


class ivmpf:
    """Shape-only interval marker for the selector fallback check."""

    __module__ = "mpmath.ctx_iv"


def _check(condition: bool, label: str) -> dict[str, Any]:
    return {"passed": bool(condition), "label": label}


def _expect_exception(
    action: Callable[[], Any], expected: type[BaseException], label: str
) -> dict[str, Any]:
    try:
        action()
    except expected as exc:
        return {"passed": True, "label": label, "exception": str(exc)}
    except Exception as exc:  # pragma: no cover - diagnostic result only
        return {
            "passed": False,
            "label": label,
            "exception": type(exc).__name__ + ": " + str(exc),
        }
    return {"passed": False, "label": label, "exception": None}


def run_fixture() -> dict[str, Any]:
    """Return adapter-only checks without entering density arithmetic."""

    graph_a = _adapter("a")
    graph_b = _adapter("b")
    jet_a = ordinary_jet("a", "a_z1", "a_z2", operations=graph_a)
    jet_b = ordinary_jet("b", "b_z1", "b_z2", operations=graph_b)
    unbound = ordinary_jet("u", "u_z1", "u_z2")
    python_jet = ordinary_jet(
        "python", "python_z1", "python_z2", operations=PYTHON_OPERATIONS
    )

    checks: list[dict[str, Any]] = []

    checks.append(
        _check(
            _select_operations(None, unbound, jet_a) is graph_a,
            "unique_explicit_adapter_is_inherited",
        )
    )
    checks.append(
        _check(
            _select_operations(graph_a, jet_a, unbound) is graph_a,
            "same_requested_adapter_is_accepted",
        )
    )
    checks.append(
        _expect_exception(
            lambda: _select_operations(None, jet_a, jet_b),
            TypeError,
            "multiple_explicit_adapters_are_rejected",
        )
    )
    checks.append(
        _expect_exception(
            lambda: _select_operations(graph_b, jet_a),
            TypeError,
            "requested_adapter_conflict_is_rejected",
        )
    )
    checks.append(
        _check(
            _select_operations(None, unbound) is PYTHON_OPERATIONS,
            "python_fallback_is_preserved_without_explicit_adapter",
        )
    )
    checks.append(
        _check(
            _select_operations(None, python_jet) is PYTHON_OPERATIONS,
            "explicit_python_adapter_identity_is_preserved",
        )
    )
    checks.append(
        _check(
            _select_operations(None, python_jet, [ivmpf()]) is PYTHON_OPERATIONS,
            "explicit_python_takes_priority_over_interval_fallback",
        )
    )
    interval_adapter = _select_operations(None, [ivmpf(), ivmpf()])
    checks.append(
        _check(
            interval_adapter is not PYTHON_OPERATIONS,
            "interval_fallback_constructs_a_non_python_adapter",
        )
    )

    checks.append(
        _expect_exception(
            lambda: _row(jet_a, 2, "row", operations=graph_b),
            TypeError,
            "row_rejects_explicit_adapter_relabel",
        )
    )
    preserved_row = _as_jet(jet_a, 2, "row", operations=None)
    checks.append(
        _check(
            preserved_row.operations is graph_a
            and preserved_row.value == jet_a.value
            and preserved_row.dz == jet_a.dz
            and preserved_row.dzz == jet_a.dzz,
            "row_unpack_preserves_outer_adapter_identity",
        )
    )
    checks.append(
        _expect_exception(
            lambda: _coordinate_jet(jet_a, operations=graph_b),
            TypeError,
            "coordinate_rejects_explicit_adapter_relabel",
        )
    )
    preserved_coordinate = _coordinate_jet(jet_a, operations=None)
    checks.append(
        _check(
            preserved_coordinate.operations is graph_a
            and preserved_coordinate.value == jet_a.value
            and preserved_coordinate.dz == jet_a.dz
            and preserved_coordinate.dzz == jet_a.dzz,
            "coordinate_preserves_outer_adapter_identity",
        )
    )

    coordinate_number_log: list[Any] = []
    coordinate_number_ops = _adapter("coordinate_number", coordinate_number_log)
    coordinate_number = _coordinate_jet(7, operations=coordinate_number_ops)
    checks.append(
        _check(
            coordinate_number_log == [7, 1, 0]
            and coordinate_number.value == ("constant", 7)
            and coordinate_number.dz == ("constant", 1)
            and coordinate_number.dzz == ("constant", 0)
            and coordinate_number.operations is coordinate_number_ops,
            "coordinate_numeric_literal_and_zero_rows_use_constant",
        )
    )
    coordinate_opaque_log: list[Any] = []
    coordinate_opaque_ops = _adapter("coordinate_opaque", coordinate_opaque_log)
    coordinate_token = object()
    coordinate_opaque = _coordinate_jet(coordinate_token, operations=coordinate_opaque_ops)
    checks.append(
        _check(
            coordinate_opaque_log == [1, 0]
            and coordinate_opaque.value is coordinate_token
            and coordinate_opaque.dz == ("constant", 1)
            and coordinate_opaque.dzz == ("constant", 0)
            and coordinate_opaque.operations is coordinate_opaque_ops,
            "coordinate_opaque_source_is_preserved_and_coordinate_rows_are_constants",
        )
    )
    constant_number_log: list[Any] = []
    constant_number_ops = _adapter("constant_number", constant_number_log)
    constant_number = _constant_jet(11, constant_number_ops)
    checks.append(
        _check(
            constant_number_log == [11, 0, 0]
            and constant_number.value == ("constant", 11)
            and constant_number.dz == ("constant", 0)
            and constant_number.dzz == ("constant", 0)
            and constant_number.operations is constant_number_ops,
            "constant_numeric_literal_and_zero_rows_use_constant",
        )
    )
    constant_opaque_log: list[Any] = []
    constant_opaque_ops = _adapter("constant_opaque", constant_opaque_log)
    constant_token = object()
    constant_opaque = _constant_jet(constant_token, constant_opaque_ops)
    checks.append(
        _check(
            constant_opaque_log == [0, 0]
            and constant_opaque.value is constant_token
            and constant_opaque.dz == ("constant", 0)
            and constant_opaque.dzz == ("constant", 0)
            and constant_opaque.operations is constant_opaque_ops,
            "constant_opaque_source_is_preserved_and_zero_rows_are_constants",
        )
    )

    row_mapping_log: list[Any] = []
    row_mapping_ops = _adapter("row_mapping", row_mapping_log)
    row_mapping = _row(
        {"value": 3, "dz": 4, "dzz": 5},
        2,
        "row_mapping",
        operations=row_mapping_ops,
    )
    checks.append(
        _check(
            row_mapping_log == [3, 4, 5]
            and row_mapping == [("constant", 3), ("constant", 4), ("constant", 5)],
            "mapping_numeric_rows_use_constant_once",
        )
    )
    row_sequence_log: list[Any] = []
    row_sequence_ops = _adapter("row_sequence", row_sequence_log)
    row_sequence = _as_jet(
        [6, 7, 8], 2, "row_sequence", operations=row_sequence_ops
    )
    checks.append(
        _check(
            row_sequence_log == [6, 7, 8]
            and row_sequence.value == ("constant", 6)
            and row_sequence.dz == ("constant", 7)
            and row_sequence.dzz == ("constant", 8)
            and row_sequence.operations is row_sequence_ops,
            "sequence_numeric_rows_use_constant_once",
        )
    )
    opaque_mapping_log: list[Any] = []
    opaque_mapping_ops = _adapter("opaque_mapping", opaque_mapping_log)
    opaque_value = object()
    opaque_dz = object()
    opaque_dzz = object()
    opaque_mapping = _as_jet(
        {"value": opaque_value, "dz": opaque_dz, "dzz": opaque_dzz},
        2,
        "opaque_mapping",
        operations=opaque_mapping_ops,
    )
    checks.append(
        _check(
            opaque_mapping_log == []
            and opaque_mapping.value is opaque_value
            and opaque_mapping.dz is opaque_dz
            and opaque_mapping.dzz is opaque_dzz,
            "mapping_opaque_rows_are_preserved",
        )
    )
    opaque_sequence_log: list[Any] = []
    opaque_sequence_ops = _adapter("opaque_sequence", opaque_sequence_log)
    opaque_sequence_values = [object(), object(), object()]
    opaque_sequence = _as_jet(
        opaque_sequence_values, 2, "opaque_sequence", operations=opaque_sequence_ops
    )
    checks.append(
        _check(
            opaque_sequence_log == []
            and opaque_sequence.value is opaque_sequence_values[0]
            and opaque_sequence.dz is opaque_sequence_values[1]
            and opaque_sequence.dzz is opaque_sequence_values[2],
            "sequence_opaque_rows_are_preserved",
        )
    )
    top_level_log: list[Any] = []
    top_level_ops = _adapter("top_level", top_level_log)
    top_level = _as_jet(
        ordinary_jet(9, 10, 11), 2, "top_level", operations=top_level_ops
    )
    checks.append(
        _check(
            top_level_log == []
            and top_level.value == 9
            and top_level.dz == 10
            and top_level.dzz == 11
            and top_level.operations is top_level_ops,
            "top_level_ordinary_jet_rows_are_not_double_converted",
        )
    )

    nested_sequence = [jet_a, "nested_z1", "nested_z2"]
    nested_mapping = {"value": jet_a, "dz": "nested_z1", "dzz": "nested_z2"}
    checks.append(
        _expect_exception(
            lambda: _select_operations(None, nested_sequence),
            ScaledDensityJetError,
            "nested_sequence_jet_is_rejected_not_borrowed",
        )
    )
    checks.append(
        _expect_exception(
            lambda: _select_operations(None, nested_mapping),
            ScaledDensityJetError,
            "nested_mapping_jet_is_rejected_not_borrowed",
        )
    )
    checks.append(
        _expect_exception(
            lambda: _as_jet(nested_sequence, 2, "nested", operations=graph_b),
            ScaledDensityJetError,
            "row_nested_sequence_is_invalid",
        )
    )
    checks.append(
        _expect_exception(
            lambda: _coordinate_jet(nested_sequence, operations=graph_b),
            ScaledDensityJetError,
            "coordinate_nested_sequence_is_invalid",
        )
    )

    passed = all(item["passed"] for item in checks)
    return {
        "schema": "scaled_density_adapter_guard_v1",
        "passed": passed,
        "source_only": True,
        "guard_only": True,
        "density_calls": 0,
        "endpoint_density_calls": 0,
        "correlation_calls": 0,
        "graph_calls": 0,
        "numeric_values_computed": False,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
        "checks": checks,
    }


__all__ = ["run_fixture"]
