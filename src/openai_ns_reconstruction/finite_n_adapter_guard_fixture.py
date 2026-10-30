"""Direct fixture for finite-N scalar-adapter identity guards.

This fixture exercises only the adapter-resolution helpers in
``original_finite_N_density_kernel``.  It deliberately does not call
``_binary``, a density or endpoint consumer, a correlation path, a source
provider, or any archived mathematical expression.  The graph callbacks
raise if arithmetic is accidentally attempted; the fixture checks only
adapter identity, ordinary-jet row preservation, and the intended rejection
or inheritance behavior.
"""

from __future__ import annotations

from typing import Any, Callable

from openai_ns_reconstruction.original_finite_N_density_kernel import (
    PYTHON_OPERATIONS,
    ScalarOperations,
    _common_operations,
    _normalize,
    coerce_ordinary_jet,
    ordinary_jet,
)


def _no_arithmetic(*_args: Any, **_kwargs: Any) -> Any:
    raise AssertionError("adapter guard fixture must not evaluate graph arithmetic")


def _adapter(constant_log: list[Any] | None = None) -> ScalarOperations:
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


def _expect_type_error(action: Callable[[], Any], label: str) -> dict[str, Any]:
    try:
        action()
    except TypeError as exc:
        return {"passed": True, "label": label, "exception": str(exc)}
    except Exception as exc:  # pragma: no cover - diagnostic result only
        return {
            "passed": False,
            "label": label,
            "exception": type(exc).__name__ + ": " + str(exc),
        }
    return {"passed": False, "label": label, "exception": None}


def _check(condition: bool, label: str) -> dict[str, Any]:
    return {"passed": bool(condition), "label": label}


def run_fixture() -> dict[str, Any]:
    """Return direct adapter-guard checks without evaluating old math."""

    graph_a = _adapter()
    graph_b = _adapter()
    jet_a = ordinary_jet("a", "a_z1", "a_z2", operations=graph_a)
    jet_b = ordinary_jet("b", "b_z1", "b_z2", operations=graph_b)
    unbound = ordinary_jet("u", "u_z1", "u_z2")
    python_jet = ordinary_jet(
        "python", "python_z1", "python_z2", operations=PYTHON_OPERATIONS
    )

    checks: list[dict[str, Any]] = []

    # An already-bound graph jet cannot be relabeled by a requested adapter.
    # The three labels document the consumer boundaries that must receive the
    # same rejection before any density, endpoint, or correlation code runs.
    for context in ("density", "endpoint", "correlation"):
        checks.append(
            _expect_type_error(
                lambda: coerce_ordinary_jet(
                    jet_a, operations=graph_b, name="existing_" + context
                ),
                "coerce_explicit_override_before_" + context,
            )
        )
    checks.append(
        _expect_type_error(
            lambda: coerce_ordinary_jet(
                jet_a, operations=PYTHON_OPERATIONS, name="explicit_python_override"
            ),
            "coerce_custom_to_explicit_python_rejected",
        )
    )
    checks.append(
        _expect_type_error(
            lambda: coerce_ordinary_jet(
                python_jet, operations=graph_a, name="explicit_graph_override"
            ),
            "coerce_explicit_python_to_custom_rejected",
        )
    )

    opaque = object()
    opaque_jet = coerce_ordinary_jet(opaque, name="opaque")
    opaque_bound = _normalize(opaque_jet, graph_a)
    checks.append(_check(opaque_bound.value is opaque and opaque_bound.operations is graph_a,
                         "opaque_unbound_source_not_converted_to_literal"))

    # Multiple explicit adapters and a conflicting requested adapter are
    # rejected before any consumer gets to normalize its fields.
    checks.append(
        _expect_type_error(
            lambda: _common_operations(jet_a, jet_b),
            "common_multiple_explicit_adapters_rejected",
        )
    )
    checks.append(
        _expect_type_error(
            lambda: _common_operations(jet_a, operations=graph_b),
            "common_requested_adapter_conflict_rejected",
        )
    )
    checks.append(
        _expect_type_error(
            lambda: _common_operations(jet_a, python_jet),
            "common_custom_and_explicit_python_conflict_rejected",
        )
    )
    checks.append(
        _expect_type_error(
            lambda: _normalize(jet_a, graph_b),
            "normalize_explicit_adapter_relabel_rejected",
        )
    )
    checks.append(
        _expect_type_error(
            lambda: _normalize(python_jet, graph_a),
            "normalize_explicit_python_relabel_rejected",
        )
    )

    # Bare scalars and unbound rows inherit only the unique requested adapter.
    bare = coerce_ordinary_jet(17, name="bare")
    checks.append(
        _check(
            bare.operations is None and bare.value == 17 and bare.dz == 0 and bare.dzz == 0,
            "bare_scalar_without_operations_is_unbound_zero_derivative",
        )
    )
    bound_bare = coerce_ordinary_jet(17, operations=graph_a, name="bound_bare")
    checks.append(
        _check(
            bound_bare.operations is graph_a
            and bound_bare.value == ("constant", 17)
            and bound_bare.dz == ("constant", 0)
            and bound_bare.dzz == ("constant", 0),
            "bare_scalar_inherits_requested_adapter_and_uses_constant_callback",
        )
    )
    inherited = _normalize(unbound, graph_a)
    checks.append(
        _check(
            inherited.operations is graph_a
            and inherited.value == unbound.value
            and inherited.dz == unbound.dz
            and inherited.dzz == unbound.dzz,
            "unbound_jet_inherits_unique_adapter_and_preserves_rows",
        )
    )
    delayed_log: list[Any] = []
    delayed_adapter = _adapter(delayed_log)
    delayed_literal = coerce_ordinary_jet(23, name="delayed_literal")
    delayed_bound = _normalize(delayed_literal, delayed_adapter)
    checks.append(
        _check(
            delayed_log == [23, 0, 0]
            and delayed_bound.operations is delayed_adapter
            and delayed_bound.value == ("constant", 23)
            and delayed_bound.dz == ("constant", 0)
            and delayed_bound.dzz == ("constant", 0),
            "unbound_literal_delays_graph_constant_conversion_until_normalize",
        )
    )
    coerce_log: list[Any] = []
    coerce_adapter = _adapter(coerce_log)
    marker = coerce_ordinary_jet(29, name="marker")
    coerced_marker = coerce_ordinary_jet(marker, operations=coerce_adapter, name="marker")
    checks.append(
        _check(
            coerce_log == [29, 0, 0]
            and coerced_marker.operations is coerce_adapter
            and coerced_marker.value == ("constant", 29)
            and coerced_marker.dz == ("constant", 0)
            and coerced_marker.dzz == ("constant", 0),
            "coerce_marker_uses_delayed_normalize_conversion",
        )
    )
    ordinary_log: list[Any] = []
    ordinary_adapter = _adapter(ordinary_log)
    ordinary_unbound = ordinary_jet("opaque", "opaque_z1", "opaque_z2")
    ordinary_bound = _normalize(ordinary_unbound, ordinary_adapter)
    checks.append(
        _check(
            ordinary_log == []
            and ordinary_bound.operations is ordinary_adapter
            and ordinary_bound.value == "opaque"
            and ordinary_bound.dz == "opaque_z1"
            and ordinary_bound.dzz == "opaque_z2",
            "ordinary_unbound_jet_does_not_trigger_literal_conversion",
        )
    )
    coerced_unbound = coerce_ordinary_jet(unbound, operations=graph_a, name="unbound")
    checks.append(
        _check(
            coerced_unbound.operations is graph_a
            and coerced_unbound.value == unbound.value
            and coerced_unbound.dz == unbound.dz
            and coerced_unbound.dzz == unbound.dzz,
            "coerce_unbound_jet_inherits_requested_adapter_and_preserves_rows",
        )
    )
    checks.append(
        _check(
            _common_operations(unbound, jet_a) is graph_a,
            "common_operations_selects_unique_explicit_adapter",
        )
    )
    checks.append(
        _check(
            _common_operations(jet_a, jet_a, operations=graph_a) is graph_a,
            "common_operations_accepts_same_requested_adapter",
        )
    )
    checks.append(
        _check(
            coerce_ordinary_jet(jet_a, operations=graph_a) is jet_a,
            "coerce_accepts_same_explicit_adapter",
        )
    )
    checks.append(
        _check(
            _normalize(jet_a, graph_a) is jet_a,
            "normalize_same_adapter_preserves_jet_identity",
        )
    )

    # Explicit Python is a real singleton adapter identity.  It is compatible
    # with unbound rows and itself, but not with a custom graph adapter.
    checks.append(
        _check(
            _common_operations(python_jet, unbound) is PYTHON_OPERATIONS,
            "explicit_python_is_singleton_adapter",
        )
    )
    checks.append(
        _check(
            _common_operations(python_jet, operations=PYTHON_OPERATIONS)
            is PYTHON_OPERATIONS,
            "explicit_python_requested_identity_is_preserved",
        )
    )
    checks.append(
        _check(
            _normalize(python_jet, PYTHON_OPERATIONS) is python_jet,
            "normalize_explicit_python_identity_is_preserved",
        )
    )
    checks.append(
        _check(
            coerce_ordinary_jet(python_jet, operations=PYTHON_OPERATIONS)
            is python_jet,
            "coerce_explicit_python_identity_is_preserved",
        )
    )
    python_inherited = _normalize(unbound, PYTHON_OPERATIONS)
    checks.append(
        _check(
            python_inherited.operations is PYTHON_OPERATIONS
            and python_inherited.value == unbound.value
            and python_inherited.dz == unbound.dz
            and python_inherited.dzz == unbound.dzz,
            "unbound_jet_can_inherit_explicit_python",
        )
    )

    # Mapping rows are also unbound until a requested adapter is supplied;
    # this checks ordinary Z1 / Z2 row preservation without doing convolution.
    row = coerce_ordinary_jet(
        {"value": "row", "dz": "row_z1", "dzz": "row_z2"},
        operations=graph_a,
        name="row",
    )
    checks.append(
        _check(
            row.operations is graph_a
            and row.value == "row"
            and row.dz == "row_z1"
            and row.dzz == "row_z2",
            "ordinary_z_rows_preserved_under_adapter_binding",
        )
    )

    passed = all(item["passed"] for item in checks)
    return {
        "schema": "finite_n_adapter_guard_v1",
        "passed": passed,
        "source_only": True,
        "guard_only": True,
        "ordinary_Z_orders": [0, 1, 2],
        "old_math_calls": 0,
        "binary_calls": 0,
        "volume_density_calls": 0,
        "endpoint_density_calls": 0,
        "correlation_calls": 0,
        "numeric_values_computed": False,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
        "checks": checks,
    }


__all__ = ["run_fixture"]
