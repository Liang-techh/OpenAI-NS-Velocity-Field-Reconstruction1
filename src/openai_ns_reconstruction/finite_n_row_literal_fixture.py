"""Direct row-literal and adapter-owned-zero fixture.

This fixture exercises only ordinary-jet coercion, delayed row binding, the
bare-number path in ``_other``, and static source checks for accumulator zeros.
It deliberately does not call ``_binary``, a density or endpoint consumer, a
correlation path, a source provider, a graph, or an older fixture.
"""

from __future__ import annotations

from fractions import Fraction
import inspect
from typing import Any, Callable

from openai_ns_reconstruction.original_finite_N_density_kernel import (
    OrdinaryJet,
    ScalarOperations,
    _normalize,
    coerce_ordinary_jet,
    endpoint_density,
    factored_endpoint_density_difference,
    factored_volume_density_difference,
    full_endpoint_density_difference,
    full_volume_density_difference,
    ordinary_jet,
    volume_density,
    _sum_jets,
)


def _no_arithmetic(*_args: Any, **_kwargs: Any) -> Any:
    raise AssertionError("row-literal fixture must not enter arithmetic")


def _adapter(constant_log: list[Any]) -> ScalarOperations:
    def constant(value: Any) -> tuple[str, Any]:
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


def _check(condition: bool, label: str) -> dict[str, Any]:
    return {"passed": bool(condition), "label": label}


def _lifted_row(values: tuple[Any, Any, Any], log: list[Any], ops: ScalarOperations) -> Any:
    row = coerce_ordinary_jet(
        {"value": values[0], "dz": values[1], "dzz": values[2]},
        operations=ops,
        name="mapping",
    )
    return row, log


def run_fixture() -> dict[str, Any]:
    """Return isolated row-literal checks without evaluating old mathematics."""

    checks: list[dict[str, Any]] = []

    mapping_log: list[Any] = []
    mapping_ops = _adapter(mapping_log)
    mapping_values = (Fraction(3, 7), Fraction(5, 11), Fraction(7, 13))
    mapping, _ = _lifted_row(mapping_values, mapping_log, mapping_ops)
    checks.append(
        _check(
            mapping_log == list(mapping_values)
            and mapping.value == ("constant", mapping_values[0])
            and mapping.dz == ("constant", mapping_values[1])
            and mapping.dzz == ("constant", mapping_values[2])
            and mapping.operations is mapping_ops,
            "explicit_mapping_numeric_rows_use_constant",
        )
    )

    sequence_log: list[Any] = []
    sequence_ops = _adapter(sequence_log)
    sequence_values = (Fraction(2, 3), 17, 19)
    sequence = coerce_ordinary_jet(sequence_values, operations=sequence_ops, name="sequence")
    checks.append(
        _check(
            sequence_log == list(sequence_values)
            and sequence.value == ("constant", sequence_values[0])
            and sequence.dz == ("constant", sequence_values[1])
            and sequence.dzz == ("constant", sequence_values[2])
            and sequence.operations is sequence_ops,
            "explicit_sequence_numeric_rows_use_constant",
        )
    )

    opaque_value = object()
    opaque_z = object()
    opaque_zz = object()
    opaque_log: list[Any] = []
    opaque_ops = _adapter(opaque_log)
    opaque = coerce_ordinary_jet(
        {"value": opaque_value, "dz": opaque_z, "dzz": opaque_zz},
        operations=opaque_ops,
        name="opaque_mapping",
    )
    checks.append(
        _check(
            opaque_log == []
            and opaque.value is opaque_value
            and opaque.dz is opaque_z
            and opaque.dzz is opaque_zz
            and opaque.operations is opaque_ops,
            "explicit_mapping_opaque_rows_preserve_identity",
        )
    )

    default_log: list[Any] = []
    default_ops = _adapter(default_log)
    default_value = object()
    default = coerce_ordinary_jet(
        {"value": default_value, "dzz": 0},
        operations=default_ops,
        name="default_mapping",
    )
    checks.append(
        _check(
            default_log == [0, 0]
            and default.value is default_value
            and default.dz == ("constant", 0)
            and default.dzz == ("constant", 0),
            "mapping_default_zero_rows_use_constant",
        )
    )

    delayed_mapping_log: list[Any] = []
    delayed_mapping_ops = _adapter(delayed_mapping_log)
    delayed_value = Fraction(4, 9)
    delayed_z = object()
    delayed_zz = Fraction(8, 13)
    delayed_mapping = coerce_ordinary_jet(
        {"value": delayed_value, "dz": delayed_z, "dzz": delayed_zz},
        name="delayed_mapping",
    )
    delayed_mapping_bound = _normalize(delayed_mapping, delayed_mapping_ops)
    checks.append(
        _check(
            delayed_mapping_log == [delayed_value, delayed_zz]
            and delayed_mapping_bound.value == ("constant", delayed_value)
            and delayed_mapping_bound.dz is delayed_z
            and delayed_mapping_bound.dzz == ("constant", delayed_zz)
            and delayed_mapping_bound.operations is delayed_mapping_ops,
            "delayed_mapping_lifts_numeric_only",
        )
    )

    delayed_sequence_log: list[Any] = []
    delayed_sequence_ops = _adapter(delayed_sequence_log)
    delayed_sequence_value = object()
    delayed_sequence_z = Fraction(10, 17)
    delayed_sequence_zz = object()
    delayed_sequence = coerce_ordinary_jet(
        (delayed_sequence_value, delayed_sequence_z, delayed_sequence_zz),
        name="delayed_sequence",
    )
    delayed_sequence_bound = _normalize(delayed_sequence, delayed_sequence_ops)
    checks.append(
        _check(
            delayed_sequence_log == [delayed_sequence_z]
            and delayed_sequence_bound.value is delayed_sequence_value
            and delayed_sequence_bound.dz == ("constant", delayed_sequence_z)
            and delayed_sequence_bound.dzz is delayed_sequence_zz
            and delayed_sequence_bound.operations is delayed_sequence_ops,
            "delayed_sequence_lifts_numeric_only",
        )
    )

    direct_log: list[Any] = []
    direct_ops = _adapter(direct_log)
    direct_value = object()
    direct_z = Fraction(11, 19)
    direct_zz = object()
    direct = ordinary_jet(direct_value, direct_z, direct_zz)
    direct_bound = _normalize(direct, direct_ops)
    checks.append(
        _check(
            direct_log == []
            and direct_bound.value is direct_value
            and direct_bound.dz == direct_z
            and direct_bound.dzz is direct_zz
            and direct_bound.operations is direct_ops,
            "direct_ordinary_jet_rows_are_not_literal_marked",
        )
    )

    scalar_log: list[Any] = []
    scalar_ops = _adapter(scalar_log)
    scalar_marker = coerce_ordinary_jet(23, name="delayed_scalar")
    scalar_bound = _normalize(scalar_marker, scalar_ops)
    checks.append(
        _check(
            scalar_log == [23, 0, 0]
            and scalar_bound.value == ("constant", 23)
            and scalar_bound.dz == ("constant", 0)
            and scalar_bound.dzz == ("constant", 0),
            "bare_number_literal_still_uses_delayed_constant",
        )
    )

    other_log: list[Any] = []
    other_ops = _adapter(other_log)
    source_value = object()
    source_jet = ordinary_jet(source_value, object(), object(), operations=other_ops)
    numeric_other = source_jet._other(Fraction(6, 17))
    checks.append(
        _check(
            other_log == [Fraction(6, 17), 0, 0]
            and numeric_other.value == ("constant", Fraction(6, 17))
            and numeric_other.dz == ("constant", 0)
            and numeric_other.dzz == ("constant", 0)
            and numeric_other.operations is other_ops,
            "other_numeric_literal_uses_constant",
        )
    )

    opaque_other_log: list[Any] = []
    opaque_other_ops = _adapter(opaque_other_log)
    opaque_other_value = object()
    opaque_source_jet = ordinary_jet("source", object(), object(), operations=opaque_other_ops)
    opaque_other = opaque_source_jet._other(opaque_other_value)
    checks.append(
        _check(
            opaque_other_log == [0, 0]
            and opaque_other.value is opaque_other_value
            and opaque_other.operations is opaque_other_ops,
            "other_opaque_value_preserves_identity",
        )
    )

    accumulator_functions: tuple[tuple[str, Callable[..., Any]], ...] = (
        ("volume_density", volume_density),
        ("endpoint_density", endpoint_density),
        ("factored_endpoint_density_difference", factored_endpoint_density_difference),
        ("factored_volume_density_difference", factored_volume_density_difference),
        ("_sum_jets", _sum_jets),
    )
    accumulator_sources = {
        name: inspect.getsource(function) for name, function in accumulator_functions
    }
    checks.append(
        _check(
            all("_numeric_constant_jet(0" in source for source in accumulator_sources.values())
            and all("OrdinaryJet.constant(0" not in source for source in accumulator_sources.values()),
            "accumulators_preconvert_zero_through_adapter",
        )
    )

    binary_source = inspect.getsource(OrdinaryJet._binary)
    checks.append(_check("a = _normalize(self, ops)" in binary_source and
                         "b = _normalize(b, ops)" in binary_source,
                         "binary_uses_guarded_normalization_static"))
    passed = all(item["passed"] for item in checks)
    return {
        "schema": "finite_n_row_literal_v1",
        "passed": passed,
        "source_only": True,
        "row_literal_only": True,
        "ordinary_Z_orders": [0, 1, 2],
        "old_math_calls": 0,
        "binary_calls": 0,
        "volume_density_calls": 0,
        "endpoint_density_calls": 0,
        "correlation_calls": 0,
        "numeric_values_computed": False,
        "runtime_model_metadata": None,
        "checks": checks,
    }


__all__ = ["run_fixture"]
