"""Pure per-source-cell assembly for the reserved ``m15`` response.

The caller supplies the already-authenticated interval objects.  This module
only preserves their interval arithmetic and source-cell association:

* ``response`` has three entries;
* each source cell has two normalized theta C0/C1 pairs;
* ``a`` is the C0 amplitude factor and ``ordinary_z1`` is its ordinary Z1
  factor used by the C1 product rule;
* ``m15_pressure`` is the global C0/C1 pressure pair.

For every cell the two theta inputs are formed first, then the exact signed
expression ``-sum(response_i * input_i)`` is evaluated.  A caller-provided
``hull`` is applied only after all cell results have been formed.

No serialized artifacts are read here.  The module does not authenticate a
source family, infer a physical endpoint, select a midpoint, or construct a
P1 term.  Those responsibilities remain with the caller.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Callable


class ShapeError(ValueError):
    """Raised when a supplied response or source-cell record has bad shape."""


_MISSING = object()
_CONTAINER_TYPES = (Mapping, list, tuple, str, bytes, bytearray)


def _atom(value: Any, path: str) -> Any:
    """Accept one arithmetic leaf without coercing or inspecting its center."""

    if isinstance(value, _CONTAINER_TYPES):
        raise ShapeError(f"{path}: expected one interval arithmetic leaf")
    return value


def _fixed_sequence(value: Any, size: int, path: str) -> Sequence[Any]:
    """Require a non-string sequence with exactly ``size`` entries."""

    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise ShapeError(f"{path}: expected a sequence of length {size}")
    if len(value) != size:
        raise ShapeError(f"{path}: expected length {size}, got {len(value)}")
    return value


def _pair(value: Any, path: str) -> tuple[Any, Any]:
    values = _fixed_sequence(value, 2, path)
    return _atom(values[0], f"{path}[0]"), _atom(values[1], f"{path}[1]")


def _triple(value: Any, path: str) -> tuple[Any, Any, Any]:
    values = _fixed_sequence(value, 3, path)
    return (
        _atom(values[0], f"{path}[0]"),
        _atom(values[1], f"{path}[1]"),
        _atom(values[2], f"{path}[2]"),
    )


def _theta_pairs(value: Any, path: str) -> tuple[tuple[Any, Any], tuple[Any, Any]]:
    values = _fixed_sequence(value, 2, path)
    return _pair(values[0], f"{path}[0]"), _pair(values[1], f"{path}[1]")


def _product_rule_pair(
    a: Any,
    ordinary_z1: Any,
    theta: tuple[Any, Any],
) -> tuple[Any, Any]:
    """Return ``(A*theta0, Z1*theta0 + A*theta1)`` for one theta input."""

    theta0, theta1 = theta
    return a * theta0, ordinary_z1 * theta0 + a * theta1


def _signed_response(
    response: tuple[Any, Any, Any],
    input_pairs: tuple[tuple[Any, Any], tuple[Any, Any], tuple[Any, Any]],
) -> tuple[Any, Any]:
    """Evaluate ``-sum(response_i * input_i)`` for C0 and C1."""

    outputs = []
    for coefficient in range(2):
        total = response[0] * input_pairs[0][coefficient]
        total = total + response[1] * input_pairs[1][coefficient]
        total = total + response[2] * input_pairs[2][coefficient]
        outputs.append(-total)
    return outputs[0], outputs[1]


def assemble_source_cell(
    *,
    response: Sequence[Any],
    normalized_theta_inputs: Sequence[Sequence[Any]],
    a: Any,
    ordinary_z1: Any,
    m15_pressure: Sequence[Any],
) -> dict[str, Any]:
    """Assemble one source cell using pure interval arithmetic.

    ``normalized_theta_inputs`` is exactly ``((theta0_C0, theta0_C1),
    (theta1_C0, theta1_C1))``.  ``a`` and ``ordinary_z1`` are scalar interval
    leaves, not midpoint values.  The returned ``inputs_C0_C1`` has the three
    pairs consumed by the signed response: theta-0, theta-1, and global m15
    pressure.
    """

    checked_response = _triple(response, "response")
    checked_theta = _theta_pairs(normalized_theta_inputs, "normalized_theta_inputs")
    checked_a = _atom(a, "a")
    checked_z1 = _atom(ordinary_z1, "ordinary_z1")
    checked_pressure = _pair(m15_pressure, "m15_pressure")
    input_zero = _product_rule_pair(checked_a, checked_z1, checked_theta[0])
    input_one = _product_rule_pair(checked_a, checked_z1, checked_theta[1])
    inputs = (input_zero, input_one, checked_pressure)
    own_swirl = (
        -(checked_response[0] * input_zero[0] + checked_response[1] * input_one[0]),
        -(checked_response[0] * input_zero[1] + checked_response[1] * input_one[1]),
    )
    pressure_term = (
        -(checked_response[2] * checked_pressure[0]),
        -(checked_response[2] * checked_pressure[1]),
    )
    coefficient = _signed_response(checked_response, inputs)
    return dict(
        inputs_C0_C1=inputs,
        own_swirl_C0_C1=own_swirl,
        pressure_C0_C1=pressure_term,
        coefficient_C0_C1=coefficient,
    )


def assemble_source_cells(
    *,
    response: Sequence[Any],
    source_cells: Sequence[Mapping[str, Any]],
    a: Any,
    ordinary_z1: Any,
    m15_pressure: Sequence[Any],
    hull: Callable[[Sequence[Any]], Any],
    error_ledger: Any = _MISSING,
) -> dict[str, Any]:
    """Assemble and finally hull a sequence of reserved-m15 source cells.

    Each source-cell mapping must contain ``cell_id`` and
    ``normalized_theta_inputs``.  If it contains ``error_ledger``, that value
    is copied by reference into the corresponding output record.  The
    top-level ``error_ledger`` argument is likewise returned unchanged when
    supplied.  No ledger is interpreted or merged.
    """

    checked_response = _triple(response, "response")
    checked_a = _atom(a, "a")
    checked_z1 = _atom(ordinary_z1, "ordinary_z1")
    checked_pressure = _pair(m15_pressure, "m15_pressure")
    if isinstance(source_cells, (str, bytes, bytearray)) or not isinstance(source_cells, Sequence):
        raise ShapeError("source_cells: expected a finite sequence")
    if not source_cells:
        raise ShapeError("source_cells: at least one source cell is required")
    if not callable(hull):
        raise ShapeError("hull: expected a callable accepting one value sequence")

    results = []
    for index, cell in enumerate(source_cells):
        path = f"source_cells[{index}]"
        if not isinstance(cell, Mapping):
            raise ShapeError(f"{path}: expected a mapping")
        if "cell_id" not in cell:
            raise ShapeError(f"{path}.cell_id: required")
        if "normalized_theta_inputs" not in cell:
            raise ShapeError(f"{path}.normalized_theta_inputs: required")
        checked_theta = _theta_pairs(cell["normalized_theta_inputs"], f"{path}.normalized_theta_inputs")
        input_zero = _product_rule_pair(checked_a, checked_z1, checked_theta[0])
        input_one = _product_rule_pair(checked_a, checked_z1, checked_theta[1])
        inputs = (input_zero, input_one, checked_pressure)
        own_swirl = (
            -(checked_response[0] * input_zero[0] + checked_response[1] * input_one[0]),
            -(checked_response[0] * input_zero[1] + checked_response[1] * input_one[1]),
        )
        pressure_term = (
            -(checked_response[2] * checked_pressure[0]),
            -(checked_response[2] * checked_pressure[1]),
        )
        coefficient = (own_swirl[0] + pressure_term[0], own_swirl[1] + pressure_term[1])
        result = dict(
            cell_id=cell["cell_id"],
            inputs_C0_C1=inputs,
            own_swirl_C0_C1=own_swirl,
            pressure_C0_C1=pressure_term,
            coefficient_C0_C1=coefficient,
        )
        if "error_ledger" in cell:
            result["error_ledger"] = cell["error_ledger"]
        results.append(result)

    hull_c0 = hull([result["coefficient_C0_C1"][0] for result in results])
    hull_c1 = hull([result["coefficient_C0_C1"][1] for result in results])
    output = dict(cell_results=tuple(results), hull_C0_C1=(hull_c0, hull_c1))
    if error_ledger is not _MISSING:
        output["error_ledger"] = error_ledger
    return output


__all__ = ["ShapeError", "assemble_source_cell", "assemble_source_cells"]
