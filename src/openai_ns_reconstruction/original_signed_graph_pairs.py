"""Small append-only graph helper for signed ordinary ``Z`` pairs.

The graph owns a private copy of an existing source prefix.  New rows are
appended after that prefix and are also exposed through :attr:`appended` so a
caller can serialize only the locally constructed suffix.  This module only
constructs graph rows; it does not execute source functions, install bounds,
or make domain assumptions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from fractions import Fraction
from typing import Any


Pair = tuple[int, int]


_REQUIRED_REFS: dict[str, tuple[str, ...]] = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "analytic_unary": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "variable", "lower", "upper"),
}


def _fraction(value: Any) -> Fraction:
    """Coerce only the exact scalar forms supported by this helper."""

    if isinstance(value, bool):
        raise TypeError("rational constants do not accept bool")
    if isinstance(value, Fraction):
        return value
    if type(value) is int or isinstance(value, str):
        return Fraction(value)
    raise TypeError("rational constants require Fraction, int, or decimal string")


class GraphPairs:
    """Build append-only scalar graph rows and first-order ordinary ``Z`` jets.

    ``base_nodes`` is copied before any rows are added.  Handles in the base
    are retained exactly, so the first locally created handle is
    ``len(base_nodes)``.  Every handle consumed by a public operation is
    checked against the current graph length; booleans are deliberately not
    accepted as integer handles.
    """

    def __init__(self, base_nodes: Sequence[Mapping[str, Any]]):
        if isinstance(base_nodes, (str, bytes)) or not isinstance(base_nodes, Sequence):
            raise TypeError("base_nodes must be a sequence of node mappings")

        copied: list[dict[str, Any]] = []
        for index, row in enumerate(base_nodes):
            if not isinstance(row, Mapping):
                raise TypeError(f"base node {index} must be a mapping")
            copied.append(deepcopy(dict(row)))

        self.nodes: list[dict[str, Any]] = copied
        self.base_nodes: tuple[dict[str, Any], ...] = tuple(deepcopy(copied))
        self.appended: list[dict[str, Any]] = []
        self._validate_all_rows()

    # ------------------------------------------------------------------
    # Handle and row validation
    # ------------------------------------------------------------------
    @staticmethod
    def _handle(value: Any, *, name: str = "node") -> int:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer graph handle")
        if value < 0:
            raise ValueError(f"{name} must be nonnegative")
        return value

    def _id(self, value: Any, *, name: str = "node") -> int:
        index = self._handle(value, name=name)
        if index >= len(self.nodes):
            raise ValueError(f"{name} {index} is outside the current graph")
        return index

    @staticmethod
    def _required_fields(operation: str) -> tuple[str, ...]:
        return _REQUIRED_REFS.get(operation, ())

    def _row_references(self, row: Mapping[str, Any], *, index: int) -> list[tuple[str, int]]:
        operation = row.get("operation")
        if not isinstance(operation, str):
            raise TypeError(f"node {index} must have a string operation")

        for field in self._required_fields(operation):
            if field not in row:
                raise ValueError(f"node {index} ({operation}) is missing {field}")

        references: list[tuple[str, int]] = []
        for field in self._required_fields(operation):
            value = row[field]
            if field == "arguments":
                if not isinstance(value, (list, tuple)):
                    raise TypeError(f"node {index} arguments must be a list or tuple")
                values = value
            elif operation == "definite_integral" and field == "variable":
                # Existing source graphs use either a graph handle or a
                # lexical variable name for this field.  The local integral
                # helper emits the latter.
                if isinstance(value, str):
                    continue
                values = (value,)
            else:
                values = (value,)

            for offset, child in enumerate(values):
                child_name = f"{field}[{offset}]" if field == "arguments" else field
                child_id = self._handle(child, name=f"node {index} {child_name}")
                references.append((child_name, child_id))
        return references

    def _validate_row(self, row: Mapping[str, Any], *, index: int, limit: int) -> None:
        references = self._row_references(row, index=index)
        for field, child in references:
            if child >= limit:
                raise ValueError(
                    f"node {index} {field} handle {child} is outside existing graph"
                )

    def _validate_all_rows(self) -> None:
        limit = len(self.nodes)
        for index, row in enumerate(self.nodes):
            self._validate_row(row, index=index, limit=limit)

    def _append_row(self, row: Mapping[str, Any]) -> int:
        index = len(self.nodes)
        copied = deepcopy(dict(row))
        self._validate_row(copied, index=index, limit=index)
        self.nodes.append(copied)
        self.appended.append(copied)
        return index

    def node(self, operation: str, **attributes: Any) -> int:
        """Append one graph row and return its new handle."""

        if not isinstance(operation, str):
            raise TypeError("operation must be a string")
        # Keep operation first in normal serialized representations while
        # avoiding mutation of the caller's attribute mapping.
        row = {"operation": operation, **attributes}
        return self._append_row(row)

    # ------------------------------------------------------------------
    # Scalar graph arithmetic
    # ------------------------------------------------------------------
    @staticmethod
    def _rational_row(row: Mapping[str, Any] | None) -> Fraction | None:
        if not row or row.get("operation") != "exact_rational":
            return None
        numerator = row.get("numerator")
        denominator = row.get("denominator")
        if isinstance(numerator, bool) or isinstance(denominator, bool):
            return None
        if not isinstance(numerator, int) or not isinstance(denominator, int) or denominator == 0:
            return None
        return Fraction(numerator, denominator)

    def _rational(self, index: int) -> Fraction | None:
        return self._rational_row(self.nodes[index])

    def _find_constant(self, value: Fraction) -> int | None:
        for index, row in enumerate(self.nodes):
            if self._rational_row(row) == value:
                return index
        return None

    def constant(self, value: Fraction | int | str) -> int:
        """Return an ``exact_rational`` handle with normalized Fraction data."""

        rational = _fraction(value)
        existing = self._find_constant(rational)
        if existing is not None:
            return existing
        return self.node(
            "exact_rational",
            numerator=rational.numerator,
            denominator=rational.denominator,
        )

    def add(self, *ids: int) -> int:
        checked = [self._id(value, name="add operand") for value in ids]
        zero = self.constant(0)
        kept = [value for value in checked if value != zero and self._rational(value) != 0]
        if not kept:
            return zero
        if len(kept) == 1:
            return kept[0]
        return self.node("sum", arguments=kept)

    def mul(self, *ids: int) -> int:
        checked = [self._id(value, name="multiply operand") for value in ids]
        one = self.constant(1)
        if any(value == one or self._rational(value) == 0 for value in checked):
            if any(self._rational(value) == 0 for value in checked):
                return self.constant(0)
        kept = [value for value in checked if value != one and self._rational(value) != 1]
        if not kept:
            return one
        if len(kept) == 1:
            return kept[0]
        return self.node("product", arguments=kept)

    def neg(self, value: int) -> int:
        checked = self._id(value, name="negative operand")
        if self._rational(checked) == 0:
            return self.constant(0)
        return self.node("negative", argument=checked)

    def div(self, numerator: int, denominator: int, certificate: Any) -> int:
        top = self._id(numerator, name="quotient numerator")
        bottom = self._id(denominator, name="quotient denominator")
        if self._rational(bottom) == 0:
            raise ValueError("quotient denominator is an exact zero")
        if self._rational(top) == 0:
            return self.constant(0)
        if self._rational(bottom) == 1:
            return top
        return self.node(
            "positive_quotient",
            numerator=top,
            denominator=bottom,
            source_positive_certificate=certificate,
        )

    # ------------------------------------------------------------------
    # C1 ordinary-Z pair arithmetic
    # ------------------------------------------------------------------
    def _pair(self, pair: Any, *, name: str = "pair") -> Pair:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError(f"{name} must contain exactly ordinary Z0 and Z1")
        return (
            self._id(pair[0], name=f"{name}[0]"),
            self._id(pair[1], name=f"{name}[1]"),
        )

    def fixed(self, value: int) -> Pair:
        checked = self._id(value, name="fixed value")
        return checked, self.constant(0)

    def padd(self, *pairs: Pair) -> Pair:
        checked = [self._pair(pair, name=f"addend {index}") for index, pair in enumerate(pairs)]
        if not checked:
            zero = self.constant(0)
            return zero, zero
        return (
            self.add(*(pair[0] for pair in checked)),
            self.add(*(pair[1] for pair in checked)),
        )

    def pneg(self, pair: Pair) -> Pair:
        value, derivative = self._pair(pair)
        return self.neg(value), self.neg(derivative)

    def pmul(self, *pairs: Pair) -> Pair:
        checked = [self._pair(pair, name=f"factor {index}") for index, pair in enumerate(pairs)]
        if not checked:
            one = self.constant(1)
            return one, self.constant(0)

        value = self.mul(*(pair[0] for pair in checked))
        terms: list[int] = []
        for derivative_index, pair in enumerate(checked):
            factors = [pair[1]]
            factors.extend(other[0] for index, other in enumerate(checked) if index != derivative_index)
            terms.append(self.mul(*factors))
        return value, self.add(*terms)

    def pscale(self, pair: Pair, scalar_id: int) -> Pair:
        value, derivative = self._pair(pair)
        scalar = self._id(scalar_id, name="scale scalar")
        return self.mul(value, scalar), self.mul(derivative, scalar)

    def pdiv(self, numerator: Pair, denominator: Pair, certificate: Any) -> Pair:
        top = self._pair(numerator, name="quotient numerator pair")
        bottom = self._pair(denominator, name="quotient denominator pair")
        value = self.div(top[0], bottom[0], certificate)
        derivative_numerator = self.add(top[1], self.neg(self.mul(value, bottom[1])))
        derivative = self.div(derivative_numerator, bottom[0], certificate)
        return value, derivative

    def at(self, pair: Pair, variable_id: int, value_id: int) -> Pair:
        checked = self._pair(pair)
        variable = self._id(variable_id, name="substitution variable")
        value = self._id(value_id, name="substitution value")
        return tuple(
            self.node(
                "function_substitution",
                expression=component,
                variable=variable,
                value=value,
                Z_independent_substitution=True,
            )
            for component in checked
        )  # type: ignore[return-value]

    def integral(
        self,
        pair: Pair,
        coordinate_id: int,
        lower_id: int,
        upper_id: int,
        name: str,
        metadata: Mapping[str, Any] | None = None,
    ) -> Pair:
        """Build C1 definite integrals over an explicitly supplied integrand.

        The caller owns all physical factors, including any ``R^2`` and
        coordinate Jacobian.  They are therefore left untouched in ``pair``.
        """

        checked = self._pair(pair)
        coordinate = self._id(coordinate_id, name="integral coordinate")
        lower = self._id(lower_id, name="integral lower endpoint")
        upper = self._id(upper_id, name="integral upper endpoint")
        if not isinstance(name, str) or not name or any(character.isspace() for character in name):
            raise ValueError("integral name must be a nonempty whitespace-free string")
        if metadata is None:
            extra: dict[str, Any] = {}
        elif isinstance(metadata, Mapping):
            extra = deepcopy(dict(metadata))
        else:
            raise TypeError("integral metadata must be a mapping or None")

        results: list[int] = []
        for order, integrand in enumerate(checked):
            binder_name = f"{name}_Z{order}"
            binder = self.node("bound_variable", name=binder_name)
            substituted = self.node(
                "function_substitution",
                expression=integrand,
                variable=coordinate,
                value=binder,
                Z_independent_substitution=True,
            )
            row = dict(extra)
            row.update(
                operation="definite_integral",
                integrand=substituted,
                variable=binder_name,
                lower=lower,
                upper=upper,
                ordinary_slow_Z_derivative_order=order,
                endpoint_independent_of_Z=True,
            )
            results.append(self.node("definite_integral", **{key: value for key, value in row.items() if key != "operation"}))
        return results[0], results[1]


__all__ = ["GraphPairs", "Pair"]
