"""Source-only identity proof for the first native radial jet at ``R_m``.

This module consumes an already materialised graph and two serialized source
binding maps.  It compares the raw ``radial_background`` rows of the
post-restore Q4 chart and the finite-N patch at their native endpoint
substitutions.  Endpoint function substitutions supplied by the caller are
protected as qualified source-function atoms; their bodies are never
expanded.  The reducer is intentionally small and fail-closed.

The proof is a leading raw native-jet statement only.  It does not execute a
saved graph, replay a producer, assign endpoint values, or establish a finite
N/Q_y/R_Q_y germ.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
RADIAL_FIELDS = ("E_y", "V_y")
ALL_FIELDS = FIELDS + RADIAL_FIELDS
ORDERS = (0, 1, 2)
_F_SUB_GUARD = "Z_independent_substitution"
_ALLOWED_OPERATIONS = frozenset(
    {
        "exact_rational",
        "sum",
        "product",
        "negative",
        "function_substitution",
        "bound_variable",
        "positive_quotient",
    }
)


class _ReductionGate(ValueError):
    """A reachable source operation that is outside the bounded reducer."""


def canonical_digest(value: Any) -> str:
    """Return the canonical JSON SHA-256 digest used by the source packets."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _module_sha256() -> str:
    digest = hashlib.sha256()
    with Path(__file__).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _handle(value: Any, path: str, node_count: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{path} must be an integer graph handle")
    if value < 0 or node_count is not None and value >= node_count:
        raise ValueError(f"{path} is outside the source graph")
    return value


def _fraction(value: Any, path: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{path} is not an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str):
        text = value.strip()
        if "/" in text:
            pieces = text.split("/")
            if len(pieces) != 2:
                raise ValueError(f"{path} is a malformed rational")
            return Fraction(int(pieces[0]), int(pieces[1]))
        return Fraction(int(text), 1)
    raise ValueError(f"{path} is not an exact rational")


def _operation(row: Mapping[str, Any], index: int) -> str:
    operation = row.get("operation", row.get("op"))
    if not isinstance(operation, str):
        raise ValueError(f"node {index} has no operation")
    return operation


def _field(row: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    raise KeyError(names[0])


def _args(row: Mapping[str, Any], index: int) -> list[Any]:
    args = row.get("arguments", row.get("operands"))
    if not _is_sequence(args) or not args:
        raise ValueError(f"node {index}.arguments must be a nonempty sequence")
    return list(args)


def _descriptor(row: Mapping[str, Any], index: int) -> tuple[Any, ...]:
    """The guarded substitution identity, excluding unrelated metadata."""

    if _operation(row, index) != "function_substitution":
        raise ValueError(f"node {index} is not a function substitution")
    return (
        row.get("expression"),
        row.get("variable"),
        row.get("value"),
    )


def _descriptor_json(row: Mapping[str, Any], index: int) -> dict[str, Any]:
    expression, variable, value = _descriptor(row, index)
    return {
        "expression": expression,
        "variable": variable,
        "endpoint_node": value,
        _F_SUB_GUARD: row.get(_F_SUB_GUARD) is True,
    }


@dataclass(frozen=True)
class _Poly:
    """A finite exact-rational polynomial in qualified source atoms."""

    terms: Mapping[tuple[str, ...], Fraction] = field(default_factory=dict)

    def __post_init__(self) -> None:
        clean: dict[tuple[str, ...], Fraction] = {}
        for monomial, coefficient in self.terms.items():
            value = Fraction(coefficient)
            if value:
                clean[tuple(str(item) for item in monomial)] = value
        object.__setattr__(self, "terms", clean)

    @classmethod
    def rational(cls, value: Fraction | int) -> "_Poly":
        value = Fraction(value)
        return cls({(): value} if value else {})

    @classmethod
    def atom(cls, name: str) -> "_Poly":
        return cls({(str(name),): Fraction(1)})

    @property
    def is_rational(self) -> bool:
        return not self.terms or set(self.terms) <= {()}

    @property
    def rational_value(self) -> Fraction:
        if not self.is_rational:
            raise ValueError("polynomial is not rational")
        return self.terms.get((), Fraction(0))

    def __add__(self, other: "_Poly") -> "_Poly":
        result = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
        return _Poly(result)

    def __neg__(self) -> "_Poly":
        return _Poly({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __mul__(self, other: "_Poly") -> "_Poly":
        if not self.terms or not other.terms:
            return _Poly.rational(0)
        result: dict[tuple[str, ...], Fraction] = {}
        for left, left_coefficient in self.terms.items():
            for right, right_coefficient in other.terms.items():
                monomial = tuple(sorted(left + right))
                result[monomial] = result.get(monomial, Fraction(0)) + left_coefficient * right_coefficient
        return _Poly(result)

    def scale(self, value: Fraction) -> "_Poly":
        return _Poly({monomial: coefficient * value for monomial, coefficient in self.terms.items()})

    def key(self) -> tuple[Any, ...]:
        return tuple(
            (tuple(monomial), coefficient.numerator, coefficient.denominator)
            for monomial, coefficient in sorted(self.terms.items(), key=lambda item: repr(item[0]))
        )

    def json(self) -> dict[str, Any]:
        return {
            "terms": [
                {
                    "monomial": list(monomial),
                    "coefficient": (
                        str(coefficient.numerator)
                        if coefficient.denominator == 1
                        else f"{coefficient.numerator}/{coefficient.denominator}"
                    ),
                }
                for monomial, coefficient in sorted(self.terms.items(), key=lambda item: repr(item[0]))
            ]
        }


@dataclass(frozen=True)
class _Alias:
    field: str
    order: int
    atom: str
    left_node: int
    right_node: int
    proof_ref: Any
    left_descriptor: tuple[Any, ...]
    right_descriptor: tuple[Any, ...]
    zero_source: bool = False


def _rows(value: Any, label: str) -> list[int]:
    if isinstance(value, Mapping):
        result: list[int] = []
        for order in ORDERS:
            candidate = value.get(str(order), value.get(order))
            if candidate is None:
                raise ValueError(f"{label} lacks ordinary order {order}")
            result.append(candidate)
        return result
    if not _is_sequence(value) or len(value) < len(ORDERS):
        raise ValueError(f"{label} lacks ordinary orders 0, 1, 2")
    return [value[order] for order in ORDERS]


def _as_rows(value: Any, label: str, node_count: int) -> list[int]:
    return [_handle(item, f"{label}[{order}]", node_count) for order, item in enumerate(_rows(value, label))]


def _unwrap(value: Mapping[str, Any], names: Sequence[str], label: str) -> tuple[Mapping[str, Any], str]:
    for name in names:
        nested = value.get(name)
        if isinstance(nested, Mapping):
            return nested, name
    return value, "direct"


def _chart(root: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    charts = root.get("charts")
    if isinstance(charts, Mapping):
        for name in ("postrestore", "post_restore"):
            chart = charts.get(name)
            if isinstance(chart, Mapping):
                return chart
    for name in ("postrestore", "post_restore", "chart"):
        chart = root.get(name)
        if isinstance(chart, Mapping):
            return chart
    raise ValueError(f"{label} postrestore/chart binding is missing")


def _radial_background(root: Mapping[str, Any], chart: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    for owner in (chart, root):
        candidate = owner.get("radial_background")
        if isinstance(candidate, Mapping):
            return candidate
    raise ValueError(f"{label} radial_background binding is missing")


def _endpoint_variable_candidates(root: Mapping[str, Any], chart: Mapping[str, Any], side: str) -> list[int]:
    keys = (
        ("native_postrestore_coordinate", "native_endpoint_variable", "coordinate", "phase")
        if side == "left"
        else ("native_endpoint_variable", "coordinate", "x_variable")
    )
    candidates: list[int] = []
    owners: list[Mapping[str, Any]] = [chart, root]
    endpoint = root.get("endpoints")
    if isinstance(endpoint, Mapping):
        rm = endpoint.get("Rm")
        if isinstance(rm, Mapping):
            owners.append(rm)
            metadata = rm.get("native_endpoint")
            if isinstance(metadata, Mapping):
                owners.append(metadata)
    for owner in owners:
        for key in keys:
            item = owner.get(key)
            if isinstance(item, int) and not isinstance(item, bool):
                candidates.append(item)
    return sorted(set(candidates))


def _exact_one(nodes: Sequence[Mapping[str, Any]], index: int) -> bool:
    try:
        row = nodes[index]
        if not isinstance(row, Mapping) or _operation(row, index) != "exact_rational":
            return False
        numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
        denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
        return denominator != 0 and numerator / denominator == 1
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def _exact_zero(nodes: Sequence[Mapping[str, Any]], index: int) -> bool:
    try:
        row = nodes[index]
        if not isinstance(row, Mapping) or _operation(row, index) != "exact_rational":
            return False
        numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
        denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
        return denominator != 0 and numerator == 0
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def _child_handles(row: Mapping[str, Any], index: int, node_count: int) -> list[tuple[str, int]]:
    operation = _operation(row, index)
    if operation not in _ALLOWED_OPERATIONS:
        raise _ReductionGate(f"node {index} uses unsupported operation {operation!r}")
    values: list[tuple[str, Any]] = []
    if operation in {"sum", "product"}:
        values.extend((f"arguments[{offset}]", item) for offset, item in enumerate(_args(row, index)))
    elif operation == "negative":
        values.append(("argument", _field(row, "argument", "operand")))
    elif operation == "positive_quotient":
        values.extend((("numerator", _field(row, "numerator")), ("denominator", _field(row, "denominator"))))
    elif operation == "function_substitution":
        values.extend(
            (
                ("expression", _field(row, "expression")),
                ("variable", _field(row, "variable")),
                ("value", _field(row, "value")),
            )
        )
    result: list[tuple[str, int]] = []
    for path, item in values:
        result.append((path, _handle(item, f"node {index}.{path}", node_count)))
    return result


class _Evaluator:
    def __init__(
        self,
        nodes: Sequence[Mapping[str, Any]],
        aliases: Mapping[int, _Alias],
        expression_aliases: Mapping[int, _Alias],
        descriptor_aliases: Mapping[tuple[Any, ...], _Alias],
        endpoint_variable: int,
        side: str,
    ) -> None:
        self.nodes = nodes
        self.aliases = dict(aliases)
        self.expression_aliases = dict(expression_aliases)
        self.descriptor_aliases = dict(descriptor_aliases)
        self.endpoint_variable = endpoint_variable
        self.side = side
        self.memo: dict[tuple[int, tuple[Any, ...]], _Poly] = {}
        self.used_nodes: set[int] = set()
        self.used_aliases: list[_Alias] = []
        self._used_alias_keys: set[tuple[str, int]] = set()

    def _row(self, index: int) -> Mapping[str, Any]:
        _handle(index, "source handle", len(self.nodes))
        row = self.nodes[index]
        if not isinstance(row, Mapping):
            raise ValueError(f"nodes[{index}] must be a mapping")
        return row

    def _remember_alias(self, alias: _Alias) -> _Poly:
        key = (alias.field, alias.order)
        if key not in self._used_alias_keys:
            self._used_alias_keys.add(key)
            self.used_aliases.append(alias)
        return _Poly.rational(0) if alias.zero_source else _Poly.atom(alias.atom)

    @staticmethod
    def _env_signature(env: Mapping[int, _Poly]) -> tuple[Any, ...]:
        return tuple((int(key), value.key()) for key, value in sorted(env.items()))

    def rewrite(self, index: int, env: Mapping[int, _Poly] | None = None) -> _Poly:
        _handle(index, "source handle", len(self.nodes))
        environment = dict(env or {})
        memo_key = (index, self._env_signature(environment))
        if memo_key in self.memo:
            return self.memo[memo_key]
        if index in environment:
            result = environment[index]
            self.memo[memo_key] = result
            return result
        alias = self.aliases.get(index)
        if alias is not None:
            self.used_nodes.add(index)
            result = self._remember_alias(alias)
            self.memo[memo_key] = result
            return result
        expression_alias = self.expression_aliases.get(index)
        if expression_alias is not None and self.endpoint_variable in environment:
            endpoint = environment[self.endpoint_variable]
            if endpoint.is_rational and endpoint.rational_value == 1:
                self.used_nodes.add(index)
                result = self._remember_alias(expression_alias)
                self.memo[memo_key] = result
                return result
        row = self._row(index)
        operation = _operation(row, index)
        if operation == "function_substitution":
            if row.get(_F_SUB_GUARD) is not True:
                raise _ReductionGate(f"node {index} lacks {_F_SUB_GUARD}=True")
            descriptor_alias = self.descriptor_aliases.get(_descriptor(row, index))
            if descriptor_alias is not None:
                self.used_nodes.add(index)
                result = self._remember_alias(descriptor_alias)
                self.memo[memo_key] = result
                return result
        for path, child in _child_handles(row, index, len(self.nodes)):
            if child >= index:
                raise _ReductionGate(f"forward reference from node {index} {path} to {child}")
        self.used_nodes.add(index)
        if operation == "exact_rational":
            numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
            denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
            if denominator == 0:
                raise _ReductionGate(f"node {index} has a zero rational denominator")
            result = _Poly.rational(numerator / denominator)
        elif operation == "bound_variable":
            name = row.get("name", row.get("symbol", f"bound[{index}]"))
            result = _Poly.atom(f"BOUND:{name}")
        elif operation == "sum":
            result = _Poly.rational(0)
            for child in _args(row, index):
                result = result + self.rewrite(_handle(child, f"node {index}.argument", len(self.nodes)), environment)
        elif operation == "product":
            result = _Poly.rational(1)
            for child in _args(row, index):
                result = result * self.rewrite(_handle(child, f"node {index}.argument", len(self.nodes)), environment)
        elif operation == "negative":
            child = _handle(_field(row, "argument", "operand"), f"node {index}.argument", len(self.nodes))
            result = -self.rewrite(child, environment)
        elif operation == "positive_quotient":
            numerator = self.rewrite(_handle(_field(row, "numerator"), f"node {index}.numerator", len(self.nodes)), environment)
            denominator = self.rewrite(_handle(_field(row, "denominator"), f"node {index}.denominator", len(self.nodes)), environment)
            if not denominator.is_rational or denominator.rational_value == 0:
                raise _ReductionGate(f"node {index} quotient denominator is not an exact nonzero rational")
            result = numerator.scale(Fraction(1, 1) / denominator.rational_value)
        elif operation == "function_substitution":
            expression = _handle(_field(row, "expression"), f"node {index}.expression", len(self.nodes))
            variable = _handle(_field(row, "variable"), f"node {index}.variable", len(self.nodes))
            value = _handle(_field(row, "value"), f"node {index}.value", len(self.nodes))
            replacement = self.rewrite(value, environment)
            scoped = dict(environment)
            scoped[variable] = replacement
            result = self.rewrite(expression, scoped)
        else:  # pragma: no cover - guarded by _child_handles
            raise _ReductionGate(f"node {index} uses unsupported operation {operation!r}")
        self.memo[memo_key] = result
        return result


def _atom_name(field: str, order: int) -> str:
    return f"FUNCTION:{field}[{order}]"


def _validate_hex_digest(value: Any, name: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"{name} must be a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueError(f"{name} must be a SHA-256 hex digest") from exc
    return value.lower()


def _register_aliases(
    nodes: Sequence[Mapping[str, Any]],
    seed_aliases: Sequence[Mapping[str, Any]],
) -> tuple[
    dict[int, _Alias],
    dict[int, _Alias],
    dict[int, _Alias],
    dict[int, _Alias],
    dict[tuple[Any, ...], _Alias],
    dict[str, Any],
]:
    if not _is_sequence(seed_aliases):
        raise ValueError("seed_aliases must be a sequence")
    expected = {(field, order) for field in ALL_FIELDS for order in ORDERS}
    seen: set[tuple[str, int]] = set()
    left_by_node: dict[int, _Alias] = {}
    right_by_node: dict[int, _Alias] = {}
    left_expression_aliases: dict[int, _Alias] = {}
    right_expression_aliases: dict[int, _Alias] = {}
    descriptor_aliases: dict[tuple[Any, ...], _Alias] = {}
    records: list[dict[str, Any]] = []
    left_variables: set[int] = set()
    right_variables: set[int] = set()
    for offset, item in enumerate(seed_aliases):
        if not isinstance(item, Mapping):
            raise ValueError(f"seed_aliases[{offset}] must be a mapping")
        field_name = item.get("field")
        order = item.get("order")
        if field_name not in ALL_FIELDS or isinstance(order, bool) or not isinstance(order, int) or order not in ORDERS:
            raise ValueError(f"seed_aliases[{offset}] has an invalid field/order")
        key = (field_name, order)
        if key in seen:
            raise ValueError(f"duplicate seed alias {field_name}[{order}]")
        seen.add(key)
        proof_ref = item.get("proof_ref")
        if proof_ref is None or proof_ref == "":
            raise ValueError(f"seed_aliases[{offset}] lacks proof_ref")
        left_node = _handle(item.get("left_node"), f"seed_aliases[{offset}].left_node", len(nodes))
        right_node = _handle(item.get("right_node"), f"seed_aliases[{offset}].right_node", len(nodes))
        left_row = nodes[left_node]
        right_row = nodes[right_node]
        if not isinstance(left_row, Mapping) or not isinstance(right_row, Mapping):
            raise ValueError(f"seed alias {field_name}[{order}] does not point to graph rows")
        if _operation(left_row, left_node) != "function_substitution" or _operation(right_row, right_node) != "function_substitution":
            raise ValueError(f"seed alias {field_name}[{order}] is not guarded by function substitution")
        for side, row, index, variables in (
            ("left", left_row, left_node, left_variables),
            ("right", right_row, right_node, right_variables),
        ):
            if row.get(_F_SUB_GUARD) is not True:
                raise ValueError(f"{side} seed alias {field_name}[{order}] lacks {_F_SUB_GUARD}=True")
            expression, variable, value = _descriptor(row, index)
            _handle(expression, f"{side} seed alias expression", len(nodes))
            variable_handle = _handle(variable, f"{side} seed alias variable", len(nodes))
            value_handle = _handle(value, f"{side} seed alias endpoint", len(nodes))
            if expression >= index or variable_handle >= index or value_handle >= index:
                raise ValueError(f"{side} seed alias {field_name}[{order}] has a forward reference")
            if not _exact_one(nodes, value_handle):
                raise ValueError(f"{side} seed alias {field_name}[{order}] endpoint is not exact rational one")
            variables.add(variable_handle)
        left_descriptor = _descriptor(left_row, left_node)
        right_descriptor = _descriptor(right_row, right_node)
        atom = _atom_name(field_name, order)
        zero_source = _exact_zero(nodes, left_descriptor[0]) or _exact_zero(nodes, right_descriptor[0])
        alias = _Alias(
            field_name,
            order,
            atom,
            left_node,
            right_node,
            deepcopy(proof_ref),
            left_descriptor,
            right_descriptor,
            zero_source,
        )
        if left_node in left_by_node or right_node in right_by_node:
            raise ValueError(f"seed alias node reused at {field_name}[{order}]")
        left_by_node[left_node] = alias
        right_by_node[right_node] = alias
        for expression, side, destination in (
            (left_descriptor[0], "left", left_expression_aliases),
            (right_descriptor[0], "right", right_expression_aliases),
        ):
            expression_handle = _handle(expression, f"{side} source expression", len(nodes))
            previous_expression = destination.get(expression_handle)
            if previous_expression is not None and not (previous_expression.zero_source and zero_source):
                raise ValueError(
                    f"source expression aliases two nonzero fields: {side} {field_name}[{order}]"
                )
            destination[expression_handle] = alias
        for descriptor, side in ((left_descriptor, "left"), (right_descriptor, "right")):
            previous = descriptor_aliases.get(descriptor)
            if previous is not None and (previous.field, previous.order) != key and not (previous.zero_source and zero_source):
                raise ValueError(f"guarded substitution descriptor aliases two source fields: {side} {field_name}[{order}]")
            descriptor_aliases[descriptor] = alias
        records.append(
            {
                "field": field_name,
                "order": order,
                "atom": atom,
                "zero_source": zero_source,
                "left_node": left_node,
                "right_node": right_node,
                "proof_ref": deepcopy(proof_ref),
                "left_guard": {
                    "expression": left_descriptor[0],
                    "variable": left_descriptor[1],
                    "endpoint_node": left_descriptor[2],
                    _F_SUB_GUARD: True,
                },
                "right_guard": {
                    "expression": right_descriptor[0],
                    "variable": right_descriptor[1],
                    "endpoint_node": right_descriptor[2],
                    _F_SUB_GUARD: True,
                },
            }
        )
    missing = sorted(expected - seen)
    if missing:
        raise ValueError(f"seed_aliases lacks fields/orders: {missing}")
    if len(left_variables) != 1 or len(right_variables) != 1:
        raise ValueError("seed aliases do not identify one native endpoint variable per side")
    return (
        left_by_node,
        right_by_node,
        left_expression_aliases,
        right_expression_aliases,
        descriptor_aliases,
        {
            "rows": records,
            "left_endpoint_variable": next(iter(left_variables)),
            "right_endpoint_variable": next(iter(right_variables)),
        },
    )


def _field_rows(radial_background: Mapping[str, Any], field_name: str, label: str, node_count: int) -> list[int]:
    if field_name not in radial_background:
        raise ValueError(f"{label} lacks radial_background.{field_name}")
    return _as_rows(radial_background[field_name], f"{label}.radial_background.{field_name}", node_count)


def _proof_ref_list(aliases: Sequence[_Alias]) -> list[Any]:
    seen: set[str] = set()
    result: list[Any] = []
    for alias in aliases:
        marker = json.dumps(alias.proof_ref, sort_keys=True, separators=(",", ":"), default=str)
        if marker not in seen:
            seen.add(marker)
            result.append(deepcopy(alias.proof_ref))
    return result


def _row_result(
    *,
    field_name: str,
    order: int,
    left_node: int,
    right_node: int,
    left: _Poly,
    right: _Poly,
    aliases: Sequence[_Alias],
) -> dict[str, Any]:
    left_json = left.json()
    right_json = right.json()
    return {
        "field": field_name,
        "order": order,
        "left_node": left_node,
        "right_node": right_node,
        "status": "proved" if left == right else "failed",
        "left_form": left_json,
        "right_form": right_json,
        "left_form_sha256": canonical_digest(left_json),
        "right_form_sha256": canonical_digest(right_json),
        "proof_refs": _proof_ref_list(aliases),
    }


def prove(
    nodes: Sequence[Mapping[str, Any]],
    modified_bindings: Mapping[str, Any],
    finite_bindings: Mapping[str, Any],
    seed_aliases: Sequence[Mapping[str, Any]],
    *,
    expected_nodes_digest: str,
    seed_proof_digest: str,
) -> dict[str, Any]:
    """Prove the leading raw native-y jet through order 2.

    ``modified_bindings`` and ``finite_bindings`` may be explicit archive
    wrappers.  The function never imports their producers, executes a graph,
    or dereferences ``proof_ref``.  Structural input errors raise; reachable
    unsupported source operations are returned as unresolved row gates.
    """

    if not _is_sequence(nodes):
        raise ValueError("nodes must be a sequence")
    _validate_hex_digest(expected_nodes_digest, "expected_nodes_digest")
    seed_proof_digest = _validate_hex_digest(seed_proof_digest, "seed_proof_digest")
    expected_digest = expected_nodes_digest.lower()
    actual_digest = canonical_digest(list(nodes))
    if actual_digest != expected_digest:
        raise ValueError("nodes digest does not match expected_nodes_digest")
    if not isinstance(modified_bindings, Mapping) or not isinstance(finite_bindings, Mapping):
        raise ValueError("modified_bindings and finite_bindings must be mappings")
    original_nodes = deepcopy(nodes)
    original_modified = deepcopy(modified_bindings)
    original_finite = deepcopy(finite_bindings)

    modified, modified_wrapper = _unwrap(
        modified_bindings,
        ("modified_Q4_pressure_C3_bindings", "modified_bindings", "q4", "postrestore"),
        "modified",
    )
    finite, finite_wrapper = _unwrap(
        finite_bindings,
        ("Rm_patch_finite_N_C5_bindings", "finite_bindings", "finite_patch", "finite"),
        "finite",
    )
    modified_chart = _chart(modified, "modified")
    finite_chart = finite.get("chart") if isinstance(finite.get("chart"), Mapping) else finite
    if not isinstance(finite_chart, Mapping):
        raise ValueError("finite chart binding is missing")
    modified_background = _radial_background(modified, modified_chart, "modified")
    finite_background = _radial_background(finite, finite_chart, "finite")
    (
        left_by_node,
        right_by_node,
        left_expression_aliases,
        right_expression_aliases,
        descriptor_aliases,
        alias_info,
    ) = _register_aliases(nodes, seed_aliases)
    left_variable = alias_info["left_endpoint_variable"]
    right_variable = alias_info["right_endpoint_variable"]
    for side, root, chart, variable in (
        ("left", modified, modified_chart, left_variable),
        ("right", finite, finite_chart, right_variable),
    ):
        declared = _endpoint_variable_candidates(root, chart, side)
        if declared and any(item != variable for item in declared):
            raise ValueError(f"{side} binding endpoint variable conflicts with guarded seed aliases")

    left_rows = {name: _field_rows(modified_background, name, "modified", len(nodes)) for name in FIELDS}
    right_rows = {name: _field_rows(finite_background, name, "finite", len(nodes)) for name in FIELDS}
    left_eval = _Evaluator(
        nodes,
        left_by_node,
        left_expression_aliases,
        descriptor_aliases,
        left_variable,
        "left",
    )
    right_eval = _Evaluator(
        nodes,
        right_by_node,
        right_expression_aliases,
        descriptor_aliases,
        right_variable,
        "right",
    )
    row_forms: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for field_name in FIELDS:
        for order in ORDERS:
            left_node = left_rows[field_name][order]
            right_node = right_rows[field_name][order]
            try:
                left_form = left_eval.rewrite(left_node, {left_variable: _Poly.rational(1)})
                right_form = right_eval.rewrite(right_node, {right_variable: _Poly.rational(1)})
            except (ValueError, _ReductionGate) as exc:
                unresolved.append(
                    {
                        "field": field_name,
                        "order": order,
                        "left_node": left_node,
                        "right_node": right_node,
                        "gate": "source_reduction_unresolved",
                        "error": str(exc),
                    }
                )
                continue
            row_forms.append(
                _row_result(
                    field_name=field_name,
                    order=order,
                    left_node=left_node,
                    right_node=right_node,
                    left=left_form,
                    right=right_form,
                    aliases=tuple(left_eval.used_aliases) + tuple(right_eval.used_aliases),
                )
            )
    failed = [row for row in row_forms if row["status"] != "proved"]
    unresolved_gates = list(unresolved)
    for row in failed:
        unresolved_gates.append(
            {
                "field": row["field"],
                "order": row["order"],
                "gate": "left_right_polynomial_mismatch",
            }
        )
    all_rows_proved = len(row_forms) == len(FIELDS) * len(ORDERS) and not failed
    result = {
        "schema": "Rm_raw_native_jet_identity_v1",
        "passed": bool(all_rows_proved and not unresolved_gates),
        "source_only": True,
        "actual_saved_graph_run": False,
        "actual_saved_graph_loaded": False,
        "leading_raw_native_jet_only": True,
        "modified_finite_N_germ": False,
        "RQy_identity_proved": False,
        "expected_nodes_digest": expected_digest,
        "seed_proof_digest": seed_proof_digest,
        "module_sha256": _module_sha256(),
        "source_graph": {
            "node_count": len(nodes),
            "nodes_digest": actual_digest,
            "modified_wrapper": modified_wrapper,
            "finite_wrapper": finite_wrapper,
            "postrestore_endpoint_variable": left_variable,
            "finite_endpoint_variable": right_variable,
            "endpoint_node_is_exact_rational_one": True,
        },
        "seed_aliases": alias_info["rows"],
        "row_forms": row_forms,
        "unresolved_gates": unresolved_gates,
        "scope": {
            "source_graph_only": True,
            "no_endpoint_function_bodies_expanded": True,
            "no_profile_or_loop_ancestor_expansion": True,
            "no_numeric_endpoint_assignment": True,
            "no_controls_or_caps_consumed": True,
            "metadata_preserved": (
                deepcopy(modified_bindings) == original_modified
                and deepcopy(finite_bindings) == original_finite
            ),
            "nodes_input_unmutated": deepcopy(nodes) == original_nodes,
        },
    }
    return result


__all__ = [
    "ALL_FIELDS",
    "FIELDS",
    "ORDERS",
    "RADIAL_FIELDS",
    "canonical_digest",
    "prove",
]
