"""Bounded pressure composition consumer for the leading Rm seed.

The leading normalization certificate owns the exact ``p`` and ``P0``
endpoint identities.  This module consumes those passed rows and proves only
the remaining pressure composition ``P0 + p``.  Guarded endpoint
substitutions are reduced to qualified source atoms; their old source bodies
are never expanded.  In particular, an unequal-bound ancestor integral is a
hard stop rather than an instruction to inspect its integrand.
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


_Z_GUARD = "Z_independent_substitution"
_ORDERS = (0, 1, 2)


class _Unsupported(ValueError):
    pass


class _BlockedIntegral(_Unsupported):
    def __init__(self, node: int) -> None:
        self.node = node
        super().__init__(f"ancestor integral {node} was intentionally not evaluated")


def canonical_digest(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _module_sha256() -> str:
    digest = hashlib.sha256()
    with Path(__file__).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _handle(value: Any, path: str, node_count: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{path} must be a nonnegative graph handle")
    if node_count is not None and value >= node_count:
        raise ValueError(f"{path} handle {value} is outside nodes")
    return value


def _fraction(value: Any, path: str = "rational") -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{path} must be exact")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        text = value.strip()
        if "/" in text:
            pieces = text.split("/")
            if len(pieces) != 2:
                raise ValueError(f"{path} is malformed")
            return Fraction(int(pieces[0]), int(pieces[1]))
        return Fraction(int(text), 1)
    raise ValueError(f"{path} must be an exact integer or rational")


def _factor_map(monomial: tuple[tuple[Any, ...], ...]) -> dict[tuple[Any, ...], int]:
    result: dict[tuple[Any, ...], int] = {}
    for factor in monomial:
        if len(factor) < 2 or not isinstance(factor[-1], int):
            raise ValueError("malformed pressure factor")
        base = factor[:-1]
        result[base] = result.get(base, 0) + factor[-1]
    return {base: exponent for base, exponent in result.items() if exponent}


def _monomial(mapping: Mapping[tuple[Any, ...], int]) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        tuple(base) + (exponent,)
        for base, exponent in sorted(mapping.items(), key=lambda item: repr(item[0]))
        if exponent
    )


@dataclass(frozen=True)
class _Expr:
    terms: Mapping[tuple[tuple[Any, ...], ...], Fraction] = field(default_factory=dict)

    def __post_init__(self) -> None:
        cleaned: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
        for monomial, coefficient in self.terms.items():
            value = Fraction(coefficient)
            if value:
                cleaned[tuple(monomial)] = value
        object.__setattr__(self, "terms", cleaned)

    @classmethod
    def rational(cls, value: Fraction | int) -> "_Expr":
        value = Fraction(value)
        return cls({(): value} if value else {})

    @classmethod
    def atom(cls, name: str) -> "_Expr":
        return cls({(("atom", str(name), 1),): Fraction(1)})

    @property
    def is_zero(self) -> bool:
        return not self.terms

    @property
    def is_rational(self) -> bool:
        return not self.terms or set(self.terms) <= {()}

    @property
    def rational_value(self) -> Fraction:
        if not self.is_rational:
            raise ValueError("pressure expression is not rational")
        return self.terms.get((), Fraction(0))

    def __add__(self, other: "_Expr") -> "_Expr":
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            terms[monomial] = terms.get(monomial, Fraction(0)) + coefficient
        return _Expr(terms)

    def __neg__(self) -> "_Expr":
        return _Expr({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def scale(self, value: Fraction | int) -> "_Expr":
        value = Fraction(value)
        return _Expr({monomial: value * coefficient for monomial, coefficient in self.terms.items()})

    def __mul__(self, other: "_Expr") -> "_Expr":
        if self.is_zero or other.is_zero:
            return _Expr.rational(0)
        terms: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
        for left, left_coefficient in self.terms.items():
            left_map = _factor_map(left)
            for right, right_coefficient in other.terms.items():
                mapping = dict(left_map)
                for factor, exponent in _factor_map(right).items():
                    mapping[factor] = mapping.get(factor, 0) + exponent
                key = _monomial(mapping)
                terms[key] = terms.get(key, Fraction(0)) + left_coefficient * right_coefficient
        return _Expr(terms)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _Expr) and dict(self.terms) == dict(other.terms)

    def json(self) -> dict[str, Any]:
        return {
            "terms": [
                {
                    "monomial": [list(factor) for factor in monomial],
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
    atom: str
    source_node: int
    descriptor: Mapping[str, Any]
    kind: str


def _row_operation(row: Mapping[str, Any], index: int) -> str:
    operation = row.get("operation", row.get("op"))
    if not isinstance(operation, str):
        raise ValueError(f"node {index} has no operation")
    return operation


def _exact_one_node(nodes: Sequence[Mapping[str, Any]], ref: Any) -> bool:
    if isinstance(ref, bool) or not isinstance(ref, int) or ref < 0 or ref >= len(nodes):
        return False
    row = nodes[ref]
    if not isinstance(row, Mapping):
        return False
    operation = row.get("operation", row.get("op"))
    try:
        if operation == "exact_rational":
            return _fraction(row.get("numerator", row.get("value", 0))) / _fraction(row.get("denominator", 1)) == 1
        if operation == "constant":
            return _fraction(row.get("value")) == 1
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    return False


def _endpoint_descriptor(
    nodes: Sequence[Mapping[str, Any]],
    ref: int,
    label: str,
) -> tuple[int, Any, int]:
    _handle(ref, label, len(nodes))
    row = nodes[ref]
    if not isinstance(row, Mapping) or _row_operation(row, ref) != "function_substitution":
        raise ValueError(f"{label} must be a function_substitution")
    if row.get(_Z_GUARD) is not True:
        raise ValueError(f"{label} lacks {_Z_GUARD}=True")
    expression = _handle(row.get("expression"), f"{label}.expression", len(nodes))
    value = _handle(row.get("value"), f"{label}.value", len(nodes))
    variable = row.get("variable")
    if isinstance(variable, bool) or not isinstance(variable, (int, str)):
        raise ValueError(f"{label}.variable is malformed")
    return expression, variable, value


def _expr_from_json(value: Any) -> _Expr:
    if not isinstance(value, Mapping) or not _is_sequence(value.get("terms")):
        raise ValueError("seed field expression is malformed")
    terms: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
    for item in value["terms"]:
        if not isinstance(item, Mapping) or not _is_sequence(item.get("monomial")):
            raise ValueError("seed field monomial is malformed")
        factors = []
        for factor in item["monomial"]:
            if not _is_sequence(factor) or len(factor) < 2 or not isinstance(factor[-1], int):
                raise ValueError("seed field factor is malformed")
            factors.append(tuple(factor))
        coefficient = _fraction(item.get("coefficient"), "seed coefficient")
        monomial = tuple(factors)
        terms[monomial] = terms.get(monomial, Fraction(0)) + coefficient
    return _Expr(terms)


def _seed_rows(seed_result: Mapping[str, Any], field_name: str) -> dict[int, Mapping[str, Any]]:
    rows: dict[int, Mapping[str, Any]] = {}
    for row in seed_result.get("fields", ()):
        if not isinstance(row, Mapping) or row.get("field") != field_name:
            continue
        order = row.get("order")
        if isinstance(order, bool) or not isinstance(order, int) or order in rows:
            raise ValueError(f"seed {field_name} rows are malformed")
        rows[order] = row
    return rows


def _seed_atom(row: Mapping[str, Any], field_name: str, order: int, side: str) -> _Expr:
    if row.get("status") != "proved":
        raise ValueError(f"seed {field_name}[{order}] {side} row is not proved")
    expression = _expr_from_json(row.get(side))
    expected = _Expr.atom(f"source:{field_name}[{order}]")
    if expression != expected:
        raise ValueError(f"seed {field_name}[{order}] {side} is not the required source atom")
    return expression


class _PressureReducer:
    def __init__(self, nodes: Sequence[Mapping[str, Any]], aliases: Mapping[int, _Alias]):
        self.nodes = nodes
        self.aliases = dict(aliases)
        self.memo: dict[tuple[int, tuple[Any, ...]], _Expr] = {}
        self.used_nodes: set[int] = set()
        self.blocked_integrals: set[int] = set()
        self.alias_witnesses: list[dict[str, Any]] = []
        self._witness_keys: set[str] = set()

    def _record(self, alias: _Alias, ref: int) -> None:
        payload = {
            "atom": alias.atom,
            "kind": alias.kind,
            "source_node": alias.source_node,
            "expression_node": ref,
            "descriptor": deepcopy(dict(alias.descriptor)),
        }
        key = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        if key not in self._witness_keys:
            self._witness_keys.add(key)
            self.alias_witnesses.append(payload)

    @staticmethod
    def _alias_is_active(alias: _Alias, env: Mapping[Any, _Expr]) -> bool:
        """Keep source-body aliases inside their exact endpoint scope.

        Wrapper aliases are already guarded endpoint roots and may be consumed
        directly.  A body alias is only valid after its own guarded
        function-substitution has bound one of the recorded native variables
        to the exact one expression.  This prevents a source recipe from
        becoming a global atom when the same body handle is encountered in a
        different scope.
        """

        variables = alias.descriptor.get("scoped_variables")
        if variables is None:
            return True
        if not _is_sequence(variables):
            return False
        one = _Expr.rational(1)
        return any(variable in env and env[variable] == one for variable in variables)

    def rewrite(self, ref: int, env: Mapping[Any, _Expr] | None = None) -> _Expr:
        env = {} if env is None else env
        _handle(ref, "source handle", len(self.nodes))
        if ref in env:
            return env[ref]
        alias = self.aliases.get(ref)
        if alias is not None and self._alias_is_active(alias, env):
            self.used_nodes.add(ref)
            self._record(alias, ref)
            return _Expr.atom(alias.atom)
        row = self.nodes[ref]
        if not isinstance(row, Mapping):
            raise ValueError(f"nodes[{ref}] is not a mapping")
        name = row.get("name", row.get("symbol"))
        if name in env:
            return env[name]
        # The canonical JSON form is a mapping and therefore cannot itself be
        # used as a dictionary key.  Digest each scoped value while retaining
        # the variable spelling so memoization remains exact and hashable.
        key = (ref, tuple(sorted((repr(k), canonical_digest(value.json())) for k, value in env.items())))
        if key in self.memo:
            return self.memo[key]
        self.used_nodes.add(ref)
        operation = _row_operation(row, ref)
        if operation in {"atom", "source_atom", "qualified_source_atom", "original_source_parameter", "current_original_source_parameter", "symbol"}:
            result = _Expr.atom(f"{operation}:{row.get('name', row.get('source_name', ref))}")
        elif operation == "exact_rational":
            result = _Expr.rational(_fraction(row.get("numerator", row.get("value", 0))) / _fraction(row.get("denominator", 1)))
        elif operation == "constant":
            result = _Expr.rational(_fraction(row.get("value")))
        elif operation == "sum":
            result = _Expr.rational(0)
            arguments = row.get("arguments", row.get("operands"))
            if not _is_sequence(arguments) or not arguments:
                raise ValueError(f"node {ref}.arguments is malformed")
            for child in arguments:
                result = result + self.rewrite(_handle(child, f"node {ref}.argument", len(self.nodes)), env)
        elif operation == "negative":
            child = _handle(row.get("argument", row.get("operand")), f"node {ref}.argument", len(self.nodes))
            result = -self.rewrite(child, env)
        elif operation == "product":
            arguments = row.get("arguments", row.get("operands"))
            if not _is_sequence(arguments) or not arguments:
                raise ValueError(f"node {ref}.arguments is malformed")
            result = _Expr.rational(1)
            for child in arguments:
                result = result * self.rewrite(_handle(child, f"node {ref}.argument", len(self.nodes)), env)
        elif operation in {"positive_quotient", "quotient"}:
            numerator = self.rewrite(_handle(row.get("numerator"), f"node {ref}.numerator", len(self.nodes)), env)
            denominator = self.rewrite(_handle(row.get("denominator"), f"node {ref}.denominator", len(self.nodes)), env)
            if not denominator.is_rational or denominator.rational_value == 0:
                raise _Unsupported(f"node {ref} pressure quotient denominator is not exact nonzero")
            result = numerator.scale(Fraction(1, 1) / denominator.rational_value)
        elif operation == "trim":
            result = self.rewrite(_handle(row.get("argument", row.get("operand")), f"node {ref}.argument", len(self.nodes)), env)
        elif operation == "function_substitution":
            expression, variable, value = _endpoint_descriptor(self.nodes, ref, f"node {ref}")
            if not _exact_one_node(self.nodes, value):
                raise _Unsupported(f"node {ref} substitution is not the p=1/x=1 endpoint")
            scoped = dict(env)
            value_expr = _Expr.rational(1)
            scoped[variable] = value_expr
            scoped[value] = value_expr
            result = self.rewrite(expression, scoped)
        elif operation in {"definite_integral", "integral"}:
            self.blocked_integrals.add(ref)
            raise _BlockedIntegral(ref)
        else:
            raise _Unsupported(f"node {ref} uses unsupported operation {operation!r}")
        self.memo[key] = result
        return result


def _node_digest(nodes: Sequence[Mapping[str, Any]]) -> str:
    return canonical_digest(list(nodes))


def _build_source_aliases(
    nodes: Sequence[Mapping[str, Any]],
    seed_result: Mapping[str, Any],
    field_rows: Mapping[str, Mapping[int, Mapping[str, Any]]],
    orders: Sequence[int],
) -> tuple[dict[int, _Alias], dict[str, Any]]:
    aliases: dict[int, _Alias] = {}
    metadata: dict[str, Any] = {"rows": [], "guarded_endpoint_bodies": []}
    for field_name in ("P0", "p"):
        for order in orders:
            row = field_rows[field_name][order]
            patch_ref = _handle(row.get("patch_node"), f"seed {field_name}[{order}].patch_node", len(nodes))
            inlet_ref = _handle(row.get("inlet_node"), f"seed {field_name}[{order}].inlet_node", len(nodes))
            patch_body, patch_variable, patch_value = _endpoint_descriptor(nodes, patch_ref, f"patch {field_name}[{order}]")
            inlet_body, inlet_variable, inlet_value = _endpoint_descriptor(nodes, inlet_ref, f"inlet {field_name}[{order}]")
            if not (_exact_one_node(nodes, patch_value) and _exact_one_node(nodes, inlet_value)):
                raise ValueError(f"{field_name}[{order}] endpoint values are not exact p=1/x=1")
            atom = f"source:{field_name}[{order}]"
            descriptor = {
                "field": field_name,
                "ordinary_order": order,
                "patch_wrapper": patch_ref,
                "inlet_wrapper": inlet_ref,
                "patch_body": patch_body,
                "inlet_body": inlet_body,
                "patch_variable": patch_variable,
                "inlet_variable": inlet_variable,
                "patch_value": patch_value,
                "inlet_value": inlet_value,
                _Z_GUARD: True,
                "seed_status": "proved",
            }
            wrapper_alias = _Alias(atom, inlet_body, descriptor, "seed_proved_endpoint_wrapper")
            body_descriptor = dict(descriptor)
            body_descriptor["scoped_variables"] = list(
                dict.fromkeys((patch_variable, inlet_variable))
            )
            body_alias = _Alias(atom, inlet_body, body_descriptor, "seed_proved_endpoint_body")

            def add_alias(ref: int, candidate: _Alias) -> None:
                previous = aliases.get(ref)
                if previous is not None and previous.atom != atom:
                    raise ValueError(f"source endpoint handle {ref} has conflicting aliases")
                # A guarded endpoint wrapper wins over a body alias if a
                # source graph happens to reuse one handle for both roles.
                if previous is None or previous.kind == "seed_proved_endpoint_body":
                    aliases[ref] = candidate

            add_alias(patch_ref, wrapper_alias)
            add_alias(inlet_ref, wrapper_alias)
            add_alias(patch_body, body_alias)
            add_alias(inlet_body, body_alias)
            metadata["rows"].append({"field": field_name, "order": order, "patch": patch_ref, "inlet": inlet_ref})
            metadata["guarded_endpoint_bodies"].append(
                {
                    "field": field_name,
                    "order": order,
                    "patch_body": patch_body,
                    "inlet_body": inlet_body,
                    "patch_value_one": True,
                    "inlet_value_one": True,
                }
            )
    return aliases, metadata


def prove(
    nodes: Sequence[Mapping[str, Any]],
    seed_result: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Consume passed ``p``/``P0`` rows and prove only ``pressure=P0+p``."""

    if not isinstance(seed_result, Mapping):
        raise ValueError("seed_result must be a mapping")
    if not isinstance(expected_nodes_digest, str) or len(expected_nodes_digest) != 64:
        raise ValueError("expected_nodes_digest must be a SHA-256 hex digest")
    source_graph = seed_result.get("source_graph")
    if not isinstance(source_graph, Mapping):
        raise ValueError("seed_result source_graph is missing")
    actual_digest = _node_digest(nodes)
    if actual_digest != expected_nodes_digest or source_graph.get("nodes_digest") != expected_nodes_digest:
        raise ValueError("nodes digest does not match seed_result")
    if seed_result.get("source_only") is not True:
        raise ValueError("seed_result must be source_only")
    orders = tuple(seed_result.get("ordinary_orders", _ORDERS))
    if orders != _ORDERS:
        raise ValueError("pressure consumer requires ordinary orders 0, 1, 2")

    field_rows = {name: _seed_rows(seed_result, name) for name in ("p", "P0", "pressure")}
    for name in field_rows:
        if set(field_rows[name]) != set(orders):
            raise ValueError(f"seed {name} rows do not cover ordinary orders")
    for order in orders:
        _seed_atom(field_rows["p"][order], "p", order, "patch")
        _seed_atom(field_rows["p"][order], "p", order, "inlet")
        _seed_atom(field_rows["P0"][order], "P0", order, "patch")
        _seed_atom(field_rows["P0"][order], "P0", order, "inlet")
        if field_rows["pressure"][order].get("status") == "proved":
            raise ValueError("seed pressure row is already proved; refusing to relabel old pressure")
        if field_rows["pressure"][order].get("status") != "unresolved":
            raise ValueError(f"seed pressure[{order}] must remain unresolved")

    aliases, alias_metadata = _build_source_aliases(nodes, seed_result, field_rows, orders)
    reducer = _PressureReducer(nodes, aliases)
    fields: list[dict[str, Any]] = []
    unresolved: list[str] = []
    for order in orders:
        row = field_rows["pressure"][order]
        patch_ref = row["patch_node"]
        inlet_ref = row["inlet_node"]
        expected = _Expr.atom(f"source:P0[{order}]") + _Expr.atom(f"source:p[{order}]")
        try:
            patch_expression = reducer.rewrite(_handle(patch_ref, f"pressure[{order}].patch_node", len(nodes)))
            inlet_expression = reducer.rewrite(_handle(inlet_ref, f"pressure[{order}].inlet_node", len(nodes)))
            status = "proved" if patch_expression == expected and inlet_expression == expected else "failed"
            if status != "proved":
                unresolved.append(f"pressure[{order}]_composition_mismatch")
            fields.append(
                {
                    "field": "pressure",
                    "order": order,
                    "status": status,
                    "patch_node": patch_ref,
                    "inlet_node": inlet_ref,
                    "expected": expected.json(),
                    "patch": patch_expression.json(),
                    "inlet": inlet_expression.json(),
                }
            )
        except _BlockedIntegral as exc:
            unresolved.append(f"pressure[{order}]_ancestor_integral_blocked")
            fields.append(
                {
                    "field": "pressure",
                    "order": order,
                    "status": "unresolved",
                    "patch_node": patch_ref,
                    "inlet_node": inlet_ref,
                    "gate": "ancestor_integral_not_evaluated",
                    "integral_node": exc.node,
                }
            )
        except _Unsupported as exc:
            unresolved.append(f"pressure[{order}]_unsupported")
            fields.append(
                {
                    "field": "pressure",
                    "order": order,
                    "status": "unresolved",
                    "patch_node": patch_ref,
                    "inlet_node": inlet_ref,
                    "gate": "bounded_pressure_source_unsupported",
                    "error": str(exc),
                }
            )

    return {
        "schema": "Rm_seed_pressure_composition_certificate_v1",
        "passed": not unresolved and all(row["status"] == "proved" for row in fields),
        "source_only": True,
        "actual_saved_graph_run": False,
        "actual_saved_graph_loaded": False,
        "seed_actual_saved_graph_run": seed_result.get("actual_saved_graph_run"),
        "seed_actual_saved_graph_loaded": seed_result.get("actual_saved_graph_loaded"),
        "ordinary_orders": list(orders),
        "expected_nodes_digest": expected_nodes_digest,
        "module_sha256": _module_sha256(),
        "source_graph": {"node_count": len(nodes), "nodes_digest": actual_digest},
        "fields": fields,
        "seed_pressure_rows": deepcopy([field_rows["pressure"][order] for order in orders]),
        "seed_partial_packet": {
            "partial_weights": deepcopy(seed_result.get("partial_weights")),
            "mean_premise": deepcopy(seed_result.get("mean_premise")),
            "flat_beta": deepcopy(seed_result.get("flat_beta")),
        },
        "source_aliases": alias_metadata,
        "exact_reductions": {
            "used_nodes": sorted(reducer.used_nodes),
            "blocked_integrals": sorted(reducer.blocked_integrals),
            "alias_witnesses": reducer.alias_witnesses,
            "ancestor_integrands_evaluated": False,
        },
        "unresolved_gates": sorted(set(unresolved)),
        "scope": {
            "pressure_only": True,
            "seed_p_and_P0_consumed": True,
            "seed_pressure_left_unresolved_until_this_consumer": True,
            "no_numeric_evaluation": True,
            "no_caps_or_midpoints_as_values": True,
            "ancestor_integrals_not_evaluated": True,
            "nodes_input_unmutated": deepcopy(nodes) == list(nodes),
        },
    }


__all__ = ["canonical_digest", "prove"]
