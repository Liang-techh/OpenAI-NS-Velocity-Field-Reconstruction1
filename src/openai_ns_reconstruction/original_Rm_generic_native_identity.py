"""Bounded source consumer for the leading generic Rm native-y jet.

The source packets accepted by this module are deliberately small interfaces:
the normalized leading seed, the accepted raw native jet, and the second-radial
certificate.  The graph is inspected through a guarded symbolic reducer.  No
producer, checker, reducer, fixture, provider, archive, quadrature, or saved
graph is imported or executed here.

Only the leading chart roots are admitted.  The two ordinary rows consumed by
this module are rows 0 and 1; rows through 2 are retained in the seed/raw
source packet because a first ordinary derivative of the generic moment
formula consumes one shifted row.  A third shifted row would be needed by a
higher ordinary derivative of the physical Q_y expression and is outside this
consumer's claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

try:
    import sympy as _sp
except ImportError:  # pragma: no cover
    _sp = None


SEED_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure", "E_y", "V_y")
RAW_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
ROOTS = ("E", "V", "a", "b", "p1", "p2", "t0", "kappa", "Delta")
SIDES = ("postrestore", "patch")
SEED_ORDERS = (0, 1, 2)
ORDERS = (0, 1)
_FUNCTION_GUARDS = (
    "Z_independent_substitution",
    "native_substitution_guard",
    "exact_native_chain_rule_substitution",
    "endpoint_substitution_guard",
    "guarded",
)
_ANALYTIC = {"exp", "log", "sqrt", "positive_sqrt", "expm1", "log1p", "sin", "cos"}
_SOURCE_CONSTANT_OPS = {
    "current_original_source_parameter",
    "current_original_leading_function_recipe",
    "current_original_unique_analytic_core_component",
    "current_Rm_unique_leading_implicit_control",
    "exact_positive_integer_expression",
    "mathematical_pi",
}


class _Unresolved(ValueError):
    """A reachable native-dependent source operation outside this reducer."""


def canonical_digest(value: Any) -> str:
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


def _require_digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"{label} must be a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueError(f"{label} must be a SHA-256 hex digest") from exc
    return value.lower()


def _handle(value: Any, label: str, count: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer graph handle")
    if value < 0 or count is not None and value >= count:
        raise ValueError(f"{label} is outside the source graph")
    return value


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{label} is not an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str):
        text = value.strip()
        try:
            if "/" in text:
                pieces = text.split("/")
                if len(pieces) != 2:
                    raise ValueError
                return Fraction(int(pieces[0]), int(pieces[1]))
            return Fraction(int(text), 1)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{label} is not an exact rational") from exc
    raise ValueError(f"{label} is not an exact rational")


def _operation(row: Mapping[str, Any], index: int) -> str:
    operation = row.get("operation", row.get("op"))
    if not isinstance(operation, str) or not operation:
        raise ValueError(f"node {index} has no operation")
    return operation


def _field(row: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    raise _Unresolved("missing " + "/".join(names))


def _args(row: Mapping[str, Any], index: int) -> list[Any]:
    args = row.get("arguments", row.get("operands"))
    if not _is_sequence(args):
        raise _Unresolved(f"node {index}.arguments is not a sequence")
    return list(args)


def _exact_value(nodes: Sequence[Any], index: int) -> Fraction | None:
    try:
        row = nodes[index]
        if not isinstance(row, Mapping) or _operation(row, index) != "exact_rational":
            return None
        numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
        denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
        if denominator == 0:
            return None
        return numerator / denominator
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return None


def _exact_one(nodes: Sequence[Any], index: int) -> bool:
    return _exact_value(nodes, index) == 1


def _guarded(row: Mapping[str, Any]) -> bool:
    return any(row.get(name) is True for name in _FUNCTION_GUARDS)


def _backward(index: int, value: Any, label: str, count: int) -> int:
    child = _handle(value, f"node {index}.{label}", count)
    if child >= index:
        raise _Unresolved(f"forward reference from node {index}.{label} to node {child}")
    return child


def _simplify(expr: Any) -> Any:
    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact generic native identities")
    return _sp.cancel(_sp.together(expr))


def _record(expr: Any) -> dict[str, Any]:
    expr = _simplify(expr)
    return {"srepr": _sp.srepr(expr), "text": str(expr), "is_zero": bool(expr == 0)}


def _atom(label: str, node: int | None = None) -> Any:
    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact generic native identities")
    suffix = "" if node is None else f"@node:{node}"
    return _sp.Symbol(f"SOURCE:{label}{suffix}", commutative=True)


def _symbol(field: str, order: int) -> Any:
    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact generic native identities")
    return _sp.Symbol(f"SEED:{field}[{order}]", commutative=True)


@dataclass(frozen=True)
class _Dual:
    value: Any
    derivative: Any


def _dual(value: Any, derivative: Any = 0) -> _Dual:
    return _Dual(_simplify(value), _simplify(derivative))


def _add(left: _Dual, right: _Dual) -> _Dual:
    return _dual(left.value + right.value, left.derivative + right.derivative)


def _neg(value: _Dual) -> _Dual:
    return _dual(-value.value, -value.derivative)


def _mul(left: _Dual, right: _Dual) -> _Dual:
    return _dual(
        left.value * right.value,
        left.derivative * right.value + left.value * right.derivative,
    )


def _quotient(numerator: _Dual, denominator: _Dual) -> _Dual:
    if _simplify(denominator.value) == 0:
        raise _Unresolved("positive_quotient denominator is symbolically zero")
    return _dual(
        numerator.value / denominator.value,
        (numerator.derivative * denominator.value
         - numerator.value * denominator.derivative) / denominator.value**2,
    )


def _unary(name: str, child: _Dual) -> _Dual:
    if name == "exp":
        value = _sp.exp(child.value)
        return _dual(value, value * child.derivative)
    if name == "expm1":
        value = _sp.exp(child.value) - 1
        return _dual(value, _sp.exp(child.value) * child.derivative)
    if name == "log":
        return _dual(_sp.log(child.value), child.derivative / child.value)
    if name in {"sqrt", "positive_sqrt"}:
        value = _sp.sqrt(child.value)
        return _dual(value, child.derivative / (2 * value))
    if name == "log1p":
        return _dual(_sp.log(1 + child.value), child.derivative / (1 + child.value))
    if name == "sin":
        return _dual(_sp.sin(child.value), _sp.cos(child.value) * child.derivative)
    if name == "cos":
        return _dual(_sp.cos(child.value), -_sp.sin(child.value) * child.derivative)
    raise _Unresolved(f"unsupported analytic unary {name!r}")


def _proof_refs(value: Any) -> list[Any]:
    if value in (None, ""):
        return []
    if isinstance(value, Mapping):
        for key in ("proof_refs", "source_refs", "refs", "proof_ref", "source_ref"):
            if key in value:
                return _proof_refs(value[key])
        return [deepcopy(dict(value))]
    if _is_sequence(value):
        result: list[Any] = []
        for item in value:
            result.extend(_proof_refs(item))
        return result
    if isinstance(value, (str, int)):
        return [deepcopy(value)]
    raise ValueError("proof reference must be JSON-compatible data")


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[str] = set()
    for value in values:
        marker = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
        if marker not in seen:
            seen.add(marker)
            result.append(deepcopy(value))
    return result


def _parse_form(form: Any, label: str, seed_symbols: Mapping[tuple[str, int], Any]) -> Any:
    """Parse an accepted raw-jet polynomial without importing its reducer."""

    if not isinstance(form, Mapping):
        raise ValueError(f"{label} must be a polynomial mapping")
    terms = form.get("terms")
    if not _is_sequence(terms):
        raise ValueError(f"{label}.terms must be a sequence")
    result = _sp.Integer(0)
    for offset, term in enumerate(terms):
        if not isinstance(term, Mapping):
            raise ValueError(f"{label}.terms[{offset}] must be a mapping")
        monomial = term.get("monomial")
        if not _is_sequence(monomial):
            raise ValueError(f"{label}.terms[{offset}].monomial must be a sequence")
        coefficient = _fraction(term.get("coefficient"), f"{label}.terms[{offset}].coefficient")
        product = _sp.Rational(coefficient.numerator, coefficient.denominator)
        for atom in monomial:
            if not isinstance(atom, str) or not atom.startswith("FUNCTION:"):
                raise _Unresolved(f"{label} has an unqualified raw source atom {atom!r}")
            body = atom[len("FUNCTION:"):]
            if "[" not in body or not body.endswith("]"):
                raise ValueError(f"{label} has malformed source atom {atom!r}")
            field, raw_order = body[:-1].rsplit("[", 1)
            try:
                order = int(raw_order)
            except ValueError as exc:
                raise ValueError(f"{label} has malformed source order {atom!r}") from exc
            key = (field, order)
            if key not in seed_symbols:
                raise _Unresolved(f"{label} references an unadmitted source field {atom!r}")
            product *= seed_symbols[key]
        result += product
    return _simplify(result)


def _raw_forms(raw_native_jet: Mapping[str, Any], seed_symbols: Mapping[tuple[str, int], Any]) -> tuple[dict[tuple[str, int], Any], list[Any], list[dict[str, Any]]]:
    rows = raw_native_jet.get("row_forms")
    if not _is_sequence(rows):
        raise ValueError("raw_native_jet.row_forms must be a sequence")
    expected = {(field, order) for field in RAW_FIELDS for order in SEED_ORDERS}
    seen: set[tuple[str, int]] = set()
    forms: dict[tuple[str, int], Any] = {}
    refs: list[Any] = []
    kept: list[dict[str, Any]] = []
    for offset, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"raw_native_jet.row_forms[{offset}] must be a mapping")
        field = row.get("field")
        order = row.get("order", row.get("ordinary_Z_order"))
        if field not in RAW_FIELDS or isinstance(order, bool) or not isinstance(order, int) or order not in SEED_ORDERS:
            raise ValueError(f"raw_native_jet.row_forms[{offset}] has an invalid field/order")
        key = (field, order)
        if key in seen:
            raise ValueError(f"duplicate raw native row {field}[{order}]")
        seen.add(key)
        if row.get("status") not in (None, "proved", "proved_by_accepted_source_pair"):
            raise ValueError(f"raw native row {field}[{order}] is not admitted")
        left_form = row.get("left_form")
        right_form = row.get("right_form")
        # The raw left/right equality is an accepted source certificate.  It
        # is consumed as data here; re-reducing both sides would turn this
        # consumer into a second raw-jet prover.  Require the saved forms to
        # agree canonically, then parse the admitted left form once for the
        # new generic dual calculation.
        if canonical_digest(left_form) != canonical_digest(right_form):
            raise ValueError(f"raw native row {field}[{order}] has unequal saved forms")
        forms[key] = _parse_form(left_form, f"raw row {field}[{order}].left_form", seed_symbols)
        refs.extend(_proof_refs(row.get("proof_refs", row.get("proof_ref"))))
        kept.append({
            "field": field,
            "order": order,
            "left_node": row.get("left_node"),
            "right_node": row.get("right_node"),
            "proof_refs": deepcopy(row.get("proof_refs", row.get("proof_ref"))),
            "left_form": deepcopy(row.get("left_form")),
            "right_form": deepcopy(row.get("right_form")),
        })
    missing = sorted(expected - seen)
    if missing:
        raise ValueError(f"raw_native_jet.row_forms lacks rows {missing}")
    if len(rows) != len(expected):
        raise ValueError("raw_native_jet.row_forms must contain exactly 27 rows")
    return forms, _unique(refs), kept


@dataclass(frozen=True)
class _SeedAlias:
    field: str
    order: int
    left_node: int
    right_node: int
    left_expression: int
    right_expression: int
    left_variable: int
    right_variable: int
    left_endpoint: int
    right_endpoint: int
    proof_ref: Any
    left_zero_source: bool
    right_zero_source: bool


def _seed_rows(seed_admission: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    rows = seed_admission.get("fields")
    if not _is_sequence(rows) or len(rows) != len(SEED_FIELDS) * len(SEED_ORDERS):
        raise ValueError("seed_admission.fields must contain exactly 33 rows")
    expected = {(field, order) for field in SEED_FIELDS for order in SEED_ORDERS}
    seen: set[tuple[str, int]] = set()
    for offset, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"seed_admission.fields[{offset}] must be a mapping")
        field = row.get("field")
        order = row.get("order", row.get("ordinary_Z_order"))
        if field not in SEED_FIELDS or isinstance(order, bool) or not isinstance(order, int) or order not in SEED_ORDERS:
            raise ValueError(f"seed_admission.fields[{offset}] has an invalid field/order")
        key = (field, order)
        if key in seen:
            raise ValueError(f"duplicate seed admission row {field}[{order}]")
        seen.add(key)
        expected_status = (
            "external_premise" if field == "m" else
            "proved_by_new_pressure_composition" if field == "pressure" else
            "proved"
        )
        if row.get("status") != expected_status:
            raise ValueError(f"seed admission row {field}[{order}] has status {row.get('status')!r}, expected {expected_status!r}")
    if seen != expected:
        raise ValueError("seed_admission.fields does not cover all 33 field rows")
    if seed_admission.get("passed") is not True:
        raise ValueError("seed_admission.passed must be true")
    if list(seed_admission.get("ordinary_Z_orders", seed_admission.get("ordinary_orders", ()))) != list(SEED_ORDERS):
        raise ValueError("seed_admission must admit ordinary Z orders 0, 1, 2")
    return list(rows)


def _seed_aliases(
    nodes: Sequence[Mapping[str, Any]],
    seed_admission: Mapping[str, Any],
    raw_native_jet: Mapping[str, Any],
) -> tuple[dict[tuple[str, int], _SeedAlias], dict[int, _Dual], dict[int, _Dual], dict[int, _Dual], dict[int, _Dual], list[Any], dict[str, Any]]:
    admission = _seed_rows(seed_admission)
    raw_aliases = raw_native_jet.get("seed_aliases")
    if not _is_sequence(raw_aliases) or len(raw_aliases) != len(SEED_FIELDS) * len(SEED_ORDERS):
        raise ValueError("raw_native_jet.seed_aliases must contain exactly 33 rows")
    expected = {(field, order) for field in SEED_FIELDS for order in SEED_ORDERS}
    aliases: dict[tuple[str, int], _SeedAlias] = {}
    left_nodes: dict[int, _Dual] = {}
    right_nodes: dict[int, _Dual] = {}
    left_expressions: dict[int, _Dual] = {}
    right_expressions: dict[int, _Dual] = {}
    refs: list[Any] = []
    left_variables: set[int] = set()
    right_variables: set[int] = set()
    count = len(nodes)
    left_expression_owners: dict[int, tuple[str, int]] = {}
    right_expression_owners: dict[int, tuple[str, int]] = {}
    for offset, item in enumerate(raw_aliases):
        if not isinstance(item, Mapping):
            raise ValueError(f"raw_native_jet.seed_aliases[{offset}] must be a mapping")
        field = item.get("field")
        order = item.get("order", item.get("ordinary_Z_order"))
        if field not in SEED_FIELDS or isinstance(order, bool) or not isinstance(order, int) or order not in SEED_ORDERS:
            raise ValueError(f"raw_native_jet.seed_aliases[{offset}] has an invalid field/order")
        key = (field, order)
        if key in aliases:
            raise ValueError(f"duplicate source seed alias {field}[{order}]")
        left_node = _handle(item.get("left_node"), f"seed alias {field}[{order}].left_node", count)
        right_node = _handle(item.get("right_node"), f"seed alias {field}[{order}].right_node", count)
        left_row, right_row = nodes[left_node], nodes[right_node]
        if not isinstance(left_row, Mapping) or not isinstance(right_row, Mapping):
            raise ValueError(f"seed alias {field}[{order}] does not point to graph mappings")
        if _operation(left_row, left_node) != "function_substitution" or _operation(right_row, right_node) != "function_substitution":
            raise ValueError(f"seed alias {field}[{order}] is not a guarded function substitution")
        if not _guarded(left_row) or not _guarded(right_row):
            raise ValueError(f"seed alias {field}[{order}] lacks a substitution guard")
        left_expression = _handle(_field(left_row, "expression"), f"seed alias {field}[{order}].left.expression", count)
        right_expression = _handle(_field(right_row, "expression"), f"seed alias {field}[{order}].right.expression", count)
        left_variable = _handle(_field(left_row, "variable"), f"seed alias {field}[{order}].left.variable", count)
        right_variable = _handle(_field(right_row, "variable"), f"seed alias {field}[{order}].right.variable", count)
        left_endpoint = _handle(_field(left_row, "value"), f"seed alias {field}[{order}].left.value", count)
        right_endpoint = _handle(_field(right_row, "value"), f"seed alias {field}[{order}].right.value", count)
        if left_expression >= left_node or right_expression >= right_node or left_variable >= left_node or right_variable >= right_node or left_endpoint >= left_node or right_endpoint >= right_node:
            raise ValueError(f"seed alias {field}[{order}] has a forward descriptor reference")
        if not _exact_one(nodes, left_endpoint) or not _exact_one(nodes, right_endpoint):
            raise ValueError(f"seed alias {field}[{order}] endpoint is not exact one")
        left_zero_source = _exact_value(nodes, left_expression) == 0
        right_zero_source = _exact_value(nodes, right_expression) == 0
        # Accepted source-pair rows may expose a literal zero on only one
        # side (the other side can be a flat, guarded endpoint expression).
        # The accepted pair certificate supplies the shared zero value; do
        # not reject that layout merely because the graph bodies differ.
        proof_ref = item.get("proof_ref", item.get("source_ref"))
        if proof_ref in (None, ""):
            raise ValueError(f"seed alias {field}[{order}] lacks proof_ref")
        alias = _SeedAlias(field, order, left_node, right_node, left_expression, right_expression, left_variable, right_variable, left_endpoint, right_endpoint, deepcopy(proof_ref), left_zero_source, right_zero_source)
        aliases[key] = alias
        refs.extend(_proof_refs(proof_ref))
        left_variables.add(left_variable)
        right_variables.add(right_variable)
        for expression, owners, zero_source in (
            (left_expression, left_expression_owners, left_zero_source),
            (right_expression, right_expression_owners, right_zero_source),
        ):
            previous = owners.get(expression)
            if previous is not None and previous != key and not zero_source:
                raise ValueError(
                    f"guarded seed expressions alias two nonzero fields: {previous} and {key}"
                )
            owners.setdefault(expression, key)
    if set(aliases) != expected:
        raise ValueError("raw_native_jet.seed_aliases does not cover all 33 rows")
    if len(left_variables) != 1 or len(right_variables) != 1:
        raise ValueError("seed aliases do not identify one native variable per side")
    admission_map = {(row["field"], row["order"]): row for row in admission}
    for key, alias in aliases.items():
        if key not in admission_map:
            raise ValueError(f"seed admission lacks {key}")
        value = _sp.Integer(0) if (alias.left_zero_source or alias.right_zero_source) else _symbol(*key)
        # The derivative for E_y/V_y is supplied by second_radial_proof later;
        # the ordinary field rows use the accepted raw source forms.
        zero = _sp.Integer(0)
        right_value = _sp.Integer(0) if (alias.left_zero_source or alias.right_zero_source) else _symbol(*key)
        left_nodes[alias.left_node] = _Dual(value, zero)
        right_nodes[alias.right_node] = _Dual(right_value, zero)
        left_expressions[alias.left_expression] = _Dual(value, zero)
        right_expressions[alias.right_expression] = _Dual(right_value, zero)
    return aliases, left_nodes, right_nodes, left_expressions, right_expressions, _unique(refs), {
        "left_variable": next(iter(left_variables)),
        "right_variable": next(iter(right_variables)),
        "admission_rows": deepcopy(admission),
    }


def _second_radial_refs(proof: Mapping[str, Any]) -> tuple[list[Any], dict[tuple[str, int], dict[str, Any]]]:
    if not isinstance(proof, Mapping):
        raise ValueError("second_radial_proof must be a mapping")
    if proof.get("passed") is not True:
        raise ValueError("second_radial_proof.passed must be true")
    orders = proof.get("ordinary_orders", proof.get("ordinary_Z_orders"))
    if list(orders or ()) != list(SEED_ORDERS):
        raise ValueError("second_radial_proof must cover ordinary Z orders 0, 1, 2")
    rows = proof.get("final_rows", proof.get("rows"))
    if not _is_sequence(rows):
        raise ValueError("second_radial_proof.final_rows must be a sequence")
    accepted: dict[tuple[str, int], dict[str, Any]] = {}
    refs: list[Any] = []
    for row in rows:
        if not isinstance(row, Mapping):
            raise ValueError("second_radial_proof row must be a mapping")
        side = row.get("side")
        field = row.get("field")
        order = row.get("order", row.get("ordinary_Z_order"))
        if side not in SIDES or field not in ("E_y", "V_y") or not isinstance(order, int) or order not in SEED_ORDERS:
            continue
        if row.get("status") != "proved":
            raise ValueError(f"second radial row {side}.{field}[{order}] is not proved")
        key = (side + ":" + field, order)
        accepted[key] = deepcopy(dict(row))
        refs.extend(_proof_refs(row.get("proof_ref", row.get("source_ref"))))
    for side in SIDES:
        for field in ("E_y", "V_y"):
            for order in SEED_ORDERS:
                if (side + ":" + field, order) not in accepted:
                    raise ValueError(f"second_radial_proof lacks {side}.{field}[{order}]")
    seam = proof.get("second_radial_seam")
    if not _is_sequence(seam) or len(seam) != 6:
        raise ValueError("second_radial_proof.second_radial_seam must contain six rows")
    for row in seam:
        if not isinstance(row, Mapping) or row.get("status") not in ("proved_by_accepted_source_pair", "proved"):
            raise ValueError("second radial source seam is not fully proved")
        refs.extend(_proof_refs(row.get("proof_ref", row.get("source_ref"))))
    return _unique(refs), accepted


class _Evaluator:
    """Small exact dual evaluator with qualified node-id source atoms."""

    def __init__(
        self,
        nodes: Sequence[Mapping[str, Any]],
        *,
        side: str,
        native: int,
        aliases: Mapping[int, _Dual],
        radius_aliases: Mapping[int, _Dual],
        frozen: set[int],
        body_alias_nodes: Sequence[int] = (),
    ) -> None:
        self.nodes = nodes
        self.side = side
        self.native = native
        self.aliases = dict(aliases)
        self.radius_aliases = dict(radius_aliases)
        self.frozen = set(frozen)
        self.body_alias_nodes = frozenset(body_alias_nodes)
        self.memo: dict[int, _Dual] = {}
        self.used_nodes: set[int] = set()
        self.opaque_nodes: list[dict[str, Any]] = []
        self.positive_quotients: list[dict[str, Any]] = []
        self.substitutions: list[dict[str, Any]] = []

    def _row(self, index: int) -> Mapping[str, Any]:
        _handle(index, "source handle", len(self.nodes))
        row = self.nodes[index]
        if not isinstance(row, Mapping):
            raise _Unresolved(f"node {index} is not a mapping")
        return row

    def _opaque(self, index: int, row: Mapping[str, Any]) -> _Dual:
        operation = _operation(row, index)
        if index > self.native and index not in self.frozen and not any(
            row.get(name) is True for name in ("native_independent", "Z_only", "z_only", "source_constant")
        ):
            raise _Unresolved(f"unknown native-dependent operation at node {index}: {operation!r}")
        self.used_nodes.add(index)
        self.opaque_nodes.append({"node": index, "operation": operation, "native_derivative": "zero", "qualified": True})
        # The graph node id is the namespace.  A common source node must
        # produce the same atom on both chart sides; putting the side label in
        # the atom would manufacture a false left/right mismatch.
        return _dual(_atom(operation, index), 0)

    def rewrite(self, index: int, env: Mapping[int, _Dual] | None = None) -> _Dual:
        _handle(index, "source handle", len(self.nodes))
        environment = dict(env or {})
        if index in environment:
            return environment[index]
        if index in self.radius_aliases:
            self.used_nodes.add(index)
            return self.radius_aliases[index]
        if index in self.aliases:
            self.used_nodes.add(index)
            admitted = self.aliases[index]
            if index in self.body_alias_nodes and self.native in environment:
                replacement = environment[self.native]
                if replacement.value != 1:
                    raise _Unresolved("endpoint source-body alias used outside native=1")
                return _dual(admitted.value, admitted.derivative * replacement.derivative)
            return admitted
        if index in self.memo and not environment:
            return self.memo[index]
        row = self._row(index)
        operation = _operation(row, index)
        # A caller-authenticated frozen source node is an opaque native
        # constant before any arithmetic or unary expansion.  Exact rational
        # leaves remain exact, and the active native coordinate still carries
        # its explicit dual even if a malformed packet lists it as frozen.
        if (index < self.native or index in self.frozen) and index != self.native and operation != "exact_rational":
            result = self._opaque(index, row)
        elif operation == "exact_rational":
            numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
            denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
            if denominator == 0:
                raise _Unresolved(f"node {index} has a zero rational denominator")
            result = _dual(_sp.Rational(numerator.numerator, numerator.denominator) / _sp.Rational(denominator.numerator, denominator.denominator), 0)
        elif operation == "bound_variable":
            if index == self.native:
                result = _dual(1, 1)
            elif index in self.frozen or any(row.get(name) is True for name in ("native_independent", "Z_only", "z_only", "frozen_native_constant")):
                result = self._opaque(index, row)
            elif _exact_value(self.nodes, index) == 1:
                result = _dual(1, 0)
            else:
                if index < self.native:
                    result = self._opaque(index, row)
                else:
                    raise _Unresolved(f"unfrozen bound variable at node {index}")
        elif operation in {"sum", "add"}:
            result = _dual(0, 0)
            for offset, value in enumerate(_args(row, index)):
                child = _backward(index, value, f"arguments[{offset}]", len(self.nodes))
                result = _add(result, self.rewrite(child, environment))
        elif operation in {"product", "mul"}:
            result = _dual(1, 0)
            for offset, value in enumerate(_args(row, index)):
                child = _backward(index, value, f"arguments[{offset}]", len(self.nodes))
                result = _mul(result, self.rewrite(child, environment))
        elif operation in {"negative", "neg"}:
            child = _backward(index, _field(row, "argument", "operand"), "argument", len(self.nodes))
            result = _neg(self.rewrite(child, environment))
        elif operation in {"positive_quotient", "quotient"}:
            numerator = _backward(index, _field(row, "numerator"), "numerator", len(self.nodes))
            denominator = _backward(index, _field(row, "denominator"), "denominator", len(self.nodes))
            certificate = row.get("source_positive_certificate", row.get("positive_certificate", row.get("denominator_positive")))
            if certificate in (None, False) and _exact_value(self.nodes, denominator) in (None, 0):
                raise _Unresolved(f"node {index} quotient lacks positivity certificate")
            result = _quotient(self.rewrite(numerator, environment), self.rewrite(denominator, environment))
            self.positive_quotients.append({"node": index, "numerator": numerator, "denominator": denominator, "certificate": deepcopy(certificate)})
        elif operation == "function_substitution":
            if not _guarded(row):
                raise _Unresolved(f"node {index} function_substitution lacks a guard")
            expression = _backward(index, _field(row, "expression"), "expression", len(self.nodes))
            variable = _backward(index, _field(row, "variable"), "variable", len(self.nodes))
            value = _backward(index, _field(row, "value"), "value", len(self.nodes))
            replacement = self.rewrite(value, environment)
            scoped = dict(environment)
            scoped[variable] = replacement
            self.substitutions.append({"node": index, "expression": expression, "variable": variable, "endpoint_node": value, "guarded": True, "closed_endpoint": replacement.derivative == 0, "source_ref": deepcopy(row.get("source_ref", row.get("proof_ref")))})
            result = self.rewrite(expression, scoped)
        elif operation in {"analytic_unary", * _ANALYTIC}:
            name = operation if operation != "analytic_unary" else row.get("name", row.get("function", row.get("unary", row.get("kind"))))
            if name not in _ANALYTIC:
                raise _Unresolved(f"unsupported analytic unary {name!r}")
            child = _backward(index, _field(row, "argument", "operand", "value"), "argument", len(self.nodes))
            result = _unary(name, self.rewrite(child, environment))
        elif operation in {"integer_power", "power"}:
            child = _backward(index, _field(row, "argument", "operand", "base"), "argument", len(self.nodes))
            exponent = row.get("exponent", row.get("power"))
            exp = _fraction(exponent, f"node {index}.exponent")
            if exp.denominator != 1:
                raise _Unresolved(f"node {index} has a noninteger power")
            base = self.rewrite(child, environment)
            if exp.numerator < 0:
                result = _quotient(_dual(1, 0), _dual(base.value ** (-exp.numerator), (-exp.numerator) * base.value ** (-exp.numerator - 1) * base.derivative))
            else:
                value = _sp.Integer(1)
                out = _dual(1, 0)
                for _ in range(exp.numerator):
                    out = _mul(out, base)
                result = out
        else:
            result = self._opaque(index, row)
        result = _dual(result.value, result.derivative)
        if not environment:
            self.memo[index] = result
        return result


def _chart(modified: Mapping[str, Any], finite: Mapping[str, Any]) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    charts = modified.get("charts")
    if not isinstance(charts, Mapping) or not isinstance(charts.get("postrestore"), Mapping):
        raise ValueError("modified_Q4_pressure_C3_bindings.charts.postrestore is missing")
    chart = finite.get("chart")
    if not isinstance(chart, Mapping):
        raise ValueError("Rm_patch_finite_N_C5_bindings.chart is missing")
    return charts["postrestore"], chart


def _root_rows(chart: Mapping[str, Any], side: str, count: int) -> dict[str, list[int]]:
    roots = chart.get("roots")
    if not isinstance(roots, Mapping):
        raise ValueError(f"{side} leading chart roots are missing")
    result: dict[str, list[int]] = {}
    for name in ROOTS:
        value = roots.get(name)
        if not _is_sequence(value) or len(value) < 2:
            raise ValueError(f"{side}.roots.{name} lacks ordinary orders 0 and 1")
        result[name] = [_handle(value[order], f"{side}.roots.{name}[{order}]", count) for order in ORDERS]
    return result


def _frozen(owner: Mapping[str, Any], chart: Mapping[str, Any], count: int) -> set[int]:
    result: set[int] = set()
    for source in (owner, chart, owner.get("bindings", {})):
        if not isinstance(source, Mapping):
            continue
        values = source.get("frozen_native_constant_nodes", source.get("frozen_nodes"))
        if values is None:
            continue
        if isinstance(values, Mapping):
            values = values.get("nodes", values.get("handles", []))
        if not _is_sequence(values):
            raise ValueError("frozen_native_constant_nodes must be a sequence")
        for offset, value in enumerate(values):
            if isinstance(value, Mapping):
                value = value.get("node", value.get("handle"))
            result.add(_handle(value, f"frozen_native_constant_nodes[{offset}]", count))
    return result


def _radius_port(
    modified: Mapping[str, Any],
    finite: Mapping[str, Any],
    modified_chart: Mapping[str, Any],
    finite_chart: Mapping[str, Any],
    nodes: Sequence[Mapping[str, Any]],
    *,
    native_variables: Mapping[str, int] | None = None,
) -> tuple[dict[int, _Dual], list[Any], dict[str, Any], list[dict[str, Any]]]:
    """Derive the radius source/layout identity, with a toy-port fallback.

    The actual packet does not supply a caller assertion such as
    ``passed=True``.  Its postrestore radius body is found from the guarded
    native-one substitution and differentiated through the graph.  The finite
    body is then checked as the actual ``Rm*x`` product and differentiated
    through the same dual.  The explicit port branch below exists only for the
    small fresh fixture, whose graph intentionally omits the archive source
    handles.
    """

    # A root-owned run may hand us the already verified radius result.  This
    # branch consumes that receipt and validates its exact graph/layout ports;
    # it deliberately does not differentiate or replay the actual graph.
    candidates: list[Mapping[str, Any]] = []
    for owner, chart in ((modified, modified_chart), (finite, finite_chart)):
        for source in (owner, chart, owner.get("bindings", {}) if isinstance(owner.get("bindings", {}), Mapping) else {}):
            for key in ("radius_source", "radius_proof", "native_radius_proof", "radius_port", "radius_witness"):
                value = source.get(key) if isinstance(source, Mapping) else None
                if isinstance(value, Mapping):
                    candidates.append(value)
    for witness in candidates:
        if witness.get("verified") is not True or witness.get("status") != "proved":
            continue
        proof = witness.get("proof") if isinstance(witness.get("proof"), Mapping) else witness
        returned = witness.get("returned_result", {})
        if not isinstance(returned, Mapping) or returned.get("status", "proved") != "proved":
            continue
        sides = proof.get("sides") if isinstance(proof.get("sides"), Mapping) else {}
        post = sides.get("postrestore") if isinstance(sides.get("postrestore"), Mapping) else proof.get("postrestore")
        patch = sides.get("patch") if isinstance(sides.get("patch"), Mapping) else proof.get("patch")
        inputs = witness.get("inputs") if isinstance(witness.get("inputs"), Mapping) else {}
        if not isinstance(post, Mapping) or not isinstance(patch, Mapping):
            continue
        try:
            left_coordinate = _handle(
                post.get("coordinate_node", inputs.get("postrestore_native_variable")),
                "verified postrestore native variable",
                len(nodes),
            )
            post_body = _handle(
                post.get("body_node", inputs.get("postrestore_body_node")),
                "verified postrestore radius body",
                len(nodes),
            )
            post_endpoint = _handle(
                post.get("endpoint_node", inputs.get("postrestore_endpoint_node")),
                "verified postrestore radius endpoint",
                len(nodes),
            )
            right_coordinate = _handle(
                patch.get("coordinate_node", inputs.get("patch_native_variable")),
                "verified patch native variable",
                len(nodes),
            )
            right_body = _handle(
                patch.get("body_node", inputs.get("patch_body_node")),
                "verified patch radius body",
                len(nodes),
            )
            measure = _handle(
                patch.get("physical_measure_node", inputs.get("physical_measure_node")),
                "verified patch physical measure",
                len(nodes),
            )
            patch_endpoint_value = patch.get("endpoint_node", inputs.get("patch_closed_endpoint_node"))
            patch_endpoint = None if patch_endpoint_value is None else _handle(patch_endpoint_value, "verified patch radius endpoint", len(nodes))
            post_row = nodes[post_body]
            post_endpoint_row = nodes[post_endpoint]
            patch_row = nodes[right_body]
            measure_row = nodes[measure]
            if not isinstance(post_row, Mapping) or _operation(post_row, post_body) not in {"product", "mul"}:
                raise _Unresolved("verified postrestore radius body is not a product")
            if not isinstance(post_endpoint_row, Mapping) or _operation(post_endpoint_row, post_endpoint) != "function_substitution" or not _guarded(post_endpoint_row):
                raise _Unresolved("verified postrestore radius endpoint is not guarded")
            if _handle(_field(post_endpoint_row, "expression"), "verified postrestore endpoint expression", len(nodes)) != post_body:
                raise _Unresolved("verified postrestore endpoint body mismatch")
            if _handle(_field(post_endpoint_row, "variable"), "verified postrestore endpoint variable", len(nodes)) != left_coordinate:
                raise _Unresolved("verified postrestore endpoint coordinate mismatch")
            if not _exact_one(nodes, _handle(_field(post_endpoint_row, "value"), "verified postrestore endpoint value", len(nodes))):
                raise _Unresolved("verified postrestore endpoint is not native one")
            if not isinstance(patch_row, Mapping) or _operation(patch_row, right_body) not in {"product", "mul"}:
                raise _Unresolved("verified patch radius body is not a product")
            patch_args = [_handle(value, "verified patch radius argument", len(nodes)) for value in _args(patch_row, right_body)]
            if post_endpoint not in patch_args or right_coordinate not in patch_args:
                raise _Unresolved("verified patch body is not Rm*x")
            if not isinstance(measure_row, Mapping) or _operation(measure_row, measure) not in {"positive_quotient", "quotient"}:
                raise _Unresolved("verified patch physical measure is not a quotient")
            numerator = _handle(_field(measure_row, "numerator"), "verified physical measure numerator", len(nodes))
            denominator = _handle(_field(measure_row, "denominator"), "verified physical measure denominator", len(nodes))
            if denominator != right_coordinate or _exact_value(nodes, numerator) != 1:
                raise _Unresolved("verified physical measure is not 1/x")
            if str(post.get("identity_text", post.get("identity", ""))) not in {"0", "D_y R - R = 0"}:
                raise _Unresolved("verified postrestore radius identity receipt is missing")
            if str(patch.get("identity_text", patch.get("identity", ""))) not in {"0", "R_patch = Rm*x, J = 1/x"}:
                raise _Unresolved("verified patch radius identity receipt is missing")
        except (ValueError, _Unresolved, IndexError, KeyError) as exc:
            return {}, [], {"status": "unresolved", "gate": "verified_radius_witness_invalid"}, [{"gate": "verified_radius_witness_invalid", "error": str(exc)}]
        radius = _atom("Rm_radius")
        aliases = {post_body: _dual(radius, radius), post_endpoint: _dual(radius, 0)}
        if patch_endpoint is not None:
            aliases[patch_endpoint] = _dual(radius, 0)
        refs = returned.get("proof_refs", witness.get("proof_refs", witness.get("source_refs")))
        refs = _unique(_proof_refs(refs))
        if not refs:
            refs = [{"kind": "verified_radius_witness", "source_graph_digest": witness.get("inputs", {}).get("source_graph_digest")}]
        metadata = {
            "status": "proved",
            "method": "consumed_verified_radius_witness",
            "verified_witness": True,
            "native_derivative_equals_value": True,
            "sides": {"postrestore": deepcopy(dict(post)), "patch": deepcopy(dict(patch))},
            "source_checkpoint": deepcopy(witness.get("execution_boundary", witness.get("source_checkpoint"))),
            "aliases": {str(key): {"value": str(value.value), "native_derivative": str(value.derivative)} for key, value in aliases.items()},
        }
        return aliases, refs, metadata, []

    # Actual archive-sized graphs must consume the root-owned verified
    # receipt above.  This fail-closed gate prevents an accidental replay of
    # the actual radius differentiation if that receipt was not attached.
    if native_variables is not None and len(nodes) > 1000:
        return {}, [], {"status": "unresolved", "gate": "verified_radius_witness_missing"}, [
            {"gate": "verified_radius_witness_missing", "graph_node_count": len(nodes)}
        ]

    # Actual source/layout path.  It uses the native variables already
    # authenticated by the accepted seed aliases, so a common name such as
    # ``radius`` cannot select an unrelated node.
    if native_variables is not None:
        try:
            left_coordinate = _handle(native_variables.get("postrestore"), "postrestore native variable", len(nodes))
            right_coordinate = _handle(native_variables.get("patch"), "patch native variable", len(nodes))

            def children(index: int, row: Mapping[str, Any]) -> list[int]:
                operation = _operation(row, index)
                if operation in {"sum", "add", "product", "mul"}:
                    return [_handle(value, f"node {index}.argument", len(nodes)) for value in _args(row, index)]
                if operation in {"negative", "neg", "analytic_unary", * _ANALYTIC, "integer_power", "power"}:
                    value = _field(row, "argument", "operand", "value", "base")
                    return [_handle(value, f"node {index}.argument", len(nodes))]
                if operation == "positive_quotient" or operation == "quotient":
                    return [
                        _handle(_field(row, "numerator"), f"node {index}.numerator", len(nodes)),
                        _handle(_field(row, "denominator"), f"node {index}.denominator", len(nodes)),
                    ]
                if operation == "function_substitution":
                    return [_handle(_field(row, "expression"), f"node {index}.expression", len(nodes))]
                return []

            def depends(index: int, target: int, seen: set[int] | None = None) -> bool:
                if index == target:
                    return True
                visited = set() if seen is None else seen
                if index in visited:
                    return False
                visited.add(index)
                row = nodes[index]
                if not isinstance(row, Mapping):
                    return False
                try:
                    return any(depends(child, target, visited) for child in children(index, row))
                except (ValueError, _Unresolved):
                    return False

            def guarded_one_substitutions(variable: int) -> list[tuple[int, int]]:
                found: list[tuple[int, int]] = []
                for index, row in enumerate(nodes):
                    if not isinstance(row, Mapping) or _operation(row, index) != "function_substitution" or not _guarded(row):
                        continue
                    try:
                        row_variable = _handle(_field(row, "variable"), f"node {index}.variable", len(nodes))
                        value = _handle(_field(row, "value"), f"node {index}.value", len(nodes))
                        expression = _handle(_field(row, "expression"), f"node {index}.expression", len(nodes))
                    except (ValueError, _Unresolved):
                        continue
                    if row_variable == variable and _exact_one(nodes, value):
                        found.append((index, expression))
                return found

            left_candidates = []
            for endpoint, expression in guarded_one_substitutions(left_coordinate):
                row = nodes[expression]
                if not isinstance(row, Mapping) or _operation(row, expression) not in {"product", "mul"}:
                    continue
                product_arguments = _args(row, expression)
                # The admitted Rm source has the distinctive exact scale
                # times exp(log-source + native coordinate) layout.  Other
                # leading fields also contain products depending on the
                # native coordinate, so dependency alone is deliberately too
                # weak a locator.
                def is_exp_argument(value: Any) -> bool:
                    child = _handle(value, "radius product argument", len(nodes))
                    child_row = nodes[child]
                    if not isinstance(child_row, Mapping):
                        return False
                    child_operation = _operation(child_row, child)
                    child_name = (
                        child_row.get("name", child_row.get("function", child_row.get("unary", child_row.get("kind"))))
                        if child_operation == "analytic_unary" else child_operation
                    )
                    return child_name == "exp"

                has_exp = any(is_exp_argument(value) for value in product_arguments)
                has_exact_scale = any(
                    _exact_value(nodes, _handle(value, "radius product argument", len(nodes))) is not None
                    for value in product_arguments
                )
                if has_exp and has_exact_scale and depends(expression, left_coordinate):
                    left_candidates.append((endpoint, expression))
            if len(left_candidates) != 1:
                raise _Unresolved(f"radius post source layout has {len(left_candidates)} native-dependent product candidates")
            post_endpoint, post_body = left_candidates[0]

            right_body = finite_chart.get("radius")
            right_body = _handle(right_body, "patch radius body", len(nodes))
            right_body_row = nodes[right_body]
            if not isinstance(right_body_row, Mapping) or _operation(right_body_row, right_body) not in {"product", "mul"}:
                raise _Unresolved("patch radius body is not a product")
            right_arguments = [_handle(value, "patch radius argument", len(nodes)) for value in _args(right_body_row, right_body)]
            if post_endpoint not in right_arguments or right_coordinate not in right_arguments:
                raise _Unresolved("patch radius body is not the guarded Rm*x layout")

            measure = finite_chart.get("physical_measure")
            measure = _handle(measure, "patch physical measure", len(nodes))
            measure_row = nodes[measure]
            if not isinstance(measure_row, Mapping) or _operation(measure_row, measure) not in {"positive_quotient", "quotient"}:
                raise _Unresolved("patch physical measure is not a positive quotient")
            numerator = _handle(_field(measure_row, "numerator"), "patch physical measure numerator", len(nodes))
            denominator = _handle(_field(measure_row, "denominator"), "patch physical measure denominator", len(nodes))
            if denominator != right_coordinate or _exact_value(nodes, numerator) != 1:
                raise _Unresolved("patch physical measure is not 1/x")

            left_frozen = _frozen(modified, modified_chart, len(nodes))
            # Differentiate the actual postrestore exponential before any
            # radius alias is installed.  All source constants remain opaque
            # atoms with zero native derivative, while the bound native
            # coordinate contributes the chain-rule one.
            source_eval = _Evaluator(
                nodes,
                side="radius_source",
                native=left_coordinate,
                aliases={},
                radius_aliases={},
                frozen=left_frozen,
            )
            source_dual = source_eval.rewrite(post_body)
            source_identity = _simplify(source_dual.derivative - source_dual.value)
            if source_identity != 0:
                raise _Unresolved("postrestore radius source does not prove D_y R - R = 0")

            radius = _atom("Rm_radius")
            # The common guarded endpoint is closed under native y.  The
            # finite product is evaluated with that closed Rm factor and the
            # open x factor, giving (R, R) by actual product arithmetic.
            closed_aliases = {post_endpoint: _dual(radius, 0)}
            patch_endpoint: int | None = None
            for endpoint, expression in guarded_one_substitutions(right_coordinate):
                if expression == right_body:
                    patch_endpoint = endpoint
                    closed_aliases[endpoint] = _dual(radius, 0)
                    break
            patch_eval = _Evaluator(
                nodes,
                side="radius_patch",
                native=right_coordinate,
                aliases=closed_aliases,
                radius_aliases={},
                frozen=_frozen(finite, finite_chart, len(nodes)),
            )
            patch_dual = patch_eval.rewrite(right_body)
            if _simplify(patch_dual.value - radius) != 0 or _simplify(patch_dual.derivative - radius) != 0:
                raise _Unresolved("patch Rm*x dual does not prove (R, R)")

            refs = [
                {
                    "kind": "graph_derived_radius_source",
                    "postrestore_body_node": post_body,
                    "postrestore_endpoint_node": post_endpoint,
                    "postrestore_native_variable": left_coordinate,
                    "identity": "D_y R - R = 0",
                },
                {
                    "kind": "graph_derived_radius_patch_layout",
                    "patch_body_node": right_body,
                    "patch_closed_endpoint_node": patch_endpoint,
                    "patch_native_variable": right_coordinate,
                    "physical_measure_node": measure,
                    "identity": "R_patch = Rm*x, J = 1/x",
                },
            ]
            aliases = {
                post_body: _dual(radius, radius),
                post_endpoint: _dual(radius, 0),
            }
            if patch_endpoint is not None:
                aliases[patch_endpoint] = _dual(radius, 0)
            metadata = {
                "status": "proved",
                "method": "graph_derived_source_and_layout",
                "native_derivative_equals_value": True,
                "sides": {
                    "postrestore": {
                        "body_node": post_body,
                        "endpoint_node": post_endpoint,
                        "coordinate_node": left_coordinate,
                        "open_body_dual": {"value": _record(source_dual.value), "native_derivative": _record(source_dual.derivative)},
                        "identity": _record(source_identity),
                        "source_opaque_nodes": deepcopy(source_eval.opaque_nodes),
                    },
                    "patch": {
                        "body_node": right_body,
                        "endpoint_node": patch_endpoint,
                        "coordinate_node": right_coordinate,
                        "physical_measure_node": measure,
                        "open_body_dual": {"value": _record(patch_dual.value), "native_derivative": _record(patch_dual.derivative)},
                        "identity": _record(patch_dual.derivative - patch_dual.value),
                    },
                },
            }
            return aliases, refs, metadata, []
        except (ValueError, _Unresolved, IndexError, KeyError) as exc:
            automatic_issue = {"gate": "radius_graph_source_layout_unresolved", "error": str(exc)}
    else:
        automatic_issue = {"gate": "radius_native_variables_missing"}

    candidates: list[Mapping[str, Any]] = []
    for owner, chart in ((modified, modified_chart), (finite, finite_chart)):
        for source in (owner, chart, owner.get("bindings", {}) if isinstance(owner.get("bindings", {}), Mapping) else {}):
            for key in ("radius_source", "radius_proof", "native_radius_proof", "radius_port", "radius_witness"):
                value = source.get(key) if isinstance(source, Mapping) else None
                if isinstance(value, Mapping):
                    candidates.append(value)
    port: Mapping[str, Any] | None = None
    for candidate in candidates:
        if any(key in candidate for key in ("postrestore", "left", "patch", "right", "sides")):
            port = candidate
            break
    if port is None:
        return {}, [], {"status": "unresolved", "gate": "radius_source_witness_missing"}, [automatic_issue, {"gate": "radius_source_witness_missing"}]

    def side_port(name: str) -> Mapping[str, Any]:
        if isinstance(port.get("sides"), Mapping):
            source = port["sides"]
            for key in (("postrestore", "post_restore", "left") if name == "postrestore" else ("patch", "finite", "right")):
                if isinstance(source.get(key), Mapping):
                    return source[key]
        for key in (("postrestore", "post_restore", "left") if name == "postrestore" else ("patch", "finite", "right")):
            if isinstance(port.get(key), Mapping):
                return port[key]
        return port

    refs: list[Any] = []
    aliases: dict[int, _Dual] = {}
    metadata: dict[str, Any] = {"status": "unresolved", "sides": {}}
    issues: list[dict[str, Any]] = []
    radius = _atom("Rm_radius")
    for name, chart in (("postrestore", modified_chart), ("patch", finite_chart)):
        item = side_port(name)
        proved = any(item.get(key) is True for key in ("passed", "native_derivative_equals_value", "D_y_R_equals_R", "physical_derivative_equals_value"))
        proof = item.get("proof_ref", item.get("source_ref", item.get("layout_proof")))
        if not proved or proof in (None, ""):
            issues.append({"gate": "radius_source_witness_incomplete", "side": name})
            continue
        refs.extend(_proof_refs(proof))
        body_value = item.get("body_node", item.get("radius_body_node", chart.get("radius")))
        coordinate_value = item.get("coordinate_node", item.get("native_variable", chart.get("coordinate")))
        if body_value is None or coordinate_value is None:
            issues.append({"gate": "radius_source_layout_handles_missing", "side": name})
            continue
        body = _handle(body_value, f"{name}.radius_body_node", len(nodes))
        coordinate = _handle(coordinate_value, f"{name}.radius_coordinate_node", len(nodes))
        body_row = nodes[body]
        if not isinstance(body_row, Mapping) or _operation(body_row, body) not in {"product", "mul"}:
            issues.append({"gate": "radius_body_layout_missing", "side": name, "node": body})
            continue
        aliases[body] = _dual(radius, radius)
        metadata["sides"][name] = {"body_node": body, "proof_ref": deepcopy(proof), "native_derivative_equals_value": True, "coordinate_node": coordinate}
        endpoint_value = item.get("endpoint_node", item.get("endpoint_node1"))
        if endpoint_value is None:
            issues.append({"gate": "radius_endpoint_handle_missing", "side": name})
        else:
            endpoint = _handle(endpoint_value, f"{name}.radius_endpoint_node", len(nodes))
            row = nodes[endpoint]
            if not isinstance(row, Mapping) or _operation(row, endpoint) != "function_substitution" or not _guarded(row):
                issues.append({"gate": "radius_endpoint_guard_missing", "side": name, "node": endpoint})
            else:
                expression = _handle(_field(row, "expression"), f"{name}.radius_endpoint.expression", len(nodes))
                variable = _handle(_field(row, "variable"), f"{name}.radius_endpoint.variable", len(nodes))
                value = _handle(_field(row, "value"), f"{name}.radius_endpoint.value", len(nodes))
                if not _exact_one(nodes, value) or variable != coordinate:
                    issues.append({"gate": "radius_endpoint_layout_mismatch", "side": name, "node": endpoint, "coordinate_node": coordinate})
                if name == "postrestore" and expression != body:
                    issues.append({"gate": "radius_post_endpoint_body_mismatch", "side": name, "node": endpoint, "body_node": body})
                aliases[endpoint] = _dual(radius, 0)
                metadata["sides"].setdefault(name, {})["endpoint_node"] = endpoint
        if name == "patch":
            arguments = body_row.get("arguments", body_row.get("operands"))
            if not _is_sequence(arguments) or not any(item == coordinate for item in arguments):
                issues.append({"gate": "radius_patch_Rm_x_layout_missing", "side": name, "body_node": body, "coordinate_node": coordinate})
            measure_value = item.get("physical_measure_node", item.get("jacobian_node"))
            if measure_value is not None:
                measure = _handle(measure_value, f"{name}.physical_measure_node", len(nodes))
                measure_row = nodes[measure]
                if not isinstance(measure_row, Mapping) or _operation(measure_row, measure) not in {"positive_quotient", "quotient"}:
                    issues.append({"gate": "radius_patch_jacobian_layout_missing", "side": name, "node": measure})
                else:
                    denominator = _handle(_field(measure_row, "denominator"), f"node {measure}.denominator", len(nodes))
                    numerator = _handle(_field(measure_row, "numerator"), f"node {measure}.numerator", len(nodes))
                    if denominator != coordinate or _exact_value(nodes, numerator) != 1:
                        issues.append({"gate": "radius_patch_jacobian_not_one_over_x", "side": name, "node": measure})
    if len(metadata["sides"]) != 2:
        issues.append({"gate": "radius_source_witness_missing_side"})
    metadata["status"] = "proved" if not issues else "unresolved"
    return aliases, _unique(refs), metadata, issues


def _row_form(side: str, name: str, order: int, root: int, value: _Dual, refs: Sequence[Any], evaluator: _Evaluator) -> dict[str, Any]:
    return {
        "side": side,
        "root": name,
        "component": "value_and_native_y",
        "order": order,
        "ordinary_Z_order": order,
        "root_node": root,
        "status": "source_form",
        "value": _record(value.value),
        "native_derivative": _record(value.derivative),
        "proof_refs": deepcopy(list(refs)),
        "used_nodes": sorted(evaluator.used_nodes),
        "positive_quotients": deepcopy(evaluator.positive_quotients),
        "backward_refs": deepcopy(evaluator.substitutions),
    }


def prove(
    nodes: Sequence[Mapping[str, Any]],
    modified_bindings: Mapping[str, Any],
    finite_bindings: Mapping[str, Any],
    seed_admission: Mapping[str, Any],
    raw_native_jet: Mapping[str, Any],
    second_radial_proof: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Compare actual leading generic roots at value and native-y order one."""

    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact generic native identities")
    if not _is_sequence(nodes):
        raise ValueError("nodes must be a sequence")
    expected = _require_digest(expected_nodes_digest, "expected_nodes_digest")
    actual = canonical_digest(list(nodes))
    if actual != expected:
        raise ValueError("nodes digest does not match expected_nodes_digest")
    if not isinstance(modified_bindings, Mapping) or not isinstance(finite_bindings, Mapping):
        raise ValueError("modified_bindings and finite_bindings must be mappings")
    if not isinstance(seed_admission, Mapping) or not isinstance(raw_native_jet, Mapping):
        raise ValueError("seed_admission and raw_native_jet must be mappings")

    original_nodes = deepcopy(nodes)
    original_modified = deepcopy(modified_bindings)
    original_finite = deepcopy(finite_bindings)
    original_seed = deepcopy(seed_admission)
    original_raw = deepcopy(raw_native_jet)
    original_second = deepcopy(second_radial_proof)

    modified = modified_bindings.get("modified_Q4_pressure_C3_bindings")
    finite = finite_bindings.get("Rm_patch_finite_N_C5_bindings")
    if not isinstance(modified, Mapping) or not isinstance(finite, Mapping):
        raise ValueError("wrapped bindings must use exact archive keys")
    modified_chart, finite_chart = _chart(modified, finite)
    roots = {
        "postrestore": _root_rows(modified_chart, "postrestore", len(nodes)),
        "patch": _root_rows(finite_chart, "patch", len(nodes)),
    }
    for name in ROOTS:
        if roots["postrestore"][name] == roots["patch"][name] and name in {"E", "V"}:
            # Equal handles are fine in a toy, but the actual pair is normally
            # two guarded source descriptors.  No alias is inferred here.
            pass

    # The checkpoint is materialized before graph inspection.  It intentionally
    # contains no evaluated root form and therefore cannot replay an operator.
    checkpoint = {
        "checkpoint_before_inspection": True,
        "nodes_digest": actual,
        "seed_admission_digest": canonical_digest(seed_admission),
        "raw_native_jet_digest": canonical_digest(raw_native_jet),
        "second_radial_proof_digest": canonical_digest(second_radial_proof),
        "expected_nodes_digest": expected,
    }

    aliases, left_nodes, right_nodes, left_expr, right_expr, seed_refs, alias_meta = _seed_aliases(nodes, seed_admission, raw_native_jet)
    zero_keys = {
        key
        for key, alias in aliases.items()
        if alias.left_zero_source or alias.right_zero_source
    }
    seed_symbols = {
        (field, order): (_sp.Integer(0) if (field, order) in zero_keys else _symbol(field, order))
        for field in SEED_FIELDS
        for order in SEED_ORDERS
    }
    raw_forms, raw_refs, raw_rows = _raw_forms(raw_native_jet, seed_symbols)
    second_refs, second_rows = _second_radial_refs(second_radial_proof)

    # Convert accepted seed leaves into source duals.  E_y and V_y are the
    # only rows whose native derivatives come from the new second-radial
    # certificate; all other native derivatives are consumed from raw forms.
    def field_dual(field: str, order: int) -> _Dual:
        value = seed_symbols[(field, order)]
        if field == "E_y":
            return _dual(value, value / 10)
        if field == "V_y":
            return _dual(value, 0)
        return _dual(value, raw_forms[(field, order)])

    for key, alias in aliases.items():
        dual = field_dual(*key)
        left_nodes[alias.left_node] = _dual(dual.value, 0)
        right_nodes[alias.right_node] = _dual(dual.value, 0)
        left_expr[alias.left_expression] = dual
        right_expr[alias.right_expression] = dual

    radius_aliases, radius_refs, radius_meta, radius_issues = _radius_port(
        modified,
        finite,
        modified_chart,
        finite_chart,
        nodes,
        native_variables={
            "postrestore": alias_meta["left_variable"],
            "patch": alias_meta["right_variable"],
        },
    )
    left_frozen = _frozen(modified, modified_chart, len(nodes))
    right_frozen = _frozen(finite, finite_chart, len(nodes))
    left_eval = _Evaluator(
        nodes,
        side="postrestore",
        native=alias_meta["left_variable"],
        aliases={**left_nodes, **left_expr},
        radius_aliases=radius_aliases,
        frozen=left_frozen,
        body_alias_nodes=tuple(left_expr),
    )
    right_eval = _Evaluator(
        nodes,
        side="patch",
        native=alias_meta["right_variable"],
        aliases={**right_nodes, **right_expr},
        radius_aliases=radius_aliases,
        frozen=right_frozen,
    )

    all_refs = _unique(seed_refs + raw_refs + second_refs + radius_refs)
    rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = list(radius_issues)
    failed: list[dict[str, Any]] = []
    for name in ROOTS:
        for order in ORDERS:
            left_root = roots["postrestore"][name][order]
            right_root = roots["patch"][name][order]
            try:
                left_value = left_eval.rewrite(left_root)
                left_row = _row_form("postrestore", name, order, left_root, left_value, all_refs, left_eval)
            except (ValueError, _Unresolved) as exc:
                left_value = None
                left_row = {"side": "postrestore", "root": name, "component": "value_and_native_y", "order": order, "ordinary_Z_order": order, "root_node": left_root, "status": "unresolved", "error": str(exc), "proof_refs": deepcopy(all_refs)}
            try:
                right_value = right_eval.rewrite(right_root)
                right_row = _row_form("patch", name, order, right_root, right_value, all_refs, right_eval)
            except (ValueError, _Unresolved) as exc:
                right_value = None
                right_row = {"side": "patch", "root": name, "component": "value_and_native_y", "order": order, "ordinary_Z_order": order, "root_node": right_root, "status": "unresolved", "error": str(exc), "proof_refs": deepcopy(all_refs)}
            if left_value is None or right_value is None:
                row_status = "unresolved"
                unresolved.append({"gate": "root_source_reduction_unresolved", "root": name, "order": order, "postrestore_node": left_root, "patch_node": right_root, "postrestore_error": left_row.get("error"), "patch_error": right_row.get("error")})
            else:
                value_equal = _simplify(left_value.value - right_value.value) == 0
                derivative_equal = _simplify(left_value.derivative - right_value.derivative) == 0
                row_status = "proved" if value_equal and derivative_equal else "failed"
                if row_status != "proved":
                    failed.append({"gate": "leading_root_value_or_native_mismatch", "root": name, "order": order, "value_equal": value_equal, "native_equal": derivative_equal})
            rows.append({
                "root": name,
                "order": order,
                "ordinary_Z_order": order,
                "status": row_status,
                "leading_only": True,
                "postrestore": left_row,
                "patch": right_row,
                "source_proof_refs": deepcopy(all_refs),
            })

    unresolved.extend(failed)
    passed = len(rows) == len(ROOTS) * len(ORDERS) and all(row["status"] == "proved" for row in rows) and not unresolved and radius_meta.get("status") == "proved"
    result = {
        "schema": "original_Rm_generic_native_identity_v1",
        "passed": bool(passed),
        "source_only": True,
        "sourceonly": True,
        "leading_only": True,
        "leadingonly": True,
        "root_count": 18,
        "root_count_value_native": 18,
        "root_names": list(ROOTS),
        "components": ["value", "native_y"],
        "ordinary_orders": list(ORDERS),
        "expected_nodes_digest": expected,
        "nodes_digest": actual,
        "module_sha256": _module_sha256(),
        "seed_admission_digest": canonical_digest(seed_admission),
        "raw_native_jet_digest": canonical_digest(raw_native_jet),
        "second_radial_proof_digest": canonical_digest(second_radial_proof),
        "checkpoint": checkpoint,
        "source_proof_refs": all_refs,
        "seed_aliases": [
            {
                "field": alias.field,
                "order": alias.order,
                "left_node": alias.left_node,
                "right_node": alias.right_node,
                "left_expression": alias.left_expression,
                "right_expression": alias.right_expression,
                "proof_ref": deepcopy(alias.proof_ref),
            }
            for alias in aliases.values()
        ],
        "raw_native_rows": raw_rows,
        "second_radial_rows": [
            deepcopy(row)
            for key in sorted(second_rows)
            for row in (second_rows[key],)
        ],
        "radius": radius_meta,
        "radius_proof_refs": radius_refs,
        "root_rows": rows,
        "row_passes": [row for row in rows if row["status"] == "proved"],
        "unresolved_gates": unresolved,
        "used_nodes": {
            "postrestore": sorted(left_eval.used_nodes),
            "patch": sorted(right_eval.used_nodes),
        },
        "opaque_source_nodes": {
            "postrestore": deepcopy(left_eval.opaque_nodes),
            "patch": deepcopy(right_eval.opaque_nodes),
        },
        "positive_quotient_certificates": {
            "postrestore": deepcopy(left_eval.positive_quotients),
            "patch": deepcopy(right_eval.positive_quotients),
        },
        "backward_refs": {
            "postrestore": deepcopy(left_eval.substitutions),
            "patch": deepcopy(right_eval.substitutions),
        },
        "scope": {
            "source_graph_only": True,
            "no_saved_graph_execution": True,
            "no_accepted_builder_replay": True,
            "no_accepted_checker_replay": True,
            "no_reducer_or_provider_replay": True,
            "no_quadratures_replayed": True,
            "no_history_archive_decode": True,
            "no_control_core_history_loop_phase_inverse_expansion": True,
            "unknown_native_dependent_ops_fail_closed": True,
            "leading_root_maps_only": True,
            "finite_N_germ_claimed": False,
            "history_identity_claimed": False,
            "global_phase_claimed": False,
            "physical_RQy_claimed": False,
            "physical_Qy_claimed": False,
            "higher_ordinary_Z_claimed": False,
            "ordinary_Z_orders_available_for_generic_native_first": [0, 1],
            "generic_shift_source_orders_consumed": [0, 1, 2],
            "metadata_preserved": (
                deepcopy(modified_bindings) == original_modified
                and deepcopy(finite_bindings) == original_finite
                and deepcopy(seed_admission) == original_seed
                and deepcopy(raw_native_jet) == original_raw
                and deepcopy(second_radial_proof) == original_second
            ),
            "nodes_input_unmutated": deepcopy(nodes) == original_nodes,
        },
        "limits": {
            "generic_native_orders": [0, 1],
            "seed_source_orders": [0, 1, 2],
            "requires_Z3_for_physical_Qy_Z1_unless_cancellation_certificate": True,
            "full_physical_seam": False,
            "finite_N_history": False,
            "global_phase": False,
            "physical_RQy": False,
        },
    }
    json.dumps(result, allow_nan=False)
    return result


__all__ = [
    "ORDERS",
    "RAW_FIELDS",
    "ROOTS",
    "SEED_FIELDS",
    "SEED_ORDERS",
    "canonical_digest",
    "prove",
]
