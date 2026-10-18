"""Bounded source-only endpoint differentiation for the second Rm jet.

The accepted first native-y packet supplies source identities for E_y and
V_y.  This module consumes those identities without importing a producer,
loading a saved report, or executing an archived graph operator.  It evaluates
only caller-selected endpoint roots with an exact symbolic dual
(value, d_native).  Reachable operations outside the small reducer become
per-row unresolved gates; no desired value is invented.
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


FIELDS = ("E_y", "V_y")
ORDERS = (0, 1, 2)
SIDES = ("postrestore", "patch")
_ANALYTIC = {"exp", "log", "sqrt", "positive_sqrt", "expm1", "log1p"}
_LAZY_FLAT = {
    "lazy_flat_branch",
    "original_lazy_flat_branch",
    "compact_flat_branch",
    "original_compact_flat_branch",
}
_FUNCTION_GUARDS = (
    "Z_independent_substitution",
    "native_substitution_guard",
    "exact_native_chain_rule_substitution",
    "endpoint_substitution_guard",
    "guarded",
)


class _Unresolved(ValueError):
    """A reachable source operation outside the bounded reducer."""


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


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _handle(value: Any, path: str, count: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{path} must be an integer graph handle")
    if value < 0 or count is not None and value >= count:
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
        try:
            if "/" in text:
                pieces = text.split("/")
                if len(pieces) != 2:
                    raise ValueError
                return Fraction(int(pieces[0]), int(pieces[1]))
            return Fraction(int(text), 1)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{path} is not an exact rational") from exc
    raise ValueError(f"{path} is not an exact rational")


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


def _validate_digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"{label} must be a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueError(f"{label} must be a SHA-256 hex digest") from exc
    return value.lower()


def _exact_value(nodes: Sequence[Any], index: int) -> Fraction | None:
    try:
        row = nodes[index]
        if not isinstance(row, Mapping) or _operation(row, index) != "exact_rational":
            return None
        num = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
        den = _fraction(row.get("denominator", 1), f"node {index}.denominator")
        if den == 0:
            return None
        return num / den
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return None


def _exact_one(nodes: Sequence[Any], index: int) -> bool:
    return _exact_value(nodes, index) == 1


def _guarded_fsub(row: Mapping[str, Any]) -> bool:
    return any(row.get(name) is True for name in _FUNCTION_GUARDS)


def _simplify(expr: Any) -> Any:
    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact endpoint differentiation")
    return _sp.simplify(expr)


def _record(expr: Any) -> dict[str, Any]:
    expr = _simplify(expr)
    return {"srepr": _sp.srepr(expr), "text": str(expr), "is_zero": bool(expr == 0)}


def _atom(name: Any) -> Any:
    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact endpoint differentiation")
    return _sp.Symbol("SOURCE:" + str(name), commutative=True)


@dataclass(frozen=True)
class _Dual:
    value: Any
    derivative: Any


def _add(left: _Dual, right: _Dual) -> _Dual:
    return _Dual(_simplify(left.value + right.value), _simplify(left.derivative + right.derivative))


def _neg(value: _Dual) -> _Dual:
    return _Dual(_simplify(-value.value), _simplify(-value.derivative))


def _mul(left: _Dual, right: _Dual) -> _Dual:
    return _Dual(
        _simplify(left.value * right.value),
        _simplify(left.derivative * right.value + left.value * right.derivative),
    )


def _quotient(numerator: _Dual, denominator: _Dual) -> _Dual:
    if _simplify(denominator.value) == 0:
        raise _Unresolved("positive_quotient denominator is symbolically zero")
    return _Dual(
        _simplify(numerator.value / denominator.value),
        _simplify(
            (numerator.derivative * denominator.value
             - numerator.value * denominator.derivative)
            / denominator.value**2
        ),
    )


def _side_name(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    key = value.lower().replace("-", "_").replace(" ", "_")
    if key in {"postrestore", "post_restore", "post"}:
        return "postrestore"
    if key in {"patch", "finite", "finite_patch", "rm_patch"}:
        return "patch"
    return None


def _side(ports: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    owners: list[Mapping[str, Any]] = []
    if isinstance(ports.get("sides"), Mapping):
        owners.append(ports["sides"])
    owners.append(ports)
    names = (
        ("postrestore", "post_restore", "post")
        if name == "postrestore"
        else ("patch", "finite", "finite_patch", "rm_patch")
    )
    for owner in owners:
        for key in names:
            candidate = owner.get(key)
            if isinstance(candidate, Mapping):
                return candidate
    raise ValueError(f"ports is missing {name} side")


def _port_handle(
    side: Mapping[str, Any],
    names: tuple[str, ...],
    label: str,
    count: int,
) -> int:
    for key in names:
        if key not in side:
            continue
        value = side[key]
        if isinstance(value, Mapping):
            value = value.get("node", value.get("handle", value.get("root")))
        return _handle(value, label, count)
    raise ValueError(f"{label} is missing")


def _rows(value: Any, label: str, count: int) -> list[int]:
    if isinstance(value, Mapping):
        values = []
        for order in ORDERS:
            item = value.get(order, value.get(str(order)))
            if item is None:
                raise ValueError(f"{label} lacks ordinary order {order}")
            values.append(_handle(item, f"{label}[{order}]", count))
        return values
    if not _is_sequence(value) or len(value) != len(ORDERS):
        raise ValueError(f"{label} must contain ordinary orders 0, 1, 2")
    return [_handle(item, f"{label}[{order}]", count) for order, item in enumerate(value)]


def _root_rows(side: Mapping[str, Any], field: str, count: int) -> list[int]:
    aliases = (field, "Ey" if field == "E_y" else "Vy")
    for key in aliases:
        if key in side:
            return _rows(side[key], f"{key} roots", count)
    roots = side.get("roots")
    if isinstance(roots, Mapping):
        for key in aliases:
            if key in roots:
                return _rows(roots[key], f"roots.{key}", count)
    raise ValueError(f"side is missing {field} roots")


def _frozen(side: Mapping[str, Any], ports: Mapping[str, Any], count: int) -> set[int]:
    values: list[Any] = []
    for owner in (ports, side):
        raw = owner.get("frozen_native_constant_nodes")
        if isinstance(raw, Mapping):
            raw = raw.get("nodes", raw.get("handles", []))
        if raw is None:
            continue
        if not _is_sequence(raw):
            raise ValueError("frozen_native_constant_nodes must be a sequence")
        values.extend(raw)
    result: set[int] = set()
    for offset, item in enumerate(values):
        if isinstance(item, Mapping):
            item = item.get("node", item.get("handle"))
        result.add(_handle(item, f"frozen_native_constant_nodes[{offset}]", count))
    return result


def _refs(value: Any, label: str) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, Mapping):
        for key in ("refs", "proof_refs", "source_refs"):
            if key in value:
                return _refs(value[key], label + "." + key)
        return [deepcopy(value)]
    if isinstance(value, (str, int)):
        return [deepcopy(value)]
    if _is_sequence(value):
        return [deepcopy(item) for item in value]
    raise ValueError(f"{label} must be proof-reference data")


def _source_refs(side: Mapping[str, Any], ports: Mapping[str, Any]) -> list[Any]:
    result: list[Any] = []
    for owner in (ports, side):
        for key in ("source_refs", "source_proof_refs", "source_refs_proof"):
            if key in owner:
                result.extend(_refs(owner[key], key))
    return result


def _pair_values(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, Mapping):
        for key in ("pairs", "rows", "identities", "accepted_pairs"):
            if key in value:
                return _pair_values(value[key])
        if any(_side_name(key) is not None for key in value):
            result: list[Any] = []
            for key in value:
                if _side_name(key) is not None:
                    result.extend(_pair_values(value[key]))
            return result
        if any(key in {"E_y", "Ey", "V_y", "Vy"} for key in value):
            result = []
            for field, rows in value.items():
                if field not in {"E_y", "Ey", "V_y", "Vy"} or not isinstance(rows, Mapping):
                    continue
                for order, row in rows.items():
                    if isinstance(row, Mapping):
                        item = dict(row)
                        item.setdefault("field", field)
                        item.setdefault("order", order)
                        result.append(item)
            return result
        return [value]
    if _is_sequence(value):
        result = []
        for item in value:
            result.extend(_pair_values(item))
        return result
    raise ValueError("accepted_Ey_Vy_source_pairs must be a mapping or sequence")


def _pair_field(value: Any) -> str | None:
    if value in {"E_y", "Ey"}:
        return "E_y"
    if value in {"V_y", "Vy"}:
        return "V_y"
    return None


def _pair_node(row: Mapping[str, Any], side: str, count: int) -> int | None:
    if side == "postrestore":
        keys = (
            "postrestore_node", "post_restore_node", "left_node",
            "postrestore_root", "post_restore_root", "left_root",
            "postrestore", "post_restore", "left",
        )
    else:
        keys = (
            "patch_node", "finite_node", "right_node",
            "patch_root", "finite_root", "right_root",
            "patch", "finite", "right",
        )
    for key in keys:
        if key not in row:
            continue
        value = row[key]
        if isinstance(value, Mapping):
            value = value.get("node", value.get("handle", value.get("root")))
        return _handle(value, f"accepted pair {side} node", count)
    return None


def _normalise_pairs(
    ports: Mapping[str, Any], count: int
) -> tuple[dict[tuple[str, int], dict[str, Any]], list[dict[str, Any]]]:
    raw: list[Any] = []
    keys = ("accepted_Ey_Vy_source_pairs", "accepted_E_y_V_y_source_pairs", "accepted_source_pairs")
    for key in keys:
        if key in ports:
            raw.extend(_pair_values(ports[key]))
    for side_name in SIDES:
        try:
            side = _side(ports, side_name)
        except ValueError:
            continue
        for key in keys:
            if key in side:
                raw.extend(_pair_values(side[key]))
    pairs: dict[tuple[str, int], dict[str, Any]] = {}
    issues: list[dict[str, Any]] = []
    for offset, item in enumerate(raw):
        if not isinstance(item, Mapping):
            issues.append({"gate": "accepted_source_pair_malformed", "index": offset})
            continue
        field = _pair_field(item.get("field", item.get("component")))
        order = item.get("order", item.get("ordinary_Z_order", item.get("Z_order")))
        if field is None or isinstance(order, bool) or not isinstance(order, int) or order not in ORDERS:
            issues.append({
                "gate": "accepted_source_pair_malformed",
                "index": offset,
                "error": "pair needs E_y/V_y and ordinary order 0, 1, 2",
            })
            continue
        key = (field, order)
        if key in pairs:
            issues.append({"gate": "accepted_source_pair_duplicate", "field": field, "order": order})
            continue
        proof = item.get("proof_ref", item.get("source_ref", item.get("identity_ref", item.get("proof"))))
        if proof in (None, ""):
            issues.append({"gate": "accepted_source_pair_missing_proof_ref", "field": field, "order": order})
        pairs[key] = {
            "field": field,
            "order": order,
            "postrestore_node": _pair_node(item, "postrestore", count),
            "patch_node": _pair_node(item, "patch", count),
            "proof_ref": deepcopy(proof),
            "metadata": deepcopy(dict(item)),
        }
    return pairs, issues


def _pair_scope(
    nodes: Sequence[Mapping[str, Any]],
    handle: int | None,
    root: int,
    native: int,
    endpoint: int,
) -> dict[str, Any] | None:
    """Validate a source root or its guarded endpoint wrapper.

    The accepted first-radial packet records endpoint function-substitution
    handles.  They are admissible pair descriptors only when their expression
    is this selected pre-endpoint root and their native replacement is exactly
    the caller's endpoint node.  The wrapper itself is not replayed here.
    """

    if handle is None:
        return None
    if handle == root:
        return {
            "kind": "selected_source_root",
            "node": handle,
            "root": root,
            "endpoint_wrapper_validated": False,
        }
    if handle < 0 or handle >= len(nodes):
        return None
    row = nodes[handle]
    if not isinstance(row, Mapping):
        return None
    try:
        if _operation(row, handle) != "function_substitution" or not _guarded_fsub(row):
            return None
        expression = _handle(row.get("expression"), f"pair wrapper {handle}.expression", len(nodes))
        variable = _handle(row.get("variable"), f"pair wrapper {handle}.variable", len(nodes))
        value = _handle(row.get("value"), f"pair wrapper {handle}.value", len(nodes))
    except (ValueError, TypeError):
        return None
    if expression != root or variable != native or value != endpoint:
        return None
    return {
        "kind": "guarded_endpoint_wrapper",
        "node": handle,
        "root": root,
        "expression": expression,
        "variable": variable,
        "endpoint_node": value,
        "guarded": True,
        "endpoint_wrapper_validated": True,
        "guard_keys": [key for key in _FUNCTION_GUARDS if row.get(key) is True],
        "source_ref": deepcopy(row.get("source_ref", row.get("proof_ref"))),
    }


class _Evaluator:
    def __init__(
        self,
        nodes: Sequence[Mapping[str, Any]],
        *,
        side: str,
        native: int,
        endpoint: int,
        frozen: set[int],
        saved_flat_premise: Any = None,
    ) -> None:
        self.nodes = nodes
        self.side = side
        self.native = native
        self.endpoint = endpoint
        self.frozen = set(frozen)
        self.saved_flat_premise = deepcopy(saved_flat_premise)
        self.memo: dict[tuple[int, tuple[Any, ...]], _Dual] = {}
        self.used_nodes: set[int] = set()
        self.used_aliases: list[dict[str, Any]] = []
        self.flat_witness: list[dict[str, Any]] = []

    def _row(self, index: int) -> Mapping[str, Any]:
        _handle(index, "source handle", len(self.nodes))
        row = self.nodes[index]
        if not isinstance(row, Mapping):
            raise _Unresolved(f"node {index} is not a mapping")
        return row

    @staticmethod
    def _signature(env: Mapping[int, _Dual]) -> tuple[Any, ...]:
        return tuple(
            (index, _sp.srepr(_simplify(pair.value)), _sp.srepr(_simplify(pair.derivative)))
            for index, pair in sorted(env.items())
        )

    def _backward(self, index: int, value: Any, label: str) -> int:
        child = _handle(value, f"node {index}.{label}", len(self.nodes))
        if child >= index:
            raise _Unresolved(f"forward reference from node {index}.{label} to node {child}")
        return child

    def _opaque(self, index: int, row: Mapping[str, Any]) -> _Dual:
        explicit_constant = any(
            row.get(name) is True
            for name in ("native_independent", "Z_only", "z_only", "frozen_native_constant")
        )
        if index > self.native and index not in self.frozen and not explicit_constant:
            raise _Unresolved(
                f"unknown native-dependent operation at node {index}: "
                f"{_operation(row, index)!r}"
            )
        # A graph handle identifies a source function.  Names and source IDs
        # are descriptive metadata only and may be reused by distinct leaves.
        # Keep every opaque source atom node-qualified so unrelated source
        # functions cannot be accidentally correlated.
        descriptor = row.get("source_atom", row.get("source_id", row.get("name", f"node_{index}")))
        label = f"SOURCE_NODE_{index}:{descriptor}"
        self.used_nodes.add(index)
        return _Dual(_atom(label), _sp.Integer(0))

    def _saved_flat_row(self, index: int) -> Mapping[str, Any] | None:
        """Find an accepted exact endpoint-flat receipt for one lazy node."""

        def visit(value: Any) -> Mapping[str, Any] | None:
            if isinstance(value, Mapping):
                if (
                    value.get("node") == index
                    and value.get("endpoint_zero") is True
                    and value.get("passed_original_beta_exact_flat_derivative") is True
                    and value.get("endpoint_dependent") is True
                ):
                    return value
                for child in value.values():
                    match = visit(child)
                    if match is not None:
                        return match
            elif _is_sequence(value):
                for child in value:
                    match = visit(child)
                    if match is not None:
                        return match
            return None

        return visit(self.saved_flat_premise)

    def _exact_scalar(self, index: int, env: Mapping[int, _Dual]) -> Any | None:
        """Evaluate only an exact scalar predicate expression.

        This deliberately does not descend through unknown source leaves before
        the selected native coordinate.  Such leaves remain opaque in the
        differentiation reducer; a saved endpoint-flat receipt can discharge
        a predicate whose source value is intentionally opaque.
        """

        _handle(index, "predicate handle", len(self.nodes))
        if index in env:
            value = _simplify(env[index].value)
            return value if value.is_Rational else None
        row = self._row(index)
        operation = _operation(row, index)
        if index < self.native and operation != "exact_rational":
            return None
        if operation == "exact_rational":
            num = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
            den = _fraction(row.get("denominator", 1), f"node {index}.denominator")
            if den == 0:
                return None
            return _sp.Rational(num.numerator, num.denominator) / _sp.Rational(den.numerator, den.denominator)
        if operation == "bound_variable":
            if index == self.native or index == self.endpoint:
                return _sp.Integer(1)
            return None
        if operation in {"sum", "add"}:
            values = [self._exact_scalar(child, env) for child in _args(row, index)]
            return None if any(value is None for value in values) else _simplify(sum(values, _sp.Integer(0)))
        if operation in {"product", "mul"}:
            result = _sp.Integer(1)
            for child in _args(row, index):
                value = self._exact_scalar(child, env)
                if value is None:
                    return None
                result = _simplify(result * value)
            return result
        if operation in {"negative", "neg"}:
            value = self._exact_scalar(_field(row, "argument", "operand"), env)
            return None if value is None else _simplify(-value)
        if operation in {"positive_quotient", "quotient"}:
            numerator = self._exact_scalar(_field(row, "numerator"), env)
            denominator = self._exact_scalar(_field(row, "denominator"), env)
            if numerator is None or denominator in (None, 0):
                return None
            return _simplify(numerator / denominator)
        return None

    def _predicate(self, index: int, env: Mapping[int, _Dual]) -> bool | None:
        """Evaluate a graph boolean at the selected endpoint when decidable."""

        _handle(index, "predicate handle", len(self.nodes))
        row = self._row(index)
        operation = _operation(row, index)
        if operation in {"logical_or", "or"}:
            values = [self._predicate(child, env) for child in _args(row, index)]
            if any(value is True for value in values):
                return True
            return None if any(value is None for value in values) else False
        if operation in {"logical_and", "and"}:
            values = [self._predicate(child, env) for child in _args(row, index)]
            if any(value is False for value in values):
                return False
            return None if any(value is None for value in values) else True
        if operation in {"logical_not", "not"}:
            child = _field(row, "argument", "operand", "value")
            value = self._predicate(child, env)
            return None if value is None else not value
        if operation in {"exact_real_comparison", "comparison", "compare"}:
            left = self._exact_scalar(_field(row, "left"), env)
            right = self._exact_scalar(_field(row, "right"), env)
            if left is None or right is None:
                return None
            operator = row.get("operator", row.get("comparator"))
            try:
                if operator == "<":
                    return bool(left < right)
                if operator == "<=":
                    return bool(left <= right)
                if operator == ">":
                    return bool(left > right)
                if operator == ">=":
                    return bool(left >= right)
                if operator in {"==", "="}:
                    return bool(left == right)
                if operator == "!=":
                    return bool(left != right)
            except TypeError:
                return None
            return None
        if operation in {"exact_boolean", "boolean", "bool"}:
            value = row.get("value", row.get("result"))
            return value if isinstance(value, bool) else None
        return None

    def _lazy(
        self,
        index: int,
        row: Mapping[str, Any],
        env: Mapping[int, _Dual],
    ) -> _Dual:
        predicate = row.get("flat_predicate", row.get("outside_predicate"))
        predicate_value: bool | None = None
        if predicate is not None:
            predicate_handle = self._backward(index, predicate, "flat_predicate")
            predicate_value = self._predicate(predicate_handle, env)
        receipt = self._saved_flat_row(index)
        receipt_matches_predicate = (
            receipt is not None
            and predicate is not None
            and receipt.get("flat_predicate") == predicate
        )
        explicit_outside = row.get(
            "outside_support",
            row.get("exact_outside_support", row.get("outside", False)),
        ) is True
        explicit_guard = any(
            row.get(name) is True
            for name in (
                "compact_Cinf_guard",
                "compact_cinf_support",
                "outside_support_proved",
                "flat_support_proof",
                "exact_zero_outside_support",
            )
        )
        # The branch flags describe the implementation, but do not by
        # themselves establish endpoint support.  Accept an exact predicate,
        # or an accepted endpoint-zero receipt tied to this node and predicate.
        outside = (explicit_outside and explicit_guard) or predicate_value is True or (
            receipt_matches_predicate
            and receipt is not None
            and receipt.get("endpoint_zero") is True
        )
        if not outside:
            if predicate_value is None and receipt is None:
                raise _Unresolved(f"flat branch node {index} has undecidable endpoint predicate")
            active = row.get("active_body", row.get("active", row.get("body")))
            if active is None:
                raise _Unresolved(f"flat branch node {index} lacks active_body inside support")
            child = self._backward(index, active, "active_body")
            self.used_nodes.add(index)
            return self.rewrite(child, env)
        flat = row.get("flat_value", row.get("outside_value", row.get("value", row.get("zero_value"))))
        self.used_nodes.add(index)
        if flat is None:
            if row.get("exact_zero_outside_support") is not True:
                raise _Unresolved(f"outside flat branch node {index} lacks flat_value")
            result = _Dual(_sp.Integer(0), _sp.Integer(0))
            flat_node = None
        else:
            flat_node = self._backward(index, flat, "flat_value")
            result = self.rewrite(flat_node, env)
        self.flat_witness.append({
            "node": index,
            "outside_support": True,
            "guarded": True,
            "outside_support_predicate_evaluated": predicate_value is True,
            "saved_endpoint_flat_receipt": deepcopy(receipt) if receipt_matches_predicate else None,
            "active_ancestry_short_circuited": True,
            "flat_value_node": flat_node,
            "flat_value_exact_zero": bool(
                _simplify(result.value) == 0 and _simplify(result.derivative) == 0
            ),
        })
        return result

    def rewrite(
        self,
        index: int,
        env: Mapping[int, _Dual] | None = None,
    ) -> _Dual:
        _handle(index, "source handle", len(self.nodes))
        environment = dict(env or {})
        key = (index, self._signature(environment))
        if key in self.memo:
            return self.memo[key]
        if index in environment:
            result = environment[index]
            self.memo[key] = result
            return result
        row = self._row(index)
        operation = _operation(row, index)
        # Backward-DAG source nodes preceding the selected native coordinate
        # cannot depend on it.  Preserve them as qualified source constants,
        # while exact rationals still get parsed exactly.
        if index < self.native and operation != "exact_rational":
            result = self._opaque(index, row)
            self.memo[key] = result
            return result

        # Accepted Z-only controls may be arithmetic nodes above the selected
        # native coordinate (for example higher implicit-control jets).  Their
        # frozen receipt is the scope boundary: preserve each as its own source
        # atom instead of expanding an unrelated recurrence.
        if (
            index in self.frozen
            and index not in {self.native, self.endpoint}
            and operation != "exact_rational"
        ):
            result = self._opaque(index, row)
            self.memo[key] = result
            return result

        if operation in _LAZY_FLAT:
            result = self._lazy(index, row, environment)
            self.memo[key] = result
            return result

        self.used_nodes.add(index)
        if operation == "exact_rational":
            num = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
            den = _fraction(row.get("denominator", 1), f"node {index}.denominator")
            if den == 0:
                raise _Unresolved(f"node {index} has a zero rational denominator")
            result = _Dual(
                _sp.Rational(num.numerator, num.denominator)
                / _sp.Rational(den.numerator, den.denominator),
                _sp.Integer(0),
            )
        elif operation == "bound_variable":
            if index == self.native:
                result = _Dual(_sp.Integer(1), _sp.Integer(1))
            elif index == self.endpoint:
                if not _exact_one(self.nodes, index):
                    raise _Unresolved("endpoint node is not exact rational one")
                result = _Dual(_sp.Integer(1), _sp.Integer(0))
            elif index in self.frozen or any(
                row.get(name) is True
                for name in ("native_independent", "Z_only", "z_only", "frozen_native_constant")
            ):
                result = self._opaque(index, row)
            else:
                raise _Unresolved(f"bound variable node {index} is not frozen")
        elif operation in {"sum", "add"}:
            result = _Dual(_sp.Integer(0), _sp.Integer(0))
            for offset, value in enumerate(_args(row, index)):
                child = self._backward(index, value, f"arguments[{offset}]")
                result = _add(result, self.rewrite(child, environment))
        elif operation in {"product", "mul"}:
            result = _Dual(_sp.Integer(1), _sp.Integer(0))
            for offset, value in enumerate(_args(row, index)):
                child = self._backward(index, value, f"arguments[{offset}]")
                result = _mul(result, self.rewrite(child, environment))
        elif operation in {"negative", "neg"}:
            child = self._backward(index, _field(row, "argument", "operand"), "argument")
            result = _neg(self.rewrite(child, environment))
        elif operation in {"positive_quotient", "quotient"}:
            numerator = self._backward(index, _field(row, "numerator"), "numerator")
            denominator = self._backward(index, _field(row, "denominator"), "denominator")
            positive = (
                row.get("source_positive_certificate") is not None
                or row.get("positive") is True
                or row.get("denominator_positive") is True
                or _exact_value(self.nodes, denominator) not in (None, 0)
            )
            if not positive:
                raise _Unresolved(f"node {index} positive_quotient lacks positivity certificate")
            result = _quotient(
                self.rewrite(numerator, environment),
                self.rewrite(denominator, environment),
            )
        elif operation in {"compact_raw_beta", "raw_compact_beta", "compact_beta"}:
            argument = self._backward(index, _field(row, "argument", "operand", "value"), "argument")
            child = self.rewrite(argument, environment)
            argument_value = _simplify(child.value)
            definition = str(row.get("exact_definition", "")).lower()
            has_compact_support = (
                row.get("compact_Cinf_guard") is True
                or row.get("compact_cinf_support") is True
                or "exact0 elsewhere" in definition
                or "zero elsewhere" in definition
            )
            if not has_compact_support:
                raise _Unresolved(f"compact_raw_beta node {index} lacks compact-support definition")
            if not getattr(argument_value, "is_Rational", False):
                raise _Unresolved(
                    f"compact_raw_beta node {index} has no exact endpoint argument outside support"
                )
            if abs(argument_value) < 1:
                raise _Unresolved(
                    f"compact_raw_beta node {index} is inside support; source body not reduced here"
                )
            # The accepted compact C-infinity definition is exactly zero,
            # including its endpoint derivatives, for |argument| >= 1.
            result = _Dual(_sp.Integer(0), _sp.Integer(0))
        elif operation == "function_substitution":
            if not _guarded_fsub(row):
                raise _Unresolved(f"node {index} function_substitution lacks a guard")
            expression = self._backward(index, _field(row, "expression"), "expression")
            variable = self._backward(index, _field(row, "variable"), "variable")
            value = self._backward(index, _field(row, "value"), "value")
            replacement = self.rewrite(value, environment)
            scoped = dict(environment)
            scoped[variable] = replacement
            self.used_aliases.append({
                "node": index,
                "expression": expression,
                "variable": variable,
                "endpoint_node": value,
                "guarded": True,
                "source_ref": deepcopy(row.get("source_ref", row.get("proof_ref"))),
            })
            result = self.rewrite(expression, scoped)
        elif operation == "analytic_unary" or operation in _ANALYTIC:
            name = operation
            if operation == "analytic_unary":
                name = row.get("name", row.get("function", row.get("unary", row.get("kind"))))
            if name not in _ANALYTIC:
                raise _Unresolved(f"unsupported analytic unary {name!r}")
            argument = self._backward(
                index, _field(row, "argument", "operand", "value"), "argument"
            )
            child = self.rewrite(argument, environment)
            if name == "exp":
                value = _sp.exp(child.value)
                derivative = value * child.derivative
            elif name == "log":
                value = _sp.log(child.value)
                derivative = child.derivative / child.value
            elif name in {"sqrt", "positive_sqrt"}:
                value = _sp.sqrt(child.value)
                derivative = child.derivative / (2 * value)
            elif name == "expm1":
                value = _sp.exp(child.value) - 1
                derivative = _sp.exp(child.value) * child.derivative
            else:
                value = _sp.log(1 + child.value)
                derivative = child.derivative / (1 + child.value)
            result = _Dual(_simplify(value), _simplify(derivative))
        else:
            explicit_constant = any(
                row.get(name) is True
                for name in ("native_independent", "Z_only", "z_only", "frozen_native_constant")
            )
            if explicit_constant or index in self.frozen:
                result = self._opaque(index, row)
            else:
                raise _Unresolved(
                    f"unsupported native-dependent operation at node {index}: {operation!r}"
                )
        result = _Dual(_simplify(result.value), _simplify(result.derivative))
        self.memo[key] = result
        return result


def _unique(values: Sequence[Any]) -> list[Any]:
    seen: set[str] = set()
    result: list[Any] = []
    for value in values:
        marker = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
        if marker not in seen:
            seen.add(marker)
            result.append(deepcopy(value))
    return result


def _row_result(
    side: str,
    field: str,
    order: int,
    root: int,
    dual: _Dual,
    factor: Any,
) -> dict[str, Any]:
    physical = _simplify(factor * dual.derivative)
    residual = _simplify(physical - dual.value / 10) if field == "E_y" else _simplify(physical)
    return {
        "side": side,
        "field": field,
        "order": order,
        "ordinary_Z_order": order,
        "root_node": root,
        "status": "proved" if residual == 0 else "failed",
        "identity": (
            "endpoint(D_y E_y - E_y/10)=0"
            if field == "E_y"
            else "endpoint(D_y V_y)=0"
        ),
        "value": _record(dual.value),
        "native_derivative": _record(dual.derivative),
        "physical_derivative": _record(physical),
        "residual": _record(residual),
        "physical_factor_at_endpoint": _record(factor),
        "physical_factor_exact_one": bool(_simplify(factor) == 1),
    }


def prove(
    nodes: Sequence[Mapping[str, Any]],
    ports: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    accepted_seed_digest: str,
) -> dict[str, Any]:
    """Consume endpoint source roots and prove the bounded second-radial rows."""

    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact endpoint differentiation")
    if not _is_sequence(nodes):
        raise ValueError("nodes must be a sequence")
    if not isinstance(ports, Mapping):
        raise ValueError("ports must be a mapping")
    expected = _validate_digest(expected_nodes_digest, "expected_nodes_digest")
    accepted = _validate_digest(accepted_seed_digest, "accepted_seed_digest")
    actual = canonical_digest(list(nodes))
    if actual != expected:
        raise ValueError("nodes digest does not match expected_nodes_digest")
    original_nodes = deepcopy(nodes)
    original_ports = deepcopy(ports)

    configs: dict[str, dict[str, Any]] = {}
    for side_name in SIDES:
        side = _side(ports, side_name)
        native = _port_handle(
            side, ("native_variable", "native_endpoint_variable"),
            f"{side_name}.native_variable", len(nodes)
        )
        endpoint = _port_handle(
            side, ("endpoint_node1", "endpoint_node", "endpoint"),
            f"{side_name}.endpoint_node1", len(nodes)
        )
        if not _exact_one(nodes, endpoint):
            raise ValueError(f"{side_name}.endpoint_node1 is not exact rational one")
        configs[side_name] = {
            "side": side,
            "native": native,
            "endpoint": endpoint,
            "roots": {field: _root_rows(side, field, len(nodes)) for field in FIELDS},
            "frozen": _frozen(side, ports, len(nodes)),
            "saved_flat_premise": deepcopy(
                side.get("saved_flat_premise", ports.get("saved_flat_premise"))
            ),
        }

    pairs, pair_issues = _normalise_pairs(ports, len(nodes))
    refs: list[Any] = []
    for config in configs.values():
        refs.extend(_source_refs(config["side"], ports))
        refs.extend(_refs(
            config["side"].get("saved_flat_premise", ports.get("saved_flat_premise")),
            "saved_flat_premise",
        ))
    refs.extend(pair["proof_ref"] for pair in pairs.values())
    refs = _unique([ref for ref in refs if ref not in (None, "")])

    evaluators: dict[str, _Evaluator] = {}
    final_rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = list(pair_issues)
    for side_name in SIDES:
        config = configs[side_name]
        evaluator = _Evaluator(
            nodes,
            side=side_name,
            native=config["native"],
            endpoint=config["endpoint"],
            frozen=config["frozen"],
            saved_flat_premise=config["saved_flat_premise"],
        )
        evaluators[side_name] = evaluator
        # Post-restore has D_y p = 1.  Patch has D_y x = x; x = 1 at
        # endpoint, so both physical factors are exactly one.
        factor = _sp.Integer(1)
        for field in FIELDS:
            for order, root in enumerate(config["roots"][field]):
                try:
                    dual = evaluator.rewrite(
                        root,
                        {config["native"]: _Dual(_sp.Integer(1), _sp.Integer(1))},
                    )
                    row = _row_result(side_name, field, order, root, dual, factor)
                    final_rows.append(row)
                    if row["status"] != "proved":
                        unresolved.append({
                            "gate": "endpoint_identity_mismatch",
                            "side": side_name,
                            "field": field,
                            "order": order,
                            "root_node": root,
                            "residual": row["residual"],
                        })
                except (ValueError, _Unresolved) as exc:
                    final_rows.append({
                        "side": side_name,
                        "field": field,
                        "order": order,
                        "ordinary_Z_order": order,
                        "root_node": root,
                        "status": "unresolved",
                        "identity": (
                            "endpoint(D_y E_y - E_y/10)=0"
                            if field == "E_y"
                            else "endpoint(D_y V_y)=0"
                        ),
                    })
                    unresolved.append({
                        "gate": "source_reduction_unresolved",
                        "side": side_name,
                        "field": field,
                        "order": order,
                        "root_node": root,
                        "error": str(exc),
                    })

    seam: list[dict[str, Any]] = []
    for field in FIELDS:
        for order in ORDERS:
            key = (field, order)
            expected_left = configs["postrestore"]["roots"][field][order]
            expected_right = configs["patch"]["roots"][field][order]
            pair = pairs.get(key)
            if pair is None:
                seam.append({
                    "field": field,
                    "order": order,
                    "status": "unresolved",
                    "gate": "missing_accepted_source_pair",
                    "postrestore_node": expected_left,
                    "patch_node": expected_right,
                })
                unresolved.append({"gate": "missing_accepted_source_pair", "field": field, "order": order})
                continue
            left, right, proof = pair["postrestore_node"], pair["patch_node"], pair["proof_ref"]
            left_scope = _pair_scope(
                nodes, left, expected_left,
                configs["postrestore"]["native"], configs["postrestore"]["endpoint"],
            )
            right_scope = _pair_scope(
                nodes, right, expected_right,
                configs["patch"]["native"], configs["patch"]["endpoint"],
            )
            if left_scope is None or right_scope is None:
                seam.append({
                    "field": field,
                    "order": order,
                    "status": "unresolved",
                    "gate": "accepted_source_pair_root_mismatch",
                    "postrestore_node": left,
                    "patch_node": right,
                    "expected_postrestore_node": expected_left,
                    "expected_patch_node": expected_right,
                    "postrestore_scope": left_scope,
                    "patch_scope": right_scope,
                    "proof_ref": deepcopy(proof),
                })
                unresolved.append({"gate": "accepted_source_pair_root_mismatch", "field": field, "order": order})
                continue
            if proof in (None, ""):
                seam.append({
                    "field": field,
                    "order": order,
                    "status": "unresolved",
                    "gate": "accepted_source_pair_missing_proof_ref",
                    "postrestore_node": left,
                    "patch_node": right,
                })
                unresolved.append({
                    "gate": "accepted_source_pair_missing_proof_ref",
                    "field": field,
                    "order": order,
                })
                continue
            seam.append({
                "field": field,
                "order": order,
                "ordinary_Z_order": order,
                "status": "proved_by_accepted_source_pair",
                "postrestore_node": left,
                "patch_node": right,
                "proof_ref": deepcopy(proof),
                "pair_identity_consumed": True,
                "independent_pair_reduction": False,
                "postrestore_scope": left_scope,
                "patch_scope": right_scope,
                "endpoint_wrappers_validated": bool(
                    left_scope["endpoint_wrapper_validated"]
                    and right_scope["endpoint_wrapper_validated"]
                ),
            })

    row_pass = len(final_rows) == 12 and all(row.get("status") == "proved" for row in final_rows)
    seam_pass = len(seam) == 6 and all(row.get("status") == "proved_by_accepted_source_pair" for row in seam)
    used = {side: sorted(ev.used_nodes) for side, ev in evaluators.items()}
    aliases = {side: deepcopy(ev.used_aliases) for side, ev in evaluators.items()}
    flat = {side: deepcopy(ev.flat_witness) for side, ev in evaluators.items()}
    port_meta = {
        side: {
            "native_variable": config["native"],
            "endpoint_node1": config["endpoint"],
            "E_y_roots": list(config["roots"]["E_y"]),
            "V_y_roots": list(config["roots"]["V_y"]),
            "frozen_native_constant_nodes": sorted(config["frozen"]),
            "saved_flat_premise": deepcopy(
                config["side"].get("saved_flat_premise", ports.get("saved_flat_premise"))
            ),
            "source_refs": _source_refs(config["side"], ports),
        }
        for side, config in configs.items()
    }
    result = {
        "schema": "original_Rm_second_radial_identity_v1",
        "passed": bool(row_pass and seam_pass and not unresolved),
        "source_only": True,
        "sourceonly": True,
        "actual_saved_graph_run": False,
        "actual_saved_graph_loaded": False,
        "accepted_builders_or_owners_replayed": False,
        "accepted_checkers_replayed": False,
        "accepted_reducers_or_providers_replayed": False,
        "quadratures_replayed": False,
        "history_archive_decoded": False,
        "physical_J_exact_one": True,
        "physicalJexact1": True,
        "ordinary_orders": list(ORDERS),
        "ordinaryorders": list(ORDERS),
        "expected_nodes_digest": expected,
        "nodes_digest": actual,
        "accepted_seed_digest": accepted,
        "module_sha256": _module_sha256(),
        "ports": port_meta,
        "accepted_seed_refs": refs,
        "acceptedseedrefs": deepcopy(refs),
        "final_rows": final_rows,
        "rows": deepcopy(final_rows),
        "second_radial_seam": seam,
        "unresolved_gates": unresolved,
        "used_nodes": used,
        "usednodes": deepcopy(used),
        "used_aliases": aliases,
        "aliases": deepcopy(aliases),
        "flat_branch_witness": flat,
        "flatbranch_witness": deepcopy(flat),
        "checkpoint": {
            "final_rows": deepcopy(final_rows),
            "unresolved_gates": deepcopy(unresolved),
            "checkpoint_before_inspection": True,
        },
        "scope": {
            "source_graph_only": True,
            "no_saved_graph_execution": True,
            "no_history_archive_decode": True,
            "no_endpoint_numeric_oracle": True,
            "no_old_normalization_replayed": True,
            "backward_DAG_enforced": True,
            "unknown_native_dependent_ops_fail_closed": True,
            "lazy_flat_guard_short_circuited_before_active_ancestry": True,
            "physical_endpoint_J_is_exact_one": True,
            "full_germ_claimed": False,
            "finite_N_germ_proved": False,
            "global_second_radial_identity_claimed": False,
            "metadata_preserved": deepcopy(nodes) == original_nodes and deepcopy(ports) == original_ports,
            "nodes_input_unmutated": deepcopy(nodes) == original_nodes,
            "ports_input_unmutated": deepcopy(ports) == original_ports,
        },
        "limits": {
            "source_endpoint_rows_only": True,
            "physical_derivative_is_slope_at_patch_x1": True,
            "accepted_pair_identities_consumed_as_premises": True,
            "full_physical_germ": False,
            "full_finite_N_history": False,
            "five_moment_closure": False,
        },
    }
    json.dumps(result, allow_nan=False)
    return result


__all__ = ["FIELDS", "ORDERS", "SIDES", "canonical_digest", "prove"]
