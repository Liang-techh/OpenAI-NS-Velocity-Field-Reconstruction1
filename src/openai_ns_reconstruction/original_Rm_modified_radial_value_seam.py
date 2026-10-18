"""Source-only VALUE congruence for the modified total radial E/V jets.

The slow native VALUE seam is an input to this module.  This module never
reconstructs or differentiates that seam.  It treats the admitted slow target
handles, generic VALUE handles, loop VALUE handles, and raw physical radial
handles as formal equality atoms and reduces only the arithmetic and
exponential constructors of ``total_radial_functions``.

The public entry point is :func:`derive`.  It is intentionally stricter than
the source graph evaluator used by older consumers: an unknown operation or a
missing source packet is an unresolved gate, and no derivative DAG is opened
to recover it.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_modified_radial_value_seam_v1"
ORDERS = (0, 1)
SIDES = ("postrestore", "patch")
GRAPH_NODE_COUNT = 205300
DERIVATIVE_PREFIX = 205147
Z_HANDLE = 20


# These are source ports, not numerical values.  Every one is checked against
# the caller's immutable graph before it can become a formal atom or target.
GENERIC_PORTS = {
    ("E", 0): (84911, 194957),
    ("E", 1): (84914, 194960),
    ("V", 0): (84945, 194993),
    ("V", 1): (84946, 194994),
}
LOOP_PORTS = {
    ("A", 0): (121070, 197890),
    ("A", 1): (121071, 197891),
    ("A_phi", 0): (121189, 198009),
    ("A_phi", 1): (121190, 198010),
    ("B_over_Pstar_phi", 0): (121291, 198123),
    ("B_over_Pstar_phi", 1): (121292, 198124),
}
SLOW_TARGET_PORTS = {
    ("A", 0): (178388, 199051),
    ("A", 1): (178710, 199421),
    ("B", 0): (182270, 203304),
    ("B", 1): (182304, 203341),
}
RAW_BODY_PORTS = {
    ("E_y", 0): (85146, 195164),
    ("E_y", 1): (85149, 195165),
    ("V_y", 0): (0, 195171),
    ("V_y", 1): (0, 195172),
}
TOTAL_TARGET_PORTS = {
    ("E_y", 0): (182725, 203774),
    ("E_y", 1): (182728, 203777),
    ("V_y", 0): (182757, 203806),
    ("V_y", 1): (182758, 203807),
}
WRAPPER_PORTS = {
    ("E_y", 0): (205147, 195850),
    ("E_y", 1): (205148, 195851),
    ("V_y", 0): (205150, 195857),
    ("V_y", 1): (205151, 195858),
}
RAW_FORM_SHA256 = {
    ("E_y", 0): "49fcff180ba65fdf49f0dac9f5cac89f01556226954a858bfe960a4db00db7ab",
    ("E_y", 1): "b764f5c67f59618459afaebe43a3a954ca81b912cee92af15637ce462b500da3",
    ("V_y", 0): "a2f2a2a871fc6c22108aeb363ce5d0ff5ddc2e4dfe0d56f5f73f2c9df22316c0",
    ("V_y", 1): "a2f2a2a871fc6c22108aeb363ce5d0ff5ddc2e4dfe0d56f5f73f2c9df22316c0",
}
WRAPPER_GUARD_KEYS = ["Z_independent_substitution"]


class ValueSeamError(ValueError):
    """Malformed admission data or a source VALUE closure failure."""


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueSeamError(message)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def canonical_digest(value: Any) -> str:
    """Return the stable digest used by source-only packets."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, f"{label} must be a SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueSeamError(f"{label} must be a SHA-256 digest") from exc
    return value.lower()


def _handle(value: Any, label: str, count: int) -> int:
    if isinstance(value, Mapping):
        value = value.get("node", value.get("handle", value.get("root")))
    _require(isinstance(value, int) and not isinstance(value, bool), f"{label} must be an integer handle")
    _require(0 <= value < count, f"{label} is outside the source graph")
    return value


def _row(nodes: Sequence[Mapping[str, Any]], index: int, label: str) -> Mapping[str, Any]:
    _require(0 <= index < len(nodes), f"{label} {index} is outside the source graph")
    value = nodes[index]
    _require(isinstance(value, Mapping), f"{label} {index} is not a mapping")
    return value


def _operation(nodes: Sequence[Mapping[str, Any]], index: int) -> str:
    value = _row(nodes, index, "source node")
    operation = value.get("operation", value.get("op"))
    _require(isinstance(operation, str) and operation, f"source node {index} has no operation")
    return operation


def _child(
    nodes: Sequence[Mapping[str, Any]], index: int, value: Any, label: str
) -> int:
    child = _handle(value, f"source node {index}.{label}", len(nodes))
    _require(child < index, f"source node {index}.{label} is a forward reference")
    _require(child < DERIVATIVE_PREFIX, f"source node {index}.{label} enters the derivative suffix")
    return child


def _exact_zero(nodes: Sequence[Mapping[str, Any]], index: int) -> bool:
    row = _row(nodes, index, "source node")
    if _operation(nodes, index) != "exact_rational":
        return False
    return row.get("numerator", row.get("value")) == 0 and row.get("denominator", 1) == 1


def _form(value: Any, label: str) -> dict[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a formal VALUE mapping")
    srepr = value.get("srepr")
    is_zero = value.get("is_zero")
    _require(isinstance(srepr, str) and srepr, f"{label}.srepr must be nonempty text")
    _require(isinstance(is_zero, bool), f"{label}.is_zero must be boolean")
    if is_zero:
        _require(srepr in {"Integer(0)", "0"}, f"{label} claims a non-exact zero")
    else:
        _require(srepr not in {"Integer(0)", "0"}, f"{label} has an unadmitted zero form")
    return {"srepr": srepr, "is_zero": is_zero}


def _owner_form(owner: Any, label: str) -> dict[str, Any]:
    _require(isinstance(owner, Mapping), f"{label} must be a mapping")
    value = owner.get("value", owner.get("form"))
    if value is None and "srepr" in owner:
        value = owner
    return _form(value, f"{label}.value")


def _side_owner(row: Mapping[str, Any], side: str, label: str) -> Mapping[str, Any]:
    owner = row.get(side)
    _require(isinstance(owner, Mapping), f"{label} lacks {side}")
    return owner


def _source_form_digest(form: Mapping[str, Any]) -> str:
    return canonical_digest({"srepr": form["srepr"], "is_zero": form["is_zero"]})


def _raw_form(owner: Any, label: str, form_sha: Any) -> tuple[dict[str, Any], str]:
    _require(isinstance(owner, Mapping), f"{label} must be a raw common-form mapping")
    terms = owner.get("terms")
    _require(_is_sequence(terms), f"{label}.terms must be a sequence")
    form_sha = _digest(form_sha, f"{label}.sha256")
    normalized_terms = deepcopy(list(terms))
    is_zero = len(normalized_terms) == 0
    if is_zero:
        form = {"srepr": "Integer(0)", "is_zero": True}
    else:
        form = {
            "srepr": json.dumps({"terms": normalized_terms}, sort_keys=True, separators=(",", ":")),
            "is_zero": False,
        }
    return form, form_sha


Formal = tuple[Any, ...]


def _zero() -> Formal:
    return ("zero",)


def _one() -> Formal:
    return ("rational", 1, 1)


def _rat(numerator: int, denominator: int = 1) -> Formal:
    _require(denominator != 0, "a formal rational has a zero denominator")
    value = Fraction(int(numerator), int(denominator))
    if value == 0:
        return _zero()
    return ("rational", value.numerator, value.denominator)


def _is_zero_expr(value: Formal) -> bool:
    return value == _zero()


def _is_one_expr(value: Formal) -> bool:
    return value == _one()


def _formal(op: str, *children: Formal) -> Formal:
    if op == "add":
        flattened: list[Formal] = []
        for child in children:
            if child[0] == "add":
                flattened.extend(child[1])
            elif not _is_zero_expr(child):
                flattened.append(child)
        if not flattened:
            return _zero()
        if len(flattened) == 1:
            return flattened[0]
        return ("add", tuple(sorted(flattened, key=repr)))
    if op == "mul":
        flattened = []
        for child in children:
            if _is_zero_expr(child):
                return _zero()
            if child[0] == "mul":
                flattened.extend(child[1])
            elif not _is_one_expr(child):
                flattened.append(child)
        if not flattened:
            return _one()
        if len(flattened) == 1:
            return flattened[0]
        return ("mul", tuple(sorted(flattened, key=repr)))
    if op == "exp":
        return ("exp", children[0])
    if op == "div":
        numerator, denominator = children
        _require(not _is_zero_expr(denominator), "a formal quotient has an exact zero denominator")
        if _is_zero_expr(numerator):
            return _zero()
        if _is_one_expr(denominator):
            return numerator
        return ("div", numerator, denominator)
    if op == "pow":
        base, exponent = children
        if exponent == _zero():
            return _one()
        return ("pow", base, exponent)
    raise ValueSeamError(f"unsupported formal constructor {op!r}")


def _formal_json(value: Formal) -> Any:
    if value[0] in {"zero", "rational", "atom"}:
        return list(value)
    if value[0] in {"exp"}:
        return [value[0], _formal_json(value[1])]
    if value[0] == "add" or value[0] == "mul":
        return [value[0], [_formal_json(item) for item in value[1]]]
    if value[0] in {"div", "pow"}:
        return [value[0], *(_formal_json(item) for item in value[1:])]
    raise ValueSeamError(f"unknown formal expression tag {value[0]!r}")


def _formal_text(value: Formal) -> str:
    tag = value[0]
    if tag == "zero":
        return "0"
    if tag == "rational":
        return str(value[1]) if value[2] == 1 else f"{value[1]}/{value[2]}"
    if tag == "atom":
        return f"{value[1]}:{value[2]}"
    if tag == "add":
        return " + ".join(_formal_text(item) for item in value[1])
    if tag == "mul":
        return " * ".join(_formal_text(item) for item in value[1])
    if tag == "div":
        return f"({_formal_text(value[1])})/({_formal_text(value[2])})"
    if tag == "pow":
        return f"({_formal_text(value[1])})^({_formal_text(value[2])})"
    if tag == "exp":
        return f"exp({_formal_text(value[1])})"
    raise ValueSeamError(f"unknown formal expression tag {tag!r}")


def _record(value: Formal) -> dict[str, Any]:
    payload = _formal_json(value)
    return {
        "srepr": json.dumps(payload, sort_keys=True, separators=(",", ":")),
        "text": _formal_text(value),
        "is_zero": _is_zero_expr(value),
        "formal_digest": canonical_digest(payload),
    }


@dataclass(frozen=True)
class _Alias:
    side: str
    node: int
    role: str
    namespace: str
    key: str
    form_digest: str
    value: Formal
    support: str
    guard: Any


class _AliasTable:
    def __init__(self, nodes: Sequence[Mapping[str, Any]]) -> None:
        self.nodes = nodes
        self.by_side: dict[str, dict[int, _Alias]] = {side: {} for side in SIDES}
        self.rows: list[dict[str, Any]] = []

    def add(self, alias: _Alias) -> None:
        _require(alias.side in SIDES, f"unknown alias side {alias.side!r}")
        _require(0 <= alias.node < DERIVATIVE_PREFIX, f"alias {alias.role} enters the derivative suffix")
        _require(isinstance(alias.guard, Mapping), f"alias {alias.role} lacks a guard mapping")
        if "ordinary_Z_order" in alias.guard:
            _require(alias.guard.get("ordinary_Z_order") in ORDERS,
                     f"alias {alias.role} has an invalid ordinary Z order guard")
        _require(alias.guard.get("Z_handle") == Z_HANDLE,
                 f"alias {alias.role} is not guarded by the admitted Z handle")
        _require(isinstance(alias.support, str) and alias.support.startswith("inherited:"),
                 f"alias {alias.role} lacks an inherited support certificate")
        if alias.node == 0:
            _require(_exact_zero(self.nodes, alias.node), "node 0 alias is not the source graph's exact zero")
            _require(_is_zero_expr(alias.value), f"node 0 alias {alias.role} is not exact zero")
        current = self.by_side[alias.side].get(alias.node)
        if current is not None:
            _require(current.value == alias.value, f"conflicting alias roles at {alias.side} node {alias.node}")
            _require(current.form_digest == alias.form_digest, f"conflicting alias forms at {alias.side} node {alias.node}")
            self.rows.append({"side": alias.side, "node": alias.node, "role": alias.role, "collision": "same_exact_form"})
            return
        self.by_side[alias.side][alias.node] = alias
        self.rows.append({
            "side": alias.side,
            "node": alias.node,
            "role": alias.role,
            "namespace": alias.namespace,
            "key": alias.key,
            "form_digest": alias.form_digest,
            "support": alias.support,
            "guard": deepcopy(alias.guard),
            "zero": _is_zero_expr(alias.value),
        })


class _FormalNormalizer:
    """Small, source-DAG-only congruence normalizer.

    Aliases are checked before the node operation.  Therefore a slow physical
    target is a formal VALUE atom and its native derivative ancestry is never
    traversed.
    """

    def __init__(self, nodes: Sequence[Mapping[str, Any]], aliases: Mapping[int, _Alias], side: str) -> None:
        self.nodes = nodes
        self.aliases = aliases
        self.side = side
        self.cache: dict[int, Formal] = {}
        self.active: set[int] = set()
        self.visited: set[int] = set()
        self.operations: Counter[str] = Counter()

    def rewrite(self, index: int) -> Formal:
        index = _handle(index, f"{self.side} target", len(self.nodes))
        alias = self.aliases.get(index)
        if alias is not None:
            return alias.value
        _require(index < DERIVATIVE_PREFIX, f"{self.side} VALUE target entered the derivative suffix at {index}")
        if index in self.cache:
            return self.cache[index]
        _require(index not in self.active, f"cyclic VALUE dependency at {self.side} node {index}")
        self.active.add(index)
        self.visited.add(index)
        operation = _operation(self.nodes, index)
        self.operations[operation] += 1
        row = _row(self.nodes, index, "source node")
        try:
            if operation == "exact_rational":
                numerator = row.get("numerator", row.get("value"))
                denominator = row.get("denominator", 1)
                _require(isinstance(numerator, int) and isinstance(denominator, int), f"source rational {index} is malformed")
                result = _rat(numerator, denominator)
            elif operation == "exact_positive_integer_expression":
                _require(row.get("exact_positive_integer") is True and row.get("Z_independent") is True,
                         f"source integer {index} lacks its exact Z-independent contract")
                result = ("atom", "SOURCE_CONSTANT", f"node_{index}")
            elif operation in {"original_source_parameter", "current_original_source_parameter"}:
                _require(row.get("Z_independent") is True, f"source parameter {index} lacks Z-independent support")
                result = ("atom", "SOURCE_PARAMETER", f"node_{index}")
            elif operation == "bound_variable":
                _require(row.get("name") == "Z", f"unadmitted bound variable at VALUE node {index}")
                result = ("atom", "BOUND_VARIABLE", "Z")
            elif operation in {"sum", "add"}:
                children = row.get("arguments", row.get("operands"))
                _require(_is_sequence(children) and children, f"source sum {index} has no arguments")
                result = _formal("add", *(self.rewrite(_child(self.nodes, index, child, "arguments")) for child in children))
            elif operation in {"product", "mul"}:
                children = row.get("arguments", row.get("operands"))
                _require(_is_sequence(children) and children, f"source product {index} has no arguments")
                result = _formal("mul", *(self.rewrite(_child(self.nodes, index, child, "arguments")) for child in children))
            elif operation in {"negative", "neg"}:
                result = _formal("mul", _rat(-1), self.rewrite(_child(self.nodes, index, row.get("argument"), "argument")))
            elif operation in {"difference", "sub", "subtract"}:
                left = self.rewrite(_child(self.nodes, index, row.get("left"), "left"))
                right = self.rewrite(_child(self.nodes, index, row.get("right"), "right"))
                result = _formal("add", left, _formal("mul", _rat(-1), right))
            elif operation in {"positive_quotient", "quotient"}:
                if operation == "positive_quotient":
                    certificate = row.get("source_positive_certificate", row.get("positive_certificate", row.get("denominator_positive")))
                    _require(certificate not in (None, "", False), f"source quotient {index} lacks its positivity guard")
                numerator = self.rewrite(_child(self.nodes, index, row.get("numerator"), "numerator"))
                denominator = self.rewrite(_child(self.nodes, index, row.get("denominator"), "denominator"))
                result = _formal("div", numerator, denominator)
            elif operation in {"integer_power", "power"}:
                base = self.rewrite(_child(self.nodes, index, row.get("base", row.get("argument")), "base"))
                exponent_raw = row.get("exponent", row.get("power"))
                if isinstance(exponent_raw, int):
                    exponent_index = _child(self.nodes, index, exponent_raw, "exponent")
                    exponent_row = _row(self.nodes, exponent_index, "power exponent")
                    _require(_operation(self.nodes, exponent_index) == "exact_rational", f"power {index} has a non-exact exponent")
                    exponent = _rat(exponent_row.get("numerator", exponent_row.get("value")), exponent_row.get("denominator", 1))
                else:
                    exponent = self.rewrite(_child(self.nodes, index, exponent_raw, "exponent"))
                result = _formal("pow", base, exponent)
            elif operation == "analytic_unary":
                _require(row.get("name") == "exp", f"unsupported analytic VALUE operation {row.get('name')!r}")
                result = _formal("exp", self.rewrite(_child(self.nodes, index, row.get("argument"), "argument")))
            else:
                raise ValueSeamError(f"unknown or unadmitted VALUE operation at node {index}: {operation!r}")
            self.cache[index] = result
            return result
        finally:
            self.active.remove(index)

    def closure(self) -> dict[str, Any]:
        used = sorted(self.visited)
        return {
            "used_node_count": len(used),
            "used_node_digest": canonical_digest(used),
            "operation_counts": {name: self.operations[name] for name in sorted(self.operations)},
            "aliases_cut_before_rewrite": True,
        }


def _packet_digest(packet: Mapping[str, Any]) -> str:
    return canonical_digest(packet)


def _validate_generic(
    packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]], aliases: _AliasTable
) -> dict[str, Any]:
    _require(packet.get("schema") == "original_Rm_generic_native_identity_v1", "generic packet has an unknown schema")
    _require(packet.get("passed") is True, "generic packet is not passed")
    _require(packet.get("source_only") is True, "generic packet is not source-only")
    _require(packet.get("leading_only") is True, "generic packet is not a leading VALUE packet")
    _require(packet.get("ordinary_orders") == [0, 1], "generic packet ordinary orders changed")
    rows = packet.get("root_rows")
    _require(_is_sequence(rows) and len(rows) == 18, "generic packet must contain 18 root rows")
    seen: set[tuple[str, int]] = set()
    selected = {key for key in GENERIC_PORTS}
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"generic row {offset} is malformed")
        root = row.get("root", row.get("name"))
        order = row.get("order", row.get("ordinary_Z_order"))
        _require(isinstance(root, str) and isinstance(order, int) and order in ORDERS, f"generic row {offset} has invalid root/order")
        key = (root, order)
        _require(key not in seen, f"duplicate generic root row {root}[{order}]")
        seen.add(key)
        _require(row.get("status") == "proved", f"generic root row {root}[{order}] is not proved")
        left = _side_owner(row, "postrestore", f"generic {root}[{order}]")
        right = _side_owner(row, "patch", f"generic {root}[{order}]")
        left_node = _handle(left.get("root_node", left.get("node")), f"generic postrestore {root}[{order}]", len(nodes))
        right_node = _handle(right.get("root_node", right.get("node")), f"generic patch {root}[{order}]", len(nodes))
        if key in selected:
            _require((left_node, right_node) == GENERIC_PORTS[key], f"generic source ports changed for {root}[{order}]")
            left_form = _owner_form(left, f"generic postrestore {root}[{order}]")
            right_form = _owner_form(right, f"generic patch {root}[{order}]")
            _require(left_form == right_form, f"generic VALUE forms differ for {root}[{order}]")
            value = _zero() if left_form["is_zero"] else ("atom", "GENERIC_VALUE", f"{root}[{order}]")
            for side, node in (("postrestore", left_node), ("patch", right_node)):
                aliases.add(_Alias(side, node, f"generic:{root}[{order}]", "GENERIC_VALUE", f"{root}[{order}]",
                                   _source_form_digest(left_form), value, "inherited:passed_generic_packet", {
                                       "ordinary_Z_order": order,
                                       "Z_handle": Z_HANDLE,
                                       "separate_ordinary_Z_atom": True,
                                       "packet_schema": packet.get("schema"),
                                   }))
    _require(seen == {(root, order) for root in packet.get("root_names", ("E", "V", "a", "b", "p1", "p2", "t0", "kappa", "Delta")) for order in ORDERS},
             "generic packet root rows do not cover the admitted roots")
    return {"schema": packet.get("schema"), "digest": _packet_digest(packet), "selected_rows": sorted(selected)}


def _validate_loop(
    packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]], aliases: _AliasTable
) -> dict[str, Any]:
    _require(packet.get("schema") == "admitted_Rm_loop_source_ports_v1", "loop packet has an unknown schema")
    _require(packet.get("passed") is True and packet.get("source_only") is True, "loop packet is not passed source data")
    _require(packet.get("mathematical_operator_executed") is False, "loop packet replayed a mathematical operator")
    _require(packet.get("numerical_values_computed") is False, "loop packet computed numerical values")
    _require(packet.get("ordinary_orders") == [0, 1], "loop packet ordinary orders changed")
    rows = packet.get("rows")
    _require(_is_sequence(rows) and len(rows) == 22, "loop packet must contain 22 VALUE rows")
    common = packet.get("common_phase_source")
    _require(isinstance(common, Mapping) and common.get("passed") is True and common.get("source_only") is True,
             "loop packet lacks its passed common phase source")
    _require(common.get("actual_phase_Z_independence_proved") is True, "loop phase is not Z-independent")
    _require(common.get("fractional_part_differentiated") is False, "loop phase fractional part was differentiated")
    _require(common.get("endpoints") == [1, 1], "loop phase endpoint contract changed")
    seen: set[tuple[str, int]] = set()
    selected = set(LOOP_PORTS)
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"loop row {offset} is malformed")
        field = row.get("field")
        order = row.get("order", row.get("ordinary_Z_order"))
        _require(isinstance(field, str) and isinstance(order, int) and order in ORDERS, f"loop row {offset} has invalid field/order")
        key = (field, order)
        _require(key not in seen, f"duplicate loop row {field}[{order}]")
        seen.add(key)
        _require(row.get("status") == "proved", f"loop row {field}[{order}] is not proved")
        if key not in selected:
            continue
        left = _side_owner(row, "postrestore", f"loop {field}[{order}]")
        right = _side_owner(row, "patch", f"loop {field}[{order}]")
        left_node = _handle(left.get("node", left.get("root_node")), f"loop postrestore {field}[{order}].node", len(nodes))
        right_node = _handle(right.get("node", right.get("root_node")), f"loop patch {field}[{order}].node", len(nodes))
        _require((left_node, right_node) == LOOP_PORTS[key], f"loop source ports changed for {field}[{order}]")
        left_form = _owner_form(left, f"loop postrestore {field}[{order}]")
        right_form = _owner_form(right, f"loop patch {field}[{order}]")
        _require(left_form == right_form, f"loop VALUE forms differ for {field}[{order}]")
        value = _zero() if left_form["is_zero"] else ("atom", "LOOP_VALUE", f"{field}[{order}]")
        for side, node in (("postrestore", left_node), ("patch", right_node)):
            aliases.add(_Alias(side, node, f"loop:{field}[{order}]", "LOOP_VALUE", f"{field}[{order}]",
                               _source_form_digest(left_form), value, "inherited:passed_loop_packet", {
                                   "ordinary_Z_order": order,
                                   "Z_handle": Z_HANDLE,
                                   "separate_ordinary_Z_atom": True,
                                   "common_phase_Z_independent": True,
                               }))
    _require(seen == {(row.get("field"), row.get("order", row.get("ordinary_Z_order"))) for row in rows},
             "loop packet row identity changed")
    _require(selected <= seen, "loop packet lacks a selected VALUE source row")
    return {"schema": packet.get("schema"), "digest": _packet_digest(packet), "selected_rows": sorted(selected)}


def _slow_value_record(row: Mapping[str, Any], side: str, label: str) -> dict[str, Any]:
    owner = _side_owner(row, side, label)
    return _owner_form(owner, label)


def _validate_slow(
    packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]], aliases: _AliasTable, expected_nodes_digest: str
) -> dict[str, Any]:
    _require(packet.get("schema") == "original_Rm_slow_native_value_seam_v1", "slow VALUE packet has an unknown schema")
    _require(packet.get("passed") is True and packet.get("source_only") is True, "slow VALUE packet is not passed source data")
    _require(packet.get("numeric_values_evaluated") is False, "slow VALUE packet computed numerical values")
    _require(packet.get("old_operators_replayed") is False, "slow VALUE packet replayed an old operator")
    if "new_native_derivative_generated" in packet:
        _require(packet.get("new_native_derivative_generated") is False, "slow VALUE packet generated a derivative")
    if "expected_nodes_digest" in packet:
        _require(_digest(packet.get("expected_nodes_digest"), "slow_value_seam.expected_nodes_digest") == expected_nodes_digest,
                 "slow VALUE packet is bound to another source graph")
    if "nodes_digest" in packet:
        _require(_digest(packet.get("nodes_digest"), "slow_value_seam.nodes_digest") == expected_nodes_digest,
                 "slow VALUE packet graph digest does not match the source graph")
    rows = packet.get("canonical_value_pairs")
    _require(_is_sequence(rows) and len(rows) == 4, "slow VALUE packet must contain four canonical value pairs")
    seen: set[tuple[str, int]] = set()
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"slow VALUE row {offset} is malformed")
        field = row.get("field")
        order = row.get("ordinary_order", row.get("order"))
        _require(field in {"A", "B"} and isinstance(order, int) and order in ORDERS,
                 f"slow VALUE row {offset} has invalid field/order")
        key = (field, order)
        _require(key not in seen, f"duplicate slow VALUE row {field}[{order}]")
        seen.add(key)
        _require(row.get("equal") is True and row.get("status") == "proved",
                 f"slow VALUE row {field}[{order}] is not an admitted equality")
        left_target = _handle(row.get("postrestore_target_handle"), f"slow postrestore target {field}[{order}]", len(nodes))
        right_target = _handle(row.get("patch_target_handle"), f"slow patch target {field}[{order}]", len(nodes))
        _require((left_target, right_target) == SLOW_TARGET_PORTS[key], f"slow target ports changed for {field}[{order}]")
        left_form = _slow_value_record(row, "postrestore", f"slow {field}[{order}]")
        right_form = _slow_value_record(row, "patch", f"slow {field}[{order}]")
        _require(left_form == right_form, f"slow VALUE forms differ for {field}[{order}]")
        _require(not left_form["is_zero"], f"slow physical target {field}[{order}] was admitted as zero")
        atom = ("atom", "SLOW_TARGET_VALUE", f"{field}[{order}]")
        form_digest = _source_form_digest(left_form)
        for side, node in (("postrestore", left_target), ("patch", right_target)):
            aliases.add(_Alias(side, node, f"slow_target:{field}[{order}]", "SLOW_TARGET_VALUE", f"{field}[{order}]",
                               form_digest, atom, "inherited:passed_slow_value_packet", {
                                   "ordinary_Z_order": order,
                                   "Z_handle": Z_HANDLE,
                                   "separate_ordinary_Z_atom": True,
                                   "value_only": True,
                                   "derivative_DAG_expanded": False,
                               }))
    _require(seen == set(SLOW_TARGET_PORTS), "slow VALUE packet does not cover A/B ordinary orders")
    return {"schema": packet.get("schema"), "digest": _packet_digest(packet), "rows": sorted(seen)}


def _validate_raw_and_wrappers(
    raw_packet: Mapping[str, Any], second_packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]], aliases: _AliasTable
) -> dict[str, Any]:
    _require(raw_packet.get("schema") == "Rm_raw_native_jet_identity_v1", "raw radial packet has an unknown schema")
    _require(raw_packet.get("passed") is True and raw_packet.get("source_only") is True, "raw radial packet is not passed source data")
    _require(raw_packet.get("leading_raw_native_jet_only") is True, "raw radial packet scope changed")
    _require(raw_packet.get("modified_finite_N_germ") is False and raw_packet.get("RQy_identity_proved") is False,
             "raw radial packet overclaims its scope")
    source_graph = raw_packet.get("source_graph")
    _require(isinstance(source_graph, Mapping), "raw radial source graph metadata is missing")
    _require(source_graph.get("modified_wrapper") == "modified_Q4_pressure_C3_bindings"
             and source_graph.get("finite_wrapper") == "Rm_patch_finite_N_C5_bindings",
             "raw radial wrapper bindings changed")
    _require(source_graph.get("postrestore_endpoint_variable") == 84528
             and source_graph.get("finite_endpoint_variable") == 185676
             and source_graph.get("endpoint_node_is_exact_rational_one") is True,
             "raw radial endpoint ports changed")
    _require(second_packet.get("schema") == "original_Rm_second_radial_identity_v1", "endpoint packet has an unknown schema")
    _require(second_packet.get("passed") is True and second_packet.get("source_only") is True,
             "endpoint packet is not passed source data")
    _require(second_packet.get("physical_J_exact_one") is True, "endpoint packet lost exact J one")
    raw_digest = raw_packet.get("expected_nodes_digest")
    second_digest = second_packet.get("expected_nodes_digest")
    if raw_digest is not None and second_digest is not None:
        _require(_digest(raw_digest, "raw_packet.expected_nodes_digest") == _digest(second_digest, "second_packet.expected_nodes_digest"),
                 "raw and endpoint packets are bound to different saved graphs")

    form_rows: dict[tuple[str, int], Mapping[str, Any]] = {}
    for row in raw_packet.get("row_forms", ()):
        _require(isinstance(row, Mapping), "raw row form is malformed")
        field = row.get("field")
        order = row.get("order", row.get("ordinary_Z_order"))
        if field in {"E", "V"} and order in ORDERS:
            form_rows[("E_y" if field == "E" else "V_y", order)] = row
    _require(set(form_rows) == set(RAW_BODY_PORTS), "raw packet lacks E_y/V_y Z 0/1 common forms")
    wrapper_rows: dict[tuple[str, int], Mapping[str, Any]] = {}
    seam = second_packet.get("second_radial_seam")
    _require(_is_sequence(seam), "endpoint packet lacks second_radial_seam rows")
    for row in seam:
        if not isinstance(row, Mapping):
            continue
        key = (row.get("field"), row.get("ordinary_Z_order", row.get("order")))
        if key in WRAPPER_PORTS:
            wrapper_rows[key] = row
    _require(set(wrapper_rows) == set(WRAPPER_PORTS), "endpoint packet lacks E_y/V_y Z 0/1 wrappers")

    for key, (left_root, right_root) in RAW_BODY_PORTS.items():
        form_row = form_rows[key]
        _require(form_row.get("status") == "proved", f"raw common form {key} is not proved")
        _require(form_row.get("left_node") == left_root and form_row.get("right_node") == right_root,
                 f"raw body ports changed for {key}")
        left_form, left_sha = _raw_form(
            form_row.get("left_form"), f"raw left {key}", form_row.get("left_form_sha256")
        )
        right_form, right_sha = _raw_form(
            form_row.get("right_form"), f"raw right {key}", form_row.get("right_form_sha256")
        )
        _require(left_sha == right_sha == RAW_FORM_SHA256[key], f"raw common form digest changed for {key}")
        _require(left_form == right_form, f"raw common forms differ for {key}")
        if key[0] == "E_y":
            _require(left_form["srepr"] == json.dumps({"terms": [{"coefficient": "1", "monomial": [f"FUNCTION:E_y[{key[1]}]"]}]}, sort_keys=True, separators=(",", ":")),
                     f"raw E_y common form changed for {key}")
            value = ("atom", "RAW_PHYSICAL_DY", f"E_y[{key[1]}]")
        else:
            _require(left_form["is_zero"] is True, f"raw V_y common form is not exact zero for {key}")
            value = _zero()
        wrapper = wrapper_rows[key]
        _require(wrapper.get("status") == "proved_by_accepted_source_pair"
                 and wrapper.get("pair_identity_consumed") is True
                 and wrapper.get("independent_pair_reduction") is False
                 and wrapper.get("endpoint_wrappers_validated") is True,
                 f"raw endpoint wrapper {key} is not an accepted pair")
        for side, root, wrapper_key, variable in (
            ("postrestore", left_root, "postrestore_scope", 84528),
            ("patch", right_root, "patch_scope", 185676),
        ):
            scope = wrapper.get(wrapper_key)
            _require(isinstance(scope, Mapping), f"raw {side} wrapper {key} is missing")
            _require(scope.get("kind") == "guarded_endpoint_wrapper"
                     and scope.get("node") == WRAPPER_PORTS[key][0 if side == "postrestore" else 1]
                     and scope.get("root") == root
                     and scope.get("expression") == root
                     and scope.get("variable") == variable
                     and scope.get("endpoint_node") == 1
                     and scope.get("guarded") is True
                     and scope.get("endpoint_wrapper_validated") is True
                     and scope.get("guard_keys") == WRAPPER_GUARD_KEYS
                     and scope.get("source_ref") is None,
                     f"raw {side} endpoint wrapper metadata changed for {key}")
            aliases.add(_Alias(side, root, f"raw_physical_D_y:{key[0]}[{key[1]}]", "RAW_PHYSICAL_DY",
                               f"{key[0]}[{key[1]}]", left_sha, value,
                               "inherited:raw_native_and_second_radial_packets", {
                                   "kind": "guarded_endpoint_wrapper",
                                   "guard_keys": list(WRAPPER_GUARD_KEYS),
                                   "Z_handle": Z_HANDLE,
                                   "ordinary_Z_order": key[1],
                                   "variable": variable,
                                   "endpoint_node": 1,
                               }))
    return {
        "raw_schema": raw_packet.get("schema"),
        "raw_digest": _packet_digest(raw_packet),
        "endpoint_schema": second_packet.get("schema"),
        "endpoint_digest": _packet_digest(second_packet),
        "wrapper_rows": len(wrapper_rows),
    }


def _validate_graph(nodes: Sequence[Mapping[str, Any]], expected_nodes_digest: str) -> tuple[str, str]:
    _require(_is_sequence(nodes), "nodes must be a sequence")
    _require(len(nodes) == GRAPH_NODE_COUNT, f"source graph must contain {GRAPH_NODE_COUNT} nodes")
    _require(all(isinstance(row, Mapping) for row in nodes), "nodes must contain mappings")
    expected = _digest(expected_nodes_digest, "expected_nodes_digest")
    actual = canonical_digest(list(nodes))
    _require(actual == expected, "nodes digest does not match expected_nodes_digest")
    prefix = canonical_digest(list(nodes[:DERIVATIVE_PREFIX]))
    return actual, prefix


def derive(
    nodes: Sequence[Mapping[str, Any]],
    generic_packet: Mapping[str, Any],
    loop_ports_admission: Mapping[str, Any],
    slow_value_seam: Mapping[str, Any],
    raw_native_identity: Mapping[str, Any],
    second_radial_identity: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Prove the four total E_y/V_y VALUE congruences.

    ``slow_value_seam`` must already be a passed value-only receipt.  The
    caller owns authentication and reconstruction of ``nodes``; this function
    only verifies the supplied graph digest, source handles, packet rows,
    guards, and exact source operations reached by the eight targets.
    """

    actual_digest, prefix_digest = _validate_graph(nodes, expected_nodes_digest)
    _require(_operation(nodes, Z_HANDLE) == "bound_variable" and _row(nodes, Z_HANDLE, "Z").get("name") == "Z",
             "graph Z handle is not the admitted Z variable")
    aliases = _AliasTable(nodes)
    packet_info = {
        "generic": _validate_generic(generic_packet, nodes, aliases),
        "loop": _validate_loop(loop_ports_admission, nodes, aliases),
        "slow": _validate_slow(slow_value_seam, nodes, aliases, actual_digest),
        "raw": _validate_raw_and_wrappers(raw_native_identity, second_radial_identity, nodes, aliases),
    }

    pair_rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    closures: dict[str, dict[str, Any]] = {}
    for field in ("E_y", "V_y"):
        for order in ORDERS:
            key = (field, order)
            left_target, right_target = TOTAL_TARGET_PORTS[key]
            left = _FormalNormalizer(nodes, aliases.by_side["postrestore"], "postrestore")
            right = _FormalNormalizer(nodes, aliases.by_side["patch"], "patch")
            left_error: str | None = None
            right_error: str | None = None
            try:
                left_value = left.rewrite(left_target)
            except ValueSeamError as exc:
                left_value = _zero()
                left_error = str(exc)
            try:
                right_value = right.rewrite(right_target)
            except ValueSeamError as exc:
                right_value = _zero()
                right_error = str(exc)
            closures[f"postrestore:{field}[{order}]"] = left.closure()
            closures[f"patch:{field}[{order}]"] = right.closure()
            if left_error is not None or right_error is not None:
                equal = False
                status = "unresolved"
                unresolved.append({
                    "gate": "value_normalization_unresolved",
                    "field": field,
                    "ordinary_order": order,
                    "postrestore_error": left_error,
                    "patch_error": right_error,
                })
            else:
                equal = left_value == right_value
                status = "proved" if equal else "failed"
                if not equal:
                    unresolved.append({
                        "gate": "canonical_value_mismatch",
                        "field": field,
                        "ordinary_order": order,
                        "postrestore": _record(left_value),
                        "patch": _record(right_value),
                    })
            pair_rows.append({
                "field": field,
                "ordinary_order": order,
                "postrestore_target_handle": left_target,
                "patch_target_handle": right_target,
                "postrestore": _record(left_value),
                "patch": _record(right_value),
                "equal": equal,
                "status": status,
            })

    passed = len(pair_rows) == 4 and not unresolved and all(row["equal"] is True for row in pair_rows)
    return {
        "schema": SCHEMA,
        "passed": bool(passed),
        "status": "proved" if passed else "unresolved",
        "source_only": True,
        "source_executed": False,
        "actual_saved_graph_run": False,
        "numeric_values_evaluated": False,
        "numeric_field_evaluated": False,
        "new_native_derivative_generated": False,
        "old_operators_replayed": False,
        "nodes_input_unmutated": canonical_digest(list(nodes)) == actual_digest,
        "expected_nodes_digest": actual_digest,
        "nodes_digest": actual_digest,
        "graph_node_count": len(nodes),
        "derivative_lookup_prefix": {"count": DERIVATIVE_PREFIX, "digest": prefix_digest},
        "canonical_value_pairs": pair_rows,
        "value_pairs": deepcopy(pair_rows),
        "alias_roles": deepcopy(aliases.rows),
        "alias_count": len(aliases.rows),
        "alias_collision_count": sum(1 for row in aliases.rows if row.get("collision") == "same_exact_form"),
        "packet_bindings": packet_info,
        "closures": closures,
        "unresolved_gates": unresolved,
        "input_data_digests": {
            "generic_packet": _packet_digest(generic_packet),
            "loop_ports_admission": _packet_digest(loop_ports_admission),
            "slow_value_seam": _packet_digest(slow_value_seam),
            "raw_native_identity": _packet_digest(raw_native_identity),
            "second_radial_identity": _packet_digest(second_radial_identity),
        },
        "target_ports": {
            "E_y": {str(order): list(TOTAL_TARGET_PORTS[("E_y", order)]) for order in ORDERS},
            "V_y": {str(order): list(TOTAL_TARGET_PORTS[("V_y", order)]) for order in ORDERS},
        },
        "source_formula": {
            "E_y": "exp(A/N) * (background_E_y + E * (slow_A_y/N + A_phi))",
            "V_y": "background_V_y + slow_B_y/N + B_over_Pstar_phi",
            "slow_targets_are_formal_atoms": True,
            "raw_D_y_atoms_are_formal_functions": True,
            "separate_ordinary_Z_atoms": True,
        },
        "scope": {
            "four_total_radial_targets_only": True,
            "ordinary_orders": [0, 1],
            "slow_value_packet_required": True,
            "slow_derivative_DAG_expanded": False,
            "derivative_generation": False,
            "physical_Q_y_claimed": False,
            "Z3_claimed": False,
            "global_C1_claimed": False,
            "Navier_Stokes_claimed": False,
        },
    }


__all__ = [
    "DERIVATIVE_PREFIX",
    "GRAPH_NODE_COUNT",
    "SCHEMA",
    "ValueSeamError",
    "canonical_digest",
    "derive",
]
