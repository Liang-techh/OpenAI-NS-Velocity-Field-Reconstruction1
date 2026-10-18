"""Exact constructor-DAG congruence for the four slow native VALUE targets.

This consumer is deliberately independent of the older slow VALUE seam.  It
consumes source-only generic VALUE/native-D admissions, the bounded common
phase admission, and an explicit paired source-to-D-to-target plan.  It does
not create a derivative, invoke a source builder, evaluate a number, or use a
symbolic algebra system.  Reachable source nodes are reduced by a small
hash-consed constructor language.  Unknown operations and missing guards are
reported as unresolved gates.

The public entry point is :func:`derive`.  The result keeps the historical
``original_Rm_slow_native_value_seam_v1`` schema so downstream VALUE consumers
can inspect the four rows, while ``proof_method`` records this independent
constructor-congruence proof.
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


SCHEMA = "original_Rm_slow_native_value_seam_v1"
PROOF_METHOD = "exact_constructor_congruence"
ORDERS = (0, 1)
SIDES = ("postrestore", "patch")
GENERIC_ROOTS = ("E", "V", "a", "b", "p1", "p2", "t0", "kappa", "Delta")
SLOW_FIELDS = ("A", "B")
GRAPH_NODE_COUNT = 205300
DERIVATIVE_PREFIX = 205147
MAX_CONSTRUCTOR_NODES = 50000
EXPECTED_DERIVATIVE_PREFIX_DIGEST = "c7ab7da8cc5c51cd782a5794e3076a6cf021a45235aae1c75bd246b01491b639"
Z_HANDLE = 20
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_PHASE_ROWS = (
    {
        "coordinate": 84528,
        "endpoint": 1,
        "radius": 84533,
        "global_offset": 119414,
        "actual_phase": 119416,
        "common_N": 18,
        "common_Ra": 65604,
        "common_phase_origin": 24954,
    },
    {
        "coordinate": 185676,
        "endpoint": 1,
        "radius": 187214,
        "global_offset": 196044,
        "actual_phase": 196046,
        "common_N": 18,
        "common_Ra": 65604,
        "common_phase_origin": 24954,
    },
)
EXPECTED_PROVENANCE_PORTS = {
    "postrestore": {
        ("A", 0): (121070, 178387, 178388),
        ("A", 1): (121071, 178709, 178710),
        ("B", 0): (121108, 182269, 182270),
        ("B", 1): (121109, 182303, 182304),
    },
    "patch": {
        ("A", 0): (197890, 199050, 199051),
        ("A", 1): (197891, 199420, 199421),
        ("B", 0): (197928, 203303, 203304),
        ("B", 1): (197929, 203340, 203341),
    },
}
EXPECTED_GENERIC_PORTS = {
    ("E", 0): (84911, 194957),
    ("E", 1): (84914, 194960),
    ("V", 0): (84945, 194993),
    ("V", 1): (84946, 194994),
    ("a", 0): (119435, 196071),
    ("a", 1): (119439, 196075),
    ("b", 0): (0, 196106),
    ("b", 1): (0, 196110),
    ("p1", 0): (119488, 196141),
    ("p1", 1): (119492, 196145),
    ("p2", 0): (97163, 195759),
    ("p2", 1): (97164, 195763),
    ("t0", 0): (0, 196211),
    ("t0", 1): (0, 196212),
    ("kappa", 0): (119435, 196278),
    ("kappa", 1): (119439, 196279),
    ("Delta", 0): (119523, 196284),
    ("Delta", 1): (119439, 196279),
}
EXPECTED_LOOP_PORTS = {
    ("A", 0): (121070, 197890),
    ("A", 1): (121071, 197891),
    ("A_phi", 0): (121189, 198009),
    ("A_phi", 1): (121190, 198010),
    ("B_over_Pstar", 0): (121108, 197928),
    ("B_over_Pstar", 1): (121109, 197929),
    ("B_over_Pstar_phi", 0): (121291, 198123),
    ("B_over_Pstar_phi", 1): (121292, 198124),
    ("E_N", 0): (121560, 198412),
    ("E_N", 1): (121562, 198414),
    ("V_N", 0): (121571, 198423),
    ("V_N", 1): (121572, 198424),
    ("m", 0): (121350, 198182),
    ("m", 1): (121351, 198183),
    ("h", 0): (121324, 198156),
    ("h", 1): (121327, 198159),
    ("k", 0): (121424, 198266),
    ("k", 1): (121425, 198267),
    ("e", 0): (121554, 198406),
    ("e", 1): (121555, 198407),
    ("p", 0): (121542, 198394),
    ("p", 1): (121543, 198395),
}
EXPECTED_PHASE_NODES = {"postrestore": 119416, "patch": 196046}
EXPECTED_COORDINATES = {"postrestore": 84528, "patch": 185676}
EXPECTED_INVERSE_J = {"postrestore": 172066, "patch": 198869}
EXPECTED_PHYSICAL_J = {"postrestore": 1, "patch": 187215}
EXPECTED_VISITED_COUNTS = {"postrestore": 365, "patch": 410}
EXPECTED_CONTROL_SYSTEM = 193470
EXPECTED_CONTROL_NODES = (193471, 193472, 193473, 193474, 193475)
EXPECTED_NATIVE_ENDPOINTS = {"postrestore": 1, "patch": 1}
EXPECTED_RADIUS_LEMMA_DIGEST = "a719562a3503be0d7b4672633adb88f7f16ed587638d58c0d5e027d9bc406b24"
EXPECTED_ENDPOINT_RADIUS_REDUCTION = {
    "left": {"closed_node": 85743, "rule": "guarded R_left(p=1)=Rm"},
    "right": {"closed_node": 187963, "body_node": 187214, "rule": "(Rm*x)|x=1=Rm*1=Rm"},
    "common_radius_source_node": 85743,
}

_HANDLE_FIELDS = (
    "argument", "operand", "base", "body", "active_body", "flat_predicate",
    "exponent", "power", "angle_variable_node",
    "phase_function", "target_phase", "inverse_angle", "query_coordinate", "free_query",
    "native_coordinate", "coordinate", "left", "right", "expression", "variable", "value",
    "integrand", "lower", "upper", "numerator", "denominator", "arguments", "operands",
    "Delta", "delta", "eta", "eta21", "original_collar_coordinate", "original_collar_upper",
    "flat_value", "rho", "Z_variable", "system", "outer", "dx", "direct", "chain",
)
_FUNCTION_GUARDS = (
    "Z_independent_substitution", "native_substitution_guard", "exact_native_chain_rule_substitution",
    "endpoint_substitution_guard", "guarded", "inverse_substitution_guard",
)
_ANALYTIC = {"exp", "expm1", "log", "log1p", "sqrt", "positive_sqrt", "sin", "cos"}
_OPAQUE_ANALYTIC = {
    "fractional_part", "gamma", "exprel", "log1prel", "original_flat_sigma",
    "original_flat_sigma_prime", "original_flat_sigma_second", "original_flat_sigma_third",
    "compact_raw_beta_derivative",
}
_COMPOUND_OPS = {
    "original_monotone_phase_inverse", "original_lazy_flat_branch", "exact_original_source_flat_sigma",
    "substitute_original_inverse_angle", "substitute_original_native_coordinate",
    "substitute_free_query_at_native_coordinate", "inverse_substitution", "native_substitution",
    "logical_or", "logical_and", "logical_not", "trim", "exact_real_comparison",
    "regular_even_radial_axis_continuation", "implicit_C1_Z_linear_system",
    "implicit_C2_ZZ_linear_system", "implicit_C3_ZZZ_linear_system", "source_bound_C1_Banach_limit",
    "C1_Banach_limit_component", "C2_Banach_limit_second_component", "C3_Banach_limit_third_component",
    "current_original_axis_pressure_function_derivative", "current_original_leading_function_recipe",
    "current_original_physical_radius", "current_original_physical_source_quantity",
    "current_original_physical_source_scale", "current_original_source_parameter",
    "current_original_unique_analytic_core_component", "current_original_unique_regular_n1_component",
    "original_Rm_unique_leading_implicit_control", "original_Rm_unique_leading_implicit_system",
    "original_compact_unique_Banach_limit_component", "original_compact_unique_Banach_limit_system",
    "proved_C1_zero_identity", "proved_C2_relative_terminal_zero_identity", "proved_C2_second_zero_identity",
    "proved_C3_relative_terminal_zero_identity", "proved_C3_third_zero_identity",
    "proved_original_compact_relative_C5_terminal_identity", "proved_original_n1_C3_transported_patch_integral_identity",
    "proved_original_quiet_linear_mass_source_identity", "compact_raw_beta", "compact_raw_beta_derivative",
}
_SOURCE_CONSTANT_OPS = {
    "exact_positive_integer_expression", "mathematical_pi", "mathematical_positive_infinity",
    "exact_original_source_heat_phi",
}


class ValueSeamError(ValueError):
    """Malformed admission data or an unresolved structural obligation."""


class _Unresolved(ValueSeamError):
    """A reachable operation outside the bounded exact constructor language."""


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueSeamError(message)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _json_safe(value: Any) -> Any:
    if isinstance(value, Fraction):
        return {"__fraction__": [value.numerator, value.denominator]}
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def canonical_digest(value: Any) -> str:
    raw = json.dumps(_json_safe(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")
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
    _require(isinstance(value, int) and not isinstance(value, bool), f"{label} must be an integer graph handle")
    _require(0 <= value < count, f"{label} is outside the source graph")
    return value


def _row(nodes: Sequence[Mapping[str, Any]], index: int, label: str = "source node") -> Mapping[str, Any]:
    _require(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < len(nodes), f"{label} {index!r} is outside the source graph")
    value = nodes[index]
    _require(isinstance(value, Mapping), f"{label} {index} is not a mapping")
    return value


def _operation(nodes: Sequence[Mapping[str, Any]], index: int) -> str:
    value = _row(nodes, index)
    operation = value.get("operation", value.get("op"))
    _require(isinstance(operation, str) and operation, f"source node {index} has no operation")
    return operation


def _child(nodes: Sequence[Mapping[str, Any]], index: int, value: Any, label: str) -> int:
    child = _handle(value, f"node {index}.{label}", len(nodes))
    _require(child < index, f"forward source reference from node {index}.{label} to node {child}")
    _require(child < DERIVATIVE_PREFIX, f"node {index}.{label} enters the derivative suffix")
    return child


def _field(row: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    raise _Unresolved("missing " + "/".join(names))


def _exact_zero(nodes: Sequence[Mapping[str, Any]], index: int) -> bool:
    row = _row(nodes, index)
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


# The canonical term language.  Terms are immutable tuples and are therefore
# hash-consed by Python's dictionaries without expanding rational polynomials.
Term = tuple[Any, ...]


def _zero() -> Term:
    return ("zero",)


def _one() -> Term:
    return ("rational", 1, 1)


def _rat(numerator: int | Fraction, denominator: int = 1) -> Term:
    _require(denominator != 0, "exact rational has a zero denominator")
    value = Fraction(numerator, denominator)
    if value == 0:
        return _zero()
    return ("rational", value.numerator, value.denominator)


def _is_zero(value: Any) -> bool:
    return value == _zero()


def _is_one(value: Any) -> bool:
    return value == _one()


def _as_fraction(value: Any) -> Fraction | None:
    if isinstance(value, tuple) and value and value[0] == "rational":
        return Fraction(value[1], value[2])
    return None


_TERM_TAGS = frozenset({
    "zero", "rational", "atom", "add", "mul", "div", "pow", "func", "compound",
    "integral", "inverse", "branch", "comparison",
})
_TERM_DIGEST_CACHE: dict[int, tuple[tuple[Any, ...], str]] = {}


def _term_digest_value(value: Any) -> Any:
    if isinstance(value, tuple):
        if value and isinstance(value[0], str) and value[0] in _TERM_TAGS:
            return ["term", _term_digest(value)]
        return ["tuple", [_term_digest_value(item) for item in value]]
    if isinstance(value, list):
        return ["list", [_term_digest_value(item) for item in value]]
    return value


def _term_digest(value: Any) -> str:
    """Return a shallow-child digest for an interned constructor DAG node.

    The cache retains each tuple strongly, so object identifiers cannot be
    reused while a digest is live.  Repeated children therefore contribute
    only their cached digest and are never serialized as repeated expression
    trees.
    """

    if not isinstance(value, tuple) or not value or not isinstance(value[0], str) or value[0] not in _TERM_TAGS:
        return canonical_digest(_term_digest_value(value))
    cached = _TERM_DIGEST_CACHE.get(id(value))
    if cached is not None and cached[0] is value:
        return cached[1]
    payload = ["TERM", value[0], *(_term_digest_value(item) for item in value[1:])]
    digest = canonical_digest(payload)
    _TERM_DIGEST_CACHE[id(value)] = (value, digest)
    return digest


def _term_key(value: Any) -> str:
    return _term_digest(value)


def _split_coefficient(value: Term) -> tuple[Fraction, Term | None]:
    direct = _as_fraction(value)
    if direct is not None:
        return direct, None
    if value and value[0] == "mul":
        coefficient = Fraction(1)
        factors: list[Term] = []
        for factor in value[1]:
            numeric = _as_fraction(factor)
            if numeric is None:
                factors.append(factor)
            else:
                coefficient *= numeric
        if not factors:
            return coefficient, None
        return coefficient, _make_mul(*factors)
    return Fraction(1), value


def _make_add(*children: Term) -> Term:
    flat: list[Term] = []
    for child in children:
        if _is_zero(child):
            continue
        if isinstance(child, tuple) and child and child[0] == "add":
            flat.extend(child[1])
        else:
            flat.append(child)
    if not flat:
        return _zero()
    coefficients: dict[str, tuple[Fraction, Term | None]] = {}
    factors: dict[str, Term | None] = {}
    for child in flat:
        coefficient, factor = _split_coefficient(child)
        key = "__constant__" if factor is None else _term_key(factor)
        factors[key] = factor
        previous = coefficients.get(key, (Fraction(0), factor))[0]
        coefficients[key] = (previous + coefficient, factor)
    terms: list[Term] = []
    for key, (coefficient, factor) in coefficients.items():
        if coefficient == 0:
            continue
        term = _rat(coefficient)
        if factor is not None:
            term = _make_mul(term, factor)
        terms.append(term)
    if not terms:
        return _zero()
    if len(terms) == 1:
        return terms[0]
    return ("add", tuple(sorted(terms, key=_term_key)))


def _make_mul(*children: Term) -> Term:
    flat: list[Term] = []
    for child in children:
        if _is_zero(child):
            return _zero()
        if isinstance(child, tuple) and child and child[0] == "mul":
            flat.extend(child[1])
        elif not _is_one(child):
            flat.append(child)
    if not flat:
        return _one()
    coefficient = Fraction(1)
    factors: list[Term] = []
    for child in flat:
        numeric = _as_fraction(child)
        if numeric is None:
            factors.append(child)
        else:
            coefficient *= numeric
    if coefficient == 0:
        return _zero()
    # Keep powers local: x*x becomes x^2, but no distributive expansion is
    # performed and no polynomial divisor computation is used.
    counts: Counter[str] = Counter(_term_key(factor) for factor in factors)
    by_key = {_term_key(factor): factor for factor in factors}
    normalized: list[Term] = []
    for key, count in counts.items():
        factor = by_key[key]
        normalized.append(factor if count == 1 else _make_pow(factor, _rat(count)))
    normalized.sort(key=_term_key)
    if coefficient != 1:
        normalized.insert(0, _rat(coefficient))
    if not normalized:
        return _rat(coefficient)
    if len(normalized) == 1:
        return normalized[0]
    return ("mul", tuple(normalized))


def _make_div(numerator: Term, denominator: Term, *, guarded: bool = False) -> Term:
    if _is_zero(denominator):
        raise _Unresolved("exact quotient has an exact zero denominator")
    if _is_zero(numerator):
        return _zero()
    if _is_one(denominator):
        return numerator
    n = _as_fraction(numerator)
    d = _as_fraction(denominator)
    if n is not None and d is not None:
        return _rat(n / d)
    if guarded and numerator == denominator:
        return _one()
    return ("div", numerator, denominator)


def _make_pow(base: Term, exponent: Term) -> Term:
    n = _as_fraction(exponent)
    if n is None or n.denominator != 1:
        raise _Unresolved("only exact integer powers are admitted")
    power = n.numerator
    if power == 0:
        return _one()
    if power == 1:
        return base
    if _is_zero(base):
        if power < 0:
            raise _Unresolved("negative power of exact zero")
        return _zero()
    numeric = _as_fraction(base)
    if numeric is not None:
        if power < 0 and numeric == 0:
            raise _Unresolved("negative power of exact zero")
        return _rat(numeric ** power)
    return ("pow", base, _rat(power))


def _make_func(name: str, argument: Term, metadata: Any = None) -> Term:
    if name == "exp" and _is_zero(argument):
        return _one()
    if name in {"log", "log1p"} and _is_one(argument if name == "log" else _make_add(_one(), argument)):
        return _zero()
    if name in {"sqrt", "positive_sqrt"} and _is_zero(argument):
        return _zero()
    if name == "sin" and _is_zero(argument):
        return _zero()
    if name == "cos" and _is_zero(argument):
        return _one()
    if name == "expm1" and _is_zero(argument):
        return _zero()
    if metadata is None:
        return ("func", name, argument)
    return ("func", name, argument, canonical_digest(metadata))


def _make_compound(name: str, arguments: Sequence[Term], metadata: Any = None) -> Term:
    digest = canonical_digest(metadata or {})
    return ("compound", name, tuple(arguments), digest)


def _term_text(value: Any) -> str:
    return f"STRUCTURAL_FUNCTION_VALUE:{_term_digest(value)}"


def _record(value: Term) -> dict[str, Any]:
    digest = _term_digest(value)
    return {
        "srepr": f"STRUCTURAL_FUNCTION_VALUE:{digest}",
        "text": f"constructor-root:{digest}",
        "is_zero": _is_zero(value),
        "formal_digest": digest,
        "constructor_root_digest": digest,
        "constructor_tag": value[0],
    }


def _unresolved_record(handle: int, error: str) -> dict[str, Any]:
    marker = f"UNRESOLVED(handle={handle}; {error})"
    return {"srepr": marker, "text": marker, "is_zero": False, "formal_digest": canonical_digest(marker)}


@dataclass(frozen=True)
class _Alias:
    side: str
    node: int
    role: str
    namespace: str
    key: str
    form_digest: str
    value: Term
    support: str
    guard: Mapping[str, Any]


class _AliasTable:
    def __init__(self, nodes: Sequence[Mapping[str, Any]]) -> None:
        self.nodes = nodes
        self.by_side: dict[str, dict[int, _Alias]] = {side: {} for side in SIDES}
        self.rows: list[dict[str, Any]] = []

    def add(self, alias: _Alias) -> None:
        _require(alias.side in SIDES, f"unknown alias side {alias.side!r}")
        _require(0 <= alias.node < DERIVATIVE_PREFIX, f"alias {alias.role} enters the derivative suffix")
        _require(isinstance(alias.guard, Mapping), f"alias {alias.role} lacks a guard mapping")
        _require(alias.guard.get("Z_handle") == Z_HANDLE, f"alias {alias.role} lacks the admitted Z guard")
        if "ordinary_Z_order" in alias.guard:
            _require(alias.guard.get("ordinary_Z_order") in ORDERS, f"alias {alias.role} has an invalid order guard")
        _require(isinstance(alias.support, str) and alias.support.startswith("inherited:"), f"alias {alias.role} lacks support provenance")
        if alias.node == 0:
            _require(_exact_zero(self.nodes, 0) and _is_zero(alias.value), f"node 0 alias {alias.role} is not exact zero")
        previous = self.by_side[alias.side].get(alias.node)
        if previous is not None:
            if previous.namespace != alias.namespace:
                # The admitted zero handle is shared by VALUE and native-D
                # source premises.  Keep the namespaces distinct in the
                # witness while allowing this one exact-zero representation
                # to serve both lookup paths; every nonzero cross-namespace
                # collision remains a hard gate.
                _require(previous.value == alias.value == _zero() and previous.form_digest == alias.form_digest, f"VALUE and native-D alias namespaces collide at {alias.side} node {alias.node}")
                self.rows.append({"side": alias.side, "node": alias.node, "role": alias.role, "namespace": alias.namespace, "collision": "same_exact_zero_distinct_namespace"})
                return
            _require(previous.value == alias.value, f"conflicting alias value at {alias.side} node {alias.node}")
            _require(previous.form_digest == alias.form_digest, f"conflicting alias form at {alias.side} node {alias.node}")
            self.rows.append({"side": alias.side, "node": alias.node, "role": alias.role, "collision": "same_exact_form"})
            return
        self.by_side[alias.side][alias.node] = alias
        self.rows.append({
            "side": alias.side, "node": alias.node, "role": alias.role, "namespace": alias.namespace,
            "key": alias.key, "form_digest": alias.form_digest, "support": alias.support,
            "guard": deepcopy(dict(alias.guard)), "zero": _is_zero(alias.value),
        })


# Graph rows carry both constructor semantics and execution/provenance
# receipts.  The latter must not become accidental function arguments merely
# because the patch chart records more receipts than the postrestore chart.
# Every key observed in the admitted reachable grammar is listed here; a new
# key fails closed in ``_constructor_metadata`` instead of being silently
# discarded.  Handle-valued fields are canonicalized as child terms by the
# caller and are therefore omitted from the metadata projection.
_KNOWN_METADATA_KEYS = frozenset({
    "name", "function", "unary", "kind", "operator", "source_node", "source_side",
    "field", "ordinary_Z_order", "ordinary_slow_Z_derivative_order", "Z_order",
    "original_chart", "native_chart", "source_formula", "guarded",
    "Z_independent", "Z_independent_substitution", "exact_native_chain_rule_substitution",
    "source_positive_certificate", "positive_certificate", "denominator_positive",
    "definition", "exact_function_not_its_interval_cap", "source_binding",
    "exact_definition", "original_radius_not_absorbed_into_raw_shape", "original_source",
    "Taylor_coefficient_factorial", "coordinate_independent", "phase_independent_leading_source",
    "defining_quantity_not_a_range_value", "leading_history_normalization_already_applied",
    "quantity", "quantity_path", "recipe", "source_family", "source_projection_binding",
    "exact_source_op_P0_function_not_enclosure_endpoint",
    "original_P0_source_registry_guard", "canonical_parameter_proof", "current_source_receipt",
    "current_source_receipt_sha256", "anchor_definition", "exact_fixed_source_parameter", "fixed_leading_input_N0",
    "parameter_enclosure_not_selected_value", "quantity_is_exact_original_parameter_not_range_endpoint",
    "source_attribute", "source_path", "source_provider", "actual_green_inverse_provider",
    "actual_source_provider", "actual_twenty_term_operator", "admitted_solution_not_a_model_or_enclosure_center",
    "axis_conditions", "component_definition", "radial_order", "radial_recurrence_provider",
    "unique_analytic_Banach_solution", "actual_anchored_gauge_integral", "actual_core_normalized_moment",
    "actual_history_name", "actual_prescribed_angular_equation", "actual_prescribed_axial_equation",
    "actual_six_history_rate", "all_hydro_pressure_swirl_sectors_retained", "anchor_is_exact_root_not_disk_center",
    "angle_endpoints_independent_of_Z", "both_original_micro_controls_retained", "comparison_not_actual_field",
    "endpoint_independent_of_Z", "exact_function_integral", "exact_native_partial_of_same_integrand",
    "exact_original_Rm_bump_weight", "exact_original_comparison_smoothing", "exact_original_positive_beta_normalization",
    "exact_power", "full_original_support", "measure", "multiplicity", "numerical_cap_not_defining_integrand",
    "numerical_value_not_installed", "original_H_root_anchor", "original_first_switch_axial_control",
    "original_first_switch_full_current_radius_direction", "original_incoming_memory_retained",
    "original_second_switch_full_current_radius_direction_and_power", "prescribed_bridge_chart",
    "range_or_midpoint_not_defining_integral", "same_unique_nonlinear_core_and_full_density_product",
    "source_measure_applied_once", "true_core_inlet_complete_history", "true_loop_C5_original_integrand",
    "ordinary_measure", "original_center", "original_radius", "Z_independent_weight",
    "exact_source", "exact_positive_integer", "J", "one_common_N_for_all_charts_and_control_denominators",
    "outer_integer_and_inner_power_not_materialized", "type", "no_additional_Jacobian",
    "exact_value_of_same_unique_original_common_box_solution", "component", "system",
    "branch_receipt_sha256", "defining_map", "exact_Z_independent_row_and_control_units",
    "exact_defect_value_nodes", "exact_source_and_weight_branch_transfer", "exact_weight_nodes", "source_bindings",
    "invAm2", "numerical_enclosures_not_values", "original_common_box_branch_receipt",
    "original_leading_fixed_point_not_Rc_finite_N_control", "original_whole_Z_base_unique_branch_and_Jacobian_admitted",
    "original_beta_exact_flat_derivative", "original_collar_definition_module", "original_collar_extra_flat_condition",
    "original_collar_proof_module", "original_collar_receipt", "original_collar_receipt_sha256",
    "original_sigma_exact_logistic_first_derivative", "original_sigma_logistic_ODE",
    "same_original_flat_extension_through_order5", "same_original_lazy_native_derivative",
    "angle_variable", "branch_inverse_is_original_not_support_cap", "derivative_positive_certificate",
    "flat_inverse", "phase_endpoints", "true_prescribed_source_C5_implicit_inverse", "angle_lower", "angle_upper",
    "lexical_alpha_renamed", "flat_predicate_preserved", "active_body_not_evaluated_on_flat_branch",
    "exact_Z_independent", "original_collar_coordinate", "original_collar_upper",
})

_ARITHMETIC_CONSTRUCTOR_OPS = frozenset({
    "sum", "add", "product", "mul", "negative", "neg", "difference", "sub", "subtract",
    "positive_quotient", "quotient", "integer_power", "power", "function_substitution",
    "exact_real_comparison", "logical_or", "logical_and", "logical_not",
    "substitute_original_inverse_angle", "inverse_substitution", "substitute_original_native_coordinate",
    "substitute_free_query_at_native_coordinate", "native_substitution",
})

_INTEGRAL_SEMANTIC_KEYS = frozenset({
    "actual_anchored_gauge_integral", "actual_core_normalized_moment", "actual_history_name",
    "actual_prescribed_angular_equation", "actual_prescribed_axial_equation", "actual_six_history_rate",
    "all_hydro_pressure_swirl_sectors_retained", "anchor_is_exact_root_not_disk_center",
    "angle_endpoints_independent_of_Z", "both_original_micro_controls_retained", "comparison_not_actual_field",
    "endpoint_independent_of_Z", "exact_function_integral", "exact_native_partial_of_same_integrand",
    "exact_original_Rm_bump_weight", "exact_original_comparison_smoothing", "exact_original_positive_beta_normalization",
    "exact_power", "full_original_support", "measure", "multiplicity", "ordinary_measure",
    "original_H_root_anchor", "original_first_switch_axial_control", "original_first_switch_full_current_radius_direction",
    "original_incoming_memory_retained", "original_second_switch_full_current_radius_direction_and_power",
    "prescribed_bridge_chart", "range_or_midpoint_not_defining_integral", "same_unique_nonlinear_core_and_full_density_product",
    "source_measure_applied_once", "true_core_inlet_complete_history", "true_loop_C5_original_integrand",
    "numerical_cap_not_defining_integrand", "numerical_value_not_installed", "original_center", "original_radius",
    "Z_independent_weight",
})


def _descriptor_without_context(value: Any) -> Any:
    """Remove chart-local coordinates from a nested source recipe only."""

    if isinstance(value, Mapping):
        return {
            str(key): _descriptor_without_context(item)
            for key, item in value.items()
            if key not in {"native_chart", "original_chart", "coordinate", "Z_variable", "angle_variable", "angle_variable_node"}
        }
    if isinstance(value, list):
        return [_descriptor_without_context(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_descriptor_without_context(item) for item in value)
    return value


def _constructor_metadata(row: Mapping[str, Any], operation: str, excluded: Sequence[str] = ()) -> dict[str, Any]:
    """Project row receipts to exact constructor identity metadata.

    Arithmetic guards and chart labels are checked by their dedicated rules;
    they are not function arguments.  Opaque source constructors retain their
    source-node identity unless an exact source descriptor is present.  This
    keeps accidental name-based merges impossible while allowing the two
    charts to share a certified source recipe.
    """

    ignored = set(excluded) | {"operation", "op", "node", "handle", "root"}
    for key in row:
        if key in ignored or key in _HANDLE_FIELDS:
            continue
        _require(key in _KNOWN_METADATA_KEYS, f"node metadata key {key!r} is outside the admitted constructor grammar")
    if operation in _ARITHMETIC_CONSTRUCTOR_OPS:
        selected: Sequence[str] = ()
    elif operation == "definite_integral":
        selected = _INTEGRAL_SEMANTIC_KEYS
    elif operation == "compact_raw_beta":
        selected = ("exact_definition", "original_radius_not_absorbed_into_raw_shape", "original_source")
    elif operation == "exact_original_source_flat_sigma":
        selected = ("definition", "exact_function_not_its_interval_cap", "source_binding")
    elif operation == "current_original_leading_function_recipe":
        selected = ("Taylor_coefficient_factorial", "Z_order", "coordinate_independent", "defining_quantity_not_a_range_value", "exact_source_op_P0_function_not_enclosure_endpoint", "leading_history_normalization_already_applied", "phase_independent_leading_source", "quantity", "quantity_path", "recipe", "source_family", "source_projection_binding")
    elif operation == "current_original_source_parameter":
        selected = ("Z_independent", "anchor_definition", "canonical_parameter_proof", "current_source_receipt", "current_source_receipt_sha256", "definition", "exact_fixed_source_parameter", "fixed_leading_input_N0", "name", "parameter_enclosure_not_selected_value", "quantity_is_exact_original_parameter_not_range_endpoint", "source_attribute", "source_path", "source_provider", "source_family")
    elif operation == "current_original_unique_analytic_core_component":
        selected = ("actual_green_inverse_provider", "actual_source_provider", "actual_twenty_term_operator", "admitted_solution_not_a_model_or_enclosure_center", "axis_conditions", "component_definition", "ordinary_slow_Z_derivative_order", "quantity", "radial_order", "radial_recurrence_provider", "source_family", "source_provider", "unique_analytic_Banach_solution")
    elif operation == "original_lazy_flat_branch":
        selected = ("active_body_not_evaluated_on_flat_branch", "original_beta_exact_flat_derivative", "original_collar_definition_module", "original_collar_extra_flat_condition", "original_collar_proof_module", "original_collar_receipt", "original_collar_receipt_sha256", "original_sigma_exact_logistic_first_derivative", "original_sigma_logistic_ODE", "same_original_flat_extension_through_order5", "same_original_lazy_native_derivative")
    elif operation == "original_monotone_phase_inverse":
        selected = ("branch_inverse_is_original_not_support_cap", "derivative_positive_certificate", "flat_inverse", "phase_endpoints", "true_prescribed_source_C5_implicit_inverse")
    elif operation == "original_Rm_unique_leading_implicit_control":
        selected = ("component", "exact_value_of_same_unique_original_common_box_solution", "source_family")
    elif operation == "original_Rm_unique_leading_implicit_system":
        selected = ("branch_receipt_sha256", "defining_map", "fixed_leading_input_N0", "invAm2", "numerical_enclosures_not_values", "original_common_box_branch_receipt", "original_leading_fixed_point_not_Rc_finite_N_control", "original_whole_Z_base_unique_branch_and_Jacobian_admitted", "source_bindings", "source_family", "source_projection_binding")
    elif operation == "exact_positive_integer_expression":
        selected = ("J", "definition", "exact_positive_integer", "one_common_N_for_all_charts_and_control_denominators", "outer_integer_and_inner_power_not_materialized", "type", "Z_independent")
    else:
        # Unknown compound operations remain source-node-qualified.  They can
        # never merge merely because two rows carry the same display label.
        selected = tuple(key for key in row if key not in ignored and key not in _HANDLE_FIELDS)
    result: dict[str, Any] = {}
    for key in selected:
        if key not in row:
            continue
        value = row[key]
        if key in {"recipe", "source_projection_binding"}:
            value = _descriptor_without_context(value)
        if key == "angle_variable" or key.endswith("_angle_variable"):
            value = "ALPHA_BINDER"
        result[key] = deepcopy(value)
    if operation in {"current_original_leading_function_recipe", "current_original_source_parameter", "current_original_unique_analytic_core_component", "original_Rm_unique_leading_implicit_control", "original_Rm_unique_leading_implicit_system"} and row.get("source_node") is not None:
        result["source_node"] = row["source_node"]
    return result


def _children(nodes: Sequence[Mapping[str, Any]], index: int, row: Mapping[str, Any]) -> list[tuple[str, int]]:
    result: list[tuple[str, int]] = []
    for field in _HANDLE_FIELDS:
        if field not in row:
            continue
        value = row[field]
        if field in {"arguments", "operands"} and _is_sequence(value):
            for offset, child in enumerate(value):
                result.append((f"{field}[{offset}]", _child(nodes, index, child, f"{field}[{offset}]")))
        elif (isinstance(value, int) and not isinstance(value, bool)) or isinstance(value, Mapping):
            result.append((field, _child(nodes, index, value, field)))
    return result


def _validate_graph(nodes: Sequence[Mapping[str, Any]], expected_nodes_digest: str) -> tuple[str, str]:
    _require(_is_sequence(nodes), "nodes must be a sequence")
    _require(len(nodes) == GRAPH_NODE_COUNT, f"source graph must contain {GRAPH_NODE_COUNT} nodes")
    _require(all(isinstance(node, Mapping) for node in nodes), "nodes must contain mappings")
    expected = _digest(expected_nodes_digest, "expected_nodes_digest")
    actual = canonical_digest(list(nodes))
    _require(actual == expected, "nodes digest does not match expected_nodes_digest")
    _require(expected == EXPECTED_NODES_DIGEST, "source graph is not the admitted 205300-node graph")
    prefix = canonical_digest(list(nodes[:DERIVATIVE_PREFIX]))
    _require(prefix == EXPECTED_DERIVATIVE_PREFIX_DIGEST, "derivative source prefix digest is not the admitted prefix")
    _require(_operation(nodes, Z_HANDLE) == "bound_variable" and nodes[Z_HANDLE].get("name") == "Z", "graph Z handle is not the admitted Z variable")
    _require(_exact_zero(nodes, 0), "graph node 0 must be exact zero")
    return actual, prefix


def _validate_generic(packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]], aliases: _AliasTable) -> dict[str, Any]:
    _require(packet.get("schema") == "original_Rm_generic_native_identity_v1", "generic packet has an unknown schema")
    _require(packet.get("passed") is True and packet.get("source_only") is True, "generic packet is not passed source data")
    _require(packet.get("leading_only", packet.get("leadingonly")) is True, "generic packet is not leading-only")
    _require(packet.get("ordinary_orders") == [0, 1], "generic packet ordinary orders changed")
    _require(packet.get("unresolved_gates", []) == [], "generic packet has unresolved gates")
    rows = packet.get("root_rows")
    _require(_is_sequence(rows) and len(rows) == 18, "generic packet must contain 18 root rows")
    seen: set[tuple[str, int]] = set()
    generic: dict[tuple[str, int], dict[str, Any]] = {}
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"generic row {offset} is malformed")
        root = row.get("root", row.get("name"))
        order = row.get("order", row.get("ordinary_Z_order"))
        _require(root in GENERIC_ROOTS and isinstance(order, int) and not isinstance(order, bool) and order in ORDERS, f"generic row {offset} has invalid root/order")
        key = (root, order)
        _require(key not in seen, f"duplicate generic row {root}[{order}]")
        seen.add(key)
        _require(row.get("status") == "proved", f"generic row {root}[{order}] is not proved")
        side_data: dict[str, dict[str, Any]] = {}
        for side in SIDES:
            owner = _side_owner(row, side, f"generic {root}[{order}]")
            node = _handle(owner.get("root_node", owner.get("node")), f"generic {side}.{root}[{order}].node", len(nodes))
            form = _owner_form(owner, f"generic {side}.{root}[{order}]")
            native = _form(owner.get("native_derivative"), f"generic {side}.{root}[{order}].native_derivative")
            side_data[side] = {"node": node, "value": form, "native_derivative": native}
        _require((side_data["postrestore"]["node"], side_data["patch"]["node"]) == EXPECTED_GENERIC_PORTS[key], f"generic source ports changed for {root}[{order}]")
        _require(side_data["postrestore"]["value"] == side_data["patch"]["value"], f"generic VALUE forms differ for {root}[{order}]")
        _require(side_data["postrestore"]["native_derivative"] == side_data["patch"]["native_derivative"], f"generic native-D forms differ for {root}[{order}]")
        value = _zero() if side_data["postrestore"]["value"]["is_zero"] else ("atom", "GENERIC_VALUE", _source_form_digest(side_data["postrestore"]["value"]))
        form_digest = _source_form_digest(side_data["postrestore"]["value"])
        for side in SIDES:
            aliases.add(_Alias(side, side_data[side]["node"], f"generic:{root}[{order}]", "GENERIC_VALUE", f"{root}[{order}]", form_digest, value, "inherited:passed_generic_packet", {"ordinary_Z_order": order, "Z_handle": Z_HANDLE, "separate_ordinary_Z_atom": True}))
        generic[key] = {"root": root, "order": order, "sides": side_data}
    _require(seen == {(root, order) for root in GENERIC_ROOTS for order in ORDERS}, "generic packet does not cover all 18 roots")
    return {"schema": packet.get("schema"), "digest": canonical_digest(packet), "rows": generic}


def _validate_common_phase(packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    common = packet.get("common_phase_source")
    _require(isinstance(common, Mapping) and common.get("schema") == "actual_Rm_common_phase_source_v1", "common phase source has an unknown schema")
    _require(common.get("passed") is True and common.get("source_only") is True, "common phase source is not passed source data")
    for key in ("actual_common_phase_proved", "smooth_periodic_lift_bound", "actual_phase_Z_independence_proved"):
        _require(common.get(key) is True, f"common phase source.{key} is not true")
    _require(common.get("fractional_part_differentiated") is False, "common phase fractional part was differentiated")
    for key in ("numerical_phase_value_computed", "native_loop_derivative_proved", "density_seam_proved", "inverse_conditioning_proved"):
        _require(common.get(key) is False, f"common phase source.{key} must be false")
    _require(common.get("phase_rows") == list(EXPECTED_PHASE_ROWS), "common phase rows changed")
    endpoints = common.get("endpoints")
    _require(_is_sequence(endpoints) and len(endpoints) == 2 and all(isinstance(value, int) and not isinstance(value, bool) for value in endpoints) and list(endpoints) == [1, 1], "common phase endpoints are not exact one")
    _require(common.get("radius_lemma_digest") == EXPECTED_RADIUS_LEMMA_DIGEST, "common phase radius lemma changed")
    _require(common.get("endpoint_radius_reduction") == EXPECTED_ENDPOINT_RADIUS_REDUCTION, "common phase endpoint radius reduction changed")
    phase_nodes = common.get("phase_nodes")
    coordinates = common.get("coordinates")
    if isinstance(phase_nodes, Mapping):
        phase_nodes = [phase_nodes["postrestore"], phase_nodes["patch"]]
    if isinstance(coordinates, Mapping):
        coordinates = [coordinates["postrestore"], coordinates["patch"]]
    _require(phase_nodes == [EXPECTED_PHASE_NODES["postrestore"], EXPECTED_PHASE_NODES["patch"]], "common phase node handles changed")
    _require(coordinates == [EXPECTED_COORDINATES["postrestore"], EXPECTED_COORDINATES["patch"]], "common phase coordinate handles changed")
    _require(all(_operation(nodes, phase_nodes[i]) == "analytic_unary" for i in range(2)), "actual phase nodes are not analytic source constructors")
    _require(all(nodes[phase_nodes[i]].get("name") == "fractional_part" for i in range(2)), "actual phase nodes are not fractional-part phase constructors")
    metadata = {
        "schema": common.get("schema"),
        "phase_rows": deepcopy(list(common["phase_rows"])),
        "endpoints": [1, 1],
        "exact_source_descriptors": deepcopy(common.get("exact_source_descriptors")),
        "endpoint_identity": deepcopy(common.get("endpoint_identity")),
        "unfolded_native_identity": deepcopy(common.get("unfolded_native_identity")),
        "proof_rules": deepcopy(common.get("proof_rules")),
        "radius_lemma_digest": common.get("radius_lemma_digest"),
        "endpoint_radius_reduction": deepcopy(common.get("endpoint_radius_reduction")),
        "source_graph_digest": common.get("source_graph_digest"),
        "same_loop_source_sha256": common.get("same_loop_source_sha256"),
    }
    return {
        "phase_nodes": {"postrestore": phase_nodes[0], "patch": phase_nodes[1]},
        "coordinates": {"postrestore": coordinates[0], "patch": coordinates[1]},
        "endpoint_certificate": {
            "native_endpoints": deepcopy(EXPECTED_NATIVE_ENDPOINTS),
            "coordinates": deepcopy(EXPECTED_COORDINATES),
            "radius_lemma_digest": EXPECTED_RADIUS_LEMMA_DIGEST,
            "endpoint_radius_reduction": deepcopy(EXPECTED_ENDPOINT_RADIUS_REDUCTION),
        },
        "metadata": metadata,
        "digest": canonical_digest(common),
        "source": common,
    }


def _validate_loop(packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    _require(packet.get("schema") == "admitted_Rm_loop_source_ports_v1", "loop admission has an unknown schema")
    _require(packet.get("passed") is True and packet.get("source_only") is True, "loop admission is not passed source data")
    _require(packet.get("mathematical_operator_executed") is False and packet.get("numerical_values_computed") is False, "loop admission executed a mathematical operator")
    _require(packet.get("graph_node_count") == 205200 and packet.get("ordinary_orders") == [0, 1], "loop admission graph/order contract changed")
    rows = packet.get("rows")
    _require(_is_sequence(rows) and len(rows) == 22, "loop admission must contain 22 VALUE rows")
    seen: set[tuple[str, int]] = set()
    compact_rows: list[dict[str, Any]] = []
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"loop row {offset} is malformed")
        field = row.get("field")
        order = row.get("order", row.get("ordinary_Z_order"))
        key = (field, order)
        _require(key in EXPECTED_LOOP_PORTS and key not in seen and isinstance(order, int) and order in ORDERS, f"loop row {offset} has an invalid field/order")
        seen.add(key)
        _require(row.get("status") == "proved", f"loop row {field}[{order}] is not proved")
        left = _side_owner(row, "postrestore", f"loop {field}[{order}]")
        right = _side_owner(row, "patch", f"loop {field}[{order}]")
        left_node = _handle(left.get("node", left.get("root_node")), f"loop postrestore {field}[{order}].node", len(nodes))
        right_node = _handle(right.get("node", right.get("root_node")), f"loop patch {field}[{order}].node", len(nodes))
        _require((left_node, right_node) == EXPECTED_LOOP_PORTS[key], f"loop source ports changed for {field}[{order}]")
        _require(_owner_form(left, f"loop postrestore {field}[{order}]") == _owner_form(right, f"loop patch {field}[{order}]") , f"loop forms differ for {field}[{order}]")
        compact_rows.append({"field": field, "order": order, "postrestore": left_node, "patch": right_node})
    _require(seen == set(EXPECTED_LOOP_PORTS), "loop admission does not cover all 22 rows")
    phase = _validate_common_phase(packet, nodes)
    loop_input = packet.get("loop_input_packet")
    _require(isinstance(loop_input, Mapping) and loop_input.get("source_only") is True and loop_input.get("no_accepted_operator_replayed") is True, "loop input packet is not admitted source data")
    return {"schema": packet.get("schema"), "digest": canonical_digest(packet), "phase": phase, "rows": compact_rows}


def _unwrap_provenance(provenance: Mapping[str, Any], side: str) -> Mapping[str, Any]:
    if isinstance(provenance.get("sides"), Mapping):
        key = "postrestore" if side == "postrestore" else "patch"
        candidate = provenance["sides"].get(key)
        _require(isinstance(candidate, Mapping), f"{side} provenance lacks its side receipt")
        return candidate
    return provenance


def _validate_provenance(provenance: Mapping[str, Any], side: str, nodes: Sequence[Mapping[str, Any]], expected: str, phase: Mapping[str, Any]) -> dict[str, Any]:
    packet = _unwrap_provenance(provenance, side)
    expected_schema = "original_Rm_native_derivative_provenance_v1" if side == "postrestore" else "Rm_control_bound_patch_native_provenance_v1"
    _require(packet.get("schema") == expected_schema, f"{side} provenance has an unknown schema")
    _require(packet.get("side") in ({"postrestore", "postslow"} if side == "postrestore" else {"patch"}), f"{side} provenance side label changed")
    _require(packet.get("passed") is True and packet.get("status") == "proved" and packet.get("source_only") is True, f"{side} provenance is not passed source data")
    for key in ("source_executed", "actual_saved_graph_run", "numeric_values_evaluated", "numeric_field_evaluated"):
        if key not in packet:
            continue
        _require(packet.get(key) is False, f"{side} provenance {key} is not false")
    if "nodes_input_unmutated" in packet:
        _require(packet.get("nodes_input_unmutated") is True, f"{side} provenance mutated the source graph")
    _require(_digest(packet.get("nodes_digest"), f"{side}.nodes_digest") == expected and _digest(packet.get("expected_nodes_digest"), f"{side}.expected_nodes_digest") == expected, f"{side} provenance graph digest does not match")
    _require(packet.get("derivative_source_prefix_count") == DERIVATIVE_PREFIX, f"{side} derivative lookup prefix changed")
    _require(_digest(packet.get("derivative_source_prefix_digest"), f"{side}.derivative_source_prefix_digest") == canonical_digest(list(nodes[:DERIVATIVE_PREFIX]),), f"{side} derivative prefix digest does not match")
    coordinate = _handle(packet.get("coordinate_handle"), f"{side}.coordinate_handle", len(nodes))
    phase_handle = _handle(packet.get("actual_phase_handle"), f"{side}.actual_phase_handle", len(nodes))
    _require(coordinate == phase["coordinates"][side] and phase_handle == phase["phase_nodes"][side], f"{side} provenance coordinate/phase is not the common pair")
    _require(packet.get("physical_J_handle") == EXPECTED_PHYSICAL_J[side] and packet.get("inverse_J_handle") == EXPECTED_INVERSE_J[side], f"{side} J handles changed")
    rows = packet.get("rows")
    _require(_is_sequence(rows) and len(rows) == 4, f"{side} provenance must contain four rows")
    compact: list[dict[str, Any]] = []
    seen: set[tuple[str, int]] = set()
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"{side} provenance row {offset} is malformed")
        field = row.get("field")
        order = row.get("ordinary_order", row.get("order"))
        key = (field, order)
        _require(field in SLOW_FIELDS and isinstance(order, int) and order in ORDERS and key not in seen, f"{side} provenance row {offset} has an invalid field/order")
        seen.add(key)
        row_side = row.get("side", side)
        allowed_row_sides = {side}
        if side == "postrestore":
            allowed_row_sides.add("postslow")
        _require(row_side in allowed_row_sides, f"{side} provenance row {field}[{order}] has a foreign side")
        _require(row.get("status") == "proved" and row.get("physical_target_matches") is True, f"{side} provenance row {field}[{order}] is not proved")
        source = _handle(row.get("source_handle", row.get("source_node")), f"{side}.{field}[{order}].source_handle", len(nodes))
        derivative = _handle(row.get("derivative_handle", row.get("derivative_root")), f"{side}.{field}[{order}].derivative_handle", len(nodes))
        target = _handle(row.get("physical_target_handle", row.get("target_handle")), f"{side}.{field}[{order}].physical_target_handle", len(nodes))
        if row.get("source_node") is not None:
            _require(row.get("source_node") == source, f"{side}.{field}[{order}] source handle mismatch")
        if row.get("derivative_root") is not None:
            _require(row.get("derivative_root") == derivative, f"{side}.{field}[{order}] derivative handle mismatch")
        if row.get("target_handle") is not None:
            _require(row.get("target_handle") == target, f"{side}.{field}[{order}] target handle mismatch")
        if row.get("expected_physical_target_handle") is not None:
            _require(row.get("expected_physical_target_handle") == target, f"{side}.{field}[{order}] expected target handle mismatch")
        _require((source, derivative, target) == EXPECTED_PROVENANCE_PORTS[side][key], f"{side}.{field}[{order}] source/D/target handles changed")
        _require(source < DERIVATIVE_PREFIX and target < DERIVATIVE_PREFIX, f"{side}.{field}[{order}] enters derivative suffix")
        compact.append({"field": field, "ordinary_order": order, "source_handle": source, "derivative_handle": derivative, "physical_target_handle": target})
    _require(seen == set(EXPECTED_PROVENANCE_PORTS[side]), f"{side} provenance rows do not cover A/B ordinary orders")
    unresolved = packet.get("unresolved_gates", [])
    _require(_is_sequence(unresolved) and not unresolved, f"{side} provenance has unresolved gates")
    visited = packet.get("visited_source_D_pairs")
    _require(_is_sequence(visited) and len(visited) == EXPECTED_VISITED_COUNTS[side], f"{side} visited source/D pair count changed")
    exact: dict[tuple[int, tuple[Any, ...]], int] = {}
    for offset, raw in enumerate(visited):
        _require(isinstance(raw, Mapping), f"{side} visited pair {offset} is malformed")
        source = _handle(raw.get("source_handle", raw.get("source_node")), f"{side} visited pair {offset}.source", len(nodes))
        derivative = _handle(raw.get("derivative_handle", raw.get("derivative_root")), f"{side} visited pair {offset}.derivative", len(nodes))
        context_raw = raw.get("derivative_context", raw.get("context"))
        _require(_is_sequence(context_raw) and len(context_raw) == 2, f"{side} visited pair {offset} context is malformed")
        context = tuple(context_raw)
        prior = exact.get((source, context))
        _require(prior is None or prior == derivative, f"{side} conflicting visited source/context derivative")
        exact[(source, context)] = derivative
    seed_handles = packet.get("seed_handles")
    _require(isinstance(seed_handles, Mapping) and seed_handles.get(str(phase_handle)) == 0, f"{side} held phase seed is not exact zero")
    _require(exact.get((phase_handle, ("node", coordinate))) == 0, f"{side} held phase native D zero binding is missing")
    if side == "patch":
        _require(packet.get("accepted_left_provenance_consumed_as_data") is True, "patch provenance did not consume the accepted left receipt as data")
        _require(packet.get("accepted_left_traversal_replayed") is False, "patch provenance replayed the left traversal")
    return {"packet": packet, "rows": sorted(compact, key=lambda item: (item["field"], item["ordinary_order"])), "coordinate": coordinate, "phase": phase_handle, "physical_J": packet["physical_J_handle"], "inverse_J": packet["inverse_J_handle"], "visited": list(visited), "exact_visited": exact, "visited_digest": canonical_digest(list(visited))}


def _validate_control_certificate(packet: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    cert = packet.get("control_certificate")
    _require(isinstance(cert, Mapping) and cert.get("passed") is True, "patch control certificate is not passed")
    if packet.get("source_family") is not None and cert.get("source_family") is not None:
        _require(cert.get("source_family") == packet.get("source_family"), "patch control source family changed")
    _require(cert.get("system") == EXPECTED_CONTROL_SYSTEM and cert.get("control_value_nodes") == list(EXPECTED_CONTROL_NODES), "patch control handles changed")
    _require(cert.get("system_input_free_variables") == ["Z"] and cert.get("control_values_not_set_to_zero") is True and cert.get("Z_dependence_retained") is True, "patch control support contract changed")
    _require(cert.get("all_higher_control_rows_checked_by_lexical_dependency") is True and cert.get("same_unique_original_common_box_solution") is True and cert.get("implicit_recurrence_replayed") is False, "patch control lexical certificate changed")
    rows = cert.get("control_support_rows")
    _require(_is_sequence(rows) and len(rows) == 35, "patch control certificate must contain 35 support rows")
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping) and row.get("free_variables") == ["Z"] and row.get("native_coordinate_independent_by_lexical_support") is True, f"patch control support row {offset} is not exactly Z-supported")
    for node in EXPECTED_CONTROL_NODES:
        _require(_operation(nodes, node) == "original_Rm_unique_leading_implicit_control", f"patch control node {node} has an unexpected operation")
    return {"system": EXPECTED_CONTROL_SYSTEM, "control_nodes": list(EXPECTED_CONTROL_NODES), "support": ["Z"], "support_row_count": len(rows), "support_rows": deepcopy(list(rows)), "support_row_digest": canonical_digest(list(rows)), "digest": canonical_digest(cert)}


def _validate_plan(plan: Mapping[str, Any], left: Mapping[str, Any], right: Mapping[str, Any], phase: Mapping[str, Any]) -> dict[str, Any]:
    _require(isinstance(plan, Mapping) and plan.get("schema") == "actual_Rm_paired_native_value_targets_v1", "paired_target_plan has an unknown schema")
    _require(plan.get("native_coordinates") == EXPECTED_COORDINATES and plan.get("held_phase") == EXPECTED_PHASE_NODES, "paired target coordinate/phase handles changed")
    native_endpoints = plan.get("native_endpoints")
    _require(isinstance(native_endpoints, Mapping) and all(isinstance(native_endpoints.get(side), int) and not isinstance(native_endpoints.get(side), bool) for side in SIDES) and native_endpoints == EXPECTED_NATIVE_ENDPOINTS, "paired target native endpoints are not the certified exact endpoint")
    endpoint_certificate = phase.get("endpoint_certificate")
    _require(isinstance(endpoint_certificate, Mapping), "common phase lacks its endpoint certificate")
    _require(endpoint_certificate.get("native_endpoints") == EXPECTED_NATIVE_ENDPOINTS, "common phase endpoint certificate changed")
    _require(endpoint_certificate.get("coordinates") == EXPECTED_COORDINATES, "common phase endpoint coordinates changed")
    _require(endpoint_certificate.get("radius_lemma_digest") == EXPECTED_RADIUS_LEMMA_DIGEST, "common phase endpoint radius lemma is not admitted")
    _require(endpoint_certificate.get("endpoint_radius_reduction") == EXPECTED_ENDPOINT_RADIUS_REDUCTION, "common phase endpoint radius reduction is not admitted")
    rows = plan.get("targets")
    expected_keys = [("A", 0), ("A", 1), ("B", 0), ("B", 1)]
    _require(_is_sequence(rows) and len(rows) == 4, "paired target plan must contain four rows")
    seen: list[tuple[str, int]] = []
    compact: list[dict[str, Any]] = []
    left_by_key = {(row["field"], row["ordinary_order"]): row for row in left["rows"]}
    right_by_key = {(row["field"], row["ordinary_order"]): row for row in right["rows"]}
    for offset, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"paired target row {offset} is malformed")
        key = (row.get("field"), row.get("ordinary_order"))
        _require(key in expected_keys and key not in seen, f"paired target row {offset} has an invalid field/order")
        seen.append(key)
        values: dict[str, Any] = {"field": key[0], "ordinary_order": key[1]}
        for side, info, source_key in (("postrestore", left, "postrestore"), ("patch", right, "patch")):
            side_row = row.get(side)
            _require(isinstance(side_row, Mapping), f"paired target row {key} lacks {side}")
            expected = (left_by_key if side == "postrestore" else right_by_key)[key]
            actual = tuple(_handle(side_row.get(name), f"paired {side}.{key}.{name}", GRAPH_NODE_COUNT) for name in ("source_handle", "derivative_handle", "physical_target_handle"))
            required = (expected["source_handle"], expected["derivative_handle"], expected["physical_target_handle"])
            _require(actual == required, f"paired target row {key} does not match {side} provenance")
            _require(info["exact_visited"].get((actual[0], ("node", info["coordinate"]))) == actual[1], f"paired target row {key} lacks exact {side} source/D context")
            values[side] = {"source_handle": actual[0], "derivative_handle": actual[1], "physical_target_handle": actual[2]}
        compact.append(values)
    _require(seen == expected_keys, "paired target plan rows are not ordered A0,A1,B0,B1")
    return {
        "schema": plan.get("schema"),
        "native_coordinates": deepcopy(EXPECTED_COORDINATES),
        "native_endpoints": deepcopy(EXPECTED_NATIVE_ENDPOINTS),
        "held_phase": deepcopy(EXPECTED_PHASE_NODES),
        "targets": compact,
        "digest": canonical_digest(plan),
    }


def _validate_loop_provenance_pairing(loop: Mapping[str, Any], left: Mapping[str, Any], right: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the four physical source ports to separately admitted loop rows."""

    loop_rows = {(row["field"], row["order"]): row for row in loop["rows"]}
    matches: list[dict[str, Any]] = []
    for field, loop_field in (("A", "A"), ("B", "B_over_Pstar")):
        for order in ORDERS:
            key = (field, order)
            loop_row = loop_rows.get((loop_field, order))
            _require(isinstance(loop_row, Mapping), f"loop admission lacks {loop_field}[{order}]")
            for side, info in (("postrestore", left), ("patch", right)):
                source = next((row["source_handle"] for row in info["rows"] if (row["field"], row["ordinary_order"]) == key), None)
                _require(source == loop_row[side], f"{side} {field}[{order}] source does not match admitted loop {loop_field}[{order}]")
                matches.append({"side": side, "field": field, "ordinary_order": order, "loop_field": loop_field, "source_handle": source})
    return {"passed": True, "matches": matches, "digest": canonical_digest(matches)}


def _d_aliases(generic: Mapping[str, Any], info: Mapping[str, Any], side: str, aliases: _AliasTable) -> tuple[dict[int, _Alias], list[dict[str, Any]], dict[int, tuple[str, dict[str, Any]]]]:
    d_aliases: dict[int, _Alias] = {}
    roles: dict[int, tuple[str, dict[str, Any]]] = {}
    bindings: list[dict[str, Any]] = []
    for (root, order), item in generic["rows"].items():
        source = item["sides"][side]["node"]
        derivative = info["exact_visited"].get((source, ("node", info["coordinate"])))
        if derivative is None:
            continue
        form = item["sides"][side]["native_derivative"]
        form_digest = _source_form_digest(form)
        prior = roles.get(derivative)
        if prior is not None:
            _require(prior[0] == form_digest, f"conflicting D alias forms at {side} node {derivative}")
        value = _zero() if form["is_zero"] else ("atom", "NATIVE_D", form_digest)
        alias = _Alias(side, derivative, f"D:{root}[{order}]", "NATIVE_D", f"{root}[{order}]", form_digest, value, "inherited:visited_source_D_provenance", {"ordinary_Z_order": order, "Z_handle": Z_HANDLE, "derivative_context": ["node", info["coordinate"]], "native_derivative_form": True})
        aliases.add(alias)
        d_aliases[derivative] = alias
        roles[derivative] = (form_digest, form)
        bindings.append({"root": root, "ordinary_order": order, "source_handle": source, "derivative_handle": derivative, "derivative_context": ["node", info["coordinate"]]})
    return d_aliases, bindings, roles


class _SupportWalker:
    """Lexical source-support walk for inherited VALUE/D aliases."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]], *, side: str, coordinate: int, control_nodes: Sequence[int] = ()) -> None:
        self.nodes = nodes
        self.side = side
        self.coordinate = coordinate
        self.control_nodes = set(control_nodes)
        self.cache: dict[int, frozenset[str]] = {}
        self.active: set[int] = set()
        self.visited: set[int] = set()
        self.operation_counts: Counter[str] = Counter()

    def free(self, index: int) -> frozenset[str]:
        _require(index not in self.active, f"cyclic lexical support dependency at node {index}")
        if index in self.cache:
            return self.cache[index]
        row = _row(self.nodes, index, "support node")
        op = _operation(self.nodes, index)
        self.active.add(index)
        self.visited.add(index)
        self.operation_counts[op] += 1
        try:
            def child(name: str) -> int:
                return _child(self.nodes, index, row[name], name)
            if op == "bound_variable":
                result = frozenset((str(row.get("name", f"node:{index}")),))
            elif op in _SOURCE_CONSTANT_OPS or op in {"exact_rational", "mathematical_pi", "mathematical_positive_infinity"}:
                result = frozenset()
            elif op in {"original_source_parameter", "current_original_source_parameter"}:
                _require(row.get("Z_independent") is True, f"source parameter {index} lacks Z-independent support")
                result = frozenset()
            elif op == "original_Rm_unique_leading_implicit_control":
                _require(self.side == "patch" and index in self.control_nodes, f"unsupported leading control {index}")
                result = frozenset(("Z",))
            elif op in {"sum", "add", "product", "mul", "logical_or", "logical_and"}:
                result = frozenset().union(*(self.free(c) for _, c in _children(self.nodes, index, row)))
            elif op in {"negative", "neg", "analytic_unary", "exact_original_source_flat_sigma", "compact_raw_beta", *_ANALYTIC, *_OPAQUE_ANALYTIC}:
                result = self.free(child("argument"))
            elif op in {"positive_quotient", "quotient", "difference", "sub", "subtract", "integer_power", "power"}:
                if op == "positive_quotient":
                    _require(row.get("source_positive_certificate", row.get("positive_certificate", row.get("denominator_positive"))) not in (None, False, ""), f"support quotient {index} lacks positivity guard")
                result = frozenset().union(*(self.free(c) for _, c in _children(self.nodes, index, row)))
            elif op == "function_substitution":
                _require(any(row.get(key) is True for key in _FUNCTION_GUARDS), f"support substitution {index} lacks a guard")
                expression = child("expression")
                variable = child("variable")
                value = child("value")
                binder = str(_row(self.nodes, variable).get("name", f"node:{variable}"))
                values = set(self.free(expression)); values.discard(binder); values.update(self.free(value)); result = frozenset(values)
            elif op in {"substitute_original_native_coordinate", "substitute_free_query_at_native_coordinate", "native_substitution"}:
                body = child("body") if "body" in row else child("expression")
                query = child("query_coordinate") if "query_coordinate" in row else child("free_query")
                variable = row.get("native_variable", row.get("variable", row.get("coordinate")))
                binder = str(_row(self.nodes, variable).get("name", f"node:{variable}")) if isinstance(variable, int) else str(variable)
                values = set(self.free(body)); values.discard(binder); values.update(self.free(query)); result = frozenset(values)
            elif op in {"substitute_original_inverse_angle", "inverse_substitution"}:
                body = child("body") if "body" in row else child("expression")
                inverse = child("inverse_angle") if "inverse_angle" in row else child("inverse")
                values = set(self.free(body)); values.discard(str(row.get("angle_variable"))); values.update(self.free(inverse)); result = frozenset(values)
            elif op == "original_monotone_phase_inverse":
                body = child("phase_function"); target = child("target_phase")
                values = set(self.free(body)); values.discard(str(row.get("angle_variable"))); values.update(self.free(target))
                for field in ("Delta", "delta", "eta", "eta21", "angle_lower", "angle_upper"):
                    if isinstance(row.get(field), int) and not isinstance(row.get(field), bool): values.update(self.free(child(field)))
                result = frozenset(values)
            elif op == "definite_integral":
                integrand = child("integrand"); lower = child("lower"); upper = child("upper")
                variable = row.get("variable")
                binder = str(_row(self.nodes, variable).get("name", f"node:{variable}")) if isinstance(variable, int) else str(variable)
                values = set(self.free(integrand)); values.discard(binder); values.update(self.free(lower)); values.update(self.free(upper)); result = frozenset(values)
            elif op == "original_lazy_flat_branch":
                result = frozenset().union(*(self.free(c) for _, c in _children(self.nodes, index, row)))
            elif op in _COMPOUND_OPS:
                result = frozenset().union(*(self.free(c) for _, c in _children(self.nodes, index, row)))
            else:
                raise ValueSeamError(f"unknown lexical support operation at node {index}: {op!r}")
            self.cache[index] = result
            return result
        finally:
            self.active.remove(index)

    def certificate(self, candidates: Mapping[int, Sequence[str]], allowed: set[str]) -> dict[str, Any]:
        rows: list[dict[str, Any]] = []
        for node, roles in sorted(candidates.items()):
            support = sorted(self.free(node))
            _require(set(support) <= allowed, f"support node {node} exposes an unadmitted variable")
            rows.append({"node": node, "roles": sorted(roles), "support": support})
        return {"passed": True, "allowed_support": sorted(allowed), "candidate_count": len(rows), "visited_node_count": len(self.visited), "visited_node_digest": canonical_digest(sorted(self.visited)), "operation_counts": dict(sorted(self.operation_counts.items())), "support_digest": canonical_digest(rows)}


class _Normalizer:
    """Bounded exact reducer over source constructor nodes."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]], *, side: str, aliases: Mapping[int, _Alias], native_coordinate: int, endpoint_endpoints: Mapping[str, int], control_nodes: Sequence[int] = ()) -> None:
        self.nodes = nodes
        self.side = side
        self.aliases = dict(aliases)
        self.native_coordinate = native_coordinate
        self.endpoint_endpoints = dict(endpoint_endpoints)
        _require(self.endpoint_endpoints.get(side) == 1, f"{side} lacks the certified native endpoint")
        self.native_name = str(nodes[native_coordinate].get("name", f"node:{native_coordinate}"))
        self.control_nodes = set(control_nodes)
        self.cache: dict[tuple[int, str], Term] = {}
        self.active: set[tuple[int, str]] = set()
        self.used_nodes: set[int] = set()
        self.operations: Counter[str] = Counter()
        self.events: list[dict[str, Any]] = []
        self.binder_depth = 0
        self.bound_name_nodes: dict[tuple[int, str], frozenset[int]] = {}

    @staticmethod
    def _name_key(name: str) -> tuple[str, str]:
        return ("__structural_name_binding__", name)

    def _name_binding(self, env: Mapping[Any, Any], name: str) -> tuple[frozenset[int], Term] | None:
        raw = env.get(self._name_key(name))
        if raw is None:
            return None
        _require(isinstance(raw, tuple) and len(raw) == 2, f"lexical name binding {name!r} is malformed")
        allowed, replacement = raw
        _require(isinstance(allowed, frozenset) and all(isinstance(item, int) and not isinstance(item, bool) for item in allowed), f"lexical name binding {name!r} has malformed node identity")
        _require(isinstance(replacement, tuple), f"lexical name binding {name!r} has a non-constructor replacement")
        return allowed, replacement

    def _collect_bound_name_nodes(self, root: int, name: str) -> frozenset[int]:
        cache_key = (root, name)
        cached = self.bound_name_nodes.get(cache_key)
        if cached is not None:
            return cached
        pending = [root]
        seen: set[int] = set()
        matches: set[int] = set()
        while pending:
            node = pending.pop()
            if node in seen:
                continue
            seen.add(node)
            row = _row(self.nodes, node)
            if _operation(self.nodes, node) == "bound_variable" and row.get("name", row.get("variable")) == name:
                matches.add(node)
            pending.extend(child for _, child in _children(self.nodes, node, row))
            _require(len(seen) <= GRAPH_NODE_COUNT, f"lexical binder walk exceeded the source graph bound at {root}")
        result = frozenset(matches)
        _require(len(result) <= 1, f"lexical name {name!r} resolves to multiple bound-variable node identities")
        self.bound_name_nodes[cache_key] = result
        return result

    def _bind_name(self, env: dict[Any, Any], name: str, replacement: Term, allowed_nodes: frozenset[int], label: str) -> None:
        _require(isinstance(name, str) and name, f"{label} has a malformed lexical name")
        _require(len(allowed_nodes) <= 1, f"{label} has colliding bound-variable node identities")
        env[self._name_key(name)] = (allowed_nodes, replacement)

    def _env_signature(self, env: Mapping[Any, Term]) -> str:
        items: list[Any] = []
        for key, value in env.items():
            if isinstance(key, tuple) and len(key) == 2 and key[0] == "__structural_name_binding__":
                _require(isinstance(value, tuple) and len(value) == 2 and isinstance(value[0], frozenset), "lexical environment contains a malformed name binding")
                items.append((repr(key), tuple(sorted(value[0])), _term_key(value[1])))
            else:
                items.append((repr(key), _term_key(value)))
        items.sort(key=lambda item: item[0])
        return canonical_digest(items)

    def _aliases_allowed(self, env: Mapping[Any, Term]) -> bool:
        if Z_HANDLE in env or self._name_binding(env, "Z") is not None:
            return False
        if self.native_coordinate in env:
            return _is_one(env[self.native_coordinate])
        native_binding = self._name_binding(env, self.native_name)
        if native_binding is not None:
            return _is_one(native_binding[1])
        return True

    def _bound(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        name = row.get("name", row.get("variable"))
        if index in env:
            return env[index]
        if isinstance(name, str):
            binding = self._name_binding(env, name)
            if binding is not None:
                allowed_nodes, replacement = binding
                _require(index in allowed_nodes, f"bound variable node {index} reuses lexical name {name!r} outside its certified scope")
                return replacement
        if index == self.native_coordinate:
            _require(self.endpoint_endpoints.get(self.side) == 1, f"native coordinate {index} lacks an explicit exact endpoint certificate")
            self.events.append({"node": index, "kind": "endpoint_native_coordinate", "value": 1})
            return _one()
        return ("atom", "BOUND_VARIABLE", str(name) if isinstance(name, str) else f"node:{index}")

    def _rewrite_child(self, index: int, value: Any, label: str, env: Mapping[Any, Term]) -> Term:
        return self.rewrite(_child(self.nodes, index, value, label), env)

    def _substitution(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        _require(any(row.get(key) is True for key in _FUNCTION_GUARDS), f"node {index} substitution lacks its guard")
        expression = _child(self.nodes, index, _field(row, "expression", "body"), "expression")
        variable = _child(self.nodes, index, _field(row, "variable", "native_variable"), "variable")
        value = _child(self.nodes, index, _field(row, "value", "endpoint", "query"), "value")
        replacement = self.rewrite(value, env)
        scoped = dict(env); scoped[variable] = replacement
        name = self.nodes[variable].get("name")
        if isinstance(name, str):
            self._bind_name(scoped, name, replacement, frozenset({variable}), f"node {index} substitution binder")
        self.events.append({"node": index, "kind": "guarded_function_substitution", "variable": variable, "value": value})
        return self.rewrite(expression, scoped)

    def _native_substitution(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        body = _child(self.nodes, index, _field(row, "body", "expression"), "body")
        query = _child(self.nodes, index, _field(row, "query_coordinate", "free_query", "value", "coordinate", "native_coordinate"), "query_coordinate")
        query_term = self.rewrite(query, env)
        variable = row.get("native_variable", row.get("variable", row.get("coordinate")))
        scoped = dict(env)
        excluded = ["body", "expression", "query_coordinate", "free_query", "value", "coordinate", "native_coordinate", "native_variable", "variable"]
        if isinstance(variable, int) and not isinstance(variable, bool):
            binder = _child(self.nodes, index, variable, "native_variable")
            scoped[binder] = query_term
            name = self.nodes[binder].get("name")
            if isinstance(name, str):
                self._bind_name(scoped, name, query_term, frozenset({binder}), f"node {index} native substitution binder")
        elif isinstance(variable, str):
            self._bind_name(scoped, variable, query_term, self._collect_bound_name_nodes(body, variable), f"node {index} native substitution binder")
        else:
            raise _Unresolved(f"node {index} native substitution lacks a variable binder")
        body_term = self.rewrite(body, scoped)
        metadata = _constructor_metadata(row, "substitute_original_native_coordinate", excluded)
        metadata["variable"] = "NATIVE_COORDINATE"
        return _make_compound(_operation(self.nodes, index), (body_term, query_term), metadata)

    def _inverse_substitution(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        angle = row.get("angle_variable")
        _require(isinstance(angle, str) and angle, f"node {index} inverse substitution lacks angle_variable")
        body = _child(self.nodes, index, _field(row, "body", "expression"), "body")
        inverse = _child(self.nodes, index, _field(row, "inverse_angle", "inverse", "value"), "inverse_angle")
        binder = ("atom", "ALPHA_BINDER", f"angle:{self.binder_depth}")
        scoped = dict(env)
        angle_node_raw = row.get("angle_variable_node")
        if angle_node_raw is not None:
            angle_node = _child(self.nodes, index, angle_node_raw, "angle_variable_node")
            scoped[angle_node] = binder
            allowed_nodes = frozenset({angle_node})
        else:
            allowed_nodes = self._collect_bound_name_nodes(body, angle)
        self._bind_name(scoped, angle, binder, allowed_nodes, f"node {index} inverse angle binder")
        self.binder_depth += 1
        try:
            body_term = self.rewrite(body, scoped)
        finally:
            self.binder_depth -= 1
        inverse_term = self.rewrite(inverse, env)
        metadata = _constructor_metadata(row, "substitute_original_inverse_angle", ("body", "expression", "inverse_angle", "inverse", "value", "angle_variable", "angle_variable_node"))
        metadata["angle_variable"] = "ALPHA_BINDER"
        return ("inverse", body_term, inverse_term, canonical_digest(metadata))

    def _phase_inverse(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        for key in ("angle_variable", "derivative_positive_certificate", "flat_inverse", "phase_endpoints"):
            _require(row.get(key) not in (None, "", False), f"node {index} phase inverse lacks {key}")
        _require(row.get("branch_inverse_is_original_not_support_cap") is True, f"node {index} phase inverse branch is not admitted")
        angle = row.get("angle_variable")
        _require(isinstance(angle, str), f"node {index} phase inverse angle binder is malformed")
        phase_function = _child(self.nodes, index, row.get("phase_function"), "phase_function")
        target_phase = _child(self.nodes, index, row.get("target_phase"), "target_phase")
        binder = ("atom", "ALPHA_BINDER", f"angle:{self.binder_depth}")
        scoped = dict(env)
        angle_node_raw = row.get("angle_variable_node")
        if angle_node_raw is not None:
            angle_node = _child(self.nodes, index, angle_node_raw, "angle_variable_node")
            scoped[angle_node] = binder
            allowed_nodes = frozenset({angle_node})
        else:
            allowed_nodes = self._collect_bound_name_nodes(phase_function, angle)
        self._bind_name(scoped, angle, binder, allowed_nodes, f"node {index} phase angle binder")
        self.binder_depth += 1
        try:
            function_term = self.rewrite(phase_function, scoped)
        finally:
            self.binder_depth -= 1
        target_term = self.rewrite(target_phase, env)
        arguments: list[Term] = [function_term, target_term]
        excluded = ["phase_function", "target_phase", "angle_variable", "angle_variable_node"]
        for key in ("Delta", "delta", "eta", "eta21", "angle_lower", "angle_upper"):
            if isinstance(row.get(key), int) and not isinstance(row.get(key), bool):
                arguments.append(self._rewrite_child(index, row[key], key, env)); excluded.append(key)
        metadata = _constructor_metadata(row, "original_monotone_phase_inverse", excluded)
        metadata["angle_variable"] = "ALPHA_BINDER"
        metadata["branch_inverse_is_original_not_support_cap"] = True
        return _make_compound("original_monotone_phase_inverse", arguments, metadata)

    def _comparison(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Any:
        left = self._rewrite_child(index, row.get("left"), "left", env)
        right = self._rewrite_child(index, row.get("right"), "right", env)
        operator = row.get("operator")
        _require(operator in {"<", "<=", "==", ">=", ">", "!="}, f"node {index} has an unsupported comparison")
        lf = _as_fraction(left); rf = _as_fraction(right)
        if lf is not None and rf is not None:
            return {"<": lf < rf, "<=": lf <= rf, "==": lf == rf, ">=": lf >= rf, ">": lf > rf, "!=": lf != rf}[operator]
        if left == right:
            return {"<": False, "<=": True, "==": True, ">=": True, ">": False, "!=": False}[operator]
        return ("comparison", operator, left, right)

    def _lazy_branch(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        _require(row.get("active_body_not_evaluated_on_flat_branch") is True, f"node {index} flat branch does not preserve laziness")
        predicate = _child(self.nodes, index, row.get("flat_predicate"), "flat_predicate")
        predicate_term = self.rewrite(predicate, env)
        _require("flat_value" in row, f"node {index} flat branch lacks an explicit flat_value handle")
        flat_value = row["flat_value"]
        _require(isinstance(flat_value, int) and not isinstance(flat_value, bool), f"node {index} flat_value must be an integer graph handle")
        flat_handle = _child(self.nodes, index, flat_value, "flat_value")
        # Delta/eta/collar handles are semantic constructor children.  Keeping
        # their raw graph numbers in branch metadata would make the identical
        # postrestore and patch branch functions look different solely because
        # the two charts have different graph offsets.
        if predicate_term is True:
            self.events.append({"node": index, "kind": "lazy_flat_branch", "selected": "flat"})
            return self.rewrite(flat_handle, env)
        active = _child(self.nodes, index, _field(row, "active_body", "body"), "active_body")
        if predicate_term is False:
            self.events.append({"node": index, "kind": "lazy_flat_branch", "selected": "active"})
            return self.rewrite(active, env)
        flat_term = self.rewrite(flat_handle, env)
        active_term = self.rewrite(active, env)
        guard_terms: list[Term] = []
        guard_fields: list[str] = []
        for field in ("Delta", "delta", "eta", "eta21", "original_collar_coordinate", "original_collar_upper"):
            value = row.get(field)
            if isinstance(value, int) and not isinstance(value, bool):
                guard_terms.append(self._rewrite_child(index, value, field, env))
                guard_fields.append(field)
        metadata = _constructor_metadata(row, "original_lazy_flat_branch", ("flat_predicate", "active_body", "body", "flat_value", *guard_fields))
        metadata["flat_predicate_preserved"] = True
        metadata["active_body_not_evaluated_on_flat_branch"] = True
        return ("branch", predicate_term, active_term, flat_term, tuple(guard_terms), canonical_digest(metadata))

    def _integral(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Term]) -> Term:
        _require(row.get("exact_function_integral") is True, f"node {index} integral lacks its exact-function guard")
        _require(row.get("numerical_value_not_installed") is True, f"node {index} integral carries a numerical value")
        integrand = _child(self.nodes, index, row.get("integrand"), "integrand")
        lower = _child(self.nodes, index, row.get("lower"), "lower")
        upper = _child(self.nodes, index, row.get("upper"), "upper")
        variable = row.get("variable")
        _require(isinstance(variable, (str, int)), f"node {index} integral lacks a lexical variable")
        binder = ("atom", "ALPHA_BINDER", f"integral:{self.binder_depth}")
        scoped = dict(env)
        if isinstance(variable, int):
            variable_node = _child(self.nodes, index, variable, "variable")
            scoped[variable_node] = binder
            name = self.nodes[variable_node].get("name")
            if isinstance(name, str):
                self._bind_name(scoped, name, binder, frozenset({variable_node}), f"node {index} integral binder")
        else:
            self._bind_name(scoped, variable, binder, self._collect_bound_name_nodes(integrand, variable), f"node {index} integral binder")
        self.binder_depth += 1
        try:
            integrand_term = self.rewrite(integrand, scoped)
        finally:
            self.binder_depth -= 1
        lower_term = self.rewrite(lower, env); upper_term = self.rewrite(upper, env)
        metadata = _constructor_metadata(row, "definite_integral", ("integrand", "lower", "upper", "variable"))
        metadata["variable"] = "ALPHA_ANGLE_BINDER"
        metadata["lexical_alpha_renamed"] = True
        return ("integral", integrand_term, lower_term, upper_term, canonical_digest(metadata))

    def _generic_compound(self, index: int, row: Mapping[str, Any], operation: str, env: Mapping[Any, Term]) -> Term:
        values: list[Term] = []
        excluded: list[str] = []
        for field, child in _children(self.nodes, index, row):
            values.append(self.rewrite(child, env)); excluded.append(field.split("[")[0])
        metadata = _constructor_metadata(row, operation, excluded)
        if operation in {"current_original_source_parameter", "original_Rm_unique_leading_implicit_control", "original_Rm_unique_leading_implicit_system"}:
            # These are source constructors, not name-based aliases.  Their
            # full exact descriptors (including source receipts) identify the
            # function; a bare label never identifies equality.
            return _make_compound(operation, values, metadata)
        if operation not in _COMPOUND_OPS:
            metadata["source_node"] = index
            metadata["source_side"] = self.side
        return _make_compound(operation, values, metadata)

    def rewrite(self, index: int, env: Mapping[Any, Term] | None = None) -> Term:
        environment = dict(env or {})
        if index in environment:
            return environment[index]
        env_signature = self._env_signature(environment)
        cache_key = (index, env_signature)
        if not environment and cache_key in self.cache:
            return self.cache[cache_key]
        alias = self.aliases.get(index)
        if alias is not None:
            if not self._aliases_allowed(environment):
                raise _Unresolved(f"endpoint alias node {index} reached under an altered native/Z substitution")
            self.used_nodes.add(index)
            _require(len(self.used_nodes) <= MAX_CONSTRUCTOR_NODES, "constructor traversal exceeded its bounded node budget")
            return alias.value
        _require(0 <= index < len(self.nodes), f"source handle {index} is outside the source graph")
        _require(index < DERIVATIVE_PREFIX, f"VALUE traversal entered derivative suffix at {index}")
        if cache_key in self.active:
            raise _Unresolved(f"cyclic VALUE dependency at node {index}")
        row = _row(self.nodes, index)
        operation = _operation(self.nodes, index)
        self.active.add(cache_key); self.used_nodes.add(index); self.operations[operation] += 1
        _require(len(self.used_nodes) <= MAX_CONSTRUCTOR_NODES, "constructor traversal exceeded its bounded node budget")
        try:
            if operation == "exact_rational":
                numerator = row.get("numerator", row.get("value")); denominator = row.get("denominator", 1)
                _require(isinstance(numerator, int) and not isinstance(numerator, bool) and isinstance(denominator, int) and not isinstance(denominator, bool), f"source rational {index} is malformed")
                result = _rat(numerator, denominator)
            elif operation == "exact_positive_integer_expression":
                _require(row.get("exact_positive_integer") is True and row.get("Z_independent") is True, f"source integer {index} lacks its exact source contract")
                result = ("atom", "SOURCE_CONSTANT", canonical_digest(_constructor_metadata(row, operation)))
            elif operation in {"original_source_parameter", "current_original_source_parameter"}:
                _require(row.get("Z_independent") is True, f"source parameter {index} lacks Z-independent support")
                result = self._generic_compound(index, row, operation, environment)
            elif operation == "bound_variable":
                result = self._bound(index, row, environment)
            elif operation in {"sum", "add"}:
                values = [self.rewrite(child, environment) for _, child in _children(self.nodes, index, row)]
                _require(values, f"source sum {index} has no arguments")
                result = _make_add(*values)
            elif operation in {"product", "mul"}:
                values = [self.rewrite(child, environment) for _, child in _children(self.nodes, index, row)]
                _require(values, f"source product {index} has no arguments")
                result = _make_mul(*values)
            elif operation in {"negative", "neg"}:
                result = _make_mul(_rat(-1), self._rewrite_child(index, _field(row, "argument", "operand"), "argument", environment))
            elif operation in {"difference", "sub", "subtract"}:
                values = [self.rewrite(child, environment) for _, child in _children(self.nodes, index, row)]
                _require(len(values) == 2, f"source difference {index} needs two arguments")
                result = _make_add(values[0], _make_mul(_rat(-1), values[1]))
            elif operation in {"positive_quotient", "quotient"}:
                numerator = self._rewrite_child(index, row.get("numerator"), "numerator", environment)
                denominator = self._rewrite_child(index, row.get("denominator"), "denominator", environment)
                guarded = operation == "positive_quotient"
                if guarded:
                    _require(row.get("source_positive_certificate", row.get("positive_certificate", row.get("denominator_positive"))) not in (None, False, ""), f"source quotient {index} lacks positivity guard")
                result = _make_div(numerator, denominator, guarded=guarded)
            elif operation in {"integer_power", "power"}:
                base = self._rewrite_child(index, row.get("base", row.get("argument", row.get("operand"))), "base", environment)
                _require("exponent" in row or "power" in row, f"power {index} lacks an explicit exponent handle")
                exponent_raw = row.get("exponent", row.get("power"))
                exponent_handle = _child(self.nodes, index, exponent_raw, "exponent")
                exponent_row = _row(self.nodes, exponent_handle, "power exponent")
                _require(_operation(self.nodes, exponent_handle) == "exact_rational", f"power {index} has a non-exact exponent")
                numerator = exponent_row.get("numerator", exponent_row.get("value"))
                denominator = exponent_row.get("denominator", 1)
                _require(isinstance(numerator, int) and not isinstance(numerator, bool) and isinstance(denominator, int) and not isinstance(denominator, bool), f"power {index} exponent rational is malformed")
                exponent = _rat(numerator, denominator)
                result = _make_pow(base, exponent)
            elif operation == "function_substitution":
                result = self._substitution(index, row, environment)
            elif operation == "analytic_unary" or operation in _ANALYTIC or operation in _OPAQUE_ANALYTIC:
                name = row.get("name", row.get("function", row.get("unary", row.get("kind"))))
                _require(isinstance(name, str), f"node {index} analytic operation lacks a function name")
                argument = self._rewrite_child(index, _field(row, "argument", "operand", "value"), "argument", environment)
                if name not in _ANALYTIC and name not in _OPAQUE_ANALYTIC:
                    raise _Unresolved(f"unsupported analytic function {name!r} at node {index}")
                result = _make_func(name, argument, None if name in _ANALYTIC else _constructor_metadata(row, "analytic_unary", ("argument", "operand", "value")))
            elif operation == "definite_integral":
                result = self._integral(index, row, environment)
            elif operation == "original_monotone_phase_inverse":
                result = self._phase_inverse(index, row, environment)
            elif operation == "original_lazy_flat_branch":
                result = self._lazy_branch(index, row, environment)
            elif operation == "exact_real_comparison":
                result = self._comparison(index, row, environment)
            elif operation in {"logical_or", "logical_and"}:
                values = [self.rewrite(child, environment) for _, child in _children(self.nodes, index, row)]
                if operation == "logical_or" and any(value is True for value in values): result = True
                elif operation == "logical_and" and any(value is False for value in values): result = False
                elif operation == "logical_or" and all(value is False for value in values): result = False
                elif operation == "logical_and" and all(value is True for value in values): result = True
                else: result = _make_compound(operation, tuple(value for value in values if isinstance(value, tuple)), _constructor_metadata(row, operation, ("arguments", "operands")))
            elif operation == "logical_not":
                value = self._rewrite_child(index, row.get("argument"), "argument", environment)
                result = (not value) if isinstance(value, bool) else _make_compound(operation, (value,), _constructor_metadata(row, operation, ("argument",)))
            elif operation in {"substitute_original_inverse_angle", "inverse_substitution"}:
                result = self._inverse_substitution(index, row, environment)
            elif operation in {"substitute_original_native_coordinate", "substitute_free_query_at_native_coordinate", "native_substitution"}:
                result = self._native_substitution(index, row, environment)
            elif operation == "exact_original_source_flat_sigma":
                argument = self._rewrite_child(index, row.get("argument"), "argument", environment)
                result = _make_func(operation, argument, _constructor_metadata(row, operation, ("argument",)))
            elif operation == "compact_raw_beta":
                _require(row.get("exact_definition", row.get("definition")) in {"exp(-1/(1-x^2)) for |x|<1; exact0 elsewhere", "exp(-1/(1-t^2)) on abs(t)<1, zero otherwise", "exp(-1/(1-t^2)) for |t|<1; exact0 elsewhere"}, f"compact beta {index} changed its exact definition")
                _require(row.get("original_radius_not_absorbed_into_raw_shape", True) is True, f"compact beta {index} changed its radius contract")
                argument = self._rewrite_child(index, row.get("argument"), "argument", environment)
                result = _make_func("compact_raw_beta", argument, _constructor_metadata(row, "compact_raw_beta", ("argument",)))
            elif operation in _SOURCE_CONSTANT_OPS:
                if operation == "exact_original_source_heat_phi":
                    _require(row.get("Z_independent", row.get("exact_source", True)) is not False, f"source constant {index} lost exact source contract")
                result = ("atom", "SOURCE_CONSTANT", canonical_digest({"operation": operation, "metadata": _constructor_metadata(row, operation)}))
            elif operation == "original_Rm_unique_leading_implicit_control":
                _require(index in self.control_nodes and row.get("system") == EXPECTED_CONTROL_SYSTEM and row.get("Z_variable") == Z_HANDLE and row.get("exact_value_of_same_unique_original_common_box_solution") is True, f"control node {index} lacks its exact source contract")
                result = _make_compound(operation, (), _constructor_metadata(row, operation, ("system", "Z_variable")))
            elif operation in _COMPOUND_OPS:
                result = self._generic_compound(index, row, operation, environment)
            else:
                raise _Unresolved(f"unsupported source operation at node {index}: {operation!r}")
            if isinstance(result, tuple):
                self.cache[cache_key] = result
            return result
        finally:
            self.active.remove(cache_key) if cache_key in self.active else None

    def closure(self) -> dict[str, Any]:
        used = sorted(self.used_nodes)
        return {"used_node_count": len(used), "used_node_digest": canonical_digest(used), "used_nodes": used, "operation_counts": dict(sorted(self.operations.items())), "endpoint_events": deepcopy(self.events), "alias_lookup_is_endpoint_guarded": True}


def _support_certificates(nodes: Sequence[Mapping[str, Any]], generic: Mapping[str, Any], d_roles: Mapping[str, Mapping[int, tuple[str, dict[str, Any]]]], infos: Mapping[str, Mapping[str, Any]], control: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for side in SIDES:
        candidates: dict[int, list[str]] = {}
        for (root, order), item in generic["rows"].items():
            candidates.setdefault(item["sides"][side]["node"], []).append(f"VALUE:{root}[{order}]")
        for node in d_roles[side]:
            candidates.setdefault(node, []).append(f"D:{node}")
        coordinate = infos[side]["coordinate"]
        native_name = str(nodes[coordinate].get("name", f"node:{coordinate}"))
        walker = _SupportWalker(nodes, side=side, coordinate=coordinate, control_nodes=control.get("control_nodes", ()) if side == "patch" else ())
        result[side] = walker.certificate(candidates, {"Z", native_name})
        result[side].update({"side": side, "coordinate": coordinate, "native_coordinate_name": native_name, "control_nodes": list(control.get("control_nodes", ())) if side == "patch" else []})
    return result


def derive(
    nodes: Sequence[Mapping[str, Any]],
    generic_packet: Mapping[str, Any],
    loop_ports_admission: Mapping[str, Any],
    left_provenance: Mapping[str, Any],
    right_provenance: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    paired_target_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Prove exact source-constructor congruence for A/B at ordinary Z 0/1.

    ``paired_target_plan`` binds each physical target to its admitted source
    and native-D context.  It is checked against the two provenance receipts
    and does not assert the equality that this function must prove.
    """

    actual, prefix = _validate_graph(nodes, expected_nodes_digest)
    _require(isinstance(generic_packet, Mapping) and isinstance(loop_ports_admission, Mapping), "generic and loop packets must be mappings")
    _require(isinstance(left_provenance, Mapping) and isinstance(right_provenance, Mapping), "both provenance inputs must be mappings")
    aliases = _AliasTable(nodes)
    generic = _validate_generic(generic_packet, nodes, aliases)
    loop = _validate_loop(loop_ports_admission, nodes)
    phase = loop["phase"]
    left = _validate_provenance(left_provenance, "postrestore", nodes, actual, phase)
    right = _validate_provenance(right_provenance, "patch", nodes, actual, phase)
    loop_pairing = _validate_loop_provenance_pairing(loop, left, right)
    control = _validate_control_certificate(right["packet"], nodes)
    plan = _validate_plan(paired_target_plan, left, right, phase)
    d_roles: dict[str, dict[int, tuple[str, dict[str, Any]]]] = {}
    d_bindings: dict[str, list[dict[str, Any]]] = {}
    for side, info in (("postrestore", left), ("patch", right)):
        _, bindings, roles = _d_aliases(generic, info, side, aliases)
        d_bindings[side] = bindings
        d_roles[side] = roles
    # The phase aliases are one common source constructor, with both actual
    # chart phase handles retained in its source metadata.
    phase_term = _make_compound("common_phase", (), phase["metadata"])
    for side in SIDES:
        _require(phase["phase_nodes"][side] not in aliases.by_side[side], f"{side} common phase collides with an admitted alias")
        _require(phase["coordinates"][side] not in aliases.by_side[side], f"{side} native coordinate collides with an admitted alias")
        aliases.add(_Alias(side, phase["phase_nodes"][side], "common_phase", "COMMON_PHASE", "actual_common_phase", phase["digest"], phase_term, "inherited:passed_common_phase_source", {"Z_handle": Z_HANDLE, "endpoint_only": True, "actual_phase_handle": phase["phase_nodes"][side]}))
    supports = _support_certificates(nodes, generic, d_roles, {"postrestore": left, "patch": right}, control)
    targets = {(row["field"], row["ordinary_order"]): row for row in plan["targets"]}
    pair_rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    closures: dict[str, Any] = {}
    normalizers = {side: _Normalizer(nodes, side=side, aliases=aliases.by_side[side], native_coordinate=info["coordinate"], endpoint_endpoints=plan["native_endpoints"], control_nodes=control["control_nodes"] if side == "patch" else ()) for side, info in (("postrestore", left), ("patch", right))}
    for field, order in (row for row in [("A", 0), ("A", 1), ("B", 0), ("B", 1)]):
        row_plan = targets[(field, order)]
        left_handle = row_plan["postrestore"]["physical_target_handle"]
        right_handle = row_plan["patch"]["physical_target_handle"]
        left_error: str | None = None; right_error: str | None = None
        try:
            left_value = normalizers["postrestore"].rewrite(left_handle)
        except (ValueSeamError, TypeError, KeyError, IndexError) as exc:
            left_value = None; left_error = str(exc)
        try:
            right_value = normalizers["patch"].rewrite(right_handle)
        except (ValueSeamError, TypeError, KeyError, IndexError) as exc:
            right_value = None; right_error = str(exc)
        closures[f"postrestore:{field}[{order}]"] = normalizers["postrestore"].closure()
        closures[f"patch:{field}[{order}]"] = normalizers["patch"].closure()
        if left_error is not None or right_error is not None:
            equal = False; status = "unresolved"
            unresolved.append({"gate": "constructor_normalization_unresolved", "field": field, "ordinary_order": order, "postrestore_target_handle": left_handle, "patch_target_handle": right_handle, "postrestore_error": left_error, "patch_error": right_error})
            left_record = _unresolved_record(left_handle, left_error or "missing")
            right_record = _unresolved_record(right_handle, right_error or "missing")
        else:
            _require(isinstance(left_value, tuple) and isinstance(right_value, tuple), f"{field}[{order}] normalized to a non-term value")
            equal = left_value == right_value
            # A nonidentical constructor root is only a failure to admit
            # congruence.  It is not a counterexample and must never be
            # reported as a functional or numerical inequality.
            status = "proved" if equal else "unresolved"
            left_record = _record(left_value); right_record = _record(right_value)
            if not equal:
                unresolved.append({"gate": "constructor_roots_not_identical", "field": field, "ordinary_order": order, "postrestore_target_handle": left_handle, "patch_target_handle": right_handle, "postrestore": left_record, "patch": right_record, "meaning": "exact constructor congruence did not admit equality; this is not an inequality claim"})
        pair_rows.append({"field": field, "ordinary_order": order, "postrestore_target_handle": left_handle, "patch_target_handle": right_handle, "postrestore": left_record, "patch": right_record, "equal": equal, "status": status})
    passed = len(pair_rows) == 4 and not unresolved and all(row["equal"] is True and row["status"] == "proved" for row in pair_rows)
    return {
        "schema": SCHEMA,
        "proof_method": PROOF_METHOD,
        "passed": bool(passed),
        "status": "proved" if passed else "unresolved",
        "source_only": True,
        "sourceonly": True,
        "numeric_values_evaluated": False,
        "no_numerical_values": True,
        "new_native_derivative_generated": False,
        "no_derivatives_generated": True,
        "old_operators_replayed": False,
        "expected_nodes_digest": actual,
        "nodes_digest": actual,
        "graph_node_count": len(nodes),
        "derivative_lookup_prefix": {"count": DERIVATIVE_PREFIX, "digest": prefix},
        "ordinary_orders": [0, 1],
        "canonical_value_pairs": pair_rows,
        "value_pairs": deepcopy(pair_rows),
        "paired_target_plan": deepcopy(plan),
        "alias_roles": deepcopy(aliases.rows),
        "alias_count": len(aliases.rows),
        "visited_bindings": d_bindings,
        "alias_support_certificates": supports,
        "loop_provenance_pairing": loop_pairing,
        "common_phase": {"phase_nodes": deepcopy(phase["phase_nodes"]), "coordinates": deepcopy(phase["coordinates"]), "endpoints": [1, 1], "endpoint_certificate": deepcopy(phase["endpoint_certificate"]), "proof_digest": phase["digest"], "semantic_metadata": deepcopy(phase["metadata"]), "constructor_alias": _record(phase_term)},
        "control_certificate": control,
        "closures": closures,
        "unresolved_gates": unresolved,
        "remaining_gates": deepcopy(unresolved),
        "kernel_api": {
            "entry_point": "derive",
            "required_keywords": ["expected_nodes_digest", "paired_target_plan"],
            "paired_plan_schema": "actual_Rm_paired_native_value_targets_v1",
            "source_to_D_is_consumed_as_data": True,
        },
        "constructor_grammar": {
            "term_tags": sorted(_TERM_TAGS),
            "supported_source_operations": sorted(
                set(_COMPOUND_OPS) | set(_SOURCE_CONSTANT_OPS) | set(_ANALYTIC) |
                set(_OPAQUE_ANALYTIC) | set(_ARITHMETIC_CONSTRUCTOR_OPS) |
                {"bound_variable", "exact_rational", "original_source_parameter", "current_original_source_parameter", "definite_integral", "original_monotone_phase_inverse", "original_lazy_flat_branch"}
            ),
            "local_exact_rules": [
                "flatten_and_sort_same_add_or_mul",
                "collect_identical_add_terms",
                "collect_identical_mul_factors_as_integer_powers",
                "guarded_same_factor_division",
                "exact_rational_arithmetic",
                "exact_zero_and_one_unary_rules",
                "exact_rational_integer_power_only",
                "lexical_alpha_renaming",
                "lazy_branch_preservation",
                "source_guarded_alias_substitution",
            ],
            "unsupported_operations": "fail_closed_as_unresolved_gate",
        },
        "semantic_metadata_rules": {
            "known_key_count": len(_KNOWN_METADATA_KEYS),
            "known_key_digest": canonical_digest(sorted(_KNOWN_METADATA_KEYS)),
            "arithmetic_receipts_are_not_function_arguments": True,
            "source_descriptors_are_retained_for_opaque_constructors": True,
            "chart_context_is_removed_only_from_source_recipe_bindings": True,
            "native_endpoint_one_is_admitted_only_by_paired_plan_and_common_phase_certificate": True,
            "bound_name_reuse_with_multiple_node_identities_fails_closed": True,
            "unknown_top_level_metadata_fails_closed": True,
        },
        "static_checks": {
            "expected_graph_digest_checked": True,
            "derivative_prefix_checked": True,
            "alias_namespaces_separated": True,
            "zero_handle_cross_namespace_exception": "same exact zero form only",
            "endpoint_aliases_guarded_against_Z_or_native_substitution": True,
            "native_endpoint_certificate_required": True,
            "common_phase_radius_endpoint_lemma_checked": True,
            "lexical_binders_require_scoped_node_identity_or_unambiguous_name": True,
            "lazy_flat_value_requires_explicit_graph_handle": True,
            "power_exponents_require_exact_rational_graph_handles": True,
            "constructor_node_budget": MAX_CONSTRUCTOR_NODES,
            "numeric_values_evaluated": False,
            "derivatives_generated": False,
            "symbolic_algebra_or_polynomial_expansion": False,
        },
        "input_data_digests": {"generic_packet": canonical_digest(generic_packet), "loop_ports_admission": canonical_digest(loop_ports_admission), "left_provenance": canonical_digest(left_provenance), "right_provenance": canonical_digest(right_provenance), "paired_target_plan": plan["digest"]},
        "scope": {"four_physical_slow_targets_only": True, "ordinary_orders": [0, 1], "source_to_D_rules_reused_from_admission": True, "constructor_dag_only": True, "derivative_generation": False, "numerical_evaluation": False, "global_C1_claimed": False, "Navier_Stokes_claimed": False, "pressure_join_claimed": False, "closure_claimed": False},
    }


__all__ = ["DERIVATIVE_PREFIX", "EXPECTED_DERIVATIVE_PREFIX_DIGEST", "GRAPH_NODE_COUNT", "MAX_CONSTRUCTOR_NODES", "PROOF_METHOD", "SCHEMA", "ValueSeamError", "canonical_digest", "derive"]
