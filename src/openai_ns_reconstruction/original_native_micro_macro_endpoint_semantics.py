"""Source-only canonicalization of the actual micro/macro endpoint ports.

This module is deliberately a leaf.  It accepts the authenticated source
packets assembled by the runner and canonicalizes only exact source DAG rows;
it never imports a legacy operator, evaluates a numerical model, or executes
one of the archived producers.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from typing import Any, Mapping


SCHEMA = "native_micro_macro_endpoint_semantics_v1"
RECORD_NAMES = [
    "lei_ren_part1_paper_compliant_current_original_native_first_micro_B_jump_source.json",
    "lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json",
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz",
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz",
]
RECORD_HASHES = {
    "lei_ren_part1_paper_compliant_current_original_native_first_micro_B_jump_source.json":
        "588c4bed6a73504a5b6af2be4fb8341ac6f3729db568678fed5693f36db07db9",
    "lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json":
        "aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28",
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz":
        "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz":
        "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz":
        "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92",
}
STAGE_FLAGS = [
    "actual_native_micro_macro_U_C1_identity_proved",
    "actual_native_micro_macro_m_Z2_identity_proved",
    "actual_native_micro_macro_Q_C1_identity_proved",
    "actual_native_micro_macro_R_identity_proved",
]
FALSE_FLAGS = [
    "all_native_B_jumps_zero_proved",
    "native_old_density_identity_proved",
    "terminal_difference_zero_proved",
    "global_five_moment_closure",
    "numeric_values_computed",
    "old_error_terms_removed",
    "old_operators_replayed",
    "time_order_recursion_claimed",
    "full_goal_complete",
]


class UnsupportedSourceOperation(ValueError):
    """Raised when a source graph contains a semantic operation this leaf lacks."""


def _plain_digest(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _jsonable(value: Any) -> Any:
    """Convert canonical tuples and exact rationals to JSON without losing type."""

    if isinstance(value, Fraction):
        return ["fraction", value.numerator, value.denominator]
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, Mapping):
        return {
            str(key): _jsonable(value[key])
            for key in sorted(value, key=lambda item: str(item))
        }
    return value


def _canonical_digest(value: Any) -> str:
    return _plain_digest(_jsonable(value))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _as_fraction(value: Any) -> Fraction | None:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return Fraction(value, 1)
    return None


def _freeze_metadata(value: Any) -> Any:
    """Freeze non-reference metadata into a deterministic, JSON-safe value."""

    if isinstance(value, Mapping):
        return tuple(
            (str(key), _freeze_metadata(value[key]))
            for key in sorted(value, key=lambda item: str(item))
        )
    if isinstance(value, list):
        return tuple(_freeze_metadata(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze_metadata(item) for item in value)
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError("non-finite metadata is not an exact source value")
    return value


def _source_node_digest(nodes: list[Mapping[str, Any]]) -> str:
    # The source graph is already decoded and authenticated by the runner.  A
    # deterministic graph digest is retained for the result without importing
    # any graph producer.
    return _plain_digest(nodes)


def _mapping_with_key(value: Any, key: str, depth: int = 0) -> Mapping[str, Any] | None:
    if depth > 8:
        return None
    if isinstance(value, Mapping):
        if key in value and isinstance(value[key], Mapping):
            return value[key]
        for child in value.values():
            found = _mapping_with_key(child, key, depth + 1)
            if found is not None:
                return found
    elif isinstance(value, list):
        for child in value[:32]:
            found = _mapping_with_key(child, key, depth + 1)
            if found is not None:
                return found
    return None


def _graph_nodes(payload: Any) -> list[Mapping[str, Any]] | None:
    """Find an exact graph without treating arbitrary metadata as a graph."""

    if not isinstance(payload, Mapping):
        return None
    direct = payload.get("exact_graph_nodes")
    if isinstance(direct, list) and all(isinstance(row, Mapping) for row in direct[:8]):
        return direct
    for key in ("definition_graph", "component_result"):
        child = payload.get(key)
        if isinstance(child, Mapping):
            found = _graph_nodes(child)
            if found is not None:
                return found
    candidate = payload.get("nodes")
    if isinstance(candidate, list) and candidate and all(
        isinstance(row, Mapping) for row in candidate[:8]
    ):
        return candidate
    return None


def _find_truth(value: Any, needle: str, depth: int = 0) -> bool:
    """Look for an explicit true source receipt, never infer one from a label."""

    if depth > 10:
        return False
    lowered = needle.lower()
    if isinstance(value, Mapping):
        for key, child in value.items():
            key_text = str(key).lower()
            if lowered in key_text and child is True:
                return True
            if _find_truth(child, needle, depth + 1):
                return True
    elif isinstance(value, list):
        for child in value[:128]:
            if _find_truth(child, needle, depth + 1):
                return True
    return False


def _parse_record(record: Mapping[str, Any]) -> Any:
    text = record.get("source_text")
    _require(isinstance(text, str), "source record lacks source_text")
    decoded = hashlib.sha256(text.encode("utf-8")).hexdigest()
    _require(decoded == record.get("decoded_sha256"), "decoded source hash mismatch")
    _require(record.get("raw_authenticated_before_decode") is True,
             "raw source was not authenticated before decode")
    stripped = text.lstrip()
    if stripped.startswith("{") or stripped.startswith("["):
        return json.loads(text)
    # Optional source-Python context is parsed only.  It is never compiled,
    # imported, or executed by this leaf.
    return {"source_ast": ast.parse(text, filename=str(record.get("source_path", "source.py")))}


def _bind_scope(scope: tuple[tuple[Any, Any], ...], variable: Any, value_id: int,
                row: Mapping[str, Any], store: "CanonicalStore") -> tuple[tuple[Any, Any], ...]:
    keys: list[Any] = []
    if _is_int(variable):
        keys.append(("node", variable))
        if 0 <= variable < len(store.raw_nodes):
            candidate = store.raw_nodes[variable]
            if isinstance(candidate, Mapping) and isinstance(candidate.get("name"), str):
                keys.append(("name", candidate["name"]))
    elif isinstance(variable, str):
        keys.append(("name", variable))
    if isinstance(row.get("variable_name"), str):
        keys.append(("name", row["variable_name"]))
    if not keys:
        raise UnsupportedSourceOperation("unbound substitution variable")
    result = scope
    for key in keys:
        result = result + ((key, value_id),)
    return result


class CanonicalStore:
    """Interned canonical DAG state shared by one authenticated source graph."""

    def __init__(self, raw_nodes: list[Mapping[str, Any]], source_name: str):
        self.raw_nodes = raw_nodes
        self.source_name = source_name
        self.index: dict[tuple[Any, ...], int] = {}
        self.keys: list[tuple[Any, ...]] = []
        self.digests: list[str] = []
        self.digest_to_id: dict[str, int] = {}
        self.memo: dict[tuple[int, tuple[Any, ...]], int] = {}
        self.active: set[tuple[int, tuple[Any, ...]]] = set()
        self.depth = 0

    def intern(self, key: tuple[Any, ...]) -> int:
        existing = self.index.get(key)
        if existing is not None:
            return existing
        node_id = len(self.keys)
        self.index[key] = node_id
        self.keys.append(key)
        digest = _canonical_digest(key)
        self.digests.append(digest)
        self.digest_to_id.setdefault(digest, node_id)
        return node_id

    def key(self, node_id: int) -> tuple[Any, ...]:
        return self.keys[node_id]

    def digest(self, node_id: int) -> str:
        return self.digests[node_id]

    def child_ref(self, node_id: int) -> str:
        return self.digest(node_id)

    def rational(self, node_id: int) -> Fraction | None:
        key = self.key(node_id)
        if key[:1] != ("rational",):
            return None
        return Fraction(int(key[1]), int(key[2]))

    def _scope_signature(self, scope: tuple[tuple[Any, Any], ...]) -> tuple[Any, ...]:
        return tuple((key, self._value_signature(value)) for key, value in scope)

    def _value_signature(self, value: Any) -> Any:
        if isinstance(value, int) and 0 <= value < len(self.keys):
            return ("canonical", self.digest(value))
        return _freeze_metadata(value)

    def lookup(self, scope: tuple[tuple[Any, Any], ...], variable: Any) -> int | None:
        keys: list[Any] = []
        if _is_int(variable):
            keys.append(("node", variable))
            if 0 <= variable < len(self.raw_nodes):
                row = self.raw_nodes[variable]
                if isinstance(row, Mapping) and isinstance(row.get("name"), str):
                    keys.append(("name", row["name"]))
        elif isinstance(variable, str):
            keys.append(("name", variable))
        for key in reversed(keys):
            for bound_key, value in reversed(scope):
                if bound_key == key:
                    return value if isinstance(value, int) else self.intern(value)
        return None

    def _ref_tuple(self, node_id: int) -> tuple[str, str]:
        return ("ref", self.digest(node_id))

    def _metadata(self, row: Mapping[str, Any], omitted: set[str]) -> tuple[Any, ...]:
        return tuple(
            (str(key), _freeze_metadata(value))
            for key, value in sorted(row.items(), key=lambda item: str(item[0]))
            if key not in omitted
        )

    def _scope_bound(self, scope: tuple[tuple[Any, Any], ...], variable: Any) -> int | None:
        return self.lookup(scope, variable)

    def _factor_form(self, node_id: int) -> tuple[Fraction, dict[str, int]]:
        key = self.key(node_id)
        kind = key[0]
        if kind == "rational":
            return Fraction(int(key[1]), int(key[2])), {}
        if kind == "product":
            coefficient = Fraction(int(key[1]), int(key[2]))
            factors = {str(base): int(power) for base, power in key[3]}
            return coefficient, factors
        return Fraction(1, 1), {self.digest(node_id): 1}

    def _from_factors(self, coefficient: Fraction, factors: Mapping[str, int]) -> int:
        cleaned = {str(base): int(power) for base, power in factors.items() if int(power) != 0}
        if coefficient == 0:
            return self.intern(("rational", 0, 1))
        if not cleaned:
            return self.intern(("rational", coefficient.numerator, coefficient.denominator))
        if coefficient == 1 and len(cleaned) == 1:
            base, power = next(iter(cleaned.items()))
            if power == 1 and base in self.digest_to_id:
                return self.digest_to_id[base]
        factor_rows = tuple(sorted((base, power) for base, power in cleaned.items()))
        return self.intern(("product", coefficient.numerator, coefficient.denominator, factor_rows))

    def _sum(self, arguments: list[int]) -> int:
        coefficient = Fraction(0, 1)
        terms: defaultdict[str, Fraction] = defaultdict(Fraction)
        for argument in arguments:
            if self.key(argument)[0] == "sum":
                for numerator, denominator, base in self.key(argument)[1]:
                    value = Fraction(int(numerator), int(denominator))
                    if str(base) == "__constant__":
                        coefficient += value
                    else:
                        terms[str(base)] += value
                continue
            factor, factors = self._factor_form(argument)
            if not factors:
                coefficient += factor
                continue
            base_id = self._from_factors(Fraction(1, 1), factors)
            terms[self.digest(base_id)] += factor
        rows: list[tuple[int, int, str]] = []
        for base, value in terms.items():
            if value:
                rows.append((value.numerator, value.denominator, base))
        rows.sort(key=lambda row: row[2])
        if coefficient:
            rows.append((coefficient.numerator, coefficient.denominator, "__constant__"))
        if not rows:
            return self.intern(("rational", 0, 1))
        if len(rows) == 1 and rows[0][2] == "__constant__":
            return self.intern(("rational", rows[0][0], rows[0][1]))
        if len(rows) == 1:
            numerator, denominator, base = rows[0]
            if base in self.digest_to_id:
                return self._product([self.intern(("rational", numerator, denominator)),
                                      self.digest_to_id[base]])
        return self.intern(("sum", tuple(rows)))

    def _product(self, arguments: list[int]) -> int:
        coefficient = Fraction(1, 1)
        factors: defaultdict[str, int] = defaultdict(int)
        for argument in arguments:
            factor, pieces = self._factor_form(argument)
            coefficient *= factor
            for base, power in pieces.items():
                factors[base] += power
        return self._from_factors(coefficient, factors)

    def _positive_quotient(self, numerator: int, denominator: int,
                           row: Mapping[str, Any]) -> int:
        certificate = row.get("source_positive_certificate")
        if not isinstance(certificate, str) or not certificate.strip():
            raise UnsupportedSourceOperation("positive quotient lacks source certificate")
        numerator_value, numerator_factors = self._factor_form(numerator)
        denominator_value, denominator_factors = self._factor_form(denominator)
        if denominator_value == 0 and not denominator_factors:
            raise UnsupportedSourceOperation("positive quotient has exact zero denominator")
        factors: defaultdict[str, int] = defaultdict(int)
        for base, power in numerator_factors.items():
            factors[base] += power
        for base, power in denominator_factors.items():
            factors[base] -= power
        coefficient = numerator_value / denominator_value
        return self._from_factors(coefficient, factors)

    def _canonical_leaf(self, node_id: int, row: Mapping[str, Any],
                        operation: str, references: tuple[str, ...]) -> int:
        reference_fields = {
            "current_original_axis_pressure_function_derivative": ("original_function", "Z_variable"),
            "current_original_leading_function_recipe": ("coordinate", "Z_variable"),
            "current_original_unique_analytic_core_component": ("rho", "Z_variable"),
        }
        omitted = {"operation", *reference_fields.get(operation, ())}
        return self.intern(("source_leaf", operation, references, self._metadata(row, omitted)))

    def _comparison(self, row: Mapping[str, Any], scope: tuple[tuple[Any, Any], ...]) -> int:
        left = self.walk(row["left"], scope)
        right = self.walk(row["right"], scope)
        left_value, right_value = self.rational(left), self.rational(right)
        operator = row.get("operator")
        if left_value is not None and right_value is not None:
            outcome = {
                "==": left_value == right_value,
                "!=": left_value != right_value,
                "<": left_value < right_value,
                "<=": left_value <= right_value,
                ">": left_value > right_value,
                ">=": left_value >= right_value,
            }.get(operator)
            if outcome is None:
                raise UnsupportedSourceOperation("unknown exact comparison operator")
            return self.intern(("rational", 1 if outcome else 0, 1))
        if operator not in {"==", "!=", "<", "<=", ">", ">="}:
            raise UnsupportedSourceOperation("unknown exact comparison operator")
        return self.intern(("comparison", str(operator), self.child_ref(left), self.child_ref(right)))

    def _logical_or(self, arguments: list[int]) -> int:
        remaining: dict[str, int] = {}
        for argument in arguments:
            value = self.rational(argument)
            if value == 1:
                return self.intern(("rational", 1, 1))
            if value == 0:
                continue
            remaining[self.digest(argument)] = argument
        if not remaining:
            return self.intern(("rational", 0, 1))
        return self.intern(("logical_or", tuple(sorted(remaining))))

    def _analytic(self, row: Mapping[str, Any], argument: int) -> int:
        name = str(row.get("name", ""))
        value = self.rational(argument)
        if value is not None:
            if name in {"exp", "exponential"} and value == 0:
                return self.intern(("rational", 1, 1))
            if name in {"expm1", "sin", "sinh", "tan", "tanh"} and value == 0:
                return self.intern(("rational", 0, 1))
            if name in {"cos", "cosh"} and value == 0:
                return self.intern(("rational", 1, 1))
            if name == "log" and value == 1:
                return self.intern(("rational", 0, 1))
        return self.intern(("analytic", name, self.child_ref(argument),
                            self._metadata(row, {"operation", "name", "argument"})))

    def _flat_sigma(self, row: Mapping[str, Any], argument: int) -> int:
        value = self.rational(argument)
        if value is not None:
            if value <= 0:
                return self.intern(("rational", 0, 1))
            if value >= 1:
                return self.intern(("rational", 1, 1))
        return self.intern(("flat_sigma", self.child_ref(argument),
                            self._metadata(row, {"operation", "argument"})))

    def _lazy_branch(self, row: Mapping[str, Any], scope: tuple[tuple[Any, Any], ...]) -> int:
        predicate = self.walk(row["flat_predicate"], scope)
        truth = self.rational(predicate)
        if truth == 1:
            return self.walk(row["flat_value"], scope)
        if truth == 0:
            return self.walk(row["active_body"], scope)
        active = self.walk(row["active_body"], scope)
        flat = self.walk(row["flat_value"], scope)
        delta = self.walk(row["Delta"], scope) if _is_int(row.get("Delta")) else None
        eta = self.walk(row["eta"], scope) if _is_int(row.get("eta")) else None
        collar_coordinate = (
            self.walk(row["original_collar_coordinate"], scope)
            if _is_int(row.get("original_collar_coordinate")) else None
        )
        collar_upper = (
            self.walk(row["original_collar_upper"], scope)
            if _is_int(row.get("original_collar_upper")) else None
        )
        return self.intern((
            "lazy_branch", self.child_ref(predicate), self.child_ref(flat),
            self.child_ref(active),
            self.child_ref(delta) if delta is not None else None,
            self.child_ref(eta) if eta is not None else None,
            self.child_ref(collar_coordinate) if collar_coordinate is not None else None,
            self.child_ref(collar_upper) if collar_upper is not None else None,
            self._metadata(row, {
                "operation", "flat_predicate", "flat_value", "active_body", "Delta", "eta",
                "original_collar_coordinate", "original_collar_upper",
                "active_body_not_evaluated_on_flat_branch",
            }),
        ))

    def _inverse(self, row: Mapping[str, Any], scope: tuple[tuple[Any, Any], ...]) -> int:
        if row.get("true_prescribed_source_C5_implicit_inverse") is not True:
            raise UnsupportedSourceOperation("inverse is not the prescribed source inverse")
        if row.get("branch_inverse_is_original_not_support_cap") is not True:
            raise UnsupportedSourceOperation("inverse original branch is not source-bound")
        if row.get("flat_inverse") != "2*pi*phase" or row.get("phase_endpoints") != [0, 1]:
            raise UnsupportedSourceOperation("inverse phase endpoints or flat branch changed")
        certificate = row.get("derivative_positive_certificate")
        if not isinstance(certificate, str) or not certificate.strip():
            certificate = row.get("positive_certificate")
        if not isinstance(certificate, str) or not certificate.strip():
            raise UnsupportedSourceOperation("inverse lacks positive derivative certificate")
        variable = row.get("angle_variable")
        if not isinstance(variable, (str, int)):
            raise UnsupportedSourceOperation("inverse lacks angle binder")
        depth = sum(1 for key, _ in scope if isinstance(key, tuple) and key[:1] == ("binder",))
        bound = ("bound", depth)
        phase_scope = scope + ((("binder", "angle"), bound),)
        phase_scope = _bind_scope(phase_scope, variable, self.intern(bound), row, self)
        phase = self.walk(row["phase_function"], phase_scope)
        target = self.walk(row["target_phase"], scope)
        lower = self.walk(row["angle_lower"], scope)
        upper = self.walk(row["angle_upper"], scope)
        delta = self.walk(row["Delta"], scope) if _is_int(row.get("Delta")) else None
        eta = self.walk(row["eta"], scope) if _is_int(row.get("eta")) else None
        # These are literal phase values, not graph handles.
        endpoint_refs = (("rational", 0, 1), ("rational", 1, 1))
        return self.intern((
            "monotone_inverse", self.child_ref(phase), self.child_ref(target),
            self.child_ref(lower), self.child_ref(upper),
            self.child_ref(delta) if delta is not None else None,
            self.child_ref(eta) if eta is not None else None,
            endpoint_refs,
            self._metadata(row, {
                "operation", "phase_function", "target_phase", "angle_lower", "angle_upper",
                "Delta", "eta", "phase_endpoints", "angle_variable",
                "derivative_positive_certificate", "positive_certificate",
                "flat_inverse", "branch_inverse_is_original_not_support_cap",
                "true_prescribed_source_C5_implicit_inverse",
            }),
        ))

    def _substitute_inverse_angle(self, row: Mapping[str, Any],
                                  scope: tuple[tuple[Any, Any], ...]) -> int:
        inverse = self.walk(row["inverse_angle"], scope)
        variable = row.get("angle_variable")
        if not isinstance(variable, (str, int)):
            raise UnsupportedSourceOperation("inverse-angle substitution lacks binder")
        bound_scope = _bind_scope(scope, variable, inverse, row, self)
        return self.walk(row["body"], bound_scope)

    def _walk_leaf_or_fail(self, node_id: int, row: Mapping[str, Any], operation: str,
                           scope: tuple[tuple[Any, Any], ...]) -> int:
        refs_by_operation = {
            "current_original_axis_pressure_function_derivative":
                ("original_function", "Z_variable"),
            "current_original_leading_function_recipe": ("coordinate", "Z_variable"),
            "current_original_unique_analytic_core_component": ("rho", "Z_variable"),
            "current_original_source_parameter": (),
            "exact_positive_integer_expression": (),
            "mathematical_pi": (),
        }
        if operation not in refs_by_operation:
            raise UnsupportedSourceOperation("unsupported source operation: " + operation)
        if operation == "current_original_unique_analytic_core_component":
            for flag in ("actual_twenty_term_operator", "unique_analytic_Banach_solution",
                         "admitted_solution_not_a_model_or_enclosure_center"):
                if row.get(flag) is not True:
                    raise UnsupportedSourceOperation("unadmitted analytic core: " + flag)
        if operation in {"current_original_leading_function_recipe",
                         "current_original_axis_pressure_function_derivative"}:
            if row.get("defining_quantity_not_a_range_value") is not True:
                raise UnsupportedSourceOperation("source recipe is not a defining function")
        refs: list[str] = []
        for field in refs_by_operation[operation]:
            reference = row.get(field)
            if not _is_int(reference):
                raise UnsupportedSourceOperation(
                    "source operation reference is not a graph node: " + field
                )
            refs.append(self.child_ref(self.walk(reference, scope)))
        return self._canonical_leaf(node_id, row, operation, tuple(refs))

    def walk(self, node_id: int, scope: tuple[tuple[Any, Any], ...] = ()) -> int:
        _require(_is_int(node_id), "source graph child is not an integer node")
        _require(0 <= node_id < len(self.raw_nodes), "source graph child is out of range")
        signature = self._scope_signature(scope)
        memo_key = (node_id, signature)
        if memo_key in self.memo:
            return self.memo[memo_key]
        if memo_key in self.active:
            raise UnsupportedSourceOperation("cyclic source graph")
        self.active.add(memo_key)
        self.depth += 1
        if self.depth > 4096:
            self.depth -= 1
            self.active.remove(memo_key)
            raise UnsupportedSourceOperation("source graph nesting exceeds canonicalizer limit")
        try:
            row = self.raw_nodes[node_id]
            operation = row.get("operation") if isinstance(row, Mapping) else None
            _require(isinstance(operation, str), "source graph row lacks operation")
            if operation == "exact_rational":
                numerator, denominator = row.get("numerator"), row.get("denominator")
                _require(_is_int(numerator) and _is_int(denominator), "malformed exact rational")
                _require(denominator != 0, "exact rational has zero denominator")
                value = Fraction(numerator, denominator)
                result = self.intern(("rational", value.numerator, value.denominator))
            elif operation == "bound_variable":
                bound = self._scope_bound(scope, node_id)
                if bound is None:
                    name = row.get("name")
                    bound = self._scope_bound(scope, name) if isinstance(name, str) else None
                if bound is not None:
                    result = bound
                else:
                    result = self.intern((
                        "free_bound_source_atom", self.source_name, node_id,
                        str(row.get("name", "")),
                    ))
            elif operation == "sum":
                arguments = row.get("arguments")
                _require(isinstance(arguments, list), "malformed sum")
                result = self._sum([self.walk(child, scope) for child in arguments])
            elif operation == "product":
                arguments = row.get("arguments")
                _require(isinstance(arguments, list), "malformed product")
                result = self._product([self.walk(child, scope) for child in arguments])
            elif operation == "negative":
                result = self._product([self.intern(("rational", -1, 1)),
                                        self.walk(row["argument"], scope)])
            elif operation == "positive_quotient":
                result = self._positive_quotient(
                    self.walk(row["numerator"], scope),
                    self.walk(row["denominator"], scope),
                    row,
                )
            elif operation == "function_substitution":
                value = self.walk(row["value"], scope)
                bound_scope = _bind_scope(scope, row.get("variable"), value, row, self)
                result = self.walk(row["expression"], bound_scope)
            elif operation == "definite_integral":
                lower = self.walk(row["lower"], scope)
                upper = self.walk(row["upper"], scope)
                if self.digest(lower) == self.digest(upper):
                    result = self.intern(("rational", 0, 1))
                else:
                    variable = row.get("variable")
                    if not isinstance(variable, (str, int)):
                        raise UnsupportedSourceOperation("integral lacks binder")
                    depth = sum(1 for key, _ in scope if isinstance(key, tuple) and key[:1] == ("binder",))
                    binder = ("bound", depth)
                    integrand_scope = scope + ((("binder", "integral"), binder),)
                    bound_scope = _bind_scope(integrand_scope, variable, self.intern(
                        ("bound", depth)
                    ), row, self)
                    integrand = self.walk(row["integrand"], bound_scope)
                    result = self.intern((
                        "definite_integral", self.child_ref(integrand),
                        self.child_ref(lower), self.child_ref(upper), ("bound", depth),
                        self._metadata(row, {
                            "operation", "integrand", "lower", "upper", "variable",
                        }),
                    ))
            elif operation == "analytic_unary":
                result = self._analytic(row, self.walk(row["argument"], scope))
            elif operation == "exact_original_source_flat_sigma":
                result = self._flat_sigma(row, self.walk(row["argument"], scope))
            elif operation in {"original_flat_sigma", "original_flat_sigma_derivative",
                               "original_flat_sigma_prime", "original_flat_sigma_second_derivative"}:
                argument = self.walk(row["argument"], scope)
                if operation != "original_flat_sigma":
                    value = self.rational(argument)
                    if value is not None and (value <= 0 or value >= 1):
                        result = self.intern(("rational", 0, 1))
                    else:
                        result = self.intern((operation, self.child_ref(argument),
                                              self._metadata(row, {"operation", "argument"})))
                else:
                    result = self._flat_sigma(row, argument)
            elif operation == "exact_real_comparison":
                result = self._comparison(row, scope)
            elif operation == "logical_or":
                arguments = row.get("arguments")
                _require(isinstance(arguments, list), "malformed logical_or")
                result = self._logical_or([self.walk(child, scope) for child in arguments])
            elif operation == "original_lazy_flat_branch":
                result = self._lazy_branch(row, scope)
            elif operation == "original_monotone_phase_inverse":
                result = self._inverse(row, scope)
            elif operation == "substitute_original_inverse_angle":
                result = self._substitute_inverse_angle(row, scope)
            else:
                result = self._walk_leaf_or_fail(node_id, row, operation, scope)
            self.memo[memo_key] = result
            return result
        finally:
            self.depth -= 1
            self.active.remove(memo_key)

    def export(self) -> dict[str, Any]:
        rows = [
            {"id": index, "key": _jsonable(key), "digest": self.digests[index]}
            for index, key in enumerate(self.keys)
        ]
        return {
            "source_name": self.source_name,
            "node_count": len(rows),
            "nodes_digest": _plain_digest(rows),
            "nodes": rows,
        }


def _record_payloads(source_records: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    payloads: dict[str, Any] = {}
    for name, record in source_records.items():
        _require(isinstance(record, Mapping), "source record must be a mapping")
        expected = RECORD_HASHES.get(name)
        if expected is not None:
            _require(record.get("raw_sha256") == expected, "authenticated source hash changed: " + name)
        payloads[name] = _parse_record(record)
    return payloads


def _binding(payload: Any, key: str) -> Mapping[str, Any] | None:
    return _mapping_with_key(payload, key)


def _physical_roots(payload: Any) -> tuple[list[Mapping[str, Any]] | None, dict[str, Any] | None]:
    binding = _binding(payload, "modified_Q4_pressure_C3_bindings")
    if binding is None:
        return None, None
    charts = binding.get("charts")
    if not isinstance(charts, Mapping):
        return None, None
    return _graph_nodes(payload), charts


def _endpoint_rows(charts: Mapping[str, Any], chart: str) -> dict[str, list[list[int]]]:
    packet = charts.get(chart)
    _require(isinstance(packet, Mapping), "missing chart packet: " + chart)
    physical = packet.get("physical")
    _require(isinstance(physical, Mapping), "missing physical chart ports: " + chart)
    result: dict[str, list[list[int]]] = {}
    for field in ("U", "m", "Q"):
        rows = physical.get(field)
        _require(isinstance(rows, list), "missing physical field: " + field)
        _require(all(_is_int(item) for item in rows), "malformed physical source jet: " + field)
        result[field] = [[item] for item in rows]
    return result


def _loop_radii(payload: Any) -> tuple[list[Mapping[str, Any]] | None, dict[str, Any] | None]:
    binding = _binding(payload, "prescribed_loop_C5_bindings")
    if binding is None:
        return None, None
    charts = binding.get("charts")
    if not isinstance(charts, Mapping):
        return None, None
    return _graph_nodes(payload), charts


def _compare_vector(store: CanonicalStore, roots: list[int], coordinate: int,
                    endpoint: int) -> dict[str, Any]:
    endpoint_value = store.walk(endpoint)
    scope = _bind_scope((), coordinate, endpoint_value, {}, store)
    handles = [store.walk(root, scope) for root in roots]
    vector = store.intern(("source_port_vector", tuple(store.child_ref(handle) for handle in handles)))
    return {
        "source_node_ids": list(roots),
        "canonical_node_ids": handles,
        "canonical_node_id": vector,
        "canonical_digest": store.digest(vector),
        "canonical_handle_digests": [store.digest(handle) for handle in handles],
    }


def _endpoint_comparisons(
    physical_store: CanonicalStore | None,
    physical_rows: dict[str, dict[str, list[list[int]]]] | None,
    radius_store: CanonicalStore | None,
    radius_rows: dict[str, int] | None,
) -> tuple[dict[str, Any], list[str]]:
    comparisons: dict[str, Any] = {}
    canonicalization_errors: list[str] = []
    scope = {
        "second_micro": {"coordinate_node": 66988, "endpoint_graph_handle": 48},
        "frozen_macro": {"coordinate_node": 70265, "endpoint_graph_handle": 0},
    }
    orders = {"U": 2, "m": 3, "Q": 2}
    for field, count in orders.items():
        for order in range(count):
            label = f"{field}_Z{order}"
            if physical_store is None or physical_rows is None:
                comparisons[label] = {
                    "field": field, "ordinary_Z_order": order, "equal": None,
                    "left": None, "right": None, "source_endpoint_scope": scope,
                }
                continue
            try:
                left = _compare_vector(physical_store, physical_rows["second_micro"][field][order], 66988, 48)
                right = _compare_vector(physical_store, physical_rows["frozen_macro"][field][order], 70265, 0)
            except (IndexError, KeyError, TypeError, ValueError, UnsupportedSourceOperation) as error:
                reason = f"{type(error).__name__}: {error}"
                canonicalization_errors.append(f"{label}: {reason}")
                comparisons[label] = {
                    "field": field,
                    "ordinary_Z_order": order,
                    "left": None,
                    "right": None,
                    "equal": None,
                    "numeric_value": None,
                    "comparison_kind": "actual_physical_source_port_canonical_digest",
                    "canonicalization_error": reason,
                }
                continue
            comparisons[label] = {
                "field": field,
                "ordinary_Z_order": order,
                "left": dict(scope["second_micro"], **left, chart="second_micro"),
                "right": dict(scope["frozen_macro"], **right, chart="frozen_macro"),
                "equal": left["canonical_digest"] == right["canonical_digest"],
                "numeric_value": None,
                "comparison_kind": "actual_physical_source_port_canonical_digest",
            }
    for order in (0,):
        label = f"R_Z{order}"
        if radius_store is None or radius_rows is None:
            comparisons[label] = {
                "field": "R", "ordinary_Z_order": order, "equal": None,
                "left": None, "right": None, "source_endpoint_scope": scope,
            }
            continue
        try:
            left = _compare_vector(radius_store, [radius_rows["second_micro"]], 66988, 48)
            right = _compare_vector(radius_store, [radius_rows["frozen_macro"]], 70265, 0)
        except (IndexError, KeyError, TypeError, ValueError, UnsupportedSourceOperation) as error:
            reason = f"{type(error).__name__}: {error}"
            canonicalization_errors.append(f"{label}: {reason}")
            comparisons[label] = {
                "field": "R",
                "ordinary_Z_order": order,
                "left": None,
                "right": None,
                "equal": None,
                "numeric_value": None,
                "comparison_kind": "actual_radius_source_port_canonical_digest",
                "canonicalization_error": reason,
                "Z_independent_source": True,
            }
            continue
        comparisons[label] = {
            "field": "R",
            "ordinary_Z_order": order,
            "left": dict(scope["second_micro"], **left, chart="second_micro"),
            "right": dict(scope["frozen_macro"], **right, chart="frozen_macro"),
            "equal": left["canonical_digest"] == right["canonical_digest"],
            "numeric_value": None,
            "comparison_kind": "actual_radius_source_port_canonical_digest",
            "Z_independent_source": True,
        }
    return comparisons, canonicalization_errors


def _source_check(payloads: Mapping[str, Any], needle: str) -> bool:
    return any(_find_truth(payload, needle) for payload in payloads.values())


def derive(source_records: Mapping[str, Mapping[str, Any]], *, expected_input_digest: str) -> dict[str, Any]:
    """Canonicalize the authenticated source endpoint ports exactly once.

    The caller owns the once-entry guard and supplies decoded, raw-authenticated
    records.  This function performs no numerical evaluation; ``passed`` means
    that all local source comparisons and prerequisite gates were accepted.
    """

    _require(isinstance(source_records, Mapping), "source_records must be a mapping")
    _require(_plain_digest(source_records) == expected_input_digest,
             "source input digest changed before canonicalization")
    _require(set(RECORD_NAMES).issubset(source_records), "new source record inventory is incomplete")
    payloads = _record_payloads(source_records)

    modified_name = next((name for name in source_records if "modified_Q4_pressure_C3_functions" in name), None)
    loop_name = next((name for name in source_records if "prescribed_loop_C5_functions" in name), None)
    _require(modified_name is not None and loop_name is not None, "required source records are missing")
    modified_nodes, modified_charts = _physical_roots(payloads[modified_name])
    loop_nodes, loop_charts = _loop_radii(payloads[loop_name])

    physical_store: CanonicalStore | None = None
    physical_rows: dict[str, dict[str, list[list[int]]]] | None = None
    physical_graph_digest: str | None = None
    canonical_stores: dict[str, CanonicalStore] = {}
    raw_graphs: dict[str, dict[str, Any]] = {}
    if modified_nodes is not None and modified_charts is not None:
        physical_rows = {
            chart: _endpoint_rows(modified_charts, chart)
            for chart in ("second_micro", "frozen_macro")
        }
        physical_store = CanonicalStore(modified_nodes, modified_name)
        canonical_stores[modified_name] = physical_store
        physical_graph_digest = _source_node_digest(modified_nodes)
        raw_graphs[modified_name] = {
            "node_count": len(modified_nodes),
            "nodes_digest": physical_graph_digest,
        }

    radius_store: CanonicalStore | None = None
    radius_rows: dict[str, int] | None = None
    radius_graph_digest: str | None = None
    if loop_nodes is not None and loop_charts is not None:
        radius_rows = {}
        for chart in ("second_micro", "frozen_macro"):
            packet = loop_charts.get(chart)
            _require(isinstance(packet, Mapping) and _is_int(packet.get("radius")),
                     "loop chart lacks actual radius source port: " + chart)
            radius_rows[chart] = packet["radius"]
        radius_store = CanonicalStore(loop_nodes, loop_name)
        canonical_stores[loop_name] = radius_store
        radius_graph_digest = _source_node_digest(loop_nodes)
        raw_graphs[loop_name] = {"node_count": len(loop_nodes), "nodes_digest": radius_graph_digest}

    comparisons, canonicalization_errors = _endpoint_comparisons(
        physical_store, physical_rows, radius_store, radius_rows
    )
    equal_values = {key: value.get("equal") is True for key, value in comparisons.items()}
    stage_flags = {
        STAGE_FLAGS[0]: all(equal_values.get(f"U_Z{order}", False) for order in (0, 1)),
        STAGE_FLAGS[1]: all(equal_values.get(f"m_Z{order}", False) for order in (0, 1, 2)),
        STAGE_FLAGS[2]: all(equal_values.get(f"Q_Z{order}", False) for order in (0, 1)),
        STAGE_FLAGS[3]: equal_values.get("R_Z0", False),
    }

    source_families = [
        payload.get("source_family")
        for payload in payloads.values()
        if isinstance(payload, Mapping) and isinstance(payload.get("source_family"), Mapping)
    ]
    same_source_family = bool(source_families) and len({
        _plain_digest(item) for item in source_families
    }) == 1
    parent = payloads[RECORD_NAMES[0]]
    weighted = payloads[RECORD_NAMES[1]]
    bridge = payloads[RECORD_NAMES[4]]
    loop_binding = payloads[loop_name]['prescribed_loop_C5_bindings']
    prior_retained = parent['retained_construction_data']
    audit = {
        "frozen_input_digest": _plain_digest(source_records) == expected_input_digest,
        "exact_five_record_inventory": set(source_records) == set(RECORD_NAMES),
        "all_records_raw_authenticated_before_decode": all(
            record.get("raw_authenticated_before_decode") is True
            for record in source_records.values()),
        "same_exact_source_family": same_source_family,
        "accepted_parent_J0_retained": parent['component_result']['passed'] is True
            and parent['component_result']['source_interface']['first_jump_zero'] is True,
        "actual_weighted_J1_source_bound": weighted['definition_graph']['nodes'][521]['arguments'] == [355, 345],
        "native_loop_is_exact_MOD_prefix": modified_nodes[:len(loop_nodes)] == loop_nodes,
        "native_bridge_is_exact_loop_prefix": loop_nodes[:len(bridge['exact_graph_nodes'])] == bridge['exact_graph_nodes'],
        "actual_native_coordinates": loop_binding['charts']['second_micro']['coordinate'] == 66988
            and loop_binding['charts']['frozen_macro']['coordinate'] == 70265,
        "actual_endpoint_graph_rationals": modified_nodes[48] == dict(operation='exact_rational', numerator=2, denominator=1)
            and modified_nodes[0] == dict(operation='exact_rational', numerator=0, denominator=1),
        "actual_nonreset_m_history": loop_binding['charts']['frozen_macro']['incoming'] == loop_binding['endpoints']['second_micro'],
        "construction_data_frozen": prior_retained['cumulative_terms'] == 212
            and prior_retained['raw_signed_terms'] == 213 and prior_retained['reduced_terms'] == 183
            and prior_retained['six_signed_channels_retained'] is True
            and prior_retained['all_old_error_data_retained'] is True
            and prior_retained['reserved_contribution_registered_once'] is True
            and prior_retained['inherited_errors_deleted'] is False,
    }
    admitted_sources = all(audit.values())
    stage_flags = {name: value and admitted_sources for name, value in stage_flags.items()}

    unresolved: list[str] = []
    if physical_store is None:
        unresolved.append("actual_physical_U_m_Q_source_ports_missing")
    if radius_store is None:
        unresolved.append("actual_physical_R_source_ports_missing")
    unresolved.extend(
        "endpoint_canonicalization_unresolved: " + error
        for error in canonicalization_errors
    )
    for field, gate in (
        ("U", "actual_native_micro_macro_U_endpoint_identity_unproved"),
        ("m", "actual_native_micro_macro_m_endpoint_identity_unproved"),
        ("Q", "actual_native_micro_macro_Q_endpoint_identity_unproved"),
        ("R", "actual_native_micro_macro_R_endpoint_identity_unproved"),
    ):
        if not stage_flags[{
            "U": STAGE_FLAGS[0], "m": STAGE_FLAGS[1],
            "Q": STAGE_FLAGS[2], "R": STAGE_FLAGS[3],
        }[field]]:
            unresolved.append(gate)
    unresolved.extend("source_prerequisite_unresolved: " + name for name, value in audit.items() if not value)
    unresolved = list(dict.fromkeys(unresolved))

    all_identity = all(stage_flags.values())
    remaining_jump_nodes = list(range(521, 536))
    source_interface = {
        "parent_jump_node": 521,
        "candidate_J1_node": 521,
        "candidate_J1_arguments": [355, 345],
        "second_micro_coordinate_node": 66988,
        "second_micro_endpoint_graph_handle": 48,
        "frozen_macro_coordinate_node": 70265,
        "frozen_macro_endpoint_graph_handle": 0,
        "actual_L_node": 25272,
        "actual_native_physical_U_m_Q_ports": physical_rows,
        "actual_native_radius_ports": radius_rows,
        "ordinary_Z_orders": [0, 1],
        "owned_Z_cell": ["0", "1/2"],
        "remaining_parent_jump_nodes": remaining_jump_nodes,
        "old_archived_16_J_terms_unchanged": True,
        "old_packet_errors_deleted": False,
        "common_U_Q_binding_used": False,
        "canonical_difference_does_not_prove_function_inequality": True,
        "first_parent_J0_proof_retained": True,
    }
    retained = dict(prior_retained)

    definition_graph = {
        "namespace": SCHEMA,
        "node_count": sum(len(store.keys) for store in canonical_stores.values()),
        "nodes_digest": _plain_digest({name: store.export()["nodes_digest"] for name, store in canonical_stores.items()}),
        "raw_graphs": raw_graphs,
        "canonical_graphs": {name: store.export() for name, store in canonical_stores.items()},
        "canonical_graph_digest": _plain_digest({
            name: store.export()["nodes_digest"] for name, store in canonical_stores.items()
        }),
        "compared_endpoint_roots": comparisons,
        "interned_dag": True,
        "memo_scope": "(source_node, lexical_scope)",
        "exact_fraction_rationals": True,
        "no_expanded_polynomial": True,
        "no_quadrature": True,
        "no_numeric_values": True,
    }

    result: dict[str, Any] = {
        "schema": SCHEMA,
        "component": "micro_macro_endpoint",
        "passed": bool(all_identity and not unresolved),
        "source_only": True,
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z_not_temporal_n": True,
        "input_digest": expected_input_digest,
        "source_family": source_families[0] if same_source_family else None,
        "source_prerequisite_audit": audit,
        "endpoint_comparisons": comparisons,
        "definition_graph": definition_graph,
        "source_interface": source_interface,
        "retained_construction_data": retained,
        "unresolved_gates": unresolved,
        "remaining_construction_gates": [
            "actual_native_micro_macro_B_jump_source_identity",
            "fourteen_other_native_B_jump_bindings",
            "native_old_active_density_identity",
            "terminal_difference_zero",
            "global_five_moment_closure",
        ],
        "source_only_scope": "actual physical source ports only; no common U/Q seed",
        "numeric_values_computed": False,
        "old_error_terms_removed": False,
        "old_operators_replayed": False,
        "runtime_model_metadata": None,
        "source_inspection_audit": {
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "old_math_operator_imports": 0,
            "numeric_evaluation_calls": 0,
            "quadrature_calls": 0,
            "canonicalization_only": True,
        },
        "stage_flags": list(STAGE_FLAGS),
        "false_flags": list(FALSE_FLAGS),
        "identity_flags": stage_flags,
    }
    result.update(stage_flags)
    result.update({flag: False for flag in FALSE_FLAGS})
    return result


__all__ = [
    "SCHEMA",
    "RECORD_NAMES",
    "RECORD_HASHES",
    "STAGE_FLAGS",
    "FALSE_FLAGS",
    "derive",
]
