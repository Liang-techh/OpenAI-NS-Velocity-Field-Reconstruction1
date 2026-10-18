"""Read-only provenance checking for the existing Rm native derivative DAG.

The verifier consumes an already saved graph and chart descriptors.  It never
imports a producer, constructs a graph, evaluates a function, decodes an
archive, or runs an old derivative implementation.  Every expected expression
is resolved by descriptor lookup in the supplied graph, so a successful row
always points at an existing source and target handle.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
import hashlib
import json
from typing import Any


DEFAULT_ORDERS = (0, 1)
_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
_FIELD_ALIASES = {
    "E": ("E", "e", "E_y"),
    "V": ("V", "v", "V_y"),
    "m": ("m", "M", "m_y"),
    "h": ("h", "H", "h_y"),
    "k": ("k", "K", "k_y"),
    "e": ("e", "epsilon", "e_y"),
    "p": ("p", "P", "p_y"),
    "pressure": ("pressure",),
    "P0": ("P0", "p0", "P_0"),
}
_F_SUB_GUARDS = (
    "Z_independent_substitution",
    "exact_native_chain_rule_substitution",
)
_VOLATILE_METADATA = frozenset(
    {
        "source_formula",
        "source_binding",
        "source_native_substitution",
        "physical_Jacobian_node",
        "physical_Jacobian_at_endpoint",
    }
)
_NATIVE_DEPENDENCY_KEYS = (
    "native_dependent",
    "native_dependency",
    "native_dependent_source",
    "depends_on_native",
    "depends_on_phase",
    "phase_dependent",
)
_SOURCE_INDEPENDENT_KEYS = (
    "source_opaque",
    "sourceopaque",
    "native_independent",
    "source_native_independent",
    "coordinate_independent",
    "Z_independent",
    "Z_independent_source",
)


class ProvenanceError(ValueError):
    """Malformed input or an explicitly unresolved structural obligation."""


class _OpenRule(ProvenanceError):
    """A reachable native source rule is outside this bounded verifier."""


def canonical_digest(value: Any) -> str:
    """Return the canonical SHA-256 digest used for saved graph bindings."""

    def default(item: Any) -> Any:
        if isinstance(item, Fraction):
            return {"__fraction__": [item.numerator, item.denominator]}
        raise TypeError(f"value is not JSON serializable: {type(item).__name__}")

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
        default=default,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProvenanceError(message)


def _sequence(value: Any, label: str) -> list[Any]:
    _require(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
        f"{label} must be a sequence",
    )
    return list(value)


def _handle(value: Any, label: str, count: int) -> int:
    if isinstance(value, Mapping):
        value = value.get("node", value.get("handle", value.get("root")))
    _require(
        isinstance(value, int) and not isinstance(value, bool),
        f"{label} must be an integer graph handle",
    )
    _require(0 <= value < count, f"{label} is outside the source graph")
    return value


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, bool):
        raise ProvenanceError(f"{label} is not an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str):
        text = value.strip()
        try:
            if "/" in text:
                pieces = text.split("/")
                _require(len(pieces) == 2, f"{label} is not an exact rational")
                return Fraction(int(pieces[0]), int(pieces[1]))
            return Fraction(int(text), 1)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ProvenanceError(f"{label} is not an exact rational") from exc
    raise ProvenanceError(f"{label} is not an exact rational")


def _freeze(value: Any) -> Any:
    if isinstance(value, Fraction):
        return ("fraction", value.numerator, value.denominator)
    if isinstance(value, Mapping):
        return tuple(sorted((str(key), _freeze(item)) for key, item in value.items()))
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return tuple(_freeze(item) for item in value)
    return value


def _op(nodes: Sequence[Any], index: int) -> str:
    row = nodes[index]
    _require(isinstance(row, Mapping), f"node {index} is not a mapping")
    operation = row.get("operation", row.get("op"))
    _require(isinstance(operation, str) and operation, f"node {index} has no operation")
    return operation


def _row(nodes: Sequence[Any], index: int, label: str = "node") -> Mapping[str, Any]:
    _handle(index, label, len(nodes))
    value = nodes[index]
    _require(isinstance(value, Mapping), f"{label} {index} is not a mapping")
    return value


def _args(nodes: Sequence[Any], index: int) -> list[int]:
    row = _row(nodes, index)
    raw = row.get("arguments", row.get("operands"))
    if raw is None:
        raise ProvenanceError(f"node {index}.arguments is missing")
    return [
        _handle(value, f"node {index}.arguments[{offset}]", len(nodes))
        for offset, value in enumerate(_sequence(raw, f"node {index}.arguments"))
    ]


def _child(nodes: Sequence[Any], row: Mapping[str, Any], key: str) -> int:
    return _handle(row.get(key), f"node.{key}", len(nodes))


def _exact(nodes: Sequence[Any], index: int) -> Fraction | None:
    if _op(nodes, index) != "exact_rational":
        return None
    row = _row(nodes, index)
    try:
        numerator = _fraction(row.get("numerator", row.get("value")), "rational numerator")
        denominator = _fraction(row.get("denominator", 1), "rational denominator")
        if denominator == 0:
            return None
        return numerator / denominator
    except (ProvenanceError, ZeroDivisionError):
        return None


def _bool_true(row: Mapping[str, Any], key: str) -> bool:
    return row.get(key) is True


def _native_dependent(row: Mapping[str, Any]) -> bool:
    for key in _NATIVE_DEPENDENCY_KEYS:
        value = row.get(key)
        if value is True:
            return True
        if isinstance(value, str) and value.lower() not in {"", "false", "none", "0", "no"}:
            return True
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and value:
            return True
    return False


def _source_independent(row: Mapping[str, Any]) -> bool:
    return any(_bool_true(row, key) for key in _SOURCE_INDEPENDENT_KEYS)


def _bound_name(nodes: Sequence[Any], value: Any, label: str = "bound variable") -> str:
    index = _handle(value, label, len(nodes))
    row = _row(nodes, index, label)
    _require(_op(nodes, index) == "bound_variable", f"{label} {index} is not a bound variable")
    name = row.get("name")
    _require(isinstance(name, str) and name, f"{label} {index} has no name")
    return name


def _field_container(chart: Mapping[str, Any]) -> Mapping[str, Any] | None:
    owners: list[Mapping[str, Any]] = [chart]
    for key in ("source_chart", "leading", "leading_chart", "source"):
        child = chart.get(key)
        if isinstance(child, Mapping):
            owners.append(child)
    for owner in owners:
        for key in ("fields", "leading_fields", "leadingFields"):
            child = owner.get(key)
            if isinstance(child, Mapping):
                return child
    if any(key in chart for key in _FIELDS):
        return chart
    return None


def _first(owner: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in owner:
            return owner[key]
    return None


def _rows(raw: Any, label: str, count: int, orders: Sequence[int]) -> dict[int, int]:
    result: dict[int, int] = {}
    if isinstance(raw, Mapping):
        for order in orders:
            value = raw.get(order, raw.get(str(order)))
            _require(value is not None, f"{label} lacks ordinary order {order}")
            result[order] = _handle(value, f"{label}[{order}]", count)
        return result
    values = _sequence(raw, label)
    for order in orders:
        _require(order < len(values), f"{label} lacks ordinary order {order}")
        result[order] = _handle(values[order], f"{label}[{order}]", count)
    return result


def _field_value(container: Mapping[str, Any], field: str) -> Any:
    for key in _FIELD_ALIASES.get(field, (field,)):
        if key in container:
            return container[key]
    return None


def _container_value(chart: Mapping[str, Any], key: str) -> Any:
    nested = chart.get("loop")
    if isinstance(nested, Mapping) and key in nested:
        return nested[key]
    return chart.get(key)


def _semantic_key(nodes: Sequence[Any], index: int) -> tuple[Any, ...]:
    row = _row(nodes, index)
    operation = _op(nodes, index)
    if operation == "exact_rational":
        value = _exact(nodes, index)
        if value is None:
            raise ProvenanceError(f"node {index} has an invalid exact rational")
        return (operation, value.numerator, value.denominator)
    if operation == "bound_variable":
        return (operation, row.get("name"))
    if operation in ("sum", "product", "logical_or"):
        return (operation, tuple(_args(nodes, index)))
    if operation == "negative":
        return (operation, _child(nodes, row, "argument"))
    if operation in ("positive_quotient",):
        return (
            operation,
            _child(nodes, row, "numerator"),
            _child(nodes, row, "denominator"),
            _freeze(row.get("source_positive_certificate")),
        )
    if operation in ("analytic_unary",):
        return (operation, row.get("name"), _child(nodes, row, "argument"))
    if operation == "exact_original_source_flat_sigma":
        return (operation, _child(nodes, row, "argument"))
    if operation in ("function_substitution",):
        return (
            operation,
            _child(nodes, row, "expression"),
            _child(nodes, row, "variable"),
            _child(nodes, row, "value"),
        )
    if operation == "definite_integral":
        variable = row.get("variable")
        if isinstance(variable, Mapping):
            variable = _freeze(variable)
        return (
            operation,
            _child(nodes, row, "integrand"),
            variable,
            _child(nodes, row, "lower"),
            _child(nodes, row, "upper"),
        )
    if operation in ("substitute_original_inverse_angle",):
        return (
            operation,
            _child(nodes, row, "body"),
            _child(nodes, row, "inverse_angle"),
            row.get("angle_variable"),
        )
    if operation == "original_monotone_phase_inverse":
        return (
            operation,
            _child(nodes, row, "phase_function"),
            _child(nodes, row, "target_phase"),
            row.get("angle_variable"),
        )
    if operation == "original_lazy_flat_branch":
        return (
            operation,
            _child(nodes, row, "active_body"),
            _child(nodes, row, "Delta"),
            _child(nodes, row, "eta"),
            _child(nodes, row, "flat_value"),
            _child(nodes, row, "flat_predicate"),
        )
    if operation == "exact_real_comparison":
        return (
            operation,
            row.get("operator"),
            _child(nodes, row, "left"),
            _child(nodes, row, "right"),
        )
    # This key is used only for diagnostics.  Expected nodes for unknown
    # operations are never synthesized or accepted by this verifier.
    return (
        operation,
        tuple(
            sorted(
                (str(key), _freeze(value))
                for key, value in row.items()
                if key not in _VOLATILE_METADATA and key not in {"operation", "op"}
            )
        ),
    )


class _DescriptorIndex:
    """Read-only semantic descriptor index over the supplied node list."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]]) -> None:
        self.nodes = nodes
        self.by_key: dict[tuple[Any, ...], list[int]] = {}
        self.errors: list[str] = []
        for index, row in enumerate(nodes):
            if not isinstance(row, Mapping):
                self.errors.append(f"node {index} is not a mapping")
                continue
            try:
                key = _semantic_key(nodes, index)
            except ProvenanceError as exc:
                self.errors.append(str(exc))
                continue
            self.by_key.setdefault(key, []).append(index)

    def _find(
        self,
        expected: Mapping[str, Any],
        *,
        required: Mapping[str, Any] | None = None,
        preserved_from: Mapping[str, Any] | None = None,
        changed: Sequence[str] = (),
    ) -> int:
        # Build a temporary descriptor key without adding it to the graph.
        operation = expected.get("operation", expected.get("op"))
        if not isinstance(operation, str):
            raise _OpenRule("expected descriptor has no operation")
        key = self._key_from_expected(expected)
        candidates = list(self.by_key.get(key, ()))
        required = {} if required is None else dict(required)
        candidates = [
            index
            for index in candidates
            if all(self.nodes[index].get(name) == value for name, value in required.items())
        ]
        if preserved_from is not None:
            ignored = {"operation", "op", *changed}
            candidates = [
                index
                for index in candidates
                if all(
                    self.nodes[index].get(name) == value
                    for name, value in preserved_from.items()
                    if name not in ignored
                )
            ]
        if not candidates:
            raise _OpenRule(f"existing descriptor not found: {self._brief(expected)}")
        if len(candidates) > 1:
            # Hash-consed source graphs normally have one descriptor.  Do not
            # silently pick a metadata variant when provenance is ambiguous.
            raise _OpenRule(
                f"existing descriptor is ambiguous ({len(candidates)} matches): {self._brief(expected)}"
            )
        return candidates[0]

    def _key_from_expected(self, row: Mapping[str, Any]) -> tuple[Any, ...]:
        operation = row.get("operation", row.get("op"))
        _require(isinstance(operation, str), "expected descriptor has no operation")
        if operation == "exact_rational":
            value = _fraction(row.get("numerator", row.get("value")), "expected rational numerator")
            denominator = _fraction(row.get("denominator", 1), "expected rational denominator")
            value /= denominator
            return (operation, value.numerator, value.denominator)
        if operation == "bound_variable":
            return (operation, row.get("name"))
        if operation in ("sum", "product", "logical_or"):
            return (operation, tuple(row.get("arguments", row.get("operands", ()))))
        if operation == "negative":
            return (operation, row.get("argument"))
        if operation == "positive_quotient":
            return (
                operation,
                row.get("numerator"),
                row.get("denominator"),
                _freeze(row.get("source_positive_certificate")),
            )
        if operation == "analytic_unary":
            return (operation, row.get("name"), row.get("argument"))
        if operation == "exact_original_source_flat_sigma":
            return (operation, row.get("argument"))
        if operation == "function_substitution":
            return (operation, row.get("expression"), row.get("variable"), row.get("value"))
        if operation == "definite_integral":
            return (
                operation,
                row.get("integrand"),
                row.get("variable"),
                row.get("lower"),
                row.get("upper"),
            )
        if operation == "substitute_original_inverse_angle":
            return (operation, row.get("body"), row.get("inverse_angle"), row.get("angle_variable"))
        if operation == "original_monotone_phase_inverse":
            return (
                operation,
                row.get("phase_function"),
                row.get("target_phase"),
                row.get("angle_variable"),
            )
        if operation == "original_lazy_flat_branch":
            return (
                operation,
                row.get("active_body"),
                row.get("Delta"),
                row.get("eta"),
                row.get("flat_value"),
                row.get("flat_predicate"),
            )
        if operation == "exact_real_comparison":
            return (operation, row.get("operator"), row.get("left"), row.get("right"))
        return (
            operation,
            tuple(
                sorted(
                    (str(key), _freeze(value))
                    for key, value in row.items()
                    if key not in _VOLATILE_METADATA and key not in {"operation", "op"}
                )
            ),
        )

    @staticmethod
    def _brief(row: Mapping[str, Any]) -> str:
        operation = row.get("operation", row.get("op"))
        return f"{operation} {dict(row)}"

    def find(self, operation: str, fields: Mapping[str, Any], **kwargs: Any) -> int:
        row = {"operation": operation, **dict(fields)}
        return self._find(row, **kwargs)

    def rational(self, value: Fraction | int) -> int:
        value = Fraction(value)
        return self.find(
            "exact_rational",
            {"numerator": value.numerator, "denominator": value.denominator},
        )

    def bound(self, name: str) -> int:
        return self.find("bound_variable", {"name": name})


class _Oracle:
    """Native derivative rules that return only handles already in ``nodes``."""

    def __init__(
        self,
        nodes: Sequence[Mapping[str, Any]],
        index: _DescriptorIndex,
        *,
        coordinate: int,
        phase: int,
        seeds: Mapping[int, int],
        expected_source_family: Any = None,
    ) -> None:
        self.nodes = nodes
        self.index = index
        self.coordinate = coordinate
        self.phase = phase
        self.main_context = ("node", coordinate)
        self.seeds = dict(seeds)
        self.expected_source_family = expected_source_family
        self.cache: dict[tuple[int, tuple[str, Any]], int] = {}
        self.free_cache: dict[int, frozenset[str]] = {}
        self.active: set[tuple[int, tuple[str, Any]]] = set()
        self.visited: dict[tuple[int, tuple[str, Any]], int] = {}

    def _ctx_name(self, context: tuple[str, Any]) -> str:
        kind, value = context
        if kind == "name":
            return str(value)
        return _bound_name(self.nodes, int(value), "derivative variable")

    def _ctx_matches(self, context: tuple[str, Any], node: int) -> bool:
        kind, value = context
        if kind == "node":
            return int(value) == node
        return self._ctx_name(context) == _bound_name(self.nodes, node, "bound variable")

    def _is_zero(self, node: int) -> bool:
        return _exact(self.nodes, node) == 0

    def _is_one(self, node: int) -> bool:
        return _exact(self.nodes, node) == 1

    def _add(self, *values: int) -> int:
        zero = self.index.rational(0)
        terms = [value for value in values if not self._is_zero(value)]
        if not terms:
            return zero
        if len(terms) == 1:
            return terms[0]
        return self.index.find("sum", {"arguments": terms})

    def _mul(self, *values: int) -> int:
        zero = self.index.rational(0)
        one = self.index.rational(1)
        if any(self._is_zero(value) for value in values):
            return zero
        factors = [value for value in values if not self._is_one(value)]
        if not factors:
            return one
        if len(factors) == 1:
            return factors[0]
        return self.index.find("product", {"arguments": factors})

    def _neg(self, value: int) -> int:
        if self._is_zero(value):
            return self.index.rational(0)
        return self.index.find("negative", {"argument": value})

    def _sub(self, left: int, right: int) -> int:
        return self._add(left, self._neg(right))

    def _quotient(self, numerator: int, denominator: int, certificate: Any) -> int:
        if self._is_zero(denominator):
            raise _OpenRule("native derivative produced a zero quotient denominator")
        if self._is_zero(numerator):
            return self.index.rational(0)
        return self.index.find(
            "positive_quotient",
            {
                "numerator": numerator,
                "denominator": denominator,
                "source_positive_certificate": certificate,
            },
        )

    def _substitution(self, expression: int, variable: int, value: int) -> int:
        return self.index.find(
            "function_substitution",
            {"expression": expression, "variable": variable, "value": value},
            required={key: True for key in _F_SUB_GUARDS},
        )

    def _inverse_at(self, body: int, inverse: int, angle: str) -> int:
        return self.index.find(
            "substitute_original_inverse_angle",
            {"body": body, "inverse_angle": inverse, "angle_variable": angle},
        )

    def _analytic(self, name: str, argument: int) -> int:
        return self.index.find("analytic_unary", {"name": name, "argument": argument})

    def _sigma_prime(self, argument: int) -> int:
        zero = self.index.rational(0)
        one = self.index.rational(1)
        two = self.index.rational(2)
        sigma = self.index.find("exact_original_source_flat_sigma", {"argument": argument}, required={
            "exact_function_not_its_interval_cap": True,
            "definition": "0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside",
            "source_binding": {
                "module": "lei_ren_part1_paper_compliant_outer_initial.py",
                "module_sha256": "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657",
                "function": "stable_sigma",
                "AST_sha256": "4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c",
            },
        })
        complement = self._sub(one, argument)
        first = self._add(
            self._quotient(
                two,
                self._mul(argument, argument, argument),
                "interior original sigma x>0",
            ),
            self._quotient(
                two,
                self._mul(complement, complement, complement),
                "interior original sigma x<1",
            ),
        )
        lower = self.index.find(
            "exact_real_comparison",
            {"operator": "<=", "left": argument, "right": zero},
        )
        upper = self.index.find(
            "exact_real_comparison",
            {"operator": ">=", "left": argument, "right": one},
        )
        predicate = self.index.find("logical_or", {"arguments": [lower, upper]})
        active = self._mul(sigma, self._sub(one, sigma), first)
        return self.index.find(
            "original_lazy_flat_branch",
            {
                "active_body": active,
                "Delta": argument,
                "eta": one,
                "flat_value": zero,
                "flat_predicate": predicate,
            },
            required={"original_sigma_exact_logistic_first_derivative": True,
                      "active_body_not_evaluated_on_flat_branch": True},
        )

    def _lazy_derivative(self, row: Mapping[str, Any], active: int, flat: int) -> int:
        expected = dict(row)
        expected["operation"] = "original_lazy_flat_branch"
        expected["active_body"] = active
        expected["flat_value"] = flat
        expected["same_original_lazy_native_derivative"] = True
        return self.index.find(
            "original_lazy_flat_branch",
            {
                "active_body": active,
                "Delta": row.get("Delta"),
                "eta": row.get("eta"),
                "flat_value": flat,
                "flat_predicate": row.get("flat_predicate"),
            },
            required={"same_original_lazy_native_derivative": True},
            preserved_from=expected,
            changed=("active_body", "flat_value", "same_original_lazy_native_derivative"),
        )

    def _integral_derivative(
        self,
        row: Mapping[str, Any],
        partial: int,
    ) -> int:
        expected = dict(row)
        expected["operation"] = "definite_integral"
        expected["integrand"] = partial
        expected["exact_native_partial_of_same_integrand"] = True
        return self.index.find(
            "definite_integral",
            {
                "integrand": partial,
                "variable": row.get("variable"),
                "lower": row.get("lower"),
                "upper": row.get("upper"),
            },
            required={"exact_native_partial_of_same_integrand": True},
            preserved_from=expected,
            changed=("integrand", "exact_native_partial_of_same_integrand"),
        )

    def _core_derivative(self, row: Mapping[str, Any], context: tuple[str, Any]) -> int:
        required = (
            "source_family",
            "unique_analytic_Banach_solution",
            "actual_twenty_term_operator",
            "admitted_solution_not_a_model_or_enclosure_center",
        )
        if row.get("source_family") is None or not all(row.get(key) is True for key in required[1:]):
            raise _OpenRule("analytic core source premises are not admitted")
        if self.expected_source_family is None or row["source_family"] != self.expected_source_family:
            raise _OpenRule("analytic core source family is not bound to the actual chart")
        terms: list[int] = []
        radial = int(row.get("radial_order", 0))
        ordinary = int(row.get("ordinary_slow_Z_derivative_order", 0))
        dr = self.derivative(_child(self.nodes, row, "rho"), context)
        dz = self.derivative(_child(self.nodes, row, "Z_variable"), context)
        if not self._is_zero(dr):
            if radial >= 2:
                raise _OpenRule("additional true radial core source is required")
            expected = dict(row)
            expected["operation"] = "current_original_unique_analytic_core_component"
            expected["radial_order"] = radial + 1
            new = self._lookup_copied_core(expected, changed=("radial_order",))
            terms.append(self._mul(new, dr))
        if not self._is_zero(dz):
            expected = dict(row)
            expected["operation"] = "current_original_unique_analytic_core_component"
            expected["ordinary_slow_Z_derivative_order"] = ordinary + 1
            new = self._lookup_copied_core(expected, changed=("ordinary_slow_Z_derivative_order",))
            terms.append(self._mul(new, dz))
        return self._add(*terms)

    def _lookup_copied_core(self, expected: Mapping[str, Any], changed: Sequence[str]) -> int:
        candidates: list[int] = []
        for index, candidate in enumerate(self.nodes):
            if not isinstance(candidate, Mapping) or _op(self.nodes, index) != expected.get("operation"):
                continue
            if all(
                candidate.get(key) == value
                for key, value in expected.items()
                if key not in {"operation", "op", *changed}
            ) and all(candidate.get(key) == expected.get(key) for key in changed):
                candidates.append(index)
        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise _OpenRule("copied analytic core component is absent from the existing DAG")
        raise _OpenRule("copied analytic core component is ambiguous")

    def _free(self, index: int) -> frozenset[str]:
        if index in self.free_cache:
            return self.free_cache[index]
        row = _row(self.nodes, index)
        operation = _op(self.nodes, index)
        if operation == "bound_variable":
            name = row.get("name")
            _require(isinstance(name, str) and name, f"bound variable {index} has no name")
            result = frozenset((name,))
        elif operation in (
            "exact_rational",
            "exact_positive_integer_expression",
            "mathematical_pi",
            "mathematical_positive_infinity",
        ):
            result = frozenset()
        elif operation in ("original_source_parameter", "current_original_source_parameter"):
            if row.get("Z_independent") is not True:
                raise _OpenRule(f"source parameter {index} lacks Z-independent certificate")
            result = frozenset()
        elif operation in ("sum", "product", "logical_or"):
            result = frozenset().union(*(self._free(child) for child in _args(self.nodes, index)))
        elif operation == "negative":
            result = self._free(_child(self.nodes, row, "argument"))
        elif operation == "positive_quotient":
            result = self._free(_child(self.nodes, row, "numerator")) | self._free(
                _child(self.nodes, row, "denominator")
            )
        elif operation in ("analytic_unary", "exact_original_source_flat_sigma"):
            result = self._free(_child(self.nodes, row, "argument"))
        elif operation == "exact_real_comparison":
            result = self._free(_child(self.nodes, row, "left")) | self._free(
                _child(self.nodes, row, "right")
            )
        elif operation == "function_substitution":
            expression = _child(self.nodes, row, "expression")
            variable = _child(self.nodes, row, "variable")
            value = _child(self.nodes, row, "value")
            variable_name = _bound_name(self.nodes, variable)
            original = set(self._free(expression))
            if variable_name in original:
                original.remove(variable_name)
                original.update(self._free(value))
            result = frozenset(original)
        elif operation == "definite_integral":
            variable = row.get("variable")
            if isinstance(variable, int) and not isinstance(variable, bool):
                variable = _bound_name(self.nodes, variable, "integral variable")
            _require(isinstance(variable, str), f"node {index} integral variable is invalid")
            result = frozenset(
                (set(self._free(_child(self.nodes, row, "integrand"))) - {variable})
                | set(self._free(_child(self.nodes, row, "lower")))
                | set(self._free(_child(self.nodes, row, "upper")))
            )
        elif operation == "substitute_original_inverse_angle":
            angle = row.get("angle_variable")
            _require(isinstance(angle, str), f"node {index} inverse angle is invalid")
            result = frozenset(
                (set(self._free(_child(self.nodes, row, "body"))) - {angle})
                | set(self._free(_child(self.nodes, row, "inverse_angle")))
            )
        elif operation == "original_monotone_phase_inverse":
            angle = row.get("angle_variable")
            _require(isinstance(angle, str), f"node {index} inverse angle is invalid")
            result = frozenset(
                (set(self._free(_child(self.nodes, row, "phase_function"))) - {angle})
                | set(self._free(_child(self.nodes, row, "target_phase")))
            )
        elif operation == "original_lazy_flat_branch":
            result = frozenset(
                set(self._free(_child(self.nodes, row, "active_body")))
                | set(self._free(_child(self.nodes, row, "flat_value")))
                | set(self._free(_child(self.nodes, row, "flat_predicate")))
            )
        elif operation == "current_original_unique_analytic_core_component":
            result = self._free(_child(self.nodes, row, "rho")) | self._free(
                _child(self.nodes, row, "Z_variable")
            )
        elif operation in ("current_original_leading_function_recipe", "current_original_axis_pressure_function_derivative"):
            if row.get("coordinate_independent") is not True:
                raise _OpenRule(f"source projection {index} lacks coordinate-independent certificate")
            result = self._free(_child(self.nodes, row, "Z_variable"))
        elif _native_dependent(row):
            raise _OpenRule(f"unknown native-dependent source operation: {operation}")
        else:
            raise _OpenRule(f"unknown source operation: {operation}")
        self.free_cache[index] = result
        return result

    def derivative(self, index: int, context: tuple[str, Any] | None = None) -> int:
        if context is None:
            context = self.main_context
        key = (int(index), context)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        if index == self.phase:
            result = self.index.rational(0)
            self.cache[key] = result
            self.visited[key] = result
            return result
        if context == self.main_context and index in self.seeds:
            result = self.seeds[index]
            self.cache[key] = result
            self.visited[key] = result
            return result
        if key in self.active:
            raise _OpenRule(f"cyclic derivative dependency at source node {index}")
        self.active.add(key)
        try:
            name = self._ctx_name(context)
            if name not in self._free(index):
                result = self.index.rational(0)
                self.cache[key] = result
                self.visited[key] = result
                return result
            row = _row(self.nodes, index)
            operation = _op(self.nodes, index)
            if operation == "bound_variable":
                result = self.index.rational(1) if self._ctx_matches(context, index) else self.index.rational(0)
            elif operation == "sum":
                result = self._add(*(self.derivative(child, context) for child in _args(self.nodes, index)))
            elif operation == "product":
                args = _args(self.nodes, index)
                result = self._add(
                    *(
                        self._mul(
                            *(self.derivative(child, context) if position == selected else child
                              for position, child in enumerate(args))
                        )
                        for selected in range(len(args))
                    )
                )
            elif operation == "negative":
                result = self._neg(self.derivative(_child(self.nodes, row, "argument"), context))
            elif operation == "positive_quotient":
                numerator = _child(self.nodes, row, "numerator")
                denominator = _child(self.nodes, row, "denominator")
                result = self._quotient(
                    self._sub(
                        self._mul(self.derivative(numerator, context), denominator),
                        self._mul(numerator, self.derivative(denominator, context)),
                    ),
                    self._mul(denominator, denominator),
                    row.get("source_positive_certificate"),
                )
            elif operation in ("analytic_unary", "exact_original_source_flat_sigma"):
                argument = _child(self.nodes, row, "argument")
                dx = self.derivative(argument, context)
                name_value = "original_flat_sigma" if operation == "exact_original_source_flat_sigma" else row.get("name")
                if name_value in ("exp", "expm1"):
                    outer = self._analytic("exp", argument)
                elif name_value == "log":
                    outer = self._quotient(
                        self.index.rational(1),
                        argument,
                        "same original positive logarithm input",
                    )
                elif name_value == "log1p":
                    outer = self._quotient(
                        self.index.rational(1),
                        self._add(self.index.rational(1), argument),
                        "same original positive 1+x",
                    )
                elif name_value in ("sqrt", "positive_sqrt"):
                    outer = self._quotient(
                        self.index.rational(1),
                        self._mul(self.index.rational(2), self._analytic("positive_sqrt", argument)),
                        "same original active positive square root",
                    )
                elif name_value == "cos":
                    outer = self._neg(self._analytic("sin", argument))
                elif name_value == "sin":
                    outer = self._analytic("cos", argument)
                elif name_value == "original_flat_sigma":
                    outer = self._sigma_prime(argument)
                elif name_value == "original_flat_sigma_prime":
                    result = self.derivative(self._sigma_prime(argument), context)
                    self.cache[key] = result
                    self.visited[key] = result
                    return result
                else:
                    raise _OpenRule(f"unissued analytic unary: {name_value}")
                result = self._mul(outer, dx)
            elif operation == "function_substitution":
                body = _child(self.nodes, row, "expression")
                variable = _child(self.nodes, row, "variable")
                argument = _child(self.nodes, row, "value")
                variable_context = ("node", variable)
                direct = self.index.rational(0) if self._ctx_matches(context, variable) else self._substitution(
                    self.derivative(body, context), variable, argument
                )
                chain = self._mul(
                    self._substitution(self.derivative(body, variable_context), variable, argument),
                    self.derivative(argument, context),
                )
                result = self._add(direct, chain)
            elif operation == "definite_integral":
                body = _child(self.nodes, row, "integrand")
                lower = _child(self.nodes, row, "lower")
                upper = _child(self.nodes, row, "upper")
                variable = row.get("variable")
                if isinstance(variable, int) and not isinstance(variable, bool):
                    variable_name = _bound_name(self.nodes, variable, "integral variable")
                    variable_handle = variable
                else:
                    _require(isinstance(variable, str), f"node {index} integral variable is invalid")
                    variable_name = variable
                    variable_handle = self.index.bound(variable)
                terms: list[int] = []
                if variable_name != name:
                    partial = self.derivative(body, context)
                    if not self._is_zero(partial):
                        terms.append(self._integral_derivative(row, partial))
                terms.append(self._mul(
                    self._substitution(body, variable_handle, upper),
                    self.derivative(upper, context),
                ))
                terms.append(self._neg(self._mul(
                    self._substitution(body, variable_handle, lower),
                    self.derivative(lower, context),
                )))
                result = self._add(*terms)
            elif operation == "substitute_original_inverse_angle":
                body = _child(self.nodes, row, "body")
                inverse = _child(self.nodes, row, "inverse_angle")
                angle = row.get("angle_variable")
                _require(isinstance(angle, str), f"node {index} inverse angle is invalid")
                direct = self.index.rational(0) if self._ctx_name(context) == angle else self._inverse_at(
                    self.derivative(body, context), inverse, angle
                )
                other = self._inverse_at(self.derivative(body, ("name", angle)), inverse, angle)
                result = self._add(direct, self._mul(other, self.derivative(inverse, context)))
            elif operation == "original_monotone_phase_inverse":
                phase_function = _child(self.nodes, row, "phase_function")
                target = _child(self.nodes, row, "target_phase")
                angle = row.get("angle_variable")
                _require(isinstance(angle, str), f"node {index} inverse angle is invalid")
                inverse = index
                held = self.index.rational(0) if self._ctx_name(context) == angle else self.derivative(
                    phase_function, context
                )
                lam = self._inverse_at(
                    self.derivative(phase_function, ("name", angle)), inverse, angle
                )
                result = self._quotient(
                    self._sub(
                        self.derivative(target, context),
                        self._inverse_at(held, inverse, angle),
                    ),
                    lam,
                    row.get(
                        "source_positive_certificate",
                        "original positive phase lambda; exact implicit native derivative",
                    ),
                )
            elif operation == "original_lazy_flat_branch":
                result = self._lazy_derivative(
                    row,
                    self.derivative(_child(self.nodes, row, "active_body"), context),
                    self.derivative(_child(self.nodes, row, "flat_value"), context),
                )
            elif operation == "current_original_unique_analytic_core_component":
                result = self._core_derivative(row, context)
            elif _native_dependent(row):
                raise _OpenRule(f"unknown native-dependent source operation: {operation}")
            else:
                raise _OpenRule(f"unissued actual native derivative operation: {operation}")
            self.cache[key] = result
            self.visited[key] = result
            return result
        finally:
            self.active.discard(key)


def _normalise_orders(ordinary_orders: Sequence[int]) -> tuple[int, ...]:
    orders = tuple(ordinary_orders)
    _require(bool(orders), "ordinary_orders must not be empty")
    _require(
        all(isinstance(item, int) and not isinstance(item, bool) and item >= 0 for item in orders),
        "ordinary_orders must contain nonnegative integers",
    )
    _require(len(set(orders)) == len(orders) and orders == tuple(sorted(orders)),
             "ordinary_orders must be sorted and unique")
    return orders


def _normalise_side(side: Any) -> str:
    _require(isinstance(side, str) and side.strip(), "side must be a nonempty string")
    value = side.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "q4": "postslow",
        "postrestore": "postslow",
        "post_slow": "postslow",
        "modified": "postslow",
        "finite": "patch",
        "rm_patch": "patch",
        "rm_patch_finite": "patch",
    }
    return aliases.get(value, value)


def _chart_handle(chart: Mapping[str, Any], *keys: str) -> int | None:
    for key in keys:
        value = chart.get(key)
        if value is None:
            continue
        if isinstance(value, Mapping):
            value = value.get("node", value.get("handle", value.get("root")))
        if isinstance(value, int) and not isinstance(value, bool):
            return value
    return None


def _find_jacobian(
    nodes: Sequence[Mapping[str, Any]],
    index: _DescriptorIndex,
    source_chart: Mapping[str, Any],
    loop_chart: Mapping[str, Any],
) -> tuple[int, int]:
    owners = [source_chart, loop_chart]
    # The saved Q4 and finite binding charts expose geometry in ``radial``;
    # the loop rows themselves live under ``loop``.  Search only these two
    # descriptor containers so a similarly named arbitrary payload cannot
    # accidentally supply a Jacobian.
    for owner in tuple(owners):
        for key in ("loop", "radial"):
            nested = owner.get(key)
            if isinstance(nested, Mapping):
                owners.append(nested)
    jacobian = next(
        (
            _chart_handle(owner, "J", "jacobian", "physical_measure", "physical_Jacobian")
            for owner in owners
            if _chart_handle(owner, "J", "jacobian", "physical_measure", "physical_Jacobian") is not None
        ),
        None,
    )
    inverse = next(
        (
            _chart_handle(owner, "invJ", "inv_j", "inverse_jacobian", "physical_measure_inverse")
            for owner in owners
            if _chart_handle(owner, "invJ", "inv_j", "inverse_jacobian", "physical_measure_inverse") is not None
        ),
        None,
    )
    if jacobian is None and inverse is not None:
        candidates = []
        for candidate, row in enumerate(nodes):
            if not isinstance(row, Mapping) or _op(nodes, candidate) != "positive_quotient":
                continue
            if row.get("numerator") == index.rational(1) and row.get("denominator") == inverse:
                candidates.append(candidate)
        if len(candidates) == 1:
            jacobian = candidates[0]
    if inverse is None and jacobian is not None:
        candidates = []
        for candidate, row in enumerate(nodes):
            if not isinstance(row, Mapping) or _op(nodes, candidate) != "positive_quotient":
                continue
            if row.get("numerator") == index.rational(1) and row.get("denominator") == jacobian:
                candidates.append(candidate)
        if len(candidates) == 1:
            inverse = candidates[0]
    if jacobian is None or inverse is None:
        raise ProvenanceError("chart Jacobian and inverse-Jacobian handles are required")
    inverse_row = _row(nodes, inverse, "inverse Jacobian")
    _require(inverse_row.get("operation") == "positive_quotient"
             and inverse_row.get("numerator") == index.rational(1)
             and inverse_row.get("denominator") == jacobian
             and inverse_row.get("source_positive_certificate") == "positive original chart log-radius Jacobian",
             "inverse Jacobian is not the original exact 1/J source")
    return jacobian, inverse


def _source_rows(
    source_chart: Mapping[str, Any],
    radial_background: Mapping[str, Any],
    nodes: Sequence[Mapping[str, Any]],
) -> tuple[dict[tuple[str, int], int], dict[str, Any]]:
    fields = _field_container(source_chart)
    _require(fields is not None, "source_chart leading fields are missing")
    seeds: dict[tuple[str, int], int] = {}
    descriptors: dict[str, Any] = {}
    for field in _FIELDS:
        source_raw = _field_value(fields, field)
        radial_raw = _field_value(radial_background, field)
        if source_raw is None and radial_raw is None:
            continue
        if source_raw is None or radial_raw is None:
            raise ProvenanceError(f"leading field/background pair is incomplete for {field}")
        source = _rows(source_raw, f"source_chart.fields.{field}", len(nodes), (0, 1, 2))
        radial = _rows(radial_raw, f"radial_background.{field}", len(nodes), (0, 1, 2))
        descriptors[field] = {"source": source, "radial": radial}
        # The actual J multiplication is installed by verify after the chart
        # Jacobian has been resolved.  Keep the handles here for diagnostics.
        for order in (0, 1, 2):
            seeds[(field, order)] = source[order]
    _require(bool(descriptors), "source_chart has no usable leading field rows")
    return seeds, descriptors


def _loop_rows(
    loop_chart: Mapping[str, Any],
    nodes: Sequence[Mapping[str, Any]],
    orders: Sequence[int],
) -> tuple[int, int, dict[str, dict[int, int]], dict[str, dict[int, int]]]:
    coordinate = _chart_handle(loop_chart, "coordinate", "coord", "native_coordinate")
    phase = _chart_handle(loop_chart, "actual_phase", "phase", "actual_phase_handle")
    _require(coordinate is not None and phase is not None, "loop_chart coordinate and actual_phase handles are required")
    source: dict[str, dict[int, int]] = {}
    target: dict[str, dict[int, int]] = {}
    for field, aliases in (("A", ("A", "a")), ("B", ("B", "B_over_Pstar", "b"))):
        raw = next((value for key in aliases if (value := _container_value(loop_chart, key)) is not None), None)
        if raw is None:
            raise ProvenanceError(f"loop_chart loop {field} rows are missing")
        source[field] = _rows(raw, f"loop_chart.loop.{field}", len(nodes), orders)
    radial = loop_chart.get("radial")
    if not isinstance(radial, Mapping):
        radial = loop_chart
    for field, aliases in (("A", ("slow_A_y", "slow_A", "A_y")), ("B", ("slow_B_y", "slow_B", "B_y"))):
        raw = next((radial.get(key) for key in aliases if key in radial), None)
        if raw is None:
            raise ProvenanceError(f"loop_chart.radial {field} target rows are missing")
        target[field] = _rows(raw, f"loop_chart.radial.{field}", len(nodes), orders)
    return coordinate, phase, source, target


def _gamma_seeds(
    value: Any,
    nodes: Sequence[Mapping[str, Any]],
) -> dict[int, int]:
    if value in (None, (), []):
        return {}
    if isinstance(value, Mapping):
        value = value.get("rows", value.get("gamma", value.get("handles")))
    rows = _sequence(value, "gamma_seed_rows")
    # Accept the original one-triple form as well as the patch chart's
    # sequence of independent (gamma0, gamma1, gamma2) triples.  Handle
    # validation happens before any bindings are installed; no global-zero
    # assumption is made for any gamma row.
    nested = bool(rows) and isinstance(rows[0], Sequence) and not isinstance(
        rows[0], (str, bytes, bytearray, Mapping)
    )
    triples = rows if nested else [rows]
    result: dict[int, int] = {}
    for triple_offset, triple in enumerate(triples):
        triple_rows = _sequence(triple, f"gamma_seed_rows[{triple_offset}]")
        _require(
            len(triple_rows) == 3,
            "each gamma_seed_rows entry must contain gamma0, gamma1, gamma2",
        )
        handles = [
            _handle(
                item,
                f"gamma_seed_rows[{triple_offset}][{offset}]",
                len(nodes),
            )
            for offset, item in enumerate(triple_rows)
        ]
        raw, first, second = handles
        beta = _row(nodes, raw, "gamma source")
        _require(beta.get("operation") == "compact_raw_beta", "gamma seed source is not original compact beta")
        for ref in (first, second):
            derivative = _row(nodes, ref, "gamma derivative")
            _require(derivative.get("operation") == "original_lazy_flat_branch"
                     and derivative.get("Delta") == beta.get("argument")
                     and derivative.get("flat_value") == 0
                     and derivative.get("original_beta_exact_flat_derivative") is True
                     and derivative.get("active_body_not_evaluated_on_flat_branch") is True,
                     "gamma derivative source contract is not authenticated")
        for source, derivative in ((handles[0], handles[1]), (handles[1], handles[2])):
            old = result.get(source)
            if old is not None and old != derivative:
                raise ProvenanceError(f"conflicting gamma seed for source node {source}")
            result[source] = derivative
    return result


def verify(
    nodes: Sequence[Mapping[str, Any]],
    source_chart: Mapping[str, Any],
    loop_chart: Mapping[str, Any],
    radial_background: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    side: str,
    gamma_seed_rows: Sequence[Any] = (),
    ordinary_orders: Sequence[int] = DEFAULT_ORDERS,
    derivative_source_prefix_count: int | None = None,
    expected_derivative_prefix_digest: str | None = None,
) -> dict[str, Any]:
    """Verify saved slow ``A`` and ``B`` rows against native chain rules.

    Only existing graph descriptors are accepted.  Unsupported reachable
    operations become unresolved row gates, and therefore can never produce a
    ``passed`` result.
    """

    base = _sequence(nodes, "nodes")
    _require(all(isinstance(row, Mapping) for row in base), "nodes must contain mappings")
    _require(isinstance(source_chart, Mapping), "source_chart must be a mapping")
    _require(isinstance(loop_chart, Mapping), "loop_chart must be a mapping")
    _require(isinstance(radial_background, Mapping), "radial_background must be a mapping")
    expected = str(expected_nodes_digest).lower()
    _require(len(expected) == 64, "expected_nodes_digest must be a SHA-256 digest")
    try:
        int(expected, 16)
    except ValueError as exc:
        raise ProvenanceError("expected_nodes_digest must be a SHA-256 digest") from exc
    before = canonical_digest(base)
    _require(before == expected, "nodes digest does not match expected_nodes_digest")
    orders = _normalise_orders(ordinary_orders)
    side_name = _normalise_side(side)
    lookup_count = len(base) if derivative_source_prefix_count is None else derivative_source_prefix_count
    _require(isinstance(lookup_count, int) and not isinstance(lookup_count, bool) and 0 < lookup_count <= len(base),
             "derivative source prefix count is invalid")
    lookup_nodes = base[:lookup_count]
    lookup_digest = canonical_digest(lookup_nodes)
    _require(expected_derivative_prefix_digest is not None and lookup_digest == expected_derivative_prefix_digest,
             "derivative source prefix digest is not authenticated")
    index = _DescriptorIndex(lookup_nodes)
    coordinate, phase, loop_sources, radial_targets = _loop_rows(loop_chart, base, orders)
    jacobian, inverse_jacobian = _find_jacobian(base, index, source_chart, loop_chart)
    leading_handles, leading = _source_rows(source_chart, radial_background, base)
    zero = index.rational(0)
    one = index.rational(1)
    seeds: dict[int, int] = {phase: zero}
    family = source_chart.get("source_family")
    _require(family is not None, "actual source family binding is required")
    _require(canonical_digest(gamma_seed_rows) == canonical_digest(source_chart.get("gamma_seed_rows", ())),
             "gamma seed rows are not bound to the source chart contract")
    seed_oracle = _Oracle(base, index, coordinate=coordinate, phase=phase, seeds={})
    for (field, order), source_handle in leading_handles.items():
        radial_handle = leading[field]["radial"][order]
        derivative_handle = seed_oracle._mul(jacobian, radial_handle)
        _require(source_handle not in (zero, one) or derivative_handle == zero,
                 "a constant leading seed cannot have nonzero native derivative")
        old = seeds.get(source_handle)
        if old is not None and old != derivative_handle:
            raise ProvenanceError(
                f"conflicting leading-field seed for source node {source_handle}"
            )
        seeds[source_handle] = derivative_handle
    for source_handle, derivative_handle in _gamma_seeds(gamma_seed_rows, base).items():
        _require(source_handle not in (zero, one) or derivative_handle == zero,
                 "a constant gamma seed cannot have nonzero native derivative")
        old = seeds.get(source_handle)
        if old is not None and old != derivative_handle:
            raise ProvenanceError(f"conflicting seed for source node {source_handle}")
        seeds[source_handle] = derivative_handle
    oracle = _Oracle(base, index, coordinate=coordinate, phase=phase, seeds=seeds,
                     expected_source_family=family)
    rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for field in ("A", "B"):
        for order in orders:
            source_handle = loop_sources[field][order]
            target_handle = radial_targets[field][order]
            row: dict[str, Any] = {
                "side": side_name,
                "field": field,
                "ordinary_order": order,
                "source_handle": source_handle,
                "source_node": source_handle,
                "physical_target_handle": target_handle,
                "target_handle": target_handle,
                "status": "unresolved",
                "formula": "physical_target = invJ * D_native(source) at held actual phase",
            }
            try:
                derivative_handle = oracle.derivative(source_handle)
                expected_target = oracle._mul(inverse_jacobian, derivative_handle)
                row.update(
                    derivative_handle=derivative_handle,
                    expected_physical_target_handle=expected_target,
                    derivative_root=derivative_handle,
                    physical_target_matches=target_handle == expected_target,
                )
                if target_handle != expected_target:
                    gate = {
                        "gate": "physical_slow_target_descriptor_mismatch",
                        "side": side_name,
                        "field": field,
                        "ordinary_order": order,
                        "source_handle": source_handle,
                        "physical_target_handle": target_handle,
                        "expected_physical_target_handle": expected_target,
                    }
                    unresolved.append(gate)
                    row["error"] = "saved physical slow target is not invJ*D_native(source)"
                else:
                    row["status"] = "proved"
            except (ProvenanceError, KeyError, TypeError, ValueError) as exc:
                row["error"] = str(exc)
                unresolved.append(
                    {
                        "gate": "native_derivative_rule_unresolved",
                        "side": side_name,
                        "field": field,
                        "ordinary_order": order,
                        "source_handle": source_handle,
                        "error": str(exc),
                    }
                )
            rows.append(row)
    after = canonical_digest(base)
    if after != before:
        unresolved.append({"gate": "nodes_input_mutated", "before": before, "after": after})
    if index.errors:
        # Unreachable malformed descriptors are irrelevant to this proof.  Do
        # not turn unrelated trailing graph data into a false failure.
        pass
    visited_pairs = [
        {
            "source_handle": source,
            "derivative_context": list(context),
            "derivative_handle": derivative,
        }
        for (source, context), derivative in sorted(oracle.visited.items(), key=lambda item: (item[0][0], repr(item[0][1])))
    ]
    visited_map: dict[str, list[dict[str, Any]]] = {}
    for pair in visited_pairs:
        visited_map.setdefault(str(pair["source_handle"]), []).append(pair)
    passed = bool(rows) and not unresolved and all(row.get("status") == "proved" for row in rows)
    return {
        "schema": "original_Rm_native_derivative_provenance_v1",
        "passed": passed,
        "status": "proved" if passed else "unresolved",
        "source_only": True,
        "source_executed": False,
        "actual_saved_graph_run": False,
        "synthetic_execution": False,
        "nodes_input_unmutated": after == before,
        "expected_nodes_digest": expected,
        "nodes_digest": after,
        "derivative_source_prefix_count": lookup_count,
        "derivative_source_prefix_digest": lookup_digest,
        "side": side_name,
        "ordinary_orders": list(orders),
        "coordinate_handle": coordinate,
        "actual_phase_handle": phase,
        "physical_J_handle": jacobian,
        "inverse_J_handle": inverse_jacobian,
        "seed_handles": {str(source): derivative for source, derivative in sorted(seeds.items())},
        "rows": rows,
        "visited_source_D_pairs": visited_pairs,
        "visited_source_to_D_pairs": visited_map,
        "unresolved_gates": unresolved,
        "closure_claim": False,
        "numeric_values_evaluated": False,
    }


__all__ = ["DEFAULT_ORDERS", "ProvenanceError", "canonical_digest", "verify"]
