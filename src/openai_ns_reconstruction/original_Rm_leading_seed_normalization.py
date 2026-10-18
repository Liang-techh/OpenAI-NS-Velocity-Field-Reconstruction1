"""Bounded exact source reducer for the actual leading ``R_m`` seed.

The reducer consumes an already materialised source graph and small binding
maps.  It deliberately knows nothing about graph producers, numerical
providers, or archived checkpoints.  Endpoint source functions are retained
as qualified atoms (with their ordinary ``Z`` rows) while the actual ``E`` and
``V`` normalisations are reduced through the source arithmetic.  Empty
equal-bound integrals and the compact flat branch are discharged before an
integrand is inspected.

The public entry point is :func:`prove`.  A proof packet can therefore be
assembled by a parent runner without importing any of the source owners.
Unsupported reachable operations produce explicit unresolved gates; they are
never assigned a desired value merely to make a row pass.
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
_PATCH_X1_SCOPE = "patch_endpoint_x1"
_INLET_X1_SCOPE = "inlet_endpoint_x1"
_MAX_EXACT_INTEGER_POWER = 16

# This is intentionally a closed list.  The two quotient spellings and the
# exact integer/symbol spellings occur in the old source graph; all are
# reduced by the same bounded implementation below.
_DECLARED_OPERATIONS = frozenset(
    {
        "exact_rational",
        "constant",
        "symbol",
        "atom",
        "source_atom",
        "qualified_source_atom",
        "original_source_parameter",
        "current_original_source_parameter",
        "sum",
        "negative",
        "product",
        "positive_quotient",
        "positive_function_quotient",
        "quotient",
        "function_substitution",
        "definite_integral",
        "integral",
        "bound_variable",
        "analytic_unary",
        "compact_raw_beta",
        "compact_raw_beta_derivative",
        "compact_raw_beta_second_derivative",
        "original_lazy_flat_branch",
        "logical_or",
        "logical_and",
        "logical_not",
        "exact_real_comparison",
        "trim",
        "exact_positive_integer_expression",
        "original_Rm_unique_leading_implicit_control",
    }
)


class _UnsupportedOperation(ValueError):
    """An operation outside the bounded reducer's declared closure."""


def canonical_digest(value: Any) -> str:
    """Return the canonical JSON SHA-256 digest of a JSON value."""

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


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _handle(value: Any, path: str, node_count: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{path} must be an integer graph handle")
    if value < 0:
        raise ValueError(f"{path} must be nonnegative")
    if node_count is not None and value >= node_count:
        raise ValueError(f"{path} handle {value} is outside nodes")
    return value


def _fraction(value: Any, path: str = "rational") -> Fraction:
    """Parse exact integer or rational literals without float conversion."""

    if isinstance(value, bool):
        raise ValueError(f"{path} must be an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        text = value.strip()
        if "/" in text:
            pieces = text.split("/")
            if len(pieces) != 2:
                raise ValueError(f"{path} has malformed rational literal")
            return Fraction(int(pieces[0]), int(pieces[1]))
        return Fraction(int(text), 1)
    raise ValueError(f"{path} must be an exact integer or rational literal")


def _fraction_text(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _row_operation(row: Mapping[str, Any], index: int) -> str:
    operation = row.get("operation", row.get("op"))
    if not isinstance(operation, str):
        raise ValueError(f"node {index} has no operation")
    return operation


def _arg_list(row: Mapping[str, Any], index: int) -> list[Any]:
    args = row.get("arguments", row.get("operands"))
    if not _is_sequence(args) or not args:
        raise ValueError(f"node {index}.arguments must be a nonempty sequence")
    return list(args)


def _field(row: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    raise KeyError(names[0])


def _factor_sort_key(value: tuple[Any, ...]) -> str:
    return repr(value)


def _factor_map(monomial: tuple[tuple[Any, ...], ...]) -> dict[tuple[Any, ...], int]:
    result: dict[tuple[Any, ...], int] = {}
    for factor in monomial:
        if not isinstance(factor, tuple) or len(factor) < 2:
            raise ValueError("malformed symbolic factor")
        base = factor[:-1]
        exponent = factor[-1]
        if isinstance(exponent, bool) or not isinstance(exponent, int):
            raise ValueError("symbolic factor exponent is not integral")
        result[base] = result.get(base, 0) + exponent
    return {key: exponent for key, exponent in result.items() if exponent}


def _monomial(mapping: Mapping[tuple[Any, ...], int]) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        tuple(base) + (exponent,)
        for base, exponent in sorted(mapping.items(), key=lambda item: _factor_sort_key(item[0]))
        if exponent
    )


def _remove_inverse_factor(expression: "_Expr", denominator_key: tuple[Any, ...]) -> "_Expr" | None:
    """Remove one matching polynomial inverse from every monomial."""

    target = ("inverse_expression", denominator_key)
    terms: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
    for monomial, coefficient in expression.terms.items():
        factors = _factor_map(monomial)
        if factors.get(target, 0) != 1:
            return None
        factors.pop(target)
        key = _monomial(factors)
        terms[key] = terms.get(key, Fraction(0)) + coefficient
    return _Expr(terms)


@dataclass(frozen=True)
class _Expr:
    """A finite exact rational Laurent expression in qualified factors."""

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
    def factor(cls, factor: tuple[Any, ...], exponent: int = 1) -> "_Expr":
        return cls({_monomial({tuple(factor): exponent}): Fraction(1)})

    @classmethod
    def atom(cls, name: str) -> "_Expr":
        return cls.factor(("atom", str(name)))

    @property
    def is_zero(self) -> bool:
        return not self.terms

    @property
    def is_rational(self) -> bool:
        return not self.terms or set(self.terms) <= {()}

    @property
    def rational_value(self) -> Fraction:
        if not self.is_rational:
            raise ValueError("expression is not rational")
        return self.terms.get((), Fraction(0))

    def scale(self, value: Fraction | int) -> "_Expr":
        value = Fraction(value)
        return _Expr({monomial: coefficient * value for monomial, coefficient in self.terms.items()})

    def __add__(self, other: "_Expr") -> "_Expr":
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            terms[monomial] = terms.get(monomial, Fraction(0)) + coefficient
        return _Expr(terms)

    def __neg__(self) -> "_Expr":
        return self.scale(-1)

    def __sub__(self, other: "_Expr") -> "_Expr":
        return self + (-other)

    def power_integer(self, exponent: int) -> "_Expr":
        """Evaluate a bounded exact integer power without numeric conversion."""

        if isinstance(exponent, bool) or not isinstance(exponent, int):
            raise _UnsupportedOperation("analytic integer power is not integral")
        if abs(exponent) > _MAX_EXACT_INTEGER_POWER:
            raise _UnsupportedOperation(
                f"analytic integer power {exponent} exceeds the exact reducer bound"
            )
        if exponent == 0:
            return _Expr.rational(1)
        base = self if exponent > 0 else self.inverse()
        remaining = abs(exponent)
        result = _Expr.rational(1)
        while remaining:
            if remaining & 1:
                result = result * base
            remaining >>= 1
            if remaining:
                base = base * base
        return result

    def __mul__(self, other: "_Expr") -> "_Expr":
        if self.is_zero or other.is_zero:
            return _Expr.rational(0)
        # Quotients retain a polynomial denominator as an inverse-expression
        # factor.  Cancel it against the exact denominator expression before
        # distributing terms.  This handles ordinary derivative rows such as
        # ``(6 Z^2-2)c/(1+Z^2)^3`` without expanding a symbolic denominator.
        if self.key() == other._inverse_denominator_key():
            removed = _remove_inverse_factor(other, self.key())
            if removed is not None:
                return removed
        if other.key() == self._inverse_denominator_key():
            removed = _remove_inverse_factor(self, other.key())
            if removed is not None:
                return removed
        # A quotient of a polynomial by itself may carry an additional
        # source factor, e.g. ``(A * H) * A^{-1}``.  The direct factor
        # cancellation above handles ``A * A^{-1}``, while this bounded exact
        # division handles the same common polynomial with a qualified source
        # atom attached.  It never guesses a value: cancellation is returned
        # only after an exact Laurent-polynomial divisibility check.
        cancelled = _cancel_inverse_operand(self, other)
        if cancelled is not None:
            return cancelled
        cancelled = _cancel_inverse_operand(other, self)
        if cancelled is not None:
            return cancelled
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

    def _inverse_denominator_key(self) -> tuple[Any, ...] | None:
        """Return a single polynomial denominator key when present."""

        key: tuple[Any, ...] | None = None
        for monomial in self.terms:
            factors = _factor_map(monomial)
            candidates = [factor[1] for factor, exponent in factors.items() if factor[0] == "inverse_expression" and exponent]
            if len(candidates) != 1:
                return None
            candidate = candidates[0]
            if key is None:
                key = candidate
            elif key != candidate:
                return None
        return key

    def inverse(self) -> "_Expr":
        if self.is_zero:
            raise ValueError("zero expression cannot be inverted")
        if len(self.terms) == 1:
            monomial, coefficient = next(iter(self.terms.items()))
            mapping = {factor: -exponent for factor, exponent in _factor_map(monomial).items()}
            return _Expr({_monomial(mapping): Fraction(1, 1) / coefficient})
        # A polynomial denominator is retained as one qualified factor.  A
        # later exp(-log(denominator)) creates this exact same factor.
        return _Expr.factor(("inverse_expression", self.key()))

    def key(self) -> tuple[Any, ...]:
        return tuple(
            (
                tuple(tuple(factor) for factor in monomial),
                coefficient.numerator,
                coefficient.denominator,
            )
            for monomial, coefficient in sorted(self.terms.items(), key=lambda item: repr(item[0]))
        )

    def json(self) -> dict[str, Any]:
        return {
            "terms": [
                {
                    "monomial": [list(factor) for factor in monomial],
                    "coefficient": _fraction_text(coefficient),
                }
                for monomial, coefficient in sorted(self.terms.items(), key=lambda item: repr(item[0]))
            ]
        }


def _expression_from_key(key: tuple[Any, ...]) -> _Expr | None:
    """Reconstruct the finite expression encoded by ``_Expr.key``."""

    terms: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
    try:
        for item in key:
            monomial_key, numerator, denominator = item
            monomial = tuple(tuple(factor) for factor in monomial_key)
            terms[monomial] = Fraction(int(numerator), int(denominator))
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    return _Expr(terms)


def _monomial_vector(
    monomial: tuple[tuple[Any, ...], ...],
    factor_order: tuple[tuple[Any, ...], ...],
) -> tuple[int, ...]:
    factors = _factor_map(monomial)
    return tuple(factors.get(base, 0) for base in factor_order)


def _divide_exact(numerator: _Expr, denominator: _Expr) -> _Expr | None:
    """Return an exact Laurent quotient when finite polynomial division works."""

    if denominator.is_zero:
        return None
    if numerator.is_zero:
        return _Expr.rational(0)
    denominator_terms = dict(denominator.terms)
    numerator_terms = dict(numerator.terms)
    if len(denominator_terms) == 1:
        denominator_monomial, denominator_coefficient = next(iter(denominator_terms.items()))
        denominator_factors = _factor_map(denominator_monomial)
        result: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
        for monomial, coefficient in numerator_terms.items():
            factors = _factor_map(monomial)
            for factor, exponent in denominator_factors.items():
                factors[factor] = factors.get(factor, 0) - exponent
            quotient_monomial = _monomial(factors)
            result[quotient_monomial] = result.get(quotient_monomial, Fraction(0)) + coefficient / denominator_coefficient
        return _Expr(result)

    factor_bases: set[tuple[Any, ...]] = set()
    for expression in (numerator, denominator):
        for monomial in expression.terms:
            factor_bases.update(_factor_map(monomial))
    factor_order = tuple(sorted(factor_bases, key=repr))
    denominator_vectors = {
        _monomial_vector(monomial, factor_order): coefficient
        for monomial, coefficient in denominator_terms.items()
    }
    remainder = {
        _monomial_vector(monomial, factor_order): coefficient
        for monomial, coefficient in numerator_terms.items()
    }
    pivot_vector, pivot_coefficient = min(denominator_vectors.items())
    quotient: dict[tuple[int, ...], Fraction] = {}
    # A divisible finite Laurent polynomial removes one leading term per
    # quotient term.  The cap keeps unsupported/nondivisible source algebra
    # fail-closed instead of allowing an unbounded symbolic loop.
    iteration_cap = max(64, (len(numerator_terms) + 1) * (len(denominator_terms) + 1) * 8)
    iterations = 0
    while remainder:
        if iterations >= iteration_cap:
            return None
        leading_vector, leading_coefficient = min(remainder.items())
        quotient_vector = tuple(left - right for left, right in zip(leading_vector, pivot_vector))
        quotient_coefficient = leading_coefficient / pivot_coefficient
        quotient[quotient_vector] = quotient.get(quotient_vector, Fraction(0)) + quotient_coefficient
        if not quotient[quotient_vector]:
            quotient.pop(quotient_vector)
        for denominator_vector, denominator_coefficient in denominator_vectors.items():
            target = tuple(left + right for left, right in zip(quotient_vector, denominator_vector))
            value = remainder.get(target, Fraction(0)) - quotient_coefficient * denominator_coefficient
            if value:
                remainder[target] = value
            else:
                remainder.pop(target, None)
        iterations += 1
    terms: dict[tuple[tuple[Any, ...], ...], Fraction] = {}
    for vector, coefficient in quotient.items():
        factors = {
            base: exponent
            for base, exponent in zip(factor_order, vector)
            if exponent
        }
        terms[_monomial(factors)] = coefficient
    candidate = _Expr(terms)
    # Verify using the ordinary exact product; this guards the leading-term
    # division against accidental cancellation in a nondivisible expression.
    if denominator * candidate != numerator:
        return None
    return candidate


def _cancel_inverse_operand(numerator: _Expr, inverse_operand: _Expr) -> _Expr | None:
    """Cancel one polynomial-inverse factor when the numerator is divisible."""

    if numerator.is_zero or len(inverse_operand.terms) != 1:
        return None
    inverse_monomial, _ = next(iter(inverse_operand.terms.items()))
    factors = _factor_map(inverse_monomial)
    inverse_factors = [
        (factor, exponent)
        for factor, exponent in factors.items()
        if factor[0] == "inverse_expression" and exponent == 1
    ]
    if len(inverse_factors) != 1:
        return None
    inverse_factor, _ = inverse_factors[0]
    denominator_key = inverse_factor[1]
    denominator = _expression_from_key(denominator_key)
    if denominator is None:
        return None
    quotient = _divide_exact(numerator, denominator)
    if quotient is None:
        return None
    remaining = _remove_inverse_factor(inverse_operand, denominator_key)
    if remaining is None:
        return None
    return quotient * remaining


@dataclass(frozen=True)
class _EnvValue:
    value_handle: int | None
    expression: _Expr


@dataclass(frozen=True)
class _Alias:
    atom: str
    source_node: int
    descriptor: Mapping[str, Any]
    kind: str


def _path_value(value: Mapping[str, Any], paths: Sequence[Sequence[str]], label: str) -> Any:
    for path in paths:
        current: Any = value
        for key in path:
            if not isinstance(current, Mapping) or key not in current:
                break
            current = current[key]
        else:
            return current
    raise ValueError(f"{label} binding is missing")


def _first_mapping(value: Any, names: Sequence[str]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("bindings must be a mapping")
    for name in names:
        nested = value.get(name)
        if isinstance(nested, Mapping):
            return nested
    return value


def _rows(value: Any, label: str, ordinary_orders: Sequence[int]) -> list[int]:
    if isinstance(value, Mapping):
        result: list[int] = []
        for order in ordinary_orders:
            candidate = value.get(str(order), value.get(order))
            result.append(_handle(candidate, f"{label}[{order}]"))
        return result
    if not _is_sequence(value):
        raise ValueError(f"{label} must be a sequence or order mapping")
    result = []
    for order in ordinary_orders:
        if order < 0 or order >= len(value):
            raise ValueError(f"{label} lacks ordinary order {order}")
        result.append(_handle(value[order], f"{label}[{order}]"))
    return result


def _source_label(field_name: str, order: int) -> str:
    return f"source:{field_name}[{order}]"


def _canonical_descriptor(row: Mapping[str, Any], index: int) -> dict[str, Any]:
    operation = _row_operation(row, index)
    if operation != "function_substitution":
        return {"node": index, "operation": operation}
    descriptor: dict[str, Any] = {"operation": operation}
    for key in ("expression", "variable", "value", _Z_GUARD):
        if key in row:
            descriptor[key] = row[key]
    descriptor[_Z_GUARD] = row.get(_Z_GUARD) is True
    return descriptor


class _Reducer:
    """Conservative exact evaluator for reachable selected source rows."""

    def __init__(
        self,
        nodes: Sequence[Mapping[str, Any]],
        *,
        protected: Mapping[int, _Alias],
        endpoint_aliases: Mapping[int, _Alias],
        ordinary_orders: Sequence[int],
        body_aliases: Mapping[int, _Alias] = (),
        zero_aliases: Sequence[int] = (),
        zero_alias_scope_guard: bool = False,
    ) -> None:
        if not _is_sequence(nodes):
            raise ValueError("nodes must be a sequence")
        self.nodes = nodes
        self.protected = dict(protected)
        self.endpoint_aliases = dict(endpoint_aliases)
        self.body_aliases = dict(body_aliases)
        self.zero_aliases = frozenset(int(ref) for ref in zero_aliases)
        self.zero_alias_scope_guard = bool(zero_alias_scope_guard)
        self._active_scope: str | None = None
        self.descriptor_aliases: dict[str, _Alias] = {}
        for alias in self.protected.values():
            if alias.descriptor.get("operation") != "function_substitution":
                continue
            key = canonical_digest(dict(alias.descriptor))
            previous = self.descriptor_aliases.get(key)
            if previous is not None and previous.atom != alias.atom:
                raise ValueError("distinct protected source ports share one substitution descriptor")
            self.descriptor_aliases[key] = alias
        self.ordinary_orders = tuple(ordinary_orders)
        self.memo: dict[tuple[str | None, int, tuple[Any, ...]], _Expr] = {}
        self.used_nodes: set[int] = set()
        self.equal_bound_integrals: set[int] = set()
        self.flat_beta_nodes: set[int] = set()
        self.alias_witnesses: list[dict[str, Any]] = []
        self._witness_keys: set[str] = set()
        self._log_bases: dict[tuple[Any, ...], _Expr] = {}
        self._exp_arguments: dict[tuple[Any, ...], _Expr] = {}

    def _row(self, index: int) -> Mapping[str, Any]:
        _handle(index, "source handle", len(self.nodes))
        row = self.nodes[index]
        if not isinstance(row, Mapping):
            raise ValueError(f"nodes[{index}] must be a mapping")
        return row

    def _operation(self, index: int, row: Mapping[str, Any]) -> str:
        operation = _row_operation(row, index)
        if operation not in _DECLARED_OPERATIONS:
            raise _UnsupportedOperation(
                f"visited node {index} uses unsupported operation {operation!r}"
            )
        return operation

    def _operand_handles(self, index: int, row: Mapping[str, Any]) -> list[tuple[str, int]]:
        operation = self._operation(index, row)
        values: list[tuple[str, Any]] = []
        if operation in {"sum", "product", "logical_or", "logical_and"}:
            values.extend((f"arguments[{i}]", item) for i, item in enumerate(_arg_list(row, index)))
        elif operation in {"negative", "trim", "exact_positive_integer_expression"}:
            values.append(("argument", _field(row, "argument", "operand")))
        elif operation in {"positive_quotient", "positive_function_quotient", "quotient"}:
            values.extend(
                (("numerator", _field(row, "numerator")), ("denominator", _field(row, "denominator")))
            )
        elif operation in {"analytic_unary", "compact_raw_beta", "compact_raw_beta_derivative", "compact_raw_beta_second_derivative"}:
            values.append(("argument", _field(row, "argument", "operand")))
        elif operation == "function_substitution":
            values.append(("expression", _field(row, "expression")))
            variable = _field(row, "variable")
            if isinstance(variable, int) and not isinstance(variable, bool):
                values.append(("variable", variable))
            elif not isinstance(variable, str):
                raise ValueError(f"node {index}.variable must be a handle or string")
            values.append(("value", _field(row, "value")))
        elif operation in {"definite_integral", "integral"}:
            values.extend(
                (
                    ("integrand", _field(row, "integrand")),
                    ("lower", _field(row, "lower")),
                    ("upper", _field(row, "upper")),
                )
            )
            variable = _field(row, "variable", "angle_variable")
            if isinstance(variable, int) and not isinstance(variable, bool):
                values.append(("variable", variable))
            elif not isinstance(variable, str):
                raise ValueError(f"node {index}.variable must be a handle or string")
        elif operation == "original_lazy_flat_branch":
            values.extend(
                (
                    ("active_body", _field(row, "active_body")),
                    ("flat_value", _field(row, "flat_value")),
                    ("flat_predicate", _field(row, "flat_predicate")),
                )
            )
        elif operation == "exact_real_comparison":
            values.extend(
                (("left", _field(row, "left")), ("right", _field(row, "right")))
            )
        elif operation == "logical_not":
            values.append(("argument", _field(row, "argument", "operand")))
        elif operation in {
            "exact_rational",
            "constant",
            "symbol",
            "atom",
            "source_atom",
            "qualified_source_atom",
            "original_source_parameter",
            "current_original_source_parameter",
            "bound_variable",
        }:
            pass
        elif operation == "original_Rm_unique_leading_implicit_control":
            # This node is a source atom unless it exposes an explicit
            # arithmetic child; metadata is intentionally left intact.
            child = row.get("value", row.get("argument"))
            if child is not None:
                values.append(("value", child))
        result: list[tuple[str, int]] = []
        for path, value in values:
            if isinstance(value, bool) or not isinstance(value, int):
                # Variable names are allowed for binders and do not need
                # graph validation.
                if path in {"variable"} and isinstance(value, str):
                    continue
                raise ValueError(f"node {index}.{path} must be a graph handle")
            result.append((path, _handle(value, f"node {index}.{path}", len(self.nodes))))
        return result

    def _validate_operands(self, index: int, row: Mapping[str, Any]) -> None:
        for path, child in self._operand_handles(index, row):
            if child >= index:
                raise ValueError(f"forward reference from node {index} {path} to {child}")

    @staticmethod
    def _env_signature(env: Mapping[Any, _EnvValue]) -> tuple[Any, ...]:
        return tuple(
            (repr(key), value.value_handle, value.expression.key())
            for key, value in sorted(env.items(), key=lambda item: repr(item[0]))
        )

    def _record_witness(
        self,
        alias: _Alias,
        *,
        via: str,
        expression_node: int,
        environment: Mapping[Any, _EnvValue],
    ) -> None:
        payload = {
            "atom": alias.atom,
            "kind": alias.kind,
            "source_node": alias.source_node,
            "expression_node": expression_node,
            "descriptor": deepcopy(dict(alias.descriptor)),
            "via": via,
            "scoped_environment": {
                str(key): {
                    "value_handle": value.value_handle,
                    "value": value.expression.json(),
                }
                for key, value in sorted(environment.items(), key=lambda item: repr(item[0]))
            },
        }
        key = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        if key not in self._witness_keys:
            self._witness_keys.add(key)
            self.alias_witnesses.append(payload)

    def _scoped_alias(self, ref: int, env: Mapping[Any, _EnvValue]) -> _Alias | None:
        for variable, binding in sorted(env.items(), key=lambda item: repr(item[0])):
            if binding.value_handle is None:
                continue
            alias = self.endpoint_aliases.get(ref)
            if alias is None:
                continue
            descriptor = alias.descriptor
            if descriptor.get("variable") == variable and descriptor.get("value") == binding.value_handle:
                self._record_witness(alias, via="scoped_function_substitution", expression_node=ref, environment=env)
                return alias
        return None

    def _rewrite_bound(self, ref: int, env: Mapping[Any, _EnvValue]) -> _Expr:
        if ref in env:
            return env[ref].expression
        row = self._row(ref)
        name = row.get("name", row.get("symbol", f"bound[{ref}]"))
        if name in env:
            return env[name].expression
        # A free symbol is a source function of Z.  It remains an atom; this
        # is preferable to guessing a polynomial body.
        return _Expr.atom(f"symbol:{name}")

    @staticmethod
    def _linear_factors(expression: _Expr) -> tuple[Fraction, list[tuple[tuple[Any, ...], Fraction]]] | None:
        constant = expression.terms.get((), Fraction(0))
        factors: list[tuple[tuple[Any, ...], Fraction]] = []
        for monomial, coefficient in expression.terms.items():
            if monomial == ():
                continue
            if len(monomial) != 1:
                return None
            factor = monomial[0][:-1]
            exponent = monomial[0][-1]
            if exponent != 1:
                return None
            factors.append((factor, coefficient))
        return constant, factors

    def _analytic(self, index: int, row: Mapping[str, Any], argument: _Expr) -> _Expr:
        name = row.get("name", row.get("function", row.get("unary")))
        if not isinstance(name, str):
            raise ValueError(f"node {index}.analytic_unary has no name")
        name = name.lower()
        if name == "exp":
            linear = self._linear_factors(argument)
            if linear is not None:
                constant, factors = linear
                result = _Expr.rational(1)
                if constant:
                    result = result * _Expr.factor(("exp_rational", constant.numerator, constant.denominator))
                for factor, coefficient in factors:
                    if factor and factor[0] == "log" and len(factor) == 2:
                        base = self._log_bases.get(factor)
                        if base is not None and coefficient.denominator == 1:
                            result = result * base.power_integer(coefficient.numerator)
                            continue
                    if coefficient.denominator == 1:
                        result = result * _Expr.factor(("exp_factor", factor, coefficient.numerator))
                    else:
                        result = result * _Expr.factor(("exp_linear_factor", factor, coefficient.numerator, coefficient.denominator))
                self._exp_arguments[result.key()] = argument
                return result
            opaque = _Expr.factor(("exp_expression", argument.key()))
            self._exp_arguments[opaque.key()] = argument
            return opaque
        if name == "log":
            if argument.is_zero:
                raise ValueError("logarithm of exact zero")
            if argument.is_rational:
                if argument.rational_value <= 0:
                    raise ValueError("logarithm requires a positive exact rational")
                if argument.rational_value == 1:
                    return _Expr.rational(0)
                factor = ("log_rational", argument.rational_value.numerator, argument.rational_value.denominator)
                return _Expr.factor(factor)
            if len(argument.terms) == 1:
                monomial, coefficient = next(iter(argument.terms.items()))
                if coefficient == 1 and len(monomial) == 1:
                    factor = monomial[0][:-1]
                    exponent = monomial[0][-1]
                    if exponent == 1 and factor and factor[0] == "exp_factor":
                        # The exponent factor is intentionally kept opaque.
                        # It still cancels by factor identity when a source
                        # graph contains the same exp/log wrapper on both
                        # sides; no inverse logarithm law is needed here.
                        original = factor[1]
                        return _Expr.factor(("log_exp_factor", original))
            factor = ("log", argument.key())
            self._log_bases[factor] = argument
            return _Expr.factor(factor)
        if name == "sqrt":
            return _Expr.factor(("sqrt_expression", argument.key()))
        raise _UnsupportedOperation(f"node {index} analytic function {name!r} is unsupported")

    def _comparison(self, index: int, row: Mapping[str, Any], left: _Expr, right: _Expr) -> bool:
        if not left.is_rational or not right.is_rational:
            raise _UnsupportedOperation(f"node {index} comparison is not exact")
        operator = row.get("operator", row.get("relation"))
        if operator == "<=":
            return left.rational_value <= right.rational_value
        if operator == "<":
            return left.rational_value < right.rational_value
        if operator == ">=":
            return left.rational_value >= right.rational_value
        if operator == ">":
            return left.rational_value > right.rational_value
        if operator in {"==", "="}:
            return left.rational_value == right.rational_value
        if operator == "!=":
            return left.rational_value != right.rational_value
        raise _UnsupportedOperation(f"node {index} comparison operator {operator!r} is unsupported")

    def _boolean(self, ref: int, env: Mapping[Any, _EnvValue]) -> bool:
        row = self._row(ref)
        operation = self._operation(ref, row)
        self._validate_operands(ref, row)
        if operation == "exact_real_comparison":
            return self._comparison(ref, row, self.rewrite(row["left"], env), self.rewrite(row["right"], env))
        if operation == "logical_or":
            return any(self._boolean(child, env) for child in _arg_list(row, ref))
        if operation == "logical_and":
            return all(self._boolean(child, env) for child in _arg_list(row, ref))
        if operation == "logical_not":
            return not self._boolean(_field(row, "argument", "operand"), env)
        value = self.rewrite(ref, env)
        if not value.is_rational:
            raise _UnsupportedOperation(f"node {ref} is not an exact boolean")
        return bool(value.rational_value)

    def rewrite(
        self,
        ref: int,
        env: Mapping[Any, _EnvValue] | None = None,
        *,
        scope: str | None = None,
    ) -> _Expr:
        """Rewrite a source node under an explicit endpoint scope."""

        if scope is None:
            return self._rewrite(ref, env)
        previous_scope = self._active_scope
        self._active_scope = scope
        try:
            return self._rewrite(ref, env)
        finally:
            self._active_scope = previous_scope

    def _rewrite(self, ref: int, env: Mapping[Any, _EnvValue] | None = None) -> _Expr:
        env = {} if env is None else env
        _handle(ref, "source handle", len(self.nodes))
        key = (self._active_scope, ref, self._env_signature(env))
        if key in self.memo:
            return self.memo[key]
        if ref in env:
            result = env[ref].expression
            self.memo[key] = result
            return result
        row_for_binding = self._row(ref)
        row_name = row_for_binding.get("name", row_for_binding.get("symbol"))
        if row_name in env:
            result = env[row_name].expression
            self.memo[key] = result
            return result
        alias = self.protected.get(ref)
        if alias is not None:
            self.used_nodes.add(ref)
            self._record_witness(alias, via="protected_source_port", expression_node=ref, environment=env)
            result = _Expr.atom(alias.atom)
            self.memo[key] = result
            return result
        body_alias = self.body_aliases.get(ref)
        if body_alias is not None:
            self.used_nodes.add(ref)
            self._record_witness(
                body_alias,
                via="guarded_common_source_body",
                expression_node=ref,
                environment=env,
            )
            result = _Expr.atom(body_alias.atom)
            self.memo[key] = result
            return result
        if self.zero_alias_scope_guard and self._active_scope == _PATCH_X1_SCOPE and ref in self.zero_aliases:
            self.used_nodes.add(ref)
            result = _Expr.rational(0)
            self.memo[key] = result
            return result
        row_for_alias = self._row(ref)
        if _row_operation(row_for_alias, ref) == "function_substitution":
            descriptor_alias = self.descriptor_aliases.get(canonical_digest(_canonical_descriptor(row_for_alias, ref)))
            if descriptor_alias is not None:
                self.used_nodes.add(ref)
                self._record_witness(
                    descriptor_alias,
                    via="matching_source_substitution_descriptor",
                    expression_node=ref,
                    environment=env,
                )
                result = _Expr.atom(descriptor_alias.atom)
                self.memo[key] = result
                return result
        scoped = self._scoped_alias(ref, env)
        if scoped is not None:
            self.used_nodes.add(ref)
            result = _Expr.atom(scoped.atom)
            self.memo[key] = result
            return result

        row = self._row(ref)
        self.used_nodes.add(ref)
        operation = self._operation(ref, row)
        self._validate_operands(ref, row)
        if operation == "exact_rational":
            numerator = row.get("numerator", row.get("value", 0))
            denominator = row.get("denominator", 1)
            result = _Expr.rational(Fraction(_fraction(numerator, f"node {ref}.numerator"), _fraction(denominator, f"node {ref}.denominator")))
        elif operation == "constant":
            result = _Expr.rational(_fraction(row.get("value"), f"node {ref}.value"))
        elif operation in {"symbol", "atom", "source_atom", "qualified_source_atom", "original_source_parameter", "current_original_source_parameter"}:
            if row.get("exact_zero") is True or row.get("zero_source") is True:
                result = _Expr.rational(0)
                self.memo[key] = result
                return result
            name = row.get("name", row.get("source_name", row.get("label", f"node[{ref}]")))
            result = _Expr.atom(f"{operation}:{name}")
        elif operation == "bound_variable":
            result = self._rewrite_bound(ref, env)
        elif operation == "sum":
            result = _Expr.rational(0)
            for child in _arg_list(row, ref):
                result = result + self.rewrite(child, env)
        elif operation == "negative":
            result = -self.rewrite(_field(row, "argument", "operand"), env)
        elif operation == "product":
            children = _arg_list(row, ref)
            # Accepted endpoint-zero weight roots are scoped source facts.
            # Discharge the product before touching its other operands so an
            # opaque leading control or recurrence body is never expanded at
            # the x=1 seam.
            if (
                self.zero_alias_scope_guard
                and self._active_scope == _PATCH_X1_SCOPE
                and any(child in self.zero_aliases for child in children)
            ):
                result = _Expr.rational(0)
                self.memo[key] = result
                return result
            result = _Expr.rational(1)
            for child in children:
                result = result * self.rewrite(child, env)
        elif operation in {"positive_quotient", "positive_function_quotient", "quotient"}:
            numerator = self.rewrite(_field(row, "numerator"), env)
            denominator = self.rewrite(_field(row, "denominator"), env)
            if denominator.is_zero:
                raise ValueError(f"node {ref} quotient denominator is zero")
            result = numerator * denominator.inverse()
        elif operation == "trim":
            result = self.rewrite(_field(row, "argument", "operand"), env)
        elif operation == "exact_positive_integer_expression":
            result = self.rewrite(_field(row, "argument", "operand"), env)
        elif operation == "analytic_unary":
            result = self._analytic(ref, row, self.rewrite(_field(row, "argument", "operand"), env))
        elif operation in {"compact_raw_beta", "compact_raw_beta_derivative", "compact_raw_beta_second_derivative"}:
            argument = self.rewrite(_field(row, "argument", "operand"), env)
            if argument.is_rational and abs(argument.rational_value) >= 1:
                self.flat_beta_nodes.add(ref)
                result = _Expr.rational(0)
            else:
                result = _Expr.factor((operation, argument.key()))
        elif operation == "function_substitution":
            if row.get(_Z_GUARD) is not True:
                raise ValueError(f"node {ref} lacks {_Z_GUARD}=True")
            expression = _handle(_field(row, "expression"), f"node {ref}.expression", len(self.nodes))
            variable = _field(row, "variable")
            value_handle = _handle(_field(row, "value"), f"node {ref}.value", len(self.nodes))
            value = self.rewrite(value_handle, env)
            scoped_env = dict(env)
            scoped_env[variable] = _EnvValue(value_handle, value)
            result = self.rewrite(expression, scoped_env)
        elif operation in {"definite_integral", "integral"}:
            lower = self.rewrite(_field(row, "lower"), env)
            upper = self.rewrite(_field(row, "upper"), env)
            if lower == upper:
                self.equal_bound_integrals.add(ref)
                result = _Expr.rational(0)
            else:
                raise _UnsupportedOperation(f"node {ref} integral has unequal bounds")
        elif operation == "original_lazy_flat_branch":
            predicate_ref = _handle(_field(row, "flat_predicate"), f"node {ref}.flat_predicate", len(self.nodes))
            if self._boolean(predicate_ref, env):
                result = self.rewrite(_field(row, "flat_value"), env)
            else:
                result = self.rewrite(_field(row, "active_body"), env)
        elif operation == "exact_real_comparison":
            # Boolean nodes are only consumed by lazy branches.  Expose exact
            # 0/1 for a direct field request without making a numeric claim.
            result = _Expr.rational(int(self._comparison(ref, row, self.rewrite(row["left"], env), self.rewrite(row["right"], env))))
        elif operation == "logical_or":
            result = _Expr.rational(int(any(self._boolean(child, env) for child in _arg_list(row, ref))))
        elif operation == "logical_and":
            result = _Expr.rational(int(all(self._boolean(child, env) for child in _arg_list(row, ref))))
        elif operation == "logical_not":
            result = _Expr.rational(int(not self._boolean(_field(row, "argument", "operand"), env)))
        elif operation == "original_Rm_unique_leading_implicit_control":
            child = row.get("value", row.get("argument"))
            result = self.rewrite(child, env) if child is not None else _Expr.atom(f"implicit_control:{ref}")
        else:  # pragma: no cover - operation validation is exhaustive
            raise _UnsupportedOperation(f"node {ref} uses unsupported operation {operation!r}")
        self.memo[key] = result
        return result


def _find_endpoint(value: Mapping[str, Any], names: Sequence[str]) -> Mapping[str, Any]:
    candidates: list[Any] = [value]
    for direct_name in ("inlet_endpoint", "actual_endpoint", "postrestore_endpoint", "endpoint_Rm"):
        direct = value.get(direct_name)
        if isinstance(direct, Mapping) and all(key in direct for key in ("fields", "radial")):
            return direct
    for name in names:
        current = value.get(name)
        if isinstance(current, Mapping):
            candidates.append(current)
    for candidate in candidates:
        for direct_name in ("inlet_endpoint", "actual_endpoint", "postrestore_endpoint", "endpoint_Rm"):
            direct = candidate.get(direct_name) if isinstance(candidate, Mapping) else None
            if isinstance(direct, Mapping):
                candidates.append(direct)
        for path in (("endpoints", "Rm"), ("Rm",), ("chart", "endpoints", "Rm")):
            try:
                endpoint = _path_value(candidate, (path,), "endpoint")
            except ValueError:
                continue
            if isinstance(endpoint, Mapping):
                return endpoint
        if all(key in candidate for key in ("fields", "radial")):
            return candidate
    raise ValueError("Rm endpoint binding is missing")


def _group_rows(endpoint: Mapping[str, Any], name: str, ordinary_orders: Sequence[int]) -> list[int]:
    candidate = endpoint.get(name)
    if candidate is None:
        raise ValueError(f"endpoint field group {name!r} is missing")
    return _rows(candidate, f"endpoint.{name}", ordinary_orders)


def _field_group(endpoint: Mapping[str, Any], names: Sequence[str], ordinary_orders: Sequence[int]) -> dict[str, list[int]]:
    for group_name in names:
        group = endpoint.get(group_name)
        if not isinstance(group, Mapping):
            continue
        result: dict[str, list[int]] = {}
        for field_name, values in group.items():
            if _is_sequence(values) or isinstance(values, Mapping):
                try:
                    result[str(field_name)] = _rows(values, f"endpoint.{group_name}.{field_name}", ordinary_orders)
                except ValueError:
                    continue
        if result:
            return result
    raise ValueError(f"endpoint field group {names[0]!r} is missing")


def _get_source_field_group(endpoint: Mapping[str, Any], ordinary_orders: Sequence[int]) -> dict[str, list[int]]:
    return _field_group(endpoint, ("fields", "field_roots", "actual_fields", "modified_fields"), ordinary_orders)


def _get_radial_group(endpoint: Mapping[str, Any], ordinary_orders: Sequence[int]) -> dict[str, list[int]]:
    return _field_group(endpoint, ("radial", "radial_fields", "actual_radial"), ordinary_orders)


def _node_is_exact_one(nodes: Sequence[Mapping[str, Any]], ref: Any) -> bool:
    if isinstance(ref, bool) or not isinstance(ref, int) or ref < 0 or ref >= len(nodes):
        return False
    row = nodes[ref]
    if not isinstance(row, Mapping):
        return False
    operation = row.get("operation", row.get("op"))
    if operation == "exact_rational":
        try:
            numerator = _fraction(row.get("numerator", row.get("value", 0)))
            denominator = _fraction(row.get("denominator", 1))
            return denominator != 0 and numerator / denominator == 1
        except (TypeError, ValueError, ZeroDivisionError):
            return False
    if operation == "constant":
        try:
            return _fraction(row.get("value")) == 1
        except (TypeError, ValueError):
            return False
    return False


def _endpoint_x1_guard(
    endpoint: Mapping[str, Any],
    nodes: Sequence[Mapping[str, Any]],
    ordinary_orders: Sequence[int],
) -> bool:
    """Confirm that an endpoint wrapper is evaluated at the exact x=1 node."""

    for key in ("endpoint_value", "x_value", "x1", "value", "at"):
        candidate = endpoint.get(key)
        if _node_is_exact_one(nodes, candidate):
            return True
        if isinstance(candidate, Mapping):
            for nested in ("node", "handle", "value", "root"):
                item = candidate.get(nested)
                if _node_is_exact_one(nodes, item):
                    return True
    for group_name in ("fields", "radial"):
        group = endpoint.get(group_name)
        if not isinstance(group, Mapping):
            continue
        for values in group.values():
            try:
                rows = _rows(values, f"endpoint.{group_name}", ordinary_orders)
            except ValueError:
                continue
            for ref in rows:
                row = nodes[ref]
                if not isinstance(row, Mapping):
                    continue
                operation = row.get("operation", row.get("op"))
                if operation != "function_substitution":
                    continue
                value = row.get("value")
                if _node_is_exact_one(nodes, value):
                    return True
    return False


def _substitution_parts(
    nodes: Sequence[Mapping[str, Any]],
    ref: int,
) -> tuple[int, Any, int] | None:
    if ref < 0 or ref >= len(nodes) or not isinstance(nodes[ref], Mapping):
        return None
    row = nodes[ref]
    if row.get("operation", row.get("op")) != "function_substitution":
        return None
    if row.get(_Z_GUARD) is not True:
        return None
    expression = row.get("expression")
    value = row.get("value")
    if isinstance(expression, bool) or not isinstance(expression, int):
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return expression, row.get("variable"), value


def _body_has_bound_variable(nodes: Sequence[Mapping[str, Any]], root: int) -> bool:
    """Reject a shared-body alias if it still contains a native bound variable."""

    seen: set[int] = set()

    def visit(ref: int) -> bool:
        if ref in seen or ref < 0 or ref >= len(nodes):
            return False
        seen.add(ref)
        row = nodes[ref]
        if not isinstance(row, Mapping):
            return True
        operation = row.get("operation", row.get("op"))
        if operation == "bound_variable":
            return True
        # Unsupported source recipes remain opaque atoms.  Equality of the
        # common body handle or an authenticated saved body pair is the guard;
        # expanding an old recipe here would defeat the bounded reducer.
        if operation not in _DECLARED_OPERATIONS:
            return False
        # Keep this structural scan independent of a reducer instance.
        children: list[int] = []
        if operation in {"sum", "product", "logical_or", "logical_and"}:
            values = row.get("arguments", row.get("operands"),)
            if _is_sequence(values):
                children.extend(item for item in values if isinstance(item, int) and not isinstance(item, bool))
        else:
            for key in (
                "argument", "operand", "numerator", "denominator", "expression",
                "integrand", "lower", "upper", "active_body", "flat_value",
                "flat_predicate", "left", "right", "value",
            ):
                item = row.get(key)
                if isinstance(item, int) and not isinstance(item, bool):
                    children.append(item)
        return any(visit(child) for child in children)

    return visit(root)


def _common_body_aliases(
    nodes: Sequence[Mapping[str, Any]],
    patch_rows: Mapping[str, Sequence[int]],
    inlet_rows: Mapping[str, Sequence[int]],
    ordinary_orders: Sequence[int],
    saved_mean: Mapping[str, Any],
) -> dict[int, _Alias]:
    """Register only guarded common source bodies for protected fields."""

    aliases: dict[int, _Alias] = {}
    protected_fields = ("m", "h", "k", "e", "p", "P0", "pressure")

    def add_pair(name: str, order: int, left: int, right: int, *, premise: str) -> None:
        if left < 0 or right < 0 or left >= len(nodes) or right >= len(nodes):
            return
        if _body_has_bound_variable(nodes, left) or _body_has_bound_variable(nodes, right):
            return
        atom = _source_label(name, order)
        descriptor = {
            "kind": "guarded_common_source_body",
            "field": name,
            "ordinary_order": order,
            "left_body": left,
            "right_body": right,
            "premise": premise,
            "Z_independent_substitution": True,
        }
        alias = _Alias(atom, right, descriptor, "guarded_common_source_body")
        previous_left = aliases.get(left)
        previous_right = aliases.get(right)
        if previous_left is not None and previous_left.atom != atom:
            return
        if previous_right is not None and previous_right.atom != atom:
            return
        aliases[left] = alias
        aliases[right] = alias

    for name in protected_fields:
        for order, patch_ref, inlet_ref in zip(
            ordinary_orders,
            patch_rows.get(name, ()),
            inlet_rows.get(name, ()),
        ):
            patch_parts = _substitution_parts(nodes, patch_ref)
            inlet_parts = _substitution_parts(nodes, inlet_ref)
            if patch_parts is None or inlet_parts is None:
                continue
            patch_body, _patch_variable, patch_value = patch_parts
            inlet_body, _inlet_variable, inlet_value = inlet_parts
            if not (_node_is_exact_one(nodes, patch_value) and _node_is_exact_one(nodes, inlet_value)):
                continue
            if patch_body == inlet_body:
                add_pair(name, order, inlet_body, patch_body, premise="shared_endpoint_body")

    certificate = saved_mean.get("certificate") if isinstance(saved_mean, Mapping) else None
    if isinstance(certificate, Mapping):
        for row in certificate.get("rows", ()):
            if not isinstance(row, Mapping) or row.get("passed") is not True:
                continue
            order = row.get("ordinary_Z_order")
            left = row.get("left_node")
            right = row.get("right_node")
            if isinstance(order, int) and order in ordinary_orders and isinstance(left, int) and isinstance(right, int):
                add_pair("m", order, left, right, premise="saved_mean_premise")
    return aliases


def _find_rows(group: Mapping[str, list[int]], names: Sequence[str], label: str) -> list[int]:
    for name in names:
        if name in group:
            return list(group[name])
    raise ValueError(f"{label} row is missing")


def _positive_parameter(node: int, nodes: Sequence[Mapping[str, Any]], label: str) -> dict[str, Any]:
    row = nodes[node]
    operation = _row_operation(row, node)
    positive = bool(
        row.get("positive", False)
        or row.get("positive_source", False)
        or row.get("positive_certificate")
        or "positive" in str(row.get("condition", "")).lower()
        or operation in {"current_original_source_parameter", "original_source_parameter"}
        or (operation == "analytic_unary" and str(row.get("name", "")).lower() == "exp")
    )
    return {"node": node, "operation": operation, "positive_guard": positive, "metadata": deepcopy(dict(row))}


def _extract_parameter(value: Mapping[str, Any], names: Sequence[str]) -> int | None:
    for name in names:
        candidate = value.get(name)
        if isinstance(candidate, int) and not isinstance(candidate, bool):
            return candidate
        if isinstance(candidate, Mapping):
            for nested in ("node", "handle", "value", "root"):
                item = candidate.get(nested)
                if isinstance(item, int) and not isinstance(item, bool):
                    return item
    return None


def _extract_logp_source(
    value: Mapping[str, Any],
    pstar_ref: int | None,
    nodes: Sequence[Mapping[str, Any]],
) -> int | None:
    """Find the real Z-independent log-P source under the frozen Pstar port."""

    explicit = _extract_parameter(
        value,
        ("logP", "log_P", "Pstar_log", "Pstar_log_argument", "logP_source"),
    )
    if explicit is not None:
        return explicit
    if pstar_ref is None or pstar_ref >= len(nodes):
        return None
    row = nodes[pstar_ref]
    operation = row.get("operation", row.get("op")) if isinstance(row, Mapping) else None
    name = row.get("name", row.get("function", row.get("unary"))) if isinstance(row, Mapping) else None
    argument = row.get("argument", row.get("operand")) if isinstance(row, Mapping) else None
    # The actual frozen Pstar port is exp(logP).  Its argument is a real
    # source parameter and is deliberately atomised before the exp rewrite;
    # this prevents ancestry expansion while preserving the exact exp/log
    # factor used by both source endpoint formulas.
    if (
        operation == "analytic_unary"
        and isinstance(name, str)
        and name.lower() == "exp"
        and isinstance(argument, int)
        and not isinstance(argument, bool)
    ):
        return argument
    return None


def _extract_am(leading: Mapping[str, Any], ordinary_orders: Sequence[int]) -> list[int] | None:
    for name in ("am", "amplitude", "normalized_amplitude", "Rm_amplitude"):
        candidate = leading.get(name)
        if candidate is None:
            continue
        try:
            return _rows(candidate, name, ordinary_orders)
        except ValueError:
            if isinstance(candidate, Mapping):
                for nested in ("rows", "ordinary", "jets", "values"):
                    if nested in candidate:
                        return _rows(candidate[nested], f"{name}.{nested}", ordinary_orders)
    normalization = leading.get("normalization")
    if isinstance(normalization, Mapping):
        return _extract_am(normalization, ordinary_orders)
    return None


def _extract_f0(value: Mapping[str, Any], ordinary_orders: Sequence[int]) -> list[int] | None:
    """Read the prescribed ``F0`` ordinary source rows without ancestry expansion."""

    for name in ("F0", "F_0", "background_F0", "frozen_F0"):
        candidate = value.get(name)
        if candidate is None:
            continue
        try:
            return _rows(candidate, name, ordinary_orders)
        except ValueError:
            if isinstance(candidate, Mapping):
                for nested in ("rows", "ordinary", "jets", "values", "source_rows"):
                    if nested in candidate:
                        return _rows(candidate[nested], f"{name}.{nested}", ordinary_orders)
    for nested_name in ("parameter_binding", "parameters", "normalization", "background"):
        nested = value.get(nested_name)
        if isinstance(nested, Mapping):
            result = _extract_f0(nested, ordinary_orders)
            if result is not None:
                return result
    return None


def _extract_control_rows(value: Mapping[str, Any]) -> list[int]:
    """Collect only explicitly bound leading-control handles.

    Controls are source functions, never numeric assignments.  The traversal
    is intentionally limited to named control containers so unrelated graph
    metadata cannot become an accidental alias.
    """

    result: list[int] = []
    seen: set[int] = set()

    def visit(candidate: Any) -> None:
        if isinstance(candidate, int) and not isinstance(candidate, bool):
            if candidate not in seen:
                seen.add(candidate)
                result.append(candidate)
            return
        if isinstance(candidate, Mapping):
            for child in candidate.values():
                visit(child)
            return
        if _is_sequence(candidate):
            for child in candidate:
                visit(child)

    for key in ("controls", "leading_controls", "implicit_controls"):
        if key in value:
            visit(value[key])
    return result


def _extract_weights(leading: Mapping[str, Any]) -> list[int]:
    values: list[int] = []
    for key in ("partial_weights", "Wpartial", "W_partial", "weights"):
        candidate = leading.get(key)
        if candidate is None:
            continue
        def visit(value: Any) -> None:
            if isinstance(value, int) and not isinstance(value, bool):
                values.append(value)
            elif isinstance(value, Mapping):
                for child in value.values():
                    visit(child)
            elif _is_sequence(value):
                for child in value:
                    visit(child)
        visit(candidate)
        if values:
            return values
    return values


def _extract_flat_roots(leading: Mapping[str, Any]) -> list[int]:
    """Return explicitly bound flat-branch probes without guessing graph rows."""

    values: list[int] = []
    for key in ("flat_beta_probe", "flat_branch", "compact_flat_branch", "beta_flat_branch"):
        candidate = leading.get(key)
        if candidate is None:
            continue
        if isinstance(candidate, int) and not isinstance(candidate, bool):
            values.append(candidate)
        elif _is_sequence(candidate):
            values.extend(
                item for item in candidate
                if isinstance(item, int) and not isinstance(item, bool)
            )
    return values


def _extract_probe_roots(leading: Mapping[str, Any]) -> list[int]:
    """Return explicitly named fresh equal-bound integral probes."""

    values: list[int] = []
    for key in ("fresh_empty_integral", "empty_integral_probe", "toy_empty_integral"):
        candidate = leading.get(key)
        if isinstance(candidate, int) and not isinstance(candidate, bool):
            values.append(candidate)
        elif _is_sequence(candidate):
            values.extend(
                item for item in candidate
                if isinstance(item, int) and not isinstance(item, bool)
            )
    return values


def _saved_mean_certificate(value: Any, ordinary_orders: Sequence[int]) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("saved_mean_premise must be a mapping")
    nested = value.get("leading_mean_certificate", value.get("mean", value))
    if not isinstance(nested, Mapping):
        raise ValueError("saved mean premise certificate is missing")
    if nested.get("passed") is not True or nested.get("source_only") is not True:
        raise ValueError("saved mean premise must be passed and source_only")
    packet_orders = nested.get("ordinary_orders", nested.get("ordinary_Z_orders", ()))
    try:
        packet_order_set = {int(order) for order in packet_orders}
    except (TypeError, ValueError):
        raise ValueError("saved mean premise ordinary orders are malformed")
    if not all(order in packet_order_set for order in ordinary_orders):
        raise ValueError("saved mean premise lacks a requested ordinary order")
    for forbidden in (
        "numeric_values",
        "caps",
        "global_opaque_alias",
        "global_opaque_aliases",
        "numeric_values_selected",
        "numeric_caps_used",
        "assigned_numeric_values",
        "global_opaque_alias_used",
        "source_values_selected",
    ):
        if nested.get(forbidden) not in (None, False, [], {}):
            raise ValueError(f"saved mean premise contains forbidden {forbidden}")
    partial = nested.get("partial_weights")
    controls = nested.get("controls")
    partial_passed = isinstance(partial, Mapping) and partial.get("passed") is True
    controls_passed = isinstance(controls, Mapping) and controls.get("passed") is True
    flat_certificate = value.get("flat", value.get("compact_flat_and_history_certificate"))
    if isinstance(flat_certificate, Mapping):
        if flat_certificate.get("passed") is not True or flat_certificate.get("source_only") is not True:
            raise ValueError("saved flat premise must be passed and source_only")
        for forbidden in (
            "numeric_values_selected",
            "numeric_caps_used",
            "assigned_numeric_values",
            "global_opaque_alias_used",
            "source_values_selected",
        ):
            if flat_certificate.get(forbidden) not in (None, False, [], {}):
                raise ValueError(f"saved flat premise contains forbidden {forbidden}")
    zero_nodes: set[int] = set()
    if isinstance(partial, Mapping):
        for row in partial.get("rows", ()):
            if isinstance(row, Mapping) and row.get("passed_exact_zero") is True:
                candidate = row.get("weight_node")
                if isinstance(candidate, int) and not isinstance(candidate, bool):
                    zero_nodes.add(candidate)
        for row in partial.get("empty_integral_rows", ()):
            if isinstance(row, Mapping) and row.get("passed") is True:
                candidate = row.get("integral_node")
                if isinstance(candidate, int) and not isinstance(candidate, bool):
                    zero_nodes.add(candidate)
    if isinstance(flat_certificate, Mapping):
        beta = flat_certificate.get("compact_raw_beta")
        if isinstance(beta, Mapping):
            for key in ("rows", "flat_derivative_rows"):
                for row in beta.get(key, ()):
                    if isinstance(row, Mapping) and row.get("endpoint_zero") is True:
                        candidate = row.get("node")
                        if isinstance(candidate, int) and not isinstance(candidate, bool):
                            zero_nodes.add(candidate)
    payload = {
        "report_sha256": value.get("report_sha256", value.get("sha256")),
        "certificate": deepcopy(dict(nested)),
    }
    return {
        "passed": True,
        "source_only": True,
        "ordinary_orders": list(ordinary_orders),
        "sha256": canonical_digest(payload),
        "partial_weights_consumed": bool(
            partial_passed
            or nested.get("partial_weight_rows")
            or nested.get("partial_weights_zero")
            or nested.get("all_partial_weights_exact_zero")
        ),
        "controls_inactivity_consumed": bool(
            controls_passed
            or nested.get("controls_inactive")
            or nested.get("controls_inactivity")
        ),
        "flat_certificate_consumed": bool(
            isinstance(flat_certificate, Mapping)
            and flat_certificate.get("passed") is True
        ),
        "zero_nodes": sorted(zero_nodes),
        "certificate_keys": sorted(str(key) for key in nested),
    }


def _node_digest(nodes: Sequence[Mapping[str, Any]]) -> str:
    return canonical_digest(list(nodes))


def _rows_digest(nodes: Sequence[Mapping[str, Any]], rows: Mapping[str, Sequence[int]]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for name, handles in rows.items():
        result[name] = [canonical_digest(nodes[index]) for index in handles]
    return result


def _json_expr(expression: _Expr) -> dict[str, Any]:
    return expression.json()


def _safe_row(
    reducer: _Reducer,
    nodes: Sequence[Mapping[str, Any]],
    patch_rows: Mapping[str, Sequence[int]],
    inlet_rows: Mapping[str, Sequence[int]],
    name: str,
    order_index: int,
    *,
    patch_scope: str = _PATCH_X1_SCOPE,
    inlet_scope: str = _INLET_X1_SCOPE,
) -> dict[str, Any]:
    patch_ref = patch_rows[name][order_index]
    inlet_ref = inlet_rows.get(name, patch_rows[name])[order_index]
    try:
        patch_expression = reducer.rewrite(patch_ref, scope=patch_scope)
        inlet_expression = reducer.rewrite(inlet_ref, scope=inlet_scope)
    except _UnsupportedOperation as exc:
        return {
            "field": name,
            "status": "unresolved",
            "patch_node": patch_ref,
            "inlet_node": inlet_ref,
            "gate": f"unsupported:{name}[{order_index}]",
            "error": str(exc),
        }
    return {
        "field": name,
        "status": "proved" if patch_expression == inlet_expression else "failed",
        "patch_node": patch_ref,
        "inlet_node": inlet_ref,
        "patch": _json_expr(patch_expression),
        "inlet": _json_expr(inlet_expression),
    }


def prove(
    nodes: Sequence[Mapping[str, Any]],
    leading_bindings: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    ordinary_orders: Sequence[int] = (0, 1, 2),
    saved_mean_premise: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Reduce selected actual leading endpoint rows.

    ``leading_bindings`` may be the raw patch snapshot or a wrapper containing
    ``true_Rm_patch_C6_bindings`` and ``prescribed_R100_Rm_C6_bindings``.  The
    function never mutates either argument.  Structural graph errors are
    rejected; unsupported reachable source operations are reported as
    unresolved row gates so already-proved rows remain available.
    """

    if not isinstance(expected_nodes_digest, str) or len(expected_nodes_digest) != 64:
        raise ValueError("expected_nodes_digest must be a SHA-256 hex digest")
    if _node_digest(nodes) != expected_nodes_digest:
        raise ValueError("nodes digest does not match expected_nodes_digest")
    orders = tuple(int(order) for order in ordinary_orders)
    if not orders or any(order < 0 for order in orders) or len(set(orders)) != len(orders):
        raise ValueError("ordinary_orders must be distinct nonnegative integers")
    original_nodes = deepcopy(nodes)
    original_bindings = deepcopy(leading_bindings)
    wrapper = _first_mapping(
        leading_bindings,
        ("true_Rm_patch_C6_bindings", "leading_bindings", "leading", "patch"),
    )
    prescribed = _first_mapping(
        leading_bindings,
        ("prescribed_R100_Rm_C6_bindings", "prescribed_bindings", "prescribed", "inlet"),
    )
    patch_endpoint = _find_endpoint(wrapper, ("endpoints", "Rm"))
    inlet_endpoint = _find_endpoint(
        prescribed,
        ("endpoints", "Rm", "inlet_endpoint", "actual_endpoint", "postrestore_endpoint"),
    )
    patch_fields = _get_source_field_group(patch_endpoint, orders)
    patch_radial = _get_radial_group(patch_endpoint, orders)
    inlet_fields = _get_source_field_group(inlet_endpoint, orders)
    inlet_radial = _get_radial_group(inlet_endpoint, orders)

    required = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
    for name in required:
        _find_rows(patch_fields, (name,), f"patch {name}")
        _find_rows(inlet_fields, (name,), f"inlet {name}")
    for name in ("E_y", "V_y"):
        _find_rows(patch_radial, (name,), f"patch {name}")
        _find_rows(inlet_radial, (name,), f"inlet {name}")

    # Protect ordinary source functions as independent atoms.  The same root
    # node may occur in defects and in the endpoint, so this preserves the
    # exact source function identity without assigning a polynomial body.
    protected: dict[int, _Alias] = {}
    endpoint_aliases: dict[int, _Alias] = {}
    protected_names = {"m", "h", "k", "e", "p", "P0"}
    for name in protected_names:
        for index, ref in enumerate(_find_rows(inlet_fields, (name,), f"inlet {name}")):
            descriptor = _canonical_descriptor(nodes[ref], ref)
            atom = _Alias(_source_label(name, orders[index]), ref, descriptor, f"inlet_{name}")
            protected.setdefault(ref, atom)
            endpoint_aliases[ref] = atom
    f0_rows = _extract_f0(prescribed, orders)
    if f0_rows is None:
        f0_rows = _extract_f0(wrapper, orders)
    if f0_rows is not None:
        for index, ref in enumerate(f0_rows):
            descriptor = _canonical_descriptor(nodes[ref], ref)
            protected.setdefault(
                ref,
                _Alias(f"source:F0[{orders[index]}]", ref, descriptor, "prescribed_F0"),
            )
    control_rows = _extract_control_rows(wrapper)
    for index, ref in enumerate(control_rows):
        descriptor = _canonical_descriptor(nodes[ref], ref)
        protected.setdefault(
            ref,
            _Alias(f"source:leading_control[{index}]", ref, descriptor, "leading_control"),
        )
    pstar_ref = _extract_parameter(wrapper, ("Pstar", "P_star", "pstar"))
    if pstar_ref is None:
        pstar_ref = _extract_parameter(prescribed, ("Pstar", "P_star", "pstar"))
    logp_ref = _extract_logp_source(wrapper, pstar_ref, nodes)
    if logp_ref is None:
        logp_ref = _extract_logp_source(prescribed, pstar_ref, nodes)
    if logp_ref is not None:
        descriptor = _canonical_descriptor(nodes[logp_ref], logp_ref)
        protected.setdefault(
            logp_ref,
            _Alias("source:logP", logp_ref, descriptor, "Pstar_log_source"),
        )

    # Validate source endpoint wrappers before reducing their arithmetic.  A
    # direct source atom is accepted; wrappers must carry the scoped guard.
    endpoint_guard_rows: list[dict[str, Any]] = []
    for group_name, group in (("patch", patch_fields), ("inlet", inlet_fields), ("patch_radial", patch_radial), ("inlet_radial", inlet_radial)):
        for field_name, rows in group.items():
            for order, ref in zip(orders, rows):
                row = nodes[ref]
                operation = _row_operation(row, ref)
                if operation == "function_substitution":
                    guarded = row.get(_Z_GUARD) is True
                    endpoint_guard_rows.append({"group": group_name, "field": field_name, "order": order, "node": ref, "guarded": guarded})
                    if not guarded:
                        raise ValueError(f"{group_name}.{field_name}[{order}] lacks {_Z_GUARD}=True")

    patch_rows = {name: _find_rows(patch_fields, (name,), f"patch {name}") for name in required}
    inlet_rows = {name: _find_rows(inlet_fields, (name,), f"inlet {name}") for name in required}
    patch_rows.update({name: _find_rows(patch_radial, (name,), f"patch {name}") for name in ("E_y", "V_y")})
    inlet_rows.update({name: _find_rows(inlet_radial, (name,), f"inlet {name}") for name in ("E_y", "V_y")})

    # Authenticate external mean/flat premises before any field or integral
    # traversal.  In particular, accepted partial weights and controls are
    # consumed as source facts rather than replayed by this reducer.
    mean_info: dict[str, Any]
    mean_gate = ""
    if saved_mean_premise is None:
        mean_info = {"passed": False, "source_only": False, "unresolved": "saved_mean_premise is required for m"}
        mean_gate = "saved_mean_premise_missing"
    else:
        try:
            mean_info = _saved_mean_certificate(saved_mean_premise, orders)
        except ValueError as exc:
            mean_info = {"passed": False, "source_only": False, "unresolved": str(exc)}
            mean_gate = "saved_mean_premise_invalid"

    patch_x1_scope_guard = _endpoint_x1_guard(patch_endpoint, nodes, orders)
    zero_scope_gate = ""
    if mean_info.get("zero_nodes") and not patch_x1_scope_guard:
        zero_scope_gate = "patch_endpoint_x1_scope_missing"
    body_aliases = _common_body_aliases(
        nodes,
        patch_rows,
        inlet_rows,
        orders,
        mean_info,
    )

    reducer = _Reducer(
        nodes,
        protected=protected,
        endpoint_aliases=endpoint_aliases,
        ordinary_orders=orders,
        body_aliases=body_aliases,
        zero_aliases=mean_info.get("zero_nodes", ()),
        zero_alias_scope_guard=patch_x1_scope_guard,
    )

    # Evaluate the positive frozen source parameter once.  An opaque atom is
    # accepted only with an explicit positive source guard; no numerical
    # value is materialised.
    pstar_info: dict[str, Any] = {"node": pstar_ref, "positive_guard": False}
    if pstar_ref is not None:
        pstar_info = _positive_parameter(pstar_ref, nodes, "Pstar")
        try:
            pstar_expr = reducer.rewrite(pstar_ref)
            pstar_info["expression"] = pstar_expr.json()
            pstar_info["exact_nonzero"] = not pstar_expr.is_zero
        except _UnsupportedOperation as exc:
            pstar_info["unresolved"] = str(exc)
    else:
        pstar_info["unresolved"] = "Pstar source parameter binding is missing"

    am_rows = _extract_am(wrapper, orders)
    am_results: list[dict[str, Any]] = []
    if am_rows is None:
        am_gate = "am_binding_missing"
    else:
        am_gate = ""
        for order, ref in zip(orders, am_rows):
            try:
                expression = reducer.rewrite(ref)
                am_results.append({"order": order, "node": ref, "expression": expression.json()})
            except _UnsupportedOperation as exc:
                am_gate = "am_normalization_unsupported"
                am_results.append({"order": order, "node": ref, "status": "unresolved", "error": str(exc)})
        # Ordinary rows are genuine ``Z`` derivatives and therefore need not
        # be equal to row 0.  Their per-row equality is checked by the E and
        # E_y endpoint comparisons below; no cross-row identity is inferred.

    f0_gate = "" if f0_rows is not None else "F0_source_rows_binding_missing"

    field_results: list[dict[str, Any]] = []
    for name in ("E", "V", "h", "k", "e", "p", "P0", "pressure", "E_y", "V_y"):
        source_rows = inlet_rows
        patch_name = name
        for index, order in enumerate(orders):
            row = _safe_row(reducer, nodes, {patch_name: patch_rows[patch_name]}, {name: source_rows[name]}, name, index)
            row["order"] = order
            field_results.append(row)

    # The accepted mean certificate owns the m identity.  Keep its endpoint
    # rows visible in this packet without traversing or re-proving them.
    for index, order in enumerate(orders):
        field_results.append(
            {
                "field": "m",
                "status": "external_premise",
                "order": order,
                "patch_node": patch_rows["m"][index],
                "inlet_node": inlet_rows["m"][index],
                "premise": "saved_mean_premise.leading_mean_certificate",
            }
        )

    # Consume selected partial-weight roots once.  If the saved premise has
    # already authenticated them, no old computation is replayed.  Otherwise
    # evaluate only rows explicitly bound by the new packet and require exact
    # zero; missing rows remain a gate.
    weight_refs = [] if mean_info.get("partial_weights_consumed") else _extract_weights(wrapper)
    weight_results: list[dict[str, Any]] = []
    for ref in weight_refs:
        try:
            expression = reducer.rewrite(ref, scope=_PATCH_X1_SCOPE)
            weight_results.append({"node": ref, "expression": expression.json(), "zero": expression.is_zero})
        except _UnsupportedOperation as exc:
            weight_results.append({"node": ref, "status": "unresolved", "error": str(exc), "zero": False})
    weights_zero = bool(mean_info.get("partial_weights_consumed")) or bool(weight_results) and all(row.get("zero") is True for row in weight_results)
    if not weights_zero:
        weight_gate = "partial_weight_zero_premise_missing_or_unproved"
    else:
        weight_gate = ""

    fresh_integral_results: list[dict[str, Any]] = []
    for ref in _extract_probe_roots(wrapper):
        try:
            expression = reducer.rewrite(ref, scope=_PATCH_X1_SCOPE)
            fresh_integral_results.append({"node": ref, "expression": expression.json(), "zero": expression.is_zero})
        except _UnsupportedOperation as exc:
            fresh_integral_results.append({"node": ref, "status": "unresolved", "error": str(exc), "zero": False})

    flat_results: list[dict[str, Any]] = []
    if not mean_info.get("flat_certificate_consumed"):
        for ref in _extract_flat_roots(wrapper):
            try:
                expression = reducer.rewrite(ref, scope=_PATCH_X1_SCOPE)
                flat_results.append({"node": ref, "expression": expression.json(), "zero": expression.is_zero})
            except _UnsupportedOperation as exc:
                flat_results.append({"node": ref, "status": "unresolved", "error": str(exc), "zero": False})
    flat_proved = bool(reducer.flat_beta_nodes) or bool(flat_results) and all(row.get("zero") is True for row in flat_results)
    flat_proved = flat_proved or bool(mean_info.get("flat_certificate_consumed"))
    flat_gate = "" if flat_proved else "compact_flat_beta_branch_not_reached"
    unresolved: list[str] = [
        gate for gate in (am_gate, f0_gate, mean_gate, weight_gate, flat_gate, zero_scope_gate) if gate
    ]
    for row in field_results:
        if row.get("status") == "unresolved":
            unresolved.append(str(row.get("gate")))
    patch_pass = all(row.get("status") == "proved" for row in field_results if row["field"] != "m")
    # The E and radial jet rows carry the actual normalisation burden.  The
    # remaining source-field rows are exact atom cancellations.
    passed = bool(
        patch_pass
        and not unresolved
        and bool(pstar_info.get("positive_guard"))
        and bool(pstar_info.get("exact_nonzero"))
        and am_rows is not None
        and f0_rows is not None
    )

    result = {
        "schema": "Rm_leading_seed_normalization_certificate_v1",
        "passed": passed,
        "conditional_partial": bool(field_results) and not passed,
        "source_only": True,
        "actual_saved_graph_run": False,
        "actual_saved_graph_loaded": False,
        "accepted_m_reproved": False,
        "ordinary_orders": list(orders),
        "expected_nodes_digest": expected_nodes_digest,
        "module_sha256": _module_sha256(),
        "source_graph": {
            "node_count": len(nodes),
            "nodes_digest": _node_digest(nodes),
            "row_digests": _rows_digest(nodes, {**patch_rows, **inlet_rows}),
        },
        "source_parameter": {
            "Pstar": pstar_info,
            "logP_source": logp_ref,
            "am": am_results,
            "F0_rows": f0_rows,
            "leading_control_rows": control_rows,
        },
        "fields": field_results,
        "mean_premise": mean_info,
        "partial_weights": {"rows": weight_results, "all_zero_or_saved": weights_zero},
        "flat_beta": {"rows": flat_results, "branch_proved": flat_proved},
        "endpoint_guards": endpoint_guard_rows,
        "exact_reductions": {
            "equal_bound_integrals": sorted(reducer.equal_bound_integrals),
            "fresh_empty_integrals": fresh_integral_results,
            "flat_beta_nodes": sorted(reducer.flat_beta_nodes),
            "accepted_zero_aliases": sorted(reducer.zero_aliases),
            "zero_alias_scope": _PATCH_X1_SCOPE,
            "zero_alias_scope_guard": patch_x1_scope_guard,
            "used_nodes": sorted(reducer.used_nodes),
            "alias_witnesses": reducer.alias_witnesses,
            "guarded_common_body_aliases": sorted(body_aliases),
        },
        "unresolved_gates": sorted(set(unresolved)),
        "scope": {
            "source_graph_only": True,
            "no_numeric_evaluation": True,
            "no_caps_or_midpoints_as_values": True,
            "no_generic_integer_replacement": True,
            "metadata_preserved": deepcopy(leading_bindings) == original_bindings,
            "nodes_input_unmutated": deepcopy(nodes) == original_nodes,
            "m_external_saved_premise_only": True,
        },
    }
    return result


__all__ = ["canonical_digest", "prove"]
