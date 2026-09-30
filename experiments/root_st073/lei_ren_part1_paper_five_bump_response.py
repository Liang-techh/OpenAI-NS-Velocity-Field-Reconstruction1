"""Formal finite response of the Section 10 five-bump moment map.

The five centered defects are retained as independent formal variables
``d1`` through ``d5``.  If ``A`` is the fixed linear five-bump matrix and
``Q`` is the quadratic contribution, this module expands

    A h + Q(h, h) + d = 0

by total degree in the five defect variables.  It is deliberately a formal
response: coefficients can be ordinary MP numbers, ``PressureWidthJet``
objects, ``AxialDual`` objects, or another compatible scalar ring, while the
defect variables remain sparse monomials until ``evaluate_terms`` is called.

Every monomial contribution is returned separately.  This keeps tiny d3/d5
terms visible even when their summed scalar response is dominated by larger
terms.  The finite response does not assert convergence, a contraction ball,
temporal recursion, global closure, or a Navier--Stokes field repair.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import json
import operator
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
except (ImportError, ValueError):
    from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap  # type: ignore


VARIABLE_COUNT = 5
COEFFICIENT_NAMES = ("c1", "c2", "xi1", "xi2", "xi3")
DEFECT_NAMES = ("d1", "d2", "d3", "d4", "d5")


def _integer(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    return int(result)


def _exact_zero(value: Any) -> bool:
    """Prune only ordinary exact scalar zeros, never tiny ring atoms."""

    if isinstance(value, (int, mp.mpf, mp.mpc)):
        return value == 0
    return False


def _multiply(left: Any, right: Any) -> Any:
    """Multiply generic coefficients with a conservative dispatch fallback."""

    try:
        return left * right
    except (TypeError, ValueError, AttributeError):
        return right * left


def _divide(left: Any, right: Any) -> Any:
    try:
        return left / right
    except (TypeError, ValueError, AttributeError):
        return (right ** -1) * left


def _sum(values: Sequence[Any], template: Any = None) -> Any:
    if not values:
        if template is not None:
            return _multiply(template, 0)
        return mp.mpf(0)
    result = values[0]
    for value in values[1:]:
        result = result + value
    return result


def _zero_like(value: Any) -> Any:
    return _multiply(value, 0)


def _normalize_index(index: Any) -> tuple[int, ...]:
    try:
        result = tuple(_integer(value, "monomial power") for value in index)
    except TypeError as exc:
        raise TypeError("a monomial index must be an iterable") from exc
    if len(result) != VARIABLE_COUNT or any(value < 0 for value in result):
        raise ValueError("a monomial index must have five nonnegative powers")
    return result


def _coerce_defects(defects: Any) -> tuple[Any, ...]:
    """Accept a five-sequence or the common row/name mappings."""

    if isinstance(defects, Mapping):
        for keys in (
            tuple(range(1, VARIABLE_COUNT + 1)),
            tuple(range(VARIABLE_COUNT)),
            DEFECT_NAMES,
            tuple(str(index) for index in range(1, VARIABLE_COUNT + 1)),
            tuple(str(index) for index in range(VARIABLE_COUNT)),
            tuple(f"d{index}" for index in range(1, VARIABLE_COUNT + 1)),
        ):
            if all(key in defects for key in keys):
                return tuple(defects[key] for key in keys)
        raise ValueError("defect mapping must provide five d1..d5 entries")
    values = tuple(defects)
    if len(values) != VARIABLE_COUNT:
        raise ValueError("defects must contain five values")
    return values


def _monomial_degree(index: tuple[int, ...]) -> int:
    return sum(index)


class SparseDefectPolynomial:
    """Sparse normalized polynomial in five independent defect variables.

    The mapping keys are five-tuples ``(p1,...,p5)`` and values are generic
    coefficients.  The optional ``max_degree`` applies total-degree
    truncation to products; no coefficient is converted to MP merely because
    it belongs to this polynomial.
    """

    variable_count = VARIABLE_COUNT

    def __init__(
        self,
        terms: Mapping[Any, Any] | None = None,
        *,
        max_degree: int | None = None,
    ) -> None:
        if max_degree is not None:
            max_degree = _integer(max_degree, "max_degree")
            if max_degree < 0:
                raise ValueError("max_degree must be nonnegative")
        self.max_degree = max_degree
        normalized: dict[tuple[int, ...], Any] = {}
        if terms is not None:
            for raw_index, coefficient in terms.items():
                index = _normalize_index(raw_index)
                if max_degree is not None and _monomial_degree(index) > max_degree:
                    continue
                if index in normalized:
                    normalized[index] = normalized[index] + coefficient
                else:
                    normalized[index] = coefficient
        self.terms = {
            index: coefficient
            for index, coefficient in normalized.items()
            if not _exact_zero(coefficient)
        }

    @classmethod
    def _from_terms(
        cls,
        terms: Mapping[tuple[int, ...], Any],
        max_degree: int | None,
    ) -> "SparseDefectPolynomial":
        result = cls.__new__(cls)
        result.max_degree = max_degree
        result.terms = {
            index: coefficient
            for index, coefficient in terms.items()
            if not _exact_zero(coefficient)
        }
        return result

    @classmethod
    def zero(cls, *, max_degree: int | None = None) -> "SparseDefectPolynomial":
        return cls._from_terms({}, max_degree)

    @classmethod
    def constant(
        cls,
        value: Any,
        *,
        max_degree: int | None = None,
    ) -> "SparseDefectPolynomial":
        return cls({(0,) * VARIABLE_COUNT: value}, max_degree=max_degree)

    @classmethod
    def monomial(
        cls,
        index: Any,
        coefficient: Any = 1,
        *,
        max_degree: int | None = None,
    ) -> "SparseDefectPolynomial":
        return cls({_normalize_index(index): coefficient}, max_degree=max_degree)

    def copy(self) -> "SparseDefectPolynomial":
        return self._from_terms(dict(self.terms), self.max_degree)

    def __bool__(self) -> bool:
        return bool(self.terms)

    def _limit_with(self, other: Any) -> int | None:
        if not isinstance(other, SparseDefectPolynomial):
            return self.max_degree
        if self.max_degree is None:
            return other.max_degree
        if other.max_degree is None:
            return self.max_degree
        return min(self.max_degree, other.max_degree)

    def _coerce(self, other: Any) -> "SparseDefectPolynomial":
        if isinstance(other, SparseDefectPolynomial):
            return other
        return SparseDefectPolynomial.constant(other, max_degree=self.max_degree)

    def truncate(self, max_degree: int) -> "SparseDefectPolynomial":
        limit = _integer(max_degree, "max_degree")
        if limit < 0:
            raise ValueError("max_degree must be nonnegative")
        return SparseDefectPolynomial._from_terms(
            {
                index: coefficient
                for index, coefficient in self.terms.items()
                if _monomial_degree(index) <= limit
            },
            limit,
        )

    def homogeneous(self, degree: int) -> "SparseDefectPolynomial":
        degree = _integer(degree, "degree")
        if degree < 0:
            raise ValueError("degree must be nonnegative")
        return SparseDefectPolynomial._from_terms(
            {
                index: coefficient
                for index, coefficient in self.terms.items()
                if _monomial_degree(index) == degree
            },
            self.max_degree,
        )

    def coefficient(self, index: Any) -> Any:
        return self.terms.get(_normalize_index(index), mp.mpf(0))

    def __add__(self, other: Any) -> "SparseDefectPolynomial":
        right = self._coerce(other)
        terms = dict(self.terms)
        for index, coefficient in right.terms.items():
            terms[index] = terms.get(index, _zero_like(coefficient)) + coefficient
        return self._from_terms(terms, self._limit_with(right))

    __radd__ = __add__

    def __neg__(self) -> "SparseDefectPolynomial":
        return self._from_terms(
            {index: -coefficient for index, coefficient in self.terms.items()},
            self.max_degree,
        )

    def __sub__(self, other: Any) -> "SparseDefectPolynomial":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "SparseDefectPolynomial":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "SparseDefectPolynomial":
        if not isinstance(other, SparseDefectPolynomial):
            return self._from_terms(
                {index: _multiply(coefficient, other) for index, coefficient in self.terms.items()},
                self.max_degree,
            )
        limit = self._limit_with(other)
        terms: dict[tuple[int, ...], Any] = {}
        for left_index, left_coefficient in self.terms.items():
            for right_index, right_coefficient in other.terms.items():
                index = tuple(
                    left_index[position] + right_index[position]
                    for position in range(VARIABLE_COUNT)
                )
                if limit is not None and _monomial_degree(index) > limit:
                    continue
                value = _multiply(left_coefficient, right_coefficient)
                if index in terms:
                    terms[index] = terms[index] + value
                else:
                    terms[index] = value
        return self._from_terms(terms, limit)

    def __rmul__(self, other: Any) -> "SparseDefectPolynomial":
        return self * other

    def __truediv__(self, other: Any) -> "SparseDefectPolynomial":
        return self._from_terms(
            {index: _divide(coefficient, other) for index, coefficient in self.terms.items()},
            self.max_degree,
        )

    def __pow__(self, exponent: Any) -> "SparseDefectPolynomial":
        exponent = _integer(exponent, "polynomial exponent")
        if exponent < 0:
            raise ValueError("negative polynomial powers are not supported")
        result = SparseDefectPolynomial.constant(1, max_degree=self.max_degree)
        factor = self
        while exponent:
            if exponent & 1:
                result = result * factor
            exponent >>= 1
            if exponent:
                factor = factor * factor
        return result

    def evaluate(self, defects: Sequence[Any]) -> Any:
        values = _coerce_defects(defects)
        contributions = []
        for index, coefficient in self.terms.items():
            monomial = 1
            for value, power in zip(values, index):
                if power:
                    monomial = _multiply(monomial, value ** power)
            contributions.append(_multiply(coefficient, monomial))
        return _sum(contributions)

    def __repr__(self) -> str:
        return (
            f"SparseDefectPolynomial({self.terms!r}, "
            f"max_degree={self.max_degree!r})"
        )


def _linear_apply(moment_map: Any, vector: Sequence[SparseDefectPolynomial]) -> tuple[SparseDefectPolynomial, ...]:
    matrix = getattr(moment_map, "linear_matrix", getattr(moment_map, "matrix", None))
    if matrix is None:
        raise AttributeError("moment_map must expose linear_matrix or matrix")
    result = []
    for row in range(VARIABLE_COUNT):
        terms = [vector[column] * matrix[row, column] for column in range(VARIABLE_COUNT)]
        result.append(_sum(terms, vector[0]))
    return tuple(result)


def _form(
    left: Sequence[SparseDefectPolynomial],
    matrix: Any,
    right: Sequence[SparseDefectPolynomial],
) -> SparseDefectPolynomial:
    terms = []
    for row in range(len(left)):
        for column in range(len(right)):
            terms.append(left[row] * matrix[row, column] * right[column])
    return _sum(terms, left[0] if left else (right[0] if right else None))


def _inverse_square(am: Any) -> Any:
    try:
        return am ** (-2)
    except (TypeError, ValueError, ZeroDivisionError):
        return 1 / (am * am)


def _symmetric_quadratic(
    moment_map: Any,
    left: Sequence[SparseDefectPolynomial],
    right: Sequence[SparseDefectPolynomial],
    am: Any,
) -> tuple[SparseDefectPolynomial, ...]:
    """Return the polarized Q(left,right), with Q(h,h)=map.quadratic(h)."""

    weights = getattr(moment_map, "quadratic_weights", None)
    if weights is None:
        raise AttributeError("moment_map must expose quadratic_weights")
    left_c, right_c = left[:2], right[:2]
    left_xi, right_xi = left[2:], right[2:]
    fg = (
        _form(left_c, weights["fg_sqrt"], right_xi)
        + _form(right_c, weights["fg_sqrt"], left_xi)
    ) / 2
    gg = _form(left_c, weights["gg"], right_c)
    ff = _form(left_xi, weights["ff"], right_xi)
    ff_over_x = _form(left_xi, weights["ff_over_x"], right_xi)
    zero = SparseDefectPolynomial.zero(max_degree=left[0].max_degree)
    return (
        zero,
        fg,
        zero,
        gg * _inverse_square(am) - ff / 2,
        ff_over_x / 2,
    )


def _quadratic_total(
    moment_map: Any,
    left: Sequence[SparseDefectPolynomial],
    right: Sequence[SparseDefectPolynomial],
    am: Any,
) -> tuple[SparseDefectPolynomial, ...]:
    return _symmetric_quadratic(moment_map, left, right, am)


def _zero_vector(max_degree: int) -> tuple[SparseDefectPolynomial, ...]:
    return tuple(
        SparseDefectPolynomial.zero(max_degree=max_degree)
        for _ in range(VARIABLE_COUNT)
    )


def _add_vectors(
    left: Sequence[SparseDefectPolynomial],
    right: Sequence[SparseDefectPolynomial],
) -> tuple[SparseDefectPolynomial, ...]:
    return tuple(left[index] + right[index] for index in range(VARIABLE_COUNT))


def _negative_vector(
    vector: Sequence[SparseDefectPolynomial],
) -> tuple[SparseDefectPolynomial, ...]:
    return tuple(-value for value in vector)


def _homogeneous_vector(
    vector: Sequence[SparseDefectPolynomial], degree: int
) -> tuple[SparseDefectPolynomial, ...]:
    return tuple(value.homogeneous(degree) for value in vector)


def _basis_defects(max_degree: int) -> tuple[SparseDefectPolynomial, ...]:
    return tuple(
        SparseDefectPolynomial.monomial(
            tuple(1 if position == variable else 0 for position in range(VARIABLE_COUNT)),
            max_degree=max_degree,
        )
        for variable in range(VARIABLE_COUNT)
    )


class FiveBumpResponse:
    """Formal response and coefficientwise residual receipts."""

    def __init__(
        self,
        moment_map: Any,
        am: Any,
        degree: int,
        defect_basis: tuple[SparseDefectPolynomial, ...],
        h_by_degree: Mapping[int, tuple[SparseDefectPolynomial, ...]],
        q_by_degree: Mapping[int, tuple[SparseDefectPolynomial, ...]],
        residual_by_degree: Mapping[int, tuple[SparseDefectPolynomial, ...]],
    ) -> None:
        self.moment_map = moment_map
        self.Am = am
        self.degree = degree
        self.defect_basis = defect_basis
        self.h_by_degree = dict(h_by_degree)
        self.q_by_degree = dict(q_by_degree)
        self.residual_by_degree = dict(residual_by_degree)
        self.h = tuple(
            _sum([self.h_by_degree[n][row] for n in self.h_by_degree], self.h_by_degree[1][row])
            for row in range(VARIABLE_COUNT)
        )
        self.residual = tuple(
            _sum(
                [self.residual_by_degree[n][row] for n in self.residual_by_degree],
                self.residual_by_degree[1][row],
            )
            for row in range(VARIABLE_COUNT)
        )
        self.metadata = {
            "degree": degree,
            "defect_variables": DEFECT_NAMES,
            "coefficient_order": COEFFICIENT_NAMES,
            "finite_response": True,
            "coefficientwise_residual_through_degree": True,
            "convergent_series_certified": False,
            "contraction_ball_certified": False,
            "temporal_recursion_certified": False,
            "global_field_installed": False,
        }
        # Flattened coefficient maps are the stable integration API.  Each
        # key is one five-tuple defect monomial and each value is the five
        # bump-coefficient response vector for that monomial.
        self.coefficients = self._flatten_by_monomial(self.h_by_degree)
        self.q_coefficients = self._flatten_by_monomial(self.q_by_degree)
        self.residual_coefficients = self._flatten_by_monomial(
            self.residual_by_degree
        )

    @staticmethod
    def _flatten_by_monomial(
        vectors_by_degree: Mapping[
            int, tuple[SparseDefectPolynomial, ...]
        ]
    ) -> dict[tuple[int, ...], tuple[Any, ...]]:
        indices: set[tuple[int, ...]] = set()
        for vector in vectors_by_degree.values():
            for polynomial in vector:
                indices.update(polynomial.terms)
        flattened: dict[tuple[int, ...], tuple[Any, ...]] = {}
        for index in sorted(indices):
            entries: list[Any] = []
            for row in range(VARIABLE_COUNT):
                coefficients = [
                    vector[row].terms[index]
                    for vector in vectors_by_degree.values()
                    if index in vector[row].terms
                ]
                entries.append(_sum(coefficients))
            flattened[index] = tuple(entries)
        return flattened

    def _evaluated_polynomial_terms(
        self,
        polynomial: SparseDefectPolynomial,
        defects: Sequence[Any],
    ) -> dict[tuple[int, ...], Any]:
        values = _coerce_defects(defects)
        result: dict[tuple[int, ...], Any] = {}
        for index, coefficient in polynomial.terms.items():
            monomial = 1
            for value, power in zip(values, index):
                if power:
                    monomial = _multiply(monomial, value ** power)
            result[index] = _multiply(coefficient, monomial)
        return result

    @staticmethod
    def _sum_term_mapping(mapping: Mapping[tuple[int, ...], Any]) -> Any:
        return _sum(list(mapping.values()))

    def _evaluate_coefficients(
        self,
        coefficients: Mapping[tuple[int, ...], Sequence[Any]],
        defects: Sequence[Any],
    ) -> dict[tuple[int, ...], tuple[Any, ...]]:
        values = _coerce_defects(defects)
        result: dict[tuple[int, ...], tuple[Any, ...]] = {}
        for index, vector in coefficients.items():
            monomial = 1
            for value, power in zip(values, index):
                if power:
                    monomial = _multiply(monomial, value ** power)
            result[index] = tuple(
                _multiply(coefficient, monomial) for coefficient in vector
            )
        return result

    @staticmethod
    def _aggregate_coefficient_terms(
        terms: Mapping[tuple[int, ...], Sequence[Any]],
    ) -> tuple[Any, ...]:
        return tuple(
            _sum([vector[row] for vector in terms.values()])
            for row in range(VARIABLE_COUNT)
        )

    def evaluate_terms(
        self, defects: Sequence[Any]
    ) -> dict[tuple[int, ...], tuple[Any, ...]]:
        """Return every evaluated h contribution by defect monomial.

        The result is ``multiindex -> five-vector``.  It intentionally keeps
        the monomial key instead of returning only the sum.  Use
        :meth:`evaluate_residual_terms` for coefficientwise residuals and
        :meth:`evaluate_receipt` when aggregated vectors and per-degree
        diagnostics are also desired.
        """

        return self._evaluate_coefficients(self.coefficients, defects)

    def evaluate_residual_terms(
        self, defects: Sequence[Any]
    ) -> dict[tuple[int, ...], tuple[Any, ...]]:
        """Evaluate the retained coefficientwise residual by monomial."""

        return self._evaluate_coefficients(self.residual_coefficients, defects)

    def evaluate_receipt(self, defects: Sequence[Any]) -> dict[str, Any]:
        """Return separated h, Q, residual, and per-degree evaluations."""

        h_terms: dict[int, dict[tuple[int, ...], Any]] = {}
        q_terms: dict[int, dict[tuple[int, ...], Any]] = {}
        residual_terms: dict[int, dict[tuple[int, ...], Any]] = {}
        h_by_degree: dict[int, tuple[Any, ...]] = {}
        q_by_degree: dict[int, tuple[Any, ...]] = {}
        residual_by_degree: dict[int, tuple[Any, ...]] = {}

        for degree, vector in self.h_by_degree.items():
            for row, polynomial in enumerate(vector):
                h_terms.setdefault(row, {}).update(
                    self._evaluated_polynomial_terms(polynomial, defects)
                )
            h_by_degree[degree] = tuple(
                self._sum_term_mapping(
                    self._evaluated_polynomial_terms(vector[row], defects)
                )
                for row in range(VARIABLE_COUNT)
            )
        for degree, vector in self.q_by_degree.items():
            for row, polynomial in enumerate(vector):
                q_terms.setdefault(row, {}).update(
                    self._evaluated_polynomial_terms(polynomial, defects)
                )
            q_by_degree[degree] = tuple(
                self._sum_term_mapping(
                    self._evaluated_polynomial_terms(vector[row], defects)
                )
                for row in range(VARIABLE_COUNT)
            )
        for degree, vector in self.residual_by_degree.items():
            for row, polynomial in enumerate(vector):
                residual_terms.setdefault(row, {}).update(
                    self._evaluated_polynomial_terms(polynomial, defects)
                )
            residual_by_degree[degree] = tuple(
                self._sum_term_mapping(
                    self._evaluated_polynomial_terms(vector[row], defects)
                )
                for row in range(VARIABLE_COUNT)
            )

        h = self._aggregate_coefficient_terms(self.evaluate_terms(defects))
        q = self._aggregate_coefficient_terms(
            self._evaluate_coefficients(self.q_coefficients, defects)
        )
        residual = self._aggregate_coefficient_terms(
            self.evaluate_residual_terms(defects)
        )
        defect_terms = {
            row: self._evaluated_polynomial_terms(self.defect_basis[row], defects)
            for row in range(VARIABLE_COUNT)
        }
        return {
            "h_terms": h_terms,
            "q_terms": q_terms,
            "residual_terms": residual_terms,
            "defect_terms": defect_terms,
            "h_by_degree": h_by_degree,
            "q_by_degree": q_by_degree,
            "residual_by_degree": residual_by_degree,
            "h": h,
            "q": q,
            "residual": residual,
            "degree": self.degree,
            "finite_response": True,
            "coefficientwise_residual_through_degree": True,
            "convergent_series_certified": False,
            "contraction_ball_certified": False,
            "temporal_recursion_certified": False,
            "global_field_installed": False,
        }


def build_response(moment_map: Any, Am: Any, degree: int = 3) -> FiveBumpResponse:
    """Build the total-degree formal response through ``degree``."""

    degree = _integer(degree, "degree")
    if degree < 1:
        raise ValueError("degree must be at least 1")
    if not hasattr(moment_map, "linear_inverse"):
        raise TypeError("moment_map must expose linear_inverse")
    if not hasattr(moment_map, "quadratic_weights"):
        raise TypeError("moment_map must expose quadratic_weights")

    defects = _basis_defects(degree)
    h_by_degree: dict[int, tuple[SparseDefectPolynomial, ...]] = {}
    q_by_degree: dict[int, tuple[SparseDefectPolynomial, ...]] = {}
    h1_rhs = _negative_vector(defects)
    h_by_degree[1] = tuple(moment_map.linear_inverse(h1_rhs))

    for current_degree in range(2, degree + 1):
        partial = _zero_vector(degree)
        for prior_degree in range(1, current_degree):
            partial = _add_vectors(partial, h_by_degree[prior_degree])
        q_full = _quadratic_total(moment_map, partial, partial, Am)
        q_current = _homogeneous_vector(q_full, current_degree)
        q_by_degree[current_degree] = q_current
        h_by_degree[current_degree] = tuple(
            moment_map.linear_inverse(_negative_vector(q_current))
        )

    zero = _zero_vector(degree)
    q_by_degree[1] = zero
    residual_by_degree: dict[int, tuple[SparseDefectPolynomial, ...]] = {}
    for current_degree in range(1, degree + 1):
        linear = _linear_apply(moment_map, h_by_degree[current_degree])
        q = q_by_degree.get(current_degree, zero)
        residual = _add_vectors(linear, q)
        if current_degree == 1:
            residual = _add_vectors(residual, defects)
        residual_by_degree[current_degree] = residual

    return FiveBumpResponse(
        moment_map,
        Am,
        degree,
        defects,
        h_by_degree,
        q_by_degree,
        residual_by_degree,
    )


def evaluate_terms(response: FiveBumpResponse, defects: Sequence[Any]) -> dict[str, Any]:
    """Functional alias for ``response.evaluate_terms(defects)``."""

    if not isinstance(response, FiveBumpResponse):
        raise TypeError("response must be a FiveBumpResponse")
    return response.evaluate_terms(defects)


def _smoke() -> dict[str, Any]:
    """Small formal smoke used by the module's direct execution path."""

    with mp.workdps(90):
        moment_map = FiveBumpMomentMap(precision=90, order=48)
        am = mp.mpf("1.7")
        response = build_response(moment_map, am, degree=3)
        defects = tuple(mp.mpf(value) for value in ("1e-3", "-2e-3", "3e-3", "-4e-3", "5e-3"))
        evaluated = response.evaluate_receipt(defects)
        coefficient_residual = max(
            (
                abs(coefficient)
                for vector in response.residual_by_degree.values()
                for polynomial in vector
                for coefficient in polynomial.terms.values()
                if isinstance(coefficient, (int, mp.mpf, mp.mpc))
            ),
            default=mp.mpf(0),
        )
        report = {
            "degree": response.degree,
            "monomial_counts_by_h_degree": {
                str(n): sum(len(polynomial.terms) for polynomial in vector)
                for n, vector in response.h_by_degree.items()
            },
            "coefficientwise_residual_max": mp.nstr(coefficient_residual, 20),
            "evaluated_residual": [mp.nstr(value, 20) for value in evaluated["residual"]],
            "finite_response": True,
            "convergent_series_certified": False,
            "global_field_installed": False,
        }
        return report


if __name__ == "__main__":
    print(json.dumps(_smoke(), indent=2))


__all__ = [
    "SparseDefectPolynomial",
    "FiveBumpResponse",
    "build_response",
    "evaluate_terms",
    "DEFECT_NAMES",
    "COEFFICIENT_NAMES",
]
