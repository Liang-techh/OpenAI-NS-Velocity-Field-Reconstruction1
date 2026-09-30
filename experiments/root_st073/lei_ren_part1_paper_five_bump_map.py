"""Section 10.2's finite five-bump nonlinear moment map (equation 10.8).

The three angular bumps are

``gamma_j(x) = beta_{1/40}(x-s_j)``, ``s=(5/4,3/2,7/4)``,

and the two axial bumps are ``b_1=gamma_1`` and ``b_2=gamma_3``.  The
coefficient order is ``h=(c_1,c_2,xi_1,xi_2,xi_3)``.  This module keeps the
fixed scalar quadrature data separate from the input coefficients, so the
same map works with ordinary MP scalars, ``PressureWidthJet`` values, and
``AxialDual`` values.

The returned rows are exactly the five integrated expressions in (10.8):

* ``integral(g)``;
* ``integral(x**(3/5)*g + sqrt(x)*f*g)``;
* ``integral(sqrt(x)*f)``;
* ``integral(-x**(1/10)*f + Am**(-2)*g**2 - f**2/2)``;
* ``integral(x**(-9/10)*f + f**2/(2*x))``.

All quadrature is finite order and is retained as a receipt.  This is a
coefficient map only: it does not claim a contraction, a cone bound, an
enclosure, or a full field closure.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import json
import operator
from pathlib import Path
from typing import Any

import mpmath as mp


SOURCE = "https://arxiv.org/html/2609.35406v2"
SOURCE_SECTION = "Part I Section 10.2, equation (10.8), with (10.9)--(10.10); Interval I.4 bump"
COEFFICIENT_NAMES = ("c1", "c2", "xi1", "xi2", "xi3")
CENTERS = ("1.25", "1.5", "1.75")
RADIUS = "0.025"
MOMENT_ROWS = ("z", "z_weighted", "theta", "z_theta", "p")


def _index(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return int(result)


def _raw_bump(t: Any) -> mp.mpf:
    """Return the paper's normalized-family raw bump ``b(t)``.

    Interval I.4 uses ``b(t)=exp(-1/(1-t^2))`` on ``|t|<1`` and zero at
    and outside the endpoints.  Keeping this expression in MP arithmetic
    avoids the float conversion used by older exploratory adapters.
    """

    coordinate = mp.mpf(t)
    if abs(coordinate) >= 1:
        return mp.mpf(0)
    return mp.exp(-1 / (1 - coordinate * coordinate))


def _zero_like(value: Any) -> Any:
    """Construct an additive zero in a scalar or coefficient ring."""

    try:
        return value * 0
    except Exception:
        return mp.mpf(0)


def _sum_terms(terms: Iterable[Any], template: Any | None = None) -> Any:
    terms = iter(terms)
    if template is None:
        try:
            first = next(terms)
        except StopIteration:
            return mp.mpf(0)
        total = first
    else:
        total = _zero_like(template)
    for term in terms:
        total = total + term
    return total


def _as_sequence(value: Any, name: str, length: int = 5) -> tuple[Any, ...]:
    """Read a five-vector while preserving arbitrary coefficient objects."""

    if hasattr(value, "rows") and hasattr(value, "cols"):
        try:
            if int(value.rows) * int(value.cols) != length:
                raise ValueError(f"{name} must have exactly {length} entries")
            result = tuple(value[index] for index in range(length))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must have exactly {length} entries") from exc
    else:
        try:
            result = tuple(value)
        except TypeError as exc:
            raise TypeError(f"{name} must be an iterable of {length} entries") from exc
        if len(result) != length:
            raise ValueError(f"{name} must have exactly {length} entries")
    return result


def _mp_matrix_to_lists(matrix: mp.matrix) -> list[list[mp.mpf]]:
    return [
        [matrix[row, col] for col in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


class FiveBumpMomentMap:
    """Arbitrary-precision five-bump map with generic coefficient support."""

    coefficient_order = COEFFICIENT_NAMES
    centers_as_strings = CENTERS
    radius_as_string = RADIUS
    finite_quadrature = True
    quadrature_enclosed = False

    def __init__(self, *, precision: int = 160, order: int = 64) -> None:
        if _index(precision, "precision") < 80:
            raise ValueError("precision must be an integer >= 80")
        if _index(order, "order") < 8:
            raise ValueError("order must be an integer >= 8")
        self.precision = int(precision)
        self.order = int(order)
        with mp.workdps(self.precision):
            self.radius = mp.mpf(1) / 40
            self.centers = tuple(mp.mpf(value) for value in CENTERS)
            self.nodes, self.weights = mp.gauss_quadrature(self.order, "legendre")
            self.beta_normalization = self._integral_raw(
                _raw_bump, mp.mpf(-1), mp.mpf(1), breaks=(mp.mpf(0),)
            )
            if self.beta_normalization <= 0:
                raise ArithmeticError("paper bump normalization must be positive")
            self._support_edges = tuple(
                edge
                for center in self.centers
                for edge in (center - self.radius, center + self.radius)
            )
            self.linear_matrix = self._build_linear_matrix(mp.mpf(1), mp.mpf(2))
            self.matrix = self.linear_matrix
            self.inverse_matrix = self._scalar_inverse(self.linear_matrix)
            self.inverse = self.inverse_matrix
            self.quadratic_weights = self._build_quadratic_weights(
                mp.mpf(1), mp.mpf(2)
            )
            # Explicit aliases make the scalar receipts easy to consume from
            # bridge code without requiring callers to know the private names.
            self.fg_weights = self.quadratic_weights["fg_sqrt"]
            self.gg_weights = self.quadratic_weights["gg"]
            self.ff_weights = self.quadratic_weights["ff"]
            self.ff_over_x_weights = self.quadratic_weights["ff_over_x"]
            self.inverse_l1_norm = max(
                (
                    sum(
                        abs(self.inverse_matrix[row, col])
                        for row in range(self.inverse_matrix.rows)
                    )
                    for col in range(self.inverse_matrix.cols)
                ),
                default=mp.mpf(0),
            )
            self.metadata = self._metadata()

    # ------------------------------------------------------------------
    # Fixed bump and MP Gauss primitives.
    # ------------------------------------------------------------------
    def _integral_raw(
        self,
        function: Any,
        left: Any,
        right: Any,
        *,
        breaks: Iterable[Any] = (),
    ) -> mp.mpf:
        left_mp = mp.mpf(left)
        right_mp = mp.mpf(right)
        if right_mp <= left_mp:
            return mp.mpf(0)
        points = [left_mp]
        points.extend(
            point
            for point in (mp.mpf(item) for item in breaks)
            if left_mp < point < right_mp
        )
        points.append(right_mp)
        points = sorted(set(points))
        total = mp.mpf(0)
        for lower, upper in zip(points, points[1:]):
            midpoint = (lower + upper) / 2
            half = (upper - lower) / 2
            for node, weight in zip(self.nodes, self.weights):
                total += half * weight * function(midpoint + half * node)
        return total

    def _interval_breaks(self, left: mp.mpf, right: mp.mpf) -> tuple[mp.mpf, ...]:
        return tuple(
            sorted(
                {
                    edge
                    for edge in self._support_edges
                    if left < edge < right
                }
            )
        )

    def _integral(self, function: Any, left: Any, right: Any) -> mp.mpf:
        left_mp = mp.mpf(left)
        right_mp = mp.mpf(right)
        return self._integral_raw(
            function,
            left_mp,
            right_mp,
            breaks=self._interval_breaks(left_mp, right_mp),
        )

    def beta(self, x: Any, center: Any) -> mp.mpf:
        """Return ``beta_{1/40}(x-center)`` with exact MP support logic."""

        t = (mp.mpf(x) - mp.mpf(center)) / self.radius
        if abs(t) >= 1:
            return mp.mpf(0)
        return _raw_bump(t) / (self.radius * self.beta_normalization)

    def _gamma(self, index: int, x: Any) -> mp.mpf:
        return self.beta(x, self.centers[index])

    def _mass_integral(self, index: int, left: Any, right: Any) -> mp.mpf:
        left_mp = mp.mpf(left)
        right_mp = mp.mpf(right)
        center = self.centers[index]
        support_left = center - self.radius
        support_right = center + self.radius
        if right_mp <= support_left or left_mp >= support_right:
            return mp.mpf(0)
        if left_mp <= support_left and right_mp >= support_right:
            # This is the defining normalization of beta_{1/40}; use the
            # identity exactly on complete supports rather than manufacture a
            # quadrature residual in row one.
            return mp.mpf(1)
        lower = max(left_mp, support_left)
        upper = min(right_mp, support_right)
        return self._integral(lambda x: self._gamma(index, x), lower, upper)

    def _weighted_gamma(
        self,
        index: int,
        power: mp.mpf,
        left: Any,
        right: Any,
    ) -> mp.mpf:
        return self._integral(
            lambda x: x**power * self._gamma(index, x), left, right
        )

    def _product_integral(
        self,
        left_index: int,
        right_index: int,
        weight_power: mp.mpf,
        left: Any,
        right: Any,
    ) -> mp.mpf:
        return self._integral(
            lambda x: x**weight_power
            * self._gamma(left_index, x)
            * self._gamma(right_index, x),
            left,
            right,
        )

    def _build_linear_matrix(self, left: Any, right: Any) -> mp.matrix:
        matrix = mp.matrix(5, 5)
        for row in range(5):
            for column in range(5):
                matrix[row, column] = mp.mpf(0)
        # b1=gamma1 and b2=gamma3, hence the axial columns are gamma 1 and 3.
        for axial_column, gamma_index in enumerate((0, 2)):
            matrix[0, axial_column] = self._mass_integral(
                gamma_index, left, right
            )
            matrix[1, axial_column] = self._weighted_gamma(
                gamma_index, mp.mpf(".6"), left, right
            )
        for angular_column in range(3):
            gamma_index = angular_column
            column = 2 + angular_column
            matrix[2, column] = self._weighted_gamma(
                gamma_index, mp.mpf(".5"), left, right
            )
            matrix[3, column] = -self._weighted_gamma(
                gamma_index, mp.mpf(".1"), left, right
            )
            matrix[4, column] = self._weighted_gamma(
                gamma_index, mp.mpf("-.9"), left, right
            )
        return matrix

    def _build_quadratic_weights(self, left: Any, right: Any) -> dict[str, mp.matrix]:
        fg = mp.matrix(2, 3)
        gg = mp.matrix(2, 2)
        ff = mp.matrix(3, 3)
        ff_over_x = mp.matrix(3, 3)
        for i in range(2):
            gamma_i = (0, 2)[i]
            for j in range(3):
                fg[i, j] = self._product_integral(
                    gamma_i, j, mp.mpf(".5"), left, right
                )
            for j in range(2):
                gamma_j = (0, 2)[j]
                gg[i, j] = self._product_integral(
                    gamma_i, gamma_j, mp.mpf(0), left, right
                )
        for i in range(3):
            for j in range(3):
                ff[i, j] = self._product_integral(
                    i, j, mp.mpf(0), left, right
                )
                ff_over_x[i, j] = self._product_integral(
                    i, j, mp.mpf(-1), left, right
                )
        return {
            "fg_sqrt": fg,
            "gg": gg,
            "ff": ff,
            "ff_over_x": ff_over_x,
        }

    @staticmethod
    def _scalar_inverse(matrix: mp.matrix) -> mp.matrix:
        """Invert only the fixed MP scalar matrix, never a coefficient ring."""

        return mp.inverse(matrix)

    def _metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_section": SOURCE_SECTION,
            "coefficient_order": self.coefficient_order,
            "centers": tuple(self.centers),
            "radius": self.radius,
            "beta_normalization": self.beta_normalization,
            "quadrature_order": self.order,
            "precision": self.precision,
            "linear_matrix": _mp_matrix_to_lists(self.linear_matrix),
            "inverse_matrix": _mp_matrix_to_lists(self.inverse_matrix),
            "quadratic_weights": {
                name: _mp_matrix_to_lists(value)
                for name, value in self.quadratic_weights.items()
            },
            "finite_quadrature": True,
            "quadrature_enclosed": False,
            "functional_closure": False,
            "cone_certified": False,
        }

    # ------------------------------------------------------------------
    # Generic coefficient map.
    # ------------------------------------------------------------------
    def _am_inverse_square(self, Am: Any) -> Any:
        try:
            return Am**(-2)
        except (TypeError, ValueError, ZeroDivisionError):
            return 1 / (Am * Am)

    @staticmethod
    def _bilinear(left: Sequence[Any], matrix: mp.matrix, right: Sequence[Any]) -> Any:
        template = left[0] if left else (right[0] if right else mp.mpf(0))
        return _sum_terms(
            (
                left[row] * matrix[row, column] * right[column]
                for row in range(matrix.rows)
                for column in range(matrix.cols)
            ),
            template,
        )

    def _linear_apply(self, h: Sequence[Any], matrix: mp.matrix | None = None) -> list[Any]:
        matrix = self.linear_matrix if matrix is None else matrix
        result: list[Any] = []
        for row in range(matrix.rows):
            result.append(
                _sum_terms(
                    (h[column] * matrix[row, column] for column in range(matrix.cols)),
                    h[0],
                )
            )
        return result

    def _quadratic_with_weights(
        self,
        h: Sequence[Any],
        Am: Any,
        weights: Mapping[str, mp.matrix],
    ) -> list[Any]:
        c = h[:2]
        xi = h[2:]
        fg = self._bilinear(c, weights["fg_sqrt"], xi)
        gg = self._bilinear(c, weights["gg"], c)
        ff = self._bilinear(xi, weights["ff"], xi)
        ff_over_x = self._bilinear(xi, weights["ff_over_x"], xi)
        zero = _zero_like(h[0])
        return [zero, fg, zero, self._am_inverse_square(Am) * gg - ff / 2, ff_over_x / 2]

    def apply(self, h: Any, Am: Any) -> tuple[Any, ...]:
        """Evaluate all five rows for scalar, jet, or dual coefficients."""

        vector = _as_sequence(h, "h")
        with mp.workdps(self.precision):
            linear = self._linear_apply(vector)
            quadratic = self._quadratic_with_weights(
                vector, Am, self.quadratic_weights
            )
            return tuple(linear[row] + quadratic[row] for row in range(5))

    evaluate = apply
    __call__ = apply

    def quadratic(self, h: Any, Am: Any) -> tuple[Any, ...]:
        """Return only the nonlinear five-row contribution."""

        vector = _as_sequence(h, "h")
        with mp.workdps(self.precision):
            return self.bilinear(vector, vector, Am)

    def bilinear(self, h: Any, k: Any, Am: Any) -> tuple[Any, ...]:
        """Return the symmetric quadratic polarization ``Q(h,k)``.

        The row-two cross term is explicitly symmetrized as
        ``(f_h*g_k + f_k*g_h)/2``.  Rows four and five use the direct
        symmetric products specified by (10.8), so formal response code can
        evaluate mixed coefficient directions without subtracting two large
        generic-ring values.
        """

        left = _as_sequence(h, "h")
        right = _as_sequence(k, "k")
        c_left, xi_left = left[:2], left[2:]
        c_right, xi_right = right[:2], right[2:]
        weights = self.quadratic_weights
        fg_left_right = self._bilinear(
            c_left, weights["fg_sqrt"], xi_right
        )
        fg_right_left = self._bilinear(
            c_right, weights["fg_sqrt"], xi_left
        )
        gg = self._bilinear(c_left, weights["gg"], c_right)
        ff = self._bilinear(xi_left, weights["ff"], xi_right)
        ff_over_x = self._bilinear(
            xi_left, weights["ff_over_x"], xi_right
        )
        template = left[0]
        return (
            _zero_like(template),
            (fg_left_right + fg_right_left) / 2,
            _zero_like(template),
            self._am_inverse_square(Am) * gg - ff / 2,
            ff_over_x / 2,
        )

    def jacobian(self, h: Any, Am: Any) -> tuple[tuple[Any, ...], ...]:
        """Return the generic coefficient Jacobian of ``apply``."""

        vector = _as_sequence(h, "h")
        c = vector[:2]
        xi = vector[2:]
        weights = self.quadratic_weights
        result: list[list[Any]] = [
            [self.linear_matrix[row, col] for col in range(5)]
            for row in range(5)
        ]
        inv_am2 = self._am_inverse_square(Am)
        for i in range(2):
            for j in range(3):
                result[1][2 + j] = result[1][2 + j] + weights["fg_sqrt"][i, j] * c[i]
                result[1][i] = result[1][i] + weights["fg_sqrt"][i, j] * xi[j]
        for i in range(2):
            result[3][i] = result[3][i] + 2 * inv_am2 * _sum_terms(
                (weights["gg"][i, k] * c[k] for k in range(2)), c[0]
            )
        for j in range(3):
            result[3][2 + j] = result[3][2 + j] - _sum_terms(
                (weights["ff"][j, k] * xi[k] for k in range(3)), xi[0]
            )
            result[4][2 + j] = result[4][2 + j] + _sum_terms(
                (weights["ff_over_x"][j, k] * xi[k] for k in range(3)), xi[0]
            )
        return tuple(tuple(row) for row in result)

    def linear_inverse(self, rhs: Any) -> tuple[Any, ...]:
        """Apply the precomputed scalar inverse to a generic five-vector."""

        vector = _as_sequence(rhs, "rhs")
        with mp.workdps(self.precision):
            result: list[Any] = []
            for row in range(5):
                result.append(
                    _sum_terms(
                        (
                            vector[column] * self.inverse_matrix[row, column]
                            for column in range(5)
                        ),
                        vector[0],
                    )
                )
            return tuple(result)

    def _weights_for_interval(self, left: Any, right: Any) -> tuple[mp.matrix, dict[str, mp.matrix]]:
        left_mp = mp.mpf(left)
        right_mp = mp.mpf(right)
        if not mp.isfinite(left_mp) or not mp.isfinite(right_mp):
            raise ValueError("partial interval endpoints must be finite")
        if left_mp < 1 or right_mp > 2:
            raise ValueError("partial interval must lie in 1 <= x <= 2")
        if right_mp < left_mp:
            raise ValueError("partial interval requires left <= right")
        return (
            self._build_linear_matrix(left_mp, right_mp),
            self._build_quadratic_weights(left_mp, right_mp),
        )

    def partial(self, h: Any, Am: Any, x1: Any = 1, x2: Any = 2) -> tuple[Any, ...]:
        """Evaluate the same map over an arbitrary subinterval of ``[1,2]``."""

        vector = _as_sequence(h, "h")
        with mp.workdps(self.precision):
            linear_matrix, quadratic_weights = self._weights_for_interval(x1, x2)
            linear = self._linear_apply(vector, linear_matrix)
            quadratic = self._quadratic_with_weights(vector, Am, quadratic_weights)
            return tuple(linear[row] + quadratic[row] for row in range(5))

    partial_increment = partial
    apply_partial = partial

    def partial_metadata(self, x1: Any = 1, x2: Any = 2) -> dict[str, Any]:
        """Return scalar receipts for a partial interval."""

        with mp.workdps(self.precision):
            matrix, weights = self._weights_for_interval(x1, x2)
            return {
                "x1": mp.mpf(x1),
                "x2": mp.mpf(x2),
                "linear_matrix": _mp_matrix_to_lists(matrix),
                "quadratic_weights": {
                    key: _mp_matrix_to_lists(value)
                    for key, value in weights.items()
                },
                "finite_quadrature": True,
                "quadrature_enclosed": False,
            }


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 50)
    if isinstance(value, mp.matrix):
        return _mp_matrix_to_lists(value)
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


__all__ = [
    "FiveBumpMomentMap",
    "COEFFICIENT_NAMES",
    "MOMENT_ROWS",
    "SOURCE",
    "SOURCE_SECTION",
]
