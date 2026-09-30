"""Section 10.2 fixed five-moment bump map.

The map is the exact finite-dimensional correction used on ``1 < x < 2``
with the source bump ``beta_{1/40}``.  The coefficient order is

    h = (c1, c2, xi1, xi2, xi3),

where ``c1,c2`` multiply the first and last axial bumps and ``xi1..xi3``
multiply the three angular bumps.  The nonlinear rows are evaluated from
the integrated products, never replaced by zero or by a target moment.

``solve(d, Am)`` solves ``evaluate(h, Am) + d = 0`` for an explicit numeric
five-vector ``d``.  ``solve_z_derivative`` uses the analytic implicit equation

    J h_Z = -d_Z - F_{Am} Am_Z.

This module is a numerical map and quadrature receipt.  It does not infer
unknown defects, certify the input contraction ball, or claim the full outer
construction.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Part I Section 10.2, equations (10.7)-(10.10)"
MOMENT_KEYS = ("z", "theta_z_centered", "theta", "z_theta_centered", "p")
COEFFICIENT_NAMES = ("c1", "c2", "xi1", "xi2", "xi3")
RADIUS = '.025'
CENTERS = ('1.25', '1.5', '1.75')


def _mp(value: Any) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(str(value))


def _vector(values: Any, name: str = "vector") -> mp.matrix:
    if isinstance(values, mp.matrix):
        if values.rows * values.cols != 5:
            raise ValueError(f"{name} must have exactly five entries")
        flat = [values[i] for i in range(values.rows * values.cols)]
    else:
        try:
            flat = list(values)
        except TypeError as exc:
            raise TypeError(f"{name} must be an iterable of five numeric values") from exc
        if len(flat) != 5:
            raise ValueError(f"{name} must have exactly five entries")
    if any(value is None for value in flat):
        raise ValueError(f"{name} must provide all five entries explicitly")
    return mp.matrix([_mp(value) for value in flat])


def _nstr(value: Any, digits: int = 50) -> str:
    return mp.nstr(_mp(value), digits)


def _max_abs(values: Iterable[Any]) -> mp.mpf:
    values = list(values)
    return max((abs(_mp(value)) for value in values), default=mp.mpf(0))


def _phi(t: mp.mpf) -> mp.mpf:
    if abs(t) >= 1:
        return mp.mpf(0)
    return mp.exp(-1 / (1 - t * t))


def _phi_derivative(t: mp.mpf, order: int) -> mp.mpf:
    if abs(t) >= 1:
        return mp.mpf(0)
    one_minus = 1 - t * t
    value = _phi(t)
    if order == 0:
        return value
    first = -2 * t / one_minus ** 2
    if order == 1:
        return value * first
    if order == 2:
        second = -2 / one_minus ** 2 - 8 * t * t / one_minus ** 3
        return value * (first * first + second)
    raise ValueError("bump derivative order must be 0, 1, or 2")


class MomentMap:
    """MP quadrature representation of the exact Section 10.2 map."""

    def __init__(self, *, precision: int = 160, order: int = 64) -> None:
        if int(precision) != precision or precision < 80:
            raise ValueError("precision must be an integer >= 80")
        if int(order) != order or order < 8:
            raise ValueError("order must be an integer >= 8")
        self.precision = int(precision)
        self.order = int(order)
        with mp.workdps(self.precision):
            self.radius = mp.mpf(RADIUS)
            self.centers = tuple(mp.mpf(center) for center in CENTERS)
            self.beta_normalization = mp.quad(_phi, [-1, 0, 1])
            if self.beta_normalization <= 0:
                raise ArithmeticError("bump normalization must be positive")
            self.nodes, self.weights = mp.gauss_quadrature(self.order, "legendre")
            self._support_edges = tuple(
                edge
                for center in self.centers
                for edge in (center - self.radius, center + self.radius)
            )
            self.matrix = self._build_matrix(1, mp.mpf(2))
            self.inverse = mp.inverse(self.matrix)
            self.inverse_l1_norm = max(
                sum(abs(self.inverse[row, col]) for row in range(5))
                for col in range(5)
            )
            self._fg = mp.matrix(2, 3)
            self._gg = mp.matrix(2, 2)
            self._ff = mp.matrix(3, 3)
            self._ff_over_x = mp.matrix(3, 3)
            for i in range(2):
                for j in range(3):
                    self._fg[i, j] = self._integral(
                        lambda x, i=i, j=j: mp.sqrt(x)
                        * self._basis_value(i, x)
                        * self._basis_value(2 + j, x),
                        1,
                        2,
                    )
            for i in range(2):
                for j in range(2):
                    self._gg[i, j] = self._integral(
                        lambda x, i=i, j=j: self._basis_value(i, x)
                        * self._basis_value(j, x),
                        1,
                        2,
                    )
            for i in range(3):
                for j in range(3):
                    product = lambda x, i=i, j=j: self._basis_value(2 + i, x) * self._basis_value(2 + j, x)
                    self._ff[i, j] = self._integral(product, 1, 2)
                    self._ff_over_x[i, j] = self._integral(
                        lambda x, product=product: product(x) / x,
                        1,
                        2,
                    )

    # ------------------------------------------------------------------
    # Fixed bump and quadrature primitives.
    # ------------------------------------------------------------------
    def _basis_value(self, index: int, x: Any, derivative: int = 0) -> mp.mpf:
        if not 0 <= int(index) < 5:
            raise IndexError("basis index must be 0..4")
        x_mp = _mp(x)
        bump_index = index if index < 2 else index - 2
        center = self.centers[0 if index == 0 else 2 if index == 1 else bump_index]
        t = (x_mp - center) / self.radius
        return _phi_derivative(t, derivative) / (
            self.radius ** (derivative + 1) * self.beta_normalization
        )

    def bump(self, index: int, x: Any, derivative: int = 0) -> mp.mpf:
        """Return bump ``gamma_{index+1}`` or axial bump value/derivative.

        Indices 0 and 1 are the axial bumps ``b1=gamma1`` and
        ``b2=gamma3``; indices 2, 3, and 4 are the angular bumps
        ``gamma1, gamma2, gamma3``.
        """

        return self._basis_value(index, x, derivative)

    def basis_values(self, x: Any, derivative: int = 0) -> tuple[mp.mpf, ...]:
        """Return all five coefficient-basis values or radial derivatives."""

        return tuple(self._basis_value(index, x, derivative) for index in range(5))

    def bump_values(self, x: Any) -> dict[str, tuple[mp.mpf, ...]]:
        """Return angular gamma values and first radial derivatives.

        The returned keys are intentionally field-adapter friendly:
        ``gamma`` has entries ``gamma1,gamma2,gamma3`` and ``gamma_x`` has
        their exact radial derivatives.  Axial values are available through
        :meth:`basis_values` at indices 0 and 1.
        """

        return {
            "gamma": tuple(self._basis_value(2 + index, x, 0) for index in range(3)),
            "gamma_x": tuple(self._basis_value(2 + index, x, 1) for index in range(3)),
        }

    def bump_derivatives(self, x: Any) -> tuple[mp.mpf, ...]:
        return tuple(self._basis_value(index, x, 1) for index in range(5))

    def _interval_breaks(self, left: mp.mpf, right: mp.mpf) -> list[mp.mpf]:
        points = [left]
        points.extend(edge for edge in self._support_edges if left < edge < right)
        points.append(right)
        return sorted(set(points))

    def _integral(self, function: Any, left: Any, right: Any) -> mp.mpf:
        left_mp = _mp(left)
        right_mp = _mp(right)
        if right_mp <= left_mp:
            return mp.mpf(0)
        total = mp.mpf(0)
        points = self._interval_breaks(left_mp, right_mp)
        for lower, upper in zip(points, points[1:]):
            midpoint = (lower + upper) / 2
            half = (upper - lower) / 2
            for node, weight in zip(self.nodes, self.weights):
                total += half * weight * function(midpoint + half * node)
        return total

    def _mass_integral(self, index: int, left: Any, right: Any) -> mp.mpf:
        """Integrate one axial bump, using its normalized mass exactly.

        The source normalization is ``integral beta = 1``.  A completed
        support therefore contributes exactly one coefficient; only an
        interval cut through a support needs quadrature.  This keeps the
        small residual axial mass from being manufactured by quadrature
        error after a bump has ended.
        """

        if index not in (0, 1):
            raise IndexError("axial mass index must be 0 or 1")
        left_mp = _mp(left)
        right_mp = _mp(right)
        if right_mp <= left_mp:
            return mp.mpf(0)
        center = self.centers[0 if index == 0 else 2]
        support_left = center - self.radius
        support_right = center + self.radius
        if right_mp <= support_left or left_mp >= support_right:
            return mp.mpf(0)
        if left_mp <= support_left and right_mp >= support_right:
            return mp.mpf(1)
        lower = max(left_mp, support_left)
        upper = min(right_mp, support_right)
        return self._integral(lambda x: self._basis_value(index, x), lower, upper)

    def _build_matrix(self, left: Any, right: Any) -> mp.matrix:
        matrix = mp.matrix(5, 5)
        for i in range(2):
            matrix[0, i] = self._mass_integral(i, left, right)
            matrix[1, i] = self._integral(
                lambda x, i=i: x ** mp.mpf(".6") * self._basis_value(i, x), left, right
            )
        for j in range(3):
            matrix[2, 2 + j] = self._integral(
                lambda x, j=j: mp.sqrt(x) * self._basis_value(2 + j, x), left, right
            )
            matrix[3, 2 + j] = -self._integral(
                lambda x, j=j: x ** mp.mpf(".1") * self._basis_value(2 + j, x), left, right
            )
            matrix[4, 2 + j] = self._integral(
                lambda x, j=j: x ** mp.mpf("-.9") * self._basis_value(2 + j, x), left, right
            )
        return matrix

    def _split_coefficients(self, h: Any) -> tuple[mp.matrix, mp.matrix]:
        vector = _vector(h, "h")
        return mp.matrix([vector[0], vector[1]]), mp.matrix([vector[2], vector[3], vector[4]])

    def _validate_Am(self, Am: Any) -> mp.mpf:
        value = _mp(Am)
        if not mp.isfinite(value) or value == 0:
            raise ValueError("Am must be finite and nonzero")
        return value

    # ------------------------------------------------------------------
    # Full and partial maps.
    # ------------------------------------------------------------------
    def _map_from_tensors(self, h: mp.matrix, Am: mp.mpf) -> mp.matrix:
        c, xi = self._split_coefficients(h)
        out = self.matrix * h
        fg = (c.T * self._fg * xi)[0]
        gg = (c.T * self._gg * c)[0]
        ff = (xi.T * self._ff * xi)[0]
        ff_x = (xi.T * self._ff_over_x * xi)[0]
        out[1] += fg
        out[3] += gg / (Am * Am) - ff / 2
        out[4] += ff_x / 2
        return out

    def evaluate(self, h: Any, Am: Any) -> mp.matrix:
        """Return the five centered moment changes in source order."""

        with mp.workdps(self.precision):
            return self._map_from_tensors(_vector(h, "h"), self._validate_Am(Am))

    def jacobian(self, h: Any, Am: Any) -> mp.matrix:
        with mp.workdps(self.precision):
            vector = _vector(h, "h")
            Am_mp = self._validate_Am(Am)
            c, xi = self._split_coefficients(vector)
            result = mp.matrix(self.matrix)
            for i in range(2):
                result[1, i] += sum(self._fg[i, j] * xi[j] for j in range(3))
                result[3, i] += 2 * sum(self._gg[i, k] * c[k] for k in range(2)) / (Am_mp * Am_mp)
            for j in range(3):
                result[1, 2 + j] += sum(self._fg[i, j] * c[i] for i in range(2))
                result[3, 2 + j] -= sum(self._ff[j, k] * xi[k] for k in range(3))
                result[4, 2 + j] += sum(self._ff_over_x[j, k] * xi[k] for k in range(3))
            return result

    def parameter_derivative(self, h: Any, Am: Any, Am_Z: Any) -> mp.matrix:
        """Return F_Am * Am_Z for the implicit Z derivative equation."""

        with mp.workdps(self.precision):
            vector = _vector(h, "h")
            Am_mp = self._validate_Am(Am)
            c, _ = self._split_coefficients(vector)
            out = mp.matrix(5, 1)
            out[3] = -2 * _mp(Am_Z) * (c.T * self._gg * c)[0] / (Am_mp ** 3)
            return out

    def quadratic_l1_bound(self, Am: Any) -> mp.mpf:
        """Conservative coefficient bound for the integrated quadratic map."""

        Am_mp = abs(self._validate_Am(Am))
        q2 = _max_abs(self._fg) / 2
        q4 = max(_max_abs(self._gg) / (Am_mp * Am_mp), _max_abs(self._ff) / 2)
        q5 = _max_abs(self._ff_over_x) / 2
        return q2 + q4 + q5

    def solve(self, d: Any, Am: Any, *, initial: Any | None = None, max_iter: int = 32) -> dict[str, Any]:
        """Solve ``evaluate(h, Am) + d = 0`` for an explicit five-vector d."""

        with mp.workdps(self.precision):
            defect = _vector(d, "d")
            Am_mp = self._validate_Am(Am)
            if initial is None:
                h = -self.inverse * defect
            else:
                h = _vector(initial, "initial")
            converged = False
            iterations = 0
            for iterations in range(1, int(max_iter) + 1):
                residual = self.evaluate(h, Am_mp) + defect
                scale = max(mp.mpf(1), max(abs(residual[i]) for i in range(5)))
                if max(abs(residual[i]) for i in range(5)) <= mp.power(10, -self.precision // 2) * scale:
                    converged = True
                    break
                h += mp.lu_solve(self.jacobian(h, Am_mp), -residual)
            residual = self.evaluate(h, Am_mp) + defect
            if not converged:
                raise ArithmeticError("MomentMap Newton solve did not converge")
            return {
                "h": h,
                "coefficients": dict(zip(COEFFICIENT_NAMES, [h[i] for i in range(5)])),
                "residual": residual,
                "jacobian": self.jacobian(h, Am_mp),
                "iterations": iterations,
                "converged": converged,
                "quadratic_l1_bound": self.quadratic_l1_bound(Am_mp),
            }

    def solve_z_derivative(
        self,
        d: Any,
        d_Z: Any,
        Am: Any,
        Am_Z: Any,
        *,
        h: Any | None = None,
    ) -> dict[str, Any]:
        """Solve ``J h_Z = -d_Z - F_Am Am_Z`` analytically."""

        with mp.workdps(self.precision):
            solution = self.solve(d, Am) if h is None else None
            h_value = solution["h"] if solution is not None else _vector(h, "h")
            defect_z = _vector(d_Z, "d_Z")
            Am_mp = self._validate_Am(Am)
            parameter = self.parameter_derivative(h_value, Am_mp, Am_Z)
            jac = self.jacobian(h_value, Am_mp)
            h_z = mp.lu_solve(jac, -defect_z - parameter)
            return {
                "h_Z": h_z,
                "coefficients_Z": dict(zip(COEFFICIENT_NAMES, [h_z[i] for i in range(5)])),
                "jacobian": jac,
                "parameter_derivative": parameter,
                "equation_residual": jac * h_z + defect_z + parameter,
                "h": h_value,
            }

    def solve_Z(self, h: Any, d_Z: Any, Am: Any, Am_Z: Any) -> mp.matrix:
        """Short adapter API when the solved ``h`` is already available."""

        with mp.workdps(self.precision):
            h_value = _vector(h, "h")
            defect_z = _vector(d_Z, "d_Z")
            Am_mp = self._validate_Am(Am)
            parameter = self.parameter_derivative(h_value, Am_mp, Am_Z)
            return mp.lu_solve(self.jacobian(h_value, Am_mp), -defect_z - parameter)

    def _direct_partial(self, x: Any, h: Any, Am: Any, *, jacobian: bool = False) -> tuple[mp.matrix, mp.matrix | None]:
        x_mp = _mp(x)
        if not 1 <= x_mp <= 2:
            raise ValueError("partial integrals require 1 <= x <= 2")
        vector = _vector(h, "h")
        Am_mp = self._validate_Am(Am)
        c, xi = self._split_coefficients(vector)

        def f_value(s):
            return sum(xi[j] * self._basis_value(2 + j, s) for j in range(3))

        def g_value(s):
            return sum(c[i] * self._basis_value(i, s) for i in range(2))

        out = mp.matrix(5, 1)
        out[0] = c[0] * self._mass_integral(0, 1, x_mp) + c[1] * self._mass_integral(1, 1, x_mp)
        out[1] = self._integral(lambda s: s ** mp.mpf(".6") * g_value(s) + mp.sqrt(s) * f_value(s) * g_value(s), 1, x_mp)
        out[2] = self._integral(lambda s: mp.sqrt(s) * f_value(s), 1, x_mp)
        out[3] = self._integral(lambda s: -s ** mp.mpf(".1") * f_value(s) + Am_mp ** -2 * g_value(s) ** 2 - f_value(s) ** 2 / 2, 1, x_mp)
        out[4] = self._integral(lambda s: s ** mp.mpf("-.9") * f_value(s) + f_value(s) ** 2 / (2 * s), 1, x_mp)
        if not jacobian:
            return out, None
        jac = self._partial_jacobian(x_mp, vector, Am_mp, f_value=f_value, g_value=g_value)
        return out, jac

    def _partial_jacobian(self, x: mp.mpf, h: mp.matrix, Am: mp.mpf, *, f_value: Any, g_value: Any) -> mp.matrix:
        jac = self._build_matrix(1, x)
        for i in range(2):
            jac[1, i] += self._integral(
                lambda s, i=i: mp.sqrt(s) * f_value(s) * self._basis_value(i, s), 1, x
            )
            c_i = h[i]
            jac[3, i] += self._integral(
                lambda s, i=i: 2 * Am ** -2 * g_value(s) * self._basis_value(i, s), 1, x
            )
        for j in range(3):
            jac[1, 2 + j] += self._integral(
                lambda s, j=j: mp.sqrt(s) * g_value(s) * self._basis_value(2 + j, s), 1, x
            )
            jac[3, 2 + j] -= self._integral(
                lambda s, j=j: f_value(s) * self._basis_value(2 + j, s), 1, x
            )
            jac[4, 2 + j] += self._integral(
                lambda s, j=j: f_value(s) * self._basis_value(2 + j, s) / s, 1, x
            )
        return jac

    def partial_integrals(self, x: Any, h: Any, Am: Any) -> mp.matrix:
        """Return the five centered moment changes accumulated over ``1<x'<x``."""

        with mp.workdps(self.precision):
            return self._direct_partial(x, h, Am)[0]

    def partial_Z(self, x: Any, h: Any, h_Z: Any, Am: Any, Am_Z: Any) -> mp.matrix:
        """Return the analytic Z derivative of ``partial_integrals``."""

        with mp.workdps(self.precision):
            value, jac = self._direct_partial(x, h, Am, jacobian=True)
            parameter = self._partial_parameter_derivative(x, h, Am, Am_Z)
            return jac * _vector(h_Z, "h_Z") + parameter

    def _partial_parameter_derivative(self, x: Any, h: Any, Am: Any, Am_Z: Any) -> mp.matrix:
        x_mp = _mp(x)
        vector = _vector(h, "h")
        Am_mp = self._validate_Am(Am)
        c, _ = self._split_coefficients(vector)
        g_value = lambda s: sum(c[i] * self._basis_value(i, s) for i in range(2))
        out = mp.matrix(5, 1)
        out[3] = self._integral(
            lambda s: -2 * _mp(Am_Z) * Am_mp ** -3 * g_value(s) ** 2,
            1,
            x_mp,
        )
        return out

    def partial_details(self, x: Any, h: Any, Am: Any) -> dict[str, Any]:
        """Return partial map plus raw square/cross primitives for field assembly."""

        x_mp = _mp(x)
        vector = _vector(h, "h")
        c, xi = self._split_coefficients(vector)
        g_value = lambda s: sum(c[i] * self._basis_value(i, s) for i in range(2))
        f_value = lambda s: sum(xi[j] * self._basis_value(2 + j, s) for j in range(3))
        value = self.partial_integrals(x_mp, vector, Am)
        return {
            "map": value,
            "axial_square": self._integral(lambda s: g_value(s) ** 2, 1, x_mp),
            "angular_square": self._integral(lambda s: f_value(s) ** 2, 1, x_mp),
            "angular_square_over_x": self._integral(lambda s: f_value(s) ** 2 / s, 1, x_mp),
            "angular_axial_cross": self._integral(lambda s: f_value(s) * g_value(s), 1, x_mp),
        }

    def axial_square_increment(self, x: Any, h: Any) -> mp.mpf:
        return self.partial_details(x, h, 1)["axial_square"]

    def axial_square_increment_Z(self, x: Any, h: Any, h_Z: Any) -> mp.mpf:
        x_mp = _mp(x)
        vector = _vector(h, "h")
        vector_z = _vector(h_Z, "h_Z")
        c, _ = self._split_coefficients(vector)
        c_z, _ = self._split_coefficients(vector_z)
        g = lambda s: sum(c[i] * self._basis_value(i, s) for i in range(2))
        g_z = lambda s: sum(c_z[i] * self._basis_value(i, s) for i in range(2))
        return 2 * self._integral(lambda s: g(s) * g_z(s), 1, x_mp)


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return _nstr(value)
    if isinstance(value, mp.matrix):
        return [
            [_json_value(value[row, col]) for col in range(value.cols)]
            for row in range(value.rows)
        ]
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def run_synthetic() -> dict[str, Any]:
    with mp.workdps(160):
        coarse = MomentMap(precision=160, order=32)
        fine = MomentMap(precision=160, order=64)
        Am0 = mp.mpf("1.7")
        AmZ = mp.mpf(".07")
        # The source contraction is a small-defect regime.  Keep every entry
        # nonzero while staying inside a numerical Newton basin.
        d0 = mp.matrix([mp.mpf("1e-8"), mp.mpf("-2e-8"), mp.mpf("3e-8"), mp.mpf("-4e-8"), mp.mpf("5e-8")])
        dZ = mp.matrix([mp.mpf("-3e-9"), mp.mpf("4e-9"), mp.mpf("-5e-9"), mp.mpf("6e-9"), mp.mpf("-7e-9")])
        solution = fine.solve(d0, Am0)
        derivative = fine.solve_z_derivative(d0, dZ, Am0, AmZ, h=solution["h"])
        z_step = mp.mpf("1e-5")
        plus = fine.solve(d0 + z_step * dZ, Am0 + z_step * AmZ)["h"]
        minus = fine.solve(d0 - z_step * dZ, Am0 - z_step * AmZ)["h"]
        fd = (plus - minus) / (2 * z_step)
        derivative_error = max(abs(fd[i] - derivative["h_Z"][i]) for i in range(5))
        map_diff = max(abs(coarse.matrix[i, j] - fine.matrix[i, j]) for i in range(5) for j in range(5))
        partial = fine.partial_details(mp.mpf("1.6"), solution["h"], Am0)
        partial_z = fine.partial_Z(mp.mpf("1.6"), solution["h"], derivative["h_Z"], Am0, AmZ)
        full_partial = fine.partial_integrals(mp.mpf("2"), solution["h"], Am0)
        full_partial_z = fine.partial_Z(mp.mpf("2"), solution["h"], derivative["h_Z"], Am0, AmZ)
        mass_identity_error = max(
            abs(fine.matrix[0, 0] - 1),
            abs(fine.matrix[0, 1] - 1),
            max(abs(fine.matrix[0, j]) for j in range(2, 5)),
        )
        partial_mass_identity_error = max(
            abs(partial["map"][0] - solution["h"][0]),
            abs(partial_z[0] - derivative["h_Z"][0]),
            abs(full_partial[0] - solution["h"][0] - solution["h"][1]),
            abs(full_partial_z[0] - derivative["h_Z"][0] - derivative["h_Z"][1]),
        )
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_section": SOURCE_SECTION,
            "radius": fine.radius,
            "centers": fine.centers,
            "beta_normalization": fine.beta_normalization,
            "coefficient_order": COEFFICIENT_NAMES,
            "matrix_order32_vs64_max_abs": map_diff,
            "inverse_l1_norm_order64": fine.inverse_l1_norm,
            "quadratic_l1_bound_Am1.7": fine.quadratic_l1_bound(Am0),
            "synthetic_defect": d0,
            "synthetic_defect_Z": dZ,
            "solution": solution,
            "Z_derivative": derivative,
            "Z_derivative_fd_max_abs_error": derivative_error,
            "exact_full_matrix_mass_row_error": mass_identity_error,
            "exact_completed_partial_mass_error": partial_mass_identity_error,
            "partial_map_x1.6": partial["map"],
            "partial_map_Z_x1.6": partial_z,
            "partial_details_x1.6": partial,
            "all_five_defects_explicit": True,
            "scope": "self-contained numerical map validation; no actual defect adoption or full cone/PDE certificate",
        }


if __name__ == "__main__":
    receipt = run_synthetic()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(_json_value(receipt), indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "matrix_order32_vs64_max_abs": _nstr(receipt["matrix_order32_vs64_max_abs"]),
                "inverse_l1_norm_order64": _nstr(receipt["inverse_l1_norm_order64"]),
                "Z_derivative_fd_max_abs_error": _nstr(receipt["Z_derivative_fd_max_abs_error"]),
                "all_five_defects_explicit": receipt["all_five_defects_explicit"],
            }
        ),
        flush=True,
    )


__all__ = ["MomentMap", "run_synthetic", "COEFFICIENT_NAMES", "MOMENT_KEYS"]
