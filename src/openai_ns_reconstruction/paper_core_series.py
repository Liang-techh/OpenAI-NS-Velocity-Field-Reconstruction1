"""Finite nonlinear near-axis radial series for the paper core equations.

This module is a numerical implementation of the triangular radial Taylor
recurrence obtained from Eqs. (4.9) and (4.13), with ``nu=1``.  It is useful
for evaluating a finite prefix of the formal profile while testing the
independent residual operator.  The prefix is not a converged fixed-point
solution, an outer matching construction, or a theorem about the paper field.

The coefficient convention is ordinary powers of ``X``::

    Phi(X, eta) = sum_n phi[n](eta) * X**n
    U(X, eta)   = sum_n u[n](eta)   * X**n
    Pi(X, eta)  = sum_n pi[n](eta)  * X**n.

The angular functions are stored on a Chebyshev--Lobatto grid.  Their eta
derivatives use the corresponding spectral differentiation matrix and values
at off-grid eta use barycentric polynomial interpolation.  Consequently all
outputs are finite floating-point approximations, including the eta
derivatives; no claim of convergence is made by this layer.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from .natural_axis import real_gradient
from .paper_core_reference import PaperCoreReference
from .profiles import LeadingProfile
from .velocity import leading_pressure_cartesian, leading_velocity_cartesian


_MAX_DEGREE = 16
_MAX_ETA_NODES = 513
_DOMAIN_X_FACTOR = 4.1


def _integer(value: Any, name: str, *, minimum: int, maximum: int) -> int:
    """Validate a bounded integer option without accepting booleans."""

    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if not minimum <= result <= maximum:
        raise ValueError(f"{name} must satisfy {minimum} <= {name} <= {maximum}")
    return result


def _finite_array(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class ChebyshevEtaGrid:
    """Chebyshev--Lobatto eta grid and finite interpolation operators.

    ``nodes`` counts the grid points, including both endpoints.  The nodes
    are stored in the conventional descending order ``cos(pi*j/(nodes-1))``.
    Values passed to :meth:`differentiate` or :meth:`interpolate` have their
    eta axis last, which makes the helpers work for one or many coefficient
    rows at once.
    """

    nodes: int = 129
    eta: np.ndarray = field(init=False, repr=False)
    differentiation_matrix: np.ndarray = field(init=False, repr=False)
    barycentric_weights: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        nodes = _integer(self.nodes, "nodes", minimum=3, maximum=_MAX_ETA_NODES)
        # A finite spectral derivative needs at least one interior node.  The
        # rounded cosine values remain in the closed physical eta interval.
        eta = np.cos(np.pi * np.arange(nodes, dtype=float) / float(nodes - 1))
        eta[0], eta[-1] = 1.0, -1.0

        # Barycentric interpolation weights are reciprocal node products.
        # The endpoint weights are half the interior weights, up to one common
        # nonzero scale which cancels in interpolation.
        weights = (-1.0) ** np.arange(nodes, dtype=float)
        weights[[0, -1]] *= 0.5

        # The standard Chebyshev--Lobatto differentiation matrix.  Construct
        # off-diagonal entries first and set the diagonal by the zero-row-sum
        # identity, avoiding a singular 0/0 calculation.
        c = (-1.0) ** np.arange(nodes, dtype=float)
        c[[0, -1]] *= 2.0
        differences = eta[:, None] - eta[None, :]
        matrix = (c[:, None] / c[None, :]) / (differences + np.eye(nodes))
        matrix = matrix - np.diag(np.sum(matrix, axis=1))
        matrix = np.asarray(matrix, dtype=float)

        object.__setattr__(self, "nodes", nodes)
        object.__setattr__(self, "eta", eta)
        object.__setattr__(self, "differentiation_matrix", matrix)
        object.__setattr__(self, "barycentric_weights", weights)

    @property
    def D(self) -> np.ndarray:
        """Compatibility alias for the eta differentiation matrix."""

        return self.differentiation_matrix

    def differentiate(self, values: Any) -> np.ndarray:
        """Differentiate values on this grid along their final axis."""

        array = _finite_array(values, "values")
        if array.ndim == 0 or array.shape[-1] != self.nodes:
            raise ValueError(
                f"values must have a final axis of length {self.nodes}"
            )
        return np.matmul(array, self.differentiation_matrix.T)

    def interpolate(self, values: Any, eta: Any) -> np.ndarray:
        """Evaluate the grid polynomial at scalar or array-valued eta.

        The final axis of ``values`` is interpreted as the grid axis.  The
        returned shape is ``values.shape[:-1] + eta.shape``; scalar eta keeps
        the coefficient-row dimensions and scalar rows therefore return a
        scalar.  Exact grid nodes are handled before forming reciprocal
        distances so endpoint evaluations stay well conditioned.
        """

        array = _finite_array(values, "values")
        if array.ndim == 0 or array.shape[-1] != self.nodes:
            raise ValueError(
                f"values must have a final axis of length {self.nodes}"
            )
        points = _finite_array(eta, "eta")
        if np.any(np.abs(points) > 1.0 + 2e-13):
            raise ValueError("eta must satisfy |eta| <= 1")
        points = np.clip(points, -1.0, 1.0)

        point_shape = points.shape
        flat_points = points.reshape(-1)
        flat_values = array.reshape((-1, self.nodes))
        output = np.empty((flat_values.shape[0], flat_points.size), dtype=float)

        # 1e-14 is below the useful resolution of this binary64 interpolation
        # and avoids an avoidable endpoint divide-by-zero for copied nodes.
        exact_tolerance = 2e-14
        for index, point in enumerate(flat_points):
            exact = np.flatnonzero(np.abs(self.eta - point) <= exact_tolerance)
            if exact.size:
                output[:, index] = flat_values[:, int(exact[0])]
                continue
            factors = self.barycentric_weights / (point - self.eta)
            output[:, index] = (flat_values @ factors) / np.sum(factors)

        return output.reshape(array.shape[:-1] + point_shape)

    def differentiate_and_interpolate(self, values: Any, eta: Any) -> np.ndarray:
        """Differentiate grid values and interpolate the derivative."""

        return self.interpolate(self.differentiate(values), eta)


def _coefficient_product(left: np.ndarray, right: np.ndarray, order: int) -> np.ndarray:
    """Return the coefficient ``order`` of two X-polynomial arrays."""

    result = np.zeros(left.shape[1], dtype=float)
    first = max(0, order - (right.shape[0] - 1))
    last = min(order, left.shape[0] - 1)
    for k in range(first, last + 1):
        result += left[k] * right[order - k]
    return result


def _polynomial_value(coefficients: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Horner evaluation for rows of X coefficients at broadcast X."""

    result = np.zeros_like(X, dtype=float)
    for coefficient in coefficients[::-1]:
        result = result * X + coefficient
    return result


def _polynomial_derivative_coefficients(coefficients: np.ndarray) -> np.ndarray:
    if coefficients.shape[0] <= 1:
        return np.zeros_like(coefficients[:1])
    factors = np.arange(1, coefficients.shape[0], dtype=float)[:, None]
    return factors * coefficients[1:]


@dataclass(frozen=True)
class PaperCoreSeries:
    """A finite, nonlinear radial Taylor prefix for the near-axis profile.

    Parameters are inherited from :class:`PaperCoreReference`; its pressure
    datum and amplitude ``g`` remain explicit independent choices in that
    reference object.  ``maxdegree`` retains coefficients through ``X`` raised
    to that degree, and ``eta_nodes`` controls the angular interpolation grid.
    Both are deliberately capped to keep this evaluator a finite diagnostic.
    """

    reference: PaperCoreReference = field(default_factory=PaperCoreReference)
    maxdegree: int = 8
    eta_nodes: int = 129

    grid: ChebyshevEtaGrid = field(init=False, repr=False)
    g_values: np.ndarray = field(init=False, repr=False)
    a_values: np.ndarray = field(init=False, repr=False)
    phi_coefficients: np.ndarray = field(init=False, repr=False)
    u_coefficients: np.ndarray = field(init=False, repr=False)
    pi_coefficients: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.reference, PaperCoreReference):
            raise TypeError("reference must be a PaperCoreReference")
        maxdegree = _integer(
            self.maxdegree, "maxdegree", minimum=0, maximum=_MAX_DEGREE
        )
        eta_nodes = _integer(
            self.eta_nodes, "eta_nodes", minimum=3, maximum=_MAX_ETA_NODES
        )
        grid = ChebyshevEtaGrid(eta_nodes)
        eta = grid.eta
        h, j, sigma, Lambda = (
            self.reference.h,
            self.reference.j,
            self.reference.sigma,
            self.reference.Lambda,
        )

        g = np.array([self.reference.amplitude(float(value)) for value in eta])
        # This is the requested a = g'/g identity, with realGradient supplied
        # by the existing axis construction.  Use the vectorized closed form
        # only as a consistency check source for the scalar helper's formula.
        a = np.array(
            [Lambda * real_gradient(h, j, sigma, float(value)) for value in eta]
        )
        if not np.all(np.isfinite(g)) or not np.all(np.isfinite(a)):
            raise ArithmeticError("axis amplitude/log derivative is non-finite")

        phi = np.zeros((maxdegree + 1, eta_nodes), dtype=float)
        u = np.zeros_like(phi)
        # The nonlinear pressure derivative is F**2.  A velocity prefix of
        # degree N therefore has an exactly integrable pressure prefix through
        # degree 2N+1.  Rows above N are filled after the triangular solve;
        # rows through N are populated during the solve because the next U
        # recurrence needs the already-known lower pressure rows.
        pi = np.zeros((2 * maxdegree + 2, eta_nodes), dtype=float)
        phi[0] = 1.0
        u[0] = 4.0 * eta + j
        pressure_scale = self.reference.pressure_scale
        pi[0] = -(pressure_scale**2) / (1.0 + eta**2) ** 2

        # Build the triangular recurrence.  At step n, only rows 0..n enter
        # the source.  Therefore phi[n+1], pi[n+1], and u[n+1] are computed in
        # that order without a coupled solve at the new radial degree.
        for n in range(maxdegree):
            rows = n + 1
            average_u = u[:rows] / np.arange(1, rows + 1, dtype=float)[:, None]
            average_u_eta = grid.differentiate(average_u)

            W = np.zeros_like(average_u)
            W[0] = 1.0 - 2.0 * (0.5 - h) * eta * average_u[0]
            W -= (1.0 - eta**2) * average_u_eta
            if rows > 1:
                W[1:] -= 2.0 * (0.5 - h) * eta * average_u[1:]

            Hc = (1.0 - eta**2)[None, :] * u[:rows]
            Hc[0] += (0.5 - h) * eta

            one_minus_2eta_u = -2.0 * eta[None, :] * u[:rows]
            one_minus_2eta_u[0] += 1.0

            phi_eta = grid.differentiate(phi[:rows])
            phi_plus = np.arange(1, rows + 1, dtype=float)[:, None] * phi[:rows]
            phi_source_eta = phi_eta + a[None, :] * phi[:rows]

            source_phi = (
                _coefficient_product(W, phi_plus, n)
                + h * _coefficient_product(one_minus_2eta_u, phi, n)
                + _coefficient_product(Hc, phi_source_eta, n)
            )
            phi[n + 1] = source_phi / (2.0 * (1.0 - 2.0 * h * eta**2) * (n + 1) * (n + 2))

            phi_square = _coefficient_product(phi, phi, n)
            pi[n + 1] = g**2 * phi_square / float(n + 1)

            u_eta = grid.differentiate(u[:rows])
            xu = np.arange(rows, dtype=float)[:, None] * u[:rows]
            xpi = np.arange(rows, dtype=float)[:, None] * pi[:rows]
            source_u = (
                _coefficient_product(W, xu, n)
                + (0.5 + h) * _coefficient_product(one_minus_2eta_u, u, n)
                + _coefficient_product(Hc, u_eta, n)
                + (1.0 - eta**2) * grid.differentiate(pi[:rows])[n]
                - 4.0 * (0.5 + h) * eta * pi[n]
                - 2.0 * eta * xpi[n]
            )
            u[n + 1] = source_u / (
                2.0 * (1.0 - 2.0 * h * eta**2) * (n + 1) ** 2
            )

        # Complete the pressure polynomial from the retained Phi prefix.  A
        # finite velocity prefix still determines every convolution row up to
        # degree 2N, so Pi_X=F**2 is represented without an additional radial
        # truncation in the evaluator.
        for n in range(2 * maxdegree + 1):
            pi[n + 1] = g**2 * _coefficient_product(phi, phi, n) / float(n + 1)

        for name, values in (
            ("g_values", g),
            ("a_values", a),
            ("phi_coefficients", phi),
            ("u_coefficients", u),
            ("pi_coefficients", pi),
        ):
            values = np.asarray(values, dtype=float)
            if not np.all(np.isfinite(values)):
                raise ArithmeticError(f"{name} contains a non-finite coefficient")
            values.setflags(write=False)
            object.__setattr__(self, name, values)
        object.__setattr__(self, "maxdegree", maxdegree)
        object.__setattr__(self, "eta_nodes", eta_nodes)
        object.__setattr__(self, "grid", grid)

    # Short aliases make the X-series convention visible to callers that use
    # the mathematical names from Eq. (4.13).
    @property
    def phi(self) -> np.ndarray:
        return self.phi_coefficients

    @property
    def Phi(self) -> np.ndarray:
        return self.phi_coefficients

    @property
    def u(self) -> np.ndarray:
        return self.u_coefficients

    @property
    def pi(self) -> np.ndarray:
        return self.pi_coefficients

    @property
    def pressure_coefficients(self) -> np.ndarray:
        return self.pi_coefficients

    @property
    def eta_derivative_matrix(self) -> np.ndarray:
        return self.grid.differentiation_matrix

    @property
    def average_u_coefficients(self) -> np.ndarray:
        return self.u_coefficients / np.arange(1, self.maxdegree + 2, dtype=float)[:, None]

    def _broadcast_domain(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_array = _finite_array(X, "X")
        eta_array = _finite_array(eta, "eta")
        if np.any(X_array < 0.0):
            raise ValueError("X must be nonnegative")
        if np.any(np.abs(eta_array) > 1.0 + 2e-13):
            raise ValueError("eta must satisfy |eta| <= 1")
        X_array, eta_array = np.broadcast_arrays(X_array, eta_array)
        if np.any(self.reference.Lambda * X_array > _DOMAIN_X_FACTOR + 2e-12):
            raise ValueError(
                f"finite core series only defined for Lambda*X <= {_DOMAIN_X_FACTOR}"
            )
        return X_array.astype(float, copy=False), np.clip(eta_array, -1.0, 1.0)

    def _evaluate_rows(self, coefficients: np.ndarray, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        grid_values = self.grid.interpolate(coefficients, eta_array)
        return _polynomial_value(grid_values, X_array)

    def _evaluate_derivative_rows(
        self, coefficients: np.ndarray, X: Any, eta: Any
    ) -> np.ndarray:
        derivative_coefficients = _polynomial_derivative_coefficients(coefficients)
        if derivative_coefficients.shape[0] == 1 and coefficients.shape[0] <= 1:
            X_array, eta_array = self._broadcast_domain(X, eta)
            return np.zeros_like(X_array, dtype=float)
        return self._evaluate_rows(derivative_coefficients, X, eta)

    def amplitude(self, eta: Any) -> np.ndarray:
        """Evaluate the cached ``PaperCoreReference`` axis amplitude ``g``.

        The coefficient recurrence uses ``g`` sampled on the eta grid.  Field
        evaluation calls the same reference amplitude at the requested eta so
        that ``F=g*Phi`` and the pressure identity ``Pi_X=F**2`` remain
        consistent between grid nodes as well.
        """

        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0 + 2e-13):
            raise ValueError("eta must satisfy |eta| <= 1")
        clipped = np.clip(eta_array, -1.0, 1.0)
        flat = np.array(
            [self.reference.amplitude(float(value)) for value in clipped.reshape(-1)],
            dtype=float,
        )
        return flat.reshape(clipped.shape)

    def log_amplitude_derivative(self, eta: Any) -> np.ndarray:
        """Return ``a(eta)=g'(eta)/g(eta)=Lambda*real_gradient(eta)``."""

        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0 + 2e-13):
            raise ValueError("eta must satisfy |eta| <= 1")
        clipped = np.clip(eta_array, -1.0, 1.0)
        ref = self.reference
        flat = np.array(
            [
                ref.Lambda
                * real_gradient(ref.h, ref.j, ref.sigma, float(value))
                for value in clipped.reshape(-1)
            ],
            dtype=float,
        )
        return flat.reshape(clipped.shape)

    def phi_value(self, X: Any, eta: Any) -> np.ndarray:
        return self._evaluate_rows(self.phi_coefficients, X, eta)

    def phi_radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        return self._evaluate_derivative_rows(self.phi_coefficients, X, eta)

    def phi_eta(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        coefficients = self.grid.differentiate(self.phi_coefficients)
        return _polynomial_value(self.grid.interpolate(coefficients, eta_array), X_array)

    def F(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        return self.amplitude(eta_array) * self.phi_value(X_array, eta_array)

    def F_radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        return self.amplitude(eta_array) * self.phi_radial_derivative(X_array, eta_array)

    def F_eta(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        phi = self.phi_value(X_array, eta_array)
        return self.amplitude(eta_array) * (
            self.phi_eta(X_array, eta_array)
            + self.log_amplitude_derivative(eta_array) * phi
        )

    def E(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        return np.sqrt(2.0 * X_array) * self.F(X_array, eta_array)

    def U(self, X: Any, eta: Any) -> np.ndarray:
        return self._evaluate_rows(self.u_coefficients, X, eta)

    def U_radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        return self._evaluate_derivative_rows(self.u_coefficients, X, eta)

    def dU_deta(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        coefficients = self.grid.differentiate(self.u_coefficients)
        return _polynomial_value(self.grid.interpolate(coefficients, eta_array), X_array)

    def Pi(self, X: Any, eta: Any) -> np.ndarray:
        # Integrate the full retained Phi**2 polynomial.  This keeps the
        # finite field's pressure derivative equal to F**2 through degree
        # 2*maxdegree, rather than introducing a second pressure truncation at
        # maxdegree.  The stored pi_coefficients remain available as the
        # grid-collocation rows used by the triangular recurrence.
        X_array, eta_array = self._broadcast_domain(X, eta)
        convolution = self._phi_square_at(eta_array)
        integral_coefficients = convolution / np.arange(
            1, convolution.shape[0] + 1, dtype=float
        ).reshape((convolution.shape[0],) + (1,) * X_array.ndim)
        integral = X_array * _polynomial_value(integral_coefficients, X_array)
        ref = self.reference
        pressure_axis = -(ref.pressure_scale**2) / (1.0 + eta_array**2) ** 2
        return pressure_axis + self.amplitude(eta_array) ** 2 * integral

    def Pi_radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        phi = _polynomial_value(
            self.grid.interpolate(self.phi_coefficients, eta_array), X_array
        )
        return self.amplitude(eta_array) ** 2 * phi**2

    def Pi_eta(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        phi_rows = self.grid.interpolate(self.phi_coefficients, eta_array)
        phi_eta_rows = self.grid.interpolate(
            self.grid.differentiate(self.phi_coefficients), eta_array
        )
        convolution = self._phi_square_at(eta_array)
        convolution_eta = np.zeros_like(convolution)
        degree = self.maxdegree
        for i in range(degree + 1):
            for j in range(degree + 1):
                convolution_eta[i + j] += (
                    phi_eta_rows[i] * phi_rows[j]
                    + phi_rows[i] * phi_eta_rows[j]
                )
        divisors = np.arange(
            1, convolution.shape[0] + 1, dtype=float
        ).reshape((convolution.shape[0],) + (1,) * X_array.ndim)
        integral = X_array * _polynomial_value(convolution / divisors, X_array)
        integral_eta = X_array * _polynomial_value(convolution_eta / divisors, X_array)
        ref = self.reference
        pressure_axis_eta = (
            4.0 * ref.pressure_scale**2 * eta_array / (1.0 + eta_array**2) ** 3
        )
        g = self.amplitude(eta_array)
        a = self.log_amplitude_derivative(eta_array)
        return pressure_axis_eta + g**2 * (2.0 * a * integral + integral_eta)

    def _phi_square_at(self, eta: np.ndarray) -> np.ndarray:
        """Convolution rows of the interpolated finite Phi prefix."""

        phi_rows = self.grid.interpolate(self.phi_coefficients, eta)
        convolution = np.zeros(
            (2 * self.maxdegree + 1,) + eta.shape, dtype=float
        )
        for i in range(self.maxdegree + 1):
            for j in range(self.maxdegree + 1):
                convolution[i + j] += phi_rows[i] * phi_rows[j]
        return convolution

    def radial_average_U(self, X: Any, eta: Any) -> np.ndarray:
        return self._evaluate_rows(self.average_u_coefficients, X, eta)

    def radial_average_dU_deta(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = self._broadcast_domain(X, eta)
        coefficients = self.grid.differentiate(self.average_u_coefficients)
        return _polynomial_value(self.grid.interpolate(coefficients, eta_array), X_array)

    def radial_pressure_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Return ``Pi_X`` for diagnostics; ``LeadingProfile`` uses ``F**2``."""

        return self.Pi_radial_derivative(X, eta)

    def pressure_radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Compatibility alias for the finite pressure's X derivative."""

        return self.Pi_radial_derivative(X, eta)

    @property
    def profile(self) -> LeadingProfile:
        """Adapt this finite prefix to the existing velocity evaluator."""

        return LeadingProfile(
            E=self.E,
            U=self.U,
            dU_deta=self.dU_deta,
            Pi=self.Pi,
            F=self.F,
            average_U=self.radial_average_U,
            average_dU_deta=self.radial_average_dU_deta,
            name="paper-core-finite-nonlinear-radial-series",
            paper_exact=False,
            provenance=(
                "Finite Chebyshev eta interpolation of the triangular nonlinear "
                "Eqs. (4.9)/(4.13) radial recurrence; independent finite pressure "
                "datum and PaperCoreReference axis amplitude; not converged."
            ),
        )

    def velocity(
        self,
        x: float,
        y: float,
        z: float,
        t: float,
        *,
        quadrature_points: int = 32,
    ) -> np.ndarray:
        """Sample the finite profile through the existing Cartesian adapter."""

        return leading_velocity_cartesian(
            x,
            y,
            z,
            t,
            self.profile,
            h=self.reference.h,
            quadrature_points=quadrature_points,
        )

    def pressure(self, x: float, y: float, z: float, t: float) -> float:
        """Sample pressure through the existing similarity-coordinate adapter."""

        return leading_pressure_cartesian(
            x, y, z, t, self.profile, h=self.reference.h
        )

    def metadata(self) -> dict[str, Any]:
        """Describe the finite numerical status and independent choices."""

        return {
            "parameters": asdict(self.reference),
            "maxdegree": self.maxdegree,
            "eta_nodes": self.eta_nodes,
            "series_variable": "X",
            "coefficient_convention": "ordinary powers X**n",
            "viscosity": 1.0,
            "paper_exact": False,
            "converged": False,
            "status": "finite formal near-axis nonlinear radial prefix",
            "domain": f"0 <= Lambda*X <= {_DOMAIN_X_FACTOR}, |eta| <= 1",
            "pressure_datum": (
                "-pressure_scale**2/(1+eta**2)**2, independent autonomous datum"
            ),
            "limitations": [
                "finite degree and Chebyshev interpolation truncation",
                "no convergence or residual theorem",
                "no exterior matching or oscillatory correction",
            ],
        }


# Concise factory aliases for callers that prefer a function over the class.
def build_paper_core_series(
    reference: PaperCoreReference | None = None,
    *,
    maxdegree: int = 8,
    eta_nodes: int = 129,
) -> PaperCoreSeries:
    """Build one bounded finite nonlinear paper-core radial prefix."""

    return PaperCoreSeries(
        reference=PaperCoreReference() if reference is None else reference,
        maxdegree=maxdegree,
        eta_nodes=eta_nodes,
    )


def paper_core_profile(
    reference: PaperCoreReference | None = None,
    *,
    maxdegree: int = 8,
    eta_nodes: int = 129,
) -> LeadingProfile:
    """Build a finite-series :class:`LeadingProfile` adapter."""

    return build_paper_core_series(
        reference, maxdegree=maxdegree, eta_nodes=eta_nodes
    ).profile


__all__ = [
    "ChebyshevEtaGrid",
    "PaperCoreSeries",
    "build_paper_core_series",
    "paper_core_profile",
]
