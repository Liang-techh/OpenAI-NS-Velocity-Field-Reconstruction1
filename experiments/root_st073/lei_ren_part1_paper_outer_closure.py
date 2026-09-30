"""Algebraic closure nuclei from Lei--Ren Part I, Sections 7.4--7.5.

This file contains only the finite-dimensional coefficient solves.  It does
not construct the source paper's outer profile, its heat tail, or the
waiting-length root.  A caller must provide the actual profile-derived
weights and moments.  In particular, a successful check in this module is
not an outer-interval or Navier--Stokes closure certificate.

The angular solver is the exact two-equation system (7.21), reduced as in
(7.23) and solved by the exact quadratic root.  The axial solver follows
(7.31) and (7.34): first ``c_j(a)=u_j+v_j*a`` is obtained with the source
sign, and then the resulting energy quadratic is solved on ``[.9,1.2]``.
No asymptotic coefficient is substituted for either exact root.

For the construction's very small ``mu`` values, use
``stable_paper_axial_weights`` and
``assemble_stable_paper_axial_affine``.  They use the finite row difference
``(W(lambda_1)-W(lambda_2))/mu``.  Below ``1e-12`` the corresponding source
RHS differences require Decimal or decimal-string input; binary floats are
rejected rather than silently treated as exact.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
import math
from typing import Any, Iterable, Sequence

import numpy as np


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTIONS = "Part I Sections 4.1, 4.12--4.13, 7.4--7.5, equations (4.1), (7.21), (7.23), (7.31)--(7.35)"
ANGULAR_EQUATION = "(7.21), normalized exactly as (7.23)"
AXIAL_LINEAR_EQUATION = "(7.31)"
AXIAL_ENERGY_EQUATION = "(7.34)"
PAPER_BUMP_DEFINITION = "(4.1), (4.13): beta_ell(t)=b(t/ell)/(ell*int_{-1}^1 b), ell=3/20"
PAPER_GAMMA_DEFINITION = "(4.12), (7.28)--(7.29): gamma_1(t)=beta(t-13/mu+3), gamma_2(t)=beta(t-13/mu+1)"

_EPS = 1.0e-14
_FLOAT_RHS_CUTOFF = 1.0e-12


@dataclass(frozen=True)
class PaperBumpWeights:
    """Quadrature-derived weights for the fixed paper bump ``beta_ell``.

    The source fixes the outer-correction width to ``ell=3/20`` in (4.13),
    while (4.1) defines the normalized family for any positive ``ell``.
    This object records the selected width and quadrature order so a caller
    cannot mistake a synthetic bump for the source bump.
    """

    mu: float
    ell: float
    quadrature_order: int
    normalization_integral: float
    beta_inf: float
    beta_l2_squared: float
    A_mu: float
    B_mu: float
    D_mu: float


@dataclass(frozen=True)
class PaperAxialWeights:
    """Normalized end-bump matrix and Gram weights from (7.31)--(7.34)."""

    mu: float
    ell: float
    quadrature_order: int
    lambda_values: tuple[float, float]
    matrix: tuple[tuple[float, float], tuple[float, float]]
    K_bump: tuple[float, float]
    gamma_centers: tuple[float, float]
    gamma_supports: tuple[tuple[float, float], tuple[float, float]]
    beta_l2_squared: float


@dataclass(frozen=True)
class StablePaperAxialWeights:
    """Finite row-difference representation of the Section 7.5 matrix.

    The original normalized rows are ``W(lambda_1)`` and ``W(lambda_2)``.
    This representation stores row zero as ``W(lambda_1)`` and row one as
    ``(W(lambda_1)-W(lambda_2))/mu``.  Thus its determinant stays finite as
    ``mu`` tends to zero, while the original determinant is exactly
    ``-mu * normalized_determinant``.  The absolute ``t=13/mu`` geometry is
    deliberately optional: stage-local offsets remain available even when
    binary floating point makes the translated support endpoints collapse.
    """

    mu: float
    ell: float
    quadrature_order: int
    lambda_values: tuple[float, float]
    matrix: tuple[tuple[float, float], tuple[float, float]]
    normalized_determinant: float
    source_determinant_estimate: float
    source_determinant_over_mu: float
    condition: float
    row_scaling: tuple[float, float]
    gamma_offset_centers: tuple[float, float]
    gamma_offset_supports: tuple[tuple[float, float], tuple[float, float]]
    absolute_geometry_supported: bool
    gamma_centers: tuple[float, float] | None
    gamma_supports: tuple[tuple[float, float], tuple[float, float]] | None
    beta_l2_squared: float


@dataclass(frozen=True)
class StablePaperAxialAssembly:
    """Stable gamma weights together with the exactly row-scaled affine solve."""

    weights: StablePaperAxialWeights
    affine: "AxialAffineCorrection"
    source_base: tuple[str, str]
    source_pulse: tuple[str, str]
    scaled_base: tuple[float, float]
    scaled_pulse: tuple[float, float]
    rhs_exact: bool
    source_rows_replay_available: bool

    def __iter__(self):
        # Keep tuple-style unpacking convenient while retaining metadata.
        yield self.weights
        yield self.affine


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _vector(value: Sequence[Any] | np.ndarray, name: str, length: int) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (length,):
        raise ValueError(f"{name} must have shape ({length},), got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array.copy()


def _matrix(value: Sequence[Sequence[Any]] | np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (2, 2):
        raise ValueError(f"{name} must have shape (2,2), got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array.copy()


def _paper_raw_bump(unit_coordinate: np.ndarray | float) -> np.ndarray | float:
    """The compactly supported raw bump ``mathfrak b`` from (4.1)."""

    values = np.asarray(unit_coordinate, dtype=float)
    result = np.zeros_like(values, dtype=float)
    inside = np.abs(values) < 1.0
    result[inside] = np.exp(-1.0 / (1.0 - values[inside] * values[inside]))
    if np.ndim(unit_coordinate) == 0:
        return float(result)
    return result


def paper_bump_weights(*, mu: Any, ell: Any = 3.0 / 20.0, quadrature_order: int = 128) -> PaperBumpWeights:
    """Construct the source ``beta_ell`` and its exact Section 7 weights.

    The integral is evaluated on the compact support with Gauss--Legendre
    quadrature.  ``ell=3/20`` is the paper's outer choice; other positive
    widths are allowed by (4.1) and are recorded explicitly.  The source
    range ``0 <= mu <= 1/60`` is enforced.
    """

    mu_value = _finite(mu, "mu")
    ell_value = _finite(ell, "ell")
    if not 0.0 <= mu_value <= 1.0 / 60.0:
        raise ValueError("mu must lie in the Part I range [0,1/60]")
    if ell_value <= 0.0:
        raise ValueError("ell must be positive")
    try:
        order = int(quadrature_order)
    except (TypeError, ValueError) as exc:
        raise ValueError("quadrature_order must be an integer at least 16") from exc
    if order < 16 or order != quadrature_order:
        raise ValueError("quadrature_order must be an integer at least 16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    raw = np.asarray(_paper_raw_bump(nodes), dtype=float)
    normalization = float(np.sum(weights * raw))
    if normalization <= 0.0 or not math.isfinite(normalization):
        raise ValueError("paper bump normalization quadrature failed")
    support_nodes = ell_value * nodes
    support_weights = ell_value * weights
    beta = raw / (ell_value * normalization)
    beta_l2_squared = float(np.sum(support_weights * beta * beta))
    # The raw bump is even and strictly increasing to its maximum at zero on
    # [-1,0], so this is its exact supremum rather than a grid estimate.
    beta_inf = math.exp(-1.0) / (ell_value * normalization)
    A_mu = float(np.sum(support_weights * np.exp((1.0 - mu_value) * support_nodes) * beta))
    B_mu = float(np.sum(support_weights * np.exp(-(1.0 + 2.0 * mu_value) * support_nodes) * beta))
    D_mu = float(
        np.sum(support_weights * np.exp(-(1.0 + 2.0 * mu_value) * support_nodes) * beta * beta)
    )
    return PaperBumpWeights(
        mu=mu_value,
        ell=ell_value,
        quadrature_order=order,
        normalization_integral=normalization,
        beta_inf=float(beta_inf),
        beta_l2_squared=beta_l2_squared,
        A_mu=A_mu,
        B_mu=B_mu,
        D_mu=D_mu,
    )


def solve_angular_bumps_from_paper_bump(
    *,
    mu: Any,
    r: Any,
    s_H: Any,
    ell: Any = 3.0 / 20.0,
    quadrature_order: int = 128,
    require_positive: bool = False,
    conditioning_limit: float = 1.0e12,
) -> AngularBumpSolution:
    """Build source bump weights, then solve the exact (7.21) system."""

    weights = paper_bump_weights(mu=mu, ell=ell, quadrature_order=quadrature_order)
    return solve_angular_bumps(
        mu=weights.mu,
        A_mu=weights.A_mu,
        B_mu=weights.B_mu,
        D_mu=weights.D_mu,
        r=r,
        s_H=s_H,
        beta_inf=weights.beta_inf,
        require_positive=require_positive,
        conditioning_limit=conditioning_limit,
    )


def paper_axial_weights(
    *, mu: Any, ell: Any = 3.0 / 20.0, quadrature_order: int = 128
) -> PaperAxialWeights:
    """Assemble the normalized gamma matrix and bump Gram weights.

    For ``gamma_1(t)=beta(t-13/mu+3)`` and
    ``gamma_2(t)=beta(t-13/mu+1)``, the source-normalized rows of (7.31) are
    ``exp(-3 lambda_i) int exp(lambda_i s) beta(s) ds`` and
    ``exp(-lambda_i) int exp(lambda_i s) beta(s) ds``.  The returned ``K_bump``
    is the unnormalized Gram factor in (7.34), evaluated in its stable
    translated form.  Pulse integrals and actual incoming targets remain
    caller data.
    """

    bump = paper_bump_weights(mu=mu, ell=ell, quadrature_order=quadrature_order)
    if bump.mu < _FLOAT_RHS_CUTOFF:
        raise ValueError(
            "ordinary absolute gamma coordinates lose float resolution for mu<1e-12; "
            "use stable_paper_axial_weights and exact Decimal/string RHS data"
        )
    nodes, weights = np.polynomial.legendre.leggauss(bump.quadrature_order)
    support_nodes = bump.ell * nodes
    support_weights = bump.ell * weights
    raw = np.asarray(_paper_raw_bump(nodes), dtype=float)
    beta = raw / (bump.ell * bump.normalization_integral)
    lambdas = (0.5 - bump.mu, 0.5 - 2.0 * bump.mu)
    matrix_rows: list[tuple[float, float]] = []
    for lam in lambdas:
        translated_integral = float(np.sum(support_weights * np.exp(lam * support_nodes) * beta))
        matrix_rows.append(
            (
                math.exp(-3.0 * lam) * translated_integral,
                math.exp(-lam) * translated_integral,
            )
        )
    gram_integral = float(
        np.sum(support_weights * np.exp(-2.0 * bump.mu * support_nodes) * beta * beta)
    )
    K_bump = (
        math.exp(-26.0 + 6.0 * bump.mu) * gram_integral,
        math.exp(-26.0 + 2.0 * bump.mu) * gram_integral,
    )
    gamma_centers = (13.0 / bump.mu - 3.0, 13.0 / bump.mu - 1.0)
    gamma_supports = tuple(
        (center - bump.ell, center + bump.ell) for center in gamma_centers
    )
    return PaperAxialWeights(
        mu=bump.mu,
        ell=bump.ell,
        quadrature_order=bump.quadrature_order,
        lambda_values=(float(lambdas[0]), float(lambdas[1])),
        matrix=(matrix_rows[0], matrix_rows[1]),
        K_bump=(float(K_bump[0]), float(K_bump[1])),
        gamma_centers=gamma_centers,
        gamma_supports=gamma_supports,
        beta_l2_squared=bump.beta_l2_squared,
    )


def assemble_paper_axial_affine(
    *,
    mu: Any,
    base: Sequence[Any] | np.ndarray,
    pulse: Sequence[Any] | np.ndarray,
    ell: Any = 3.0 / 20.0,
    quadrature_order: int = 128,
    conditioning_limit: float = 1.0e12,
) -> tuple[PaperAxialWeights, AxialAffineCorrection]:
    """Build the source gamma matrix, then apply the exact (7.31) signs."""

    weights = paper_axial_weights(mu=mu, ell=ell, quadrature_order=quadrature_order)
    affine = solve_axial_affine(
        matrix=weights.matrix,
        base=base,
        pulse=pulse,
        conditioning_limit=conditioning_limit,
    )
    return weights, affine


def _paper_bump_arrays(bump: PaperBumpWeights) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(bump.quadrature_order)
    support_nodes = bump.ell * nodes
    support_weights = bump.ell * weights
    raw = np.asarray(_paper_raw_bump(nodes), dtype=float)
    beta = raw / (bump.ell * bump.normalization_integral)
    return support_nodes, support_weights, beta


def stable_paper_axial_weights(
    *,
    mu: Any,
    ell: Any = 3.0 / 20.0,
    quadrature_order: int = 128,
) -> StablePaperAxialWeights:
    """Assemble the finite row-difference matrix for tiny ``mu``.

    With ``W_k(lambda)=int exp(lambda*(s-shift_k))*beta(s) ds``, the matrix
    rows are ``W(lambda_1)`` and ``(W(lambda_1)-W(lambda_2))/mu``.  The
    latter is evaluated as ``expm1`` under the integral, so it remains
    meaningful when ``lambda_1`` and ``lambda_2`` round to the same float.
    This only stabilizes the matrix.  The corresponding RHS difference is
    handled by :func:`assemble_stable_paper_axial_affine`, which requires
    exact Decimal/string input below ``1e-12``.
    """

    bump = paper_bump_weights(mu=mu, ell=ell, quadrature_order=quadrature_order)
    if bump.mu <= 0.0:
        raise ValueError("stable axial row-difference assembly requires mu>0")
    support_nodes, support_weights, beta = _paper_bump_arrays(bump)
    lambda_1 = 0.5 - bump.mu
    lambda_2 = 0.5 - 2.0 * bump.mu
    shifts = (3.0, 1.0)
    row_zero: list[float] = []
    row_difference: list[float] = []
    for shift in shifts:
        shifted = support_nodes - shift
        kernel = np.exp(lambda_1 * shifted) * beta
        row_zero.append(float(np.sum(support_weights * kernel)))
        # -expm1(-mu*shifted)/mu is the divided difference of the two
        # exponentials and avoids subtracting equal W rows.
        divided_difference = -np.expm1(-bump.mu * shifted) / bump.mu
        row_difference.append(float(np.sum(support_weights * kernel * divided_difference)))
    matrix = np.asarray([row_zero, row_difference], dtype=float)
    normalized_determinant = float(np.linalg.det(matrix))
    if not math.isfinite(normalized_determinant) or abs(normalized_determinant) <= 1.0e-14:
        raise ValueError(
            "stable axial row-difference matrix is degenerate: "
            f"normalized_determinant={normalized_determinant:.17g}"
        )
    condition = float(np.linalg.cond(matrix))
    if not math.isfinite(condition):
        raise ValueError("stable axial row-difference matrix has nonfinite condition")
    gamma_offset_centers = (-3.0, -1.0)
    gamma_offset_supports = (
        (-3.0 - bump.ell, -3.0 + bump.ell),
        (-1.0 - bump.ell, -1.0 + bump.ell),
    )
    absolute_geometry_supported = bump.mu >= _FLOAT_RHS_CUTOFF
    if absolute_geometry_supported:
        gamma_centers = tuple(13.0 / bump.mu + offset for offset in gamma_offset_centers)
        gamma_supports = tuple(
            (center - bump.ell, center + bump.ell) for center in gamma_centers
        )
        absolute_geometry_supported = all(
            low < high and low != high
            for low, high in gamma_supports
        )
        if not absolute_geometry_supported:
            gamma_centers = None
            gamma_supports = None
    else:
        gamma_centers = None
        gamma_supports = None
    return StablePaperAxialWeights(
        mu=bump.mu,
        ell=bump.ell,
        quadrature_order=bump.quadrature_order,
        lambda_values=(float(lambda_1), float(lambda_2)),
        matrix=(tuple(float(value) for value in matrix[0]), tuple(float(value) for value in matrix[1])),
        normalized_determinant=normalized_determinant,
        source_determinant_estimate=-bump.mu * normalized_determinant,
        source_determinant_over_mu=-normalized_determinant,
        condition=condition,
        row_scaling=(1.0, 1.0 / bump.mu),
        gamma_offset_centers=gamma_offset_centers,
        gamma_offset_supports=gamma_offset_supports,
        absolute_geometry_supported=absolute_geometry_supported,
        gamma_centers=gamma_centers,
        gamma_supports=gamma_supports,
        beta_l2_squared=bump.beta_l2_squared,
    )


def _decimal_exact(value: Any, name: str) -> Decimal:
    if isinstance(value, Decimal):
        result = value
    elif isinstance(value, str):
        try:
            result = Decimal(value)
        except InvalidOperation as exc:
            raise ValueError(f"{name} must be a finite Decimal string") from exc
    else:
        raise ValueError(
            f"{name} must be supplied as Decimal or decimal string when mu<1e-12; "
            "a binary float cannot recover the divided RHS"
        )
    if not result.is_finite():
        raise ValueError(f"{name} must be finite")
    return result


def _stable_rhs_pair(
    values: Sequence[Any] | np.ndarray,
    *,
    mu: float,
    name: str,
) -> tuple[np.ndarray, tuple[str, str], bool]:
    try:
        if len(values) != 2:
            raise ValueError(f"{name} must have two entries")
    except TypeError as exc:
        raise ValueError(f"{name} must have two entries") from exc
    exact_required = mu < _FLOAT_RHS_CUTOFF
    if exact_required:
        first = _decimal_exact(values[0], f"{name}[0]")
        second = _decimal_exact(values[1], f"{name}[1]")
        with localcontext() as context:
            context.prec = 100
            mu_decimal = Decimal(str(mu))
            if mu_decimal <= 0:
                raise ValueError("mu must be positive for divided RHS")
            difference = (first - second) / mu_decimal
        scaled = np.asarray([float(first), float(difference)], dtype=float)
        return scaled, (str(first), str(second)), True
    array = _vector(values, name, 2)
    if mu <= 0.0:
        raise ValueError("mu must be positive for divided RHS")
    scaled = np.asarray([array[0], (array[0] - array[1]) / mu], dtype=float)
    return scaled, (str(values[0]), str(values[1])), False


def assemble_stable_paper_axial_affine(
    *,
    mu: Any,
    base: Sequence[Any] | np.ndarray,
    pulse: Sequence[Any] | np.ndarray,
    ell: Any = 3.0 / 20.0,
    quadrature_order: int = 128,
    conditioning_limit: float = 1.0e12,
) -> StablePaperAxialAssembly:
    """Solve (7.31) after exact row-difference scaling for tiny ``mu``.

    The transformed RHS is ``[base_1,(base_1-base_2)/mu]`` and likewise for
    ``pulse``.  Below ``1e-12`` both original pairs must be Decimal objects or
    decimal strings; float inputs are rejected rather than numerically
    reconstructed.  The returned affine coefficients solve the scaled system
    with the same source sign ``A c=-base-a*pulse``.
    """

    weights = stable_paper_axial_weights(mu=mu, ell=ell, quadrature_order=quadrature_order)
    scaled_base, source_base, base_exact = _stable_rhs_pair(
        base, mu=weights.mu, name="base"
    )
    scaled_pulse, source_pulse, pulse_exact = _stable_rhs_pair(
        pulse, mu=weights.mu, name="pulse"
    )
    affine = solve_axial_affine(
        matrix=weights.matrix,
        base=scaled_base,
        pulse=scaled_pulse,
        conditioning_limit=conditioning_limit,
    )
    return StablePaperAxialAssembly(
        weights=weights,
        affine=affine,
        source_base=source_base,
        source_pulse=source_pulse,
        scaled_base=(float(scaled_base[0]), float(scaled_base[1])),
        scaled_pulse=(float(scaled_pulse[0]), float(scaled_pulse[1])),
        rhs_exact=bool(base_exact and pulse_exact),
        source_rows_replay_available=weights.absolute_geometry_supported,
    )


def _real_quadratic_roots(a: float, b: float, c: float, *, name: str) -> tuple[float, ...]:
    """Return all finite real roots without replacing a quadratic by a series."""

    a = _finite(a, f"{name}.a")
    b = _finite(b, f"{name}.b")
    c = _finite(c, f"{name}.c")
    scale = max(abs(b * b), abs(4.0 * a * c), 1.0)
    discriminant = b * b - 4.0 * a * c
    if discriminant < -1.0e-13 * scale:
        raise ValueError(f"{name} has no real root: discriminant={discriminant:.17g}")
    if abs(discriminant) <= 1.0e-13 * scale:
        discriminant = 0.0
    if abs(a) <= _EPS * max(abs(b), abs(c), 1.0):
        if abs(b) <= _EPS * max(abs(c), 1.0):
            raise ValueError(f"{name} is degenerate: both quadratic coefficients vanish")
        return (-c / b,)
    root_discriminant = math.sqrt(max(discriminant, 0.0))
    # The two direct roots are adequate in the source parameter range and
    # retain both branches for the explicit small-branch selection below.
    roots = [(-b - root_discriminant) / (2.0 * a), (-b + root_discriminant) / (2.0 * a)]
    unique: list[float] = []
    for root in roots:
        if math.isfinite(root) and not any(abs(root - old) <= 1.0e-13 * max(1.0, abs(root), abs(old)) for old in unique):
            unique.append(root)
    if not unique:
        raise ValueError(f"{name} has no finite real root")
    return tuple(unique)


@dataclass(frozen=True)
class AngularBumpSolution:
    """Exact small-branch solution of the normalized angular system."""

    mu: float
    A_mu: float
    B_mu: float
    D_mu: float
    r: float
    s_H: float
    d1: float
    d2: float
    alpha: float
    beta2: float
    rho: float
    sigma: float
    q: float
    discriminant: float
    scalar_derivative: float
    jacobian_condition: float
    linear_determinant: float
    residual_equations: tuple[float, float]
    branch: str
    reachable: bool
    well_conditioned: bool
    positivity_checked: bool
    positive: bool | None
    beta_inf: float | None

    @property
    def residual_max(self) -> float:
        return max(abs(value) for value in self.residual_equations)

    def implicit_derivative(
        self,
        *,
        dr: float = 0.0,
        ds_H: float = 0.0,
        dmu: float = 0.0,
        dA_mu: float = 0.0,
        dB_mu: float = 0.0,
        dD_mu: float = 0.0,
    ) -> "AngularBumpDerivative":
        return angular_implicit_derivative(
            self,
            dr=dr,
            ds_H=ds_H,
            dmu=dmu,
            dA_mu=dA_mu,
            dB_mu=dB_mu,
            dD_mu=dD_mu,
        )


@dataclass(frozen=True)
class AngularBumpDerivative:
    """Implicit derivative ``(dd1, dd2)`` for a supplied input direction."""

    dd1: float
    dd2: float
    drho: float
    dsigma: float
    jacobian_condition: float


def _angular_normalized_parameters(
    mu: float, A_mu: float, B_mu: float, D_mu: float, r: float, s_H: float
) -> tuple[float, float, float, float, float, float]:
    mu = _finite(mu, "mu")
    A_mu = _finite(A_mu, "A_mu")
    B_mu = _finite(B_mu, "B_mu")
    D_mu = _finite(D_mu, "D_mu")
    r = _finite(r, "r")
    s_H = _finite(s_H, "s_H")
    if not 0.0 <= mu <= 1.0 / 60.0:
        raise ValueError("mu must lie in the Part I range [0,1/60]")
    if A_mu <= 0.0 or B_mu <= 0.0 or D_mu <= 0.0:
        raise ValueError("A_mu, B_mu and D_mu must be strictly positive")
    alpha = math.exp(-2.0 * (1.0 - mu))
    beta2 = math.exp(-2.0 * (1.0 + 2.0 * mu))
    rho = math.exp(1.0 - mu) * r / A_mu
    sigma = math.exp(-3.0 * (1.0 + 2.0 * mu)) * s_H / B_mu
    q = D_mu / (2.0 * B_mu)
    return mu, alpha, beta2, rho, sigma, q


def _angular_equations(
    *, alpha: float, beta2: float, rho: float, sigma: float, q: float, d1: float, d2: float
) -> tuple[float, float]:
    return (
        alpha * d1 + d2 - rho,
        d1 + beta2 * d2 + q * (d1 * d1 + beta2 * d2 * d2) - sigma,
    )


def solve_angular_bumps(
    *,
    mu: Any,
    A_mu: Any,
    B_mu: Any,
    D_mu: Any,
    r: Any,
    s_H: Any,
    beta_inf: Any | None = None,
    require_positive: bool = False,
    conditioning_limit: float = 1.0e12,
) -> AngularBumpSolution:
    """Solve the exact angular bump equations (7.21).

    ``r`` and ``s_H`` are the *actual* profile discrepancies at the supplied
    axial coordinate.  ``beta_inf`` is an optional independently established
    ``||beta||_infinity``.  If it is supplied, positivity checks
    ``1+d_j beta >= 0`` using the sufficient bound
    ``1+min(d_j) beta_inf >= 0``.  The solver does not silently assert
    positivity when the bump norm was not supplied.
    """

    mu_value, alpha, beta2, rho, sigma, q = _angular_normalized_parameters(
        mu, A_mu, B_mu, D_mu, r, s_H
    )
    A_value = _finite(A_mu, "A_mu")
    B_value = _finite(B_mu, "B_mu")
    D_value = _finite(D_mu, "D_mu")
    r_value = _finite(r, "r")
    s_value = _finite(s_H, "s_H")
    linear_determinant = alpha * beta2 - 1.0
    if abs(linear_determinant) <= 1.0e-12:
        raise ValueError("angular linear system is degenerate")

    # Substitute d2=rho-alpha*d1 into the second normalized equation.
    c2 = q * (1.0 + beta2 * alpha * alpha)
    c1 = 1.0 - beta2 * alpha - 2.0 * q * beta2 * alpha * rho
    c0 = beta2 * rho + q * beta2 * rho * rho - sigma
    scale = max(abs(c1 * c1), abs(4.0 * c2 * c0), 1.0)
    discriminant = c1 * c1 - 4.0 * c2 * c0
    if discriminant < -1.0e-13 * scale:
        raise ValueError(
            "angular data are unreachable on the real small-branch solve: "
            f"discriminant={discriminant:.17g}"
        )
    if abs(discriminant) <= 1.0e-13 * scale:
        discriminant = 0.0
    roots = _real_quadratic_roots(c2, c1, c0, name="angular scalar equation")
    # The small branch is the one continuous with the q=0 linear solve.
    if abs(c1) <= _EPS:
        linear_d1 = -c0 / max(abs(c1), _EPS)
    else:
        linear_d1 = -c0 / c1
    d1 = min(roots, key=lambda value: abs(value - linear_d1))
    d2 = rho - alpha * d1
    residuals = _angular_equations(
        alpha=alpha, beta2=beta2, rho=rho, sigma=sigma, q=q, d1=d1, d2=d2
    )
    scalar_derivative = c1 + 2.0 * c2 * d1
    jacobian = np.array(
        [[alpha, 1.0], [1.0 + 2.0 * q * d1, beta2 * (1.0 + 2.0 * q * d2)]],
        dtype=float,
    )
    jacobian_condition = float(np.linalg.cond(jacobian))
    well_conditioned = (
        math.isfinite(jacobian_condition)
        and abs(scalar_derivative) > 1.0e-12
        and jacobian_condition <= float(conditioning_limit)
    )
    if not well_conditioned:
        raise ValueError(
            "angular small-branch Jacobian is ill-conditioned: "
            f"cond={jacobian_condition:.17g}, scalar_derivative={scalar_derivative:.17g}"
        )

    if beta_inf is None:
        beta_bound = None
        positivity_checked = False
        positive: bool | None = None
    else:
        beta_bound = _finite(beta_inf, "beta_inf")
        if beta_bound < 0.0:
            raise ValueError("beta_inf must be nonnegative")
        positivity_checked = True
        positive = min(1.0 + d1 * beta_bound, 1.0 + d2 * beta_bound) >= -1.0e-12
        if require_positive and not positive:
            raise ValueError(
                "angular bump multiplier fails the supplied positivity bound: "
                f"d1={d1:.17g}, d2={d2:.17g}, beta_inf={beta_bound:.17g}"
            )
    return AngularBumpSolution(
        mu=mu_value,
        A_mu=A_value,
        B_mu=B_value,
        D_mu=D_value,
        r=r_value,
        s_H=s_value,
        d1=float(d1),
        d2=float(d2),
        alpha=float(alpha),
        beta2=float(beta2),
        rho=float(rho),
        sigma=float(sigma),
        q=float(q),
        discriminant=float(discriminant),
        scalar_derivative=float(scalar_derivative),
        jacobian_condition=jacobian_condition,
        linear_determinant=float(linear_determinant),
        residual_equations=(float(residuals[0]), float(residuals[1])),
        branch="small_q_continuation",
        reachable=True,
        well_conditioned=well_conditioned,
        positivity_checked=positivity_checked,
        positive=positive,
        beta_inf=beta_bound,
    )


def angular_implicit_derivative(
    solution: AngularBumpSolution,
    *,
    dr: Any = 0.0,
    ds_H: Any = 0.0,
    dmu: Any = 0.0,
    dA_mu: Any = 0.0,
    dB_mu: Any = 0.0,
    dD_mu: Any = 0.0,
) -> AngularBumpDerivative:
    """Differentiate the exact normalized angular equations implicitly.

    All directions are optional.  Setting only ``dr`` or ``ds_H`` gives the
    common profile-discrepancy derivatives; the remaining arguments permit a
    full derivative when the computed weights vary with the parameter.
    """

    dr = _finite(dr, "dr")
    ds_H = _finite(ds_H, "ds_H")
    dmu = _finite(dmu, "dmu")
    dA_mu = _finite(dA_mu, "dA_mu")
    dB_mu = _finite(dB_mu, "dB_mu")
    dD_mu = _finite(dD_mu, "dD_mu")
    d_alpha = 2.0 * solution.alpha * dmu
    d_beta2 = -4.0 * solution.beta2 * dmu
    d_rho = (
        math.exp(1.0 - solution.mu) / solution.A_mu * dr
        - solution.rho * dmu
        - solution.rho * dA_mu / solution.A_mu
    )
    d_sigma = (
        math.exp(-3.0 * (1.0 + 2.0 * solution.mu)) / solution.B_mu * ds_H
        - 6.0 * solution.sigma * dmu
        - solution.sigma * dB_mu / solution.B_mu
    )
    d_q = solution.q * (dD_mu / solution.D_mu - dB_mu / solution.B_mu)
    Q = solution.d1 * solution.d1 + solution.beta2 * solution.d2 * solution.d2
    jacobian = np.array(
        [[solution.alpha, 1.0], [1.0 + 2.0 * solution.q * solution.d1, solution.beta2 * (1.0 + 2.0 * solution.q * solution.d2)]],
        dtype=float,
    )
    rhs = np.array(
        [
            d_rho - d_alpha * solution.d1,
            d_sigma
            - d_beta2 * solution.d2
            - d_q * Q
            - solution.q * d_beta2 * solution.d2 * solution.d2,
        ],
        dtype=float,
    )
    try:
        derivative = np.linalg.solve(jacobian, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("angular implicit derivative has a singular Jacobian") from exc
    return AngularBumpDerivative(
        dd1=float(derivative[0]),
        dd2=float(derivative[1]),
        drho=float(d_rho),
        dsigma=float(d_sigma),
        jacobian_condition=float(np.linalg.cond(jacobian)),
    )


@dataclass(frozen=True)
class AxialAffineCorrection:
    """The exact affine coefficients ``c(a)=u+v*a`` from (7.31)."""

    matrix: tuple[tuple[float, float], tuple[float, float]]
    base: tuple[float, float]
    pulse: tuple[float, float]
    u: tuple[float, float]
    v: tuple[float, float]
    determinant: float
    condition: float
    source_sign: str

    def coefficients(self, a: Any) -> tuple[float, float]:
        amplitude = _finite(a, "a")
        u = np.asarray(self.u, dtype=float)
        v = np.asarray(self.v, dtype=float)
        result = u + amplitude * v
        return float(result[0]), float(result[1])

    def residual(self, a: Any) -> tuple[float, float]:
        amplitude = _finite(a, "a")
        matrix = np.asarray(self.matrix, dtype=float)
        base = np.asarray(self.base, dtype=float)
        pulse = np.asarray(self.pulse, dtype=float)
        coeff = np.asarray(self.coefficients(amplitude), dtype=float)
        residual = matrix @ coeff + base + amplitude * pulse
        return float(residual[0]), float(residual[1])


def solve_axial_affine(
    *,
    matrix: Sequence[Sequence[Any]] | np.ndarray,
    base: Sequence[Any] | np.ndarray,
    pulse: Sequence[Any] | np.ndarray,
    conditioning_limit: float = 1.0e12,
) -> AxialAffineCorrection:
    """Solve ``A c(a) = -base - a pulse`` exactly as in (7.31)."""

    A = _matrix(matrix, "matrix")
    m = _vector(base, "base", 2)
    p = _vector(pulse, "pulse", 2)
    determinant = float(np.linalg.det(A))
    condition = float(np.linalg.cond(A))
    if abs(determinant) <= 1.0e-14 * max(float(np.linalg.norm(A, ord=2) ** 2), 1.0):
        raise ValueError(f"axial affine matrix is degenerate: determinant={determinant:.17g}")
    if not math.isfinite(condition) or condition > float(conditioning_limit):
        raise ValueError(f"axial affine matrix is ill-conditioned: cond={condition:.17g}")
    u = np.linalg.solve(A, -m)
    v = np.linalg.solve(A, -p)
    return AxialAffineCorrection(
        matrix=tuple(tuple(float(value) for value in row) for row in A),
        base=(float(m[0]), float(m[1])),
        pulse=(float(p[0]), float(p[1])),
        u=(float(u[0]), float(u[1])),
        v=(float(v[0]), float(v[1])),
        determinant=determinant,
        condition=condition,
        source_sign="A c = -base - a*pulse",
    )


@dataclass(frozen=True)
class AxialEnergySolution:
    """Exact root of the source energy equation (7.34)."""

    affine: AxialAffineCorrection
    K_p: float
    K_bump: tuple[float, float]
    mu: float
    energy_target: float
    interval: tuple[float, float]
    quadratic_coefficients: tuple[float, float, float]
    roots: tuple[float, ...]
    a_p: float
    c_at_root: tuple[float, float]
    endpoint_values: tuple[float, float]
    derivative_at_root: float
    discriminant: float
    reachable: bool
    well_conditioned: bool

    def residual(self, a: Any | None = None) -> float:
        amplitude = self.a_p if a is None else _finite(a, "a")
        qa, qb, qc = self.quadratic_coefficients
        return float(qa * amplitude * amplitude + qb * amplitude + qc)

    def implicit_derivative(
        self,
        *,
        d_energy_target: float = 0.0,
        d_K_p: float = 0.0,
        d_K_bump: Iterable[float] = (0.0, 0.0),
        d_mu: float = 0.0,
        d_u: Iterable[float] = (0.0, 0.0),
        d_v: Iterable[float] = (0.0, 0.0),
    ) -> float:
        return axial_energy_implicit_derivative(
            self,
            d_energy_target=d_energy_target,
            d_K_p=d_K_p,
            d_K_bump=d_K_bump,
            d_mu=d_mu,
            d_u=d_u,
            d_v=d_v,
        )


def solve_axial_energy_root(
    affine: AxialAffineCorrection,
    *,
    K_p: Any,
    K_bump: Sequence[Any] | np.ndarray,
    mu: Any,
    energy_target: Any,
    interval: tuple[Any, Any] = (0.9, 1.2),
    conditioning_limit: float = 1.0e12,
) -> AxialEnergySolution:
    """Solve the exact quadratic obtained from (7.34) on ``[.9,1.2]``."""

    K_p_value = _finite(K_p, "K_p")
    K = _vector(K_bump, "K_bump", 2)
    mu_value = _finite(mu, "mu")
    target = _finite(energy_target, "energy_target")
    lower = _finite(interval[0], "interval[0]")
    upper = _finite(interval[1], "interval[1]")
    if not lower < upper or lower < 0.0:
        raise ValueError("energy interval must satisfy 0 <= lower < upper")
    if K_p_value <= 0.0 or np.any(K < 0.0) or mu_value < 0.0:
        raise ValueError("K_p, K_bump and mu must obey the source nonnegative energy ranges")
    u = np.asarray(affine.u, dtype=float)
    v = np.asarray(affine.v, dtype=float)
    quadratic = K_p_value + mu_value * float(np.sum(K * v * v))
    linear = 2.0 * mu_value * float(np.sum(K * u * v))
    constant = mu_value * float(np.sum(K * u * u)) - target
    if quadratic <= 0.0:
        raise ValueError("energy quadratic has nonpositive leading coefficient")
    values = (quadratic * lower * lower + linear * lower + constant, quadratic * upper * upper + linear * upper + constant)
    if values[0] > 1.0e-12 * max(abs(target), 1.0) or values[1] < -1.0e-12 * max(abs(target), 1.0):
        raise ValueError(
            "energy target is unreachable on the requested interval: "
            f"endpoint_values={values}"
        )
    roots = _real_quadratic_roots(quadratic, linear, constant, name="axial energy equation")
    candidates = [root for root in roots if lower - 1.0e-12 <= root <= upper + 1.0e-12]
    positive_slope = [root for root in candidates if 2.0 * quadratic * root + linear > 0.0]
    if len(positive_slope) != 1:
        raise ValueError(
            "energy root is not uniquely reachable with positive source slope: "
            f"roots={roots}, interval={interval}"
        )
    root = float(positive_slope[0])
    derivative = 2.0 * quadratic * root + linear
    condition = abs(derivative)
    well_conditioned = math.isfinite(condition) and condition > 1.0e-12 and condition <= float(conditioning_limit)
    if not well_conditioned:
        raise ValueError(f"energy root is ill-conditioned: derivative={derivative:.17g}")
    return AxialEnergySolution(
        affine=affine,
        K_p=K_p_value,
        K_bump=(float(K[0]), float(K[1])),
        mu=mu_value,
        energy_target=target,
        interval=(lower, upper),
        quadratic_coefficients=(float(quadratic), float(linear), float(constant)),
        roots=tuple(float(value) for value in roots),
        a_p=root,
        c_at_root=affine.coefficients(root),
        endpoint_values=(float(values[0]), float(values[1])),
        derivative_at_root=float(derivative),
        discriminant=float(linear * linear - 4.0 * quadratic * constant),
        reachable=True,
        well_conditioned=True,
    )


def axial_energy_implicit_derivative(
    solution: AxialEnergySolution,
    *,
    d_energy_target: Any = 0.0,
    d_K_p: Any = 0.0,
    d_K_bump: Iterable[Any] = (0.0, 0.0),
    d_mu: Any = 0.0,
    d_u: Iterable[Any] = (0.0, 0.0),
    d_v: Iterable[Any] = (0.0, 0.0),
) -> float:
    """Return the exact implicit derivative of ``a_p`` from (7.34)."""

    d_energy_target = _finite(d_energy_target, "d_energy_target")
    d_K_p = _finite(d_K_p, "d_K_p")
    d_K = _vector(tuple(d_K_bump), "d_K_bump", 2)
    d_u_array = _vector(tuple(d_u), "d_u", 2)
    d_v_array = _vector(tuple(d_v), "d_v", 2)
    d_mu = _finite(d_mu, "d_mu")
    K = np.asarray(solution.K_bump, dtype=float)
    u = np.asarray(solution.affine.u, dtype=float)
    v = np.asarray(solution.affine.v, dtype=float)
    du = d_u_array
    dv = d_v_array
    dqa = d_K_p + d_mu * float(np.sum(K * v * v)) + solution.mu * float(np.sum(d_K * v * v + 2.0 * K * v * dv))
    dqb = 2.0 * (
        d_mu * float(np.sum(K * u * v))
        + solution.mu * float(np.sum(d_K * u * v + K * du * v + K * u * dv))
    )
    dqc = d_mu * float(np.sum(K * u * u)) + solution.mu * float(np.sum(d_K * u * u + 2.0 * K * u * du)) - d_energy_target
    qa, qb, _ = solution.quadratic_coefficients
    numerator = solution.a_p * solution.a_p * dqa + solution.a_p * dqb + dqc
    return float(-numerator / solution.derivative_at_root)


def angular_equation_residuals(solution: AngularBumpSolution) -> tuple[float, float]:
    """Recompute (7.23) independently from a returned solution."""

    return _angular_equations(
        alpha=solution.alpha,
        beta2=solution.beta2,
        rho=solution.rho,
        sigma=solution.sigma,
        q=solution.q,
        d1=solution.d1,
        d2=solution.d2,
    )


def axial_energy_coefficients(solution: AxialEnergySolution) -> tuple[float, float, float]:
    """Expose ``(q_2,q_1,q_0)`` for independent equation replay."""

    return solution.quadratic_coefficients


__all__ = [
    "SOURCE",
    "SOURCE_VERSION",
    "SOURCE_SECTIONS",
    "ANGULAR_EQUATION",
    "AXIAL_LINEAR_EQUATION",
    "AXIAL_ENERGY_EQUATION",
    "PAPER_BUMP_DEFINITION",
    "PAPER_GAMMA_DEFINITION",
    "PaperBumpWeights",
    "PaperAxialWeights",
    "StablePaperAxialWeights",
    "StablePaperAxialAssembly",
    "paper_bump_weights",
    "solve_angular_bumps_from_paper_bump",
    "paper_axial_weights",
    "assemble_paper_axial_affine",
    "stable_paper_axial_weights",
    "assemble_stable_paper_axial_affine",
    "AngularBumpSolution",
    "AngularBumpDerivative",
    "solve_angular_bumps",
    "angular_implicit_derivative",
    "angular_equation_residuals",
    "AxialAffineCorrection",
    "solve_axial_affine",
    "AxialEnergySolution",
    "solve_axial_energy_root",
    "axial_energy_implicit_derivative",
    "axial_energy_coefficients",
]
