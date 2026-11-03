"""Self-contained exact-rational whole-``Z`` enclosure for the core gauge.

This module is deliberately independent of the archived source owners and of
the saved graph fixtures.  It implements only the mathematical leaf

    G(Z) = integral_a^Z (1 - delta*v**2) H(v)
                     / (H(v)**2 + sigma**2) dv,
    H(v) = (1-delta)*v/2 + (1-v**2)*(4*v+j),

with the original exact ``j`` and the directed enclosure
``0 <= delta <= 10**-200``.  Every arithmetic operation used for the
enclosure is a ``Fraction`` operation.  No decimal midpoint is selected for
``delta`` or for the anchor.  Bisection midpoints are partition boundaries,
not source values.

The interval integral is a finite, fixed geometry construction.  Panels are
placed geometrically from the *lower* anchor endpoint at scales of ``sigma``;
each panel uses a rational rectangle bound for the full integrand interval and
is immediately rounded outward to a common ``2**(-384)`` dyadic grid before
accumulation.  Per-panel records are discarded; only cell totals and error
budgets are retained.
The anchor uncertainty is propagated as a separate signed correction.  The
function is proved nonnegative and monotone on each side of the unique
anchor, so endpoint enclosures cover any requested rational ``Z`` cell.

This is a mathematical leaf only.  Source-owner, physical-m15 and closure
claims are intentionally false.

The admitted Z domain is exactly ``[-1,1]``.  "Whole-Z" means the complete
function on that admitted interval; this module makes no real-axis extension.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any


ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)

# The admitted original minimum branch.  It is an exact rational singleton,
# not an interval endpoint or a value recovered from a saved calculation.
ACTUAL_J = Fraction(1, 5721363714816000000000)
ACTUAL_SIGMA = ACTUAL_J / 500

# The upper endpoint is a strict proved cap.  The lower endpoint is zero only
# as an enclosure endpoint; the actual exponential delta is strictly positive.
DELTA_UPPER = Fraction(1, 10**200)
# A fixed dyadic root tolerance keeps the bisection certificate and every
# anchor-panel denominator bounded.  With sigma=j/500, the global derivative
# cap turns this into an explicit anchor-shift error of order 10**(-52).
ROOT_TOLERANCE = Fraction(1, 1 << 256)

# These are deliberately finite and fixed before any cell is evaluated.
ROOT_GEOMETRIC_LEVELS = 96
MAX_ROOT_BISECTION_STEPS = 1024
MAX_PANELS = 512
DYADIC_BITS = 384
DYADIC_DENOMINATOR = 1 << DYADIC_BITS
DYADIC_STEP = Fraction(1, DYADIC_DENOMINATOR)

# Fixed input/result limits prevent hostile Fraction values or a long panel
# chain from producing unbounded arithmetic or JSON integers.
MAX_INPUT_FRACTION_BITS = 4096
# 12000 binary bits are below Python's default 4300-decimal-digit JSON
# guard, while allowing one exact rational quotient before its dyadic round.
MAX_RESULT_FRACTION_BITS = 12000
MAX_DECIMAL_INPUT_DIGITS = 2500


def _bit_length(value: int) -> int:
    return abs(value).bit_length()


def _check_fraction_bits(value: Fraction, name: str, *, input_value: bool = False) -> Fraction:
    limit = MAX_INPUT_FRACTION_BITS if input_value else MAX_RESULT_FRACTION_BITS
    if (_bit_length(value.numerator) > limit
            or _bit_length(value.denominator) > limit):
        raise ArithmeticError(f"{name} exceeds the fixed Fraction bit budget")
    return value


def _coerce_fraction(value: Any, name: str = "value") -> Fraction:
    """Coerce an internal operand under the result (not input) bit budget."""

    if isinstance(value, bool):
        raise TypeError(f"{name} must not be bool")
    if isinstance(value, Fraction):
        return _check_fraction_bits(value, name)
    if isinstance(value, int):
        if _bit_length(value) > MAX_INPUT_FRACTION_BITS:
            raise ArithmeticError(f"{name} exceeds the fixed integer bit budget")
        return _check_fraction_bits(Fraction(value), name)
    if isinstance(value, str):
        if len(value) > MAX_DECIMAL_INPUT_DIGITS:
            raise ArithmeticError(f"{name} exceeds the fixed decimal text budget")
        return _check_fraction_bits(Fraction(value), name)
    raise TypeError(f"{name} must be an exact rational value")


def _fraction(value: Any, name: str = "value") -> Fraction:
    """Convert only exact integer/string/Fraction inputs."""

    if isinstance(value, bool):
        raise TypeError(f"{name} must not be bool")
    if isinstance(value, Fraction):
        return _check_fraction_bits(value, name, input_value=True)
    if isinstance(value, int):
        if _bit_length(value) > MAX_INPUT_FRACTION_BITS:
            raise ArithmeticError(f"{name} exceeds the fixed integer bit budget")
        return Fraction(value)
    if isinstance(value, str):
        if len(value) > MAX_DECIMAL_INPUT_DIGITS:
            raise ArithmeticError(f"{name} exceeds the fixed decimal text budget")
        return _check_fraction_bits(Fraction(value), name, input_value=True)
    raise TypeError(f"{name} must be an exact Fraction, integer, or string")


@dataclass(frozen=True)
class RationalInterval:
    """Closed rational interval with exact outward arithmetic."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        lo = _coerce_fraction(self.lower, "lower")
        hi = _coerce_fraction(self.upper, "upper")
        if lo > hi:
            raise ValueError("interval lower endpoint exceeds upper endpoint")
        object.__setattr__(self, "lower", lo)
        object.__setattr__(self, "upper", hi)

    @classmethod
    def point(cls, value: Any) -> "RationalInterval":
        value = _coerce_fraction(value)
        return cls(value, value)

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower

    @property
    def radius(self) -> Fraction:
        return self.width / 2

    def __add__(self, other: Any) -> "RationalInterval":
        other = _interval(other)
        return type(self)(self.lower + other.lower, self.upper + other.upper)

    __radd__ = __add__

    def __neg__(self) -> "RationalInterval":
        return type(self)(-self.upper, -self.lower)

    def __sub__(self, other: Any) -> "RationalInterval":
        return self + (-_interval(other))

    def __rsub__(self, other: Any) -> "RationalInterval":
        return _interval(other) - self

    def __mul__(self, other: Any) -> "RationalInterval":
        other = _interval(other)
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return type(self)(min(products), max(products))

    __rmul__ = __mul__

    def scale(self, scalar: Any) -> "RationalInterval":
        return self * RationalInterval.point(scalar)

    def square(self) -> "RationalInterval":
        if self.lower <= 0 <= self.upper:
            return type(self)(ZERO, max(self.lower * self.lower,
                                        self.upper * self.upper))
        values = (self.lower * self.lower, self.upper * self.upper)
        return type(self)(min(values), max(values))

    def divide_positive(self, denominator: Any) -> "RationalInterval":
        denominator = _interval(denominator)
        if denominator.lower <= 0:
            raise ArithmeticError("positive denominator proof failed")
        return self * type(self)(ONE / denominator.upper,
                                 ONE / denominator.lower)

    def clamp_lower(self, lower: Any) -> "RationalInterval":
        lower = _fraction(lower)
        # The sign theorem supplies the lower bound.  If the computed upper
        # endpoint misses that theorem, the interval proof has failed; do not
        # silently collapse a disjoint interval to the theorem boundary.
        if self.upper < lower:
            raise ArithmeticError("interval misses the asserted lower bound")
        return type(self)(max(self.lower, lower), self.upper)

    def hull(self, other: Any) -> "RationalInterval":
        other = _interval(other)
        return type(self)(min(self.lower, other.lower),
                          max(self.upper, other.upper))

    def absolute_upper(self) -> Fraction:
        return max(abs(self.lower), abs(self.upper))

    def record(self) -> dict[str, Any]:
        return {
            "lower": {"numerator": self.lower.numerator,
                      "denominator": self.lower.denominator},
            "upper": {"numerator": self.upper.numerator,
                      "denominator": self.upper.denominator},
            "width": {"numerator": self.width.numerator,
                      "denominator": self.width.denominator},
        }


def _interval(value: Any) -> RationalInterval:
    if isinstance(value, RationalInterval):
        return value
    return RationalInterval.point(_fraction(value))


def _floor_dyadic(value: Fraction) -> Fraction:
    """Largest multiple of 2**(-DYADIC_BITS) not exceeding value."""

    value = _check_fraction_bits(value, "dyadic input")
    scaled = value * DYADIC_DENOMINATOR
    _check_fraction_bits(scaled, "dyadic scaled input")
    quotient = scaled.numerator // scaled.denominator
    return _check_fraction_bits(Fraction(quotient, DYADIC_DENOMINATOR),
                                "dyadic floor")


def _ceil_dyadic(value: Fraction) -> Fraction:
    """Smallest multiple of 2**(-DYADIC_BITS) not below value."""

    value = _check_fraction_bits(value, "dyadic input")
    scaled = value * DYADIC_DENOMINATOR
    _check_fraction_bits(scaled, "dyadic scaled input")
    quotient = -((-scaled.numerator) // scaled.denominator)
    return _check_fraction_bits(Fraction(quotient, DYADIC_DENOMINATOR),
                                "dyadic ceiling")


def _dyadic_width_upper(value: RationalInterval) -> Fraction:
    """Outward dyadic upper bound for an interval width without an LCM."""

    return _ceil_dyadic(value.upper) - _floor_dyadic(value.lower)


def _round_interval_dyadic(value: RationalInterval) -> tuple[RationalInterval, Fraction]:
    """Round both endpoints outward and return the two-step budget."""

    return (RationalInterval(_floor_dyadic(value.lower),
                             _ceil_dyadic(value.upper)),
            2 * DYADIC_STEP)


@dataclass(frozen=True)
class GaugeSource:
    """Exact source recipes and the directed delta enclosure."""

    j: Fraction = ACTUAL_J
    delta: RationalInterval = RationalInterval(ZERO, DELTA_UPPER)
    sigma: Fraction = ACTUAL_SIGMA

    def __post_init__(self) -> None:
        j = _fraction(self.j, "j")
        sigma = _fraction(self.sigma, "sigma")
        delta = _interval(self.delta)
        if j <= 0 or sigma <= 0:
            raise ValueError("j and sigma must be positive")
        if delta.lower < 0 or delta.upper >= ONE:
            raise ValueError("delta enclosure must lie in [0,1)")
        if sigma != j / 500:
            raise ValueError("sigma must be exactly j/500")
        object.__setattr__(self, "j", j)
        object.__setattr__(self, "delta", delta)
        object.__setattr__(self, "sigma", sigma)

    def record(self) -> dict[str, Any]:
        return {
            "j": {"numerator": self.j.numerator,
                  "denominator": self.j.denominator},
            "sigma": {"numerator": self.sigma.numerator,
                      "denominator": self.sigma.denominator},
            "delta_interval": self.delta.record(),
            "j_is_actual_exact_fraction": self.j == ACTUAL_J,
            "sigma_is_exact_j_over_500": self.sigma == ACTUAL_SIGMA,
            "delta_lower_is_not_selected_value": True,
            "delta_upper_is_not_selected_value": True,
        }


def actual_source() -> GaugeSource:
    """Return the actual source recipe without selecting a delta value."""

    return GaugeSource()


def delta_bound_proof() -> dict[str, Any]:
    """Return the exact elementary proof used for the delta cap."""

    exp40_lower = 1 + 40 + 40**2 // 2
    exponent_upper = -4 * exp40_lower - 74
    return {
        "exp40_strict_lower_bound": exp40_lower,
        "e_cubed_series_lower_bound": "1 + 3 + 9/2 + 27/6 = 13 > 10",
        "e_cubed_strict_lower_bound": 10,
        "actual_delta_formula": "exp(-4*(exp(40)+11)-30)",
        "delta_exponent_form": "exp(-4*exp(40)-74)",
        "exponent_strict_upper_bound": exponent_upper,
        "comparison": (
            "exp(-4*exp(40)-74) < exp(-3438) "
            "= exp(-3*1146) < 10**(-1146) < 10**(-200)"
        ),
        "delta_interval": RationalInterval(ZERO, DELTA_UPPER).record(),
        "actual_delta_strictly_positive": True,
        "cap_is_not_a_selected_delta": True,
    }


def _h_point(v: Fraction, delta: Fraction, j: Fraction) -> Fraction:
    return ((ONE - delta) * v) / 2 + (ONE - v * v) * (4 * v + j)


def _h_interval(v: RationalInterval, source: GaugeSource) -> RationalInterval:
    one_minus_delta = RationalInterval(ONE - source.delta.upper,
                                        ONE - source.delta.lower)
    first = (one_minus_delta * v).scale(HALF)
    second = (RationalInterval.point(ONE) - v.square()) * (
        v.scale(4) + RationalInterval.point(source.j)
    )
    return first + second


def _l_interval(v: RationalInterval, source: GaugeSource) -> RationalInterval:
    return RationalInterval.point(ONE) - source.delta * v.square()


def gauge_derivative_interval(v: RationalInterval, source: GaugeSource | None = None) -> RationalInterval:
    """Enclose ``G_Z`` on the entire rational input cell ``v``."""

    source = source or actual_source()
    v = _interval(v)
    if v.lower < -ONE or v.upper > ONE:
        raise ValueError("Z cell must lie in [-1,1]")
    h = _h_interval(v, source)
    l = _l_interval(v, source)
    if l.lower < ZERO or l.upper > ONE:
        raise ArithmeticError("0 <= L <= 1 proof failed")
    numerator = l * h
    denominator = h.square() + RationalInterval.point(source.sigma * source.sigma)
    if denominator.lower <= 0:
        raise ArithmeticError("H^2+sigma^2 positivity proof failed")
    raw = numerator.divide_positive(denominator)
    global_bound = ONE / (2 * source.sigma)
    bounded = RationalInterval(max(raw.lower, -global_bound),
                               min(raw.upper, global_bound))
    if bounded.lower > bounded.upper:
        raise ArithmeticError("global |G_Z| <= 1/(2*sigma) intersection failed")
    return bounded


@dataclass(frozen=True)
class RootBracket:
    interval: RationalInterval
    delta: Fraction
    iterations: int

    def record(self) -> dict[str, Any]:
        return {
            "delta": {"numerator": self.delta.numerator,
                      "denominator": self.delta.denominator},
            "interval": self.interval.record(),
            "iterations": self.iterations,
            "bisection_midpoints_are_partition_boundaries_only": True,
        }


@dataclass(frozen=True)
class AnchorCertificate:
    interval: RationalInterval
    delta_lower_root: RootBracket
    delta_upper_root: RootBracket
    root_tolerance: Fraction
    derivative_lower_on_bracket: Fraction

    def record(self) -> dict[str, Any]:
        return {
            "interval": self.interval.record(),
            "delta_lower_root": self.delta_lower_root.record(),
            "delta_upper_root": self.delta_upper_root.record(),
            "root_tolerance": {
                "numerator": self.root_tolerance.numerator,
                "denominator": self.root_tolerance.denominator,
            },
            "derivative_lower_on_[-j,0]": {
                "numerator": self.derivative_lower_on_bracket.numerator,
                "denominator": self.derivative_lower_on_bracket.denominator,
            },
            "unique_root_proved_for_every_delta_in_interval": True,
            "delta_root_ordering_proved": (
                "root(delta_upper) <= root(delta) <= root(delta_lower)"
            ),
            "root_center_not_selected": True,
        }


def _root_bracket(
    delta: Fraction,
    j: Fraction,
    tolerance: Fraction,
    max_steps: int,
) -> RootBracket:
    """Monotone exact bisection for one fixed *endpoint* delta."""

    lo, hi = -j, ZERO
    if _h_point(lo, delta, j) >= 0 or _h_point(hi, delta, j) <= 0:
        raise ArithmeticError("root endpoint sign bracket failed")
    for iteration in range(1, max_steps + 1):
        if hi - lo <= tolerance:
            return RootBracket(RationalInterval(lo, hi), delta, iteration - 1)
        mid = (lo + hi) / 2
        value = _h_point(mid, delta, j)
        if value < 0:
            lo = mid
        elif value > 0:
            hi = mid
        else:
            return RootBracket(RationalInterval(mid, mid), delta, iteration)
    raise ArithmeticError("root bisection budget exhausted")


def anchor_certificate(
    source: GaugeSource | None = None,
    *,
    tolerance: Fraction = ROOT_TOLERANCE,
    max_steps: int = MAX_ROOT_BISECTION_STEPS,
) -> AnchorCertificate:
    """Enclose all anchors induced by the directed delta interval."""

    source = source or actual_source()
    tolerance = _fraction(tolerance, "tolerance")
    if tolerance <= 0:
        raise ValueError("root tolerance must be positive")
    if not isinstance(max_steps, int) or isinstance(max_steps, bool) or max_steps <= 0:
        raise ValueError("max_steps must be a positive integer")
    if max_steps > MAX_ROOT_BISECTION_STEPS:
        raise ValueError("max_steps exceeds the fixed bisection budget")
    derivative_lower = (ONE - source.delta.upper) / 2 + 4 - 12 * source.j * source.j
    if derivative_lower <= 0:
        raise ArithmeticError("H is not certified increasing on [-j,0]")
    lower_root = _root_bracket(source.delta.upper, source.j, tolerance, max_steps)
    upper_root = _root_bracket(source.delta.lower, source.j, tolerance, max_steps)
    interval = RationalInterval(
        min(lower_root.interval.lower, upper_root.interval.lower),
        max(lower_root.interval.upper, upper_root.interval.upper),
    )
    return AnchorCertificate(interval, lower_root, upper_root, tolerance,
                             derivative_lower)


def _validate_anchor_certificate(
    anchor: AnchorCertificate,
    source: GaugeSource,
) -> None:
    """Reject an anchor certificate made for another source or root proof."""

    if not isinstance(anchor, AnchorCertificate):
        raise TypeError("anchor must be an AnchorCertificate")
    if not isinstance(anchor.interval, RationalInterval):
        raise TypeError("anchor interval must be a RationalInterval")
    tolerance = _fraction(anchor.root_tolerance, "anchor.root_tolerance")
    if tolerance <= 0:
        raise ArithmeticError("anchor root tolerance must be positive")
    j = source.j
    if anchor.interval.lower < -j or anchor.interval.upper > ZERO:
        raise ArithmeticError("anchor envelope is outside [-j,0]")

    # Recompute the monotonicity lower bound from this source.  A caller may
    # not smuggle in a positive-looking certificate whose stored lower bound
    # is unrelated to the current j or delta enclosure.
    expected_derivative_lower = (
        (ONE - source.delta.upper) / 2 + 4 - 12 * j * j
    )
    stored_derivative_lower = _fraction(
        anchor.derivative_lower_on_bracket,
        "anchor.derivative_lower_on_bracket",
    )
    if expected_derivative_lower <= ZERO:
        raise ArithmeticError("source does not prove H increasing on [-j,0]")
    if stored_derivative_lower != expected_derivative_lower:
        raise ArithmeticError("anchor monotonicity certificate does not match source")

    endpoint_roots = (
        (anchor.delta_lower_root, source.delta.upper),
        (anchor.delta_upper_root, source.delta.lower),
    )
    root_intervals: list[RationalInterval] = []
    for root, expected_delta in endpoint_roots:
        if not isinstance(root, RootBracket):
            raise TypeError("anchor endpoint root must be a RootBracket")
        if root.delta != expected_delta:
            raise ArithmeticError("anchor endpoint delta is bound to another source")
        if not isinstance(root.interval, RationalInterval):
            raise TypeError("anchor endpoint interval must be a RationalInterval")
        if (not isinstance(root.iterations, int)
                or isinstance(root.iterations, bool)
                or root.iterations < 0
                or root.iterations > MAX_ROOT_BISECTION_STEPS):
            raise ArithmeticError("anchor endpoint bisection metadata is invalid")
        if root.interval.lower < -j or root.interval.upper > ZERO:
            raise ArithmeticError("anchor endpoint root lies outside [-j,0]")
        # Each endpoint bracket must independently prove H(lo) <= 0 <= H(hi).
        h_lower = _h_point(root.interval.lower, root.delta, j)
        h_upper = _h_point(root.interval.upper, root.delta, j)
        if h_lower > ZERO or h_upper < ZERO:
            raise ArithmeticError("anchor endpoint H sign bracket failed")
        root_intervals.append(root.interval)

    root_lower = min(interval.lower for interval in root_intervals)
    root_upper = max(interval.upper for interval in root_intervals)
    if anchor.interval.lower > root_lower or anchor.interval.upper < root_upper:
        raise ArithmeticError("anchor envelope does not contain endpoint roots")


@dataclass(frozen=True)
class SegmentIntegral:
    interval: RationalInterval
    rectangle_spread: Fraction
    absolute_integral_upper: Fraction
    panel_count: int
    rounding_bound: Fraction

    def record(self) -> dict[str, Any]:
        return {
            "interval": self.interval.record(),
            "rectangle_spread": {
                "numerator": self.rectangle_spread.numerator,
                "denominator": self.rectangle_spread.denominator,
            },
            "absolute_integral_upper": {
                "numerator": self.absolute_integral_upper.numerator,
                "denominator": self.absolute_integral_upper.denominator,
            },
            "panel_count": self.panel_count,
            "rounding_bound": {
                "numerator": self.rounding_bound.numerator,
                "denominator": self.rounding_bound.denominator,
            },
        }


def _geometry_edges(
    left: Fraction,
    right: Fraction,
    anchor_lower: Fraction,
    anchor_upper: Fraction,
    sigma: Fraction,
    levels: int,
) -> tuple[Fraction, ...]:
    """Fixed finite geometry, anchored at ``anchor_lower``."""

    values = {left, right, anchor_lower, anchor_upper, ZERO}
    for level in range(levels + 1):
        step = sigma * (1 << level)
        values.add(anchor_lower - step)
        values.add(anchor_lower + step)
    values = sorted(value for value in values if left <= value <= right)
    if not values or values[0] != left:
        values.insert(0, left)
    if values[-1] != right:
        values.append(right)
    return tuple(dict.fromkeys(values))


def _panel_contribution(
    left: Fraction,
    right: Fraction,
    source: GaugeSource,
) -> tuple[RationalInterval, RationalInterval]:
    """Return exact derivative and exact integral contribution for one panel."""

    if left >= right:
        raise ValueError("panel endpoints must be strictly ordered")
    derivative = gauge_derivative_interval(RationalInterval(left, right), source)
    width = right - left
    contribution = RationalInterval(derivative.lower * width,
                                    derivative.upper * width)
    return derivative, contribution


def _integrate_between(
    left: Fraction,
    right: Fraction,
    source: GaugeSource,
    anchor: AnchorCertificate,
    levels: int,
    max_panels: int,
) -> SegmentIntegral:
    if left == right:
        return SegmentIntegral(RationalInterval.point(ZERO), ZERO, ZERO, 0, ZERO)
    orientation = 1
    if left > right:
        left, right = right, left
        orientation = -1
    edges = _geometry_edges(left, right, anchor.interval.lower,
                            anchor.interval.upper, source.sigma, levels)
    if len(edges) - 1 > max_panels:
        raise ArithmeticError("fixed panel budget exceeded")
    panel_count = len(edges) - 1
    total = RationalInterval.point(ZERO)
    spread = ZERO
    absolute_integral_upper = ZERO
    for a, b in zip(edges, edges[1:]):
        derivative, exact_contribution = _panel_contribution(a, b, source)
        rounded_contribution, _ = _round_interval_dyadic(exact_contribution)
        total = total + rounded_contribution
        # These are already outward-rounded dyadic diagnostics, so their
        # denominators remain exactly the fixed common dyadic denominator.
        spread += _dyadic_width_upper(exact_contribution)
        absolute_integral_upper += _ceil_dyadic(
            derivative.absolute_upper() * (b - a)
        )
    if orientation < 0:
        total = -total
    return SegmentIntegral(
        total,
        spread,
        absolute_integral_upper,
        panel_count,
        panel_count * 2 * DYADIC_STEP,
    )


def _derivative_hull_between(
    left: Fraction,
    right: Fraction,
    source: GaugeSource,
    anchor: AnchorCertificate,
    levels: int,
    max_panels: int,
) -> RationalInterval:
    """Hull panelwise G_Z intervals instead of one coarse H rectangle."""

    if left == right:
        return gauge_derivative_interval(RationalInterval.point(left), source)
    if left > right:
        left, right = right, left
    edges = _geometry_edges(left, right, anchor.interval.lower,
                            anchor.interval.upper, source.sigma, levels)
    if len(edges) - 1 > max_panels:
        raise ArithmeticError("fixed derivative panel budget exceeded")
    hull = None
    for a, b in zip(edges, edges[1:]):
        derivative, _ = _panel_contribution(a, b, source)
        hull = derivative if hull is None else hull.hull(derivative)
    if hull is None:
        raise ArithmeticError("empty derivative panel partition")
    return hull


def _anchor_shift(
    source: GaugeSource,
    anchor: AnchorCertificate,
    levels: int,
    max_panels: int,
) -> tuple[RationalInterval, Fraction, SegmentIntegral]:
    """Bound ``integral_{a_lo}^a g`` without selecting ``a``."""

    lower = anchor.interval.lower
    upper = anchor.interval.upper
    segment = _integrate_between(lower, upper, source, anchor, levels, max_panels)
    absolute_bound = segment.absolute_integral_upper
    # Every admissible actual root lies to the right of a_lo; on [a_lo,a]
    # H<0, hence g<0.  The interval below retains only this proved sign.
    return RationalInterval(-absolute_bound, ZERO), absolute_bound, segment


def _endpoint_value(
    z: Fraction,
    source: GaugeSource,
    anchor: AnchorCertificate,
    levels: int,
    max_panels: int,
    anchor_shift: RationalInterval,
) -> tuple[RationalInterval, SegmentIntegral]:
    base = _integrate_between(anchor.interval.lower, z, source, anchor,
                              levels, max_panels)
    return (base.interval - anchor_shift).clamp_lower(ZERO), base


def _cell_value_interval(
    left: Fraction,
    right: Fraction,
    left_value: RationalInterval,
    right_value: RationalInterval,
    anchor: AnchorCertificate,
) -> RationalInterval:
    if left <= anchor.interval.upper and right >= anchor.interval.lower:
        # The unique root is in, or can touch, this Z cell.  G(a)=0 and G is
        # nonnegative, so the lower endpoint is exact and the upper endpoint
        # is the larger of the two monotone-side endpoint bounds.
        return RationalInterval(ZERO,
                                max(left_value.upper, right_value.upper))
    lower = max(ZERO, min(left_value.lower, right_value.lower))
    upper = max(left_value.upper, right_value.upper)
    return RationalInterval(lower, upper)


def evaluate_wholez_cell(
    z_left: Any,
    z_right: Any,
    *,
    source: GaugeSource | None = None,
    anchor: AnchorCertificate | None = None,
    root_tolerance: Fraction = ROOT_TOLERANCE,
    geometry_levels: int = ROOT_GEOMETRIC_LEVELS,
    max_panels: int = MAX_PANELS,
) -> dict[str, Any]:
    """Return exact rational enclosures for ``G`` and ``G_Z`` on a Z cell.

    ``z_left`` and ``z_right`` are rational cell endpoints, not sampling
    points.  The returned ``G_interval`` encloses every function value on the
    complete cell in the admitted domain ``Z in [-1,1]``, while
    ``G_Z_interval`` encloses the derivative on that same cell.  No real-axis
    extension is claimed.
    """

    source = source or actual_source()
    root_tolerance = _fraction(root_tolerance, "root_tolerance")
    if root_tolerance <= 0:
        raise ValueError("root tolerance must be positive")
    left, right = _fraction(z_left, "z_left"), _fraction(z_right, "z_right")
    if left > right:
        raise ValueError("Z cell endpoints must be ordered")
    if left < -ONE or right > ONE:
        raise ValueError("Z cell must lie in [-1,1]")
    if not isinstance(geometry_levels, int) or isinstance(geometry_levels, bool) or geometry_levels < 0:
        raise ValueError("geometry_levels must be a nonnegative integer")
    if geometry_levels > ROOT_GEOMETRIC_LEVELS:
        raise ValueError("geometry_levels exceeds the fixed geometry budget")
    if (not isinstance(max_panels, int) or isinstance(max_panels, bool)
            or max_panels <= 0 or max_panels > MAX_PANELS):
        raise ValueError("max_panels exceeds the fixed panel budget")
    if anchor is None:
        anchor = anchor_certificate(source, tolerance=root_tolerance)
    else:
        _validate_anchor_certificate(anchor, source)
    shift, anchor_bound, anchor_segment = _anchor_shift(
        source, anchor, geometry_levels, max_panels
    )
    left_value, left_segment = _endpoint_value(
        left, source, anchor, geometry_levels, max_panels, shift
    )
    right_value, right_segment = _endpoint_value(
        right, source, anchor, geometry_levels, max_panels, shift
    )
    value = _cell_value_interval(left, right, left_value, right_value, anchor)
    derivative_raw = _derivative_hull_between(
        left, right, source, anchor, geometry_levels, max_panels
    )
    # Keep the externally returned derivative on the same bounded dyadic
    # grid as the integral totals.  The exact panel hull is used for the
    # proof; only its serialized cell value is rounded outward here.  This
    # avoids retaining a cell-wide Fraction with a fresh denominator while
    # keeping the output enclosure honest.
    derivative, derivative_rounding = _round_interval_dyadic(derivative_raw)
    integration_spread = (left_segment.rectangle_spread
                          + right_segment.rectangle_spread
                          + anchor_segment.rectangle_spread)
    total_panels = (left_segment.panel_count + right_segment.panel_count
                    + anchor_segment.panel_count)
    rounding_bound = (left_segment.rounding_bound
                      + right_segment.rounding_bound
                      + anchor_segment.rounding_bound)
    return {
        "schema": "original_core_gauge_wholez_v1",
        "admitted_Z_domain": RationalInterval(-ONE, ONE).record(),
        "Z_cell": RationalInterval(left, right).record(),
        "G_interval": value.record(),
        "G_Z_interval": derivative.record(),
        "endpoint_values": {"left": left_value.record(),
                             "right": right_value.record()},
        "anchor": anchor.record(),
        "anchor_shift_interval": shift.record(),
        "panel_budget": {
            "geometry_anchor": "anchor_interval.lower",
            "geometry_levels": geometry_levels,
            "max_panels": max_panels,
            "finite_geometry_declared_before_evaluation": True,
            "anchor_panel_count": anchor_segment.panel_count,
            "left_panel_count": left_segment.panel_count,
            "right_panel_count": right_segment.panel_count,
            "total_panel_count": total_panels,
        },
        "error_budget": {
            "integration_truncation": RationalInterval(ZERO, integration_spread).record(),
            "anchor_uncertainty": RationalInterval(ZERO, anchor_bound).record(),
            "root_interval_width": {
                "numerator": anchor.interval.width.numerator,
                "denominator": anchor.interval.width.denominator,
            },
            "rounding": RationalInterval(ZERO, rounding_bound).record(),
            "G_Z_rounding": RationalInterval(ZERO, derivative_rounding).record(),
            "dyadic_bits": DYADIC_BITS,
            "dyadic_step": {
                "numerator": DYADIC_STEP.numerator,
                "denominator": DYADIC_STEP.denominator,
            },
            "rounding_bound_formula": "2 * total_panel_count * 2**(-384)",
            "rounding_bound_is_at_most_two_steps_per_panel": True,
            "G_Z_rounding_bound_formula": "2 * 2**(-384), separate from integral panel rounding",
            "all_endpoints_exact_fraction": True,
        },
        "proofs": {
            "delta_bound": delta_bound_proof(),
            "H_increasing_on_anchor_bracket": True,
            "unique_anchor_in_[-j,0]": True,
            "H_negative_left_of_anchor_and_positive_right": True,
            "L_strictly_positive": True,
            "H_squared_plus_sigma_squared_strictly_positive": True,
            "G_nonnegative": True,
            "G_monotone_on_each_side_of_anchor": True,
            "endpoint_hull_covers_complete_Z_cell": True,
            "global_G_Z_bound": "abs(G_Z) <= 1/(2*sigma)",
            "global_G_Z_bound_intersected_with_cell_interval": True,
            "admitted_domain_only_Z_in_[-1,1]": True,
            "no_real_axis_extension": True,
        },
        "source": source.record(),
        "source_owner_proved": False,
        "physical_m15": False,
        "closure": False,
        "old_fixture_imported": False,
        "ancestor_recomputed": False,
        "delta_center_selected": False,
        "anchor_center_selected": False,
        "midpoint_sample_used": False,
    }


class ActualCoreGaugeWholeZ:
    """Reusable owner for multiple rational Z cells."""

    def __init__(
        self,
        *,
        source: GaugeSource | None = None,
        root_tolerance: Fraction = ROOT_TOLERANCE,
        geometry_levels: int = ROOT_GEOMETRIC_LEVELS,
        max_panels: int = MAX_PANELS,
    ) -> None:
        self.source = source or actual_source()
        self.root_tolerance = _fraction(root_tolerance, "root_tolerance")
        if not isinstance(geometry_levels, int) or isinstance(geometry_levels, bool):
            raise ValueError("geometry_levels must be an integer")
        if geometry_levels < 0 or geometry_levels > ROOT_GEOMETRIC_LEVELS:
            raise ValueError("geometry_levels exceeds the fixed geometry budget")
        if (not isinstance(max_panels, int) or isinstance(max_panels, bool)
                or max_panels <= 0 or max_panels > MAX_PANELS):
            raise ValueError("max_panels exceeds the fixed panel budget")
        self.geometry_levels = geometry_levels
        self.max_panels = max_panels
        self.anchor = anchor_certificate(self.source, tolerance=self.root_tolerance)

    def evaluate(self, z_left: Any, z_right: Any | None = None) -> dict[str, Any]:
        if z_right is None:
            if not isinstance(z_left, (tuple, list)) or len(z_left) != 2:
                raise TypeError("evaluate expects (z_left,z_right) or two endpoints")
            z_left, z_right = z_left
        return evaluate_wholez_cell(
            z_left,
            z_right,
            source=self.source,
            anchor=self.anchor,
            root_tolerance=self.root_tolerance,
            geometry_levels=self.geometry_levels,
            max_panels=self.max_panels,
        )

    def evaluate_cell(self, z_left: Any, z_right: Any | None = None) -> dict[str, Any]:
        return self.evaluate(z_left, z_right)

    def derivative(self, z_left: Any, z_right: Any | None = None) -> RationalInterval:
        if z_right is None:
            if not isinstance(z_left, (tuple, list)) or len(z_left) != 2:
                raise TypeError("derivative expects (z_left,z_right) or two endpoints")
            z_left, z_right = z_left
        left, right = _fraction(z_left, "z_left"), _fraction(z_right, "z_right")
        if left > right:
            raise ValueError("Z cell endpoints must be ordered")
        return gauge_derivative_interval(RationalInterval(left, right), self.source)


def source_contract() -> dict[str, Any]:
    """Static contract for callers before any real cell is evaluated."""

    source = actual_source()
    return {
        "schema": "original_core_gauge_wholez_v1_source_contract",
        "source": source.record(),
        "delta_bound_proof": delta_bound_proof(),
        "domain": RationalInterval(-ONE, ONE).record(),
        "domain_statement": "admitted Z in [-1,1]; no real-axis extension",
        "formula": {
            "H": "((1-delta)/2)*v + (1-v^2)*(4*v+j)",
            "L": "1-delta*v^2",
            "G_Z": "L*H/(H^2+(j/500)^2)",
            "G": "integral from the unique H-root a to Z",
        },
        "source_owner_proved": False,
        "physical_m15": False,
        "closure": False,
    }


__all__ = [
    "ACTUAL_J",
    "ACTUAL_SIGMA",
    "DELTA_UPPER",
    "ROOT_TOLERANCE",
    "RationalInterval",
    "GaugeSource",
    "RootBracket",
    "AnchorCertificate",
    "ActualCoreGaugeWholeZ",
    "actual_source",
    "anchor_certificate",
    "delta_bound_proof",
    "gauge_derivative_interval",
    "evaluate_wholez_cell",
    "source_contract",
]
