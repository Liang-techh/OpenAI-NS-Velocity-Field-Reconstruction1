"""Independent rational interval leaf for the original flat ``sigma``.

This module is deliberately a new interval family.  It does not import the
point-valued sigma implementation, a graph/fixture producer, or any of the
older project modules.  Its public result encloses only

    sigma(x) and d sigma(x) / d x,

where ``x`` is an ordinary real input coordinate.  It is not a derivative in
the slow ``Z`` coordinate and it does not establish an owner or source
binding.

The source function is

    sigma(x) = 0,                              x <= 0,
               a / (a + b),                    0 < x < 1,
               1,                              x >= 1,

    a = exp(-1/x**2),       b = exp(-1/(1-x)**2).

For the open interval, its ordinary derivative is

    sigma'(x) = sigma(x) * (1 - sigma(x))
                * (2/x**3 + 2/(1-x)**3).

All interval endpoints are ``fractions.Fraction`` instances.  The only
transcendental operation is enclosed by a proved rational construction:
for ``0 <= y <= 1`` the alternating Taylor series for ``exp(-y)`` is used;
for a larger exponent it is range-reduced by powers of two.  Very large
exponents use the elementary positive-series inequality

    exp(t) >= t**p / p!,

which gives the rational upper bound ``exp(-t) <= p! / t**p``.  The fallback
lower bound is zero; the bounded exact power is constructed only after the
same numerator/denominator budget check as every other intermediate.  No
midpoint or floating-point value is used as a function value.

The derivative enclosure is intentionally conservative near the flat
endpoints.  On ``0 < x <= 1/4`` the inequalities

    t = x**(-2) - (1-x)**(-2) >= (8/9) x**(-2),
    2 x**(-3) + 2 (1-x)**(-3) <= (56/27) x**(-3)

and ``exp(-t) <= 8!/t**8`` give a finite rational cap.  The right endpoint
uses reflection.  The central pieces use monotonic endpoint bounds.  Thus
an interval containing 0, 1, or 1/2 never divides by an endpoint distance.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


SCHEMA = "n1_sigma_interval_rational_v1"
COMPONENT = "original_flat_sigma_interval_rational"
SOURCE_DEFINITION = (
    "0 for x<=0; 1 for x>=1; "
    "exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside"
)

# These constants are algorithm controls, not numerical data from an earlier
# computation.  The defaults give exact rational endpoints and a useful
# enclosure while keeping repeated interval squaring bounded.  The hard
# maxima are part of the fail-closed contract: a caller cannot silently ask
# this leaf to allocate an unbounded rational, Taylor series, reduction chain,
# or fallback exponent.
DEFAULT_SERIES_TERMS = 24
DEFAULT_MAX_REDUCTION_STEPS = 8
DEFAULT_FALLBACK_POWER = 8
_ENDPOINT_CAP_POWER = 8

MAX_INPUT_NUMERATOR_BITS = 4096
MAX_INPUT_DENOMINATOR_BITS = 4096
MAX_INTERMEDIATE_NUMERATOR_BITS = 65536
MAX_INTERMEDIATE_DENOMINATOR_BITS = 65536
MAX_SERIES_TERMS = 128
MAX_REDUCTION_STEPS = 16
MAX_FALLBACK_POWER = 32


def _factorial(n: int) -> int:
    result = 1
    for index in range(2, n + 1):
        result *= index
    return result


def _require_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    return value


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


class ResourceLimitError(ArithmeticError):
    """Raised when exact rational work exceeds the bounded leaf budget."""


@dataclass(frozen=True)
class ResourceLimits:
    """Hard resource ceilings for exact rational evaluation."""

    max_input_numerator_bits: int = MAX_INPUT_NUMERATOR_BITS
    max_input_denominator_bits: int = MAX_INPUT_DENOMINATOR_BITS
    max_intermediate_numerator_bits: int = MAX_INTERMEDIATE_NUMERATOR_BITS
    max_intermediate_denominator_bits: int = MAX_INTERMEDIATE_DENOMINATOR_BITS
    max_series_terms: int = MAX_SERIES_TERMS
    max_reduction_steps: int = MAX_REDUCTION_STEPS
    max_fallback_power: int = MAX_FALLBACK_POWER

    def __post_init__(self) -> None:
        fields = (
            ("max_input_numerator_bits", self.max_input_numerator_bits,
             MAX_INPUT_NUMERATOR_BITS),
            ("max_input_denominator_bits", self.max_input_denominator_bits,
             MAX_INPUT_DENOMINATOR_BITS),
            ("max_intermediate_numerator_bits",
             self.max_intermediate_numerator_bits,
             MAX_INTERMEDIATE_NUMERATOR_BITS),
            ("max_intermediate_denominator_bits",
             self.max_intermediate_denominator_bits,
             MAX_INTERMEDIATE_DENOMINATOR_BITS),
            ("max_series_terms", self.max_series_terms, MAX_SERIES_TERMS),
            ("max_reduction_steps", self.max_reduction_steps,
             MAX_REDUCTION_STEPS),
            ("max_fallback_power", self.max_fallback_power,
             MAX_FALLBACK_POWER),
        )
        for name, value, hard_max in fields:
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
            if value > hard_max:
                raise ResourceLimitError(f"{name} exceeds hard resource maximum")
        if self.max_input_numerator_bits > self.max_intermediate_numerator_bits:
            raise ValueError("input numerator budget exceeds intermediate budget")
        if self.max_input_denominator_bits > self.max_intermediate_denominator_bits:
            raise ValueError("input denominator budget exceeds intermediate budget")


DEFAULT_RESOURCE_LIMITS = ResourceLimits()


def _check_fraction_bits(
    value: Fraction,
    limits: ResourceLimits,
    name: str,
    *,
    input_value: bool = False,
) -> Fraction:
    numerator_bits = abs(value.numerator).bit_length()
    denominator_bits = value.denominator.bit_length()
    numerator_limit = (
        limits.max_input_numerator_bits
        if input_value
        else limits.max_intermediate_numerator_bits
    )
    denominator_limit = (
        limits.max_input_denominator_bits
        if input_value
        else limits.max_intermediate_denominator_bits
    )
    if numerator_bits > numerator_limit:
        raise ResourceLimitError(f"{name} numerator exceeds rational bit budget")
    if denominator_bits > denominator_limit:
        raise ResourceLimitError(f"{name} denominator exceeds rational bit budget")
    return value


def _checked_add(
    left: Fraction, right: Fraction, limits: ResourceLimits, name: str
) -> Fraction:
    if (left.numerator.bit_length() + right.denominator.bit_length() + 1
            > limits.max_intermediate_numerator_bits):
        raise ResourceLimitError(f"{name} numerator addition exceeds bit budget")
    if (right.numerator.bit_length() + left.denominator.bit_length() + 1
            > limits.max_intermediate_numerator_bits):
        raise ResourceLimitError(f"{name} numerator addition exceeds bit budget")
    if (left.denominator.bit_length() + right.denominator.bit_length()
            > limits.max_intermediate_denominator_bits):
        raise ResourceLimitError(f"{name} denominator addition exceeds bit budget")
    return _check_fraction_bits(left + right, limits, name)


def _checked_sub(
    left: Fraction, right: Fraction, limits: ResourceLimits, name: str
) -> Fraction:
    return _checked_add(left, -right, limits, name)


def _checked_mul(
    left: Fraction, right: Fraction, limits: ResourceLimits, name: str
) -> Fraction:
    if abs(left) == 1:
        numerator_bound = right.numerator.bit_length()
    elif abs(right) == 1:
        numerator_bound = left.numerator.bit_length()
    else:
        numerator_bound = left.numerator.bit_length() + right.numerator.bit_length()
    if numerator_bound > limits.max_intermediate_numerator_bits:
        raise ResourceLimitError(f"{name} numerator multiplication exceeds bit budget")
    if left.denominator == 1 or right.denominator == 1:
        denominator_bound = max(left.denominator.bit_length(), right.denominator.bit_length())
    else:
        denominator_bound = left.denominator.bit_length() + right.denominator.bit_length()
    if denominator_bound > limits.max_intermediate_denominator_bits:
        raise ResourceLimitError(f"{name} denominator multiplication exceeds bit budget")
    return _check_fraction_bits(left * right, limits, name)


def _checked_div(
    left: Fraction, right: Fraction, limits: ResourceLimits, name: str
) -> Fraction:
    if right == 0:
        raise ZeroDivisionError(f"{name} divides by zero")
    if abs(left) == 1 or right.denominator == 1:
        numerator_bound = max(left.numerator.bit_length(), right.denominator.bit_length())
    else:
        numerator_bound = left.numerator.bit_length() + right.denominator.bit_length()
    if numerator_bound > limits.max_intermediate_numerator_bits:
        raise ResourceLimitError(f"{name} numerator division exceeds bit budget")
    if left.denominator == 1 or abs(right) == 1:
        denominator_bound = max(left.denominator.bit_length(), right.numerator.bit_length())
    else:
        denominator_bound = left.denominator.bit_length() + right.numerator.bit_length()
    if denominator_bound > limits.max_intermediate_denominator_bits:
        raise ResourceLimitError(f"{name} denominator division exceeds bit budget")
    return _check_fraction_bits(left / right, limits, name)


def _checked_pow(
    value: Fraction, exponent: int, limits: ResourceLimits, name: str
) -> Fraction:
    if isinstance(exponent, bool) or not isinstance(exponent, int) or exponent < 0:
        raise ValueError(f"{name} exponent must be a nonnegative integer")
    result = _check_fraction_bits(Fraction(1), limits, f"{name} identity")
    base = _check_fraction_bits(value, limits, f"{name} base")
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _checked_mul(result, base, limits, name)
        remaining //= 2
        if remaining:
            base = _checked_mul(base, base, limits, name)
    return result


@dataclass(frozen=True)
class RationalInterval:
    """A closed interval with exact rational endpoints."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction):
            raise TypeError("lower must be a Fraction")
        if not isinstance(self.upper, Fraction):
            raise TypeError("upper must be a Fraction")
        if abs(self.lower.numerator).bit_length() > MAX_INTERMEDIATE_NUMERATOR_BITS:
            raise ResourceLimitError("interval lower numerator exceeds hard bit budget")
        if self.lower.denominator.bit_length() > MAX_INTERMEDIATE_DENOMINATOR_BITS:
            raise ResourceLimitError("interval lower denominator exceeds hard bit budget")
        if abs(self.upper.numerator).bit_length() > MAX_INTERMEDIATE_NUMERATOR_BITS:
            raise ResourceLimitError("interval upper numerator exceeds hard bit budget")
        if self.upper.denominator.bit_length() > MAX_INTERMEDIATE_DENOMINATOR_BITS:
            raise ResourceLimitError("interval upper denominator exceeds hard bit budget")
        if self.lower > self.upper:
            raise ValueError("interval bounds must be ordered")

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower

    def contains(self, value: Fraction) -> bool:
        value = _require_fraction(value, "value")
        return self.lower <= value <= self.upper


@dataclass(frozen=True)
class SigmaIntervalResult:
    """Certified ordinary-x sigma and derivative enclosure.

    ``sigma`` encloses ``{sigma(x): x in input_interval}`` and
    ``d_sigma_dx`` encloses ``{d sigma(x) / d x: x in input_interval}``.
    ``ordinary_x_derivative`` is explicit to prevent a caller from treating
    this as an ordinary-Z jet.  ``source_identity_verified`` is always false
    here, and ``owner_binding`` is always false here.
    """

    input_interval: RationalInterval
    sigma: RationalInterval
    d_sigma_dx: RationalInterval
    certified: bool = True
    ordinary_x_derivative: bool = True
    z_derivative: bool = False
    owner_binding: bool = False
    source_definition: str = SOURCE_DEFINITION
    component: str = COMPONENT
    algorithm_family: str = SCHEMA
    source_identity_verified: bool = False
    used_legacy_modules: bool = False

    @property
    def derivative(self) -> RationalInterval:
        """Compatibility spelling for the ordinary-x derivative only."""

        return self.d_sigma_dx


def _interval_square_nonnegative(
    interval: RationalInterval, limits: ResourceLimits
) -> RationalInterval:
    if interval.lower < 0:
        raise ValueError("nonnegative interval required")
    lower = _checked_mul(interval.lower, interval.lower, limits, "interval square")
    upper = _checked_mul(interval.upper, interval.upper, limits, "interval square")
    return RationalInterval(lower, upper)


def _exp_unit_bounds(
    value: Fraction, terms: int, limits: ResourceLimits
) -> RationalInterval:
    """Enclose ``exp(-value)`` for ``0 <= value <= 1``.

    If ``S_N`` is the alternating partial sum through degree ``N``, the
    alternating-series theorem places the exact value between ``S_N`` and
    ``S_N + (-1) ** (N + 1) * value ** (N + 1) / (N + 1)!``.
    """

    if not Fraction(0) <= value <= Fraction(1):
        raise ValueError("unit exponential argument must lie in [0, 1]")
    value = _check_fraction_bits(value, limits, "unit exponential argument")
    total = Fraction(1)
    term = Fraction(1)
    for degree in range(1, terms + 1):
        term = _checked_mul(term, value, limits, "exponential term")
        term = _checked_div(term, Fraction(degree), limits, "exponential term")
        if degree % 2:
            total = _checked_sub(total, term, limits, "exponential sum")
        else:
            total = _checked_add(total, term, limits, "exponential sum")
    next_term = _checked_mul(term, value, limits, "exponential next term")
    next_term = _checked_div(
        next_term, Fraction(terms + 1), limits, "exponential next term"
    )
    if terms % 2:
        lower, upper = total, _checked_add(
            total, next_term, limits, "exponential upper"
        )
    else:
        lower, upper = _checked_sub(
            total, next_term, limits, "exponential lower"
        ), total
    if lower < 0:
        lower = Fraction(0)
    if upper > 1:
        upper = Fraction(1)
    return RationalInterval(lower, upper)


def _rational_exp_negative_bounded(
    value: Fraction,
    limits: ResourceLimits,
    series_terms: int,
    max_reduction_steps: int,
    fallback_power: int,
    *,
    input_value: bool,
) -> RationalInterval:
    """Bound ``exp(-value)`` after the public resource checks."""

    value = _check_fraction_bits(
        value, limits, "exponential input", input_value=input_value
    )
    if value < 0:
        raise ValueError("value must be nonnegative")
    if value == 0:
        return RationalInterval(Fraction(1), Fraction(1))

    reduced = value
    squarings = 0
    while reduced > 1 and squarings < max_reduction_steps:
        reduced = _checked_div(reduced, Fraction(2), limits, "exponential reduction")
        squarings += 1
    if reduced > 1:
        powered = _checked_pow(value, fallback_power, limits, "exponential fallback")
        cap = _checked_div(
            Fraction(_factorial(fallback_power)),
            powered,
            limits,
            "exponential fallback",
        )
        return RationalInterval(Fraction(0), min(Fraction(1), cap))

    interval = _exp_unit_bounds(reduced, series_terms, limits)
    for _ in range(squarings):
        interval = _interval_square_nonnegative(interval, limits)
    return interval


def rational_exp_negative(
    value: Fraction,
    *,
    series_terms: int = DEFAULT_SERIES_TERMS,
    max_reduction_steps: int = DEFAULT_MAX_REDUCTION_STEPS,
    fallback_power: int = DEFAULT_FALLBACK_POWER,
    resource_limits: ResourceLimits = DEFAULT_RESOURCE_LIMITS,
) -> RationalInterval:
    """Return a proved rational enclosure for ``exp(-value)``.

    ``value`` must be a nonnegative exact rational.  For a reduced argument
    ``u = value / 2**s <= 1``, the alternating interval for ``exp(-u)`` is
    squared ``s`` times.  The fallback uses the positive Taylor inequality
    ``exp(value) >= value**p / p!``.  Every power, series term, and reduction
    is checked against ``resource_limits`` before and after construction;
    exceeding a hard bound raises :class:`ResourceLimitError`.
    """

    value = _require_fraction(value, "value")
    if not isinstance(resource_limits, ResourceLimits):
        raise TypeError("resource_limits must be ResourceLimits")
    series_terms = _positive_integer(series_terms, "series_terms")
    max_reduction_steps = _positive_integer(
        max_reduction_steps, "max_reduction_steps"
    )
    fallback_power = _positive_integer(fallback_power, "fallback_power")
    if series_terms > resource_limits.max_series_terms:
        raise ResourceLimitError("series_terms exceeds resource budget")
    if max_reduction_steps > resource_limits.max_reduction_steps:
        raise ResourceLimitError("max_reduction_steps exceeds resource budget")
    if fallback_power > resource_limits.max_fallback_power:
        raise ResourceLimitError("fallback_power exceeds resource budget")
    return _rational_exp_negative_bounded(
        value,
        resource_limits,
        series_terms,
        max_reduction_steps,
        fallback_power,
        input_value=True,
    )


def _odds_from_exp(
    value: Fraction,
    limits: ResourceLimits,
    *,
    series_terms: int,
    max_reduction_steps: int,
    fallback_power: int,
) -> RationalInterval:
    """Enclose ``q = exp(-value)`` for a left-half point ``x <= 1/2``."""

    return _rational_exp_negative_bounded(
        value,
        limits,
        series_terms,
        max_reduction_steps,
        fallback_power,
        input_value=False,
    )


def _sigma_point(
    x: Fraction,
    *,
    limits: ResourceLimits,
    series_terms: int,
    max_reduction_steps: int,
    fallback_power: int,
) -> tuple[RationalInterval, RationalInterval]:
    """Point-input enclosure used only at rational interval endpoints."""

    if x <= 0:
        zero = RationalInterval(Fraction(0), Fraction(0))
        return zero, zero
    if x >= 1:
        one = RationalInterval(Fraction(1), Fraction(1))
        zero = RationalInterval(Fraction(0), Fraction(0))
        return one, zero
    if x == Fraction(1, 2):
        return (
            RationalInterval(Fraction(1, 2), Fraction(1, 2)),
            RationalInterval(Fraction(8), Fraction(8)),
        )
    if x > Fraction(1, 2):
        reflected_sigma, reflected_derivative = _sigma_point(
            _checked_sub(Fraction(1), x, limits, "sigma reflection"),
            limits=limits,
            series_terms=series_terms,
            max_reduction_steps=max_reduction_steps,
            fallback_power=fallback_power,
        )
        return (
            RationalInterval(
                _checked_sub(
                    Fraction(1), reflected_sigma.upper, limits, "sigma reflection"
                ),
                _checked_sub(
                    Fraction(1), reflected_sigma.lower, limits, "sigma reflection"
                ),
            ),
            reflected_derivative,
        )

    x_square = _checked_pow(x, 2, limits, "sigma x square")
    reflected_x = _checked_sub(Fraction(1), x, limits, "sigma complement")
    reflected_square = _checked_pow(reflected_x, 2, limits, "sigma complement square")
    left_inverse = _checked_div(Fraction(1), x_square, limits, "sigma exponent")
    right_inverse = _checked_div(
        Fraction(1), reflected_square, limits, "sigma exponent"
    )
    t = _checked_sub(left_inverse, right_inverse, limits, "sigma exponent")
    odds = _odds_from_exp(
        t,
        limits,
        series_terms=series_terms,
        max_reduction_steps=max_reduction_steps,
        fallback_power=fallback_power,
    )
    one_plus_lower = _checked_add(
        Fraction(1), odds.lower, limits, "sigma denominator"
    )
    one_plus_upper = _checked_add(
        Fraction(1), odds.upper, limits, "sigma denominator"
    )
    sigma = RationalInterval(
        _checked_div(odds.lower, one_plus_lower, limits, "sigma value"),
        _checked_div(odds.upper, one_plus_upper, limits, "sigma value"),
    )
    # q/(1+q)^2 is increasing for 0 <= q <= 1.  The source derivative is
    # q/(1+q)^2 * (2/x^3 + 2/(1-x)^3), which avoids an interval division by a
    # vanishing endpoint distance.
    lower_square = _checked_pow(one_plus_lower, 2, limits, "sigma derivative factor")
    upper_square = _checked_pow(one_plus_upper, 2, limits, "sigma derivative factor")
    derivative_factor = RationalInterval(
        _checked_div(odds.lower, lower_square, limits, "sigma derivative factor"),
        _checked_div(odds.upper, upper_square, limits, "sigma derivative factor"),
    )
    x_cube = _checked_pow(x, 3, limits, "sigma derivative scale")
    reflected_cube = _checked_pow(reflected_x, 3, limits, "sigma derivative scale")
    derivative_scale = _checked_add(
        _checked_div(Fraction(2), x_cube, limits, "sigma derivative scale"),
        _checked_div(Fraction(2), reflected_cube, limits, "sigma derivative scale"),
        limits,
        "sigma derivative scale",
    )
    derivative = RationalInterval(
        _checked_mul(
            derivative_factor.lower, derivative_scale, limits, "sigma derivative"
        ),
        _checked_mul(
            derivative_factor.upper, derivative_scale, limits, "sigma derivative"
        ),
    )
    return sigma, derivative


def _endpoint_derivative_cap() -> Fraction:
    """A rational cap valid on ``0 <= x <= 1/4`` and by reflection on the right."""

    p = _ENDPOINT_CAP_POWER
    return (
        Fraction(_factorial(p), 1)
        * Fraction(9, 8) ** p
        * Fraction(56, 27)
        * Fraction(1, 4) ** (2 * p - 3)
    )


_ENDPOINT_DERIVATIVE_CAP = _endpoint_derivative_cap()
ENDPOINT_DERIVATIVE_CAP = _ENDPOINT_DERIVATIVE_CAP


def _left_derivative_upper(
    lower: Fraction,
    upper: Fraction,
    *,
    limits: ResourceLimits,
    series_terms: int,
    max_reduction_steps: int,
    fallback_power: int,
) -> Fraction:
    """Upper-bound sigma' on ``0 <= lower < upper <= 1/2``."""

    if lower < 0 or upper > Fraction(1, 2) or lower >= upper:
        raise ValueError("left derivative interval must lie in (0, 1/2]")
    if upper <= Fraction(1, 4):
        return _ENDPOINT_DERIVATIVE_CAP
    cap = _ENDPOINT_DERIVATIVE_CAP if lower < Fraction(1, 4) else Fraction(0)
    central_lower = max(lower, Fraction(1, 4))
    if central_lower >= upper:
        return cap

    if upper == Fraction(1, 2):
        odds_upper = Fraction(1)
    else:
        upper_square = _checked_pow(upper, 2, limits, "derivative endpoint")
        reflected_upper = _checked_sub(
            Fraction(1), upper, limits, "derivative endpoint"
        )
        reflected_square = _checked_pow(
            reflected_upper, 2, limits, "derivative endpoint"
        )
        t_upper = _checked_sub(
            _checked_div(Fraction(1), upper_square, limits, "derivative endpoint"),
            _checked_div(
                Fraction(1), reflected_square, limits, "derivative endpoint"
            ),
            limits,
            "derivative endpoint",
        )
        odds_upper = _odds_from_exp(
            t_upper,
            limits,
            series_terms=series_terms,
            max_reduction_steps=max_reduction_steps,
            fallback_power=fallback_power,
        ).upper
    one_plus_odds = _checked_add(
        Fraction(1), odds_upper, limits, "derivative factor"
    )
    odds_factor_upper = _checked_div(
        odds_upper,
        _checked_pow(one_plus_odds, 2, limits, "derivative factor"),
        limits,
        "derivative factor",
    )
    central_cube = _checked_pow(central_lower, 3, limits, "derivative scale")
    reflected_upper = _checked_sub(
        Fraction(1), upper, limits, "derivative scale"
    )
    reflected_cube = _checked_pow(reflected_upper, 3, limits, "derivative scale")
    derivative_scale_upper = _checked_add(
        _checked_div(Fraction(2), central_cube, limits, "derivative scale"),
        _checked_div(Fraction(2), reflected_cube, limits, "derivative scale"),
        limits,
        "derivative scale",
    )
    return max(
        cap,
        _checked_mul(
            odds_factor_upper, derivative_scale_upper, limits, "derivative upper"
        ),
    )


def _derivative_upper_on_unit_interval(
    lower: Fraction,
    upper: Fraction,
    *,
    limits: ResourceLimits,
    series_terms: int,
    max_reduction_steps: int,
    fallback_power: int,
) -> Fraction:
    """Bound the derivative on the portion of an input interval in [0, 1]."""

    clipped_lower = max(lower, Fraction(0))
    clipped_upper = min(upper, Fraction(1))
    if clipped_lower >= clipped_upper:
        return Fraction(0)

    cap = Fraction(0)
    left_upper = min(clipped_upper, Fraction(1, 2))
    if clipped_lower < left_upper:
        cap = max(
            cap,
            _left_derivative_upper(
                clipped_lower,
                left_upper,
                limits=limits,
                series_terms=series_terms,
                max_reduction_steps=max_reduction_steps,
                fallback_power=fallback_power,
            ),
        )

    right_lower = max(clipped_lower, Fraction(1, 2))
    if right_lower < clipped_upper:
        # y = 1 - x maps this portion to the left half and preserves sigma'.
        reflected_lower = _checked_sub(
            Fraction(1), clipped_upper, limits, "derivative reflection"
        )
        reflected_upper = _checked_sub(
            Fraction(1), right_lower, limits, "derivative reflection"
        )
        cap = max(
            cap,
            _left_derivative_upper(
                reflected_lower,
                reflected_upper,
                limits=limits,
                series_terms=series_terms,
                max_reduction_steps=max_reduction_steps,
                fallback_power=fallback_power,
            ),
        )
    return cap


def evaluate_sigma_interval(
    lower: Fraction,
    upper: Fraction,
    *,
    series_terms: int = DEFAULT_SERIES_TERMS,
    max_reduction_steps: int = DEFAULT_MAX_REDUCTION_STEPS,
    fallback_power: int = DEFAULT_FALLBACK_POWER,
    resource_limits: ResourceLimits = DEFAULT_RESOURCE_LIMITS,
) -> SigmaIntervalResult:
    """Enclose sigma and its ordinary input-x derivative on a rational cell.

    The cell is closed and may lie outside [0, 1] or cross any of 0, 1/2,
    and 1.  Sigma is monotone, so its enclosure uses the two rational
    endpoint enclosures.  The derivative lower endpoint is zero for a
    nondegenerate cell (a safe full-cell bound); the upper endpoint is the
    union of the flat endpoint caps and central monotone bounds.  A
    degenerate cell uses the direct rational-point enclosure, retaining a
    nonzero derivative lower endpoint when one is exactly certified.
    """

    lower = _require_fraction(lower, "lower")
    upper = _require_fraction(upper, "upper")
    if not isinstance(resource_limits, ResourceLimits):
        raise TypeError("resource_limits must be ResourceLimits")
    _check_fraction_bits(lower, resource_limits, "lower input", input_value=True)
    _check_fraction_bits(upper, resource_limits, "upper input", input_value=True)
    if lower > upper:
        raise ValueError("lower must not exceed upper")
    series_terms = _positive_integer(series_terms, "series_terms")
    max_reduction_steps = _positive_integer(
        max_reduction_steps, "max_reduction_steps"
    )
    fallback_power = _positive_integer(fallback_power, "fallback_power")
    if series_terms > resource_limits.max_series_terms:
        raise ResourceLimitError("series_terms exceeds resource budget")
    if max_reduction_steps > resource_limits.max_reduction_steps:
        raise ResourceLimitError("max_reduction_steps exceeds resource budget")
    if fallback_power > resource_limits.max_fallback_power:
        raise ResourceLimitError("fallback_power exceeds resource budget")

    input_interval = RationalInterval(lower, upper)
    if lower == upper:
        sigma, derivative = _sigma_point(
            lower,
            limits=resource_limits,
            series_terms=series_terms,
            max_reduction_steps=max_reduction_steps,
            fallback_power=fallback_power,
        )
        return SigmaIntervalResult(input_interval, sigma, derivative)

    sigma_lower, _ = _sigma_point(
        lower,
        limits=resource_limits,
        series_terms=series_terms,
        max_reduction_steps=max_reduction_steps,
        fallback_power=fallback_power,
    )
    sigma_upper, _ = _sigma_point(
        upper,
        limits=resource_limits,
        series_terms=series_terms,
        max_reduction_steps=max_reduction_steps,
        fallback_power=fallback_power,
    )
    sigma = RationalInterval(sigma_lower.lower, sigma_upper.upper)
    derivative_upper = _derivative_upper_on_unit_interval(
        lower,
        upper,
        limits=resource_limits,
        series_terms=series_terms,
        max_reduction_steps=max_reduction_steps,
        fallback_power=fallback_power,
    )
    derivative = RationalInterval(Fraction(0), derivative_upper)
    return SigmaIntervalResult(input_interval, sigma, derivative)


def sigma_interval(
    lower: Fraction,
    upper: Fraction,
    **controls: object,
) -> SigmaIntervalResult:
    """Short alias for :func:`evaluate_sigma_interval`."""

    return evaluate_sigma_interval(lower, upper, **controls)


__all__ = [
    "COMPONENT",
    "DEFAULT_FALLBACK_POWER",
    "DEFAULT_MAX_REDUCTION_STEPS",
    "DEFAULT_RESOURCE_LIMITS",
    "DEFAULT_SERIES_TERMS",
    "ENDPOINT_DERIVATIVE_CAP",
    "MAX_FALLBACK_POWER",
    "MAX_INPUT_DENOMINATOR_BITS",
    "MAX_INPUT_NUMERATOR_BITS",
    "MAX_INTERMEDIATE_DENOMINATOR_BITS",
    "MAX_INTERMEDIATE_NUMERATOR_BITS",
    "MAX_REDUCTION_STEPS",
    "MAX_SERIES_TERMS",
    "RationalInterval",
    "ResourceLimitError",
    "ResourceLimits",
    "SCHEMA",
    "SOURCE_DEFINITION",
    "SigmaIntervalResult",
    "evaluate_sigma_interval",
    "rational_exp_negative",
    "sigma_interval",
]
