"""Directed rational enclosures for the independent gap factor.

This module certifies only the scalar function

    ``g(D) = 1 - exp(-D)``

for a non-negative exact ``Fraction`` interval ``[D_lo, D_hi]``.  It is
deliberately independent of the F source function, endpoint identities,
transport formulae, and final error estimates.  In particular, a returned
endpoint is an enclosure bound, never a midpoint or a numerical estimate
presented as a value.

Only integer and :class:`fractions.Fraction` arithmetic is used.  There is no
floating-point conversion and no call to a floating-point exponential routine
(or to any legacy reconstruction module).

The enclosure has two exact paths.

* On a reduced argument ``t`` in ``[0, 1]``, the alternating Taylor partial
  sums

      S_n(t) = sum((-t)**j / j!, j=0..n)

  satisfy ``S_(2r+1) <= exp(-t) <= S_(2r)``.  The terms decrease because
  ``0 <= t <= 1`` and ``t/(j+1) <= 1``.  Thus every bound is rational and the
  alternating-series remainder has the sign of its first omitted term.  For
  ``x > 1`` we may set ``t = x / 2**k`` with the least ``k`` for which
  ``t <= 1``.  Monotonicity of multiplication on non-negative numbers gives
  ``S_lower(t)**(2**k) <= exp(-x) <= S_upper(t)**(2**k)``.

* For a large argument, let ``m = floor(x)``.  The defining positive series
  ``e = sum(1/j!, j>=0)`` has ``e > 1 + 1 = 2`` (the remaining terms are
  positive).  Consequently, for every integer ``m >= 1``,

      exp(-x) <= exp(-m) = e**(-m) < 2**(-m).

  We use ``2**(-k)`` with ``k = min(m, max_bits, tail_exponent_cap)``.
  This is a deliberately capped, strict rational upper bound: its denominator
  is never larger than the caller's explicit bit budget.  When the budget is
  zero the safe result is the full interval ``[0, 1]``; no large integer is
  allocated.

The Taylor route also has a fixed conservative ``65_536``-bit intermediate
work budget.  Its estimate uses the input numerator/denominator bit lengths,
the term count and its binary length, the dyadic reduction exponent, and the
final powering scale.  If that estimate is exceeded, the code takes the
strict dyadic tail route before constructing the reduced series.  Thus
``max_bits`` controls only the denominator budget of that tail route; it is
not a promise of Taylor accuracy, and the explicit ``max_terms`` knob is
still subject to the work guard.

Finally, ``g`` is monotone on the non-negative half-line: if ``a <= b``, then
``exp(-b) = exp(-a)*exp(-(b-a)) <= exp(-a)`` because the latter factor lies in
``(0, 1]``.  Thus, if ``[E_lo(a), E_hi(a)]`` and
``[E_lo(b), E_hi(b)]`` are the endpoint exponential enclosures, the returned
cross-endpoint pair is
``1 - E_hi(D_lo) <= g(D) <= 1 - E_lo(D_hi)`` for every
``D_lo <= D <= D_hi``.  This is the complete directed-interval composition;
no midpoint or unproved representative is introduced.

The function ``canonical_native_decay_enclosure`` is a separate conditional
static theorem for the historical canonical relation
``mu = (1/1000)*exp(-4*(exp(40)+11))`` and ``Delta_y = 2/mu - 4``.  It does
not evaluate that expression.  If and only if the caller has independently
proved that relation for the object in hand, then ``0 < mu < 1/1000`` gives
``Delta_y > 1996`` and hence the returned ``[1 - 2**(-bits), 1]`` encloses
``g(Delta_y)`` for ``0 <= bits <= 1996``.  This conditional statement does
not establish F ownership, source identity, or an endpoint identity.
"""

from fractions import Fraction


_DEFAULT_MAX_TERMS = 32
_DEFAULT_MAX_BITS = 128
_DEFAULT_TAIL_EXPONENT_CAP = 128
_HARD_MAX_TERMS = 4096
_HARD_MAX_BITS = 4096
_CANONICAL_MAX_BITS = 1996
_SERIES_WORK_BIT_LIMIT = 65_536

# Range reduction is intentionally kept small.  It bounds the number of exact
# squarings and prevents an apparently harmless caller input from asking for
# an enormous rational intermediate.  Larger arguments use the independently
# proved dyadic tail instead.
_SERIES_ARGUMENT_LIMIT = Fraction(32, 1)


def _require_fraction(value: object, name: str) -> Fraction:
    """Require an exact Fraction endpoint without implicit coercion."""

    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value < 0:
        raise ValueError(f"{name} must be non-negative")
    return value


def _require_budget(value: object, name: str, *, allow_zero: bool) -> int:
    """Validate an integer resource budget before any arithmetic uses it."""

    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(f"{name} must be an int")
    minimum = 0 if allow_zero else 1
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _validate_budgets(
    max_terms: object,
    max_bits: object,
    tail_exponent_cap: object | None,
) -> tuple[int, int, int]:
    """Validate and cap resource knobs without allocating caller-sized powers."""

    terms = _require_budget(max_terms, "max_terms", allow_zero=False)
    bits = _require_budget(max_bits, "max_bits", allow_zero=True)
    if terms > _HARD_MAX_TERMS:
        raise ValueError(f"max_terms must be <= {_HARD_MAX_TERMS}")
    if bits > _HARD_MAX_BITS:
        raise ValueError(f"max_bits must be <= {_HARD_MAX_BITS}")

    if tail_exponent_cap is None:
        cap = bits
    else:
        cap = _require_budget(
            tail_exponent_cap, "tail_exponent_cap", allow_zero=True
        )
        if cap > _HARD_MAX_BITS:
            raise ValueError(
                f"tail_exponent_cap must be <= {_HARD_MAX_BITS}"
            )
        cap = min(cap, bits)
    return terms, bits, cap


def _unit_taylor_interval(
    t: Fraction, max_terms: int
) -> tuple[Fraction, Fraction]:
    """Return rational lower/upper bounds for exp(-t), ``0 <= t <= 1``."""

    if not Fraction(0, 1) <= t <= Fraction(1, 1):
        raise ValueError("Taylor argument must lie in [0, 1]")
    if t == 0:
        return Fraction(1, 1), Fraction(1, 1)

    # A pair of adjacent partial sums is needed for a directed enclosure.  A
    # one-term budget cannot provide both signs, so fail closed.
    if max_terms < 2:
        return Fraction(0, 1), Fraction(1, 1)

    highest_index = max_terms - 1
    if highest_index % 2 == 0:
        upper_index = highest_index
        lower_index = highest_index - 1
    else:
        lower_index = highest_index
        upper_index = highest_index - 1

    term = Fraction(1, 1)
    partial = Fraction(1, 1)  # S_0
    lower: Fraction | None = None
    upper: Fraction | None = partial if upper_index == 0 else None

    for index in range(1, highest_index + 1):
        # term is (-t)**(index-1)/(index-1)! before this update.
        term = -term * t / index
        partial += term
        if index == lower_index:
            lower = partial
        if index == upper_index:
            upper = partial

    # For t in [0, 1], S_(odd) is non-negative; the min/max guards are
    # harmless exact clipping and keep the public invariant explicit.
    if lower is None or upper is None:
        return Fraction(0, 1), Fraction(1, 1)
    lower = max(Fraction(0, 1), lower)
    upper = min(Fraction(1, 1), upper)
    if lower > upper:
        # This is unreachable under the alternating-series proof, but a broad
        # interval is safer than exposing an invalid enclosure if the helper
        # is changed in the future.
        return Fraction(0, 1), Fraction(1, 1)
    return lower, upper


def _pow_fraction(base: Fraction, exponent: int) -> Fraction:
    """Exact non-negative integer power by bounded exponentiation by squaring."""

    result = Fraction(1, 1)
    factor = base
    remaining = exponent
    while remaining:
        if remaining & 1:
            result *= factor
        remaining >>= 1
        if remaining:
            factor *= factor
    return result


def _ceil_power_of_two_scale(value: Fraction) -> tuple[int, int]:
    """Return ``(k, 2**k)`` with ``value <= 2**k`` for ``value > 1``."""

    scale = 1
    exponent = 0
    while value > scale:
        scale <<= 1
        exponent += 1
    return exponent, scale


def _tail_exp_interval(
    value: Fraction, max_bits: int, tail_exponent_cap: int
) -> tuple[Fraction, Fraction]:
    """Return ``0 <= exp(-value) <= 2**(-k)`` with a capped denominator."""

    if value == 0:
        return Fraction(1, 1), Fraction(1, 1)
    capped_exponent = min(max_bits, tail_exponent_cap)
    if capped_exponent <= 0:
        return Fraction(0, 1), Fraction(1, 1)

    # Compare against the bounded cap before taking a quotient.  For a huge
    # numerator this short-circuits without materializing its floor; the
    # multiplication is by at most the caller's hard-capped budget.
    if value.numerator >= capped_exponent * value.denominator:
        exponent = capped_exponent
    else:
        # Here value < capped_exponent, so the quotient is itself bounded.
        integer_part = value.numerator // value.denominator
        exponent = integer_part
    if exponent <= 0:
        return Fraction(0, 1), Fraction(1, 1)
    # exponent is budget-capped before the shift, so this allocation is
    # bounded by the explicit caller budget (and the hard module cap).
    return Fraction(0, 1), Fraction(1, 1 << exponent)


def _series_work_within_budget(
    value: Fraction, max_terms: int, reduction_exponent: int, scale: int
) -> bool:
    """Conservatively bound Taylor/powering intermediate numerator bits.

    The input bit length covers the numerator and denominator of ``value``.
    Multiplication by ``max_terms`` covers powers and common denominators; the
    ``max_terms.bit_length()`` term covers factorial/index growth; the
    reduction exponent covers dividing by the dyadic scale; and the final
    ``scale`` factor covers exact powering.  This is intentionally an
    overestimate, used only to fail closed before expensive Fraction work.
    """

    input_bits = max(
        abs(value.numerator).bit_length(),
        value.denominator.bit_length(),
        1,
    )
    if input_bits > _SERIES_WORK_BIT_LIMIT:
        return False
    per_term_bits = (
        input_bits + max_terms.bit_length() + reduction_exponent + 8
    )
    estimated_bits = max_terms * per_term_bits * scale
    return estimated_bits <= _SERIES_WORK_BIT_LIMIT


def _exp_neg_interval(
    value: Fraction,
    *,
    max_terms: int,
    max_bits: int,
    tail_exponent_cap: int,
) -> tuple[Fraction, Fraction]:
    """Enclose ``exp(-value)`` using only exact rational operations."""

    if value == 0:
        return Fraction(1, 1), Fraction(1, 1)

    if value <= _SERIES_ARGUMENT_LIMIT and max_terms >= 2:
        reduction_exponent, scale = _ceil_power_of_two_scale(value)
        # ``scale`` is at most 32 on this path.  The work estimate is checked
        # before constructing the reduced Fraction or any Taylor term.
        if _series_work_within_budget(
            value, max_terms, reduction_exponent, scale
        ):
            reduced = value / scale
            reduced_lower, reduced_upper = _unit_taylor_interval(
                reduced, max_terms
            )
            power = scale  # 2**reduction_exponent
            return (
                _pow_fraction(reduced_lower, power),
                _pow_fraction(reduced_upper, power),
            )

    return _tail_exp_interval(value, max_bits, tail_exponent_cap)


def _gap_factor_at(
    value: Fraction,
    *,
    max_terms: int,
    max_bits: int,
    tail_exponent_cap: int,
) -> tuple[Fraction, Fraction]:
    """Enclose ``g(value)`` at one non-negative exact endpoint."""

    exp_lower, exp_upper = _exp_neg_interval(
        value,
        max_terms=max_terms,
        max_bits=max_bits,
        tail_exponent_cap=tail_exponent_cap,
    )
    # g = 1 - exp(-D) reverses the enclosure endpoints.  No representative
    # value is formed.
    gap_lower = Fraction(1, 1) - exp_upper
    gap_upper = Fraction(1, 1) - exp_lower
    gap_lower = max(Fraction(0, 1), gap_lower)
    gap_upper = min(Fraction(1, 1), gap_upper)
    if gap_lower > gap_upper:
        return Fraction(0, 1), Fraction(1, 1)
    return gap_lower, gap_upper


def gap_decay_interval(
    D_lo: Fraction,
    D_hi: Fraction,
    *,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_bits: int = _DEFAULT_MAX_BITS,
    tail_exponent_cap: int | None = _DEFAULT_TAIL_EXPONENT_CAP,
) -> tuple[Fraction, Fraction]:
    """Return a directed closed Fraction interval for ``1 - exp(-D)``.

    Parameters
    ----------
    D_lo, D_hi:
        Exact, non-negative ``Fraction`` endpoints with ``D_lo <= D_hi``.
    max_terms:
        Number of alternating Taylor terms available after range reduction.
        A value of one deliberately fails closed to a wide interval; larger
        values are bounded by a hard limit to prevent an accidental runaway.
    max_bits:
        Maximum dyadic exponent that this call may materialize in the large-D
        tail path.  It does not promise a Taylor precision: the Taylor path
        is instead subject to a fixed conservative intermediate-bit budget,
        and an over-budget call falls back to the strict tail interval.
        ``max_bits=0`` therefore gives ``[0, 1]`` whenever the tail path is
        needed, but does not disable an affordable Taylor enclosure.
    tail_exponent_cap:
        Optional additional cap on the dyadic exponent.  The effective
        exponent is ``min(floor(D), max_bits, tail_exponent_cap)``.  The
        default cap of 128 means even very large D only materializes a
        129-bit denominator.

    Returns
    -------
    tuple[Fraction, Fraction]
        ``(g_lo, g_hi)`` such that every ``D`` in ``[D_lo, D_hi]`` satisfies
        ``g_lo <= 1 - exp(-D) <= g_hi``.  The bounds are closed, directed,
        and independent; neither is a midpoint or a claimed point value.

    Raises
    ------
    TypeError, ValueError
        If endpoints are not exact non-negative Fractions, are reversed, or a
        resource budget is invalid.
    """

    lower_endpoint = _require_fraction(D_lo, "D_lo")
    upper_endpoint = _require_fraction(D_hi, "D_hi")
    if lower_endpoint > upper_endpoint:
        raise ValueError("D_lo must be <= D_hi")
    terms, bits, cap = _validate_budgets(
        max_terms, max_bits, tail_exponent_cap
    )

    lower_gap, _ = _gap_factor_at(
        lower_endpoint,
        max_terms=terms,
        max_bits=bits,
        tail_exponent_cap=cap,
    )
    _, upper_gap = _gap_factor_at(
        upper_endpoint,
        max_terms=terms,
        max_bits=bits,
        tail_exponent_cap=cap,
    )
    # Monotonicity of 1-exp(-D) proves this cross-endpoint composition.  The
    # defensive check is exact and returns a valid broad interval on any
    # future helper regression.
    if lower_gap > upper_gap:
        return Fraction(0, 1), Fraction(1, 1)
    return lower_gap, upper_gap


def canonical_native_decay_enclosure(
    bits: int = _DEFAULT_TAIL_EXPONENT_CAP,
) -> tuple[Fraction, Fraction]:
    """Return the conditional canonical ``Delta_y`` gap enclosure.

    This is a static theorem, not an evaluation.  Under the independently
    supplied premise

        ``mu = (1/1000) * exp(-4 * (exp(40) + 11))`` and
        ``Delta_y = 2/mu - 4``,

    positivity of the exponential implies ``0 < mu < 1/1000``.  Therefore
    ``Delta_y > 1996``.  The proof ``e > 2`` in this module then gives
    ``exp(-Delta_y) < 2**(-1996) <= 2**(-bits)`` for every
    ``0 <= bits <= 1996``.  Consequently the returned interval encloses
    ``1 - exp(-Delta_y)``.  It intentionally does not assert that any F field
    owns this ``mu``, nor that an endpoint identity or final error follows.

    The only materialized denominator is ``2**bits`` with ``bits <= 1996``.
    """

    checked_bits = _require_budget(bits, "bits", allow_zero=True)
    if checked_bits > _CANONICAL_MAX_BITS:
        raise ValueError(f"bits must be <= {_CANONICAL_MAX_BITS}")
    return Fraction(1, 1) - Fraction(1, 1 << checked_bits), Fraction(1, 1)


# Descriptive aliases keep the tuple API discoverable without creating another
# implementation or changing the proof obligations.
directed_gap_decay_interval = gap_decay_interval
gap_factor_interval = gap_decay_interval


__all__ = [
    "canonical_native_decay_enclosure",
    "directed_gap_decay_interval",
    "gap_decay_interval",
    "gap_factor_interval",
]
