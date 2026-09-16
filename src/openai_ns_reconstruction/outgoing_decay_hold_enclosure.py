"""Certified rational enclosure of the pinned outgoing ``decayHold`` scalar.

For one actual :class:`~openai_ns_reconstruction.outgoing_tail.TailData`, let
``L = releaseLag(rampEnd)``, ``D = tailDebt`` and ``k = 1 - h``.  The pinned
construction defines

    decayHold = log(L / D) / k.

This module composes the independently certified actual-source release-lag and
tail-debt intervals.  It first proves ``L.lower > D.upper > 0`` and then uses
the monotonicity of ``log`` and exact rational logarithm enclosures.  Source
accuracy is budgeted from the requested final hold width: lag and debt receive
different absolute tolerances because their logarithmic sensitivities scale as
``1/L`` and ``1/D``.  The complete hold interval is checked directly and the
sources are refined further if the first budget is insufficient.

The result is a certified finite scalar enclosure only.  It does not certify
the later transition geometry, pressure integral, a global AxisSpace norm, or
the paper-exact reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .axis_amplitude_derivative_enclosure import rational_log_enclosure
from .outgoing_release_lag_enclosure import (
    ReleaseLagEnclosure,
    validated_release_lag_enclosure,
)
from .outgoing_sigma_enclosure import RationalInterval
from .outgoing_tail import TailData
from .outgoing_tail_debt_enclosure import (
    TailDebtEnclosure,
    validated_tail_debt_enclosure,
)


_DEFAULT_TOLERANCE = Fraction(1, 64)
_DEFAULT_MAX_CELLS = 16384
_DEFAULT_MAX_TERMS = 1024
_DEFAULT_MAX_SQUARINGS = 4096
_DEFAULT_MAX_REFINEMENTS = 8
_PILOT_LAG_TOLERANCE = Fraction(1, 64)
_PILOT_DEBT_TOLERANCE = Fraction(1, 10**8)
_LAG_HOLD_SHARE = Fraction(1, 8)
_DEBT_HOLD_SHARE = Fraction(3, 4)
_LOG_ENDPOINT_SHARE = Fraction(1, 16)


def _positive_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _positive_cap(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _require_actual_data(data: object) -> tuple[TailData, Fraction, Fraction]:
    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    h = Fraction.from_float(data.h)
    lam = Fraction.from_float(data.core.lam)
    if not Fraction(0) < 2 * h < lam < Fraction(1, 10):
        raise ValueError("TailData must satisfy 0 < 2*h < lambda < 1/10")
    return data, h, lam


def _validate_source_identity(
    *,
    h: Fraction,
    lam: Fraction,
    release_lag: ReleaseLagEnclosure,
    tail_debt: TailDebtEnclosure,
) -> None:
    if not isinstance(release_lag, ReleaseLagEnclosure):
        raise TypeError("release_lag must be a ReleaseLagEnclosure")
    if not isinstance(tail_debt, TailDebtEnclosure):
        raise TypeError("tail_debt must be a TailDebtEnclosure")
    if release_lag.h != h or release_lag.lam != lam:
        raise ValueError("release lag parameters do not match TailData")
    if tail_debt.h != h:
        raise ValueError("tail debt parameters do not match TailData")


def _initial_source_tolerances(
    *,
    h: Fraction,
    hold_tolerance: Fraction,
    release_lag: ReleaseLagEnclosure,
    tail_debt: TailDebtEnclosure,
) -> tuple[Fraction, Fraction, Fraction]:
    """Allocate source/log tolerances from the requested final hold width.

    The first-order logarithmic sensitivities are ``1 / ((1-h)L)`` and
    ``1 / ((1-h)D)``.  We therefore scale the two absolute source budgets by
    certified positive pilot lower bounds rather than giving them the same
    tolerance.  Debt receives more of the final hold-width budget because its
    absolute scale is much smaller; the final enclosure check remains the
    authority if this sensitivity allocation is not yet tight enough.

    The pilot enclosures are used only to obtain certified positive scale
    information.  They do not cap these budgets: otherwise a fixed pilot
    tolerance can dominate two different requested final hold widths and make
    the reported source budget insensitive to the caller's requested accuracy.
    """

    k = 1 - h
    lag_tolerance = (
        k * hold_tolerance * _LAG_HOLD_SHARE * release_lag.lag.lower
    )
    debt_tolerance = (
        k * hold_tolerance * _DEBT_HOLD_SHARE * tail_debt.debt.lower
    )
    log_tolerance = k * hold_tolerance * _LOG_ENDPOINT_SHARE
    if lag_tolerance <= 0 or debt_tolerance <= 0 or log_tolerance <= 0:
        raise ArithmeticError("decay-hold tolerance allocation lost positivity")
    return lag_tolerance, debt_tolerance, log_tolerance


@dataclass(frozen=True)
class DecayHoldEnclosure:
    """Exact rational bounds for ``OutgoingTail.decayHold``."""

    h: Fraction
    lam: Fraction
    release_lag_source: ReleaseLagEnclosure
    tail_debt_source: TailDebtEnclosure
    log_ratio: RationalInterval
    hold: RationalInterval
    requested_tolerance: Fraction
    lag_tolerance: Fraction
    debt_tolerance: Fraction
    log_tolerance: Fraction
    refinements: int

    def __post_init__(self) -> None:
        if not isinstance(self.h, Fraction) or not isinstance(self.lam, Fraction):
            raise TypeError("h and lam must be Fractions")
        if not Fraction(0) < 2 * self.h < self.lam < Fraction(1, 10):
            raise ValueError("h and lam must satisfy 0 < 2*h < lam < 1/10")
        _validate_source_identity(
            h=self.h,
            lam=self.lam,
            release_lag=self.release_lag_source,
            tail_debt=self.tail_debt_source,
        )
        if self.release_lag_source.lag.lower <= self.tail_debt_source.debt.upper:
            raise ValueError("release lag lower bound must exceed tail debt upper bound")
        if not isinstance(self.log_ratio, RationalInterval) or self.log_ratio.lower <= 0:
            raise ValueError("log_ratio must be a positive RationalInterval")
        if not isinstance(self.hold, RationalInterval) or self.hold.lower <= 0:
            raise ValueError("hold must be a positive RationalInterval")
        _positive_fraction(self.requested_tolerance, "requested_tolerance")
        _positive_fraction(self.lag_tolerance, "lag_tolerance")
        _positive_fraction(self.debt_tolerance, "debt_tolerance")
        _positive_fraction(self.log_tolerance, "log_tolerance")
        if isinstance(self.refinements, bool) or not isinstance(self.refinements, int):
            raise ValueError("refinements must be a nonnegative integer")
        if self.refinements < 0:
            raise ValueError("refinements must be a nonnegative integer")

    @property
    def decay_hold(self) -> RationalInterval:
        return self.hold

    @property
    def lag(self) -> RationalInterval:
        return self.release_lag_source.lag

    @property
    def debt(self) -> RationalInterval:
        return self.tail_debt_source.debt

    @property
    def lag_cells(self) -> int:
        return self.release_lag_source.cells

    @property
    def debt_cells(self) -> int:
        return self.tail_debt_source.cells

    @property
    def width(self) -> Fraction:
        return self.hold.width

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False


def _compose_decay_hold(
    *,
    h: Fraction,
    lam: Fraction,
    release_lag: ReleaseLagEnclosure,
    tail_debt: TailDebtEnclosure,
    hold_tolerance: Fraction,
    lag_tolerance: Fraction,
    debt_tolerance: Fraction,
    log_tolerance: Fraction,
    max_terms: int,
    refinements: int,
) -> DecayHoldEnclosure:
    _validate_source_identity(
        h=h,
        lam=lam,
        release_lag=release_lag,
        tail_debt=tail_debt,
    )
    if release_lag.lag.lower <= tail_debt.debt.upper:
        raise ArithmeticError("release lag lower bound must exceed tail debt upper bound")
    if tail_debt.debt.lower <= 0:
        raise ArithmeticError("tail debt lower bound must be positive")

    ratio_lower = release_lag.lag.lower / tail_debt.debt.upper
    ratio_upper = release_lag.lag.upper / tail_debt.debt.lower
    if not Fraction(1) < ratio_lower <= ratio_upper:
        raise ArithmeticError("decay-hold ratio interval must lie strictly above one")

    lower_log = rational_log_enclosure(
        ratio_lower,
        absolute_tolerance=log_tolerance,
        max_terms=max_terms,
    )
    upper_log = rational_log_enclosure(
        ratio_upper,
        absolute_tolerance=log_tolerance,
        max_terms=max_terms,
    )
    log_ratio = RationalInterval(lower_log.lower, upper_log.upper)
    k = 1 - h
    hold = RationalInterval(log_ratio.lower / k, log_ratio.upper / k)
    if hold.lower <= 0:
        raise ArithmeticError("decay hold enclosure lost positivity")

    return DecayHoldEnclosure(
        h=h,
        lam=lam,
        release_lag_source=release_lag,
        tail_debt_source=tail_debt,
        log_ratio=log_ratio,
        hold=hold,
        requested_tolerance=hold_tolerance,
        lag_tolerance=lag_tolerance,
        debt_tolerance=debt_tolerance,
        log_tolerance=log_tolerance,
        refinements=refinements,
    )


def validated_decay_hold_enclosure(
    data: TailData,
    *,
    absolute_tolerance: Fraction = _DEFAULT_TOLERANCE,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_squarings: int = _DEFAULT_MAX_SQUARINGS,
    max_refinements: int = _DEFAULT_MAX_REFINEMENTS,
    release_lag: ReleaseLagEnclosure | None = None,
    tail_debt: TailDebtEnclosure | None = None,
) -> DecayHoldEnclosure:
    """Return a certified rational interval for the actual ``decayHold``.

    Passing both source enclosures is supported for provenance-preserving reuse;
    they must carry the exact ``h``/``lambda`` values of ``data`` and already be
    narrow enough for the requested final tolerance.  The normal path computes
    both sources from the same ``TailData`` and adaptively refines them.
    """

    data, h, lam = _require_actual_data(data)
    tolerance = _positive_fraction(absolute_tolerance, "absolute_tolerance")
    max_cells = _positive_cap(max_cells, "max_cells")
    max_terms = _positive_cap(max_terms, "max_terms")
    max_squarings = _positive_cap(max_squarings, "max_squarings")
    max_refinements = _positive_cap(max_refinements, "max_refinements")

    if (release_lag is None) != (tail_debt is None):
        raise ValueError("release_lag and tail_debt must be supplied together")

    if release_lag is not None and tail_debt is not None:
        _validate_source_identity(
            h=h,
            lam=lam,
            release_lag=release_lag,
            tail_debt=tail_debt,
        )
        log_tolerance = (1 - h) * tolerance * _LOG_ENDPOINT_SHARE
        result = _compose_decay_hold(
            h=h,
            lam=lam,
            release_lag=release_lag,
            tail_debt=tail_debt,
            hold_tolerance=tolerance,
            lag_tolerance=max(release_lag.width, Fraction(1, 10**120)),
            debt_tolerance=max(tail_debt.width, Fraction(1, 10**120)),
            log_tolerance=log_tolerance,
            max_terms=max_terms,
            refinements=0,
        )
        if result.width > tolerance:
            raise ArithmeticError("supplied decay-hold source enclosures are too wide")
        return result

    try:
        pilot_lag = validated_release_lag_enclosure(
            data,
            absolute_tolerance=_PILOT_LAG_TOLERANCE,
            max_cells=max_cells,
            max_terms=max_terms,
            max_squarings=max_squarings,
        )
        pilot_debt = validated_tail_debt_enclosure(
            data,
            absolute_tolerance=_PILOT_DEBT_TOLERANCE,
            max_cells=max_cells,
            max_terms=max_terms,
            max_squarings=max_squarings,
        )
    except ArithmeticError as exc:
        raise ArithmeticError(f"decay hold source enclosure failed: {exc}") from exc

    lag_tolerance, debt_tolerance, log_tolerance = _initial_source_tolerances(
        h=h,
        hold_tolerance=tolerance,
        release_lag=pilot_lag,
        tail_debt=pilot_debt,
    )

    for refinement in range(max_refinements):
        try:
            release_lag = (
                pilot_lag
                if pilot_lag.width <= lag_tolerance
                else validated_release_lag_enclosure(
                    data,
                    absolute_tolerance=lag_tolerance,
                    max_cells=max_cells,
                    max_terms=max_terms,
                    max_squarings=max_squarings,
                )
            )
            tail_debt = (
                pilot_debt
                if pilot_debt.width <= debt_tolerance
                else validated_tail_debt_enclosure(
                    data,
                    absolute_tolerance=debt_tolerance,
                    max_cells=max_cells,
                    max_terms=max_terms,
                    max_squarings=max_squarings,
                )
            )
        except ArithmeticError as exc:
            raise ArithmeticError(f"decay hold source enclosure failed: {exc}") from exc

        result = _compose_decay_hold(
            h=h,
            lam=lam,
            release_lag=release_lag,
            tail_debt=tail_debt,
            hold_tolerance=tolerance,
            lag_tolerance=lag_tolerance,
            debt_tolerance=debt_tolerance,
            log_tolerance=log_tolerance,
            max_terms=max_terms,
            refinements=refinement,
        )
        if result.width <= tolerance:
            return result

        lag_tolerance /= 2
        debt_tolerance /= 2
        log_tolerance /= 2

    raise ArithmeticError(
        "decay hold enclosure did not reach requested tolerance before max_refinements"
    )


decay_hold_enclosure = validated_decay_hold_enclosure
validated_decay_hold = validated_decay_hold_enclosure


__all__ = [
    "DecayHoldEnclosure",
    "decay_hold_enclosure",
    "validated_decay_hold",
    "validated_decay_hold_enclosure",
]
