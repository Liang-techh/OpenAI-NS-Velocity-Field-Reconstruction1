"""Cancellation-free continuous waiting-length matching.

For the zero-``Z`` angular moment, let ``X_t`` be the normalized incoming
moment at ``y_t`` (the start of the waiting/terminal collar) and put
``a = 1 - delta/2``.  With the heat factor replaced by ``H_delta = 1``,

``K0(t) = 1 - epsilon*k(t)``,

where

``k(t) = 1 - sigma(t) + sigma(t)*f((3-t)/2)``,
``f(s) = exp(-1/s**2)`` for ``s > 0``.

The finite collar quantity needed by (7.9)--(7.10) is

``J = integral_0^3 exp(a*t)*k(t) dt``.

Keeping the positive terms separate gives the exact matching equation

``0 = (X_t - 1/a) - exp(a*tau)*epsilon/(1-epsilon)*(1/a + J)``.

Therefore

``tau = log((X_t - 1/a)*(1-epsilon) /
           (epsilon*(1/a + J))) / a``.

The adapter evaluates ``X_t`` from an existing
``ContinuousAngularMoments`` engine and uses its MP nodes when available.
It does not rebuild a field, install a schedule, or subtract the whole
terminal bracket ``K0`` from one.  The returned values are raw mpmath
numbers so an exponentially small ``epsilon`` remains observable.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_continuous_axial_pulse import (  # noqa: E402
    ContinuousAxialPulse,
)


SOURCE = "https://arxiv.org/html/2609.35406v2"
SOURCE_VERSION = "2609.35406v2"


def _mp(value: Any) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(str(value))


def _n(value: Any, digits: int) -> str:
    return mp.nstr(_mp(value), int(digits))


def _flat_edge(value: mp.mpf) -> mp.mpf:
    value = _mp(value)
    if value <= 0:
        return mp.mpf(0)
    return mp.exp(-1 / (value * value))


def _k_factor(t: mp.mpf) -> mp.mpf:
    """Return ``k(t)=1-sigma(t)+sigma(t)f((3-t)/2)``."""

    t = _mp(t)
    if t <= 0:
        return mp.mpf(1)
    if t >= 3:
        return mp.mpf(0)
    sigma, _ = ContinuousAxialPulse.sigma_pair(t)
    return 1 - sigma + sigma * _flat_edge((3 - t) / 2)


def _k0_factor(t: mp.mpf, epsilon: mp.mpf) -> mp.mpf:
    """Return the heat-collar bracket ``K0=1-epsilon*k``."""

    return 1 - _mp(epsilon) * _k_factor(_mp(t))


class ContinuousWaitingMatch:
    """Solve the continuous zero-``Z`` pre-heat waiting equation.

    ``engine`` is an existing continuous angular moment engine.  The adapter
    only calls its private normalized baseline ``_baseline`` at ``y_t``;
    this is intentional because that baseline is the source of the actual
    incoming moment used by the surrounding continuous construction.
    """

    def __init__(
        self,
        engine: Any,
        schedule: Any | None = None,
        *,
        quadrature_order: int | None = None,
        precision: int | None = None,
    ) -> None:
        if not hasattr(engine, "_baseline"):
            raise TypeError("engine must expose ContinuousAngularMoments._baseline")
        if schedule is None:
            schedule = getattr(engine, "schedule", None)
        if schedule is None:
            raise TypeError("schedule must be supplied or available on engine")
        order = getattr(engine, "order", None) if quadrature_order is None else quadrature_order
        if order is None:
            order = 192
        if int(order) != order or int(order) < 16:
            raise ValueError("quadrature_order must be an integer at least 16")
        engine_precision = int(getattr(engine, "precision", 80))
        self.precision = max(32, engine_precision if precision is None else int(precision))
        self.engine = engine
        self.schedule = schedule
        self.quadrature_order = int(order)
        with mp.workdps(self.precision):
            self.delta = _mp(schedule.delta)
            self.epsilon = _mp(schedule.epsilon)
            self.a = 1 - self.delta / 2
            self.y_t = _mp(schedule.y_t)
        self._nodes_cache: tuple[tuple[mp.mpf, mp.mpf], ...] | None = None

    def nodes(self) -> tuple[tuple[mp.mpf, mp.mpf], ...]:
        """Return [0,1] nodes shared with the angular engine when present."""

        if self._nodes_cache is not None:
            return self._nodes_cache
        supplied = getattr(self.engine, "nodes", None)
        if supplied is not None and len(supplied) == self.quadrature_order:
            self._nodes_cache = tuple((_mp(node), _mp(weight)) for node, weight in supplied)
            return self._nodes_cache
        with mp.workdps(self.precision):
            nodes, weights = mp.gauss_quadrature(self.quadrature_order, "legendre")
            self._nodes_cache = tuple(
                ((_mp(node) + 1) / 2, _mp(weight) / 2)
                for node, weight in zip(nodes, weights)
            )
        return self._nodes_cache

    @property
    def nodes_shared_with_engine(self) -> bool:
        supplied = getattr(self.engine, "nodes", None)
        return supplied is not None and len(supplied) == self.quadrature_order

    def incoming_X_t(self) -> mp.mpf:
        """Read the actual normalized incoming moment at ``y_t``."""

        with mp.workdps(self.precision):
            # Keep the Decimal checkpoint spelling intact.  The engine's
            # baseline converts it through its own exact-radius helper.
            row = self.engine._baseline(str(self.schedule.y_t), "0")
            if not row:
                raise ArithmeticError("continuous angular baseline returned no X_t")
            return _mp(row[0])

    def collar_k_integral(self) -> mp.mpf:
        """Return ``J = integral_0^3 exp(a*t) k(t) dt``."""

        with mp.workdps(self.precision):
            three = mp.mpf(3)
            total = mp.fsum(
                weight
                * mp.exp(self.a * three * node)
                * _k_factor(three * node)
                for node, weight in self.nodes()
            )
            return three * total

    def collar_k0_integral(self) -> mp.mpf:
        """Return ``integral_0^3 exp(a*t) K0(t) dt`` independently."""

        with mp.workdps(self.precision):
            three = mp.mpf(3)
            total = mp.fsum(
                weight
                * mp.exp(self.a * three * node)
                * _k0_factor(three * node, self.epsilon)
                for node, weight in self.nodes()
            )
            return three * total

    def solve(self, *, X_t: Any | None = None, J: Any | None = None) -> dict[str, Any]:
        """Solve the nonnegative waiting root and retain all positive atoms."""

        with mp.workdps(self.precision):
            incoming = self.incoming_X_t() if X_t is None else _mp(X_t)
            collar = self.collar_k_integral() if J is None else _mp(J)
            one_over_a = 1 / self.a
            incoming_excess = incoming - one_over_a
            epsilon_fraction = self.epsilon / (1 - self.epsilon)
            collar_positive_term = one_over_a + collar
            growth_coefficient = epsilon_fraction * collar_positive_term
            if incoming_excess <= 0:
                raise ValueError(
                    "no positive continuous waiting root: X_t must exceed 1/a"
                )
            if growth_coefficient <= 0:
                raise ValueError("waiting growth coefficient must be positive")
            exp_a_tau = incoming_excess / growth_coefficient
            if exp_a_tau < 1:
                raise ValueError(
                    "no nonnegative continuous waiting root: exp(a*tau) < 1"
                )
            tau = mp.log(exp_a_tau) / self.a
            replay_positive = exp_a_tau * growth_coefficient
            residual = incoming_excess - replay_positive
            scale = max(abs(incoming_excess), abs(replay_positive), mp.mpf(1))
            normalized_residual = residual / scale
            # This is the displayed (7.10) left side, retained separately
            # from all pressure and angular-bump corrections.
            matching_left_side = (
                incoming - one_over_a
                - mp.exp(self.a * tau)
                * self.epsilon
                / (1 - self.epsilon)
                * (one_over_a + collar)
            )
            result = {
                "source": SOURCE,
                "source_version": SOURCE_VERSION,
                "stage": "continuous_zero_Z_H1_waiting_match",
                "quadrature_order": self.quadrature_order,
                "precision": self.precision,
                "nodes_shared_with_engine": self.nodes_shared_with_engine,
                "X_t": incoming,
                "J": collar,
                "a": self.a,
                "delta": self.delta,
                "epsilon": self.epsilon,
                "waiting_length": tau,
                "tau": tau,
                "exp_a_tau": exp_a_tau,
                "one_over_a": one_over_a,
                "incoming_excess_X_t_minus_one_over_a": incoming_excess,
                "epsilon_over_one_minus_epsilon": epsilon_fraction,
                "one_over_a_plus_J": collar_positive_term,
                "positive_growth_coefficient": growth_coefficient,
                "replayed_positive_term": replay_positive,
                "matching_left_side": matching_left_side,
                "normalized_residual": normalized_residual,
                "normalized_match_residual": matching_left_side / scale,
                "collar_k0_integral": self.collar_k0_integral(),
                "formula": (
                    "0=(X_t-1/a)-exp(a*tau)*epsilon/(1-epsilon)*(1/a+J)"
                ),
                "k_formula": "1-sigma(t)+sigma(t)*exp(-1/((3-t)/2)^2)",
                "K0_formula": "1-epsilon*k(t)",
                "H_replaced_by_one": True,
                "whole_K0_subtraction_used": False,
                "quadrature_error_enclosed": False,
                "actual_heat_moment_closed": False,
                "other_four_moments_closed": False,
                "field_installed": False,
            }
            return result

    def receipt(self, *, X_t: Any | None = None, J: Any | None = None) -> dict[str, Any]:
        """JSON-friendly receipt retaining the positive factors explicitly."""

        row = self.solve(X_t=X_t, J=J)
        digits = self.precision
        keys = (
            "X_t",
            "J",
            "a",
            "delta",
            "epsilon",
            "waiting_length",
            "tau",
            "exp_a_tau",
            "one_over_a",
            "incoming_excess_X_t_minus_one_over_a",
            "epsilon_over_one_minus_epsilon",
            "one_over_a_plus_J",
            "positive_growth_coefficient",
            "replayed_positive_term",
            "matching_left_side",
            "normalized_residual",
            "normalized_match_residual",
            "collar_k0_integral",
        )
        result = {
            key: _n(row[key], digits)
            for key in keys
        }
        result.update(
            {
                "source": row["source"],
                "source_version": row["source_version"],
                "stage": row["stage"],
                "quadrature_order": row["quadrature_order"],
                "precision": row["precision"],
                "nodes_shared_with_engine": row["nodes_shared_with_engine"],
                "formula": row["formula"],
                "k_formula": row["k_formula"],
                "K0_formula": row["K0_formula"],
                "H_replaced_by_one": row["H_replaced_by_one"],
                "whole_K0_subtraction_used": row["whole_K0_subtraction_used"],
                "quadrature_error_enclosed": row["quadrature_error_enclosed"],
                "actual_heat_moment_closed": row["actual_heat_moment_closed"],
                "other_four_moments_closed": row["other_four_moments_closed"],
                "field_installed": row["field_installed"],
            }
        )
        return result


def solve_continuous_waiting_match(
    engine: Any,
    schedule: Any | None = None,
    *,
    quadrature_order: int | None = None,
    precision: int | None = None,
    X_t: Any | None = None,
    J: Any | None = None,
) -> dict[str, Any]:
    """Functional wrapper around :class:`ContinuousWaitingMatch`."""

    solver = ContinuousWaitingMatch(
        engine,
        schedule,
        quadrature_order=quadrature_order,
        precision=precision,
    )
    return solver.solve(X_t=X_t, J=J)


def coherent_waiting_length(
    engine: Any,
    schedule: Any | None = None,
    *,
    quadrature_order: int | None = None,
    precision: int | None = None,
) -> dict[str, Any]:
    """Parent-integration spelling for the coherent zero-``Z`` root."""

    return solve_continuous_waiting_match(
        engine,
        schedule,
        quadrature_order=quadrature_order,
        precision=precision,
    )


# Descriptive aliases for integration code using the paper's terminology.
ContinuousWaitingLengthMatch = ContinuousWaitingMatch
match_waiting_length = solve_continuous_waiting_match


def independent_k0_fixture(*, order: int = 192, precision: int = 100) -> dict[str, Any]:
    """Moderate-epsilon fixture checking the sign and collar factors.

    The fixture does not construct a schedule or field.  It compares the
    direct ``K0=1-epsilon*k`` integral against
    ``integral(exp(a*t)) - epsilon*J`` using independent adaptive quadrature.
    """

    if int(order) != order or int(order) < 16:
        raise ValueError("order must be an integer at least 16")
    with mp.workdps(int(precision)):
        a = mp.mpf("0.75")
        epsilon = mp.mpf("0.2")
        nodes, weights = mp.gauss_quadrature(int(order), "legendre")
        half_nodes = [((node + 1) / 2, weight / 2) for node, weight in zip(nodes, weights)]
        J_gauss = 3 * mp.fsum(
            weight * mp.exp(a * 3 * node) * _k_factor(3 * node)
            for node, weight in half_nodes
        )
        K0_gauss = 3 * mp.fsum(
            weight * mp.exp(a * 3 * node) * _k0_factor(3 * node, epsilon)
            for node, weight in half_nodes
        )
        whole_exp = (mp.exp(3 * a) - 1) / a
        K0_decomposed = whole_exp - epsilon * J_gauss
        K0_adaptive = mp.quad(
            lambda t: mp.exp(a * t) * _k0_factor(t, epsilon), [0, 1, 2, 3]
        )
        decomposition_residual = K0_decomposed - K0_adaptive
        return {
            "stage": "independent_moderate_epsilon_K0_fixture",
            "order": int(order),
            "precision": int(precision),
            "a": a,
            "epsilon": epsilon,
            "J": J_gauss,
            "K0_direct_gauss": K0_gauss,
            "K0_decomposed": K0_decomposed,
            "K0_adaptive": K0_adaptive,
            "decomposition_residual": decomposition_residual,
            "positive_J": bool(J_gauss > 0),
            "positive_K0": bool(K0_gauss > 0),
            # The independent adaptive replay sees the flat sigma edge as a
            # numerically narrow transition.  This tolerance checks the
            # sign/factor decomposition at the fixture's finite quadrature
            # order; it is deliberately not a quadrature enclosure claim.
            "sign_and_factor_check": bool(abs(decomposition_residual) < mp.mpf("1e-14")),
            "whole_K0_subtraction_used": False,
        }


__all__ = [
    "ContinuousWaitingMatch",
    "ContinuousWaitingLengthMatch",
    "solve_continuous_waiting_match",
    "coherent_waiting_length",
    "match_waiting_length",
    "independent_k0_fixture",
]
