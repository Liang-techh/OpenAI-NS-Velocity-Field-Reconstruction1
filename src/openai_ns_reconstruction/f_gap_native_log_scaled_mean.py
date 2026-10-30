"""Caller-owned log-scaled ordinary-``Z`` arithmetic for the F-gap mean.

The static source formula in
``experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py``
has, on the gap chart,

``m = -Pstar*u*(C1*W11 + C2*W12) * exp(phi)``

with

``phi = log_end_scale + (13/2 - xi)/mu - 13``.

This module implements only that formula as ordinary derivative jets through
order three.  It never imports a producer, follows a graph, evaluates a source
function, calls the older Q adapter, converts to ``float``, extracts a
midpoint, or materializes ``exp(phi0)``.  Every source quantity is supplied by
the caller as four explicit ordinary-``Z`` rows.  A non-primitive caller-
directed interval object is retained by the arithmetic operators; an opaque
integer cannot accidentally be interpreted as a source value.

The fixed-``xi`` coordinate is a caller premise.  Consequently ``xi`` has no
ordinary-``Z`` rows in this module.  That is distinct from a moving-endpoint
or radius transport calculation.

The output convention is

``m^(k) = exp(log_scale) * scaled_m[k]`` for ``k = 0, 1, 2, 3``.

The optional ``scaled_my`` rows are the conditional gap identity
``m_y = -m`` through rows 0, 1, and 2.  They are explicitly marked as
conditioned output and are not an independent source proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Any, Callable, Generic, Mapping, Sequence, TypeVar


Scalar = TypeVar("Scalar")
SCHEMA = "f_gap_native_log_scaled_mean_v1"
FIXED_COORDINATE = "fixed_xi"


def _reject_primitive(name: str, value: Any) -> None:
    """Reject values that could hide a graph id or force binary rounding.

    Fractions are the exact built-in scalar accepted by the contract.  Other
    non-primitive objects are intentionally opaque caller-directed interval
    values; their arithmetic is never inspected, narrowed, or converted.
    """

    if isinstance(value, bool) or isinstance(value, int):
        raise TypeError(
            f"{name} must be Fraction or a caller-directed interval object; "
            "plain integers are not accepted as source values or node ids"
        )
    if isinstance(value, float):
        raise TypeError(f"{name} must not be float; pass Fraction or an interval object")


@dataclass(frozen=True, slots=True)
class OrdinaryZJet3(Generic[Scalar]):
    """Four explicit ordinary derivatives at fixed ``xi``.

    Rows are ordinary derivatives, not divided Taylor coefficients.  The
    constructor intentionally has no omitted-row defaults.
    """

    rows: tuple[Scalar, Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 4:
            raise ValueError("OrdinaryZJet3 requires exactly four explicit rows")
        for order, value in enumerate(raw):
            _reject_primitive(f"jet row {order}", value)
        object.__setattr__(self, "rows", raw)

    @property
    def z0(self) -> Scalar:
        return self.rows[0]

    @property
    def z1(self) -> Scalar:
        return self.rows[1]

    @property
    def z2(self) -> Scalar:
        return self.rows[2]

    @property
    def z3(self) -> Scalar:
        return self.rows[3]


def make_ordinary_z_jet3(z0: Scalar, z1: Scalar, z2: Scalar, z3: Scalar) -> OrdinaryZJet3[Scalar]:
    """Build a jet while requiring all four rows explicitly."""

    return OrdinaryZJet3((z0, z1, z2, z3))


def _coerce_jet(value: OrdinaryZJet3[Scalar] | Sequence[Scalar], name: str) -> OrdinaryZJet3[Scalar]:
    if isinstance(value, OrdinaryZJet3):
        return value
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must provide four ordinary Z rows")
    try:
        rows = tuple(value)
    except TypeError as exc:
        raise TypeError(f"{name} must provide four ordinary Z rows") from exc
    try:
        return OrdinaryZJet3(rows)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"{name}: {exc}") from exc


def _binomial(n: int, k: int) -> Fraction:
    return Fraction(comb(n, k), 1)


def _add_values(left: Scalar, right: Scalar) -> Scalar:
    return left + right  # type: ignore[operator,return-value]


def _sub_values(left: Scalar, right: Scalar) -> Scalar:
    return left - right  # type: ignore[operator,return-value]


def _neg_value(value: Scalar) -> Scalar:
    return -value  # type: ignore[operator,return-value]


def _mul_values(left: Scalar, right: Scalar) -> Scalar:
    return left * right  # type: ignore[operator,return-value]


def _div_values(left: Scalar, right: Scalar) -> Scalar:
    return left / right  # type: ignore[operator,return-value]


def _sum_terms(terms: Sequence[Scalar]) -> Scalar:
    if not terms:
        raise ValueError("at least one ordinary-jet term is required")
    total = terms[0]
    for term in terms[1:]:
        total = _add_values(total, term)
    return total


def _jet_add(*jets: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    if not jets:
        raise ValueError("jet addition requires at least one operand")
    return OrdinaryZJet3(
        tuple(_sum_terms([jet.rows[order] for jet in jets]) for order in range(4))  # type: ignore[arg-type]
    )


def _jet_sub(left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(
        tuple(_sub_values(left.rows[order], right.rows[order]) for order in range(4))  # type: ignore[arg-type]
    )


def _jet_neg(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(tuple(_neg_value(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet_product(left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    """Ordinary Leibniz product with the required binomial factors."""

    rows: list[Scalar] = []
    for order in range(4):
        terms = [
            _mul_values(
                _mul_values(_binomial(order, split), left.rows[split]),
                right.rows[order - split],
            )
            for split in range(order + 1)
        ]
        rows.append(_sum_terms(terms))
    return OrdinaryZJet3(tuple(rows))  # type: ignore[arg-type]


def _jet_reciprocal(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    """Return the ordinary derivative jet of ``1 / jet``.

    The recurrence differentiates ``jet * reciprocal = 1`` and therefore keeps
    the ordinary-derivative binomial factors without using a Taylor
    coefficient convention.
    """

    rows: list[Scalar] = [_div_values(Fraction(1, 1), jet.rows[0])]
    for order in range(1, 4):
        terms = [
            _mul_values(
                _mul_values(_binomial(order, split), jet.rows[split]),
                rows[order - split],
            )
            for split in range(1, order + 1)
        ]
        rows.append(_div_values(_neg_value(_sum_terms(terms)), jet.rows[0]))
    return OrdinaryZJet3(tuple(rows))  # type: ignore[arg-type]


def _jet_constant(value: Scalar) -> OrdinaryZJet3[Scalar]:
    _reject_primitive("fixed scalar", value)
    zero = Fraction(0, 1)
    return OrdinaryZJet3((value, zero, zero, zero))


@dataclass(frozen=True, slots=True)
class FGapNativeLogScaledMean(Generic[Scalar]):
    """Unmaterialized source-formula result for the gap own mean."""

    amplitude: OrdinaryZJet3[Scalar]
    phi: OrdinaryZJet3[Scalar]
    scaled_m: OrdinaryZJet3[Scalar]
    log_scale: Scalar
    xi: Scalar
    scaled_my: tuple[Scalar, Scalar, Scalar]
    coordinate: str = FIXED_COORDINATE
    schema: str = SCHEMA
    pstar_applied_once: bool = True
    scale_materialized: bool = False
    source_formula_conditioned: bool = True
    source_identity_proved: bool = False
    source_values_authenticated: bool = False
    gap_U_zero_conditioned: bool = True
    gap_Uy_zero_conditioned: bool = True
    my_identity_conditioned: bool = True
    own_history_not_zeroed: bool = True

    @property
    def scaled_m_y(self) -> tuple[Scalar, Scalar, Scalar]:
        """Alias naming the conditional mixed-y rows explicitly."""

        return self.scaled_my

    @property
    def actual_m_statement(self) -> str:
        return "m^(k) = exp(log_scale) * scaled_m[k], k=0..3; exp is not materialized"

    def as_dict(self) -> dict[str, Any]:
        """Return metadata and unmaterialized rows without coercion."""

        return {
            "schema": self.schema,
            "coordinate": self.coordinate,
            "xi": self.xi,
            "phi_rows": self.phi.rows,
            "log_scale": self.log_scale,
            "amplitude_rows": self.amplitude.rows,
            "scaled_m_rows": self.scaled_m.rows,
            "scaled_my_rows": self.scaled_my,
            "actual_m": self.actual_m_statement,
            "pstar_applied_once": self.pstar_applied_once,
            "scale_materialized": self.scale_materialized,
            "source_formula_conditioned": self.source_formula_conditioned,
            "source_identity_proved": self.source_identity_proved,
            "source_values_authenticated": self.source_values_authenticated,
            "gap_U_zero_conditioned": self.gap_U_zero_conditioned,
            "gap_Uy_zero_conditioned": self.gap_Uy_zero_conditioned,
            "my_identity_conditioned": self.my_identity_conditioned,
            "own_history_not_zeroed": self.own_history_not_zeroed,
        }


def build_f_gap_native_log_scaled_mean(
    Pstar: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    u: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    C1: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    C2: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    W11: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    W12: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    log_end_scale: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    mu: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    xi: Scalar,
    *,
    denominator_nonzero_gate: Callable[[Scalar], bool],
) -> FGapNativeLogScaledMean[Scalar]:
    """Build the caller-owned fixed-``xi`` log-scaled gap mean.

    ``Pstar`` is included exactly once in the signed amplitude.  All eight
    source jets require four explicit rows.  ``xi`` is a fixed coordinate
    value, so its ordinary-``Z`` rows are intentionally not an input.

    The gate receives ``L = mu[0]``, the actual reciprocal denominator, and
    must return the built-in singleton ``True``.  This is a caller premise; it
    is not a source or interval proof.
    """

    source_jets = {
        "Pstar": _coerce_jet(Pstar, "Pstar"),
        "u": _coerce_jet(u, "u"),
        "C1": _coerce_jet(C1, "C1"),
        "C2": _coerce_jet(C2, "C2"),
        "W11": _coerce_jet(W11, "W11"),
        "W12": _coerce_jet(W12, "W12"),
        "log_end_scale": _coerce_jet(log_end_scale, "log_end_scale"),
        "mu": _coerce_jet(mu, "mu"),
    }
    _reject_primitive("xi", xi)
    if not callable(denominator_nonzero_gate):
        raise TypeError("denominator_nonzero_gate must be callable")

    L = source_jets["mu"].rows[0]
    if denominator_nonzero_gate(L) is not True:
        raise ValueError("denominator_nonzero_gate(mu[0]) must return built-in True")

    coefficient_row = _jet_add(
        _jet_product(source_jets["C1"], source_jets["W11"]),
        _jet_product(source_jets["C2"], source_jets["W12"]),
    )
    amplitude = _jet_neg(
        _jet_product(
            _jet_product(source_jets["Pstar"], source_jets["u"]),
            coefficient_row,
        )
    )

    thirteen_half = Fraction(13, 2)
    numerator = _jet_sub(_jet_constant(thirteen_half), _jet_constant(xi))
    reciprocal_mu = _jet_reciprocal(source_jets["mu"])
    phi = _jet_add(
        source_jets["log_end_scale"],
        _jet_product(numerator, reciprocal_mu),
        _jet_constant(Fraction(-13, 1)),
    )

    phi1 = phi.rows[1]
    phi2 = phi.rows[2]
    phi3 = phi.rows[3]
    bell2 = _add_values(phi2, _mul_values(phi1, phi1))
    bell3 = _add_values(
        _add_values(phi3, _mul_values(Fraction(3, 1), _mul_values(phi1, phi2))),
        _mul_values(_mul_values(phi1, phi1), phi1),
    )
    scaled_m = OrdinaryZJet3(
        (
            amplitude.rows[0],
            _add_values(amplitude.rows[1], _mul_values(amplitude.rows[0], phi1)),
            _add_values(
                _add_values(amplitude.rows[2], _mul_values(Fraction(2, 1), _mul_values(amplitude.rows[1], phi1))),
                _mul_values(amplitude.rows[0], bell2),
            ),
            _add_values(
                _add_values(
                    _add_values(
                        amplitude.rows[3],
                        _mul_values(Fraction(3, 1), _mul_values(amplitude.rows[2], phi1)),
                    ),
                    _mul_values(Fraction(3, 1), _mul_values(amplitude.rows[1], bell2)),
                ),
                _mul_values(amplitude.rows[0], bell3),
            ),
        )
    )

    scaled_my = tuple(_neg_value(scaled_m.rows[order]) for order in range(3))
    return FGapNativeLogScaledMean(
        amplitude=amplitude,
        phi=phi,
        scaled_m=scaled_m,
        log_scale=phi.rows[0],
        xi=xi,
        scaled_my=scaled_my,  # type: ignore[arg-type]
    )


def hand_check_denominator_nonzero(value: Fraction) -> bool:
    """Gate used by the independent exact-rational hand fixture."""

    return value != 0


def hand_check_fixture() -> Mapping[str, Any]:
    """Return a nonzero-derivative rational fixture and independent targets.

    The fixture is caller-owned arithmetic data only; it is not an F source
    identity, producer receipt, or numerical result.  Every supplied source
    jet row is nonzero.  The selected ``log_end_scale`` rows make the exact
    quotient produce ``phi=(1,2,3,4)`` for ``mu=(1/2,1/4,1/8,1/16)`` and
    fixed ``xi=11``.
    """

    f = Fraction
    ones = OrdinaryZJet3((f(1), f(1), f(1), f(1)))
    mu = OrdinaryZJet3((f(1, 2), f(1, 4), f(1, 8), f(1, 16)))
    log_scale = OrdinaryZJet3((f(23), f(-5, 2), f(21, 4), f(23, 8)))
    inputs = {
        "Pstar": ones,
        "u": ones,
        "C1": ones,
        "C2": ones,
        "W11": ones,
        "W12": ones,
        "log_end_scale": log_scale,
        "mu": mu,
        "xi": f(11),
        "denominator_nonzero_gate": hand_check_denominator_nonzero,
    }
    return {
        "inputs": inputs,
        "expected": {
            "amplitude": (f(-2), f(-8), f(-32), f(-128)),
            "phi": (f(1), f(2), f(3), f(4)),
            "log_scale": f(1),
            "scaled_m": (f(-2), f(-12), f(-78), f(-548)),
            "scaled_my": (f(2), f(12), f(78)),
            "actual_m_statement": "m^(k) = exp(log_scale) * scaled_m[k], k=0..3; exp is not materialized",
        },
        "source_identity_proved": False,
        "source_values_authenticated": False,
        "coordinate": FIXED_COORDINATE,
    }


__all__ = [
    "FIXED_COORDINATE",
    "FGapNativeLogScaledMean",
    "OrdinaryZJet3",
    "SCHEMA",
    "build_f_gap_native_log_scaled_mean",
    "hand_check_denominator_nonzero",
    "hand_check_fixture",
    "make_ordinary_z_jet3",
]
