"""Signed native ``p2`` source with ordinary :math:`Z` jets.

This module is a deliberately small source boundary.  It evaluates the
original pressure numerator from signed coefficient jets and keeps the
``B=phi**2*b`` and ``C=phi**2*c`` normalization until after the exact
algebraic cancellation

``R*F0**2*phi**2 = Utheta**2/2``.

The audited native numerator is

``N = -U + (1-delta)*(M-Z*MZ)/2 + U*((1-delta)*Z*M+d*MZ)``
``    + 2*delta*Z*A-d*AZ + Pstar**2*(2*(1+delta)*Z*P0-d*P0Z)``
``    + R*(-2*delta*Z*F0**2*B``
``       + d*(F0**2_Z*B+F0**2*B_Z)``
``       + 2*(1+delta)*Z*F0**2*C``
``       - d*(F0**2_Z*C+F0**2*C_Z))``.

The public numerical entry point accepts ordinary derivatives (not Taylor
coefficients).  Histories are required through order two because the first
ordinary ``Z`` derivative of ``p2`` consumes ``M_ZZ``, ``A_ZZ``, ``P0_ZZ``,
``b_ZZ`` and ``c_ZZ``.  The result exposes the two signed, independently
scaled pieces

``p2 = (R/Utheta)*H/L + (R*Utheta/2)*S/L``

Here ``H`` is the audited native hydro plus ``Pstar**2*P0`` part and ``S`` is
the normalized signed swirl part.  Large positive scales are never replaced
by a selected signed centre or by an absolute logarithmic cap.  When
materializing their product would be unsafe, ``p2`` is left as ``None`` while
the signed coefficient jets and both positive scale jets remain available.

No accepted producer, owner, checker, recurrence, quadrature, or fixture is
imported here.  The native-control report hash below is provenance only; it is
never hydrated as a numerical value.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MODULE_RELATIVE_PATH = (
    "src/openai_ns_reconstruction/original_signed_native_p2_kernel.py"
)
NATIVE_V2_REPORT_RELATIVE_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
NATIVE_V2_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
GENERIC_PRESSURE_SOURCE_RELATIVE_PATH = (
    "src/openai_ns_reconstruction/original_long_reference_native_controls.py"
)
SAFE_MATERIALIZATION_LOG_LIMIT = 10000.0


class OrdinaryZJet(tuple):
    """A tuple of ordinary ``Z`` derivatives.

    Entries are derivatives, so the second row is ``f_ZZ`` rather than the
    Taylor coefficient ``f_ZZ/2``.  The class intentionally remains a tuple
    so existing source consumers can use indexing and tuple comparisons.
    """

    def __new__(cls, values: Iterable[Any] = ()):
        return super().__new__(cls, tuple(values))

    @property
    def values(self) -> tuple[Any, ...]:
        return tuple(self)

    @property
    def value(self) -> Any:
        return self[0]

    @property
    def first(self) -> Any:
        return self[1]

    @property
    def second(self) -> Any:
        return self[2]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_")


def _zero_like(value: Any) -> Any:
    try:
        return value * 0
    except Exception:
        return 0


def _one_like(value: Any) -> Any:
    try:
        return _zero_like(value) + 1
    except Exception:
        return 1


_SIGNED_RECORD_CAP_KEYS = frozenset(
    {
        "cap",
        "caps",
        "logcap",
        "log_cap",
        "log_absolute_upper_bound",
        "absolute_upper_bound",
        "bound",
        "bounds",
        "majorant",
        "majorants",
        "upper_bound",
        "error_cap",
    }
)


def _record_endpoint(context: Any, value: Any) -> Any:
    if hasattr(value, "_mpf_"):
        return value
    if isinstance(value, (tuple, list)) and len(value) == 4 and all(
        type(item) is int for item in value
    ):
        from mpmath import make_mpf
        return make_mpf(tuple(value))
    return context.mpf(value)


def _decode_signed_record(context: Any, value: Mapping[str, Any], *, name: str) -> Any:
    keys = {str(key).lower() for key in value}
    if _SIGNED_RECORD_CAP_KEYS.intersection(keys):
        raise TypeError(f"{name} is a logarithmic cap, not a signed finite value")
    if "ordinary_z_order" in keys:
        value = value.get("value", value.get("signed_value", value))
        if not isinstance(value, Mapping):
            return _coerce_scalar(value, context=context)
        keys = {str(key).lower() for key in value}
    if "signed_value" in keys:
        marker = value.get("signed_value")
        if marker is True:
            if "value" in value:
                return _coerce_scalar(value["value"], context=context)
            if {
                "lower_exact_mpf_tuple", "upper_exact_mpf_tuple"
            }.issubset(keys):
                low = _record_endpoint(context, value["lower_exact_mpf_tuple"])
                high = _record_endpoint(context, value["upper_exact_mpf_tuple"])
                return context.mpf([low, high])
            raise TypeError(f"{name} signed record has no finite endpoints")
        return _coerce_scalar(marker, context=context)
    if "value" in keys and len(keys) <= 5:
        return _coerce_scalar(value["value"], context=context)
    if {"lower_exact_mpf_tuple", "upper_exact_mpf_tuple"}.issubset(keys):
        low = _record_endpoint(context, value["lower_exact_mpf_tuple"])
        high = _record_endpoint(context, value["upper_exact_mpf_tuple"])
        return context.mpf([low, high])
    if {"lower", "upper"}.issubset(keys):
        low = _record_endpoint(context, value["lower"])
        high = _record_endpoint(context, value["upper"])
        return context.mpf([low, high])
    raise TypeError(f"{name} must be an explicit signed finite value record")


def _coerce_scalar(value: Any, *, context: Any = None) -> Any:
    """Convert a scalar to a directed interval when a context is supplied."""
    if value is None:
        raise TypeError("signed source values cannot be unknown")
    if context is not None and isinstance(value, Mapping):
        return _decode_signed_record(context, value, name="signed value")
    if context is None or _is_interval(value):
        return value
    # MPIntervalContext does not accept Fraction directly.  Division of the
    # two integer point intervals is directed and contains the rational.
    if isinstance(value, Fraction):
        return context.mpf(value.numerator) / context.mpf(value.denominator)
    if isinstance(value, int):
        return context.mpf(value)
    if isinstance(value, float):
        return context.mpf(repr(value))
    if isinstance(value, str):
        return context.mpf(value)
    return value


def _mapping_value(value: Mapping[str, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in value:
            return value[name]
    return None


def _jet(value: Any, *, name: str, context: Any = None,
         variable: bool = False, minimum: int = 3) -> OrdinaryZJet:
    """Normalize a scalar, sequence, or named row mapping to a jet."""
    if isinstance(value, OrdinaryZJet):
        rows = list(value)
    elif isinstance(value, Mapping):
        if "rows" in value:
            rows = list(value["rows"])
        elif "values" in value:
            rows = list(value["values"])
        elif "ordinary_Z" in value:
            rows = list(value["ordinary_Z"])
        elif any(
            key in value
            for key in (
                "signed_value", "value", "ordinary_Z_order",
                "lower_exact_mpf_tuple", "upper_exact_mpf_tuple",
                "lower", "upper",
            )
        ):
            rows = [value]
        else:
            rows = []
            for order in range(minimum):
                if order in value:
                    rows.append(value[order])
                elif str(order) in value:
                    rows.append(value[str(order)])
                else:
                    break
            if len(rows) == 0:
                raise TypeError(f"{name} must be a scalar or ordinary-Z rows")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        rows = list(value)
    else:
        rows = [value]
    if len(rows) == 1:
        if variable:
            rows = [rows[0], _one_like(rows[0]), _zero_like(rows[0])]
        else:
            rows = [rows[0]] + [_zero_like(rows[0]) for _ in range(minimum - 1)]
    if len(rows) < minimum:
        raise ValueError(
            f"{name} requires ordinary-Z rows 0..{minimum - 1}; got {len(rows)}"
        )
    rows = [_coerce_scalar(row, context=context) for row in rows[:minimum]]
    return OrdinaryZJet(rows)


def _constant(value: Any, *, context: Any = None, order: int = 3) -> OrdinaryZJet:
    value = _coerce_scalar(value, context=context)
    return OrdinaryZJet([value] + [_zero_like(value) for _ in range(order - 1)])


def _add(left: Sequence[Any], right: Sequence[Any]) -> OrdinaryZJet:
    return OrdinaryZJet(a + b for a, b in zip(left, right))


def _neg(value: Sequence[Any]) -> OrdinaryZJet:
    return OrdinaryZJet(-a for a in value)


def _sub(left: Sequence[Any], right: Sequence[Any]) -> OrdinaryZJet:
    return _add(left, _neg(right))


def _safe_div(left: Any, right: Any) -> Any:
    """Divide without coercing exact large integers to an overflowing float."""
    if isinstance(left, Fraction) or isinstance(right, Fraction):
        return Fraction(left) / Fraction(right)
    if isinstance(left, int) and isinstance(right, int):
        return Fraction(left, right)
    return left / right


def _mul(left: Sequence[Any], right: Sequence[Any]) -> OrdinaryZJet:
    """Ordinary derivative product through the available order."""
    count = min(len(left), len(right))
    rows: list[Any] = []
    for order in range(count):
        # The binomial factors are for ordinary derivatives, not Taylor rows.
        value = _zero_like(left[0] * right[0])
        for index in range(order + 1):
            value = value + math.comb(order, index) * left[index] * right[order - index]
        rows.append(value)
    return OrdinaryZJet(rows)


def _scale(left: Sequence[Any], scalar: Any) -> OrdinaryZJet:
    if isinstance(scalar, Fraction):
        return OrdinaryZJet(
            _safe_div(value * scalar.numerator, scalar.denominator) for value in left
        )
    return OrdinaryZJet(value * scalar for value in left)


def _inv(value: Sequence[Any], *, name: str = "denominator") -> OrdinaryZJet:
    if not value:
        raise ValueError(f"{name} has no rows")
    if _contains_zero(value[0]):
        raise ZeroDivisionError(f"{name} ordinary-Z0 interval contains zero")
    a0 = value[0]
    one = _one_like(a0)
    rows = [_safe_div(one, a0)]
    if len(value) >= 2:
        rows.append(_safe_div(-value[1], a0 * a0))
    if len(value) >= 3:
        rows.append(
            _safe_div(2 * value[1] * value[1], a0 * a0 * a0)
            - _safe_div(value[2], a0 * a0)
        )
    return OrdinaryZJet(rows)


def _quotient(left: Sequence[Any], right: Sequence[Any], *, name: str = "denominator") -> OrdinaryZJet:
    return _mul(left, _inv(right, name=name))


def _shift(value: Sequence[Any]) -> OrdinaryZJet:
    if len(value) < 3:
        raise ValueError("A history through ordinary Z2 is required before shifting")
    return OrdinaryZJet((value[1], value[2], _zero_like(value[0])))


def _contains_zero(value: Any) -> bool:
    if _is_interval(value):
        lower, upper = value._mpi_
        try:
            from mpmath.libmp import fzero, mpf_cmp
            return mpf_cmp(lower, fzero) <= 0 and mpf_cmp(upper, fzero) >= 0
        except Exception:
            return float(lower) <= 0.0 <= float(upper)
    try:
        return value == 0
    except Exception:
        return False


def _strict_positive(value: Any, *, name: str) -> None:
    if _is_interval(value):
        lower = value._mpi_[0]
        try:
            from mpmath.libmp import fzero, mpf_cmp
            if mpf_cmp(lower, fzero) <= 0:
                raise ValueError(f"{name} must be strictly positive")
        except ImportError:
            if float(lower) <= 0.0:
                raise ValueError(f"{name} must be strictly positive")
    elif value <= 0:
        raise ValueError(f"{name} must be strictly positive")


def _log_magnitude(value: Any) -> float:
    """Return a finite conservative log magnitude for materialization gating."""
    if _is_interval(value):
        endpoints = value._mpi_
        endpoint = endpoints[0]
        try:
            from mpmath import log, make_mpf
            magnitude = max(abs(make_mpf(endpoints[0])), abs(make_mpf(endpoints[1])))
            if magnitude == 0:
                return float("-inf")
            return float(log(magnitude))
        except Exception:
            try:
                return math.log(max(abs(float(endpoint)), 1e-300))
            except Exception:
                return float("inf")
    if isinstance(value, Fraction):
        if value == 0:
            return float("-inf")
        # Avoid converting a huge Fraction to a binary float.
        return math.log(abs(value.numerator)) - math.log(abs(value.denominator))
    try:
        magnitude = abs(value)
        if magnitude == 0:
            return float("-inf")
        return math.log(float(magnitude))
    except Exception:
        try:
            from mpmath import log
            return float(log(abs(value)))
        except Exception:
            return float("inf")


def _positive_log(value: Any) -> Any:
    """Retain a positive scale's logarithm without selecting a coefficient."""
    if _is_interval(value):
        try:
            from mpmath import iv
            return iv.log(value)
        except Exception:
            return None
    if isinstance(value, Fraction):
        try:
            from mpmath import mp
            return mp.log(mp.mpf(value.numerator) / value.denominator)
        except Exception:
            return math.log(abs(value.numerator)) - math.log(abs(value.denominator))
    try:
        from mpmath import mp
        return mp.log(value)
    except Exception:
        return math.log(value)


def _relative_inverse(jet: Sequence[Any]) -> OrdinaryZJet:
    """Return ``jet/jet[0]`` as ordinary rows, with row zero equal to one."""
    a0 = jet[0]
    return OrdinaryZJet(
        [
            _one_like(a0),
            -_safe_div(jet[1], a0),
            _safe_div(2 * jet[1] * jet[1], a0 * a0)
            - _safe_div(jet[2], a0),
        ]
    )


@dataclass(frozen=True)
class PositiveScaleJet:
    """A positive base scale and an independent signed relative coefficient jet."""

    name: str
    base: Any
    log_base: Any
    relative: OrdinaryZJet
    signed_jet: OrdinaryZJet

    @property
    def positive_base_scale(self) -> Any:
        return self.base

    @property
    def signed_coefficient_jet(self) -> OrdinaryZJet:
        return self.relative


def _scale_jet(name: str, base: Any, relative: Sequence[Any]) -> PositiveScaleJet:
    _strict_positive(base, name=name)
    return PositiveScaleJet(
        name=name,
        base=base,
        log_base=_positive_log(base),
        relative=OrdinaryZJet(relative),
        signed_jet=_scale(relative, base),
    )


@dataclass(frozen=True)
class SignedP2SourceInputs:
    """Inputs to the pure signed source evaluator.

    ``U``, ``M``, ``A``, ``P0``, ``b`` and ``c`` are the signed source values
    in the native formula.  ``Pstar`` is applied exactly once as the physical
    ``Pstar**2*P0`` source factor.  Omitting it uses the neutral exact scale
    ``Pstar=1``.
    """

    Z: Any
    delta: Any
    R: Any
    Utheta: Any
    F0: Any
    U: Any
    M: Any
    A: Any
    P0: Any
    b: Any
    c: Any
    Pstar: Any | None = None
    loggrad_Utheta: Any | None = None
    loggrad_F0: Any | None = None

    @classmethod
    def from_mapping(cls, source: Mapping[str, Any]) -> "SignedP2SourceInputs":
        def required(*names: str) -> Any:
            value = _mapping_value(source, names)
            if value is None:
                raise KeyError("Missing signed p2 source field: " + "/".join(names))
            return value

        def optional(*names: str) -> Any:
            return _mapping_value(source, names)

        return cls(
            Z=required("Z", "z"),
            delta=required("delta", "Delta", "δ"),
            R=required("R", "radius"),
            Utheta=required("Utheta", "U_theta", "utheta", "Uθ"),
            F0=required("F0", "f0", "F_0"),
            U=required("U", "u"),
            M=required("M", "m"),
            A=required("A", "a"),
            P0=required("P0", "p0", "pressure"),
            b=required("b", "B_normalized", "b_history"),
            c=required("c", "C_normalized", "c_history"),
            Pstar=optional("Pstar", "P_star", "pstar"),
            loggrad_Utheta=optional("loggrad_Utheta", "Utheta_loggrad", "log_Utheta_Z"),
            loggrad_F0=optional("loggrad_F0", "F0_loggrad", "log_F0_Z"),
        )


@dataclass(frozen=True)
class P2Decomposition:
    """The signed factorized pressure source and its independent scale pieces."""

    H: OrdinaryZJet
    S: OrdinaryZJet
    L: OrdinaryZJet
    d: OrdinaryZJet
    g: OrdinaryZJet
    phi_squared: OrdinaryZJet
    R_over_Utheta: PositiveScaleJet
    R_times_Utheta_half: PositiveScaleJet
    hydro_pressure_piece: OrdinaryZJet
    normalized_swirl_piece: OrdinaryZJet
    factorized_p2: OrdinaryZJet
    direct_original_p2: OrdinaryZJet | None
    original_N: OrdinaryZJet
    normalized_original_N: OrdinaryZJet


@dataclass(frozen=True)
class SignedP2SourceResult:
    """Output of :func:`evaluate_signed_p2_source`."""

    p2: OrdinaryZJet | None
    decomposition: P2Decomposition
    ordinary_Z_orders: tuple[int, int]
    materialized: bool
    materialization_reason: str
    pstar_applied: bool
    source_phi_identity: str
    normalized_swirl_identity: str
    source_values_are_signed: bool = True
    cap_values_promoted: bool = False
    dominant_signed_centre_selected: bool = False

    @property
    def H(self) -> OrdinaryZJet:
        return self.decomposition.H

    @property
    def S(self) -> OrdinaryZJet:
        return self.decomposition.S

    def __getitem__(self, key: str) -> Any:
        if key == "p2":
            return self.p2
        if hasattr(self, key):
            return getattr(self, key)
        if hasattr(self.decomposition, key):
            return getattr(self.decomposition, key)
        raise KeyError(key)

    def as_dict(self) -> dict[str, Any]:
        return {
            "p2": None if self.p2 is None else list(self.p2),
            "ordinary_Z_orders": list(self.ordinary_Z_orders),
            "materialized": self.materialized,
            "materialization_reason": self.materialization_reason,
            "pstar_applied": self.pstar_applied,
            "source_phi_identity": self.source_phi_identity,
            "normalized_swirl_identity": self.normalized_swirl_identity,
            "source_values_are_signed": self.source_values_are_signed,
            "cap_values_promoted": self.cap_values_promoted,
            "dominant_signed_centre_selected": self.dominant_signed_centre_selected,
            "decomposition": {
                "H": list(self.H),
                "S": list(self.S),
                "L": list(self.decomposition.L),
                "d": list(self.decomposition.d),
                "g": list(self.decomposition.g),
                "phi_squared": list(self.decomposition.phi_squared),
                "R_over_Utheta": list(self.decomposition.R_over_Utheta.signed_jet),
                "R_times_Utheta_half": list(self.decomposition.R_times_Utheta_half.signed_jet),
                "hydro_pressure_piece": list(self.decomposition.hydro_pressure_piece),
                "normalized_swirl_piece": list(self.decomposition.normalized_swirl_piece),
                "factorized_p2": list(self.decomposition.factorized_p2),
                "direct_original_p2": (
                    None if self.decomposition.direct_original_p2 is None
                    else list(self.decomposition.direct_original_p2)
                ),
                "original_N": list(self.decomposition.original_N),
                "normalized_original_N": list(self.decomposition.normalized_original_N),
            },
        }


@dataclass(frozen=True)
class SignedScaledP2SourceInputs:
    """Finite signed inputs for the scale-free native coefficient ``K``.

    ``hydro`` is the signed native hydro contribution before division by
    ``Utheta**2``.  ``P0`` is the raw signed pressure history; the native
    operator ``2*(1+delta)*Z*P0-d*P0_Z`` is formed inside the API.  The two
    explicit finite ratios are ``1/Utheta**2`` and the ratio multiplying that
    pressure operator.  For raw normalized pressure use
    ``Pstar_squared_over_Utheta_squared`` as ``pressure_ratio``; when ``P0``
    already contains the physical
    ``Pstar**2`` scale, pass ``inverse_Utheta_squared`` as that ratio and set
    ``pressure_is_physical=True``.  This prevents a second Pstar scaling.
    """

    Z: Any
    delta: Any
    hydro: Any
    P0: Any
    b: Any
    c: Any
    inverse_Utheta_squared: Any | None = None
    pressure_ratio: Any | None = None
    Utheta: Any | None = None
    g: Any | None = None
    loggrad_Utheta: Any | None = None
    pressure_is_physical: bool = False

    @classmethod
    def from_mapping(cls, source: Mapping[str, Any]) -> "SignedScaledP2SourceInputs":
        def required(*names: str) -> Any:
            value = _mapping_value(source, names)
            if value is None:
                raise KeyError("Missing scaled signed p2 field: " + "/".join(names))
            return value

        def optional(*names: str) -> Any:
            return _mapping_value(source, names)

        return cls(
            Z=required("Z", "z"),
            delta=required("delta", "Delta", "δ"),
            hydro=required("hydro", "H_hydro", "hydro_source", "hydro_jets"),
            P0=required("P0", "pressure", "P0_scaled", "pressure_jets"),
            b=required("b", "B_normalized", "b_history"),
            c=required("c", "C_normalized", "c_history"),
            inverse_Utheta_squared=optional(
                "inverse_Utheta_squared", "Utheta_inv_sq", "inv_Utheta2",
                "inverse_utheta_squared",
            ),
            pressure_ratio=optional(
                "pressure_ratio", "Pstar_squared_over_Utheta_squared",
                "Pstar2_over_Utheta2", "Pstar2_over_utheta2",
            ),
            Utheta=optional("Utheta", "U_theta", "utheta", "Uθ"),
            g=optional("g", "Utheta_squared_loggrad", "phi_squared_loggrad"),
            loggrad_Utheta=optional(
                "loggrad_Utheta", "Utheta_loggrad", "log_Utheta_Z"
            ),
            pressure_is_physical=bool(optional("pressure_is_physical")),
        )


@dataclass(frozen=True)
class SignedScaledP2SourceResult:
    """Signed finite coefficient ``K=p2/(R*Utheta/2)`` and its pieces."""

    K: OrdinaryZJet | None
    K_numerator: OrdinaryZJet
    L: OrdinaryZJet
    g: OrdinaryZJet
    S: OrdinaryZJet
    hydro: OrdinaryZJet
    hydro_scaled: OrdinaryZJet
    P0: OrdinaryZJet
    pressure_operator: OrdinaryZJet
    pressure_scaled: OrdinaryZJet
    inverse_Utheta_squared: OrdinaryZJet
    pressure_ratio: OrdinaryZJet
    ordinary_Z_orders: tuple[int, int]
    materialized: bool
    materialization_reason: str
    pressure_is_physical: bool
    hydro_retained: bool = True
    source_values_are_signed: bool = True
    cap_values_promoted: bool = False
    dominant_signed_centre_selected: bool = False

    @property
    def k(self) -> OrdinaryZJet | None:
        return self.K

    @property
    def coefficient(self) -> OrdinaryZJet | None:
        return self.K

    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        raise KeyError(key)

    def as_dict(self) -> dict[str, Any]:
        return {
            "K": None if self.K is None else list(self.K),
            "K_numerator": list(self.K_numerator),
            "L": list(self.L),
            "g": list(self.g),
            "S": list(self.S),
            "hydro": list(self.hydro),
            "hydro_scaled": list(self.hydro_scaled),
            "P0": list(self.P0),
            "pressure_operator": list(self.pressure_operator),
            "pressure_scaled": list(self.pressure_scaled),
            "inverse_Utheta_squared": list(self.inverse_Utheta_squared),
            "pressure_ratio": list(self.pressure_ratio),
            "ordinary_Z_orders": list(self.ordinary_Z_orders),
            "materialized": self.materialized,
            "materialization_reason": self.materialization_reason,
            "pressure_is_physical": self.pressure_is_physical,
            "hydro_retained": self.hydro_retained,
            "source_values_are_signed": self.source_values_are_signed,
            "cap_values_promoted": self.cap_values_promoted,
            "dominant_signed_centre_selected": self.dominant_signed_centre_selected,
        }


def _prepare_inputs(source: SignedP2SourceInputs | Mapping[str, Any], *, context: Any = None) -> dict[str, OrdinaryZJet]:
    if isinstance(source, Mapping):
        source = SignedP2SourceInputs.from_mapping(source)
    if not isinstance(source, SignedP2SourceInputs):
        raise TypeError("signed p2 inputs must be SignedP2SourceInputs or a mapping")
    # A scalar Z denotes the evaluation point of the independent variable and
    # therefore receives ordinary rows Z_Z=1, Z_ZZ=0.
    Z = _jet(source.Z, name="Z", context=context, variable=True)
    delta = _jet(source.delta, name="delta", context=context)
    R = _jet(source.R, name="R", context=context)
    Pstar = _jet(1 if source.Pstar is None else source.Pstar,
                 name="Pstar", context=context)
    Utheta = _jet(source.Utheta, name="Utheta", context=context)
    F0 = _jet(source.F0, name="F0", context=context)
    U = _jet(source.U, name="U", context=context)
    M = _jet(source.M, name="M", context=context)
    A = _jet(source.A, name="A", context=context)
    P0 = _jet(source.P0, name="P0", context=context)
    b = _jet(source.b, name="b normalized history", context=context)
    c = _jet(source.c, name="c normalized history", context=context)
    if any(not _is_zero(row) for row in R[1:]):
        raise ValueError("R must be Z-independent in the signed native p2 source")
    if any(not _is_zero(row) for row in Pstar[1:]):
        raise ValueError("Pstar must be Z-independent in the signed native p2 source")
    _strict_positive(R[0], name="R")
    _strict_positive(Pstar[0], name="Pstar")
    _strict_positive(Utheta[0], name="Utheta")
    _strict_positive(F0[0], name="F0")
    return dict(Z=Z, delta=delta, R=R, Pstar=Pstar, Utheta=Utheta, F0=F0,
                U=U, M=M, A=A, P0=P0, b=b, c=c,
                _loggrad_Utheta=source.loggrad_Utheta,
                _loggrad_F0=source.loggrad_F0,
                _Pstar=source.Pstar)


def _is_zero(value: Any) -> bool:
    if _is_interval(value):
        return value._mpi_[0] == value._mpi_[1] == (0, 0, 0, 0)
    return value == 0


def _finite_scalar(value: Any) -> bool:
    if _is_interval(value):
        try:
            from mpmath import isfinite, make_mpf
            lower, upper = value._mpi_
            return bool(isfinite(make_mpf(lower)) and isfinite(make_mpf(upper)))
        except Exception:
            return True
    if isinstance(value, float):
        return math.isfinite(value)
    try:
        from mpmath import isfinite
        return bool(isfinite(value))
    except Exception:
        return True


def _require_finite_jet(value: Sequence[Any], *, name: str) -> None:
    if not all(_finite_scalar(row) for row in value):
        raise ValueError(f"{name} must contain finite signed values")


def _loggrad_jet(value: Sequence[Any], *, supplied: Any = None,
                 name: str) -> OrdinaryZJet:
    if supplied is not None:
        return _jet(supplied, name=name, minimum=2)
    first = _safe_div(value[1], value[0])
    second = _safe_div(value[2], value[0]) - first * first
    return OrdinaryZJet((first, second, _zero_like(first)))


def _direct_original_p2(data: Mapping[str, OrdinaryZJet], *, phi_squared: Sequence[Any],
                        H: Sequence[Any], d: Sequence[Any], L: Sequence[Any]) -> tuple[OrdinaryZJet, OrdinaryZJet]:
    """Evaluate the uncontracted formula for audit and fixture comparison."""
    R, F0, Utheta = data["R"], data["F0"], data["Utheta"]
    B = _mul(phi_squared, data["b"])
    C = _mul(phi_squared, data["c"])
    BZ = _shift(B)
    CZ = _shift(C)
    F2 = _mul(F0, F0)
    F2Z = _shift(F2)
    Z, delta = data["Z"], data["delta"]
    one = _constant(1, context=None)
    one_plus_delta = _add(one, delta)
    # Keep the audited signed source sectors in their original source shape.
    # In particular, F0^2 is differentiated together with B and C before
    # the exact F0^2*phi^2 normalization is applied.
    swirl = _add(
        _add(
            _neg(_scale(_mul(_mul(delta, Z), _mul(F2, B)), 2)),
            _mul(d, _add(_mul(F2Z, B), _mul(F2, BZ))),
        ),
        _sub(
            _scale(_mul(_mul(one_plus_delta, Z), _mul(F2, C)), 2),
            _mul(d, _add(_mul(F2Z, C), _mul(F2, CZ))),
        ),
    )
    # The first addend above is -2*delta*Z*F0^2*B + d*(F0^2_Z*B+F0^2*B_Z);
    # the second is the corresponding positive C sector.
    original_N = _add(H, _mul(R, swirl))
    p2 = _quotient(_mul(R, original_N), _mul(Utheta, L), name="Utheta*L")
    return OrdinaryZJet(p2), OrdinaryZJet(original_N)


def evaluate_signed_p2_source(
    source: SignedP2SourceInputs | Mapping[str, Any],
    *,
    context: Any = None,
    materialize: bool | None = None,
    audit_uncontracted: bool = False,
) -> SignedP2SourceResult:
    """Evaluate signed ordinary ``p2`` rows 0 and 1.

    ``materialize=None`` applies a conservative logarithmic safety gate;
    ``False`` always returns the decomposition only.  Even with
    ``materialize=True`` a scale beyond the safety gate is kept factored.
    ``context`` may be an ``mpmath.ctx_iv.MPIntervalContext``; raw rational
    inputs then become directed interval points and the output remains an
    enclosure of the signed rational source.
    """
    data = _prepare_inputs(source, context=context)
    Z, delta, R, Pstar, Utheta, F0 = (
        data[key] for key in ("Z", "delta", "R", "Pstar", "Utheta", "F0")
    )
    U, M, A, P0, b, c = (data[key] for key in ("U", "M", "A", "P0", "b", "c"))
    one = _constant(1, context=context)
    d = _sub(one, _mul(Z, Z))
    L = _sub(one, _mul(delta, _mul(Z, Z)))
    MZ, AZ, P0Z, bZ, cZ = (_shift(value) for value in (M, A, P0, b, c))
    one_minus_delta = _sub(one, delta)
    one_plus_delta = _add(one, delta)

    # H is the audited native hydro plus Pstar^2*P0 part of N.  The signs in
    # the M sector and the leading -U are part of the source identity.
    H = _neg(U)
    H = _add(
        H,
        _scale(_mul(one_minus_delta, _sub(M, _mul(Z, MZ))), Fraction(1, 2)),
    )
    H = _add(H, _mul(U, _add(_mul(_mul(one_minus_delta, Z), M), _mul(d, MZ))))
    H = _add(H, _sub(_scale(_mul(_mul(delta, Z), A), 2), _mul(d, AZ)))
    pressure = _sub(
        _scale(_mul(_mul(one_plus_delta, Z), P0), 2),
        _mul(d, P0Z),
    )
    H = _add(H, _mul(_mul(Pstar, Pstar), pressure))

    # g is the ordinary logarithmic derivative of F0^2*phi^2 = Utheta^2/(2R).
    # The F0 derivative cancels before the signed swirl source is formed.
    # The optional Utheta log-gradient port therefore determines g directly;
    # loggrad_F0 is retained only as an input compatibility port.
    log_ut = _loggrad_jet(Utheta, supplied=data["_loggrad_Utheta"], name="loggrad_Utheta")
    g = _scale(log_ut, 2)

    # Normalized swirl S.  Contracting F0^2*phi^2 before any interval boxing
    # is the key source identity; B and C are never independently scaled.
    b_plus_gb = _add(bZ, _mul(g, b))
    c_plus_gc = _add(cZ, _mul(g, c))
    S = _add(
        _add(_neg(_scale(_mul(_mul(delta, Z), b), 2)), _mul(d, b_plus_gb)),
        _sub(_scale(_mul(_mul(one_plus_delta, Z), c), 2), _mul(d, c_plus_gc)),
    )

    # phi^2 is retained as a same-source jet for provenance/direct formula
    # inspection, while the factorized path uses S and never multiplies its
    # large positive scale into the signed history prematurely.
    phi_squared = _quotient(
        _mul(Utheta, Utheta), _scale(_mul(_mul(R, F0), F0), 2), name="2*R*F0^2"
    )

    # Positive scale jets: the relative rows are signed coefficient jets,
    # while base/log_base remain independent positive scale objects.
    inv_ut_relative = _relative_inverse(Utheta)
    ut_relative = OrdinaryZJet((
        _one_like(Utheta[0]),
        _safe_div(Utheta[1], Utheta[0]),
        _safe_div(Utheta[2], Utheta[0]),
    ))
    scale_R_over_ut = _scale_jet(
        "R/Utheta", _safe_div(R[0], Utheta[0]), inv_ut_relative
    )
    scale_R_ut_half = _scale_jet(
        "R*Utheta/2", _safe_div(R[0] * Utheta[0], 2), ut_relative
    )
    inv_L = _inv(L, name="L")
    H_over_L = _mul(H, inv_L)
    S_over_L = _mul(S, inv_L)
    hydro_piece = _mul(scale_R_over_ut.signed_jet, H_over_L)
    swirl_piece = _mul(scale_R_ut_half.signed_jet, S_over_L)
    factorized_p2 = _add(hydro_piece, swirl_piece)

    normalized_original_N = _add(H, _scale(_mul(_mul(Utheta, Utheta), S), Fraction(1, 2)))
    # The production path never independently materializes R*F0^2*B/C.  The
    # direct uncontracted expression is available only as an explicit small
    # fixture audit; the normalized identity is the source N used below.
    direct_p2: OrdinaryZJet | None = None
    original_N = normalized_original_N

    # The signed result itself is materialized only when the two independent
    # positive scales are within the safe log budget.  This gate is semantic,
    # not an arithmetic overflow workaround: no dominant signed term is used.
    scale_logs = [
        _log_magnitude(scale_R_over_ut.base),
        _log_magnitude(scale_R_ut_half.base),
        _log_magnitude(H[0]),
        _log_magnitude(S[0]),
        _log_magnitude(inv_L[0]),
    ]
    safe = all(value <= SAFE_MATERIALIZATION_LOG_LIMIT for value in scale_logs)
    requested = True if materialize is None else bool(materialize)
    do_materialize = requested and safe
    if do_materialize:
        p2 = OrdinaryZJet(factorized_p2[:2])
        reason = "signed factorized coefficients materialized within positive-scale log budget"
    else:
        p2 = None
        reason = (
            "kept factorized because an independent positive scale exceeds the "
            "materialization budget"
            if not safe else "materialize=False requested"
        )
    if audit_uncontracted and safe:
        direct_p2, original_N = _direct_original_p2(
            data, phi_squared=phi_squared, H=H, d=d, L=L
        )
    return SignedP2SourceResult(
        p2=p2,
        decomposition=P2Decomposition(
            H=H, S=S, L=L, d=d, g=g, phi_squared=phi_squared,
            R_over_Utheta=scale_R_over_ut,
            R_times_Utheta_half=scale_R_ut_half,
            hydro_pressure_piece=hydro_piece,
            normalized_swirl_piece=swirl_piece,
            factorized_p2=factorized_p2,
            direct_original_p2=direct_p2,
            original_N=original_N,
            normalized_original_N=normalized_original_N,
        ),
        ordinary_Z_orders=(0, 1),
        materialized=do_materialize,
        materialization_reason=reason,
        pstar_applied=True,
        source_phi_identity="phi=Utheta/(F0*sqrt(2*R)); R*F0^2*phi^2=Utheta^2/2",
        normalized_swirl_identity=(
            "R*[-2*delta*Z*F0^2*B+d*(F0^2_Z*B+F0^2*B_Z) "
            "+2*(1+delta)*Z*F0^2*C-d*(F0^2_Z*C+F0^2*C_Z)] "
            "= (Utheta^2/2)*[-2*delta*Z*b+d*(b_Z+g*b) "
            "+2*(1+delta)*Z*c-d*(c_Z+g*c)]"
        ),
    )


def _prepare_scaled_inputs(
    source: SignedScaledP2SourceInputs | Mapping[str, Any], *, context: Any = None
) -> dict[str, Any]:
    if isinstance(source, Mapping):
        source = SignedScaledP2SourceInputs.from_mapping(source)
    if not isinstance(source, SignedScaledP2SourceInputs):
        raise TypeError(
            "scaled signed p2 inputs must be SignedScaledP2SourceInputs or a mapping"
        )
    Z = _jet(source.Z, name="Z", context=context, variable=True)
    delta = _jet(source.delta, name="delta", context=context)
    # The scaled coefficient only needs hydro Z0/Z1.  It is assembled by the
    # caller from source fields through Z2, but the resulting hydro term has
    # no claimed Z2 row and must not be padded with an invented value.
    hydro = _jet(source.hydro, name="hydro", context=context, minimum=2)
    P0 = _jet(source.P0, name="P0 signed", context=context)
    b = _jet(source.b, name="b normalized history", context=context)
    c = _jet(source.c, name="c normalized history", context=context)
    Utheta = (
        None
        if source.Utheta is None
        else _jet(source.Utheta, name="Utheta", context=context)
    )
    inverse = (
        None
        if source.inverse_Utheta_squared is None
        else _jet(
            source.inverse_Utheta_squared,
            name="inverse_Utheta_squared",
            context=context,
        )
    )
    pressure_ratio = (
        None
        if source.pressure_ratio is None
        else _jet(source.pressure_ratio, name="pressure_ratio", context=context)
    )
    if Utheta is not None:
        _strict_positive(Utheta[0], name="Utheta")
    if inverse is None:
        if Utheta is None:
            raise ValueError(
                "inverse_Utheta_squared is required when Utheta is omitted"
            )
        inverse = _quotient(
            _constant(1, context=context), _mul(Utheta, Utheta),
            name="Utheta^2",
        )
    if pressure_ratio is None:
        if source.pressure_is_physical:
            pressure_ratio = inverse
        else:
            raise ValueError(
                "pressure_ratio is required for normalized P0; pass the finite "
                "Pstar^2/Utheta^2 ratio or mark P0 physical"
            )
    if source.g is not None:
        g = _jet(source.g, name="g", context=context, minimum=2)
    elif Utheta is not None:
        log_ut = _loggrad_jet(
            Utheta, supplied=source.loggrad_Utheta, name="loggrad_Utheta"
        )
        g = _scale(log_ut, 2)
    else:
        raise ValueError("g or Utheta/loggrad_Utheta is required for normalized swirl")
    for name, value in (
        ("hydro", hydro),
        ("P0 signed", P0),
        ("b normalized history", b),
        ("c normalized history", c),
        ("inverse_Utheta_squared", inverse),
        ("pressure_ratio", pressure_ratio),
        ("g", g),
    ):
        _require_finite_jet(value, name=name)
    return {
        "Z": Z,
        "delta": delta,
        "hydro": hydro,
        "P0": P0,
        "b": b,
        "c": c,
        "inverse_Utheta_squared": inverse,
        "pressure_ratio": pressure_ratio,
        "g": g,
        "pressure_is_physical": source.pressure_is_physical,
    }


def evaluate_scaled_signed_p2_source(
    source: SignedScaledP2SourceInputs | Mapping[str, Any] | None = None,
    *,
    context: Any = None,
    materialize: bool | None = None,
    Z: Any | None = None,
    delta: Any | None = None,
    hydro: Any | None = None,
    P0: Any | None = None,
    b: Any | None = None,
    c: Any | None = None,
    inverse_Utheta_squared: Any | None = None,
    pressure_ratio: Any | None = None,
    Utheta: Any | None = None,
    g: Any | None = None,
    loggrad_Utheta: Any | None = None,
    pressure_is_physical: bool = False,
) -> SignedScaledP2SourceResult:
    """Evaluate finite signed ``K=p2/(R*Utheta/2)`` rows 0 and 1.

    This API is for the actual saved source ports where ``P0`` may already
    carry the physical ``Pstar**2`` scale.  It accepts that scale only through
    the explicit finite ``pressure_ratio`` and never applies Pstar again.
    ``hydro`` is retained even when its scaled contribution is tiny.
    """
    if source is None:
        required = {"Z": Z, "delta": delta, "hydro": hydro,
                    "P0": P0, "b": b, "c": c}
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise TypeError("missing scaled signed p2 fields: " + ", ".join(missing))
        source = SignedScaledP2SourceInputs(
            Z=Z,
            delta=delta,
            hydro=hydro,
            P0=P0,
            b=b,
            c=c,
            inverse_Utheta_squared=inverse_Utheta_squared,
            pressure_ratio=pressure_ratio,
            Utheta=Utheta,
            g=g,
            loggrad_Utheta=loggrad_Utheta,
            pressure_is_physical=pressure_is_physical,
        )
    elif any(value is not None for value in (Z, delta, hydro, P0, b, c,
                                              inverse_Utheta_squared,
                                              pressure_ratio, Utheta, g,
                                              loggrad_Utheta)):
        raise TypeError("pass either source or explicit scaled signed p2 fields")
    data = _prepare_scaled_inputs(source, context=context)
    Z, delta = data["Z"], data["delta"]
    hydro, P0 = data["hydro"], data["P0"]
    b, c = data["b"], data["c"]
    inverse = data["inverse_Utheta_squared"]
    pressure_ratio = data["pressure_ratio"]
    g = data["g"]
    one = _constant(1, context=context)
    d = _sub(one, _mul(Z, Z))
    L = _sub(one, _mul(delta, _mul(Z, Z)))
    bZ, cZ = _shift(b), _shift(c)
    P0Z = _shift(P0)
    one_plus_delta = _add(one, delta)
    S = _add(
        _add(
            _neg(_scale(_mul(_mul(delta, Z), b), 2)),
            _mul(d, _add(bZ, _mul(g, b))),
        ),
        _sub(
            _scale(_mul(_mul(one_plus_delta, Z), c), 2),
            _mul(d, _add(cZ, _mul(g, c))),
        ),
    )
    hydro_scaled = _scale(_mul(inverse, hydro), 2)
    pressure_operator = _sub(
        _scale(_mul(_mul(one_plus_delta, Z), P0), 2),
        _mul(d, P0Z),
    )
    pressure_scaled = _scale(_mul(pressure_ratio, pressure_operator), 2)
    K_numerator = _add(S, _add(hydro_scaled, pressure_scaled))
    K_full = _quotient(K_numerator, L, name="L")
    requested = True if materialize is None else bool(materialize)
    if requested:
        K = OrdinaryZJet(K_full[:2])
        reason = "finite signed scaled coefficient materialized"
    else:
        K = None
        reason = "materialize=False requested"
    return SignedScaledP2SourceResult(
        K=K,
        K_numerator=K_numerator,
        L=L,
        g=g,
        S=S,
        hydro=hydro,
        hydro_scaled=hydro_scaled,
        P0=P0,
        pressure_operator=pressure_operator,
        pressure_scaled=pressure_scaled,
        inverse_Utheta_squared=inverse,
        pressure_ratio=pressure_ratio,
        ordinary_Z_orders=(0, 1),
        materialized=requested,
        materialization_reason=reason,
        pressure_is_physical=data["pressure_is_physical"],
    )


# Explicit aliases keep the source boundary easy to find for downstream
# consumers without introducing a second implementation.
finite_directed_iv_p2_source = evaluate_signed_p2_source
directed_iv_p2_source = evaluate_signed_p2_source
evaluate_finite_directed_p2_source = evaluate_signed_p2_source
finite_directed_iv_scaled_p2_source = evaluate_scaled_signed_p2_source
directed_iv_scaled_p2_source = evaluate_scaled_signed_p2_source
evaluate_finite_directed_scaled_p2_source = evaluate_scaled_signed_p2_source


def p2_ordinary_z01(
    source: SignedP2SourceInputs | Mapping[str, Any],
    *, context: Any = None,
    materialize: bool | None = None,
) -> SignedP2SourceResult:
    """Named ordinary-Z0/Z1 facade returning the signed decomposition result."""
    return evaluate_signed_p2_source(source, context=context, materialize=materialize)


def _sympy_symbols() -> dict[str, Any]:
    import sympy as sp

    Z = sp.Symbol("Z")
    delta = sp.Symbol("delta")
    R = sp.Symbol("R", positive=True)
    Utheta = sp.Function("Utheta")(Z)
    F0 = sp.Function("F0")(Z)
    U = sp.Function("U")(Z)
    M = sp.Function("M")(Z)
    A = sp.Function("A")(Z)
    P0 = sp.Function("P0")(Z)
    b = sp.Function("b")(Z)
    c = sp.Function("c")(Z)
    return dict(sp=sp, Z=Z, delta=delta, R=R, Utheta=Utheta, F0=F0,
                U=U, M=M, A=A, P0=P0, b=b, c=c)


def sympy_signed_p2_source() -> dict[str, Any]:
    """Return the full and exact normalized SymPy source expressions."""
    symbols = _sympy_symbols()
    sp = symbols.pop("sp")
    Z, delta, R = symbols["Z"], symbols["delta"], symbols["R"]
    Pstar = sp.Symbol("Pstar", positive=True)
    symbols["Pstar"] = Pstar
    Utheta, F0 = symbols["Utheta"], symbols["F0"]
    U, M, A, P0, b, c = (symbols[name] for name in ("U", "M", "A", "P0", "b", "c"))
    d = 1 - Z**2
    L = 1 - delta * Z**2
    MZ, AZ, P0Z = sp.diff(M, Z), sp.diff(A, Z), sp.diff(P0, Z)
    bZ, cZ = sp.diff(b, Z), sp.diff(c, Z)
    phi_squared = Utheta**2 / (2 * R * F0**2)
    B, C = phi_squared * b, phi_squared * c
    BZ, CZ = sp.diff(B, Z), sp.diff(C, Z)
    H = (
        -U + (1 - delta) * (M - Z * MZ) / 2
        + U * ((1 - delta) * Z * M + d * MZ)
        + 2 * delta * Z * A - d * AZ
        + Pstar**2 * (2 * (1 + delta) * Z * P0 - d * P0Z)
    )
    original_swirl = (
        -2 * delta * Z * F0**2 * B
        + d * (sp.diff(F0**2, Z) * B + F0**2 * BZ)
        + 2 * (1 + delta) * Z * F0**2 * C
        - d * (sp.diff(F0**2, Z) * C + F0**2 * CZ)
    )
    original_N = H + R * original_swirl
    original_p2 = R * original_N / (Utheta * L)
    g = sp.diff(Utheta**2, Z) / Utheta**2
    S = (
        -2 * delta * Z * b + d * (bZ + g * b)
        + 2 * (1 + delta) * Z * c - d * (cZ + g * c)
    )
    normalized_N = H + Utheta**2 * S / 2
    factorized_p2 = R * H / (Utheta * L) + R * Utheta * S / (2 * L)
    return {
        "symbols": symbols,
        "Pstar": Pstar,
        "d": d,
        "L": L,
        "phi_squared": phi_squared,
        "g": g,
        "B": B,
        "C": C,
        "B_Z": BZ,
        "C_Z": CZ,
        "H": H,
        "S": S,
        "original_swirl": original_swirl,
        "original_N": original_N,
        "normalized_original_N": normalized_N,
        "original_p2": original_p2,
        "factorized_p2": factorized_p2,
        "normalization_identity": sp.simplify(
            R * original_swirl - Utheta**2 * S / 2
        ),
        "factorization_identity": sp.simplify(original_p2 - factorized_p2),
        "ordinary_Z_orders": [0, 1],
        "R_is_Z_constant": True,
        "source_phi_identity": "phi=Utheta/(F0*sqrt(2*R))",
    }


sympy_p2_source = sympy_signed_p2_source
build_sympy_signed_p2_source = sympy_signed_p2_source


def sympy_signed_p2_ordinary_z01(
    expression: Any | None = None,
    *,
    z_value: Any | None = None,
) -> tuple[Any, Any]:
    """Return an expression and its ordinary ``Z`` first derivative."""
    import sympy as sp

    source = sympy_signed_p2_source()
    if expression is None:
        expression = source["original_p2"]
    Z = source["symbols"]["Z"]
    rows = (expression, sp.diff(expression, Z))
    if z_value is not None:
        rows = tuple(sp.simplify(row.subs(Z, z_value)) for row in rows)
    return rows


def source_preflight(path: str | Path | None = None) -> dict[str, Any]:
    """Check the source file's LF and single-EOF contract before fixtures."""
    target = Path(path) if path is not None else Path(__file__)
    raw = target.read_bytes()
    lf_only = b"\r" not in raw
    single_eof = raw.endswith(b"\n") and not raw.endswith(b"\n\n")
    return {
        "path": target.as_posix(),
        "lf_only": lf_only,
        "single_eof_newline": single_eof,
        "passed": lf_only and single_eof,
    }


def _rational_interval(context: Any, value: Any) -> Any:
    value = Fraction(value)
    return context.mpf(value.numerator) / context.mpf(value.denominator)


def _sympy_jet_at(sp: Any, expression: Any, Z: Any, point: Any) -> list[Any]:
    return [sp.cancel(sp.diff(expression, Z, order).subs(Z, point)) for order in range(3)]


def _interval_contains(context: Any, interval: Any, rational: Any) -> bool:
    target = _rational_interval(context, rational)
    if not _is_interval(interval):
        return interval == rational
    lower, upper = interval._mpi_
    from mpmath.libmp import mpf_cmp
    return mpf_cmp(lower, target._mpi_[0]) <= 0 and mpf_cmp(target._mpi_[1], upper) <= 0


def rational_sympy_fixture() -> dict[str, Any]:
    """Run a fresh rational source/factorization and directed-IV containment fixture.

    The preflight is intentionally the first operation.  The fixture uses
    signed ``b`` and ``c`` histories with nonzero second derivatives and a
    nonconstant ``F0``.  It does not read or execute any accepted report.
    """
    preflight = source_preflight()
    if not preflight["passed"]:
        return {"passed": False, "preflight": preflight,
                "module_sha256": _module_sha256()}
    import sympy as sp
    from mpmath.ctx_iv import MPIntervalContext

    Z = sp.Symbol("Z")
    z0 = sp.Rational(2, 7)
    delta = sp.Rational(2, 9)
    R = sp.Rational(5, 3)
    Pstar = sp.Rational(7, 5)
    Utheta = 3 + 2 * Z + Z**3
    F0 = 2 + Z + Z**2
    U = -1 + 3 * Z - Z**2
    M = 2 - Z + 2 * Z**2 - Z**3
    A = -3 + Z + Z**2 + 2 * Z**3
    P0 = 1 - 2 * Z + Z**2 + Z**3
    b = -2 + 3 * Z + 5 * Z**2 - Z**3
    c = 1 - 4 * Z + 2 * Z**2 + 3 * Z**3
    d = 1 - Z**2
    L = 1 - delta * Z**2
    phi2 = Utheta**2 / (2 * R * F0**2)
    B, C = phi2 * b, phi2 * c
    H = (
        -U + (1 - delta) * (M - Z * sp.diff(M, Z)) / 2
        + U * ((1 - delta) * Z * M + d * sp.diff(M, Z))
        + 2 * delta * Z * A - d * sp.diff(A, Z)
        + Pstar**2 * (2 * (1 + delta) * Z * P0 - d * sp.diff(P0, Z))
    )
    pressure_source = 2 * (1 + delta) * Z * P0 - d * sp.diff(P0, Z)
    hydro_source = H - Pstar**2 * pressure_source
    original_N = H + R * (
        -2 * delta * Z * F0**2 * B
        + d * (sp.diff(F0**2, Z) * B + F0**2 * sp.diff(B, Z))
        + 2 * (1 + delta) * Z * F0**2 * C
        - d * (sp.diff(F0**2, Z) * C + F0**2 * sp.diff(C, Z))
    )
    original_p2 = R * original_N / (Utheta * L)
    g = sp.diff(Utheta**2, Z) / Utheta**2
    S = (
        -2 * delta * Z * b + d * (sp.diff(b, Z) + g * b)
        + 2 * (1 + delta) * Z * c - d * (sp.diff(c, Z) + g * c)
    )
    factorized_p2 = R * H / (Utheta * L) + R * Utheta * S / (2 * L)
    exact_factorization = sp.simplify(original_p2 - factorized_p2) == 0
    exact_normalization = sp.simplify(
        R * (
            -2 * delta * Z * F0**2 * B
            + d * (sp.diff(F0**2, Z) * B + F0**2 * sp.diff(B, Z))
            + 2 * (1 + delta) * Z * F0**2 * C
            - d * (sp.diff(F0**2, Z) * C + F0**2 * sp.diff(C, Z))
        ) - Utheta**2 * S / 2
    ) == 0
    original_rows = [sp.factor(sp.diff(original_p2, Z, order).subs(Z, z0)) for order in (0, 1)]
    fact_rows = [sp.factor(sp.diff(factorized_p2, Z, order).subs(Z, z0)) for order in (0, 1)]
    exact_rows = all(sp.simplify(a - b_) == 0 for a, b_ in zip(original_rows, fact_rows))
    f0_z_nonzero = sp.diff(F0, Z).subs(Z, z0) != 0
    b_c_z2_nonzero = sp.diff(b, Z, 2).subs(Z, z0) != 0 and sp.diff(c, Z, 2).subs(Z, z0) != 0

    context = MPIntervalContext()
    context.prec = 256
    inputs = SignedP2SourceInputs(
        Z=[_rational_interval(context, z0), _rational_interval(context, 1), _rational_interval(context, 0)],
        delta=_rational_interval(context, delta),
        R=_rational_interval(context, R),
        Utheta=[_rational_interval(context, x) for x in _sympy_jet_at(sp, Utheta, Z, z0)],
        F0=[_rational_interval(context, x) for x in _sympy_jet_at(sp, F0, Z, z0)],
        U=[_rational_interval(context, x) for x in _sympy_jet_at(sp, U, Z, z0)],
        M=[_rational_interval(context, x) for x in _sympy_jet_at(sp, M, Z, z0)],
        A=[_rational_interval(context, x) for x in _sympy_jet_at(sp, A, Z, z0)],
        P0=[_rational_interval(context, x) for x in _sympy_jet_at(sp, P0, Z, z0)],
        Pstar=_rational_interval(context, Pstar),
        b=[_rational_interval(context, x) for x in _sympy_jet_at(sp, b, Z, z0)],
        c=[_rational_interval(context, x) for x in _sympy_jet_at(sp, c, Z, z0)],
    )
    result = evaluate_signed_p2_source(inputs, context=context, materialize=True)
    interval_containment = result.p2 is not None and all(
        _interval_contains(context, result.p2[index], value)
        for index, value in enumerate(original_rows)
    )
    def rational_jet(expression: Any) -> list[Fraction]:
        return [
            Fraction(sp.cancel(sp.diff(expression, Z, order).subs(Z, z0)))
            for order in range(3)
        ]

    inverse_u2 = 1 / Utheta**2
    pressure_ratio = Pstar**2 / Utheta**2
    scaled_K = sp.cancel((S + 2 * hydro_source * inverse_u2
                          + 2 * pressure_ratio * pressure_source) / L)
    scaled_rows = [
        sp.cancel(sp.diff(scaled_K, Z, order).subs(Z, z0))
        for order in (0, 1)
    ]
    hydro_rows = rational_jet(hydro_source)[:2]
    scaled_inputs = SignedScaledP2SourceInputs(
        Z=[Fraction(z0), Fraction(1), Fraction(0)],
        delta=Fraction(delta),
        hydro=hydro_rows,
        P0=rational_jet(P0),
        b=rational_jet(b),
        c=rational_jet(c),
        inverse_Utheta_squared=rational_jet(inverse_u2),
        pressure_ratio=rational_jet(pressure_ratio),
        g=rational_jet(g),
    )
    scaled_result = evaluate_scaled_signed_p2_source(scaled_inputs)
    scaled_exact = scaled_result.K is not None and all(
        sp.simplify(scaled_result.K[index] - scaled_rows[index]) == 0
        for index in (0, 1)
    )
    scaled_physical_inputs = SignedScaledP2SourceInputs(
        Z=[Fraction(z0), Fraction(1), Fraction(0)],
        delta=Fraction(delta),
        hydro=hydro_rows,
        P0=rational_jet(Pstar**2 * P0),
        b=rational_jet(b),
        c=rational_jet(c),
        inverse_Utheta_squared=rational_jet(inverse_u2),
        pressure_is_physical=True,
        g=rational_jet(g),
    )
    scaled_physical_result = evaluate_scaled_signed_p2_source(
        scaled_physical_inputs
    )
    scaled_physical_exact = scaled_physical_result.K is not None and all(
        sp.simplify(scaled_physical_result.K[index] - scaled_rows[index]) == 0
        for index in (0, 1)
    )
    # The final field is intentionally a receipt, not a production report.
    checks = {
        "preflight": preflight["passed"],
        "exact_factorization": exact_factorization,
        "exact_normalization": exact_normalization,
        "exact_Z0_Z1_rows": exact_rows,
        "nonzero_F0_Z": f0_z_nonzero,
        "signed_b_c_Z2": b_c_z2_nonzero,
        "directed_iv_contains_signed_rational": interval_containment,
        "scaled_coefficient_exact": scaled_exact,
        "scaled_physical_P0_no_double_Pstar": scaled_physical_exact,
        "scaled_hydro_uses_only_Z0_Z1": scaled_result.hydro == tuple(hydro_rows),
        "ordinary_Z_orders": result.ordinary_Z_orders == (0, 1),
        "R_constant": True,
        "Pstar_applied_once": result.pstar_applied is True,
        "source_values_are_signed": result.source_values_are_signed,
        "cap_values_not_promoted": result.cap_values_promoted is False,
        "dominant_signed_centre_not_selected": result.dominant_signed_centre_selected is False,
    }
    return {
        "passed": all(checks.values()),
        "module_sha256": _module_sha256(),
        "ordinary_Z_orders": [0, 1],
        "checks": checks,
        "preflight": preflight,
        "source_phi_identity": "phi=Utheta/(F0*sqrt(2*R)); R*F0^2*phi^2=Utheta^2/2",
        "materialized": result.materialized,
        "decomposition_preserved": result.p2 is not None,
        "accepted_report_read": False,
        "accepted_owner_checker_fixture_recurrence_quadrature_called": False,
    }


def _module_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def write_rational_sympy_fixture(path: str | Path = "tmp/signed_native_p2_fixture.json") -> dict[str, Any]:
    """Write the optional scratch fixture; no production report is changed."""
    result = rational_sympy_fixture()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return result


run_rational_sympy_fixture = rational_sympy_fixture


__all__ = [
    "MODULE_RELATIVE_PATH",
    "NATIVE_V2_REPORT_RELATIVE_PATH",
    "NATIVE_V2_REPORT_SHA256",
    "GENERIC_PRESSURE_SOURCE_RELATIVE_PATH",
    "OrdinaryZJet",
    "PositiveScaleJet",
    "SignedP2SourceInputs",
    "P2Decomposition",
    "SignedP2SourceResult",
    "SignedScaledP2SourceInputs",
    "SignedScaledP2SourceResult",
    "evaluate_signed_p2_source",
    "finite_directed_iv_p2_source",
    "directed_iv_p2_source",
    "evaluate_finite_directed_p2_source",
    "evaluate_scaled_signed_p2_source",
    "finite_directed_iv_scaled_p2_source",
    "directed_iv_scaled_p2_source",
    "evaluate_finite_directed_scaled_p2_source",
    "p2_ordinary_z01",
    "sympy_signed_p2_source",
    "sympy_p2_source",
    "build_sympy_signed_p2_source",
    "sympy_signed_p2_ordinary_z01",
    "source_preflight",
    "rational_sympy_fixture",
    "run_rational_sympy_fixture",
    "write_rational_sympy_fixture",
]


if __name__ == "__main__":
    print(json.dumps(rational_sympy_fixture(), indent=2, sort_keys=True))
