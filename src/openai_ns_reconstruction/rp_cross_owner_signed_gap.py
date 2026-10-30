"""Conditional signed gap between the old and current ``Rp`` owners.

The caller supplies two opaque tokens for the separate factors
``C = exp(logRc_old)`` and ``B = exp(Tw_old-2)``; their meaning and positivity
are caller premises. This leaf never evaluates either exponential or forms a subtraction of
the two physical radii.  The only evaluated expressions are the caller-owned
correlated ports ``d``, ``t``, ``f``, ``g``, ``h`` and ``dZ``::

    logscale = logRc_old + Tw_old - 2
    delta = d + t + d*t
    delta_Z = (h + f)*delta + dZ*(1 + t) + (1 + d)*(1 + t)*g

Here ``d`` is a direct correlated mantissa, ``t`` is the caller's ``expm1``
port, and ``dZ`` is the direct ordinary-``Z1`` mantissa port.  They are not
reconstructed from independent values.  Source ownership, physical units,
endpoint identity, interval/error certification and global claims therefore
remain false.  The supplied historical error ledger is returned unchanged.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any


SCHEMA = "rp_cross_owner_signed_gap_v1"
COMPONENT = SCHEMA
P1_IS_M15 = False
P1_NOT_M15 = True
ROW_CONVENTION = "ordinary Z1 derivative"
SOURCE_RELATIONS = {
    "old_logRp": "logRc_old - 2 + Tw_old",
    "current_Rp": "Rc_current * exp(Tw_current - 2)",
    "compact_quiet_R2": "R2 = 2*Rc",
    "compact_quiet_length": "length = Tw - 2 - ln(2)",
}

SOURCE_OWNER = False
SOURCE_OWNERSHIP = False
PHYSICAL = False
CERTIFIED = False
ENDPOINT = False
GLOBAL = False

_MISSING = object()
_ONE = Fraction(1, 1)
_TWO = Fraction(2, 1)


class RPCrossOwnerSignedGapInputError(TypeError):
    """A caller port is missing or is not an admitted scalar."""


CrossOwnerSignedGapInputError = RPCrossOwnerSignedGapInputError


def _is_midpoint(value: Any) -> bool:
    """Reject explicit midpoint wrappers without inspecting their value."""

    if "midpoint" in type(value).__name__.lower():
        return True
    for name in ("is_midpoint", "_is_midpoint"):
        try:
            if getattr(value, name, False) is True:
                return True
        except Exception:  # pragma: no cover - opaque caller object
            return False
    return isinstance(value, Mapping) and "midpoint" in value


def _scalar(value: Any, name: str) -> Any:
    if value is _MISSING or value is None:
        raise RPCrossOwnerSignedGapInputError(f"missing caller scalar: {name}")
    if isinstance(value, (bool, int, float, complex, str, bytes, bytearray)):
        raise RPCrossOwnerSignedGapInputError(
            f"{name} must be Fraction or a caller-directed scalar"
        )
    if isinstance(value, (list, tuple, set, frozenset, dict)):
        raise RPCrossOwnerSignedGapInputError(f"{name} must be one scalar")
    if _is_midpoint(value):
        raise RPCrossOwnerSignedGapInputError(
            f"{name} must not be a midpoint-selected wrapper"
        )
    return value


def _add(left: Any, right: Any) -> Any:
    try:
        return left + right
    except TypeError:
        return right + left


def _sub(left: Any, right: Any) -> Any:
    try:
        return left - right
    except TypeError:
        return -(right - left)


def _mul(left: Any, right: Any) -> Any:
    try:
        return left * right
    except TypeError:
        return right * left


def _pick(value: Any, names: tuple[str, ...], label: str) -> Any:
    if isinstance(value, Mapping):
        for name in names:
            if name in value:
                return value[name]
    else:
        for name in names:
            try:
                return getattr(value, name)
            except AttributeError:
                continue
    raise RPCrossOwnerSignedGapInputError(
        f"{label} is missing; expected one of {', '.join(names)}"
    )


@dataclass(frozen=True, slots=True)
class RPCrossOwnerSignedGapInput:
    """Caller-owned scale tokens and correlated ordinary-``Z1`` ports."""

    C: Any
    B: Any
    logRc_old: Any
    Tw_old: Any
    d: Any
    t: Any
    f: Any
    g: Any
    h: Any
    dZ: Any
    error_ledger: Any = _MISSING

    def __post_init__(self) -> None:
        for name in ("C", "B", "logRc_old", "Tw_old", "d", "t", "f", "g", "h", "dZ"):
            _scalar(getattr(self, name), name)
        if self.error_ledger is _MISSING or self.error_ledger is None:
            raise RPCrossOwnerSignedGapInputError("missing caller error ledger")

    @property
    def log_rc_old(self) -> Any:
        return self.logRc_old

    @property
    def tw_old(self) -> Any:
        return self.Tw_old

    @property
    def d_Z(self) -> Any:
        return self.dZ


CrossOwnerRadiusInput = RPCrossOwnerSignedGapInput
RPSignedGapInput = RPCrossOwnerSignedGapInput


def _input(value: Any = _MISSING, **fields: Any) -> RPCrossOwnerSignedGapInput:
    if value is not _MISSING and fields:
        raise RPCrossOwnerSignedGapInputError("pass one input object or named ports, not both")
    if isinstance(value, RPCrossOwnerSignedGapInput):
        return value
    raw = value if value is not _MISSING else fields
    aliases = {
        "C": ("C", "c", "C_token", "old_C"),
        "B": ("B", "b", "B_token", "old_B"),
        "logRc_old": ("logRc_old", "log_rc_old", "old_logRc"),
        "Tw_old": ("Tw_old", "tw_old", "old_Tw"),
        "d": ("d", "direct_d", "delta_Rc_over_C"),
        "t": ("t", "expm1_gap", "expm1_Tw_delta", "tw_expm1"),
        "f": ("f", "f_Z", "Tw_old_Z", "tw_old_Z"),
        "g": ("g", "g_Z", "Tw_current_Z_minus_Tw_old_Z", "tw_delta_Z"),
        "h": ("h", "h_Z", "logRc_old_Z", "log_rc_old_Z"),
        "dZ": ("dZ", "d_Z", "delta_Z", "direct_dZ"),
        "error_ledger": ("error_ledger", "ledger", "old_error_ledger"),
    }
    return RPCrossOwnerSignedGapInput(
        **{name: _pick(raw, names, "cross-owner input") for name, names in aliases.items()}
    )


@dataclass(frozen=True, slots=True)
class RPCrossOwnerSignedGapResult:
    """Unmaterialized log-scale plus the signed value and ordinary ``Z1``."""

    logscale: Any
    signed_mantissa_delta: Any
    signed_mantissa_delta_Z1: Any
    position_mantissa: Any
    duration_mantissa: Any
    position_mantissa_Z1: Any
    duration_mantissa_Z1: Any
    C: Any
    B: Any
    d: Any
    t: Any
    f: Any
    g: Any
    h: Any
    dZ: Any
    old_logRp: Any
    error_ledger: Any

    schema: str = SCHEMA
    row_convention: str = ROW_CONVENTION
    ordinary_Z1_is_derivative: bool | None = None
    P1_is_m15: bool = P1_IS_M15
    P1_not_m15: bool = P1_NOT_M15
    conditional_algebra_only: bool = True
    exponential_evaluated: bool = False
    raw_radius_difference_formed: bool = False
    independent_Z_assumed: bool = False
    direct_mantissa_is_bound_or_overlap: bool | None = None
    d_is_function_bound: bool | None = None
    d_is_center: bool | None = None
    interval_overlap_used: bool | None = None
    same_source_callback_proved: bool = False
    same_source_callback_required: bool = True
    cross_owner_function_identity_proved: bool = False
    source_owner: bool = SOURCE_OWNER
    source_ownership: bool = SOURCE_OWNERSHIP
    source_identity_proved: bool = False
    physical: bool = PHYSICAL
    physical_scale_proved: bool = False
    physical_equivalence_proved: bool = False
    certified: bool = CERTIFIED
    interval_certified: bool = False
    error_certified: bool = False
    quadrature_error_certified: bool = False
    endpoint: bool = ENDPOINT
    endpoint_owner_bound: bool = False
    endpoint_identity_proved: bool = False
    moving_endpoint_chain_rule_proved: bool = False
    global_scope: bool = GLOBAL
    global_certified: bool = False
    full_goal_complete: bool = False
    source_relations: Mapping[str, str] = field(
        default_factory=lambda: SOURCE_RELATIONS
    )
    old_error_terms_removed: bool = False
    error_ledger_changed: bool = False
    extra_ledger_terms_added: bool = False
    historical_error_ledger_preserved: bool = True

    @property
    def log_scale(self) -> Any:
        return self.logscale

    @property
    def signed_delta(self) -> Any:
        return self.signed_mantissa_delta

    @property
    def delta(self) -> Any:
        return self.signed_mantissa_delta

    @property
    def delta_Z(self) -> Any:
        return self.signed_mantissa_delta_Z1

    @property
    def ordinary_Z1(self) -> Any:
        return self.signed_mantissa_delta_Z1

    @property
    def Z1(self) -> Any:
        return self.signed_mantissa_delta_Z1

    @property
    def old_log_scale(self) -> Any:
        return self.logscale

    @property
    def old_error_ledger(self) -> Any:
        return self.error_ledger

    def as_dict(self) -> dict[str, Any]:
        """Expose formulas and flags without converting caller objects."""

        data = {
            "schema": self.schema,
            "logscale": self.logscale,
            "log_scale": self.logscale,
            "signed_mantissa_delta": self.signed_mantissa_delta,
            "signed_mantissa_delta_Z1": self.signed_mantissa_delta_Z1,
            "position_mantissa": self.position_mantissa,
            "duration_mantissa": self.duration_mantissa,
            "position_mantissa_Z1": self.position_mantissa_Z1,
            "duration_mantissa_Z1": self.duration_mantissa_Z1,
            "delta": self.delta,
            "delta_Z": self.delta_Z,
            "ordinary_Z1": self.ordinary_Z1,
            "old_logRp": self.old_logRp,
            "scale_tokens": {"C": self.C, "B": self.B},
            "source_relations": self.source_relations,
            "C": self.C,
            "B": self.B,
            "error_ledger": self.error_ledger,
            "row_convention": self.row_convention,
            "ordinary_Z1_is_derivative": self.ordinary_Z1_is_derivative,
            "P1_is_m15": self.P1_is_m15,
            "P1_not_m15": self.P1_not_m15,
        }
        for name in (
            "conditional_algebra_only", "exponential_evaluated",
            "raw_radius_difference_formed", "independent_Z_assumed",
            "direct_mantissa_is_bound_or_overlap", "same_source_callback_proved",
            "d_is_function_bound", "d_is_center", "interval_overlap_used",
            "same_source_callback_required",
            "cross_owner_function_identity_proved", "source_owner", "source_ownership",
            "source_identity_proved", "physical", "physical_scale_proved",
            "physical_equivalence_proved", "certified", "interval_certified",
            "error_certified", "quadrature_error_certified", "endpoint",
            "endpoint_owner_bound", "endpoint_identity_proved",
            "moving_endpoint_chain_rule_proved", "global_scope", "global_certified",
            "full_goal_complete", "old_error_terms_removed", "error_ledger_changed",
            "extra_ledger_terms_added", "historical_error_ledger_preserved",
        ):
            data[name] = getattr(self, name)
        return data


def build_rp_cross_owner_signed_gap(
    inputs: Any = _MISSING, **ports: Any
) -> RPCrossOwnerSignedGapResult:
    """Build the conditional cross-owner signed ``Rp`` gap.

    ``inputs`` may be :class:`RPCrossOwnerSignedGapInput` or a mapping/object
    exposing its fields.  Named ports use the same names and aliases accepted
    by :func:`_input`.  All arithmetic inputs must be ``Fraction`` or a
    caller-directed scalar; the ledger is retained by identity.
    """

    source = _input(inputs, **ports)
    try:
        logscale = _sub(_add(source.logRc_old, source.Tw_old), _TWO)
        old_log_rp = _add(_sub(source.logRc_old, _TWO), source.Tw_old)
        one_plus_t = _add(_ONE, source.t)
        one_plus_d = _add(_ONE, source.d)
        position = _mul(source.d, one_plus_t)
        duration = source.t
        delta = _add(position, duration)
        position_z = _add(_add(_mul(_add(source.h, source.f), position),
                               _mul(source.dZ, one_plus_t)),
                          _mul(_mul(source.d, one_plus_t), source.g))
        duration_z = _add(_mul(_add(source.h, source.f), duration),
                          _mul(one_plus_t, source.g))
        delta_z = _add(
            _add(
                _mul(_add(source.h, source.f), delta),
                _mul(source.dZ, one_plus_t),
            ),
            _mul(_mul(one_plus_d, one_plus_t), source.g),
        )
    except (ArithmeticError, OverflowError, TypeError, ValueError) as exc:
        raise RPCrossOwnerSignedGapInputError(
            f"caller-directed cross-owner arithmetic failed: {exc}"
        ) from exc
    return RPCrossOwnerSignedGapResult(
        logscale=logscale,
        signed_mantissa_delta=delta,
        signed_mantissa_delta_Z1=delta_z,
        position_mantissa=position,
        duration_mantissa=duration,
        position_mantissa_Z1=position_z,
        duration_mantissa_Z1=duration_z,
        C=source.C,
        B=source.B,
        d=source.d,
        t=source.t,
        f=source.f,
        g=source.g,
        h=source.h,
        dZ=source.dZ,
        old_logRp=old_log_rp,
        error_ledger=source.error_ledger,
    )


def evaluate_rp_cross_owner_signed_gap(inputs: Any = _MISSING, **ports: Any) -> RPCrossOwnerSignedGapResult:
    """Alias for :func:`build_rp_cross_owner_signed_gap`."""

    return build_rp_cross_owner_signed_gap(inputs, **ports)


build_cross_owner_radius_signed_gap = build_rp_cross_owner_signed_gap
compute_rp_cross_owner_signed_gap = build_rp_cross_owner_signed_gap
reconstruct_rp_cross_owner_signed_gap = build_rp_cross_owner_signed_gap
RPSignedGapResult = RPCrossOwnerSignedGapResult
CrossOwnerRadiusSignedGapResult = RPCrossOwnerSignedGapResult


__all__ = [
    "SCHEMA", "COMPONENT", "P1_IS_M15", "P1_NOT_M15", "ROW_CONVENTION", "SOURCE_RELATIONS",
    "RPCrossOwnerSignedGapInputError", "CrossOwnerSignedGapInputError",
    "RPCrossOwnerSignedGapInput", "CrossOwnerRadiusInput", "RPSignedGapInput",
    "RPCrossOwnerSignedGapResult", "CrossOwnerRadiusSignedGapResult", "RPSignedGapResult",
    "build_rp_cross_owner_signed_gap", "evaluate_rp_cross_owner_signed_gap",
    "build_cross_owner_radius_signed_gap", "compute_rp_cross_owner_signed_gap",
    "reconstruct_rp_cross_owner_signed_gap",
]
