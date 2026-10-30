"""Direct signed ``P-F`` endpoint algebra for an ordinary ``Z1`` jet.

The endpoint primitive is supplied as a frame value plus a correlated signed
difference.  With ``t = expm1(d)`` and

``H = dQ_y + t * (Q_yF + dQ_y)``,

the endpoint difference is ``E = R_F * H``.  Its ordinary slow-``Z``
derivative is

``H_Z = dQ_y_Z + (1 + t) * d_Z * (Q_yF + dQ_y)``
``      + t * (Q_yF_Z + dQ_y_Z)``

For a fixed endpoint coordinate ``xi``, use the four caller-supplied partials
``R_F_Z_at_xi``, ``d_Z_at_xi``, ``Q_y_F_Z_at_xi`` and ``dQ_y_Z_at_xi``:

``E_Z|xi = R_F_Z_at_xi * H + R_F * H_Z|xi``.

If the endpoint moves with speed ``xi_e_Z``, its additional correction is
``xi_e_Z * E_xi``, where ``E_xi = R_F_xi * H + R_F * H_xi`` and ``H_xi``
retains all four caller-supplied ``xi`` derivatives.  A total derivative is
returned only when ``xi_e_Z`` and all four ``xi`` derivatives are supplied;
missing moving data is never treated as zero.

``d`` and ``expm1_d`` are source-owned inputs.  This module deliberately
does not call an exponential routine: a caller must provide a high-
precision or interval ``expm1(d)`` value, and may provide a safety bound.
Arithmetic stays in the caller's exact or interval domain; there is no
midpoint conversion and no conversion through ``float``.

The input contract is intentionally narrow.  ``R_F``, ``Q_yF`` and all
ordinary ``Z`` derivatives come from the caller's frame/source adapter,
while ``dQ_y`` and ``dQ_y_Z`` are supplied signed differences.  No ``P``
value is accepted and no large ``P-F`` subtraction is performed here.
Missing ``y``-derivatives are rejected.  ``interval_rounding_certified`` is
recorded only as a caller assertion; this module's ``certified`` field is
always false.  No field proves source identity, authentication, or
physical-source binding.

The hand checks in :data:`HAND_CHECK_EXAMPLES` use ``d = 0`` so the caller-
provided ``expm1_d = 0`` is exact.  One has nonzero ``d_Z_at_xi``; another
supplies complete moving-coordinate data.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from decimal import Decimal
from numbers import Number, Rational
from typing import Any


COMPONENT = "pulse_signed_endpoint_interval_v2"
SCHEMA = COMPONENT
_MISSING = object()


class EndpointInputError(ValueError):
    """Caller data cannot be evaluated safely in one numeric domain."""


def _interval_like(value: Any) -> bool:
    """Recognize common interval-shaped values without importing a backend."""

    if getattr(value, "_mpi_", None) is not None:
        return True
    for name in ("lower", "upper", "inf", "sup", "lo", "hi"):
        if hasattr(value, name):
            return True
    return False


def _call_bound(value: Any, names: tuple[str, ...]) -> Any | None:
    """Return an exposed interval endpoint when the backend has one."""

    for name in names:
        if not hasattr(value, name):
            continue
        endpoint = getattr(value, name)
        try:
            endpoint = endpoint() if callable(endpoint) else endpoint
        except (TypeError, ValueError, ArithmeticError):
            continue
        return endpoint
    return None


def _validate_value(value: Any, label: str) -> bool:
    """Validate one exact scalar or interval without changing its type.

    The return value records whether ``value`` looks interval-valued.  Built
    in binary floats are rejected because silently accepting them would make
    this backend less precise than the caller's declared contract.
    """

    if value is _MISSING:
        raise EndpointInputError(f"{label} is required")
    if isinstance(value, bool) or isinstance(value, complex):
        raise EndpointInputError(f"{label} must be a finite real exact/interval value")
    if isinstance(value, float):
        raise EndpointInputError(
            f"{label} must use an exact or interval backend; binary floats are unsupported"
        )

    is_interval = _interval_like(value)
    if is_interval:
        finite = getattr(value, "is_finite", None)
        if finite is not None:
            try:
                finite = finite() if callable(finite) else finite
            except (TypeError, ValueError, ArithmeticError) as exc:
                raise EndpointInputError(f"{label} finiteness is undecidable") from exc
            if isinstance(finite, bool) and not finite:
                raise EndpointInputError(f"{label} must be finite")
        return True

    if isinstance(value, Decimal):
        if not value.is_finite():
            raise EndpointInputError(f"{label} must be finite")
        return False
    if isinstance(value, (int, Rational)):
        return False

    # mpmath-like scalar values expose ``_mpf_``.  Avoid importing mpmath or
    # any repository module; reject its obvious non-finite spellings.
    if getattr(value, "_mpf_", None) is not None:
        text = str(value).strip().lower()
        if text in {
            "nan", "+nan", "-nan", "inf", "+inf", "-inf",
            "infinity", "+infinity", "-infinity",
        }:
            raise EndpointInputError(f"{label} must be finite")
        return False

    finite = getattr(value, "is_finite", None)
    if finite is not None:
        try:
            finite = finite() if callable(finite) else finite
        except (TypeError, ValueError, ArithmeticError) as exc:
            raise EndpointInputError(f"{label} finiteness is undecidable") from exc
        if isinstance(finite, bool) and not finite:
            raise EndpointInputError(f"{label} must be finite")
        return False

    if isinstance(value, Number):
        # Exact third-party real types may register as Number.  A reflexive
        # equality check catches NaN without coercing the value.
        try:
            if not bool(value == value):
                raise EndpointInputError(f"{label} must be finite")
        except EndpointInputError:
            raise
        except (TypeError, ValueError, ArithmeticError) as exc:
            raise EndpointInputError(f"{label} finiteness is undecidable") from exc
        return False

    raise EndpointInputError(
        f"{label} must be a finite exact scalar or caller-owned interval value"
    )


def _same_value(left: Any, right: Any) -> bool:
    """Return a definite equality result, never an interval truth value."""

    if left is right:
        return True
    if _interval_like(left) and _interval_like(right):
        left_lower = _call_bound(left, ("lower", "inf", "lo", "a"))
        left_upper = _call_bound(left, ("upper", "sup", "hi", "b"))
        right_lower = _call_bound(right, ("lower", "inf", "lo", "a"))
        right_upper = _call_bound(right, ("upper", "sup", "hi", "b"))
        if None not in (left_lower, left_upper, right_lower, right_upper):
            return _same_value(left_lower, right_lower) and _same_value(left_upper, right_upper)
    try:
        result = left == right
    except (TypeError, ValueError, ArithmeticError):
        return False
    return isinstance(result, bool) and result


def _exact_zero(value: Any) -> bool:
    try:
        result = value == 0
    except (TypeError, ValueError, ArithmeticError):
        return False
    return isinstance(result, bool) and result


def _require_positive(value: Any, label: str) -> None:
    """Require a radius interval to lie strictly above zero."""

    if _interval_like(value):
        lower = _call_bound(value, ("lower", "inf", "lo", "a"))
        if lower is not None:
            try:
                if not bool(lower > 0):
                    raise EndpointInputError(f"{label} must be strictly positive")
                return
            except EndpointInputError:
                raise
            except (TypeError, ValueError, ArithmeticError) as exc:
                raise EndpointInputError(f"{label} positivity is undecidable") from exc
    try:
        positive = value > 0
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise EndpointInputError(f"{label} positivity is undecidable") from exc
    if not isinstance(positive, bool) or not positive:
        raise EndpointInputError(f"{label} must be strictly positive")


def _require_nonnegative(value: Any, label: str) -> None:
    if _interval_like(value):
        lower = _call_bound(value, ("lower", "inf", "lo", "a"))
        if lower is not None:
            try:
                if not bool(lower >= 0):
                    raise EndpointInputError(f"{label} must be nonnegative")
                return
            except EndpointInputError:
                raise
            except (TypeError, ValueError, ArithmeticError) as exc:
                raise EndpointInputError(f"{label} sign is undecidable") from exc
    try:
        nonnegative = value >= 0
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise EndpointInputError(f"{label} sign is undecidable") from exc
    if not isinstance(nonnegative, bool) or not nonnegative:
        raise EndpointInputError(f"{label} must be nonnegative")


def _resolve_alias(
    value: Any,
    aliases: dict[str, Any],
    canonical: str,
    names: tuple[str, ...],
) -> Any:
    """Resolve one canonical input and optional spelling aliases."""

    chosen = value
    for name in names:
        if name not in aliases:
            continue
        candidate = aliases.pop(name)
        if chosen is not _MISSING and not _same_value(chosen, candidate):
            raise EndpointInputError(f"{canonical} and {name} disagree")
        chosen = candidate
    return chosen


def _resolve_expm1(expm1_d: Any, expm1: Any, aliases: dict[str, Any]) -> Any:
    chosen = expm1_d
    if expm1 is not _MISSING:
        if chosen is not _MISSING and not _same_value(chosen, expm1):
            raise EndpointInputError("expm1_d and expm1 disagree")
        chosen = expm1
    for name in ("expm1_value", "expm1_backend_value", "expm1_d_value"):
        if name not in aliases:
            continue
        candidate = aliases.pop(name)
        if chosen is not _MISSING and not _same_value(chosen, candidate):
            raise EndpointInputError(f"expm1_d and {name} disagree")
        chosen = candidate
    if chosen is _MISSING:
        raise EndpointInputError(
            "expm1_d is required; provide a caller-evaluated high-precision or interval expm1(d)"
        )
    return chosen


def _bound_check(expm1_d: Any, bound: Any) -> bool:
    """Check a supplied absolute safety bound when the backend permits it."""

    if _interval_like(expm1_d) or _interval_like(bound):
        lower = _call_bound(expm1_d, ("lower", "inf", "lo", "a"))
        upper = _call_bound(expm1_d, ("upper", "sup", "hi", "b"))
        if lower is None or upper is None:
            return False
        try:
            magnitude = max(abs(lower), abs(upper))
            checked = magnitude <= bound
        except (TypeError, ValueError, ArithmeticError):
            return False
        return isinstance(checked, bool) and checked
    try:
        checked = abs(expm1_d) <= bound
    except (TypeError, ValueError, ArithmeticError):
        return False
    if not isinstance(checked, bool):
        return False
    if not checked:
        raise EndpointInputError("supplied expm1_d exceeds expm1_bound")
    return True


@dataclass(frozen=True)
class SignedEndpointResult:
    """Direct endpoint value and ordinary ``Z1`` algebra.

    ``fixed_coordinate_partial`` is the complete fixed-``xi`` derivative
    ``R_F_Z_at_xi * H + R_F * H_Z_at_xi``.  ``moving_endpoint_correction``
    and ``moving_endpoint_total`` are ``None`` unless the caller supplies
    ``xi_e_Z`` and all four ``xi`` derivatives.  ``certified`` is always
    false; ``interval_rounding_certified`` records only caller input.
    """

    R_F: Any
    R_F_Z_at_xi: Any
    Q_y_F: Any
    Q_y_F_Z_at_xi: Any
    dQ_y: Any
    dQ_y_Z_at_xi: Any
    d: Any
    d_Z_at_xi: Any
    log_R_F: Any | None
    expm1_d: Any
    exp_d: Any
    frame_plus_signed_Q_y: Any
    frame_plus_signed_Q_y_Z: Any
    H: Any
    H_Z_at_xi: Any
    R_F_xi: Any | None
    d_xi: Any | None
    Q_y_F_xi: Any | None
    dQ_y_xi: Any | None
    xi_e_Z: Any | None
    H_xi: Any | None
    endpoint: Any
    endpoint_Z: Any
    fixed_coordinate_partial: Any
    fixed_coordinate_radius_term: Any
    endpoint_xi: Any | None
    moving_endpoint_correction: Any | None
    moving_endpoint_total: Any | None
    moving_data_complete: bool
    interval_inputs_present: bool
    expm1_bound_supplied: bool
    expm1_bound_checked: bool
    interval_rounding_certified: bool
    certified: bool
    source_identity_proved: bool
    source_values_bound: bool
    authentication_proved: bool
    status: str

    @property
    def E(self) -> Any:
        return self.endpoint

    @property
    def E_Z(self) -> Any:
        return self.endpoint_Z

    @property
    def R_F_Z(self) -> Any:
        return self.R_F_Z_at_xi

    @property
    def Q_y_F_Z(self) -> Any:
        return self.Q_y_F_Z_at_xi

    @property
    def dQ_y_Z(self) -> Any:
        return self.dQ_y_Z_at_xi

    @property
    def d_Z(self) -> Any:
        return self.d_Z_at_xi

    @property
    def H_Z(self) -> Any:
        return self.H_Z_at_xi

    @property
    def endpoint_Z_fixed_coordinate(self) -> Any:
        return self.endpoint_Z

    @property
    def endpoint_Z_total(self) -> Any | None:
        return self.moving_endpoint_total

    @property
    def caller_interval_rounding_asserted(self) -> bool:
        return self.interval_rounding_certified

    @property
    def fixed_coordinate_Z1(self) -> Any:
        return self.fixed_coordinate_partial

    @property
    def moving_endpoint_Z1(self) -> Any:
        return self.moving_endpoint_total

    @property
    def endpoint_fixed_coordinate_partial(self) -> Any:
        return self.fixed_coordinate_partial

    @property
    def endpoint_moving_total(self) -> Any:
        return self.moving_endpoint_total

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "E": "endpoint",
            "E_Z": "endpoint_Z",
            "fixed_coordinate_Z1": "fixed_coordinate_partial",
            "endpoint_Z_fixed_coordinate": "endpoint_Z",
            "endpoint_Z_total": "moving_endpoint_total",
            "moving_endpoint_Z1": "moving_endpoint_total",
        }
        return getattr(self, aliases.get(key, key))

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except AttributeError:
            return default

    def keys(self) -> tuple[str, ...]:
        return tuple(field.name for field in fields(self))

    def items(self) -> tuple[tuple[str, Any], ...]:
        return tuple((key, self[key]) for key in self.keys())

    def as_dict(self) -> dict[str, Any]:
        """Return fields without converting any backend values."""

        return {field.name: getattr(self, field.name) for field in fields(self)}


def evaluate_signed_endpoint(
    R_F: Any = _MISSING,
    R_F_Z_at_xi: Any = _MISSING,
    Q_y_F: Any = _MISSING,
    Q_y_F_Z_at_xi: Any = _MISSING,
    dQ_y: Any = _MISSING,
    dQ_y_Z_at_xi: Any = _MISSING,
    d: Any = _MISSING,
    d_Z_at_xi: Any = _MISSING,
    *,
    expm1_d: Any = _MISSING,
    expm1: Any = _MISSING,
    expm1_bound: Any | None = None,
    log_R_F: Any | None = None,
    xi_e_Z: Any = _MISSING,
    R_F_xi: Any = _MISSING,
    d_xi: Any = _MISSING,
    Q_y_F_xi: Any = _MISSING,
    dQ_y_xi: Any = _MISSING,
    interval_rounding_certified: bool = False,
    **aliases: Any,
) -> SignedEndpointResult:
    """Evaluate endpoint value and fixed-``xi`` ordinary ``Z1`` algebra.

    Positional order is ``R_F, R_F_Z_at_xi, Q_y_F, Q_y_F_Z_at_xi,
    dQ_y, dQ_y_Z_at_xi, d, d_Z_at_xi``.  The four ``*_at_xi`` values are
    fixed-coordinate partials, not moving-endpoint totals.  Common compact
    spellings (``Qy_F``, ``dQy``, ``dZ``, and ``RF``) are accepted as aliases.

    Optional moving data consists of ``xi_e_Z`` and all four of
    ``R_F_xi``, ``d_xi``, ``Q_y_F_xi`` and ``dQ_y_xi``.  If any member is
    supplied, all five are required; otherwise ``moving_endpoint_total`` is
    ``None``.  This deliberately does not infer a zero endpoint speed.

    ``expm1_d`` (or ``expm1``) is mandatory and must be evaluated by the
    caller's high-precision/interval backend.  ``exp_d`` in the result is
    only the algebraic ``1 + expm1_d``; this routine never calls ``exp``.
    ``expm1_bound`` is an optional caller-supplied absolute safety bound.
    ``log_R_F`` is an optional source-owned scale label and is retained
    formally, never exponentiated.  ``H`` and ``H_Z_at_xi`` are dimensionless
    brackets.  ``interval_rounding_certified`` is recorded only as a caller
    assertion; this result's ``certified`` field is always false.
    """

    R_F = _resolve_alias(R_F, aliases, "R_F", ("RF", "radius_f", "R_f"))
    R_F_Z_at_xi = _resolve_alias(
        R_F_Z_at_xi,
        aliases,
        "R_F_Z_at_xi",
        ("R_F_Z", "RF_Z", "R_F_Z_xi", "radius_f_Z", "R_f_Z"),
    )
    Q_y_F = _resolve_alias(
        Q_y_F, aliases, "Q_y_F", ("Qy_F", "Q_yF", "QyF", "frame_Q_y")
    )
    Q_y_F_Z_at_xi = _resolve_alias(
        Q_y_F_Z_at_xi,
        aliases,
        "Q_y_F_Z_at_xi",
        (
            "Q_y_F_Z",
            "Qy_F_Z",
            "Q_yF_Z",
            "QyF_Z",
            "Q_yFZ",
            "Q_y_F_Z_xi",
            "frame_Q_y_Z",
        ),
    )
    dQ_y = _resolve_alias(
        dQ_y, aliases, "dQ_y", ("dQy", "DeltaQ_y", "delta_Q_y")
    )
    dQ_y_Z_at_xi = _resolve_alias(
        dQ_y_Z_at_xi,
        aliases,
        "dQ_y_Z_at_xi",
        (
            "dQ_y_Z",
            "dQy_Z",
            "dQ_yZ",
            "DeltaQ_y_Z",
            "delta_Q_y_Z",
            "dQ_y_Z_xi",
        ),
    )
    d = _resolve_alias(d, aliases, "d", ("log_radius_ratio",))
    d_Z_at_xi = _resolve_alias(
        d_Z_at_xi,
        aliases,
        "d_Z_at_xi",
        ("d_Z", "dZ", "log_radius_ratio_Z", "d_Z_xi"),
    )
    expm1_d = _resolve_expm1(expm1_d, expm1, aliases)
    if "safe_expm1_bound" in aliases:
        candidate = aliases.pop("safe_expm1_bound")
        if expm1_bound is not None and not _same_value(expm1_bound, candidate):
            raise EndpointInputError("expm1_bound and safe_expm1_bound disagree")
        expm1_bound = candidate
    if "log_radius_f" in aliases:
        candidate = aliases.pop("log_radius_f")
        if log_R_F is not None and not _same_value(log_R_F, candidate):
            raise EndpointInputError("log_R_F and log_radius_f disagree")
        log_R_F = candidate
    xi_e_Z = _resolve_alias(
        xi_e_Z,
        aliases,
        "xi_e_Z",
        ("xi_Z", "endpoint_xi_Z", "xi_endpoint_Z"),
    )
    R_F_xi = _resolve_alias(
        R_F_xi, aliases, "R_F_xi", ("RF_xi", "radius_f_xi")
    )
    d_xi = _resolve_alias(d_xi, aliases, "d_xi", ("dXi", "log_radius_ratio_xi"))
    Q_y_F_xi = _resolve_alias(
        Q_y_F_xi,
        aliases,
        "Q_y_F_xi",
        ("Qy_F_xi", "Q_yF_xi", "QyFxi"),
    )
    dQ_y_xi = _resolve_alias(
        dQ_y_xi,
        aliases,
        "dQ_y_xi",
        ("dQy_xi", "dQ_yxi", "DeltaQ_y_xi"),
    )
    if aliases:
        unknown = ", ".join(sorted(aliases))
        raise EndpointInputError(f"unknown endpoint input(s): {unknown}")
    if not isinstance(interval_rounding_certified, bool):
        raise EndpointInputError("interval_rounding_certified must be boolean")

    values = {
        "R_F": R_F,
        "R_F_Z_at_xi": R_F_Z_at_xi,
        "Q_y_F": Q_y_F,
        "Q_y_F_Z_at_xi": Q_y_F_Z_at_xi,
        "dQ_y": dQ_y,
        "dQ_y_Z_at_xi": dQ_y_Z_at_xi,
        "d": d,
        "d_Z_at_xi": d_Z_at_xi,
        "expm1_d": expm1_d,
    }
    if log_R_F is not None:
        values["log_R_F"] = log_R_F
    interval_inputs_present = False
    for label, value in values.items():
        interval_inputs_present = _validate_value(value, label) or interval_inputs_present
    _require_positive(R_F, "R_F")
    if _exact_zero(d) and not _exact_zero(expm1_d):
        raise EndpointInputError("d = 0 requires caller expm1_d = 0")

    moving_values = {
        "xi_e_Z": xi_e_Z,
        "R_F_xi": R_F_xi,
        "d_xi": d_xi,
        "Q_y_F_xi": Q_y_F_xi,
        "dQ_y_xi": dQ_y_xi,
    }
    moving_supplied = tuple(
        name
        for name, value in moving_values.items()
        if value is not _MISSING and value is not None
    )
    moving_data_complete = len(moving_supplied) == len(moving_values)
    if moving_supplied and not moving_data_complete:
        missing = ", ".join(
            name
            for name, value in moving_values.items()
            if value is _MISSING or value is None
        )
        raise EndpointInputError(
            "moving endpoint data is incomplete; missing " + missing
        )
    if moving_data_complete:
        for label, value in moving_values.items():
            interval_inputs_present = _validate_value(value, label) or interval_inputs_present

    bound_supplied = expm1_bound is not None
    bound_checked = False
    if bound_supplied:
        interval_inputs_present = (
            _validate_value(expm1_bound, "expm1_bound")
            or interval_inputs_present
        )
        _require_nonnegative(expm1_bound, "expm1_bound")
        bound_checked = _bound_check(expm1_d, expm1_bound)

    try:
        # This is deliberately the stable signed form.  It never forms a
        # large physical P value and never computes P - F.
        exp_d = 1 + expm1_d
        frame_plus_signed_Q_y = Q_y_F + dQ_y
        frame_plus_signed_Q_y_Z = Q_y_F_Z_at_xi + dQ_y_Z_at_xi
        H = dQ_y + expm1_d * frame_plus_signed_Q_y
        H_Z_at_xi = (
            dQ_y_Z_at_xi
            + exp_d * d_Z_at_xi * frame_plus_signed_Q_y
            + expm1_d * frame_plus_signed_Q_y_Z
        )
        endpoint = R_F * H
        fixed_coordinate_radius_term = R_F_Z_at_xi * H
        fixed_coordinate_partial = fixed_coordinate_radius_term + R_F * H_Z_at_xi
        H_xi = endpoint_xi = moving_endpoint_correction = moving_endpoint_total = None
        if moving_data_complete:
            frame_plus_signed_Q_y_xi = Q_y_F_xi + dQ_y_xi
            H_xi = (
                dQ_y_xi
                + exp_d * d_xi * frame_plus_signed_Q_y
                + expm1_d * frame_plus_signed_Q_y_xi
            )
            endpoint_xi = R_F_xi * H + R_F * H_xi
            moving_endpoint_correction = xi_e_Z * endpoint_xi
            moving_endpoint_total = fixed_coordinate_partial + moving_endpoint_correction
    except (TypeError, ValueError, ArithmeticError, OverflowError) as exc:
        raise EndpointInputError(
            "endpoint arithmetic failed in the caller's exact/interval domain"
        ) from exc

    for label, value in {
        "exp_d": exp_d,
        "H": H,
        "H_Z_at_xi": H_Z_at_xi,
        "endpoint": endpoint,
        "fixed_coordinate_partial": fixed_coordinate_partial,
    }.items():
        _validate_value(value, label)
    for label, value in {
        "H_xi": H_xi,
        "endpoint_xi": endpoint_xi,
        "moving_endpoint_correction": moving_endpoint_correction,
        "moving_endpoint_total": moving_endpoint_total,
    }.items():
        if value is not None:
            _validate_value(value, label)

    return SignedEndpointResult(
        R_F=R_F,
        R_F_Z_at_xi=R_F_Z_at_xi,
        Q_y_F=Q_y_F,
        Q_y_F_Z_at_xi=Q_y_F_Z_at_xi,
        dQ_y=dQ_y,
        dQ_y_Z_at_xi=dQ_y_Z_at_xi,
        d=d,
        d_Z_at_xi=d_Z_at_xi,
        log_R_F=log_R_F,
        expm1_d=expm1_d,
        exp_d=exp_d,
        frame_plus_signed_Q_y=frame_plus_signed_Q_y,
        frame_plus_signed_Q_y_Z=frame_plus_signed_Q_y_Z,
        H=H,
        H_Z_at_xi=H_Z_at_xi,
        R_F_xi=None if not moving_data_complete else R_F_xi,
        d_xi=None if not moving_data_complete else d_xi,
        Q_y_F_xi=None if not moving_data_complete else Q_y_F_xi,
        dQ_y_xi=None if not moving_data_complete else dQ_y_xi,
        xi_e_Z=None if not moving_data_complete else xi_e_Z,
        H_xi=H_xi,
        endpoint=endpoint,
        endpoint_Z=fixed_coordinate_partial,
        fixed_coordinate_partial=fixed_coordinate_partial,
        fixed_coordinate_radius_term=fixed_coordinate_radius_term,
        endpoint_xi=endpoint_xi,
        moving_endpoint_correction=moving_endpoint_correction,
        moving_endpoint_total=moving_endpoint_total,
        moving_data_complete=moving_data_complete,
        interval_inputs_present=interval_inputs_present,
        expm1_bound_supplied=bound_supplied,
        expm1_bound_checked=bound_checked,
        interval_rounding_certified=interval_rounding_certified,
        certified=False,
        source_identity_proved=False,
        source_values_bound=False,
        authentication_proved=False,
        status="numerical_returned",
    )


def evaluate_signed_endpoint_from_mapping(
    frame: Mapping[str, Any],
    signed: Mapping[str, Any],
    *,
    d: Any = _MISSING,
    d_Z_at_xi: Any = _MISSING,
    expm1_d: Any = _MISSING,
    expm1: Any = _MISSING,
    expm1_bound: Any | None = None,
    log_R_F: Any | None = None,
    xi_e_Z: Any = _MISSING,
    R_F_xi: Any = _MISSING,
    d_xi: Any = _MISSING,
    Q_y_F_xi: Any = _MISSING,
    dQ_y_xi: Any = _MISSING,
    interval_rounding_certified: bool = False,
    **aliases: Any,
) -> SignedEndpointResult:
    """Adapt separate frame and signed-difference mappings.

    Keeping the mappings separate prevents accidental reconstruction of ``P``
    and ``F``.  Every required ``y`` derivative must be present.
    """

    if not isinstance(frame, Mapping) or not isinstance(signed, Mapping):
        raise EndpointInputError("frame and signed must be mappings")

    def pick(source: Mapping[str, Any], canonical: str, names: tuple[str, ...]) -> Any:
        found = [name for name in (canonical,) + names if name in source]
        if not found:
            raise EndpointInputError(f"mapping is missing required {canonical}")
        value = source[found[0]]
        for name in found[1:]:
            if not _same_value(value, source[name]):
                raise EndpointInputError(f"mapping aliases for {canonical} disagree")
        return value

    frame_names = {
        "R_F": ("RF", "radius_f", "R_f"),
        "R_F_Z_at_xi": ("R_F_Z", "RF_Z", "radius_f_Z", "R_f_Z"),
        "Q_y_F": ("Qy_F", "Q_yF", "QyF", "frame_Q_y"),
        "Q_y_F_Z_at_xi": (
            "Q_y_F_Z",
            "Qy_F_Z",
            "Q_yF_Z",
            "QyF_Z",
            "frame_Q_y_Z",
            "Q_yFZ",
        ),
    }
    signed_names = {
        "dQ_y": ("dQy", "DeltaQ_y", "delta_Q_y"),
        "dQ_y_Z_at_xi": (
            "dQ_y_Z",
            "dQy_Z",
            "dQ_yZ",
            "DeltaQ_y_Z",
            "delta_Q_y_Z",
        ),
    }
    frame_values = {name: pick(frame, name, names) for name, names in frame_names.items()}
    signed_values = {name: pick(signed, name, names) for name, names in signed_names.items()}
    return evaluate_signed_endpoint(
        **frame_values,
        **signed_values,
        d=d,
        d_Z_at_xi=d_Z_at_xi,
        expm1_d=expm1_d,
        expm1=expm1,
        expm1_bound=expm1_bound,
        log_R_F=log_R_F,
        xi_e_Z=xi_e_Z,
        R_F_xi=R_F_xi,
        d_xi=d_xi,
        Q_y_F_xi=Q_y_F_xi,
        dQ_y_xi=dQ_y_xi,
        interval_rounding_certified=interval_rounding_certified,
        **aliases,
    )


# Friendly names retain the same direct signed-difference contract.  There
# is intentionally no raw P/F subtraction path.
signed_endpoint_and_Z1 = evaluate_signed_endpoint
evaluate_endpoint = evaluate_signed_endpoint
ordinary_z1_endpoint = evaluate_signed_endpoint
compute_signed_endpoint = evaluate_signed_endpoint
pulse_signed_endpoint = evaluate_signed_endpoint
evaluate_endpoint_from_mapping = evaluate_signed_endpoint_from_mapping


HAND_CHECK_EXAMPLES = (
    {
        "name": "d_zero_fixed_coordinate_only",
        "inputs": {
            "R_F": 2,
            "R_F_Z_at_xi": 3,
            "Q_y_F": 5,
            "Q_y_F_Z_at_xi": 7,
            "dQ_y": 11,
            "dQ_y_Z_at_xi": 13,
            "d": 0,
            "d_Z_at_xi": 0,
            "expm1_d": 0,
        },
        "expected": {
            "endpoint": 22,
            "endpoint_Z": 59,
            "fixed_coordinate_partial": 59,
            "moving_endpoint_correction": None,
            "moving_endpoint_total": None,
        },
    },
    {
        "name": "d_zero_nonzero_dZ_chain",
        "inputs": {
            "R_F": 2,
            "R_F_Z_at_xi": 3,
            "Q_y_F": 5,
            "Q_y_F_Z_at_xi": 7,
            "dQ_y": 11,
            "dQ_y_Z_at_xi": 13,
            "d": 0,
            "d_Z_at_xi": 2,
            "expm1_d": 0,
        },
        "expected": {
            "endpoint": 22,
            "endpoint_Z": 123,
            "fixed_coordinate_partial": 123,
            "moving_endpoint_correction": None,
            "moving_endpoint_total": None,
        },
    },
    {
        "name": "complete_moving_xi_chain",
        "inputs": {
            "R_F": 2,
            "R_F_Z_at_xi": 3,
            "Q_y_F": 5,
            "Q_y_F_Z_at_xi": 7,
            "dQ_y": 11,
            "dQ_y_Z_at_xi": 13,
            "d": 0,
            "d_Z_at_xi": 0,
            "expm1_d": 0,
            "xi_e_Z": 2,
            "R_F_xi": 1,
            "d_xi": 0,
            "Q_y_F_xi": 0,
            "dQ_y_xi": 0,
        },
        "expected": {
            "endpoint": 22,
            "endpoint_Z": 59,
            "fixed_coordinate_partial": 59,
            "endpoint_xi": 11,
            "moving_endpoint_correction": 22,
            "moving_endpoint_total": 81,
        },
    },
)


def hand_check_examples() -> tuple[dict[str, Any], ...]:
    """Return independent integer hand checks without evaluating the module."""

    return HAND_CHECK_EXAMPLES


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "EndpointInputError",
    "SignedEndpointResult",
    "evaluate_signed_endpoint",
    "evaluate_signed_endpoint_from_mapping",
    "evaluate_endpoint_from_mapping",
    "signed_endpoint_and_Z1",
    "evaluate_endpoint",
    "ordinary_z1_endpoint",
    "compute_signed_endpoint",
    "pulse_signed_endpoint",
    "HAND_CHECK_EXAMPLES",
    "hand_check_examples",
]




