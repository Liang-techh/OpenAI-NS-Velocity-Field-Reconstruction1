"""Caller-owned signed P--F scaled transport algebra.

This module is deliberately a small algebra adapter.  It does not read a
source implementation, bind a physical radius or chart measure, and does not
certify a source identity.  Every value and every derivative is supplied by
the caller.  In particular, ``t`` is already the caller's value of
``expm1(kappa)``; this module never evaluates an exponential or reconstructs
``t`` from another quantity.

The fixed-chart formulas are

``t2 = t * (2 + t)``

``B = (1 + t2) * dS + t2 * SF``

``B_Z = (1 + t2) * (dS_Z + 2 * kappa_Z * (dS + SF)) + t2 * SF_Z``

``H = (1 + t) * dQy + t * QyF``

``H_Z = (1 + t) * (dQy_Z + kappa_Z * (dQy + QyF)) + t * QyF_Z``

Thus the density-scale mantissas are ``omega * B / 2`` and
``(omega_Z * B + omega * (2 * ell_Z * B + B_Z)) / 2``.  The endpoint-scale
mantissas are ``H`` and ``ell_Z * H + H_Z``.  The omitted exponential scale
factors are part of the caller's physical interpretation, not a claim made by
this module.

Coordinate convention: all fixed-chart quantities are evaluated at the same
caller-owned position coordinate (the module does not choose or relocate that
coordinate).  ``dS`` and ``dQy`` denote the supplied P--F difference terms,
not derivatives that this module infers.  A ``_Z`` suffix means an ordinary
Z derivative at fixed chart position; ``ell_Z`` is the ordinary derivative of
``ell = log(R_F)`` and ``omega_Z`` is the ordinary derivative of the supplied
chart coefficient.  Optional ``_xi`` values are ordinary derivatives with
respect to the caller's position coordinate at fixed Z.  ``xi_endpoint_Z``
is the moving endpoint coordinate derivative.  These meanings are premises,
not source-level conclusions.

The optional moving-endpoint block is all-or-nothing.  If present, it also
forms ``B_xi`` and the fixed-chart density ``xi`` mantissa so that the
``SF_xi`` and ``dS_xi`` inputs are not silently discarded.  The moving density
formula assumes that the supplied chart coefficient ``omega`` has zero
``xi`` derivative; no ``omega_xi`` is available here.  Under that premise its
mantissa is ``omega * (2 * ell_xi * B + B_xi) / 2``, where

``B_xi = (1 + t2) * (dS_xi + 2 * kappa_xi * (dS + SF)) + t2 * SF_xi``

and

``H_xi = (1 + t) * (dQy_xi + kappa_xi * (dQy + QyF)) + t * QyF_xi``.
The endpoint total derivative mantissa is

``(ell_Z * H + H_Z) + xi_endpoint_Z * (ell_xi * H + H_xi)``.

This is a caller-owned formula result only.  ``source_ownership``,
``physical_scale`` and ``certified`` are explicitly false, and no historical
error is declared to have been removed.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from numbers import Real
from typing import Any


_MISSING = object()


class ScaledTransportInputError(ValueError):
    """Raised when a caller-owned scalar or derivative is missing or unsafe."""


def _is_interval_like(value: Any) -> bool:
    """Recognize common caller interval objects without coercing them."""

    value_type = type(value)
    type_name = value_type.__name__.lower()
    module_name = value_type.__module__.lower()
    if "interval" in type_name or "interval" in module_name:
        return True
    if hasattr(value, "_mpi_"):
        return True
    endpoint_pairs = (("lower", "upper"), ("lo", "hi"), ("inf", "sup"), ("a", "b"))
    return any(hasattr(value, lower) and hasattr(value, upper) for lower, upper in endpoint_pairs)


def _require_scalar(value: Any, name: str) -> Any:
    """Require a non-coerced caller scalar, Fraction, Decimal, or interval."""

    if value is _MISSING:
        raise ScaledTransportInputError(f"missing required caller input: {name}")
    if value is None:
        raise ScaledTransportInputError(f"caller input {name} cannot be None")
    if isinstance(value, bool):
        raise ScaledTransportInputError(f"caller input {name} cannot be bool")
    if isinstance(value, (int, float)):
        raise ScaledTransportInputError(
            f"caller input {name} must retain Fraction/interval semantics; plain int/float is rejected"
        )
    if isinstance(value, complex):
        raise ScaledTransportInputError(f"caller input {name} must be real or interval-valued")
    if isinstance(value, Decimal) and not value.is_finite():
        raise ScaledTransportInputError(f"caller input {name} must be finite")
    if isinstance(value, Fraction):
        return value
    if _is_interval_like(value):
        return value
    if isinstance(value, Real):
        raise ScaledTransportInputError(
            f"caller input {name} uses an unsupported Real type; pass Fraction or an interval"
        )
    if isinstance(value, (str, bytes, bytearray)):
        raise ScaledTransportInputError(f"caller input {name} must be numeric, not text")
    # Decimal and caller-owned symbolic/interval scalar types are intentionally
    # returned unchanged.  No conversion to float or int is performed.
    return value


def _arithmetic_error(name: str, error: Exception) -> ScaledTransportInputError:
    return ScaledTransportInputError(f"caller-owned arithmetic failed while forming {name}: {error}")


@dataclass(frozen=True)
class SignedScaledTransportZ1Result:
    """Formula outputs with explicit ownership and certification metadata."""

    t: Any
    t2: Any
    B: Any
    B_Z: Any
    H: Any
    H_Z: Any
    density_mantissa: Any
    density_Z_mantissa: Any
    endpoint_mantissa: Any
    endpoint_Z_mantissa: Any
    moving_inputs_supplied: bool
    B_xi: Any = None
    H_xi: Any = None
    density_xi_mantissa: Any = None
    density_Z_total_mantissa: Any = None
    endpoint_xi_mantissa: Any = None
    endpoint_Z_total_mantissa: Any = None
    # These are intentionally data fields, not inferred claims.
    source_ownership: bool = False
    physical_scale: bool = False
    certified: bool = False
    source_ownership_proved: bool = False
    physical_scale_proved: bool = False
    certification: bool = False
    source_identity_proved: bool = False
    physical_source_bound: bool = False
    status: str = "caller_formula_only"

    @property
    def density(self) -> Any:
        """Alias for the density-scale mantissa."""

        return self.density_mantissa

    @property
    def density_Z(self) -> Any:
        """Alias for the fixed-chart physical Z1 density mantissa."""

        return self.density_Z_mantissa

    @property
    def endpoint(self) -> Any:
        """Alias for the endpoint-scale mantissa."""

        return self.endpoint_mantissa

    @property
    def endpoint_Z(self) -> Any:
        """Alias for the fixed-chart endpoint Z1 mantissa."""

        return self.endpoint_Z_mantissa

    @property
    def endpoint_Z_total(self) -> Any:
        """Alias for the optional moving-endpoint total derivative."""

        return self.endpoint_Z_total_mantissa


def evaluate_signed_scaled_transport_z1(
    *,
    omega: Any = _MISSING,
    omega_Z: Any = _MISSING,
    ell_Z: Any = _MISSING,
    kappa_Z: Any = _MISSING,
    t: Any = _MISSING,
    SF: Any = _MISSING,
    SF_Z: Any = _MISSING,
    dS: Any = _MISSING,
    dS_Z: Any = _MISSING,
    QyF: Any = _MISSING,
    QyF_Z: Any = _MISSING,
    dQy: Any = _MISSING,
    dQy_Z: Any = _MISSING,
    ell_xi: Any = _MISSING,
    kappa_xi: Any = _MISSING,
    SF_xi: Any = _MISSING,
    dS_xi: Any = _MISSING,
    QyF_xi: Any = _MISSING,
    dQy_xi: Any = _MISSING,
    xi_endpoint_Z: Any = _MISSING,
) -> SignedScaledTransportZ1Result:
    """Evaluate fixed-chart signed P--F density and endpoint Z1 algebra.

    All thirteen fixed-chart inputs are mandatory.  The seven optional
    moving-endpoint inputs are an all-or-nothing block.  ``t`` must already be
    the caller's ``expm1(kappa)`` value; no relation can be checked because
    ``kappa`` itself is deliberately not reconstructed here.
    """

    required = {
        "omega": omega,
        "omega_Z": omega_Z,
        "ell_Z": ell_Z,
        "kappa_Z": kappa_Z,
        "t": t,
        "SF": SF,
        "SF_Z": SF_Z,
        "dS": dS,
        "dS_Z": dS_Z,
        "QyF": QyF,
        "QyF_Z": QyF_Z,
        "dQy": dQy,
        "dQy_Z": dQy_Z,
    }
    values = {name: _require_scalar(value, name) for name, value in required.items()}

    moving = {
        "ell_xi": ell_xi,
        "kappa_xi": kappa_xi,
        "SF_xi": SF_xi,
        "dS_xi": dS_xi,
        "QyF_xi": QyF_xi,
        "dQy_xi": dQy_xi,
        "xi_endpoint_Z": xi_endpoint_Z,
    }
    supplied = {name: value is not _MISSING for name, value in moving.items()}
    if any(supplied.values()) and not all(supplied.values()):
        missing = ", ".join(name for name, present in supplied.items() if not present)
        raise ScaledTransportInputError(
            "moving-endpoint inputs are all-or-nothing; missing values cannot be filled with zero: "
            + missing
        )
    moving_values = {
        name: _require_scalar(value, name) for name, value in moving.items()
    } if all(supplied.values()) else {}

    try:
        one_plus_t = 1 + values["t"]
        values_t2 = values["t"] * (2 + values["t"])
        one_plus_t2 = 1 + values_t2
        B = one_plus_t2 * values["dS"] + values_t2 * values["SF"]
        B_Z = one_plus_t2 * (
            values["dS_Z"] + 2 * values["kappa_Z"] * (values["dS"] + values["SF"])
        ) + values_t2 * values["SF_Z"]
        H = one_plus_t * values["dQy"] + values["t"] * values["QyF"]
        H_Z = one_plus_t * (
            values["dQy_Z"] + values["kappa_Z"] * (values["dQy"] + values["QyF"])
        ) + values["t"] * values["QyF_Z"]
        density_mantissa = values["omega"] * B / 2
        density_Z_mantissa = (
            values["omega_Z"] * B
            + values["omega"] * (2 * values["ell_Z"] * B + B_Z)
        ) / 2
        endpoint_mantissa = H
        endpoint_Z_mantissa = values["ell_Z"] * H + H_Z
    except (ArithmeticError, OverflowError, TypeError, ValueError) as error:
        raise _arithmetic_error("fixed-chart transport terms", error) from error

    if not moving_values:
        return SignedScaledTransportZ1Result(
            t=values["t"],
            t2=values_t2,
            B=B,
            B_Z=B_Z,
            H=H,
            H_Z=H_Z,
            density_mantissa=density_mantissa,
            density_Z_mantissa=density_Z_mantissa,
            endpoint_mantissa=endpoint_mantissa,
            endpoint_Z_mantissa=endpoint_Z_mantissa,
            moving_inputs_supplied=False,
        )

    try:
        B_xi = (1 + values_t2) * (
            moving_values["dS_xi"]
            + 2 * moving_values["kappa_xi"] * (values["dS"] + values["SF"])
        ) + values_t2 * moving_values["SF_xi"]
        H_xi = (1 + values["t"]) * (
            moving_values["dQy_xi"]
            + moving_values["kappa_xi"] * (values["dQy"] + values["QyF"])
        ) + values["t"] * moving_values["QyF_xi"]
        # The chart coefficient omega is caller-declared fixed in xi because
        # no omega_xi input exists.  This is a diagnostic, not a source claim.
        density_xi_mantissa = values["omega"] * (
            2 * moving_values["ell_xi"] * B + B_xi
        ) / 2
        endpoint_xi_mantissa = moving_values["ell_xi"] * H + H_xi
        density_Z_total_mantissa = (
            density_Z_mantissa
            + moving_values["xi_endpoint_Z"] * density_xi_mantissa
        )
        endpoint_Z_total_mantissa = (
            endpoint_Z_mantissa
            + moving_values["xi_endpoint_Z"] * endpoint_xi_mantissa
        )
    except (ArithmeticError, OverflowError, TypeError, ValueError) as error:
        raise _arithmetic_error("moving-endpoint transport terms", error) from error

    return SignedScaledTransportZ1Result(
        t=values["t"],
        t2=values_t2,
        B=B,
        B_Z=B_Z,
        H=H,
        H_Z=H_Z,
        density_mantissa=density_mantissa,
        density_Z_mantissa=density_Z_mantissa,
        endpoint_mantissa=endpoint_mantissa,
        endpoint_Z_mantissa=endpoint_Z_mantissa,
        moving_inputs_supplied=True,
        B_xi=B_xi,
        H_xi=H_xi,
        density_xi_mantissa=density_xi_mantissa,
        density_Z_total_mantissa=density_Z_total_mantissa,
        endpoint_xi_mantissa=endpoint_xi_mantissa,
        endpoint_Z_total_mantissa=endpoint_Z_total_mantissa,
    )


# Descriptive aliases keep the standalone formula easy to locate without
# introducing another implementation or importing any existing module.
signed_scaled_transport_z1 = evaluate_signed_scaled_transport_z1
compute_signed_scaled_transport_z1 = evaluate_signed_scaled_transport_z1


__all__ = [
    "ScaledTransportInputError",
    "SignedScaledTransportZ1Result",
    "evaluate_signed_scaled_transport_z1",
    "signed_scaled_transport_z1",
    "compute_signed_scaled_transport_z1",
]
