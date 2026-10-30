"""Source-shaped correlated difference for the complete Gamma heat term.

This leaf deliberately stops one layer above the historical Gamma evaluator.
It carries the exact algebra of

    S * square,   square = D * (2 * pre - a * S * D),

and of the complete source energy integral, but it never imports an older
Gamma kernel, quadrature routine, provider, or source graph.  ``K`` is a
caller-owned total correlation and is recorded as ``A*K0``; ``K0`` is a
reference to the previously published source-normalized correlation, not an
implementation dependency here.

For a common ``a``, ``pre``, and ``A`` at P and F, put

    S_F = SF,  DeltaS = ds,  D_F = DF,  D_P = DF + ds*K,
    X = SF*DF,  Y0 = DF + SF*K.

The returned coefficients are kept separate:

    G0 = 2*(pre-a*X)*Y0
    G1 = 2*(pre-a*X)*K-a*Y0**2
    G2 = -2*a*Y0*K
    G3 = -a*K**2

and the exact difference is ``ds*(G0 + ds*G1 + ds**2*G2 + ds**3*G3)``.
No source value, double-exponential scale, integral, or physical error is
materialized by this module.  The point interfaces accept caller-supplied
ordinary arithmetic, ``Fraction``, SymPy, or mpmath-like values and retain
the scale factors instead of silently summing or cancelling small terms.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "gamma_weighted_energy_delta_v1"
COMPONENT = SCHEMA
PHASE = "gamma_weighted_energy_delta_pending"

SOURCE_VALUE_FLAGS = {
    "source_values_computed": False,
    "source_value_computation": False,
    "historical_source_values_computed": False,
    "actual_source_value_computed": False,
    "certified": False,
}

SOURCE_ASSUMPTIONS = {
    "source_function": (
        "original_n1_Gamma_future_functions.py lines 121-177; "
        "source definition is retained as references only"
    ),
    "collar": "t in [0,3], A=sigma(t)*C(t), pre=1-epsilon*W(t)",
    "outer": "t in [3,infinity), A=1, pre=1",
    "D": "D=A*Dhat",
    "square": "square=D*(2*pre-a*S*D)",
    "energy_hat": "Integral_0^infinity exp(-delta*t)*square dt",
    "infinity": "mathematical positive infinity; full improper integral retained",
    "K": "K=A*K0, with K0 a caller-owned source-normalized correlation reference",
    "S": "original source S argument; no source value is reconstructed here",
    "complete_future_roles": (
        "complete archive direct function_bindings.energy.complete/weighted and "
        "function_bindings.pieces.Nrel/tail_multiplier/Gamma_heat_difference"
    ),
    "gamma_piece_roles": (
        "Gamma archive direct function_bindings.pieces.collar.energy, "
        "exterior.energy, actual_scaled_deficit, preheat_bracket, and a"
    ),
    "P_F": "common a, pre, and A are a gate, not a conclusion of this leaf",
}

SOURCE_POINTER_NAMES = (
    "source_S",
    "Dhat",
    "collar_energy",
    "outer_energy",
    "actual_scaled_deficit",
    "preheat_bracket",
    "a",
)

COMPLETE_FUTURE_POINTER_NAMES = (
    "energy_complete",
    "energy_weighted",
    "Nrel",
    "tail_multiplier",
    "Gamma_heat_difference",
)


class GammaWeightedEnergyDeltaInputError(ValueError):
    """Raised when a formal/point input or source manifest is malformed."""


class GammaWeightedEnergyDeltaSourceError(RuntimeError):
    """Raised when authenticated source metadata cannot be bound safely."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GammaWeightedEnergyDeltaSourceError(message)


def _json_default(value: Any) -> Any:
    """Serialize exact and arbitrary-precision point values without rounding."""

    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "_mpf_"):
        return str(value)
    # SymPy Basic objects intentionally remain textual in receipts.  Importing
    # SymPy at module import would make the formal leaf depend on the optional
    # symbolic extra, so use its stable textual representation instead.
    if type(value).__module__.split(".", 1)[0] == "sympy":
        return str(value)
    if hasattr(value, "as_dict") and callable(value.as_dict):
        return value.as_dict()
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _display(value: Any) -> Any:
    """Make a point value suitable for a JSON-like result packet."""

    try:
        json.dumps(value, allow_nan=False)
    except (TypeError, ValueError):
        return _json_default(value)
    return value


def _allowed_source_kind(source_kind: str) -> str:
    aliases = {
        "archived_exact": "archived_exact_source_function",
        "archived": "archived_exact_source_function",
        "synthetic": "synthetic_point",
        "rounded": "uncertified_rounding",
    }
    source_kind = aliases.get(source_kind, source_kind)
    allowed = {
        "archived_exact_source_function",
        "synthetic_point",
        "uncertified_rounding",
        "formal_symbolic",
    }
    if source_kind not in allowed:
        raise GammaWeightedEnergyDeltaInputError(
            "source_kind must be archived_exact_source_function, "
            "synthetic_point, uncertified_rounding, or formal_symbolic"
        )
    return source_kind


def _materialization_flags(
    *,
    source_kind: str,
    rounded: bool,
    scaled_terms_materialized: bool,
    total_sum_materialized: bool,
) -> dict[str, Any]:
    source_kind = _allowed_source_kind(source_kind)
    if rounded:
        source_kind = "uncertified_rounding"
    return {
        "source_kind": source_kind,
        "source_kind_is_caller_declared": True,
        "caller_declared_archived_exact_source_function": source_kind == "archived_exact_source_function",
        "archived_exact_source_function": False,
        "source_identity_proved": False,
        "synthetic_point": source_kind == "synthetic_point",
        "uncertified_rounding": source_kind == "uncertified_rounding",
        "formal_symbolic": source_kind == "formal_symbolic",
        "deltaS_factored": True,
        "SF_retained_separately": True,
        "scaled_terms_materialized": bool(scaled_terms_materialized),
        "scaled_polynomial_sum_materialized": bool(total_sum_materialized),
        "small_terms_structurally_retained": True,
        "floating_point_cancellation_not_excluded": True,
        "double_exponential_materialized": False,
        "source_values_computed": False,
        "source_value_computation": False,
        "historical_source_values_computed": False,
        "certified": False,
    }


def _jet_pair(value: Any, derivative: Any, label: str) -> tuple[Any, Any]:
    """Accept a scalar plus derivative or a small caller-owned total Z jet."""

    if derivative is not None:
        return value, derivative
    if isinstance(value, Z1Jet):
        return value.value, value.Z
    if isinstance(value, Mapping):
        for value_key, derivative_key in (
            ("value", "Z"),
            ("value", "z"),
            ("base", "derivative"),
            ("0", "1"),
        ):
            if value_key in value and derivative_key in value:
                return value[value_key], value[derivative_key]
        raise GammaWeightedEnergyDeltaInputError(
            f"{label} mapping must contain value and Z entries"
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return value[0], value[1]
        raise GammaWeightedEnergyDeltaInputError(
            f"{label} sequence jet must have exactly two entries"
        )
    raise GammaWeightedEnergyDeltaInputError(
        f"{label} requires an explicit ordinary-Z derivative"
    )


@dataclass(frozen=True)
class Z1Jet:
    """A caller-owned total ordinary-Z one-jet."""

    value: Any
    Z: Any

    @property
    def derivative(self) -> Any:
        return self.Z

    def as_dict(self) -> dict[str, Any]:
        return {"value": _display(self.value), "Z": _display(self.Z)}


@dataclass(frozen=True)
class GammaEnergyDeltaCoefficients:
    """Four signed polynomial coefficients for the common-source delta."""

    a: Any
    pre: Any
    SF: Any
    DeltaS: Any
    DF: Any
    K: Any
    X: Any
    Y0: Any
    G0: Any
    G1: Any
    G2: Any
    G3: Any
    materialization: Mapping[str, Any]

    @property
    def ds(self) -> Any:
        return self.DeltaS

    @property
    def G(self) -> tuple[Any, Any, Any, Any]:
        return self.G0, self.G1, self.G2, self.G3

    @property
    def coefficients(self) -> tuple[Any, Any, Any, Any]:
        return self.G

    @property
    def scaled_terms(self) -> tuple[Any, Any, Any, Any]:
        """Return ds*G0, ds^2*G1, ds^3*G2, ds^4*G3 separately."""

        ds = self.DeltaS
        return ds * self.G0, ds**2 * self.G1, ds**3 * self.G2, ds**4 * self.G3

    def materialize_total(self) -> Any:
        """Explicitly sum the four retained scaled terms on caller request."""

        return sum(self.scaled_terms, 0)

    def as_dict(self, *, include_scaled_terms: bool = True) -> dict[str, Any]:
        result: dict[str, Any] = {
            "a": _display(self.a),
            "pre": _display(self.pre),
            "SF": _display(self.SF),
            "DeltaS": _display(self.DeltaS),
            "DF": _display(self.DF),
            "K": _display(self.K),
            "K_definition": "A*K0",
            "X": _display(self.X),
            "Y0": _display(self.Y0),
            "G": {
                "G0": _display(self.G0),
                "G1": _display(self.G1),
                "G2": _display(self.G2),
                "G3": _display(self.G3),
            },
            "polynomial": "DeltaS*(G0 + DeltaS*G1 + DeltaS^2*G2 + DeltaS^3*G3)",
            "materialization": dict(self.materialization),
        }
        if include_scaled_terms:
            scaled = self.scaled_terms
            result["scaled_terms"] = {
                "DeltaS_G0": _display(scaled[0]),
                "DeltaS2_G1": _display(scaled[1]),
                "DeltaS3_G2": _display(scaled[2]),
                "DeltaS4_G3": _display(scaled[3]),
            }
        return result

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "ds": "DeltaS",
            "delta_s": "DeltaS",
            "g0": "G0",
            "g1": "G1",
            "g2": "G2",
            "g3": "G3",
        }
        return getattr(self, aliases.get(key, key))


@dataclass(frozen=True)
class GammaEnergyDeltaZ1:
    """The complete ordinary-Z1 derivative, split into its two ds routes."""

    base: GammaEnergyDeltaCoefficients
    a_Z: Any
    pre_Z: Any
    SF_Z: Any
    DeltaS_Z: Any
    DF_Z: Any
    K_Z: Any
    X_Z: Any
    Y0_Z: Any
    G0_Z: Any
    G1_Z: Any
    G2_Z: Any
    G3_Z: Any
    G: Any
    G_Z: Any
    DeltaS_Z_times_G: Any
    DeltaS_times_G_Z: Any
    total_Z: Any
    materialization: Mapping[str, Any]

    @property
    def ds_Z(self) -> Any:
        return self.DeltaS_Z

    @property
    def coefficients_Z(self) -> tuple[Any, Any, Any, Any]:
        return self.G0_Z, self.G1_Z, self.G2_Z, self.G3_Z

    def as_dict(self) -> dict[str, Any]:
        return {
            "base": self.base.as_dict(include_scaled_terms=False),
            "jets": {
                "a_Z": _display(self.a_Z),
                "pre_Z": _display(self.pre_Z),
                "SF_Z": _display(self.SF_Z),
                "DeltaS_Z": _display(self.DeltaS_Z),
                "DF_Z": _display(self.DF_Z),
                "K_Z": _display(self.K_Z),
                "X_Z": _display(self.X_Z),
                "Y0_Z": _display(self.Y0_Z),
            },
            "G_Z": {
                "G0_Z": _display(self.G0_Z),
                "G1_Z": _display(self.G1_Z),
                "G2_Z": _display(self.G2_Z),
                "G3_Z": _display(self.G3_Z),
            },
            "G_polynomial": _display(self.G),
            "G_polynomial_Z": _display(self.G_Z),
            "retained_routes": {
                "DeltaS_Z_times_G": _display(self.DeltaS_Z_times_G),
                "DeltaS_times_G_Z": _display(self.DeltaS_times_G_Z),
            },
            "total_Z": _display(self.total_Z),
            "formula": (
                "DeltaS_Z*G + DeltaS*G_Z; "
                "G_Z=sum_j DeltaS^j*G_j_Z + "
                "DeltaS_Z*(G1 + 2*DeltaS*G2 + 3*DeltaS^2*G3)"
            ),
            "materialization": dict(self.materialization),
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "ds_Z": "DeltaS_Z",
            "delta_s_Z": "DeltaS_Z",
            "total": "total_Z",
        }
        return getattr(self, aliases.get(key, key))


def point_coefficients(
    a: Any,
    pre: Any,
    SF: Any,
    DeltaS: Any,
    DF: Any,
    K: Any,
    *,
    source_kind: str = "synthetic_point",
    rounded: bool = False,
    materialize_scaled_terms: bool = True,
    materialize_total: bool = False,
) -> GammaEnergyDeltaCoefficients:
    """Build exact/synthetic signed coefficients without any source evaluation.

    Values are intentionally not converted to float or mpmath.  This keeps
    ``Fraction``, SymPy, and caller-owned arbitrary-precision values in their
    native arithmetic domain.  ``DF`` is the full F-side ``D=A*Dhat`` and
    ``K`` is the full correlation ``A*K0``; callers must not pass the old
    ``K1`` as a total ordinary-Z derivative.
    """

    X = SF * DF
    Y0 = DF + SF * K
    common = pre - a * X
    G0 = 2 * common * Y0
    G1 = 2 * common * K - a * Y0**2
    G2 = -2 * a * Y0 * K
    G3 = -a * K**2
    materialization = _materialization_flags(
        source_kind=source_kind,
        rounded=rounded,
        scaled_terms_materialized=materialize_scaled_terms,
        total_sum_materialized=materialize_total,
    )
    return GammaEnergyDeltaCoefficients(
        a,
        pre,
        SF,
        DeltaS,
        DF,
        K,
        X,
        Y0,
        G0,
        G1,
        G2,
        G3,
        materialization,
    )


def ordinary_z1_coefficients(
    a: Any,
    pre: Any,
    SF: Any,
    DeltaS: Any,
    DF: Any,
    K: Any,
    *,
    a_Z: Any = None,
    pre_Z: Any = None,
    SF_Z: Any = None,
    DeltaS_Z: Any = None,
    DF_Z: Any = None,
    K_Z: Any = None,
    source_kind: str = "synthetic_point",
    rounded: bool = False,
) -> GammaEnergyDeltaZ1:
    """Build the complete ordinary-Z1 polynomial jet by Leibniz.

    Each core argument may itself be a ``Z1Jet``, a two-entry sequence, or a
    mapping with ``value``/``Z`` keys.  Explicit derivative keywords take
    precedence.  Missing derivatives are rejected; a zero derivative must
    be supplied explicitly by the caller.
    """

    a0, az = _jet_pair(a, a_Z, "a")
    pre0, prez = _jet_pair(pre, pre_Z, "pre")
    sf0, sfz = _jet_pair(SF, SF_Z, "SF")
    ds0, dsz = _jet_pair(DeltaS, DeltaS_Z, "DeltaS")
    df0, dfz = _jet_pair(DF, DF_Z, "DF")
    k0, kz = _jet_pair(K, K_Z, "K")
    base = point_coefficients(
        a0,
        pre0,
        sf0,
        ds0,
        df0,
        k0,
        source_kind=source_kind,
        rounded=rounded,
        materialize_scaled_terms=True,
    )
    Xz = sfz * df0 + sf0 * dfz
    Y0z = dfz + sfz * k0 + sf0 * kz
    common = pre0 - a0 * base.X
    common_z = prez - az * base.X - a0 * Xz
    G0z = 2 * (common_z * base.Y0 + common * Y0z)
    G1z = (
        2 * (common_z * k0 + common * kz)
        - az * base.Y0**2
        - 2 * a0 * base.Y0 * Y0z
    )
    G2z = -2 * (az * base.Y0 * k0 + a0 * Y0z * k0 + a0 * base.Y0 * kz)
    G3z = -(az * k0**2 + 2 * a0 * k0 * kz)
    ds = ds0
    G = base.G0 + ds * base.G1 + ds**2 * base.G2 + ds**3 * base.G3
    Gz = (
        G0z
        + ds * G1z
        + ds**2 * G2z
        + ds**3 * G3z
        + dsz * (base.G1 + 2 * ds * base.G2 + 3 * ds**2 * base.G3)
    )
    route_ds_z = dsz * G
    route_ds = ds * Gz
    flags = dict(base.materialization)
    flags.update(
        {
            "total_Z_jet": True,
            "SF_Z_retained": True,
            "DeltaS_Z_retained": True,
            "a_Z_retained": True,
            "pre_Z_retained": True,
            "DF_Z_retained": True,
            "K_Z_retained": True,
            "DeltaS_Z_times_G_materialized": True,
            "DeltaS_times_G_Z_materialized": True,
            "K1_used_as_total_K_Z": False,
        }
    )
    return GammaEnergyDeltaZ1(
        base,
        az,
        prez,
        sfz,
        dsz,
        dfz,
        kz,
        Xz,
        Y0z,
        G0z,
        G1z,
        G2z,
        G3z,
        G,
        Gz,
        route_ds_z,
        route_ds,
        route_ds_z + route_ds,
        flags,
    )


def compose_total_K(
    A: Any,
    K0: Any,
    *,
    A_Z: Any,
    K0_Z: Any,
) -> Z1Jet:
    """Compose ``K=A*K0`` and retain both ordinary-Z chain routes."""

    return Z1Jet(A * K0, A_Z * K0 + A * K0_Z)


def general_source_difference(
    *,
    SP: Any,
    SF: Any,
    DP: Any,
    DF: Any,
    aP: Any,
    aF: Any,
    preP: Any,
    preF: Any,
    DeltaS: Any = None,
    DeltaD: Any = None,
) -> dict[str, Any]:
    """Return the general unevaluated identity without common-a/pre gates.

    If ``a``, ``pre``, or ``A`` differ between P and F, ``DeltaD`` must be an
    independently source-bound difference.  This function does not replace it
    with ``DeltaS*K``.  If omitted, ``DP-DF`` is merely the caller's supplied
    algebraic difference, not a source-provider call.
    """

    ds = SP - SF if DeltaS is None else DeltaS
    dd = DP - DF if DeltaD is None else DeltaD
    XF = SF * DF
    dX = SP * dd + ds * DF
    first = dX * (2 * preP - aP * (2 * XF + dX))
    second = 2 * (preP - preF) * XF
    third = -(aP - aF) * XF**2
    value = first + second + third
    return {
        "SP": SP,
        "SF": SF,
        "DP": DP,
        "DF": DF,
        "aP": aP,
        "aF": aF,
        "preP": preP,
        "preF": preF,
        "DeltaS": ds,
        "DeltaD": dd,
        "XF": XF,
        "DeltaX": dX,
        "retained_terms": {
            "DeltaX_bracket": first,
            "Delta_pre_times_XF": second,
            "Delta_a_times_XF_squared": third,
        },
        "value": value,
        "formula": (
            "DeltaX*[2*preP-aP*(2*XF+DeltaX)] + "
            "2*(preP-preF)*XF - (aP-aF)*XF^2"
        ),
        "common_specialization_gate": False,
        "DeltaD_independent_source_difference_required": True,
        **SOURCE_VALUE_FLAGS,
    }


def assemble_heat_difference(
    *,
    DeltaSE: Any,
    SE_F: Any,
    M_P: Any,
    M_F: Any,
    DeltaSE_Z: Any = None,
    SE_F_Z: Any = None,
    M_P_Z: Any = None,
    M_F_Z: Any = None,
) -> dict[str, Any]:
    """Keep the complete ``M=Nrel*tail_multiplier*a`` P/F difference.

    ``DeltaSE`` means the already-integrated difference in ``S*energy_hat``;
    this routine never computes an integral.  All four Z1 routes are returned
    individually so a common-prefactor shortcut cannot erase ``DeltaM``.
    """

    dM = M_P - M_F
    delta_heat_first = M_P * DeltaSE
    delta_heat_second = dM * SE_F
    delta_heat = delta_heat_first + delta_heat_second
    result: dict[str, Any] = {
        "M_P": M_P,
        "M_F": M_F,
        "DeltaM": dM,
        "DeltaSE": DeltaSE,
        "SE_F": SE_F,
        "retained_terms": {
            "M_P_times_DeltaSE": delta_heat_first,
            "DeltaM_times_SE_F": delta_heat_second,
        },
        "Delta_heat": delta_heat,
        "formula": "M_P*Delta(SE) + DeltaM*(SE)_F",
        "common_prefactor_gate": False,
        "other_factor_same_source_gate": False,
        **SOURCE_VALUE_FLAGS,
    }
    if any(item is not None for item in (DeltaSE_Z, SE_F_Z, M_P_Z, M_F_Z)):
        if any(item is None for item in (DeltaSE_Z, SE_F_Z, M_P_Z, M_F_Z)):
            raise GammaWeightedEnergyDeltaInputError(
                "all four heat Z1 inputs are required together"
            )
        dM_Z = M_P_Z - M_F_Z
        z_first = M_P_Z * DeltaSE
        z_second = M_P * DeltaSE_Z
        z_third = dM_Z * SE_F
        z_fourth = dM * SE_F_Z
        result["Z1"] = {
            "M_P_Z_times_DeltaSE": z_first,
            "M_P_times_DeltaSE_Z": z_second,
            "DeltaM_Z_times_SE_F": z_third,
            "DeltaM_times_SE_F_Z": z_fourth,
            "DeltaM_Z": dM_Z,
            "Delta_heat_Z": z_first + z_second + z_third + z_fourth,
            "formula": (
                "M_P_Z*Delta(SE) + M_P*Delta(SE)_Z + "
                "DeltaM_Z*(SE)_F + DeltaM*(SE)_F,Z"
            ),
            "all_four_routes_retained": True,
        }
    else:
        result["Z1"] = None
    return result


def assemble_weighted_future_difference(
    *,
    B_P: Any,
    B_F: Any,
    DeltaHeat: Any,
    Heat_F: Any,
    DeltaHeat_Z: Any = None,
    Heat_F_Z: Any = None,
    B_P_Z: Any = None,
    B_F_Z: Any = None,
) -> dict[str, Any]:
    """Retain the weighted heat contribution's P/F and ordinary-Z routes.

    The source factor is ``B=-mu*exp(-26)/2``.  Directly multiplying a heat
    difference by one B is allowed only after the caller proves the P/F B and
    common-source gate.  This function therefore returns the general product
    identity first.
    """

    dB = B_P - B_F
    first = B_P * DeltaHeat
    second = dB * Heat_F
    result: dict[str, Any] = {
        "B_P": B_P,
        "B_F": B_F,
        "DeltaB": dB,
        "source_factor": "-mu*exp(-26)/2",
        "retained_terms": {
            "B_P_times_DeltaHeat": first,
            "DeltaB_times_Heat_F": second,
        },
        "Delta_weighted_heat_contribution": first + second,
        "formula": "B_P*DeltaHeat + DeltaB*Heat_F",
        "common_mu_gate": False,
        "common_source_factor_gate": False,
        "direct_common_multiply_allowed": False,
        **SOURCE_VALUE_FLAGS,
    }
    if any(item is not None for item in (DeltaHeat_Z, Heat_F_Z, B_P_Z, B_F_Z)):
        if any(item is None for item in (DeltaHeat_Z, Heat_F_Z, B_P_Z, B_F_Z)):
            raise GammaWeightedEnergyDeltaInputError(
                "all four weighted-future Z1 inputs are required together"
            )
        dB_Z = B_P_Z - B_F_Z
        routes = {
            "B_P_Z_times_DeltaHeat": B_P_Z * DeltaHeat,
            "B_P_times_DeltaHeat_Z": B_P * DeltaHeat_Z,
            "DeltaB_Z_times_Heat_F": dB_Z * Heat_F,
            "DeltaB_times_Heat_F_Z": dB * Heat_F_Z,
        }
        result["Z1"] = {
            **routes,
            "DeltaB_Z": dB_Z,
            "Delta_weighted_heat_contribution_Z": sum(routes.values(), 0),
            "formula": (
                "B_P_Z*DeltaHeat + B_P*DeltaHeat_Z + "
                "DeltaB_Z*Heat_F + DeltaB*Heat_F_Z"
            ),
            "all_four_routes_retained": True,
        }
    else:
        result["Z1"] = None
    return result


def time_weight_difference_structure(
    *,
    delta_P: Any = "delta_P",
    delta_F: Any = "delta_F",
    DeltaH: Any = "DeltaH",
    H_F: Any = "H_F",
    delta_Z: Any = "delta_Z",
    DeltaH_Z: Any = "DeltaH_Z",
    H_F_Z: Any = "H_F_Z",
) -> dict[str, Any]:
    """Keep the unevaluated P/F time-weight difference and its Z1 routes.

    This helper is structural: callers may pass symbols, exact values, or
    textual source references.  It never evaluates an exponential integral.
    ``w_i=exp(-delta_i*t)`` and ``Delta w=w_F*expm1(-t*Delta delta)`` are
    retained so a non-common ``delta`` cannot be silently replaced by a
    common weight.
    """

    return {
        "w_P": "exp(-delta_P*t)",
        "w_F": "exp(-delta_F*t)",
        "Delta_delta": "delta_P-delta_F",
        "Delta_w": "exp(-delta_F*t)*expm1(-t*(delta_P-delta_F))",
        "Delta_wH": "w_P*DeltaH + Delta_w*H_F",
        "point_inputs": {
            "delta_P": delta_P,
            "delta_F": delta_F,
            "DeltaH": DeltaH,
            "H_F": H_F,
        },
        "retained_terms": {
            "w_P_times_DeltaH": "exp(-delta_P*t)*DeltaH",
            "Delta_w_times_H_F": "exp(-delta_F*t)*expm1(-t*(delta_P-delta_F))*H_F",
        },
        "common_delta_direct_weight_gate": False,
        "ordinary_Z1": {
            "w_P_Z": "-t*delta_P_Z*exp(-delta_P*t)",
            "w_F_Z": "-t*delta_F_Z*exp(-delta_F*t)",
            "Delta_w_Z": "w_F_Z*expm1(-t*Delta_delta) + w_F*exp(-t*Delta_delta)*(-t*Delta_delta_Z)",
            "Delta_wH_Z": "w_P_Z*DeltaH + w_P*DeltaH_Z + Delta_w_Z*H_F + Delta_w*H_F_Z",
            "point_inputs": {
                "delta_Z": delta_Z,
                "DeltaH_Z": DeltaH_Z,
                "H_F_Z": H_F_Z,
            },
        },
        "source_values_computed": False,
        "certified": False,
    }


def future_prefactor(
    Nrel: Any,
    tail_multiplier: Any,
    a: Any,
    *,
    Nrel_Z: Any = None,
    tail_multiplier_Z: Any = None,
    a_Z: Any = None,
) -> dict[str, Any]:
    """Materialize M and, when requested, its complete ordinary-Z product jet."""

    M = Nrel * tail_multiplier * a
    result: dict[str, Any] = {
        "Nrel": Nrel,
        "tail_multiplier": tail_multiplier,
        "a": a,
        "M": M,
        "formula": "M=Nrel*tail_multiplier*a",
        "Nrel_Z_must_be_retained": True,
        **SOURCE_VALUE_FLAGS,
    }
    if any(item is not None for item in (Nrel_Z, tail_multiplier_Z, a_Z)):
        if any(item is None for item in (Nrel_Z, tail_multiplier_Z, a_Z)):
            raise GammaWeightedEnergyDeltaInputError(
                "all three future prefactor Z1 inputs are required together"
            )
        M_Z = Nrel_Z * tail_multiplier * a
        M_Z += Nrel * tail_multiplier_Z * a
        M_Z += Nrel * tail_multiplier * a_Z
        result["M_Z"] = M_Z
        result["Z1_routes"] = {
            "Nrel_Z_times_tail_multiplier_times_a": Nrel_Z * tail_multiplier * a,
            "Nrel_times_tail_multiplier_Z_times_a": Nrel * tail_multiplier_Z * a,
            "Nrel_times_tail_multiplier_times_a_Z": Nrel * tail_multiplier * a_Z,
        }
    else:
        result["M_Z"] = None
        result["Z1_routes"] = None
    return result


def formal_source_expression(
    *,
    common_S_t_independent: bool = False,
    common_delta: bool = False,
    delta_Z_zero: bool = False,
    common_a_pre_A: bool = False,
    common_M: bool = False,
    common_mu_weight: bool = False,
    source_manifest_schema: str = "gamma_weighted_energy_delta_checkpoint_v1",
) -> dict[str, Any]:
    """Return the unevaluated source and full physical-weight identities."""

    for name, value in (
        ("common_S_t_independent", common_S_t_independent),
        ("common_delta", common_delta),
        ("delta_Z_zero", delta_Z_zero),
        ("common_a_pre_A", common_a_pre_A),
        ("common_M", common_M),
        ("common_mu_weight", common_mu_weight),
    ):
        if not isinstance(value, bool):
            raise GammaWeightedEnergyDeltaInputError(f"{name} must be bool")
    return {
        "schema": SCHEMA,
        "kind": "exact_unevaluated_complete_Gamma_heat_difference",
        "source_assumptions": dict(SOURCE_ASSUMPTIONS),
        "source_definition": {
            "function": "Gamma_future_functions(g, delta, epsilon, S)",
            "collar": {
                "t_interval": "[0,3]",
                "A": "sigma(t)*C(t)",
                "C": "1-phi(t)",
                "W": "(1-sigma(t)) + sigma(t)*phi(t)",
                "pre": "1-epsilon*W(t)",
                "D": "A*Dhat",
                "square": "D*(2*pre-a*S*D)",
                "energy_integral": "Integral_t=0^3 exp(-delta*t)*square dt",
            },
            "outer": {
                "t_interval": "[3,infinity)",
                "A": "1",
                "pre": "1",
                "D": "Dhat",
                "square": "Dhat*(2-a*S*Dhat)",
                "energy_integral": "Integral_t=3^infinity exp(-delta*t)*square dt",
                "infinity": "mathematical_positive_infinity",
            },
            "energy_hat": (
                "Integral_t=0^infinity exp(-delta*t)*square dt = collar + outer"
            ),
            "Dhat": "original source Dhat reference; no provider call",
            "S": "original source S reference; no double-exponential materialization",
        },
        "common_specialization": {
            "gate": common_a_pre_A,
            "D_F": "DF=A*Dhat_F",
            "D_P": "DF + DeltaS*K",
            "K": "A*K0",
            "SF": "S_F",
            "DeltaS": "S_P-S_F",
            "X": "SF*DF",
            "Y0": "DF+SF*K",
            "G0": "2*(pre-a*X)*Y0",
            "G1": "2*(pre-a*X)*K-a*Y0^2",
            "G2": "-2*a*Y0*K",
            "G3": "-a*K^2",
            "Delta_S_square": "DeltaS*(G0+DeltaS*G1+DeltaS^2*G2+DeltaS^3*G3)",
            "SF_and_DeltaS_not_combined": True,
            "no_K1_as_total_K_Z": True,
        },
        "general_difference": {
            "XF": "SF*DF",
            "DeltaX": "S_P*DeltaD + DeltaS*DF",
            "identity": (
                "DeltaX*[2*pre_P-a_P*(2*XF+DeltaX)] + "
                "2*(pre_P-pre_F)*XF - (a_P-a_F)*XF^2"
            ),
            "DeltaD": "independent source difference when a/pre/A are not common",
            "not_replaced_by_DeltaS_K_without_gate": True,
        },
        "ordinary_Z1": {
            "inputs": "total jets a_Z, pre_Z, SF_Z, DeltaS_Z, DF_Z, K_Z",
            "X_Z": "SF_Z*DF + SF*DF_Z",
            "Y0_Z": "DF_Z + SF_Z*K + SF*K_Z",
            "G0_Z": "2*((pre_Z-a_Z*X-a*X_Z)*Y0 + (pre-a*X)*Y0_Z)",
            "G1_Z": "2*((pre_Z-a_Z*X-a*X_Z)*K + (pre-a*X)*K_Z) - a_Z*Y0^2 - 2*a*Y0*Y0_Z",
            "G2_Z": "-2*(a_Z*Y0*K+a*Y0_Z*K+a*Y0*K_Z)",
            "G3_Z": "-(a_Z*K^2+2*a*K*K_Z)",
            "G_Z": (
                "sum_j DeltaS^j*G_j_Z + DeltaS_Z*(G1 + "
                "2*DeltaS*G2 + 3*DeltaS^2*G3)"
            ),
            "Delta_S_square_Z": "DeltaS_Z*G + DeltaS*G_Z",
            "retained_routes": ["DeltaS_Z*G", "DeltaS*G_Z"],
            "K_Z": "A_Z*K0 + A*K0_Z; complete K0 total Z is required",
            "common_delta_source": "A_Z=0 only if the source's common delta makes A Z-independent",
            "nonzero_SF_Z_and_DeltaS_Z_retained": True,
        },
        "complete_energy": {
            "S_energy_hat": "Integral_0^infinity exp(-delta*t)*(S*square) dt",
            "common_S_t_independent_gate": common_S_t_independent,
            "common_delta_gate": common_delta,
            "Delta_S_energy_hat": (
                "Integral_0^infinity exp(-delta*t)*Delta(S*square) dt"
                if common_S_t_independent and common_delta
                else "unevaluated P/F weighted difference; common S and common delta are not both asserted"
            ),
            "noncommon_delta": {
                "w_P": "exp(-delta_P*t)",
                "w_F": "exp(-delta_F*t)",
                "Delta_w": "w_F*expm1(-t*(delta_P-delta_F))",
                "Delta_w_times_H": "w_P*DeltaH + Delta_w*H_F",
                "direct_common_weight_forbidden": True,
            },
            "delta_Z_zero_gate": delta_Z_zero,
            "weight_Z_route": (
                "no -delta_Z*t weight route"
                if delta_Z_zero
                else "add -delta_Z*Integral_0^infinity t*exp(-delta*t)*Delta(S*square) dt"
            ),
            "weight_Z_general": (
                "w_P_Z*DeltaH + w_P*DeltaH_Z + Delta_w_Z*H_F + Delta_w*H_F_Z"
            ),
            "collar_and_outer_retained": True,
            "full_improper_infinity_retained": True,
        },
        "future_heat": {
            "M": "Nrel*tail_multiplier*a",
            "Delta_heat": "M_P*Delta(S*energy_hat) + DeltaM*(S_F*energy_hat_F)",
            "DeltaM": "M_P-M_F",
            "common_M_gate": common_M,
            "Z1": (
                "M_P_Z*Delta(SE) + M_P*Delta(SE)_Z + "
                "DeltaM_Z*(SE)_F + DeltaM*(SE)_F,Z"
            ),
            "Nrel": "(1/4)*exp(-2*mu*(100+Lrel))*(1+Z^2)^2",
            "Nrel_Z_over_Nrel": "4*Z/(1+Z^2) when the displayed Nrel source form applies",
            "Nrel_Z_must_be_retained": True,
            "complete_factor_gate": "direct M*Delta(SE) only after P/F M equality is proved",
        },
        "weighted_future": {
            "factor": "B=-mu*exp(-26)/2",
            "Delta_weighted_future": "B_P*Delta_heat + DeltaB*heat_F",
            "common_mu_weight_gate": common_mu_weight,
            "direct_multiply_gate": common_mu_weight,
            "Z1": "B_P_Z*Delta_heat + B_P*Delta_heat_Z + DeltaB_Z*heat_F + DeltaB*heat_F_Z",
        },
        "source_binding": {
            "checkpoint_schema": source_manifest_schema,
            "raw_inputs_authenticated_before_decompression": True,
            "direct_json_pointer_only": True,
            "no_source_graph_traversal": True,
            "previous_result": "tmp/gamma_correlated_scale_delta_v1/result.json",
            "required_gamma_roles": list(SOURCE_POINTER_NAMES),
            "required_complete_future_roles": list(COMPLETE_FUTURE_POINTER_NAMES),
        },
        "classification": {
            "archived_exact_source_function": "role references and raw SHA only; not values",
            "synthetic_point": "caller-owned arithmetic point; not historical source",
            "uncertified_rounding": "rounded point estimate; no certificate",
            "source_values_computed": False,
            "full_future_P_minus_F_claimed": False,
            "old_error_deleted": False,
        },
        "gates": {
            "common_a_pre_A_required_for_G_polynomial": True,
            "common_S_t_independent_required_to_exchange_delta_and_integral": True,
            "common_delta_required_to_exchange_one_weight": True,
            "delta_Z_zero_required_for_weightless_Z1": True,
            "same_M_required_for_direct_heat_multiply": True,
            "same_mu_and_weight_required_for_direct_weighted_multiply": True,
            "A_Z_and_complete_K0_Z_required": True,
            "source_S_and_Dhat_values_not_materialized": True,
        },
        **SOURCE_VALUE_FLAGS,
    }


def formal_energy_expression(**kwargs: Any) -> dict[str, Any]:
    """Alias emphasizing the complete energy source expression."""

    return formal_source_expression(**kwargs)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _relative_path(root: Path, raw: Any, label: str) -> str:
    _require(isinstance(raw, str) and bool(raw), f"{label} path is not text")
    path = Path(raw)
    if path.is_absolute():
        try:
            return path.resolve().relative_to(root.resolve()).as_posix()
        except ValueError as exc:
            raise GammaWeightedEnergyDeltaSourceError(
                f"{label} path is outside repository"
            ) from exc
    _require(".." not in path.parts, f"{label} path escapes repository")
    return path.as_posix()


def _manifest_hashes(root: Path, source: Mapping[str, Any]) -> dict[str, str]:
    """Collect direct path/hash entries; never traverse a source graph."""

    result: dict[str, str] = {}

    def add(raw_path: Any, raw_digest: Any, label: str) -> None:
        relative = _relative_path(root, raw_path, label)
        _require(
            isinstance(raw_digest, str) and len(raw_digest) == 64,
            f"{label} has no SHA-256",
        )
        digest = raw_digest.lower()
        _require(
            all(char in "0123456789abcdef" for char in digest),
            f"{label} has invalid SHA-256",
        )
        prior = result.get(relative)
        _require(prior is None or prior == digest, f"conflicting SHA-256 for {relative}")
        result[relative] = digest

    for key in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256", "code_sha256"):
        mapping = source.get(key)
        if isinstance(mapping, Mapping):
            for raw_path, raw_digest in mapping.items():
                add(raw_path, raw_digest, f"source.{key}[{raw_path!r}]")
    files = source.get("files")
    _require(isinstance(files, list) and bool(files), "source files manifest missing")
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        add(item.get("path"), item.get("sha256"), f"source.files[{index}]")
    return result


def _json_pointer(document: Any, pointer: str) -> Any:
    """Resolve one RFC 6901 pointer only; no graph or reference traversal."""

    if not isinstance(pointer, str) or not pointer.startswith("/"):
        raise GammaWeightedEnergyDeltaSourceError("JSON pointer must begin with '/'")
    current = document
    for raw_token in pointer.split("/")[1:]:
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, Mapping):
            if token not in current:
                raise GammaWeightedEnergyDeltaSourceError(
                    f"JSON pointer component is absent: {token}"
                )
            current = current[token]
        elif isinstance(current, list):
            try:
                index = int(token)
            except ValueError as exc:
                raise GammaWeightedEnergyDeltaSourceError(
                    f"JSON pointer list component is not an index: {token}"
                ) from exc
            if index < 0 or index >= len(current):
                raise GammaWeightedEnergyDeltaSourceError(
                    f"JSON pointer list index is absent: {token}"
                )
            current = current[index]
        else:
            raise GammaWeightedEnergyDeltaSourceError(
                "JSON pointer attempted to descend through a scalar"
            )
    return current


def _pointer_groups(source: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    groups = source.get("source_pointers")
    if groups is None:
        groups = source.get("json_pointers")
    _require(isinstance(groups, Mapping), "source_pointers mapping is missing")
    result: dict[str, Mapping[str, Any]] = {}
    for name, value in groups.items():
        if isinstance(value, Mapping):
            result[str(name)] = value
    _require("gamma" in result, "source_pointers.gamma is missing")
    _require("complete_future" in result, "source_pointers.complete_future is missing")
    return result


def _read_pointer_group(
    root: Path,
    group_name: str,
    group: Mapping[str, Any],
    raw_hashes: Mapping[str, str],
) -> dict[str, Any]:
    raw_path = group.get("path", group.get("archive"))
    if raw_path is None:
        return {
            "path": None,
            "sha256": None,
            "encoding": None,
            "roles": {},
            "raw_authenticated_before_decompression": False,
            "source_graph_traversed": False,
            "available": False,
            "missing_reason": "pointer source path was not declared",
            "receipt_pin": {"declared": False, "checked": False, "accepted": False},
        }
    relative = _relative_path(root, raw_path, f"source_pointers.{group_name}")
    path = root / relative
    expected = raw_hashes.get(relative)
    _require(expected is not None, f"pointer source is not hash-pinned: {relative}")
    _require(path.is_file(), f"pointer source is missing: {relative}")
    _require(_sha256(path) == expected, f"pointer source changed: {relative}")
    encoding = group.get("encoding", "gzip_json" if relative.endswith(".gz") else "json")
    # This branch is reached only after the raw file has been authenticated.
    if encoding in ("gzip_json", "json.gz", "gzip") or relative.endswith(".gz"):
        document = json.loads(gzip.decompress(path.read_bytes()))
    elif encoding in ("json", "plain_json"):
        document = json.loads(path.read_bytes())
    else:
        raise GammaWeightedEnergyDeltaSourceError(
            f"unsupported pointer encoding for {group_name}: {encoding}"
        )
    pointers = group.get("pointers")
    if not isinstance(pointers, Mapping):
        pointers = {}
    selected: dict[str, Any] = {}
    for role, pointer in pointers.items():
        role_name = str(role)
        try:
            reference = _json_pointer(document, pointer)
        except GammaWeightedEnergyDeltaSourceError as exc:
            # A missing saved role is evidence for a false gate, not permission
            # to call a nearby node or claim that a value was found.
            selected[role_name] = {
                "path": relative,
                "json_pointer": pointer,
                "available": False,
                "missing_reason": str(exc),
                "reference_is_value": False,
            }
        else:
            selected[role_name] = {
                "path": relative,
                "json_pointer": pointer,
                "reference": reference,
                "available": True,
                "reference_is_value": False,
            }

    receipt_spec = group.get("receipt")
    receipt_pin: dict[str, Any] = {
        "declared": isinstance(receipt_spec, Mapping),
        "checked": False,
        "accepted": False,
    }
    if isinstance(receipt_spec, Mapping):
        receipt_raw = receipt_spec.get("path")
        if receipt_raw is not None:
            receipt_relative = _relative_path(root, receipt_raw, f"{group_name} receipt")
            receipt_path = root / receipt_relative
            receipt_expected = raw_hashes.get(receipt_relative)
            _require(receipt_expected is not None, f"receipt is not hash-pinned: {receipt_relative}")
            _require(receipt_path.is_file(), f"receipt is missing: {receipt_relative}")
            _require(_sha256(receipt_path) == receipt_expected, f"receipt changed: {receipt_relative}")
            receipt_document = json.loads(receipt_path.read_bytes())
            receipt_pin.update({"path": receipt_relative, "sha256": receipt_expected})
            receipt_pointers = receipt_spec.get("pointers", {})
            if isinstance(receipt_pointers, Mapping):
                receipt_pin["roles"] = {}
                for role, pointer in receipt_pointers.items():
                    role_name = str(role)
                    try:
                        reference = _json_pointer(receipt_document, pointer)
                    except GammaWeightedEnergyDeltaSourceError as exc:
                        receipt_pin["roles"][role_name] = {
                            "json_pointer": pointer,
                            "available": False,
                            "missing_reason": str(exc),
                        }
                    else:
                        receipt_pin["roles"][role_name] = {
                            "json_pointer": pointer,
                            "reference": reference,
                            "available": True,
                            "reference_is_value": False,
                        }
                archive_pin = receipt_pin["roles"].get("archive_sha256")
                receipt_pin["checked"] = True
                receipt_pin["accepted"] = bool(
                    archive_pin
                    and archive_pin.get("available") is True
                    and archive_pin.get("reference") == expected
                )
            else:
                receipt_pin["checked"] = False
                receipt_pin["accepted"] = False
    return {
        "path": relative,
        "sha256": expected,
        "encoding": encoding,
        "roles": selected,
        "raw_authenticated_before_decompression": True,
        "source_graph_traversed": False,
        "available": True,
        "receipt_pin": receipt_pin,
    }


def bind_archived_source(
    root: Any,
    source: Mapping[str, Any],
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Authenticate raw source bytes, then read only declared JSON pointers.

    The historical Python files are hashed but never imported.  Gzip decoding
    occurs only after every manifest path has been hashed.  The returned role
    references are provenance metadata, not source values.
    """

    repository = Path(root).resolve()
    if hashes is None:
        hashes = _manifest_hashes(repository, source)
    normalized = {
        _relative_path(repository, path, "hash path"): str(digest).lower()
        for path, digest in hashes.items()
    }
    for relative, expected in normalized.items():
        path = repository / relative
        _require(path.is_file(), f"frozen input is missing: {relative}")
        _require(_sha256(path) == expected, f"frozen input changed: {relative}")

    previous_path = "tmp/gamma_correlated_scale_delta_v1/result.json"
    previous_expected = normalized.get(previous_path)
    _require(previous_expected is not None, "previous Gamma delta result is not hash-pinned")
    previous = json.loads((repository / previous_path).read_bytes())
    _require(
        previous.get("component") == "gamma_correlated_scale_delta_v1"
        and previous.get("status") == "returned",
        "previous Gamma delta result is not a successful returned receipt",
    )
    previous_result = previous.get("result")
    _require(isinstance(previous_result, Mapping), "previous Gamma delta result body is missing")
    _require(
        previous_result.get("checks_passed") is True
        and previous.get("source_values_computed") is False
        and previous.get("source_value_computation") is False,
        "previous Gamma delta receipt flags are not accepted",
    )

    groups = _pointer_groups(source)
    gamma = _read_pointer_group(repository, "gamma", groups["gamma"], normalized)
    complete = _read_pointer_group(
        repository, "complete_future", groups["complete_future"], normalized
    )
    gamma_roles = gamma["roles"]
    complete_roles = complete["roles"]
    gamma_aliases = {
        "collar_energy": ("collar_energy", "collar"),
        "outer_energy": ("outer_energy", "outer", "exterior_energy"),
    }
    gamma_availability = {
        name: any(gamma_roles.get(alias, {}).get("available") is True for alias in aliases)
        for name, aliases in {
            **{name: (name,) for name in SOURCE_POINTER_NAMES if name not in gamma_aliases},
            **gamma_aliases,
        }.items()
    }
    complete_availability = {
        name: complete_roles.get(name, {}).get("available") is True
        for name in COMPLETE_FUTURE_POINTER_NAMES
    }
    return {
        "schema": SCHEMA,
        "source_namespace": "historical_Gamma_and_complete_future_archives",
        "gamma": gamma,
        "complete_future": complete,
        "previous_result": {
            "path": previous_path,
            "sha256": previous_expected,
            "component": previous.get("component"),
            "status": previous.get("status"),
            "source_values_computed": False,
            "source_value_computation": False,
        },
        "raw_inputs_authenticated_before_decompression": True,
        "source_graph_traversed": False,
        "source_refs_are_not_values": True,
        "gamma_role_gates": gamma_availability,
        "complete_future_role_gates": complete_availability,
        "missing_gamma_roles": [name for name, ok in gamma_availability.items() if not ok],
        "missing_complete_future_roles": [
            name for name, ok in complete_availability.items() if not ok
        ],
        "receipt_pins": {
            "gamma": dict(gamma.get("receipt_pin", {})),
            "complete_future": dict(complete.get("receipt_pin", {})),
        },
        "P_F_owner_correspondence_proved": False,
        "common_a_pre_A_proved": False,
        "full_future_P_F_proved": False,
        **SOURCE_VALUE_FLAGS,
    }


def _fixture_exact_sympy() -> dict[str, Any]:
    """Independent exact SymPy fixture; never calls a source or quadrature."""

    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - environment-dependent optional extra
        return {
            "fixture": "independent_exact_sympy_common_and_Z1",
            "checks_passed": False,
            "sympy_available": False,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "source_values_computed": False,
            "certified": False,
        }

    q = sp.Rational
    a = q(3, 5)
    pre = q(7, 6)
    sf = q(4, 5)
    ds = q(1, 7)
    df = q(2, 3)
    k = q(5, 8)
    z = sp.symbols("z")
    az, prez, sfz, dsz, dfz, kz = q(2, 11), -q(1, 13), q(3, 17), -q(2, 19), q(1, 23), q(4, 29)

    def square_product(S: Any, D: Any, aa: Any, pp: Any) -> Any:
        return S * D * (2 * pp - aa * S * D)

    point = point_coefficients(a, pre, sf, ds, df, k, source_kind="synthetic_point")
    direct = sp.expand(
        square_product(sf + ds, df + ds * k, a, pre)
        - square_product(sf, df, a, pre)
    )
    factorized = sp.expand(point.materialize_total())
    common_identity = sp.simplify(direct - factorized) == 0
    zero_point = point_coefficients(a, pre, sf, 0, df, k)
    zero_identity = sp.simplify(zero_point.materialize_total()) == 0

    affine = lambda value, derivative: value + derivative * z
    z1 = ordinary_z1_coefficients(
        a,
        pre,
        sf,
        ds,
        df,
        k,
        a_Z=az,
        pre_Z=prez,
        SF_Z=sfz,
        DeltaS_Z=dsz,
        DF_Z=dfz,
        K_Z=kz,
        source_kind="synthetic_point",
    )
    direct_z = sp.diff(
        square_product(
            affine(sf, sfz) + affine(ds, dsz),
            affine(df, dfz) + affine(ds, dsz) * affine(k, kz),
            affine(a, az),
            affine(pre, prez),
        )
        - square_product(affine(sf, sfz), affine(df, dfz), affine(a, az), affine(pre, prez)),
        z,
    ).subs(z, 0)
    z_identity = sp.simplify(direct_z - z1.total_Z) == 0
    route_ds_z = sp.expand(z1.DeltaS_Z_times_G)
    route_ds = sp.expand(z1.DeltaS_times_G_Z)
    route_nonzero = route_ds_z != 0 and route_ds != 0
    omission_detected = sp.simplify(direct_z - route_ds) != 0 and sp.simplify(direct_z - route_ds_z) != 0

    ds_zero_with_z = ordinary_z1_coefficients(
        a,
        pre,
        sf,
        0,
        df,
        k,
        a_Z=az,
        pre_Z=prez,
        SF_Z=sfz,
        DeltaS_Z=dsz,
        DF_Z=dfz,
        K_Z=kz,
        source_kind="synthetic_point",
    )
    ds_zero_derivative_identity = sp.simplify(ds_zero_with_z.total_Z - ds_zero_with_z.DeltaS_Z_times_G) == 0

    general_sp = q(9, 10)
    general_sf = q(2, 5)
    general_dp = q(7, 9)
    general_df = q(1, 3)
    general_aP = q(11, 13)
    general_aF = q(4, 9)
    general_preP = q(8, 7)
    general_preF = q(5, 6)
    general_ds = general_sp - general_sf
    general_dd = general_dp - general_df
    general = general_source_difference(
        SP=general_sp,
        SF=general_sf,
        DP=general_dp,
        DF=general_df,
        aP=general_aP,
        aF=general_aF,
        preP=general_preP,
        preF=general_preF,
        DeltaS=general_ds,
        DeltaD=general_dd,
    )

    general_direct = sp.expand(
        square_product(general_sp, general_dp, general_aP, general_preP)
        - square_product(general_sf, general_df, general_aF, general_preF)
    )
    general_identity = sp.simplify(general["value"] - general_direct) == 0
    general_delta_x_identity = sp.simplify(
        general["DeltaX"] - (general_sp * general_dp - general_sf * general_df)
    ) == 0
    general_differing_inputs = bool(
        general_aP != general_aF
        and general_preP != general_preF
        and general_dp != general_df
        and general["DeltaD"] == general_dd
        and general["DeltaD"] != general_ds * k
    )

    heat_dse = q(3, 8)
    heat_dse_z = q(5, 13)
    heat_sef = q(4, 9)
    heat_sef_z = -q(2, 15)
    heat_mp = q(7, 10)
    heat_mp_z = q(3, 11)
    heat_mf = q(5, 12)
    heat_mf_z = -q(1, 14)
    heat = assemble_heat_difference(
        DeltaSE=heat_dse,
        SE_F=heat_sef,
        M_P=heat_mp,
        M_F=heat_mf,
        DeltaSE_Z=heat_dse_z,
        SE_F_Z=heat_sef_z,
        M_P_Z=heat_mp_z,
        M_F_Z=heat_mf_z,
    )
    heat_z1 = heat["Z1"]
    heat_route_expected = {
        "M_P_Z_times_DeltaSE": heat_mp_z * heat_dse,
        "M_P_times_DeltaSE_Z": heat_mp * heat_dse_z,
        "DeltaM_Z_times_SE_F": (heat_mp_z - heat_mf_z) * heat_sef,
        "DeltaM_times_SE_F_Z": (heat_mp - heat_mf) * heat_sef_z,
    }
    heat_direct_z = sp.diff(
        affine(heat_mp, heat_mp_z)
        * (affine(heat_sef, heat_sef_z) + affine(heat_dse, heat_dse_z))
        - affine(heat_mf, heat_mf_z) * affine(heat_sef, heat_sef_z),
        z,
    ).subs(z, 0)
    heat_routes_exact = bool(
        all(
            sp.simplify(heat_z1[name] - expected) == 0
            for name, expected in heat_route_expected.items()
        )
    )
    heat_routes_nonzero = bool(
        all(sp.simplify(heat_z1[name]) != 0 for name in heat_route_expected)
    )
    heat_z_identity = sp.simplify(heat_direct_z - heat_z1["Delta_heat_Z"]) == 0

    weighted_dheat = q(2, 9)
    weighted_dheat_z = -q(4, 17)
    weighted_hf = q(5, 7)
    weighted_hf_z = q(3, 16)
    weighted_bp = q(7, 8)
    weighted_bp_z = q(2, 13)
    weighted_bf = q(1, 6)
    weighted_bf_z = -q(5, 19)
    weighted = assemble_weighted_future_difference(
        B_P=weighted_bp,
        B_F=weighted_bf,
        DeltaHeat=weighted_dheat,
        Heat_F=weighted_hf,
        DeltaHeat_Z=weighted_dheat_z,
        Heat_F_Z=weighted_hf_z,
        B_P_Z=weighted_bp_z,
        B_F_Z=weighted_bf_z,
    )
    weighted_z1 = weighted["Z1"]
    weighted_route_expected = {
        "B_P_Z_times_DeltaHeat": weighted_bp_z * weighted_dheat,
        "B_P_times_DeltaHeat_Z": weighted_bp * weighted_dheat_z,
        "DeltaB_Z_times_Heat_F": (weighted_bp_z - weighted_bf_z) * weighted_hf,
        "DeltaB_times_Heat_F_Z": (weighted_bp - weighted_bf) * weighted_hf_z,
    }
    weighted_direct_z = sp.diff(
        affine(weighted_bp, weighted_bp_z) * affine(weighted_dheat, weighted_dheat_z)
        + (
            affine(weighted_bp, weighted_bp_z)
            - affine(weighted_bf, weighted_bf_z)
        )
        * affine(weighted_hf, weighted_hf_z),
        z,
    ).subs(z, 0)
    weighted_routes_exact = bool(
        all(
            sp.simplify(weighted_z1[name] - expected) == 0
            for name, expected in weighted_route_expected.items()
        )
    )
    weighted_routes_nonzero = bool(
        all(sp.simplify(weighted_z1[name]) != 0 for name in weighted_route_expected)
    )
    weighted_z_identity = sp.simplify(
        weighted_direct_z - weighted_z1["Delta_weighted_heat_contribution_Z"]
    ) == 0

    delta_P, delta_F, DeltaH, H_F = sp.symbols("delta_P delta_F DeltaH H_F")
    delta_P_Z, delta_F_Z, DeltaH_Z, H_F_Z = sp.symbols(
        "delta_P_Z delta_F_Z DeltaH_Z H_F_Z"
    )
    t = sp.symbols("t")
    time_weight = time_weight_difference_structure(
        delta_P=delta_P,
        delta_F=delta_F,
        DeltaH=DeltaH,
        H_F=H_F,
        delta_Z=delta_P_Z - delta_F_Z,
        DeltaH_Z=DeltaH_Z,
        H_F_Z=H_F_Z,
    )
    delta_delta = delta_P - delta_F
    delta_w_chart = sp.exp(-(delta_F + delta_F_Z * z) * t) * (
        sp.exp(
            -t
            * ((delta_P + delta_P_Z * z) - (delta_F + delta_F_Z * z))
        )
        - 1
    )
    delta_w_direct_z = sp.diff(delta_w_chart, z).subs(z, 0)
    delta_w_expected_z = (
        (-t * delta_F_Z * sp.exp(-delta_F * t))
        * (sp.exp(-t * delta_delta) - 1)
        + sp.exp(-delta_F * t)
        * sp.exp(-t * delta_delta)
        * (-t * (delta_P_Z - delta_F_Z))
    )
    time_delta_w_formula_exact = time_weight["ordinary_Z1"]["Delta_w_Z"] == (
        "w_F_Z*expm1(-t*Delta_delta) + "
        "w_F*exp(-t*Delta_delta)*(-t*Delta_delta_Z)"
    )
    time_delta_w_derivative_exact = sp.simplify(
        delta_w_direct_z - delta_w_expected_z
    ) == 0

    missing_heat_z_raises = False
    try:
        assemble_heat_difference(
            DeltaSE=heat_dse,
            SE_F=heat_sef,
            M_P=heat_mp,
            M_F=heat_mf,
            DeltaSE_Z=heat_dse_z,
        )
    except GammaWeightedEnergyDeltaInputError:
        missing_heat_z_raises = True
    missing_weighted_z_raises = False
    try:
        assemble_weighted_future_difference(
            B_P=weighted_bp,
            B_F=weighted_bf,
            DeltaHeat=weighted_dheat,
            Heat_F=weighted_hf,
            DeltaHeat_Z=weighted_dheat_z,
        )
    except GammaWeightedEnergyDeltaInputError:
        missing_weighted_z_raises = True

    missing_jet_z_raises = False
    try:
        ordinary_z1_coefficients(a, pre, sf, ds, df, k)
    except GammaWeightedEnergyDeltaInputError:
        missing_jet_z_raises = True
    missing_z_derivative_raises = (
        missing_heat_z_raises and missing_weighted_z_raises and missing_jet_z_raises
    )
    checks = {
        "common_factorization_exact": common_identity,
        "ds_zero_exact": zero_identity,
        "general_identity_exact": general_identity,
        "general_differing_a_pre_D_inputs": general_differing_inputs,
        "general_delta_X_exact": general_delta_x_identity,
        "ordinary_Z1_exact": z_identity,
        "both_Z1_routes_nonzero": route_nonzero,
        "nonzero_ds_Z_route_not_omitted": omission_detected,
        "ds_zero_keeps_ds_Z_times_G": ds_zero_derivative_identity,
        "assemble_heat_Z1_four_routes_exact": heat_routes_exact,
        "assemble_heat_Z1_four_routes_nonzero": heat_routes_nonzero,
        "assemble_heat_Z1_exact": heat_z_identity,
        "assemble_weighted_future_Z1_four_routes_exact": weighted_routes_exact,
        "assemble_weighted_future_Z1_four_routes_nonzero": weighted_routes_nonzero,
        "assemble_weighted_future_Z1_exact": weighted_z_identity,
        "time_weight_Delta_w_Z_formula_exact": time_delta_w_formula_exact,
        "time_weight_Delta_w_Z_derivative_exact": time_delta_w_derivative_exact,
        "missing_Z_derivative_raises": missing_z_derivative_raises,
    }
    return {
        "schema": SCHEMA,
        "fixture": "independent_exact_sympy_common_and_Z1",
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "sympy_available": True,
        "exact_arithmetic": True,
        "nonzero_inputs": {
            "DeltaS_Z": str(dsz),
            "SF_Z": str(sfz),
            "a_Z": str(az),
            "pre_Z": str(prez),
            "DF_Z": str(dfz),
            "K_Z": str(kz),
        },
        "routes": {
            "DeltaS_Z_times_G": str(route_ds_z),
            "DeltaS_times_G_Z": str(route_ds),
        },
        "general_differing_inputs": {
            "SP": str(general_sp),
            "SF": str(general_sf),
            "DP": str(general_dp),
            "DF": str(general_df),
            "aP": str(general_aP),
            "aF": str(general_aF),
            "preP": str(general_preP),
            "preF": str(general_preF),
            "DeltaS": str(general_ds),
            "DeltaD": str(general_dd),
        },
        "assembled_Z1_routes": {
            "heat": {name: str(heat_z1[name]) for name in heat_route_expected},
            "weighted_future": {
                name: str(weighted_z1[name]) for name in weighted_route_expected
            },
        },
        "time_weight": {
            "Delta_w_Z_formula": time_weight["ordinary_Z1"]["Delta_w_Z"],
            "Delta_w_Z_direct": str(delta_w_direct_z),
            "Delta_w_Z_expected": str(delta_w_expected_z),
        },
        "missing_Z_derivative_errors": {
            "ordinary_z1_coefficients": missing_jet_z_raises,
            "assemble_heat_difference": missing_heat_z_raises,
            "assemble_weighted_future_difference": missing_weighted_z_raises,
        },
        "source_values_computed": False,
        "source_value_computation": False,
        "certified": False,
        "uncertified_rounding": False,
    }


def verify_synthetic_fixture(
    root: Any = None,
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Bind authenticated role metadata (when supplied) and run one fixture."""

    binding = None
    if root is not None:
        if source is None:
            source_path = Path(root) / "tmp" / "gamma_weighted_energy_delta_v1" / "source.json"
            source = json.loads(source_path.read_bytes())
        binding = bind_archived_source(root, source, hashes)
    fixture = _fixture_exact_sympy()
    result: dict[str, Any] = {
        "schema": SCHEMA,
        "fixture": fixture,
        "formal_source_expression": formal_source_expression(),
        "source_values_computed": False,
        "source_value_computation": False,
        "historical_source_values_computed": False,
        "certified": False,
        "source_binding_performed": binding is not None,
        "source_binding": binding,
        "claims": {
            "archived_exact_source_function_referenced": binding is not None,
            "synthetic_points_only": True,
            "uncertified_rounding": False,
            "full_future_P_F_claimed": False,
            "old_error_deleted": False,
        },
    }
    result["checks_passed"] = bool(fixture.get("checks_passed") is True)
    return result


# Deliberately descriptive aliases for callers that use the source notation.
build_coefficients = point_coefficients
evaluate_coefficients = point_coefficients
point_energy_delta = point_coefficients
build_z1_coefficients = ordinary_z1_coefficients
point_z1_energy_delta = ordinary_z1_coefficients
general_delta_identity = general_source_difference
assemble_future_heat_difference = assemble_heat_difference
assemble_weighted_future = assemble_weighted_future_difference
build_formal_source_expression = formal_source_expression


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "PHASE",
    "SOURCE_VALUE_FLAGS",
    "SOURCE_ASSUMPTIONS",
    "SOURCE_POINTER_NAMES",
    "COMPLETE_FUTURE_POINTER_NAMES",
    "GammaWeightedEnergyDeltaInputError",
    "GammaWeightedEnergyDeltaSourceError",
    "Z1Jet",
    "GammaEnergyDeltaCoefficients",
    "GammaEnergyDeltaZ1",
    "point_coefficients",
    "build_coefficients",
    "evaluate_coefficients",
    "point_energy_delta",
    "ordinary_z1_coefficients",
    "build_z1_coefficients",
    "point_z1_energy_delta",
    "compose_total_K",
    "general_source_difference",
    "general_delta_identity",
    "assemble_heat_difference",
    "assemble_future_heat_difference",
    "assemble_weighted_future_difference",
    "assemble_weighted_future",
    "time_weight_difference_structure",
    "future_prefactor",
    "formal_source_expression",
    "formal_energy_expression",
    "build_formal_source_expression",
    "bind_archived_source",
    "verify_synthetic_fixture",
]
