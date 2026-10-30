"""Non-common P/F ordinary-Z1 algebra for the Gamma heat integrand.

This is a deliberately standalone algebra leaf.  It does not import a
Gamma/source provider, an archived module, a quadrature routine, or any
older calculation.  The only quantity represented here is the pointwise
integrand

    q(S, D, a, pre) = S * D * (2 * pre - a * S * D).

P and F are allowed to have different ``S``, ``D``, ``a`` and ``pre``.  The
returned packets keep the independent differences ``DeltaS``, ``DeltaD``,
``DeltaA`` and ``DeltaPre`` visible.  They also retain the three-factor
decomposition of ``DeltaX`` for ``X = S*D`` instead of replacing it with a
common-source shortcut.

``ordinary_z1_difference`` accepts caller-owned total ordinary-Z one-jets
for every P/F input and returns both branch Z1 values and the complete
P-minus-F Z1, with every product route retained.  A missing derivative is
an input error; an explicitly supplied zero derivative is valid.

The optional SymPy fixture is synthetic and is intentionally not called at
module import.  It is a diagnostic identity check only: it does not bind
source owners, evaluate a physical integral, establish an error bound, or
certify any Navier--Stokes statement.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Mapping, Sequence


SCHEMA = "gamma_noncommon_z1_v1"
COMPONENT = SCHEMA

SOURCE_VALUE_FLAGS = {
    "source_values_computed": False,
    "source_value_computation": False,
    "historical_source_values_computed": False,
    "actual_source_value_computed": False,
    "source_identity_proved": False,
    "certified": False,
}

FORMULA = "q(S,D,a,pre)=S*D*(2*pre-a*S*D)"
DIFFERENCE_FORMULA = (
    "DeltaX*(2*preP-aP*(2*XF+DeltaX)) "
    "+ 2*DeltaPre*XF - DeltaA*XF^2"
)
Z1_FORMULA = "qP_Z-qF_Z, with all P/F and difference product routes retained"


class GammaNoncommonInputError(ValueError):
    """Raised when a non-common point or ordinary-Z1 input is incomplete."""


@dataclass(frozen=True)
class Z1Jet:
    """A caller-owned scalar and its total ordinary-Z derivative."""

    value: Any
    Z: Any

    @property
    def derivative(self) -> Any:
        return self.Z

    def as_dict(self) -> dict[str, Any]:
        return {"value": _display(self.value), "Z": _display(self.Z)}


def _display(value: Any) -> Any:
    """Keep serialization optional and never coerce arithmetic values."""

    if isinstance(value, Mapping):
        return {key: _display(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_display(item) for item in value]
    try:
        json.dumps(value, allow_nan=False)
    except (TypeError, ValueError):
        if hasattr(value, "_mpf_"):
            return str(value)
        if type(value).__module__.split(".", 1)[0] == "sympy":
            return str(value)
        if hasattr(value, "as_dict") and callable(value.as_dict):
            return value.as_dict()
        return str(value)
    return value


def _jet_pair(value: Any, derivative: Any, label: str) -> tuple[Any, Any]:
    """Resolve one required jet without silently inventing a derivative."""

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
        raise GammaNoncommonInputError(
            f"{label} mapping must contain value and Z entries"
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return value[0], value[1]
        raise GammaNoncommonInputError(
            f"{label} sequence jet must have exactly two entries"
        )
    raise GammaNoncommonInputError(
        f"{label} requires an explicit ordinary-Z derivative"
    )


def _optional_delta_pair(
    value: Any,
    derivative: Any,
    label: str,
) -> tuple[Any | None, Any | None, bool]:
    """Resolve an optional declared difference and report whether it was set."""

    if value is None:
        # The value difference can be inferred from the two branch values
        # while its derivative remains an independent caller declaration.
        return None, derivative, False
    if isinstance(value, Z1Jet):
        if derivative is not None:
            return value.value, derivative, True
        return value.value, value.Z, True
    if isinstance(value, Mapping):
        if derivative is not None:
            if "value" in value:
                return value["value"], derivative, True
            if "base" in value:
                return value["base"], derivative, True
            raise GammaNoncommonInputError(
                f"{label} mapping must contain value or base"
            )
        return (*_jet_pair(value, None, label), True)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 2:
            raise GammaNoncommonInputError(
                f"{label} sequence jet must have exactly two entries"
            )
        return value[0], value[1] if derivative is None else derivative, True
    if derivative is None:
        return value, None, True
    return value, derivative, True


def _materialization_flags(
    *, source_kind: str, rounded: bool, delta_factorized: bool
) -> dict[str, Any]:
    allowed = {
        "archived_exact_source_function",
        "synthetic_point",
        "uncertified_rounding",
        "formal_symbolic",
    }
    aliases = {
        "archived": "archived_exact_source_function",
        "archived_exact": "archived_exact_source_function",
        "synthetic": "synthetic_point",
        "rounded": "uncertified_rounding",
    }
    source_kind = aliases.get(source_kind, source_kind)
    if rounded:
        source_kind = "uncertified_rounding"
    if source_kind not in allowed:
        raise GammaNoncommonInputError(
            "source_kind must be archived_exact_source_function, "
            "synthetic_point, uncertified_rounding, or formal_symbolic"
        )
    return {
        "source_kind": source_kind,
        "source_kind_is_caller_declared": True,
        "caller_declared_archived_exact_source_function": (
            source_kind == "archived_exact_source_function"
        ),
        "archived_exact_source_function": False,
        "source_identity_proved": False,
        "synthetic_point": source_kind == "synthetic_point",
        "formal_symbolic": source_kind == "formal_symbolic",
        "uncertified_rounding": source_kind == "uncertified_rounding",
        "DeltaS_retained": True,
        "DeltaD_retained": True,
        "DeltaA_retained": True,
        "DeltaPre_retained": True,
        "DeltaX_factorized": bool(delta_factorized),
        "small_terms_structurally_retained": True,
        "floating_point_cancellation_not_excluded": True,
        "source_values_computed": False,
        "source_value_computation": False,
        "historical_source_values_computed": False,
        "actual_source_value_computed": False,
        "certified": False,
    }


@dataclass(frozen=True)
class NoncommonGammaDifference:
    """Point values and the exact P-minus-F difference packet."""

    SP: Any
    SF: Any
    DP: Any
    DF: Any
    aP: Any
    aF: Any
    preP: Any
    preF: Any
    DeltaS: Any
    DeltaD: Any
    DeltaA: Any
    DeltaPre: Any
    XP: Any
    XF: Any
    DeltaX: Any
    DeltaX_factorized: Any
    qP: Any
    qF: Any
    value: Any
    value_factorized: Any
    direct_factorized_residual: Any
    retained_terms: Mapping[str, Any]
    factorized_terms: Mapping[str, Any]
    difference_consistency_gate: Mapping[str, Any]
    materialization: Mapping[str, Any]

    @property
    def ds(self) -> Any:
        return self.DeltaS

    @property
    def dd(self) -> Any:
        return self.DeltaD

    @property
    def da(self) -> Any:
        return self.DeltaA

    @property
    def dpre(self) -> Any:
        return self.DeltaPre

    @property
    def conditional_factorized_value(self) -> Any:
        return self.value_factorized

    @property
    def values(self) -> dict[str, Any]:
        return {
            "SP": self.SP,
            "SF": self.SF,
            "DP": self.DP,
            "DF": self.DF,
            "aP": self.aP,
            "aF": self.aF,
            "preP": self.preP,
            "preF": self.preF,
            "DeltaS": self.DeltaS,
            "DeltaD": self.DeltaD,
            "DeltaA": self.DeltaA,
            "DeltaPre": self.DeltaPre,
            "XP": self.XP,
            "XF": self.XF,
            "DeltaX": self.DeltaX,
            "DeltaX_factorized": self.DeltaX_factorized,
            "qP": self.qP,
            "qF": self.qF,
            "value": self.value,
            "value_factorized": self.value_factorized,
            "direct_factorized_residual": self.direct_factorized_residual,
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "values": {name: _display(value) for name, value in self.values.items()},
            "retained_terms": {
                name: _display(value) for name, value in self.retained_terms.items()
            },
            "conditional_factorized": {
                "value": _display(self.value_factorized),
                "terms": {
                    name: _display(value)
                    for name, value in self.factorized_terms.items()
                },
                "direct_minus_factorized_residual": _display(
                    self.direct_factorized_residual
                ),
            },
            "difference_consistency_gate": _display(
                self.difference_consistency_gate
            ),
            "formula": DIFFERENCE_FORMULA,
            "integrand_formula": FORMULA,
            "materialization": dict(self.materialization),
            **SOURCE_VALUE_FLAGS,
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "ds": "DeltaS",
            "dd": "DeltaD",
            "da": "DeltaA",
            "dpre": "DeltaPre",
            "delta_s": "DeltaS",
            "delta_d": "DeltaD",
            "delta_a": "DeltaA",
            "delta_pre": "DeltaPre",
            "delta_x": "DeltaX",
            "total": "value",
            "factorized_value": "value_factorized",
            "conditional_factorized_value": "value_factorized",
            "direct_minus_factorized": "direct_factorized_residual",
        }
        if key == "values":
            return self.values
        if key == "Z1":
            return None
        if key == "factorized":
            return self.value_factorized
        if key == "difference_consistency_gate":
            return self.difference_consistency_gate
        return getattr(self, aliases.get(key, key))


def _q(S: Any, D: Any, a: Any, pre: Any) -> Any:
    """The source-shaped point integrand, with no source evaluation."""

    X = S * D
    return X * (2 * pre - a * X)


def _difference_value_from_delta_x(
    DeltaX: Any,
    XF: Any,
    aP: Any,
    DeltaA: Any,
    preP: Any,
    DeltaPre: Any,
) -> tuple[Any, dict[str, Any]]:
    """Evaluate the five signed value routes for a supplied DeltaX."""

    terms = {
        "two_preP_times_DeltaX": 2 * preP * DeltaX,
        "two_DeltaPre_times_XF": 2 * DeltaPre * XF,
        "minus_aP_times_two_XF_DeltaX": -aP * (2 * XF * DeltaX),
        "minus_aP_times_DeltaX_squared": -aP * DeltaX**2,
        "minus_DeltaA_times_XF_squared": -DeltaA * XF**2,
        "DeltaX_bracket": DeltaX * (2 * preP - aP * (2 * XF + DeltaX)),
    }
    value = (
        terms["two_preP_times_DeltaX"]
        + terms["two_DeltaPre_times_XF"]
        + terms["minus_aP_times_two_XF_DeltaX"]
        + terms["minus_aP_times_DeltaX_squared"]
        + terms["minus_DeltaA_times_XF_squared"]
    )
    return value, terms


def _difference_z1_from_delta_x(
    DeltaX: Any,
    DeltaX_Z: Any,
    XF: Any,
    XF_Z: Any,
    aP: Any,
    aP_Z: Any,
    DeltaA: Any,
    DeltaA_Z: Any,
    preP: Any,
    preP_Z: Any,
    DeltaPre: Any,
    DeltaPre_Z: Any,
) -> tuple[Any, dict[str, Any]]:
    """Evaluate all ordinary-Z1 routes for a supplied DeltaX one-jet."""

    bracket = 2 * XF * DeltaX + DeltaX**2
    bracket_Z = (
        2 * XF_Z * DeltaX
        + 2 * XF * DeltaX_Z
        + 2 * DeltaX * DeltaX_Z
    )
    terms = {
        "two_preP_Z_times_DeltaX": 2 * preP_Z * DeltaX,
        "two_preP_times_DeltaX_Z": 2 * preP * DeltaX_Z,
        "two_DeltaPre_Z_times_XF": 2 * DeltaPre_Z * XF,
        "two_DeltaPre_times_XF_Z": 2 * DeltaPre * XF_Z,
        "minus_aP_Z_times_bracket": -aP_Z * bracket,
        "minus_aP_times_bracket_Z": -aP * bracket_Z,
        "minus_DeltaA_Z_times_XF_squared": -DeltaA_Z * XF**2,
        "minus_DeltaA_times_two_XF_XF_Z": -DeltaA * (2 * XF * XF_Z),
        "bracket": bracket,
        "bracket_Z": bracket_Z,
    }
    total_Z = (
        terms["two_preP_Z_times_DeltaX"]
        + terms["two_preP_times_DeltaX_Z"]
        + terms["two_DeltaPre_Z_times_XF"]
        + terms["two_DeltaPre_times_XF_Z"]
        + terms["minus_aP_Z_times_bracket"]
        + terms["minus_aP_times_bracket_Z"]
        + terms["minus_DeltaA_Z_times_XF_squared"]
        + terms["minus_DeltaA_times_two_XF_XF_Z"]
    )
    return total_Z, terms


def point_difference(
    SP: Any,
    SF: Any,
    DP: Any,
    DF: Any,
    aP: Any,
    aF: Any,
    preP: Any,
    preF: Any,
    *,
    DeltaS: Any = None,
    DeltaD: Any = None,
    source_kind: str = "synthetic_point",
    rounded: bool = False,
    _deltaS_declared: bool | None = None,
    _deltaD_declared: bool | None = None,
) -> NoncommonGammaDifference:
    """Return the general non-common P/F value difference.

    ``DeltaS`` and ``DeltaD`` default to the branch differences but may be
    supplied as caller-owned declarations.  ``DeltaX`` is always formed from
    the branch products ``SP*DP-SF*DF``; the separately exposed
    ``DeltaX_factorized`` records the declared ``DeltaS``/``DeltaD`` routes.
    ``value_factorized`` applies the same exact difference formula to that
    declared factorization, but its consistency gate is deliberately
    unpassed: a declaration is not evidence that it equals the branch
    product difference.  Thus even a malformed independent declaration
    cannot make the reported direct branch value silently change.
    """

    ds = SP - SF if DeltaS is None else DeltaS
    dd = DP - DF if DeltaD is None else DeltaD
    da = aP - aF
    dpre = preP - preF
    XP = SP * DP
    XF = SF * DF
    dX = XP - XF
    dX_factorized = SF * dd + DF * ds + ds * dd

    value, terms = _difference_value_from_delta_x(
        dX, XF, aP, da, preP, dpre
    )
    value_factorized, factorized_terms = _difference_value_from_delta_x(
        dX_factorized, XF, aP, da, preP, dpre
    )
    factorized_terms = {
        **factorized_terms,
        "SF_times_DeltaD": SF * dd,
        "DF_times_DeltaS": DF * ds,
        "DeltaS_times_DeltaD": ds * dd,
    }
    terms = {
        **terms,
        "SF_times_DeltaD": SF * dd,
        "DF_times_DeltaS": DF * ds,
        "DeltaS_times_DeltaD": ds * dd,
    }
    direct_factorized_residual = value - value_factorized
    difference_consistency_gate = {
        "required_for_factorized_identity": True,
        "passed": False,
        "consistency_proved": False,
        "status": "unpassed_caller_gate",
        "caller_must_prove": True,
        "declarations_are_not_evidence": True,
        "branch_delta_equality_is_not_checked_here": True,
        "direct_minus_factorized_residual": direct_factorized_residual,
    }
    flags = _materialization_flags(
        source_kind=source_kind,
        rounded=rounded,
        delta_factorized=True,
    )
    flags.update(
        {
            "DeltaS_declared_independently": (
                DeltaS is not None
                if _deltaS_declared is None
                else bool(_deltaS_declared)
            ),
            "DeltaD_declared_independently": (
                DeltaD is not None
                if _deltaD_declared is None
                else bool(_deltaD_declared)
            ),
            "DeltaX_direct_from_branch_products": True,
            "DeltaX_factorized_from_declared_differences": True,
            "conditional_factorized_value_returned": True,
            "difference_consistency_gate_passed": False,
            "difference_consistency_caller_must_prove": True,
            "common_a_pre_D_S_assumed": False,
        }
    )
    return NoncommonGammaDifference(
        SP,
        SF,
        DP,
        DF,
        aP,
        aF,
        preP,
        preF,
        ds,
        dd,
        da,
        dpre,
        XP,
        XF,
        dX,
        dX_factorized,
        _q(SP, DP, aP, preP),
        _q(SF, DF, aF, preF),
        value,
        value_factorized,
        direct_factorized_residual,
        terms,
        factorized_terms,
        difference_consistency_gate,
        flags,
    )


@dataclass(frozen=True)
class NoncommonGammaZ1:
    """Complete ordinary-Z1 packet for a non-common P/F difference."""

    base: NoncommonGammaDifference
    SP_Z: Any
    SF_Z: Any
    DP_Z: Any
    DF_Z: Any
    aP_Z: Any
    aF_Z: Any
    preP_Z: Any
    preF_Z: Any
    DeltaS_Z: Any
    DeltaD_Z: Any
    DeltaA_Z: Any
    DeltaPre_Z: Any
    XP_Z: Any
    XF_Z: Any
    DeltaX_Z: Any
    DeltaX_Z_factorized: Any
    qP_Z: Any
    qF_Z: Any
    total_Z: Any
    total_Z_factorized: Any
    direct_factorized_Z1_residual: Any
    retained_Z_terms: Mapping[str, Any]
    factorized_Z_terms: Mapping[str, Any]
    difference_consistency_gate: Mapping[str, Any]
    materialization: Mapping[str, Any]

    @property
    def value(self) -> Any:
        return self.base.value

    @property
    def value_Z(self) -> Any:
        return self.total_Z

    @property
    def Z1(self) -> Any:
        return self.total_Z

    @property
    def conditional_factorized_Z1(self) -> Any:
        return self.total_Z_factorized

    @property
    def ds_Z(self) -> Any:
        return self.DeltaS_Z

    @property
    def dd_Z(self) -> Any:
        return self.DeltaD_Z

    @property
    def da_Z(self) -> Any:
        return self.DeltaA_Z

    @property
    def dpre_Z(self) -> Any:
        return self.DeltaPre_Z

    @property
    def values(self) -> dict[str, Any]:
        return self.base.values

    @property
    def jets(self) -> dict[str, Any]:
        return {
            "SP_Z": self.SP_Z,
            "SF_Z": self.SF_Z,
            "DP_Z": self.DP_Z,
            "DF_Z": self.DF_Z,
            "aP_Z": self.aP_Z,
            "aF_Z": self.aF_Z,
            "preP_Z": self.preP_Z,
            "preF_Z": self.preF_Z,
            "DeltaS_Z": self.DeltaS_Z,
            "DeltaD_Z": self.DeltaD_Z,
            "DeltaA_Z": self.DeltaA_Z,
            "DeltaPre_Z": self.DeltaPre_Z,
            "XP_Z": self.XP_Z,
            "XF_Z": self.XF_Z,
            "DeltaX_Z": self.DeltaX_Z,
            "DeltaX_Z_factorized": self.DeltaX_Z_factorized,
            "qP_Z": self.qP_Z,
            "qF_Z": self.qF_Z,
            "total_Z": self.total_Z,
            "total_Z_factorized": self.total_Z_factorized,
            "direct_factorized_Z1_residual": self.direct_factorized_Z1_residual,
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "values": {name: _display(value) for name, value in self.values.items()},
            "Z1": {name: _display(value) for name, value in self.jets.items()},
            "retained_Z_terms": {
                name: _display(value)
                for name, value in self.retained_Z_terms.items()
            },
            "conditional_factorized_Z1": {
                "total_Z": _display(self.total_Z_factorized),
                "terms": {
                    name: _display(value)
                    for name, value in self.factorized_Z_terms.items()
                },
                "direct_minus_factorized_residual": _display(
                    self.direct_factorized_Z1_residual
                ),
            },
            "difference_consistency_gate": _display(
                self.difference_consistency_gate
            ),
            "formula": Z1_FORMULA,
            "difference_formula": DIFFERENCE_FORMULA,
            "materialization": dict(self.materialization),
            **SOURCE_VALUE_FLAGS,
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "total": "total_Z",
            "value_Z": "total_Z",
            "delta_s_Z": "DeltaS_Z",
            "delta_d_Z": "DeltaD_Z",
            "delta_a_Z": "DeltaA_Z",
            "delta_pre_Z": "DeltaPre_Z",
            "factorized_Z1": "total_Z_factorized",
            "conditional_factorized_Z1": "total_Z_factorized",
            "direct_minus_factorized_Z1": "direct_factorized_Z1_residual",
        }
        if key == "values":
            return self.values
        if key in {"Z1", "jets"}:
            return self.jets
        if key == "factorized":
            return self.total_Z_factorized
        if key == "difference_consistency_gate":
            return self.difference_consistency_gate
        return getattr(self, aliases.get(key, key))


def _time_formal_fields(
    delta_P: Any,
    delta_F: Any,
    delta_P_Z: Any,
    delta_F_Z: Any,
    t: Any,
) -> dict[str, str]:
    """Return unevaluated exp/expm1 strings for the time-weight one-jet."""

    dp = str(delta_P)
    df = str(delta_F)
    dpz = str(delta_P_Z)
    dfz = str(delta_F_Z)
    tt = str(t)
    delta_delta = f"({dp}-{df})"
    delta_delta_Z = f"({dpz}-{dfz})"
    wP = f"exp(-({dp})*({tt}))"
    wF = f"exp(-({df})*({tt}))"
    wP_Z = f"-({tt})*({dpz})*{wP}"
    wF_Z = f"-({tt})*({dfz})*{wF}"
    Delta_w = f"{wF}*expm1(-({tt})*{delta_delta})"
    Delta_w_Z = (
        f"({wF_Z})*expm1(-({tt})*{delta_delta}) + "
        f"{wF}*exp(-({tt})*{delta_delta})*(-({tt})*{delta_delta_Z})"
    )
    return {
        "wP": wP,
        "wF": wF,
        "Delta_delta": delta_delta,
        "Delta_delta_Z": delta_delta_Z,
        "Delta_w": Delta_w,
        "wP_Z": wP_Z,
        "wF_Z": wF_Z,
        "Delta_w_Z": Delta_w_Z,
    }


@dataclass(frozen=True)
class NoncommonTimeWeightZ1:
    """Formal non-common time-weight values and their ordinary-Z1 jet."""

    delta_P: Any
    delta_F: Any
    delta_P_Z: Any
    delta_F_Z: Any
    t: Any
    Delta_delta: Any
    Delta_delta_Z: Any
    wP: Any
    wF: Any
    Delta_w: Any
    wP_Z: Any
    wF_Z: Any
    Delta_w_Z: Any
    DeltaH: Any
    H_F: Any
    DeltaH_Z: Any
    H_F_Z: Any
    Delta_wH: Any
    Delta_wH_Z: Any
    product_Z_routes: Mapping[str, Any]
    formal: Mapping[str, str]
    materialization: Mapping[str, Any]

    @property
    def values(self) -> dict[str, Any]:
        return {
            "delta_P": self.delta_P,
            "delta_F": self.delta_F,
            "delta_P_Z": self.delta_P_Z,
            "delta_F_Z": self.delta_F_Z,
            "t": self.t,
            "Delta_delta": self.Delta_delta,
            "Delta_delta_Z": self.Delta_delta_Z,
            "wP": self.wP,
            "wF": self.wF,
            "Delta_w": self.Delta_w,
            "wP_Z": self.wP_Z,
            "wF_Z": self.wF_Z,
            "Delta_w_Z": self.Delta_w_Z,
            "DeltaH": self.DeltaH,
            "H_F": self.H_F,
            "DeltaH_Z": self.DeltaH_Z,
            "H_F_Z": self.H_F_Z,
            "Delta_wH": self.Delta_wH,
            "Delta_wH_Z": self.Delta_wH_Z,
        }

    @property
    def Z1(self) -> Any:
        return self.Delta_w_Z

    @property
    def weight_Z1(self) -> Any:
        return self.Delta_w_Z

    @property
    def product_Z1(self) -> Any:
        return self.Delta_wH_Z

    @property
    def w_P(self) -> Any:
        return self.wP

    @property
    def w_F(self) -> Any:
        return self.wF

    def as_dict(self) -> dict[str, Any]:
        return {
            "values": _display(self.values),
            "formal": dict(self.formal),
            "product_Z_routes": _display(self.product_Z_routes),
            "materialization": _display(self.materialization),
            **SOURCE_VALUE_FLAGS,
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "delta_delta": "Delta_delta",
            "delta_delta_Z": "Delta_delta_Z",
            "delta_w": "Delta_w",
            "delta_w_Z": "Delta_w_Z",
            "wP_z": "wP_Z",
            "wF_z": "wF_Z",
            "delta_wH": "Delta_wH",
            "delta_wH_Z": "Delta_wH_Z",
            "product_routes": "product_Z_routes",
            "w_P": "wP",
            "w_F": "wF",
            "w_P_Z": "wP_Z",
            "w_F_Z": "wF_Z",
        }
        if key == "values":
            return self.values
        if key == "formal":
            return self.formal
        if key in {"Z1", "weight_Z1"}:
            return self.Delta_w_Z
        if key == "product_Z1":
            return self.Delta_wH_Z
        return getattr(self, aliases.get(key, key))


def noncommon_time_weight_z1(
    delta_P: Any,
    delta_F: Any,
    *,
    delta_P_Z: Any = None,
    delta_F_Z: Any = None,
    t: Any = "t",
    DeltaH: Any = None,
    H_F: Any = None,
    DeltaH_Z: Any = None,
    H_F_Z: Any = None,
    exp_fn: Any = None,
) -> NoncommonTimeWeightZ1:
    """Return non-common time weights and complete ordinary-Z1 routes.

    ``delta_P_Z`` and ``delta_F_Z`` are explicit required inputs; omitted
    derivatives are rejected, while zero is valid.  With the default
    ``exp_fn=None`` all weights remain unevaluated formal ``exp``/``expm1``
    expressions.  A caller may provide an exact symbolic exponential such
    as ``sympy.exp``; this is still caller-owned algebra and never a source
    evaluation.  The function does not accept a caller-supplied ``Delta_w``.

    If any of ``DeltaH``, ``H_F``, ``DeltaH_Z`` or ``H_F_Z`` is supplied, all
    four are required and the product
    ``Delta(wH)=wP*DeltaH+Delta_w*H_F`` is returned with its four Z1 routes.
    """

    if delta_P_Z is None or delta_F_Z is None:
        raise GammaNoncommonInputError(
            "delta_P_Z and delta_F_Z are required explicitly"
        )
    product_inputs = (DeltaH, H_F, DeltaH_Z, H_F_Z)
    if any(item is not None for item in product_inputs) and any(
        item is None for item in product_inputs
    ):
        raise GammaNoncommonInputError(
            "DeltaH, H_F, DeltaH_Z, and H_F_Z are required together"
        )

    Delta_delta = delta_P - delta_F
    Delta_delta_Z = delta_P_Z - delta_F_Z
    formal = _time_formal_fields(
        delta_P,
        delta_F,
        delta_P_Z,
        delta_F_Z,
        t,
    )

    if exp_fn is None:
        wP = formal["wP"]
        wF = formal["wF"]
        Delta_w = formal["Delta_w"]
        wP_Z = formal["wP_Z"]
        wF_Z = formal["wF_Z"]
        Delta_w_Z = formal["Delta_w_Z"]
    else:
        wP = exp_fn(-delta_P * t)
        wF = exp_fn(-delta_F * t)
        exp_delta = exp_fn(-t * Delta_delta)
        Delta_w = wF * (exp_delta - 1)
        wP_Z = -t * delta_P_Z * wP
        wF_Z = -t * delta_F_Z * wF
        Delta_w_Z = (
            wF_Z * (exp_delta - 1)
            + wF * exp_delta * (-t * Delta_delta_Z)
        )

    Delta_wH = None
    Delta_wH_Z = None
    product_Z_routes: Mapping[str, Any] = {}
    if all(item is not None for item in product_inputs):
        if exp_fn is None:
            Delta_wH = f"({wP})*({DeltaH}) + ({Delta_w})*({H_F})"
            product_Z_routes = {
                "wP_Z_times_DeltaH": f"({wP_Z})*({DeltaH})",
                "wP_times_DeltaH_Z": f"({wP})*({DeltaH_Z})",
                "Delta_w_Z_times_H_F": f"({Delta_w_Z})*({H_F})",
                "Delta_w_times_H_F_Z": f"({Delta_w})*({H_F_Z})",
            }
            Delta_wH_Z = " + ".join(product_Z_routes.values())
        else:
            Delta_wH = wP * DeltaH + Delta_w * H_F
            product_Z_routes = {
                "wP_Z_times_DeltaH": wP_Z * DeltaH,
                "wP_times_DeltaH_Z": wP * DeltaH_Z,
                "Delta_w_Z_times_H_F": Delta_w_Z * H_F,
                "Delta_w_times_H_F_Z": Delta_w * H_F_Z,
            }
            Delta_wH_Z = sum(product_Z_routes.values(), 0)

    materialization = {
        "delta_P_and_delta_F_explicit": True,
        "delta_P_Z_and_delta_F_Z_explicit": True,
        "both_delta_Z_nonzero_not_assumed": True,
        "Delta_delta_retained": True,
        "Delta_delta_Z_retained": True,
        "wP_wF_Delta_w_retained": True,
        "expm1_formal_representation_retained": True,
        "numeric_expm1_used": False,
        "Delta_w_Z_uses_exp_not_exprel": True,
        "exp_evaluated": exp_fn is not None,
        "caller_weight_difference_accepted": False,
        "caller_weight_difference_proved": False,
        "product_four_Z1_routes_retained": bool(product_Z_routes),
        "source_values_computed": False,
        "source_value_computation": False,
        "historical_source_values_computed": False,
        "certified": False,
    }
    return NoncommonTimeWeightZ1(
        delta_P,
        delta_F,
        delta_P_Z,
        delta_F_Z,
        t,
        Delta_delta,
        Delta_delta_Z,
        wP,
        wF,
        Delta_w,
        wP_Z,
        wF_Z,
        Delta_w_Z,
        DeltaH,
        H_F,
        DeltaH_Z,
        H_F_Z,
        Delta_wH,
        Delta_wH_Z,
        product_Z_routes,
        formal,
        materialization,
    )


def ordinary_z1_difference(
    SP: Any,
    SF: Any,
    DP: Any,
    DF: Any,
    aP: Any,
    aF: Any,
    preP: Any,
    preF: Any,
    *,
    SP_Z: Any = None,
    SF_Z: Any = None,
    DP_Z: Any = None,
    DF_Z: Any = None,
    aP_Z: Any = None,
    aF_Z: Any = None,
    preP_Z: Any = None,
    preF_Z: Any = None,
    DeltaS: Any = None,
    DeltaD: Any = None,
    DeltaS_Z: Any = None,
    DeltaD_Z: Any = None,
    source_kind: str = "synthetic_point",
    rounded: bool = False,
) -> NoncommonGammaZ1:
    """Return branch values and the complete ordinary-Z1 difference.

    All eight P/F input derivatives are required, either as explicit ``*_Z``
    keywords or inside a ``Z1Jet``/two-entry jet.  ``DeltaS_Z`` and
    ``DeltaD_Z`` are optional independent declarations; when omitted they are
    formed from the separately supplied branch derivatives.  The returned
    ``total_Z_factorized`` is a conditional route using those declarations;
    its consistency gate remains unpassed until the caller proves the
    declared differences equal the branch differences.
    """

    sp0, spz = _jet_pair(SP, SP_Z, "SP")
    sf0, sfz = _jet_pair(SF, SF_Z, "SF")
    dp0, dpz = _jet_pair(DP, DP_Z, "DP")
    df0, dfz = _jet_pair(DF, DF_Z, "DF")
    ap0, apz = _jet_pair(aP, aP_Z, "aP")
    af0, afz = _jet_pair(aF, aF_Z, "aF")
    pp0, ppz = _jet_pair(preP, preP_Z, "preP")
    pf0, pfz = _jet_pair(preF, preF_Z, "preF")

    ds_decl, dsz_decl, ds_supplied = _optional_delta_pair(
        DeltaS, DeltaS_Z, "DeltaS"
    )
    dd_decl, ddz_decl, dd_supplied = _optional_delta_pair(
        DeltaD, DeltaD_Z, "DeltaD"
    )
    ds = sp0 - sf0 if not ds_supplied else ds_decl
    dd = dp0 - df0 if not dd_supplied else dd_decl
    dsz = spz - sfz if dsz_decl is None else dsz_decl
    ddz = dpz - dfz if ddz_decl is None else ddz_decl

    base = point_difference(
        sp0,
        sf0,
        dp0,
        df0,
        ap0,
        af0,
        pp0,
        pf0,
        DeltaS=ds,
        DeltaD=dd,
        source_kind=source_kind,
        rounded=rounded,
        _deltaS_declared=ds_supplied,
        _deltaD_declared=dd_supplied,
    )

    XPz = spz * dp0 + sp0 * dpz
    XFz = sfz * df0 + sf0 * dfz
    dXz = XPz - XFz
    dXz_factorized = (
        sfz * dd
        + sf0 * ddz
        + dfz * ds
        + df0 * dsz
        + dsz * dd
        + ds * ddz
    )
    daz = apz - afz
    dprez = ppz - pfz
    qPz = (
        2 * ppz * base.XP
        + 2 * pp0 * XPz
        - apz * base.XP**2
        - 2 * ap0 * base.XP * XPz
    )
    qFz = (
        2 * pfz * base.XF
        + 2 * pf0 * XFz
        - afz * base.XF**2
        - 2 * af0 * base.XF * XFz
    )

    bracket = 2 * base.XF * base.DeltaX + base.DeltaX**2
    bracket_z = (
        2 * XFz * base.DeltaX
        + 2 * base.XF * dXz
        + 2 * base.DeltaX * dXz
    )
    z_terms = {
        "two_preP_Z_times_DeltaX": 2 * ppz * base.DeltaX,
        "two_preP_times_DeltaX_Z": 2 * pp0 * dXz,
        "two_DeltaPre_Z_times_XF": 2 * dprez * base.XF,
        "two_DeltaPre_times_XF_Z": 2 * base.DeltaPre * XFz,
        "minus_aP_Z_times_bracket": -apz * bracket,
        "minus_aP_times_bracket_Z": -ap0 * bracket_z,
        "minus_DeltaA_Z_times_XF_squared": -daz * base.XF**2,
        "minus_DeltaA_times_two_XF_XF_Z": -base.DeltaA * (2 * base.XF * XFz),
        "DeltaX_Z_factorized_routes": dXz_factorized,
    }
    total_Z = qPz - qFz
    retained_total_Z = (
        z_terms["two_preP_Z_times_DeltaX"]
        + z_terms["two_preP_times_DeltaX_Z"]
        + z_terms["two_DeltaPre_Z_times_XF"]
        + z_terms["two_DeltaPre_times_XF_Z"]
        + z_terms["minus_aP_Z_times_bracket"]
        + z_terms["minus_aP_times_bracket_Z"]
        + z_terms["minus_DeltaA_Z_times_XF_squared"]
        + z_terms["minus_DeltaA_times_two_XF_XF_Z"]
    )
    total_Z_factorized, factorized_Z_terms = _difference_z1_from_delta_x(
        base.DeltaX_factorized,
        dXz_factorized,
        base.XF,
        XFz,
        ap0,
        apz,
        base.DeltaA,
        daz,
        pp0,
        ppz,
        base.DeltaPre,
        dprez,
    )
    factorized_Z_terms = {
        **factorized_Z_terms,
        "DeltaX_Z_factorized_routes": dXz_factorized,
    }
    direct_factorized_Z1_residual = total_Z - total_Z_factorized
    difference_consistency_gate = {
        "required_for_factorized_identity": True,
        "passed": False,
        "consistency_proved": False,
        "status": "unpassed_caller_gate",
        "caller_must_prove": True,
        "declarations_are_not_evidence": True,
        "branch_delta_equality_is_not_checked_here": True,
        "direct_minus_factorized_value": base.direct_factorized_residual,
        "direct_minus_factorized_Z1": direct_factorized_Z1_residual,
        "value_gate": dict(base.difference_consistency_gate),
    }
    flags = dict(base.materialization)
    flags.update(
        {
            "ordinary_Z1_complete": True,
            "all_P_F_value_jets_required": True,
            "DeltaS_Z_retained": True,
            "DeltaD_Z_retained": True,
            "DeltaS_Z_declared_independently": dsz_decl is not None,
            "DeltaD_Z_declared_independently": ddz_decl is not None,
            "DeltaA_Z_retained": True,
            "DeltaPre_Z_retained": True,
            "DeltaX_Z_direct_from_branch_products": True,
            "DeltaX_Z_factorized_routes_retained": True,
            "conditional_factorized_Z1_returned": True,
            "difference_consistency_gate_passed": False,
            "difference_consistency_caller_must_prove": True,
            "retained_Z_sum_matches_direct_branch_difference": (
                retained_total_Z == total_Z
            ),
            "common_a_pre_D_S_assumed": False,
        }
    )
    return NoncommonGammaZ1(
        base,
        spz,
        sfz,
        dpz,
        dfz,
        apz,
        afz,
        ppz,
        pfz,
        dsz,
        ddz,
        daz,
        dprez,
        XPz,
        XFz,
        dXz,
        dXz_factorized,
        qPz,
        qFz,
        total_Z,
        total_Z_factorized,
        direct_factorized_Z1_residual,
        {
            **z_terms,
            "retained_total_Z": retained_total_Z,
            "direct_qP_minus_qF_Z": total_Z,
        },
        factorized_Z_terms,
        difference_consistency_gate,
        flags,
    )


def _fixture_exact_sympy() -> dict[str, Any]:
    """Build an independent exact SymPy fixture; this function is uncalled."""

    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - optional diagnostic dependency
        return {
            "fixture": "independent_exact_sympy_noncommon_z1",
            "checks_passed": False,
            "sympy_available": False,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            **SOURCE_VALUE_FLAGS,
        }

    q = sp.Rational
    z = sp.symbols("z")

    def integrand(S: Any, D: Any, a: Any, pre: Any) -> Any:
        X = S * D
        return X * (2 * pre - a * X)

    vals = {
        "SP": q(13, 17),
        "SF": q(5, 11),
        "DP": q(7, 13),
        "DF": q(2, 9),
        "aP": q(11, 19),
        "aF": q(3, 8),
        "preP": q(17, 12),
        "preF": q(9, 10),
    }
    jets = {
        "SP_Z": q(2, 23),
        "SF_Z": -q(3, 29),
        "DP_Z": q(5, 31),
        "DF_Z": q(7, 37),
        "aP_Z": -q(11, 41),
        "aF_Z": q(13, 43),
        "preP_Z": q(17, 47),
        "preF_Z": -q(19, 53),
    }
    ds = vals["SP"] - vals["SF"]
    dd = vals["DP"] - vals["DF"]
    dsz = jets["SP_Z"] - jets["SF_Z"]
    ddz = jets["DP_Z"] - jets["DF_Z"]
    point = point_difference(**vals, DeltaS=ds, DeltaD=dd)
    direct_value = sp.expand(
        integrand(vals["SP"], vals["DP"], vals["aP"], vals["preP"])
        - integrand(vals["SF"], vals["DF"], vals["aF"], vals["preF"])
    )
    affine = lambda name: vals[name] + jets[f"{name}_Z"] * z
    z1 = ordinary_z1_difference(
        **vals,
        **jets,
        DeltaS=ds,
        DeltaD=dd,
        DeltaS_Z=dsz,
        DeltaD_Z=ddz,
    )
    z1_derived = ordinary_z1_difference(**vals, **jets)
    direct_z = sp.diff(
        integrand(affine("SP"), affine("DP"), affine("aP"), affine("preP"))
        - integrand(affine("SF"), affine("DF"), affine("aF"), affine("preF")),
        z,
    ).subs(z, 0)

    zero = point_difference(
        vals["SF"],
        vals["SF"],
        vals["DF"],
        vals["DF"],
        vals["aP"],
        vals["aF"],
        vals["preP"],
        vals["preF"],
    )
    zero_expected = (
        2 * (vals["preP"] - vals["preF"]) * vals["SF"] * vals["DF"]
        - (vals["aP"] - vals["aF"]) * (vals["SF"] * vals["DF"]) ** 2
    )
    eps = q(1, 10**40)
    near = point_difference(
        vals["SF"] + eps,
        vals["SF"],
        vals["DF"] - eps,
        vals["DF"],
        vals["aP"],
        vals["aF"],
        vals["preP"],
        vals["preF"],
    )
    near_direct = sp.expand(
        integrand(vals["SF"] + eps, vals["DF"] - eps, vals["aP"], vals["preP"])
        - integrand(vals["SF"], vals["DF"], vals["aF"], vals["preF"])
    )
    inconsistent = point_difference(
        **vals,
        DeltaS=ds + eps,
        DeltaD=dd,
    )
    inconsistent_z1 = ordinary_z1_difference(
        **vals,
        **jets,
        DeltaS=ds + eps,
        DeltaD=dd,
        DeltaS_Z=dsz + eps,
        DeltaD_Z=ddz,
    )
    only_dsz_inconsistent = ordinary_z1_difference(
        **vals,
        **jets,
        DeltaS_Z=dsz + eps,
    )
    only_ddz_inconsistent = ordinary_z1_difference(
        **vals,
        **jets,
        DeltaD_Z=ddz + eps,
    )

    retained_z_names = (
        "two_preP_Z_times_DeltaX",
        "two_preP_times_DeltaX_Z",
        "two_DeltaPre_Z_times_XF",
        "two_DeltaPre_times_XF_Z",
        "minus_aP_Z_times_bracket",
        "minus_aP_times_bracket_Z",
        "minus_DeltaA_Z_times_XF_squared",
        "minus_DeltaA_times_two_XF_XF_Z",
    )
    retained_direct_sum = sum(
        (z1.retained_Z_terms[name] for name in retained_z_names), sp.Integer(0)
    )
    retained_factorized_sum = sum(
        (z1.factorized_Z_terms[name] for name in retained_z_names), sp.Integer(0)
    )

    missing_raises = False
    missing_jets = dict(jets)
    missing_jets.pop("DF_Z")
    try:
        ordinary_z1_difference(**vals, **missing_jets)
    except GammaNoncommonInputError:
        missing_raises = True

    checks = {
        "noncommon_value_exact": sp.simplify(point.value - direct_value) == 0,
        "correct_declared_value_routes_equal": sp.simplify(
            point.value - point.value_factorized
        )
        == 0,
        "deltaS_retained": point.DeltaS == ds,
        "deltaD_retained": point.DeltaD == dd,
        "deltaX_factorization_exact": sp.simplify(
            point.DeltaX - point.DeltaX_factorized
        )
        == 0,
        "ordinary_Z1_exact": sp.simplify(z1.total_Z - direct_z) == 0,
        "correct_declared_Z1_routes_equal": sp.simplify(
            z1.total_Z - z1.total_Z_factorized
        )
        == 0,
        "retained_Z_routes_sum_to_direct": sp.simplify(
            retained_direct_sum - z1.total_Z
        )
        == 0,
        "conditional_Z_routes_sum_to_factorized": sp.simplify(
            retained_factorized_sum - z1.total_Z_factorized
        )
        == 0,
        "derived_delta_flags_not_independent": (
            z1_derived.base.materialization["DeltaS_declared_independently"]
            is False
            and z1_derived.base.materialization["DeltaD_declared_independently"]
            is False
            and z1_derived.materialization["DeltaS_Z_declared_independently"]
            is False
            and z1_derived.materialization["DeltaD_Z_declared_independently"]
            is False
        ),
        "explicit_delta_flags_independent": (
            z1.base.materialization["DeltaS_declared_independently"] is True
            and z1.base.materialization["DeltaD_declared_independently"] is True
            and z1.materialization["DeltaS_Z_declared_independently"] is True
            and z1.materialization["DeltaD_Z_declared_independently"] is True
        ),
        "all_branch_Z1_values_exposed": all(
            value is not None
            for value in (
                z1.qP_Z,
                z1.qF_Z,
                z1.DeltaS_Z,
                z1.DeltaD_Z,
                z1.DeltaA_Z,
                z1.DeltaPre_Z,
            )
        ),
        "zero_differences_exact": sp.simplify(zero.value - zero_expected) == 0
        and zero.DeltaS == 0
        and zero.DeltaD == 0
        and sp.simplify(zero.value - zero.value_factorized) == 0,
        "near_zero_differences_exact": sp.simplify(near.value - near_direct) == 0
        and near.DeltaS == eps
        and near.DeltaD == -eps
        and near.DeltaX_factorized != 0
        and sp.simplify(near.value - near.value_factorized) == 0,
        "inconsistent_value_declaration_not_identity": (
            sp.simplify(inconsistent.value - inconsistent.value_factorized) != 0
            and inconsistent.difference_consistency_gate["passed"] is False
        ),
        "inconsistent_Z1_declaration_not_identity": (
            sp.simplify(
                inconsistent_z1.total_Z - inconsistent_z1.total_Z_factorized
            )
            != 0
            and inconsistent_z1.difference_consistency_gate["passed"] is False
        ),
        "only_DeltaS_Z_inconsistency_retained": (
            sp.simplify(
                only_dsz_inconsistent.total_Z
                - only_dsz_inconsistent.total_Z_factorized
            )
            != 0
            and only_dsz_inconsistent.difference_consistency_gate["passed"]
            is False
        ),
        "only_DeltaD_Z_inconsistency_retained": (
            sp.simplify(
                only_ddz_inconsistent.total_Z
                - only_ddz_inconsistent.total_Z_factorized
            )
            != 0
            and only_ddz_inconsistent.difference_consistency_gate["passed"]
            is False
        ),
        "consistency_gate_never_auto_passes": (
            point.difference_consistency_gate["passed"] is False
            and z1.difference_consistency_gate["passed"] is False
        ),
        "missing_derivative_rejected": missing_raises,
    }
    return {
        "schema": SCHEMA,
        "fixture": "independent_exact_sympy_noncommon_z1",
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "sympy_available": True,
        "exact_arithmetic": True,
        "differing_inputs": {name: str(value) for name, value in vals.items()},
        "nonzero_Z_jets": {name: str(value) for name, value in jets.items()},
        "zero_and_near_zero": {"zero": "DeltaS=DeltaD=0", "near": str(eps)},
        **SOURCE_VALUE_FLAGS,
    }


def _fixture_exact_sympy_time_weight() -> dict[str, Any]:
    """Build an independent exact SymPy time-weight fixture; leave it uncalled."""

    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - optional diagnostic dependency
        return {
            "fixture": "independent_exact_sympy_noncommon_time_weight_z1",
            "checks_passed": False,
            "sympy_available": False,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            **SOURCE_VALUE_FLAGS,
        }

    q = sp.Rational
    t = sp.symbols("t")
    z = sp.symbols("z")
    delta = q(2, 5)
    delta_P_Z = q(3, 7)
    delta_F_Z = -q(5, 11)
    DeltaH = q(4, 9)
    H_F = q(7, 13)
    DeltaH_Z = q(2, 17)
    H_F_Z = -q(3, 19)

    packet = noncommon_time_weight_z1(
        delta,
        delta,
        delta_P_Z=delta_P_Z,
        delta_F_Z=delta_F_Z,
        t=t,
        DeltaH=DeltaH,
        H_F=H_F,
        DeltaH_Z=DeltaH_Z,
        H_F_Z=H_F_Z,
        exp_fn=sp.exp,
    )
    wP = sp.exp(-delta * t)
    wF = sp.exp(-delta * t)
    Delta_delta = sp.Integer(0)
    Delta_delta_Z = delta_P_Z - delta_F_Z
    Delta_w = wF * (sp.exp(-t * Delta_delta) - 1)
    wP_Z = -t * delta_P_Z * wP
    wF_Z = -t * delta_F_Z * wF
    Delta_w_Z = (
        wF_Z * (sp.exp(-t * Delta_delta) - 1)
        + wF * sp.exp(-t * Delta_delta) * (-t * Delta_delta_Z)
    )
    chart_delta_P = delta + delta_P_Z * z
    chart_delta_F = delta + delta_F_Z * z
    chart_DeltaH = DeltaH + DeltaH_Z * z
    chart_H_F = H_F + H_F_Z * z
    chart_wP = sp.exp(-chart_delta_P * t)
    chart_wF = sp.exp(-chart_delta_F * t)
    chart_Delta_w = chart_wF * (
        sp.exp(-t * (chart_delta_P - chart_delta_F)) - 1
    )
    direct_product_Z = sp.diff(
        chart_wP * chart_DeltaH + chart_Delta_w * chart_H_F,
        z,
    ).subs(z, 0)

    expected_routes = {
        "wP_Z_times_DeltaH": wP_Z * DeltaH,
        "wP_times_DeltaH_Z": wP * DeltaH_Z,
        "Delta_w_Z_times_H_F": Delta_w_Z * H_F,
        "Delta_w_times_H_F_Z": Delta_w * H_F_Z,
    }
    missing_derivative_raises = False
    try:
        noncommon_time_weight_z1(delta, delta, delta_P_Z=delta_P_Z, t=t)
    except GammaNoncommonInputError:
        missing_derivative_raises = True

    formal = noncommon_time_weight_z1(
        delta,
        delta,
        delta_P_Z=delta_P_Z,
        delta_F_Z=delta_F_Z,
        t="t",
    )
    unequal_delta_P = q(3, 5)
    unequal = noncommon_time_weight_z1(
        unequal_delta_P, delta, delta_P_Z=delta_P_Z, delta_F_Z=delta_F_Z,
        t=t, DeltaH=DeltaH, H_F=H_F, DeltaH_Z=DeltaH_Z,
        H_F_Z=H_F_Z, exp_fn=sp.exp,
    )
    unequal_direct_Z = sp.diff(
        sp.exp(-(unequal_delta_P + delta_P_Z*z)*t)
        * (chart_H_F + chart_DeltaH)
        - sp.exp(-(delta + delta_F_Z*z)*t) * chart_H_F, z
    ).subs(z, 0)
    checks = {
        "both_delta_Z_nonzero_and_different": (
            delta_P_Z != 0
            and delta_F_Z != 0
            and delta_P_Z != delta_F_Z
        ),
        "Delta_delta_zero": packet.Delta_delta == 0,
        "Delta_delta_Z_nonzero": packet.Delta_delta_Z != 0,
        "wP_exact": sp.simplify(packet.wP - wP) == 0,
        "wF_exact": sp.simplify(packet.wF - wF) == 0,
        "Delta_w_exact": sp.simplify(packet.Delta_w - Delta_w) == 0,
        "wP_Z_exact": sp.simplify(packet.wP_Z - wP_Z) == 0,
        "wF_Z_exact": sp.simplify(packet.wF_Z - wF_Z) == 0,
        "Delta_w_Z_exp_exact": sp.simplify(packet.Delta_w_Z - Delta_w_Z) == 0,
        "four_product_Z1_routes_exact": all(
            sp.simplify(packet.product_Z_routes[name] - expected) == 0
            for name, expected in expected_routes.items()
        ),
        "Delta_wH_Z_exact": sp.simplify(
            packet.Delta_wH_Z - direct_product_Z
        )
        == 0,
        "unequal_delta_Z1_exact": sp.simplify(unequal.Delta_wH_Z - unequal_direct_Z) == 0,
        "unequal_delta_all_four_routes_nonzero": all(
            sp.simplify(route) != 0 for route in unequal.product_Z_routes.values()
        ),
        "Delta_w_Z_uses_exp_not_exprel": (
            "exp(" in packet.formal["Delta_w_Z"]
            and "exprel" not in packet.formal["Delta_w_Z"]
        ),
        "default_form_is_not_evaluated": (
            isinstance(formal.wP, str)
            and isinstance(formal.Delta_w, str)
            and formal.materialization["exp_evaluated"] is False
        ),
        "missing_delta_F_Z_rejected": missing_derivative_raises,
    }
    return {
        "schema": SCHEMA,
        "fixture": "independent_exact_sympy_noncommon_time_weight_z1",
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "sympy_available": True,
        "exact_arithmetic": True,
        "delta_values": {"delta_P": str(delta), "delta_F": str(delta)},
        "delta_Z_values": {
            "delta_P_Z": str(delta_P_Z),
            "delta_F_Z": str(delta_F_Z),
            "Delta_delta_Z": str(Delta_delta_Z),
        },
        **SOURCE_VALUE_FLAGS,
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run the optional synthetic fixture when a caller explicitly requests it."""

    result = _fixture_exact_sympy()
    time_weight = _fixture_exact_sympy_time_weight()
    result["time_weight_fixture"] = time_weight
    result["checks_passed"] = bool(
        result.get("checks_passed") is True
        and time_weight.get("checks_passed") is True
    )
    return result


fixture_exact_sympy = _fixture_exact_sympy
fixture_exact_sympy_time_weight = _fixture_exact_sympy_time_weight


def verify_time_weight_fixture() -> dict[str, Any]:
    """Run only the optional synthetic time-weight fixture on explicit request."""

    return _fixture_exact_sympy_time_weight()


# Descriptive aliases keep the source notation discoverable without importing
# the preceding common-source leaf.
general_noncommon_difference = point_difference
general_source_difference = point_difference
noncommon_z1_difference = ordinary_z1_difference
general_noncommon_z1_difference = ordinary_z1_difference
general_noncommon_z1 = ordinary_z1_difference
ordinary_Z1_difference = ordinary_z1_difference
build_noncommon_difference = point_difference
build_noncommon_z1 = ordinary_z1_difference
GammaNoncommonDifference = NoncommonGammaDifference
GammaNoncommonZ1 = NoncommonGammaZ1
time_weight_difference_z1 = noncommon_time_weight_z1
time_weight_z1 = noncommon_time_weight_z1
build_time_weight_z1 = noncommon_time_weight_z1
GammaNoncommonTimeWeightZ1 = NoncommonTimeWeightZ1


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "SOURCE_VALUE_FLAGS",
    "FORMULA",
    "DIFFERENCE_FORMULA",
    "Z1_FORMULA",
    "GammaNoncommonInputError",
    "Z1Jet",
    "NoncommonGammaDifference",
    "NoncommonGammaZ1",
    "NoncommonTimeWeightZ1",
    "GammaNoncommonDifference",
    "GammaNoncommonZ1",
    "GammaNoncommonTimeWeightZ1",
    "point_difference",
    "general_noncommon_difference",
    "general_source_difference",
    "ordinary_z1_difference",
    "noncommon_z1_difference",
    "general_noncommon_z1_difference",
    "general_noncommon_z1",
    "ordinary_Z1_difference",
    "build_noncommon_difference",
    "build_noncommon_z1",
    "noncommon_time_weight_z1",
    "time_weight_difference_z1",
    "time_weight_z1",
    "build_time_weight_z1",
    "fixture_exact_sympy",
    "fixture_exact_sympy_time_weight",
    "verify_time_weight_fixture",
    "verify_synthetic_fixture",
]
