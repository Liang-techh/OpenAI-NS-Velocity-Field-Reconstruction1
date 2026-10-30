"""Conditional P-minus-F difference for a non-common quadratic selector.

The selector on branch ``b`` is the source-shaped relation

    A2_b*a_b**2 + A1_b*a_b + A0_b = 0.

For separate P and F branches the cancellation-stable, source-oriented
identity used here is

    Delta_a = -(DeltaA2*aF**2 + DeltaA1*aF + DeltaA0) / D,
    D = A2P*(aP + aF) + A1P.

``DeltaA2``, ``DeltaA1`` and ``DeltaA0`` are required caller-owned signed
inputs; their ordinary-Z derivatives are required for Z1. No branch
subtraction is substituted for a missing signed difference. No common-source, midpoint,
interval-endpoint, precision, rounding, or source-value claim is made here.

The optional SymPy fixture is independent and exact.  It is defined but never
called at import time.  It is a synthetic algebra check only; it does not
bind Q/P/F owners, prove source identity, certify a positive root, or certify
the denominator or any numerical error.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Mapping, Sequence


SCHEMA = "pf_selector_correlated_root_v1"
COMPONENT = SCHEMA

FORMULA = (
    "Delta_a=-(DeltaA2*aF**2+DeltaA1*aF+DeltaA0)/"
    "(A2P*(aP+aF)+A1P)"
)
Z1_FORMULA = "ordinary-Z1 chain rule for the signed numerator/denominator quotient"

SOURCE_VALUE_FLAGS = {
    "source_values_computed": False,
    "source_value_computation": False,
    "historical_source_values_computed": False,
    "actual_source_value_computed": False,
    "source_identity_proved": False,
    "certified": False,
    "numerical_precision_certified": False,
    "rounding_certified": False,
    "nonzero_certified": False,
}


class PFSelectorInputError(ValueError):
    """Raised when a selector value or ordinary-Z1 input is incomplete."""


class OpaqueSourceReference:
    """Marker for a source-owned opaque id, never an arithmetic scalar.

    A naked integer cannot be distinguished from a numeric integer at this
    layer.  Callers must therefore wrap saved source handles in this marker,
    or provide an equivalent mapping with ``opaque``/``source_owned`` set.
    """

    __slots__ = ("identifier", "owner")

    def __init__(self, identifier: Any, owner: Any = None) -> None:
        self.identifier = identifier
        self.owner = owner

    def __repr__(self) -> str:  # pragma: no cover - diagnostic only
        return f"OpaqueSourceReference({self.identifier!r}, owner={self.owner!r})"


SourceOpaqueId = OpaqueSourceReference


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
    """Serialize diagnostics without coercing caller-owned arithmetic."""

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


def _is_opaque(value: Any) -> bool:
    """Recognize explicit opaque source references without rejecting integers."""

    if isinstance(value, OpaqueSourceReference):
        return True
    if isinstance(value, Mapping):
        if any(key in value for key in ("source_node_id", "node_id", "source_role", "namespace", "owner", "kind")):
            return True
        for key in (
            "opaque",
            "is_opaque",
            "source_owned",
            "opaque_id",
            "node_id_is_not_a_numeric_value",
            "reference_is_opaque_metadata",
            "roles_are_opaque_metadata",
        ):
            if value.get(key) is True:
                return True
        if value.get("kind") in {"opaque_saved_role", "saved_opaque_ref"}:
            return True
    for attribute in ("is_opaque", "opaque", "source_owned"):
        try:
            if getattr(value, attribute) is True:
                return True
        except Exception:
            pass
    return False


def _reject_opaque(value: Any, label: str) -> None:
    if _is_opaque(value):
        raise PFSelectorInputError(
            f"{label} is a source-owned opaque reference, not an arithmetic value"
        )


def _jet_pair(value: Any, derivative: Any, label: str) -> tuple[Any, Any]:
    """Resolve one required scalar/ordinary-Z1 jet without inventing Z."""

    _reject_opaque(value, label)
    _reject_opaque(derivative, label + "_Z")
    if derivative is not None:
        value0 = value.value if isinstance(value, Z1Jet) else value
        if isinstance(value, Mapping):
            if "value" in value:
                value0 = value["value"]
            elif "base" in value:
                value0 = value["base"]
        _reject_opaque(value0, label)
        _reject_opaque(derivative, label + "_Z")
        return value0, derivative
    if isinstance(value, Z1Jet):
        _reject_opaque(value.value, label)
        _reject_opaque(value.Z, label + "_Z")
        return value.value, value.Z
    if isinstance(value, Mapping):
        for value_key, derivative_key in (
            ("value", "Z"),
            ("value", "z"),
            ("base", "derivative"),
            ("0", "1"),
        ):
            if value_key in value and derivative_key in value:
                value0, derivative0 = value[value_key], value[derivative_key]
                _reject_opaque(value0, label)
                _reject_opaque(derivative0, label + "_Z")
                return value0, derivative0
        raise PFSelectorInputError(f"{label} mapping must contain value and Z entries")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            _reject_opaque(value[0], label)
            _reject_opaque(value[1], label + "_Z")
            return value[0], value[1]
        raise PFSelectorInputError(f"{label} sequence jet must have exactly two entries")
    raise PFSelectorInputError(f"{label} requires an explicit ordinary-Z derivative")


def _value_only(value: Any, label: str) -> Any:
    """Resolve a point value while refusing explicit opaque source handles."""

    _reject_opaque(value, label)
    if isinstance(value, Z1Jet):
        value = value.value
    elif isinstance(value, Mapping):
        if "value" in value:
            value = value["value"]
        elif "base" in value:
            value = value["base"]
    _reject_opaque(value, label)
    return value


def _optional_delta(
    value: Any,
    derivative: Any,
    branch_value: Any,
    branch_derivative: Any | None,
    label: str,
) -> tuple[Any | None, Any | None, bool, bool]:
    """Resolve an independently declared delta and report both declarations."""

    _reject_opaque(value, label)
    _reject_opaque(derivative, label + "_Z")
    supplied_value = value is not None
    supplied_derivative = derivative is not None
    if isinstance(value, Z1Jet):
        value0 = value.value
        if derivative is None:
            derivative = value.Z
            supplied_derivative = True
    elif isinstance(value, Mapping):
        if "value" in value:
            value0 = value["value"]
        elif "base" in value:
            value0 = value["base"]
        else:
            raise PFSelectorInputError(f"{label} mapping must contain value or base")
        if derivative is None and "Z" in value:
            derivative = value["Z"]
            supplied_derivative = True
    else:
        value0 = value
    if supplied_value:
        _reject_opaque(value0, label)
    if supplied_derivative:
        _reject_opaque(derivative, label + "_Z")
    if not supplied_value:
        raise PFSelectorInputError(f"{label} requires an explicit signed difference")
    if branch_derivative is not None and not supplied_derivative:
        raise PFSelectorInputError(f"{label}_Z requires an explicit signed difference derivative")
    return value0, derivative, supplied_value, supplied_derivative


def _definitely_zero(value: Any) -> bool:
    """Return True only for a value that is unambiguously equal to zero."""

    if value is None:
        return False
    try:
        answer = value == 0
    except Exception:
        return False
    return answer is True or type(answer).__name__ == "BooleanTrue"


def _safe_divide(numerator: Any, denominator: Any, label: str) -> Any:
    if _definitely_zero(denominator):
        raise PFSelectorInputError(f"{label} is zero; the signed root relation is undefined")
    return numerator / denominator


def _materialization_flags(
    *,
    source_kind: str,
    rounded: bool,
    delta_values_independent: Mapping[str, bool],
    delta_derivatives_independent: Mapping[str, bool],
) -> dict[str, Any]:
    return {
        **SOURCE_VALUE_FLAGS,
        "source_kind": source_kind,
        "source_kind_is_caller_declared": True,
        "synthetic_point": source_kind == "synthetic_point",
        "formal_symbolic": source_kind == "formal_symbolic",
        "uncertified_rounding": bool(rounded),
        "rounded_input_declared": bool(rounded),
        "parameter_interval_midpoint_used": False,
        "P_and_F_common_input_assumed": False,
        "all_P_F_coefficients_explicit": True,
        "all_P_F_root_values_explicit": True,
        "DeltaA_values_caller_owned": True,
        "DeltaA_Z_caller_owned": True,
        "DeltaA_values_independent": dict(delta_values_independent),
        "DeltaA_Z_independent": dict(delta_derivatives_independent),
        "signed_three_coefficient_terms_retained": True,
        "ordinary_Z1_chain_terms_retained": False,
        "conditional_formula_returned": True,
    }


def _gate(required: str, *, detail: Any = None) -> dict[str, Any]:
    result = {
        "passed": False,
        "status": "unpassed_caller_gate",
        "required_for_claim": required,
        "caller_must_prove": True,
        "proof_claimed": False,
    }
    if detail is not None:
        result["detail"] = detail
    return result


@dataclass(frozen=True)
class NoncommonSelectorDifference:
    """Point selector difference with signed numerator/denominator routes."""

    A2P: Any
    A2F: Any
    A1P: Any
    A1F: Any
    A0P: Any
    A0F: Any
    aP: Any
    aF: Any
    DeltaA2: Any
    DeltaA1: Any
    DeltaA0: Any
    numerator: Any
    denominator: Any
    Delta_a: Any
    numerator_terms: Mapping[str, Any]
    denominator_terms: Mapping[str, Any]
    branch_equation_residuals: Mapping[str, Any]
    difference_residuals: Mapping[str, Any]
    gates: Mapping[str, Any]
    materialization: Mapping[str, Any]

    @property
    def value(self) -> Any:
        return self.Delta_a

    @property
    def delta_a(self) -> Any:
        return self.Delta_a

    @property
    def da(self) -> Any:
        return self.Delta_a

    @property
    def signed_terms(self) -> Mapping[str, Any]:
        return self.numerator_terms

    @property
    def conditional_formula(self) -> Any:
        return self.Delta_a

    def as_dict(self) -> dict[str, Any]:
        values = {
            "A2P": self.A2P,
            "A2F": self.A2F,
            "A1P": self.A1P,
            "A1F": self.A1F,
            "A0P": self.A0P,
            "A0F": self.A0F,
            "aP": self.aP,
            "aF": self.aF,
            "DeltaA2": self.DeltaA2,
            "DeltaA1": self.DeltaA1,
            "DeltaA0": self.DeltaA0,
            "numerator": self.numerator,
            "denominator": self.denominator,
            "Delta_a": self.Delta_a,
        }
        return {
            "schema": SCHEMA,
            "values": {key: _display(value) for key, value in values.items()},
            "signed_numerator_terms": {
                key: _display(value) for key, value in self.numerator_terms.items()
            },
            "denominator_terms": {
                key: _display(value) for key, value in self.denominator_terms.items()
            },
            "branch_equation_residuals": {
                key: _display(value)
                for key, value in self.branch_equation_residuals.items()
            },
            "difference_residuals": {
                key: _display(value) for key, value in self.difference_residuals.items()
            },
            "gates": _display(self.gates),
            "formula": FORMULA,
            "materialization": dict(self.materialization),
            **SOURCE_VALUE_FLAGS,
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "delta_a": "Delta_a",
            "Da": "Delta_a",
            "dA2": "DeltaA2",
            "dA1": "DeltaA1",
            "dA0": "DeltaA0",
            "N": "numerator",
            "D": "denominator",
            "factorized": "Delta_a",
            "conditional": "Delta_a",
        }
        if key == "values":
            return {
                name: getattr(self, name)
                for name in (
                    "A2P",
                    "A2F",
                    "A1P",
                    "A1F",
                    "A0P",
                    "A0F",
                    "aP",
                    "aF",
                    "DeltaA2",
                    "DeltaA1",
                    "DeltaA0",
                    "numerator",
                    "denominator",
                    "Delta_a",
                )
            }
        if key in ("gates", "gate"):
            return self.gates
        if key in ("terms", "signed_terms"):
            return self.numerator_terms
        return getattr(self, aliases.get(key, key))


def point_difference(
    A2P: Any,
    A2F: Any,
    A1P: Any,
    A1F: Any,
    A0P: Any,
    A0F: Any,
    aP: Any,
    aF: Any,
    *,
    DeltaA2: Any = None,
    DeltaA1: Any = None,
    DeltaA0: Any = None,
    source_kind: str = "synthetic_point",
    rounded: bool = False,
    q_owner: Any = None,
    p_owner: Any = None,
    f_owner: Any = None,
    derivative_semantics: Any = None,
) -> NoncommonSelectorDifference:
    """Build the conditional signed P-minus-F selector root difference.

    All three delta coefficients must be supplied explicitly. The identity
    gate remains false and no source correlation is inferred.
    """

    A2P0, A2F0, A1P0, A1F0, A0P0, A0F0 = (
        _value_only(A2P, "A2P"),
        _value_only(A2F, "A2F"),
        _value_only(A1P, "A1P"),
        _value_only(A1F, "A1F"),
        _value_only(A0P, "A0P"),
        _value_only(A0F, "A0F"),
    )
    aP0, aF0 = _value_only(aP, "aP"), _value_only(aF, "aF")

    dA2, _, dA2_ind, _ = _optional_delta(
        DeltaA2, None, A2P0 - A2F0, None, "DeltaA2"
    )
    dA1, _, dA1_ind, _ = _optional_delta(
        DeltaA1, None, A1P0 - A1F0, None, "DeltaA1"
    )
    dA0, _, dA0_ind, _ = _optional_delta(
        DeltaA0, None, A0P0 - A0F0, None, "DeltaA0"
    )

    numerator_terms = {
        "DeltaA2_times_aF_squared": dA2 * aF0**2,
        "DeltaA1_times_aF": dA1 * aF0,
        "DeltaA0": dA0,
    }
    numerator = sum(numerator_terms.values(), 0)
    denominator_terms = {
        "A2P_times_aP_plus_aF": A2P0 * (aP0 + aF0),
        "A1P": A1P0,
    }
    denominator = sum(denominator_terms.values(), 0)
    if _definitely_zero(denominator):
        raise PFSelectorInputError(
            "A2P*(aP+aF)+A1P is zero; nonzero denominator is required"
        )
    Delta_a = -numerator / denominator

    branch_residuals = {
        "P": A2P0 * aP0**2 + A1P0 * aP0 + A0P0,
        "F": A2F0 * aF0**2 + A1F0 * aF0 + A0F0,
    }
    difference_residuals = {
        "DeltaA2_minus_branch_difference": dA2 - (A2P0 - A2F0),
        "DeltaA1_minus_branch_difference": dA1 - (A1P0 - A1F0),
        "DeltaA0_minus_branch_difference": dA0 - (A0P0 - A0F0),
    }
    delta_values_independent = {
        "DeltaA2": dA2_ind,
        "DeltaA1": dA1_ind,
        "DeltaA0": dA0_ind,
    }
    delta_derivatives_independent = {
        "DeltaA2_Z": False,
        "DeltaA1_Z": False,
        "DeltaA0_Z": False,
    }
    gates = {
        "denominator_nonzero": _gate(
            "caller proof that A2P*(aP+aF)+A1P is nonzero"
        ),
        "positive_root_P": _gate("caller proof that aP is the selected positive root"),
        "positive_root_F": _gate("caller proof that aF is the selected positive root"),
        "A2_positive_both_branches": _gate("source-bound A2P and A2F strict positivity"),
        "A0_negative_both_branches": _gate("source-bound A0P and A0F strict negativity"),
        "discriminants_positive": _gate("source-bound strict positive discriminants"),
        "cross_denominator_positive_margin": _gate("certified positive lower bound for cross denominator"),
        "implicit_derivative_denominators_positive_margin": _gate("certified positive lower bounds for both branch derivatives"),
        "coefficient_difference_consistency": _gate(
            "caller proof that each DeltaAi equals AiP-AiF",
            detail=difference_residuals,
        ),
        "Q_owner": _gate("caller-owned Q/selector owner and source identity"),
        "P_owner": _gate("caller-owned P selector owner and source identity"),
        "F_owner": _gate("caller-owned F selector owner and source identity"),
        "Q_P_F_owner": _gate("caller proof of Q/P/F owner alignment"),
        "derivative_semantics": _gate(
            "caller proof that ordinary-Z derivatives have the claimed source semantics"
        ),
        "source_identity": _gate("caller proof of the P/F source-function identity"),
        "numeric_precision": _gate("caller precision and arithmetic error analysis"),
        "rounding": _gate("caller rounding/interval error analysis"),
        "nonzero": _gate("caller proof of all required nonzero quantities"),
    }
    materialization = _materialization_flags(
        source_kind=source_kind,
        rounded=rounded,
        delta_values_independent=delta_values_independent,
        delta_derivatives_independent=delta_derivatives_independent,
    )
    materialization.update(
        {
            "ordinary_Z1_chain_terms_retained": False,
            "denominator_zero_rejected_if_definite": True,
            "positive_root_checked": False,
            "owner_arguments_recorded_only": {
                "q_owner": q_owner is not None,
                "p_owner": p_owner is not None,
                "f_owner": f_owner is not None,
            },
            "derivative_semantics_argument_recorded_only": derivative_semantics
            is not None,
        }
    )
    return NoncommonSelectorDifference(
        A2P0,
        A2F0,
        A1P0,
        A1F0,
        A0P0,
        A0F0,
        aP0,
        aF0,
        dA2,
        dA1,
        dA0,
        numerator,
        denominator,
        Delta_a,
        numerator_terms,
        denominator_terms,
        branch_residuals,
        difference_residuals,
        gates,
        materialization,
    )


@dataclass(frozen=True)
class NoncommonSelectorZ1:
    """Complete ordinary-Z1 selector difference packet."""

    base: NoncommonSelectorDifference
    A2P_Z: Any
    A2F_Z: Any
    A1P_Z: Any
    A1F_Z: Any
    A0P_Z: Any
    A0F_Z: Any
    aP_Z: Any
    aF_Z: Any
    DeltaA2_Z: Any
    DeltaA1_Z: Any
    DeltaA0_Z: Any
    numerator_Z: Any
    denominator_Z: Any
    Delta_a_Z: Any
    direct_aP_Z: Any
    direct_aF_Z: Any
    direct_Delta_a_Z: Any
    numerator_Z_terms: Mapping[str, Any]
    denominator_Z_terms: Mapping[str, Any]
    direct_root_terms: Mapping[str, Any]
    derivative_residuals: Mapping[str, Any]
    gates: Mapping[str, Any]
    materialization: Mapping[str, Any]

    @property
    def value(self) -> Any:
        return self.base.value

    @property
    def Delta_a(self) -> Any:
        return self.base.Delta_a

    @property
    def delta_a(self) -> Any:
        return self.base.Delta_a

    @property
    def Z1(self) -> Any:
        return self.Delta_a_Z

    @property
    def ordinary_Z1(self) -> Any:
        return self.Delta_a_Z

    @property
    def numerator(self) -> Any:
        return self.base.numerator

    @property
    def denominator(self) -> Any:
        return self.base.denominator

    @property
    def signed_terms(self) -> Mapping[str, Any]:
        return self.base.numerator_terms

    @property
    def gates_all(self) -> Mapping[str, Any]:
        return self.gates

    @property
    def values(self) -> dict[str, Any]:
        return {
            "Delta_a": self.base.Delta_a,
            "numerator": self.base.numerator,
            "denominator": self.base.denominator,
            "DeltaA2": self.base.DeltaA2,
            "DeltaA1": self.base.DeltaA1,
            "DeltaA0": self.base.DeltaA0,
        }

    @property
    def jets(self) -> dict[str, Any]:
        return {
            "A2P_Z": self.A2P_Z,
            "A2F_Z": self.A2F_Z,
            "A1P_Z": self.A1P_Z,
            "A1F_Z": self.A1F_Z,
            "A0P_Z": self.A0P_Z,
            "A0F_Z": self.A0F_Z,
            "aP_Z": self.aP_Z,
            "aF_Z": self.aF_Z,
            "DeltaA2_Z": self.DeltaA2_Z,
            "DeltaA1_Z": self.DeltaA1_Z,
            "DeltaA0_Z": self.DeltaA0_Z,
            "numerator_Z": self.numerator_Z,
            "denominator_Z": self.denominator_Z,
            "Delta_a_Z": self.Delta_a_Z,
            "direct_aP_Z": self.direct_aP_Z,
            "direct_aF_Z": self.direct_aF_Z,
            "direct_Delta_a_Z": self.direct_Delta_a_Z,
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "values": {key: _display(value) for key, value in self.values.items()},
            "Z1": {key: _display(value) for key, value in self.jets.items()},
            "signed_numerator_terms": {
                key: _display(value) for key, value in self.numerator_Z_terms.items()
            },
            "denominator_Z_terms": {
                key: _display(value) for key, value in self.denominator_Z_terms.items()
            },
            "direct_root_terms": {
                key: _display(value) for key, value in self.direct_root_terms.items()
            },
            "derivative_residuals": {
                key: _display(value) for key, value in self.derivative_residuals.items()
            },
            "gates": _display(self.gates),
            "formula": FORMULA,
            "Z1_formula": Z1_FORMULA,
            "materialization": dict(self.materialization),
            **SOURCE_VALUE_FLAGS,
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "delta_a": "Delta_a",
            "Da": "Delta_a",
            "dA2": "DeltaA2",
            "dA1": "DeltaA1",
            "dA0": "DeltaA0",
            "N": "numerator",
            "D": "denominator",
            "D_Z": "denominator_Z",
            "delta_a_Z": "Delta_a_Z",
            "Da_Z": "Delta_a_Z",
        }
        if key in ("gates", "gate"):
            return self.gates
        if key == "values":
            return self.values
        if key in ("jets", "Z1"):
            return self.jets
        if key in ("signed_terms", "terms"):
            return self.numerator_Z_terms
        return getattr(self, aliases.get(key, key))

    def __getattr__(self, name: str) -> Any:
        """Expose point fields on the Z1 packet without copying ownership."""

        if name in {
            "A2P",
            "A2F",
            "A1P",
            "A1F",
            "A0P",
            "A0F",
            "aP",
            "aF",
            "DeltaA2",
            "DeltaA1",
            "DeltaA0",
            "branch_equation_residuals",
            "difference_residuals",
            "numerator_terms",
            "denominator_terms",
        }:
            return getattr(self.base, name)
        raise AttributeError(name)


def ordinary_z1_difference(
    A2P: Any,
    A2F: Any,
    A1P: Any,
    A1F: Any,
    A0P: Any,
    A0F: Any,
    aP: Any,
    aF: Any,
    *,
    A2P_Z: Any = None,
    A2F_Z: Any = None,
    A1P_Z: Any = None,
    A1F_Z: Any = None,
    A0P_Z: Any = None,
    A0F_Z: Any = None,
    aP_Z: Any = None,
    aF_Z: Any = None,
    DeltaA2: Any = None,
    DeltaA1: Any = None,
    DeltaA0: Any = None,
    DeltaA2_Z: Any = None,
    DeltaA1_Z: Any = None,
    DeltaA0_Z: Any = None,
    source_kind: str = "synthetic_point",
    rounded: bool = False,
    q_owner: Any = None,
    p_owner: Any = None,
    f_owner: Any = None,
    derivative_semantics: Any = None,
) -> NoncommonSelectorZ1:
    """Return the complete P-minus-F selector value and ordinary-Z1 routes.

    Every P/F coefficient and root derivative is required.  A missing
    derivative raises ``PFSelectorInputError``; an explicit zero derivative is
    valid.  Independent delta derivatives are accepted and never overwritten
    by branch subtraction. All delta values and ordinary-Z derivatives are
    required explicitly and remain subject to the unpassed consistency gate.
    """

    A2P0, A2Pz = _jet_pair(A2P, A2P_Z, "A2P")
    A2F0, A2Fz = _jet_pair(A2F, A2F_Z, "A2F")
    A1P0, A1Pz = _jet_pair(A1P, A1P_Z, "A1P")
    A1F0, A1Fz = _jet_pair(A1F, A1F_Z, "A1F")
    A0P0, A0Pz = _jet_pair(A0P, A0P_Z, "A0P")
    A0F0, A0Fz = _jet_pair(A0F, A0F_Z, "A0F")
    aP0, aPz = _jet_pair(aP, aP_Z, "aP")
    aF0, aFz = _jet_pair(aF, aF_Z, "aF")

    dA2, dA2z, dA2_ind, dA2z_ind = _optional_delta(
        DeltaA2, DeltaA2_Z, A2P0 - A2F0, A2Pz - A2Fz, "DeltaA2"
    )
    dA1, dA1z, dA1_ind, dA1z_ind = _optional_delta(
        DeltaA1, DeltaA1_Z, A1P0 - A1F0, A1Pz - A1Fz, "DeltaA1"
    )
    dA0, dA0z, dA0_ind, dA0z_ind = _optional_delta(
        DeltaA0, DeltaA0_Z, A0P0 - A0F0, A0Pz - A0Fz, "DeltaA0"
    )
    base = point_difference(
        A2P0,
        A2F0,
        A1P0,
        A1F0,
        A0P0,
        A0F0,
        aP0,
        aF0,
        DeltaA2=dA2,
        DeltaA1=dA1,
        DeltaA0=dA0,
        source_kind=source_kind,
        rounded=rounded,
        q_owner=q_owner,
        p_owner=p_owner,
        f_owner=f_owner,
        derivative_semantics=derivative_semantics,
    )

    numerator_Z_terms = {
        "DeltaA2_Z_times_aF_squared": dA2z * aF0**2,
        "DeltaA2_times_two_aF_times_aF_Z": dA2 * (2 * aF0 * aFz),
        "DeltaA1_Z_times_aF": dA1z * aF0,
        "DeltaA1_times_aF_Z": dA1 * aFz,
        "DeltaA0_Z": dA0z,
    }
    numerator_Z = sum(numerator_Z_terms.values(), 0)
    denominator_Z_terms = {
        "A2P_Z_times_aP_plus_aF": A2Pz * (aP0 + aF0),
        "A2P_times_aP_Z_plus_aF_Z": A2P0 * (aPz + aFz),
        "A1P_Z": A1Pz,
    }
    denominator_Z = sum(denominator_Z_terms.values(), 0)
    Delta_a_Z = -numerator_Z / base.denominator + (
        base.numerator * denominator_Z / base.denominator**2
    )

    direct_root_terms = {
        "P_A2_Z_times_aP_squared": A2Pz * aP0**2,
        "P_A1_Z_times_aP": A1Pz * aP0,
        "P_A0_Z": A0Pz,
        "F_A2_Z_times_aF_squared": A2Fz * aF0**2,
        "F_A1_Z_times_aF": A1Fz * aF0,
        "F_A0_Z": A0Fz,
        "P_implicit_denominator": 2 * A2P0 * aP0 + A1P0,
        "F_implicit_denominator": 2 * A2F0 * aF0 + A1F0,
    }
    direct_P_numerator = (
        direct_root_terms["P_A2_Z_times_aP_squared"]
        + direct_root_terms["P_A1_Z_times_aP"]
        + direct_root_terms["P_A0_Z"]
    )
    direct_F_numerator = (
        direct_root_terms["F_A2_Z_times_aF_squared"]
        + direct_root_terms["F_A1_Z_times_aF"]
        + direct_root_terms["F_A0_Z"]
    )
    direct_aP_Z = _safe_divide(
        -direct_P_numerator,
        direct_root_terms["P_implicit_denominator"],
        "P implicit root derivative denominator",
    )
    direct_aF_Z = _safe_divide(
        -direct_F_numerator,
        direct_root_terms["F_implicit_denominator"],
        "F implicit root derivative denominator",
    )
    direct_Delta_a_Z = aPz - aFz
    derivative_residuals = {
        "aP_Z_minus_implicit_root_derivative": aPz - direct_aP_Z,
        "aF_Z_minus_implicit_root_derivative": aFz - direct_aF_Z,
        "Delta_a_Z_minus_direct_branch_difference": Delta_a_Z - direct_Delta_a_Z,
        "DeltaA2_Z_minus_branch_difference": dA2z - (A2Pz - A2Fz),
        "DeltaA1_Z_minus_branch_difference": dA1z - (A1Pz - A1Fz),
        "DeltaA0_Z_minus_branch_difference": dA0z - (A0Pz - A0Fz),
    }

    gates = dict(base.gates)
    gates.update(
        {
            "coefficient_Z1_difference_consistency": _gate(
                "caller proof that each DeltaAi_Z equals AiP_Z-AiF_Z",
                detail={
                    "DeltaA2_Z_residual": derivative_residuals[
                        "DeltaA2_Z_minus_branch_difference"
                    ],
                    "DeltaA1_Z_residual": derivative_residuals[
                        "DeltaA1_Z_minus_branch_difference"
                    ],
                    "DeltaA0_Z_residual": derivative_residuals[
                        "DeltaA0_Z_minus_branch_difference"
                    ],
                },
            ),
            "ordinary_Z1_transport": _gate(
                "caller proof of ordinary-Z transport and derivative commutation"
            ),
            "direct_root_derivative": _gate(
                "caller proof that aP/aF are the selected differentiable roots"
            ),
        }
    )
    materialization = dict(base.materialization)
    materialization.update(
        {
            "ordinary_Z1_chain_terms_retained": True,
            "all_P_F_value_jets_required": True,
            "all_P_F_root_Z1_required": True,
            "DeltaA_values_independent": {
                "DeltaA2": dA2_ind,
                "DeltaA1": dA1_ind,
                "DeltaA0": dA0_ind,
            },
            "DeltaA_Z_independent": {
                "DeltaA2_Z": dA2z_ind,
                "DeltaA1_Z": dA1z_ind,
                "DeltaA0_Z": dA0z_ind,
            },
            "DeltaA2_Z_independent": dA2z_ind,
            "DeltaA1_Z_independent": dA1z_ind,
            "DeltaA0_Z_independent": dA0z_ind,
            "DeltaA_Z_not_forced_by_subtraction": True,
            "direct_implicit_root_derivatives_exposed": True,
            "conditional_Z1_formula_returned": True,
        }
    )
    return NoncommonSelectorZ1(
        base,
        A2Pz,
        A2Fz,
        A1Pz,
        A1Fz,
        A0Pz,
        A0Fz,
        aPz,
        aFz,
        dA2z,
        dA1z,
        dA0z,
        numerator_Z,
        denominator_Z,
        Delta_a_Z,
        direct_aP_Z,
        direct_aF_Z,
        direct_Delta_a_Z,
        numerator_Z_terms,
        denominator_Z_terms,
        direct_root_terms,
        derivative_residuals,
        gates,
        materialization,
    )


def selector_difference_from_coefficients(
    AP: Mapping[str, Any] | Sequence[Any],
    AF: Mapping[str, Any] | Sequence[Any],
    aP: Any,
    aF: Any,
    **kwargs: Any,
) -> NoncommonSelectorDifference:
    """Convenience wrapper accepting ``AP``/``AF`` triples.

    Mapping keys may be ``A2/A1/A0`` or ``2/1/0``; a three-entry sequence is
    interpreted in ``(A2,A1,A0)`` order.  This wrapper is value-only; callers
    needing ordinary-Z1 should use :func:`ordinary_z1_from_coefficients`.
    """

    A2P, A1P, A0P = _coefficient_triplet(AP, "AP")
    A2F, A1F, A0F = _coefficient_triplet(AF, "AF")
    return point_difference(A2P, A2F, A1P, A1F, A0P, A0F, aP, aF, **kwargs)


def ordinary_z1_from_coefficients(
    AP: Mapping[str, Any] | Sequence[Any],
    AF: Mapping[str, Any] | Sequence[Any],
    aP: Any,
    aF: Any,
    **kwargs: Any,
) -> NoncommonSelectorZ1:
    """Convenience wrapper for coefficient triples carrying value/Z jets."""

    A2P, A1P, A0P = _coefficient_triplet(AP, "AP")
    A2F, A1F, A0F = _coefficient_triplet(AF, "AF")
    return ordinary_z1_difference(A2P, A2F, A1P, A1F, A0P, A0F, aP, aF, **kwargs)


def _coefficient_triplet(value: Mapping[str, Any] | Sequence[Any], label: str) -> tuple[Any, Any, Any]:
    if isinstance(value, Mapping):
        keys = (("A2", "A1", "A0"), ("2", "1", "0"))
        for candidates in keys:
            if all(key in value for key in candidates):
                return tuple(value[key] for key in candidates)  # type: ignore[return-value]
        raise PFSelectorInputError(f"{label} must contain A2, A1, and A0 entries")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 3:
            return value[0], value[1], value[2]
    raise PFSelectorInputError(f"{label} must be a coefficient triple")


def _fixture_exact_sympy() -> dict[str, Any]:
    """Build an independent exact SymPy fixture; this function is uncalled."""

    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - optional dependency
        return {
            "schema": SCHEMA,
            "fixture": "independent_exact_sympy_pf_selector_correlated_root",
            "checks_passed": False,
            "sympy_available": False,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            **SOURCE_VALUE_FLAGS,
        }

    z = sp.symbols("z")
    tiny = sp.Rational(1, 10**80)

    # Two independently chosen positive roots and distinct negative companion
    # roots generate exact coefficient functions without importing a source.
    aP_expr = 2 + z / 7
    aF_expr = 1 + z / 5
    otherP_expr = -3 + z / 11
    otherF_expr = -2 + z / 13
    A2P_expr = 3 + z / 17
    A2F_expr = 5 + z / 19
    A1P_expr = -A2P_expr * (aP_expr + otherP_expr)
    A0P_expr = A2P_expr * aP_expr * otherP_expr
    A1F_expr = -A2F_expr * (aF_expr + otherF_expr)
    A0F_expr = A2F_expr * aF_expr * otherF_expr
    coefficient_exprs = {
        "A2P": A2P_expr,
        "A2F": A2F_expr,
        "A1P": A1P_expr,
        "A1F": A1F_expr,
        "A0P": A0P_expr,
        "A0F": A0F_expr,
    }
    root_exprs = {"aP": aP_expr, "aF": aF_expr}
    values = {
        name: sp.expand(expr).subs(z, 0)
        for name, expr in {**coefficient_exprs, **root_exprs}.items()
    }
    jets = {
        name + "_Z": sp.diff(expr, z).subs(z, 0)
        for name, expr in {**coefficient_exprs, **root_exprs}.items()
    }
    delta_values = {
        "DeltaA2": values["A2P"] - values["A2F"],
        "DeltaA1": values["A1P"] - values["A1F"],
        "DeltaA0": values["A0P"] - values["A0F"],
    }
    delta_jets = {
        "DeltaA2_Z": jets["A2P_Z"] - jets["A2F_Z"],
        "DeltaA1_Z": jets["A1P_Z"] - jets["A1F_Z"],
        "DeltaA0_Z": jets["A0P_Z"] - jets["A0F_Z"],
    }
    point = point_difference(
        values["A2P"],
        values["A2F"],
        values["A1P"],
        values["A1F"],
        values["A0P"],
        values["A0F"],
        values["aP"],
        values["aF"],
        **delta_values,
        source_kind="formal_symbolic",
    )
    z1 = ordinary_z1_difference(
        **values,
        **jets,
        **delta_values,
        **delta_jets,
        source_kind="formal_symbolic",
    )
    direct_delta_a = sp.diff(aP_expr - aF_expr, z).subs(z, 0)
    direct_root_P = sp.diff(aP_expr, z).subs(z, 0)
    direct_root_F = sp.diff(aF_expr, z).subs(z, 0)
    uP_expr, uF_expr = 3+z/11, 2+z/13
    vP_expr, vF_expr = 4+z/17, 1+z/19
    delta_C_factored = (uP_expr-uF_expr) + aP_expr*(vP_expr-vF_expr) + (aP_expr-aF_expr)*vF_expr
    delta_C_direct = (uP_expr+aP_expr*vP_expr)-(uF_expr+aF_expr*vF_expr)

    near_A2F = values["A2F"]
    near_A1F = values["A1F"]
    near_A0F = values["A0F"]
    near_aP = values["aF"] + tiny
    near_otherP = -2 - tiny
    near_A2P = near_A2F + tiny
    near_A1P = -near_A2P * (near_aP + near_otherP)
    near_A0P = near_A2P * near_aP * near_otherP
    near = point_difference(
        near_A2P,
        near_A2F,
        near_A1P,
        near_A1F,
        near_A0P,
        near_A0F,
        near_aP,
        values["aF"],
        DeltaA2=tiny,
        DeltaA1=near_A1P-near_A1F,
        DeltaA0=near_A0P-near_A0F,
        source_kind="formal_symbolic",
    )
    inconsistent = point_difference(
        values["A2P"],
        values["A2F"],
        values["A1P"],
        values["A1F"],
        values["A0P"],
        values["A0F"],
        values["aP"],
        values["aF"],
        DeltaA2=delta_values["DeltaA2"] + tiny,
        DeltaA1=delta_values["DeltaA1"],
        DeltaA0=delta_values["DeltaA0"],
        source_kind="formal_symbolic",
    )

    denominator_zero_raises = False
    try:
        point_difference(1, 1, -2, 0, 1, 1, 1, 1,
                         DeltaA2=0, DeltaA1=-2, DeltaA0=0)
    except PFSelectorInputError:
        denominator_zero_raises = True

    missing_derivative_raises = False
    missing_jets = dict(jets)
    missing_jets.pop("A0F_Z")
    try:
        ordinary_z1_difference(
            **values,
            **missing_jets,
            **delta_values,
            **delta_jets,
            source_kind="formal_symbolic",
        )
    except PFSelectorInputError:
        missing_derivative_raises = True

    missing_delta_raises = False
    try:
        point_difference(**values)
    except PFSelectorInputError:
        missing_delta_raises = True
    missing_delta_Z_raises = False
    try:
        incomplete_delta_Z = dict(delta_jets)
        incomplete_delta_Z.pop("DeltaA0_Z")
        ordinary_z1_difference(**values, **jets, **delta_values, **incomplete_delta_Z)
    except PFSelectorInputError:
        missing_delta_Z_raises = True
    opaque_mapping_raises = False
    try:
        opaque = dict(values)
        opaque["A2P"] = {"value": values["A2P"], "node_id": 123}
        point_difference(**opaque, **delta_values)
    except PFSelectorInputError:
        opaque_mapping_raises = True

    checks = {
        "all_three_coefficient_differences_nonzero": all(
            value != 0 for value in delta_values.values()
        ),
        "all_input_Z1_nonzero": all(value != 0 for value in jets.values()),
        "P_branch_equation_exact": sp.simplify(point.branch_equation_residuals["P"]) == 0,
        "F_branch_equation_exact": sp.simplify(point.branch_equation_residuals["F"]) == 0,
        "signed_value_exact": sp.simplify(
            point.Delta_a
            + (
                delta_values["DeltaA2"] * values["aF"] ** 2
                + delta_values["DeltaA1"] * values["aF"]
                + delta_values["DeltaA0"]
            )
            / point.denominator
        )
        == 0,
        "all_Z1_chain_terms_nonzero": all(
            value != 0 for value in z1.numerator_Z_terms.values()
        )
        and all(value != 0 for value in z1.denominator_Z_terms.values()),
        "P_direct_root_derivative_exact": sp.simplify(
            z1.direct_aP_Z - direct_root_P
        )
        == 0,
        "F_direct_root_derivative_exact": sp.simplify(
            z1.direct_aF_Z - direct_root_F
        )
        == 0,
        "Delta_direct_root_derivative_exact": sp.simplify(
            z1.direct_Delta_a_Z - direct_delta_a
        )
        == 0,
        "conditional_Z1_exact_for_consistent_deltas": sp.simplify(
            z1.Delta_a_Z - direct_delta_a
        )
        == 0,
        "coefficient_transport_Z0_Z1_exact": all(
            sp.simplify(sp.diff(delta_C_factored-delta_C_direct, z, j).subs(z, 0)) == 0
            for j in (0, 1)
        ),
        "tiny_nonzero_exact_difference": (
            near.DeltaA2 != 0
            and near.DeltaA1 != 0
            and near.DeltaA0 != 0
            and near.DeltaA2 == tiny
            and near.branch_equation_residuals["P"] == 0
            and near.branch_equation_residuals["F"] == 0
            and near.Delta_a == tiny
        ),
        "inconsistent_delta_retained_and_gate_unpassed": (
            inconsistent.difference_residuals["DeltaA2_minus_branch_difference"]
            == tiny
            and inconsistent.gates["coefficient_difference_consistency"]["passed"]
            is False
        ),
        "denominator_zero_rejected": denominator_zero_raises,
        "missing_Z_derivative_rejected": missing_derivative_raises,
        "missing_signed_delta_rejected": missing_delta_raises,
        "missing_signed_delta_Z_rejected": missing_delta_Z_raises,
        "opaque_metadata_mapping_rejected": opaque_mapping_raises,
        "gates_remain_unpassed": all(
            gate["passed"] is False for gate in z1.gates.values()
        ),
        "no_midpoint_or_rounding_claim": (
            point.materialization["parameter_interval_midpoint_used"] is False
            and point.materialization["rounding_certified"] is False
        ),
    }
    return {
        "schema": SCHEMA,
        "fixture": "independent_exact_sympy_pf_selector_correlated_root",
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "sympy_available": True,
        "exact_arithmetic": True,
        "all_three_differences": {
            key: str(value) for key, value in delta_values.items()
        },
        "all_Z1_inputs": {key: str(value) for key, value in jets.items()},
        "tiny_exact_difference": str(tiny),
        **SOURCE_VALUE_FLAGS,
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run the optional synthetic fixture only when explicitly requested."""

    return _fixture_exact_sympy()


# Source-notation aliases keep the component discoverable without importing an
# earlier common-source selector implementation.
selector_root_difference = point_difference
correlated_root_difference = point_difference
correlated_selector_difference = point_difference
signed_selector_difference = point_difference
signed_root_difference = point_difference
build_selector_difference = point_difference
build_root_difference = point_difference
selector_root_z1_difference = ordinary_z1_difference
correlated_root_z1_difference = ordinary_z1_difference
correlated_selector_z1_difference = ordinary_z1_difference
signed_root_z1_difference = ordinary_z1_difference
ordinary_Z1_difference = ordinary_z1_difference
build_selector_z1 = ordinary_z1_difference
build_root_z1 = ordinary_z1_difference
fixture_exact_sympy = _fixture_exact_sympy
PFSelectorDifference = NoncommonSelectorDifference
PFSelectorZ1 = NoncommonSelectorZ1


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "FORMULA",
    "Z1_FORMULA",
    "SOURCE_VALUE_FLAGS",
    "PFSelectorInputError",
    "OpaqueSourceReference",
    "SourceOpaqueId",
    "Z1Jet",
    "NoncommonSelectorDifference",
    "NoncommonSelectorZ1",
    "point_difference",
    "selector_root_difference",
    "correlated_root_difference",
    "correlated_selector_difference",
    "signed_selector_difference",
    "signed_root_difference",
    "build_selector_difference",
    "build_root_difference",
    "ordinary_z1_difference",
    "selector_root_z1_difference",
    "correlated_root_z1_difference",
    "correlated_selector_z1_difference",
    "signed_root_z1_difference",
    "ordinary_Z1_difference",
    "build_selector_z1",
    "build_root_z1",
    "PFSelectorDifference",
    "PFSelectorZ1",
    "selector_difference_from_coefficients",
    "ordinary_z1_from_coefficients",
    "fixture_exact_sympy",
    "verify_synthetic_fixture",
]
