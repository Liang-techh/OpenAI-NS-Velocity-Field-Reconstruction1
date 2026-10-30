"""Source-shaped P/F selector algebra with caller-owned signed differences.

The original selector first forms ``R_i=k_i I_i`` and the two ``r``/``t``
right-hand sides, solves the two by two system, and then builds

``A2=K+sum(nu*v**2)``, ``A1=2*sum(nu*u*v)`` and
``A0=-future+mu*E-base+sum(nu*u**2)``.

This module ports only that algebra.  Every primitive is supplied by the
caller as separate P and F values, an explicit signed ``Delta`` value, and
the corresponding ordinary-Z1 rows.  In particular, no missing difference
is replaced by ``P-F``.  Product differences use

``Delta(x*y)=xP*Delta(y)+Delta(x)*yF``

and quotients use

``Delta(N/det)=(Delta(N)*detF-NF*Delta(det))/(detP*detF)``.

The module has no source graph, old selector import, value provider,
precision policy, or certification claim.  The optional exact SymPy fixture is
defined but is never called at import time.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any


SCHEMA = "pf_selector_primitive_delta_v2"
COMPONENT = SCHEMA
FORMULA = "Delta(N/det)=(DeltaN*detF-NF*DeltaDet)/(detP*detF)"
PRODUCT_FORMULA = "Delta(x*y)=xP*DeltaY+DeltaX*yF"
Z1_FORMULA = "ordinary-Z1 chain rule applied to every signed product and quotient route"

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


class PrimitiveDeltaInputError(ValueError):
    """Base error for incomplete or non-arithmetic caller input."""


class MissingSignedDifferenceError(PrimitiveDeltaInputError):
    """Raised when a primitive omits its explicit signed difference or Z1."""


class OpaqueSourceValueError(PrimitiveDeltaInputError):
    """Raised when saved source metadata is supplied as an arithmetic value."""


class OpaqueSourceReference:
    """Explicit marker for a source node/role, never an arithmetic scalar.

    Plain integers remain caller-owned arithmetic values at this boundary: a
    bare integer cannot be distinguished from an integer source value.  Saved
    handles must therefore be wrapped here or supplied with an explicit
    ``node_id``/``source_role`` marker, which is rejected before arithmetic.
    """

    __slots__ = ("identifier", "owner")

    def __init__(self, identifier: Any, owner: Any = None) -> None:
        self.identifier = identifier
        self.owner = owner

    def __repr__(self) -> str:  # pragma: no cover - diagnostic only
        return f"OpaqueSourceReference({self.identifier!r}, owner={self.owner!r})"


SourceOpaqueId = OpaqueSourceReference


_MISSING = object()
_OPAQUE_KEYS = {
    "source_node_id",
    "node_id",
    "source_role",
    "namespace",
    "kind",
    "reference",
    "saved_reference",
    "mapped_reference",
}


def _display(value: Any) -> Any:
    """Serialize diagnostics without coercing caller-owned arithmetic."""

    if isinstance(value, Mapping):
        return {key: _display(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_display(item) for item in value]
    try:
        import json

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
    if isinstance(value, OpaqueSourceReference):
        return True
    if isinstance(value, Mapping):
        if set(value).intersection(_OPAQUE_KEYS):
            return True
        for marker in (
            "opaque",
            "is_opaque",
            "source_owned",
            "reference_is_opaque_metadata",
            "node_id_is_not_a_numeric_value",
            "roles_are_opaque_metadata",
        ):
            if value.get(marker) is True:
                return True
    for attribute in (
        "source_node_id",
        "node_id",
        "source_role",
        "is_opaque",
        "opaque",
    ):
        try:
            if getattr(value, attribute) is not None:
                return True
        except Exception:
            pass
    return False


def _reject_opaque(value: Any, label: str) -> None:
    if _is_opaque(value):
        raise OpaqueSourceValueError(
            f"{label} is source-owned metadata, not an arithmetic value"
        )
    if isinstance(value, (str, bytes, bytearray)):
        raise OpaqueSourceValueError(f"{label} is text metadata, not an arithmetic value")
    if value is None:
        raise PrimitiveDeltaInputError(f"{label} is missing")


def _pick(mapping: Mapping[str, Any], keys: Sequence[str], label: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    raise PrimitiveDeltaInputError(f"{label} is missing")


def _maybe_pick(mapping: Mapping[str, Any], keys: Sequence[str]) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return _MISSING


def _branch_jet(value: Any, derivative: Any, label: str) -> tuple[Any, Any]:
    """Read one branch as ``(value,Z1)`` without inventing a derivative."""

    if isinstance(value, Mapping) and _is_opaque(value):
        _reject_opaque(value, label)
    if derivative is not _MISSING:
        _reject_opaque(value, label)
        _reject_opaque(derivative, label + "_Z")
        return value, derivative
    if isinstance(value, Z1Jet):
        _reject_opaque(value.value, label)
        _reject_opaque(value.Z, label + "_Z")
        return value.value, value.Z
    if isinstance(value, Mapping):
        value_key = _maybe_pick(value, ("value", "base", "Z0", "z0", "0"))
        z_key = _maybe_pick(value, ("Z", "z", "derivative", "Z1", "z1", "1"))
        if value_key is not _MISSING and z_key is not _MISSING:
            _reject_opaque(value_key, label)
            _reject_opaque(z_key, label + "_Z")
            return value_key, z_key
        raise PrimitiveDeltaInputError(f"{label} requires value and ordinary-Z1 entries")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            _reject_opaque(value[0], label)
            _reject_opaque(value[1], label + "_Z")
            return value[0], value[1]
        raise PrimitiveDeltaInputError(f"{label} sequence jet must have two entries")
    raise PrimitiveDeltaInputError(f"{label} requires an explicit ordinary-Z1 derivative")


@dataclass(frozen=True)
class Z1Jet:
    """One caller-owned scalar and its ordinary-Z1 row."""

    value: Any
    Z: Any

    @property
    def derivative(self) -> Any:
        return self.Z

    def __post_init__(self) -> None:
        _reject_opaque(self.value, "Z1Jet.value")
        _reject_opaque(self.Z, "Z1Jet.Z")


@dataclass(frozen=True, init=False)
class PrimitiveDelta:
    """P/F value and explicit signed ordinary-Z1 difference.

    ``P``/``F`` and their derivatives are retained even when the supplied
    ``Delta`` is intentionally inconsistent with ``P-F``.  The inconsistency
    is diagnostic evidence for an unpassed caller gate, never a reason to
    overwrite the caller's signed input.
    """

    P: Any
    F: Any
    Delta: Any
    P_Z: Any
    F_Z: Any
    Delta_Z: Any
    name: str
    owner: Any
    delta_terms: Mapping[str, Any]
    delta_Z_terms: Mapping[str, Any]

    def __init__(
        self,
        P: Any = _MISSING,
        F: Any = _MISSING,
        Delta: Any = _MISSING,
        P_Z: Any = _MISSING,
        F_Z: Any = _MISSING,
        Delta_Z: Any = _MISSING,
        *,
        p: Any = _MISSING,
        f: Any = _MISSING,
        delta: Any = _MISSING,
        p_Z: Any = _MISSING,
        f_Z: Any = _MISSING,
        delta_Z: Any = _MISSING,
        name: str = "primitive",
        owner: Any = None,
        delta_terms: Mapping[str, Any] | None = None,
        delta_Z_terms: Mapping[str, Any] | None = None,
        **aliases: Any,
    ) -> None:
        P_Z = P_Z if P_Z is not _MISSING else aliases.pop("PZ", _MISSING)
        F_Z = F_Z if F_Z is not _MISSING else aliases.pop("FZ", _MISSING)
        Delta_Z = Delta_Z if Delta_Z is not _MISSING else aliases.pop("DeltaZ", _MISSING)
        p_Z = p_Z if p_Z is not _MISSING else aliases.pop("pZ", _MISSING)
        f_Z = f_Z if f_Z is not _MISSING else aliases.pop("fZ", _MISSING)
        delta_Z = delta_Z if delta_Z is not _MISSING else aliases.pop("deltaZ", _MISSING)
        if aliases:
            unexpected = next(iter(aliases))
            raise TypeError(f"unexpected PrimitiveDelta keyword: {unexpected}")
        pairs = (
            ("P", P, p),
            ("F", F, f),
            ("Delta", Delta, delta),
            ("P_Z", P_Z, p_Z),
            ("F_Z", F_Z, f_Z),
            ("Delta_Z", Delta_Z, delta_Z),
        )
        values: dict[str, Any] = {}
        for field_name, upper, lower in pairs:
            value = upper if upper is not _MISSING else lower
            if value is _MISSING:
                raise MissingSignedDifferenceError(
                    f"{field_name} is required for explicit P/F Delta algebra"
                )
            _reject_opaque(value, field_name)
            values[field_name] = value
        object.__setattr__(self, "P", values["P"])
        object.__setattr__(self, "F", values["F"])
        object.__setattr__(self, "Delta", values["Delta"])
        object.__setattr__(self, "P_Z", values["P_Z"])
        object.__setattr__(self, "F_Z", values["F_Z"])
        object.__setattr__(self, "Delta_Z", values["Delta_Z"])
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "owner", owner)
        object.__setattr__(self, "delta_terms", dict(delta_terms or {"explicit_Delta": values["Delta"]}))
        object.__setattr__(self, "delta_Z_terms", dict(delta_Z_terms or {"explicit_Delta_Z": values["Delta_Z"]}))

    @property
    def p(self) -> Any:
        return self.P

    @property
    def f(self) -> Any:
        return self.F

    @property
    def delta(self) -> Any:
        return self.Delta

    @property
    def p_z(self) -> Any:
        return self.P_Z

    @property
    def f_z(self) -> Any:
        return self.F_Z

    @property
    def delta_z(self) -> Any:
        return self.Delta_Z

    @property
    def Z1(self) -> Any:
        return self.Delta_Z

    @property
    def value(self) -> Any:
        return self.Delta

    @property
    def difference_residual(self) -> Any:
        return self.Delta - (self.P - self.F)

    @property
    def difference_Z1_residual(self) -> Any:
        return self.Delta_Z - (self.P_Z - self.F_Z)

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "P": _display(self.P),
            "F": _display(self.F),
            "Delta": _display(self.Delta),
            "P_Z": _display(self.P_Z),
            "F_Z": _display(self.F_Z),
            "Delta_Z": _display(self.Delta_Z),
            "delta_terms": _display(self.delta_terms),
            "delta_Z_terms": _display(self.delta_Z_terms),
            "difference_residual": _display(self.difference_residual),
            "difference_Z1_residual": _display(self.difference_Z1_residual),
            "owner": _display(self.owner),
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "p": "P",
            "f": "F",
            "delta": "Delta",
            "Pz": "P_Z",
            "Fz": "F_Z",
            "Deltaz": "Delta_Z",
            "d": "Delta",
            "dZ": "Delta_Z",
            "Z1": "Delta_Z",
        }
        return getattr(self, aliases.get(key, key))


PFPrimitive = PrimitiveDelta
SignedPrimitive = PrimitiveDelta
PrimitiveJet = PrimitiveDelta


def primitive_delta(
    P: Any = _MISSING,
    F: Any = _MISSING,
    Delta: Any = _MISSING,
    P_Z: Any = _MISSING,
    F_Z: Any = _MISSING,
    Delta_Z: Any = _MISSING,
    **kwargs: Any,
) -> PrimitiveDelta:
    """Construct one explicit caller-owned primitive."""

    return PrimitiveDelta(P, F, Delta, P_Z, F_Z, Delta_Z, **kwargs)


def _primitive(value: Any, label: str) -> PrimitiveDelta:
    """Normalize the accepted explicit primitive spellings."""

    if isinstance(value, PrimitiveDelta):
        return value
    if isinstance(value, Mapping):
        if set(value).intersection(_OPAQUE_KEYS) and not set(value).intersection(
            {"P", "p", "F", "f", "Delta", "delta", "D", "difference"}
        ):
            _reject_opaque(value, label)
        # A top-level primitive may carry branch jets or scalar branches plus
        # explicit derivative keys.  Delta is deliberately mandatory.
        p_raw = _pick(value, ("P", "p", "P_value", "value_P", "p_value"), label + ".P")
        f_raw = _pick(value, ("F", "f", "F_value", "value_F", "f_value"), label + ".F")
        d_raw = _maybe_pick(value, ("Delta", "delta", "D", "signed_difference", "difference"))
        if d_raw is _MISSING:
            raise MissingSignedDifferenceError(label + ".Delta is missing")
        z1 = value.get("Z1", value.get("z1", _MISSING))
        pz = _maybe_pick(value, ("P_Z", "Pz", "PZ", "p_Z", "pZ", "Z_P", "P_Z1", "Z1_P"))
        fz = _maybe_pick(value, ("F_Z", "Fz", "FZ", "f_Z", "fZ", "Z_F", "F_Z1", "Z1_F"))
        dz = _maybe_pick(value, ("Delta_Z", "DeltaZ", "Deltaz", "delta_Z", "deltaZ", "Z_Delta", "Delta_Z1", "Z1_Delta"))
        if z1 is not _MISSING and isinstance(z1, Mapping):
            if pz is _MISSING:
                pz = _maybe_pick(z1, ("P", "p", "P_Z", "p_Z"))
            if fz is _MISSING:
                fz = _maybe_pick(z1, ("F", "f", "F_Z", "f_Z"))
            if dz is _MISSING:
                dz = _maybe_pick(z1, ("Delta", "delta", "Delta_Z", "delta_Z"))
        if pz is _MISSING:
            pz = _MISSING
        if fz is _MISSING:
            fz = _MISSING
        if dz is _MISSING:
            dz = _MISSING
        p, pz = _branch_jet(p_raw, pz, label + ".P")
        f, fz = _branch_jet(f_raw, fz, label + ".F")
        try:
            d, dz = _branch_jet(d_raw, dz, label + ".Delta")
        except OpaqueSourceValueError:
            raise
        except PrimitiveDeltaInputError as exc:
            raise MissingSignedDifferenceError(label + ".Delta_Z is missing or malformed") from exc
        owner = value.get("owner", value.get("source_owner"))
        return PrimitiveDelta(
            p,
            f,
            d,
            pz,
            fz,
            dz,
            name=str(value.get("name", label)),
            owner=owner,
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 6:
            return PrimitiveDelta(*value, name=label)
        if len(value) == 3:
            p, pz = _branch_jet(value[0], _MISSING, label + ".P")
            f, fz = _branch_jet(value[1], _MISSING, label + ".F")
            try:
                d, dz = _branch_jet(value[2], _MISSING, label + ".Delta")
            except OpaqueSourceValueError:
                raise
            except PrimitiveDeltaInputError as exc:
                raise MissingSignedDifferenceError(label + ".Delta_Z is missing or malformed") from exc
            return PrimitiveDelta(p, f, d, pz, fz, dz, name=label)
    raise PrimitiveDeltaInputError(
        f"{label} must provide P, F, Delta and all three ordinary-Z1 rows"
    )


def _vector(value: Any, label: str) -> tuple[PrimitiveDelta, PrimitiveDelta]:
    if isinstance(value, Mapping):
        row_aliases = {
            "first_row": (("A0", "A_0", "first_row_0"), ("A1", "A_1", "first_row_1")),
            "divided_row": (("D0", "D_0", "divided_row_0"), ("D1", "D_1", "divided_row_1")),
            "incoming_factors": (("k0", "k_0", "factor0"), ("k1", "k_1", "factor1")),
            "incoming": (("I0", "I_0", "incoming0"), ("I1", "I_1", "incoming1")),
            "scaled_main_rows": (("M0", "M_0", "main0"), ("M1", "M_1", "main1")),
            "end_energy": (("nu0", "ν0", "nu_0", "end0"), ("nu1", "ν1", "nu_1", "end1")),
        }
        aliases = row_aliases.get(label, ((), ()))
        first = _maybe_pick(value, (0, "0", "first", "row0", label + "0", *aliases[0]))
        second = _maybe_pick(value, (1, "1", "second", "row1", label + "1", *aliases[1]))
        if first is _MISSING or second is _MISSING:
            raise PrimitiveDeltaInputError(f"{label} must contain rows 0 and 1")
        return _primitive(first, f"{label}[0]"), _primitive(second, f"{label}[1]")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) < 2:
            raise PrimitiveDeltaInputError(f"{label} must contain rows 0 and 1")
        return _primitive(value[0], f"{label}[0]"), _primitive(value[1], f"{label}[1]")
    raise PrimitiveDeltaInputError(f"{label} must be a two-row sequence or mapping")


def _zero_like(value: Any) -> Any:
    return value * 0


def _sum(values: Sequence[Any]) -> Any:
    if not values:
        return 0
    total = values[0]
    for value in values[1:]:
        total = total + value
    return total


def _is_definitely_zero(value: Any) -> bool:
    if value is None:
        return False
    try:
        equal = value == 0
    except Exception:
        return False
    return equal is True or type(equal).__name__ == "BooleanTrue"


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


def _child_terms(prefix: str, terms: Mapping[str, Any]) -> dict[str, Any]:
    return {prefix: dict(terms)}


def signed_add(left: PrimitiveDelta, right: PrimitiveDelta, *, name: str = "sum") -> PrimitiveDelta:
    left = _primitive(left, name + ".left")
    right = _primitive(right, name + ".right")
    return PrimitiveDelta(
        left.P + right.P,
        left.F + right.F,
        left.Delta + right.Delta,
        left.P_Z + right.P_Z,
        left.F_Z + right.F_Z,
        left.Delta_Z + right.Delta_Z,
        name=name,
        delta_terms={
            "left_Delta": left.Delta,
            "right_Delta": right.Delta,
            **_child_terms("left_delta_terms", left.delta_terms),
            **_child_terms("right_delta_terms", right.delta_terms),
        },
        delta_Z_terms={
            "left_Delta_Z": left.Delta_Z,
            "right_Delta_Z": right.Delta_Z,
            **_child_terms("left_delta_Z_terms", left.delta_Z_terms),
            **_child_terms("right_delta_Z_terms", right.delta_Z_terms),
        },
    )


def signed_negate(value: PrimitiveDelta, *, name: str = "negative") -> PrimitiveDelta:
    value = _primitive(value, name + ".argument")
    return PrimitiveDelta(
        -value.P,
        -value.F,
        -value.Delta,
        -value.P_Z,
        -value.F_Z,
        -value.Delta_Z,
        name=name,
        delta_terms={"negative_Delta": -value.Delta, **_child_terms("child_delta_terms", value.delta_terms)},
        delta_Z_terms={"negative_Delta_Z": -value.Delta_Z, **_child_terms("child_delta_Z_terms", value.delta_Z_terms)},
    )


def signed_subtract(left: PrimitiveDelta, right: PrimitiveDelta, *, name: str = "difference") -> PrimitiveDelta:
    return signed_add(left, signed_negate(right, name=name + ".negative_right"), name=name)


def signed_scale(factor: Any, value: PrimitiveDelta, *, name: str = "scale") -> PrimitiveDelta:
    _reject_opaque(factor, name + ".factor")
    value = _primitive(value, name + ".argument")
    return PrimitiveDelta(
        factor * value.P,
        factor * value.F,
        factor * value.Delta,
        factor * value.P_Z,
        factor * value.F_Z,
        factor * value.Delta_Z,
        name=name,
        delta_terms={
            "factor_times_Delta": factor * value.Delta,
            **_child_terms("child_delta_terms", value.delta_terms),
        },
        delta_Z_terms={
            "factor_times_Delta_Z": factor * value.Delta_Z,
            **_child_terms("child_delta_Z_terms", value.delta_Z_terms),
        },
    )


def signed_product(left: PrimitiveDelta, right: PrimitiveDelta, *, name: str = "product") -> PrimitiveDelta:
    """Use the asymmetric source-shaped ``Delta(xy)`` identity."""

    left = _primitive(left, name + ".left")
    right = _primitive(right, name + ".right")
    delta_terms = {
        "left_P_times_right_Delta": left.P * right.Delta,
        "left_Delta_times_right_F": left.Delta * right.F,
        **_child_terms("left_delta_terms", left.delta_terms),
        **_child_terms("right_delta_terms", right.delta_terms),
    }
    delta_Z_terms = {
        "left_P_Z_times_right_Delta": left.P_Z * right.Delta,
        "left_P_times_right_Delta_Z": left.P * right.Delta_Z,
        "left_Delta_Z_times_right_F": left.Delta_Z * right.F,
        "left_Delta_times_right_F_Z": left.Delta * right.F_Z,
        **_child_terms("left_delta_Z_terms", left.delta_Z_terms),
        **_child_terms("right_delta_Z_terms", right.delta_Z_terms),
    }
    return PrimitiveDelta(
        left.P * right.P,
        left.F * right.F,
        _sum((delta_terms["left_P_times_right_Delta"], delta_terms["left_Delta_times_right_F"])),
        left.P_Z * right.P + left.P * right.P_Z,
        left.F_Z * right.F + left.F * right.F_Z,
        _sum(
            (
                delta_Z_terms["left_P_Z_times_right_Delta"],
                delta_Z_terms["left_P_times_right_Delta_Z"],
                delta_Z_terms["left_Delta_Z_times_right_F"],
                delta_Z_terms["left_Delta_times_right_F_Z"],
            )
        ),
        name=name,
        delta_terms=delta_terms,
        delta_Z_terms=delta_Z_terms,
    )


def signed_quotient(
    numerator: PrimitiveDelta,
    denominator: PrimitiveDelta,
    *,
    name: str = "quotient",
) -> PrimitiveDelta:
    """Use the signed quotient identity, retaining all numerator routes."""

    numerator = _primitive(numerator, name + ".numerator")
    denominator = _primitive(denominator, name + ".denominator")
    if _is_definitely_zero(denominator.P) or _is_definitely_zero(denominator.F):
        raise PrimitiveDeltaInputError(f"{name} requires nonzero P and F denominator values")

    numerator_terms = {
        "DeltaN_times_detF": numerator.Delta * denominator.F,
        "minus_NF_times_DeltaDet": -(numerator.F * denominator.Delta),
    }
    signed_numerator = _sum(tuple(numerator_terms.values()))
    denominator_terms = {
        "detP_times_detF": denominator.P * denominator.F,
    }
    product_denominator = denominator_terms["detP_times_detF"]
    numerator_Z_terms = {
        "DeltaN_Z_times_detF": numerator.Delta_Z * denominator.F,
        "DeltaN_times_detF_Z": numerator.Delta * denominator.F_Z,
        "minus_NF_Z_times_DeltaDet": -(numerator.F_Z * denominator.Delta),
        "minus_NF_times_DeltaDet_Z": -(numerator.F * denominator.Delta_Z),
    }
    signed_numerator_Z = _sum(tuple(numerator_Z_terms.values()))
    denominator_Z_terms = {
        "detP_Z_times_detF": denominator.P_Z * denominator.F,
        "detP_times_detF_Z": denominator.P * denominator.F_Z,
    }
    product_denominator_Z = _sum(tuple(denominator_Z_terms.values()))
    delta_Z = (
        signed_numerator_Z * product_denominator
        - signed_numerator * product_denominator_Z
    ) / product_denominator**2
    return PrimitiveDelta(
        numerator.P / denominator.P,
        numerator.F / denominator.F,
        signed_numerator / product_denominator,
        (numerator.P_Z * denominator.P - numerator.P * denominator.P_Z) / denominator.P**2,
        (numerator.F_Z * denominator.F - numerator.F * denominator.F_Z) / denominator.F**2,
        delta_Z,
        name=name,
        delta_terms={
            **numerator_terms,
            "signed_numerator": signed_numerator,
            **denominator_terms,
            "signed_quotient": signed_numerator / product_denominator,
            **_child_terms("numerator_delta_terms", numerator.delta_terms),
            **_child_terms("denominator_delta_terms", denominator.delta_terms),
        },
        delta_Z_terms={
            **numerator_Z_terms,
            "signed_numerator_Z": signed_numerator_Z,
            **denominator_Z_terms,
            "signed_denominator_Z": product_denominator_Z,
            "quotient_rule_numerator": signed_numerator_Z * product_denominator,
            "quotient_rule_denominator": -signed_numerator * product_denominator_Z,
            "signed_quotient_Z": delta_Z,
            **_child_terms("numerator_delta_Z_terms", numerator.delta_Z_terms),
            **_child_terms("denominator_delta_Z_terms", denominator.delta_Z_terms),
        },
    )


delta_product = signed_product
delta_quotient = signed_quotient
product_delta = signed_product
quotient_delta = signed_quotient


def _add_many(values: Sequence[PrimitiveDelta], *, name: str) -> PrimitiveDelta:
    if not values:
        raise PrimitiveDeltaInputError(f"{name} cannot be empty")
    result = values[0]
    for index, value in enumerate(values[1:], start=1):
        result = signed_add(result, value, name=f"{name}[0:{index + 1}]")
    return result


@dataclass(frozen=True)
class SelectorPrimitiveDeltaResult:
    """Complete source-shaped selector primitive and coefficient packet."""

    inputs: Mapping[str, Any]
    R: tuple[PrimitiveDelta, PrimitiveDelta]
    r: tuple[PrimitiveDelta, PrimitiveDelta]
    t: tuple[PrimitiveDelta, PrimitiveDelta]
    u: tuple[PrimitiveDelta, PrimitiveDelta]
    v: tuple[PrimitiveDelta, PrimitiveDelta]
    coefficients: Mapping[str, PrimitiveDelta]
    coefficient_delta_terms: Mapping[str, Any]
    input_difference_residuals: Mapping[str, Any]
    derived_difference_residuals: Mapping[str, Any]
    owners: Mapping[str, Any]
    gates: Mapping[str, Any]
    materialization: Mapping[str, Any]

    @property
    def A2(self) -> PrimitiveDelta:
        return self.coefficients["A2"]

    @property
    def A1(self) -> PrimitiveDelta:
        return self.coefficients["A1"]

    @property
    def A0(self) -> PrimitiveDelta:
        return self.coefficients["A0"]

    @property
    def DeltaA2(self) -> Any:
        return self.A2.Delta

    @property
    def DeltaA1(self) -> Any:
        return self.A1.Delta

    @property
    def DeltaA0(self) -> Any:
        return self.A0.Delta

    @property
    def DeltaA2_Z(self) -> Any:
        return self.A2.Delta_Z

    @property
    def DeltaA1_Z(self) -> Any:
        return self.A1.Delta_Z

    @property
    def DeltaA0_Z(self) -> Any:
        return self.A0.Delta_Z

    @property
    def dA2(self) -> Any:
        return self.DeltaA2

    @property
    def dA1(self) -> Any:
        return self.DeltaA1

    @property
    def dA0(self) -> Any:
        return self.DeltaA0

    @property
    def Delta_A2(self) -> Any:
        return self.DeltaA2

    @property
    def Delta_A1(self) -> Any:
        return self.DeltaA1

    @property
    def Delta_A0(self) -> Any:
        return self.DeltaA0

    @property
    def A2P(self) -> Any:
        return self.A2.P

    @property
    def A2F(self) -> Any:
        return self.A2.F

    @property
    def A1P(self) -> Any:
        return self.A1.P

    @property
    def A1F(self) -> Any:
        return self.A1.F

    @property
    def A0P(self) -> Any:
        return self.A0.P

    @property
    def A0F(self) -> Any:
        return self.A0.F

    @property
    def delta_terms(self) -> Mapping[str, Any]:
        return self.coefficient_delta_terms

    @property
    def z1(self) -> Mapping[str, Any]:
        return {
            "DeltaA2_Z": self.DeltaA2_Z,
            "DeltaA1_Z": self.DeltaA1_Z,
            "DeltaA0_Z": self.DeltaA0_Z,
        }

    @property
    def delta_coefficients(self) -> Mapping[str, Any]:
        return {key: value.Delta for key, value in self.coefficients.items()}

    @property
    def coefficient_primitives(self) -> Mapping[str, PrimitiveDelta]:
        return self.coefficients

    @property
    def root_inputs(self) -> dict[str, Any]:
        """Keyword packet consumable by a caller-owned root interface.

        Root values are intentionally not included: the caller must provide
        separate ``aP``/``aF`` values and their Z1 rows to its root routine.
        """

        packet: dict[str, Any] = {}
        for key in ("A2", "A1", "A0"):
            primitive = self.coefficients[key]
            packet[key + "P"] = primitive.P
            packet[key + "F"] = primitive.F
            packet[key + "P_Z"] = primitive.P_Z
            packet[key + "F_Z"] = primitive.F_Z
            packet["Delta" + key] = primitive.Delta
            packet["Delta" + key + "_Z"] = primitive.Delta_Z
        return packet

    @property
    def root_interface_inputs(self) -> dict[str, Any]:
        return self.root_inputs

    def root_kwargs(
        self,
        aP: Any,
        aF: Any,
        *,
        aP_Z: Any = _MISSING,
        aF_Z: Any = _MISSING,
    ) -> dict[str, Any]:
        """Add caller-owned selector roots to :attr:`root_inputs`."""

        aP0, aPz = _branch_jet(aP, aP_Z, "aP")
        aF0, aFz = _branch_jet(aF, aF_Z, "aF")
        packet = dict(self.root_inputs)
        packet.update({"aP": aP0, "aF": aF0, "aP_Z": aPz, "aF_Z": aFz})
        return packet

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "inputs": {
                key: [item.as_dict() for item in value]
                if isinstance(value, tuple)
                else value.as_dict()
                if isinstance(value, PrimitiveDelta)
                else _display(value)
                for key, value in self.inputs.items()
            },
            "R": [item.as_dict() for item in self.R],
            "r": [item.as_dict() for item in self.r],
            "t": [item.as_dict() for item in self.t],
            "u": [item.as_dict() for item in self.u],
            "v": [item.as_dict() for item in self.v],
            "coefficients": {key: value.as_dict() for key, value in self.coefficients.items()},
            "DeltaA": {
                key: _display(value) for key, value in self.delta_coefficients.items()
            },
            "DeltaA_Z": {key: _display(value) for key, value in self.z1.items()},
            "coefficient_delta_terms": _display(self.coefficient_delta_terms),
            "input_difference_residuals": _display(self.input_difference_residuals),
            "derived_difference_residuals": _display(self.derived_difference_residuals),
            "owners": _display(self.owners),
            "gates": _display(self.gates),
            "materialization": dict(self.materialization),
            "formula": FORMULA,
            "product_formula": PRODUCT_FORMULA,
            "Z1_formula": Z1_FORMULA,
            **SOURCE_VALUE_FLAGS,
        }

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "dA2": "DeltaA2",
            "dA1": "DeltaA1",
            "dA0": "DeltaA0",
            "Delta_A2": "DeltaA2",
            "Delta_A1": "DeltaA1",
            "Delta_A0": "DeltaA0",
            "DeltaA2Z": "DeltaA2_Z",
            "DeltaA1Z": "DeltaA1_Z",
            "DeltaA0Z": "DeltaA0_Z",
        }
        key = aliases.get(key, key)
        if key in ("gates", "gate"):
            return self.gates
        if key in ("root_inputs", "root_interface_inputs"):
            return self.root_inputs
        if key in self.coefficients:
            return self.coefficients[key]
        if key in ("DeltaA2", "DeltaA1", "DeltaA0", "DeltaA2_Z", "DeltaA1_Z", "DeltaA0_Z"):
            return getattr(self, key)
        return getattr(self, key)


PrimitiveSelectorDelta = SelectorPrimitiveDeltaResult
PFSelectorPrimitiveDelta = SelectorPrimitiveDeltaResult


def _field_from_W(W: Mapping[str, Any] | None, current: Any, keys: Sequence[str], label: str) -> Any:
    if current is not None:
        return current
    if W is not None:
        found = _maybe_pick(W, keys)
        if found is not _MISSING:
            return found
    raise PrimitiveDeltaInputError(f"{label} is required")


def _owner_map(owners: Mapping[str, Any] | None) -> dict[str, Any]:
    if owners is None:
        return {}
    if not isinstance(owners, Mapping):
        raise PrimitiveDeltaInputError("owners must be a mapping and are metadata only")
    return dict(owners)


def build_selector_primitive_delta(
    mu: Any = None,
    first_row: Any = None,
    divided_row: Any = None,
    determinant: Any = None,
    incoming_factors: Any = None,
    incoming: Any = None,
    scaled_main_rows: Any = None,
    fixed_pulse_energy: Any = None,
    end_energy: Any = None,
    incoming_energy: Any = None,
    future: Any = None,
    base: Any = None,
    *,
    W: Mapping[str, Any] | None = None,
    det: Any = None,
    nu: Any = None,
    K: Any = None,
    E: Any = None,
    owners: Mapping[str, Any] | None = None,
    source_kind: str = "caller_owned",
    rounded: bool = False,
) -> SelectorPrimitiveDeltaResult:
    """Build the independent signed P/F selector primitive algebra.

    ``W`` is an optional source-shaped mapping equivalent to the named
    arguments.  It is never traversed as a graph and is only read at the
    direct field keys listed below.  Every field must still be a primitive or
    two-row collection of primitives with explicit Delta and Delta_Z rows.
    """

    # Permit the natural ``build_selector_primitive_delta(W)`` spelling while
    # retaining the explicit positional primitive form below.
    if (
        W is None
        and isinstance(mu, Mapping)
        and first_row is None
        and divided_row is None
        and determinant is None
        and incoming_factors is None
    ):
        W, mu = mu, None
    if W is not None and not isinstance(W, Mapping):
        raise PrimitiveDeltaInputError("W must be a direct mapping of caller-owned primitives")
    determinant = determinant if determinant is not None else det
    fixed_pulse_energy = fixed_pulse_energy if fixed_pulse_energy is not None else K
    end_energy = end_energy if end_energy is not None else nu
    incoming_energy = incoming_energy if incoming_energy is not None else E

    raw = {
        "mu": _field_from_W(W, mu, ("mu", "MU"), "mu"),
        "first_row": _field_from_W(W, first_row, ("first_row", "A", "A_rows"), "first_row"),
        "divided_row": _field_from_W(W, divided_row, ("divided_row", "D", "D_rows"), "divided_row"),
        "determinant": _field_from_W(
            W,
            determinant,
            ("det", "determinant", "positive_determinant_absolute", "positive_determinant"),
            "determinant",
        ),
        "incoming_factors": _field_from_W(
            W, incoming_factors, ("incoming_factors", "k", "factors"), "incoming_factors"
        ),
        "incoming": _field_from_W(W, incoming, ("incoming", "I", "inlet"), "incoming"),
        "scaled_main_rows": _field_from_W(
            W, scaled_main_rows, ("scaled_main_rows", "M", "main_rows"), "scaled_main_rows"
        ),
        "fixed_pulse_energy": _field_from_W(
            W, fixed_pulse_energy, ("fixed_pulse_energy", "K", "fixed_energy"), "fixed_pulse_energy"
        ),
        "end_energy": _field_from_W(
            W, end_energy, ("exact_nonzero_end_energy", "end_energy", "nu", "ν"), "end_energy"
        ),
        "incoming_energy": _field_from_W(
            W, incoming_energy, ("incoming_energy", "E"), "incoming_energy"
        ),
        "future": _field_from_W(W, future, ("future", "future_energy"), "future"),
        "base": _field_from_W(W, base, ("base", "base_term"), "base"),
    }

    mu_p = _primitive(raw["mu"], "mu")
    A0, A1 = _vector(raw["first_row"], "first_row")
    D0, D1 = _vector(raw["divided_row"], "divided_row")
    determinant_p = _primitive(raw["determinant"], "determinant")
    k0, k1 = _vector(raw["incoming_factors"], "incoming_factors")
    I0, I1 = _vector(raw["incoming"], "incoming")
    M0, M1 = _vector(raw["scaled_main_rows"], "scaled_main_rows")
    K_p = _primitive(raw["fixed_pulse_energy"], "fixed_pulse_energy")
    nu0, nu1 = _vector(raw["end_energy"], "end_energy")
    E_p = _primitive(raw["incoming_energy"], "incoming_energy")
    future_p = _primitive(raw["future"], "future")
    base_p = _primitive(raw["base"], "base")

    R0 = signed_product(k0, I0, name="R0=k0*I0")
    R1 = signed_product(k1, I1, name="R1=k1*I1")
    r0 = signed_negate(signed_product(mu_p, R0, name="mu*R0"), name="r0=-mu*R0")
    r1 = signed_negate(signed_subtract(R1, R0, name="R1-R0"), name="r1=-(R1-R0)")
    t0 = signed_negate(signed_product(mu_p, M0, name="mu*M0"), name="t0=-mu*M0")
    t1 = signed_negate(signed_subtract(M1, M0, name="M1-M0"), name="t1=-(M1-M0)")

    def solve(
        first: tuple[PrimitiveDelta, PrimitiveDelta],
        divided: tuple[PrimitiveDelta, PrimitiveDelta],
        rhs: tuple[PrimitiveDelta, PrimitiveDelta],
        label: str,
    ) -> tuple[PrimitiveDelta, PrimitiveDelta]:
        rhs0, rhs1 = rhs
        n1 = signed_subtract(
            signed_product(first[1], rhs1, name=label + ".A1*rhs1"),
            signed_product(divided[1], rhs0, name=label + ".D1*rhs0"),
            name=label + ".N1",
        )
        n2 = signed_subtract(
            signed_product(divided[0], rhs0, name=label + ".D0*rhs0"),
            signed_product(first[0], rhs1, name=label + ".A0*rhs1"),
            name=label + ".N2",
        )
        return (
            signed_quotient(n1, determinant_p, name=label + ".u1"),
            signed_quotient(n2, determinant_p, name=label + ".u2"),
        )

    # The tuple arguments keep the displayed source equations and row order
    # explicit: u1=(A1*r1-D1*r0)/det and u2=(D0*r0-A0*r1)/det.
    u = solve((A0, A1), (D0, D1), (r0, r1), "u")
    v = solve((A0, A1), (D0, D1), (t0, t1), "v")

    v0_sq = signed_product(v[0], v[0], name="v0^2")
    v1_sq = signed_product(v[1], v[1], name="v1^2")
    nu0_v0_sq = signed_product(nu0, v0_sq, name="nu0*v0^2")
    nu1_v1_sq = signed_product(nu1, v1_sq, name="nu1*v1^2")
    A2_p = _add_many((K_p, nu0_v0_sq, nu1_v1_sq), name="A2")

    uv0 = signed_product(u[0], v[0], name="u0*v0")
    uv1 = signed_product(u[1], v[1], name="u1*v1")
    nu0_uv0 = signed_product(nu0, uv0, name="nu0*u0*v0")
    nu1_uv1 = signed_product(nu1, uv1, name="nu1*u1*v1")
    A1_p = signed_scale(2, _add_many((nu0_uv0, nu1_uv1), name="sum_nu_uv"), name="A1=2*sum_nu_uv")

    mu_E = signed_product(mu_p, E_p, name="mu*E")
    u0_sq = signed_product(u[0], u[0], name="u0^2")
    u1_sq = signed_product(u[1], u[1], name="u1^2")
    nu0_u0_sq = signed_product(nu0, u0_sq, name="nu0*u0^2")
    nu1_u1_sq = signed_product(nu1, u1_sq, name="nu1*u1^2")
    A0_p = _add_many(
        (
            signed_negate(future_p, name="-future"),
            mu_E,
            signed_negate(base_p, name="-base"),
            nu0_u0_sq,
            nu1_u1_sq,
        ),
        name="A0",
    )

    all_inputs: dict[str, Any] = {
        "mu": mu_p,
        "first_row": (A0, A1),
        "divided_row": (D0, D1),
        "determinant": determinant_p,
        "incoming_factors": (k0, k1),
        "incoming": (I0, I1),
        "scaled_main_rows": (M0, M1),
        "fixed_pulse_energy": K_p,
        "end_energy": (nu0, nu1),
        "incoming_energy": E_p,
        "future": future_p,
        "base": base_p,
    }
    owner_map = _owner_map(owners)
    input_residuals: dict[str, Any] = {}
    for key, item in all_inputs.items():
        if isinstance(item, tuple):
            input_residuals[key] = [
                {"Delta": value.difference_residual, "Delta_Z": value.difference_Z1_residual}
                for value in item
            ]
        else:
            input_residuals[key] = {
                "Delta": item.difference_residual,
                "Delta_Z": item.difference_Z1_residual,
            }
    coefficients = {"A2": A2_p, "A1": A1_p, "A0": A0_p}
    derived_residuals = {
        key: {"Delta": value.difference_residual, "Delta_Z": value.difference_Z1_residual}
        for key, value in coefficients.items()
    }

    determinant_nonzero = not (_is_definitely_zero(determinant_p.P) or _is_definitely_zero(determinant_p.F))
    gates = {
        "determinant_nonzero": _gate(
            "caller proof that determinant P and F values are nonzero",
            detail={"algebraic_nonzero_check": determinant_nonzero},
        ),
        "determinant_positive": _gate("caller proof of a positive determinant on both branches"),
        "positive_determinant": _gate("caller proof of a positive determinant margin"),
        "source_owners": _gate("caller-owned source owner binding for every primitive", detail=owner_map),
        "signed_difference_identities": _gate(
            "caller proof that every supplied Delta equals P-F",
            detail=input_residuals,
        ),
        "signed_difference_Z1_identities": _gate(
            "caller proof that every supplied Delta_Z equals P_Z-F_Z",
            detail=input_residuals,
        ),
        "derived_coefficient_difference_identities": _gate(
            "caller proof of coefficient Delta identities", detail=derived_residuals
        ),
        "ordinary_Z1_transport": _gate("caller proof of ordinary-Z1 source transport and commutation"),
        "precision": _gate("caller precision and arithmetic error analysis"),
        "rounding": _gate("caller rounding or interval error analysis"),
        "nonzero": _gate("caller proof of every required nonzero quantity"),
    }
    materialization = {
        **SOURCE_VALUE_FLAGS,
        "source_kind": source_kind,
        "source_kind_is_caller_declared": True,
        "rounded_input_declared": bool(rounded),
        "parameter_interval_midpoint_used": False,
        "P_and_F_common_input_assumed": False,
        "automatic_P_minus_F_fallback": False,
        "missing_signed_difference_rejected": True,
        "missing_signed_difference_Z1_rejected": True,
        "all_primitive_P_F_Delta_explicit": True,
        "all_primitive_Z1_explicit": True,
        "signed_product_identity_used": PRODUCT_FORMULA,
        "signed_quotient_identity_used": FORMULA,
        "all_coefficient_delta_terms_retained": True,
        "root_interface_delta_values_exposed": True,
        "root_interface_delta_Z1_exposed": True,
        "owner_arguments_recorded_only": bool(owner_map),
        "determinant_nonzero_checked_for_undefined_quotient": determinant_nonzero,
        "determinant_positive_certified": False,
        "precision_certified": False,
        "rounding_certified": False,
    }
    coefficient_delta_terms = {
        "A2": {
            "DeltaK": K_p.Delta,
            "Delta_nu0_v0_squared": nu0_v0_sq.Delta,
            "Delta_nu1_v1_squared": nu1_v1_sq.Delta,
            "expanded": A2_p.delta_terms,
        },
        "A1": {
            "factor_2_times_Delta_nu0_u0_v0": 2 * nu0_uv0.Delta,
            "factor_2_times_Delta_nu1_u1_v1": 2 * nu1_uv1.Delta,
            "expanded": A1_p.delta_terms,
        },
        "A0": {
            "minus_Delta_future": -future_p.Delta,
            "Delta_mu_E": mu_E.Delta,
            "minus_Delta_base": -base_p.Delta,
            "Delta_nu0_u0_squared": nu0_u0_sq.Delta,
            "Delta_nu1_u1_squared": nu1_u1_sq.Delta,
            "expanded": A0_p.delta_terms,
        },
    }
    return SelectorPrimitiveDeltaResult(
        inputs=all_inputs,
        R=(R0, R1),
        r=(r0, r1),
        t=(t0, t1),
        u=u,
        v=v,
        coefficients=coefficients,
        coefficient_delta_terms=coefficient_delta_terms,
        input_difference_residuals=input_residuals,
        derived_difference_residuals=derived_residuals,
        owners=owner_map,
        gates=gates,
        materialization=materialization,
    )


def build_source_shaped_selector_delta(*args: Any, **kwargs: Any) -> SelectorPrimitiveDeltaResult:
    return build_selector_primitive_delta(*args, **kwargs)


def source_shaped_selector_delta(*args: Any, **kwargs: Any) -> SelectorPrimitiveDeltaResult:
    return build_selector_primitive_delta(*args, **kwargs)


def compute_selector_primitive_delta(*args: Any, **kwargs: Any) -> SelectorPrimitiveDeltaResult:
    return build_selector_primitive_delta(*args, **kwargs)


def root_interface_inputs(result: SelectorPrimitiveDeltaResult) -> dict[str, Any]:
    """Return explicit coefficient P/F/Delta/Z1 routes for a root caller."""

    if not isinstance(result, SelectorPrimitiveDeltaResult):
        raise PrimitiveDeltaInputError("root interface requires SelectorPrimitiveDeltaResult")
    return result.root_inputs


def selector_root_inputs(result: SelectorPrimitiveDeltaResult) -> dict[str, Any]:
    return root_interface_inputs(result)


def _fixture_exact_sympy() -> dict[str, Any]:
    """Build an independent exact SymPy fixture; this function is uncalled."""

    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - optional dependency
        return {
            "schema": SCHEMA,
            "fixture": "independent_exact_sympy_pf_selector_primitive_delta",
            "checks_passed": False,
            "sympy_available": False,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            **SOURCE_VALUE_FLAGS,
        }

    z = sp.symbols("z")
    tiny = sp.Rational(1, 10**80)

    def make(name: str, f0: int, fs: int, d0: Any, ds: int) -> tuple[PrimitiveDelta, Any, Any, Any]:
        f_expr = f0 + fs * z
        d_expr = d0 + ds * z
        p_expr = f_expr + d_expr
        primitive = PrimitiveDelta(
            p_expr.subs(z, 0),
            f_expr.subs(z, 0),
            d_expr.subs(z, 0),
            sp.diff(p_expr, z).subs(z, 0),
            sp.diff(f_expr, z).subs(z, 0),
            sp.diff(d_expr, z).subs(z, 0),
            name=name,
        )
        return primitive, p_expr, f_expr, d_expr

    # Every P/F pair is distinct.  Determinant P-F is tiny but exact and all
    # primitive Z1 rows are nonzero.
    mu, muP, muF, muD = make("mu", 3, 2, tiny * 2, 3)
    A0, A0P, A0F, A0D = make("A0", 5, 2, tiny * 3, 2)
    A1, A1P, A1F, A1D = make("A1", 7, 3, tiny * 5, 3)
    D0, D0P, D0F, D0D = make("D0", 11, 2, tiny * 7, 2)
    D1, D1P, D1F, D1D = make("D1", 13, 3, tiny * 11, 3)
    det, detP, detF, detD = make("det", 17, 2, tiny, 1)
    k0, k0P, k0F, k0D = make("k0", 19, 2, tiny * 13, 2)
    k1, k1P, k1F, k1D = make("k1", 23, 3, tiny * 17, 3)
    I0, I0P, I0F, I0D = make("I0", 29, 2, tiny * 19, 2)
    I1, I1P, I1F, I1D = make("I1", 31, 3, tiny * 23, 3)
    M0, M0P, M0F, M0D = make("M0", 37, 2, tiny * 29, 2)
    M1, M1P, M1F, M1D = make("M1", 41, 3, tiny * 31, 3)
    K, KP, KF, KD = make("K", 43, 2, tiny * 37, 2)
    nu0, nu0P, nu0F, nu0D = make("nu0", 47, 3, tiny * 41, 3)
    nu1, nu1P, nu1F, nu1D = make("nu1", 53, 4, tiny * 43, 4)
    E, EP, EF, ED = make("E", 59, 2, tiny * 47, 2)
    future, futureP, futureF, futureD = make("future", 61, 3, tiny * 53, 3)
    base, baseP, baseF, baseD = make("base", 67, 2, tiny * 59, 2)

    packet = build_selector_primitive_delta(
        mu,
        (A0, A1),
        (D0, D1),
        det,
        (k0, k1),
        (I0, I1),
        (M0, M1),
        K,
        (nu0, nu1),
        E,
        future,
        base,
        source_kind="formal_symbolic",
    )

    # Reconstruct direct symbolic P/F expressions independently of the packet
    # operations, then compare both Z0 and ordinary Z1 rows.
    R0P, R0F = k0P * I0P, k0F * I0F
    R1P, R1F = k1P * I1P, k1F * I1F
    r0P, r0F = -muP * R0P, -muF * R0F
    r1P, r1F = -(R1P - R0P), -(R1F - R0F)
    t0P, t0F = -muP * M0P, -muF * M0F
    t1P, t1F = -(M1P - M0P), -(M1F - M0F)
    u0P = (A1P * r1P - D1P * r0P) / detP
    u0F = (A1F * r1F - D1F * r0F) / detF
    u1P = (D0P * r0P - A0P * r1P) / detP
    u1F = (D0F * r0F - A0F * r1F) / detF
    v0P = (A1P * t1P - D1P * t0P) / detP
    v0F = (A1F * t1F - D1F * t0F) / detF
    v1P = (D0P * t0P - A0P * t1P) / detP
    v1F = (D0F * t0F - A0F * t1F) / detF
    A2P_direct = KP + nu0P * v0P**2 + nu1P * v1P**2
    A2F_direct = KF + nu0F * v0F**2 + nu1F * v1F**2
    A1P_direct = 2 * (nu0P * u0P * v0P + nu1P * u1P * v1P)
    A1F_direct = 2 * (nu0F * u0F * v0F + nu1F * u1F * v1F)
    A0P_direct = -futureP + muP * EP - baseP + nu0P * u0P**2 + nu1P * u1P**2
    A0F_direct = -futureF + muF * EF - baseF + nu0F * u0F**2 + nu1F * u1F**2
    direct = {"A2": (A2P_direct, A2F_direct), "A1": (A1P_direct, A1F_direct), "A0": (A0P_direct, A0F_direct)}

    def at_zero(expr: Any) -> Any:
        return sp.simplify(expr.subs(z, 0))

    input_primitives = [
        item
        for value in packet.inputs.values()
        for item in (value if isinstance(value, tuple) else (value,))
        if isinstance(item, PrimitiveDelta)
    ]
    checks = {
        "all_branch_primitives_different": all(
            item.P != item.F for item in input_primitives
        ),
        "all_input_Z1_nonzero": all(
            item.P_Z != 0 and item.F_Z != 0 and item.Delta_Z != 0
            for item in input_primitives
        ),
        "tiny_endpoint_difference_retained": packet.inputs["determinant"].Delta == tiny and tiny != 0,
        "A2_symbolic_and_direct_Z0": sp.simplify(packet.A2.P - at_zero(A2P_direct)) == 0 and sp.simplify(packet.A2.F - at_zero(A2F_direct)) == 0,
        "A1_symbolic_and_direct_Z0": sp.simplify(packet.A1.P - at_zero(A1P_direct)) == 0 and sp.simplify(packet.A1.F - at_zero(A1F_direct)) == 0,
        "A0_symbolic_and_direct_Z0": sp.simplify(packet.A0.P - at_zero(A0P_direct)) == 0 and sp.simplify(packet.A0.F - at_zero(A0F_direct)) == 0,
        "A2_symbolic_and_direct_Z1": sp.simplify(packet.A2.P_Z - sp.diff(A2P_direct, z).subs(z, 0)) == 0 and sp.simplify(packet.A2.F_Z - sp.diff(A2F_direct, z).subs(z, 0)) == 0,
        "A1_symbolic_and_direct_Z1": sp.simplify(packet.A1.P_Z - sp.diff(A1P_direct, z).subs(z, 0)) == 0 and sp.simplify(packet.A1.F_Z - sp.diff(A1F_direct, z).subs(z, 0)) == 0,
        "A0_symbolic_and_direct_Z1": sp.simplify(packet.A0.P_Z - sp.diff(A0P_direct, z).subs(z, 0)) == 0 and sp.simplify(packet.A0.F_Z - sp.diff(A0F_direct, z).subs(z, 0)) == 0,
        "signed_delta_matches_direct_difference": all(
            sp.simplify(getattr(packet, key).Delta - (getattr(packet, key).P - getattr(packet, key).F)) == 0
            for key in ("A2", "A1", "A0")
        ),
        "signed_delta_Z1_matches_independent_derivative": all(
            sp.simplify(getattr(packet, key).Delta_Z -
                        sp.diff(direct[key][0] - direct[key][1], z).subs(z, 0)) == 0
            for key in ("A2", "A1", "A0")
        ),
    }

    missing_delta_rejected = False
    try:
        missing = {"P": 2, "F": 1, "P_Z": 1, "F_Z": 1, "Delta_Z": 1}
        _primitive(missing, "missing")
    except MissingSignedDifferenceError:
        missing_delta_rejected = True
    checks["missing_Delta_rejected"] = missing_delta_rejected
    missing_delta_Z_rejected = False
    try:
        _primitive({"P": 2, "F": 1, "Delta": 1, "P_Z": 1, "F_Z": 1}, "missing_Z")
    except MissingSignedDifferenceError:
        missing_delta_Z_rejected = True
    checks["missing_Delta_Z_rejected"] = missing_delta_Z_rejected
    sequence_delta_Z_rejected = False
    try:
        _primitive(((2, 1), (1, 1), 1), "missing_sequence_Z")
    except MissingSignedDifferenceError:
        sequence_delta_Z_rejected = True
    checks["missing_sequence_Delta_Z_rejected"] = sequence_delta_Z_rejected
    opaque_rejected = False
    try:
        PrimitiveDelta({"value": 2, "source_node_id": 17}, 1, 1, 1, 1, 1)
    except OpaqueSourceValueError:
        opaque_rejected = True
    checks["opaque_metadata_rejected"] = opaque_rejected
    zero_determinant_rejected = False
    try:
        signed_quotient(PrimitiveDelta(2, 1, 1, 1, 1, 0), PrimitiveDelta(0, 1, -1, 1, 1, 0))
    except PrimitiveDeltaInputError:
        zero_determinant_rejected = True
    checks["zero_determinant_rejected"] = zero_determinant_rejected

    mismatch = PrimitiveDelta(3, 1, tiny, 4, 2, tiny)
    checks["explicit_mismatch_not_replaced"] = mismatch.Delta == tiny and mismatch.difference_residual != 0
    checks["explicit_Z1_mismatch_not_replaced"] = mismatch.Delta_Z == tiny and mismatch.difference_Z1_residual != 0
    checks["gates_remain_unpassed"] = all(gate["passed"] is False for gate in packet.gates.values())
    checks["no_precision_or_rounding_certification"] = (
        packet.materialization["precision_certified"] is False
        and packet.materialization["rounding_certified"] is False
    )

    return {
        "schema": SCHEMA,
        "fixture": "independent_exact_sympy_pf_selector_primitive_delta",
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "sympy_available": True,
        "exact_arithmetic": True,
        "tiny_endpoint_difference": str(tiny),
        "DeltaA": {key: str(value) for key, value in packet.delta_coefficients.items()},
        "DeltaA_Z": {key: str(value) for key, value in packet.z1.items()},
        "error_is_certification": False,
        **SOURCE_VALUE_FLAGS,
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run the exact fixture only when explicitly requested by the caller."""

    return _fixture_exact_sympy()


fixture_exact_sympy = _fixture_exact_sympy
build_primitive_delta = primitive_delta
signed_difference_product = signed_product
signed_difference_quotient = signed_quotient
selector_primitive_delta = build_selector_primitive_delta
pf_selector_primitive_delta = build_selector_primitive_delta


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "FORMULA",
    "PRODUCT_FORMULA",
    "Z1_FORMULA",
    "SOURCE_VALUE_FLAGS",
    "PrimitiveDeltaInputError",
    "MissingSignedDifferenceError",
    "OpaqueSourceValueError",
    "OpaqueSourceReference",
    "SourceOpaqueId",
    "Z1Jet",
    "PrimitiveDelta",
    "PFPrimitive",
    "SignedPrimitive",
    "PrimitiveJet",
    "primitive_delta",
    "build_primitive_delta",
    "signed_add",
    "signed_negate",
    "signed_subtract",
    "signed_scale",
    "signed_product",
    "delta_product",
    "product_delta",
    "signed_quotient",
    "delta_quotient",
    "quotient_delta",
    "SelectorPrimitiveDeltaResult",
    "PrimitiveSelectorDelta",
    "PFSelectorPrimitiveDelta",
    "build_selector_primitive_delta",
    "build_source_shaped_selector_delta",
    "source_shaped_selector_delta",
    "compute_selector_primitive_delta",
    "selector_primitive_delta",
    "pf_selector_primitive_delta",
    "root_interface_inputs",
    "selector_root_inputs",
    "fixture_exact_sympy",
    "verify_synthetic_fixture",
]
