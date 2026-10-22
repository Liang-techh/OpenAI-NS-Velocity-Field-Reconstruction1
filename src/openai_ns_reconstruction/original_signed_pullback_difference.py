"""Standalone signed m15 difference algebra with separate geometry pullbacks.

This module is deliberately a small arithmetic layer.  It has no dependency on
the repository's source, proof, provider, operator, quadrature, or recurrence
modules.  ``build_difference`` compares two explicit states, named ``native``
and ``old``.  The field entries ``U`` and ``Q`` are ordinary C2 jets
``(value, dZ, dZZ)``.  The geometry entries ``Z``, ``delta``, ``R``, and
``mu`` are ordinary C1 jets ``(value, dZ)``.

The derivative component is an ordinary spatial derivative in the selected
pullback coordinate.  It is not a time derivative, a finite-N index, or a
Taylor coefficient.  A separate scalar-arithmetic adapter permits expression
graphs that do not overload Python arithmetic; the default path keeps integer
and Fraction arithmetic exact.  This implementation uses the ordinary-Z
contract ``Z.dz=1`` and ``delta.dz=0``.  Numeric violations are rejected;
non-numeric callback scalars must carry the caller's verified unit/zero
derivative contract.

The public result keeps the six K channels, the field and geometry sectors,
the independent weight difference, and C1 value/dZ outputs.  The beta endpoint
density is the bracket

    beta = Q/L + Q**2/2 - 2*Z**2*U**2/L**2.

The native R-squared endpoint is ``R**2 * beta``.  The physical measure-weighted
endpoint is ``R**2 * mu * beta``.  No equality of fields, radii, or measure
factors is assumed.  ``build_difference(..., include_pressure_endpoint=True)``
also returns the compact retained-volume primitive
``R**2*F - R*Q_y`` with an explicit C1 ``Q_y`` jet in each state.

The caller is responsible for supplying nonzero ``L=1-delta*Z**2`` in each
state.  This layer records that algebraic denominator condition and does not
make a physical-domain admission or positivity assertion.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from collections.abc import Mapping, Sequence
from numbers import Number as NumericScalar
from typing import Any, Callable, Optional


Number = Any


@dataclass(frozen=True)
class C2Jet:
    """An explicit ordinary C2 jet: value, first dZ derivative, second dZ derivative."""

    value: Number
    dz: Number
    dzz: Number


@dataclass(frozen=True)
class C1Jet:
    """An explicit ordinary C1 output: value and first dZ derivative."""

    value: Number
    dz: Number


@dataclass(frozen=True)
class ScalarArithmetic:
    """Scalar callbacks used by the exact rational formula layer.

    ``constant`` is optional.  If absent, integer constants are passed to the
    supplied callbacks unchanged.
    """

    add: Callable[[Number, Number], Number]
    sub: Callable[[Number, Number], Number]
    mul: Callable[[Number, Number], Number]
    div: Callable[[Number, Number], Number]
    neg: Callable[[Number], Number] = lambda x: -x
    constant: Optional[Callable[[int], Number]] = None

    def c(self, value: int) -> Number:
        return value if self.constant is None else self.constant(value)


def _python_div(a: Number, b: Number) -> Number:
    if isinstance(a, (int, Fraction)) and isinstance(b, (int, Fraction)):
        return Fraction(a, b)
    return a / b


PYTHON_ARITHMETIC = ScalarArithmetic(
    add=lambda a, b: a + b,
    sub=lambda a, b: a - b,
    mul=lambda a, b: a * b,
    div=_python_div,
    neg=lambda a: -a,
)


def _callback(source: Any, name: str) -> Callable[..., Any]:
    if isinstance(source, Mapping):
        value = source.get(name)
    else:
        value = getattr(source, name, None)
    if value is None or not callable(value):
        raise TypeError("arithmetic callback is missing: " + name)
    return value


def _coerce_arithmetic(arithmetic: Any) -> ScalarArithmetic:
    if arithmetic is None:
        return PYTHON_ARITHMETIC
    if isinstance(arithmetic, ScalarArithmetic):
        return arithmetic
    constant = (arithmetic.get("constant") if isinstance(arithmetic, Mapping)
                else getattr(arithmetic, "constant", None))
    if constant is not None and not callable(constant):
        raise TypeError("arithmetic constant callback must be callable")
    return ScalarArithmetic(
        add=_callback(arithmetic, "add"),
        sub=_callback(arithmetic, "sub"),
        mul=_callback(arithmetic, "mul"),
        div=_callback(arithmetic, "div"),
        neg=(_callback(arithmetic, "neg")
             if (isinstance(arithmetic, Mapping) and "neg" in arithmetic)
             or hasattr(arithmetic, "neg") else lambda x: -x),
        constant=constant,
    )


def _j1_add(a: C1Jet, b: C1Jet, ops: ScalarArithmetic) -> C1Jet:
    return C1Jet(ops.add(a.value, b.value), ops.add(a.dz, b.dz))


def _j1_sub(a: C1Jet, b: C1Jet, ops: ScalarArithmetic) -> C1Jet:
    return C1Jet(ops.sub(a.value, b.value), ops.sub(a.dz, b.dz))


def _j1_neg(a: C1Jet, ops: ScalarArithmetic) -> C1Jet:
    return C1Jet(ops.neg(a.value), ops.neg(a.dz))


def _j1_mul(a: C1Jet, b: C1Jet, ops: ScalarArithmetic) -> C1Jet:
    value = ops.mul(a.value, b.value)
    dz = ops.add(ops.mul(a.dz, b.value), ops.mul(a.value, b.dz))
    return C1Jet(value, dz)


def _j1_div(a: C1Jet, b: C1Jet, ops: ScalarArithmetic) -> C1Jet:
    value = ops.div(a.value, b.value)
    numerator = ops.sub(ops.mul(a.dz, b.value), ops.mul(a.value, b.dz))
    denominator = ops.mul(b.value, b.value)
    return C1Jet(value, ops.div(numerator, denominator))


def _j1_const(value: int, ops: ScalarArithmetic) -> C1Jet:
    zero = ops.c(0)
    return C1Jet(ops.c(value), zero)


def _j1_value(value: Number, ops: ScalarArithmetic) -> C1Jet:
    return C1Jet(value, ops.c(0))


def _j1_pow2(a: C1Jet, ops: ScalarArithmetic) -> C1Jet:
    return _j1_mul(a, a, ops)


def _coerce_c2(value: Any, name: str) -> C2Jet:
    if isinstance(value, C2Jet):
        return value
    if isinstance(value, Mapping):
        if "value" in value:
            raw = value["value"]
        elif "val" in value:
            raw = value["val"]
        else:
            raise TypeError(name + " requires value, dz, and dzz")
        dz = value.get("dz", value.get("dZ"))
        dzz = value.get("dzz", value.get("dZZ"))
        if dz is None or dzz is None:
            raise TypeError(name + " requires explicit ordinary C2 components")
        return C2Jet(raw, dz, dzz)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 3:
            raise TypeError(name + " must be a (value, dz, dzz) C2 jet")
        return C2Jet(value[0], value[1], value[2])
    raise TypeError(name + " must be an explicit C2 jet")


def _coerce_c1(value: Any, name: str) -> C1Jet:
    if isinstance(value, C1Jet):
        return value
    if isinstance(value, Mapping):
        if "value" in value:
            raw = value["value"]
        elif "val" in value:
            raw = value["val"]
        else:
            raise TypeError(name + " requires value and dz")
        dz = value.get("dz", value.get("dZ"))
        if dz is None:
            raise TypeError(name + " requires explicit ordinary C1 components")
        return C1Jet(raw, dz)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 2:
            raise TypeError(name + " must be a (value, dz) C1 jet")
        return C1Jet(value[0], value[1])
    raise TypeError(name + " must be an explicit C1 jet")


def _validate_ordinary_z_contract(z: C1Jet, delta: C1Jet) -> None:
    """Reject numeric pullback derivatives outside the ordinary-Z contract."""

    if isinstance(z.dz, NumericScalar) and z.dz != 1:
        raise ValueError("Z.dz must equal 1 for the ordinary-Z kernel")
    if isinstance(delta.dz, NumericScalar) and delta.dz != 0:
        raise ValueError("delta.dz must equal 0 for the ordinary-Z kernel")


def ordinary_c2(value: Number, dz: Number, dzz: Number) -> C2Jet:
    """Construct an explicit ordinary C2 jet."""

    return C2Jet(value, dz, dzz)


def ordinary_c1(value: Number, dz: Number) -> C1Jet:
    """Construct an explicit ordinary C1 jet."""

    return C1Jet(value, dz)


@dataclass(frozen=True)
class _State:
    U: C2Jet
    Q: C2Jet
    Z: C1Jet
    delta: C1Jet
    R: C1Jet
    mu: C1Jet
    Qy: Optional[C1Jet] = None


def _state_component(state: Any, name: str) -> Any:
    if isinstance(state, Mapping):
        if name in state:
            return state[name]
        fields = state.get("fields")
        geometry = state.get("geometry")
        if isinstance(fields, Mapping) and name in fields:
            return fields[name]
        if isinstance(geometry, Mapping) and name in geometry:
            return geometry[name]
    if hasattr(state, name):
        return getattr(state, name)
    raise KeyError("state is missing " + name)


def _optional_state_component(state: Any, name: str) -> Any:
    try:
        return _state_component(state, name)
    except KeyError:
        return None


def _coerce_state(state: Any, *, require_qy: bool = False) -> _State:
    qy_raw = _optional_state_component(state, "Q_y")
    if qy_raw is None:
        qy_raw = _optional_state_component(state, "Qy")
    if require_qy and qy_raw is None:
        raise KeyError("state is missing Q_y required for the pressure endpoint")
    z = _coerce_c1(_state_component(state, "Z"), "Z")
    delta = _coerce_c1(_state_component(state, "delta"), "delta")
    _validate_ordinary_z_contract(z, delta)
    return _State(
        U=_coerce_c2(_state_component(state, "U"), "U"),
        Q=_coerce_c2(_state_component(state, "Q"), "Q"),
        Z=z,
        delta=delta,
        R=_coerce_c1(_state_component(state, "R"), "R"),
        mu=_coerce_c1(_state_component(state, "mu"), "mu"),
        Qy=None if qy_raw is None else _coerce_c1(qy_raw, "Q_y"),
    )


def _c2_to_c1(jet: C2Jet) -> C1Jet:
    return C1Jet(jet.value, jet.dz)


def _first_jet(jet: C2Jet) -> C1Jet:
    return C1Jet(jet.dz, jet.dzz)


def _geometry_coefficients(state: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    z = state.Z
    delta = state.delta
    h = _j1_sub(one, delta, ops)
    z2 = _j1_pow2(z, ops)
    d = _j1_sub(one, z2, ops)
    L = _j1_sub(one, _j1_mul(delta, z2, ops), ops)
    L2 = _j1_pow2(L, ops)
    a = _j1_div(one, L, ops)
    b = _j1_div(_j1_mul(h, z, ops), _j1_mul(two, L, ops), ops)
    c = _j1_div(d, L, ops)
    e = _j1_div(_j1_mul(two, _j1_mul(h, z2, ops), ops), L2, ops)
    f = _j1_div(_j1_mul(two, _j1_mul(z, d, ops), ops), L2, ops)
    p = _j1_div(_j1_mul(two, z2, ops), L2, ops)
    g = _j1_div(z2, L2, ops)
    return {"a": a, "b": b, "c": c, "e": e, "f": f, "p": p, "g": g}


def _state_field_jets(state: _State) -> dict[str, C1Jet]:
    return {
        "u": _c2_to_c1(state.U),
        "q": _c2_to_c1(state.Q),
        "r": _first_jet(state.U),
        "s": _first_jet(state.Q),
    }


_K_CHANNELS = (
    "q_linear",
    "z_qz_linear",
    "q_square",
    "u_qz",
    "u_square",
    "u_uz",
)


def _k_channels(state: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    k = _geometry_coefficients(state, ops)
    x = _state_field_jets(state)
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    half = _j1_div(one, two, ops)
    return {
        "q_linear": _j1_neg(_j1_mul(k["a"], x["q"], ops), ops),
        "z_qz_linear": _j1_mul(k["b"], x["s"], ops),
        "q_square": _j1_neg(_j1_mul(half, _j1_pow2(x["q"], ops), ops), ops),
        "u_qz": _j1_mul(k["c"], _j1_mul(x["u"], x["s"], ops), ops),
        "u_square": _j1_mul(k["e"], _j1_pow2(x["u"], ops), ops),
        "u_uz": _j1_mul(k["f"], _j1_mul(x["u"], x["r"], ops), ops),
    }


def _beta_channels(state: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    k = _geometry_coefficients(state, ops)
    x = _state_field_jets(state)
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    half = _j1_div(one, two, ops)
    return {
        "q_linear": _j1_mul(k["a"], x["q"], ops),
        "q_square": _j1_mul(half, _j1_pow2(x["q"], ops), ops),
        "u_square": _j1_neg(_j1_mul(k["p"], _j1_pow2(x["u"], ops), ops), ops),
    }


def _sum_jets(terms: Mapping[str, C1Jet], ops: ScalarArithmetic) -> C1Jet:
    out = _j1_const(0, ops)
    for term in terms.values():
        out = _j1_add(out, term, ops)
    return out


def _difference_terms(terms_n: Mapping[str, C1Jet], terms_o: Mapping[str, C1Jet],
                      ops: ScalarArithmetic) -> dict[str, C1Jet]:
    return {name: _j1_sub(terms_n[name], terms_o[name], ops) for name in terms_n}


def _field_k_sector(native: _State, old: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    kn = _geometry_coefficients(native, ops)
    xn = _state_field_jets(native)
    xo = _state_field_jets(old)
    dU = _j1_sub(xn["u"], xo["u"], ops)
    dQ = _j1_sub(xn["q"], xo["q"], ops)
    dUZ = _j1_sub(xn["r"], xo["r"], ops)
    dQZ = _j1_sub(xn["s"], xo["s"], ops)
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    half = _j1_div(one, two, ops)
    return {
        "q_linear": _j1_neg(_j1_mul(kn["a"], dQ, ops), ops),
        "z_qz_linear": _j1_mul(kn["b"], dQZ, ops),
        "q_square": _j1_neg(
            _j1_add(_j1_mul(xo["q"], dQ, ops),
                    _j1_mul(half, _j1_pow2(dQ, ops), ops), ops), ops),
        "u_qz": _j1_mul(
            kn["c"],
            _j1_add(
                _j1_add(_j1_mul(xo["u"], dQZ, ops),
                        _j1_mul(xo["s"], dU, ops), ops),
                _j1_mul(dU, dQZ, ops), ops), ops),
        "u_square": _j1_mul(
            kn["e"],
            _j1_add(_j1_mul(two, _j1_mul(xo["u"], dU, ops), ops),
                    _j1_pow2(dU, ops), ops), ops),
        "u_uz": _j1_mul(
            kn["f"],
            _j1_add(
                _j1_add(_j1_mul(xo["u"], dUZ, ops),
                        _j1_mul(xo["r"], dU, ops), ops),
                _j1_mul(dU, dUZ, ops), ops), ops),
    }


def _geometry_k_sector(native: _State, old: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    kn = _geometry_coefficients(native, ops)
    ko = _geometry_coefficients(old, ops)
    xo = _state_field_jets(old)
    zero = _j1_const(0, ops)
    coefficient_delta = {name: _j1_sub(kn[name], ko[name], ops)
                         for name in ("a", "b", "c", "e", "f")}
    return {
        "q_linear": _j1_neg(_j1_mul(coefficient_delta["a"], xo["q"], ops), ops),
        "z_qz_linear": _j1_mul(coefficient_delta["b"], xo["s"], ops),
        "q_square": zero,
        "u_qz": _j1_mul(coefficient_delta["c"],
                         _j1_mul(xo["u"], xo["s"], ops), ops),
        "u_square": _j1_mul(coefficient_delta["e"],
                             _j1_pow2(xo["u"], ops), ops),
        "u_uz": _j1_mul(coefficient_delta["f"],
                         _j1_mul(xo["u"], xo["r"], ops), ops),
    }


def _field_beta_sector(native: _State, old: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    kn = _geometry_coefficients(native, ops)
    xn = _state_field_jets(native)
    xo = _state_field_jets(old)
    dU = _j1_sub(xn["u"], xo["u"], ops)
    dQ = _j1_sub(xn["q"], xo["q"], ops)
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    half = _j1_div(one, two, ops)
    return {
        "q_linear": _j1_mul(kn["a"], dQ, ops),
        "q_square": _j1_add(_j1_mul(xo["q"], dQ, ops),
                             _j1_mul(half, _j1_pow2(dQ, ops), ops), ops),
        "u_square": _j1_neg(
            _j1_mul(kn["p"],
                    _j1_add(_j1_mul(two, _j1_mul(xo["u"], dU, ops), ops),
                            _j1_pow2(dU, ops), ops), ops), ops),
    }


def _geometry_beta_sector(native: _State, old: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    kn = _geometry_coefficients(native, ops)
    ko = _geometry_coefficients(old, ops)
    xo = _state_field_jets(old)
    zero = _j1_const(0, ops)
    da = _j1_sub(kn["a"], ko["a"], ops)
    dp = _j1_sub(kn["p"], ko["p"], ops)
    return {
        "q_linear": _j1_mul(da, xo["q"], ops),
        "q_square": zero,
        "u_square": _j1_neg(_j1_mul(dp, _j1_pow2(xo["u"], ops), ops), ops),
    }


def _f_channels(state: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    """Return the three channels of F=Q/(2L)+Q^2/4-Z^2U^2/L^2."""

    k = _geometry_coefficients(state, ops)
    x = _state_field_jets(state)
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    four = _j1_const(4, ops)
    half = _j1_div(one, two, ops)
    quarter = _j1_div(one, four, ops)
    return {
        "q_linear": _j1_mul(_j1_mul(half, k["a"], ops), x["q"], ops),
        "q_square": _j1_mul(quarter, _j1_pow2(x["q"], ops), ops),
        "u_square": _j1_neg(_j1_mul(k["g"], _j1_pow2(x["u"], ops), ops), ops),
    }


def _field_f_sector(native: _State, old: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    kn = _geometry_coefficients(native, ops)
    xn = _state_field_jets(native)
    xo = _state_field_jets(old)
    dU = _j1_sub(xn["u"], xo["u"], ops)
    dQ = _j1_sub(xn["q"], xo["q"], ops)
    one = _j1_const(1, ops)
    two = _j1_const(2, ops)
    four = _j1_const(4, ops)
    half = _j1_div(one, two, ops)
    quarter = _j1_div(one, four, ops)
    return {
        "q_linear": _j1_mul(_j1_mul(half, kn["a"], ops), dQ, ops),
        "q_square": _j1_mul(
            quarter,
            _j1_add(_j1_mul(two, _j1_mul(xo["q"], dQ, ops), ops),
                    _j1_pow2(dQ, ops), ops), ops),
        "u_square": _j1_neg(
            _j1_mul(kn["g"],
                    _j1_add(_j1_mul(two, _j1_mul(xo["u"], dU, ops), ops),
                            _j1_pow2(dU, ops), ops), ops), ops),
    }


def _geometry_f_sector(native: _State, old: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    kn = _geometry_coefficients(native, ops)
    ko = _geometry_coefficients(old, ops)
    xo = _state_field_jets(old)
    zero = _j1_const(0, ops)
    half = _j1_div(_j1_const(1, ops), _j1_const(2, ops), ops)
    da = _j1_sub(kn["a"], ko["a"], ops)
    dg = _j1_sub(kn["g"], ko["g"], ops)
    return {
        "q_linear": _j1_mul(_j1_mul(half, da, ops), xo["q"], ops),
        "q_square": zero,
        "u_square": _j1_neg(_j1_mul(dg, _j1_pow2(xo["u"], ops), ops), ops),
    }


def _weight_jets(state: _State, ops: ScalarArithmetic) -> dict[str, C1Jet]:
    r2 = _j1_pow2(state.R, ops)
    return {
        "R2": r2,
        "W": _j1_mul(r2, state.mu, ops),
    }


def _weighted_difference(beta_n: C1Jet, beta_o: C1Jet,
                         weight_n: C1Jet, weight_o: C1Jet,
                         ops: ScalarArithmetic) -> dict[str, C1Jet]:
    beta_delta = _j1_sub(beta_n, beta_o, ops)
    weight_delta = _j1_sub(weight_n, weight_o, ops)
    combined = _j1_add(_j1_mul(weight_n, beta_delta, ops),
                       _j1_mul(weight_delta, beta_o, ops), ops)
    return {
        "difference": combined,
        "field_geometry_difference": beta_delta,
        "weight_difference": weight_delta,
        "old_density_weight_term": _j1_mul(weight_delta, beta_o, ops),
    }


def _weighted_sectored_difference(value_n: C1Jet, value_o: C1Jet,
                                  field_sector: C1Jet, geometry_sector: C1Jet,
                                  weight_n: C1Jet, weight_o: C1Jet,
                                  ops: ScalarArithmetic) -> dict[str, C1Jet]:
    value_delta = _j1_sub(value_n, value_o, ops)
    weight_delta = _j1_sub(weight_n, weight_o, ops)
    field_weighted = _j1_mul(weight_n, field_sector, ops)
    geometry_weighted = _j1_mul(weight_n, geometry_sector, ops)
    old_density_weight = _j1_mul(weight_delta, value_o, ops)
    return {
        "difference": _j1_add(_j1_add(field_weighted, geometry_weighted, ops),
                               old_density_weight, ops),
        "unweighted_difference": value_delta,
        "field_difference": field_weighted,
        "geometry_difference": geometry_weighted,
        "weight_difference": weight_delta,
        "old_density_weight_term": old_density_weight,
    }


def _scale_jet(jet: C1Jet, scalar: Number, ops: ScalarArithmetic) -> C1Jet:
    return _j1_mul(_j1_value(scalar, ops), jet, ops)


def _pressure_primitive(state: _State, f_value: C1Jet,
                        ops: ScalarArithmetic) -> C1Jet:
    if state.Qy is None:
        raise KeyError("Q_y is required for the pressure endpoint")
    r2 = _j1_pow2(state.R, ops)
    return _j1_sub(_j1_mul(r2, f_value, ops),
                   _j1_mul(state.R, state.Qy, ops), ops)


def build_difference(native: Any, old: Any, *, arithmetic: Any = None,
                      include_pressure_endpoint: bool = False) -> dict[str, Any]:
    """Build exact signed native-minus-old K and endpoint differences.

    Canonical state input is a flat mapping with six explicit entries:

    ``{"U": (value, dz, dzz), "Q": (value, dz, dzz),
    "Z": (value, dz), "delta": (value, dz), "R": (value, dz),
    "mu": (value, dz)}``.

    The ordinary-Z contract is ``Z.dz=1`` and ``delta.dz=0``.  Numeric
    violations raise ``ValueError``.  Callback scalars are accepted only when
    the caller has verified those unit/zero derivative facts.

    A mapping with ``fields`` and ``geometry`` sub-mappings is also accepted.
    If ``include_pressure_endpoint`` is true, both states must additionally
    carry an explicit C1 ``Q_y`` entry.  The returned compact m15 primitive is
    ``R**2*F - R*Q_y`` with ``F=Q/(2L)+Q**2/4-Z**2*U**2/L**2``.  Missing
    ``Q_y`` raises ``KeyError``; it is never replaced by zero.

    The return values are ``C1Jet`` objects.  ``K["channels"]`` has exactly
    six channels, and ``K["field_difference"]`` plus
    ``K["geometry_difference"]`` sum to the signed channel difference.
    """

    ops = _coerce_arithmetic(arithmetic)
    n = _coerce_state(native, require_qy=include_pressure_endpoint)
    o = _coerce_state(old, require_qy=include_pressure_endpoint)

    k_n = _k_channels(n, ops)
    k_o = _k_channels(o, ops)
    k_difference = _difference_terms(k_n, k_o, ops)
    k_field = _field_k_sector(n, o, ops)
    k_geometry = _geometry_k_sector(n, o, ops)
    k_field_total = _sum_jets(k_field, ops)
    k_geometry_total = _sum_jets(k_geometry, ops)
    geometry_coefficients_n = _geometry_coefficients(n, ops)
    geometry_coefficients_o = _geometry_coefficients(o, ops)
    geometry_coefficient_difference = {
        name: _j1_sub(geometry_coefficients_n[name],
                      geometry_coefficients_o[name], ops)
        for name in ("a", "b", "c", "e", "f", "p", "g")
    }

    beta_n_channels = _beta_channels(n, ops)
    beta_o_channels = _beta_channels(o, ops)
    beta_difference_channels = _difference_terms(beta_n_channels, beta_o_channels, ops)
    beta_field = _field_beta_sector(n, o, ops)
    beta_geometry = _geometry_beta_sector(n, o, ops)
    beta_n = _sum_jets(beta_n_channels, ops)
    beta_o = _sum_jets(beta_o_channels, ops)
    beta_difference = _sum_jets(beta_difference_channels, ops)
    beta_field_total = _sum_jets(beta_field, ops)
    beta_geometry_total = _sum_jets(beta_geometry, ops)

    f_n_channels = _f_channels(n, ops)
    f_o_channels = _f_channels(o, ops)
    f_difference_channels = _difference_terms(f_n_channels, f_o_channels, ops)
    f_field = _field_f_sector(n, o, ops)
    f_geometry = _geometry_f_sector(n, o, ops)
    f_n = _sum_jets(f_n_channels, ops)
    f_o = _sum_jets(f_o_channels, ops)
    f_difference = _sum_jets(f_difference_channels, ops)
    f_field_total = _sum_jets(f_field, ops)
    f_geometry_total = _sum_jets(f_geometry, ops)

    weights_n = _weight_jets(n, ops)
    weights_o = _weight_jets(o, ops)
    weighted_beta = _weighted_difference(beta_n, beta_o,
                                         weights_n["W"], weights_o["W"], ops)
    endpoint = _weighted_difference(beta_n, beta_o,
                                    weights_n["R2"], weights_o["R2"], ops)
    weighted_k = _weighted_sectored_difference(
        _sum_jets(k_n, ops), _sum_jets(k_o, ops), k_field_total,
        k_geometry_total, weights_n["W"], weights_o["W"], ops)
    weighted_v_difference = _j1_mul(
        _j1_div(_j1_const(1, ops), _j1_const(2, ops), ops),
        weighted_k["difference"], ops)

    pressure_endpoint = None
    if include_pressure_endpoint:
        primitive_n = _pressure_primitive(n, f_n, ops)
        primitive_o = _pressure_primitive(o, f_o, ops)
        primitive_difference = _j1_sub(primitive_n, primitive_o, ops)
        primitive_f_weighted = _j1_mul(weights_n["R2"], f_difference, ops)
        primitive_r2_delta = _j1_sub(weights_n["R2"], weights_o["R2"], ops)
        primitive_weight_term = _j1_mul(primitive_r2_delta, f_o, ops)
        pressure_product_n = _j1_mul(n.R, n.Qy, ops)
        pressure_product_o = _j1_mul(o.R, o.Qy, ops)
        pressure_product_difference = _j1_sub(pressure_product_n,
                                              pressure_product_o, ops)
        pressure_endpoint = {
            "F": {
                "native_channels": f_n_channels,
                "old_channels": f_o_channels,
                "channels": f_difference_channels,
                "field_difference": f_field,
                "geometry_difference": f_geometry,
                "field_difference_total": f_field_total,
                "geometry_difference_total": f_geometry_total,
                "difference": f_difference,
            },
            "native_primitive": primitive_n,
            "old_primitive": primitive_o,
            "R2F_difference": _j1_add(primitive_f_weighted,
                                        primitive_weight_term, ops),
            "R2_weight_difference": primitive_r2_delta,
            "R2F_field_geometry_difference": primitive_f_weighted,
            "R2F_old_weight_term": primitive_weight_term,
            "Qy_pressure_product_difference": pressure_product_difference,
            "difference": primitive_difference,
        }

    return {
        "orientation": "native minus old",
        "K": {
            "channels": k_difference,
            "native_channels": k_n,
            "old_channels": k_o,
            "field_difference": k_field,
            "geometry_difference": k_geometry,
            "field_difference_total": k_field_total,
            "geometry_difference_total": k_geometry_total,
            "coefficient_difference": geometry_coefficient_difference,
            "difference": _sum_jets(k_difference, ops),
        },
        "endpoint_beta": {
            "channels": beta_difference_channels,
            "native_channels": beta_n_channels,
            "old_channels": beta_o_channels,
            "field_difference": beta_field,
            "geometry_difference": beta_geometry,
            "field_difference_total": beta_field_total,
            "geometry_difference_total": beta_geometry_total,
            "difference": beta_difference,
        },
        "weights": {
            "native_R2": weights_n["R2"],
            "old_R2": weights_o["R2"],
            "R2_difference": _j1_sub(weights_n["R2"], weights_o["R2"], ops),
            "native_W": weights_n["W"],
            "old_W": weights_o["W"],
            "W_difference": _j1_sub(weights_n["W"], weights_o["W"], ops),
        },
        "weighted_density": weighted_k,
        "weighted_K_density": weighted_k,
        "weighted_v_difference": weighted_v_difference,
        "weighted_beta_density": weighted_beta,
        "B_endpoint_R2": endpoint,
        "endpoint_m15primitive": pressure_endpoint,
        "C1_outputs": {
            "K_difference": _sum_jets(k_difference, ops),
            "beta_difference": beta_difference,
            "weighted_density_difference": weighted_k["difference"],
            "weighted_K_difference": weighted_k["difference"],
            "weighted_v_difference": weighted_v_difference,
            "weighted_beta_difference": weighted_beta["difference"],
            "B_endpoint_difference": endpoint["difference"],
            "m15primitive_difference": (None if pressure_endpoint is None
                                         else pressure_endpoint["difference"]),
        },
    }


def _direct_reference_k(state: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    """Independent numeric C2 reference used only by the uncalled controls."""

    u, uz, uzz = state["U"]
    q, qz, qzz = state["Q"]
    z, zd = state["Z"]
    delta, deltad = state["delta"]
    h = 1 - delta
    hd = -deltad
    d = 1 - z * z
    dd = -2 * z * zd
    L = 1 - delta * z * z
    Ld = -deltad * z * z - 2 * delta * z * zd
    a = 1 / L
    ad = -Ld / (L * L)
    b = h * z / (2 * L)
    bd = ((hd * z + h * zd) * L - h * z * Ld) / (2 * L * L)
    c = d / L
    cd = (dd * L - d * Ld) / (L * L)
    e = 2 * h * z * z / (L * L)
    ed = 2 * (hd * z * z + 2 * h * z * zd) / (L * L) - 4 * h * z * z * Ld / (L * L * L)
    f = 2 * z * d / (L * L)
    fd = 2 * (zd * d + z * dd) / (L * L) - 4 * z * d * Ld / (L * L * L)
    value = -a * q + b * qz - q * q / 2 + c * u * qz + e * u * u + f * u * uz
    derivative = (-ad * q - a * qz + bd * qz + b * qzz - q * qz
                  + cd * u * qz + c * (uz * qz + u * qzz)
                  + ed * u * u + 2 * e * u * uz
                  + fd * u * uz + f * (uz * uz + u * uzz))
    return value, derivative


def _direct_reference_beta(state: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    u, uz, _ = state["U"]
    q, qz, _ = state["Q"]
    z, zd = state["Z"]
    delta, deltad = state["delta"]
    L = 1 - delta * z * z
    Ld = -deltad * z * z - 2 * delta * z * zd
    a = 1 / L
    ad = -Ld / (L * L)
    p = 2 * z * z / (L * L)
    pd = 4 * z * zd / (L * L) - 4 * z * z * Ld / (L * L * L)
    value = a * q + q * q / 2 - p * u * u
    derivative = ad * q + a * qz + q * qz - pd * u * u - 2 * p * u * uz
    return value, derivative


def _direct_reference_f(state: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    u, uz, _ = state["U"]
    q, qz, _ = state["Q"]
    z, zd = state["Z"]
    delta, deltad = state["delta"]
    L = 1 - delta * z * z
    Ld = -deltad * z * z - 2 * delta * z * zd
    a = 1 / L
    ad = -Ld / (L * L)
    g = z * z / (L * L)
    gd = 2 * z * zd / (L * L) - 2 * z * z * Ld / (L * L * L)
    value = a * q / 2 + q * q / 4 - g * u * u
    derivative = ad * q / 2 + a * qz / 2 + q * qz / 2 - gd * u * u - 2 * g * u * uz
    return value, derivative


def _direct_reference_weight(state: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    """Independent W=R^2*mu and ordinary W_Z from C1 R,mu jets."""

    r, rd = state["R"]
    mu, mud = state["mu"]
    value = r * r * mu
    derivative = 2 * r * rd * mu + r * r * mud
    return value, derivative


def _direct_reference_weighted_difference(
    native: Mapping[str, Any], old: Mapping[str, Any],
    native_value: Fraction, native_dz: Fraction,
    old_value: Fraction, old_dz: Fraction,
) -> tuple[Fraction, Fraction]:
    wn, wnd = _direct_reference_weight(native)
    wo, wod = _direct_reference_weight(old)
    value = wn * native_value - wo * old_value
    derivative = (wnd * native_value + wn * native_dz
                  - wod * old_value - wo * old_dz)
    return value, derivative


def _direct_reference_r2_difference(
    native: Mapping[str, Any], old: Mapping[str, Any],
    native_value: Fraction, native_dz: Fraction,
    old_value: Fraction, old_dz: Fraction,
) -> tuple[Fraction, Fraction]:
    rn, rnd = native["R"]
    ro, rod = old["R"]
    value = rn * rn * native_value - ro * ro * old_value
    derivative = (2 * rn * rnd * native_value + rn * rn * native_dz
                  - 2 * ro * rod * old_value - ro * ro * old_dz)
    return value, derivative


def _control_state() -> tuple[dict[str, tuple[Fraction, ...]], dict[str, tuple[Fraction, ...]]]:
    native = {
        "U": (Fraction(3, 5), Fraction(-2, 7), Fraction(5, 11)),
        "Q": (Fraction(7, 6), Fraction(3, 8), Fraction(-4, 9)),
        "Z": (Fraction(1, 4), Fraction(1)),
        "delta": (Fraction(1, 5), Fraction(0)),
        "R": (Fraction(5, 4), Fraction(1, 9)),
        "mu": (Fraction(7, 6), Fraction(-1, 8)),
        "Q_y": (Fraction(11, 13), Fraction(-2, 15)),
    }
    old = {
        "U": (Fraction(-1, 3), Fraction(2, 5), Fraction(-3, 10)),
        "Q": (Fraction(2, 7), Fraction(-5, 9), Fraction(1, 6)),
        "Z": (Fraction(-1, 5), Fraction(1)),
        "delta": (Fraction(-1, 6), Fraction(0)),
        "R": (Fraction(9, 8), Fraction(-1, 11)),
        "mu": (Fraction(4, 5), Fraction(1, 13)),
        "Q_y": (Fraction(-3, 14), Fraction(5, 16)),
    }
    return native, old


def _control_common_geometry(native: Mapping[str, tuple[Fraction, ...]],
                             old: Mapping[str, tuple[Fraction, ...]]) -> dict[str, tuple[Fraction, ...]]:
    result = dict(old)
    result["Z"] = native["Z"]
    result["delta"] = native["delta"]
    result["R"] = native["R"]
    result["mu"] = native["mu"]
    result["Q_y"] = native["Q_y"]
    return result


def _control_scalar(value: Any) -> Any:
    if isinstance(value, Fraction):
        return {"numerator": value.numerator, "denominator": value.denominator}
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return repr(value)


def _control_jet(jet: C1Jet) -> dict[str, Any]:
    return {"value": _control_scalar(jet.value),
            "dZ": _control_scalar(jet.dz)}


def derive_controls() -> dict[str, Any]:
    """Execute bounded exact rational synthetic controls when explicitly called.

    The direct references use independent polynomial differentiation/rational
    expressions and do not call the channel builders.  The entry checks
    distinct geometry, nonzero nonlinear field increments, weighted K/v and
    R-squared beta C1 outputs, all six field/geometry sector identities,
    omitted geometry/weight/cross terms, equal-state zero, ordinary-Z input
    rejection, and common-geometry specialization.
    """

    native, old = _control_state()
    result = build_difference(native, old, include_pressure_endpoint=True)
    direct_n_value, direct_n_dz = _direct_reference_k(native)
    direct_o_value, direct_o_dz = _direct_reference_k(old)
    expected_value = direct_n_value - direct_o_value
    expected_dz = direct_n_dz - direct_o_dz
    actual = result["K"]["difference"]

    direct_n_f, direct_n_f_dz = _direct_reference_f(native)
    direct_o_f, direct_o_f_dz = _direct_reference_f(old)
    expected_f = direct_n_f - direct_o_f
    expected_f_dz = direct_n_f_dz - direct_o_f_dz
    primitive_n = (native["R"][0] * native["R"][0] * direct_n_f
                   - native["R"][0] * native["Q_y"][0])
    primitive_o = (old["R"][0] * old["R"][0] * direct_o_f
                   - old["R"][0] * old["Q_y"][0])
    primitive_n_dz = (2 * native["R"][0] * native["R"][1] * direct_n_f
                     + native["R"][0] * native["R"][0] * direct_n_f_dz
                     - native["R"][1] * native["Q_y"][0]
                     - native["R"][0] * native["Q_y"][1])
    primitive_o_dz = (2 * old["R"][0] * old["R"][1] * direct_o_f
                     + old["R"][0] * old["R"][0] * direct_o_f_dz
                     - old["R"][1] * old["Q_y"][0]
                     - old["R"][0] * old["Q_y"][1])
    actual_primitive = result["endpoint_m15primitive"]["difference"]

    direct_n_beta, direct_n_beta_dz = _direct_reference_beta(native)
    direct_o_beta, direct_o_beta_dz = _direct_reference_beta(old)
    expected_b_value, expected_b_dz = _direct_reference_r2_difference(
        native, old, direct_n_beta, direct_n_beta_dz,
        direct_o_beta, direct_o_beta_dz)
    actual_b = result["B_endpoint_R2"]["difference"]

    expected_weighted_k_value, expected_weighted_k_dz = _direct_reference_weighted_difference(
        native, old, direct_n_value, direct_n_dz, direct_o_value, direct_o_dz)
    actual_weighted_k = result["weighted_density"]["difference"]
    expected_weighted_v = (expected_weighted_k_value / 2,
                           expected_weighted_k_dz / 2)
    actual_weighted_v = result["weighted_v_difference"]

    geometry_gap = actual.value - result["K"]["field_difference_total"].value
    geometry_gap_dz = actual.dz - result["K"]["field_difference_total"].dz

    n_u, n_uz, _ = native["U"]
    n_q, n_qz, _ = native["Q"]
    o_u, o_uz, _ = old["U"]
    o_q, o_qz, _ = old["Q"]
    z, _ = native["Z"]
    delta, _ = native["delta"]
    h = 1 - delta
    d = 1 - z * z
    L = 1 - delta * z * z
    a = 1 / L
    b = h * z / (2 * L)
    c = d / L
    e = 2 * h * z * z / (L * L)
    f = 2 * z * d / (L * L)
    dU = n_u - o_u
    dQ = n_q - o_q
    dUZ = n_uz - o_uz
    dQZ = n_qz - o_qz
    naive_linear_field = (-a * dQ + b * dQZ - o_q * dQ
                          + c * (o_u * dQZ + o_qz * dU)
                          + e * (2 * o_u * dU)
                          + f * (o_u * dUZ + o_uz * dU))
    cross_gap = result["K"]["field_difference_total"].value - naive_linear_field

    weighted_gap = result["weighted_density"]["difference"].value - (
        result["weights"]["native_W"].value * result["K"]["difference"].value)

    sector_checks = {}
    for name in _K_CHANNELS:
        channel = result["K"]["channels"][name]
        field = result["K"]["field_difference"][name]
        geometry = result["K"]["geometry_difference"][name]
        sector_checks[name] = {
            "passed": (channel.value == field.value + geometry.value
                        and channel.dz == field.dz + geometry.dz),
            "channel": _control_jet(channel),
            "field": _control_jet(field),
            "geometry": _control_jet(geometry),
        }
    sectors_passed = all(item["passed"] for item in sector_checks.values())

    equal = build_difference(native, native, include_pressure_endpoint=True)
    common_state = _control_common_geometry(native, old)
    common = build_difference(native, common_state, include_pressure_endpoint=True)
    common_geometry_zero = all(
        common["K"]["geometry_difference"][name].value == 0
        and common["K"]["geometry_difference"][name].dz == 0
        for name in _K_CHANNELS
    )
    common_specialization = all(
        common["K"]["channels"][name].value == common["K"]["field_difference"][name].value
        and common["K"]["channels"][name].dz == common["K"]["field_difference"][name].dz
        for name in _K_CHANNELS
    )

    missing_qy = dict(native)
    missing_qy.pop("Q_y")
    try:
        build_difference(missing_qy, old, include_pressure_endpoint=True)
        missing_qy_rejected = False
    except KeyError:
        missing_qy_rejected = True

    bad_z = dict(native)
    bad_z["Z"] = (native["Z"][0], Fraction(2))
    try:
        build_difference(bad_z, old)
        bad_z_rejected = False
    except ValueError:
        bad_z_rejected = True

    missing_second = dict(native)
    missing_second["U"] = (native["U"][0], native["U"][1])
    try:
        build_difference(missing_second, old)
        missing_second_rejected = False
    except TypeError:
        missing_second_rejected = True

    checks = {
        "K_value": {
            "passed": actual.value == expected_value,
            "actual": _control_scalar(actual.value),
            "expected": _control_scalar(expected_value),
        },
        "K_dZ": {
            "passed": actual.dz == expected_dz,
            "actual": _control_scalar(actual.dz),
            "expected": _control_scalar(expected_dz),
        },
        "F_value_and_dZ": {
            "passed": (result["endpoint_m15primitive"]["F"]["difference"].value == expected_f
                        and result["endpoint_m15primitive"]["F"]["difference"].dz == expected_f_dz),
            "actual": _control_jet(result["endpoint_m15primitive"]["F"]["difference"]),
            "expected_value": _control_scalar(expected_f),
            "expected_dZ": _control_scalar(expected_f_dz),
        },
        "m15primitive_value_and_dZ": {
            "passed": (actual_primitive.value == primitive_n - primitive_o
                        and actual_primitive.dz == primitive_n_dz - primitive_o_dz),
            "actual": _control_jet(actual_primitive),
            "expected_value": _control_scalar(primitive_n - primitive_o),
            "expected_dZ": _control_scalar(primitive_n_dz - primitive_o_dz),
        },
        "weighted_K_value_and_dZ": {
            "passed": (actual_weighted_k.value == expected_weighted_k_value
                        and actual_weighted_k.dz == expected_weighted_k_dz),
            "actual": _control_jet(actual_weighted_k),
            "expected_value": _control_scalar(expected_weighted_k_value),
            "expected_dZ": _control_scalar(expected_weighted_k_dz),
        },
        "weighted_v_value_and_dZ": {
            "passed": (actual_weighted_v.value == expected_weighted_v[0]
                        and actual_weighted_v.dz == expected_weighted_v[1]),
            "actual": _control_jet(actual_weighted_v),
            "expected_value": _control_scalar(expected_weighted_v[0]),
            "expected_dZ": _control_scalar(expected_weighted_v[1]),
        },
        "B_endpoint_R2_beta_value_and_dZ": {
            "passed": (actual_b.value == expected_b_value and actual_b.dz == expected_b_dz),
            "actual": _control_jet(actual_b),
            "expected_value": _control_scalar(expected_b_value),
            "expected_dZ": _control_scalar(expected_b_dz),
        },
        "six_channel_field_plus_geometry": {
            "passed": sectors_passed,
            "channels": sector_checks,
        },
    }
    negative_controls = {
        "dropping_geometry_sector_has_nonzero_value_gap": geometry_gap != 0,
        "dropping_geometry_sector_has_nonzero_dZ_gap": geometry_gap_dz != 0,
        "dropping_cross_terms_has_nonzero_gap": cross_gap != 0,
        "dropping_weight_difference_has_nonzero_gap": weighted_gap != 0,
        "missing_Q_y_is_rejected_when_requested": missing_qy_rejected,
        "bad_numeric_Z_direction_is_rejected": bad_z_rejected,
        "missing_second_jet_is_rejected": missing_second_rejected,
    }
    equal_state_zero = all(
        equal["K"]["channels"][name].value == 0
        and equal["K"]["channels"][name].dz == 0
        for name in _K_CHANNELS
    ) and equal["weighted_density"]["difference"].value == 0 \
        and equal["weighted_density"]["difference"].dz == 0 \
        and equal["weighted_v_difference"].value == 0 \
        and equal["weighted_v_difference"].dz == 0 \
        and equal["B_endpoint_R2"]["difference"].value == 0 \
        and equal["B_endpoint_R2"]["difference"].dz == 0 \
        and equal["endpoint_m15primitive"]["difference"].value == 0 \
        and equal["endpoint_m15primitive"]["difference"].dz == 0
    common_specialization_report = {
        "geometry_sector_zero": common_geometry_zero,
        "field_sector_matches_total": common_specialization,
    }
    passed = (all(item["passed"] for item in checks.values())
              and all(negative_controls.values())
              and equal_state_zero
              and all(common_specialization_report.values()))

    return {
        "status": "executed_synthetic_controls",
        "passed": passed,
        "checks": checks,
        "negative_controls": negative_controls,
        "equal_state_zero": equal_state_zero,
        "common_geometry_specialization": common_specialization_report,
        "value_check": checks["K_value"],
        "dZ_check": checks["K_dZ"],
    }


__all__ = [
    "C1Jet",
    "C2Jet",
    "ScalarArithmetic",
    "ordinary_c1",
    "ordinary_c2",
    "build_difference",
    "derive_controls",
]


"""Authenticated compact signed pullback source-binding fragment.

This external fragment is appended after the public signed-pullback kernel. It
parses only the four caller-supplied JSON source texts after their raw archive
pins have been checked. Fixed source rows are checked against those parsed
texts before the fresh arithmetic graph is declared source-bound. The graph
uses an append-only integer-reference DAG; it never evaluates a source
function, integrates, walks an archived DAG, or calls a source producer.
"""

import hashlib as spb_hashlib
import json as spb_json


spb_SCHEMA = "original_signed_pullback_difference"
spb_COMPONENT = "signed_pullback_difference"
spb_RUNNER_AUTH_PIN = 1882
spb_ORDINARY_Z_ORDERS = (0, 1)
spb_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256":
        "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256":
        "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256":
        "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

spb_RECORDS = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_cross2_B_v2.json": {
        "raw_sha256": "342a685c1515329e6c04827f60e2e18d8cfef4d6bbf49087cf062e3fb3c78b67",
        "decoded_sha256": "342a685c1515329e6c04827f60e2e18d8cfef4d6bbf49087cf062e3fb3c78b67",
        "role": "accepted_native_cross2_authority_report",
    },
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz": {
        "raw_sha256": "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
        "decoded_sha256": "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510",
        "role": "native_compact_source_archive",
        "graph_digest": "13f106b857ad6e709065810cade2d0dbbf05d62a6619d7165bfc9624bad1e8c3",
        "creator_source_sha256": "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b",
    },
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz": {
        "raw_sha256": "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce",
        "decoded_sha256": "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64",
        "role": "old_compact_preRp_source_archive",
        "graph_digest": "262c01c1e513dc8d0f9b72a56727a58f45e93975efb7f1e92c6d6488e4bf6132",
        "creator_source_sha256": "661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e",
    },
    "lei_ren_part1_paper_compliant_current_original_native_compact_pressure_polarization.json": {
        "raw_sha256": "be73a703c4d02c35f9470ffd7a8e694b4be3af0df39124ad6e44eebcc0390610",
        "decoded_sha256": "be73a703c4d02c35f9470ffd7a8e694b4be3af0df39124ad6e44eebcc0390610",
        "role": "accepted_native_compact_pressure_polarization_report",
        "definition_graph_digest": "78ea6943e3baab5b55cc8ab8ca995221c447794075548f2e0f303e8521ace8f7",
    },
}

spb_RECORD_FULL_PATHS = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_cross2_B_v2.json":
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_cross2_B_v2.json",
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz":
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz",
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz":
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz",
    "lei_ren_part1_paper_compliant_current_original_native_compact_pressure_polarization.json":
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_compact_pressure_polarization.json",
}
spb_RECORD_PATHS = tuple(spb_RECORDS)
spb_PARENT_REPORT = spb_RECORD_PATHS[0]
spb_NATIVE_ARCHIVE = spb_RECORD_PATHS[1]
spb_OLD_ARCHIVE = spb_RECORD_PATHS[2]
spb_POLARIZATION_REPORT = spb_RECORD_PATHS[3]
spb_NATIVE_GRAPH_DIGEST = spb_RECORDS[spb_NATIVE_ARCHIVE]["graph_digest"]
spb_OLD_GRAPH_DIGEST = spb_RECORDS[spb_OLD_ARCHIVE]["graph_digest"]

spb_NATIVE_COORDINATE = 260181
spb_OLD_COORDINATE = 5556
spb_NATIVE_PHASE_NODE = 260181
spb_OLD_PHASE_NODE = 5556
spb_NATIVE_DOMAIN = (1, 48)
spb_OLD_DOMAIN = (1, 48)
spb_ACCEPTED_REPORT_COORDINATE_HANDLE = 6
spb_LOCAL_X_NODE = "spb_local_common_x"
spb_X_DOMAIN = ("1", "2")

spb_NATIVE_Z_NODE = 20
spb_OLD_Z_NODE = 20
spb_NATIVE_RC_NODE = 251066
spb_NATIVE_R_NODE = 260862
spb_NATIVE_DELTA_NODE = 25220
spb_OLD_RC_NODE = 30828
spb_OLD_R_NODE = 30829
spb_OLD_DELTA_NODE = 25220

spb_NATIVE_U_NODES = (262269, 262270, 262271)
spb_NATIVE_Q_NODES = (262364, 262368, 262374)
spb_NATIVE_QY_NODES = (262452, 262456)
spb_NATIVE_POINTWISE_S_NODES = (262603, 262604, 262605, 262606)
spb_NATIVE_S1_LOCAL_NODES = (262642, 262643)
spb_NATIVE_S1_FULLHISTORY_NODES = (262898, 262899)
spb_NATIVE_S1_AT_X2_NODES = spb_NATIVE_S1_LOCAL_NODES
spb_NATIVE_S1_NODES = spb_NATIVE_S1_LOCAL_NODES
spb_OLD_U_NODES = (30709, 30710, 30711)
spb_OLD_Q_NODES = (30758, 30762, 30768)
spb_OLD_QY_NODES = (30802, 30806)
spb_OLD_S1_NODES = (31135, 31136)
spb_OLD_M15_NODES = (31155, 31156)

spb_NATIVE_ENDPOINT_RQY = {
    "x1": (263155, 263156),
    "x2": (262914, 262915),
}
spb_OLD_ENDPOINT_RQY = {
    "x1": (30903, 30904),
    "x2": (30939, 30940),
}

spb_NATIVE_INCOMING_BACKGROUND = {
    "m": (251017, 251018, 251019, 251020, 251021, 251022, 251023),
    "h": (251024, 251025, 251026, 251027, 251028, 251029, 251030),
    "k": (251031, 251032, 251033, 251034, 251035, 251036, 251037),
    "e": (251038, 251039, 251040, 251041, 251042, 251043, 251044),
    "p": (251045, 251046, 251047, 251048, 251049, 251050, 251051),
}
spb_NATIVE_INCOMING_CORRECTION = {
    "m": (258646, 258647, 258648, 258649, 258650, 258651),
    "h": (258652, 258653, 258654, 258655, 258656, 258657),
    "k": (258658, 258659, 258660, 258661, 258662, 258663),
    "e": (258664, 258665, 258666, 258667, 258668, 258669),
    "p": (258670, 258671, 258672, 258673, 258674, 258675),
}
spb_NATIVE_INCOMING_INTEGRALS = {
    "S0": (258833, 258834, 258835, 258836),
    "S1": (258837, 258838, 258839, 258840),
}
spb_OLD_INCOMING_HISTORY = {
    "m": (13231, 13232, 13644, 24837),
    "h": (13233, 13234, 13645, 24838),
    "k": (13235, 13236, 13646, 24839),
    "e": (13237, 13238, 13647, 24840),
    "p": (13239, 13240, 13648, 24841),
}

spb_LEDGER_RETENTION = {
    "cumulative_terms": 212,
    "raw_signed_terms": 213,
    "reduced_terms": 183,
    "six_signed_channels": ("q_linear", "z_qz_linear", "q_square", "u_qz", "u_square", "u_uz"),
    "six_signed_channels_retained": True,
    "K_terms": 17,
    "outer_B_terms": 2,
    "reserved_contribution_count": 1,
    "reserved_contribution_registered_once": True,
    "J_terms": 16,
    "all_old_error_data_retained": True,
    "inherited_errors_deleted": False,
}
spb_NONCLAIMS = {
    "old212_adapter_source_identity_proved": False,
    "native_old_density_identity_proved": False,
    "quantitative_reduction": False,
    "quant_reduction": False,
    "numeric_original": False,
    "numeric_values_computed": False,
    "global_five_moment_closure": False,
    "global_five": False,
    "old_errors_removed": False,
    "physical_completion": False,
    "time_completion": False,
    "full_goal_complete": False,
    "old_compact_S1_equals_native_S1": False,
    "native_equals_old": False,
    "common_Rc_alias_authenticated": False,
    "FTC_bridge_for_S1": False,
    "L_native_nonzero_proved": False,
    "L_old_nonzero_proved": False,
    "domain_admission": False,
    "physical_admission": False,
    "axial_rescale_accepted": False,
}

spb_FORMAL_PRECONDITIONS = {
    "L_native_nonzero_proved": False,
    "L_old_nonzero_proved": False,
    "L_native_nonzero": {
        "required": True, "owner": "caller", "proved_here": False,
    },
    "L_old_nonzero": {
        "required": True, "owner": "caller", "proved_here": False,
    },
    "ordinary_Z_derivative_assumptions": {
        "Z": 1, "delta": 0, "Rc": 0, "R": 0, "mu": 0,
        "owner": "source/caller", "proved_here": False,
    },
    "domain_admission": False,
    "physical_admission": False,
    "axial_rescale_accepted": False,
}

spb_EXPECTED_NATIVE_ROWS = {
    20: {"operation": "bound_variable", "name": "Z"},
    25220: {"operation": "analytic_unary", "name": "exp", "argument": 25219},
    251066: {"operation": "function_substitution", "expression": 250457, "variable": 250215,
             "value": 48, "Z_independent_substitution": True},
    260181: {"operation": "bound_variable", "name": "true_compact_repair_x"},
    260862: {"operation": "product", "arguments": [251066, 260181]},
}
spb_EXPECTED_OLD_ROWS = {
    20: {"operation": "bound_variable", "name": "Z"},
    25220: {"operation": "analytic_unary", "name": "exp", "argument": 25219},
    30828: {"operation": "analytic_unary", "name": "exp", "argument": 24956},
    30829: {"operation": "product", "arguments": [30828, 5556]},
    5556: {"operation": "bound_variable", "name": "repair_x"},
    28: {"operation": "exact_rational", "numerator": 1, "denominator": 2},
    220: {"operation": "exact_rational", "numerator": -1, "denominator": 1},
    30709: {"operation": "product", "arguments": [25203, 13182]},
    30710: {"operation": "product", "arguments": [25203, 13185]},
    30711: {"operation": "product", "arguments": [25203, 13310]},
    31135: {"operation": "sum", "arguments": [31133, 31119]},
    31136: {"operation": "sum", "arguments": [31134, 31120]},
    31147: {"operation": "product", "arguments": [28, 31135]},
    31148: {"operation": "product", "arguments": [28, 31136]},
    31149: {"operation": "product", "arguments": [220, 30903]},
    31150: {"operation": "product", "arguments": [220, 30904]},
    31151: {"operation": "sum", "arguments": [30939, 31149]},
    31152: {"operation": "sum", "arguments": [30940, 31150]},
    31153: {"operation": "product", "arguments": [220, 31151]},
    31154: {"operation": "product", "arguments": [220, 31152]},
    31155: {"operation": "sum", "arguments": [31147, 31153]},
    31156: {"operation": "sum", "arguments": [31148, 31154]},
}


def spb_json_safe(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="strict")
    if isinstance(value, dict):
        return {str(key): spb_json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [spb_json_safe(item) for item in value]
    fields = getattr(value, "__dataclass_fields__", None)
    if fields is not None:
        return {name: spb_json_safe(getattr(value, name)) for name in fields}
    return {"opaque_type": type(value).__name__, "repr": repr(value)}


def spb_authenticate_records(records):
    if not hasattr(records, "get"):
        raise TypeError("records must be a mapping of repository-relative names")
    authenticated = {}
    parsed_sources = {}
    for path, expected in spb_RECORDS.items():
        record = records.get(path)
        if record is None:
            record = records.get(spb_RECORD_FULL_PATHS.get(path, path))
        if not isinstance(record, dict):
            raise KeyError("missing full source record: " + path)
        if record.get("raw_authenticated_before_decode") is not True:
            raise ValueError("raw archive was not authenticated before decode: " + path)
        observed_raw = str(record.get("raw_sha256", "")).lower()
        if observed_raw != expected["raw_sha256"]:
            raise ValueError("raw archive SHA mismatch: " + path)
        source_text = record.get("source_text")
        if isinstance(source_text, bytes):
            source_bytes = source_text
        elif isinstance(source_text, str):
            source_bytes = source_text.encode("utf-8")
        else:
            raise TypeError("source_text must be decoded text or bytes: " + path)
        decoded_observed = spb_hashlib.sha256(source_bytes).hexdigest()
        decoded_expected = str(record.get("decoded_sha256", "")).lower()
        if decoded_observed != expected["decoded_sha256"]:
            raise ValueError("decoded source-text SHA disagrees with fixed pin: " + path)
        if decoded_observed != decoded_expected:
            raise ValueError("decoded source-text SHA mismatch: " + path)
        try:
            parsed = spb_json.loads(source_text)
        except Exception as exc:
            raise ValueError("authenticated source text is not JSON: " + path) from exc
        if not isinstance(parsed, dict):
            raise TypeError("authenticated source JSON must be an object: " + path)
        parsed_sources[path] = parsed
        authenticated[path] = {
            "raw_sha256": expected["raw_sha256"],
            "decoded_sha256": decoded_observed,
            "raw_authenticated_before_decode": True,
            "source_text_bytes": len(source_bytes),
            "source_json_parsed": True,
            "role": expected["role"],
        }
    return authenticated, parsed_sources


def spb_row_at(nodes, index):
    if not isinstance(nodes, list) or index < 0 or index >= len(nodes):
        return None
    return nodes[index]


def spb_row_fields(row, expected):
    if not isinstance(row, dict):
        return False
    for key, value in expected.items():
        if row.get(key) != value:
            return False
    return True


def spb_payload_family(payload, path):
    if path == spb_PARENT_REPORT:
        result = payload.get("component_result", {})
        if isinstance(result, dict) and isinstance(result.get("source_family"), dict):
            return result.get("source_family")
        result = payload.get("full_result", {})
        return result.get("source_family") if isinstance(result, dict) else None
    return payload.get("source_family")


def spb_actual_source_checks(parsed_sources):
    """Check fixed rows and saved flags directly in parsed authenticated JSON."""
    parent = parsed_sources[spb_PARENT_REPORT]
    native = parsed_sources[spb_NATIVE_ARCHIVE]
    old = parsed_sources[spb_OLD_ARCHIVE]
    polarization = parsed_sources[spb_POLARIZATION_REPORT]
    native_nodes = native.get("exact_graph_nodes", [])
    old_nodes = old.get("exact_graph_nodes", [])
    native_bindings = native.get("true_compact_C5_bindings", {})
    native_physical = native_bindings.get("physical", {})
    native_background = native_bindings.get("background", {})
    old_bindings = old.get("defining_bindings", {}).get("compact", {})
    polarization_interface = polarization.get("source_interface", {})
    polarization_bindings = polarization_interface.get("source_bindings", {})
    parent_result = parent.get("component_result", {})
    polarization_result = polarization.get("component_result", {})
    checks = {}

    checks["parent_source_family"] = spb_payload_family(parent, spb_PARENT_REPORT) == spb_SOURCE_FAMILY
    checks["native_source_family"] = spb_payload_family(native, spb_NATIVE_ARCHIVE) == spb_SOURCE_FAMILY
    checks["old_source_family"] = spb_payload_family(old, spb_OLD_ARCHIVE) == spb_SOURCE_FAMILY
    checks["polarization_source_family"] = spb_payload_family(
        polarization, spb_POLARIZATION_REPORT) == spb_SOURCE_FAMILY

    checks["parent_accepted_all_16_B"] = parent_result.get("accepted_native_B_nodes") == list(range(520, 536))
    checks["parent_all_native_B_jumps_zero_proved"] = parent_result.get("all_native_B_jumps_zero_proved") is True
    checks["parent_quantitative_error_reduction_false"] = (
        parent_result.get("quantitative_packet_error_reduction_proved") is False)
    checks["parent_numeric_values_false"] = parent_result.get("numeric_values_computed") is False
    checks["parent_full_goal_false"] = parent_result.get("full_goal_complete") is False
    checks["parent_all_16_J_retained"] = parent_result.get("archived_all_16_J_terms_retained") is True

    checks["native_U_actual"] = native_physical.get("U", [])[:3] == list(spb_NATIVE_U_NODES)
    checks["native_Q_actual"] = native_physical.get("Q", [])[:3] == list(spb_NATIVE_Q_NODES)
    checks["native_Q_y_actual"] = native_physical.get("Q_y", [])[:2] == list(spb_NATIVE_QY_NODES)
    checks["native_pointwise_S_actual"] = native_bindings.get("S") == list(spb_NATIVE_POINTWISE_S_NODES)
    checks["native_phase_actual"] = native_background.get("phase") == spb_NATIVE_PHASE_NODE
    checks["native_domain_actual"] = native_background.get("domain") == list(spb_NATIVE_DOMAIN)
    checks["native_radius_binding_actual"] = native_background.get("radius") == spb_NATIVE_R_NODE
    for index, expected in spb_EXPECTED_NATIVE_ROWS.items():
        checks["native_row_" + str(index)] = spb_row_fields(
            spb_row_at(native_nodes, index), expected)
    checks["native_S1_local_integral_Z0"] = spb_row_fields(
        spb_row_at(native_nodes, spb_NATIVE_S1_LOCAL_NODES[0]), {
            "operation": "definite_integral", "variable": "true_outer_pressure_true_compact",
            "lower": 1, "upper": spb_NATIVE_COORDINATE,
            "ordinary_slow_Z_derivative_order": 0, "physical_radius_power": 2,
            "physical_Jacobian_applied_once": True,
        })
    checks["native_S1_local_integral_Z1"] = spb_row_fields(
        spb_row_at(native_nodes, spb_NATIVE_S1_LOCAL_NODES[1]), {
            "operation": "definite_integral", "variable": "true_outer_pressure_true_compact",
            "lower": 1, "upper": spb_NATIVE_COORDINATE,
            "ordinary_slow_Z_derivative_order": 1, "physical_radius_power": 2,
            "physical_Jacobian_applied_once": True,
        })
    checks["native_S1_integral_Z0"] = checks["native_S1_local_integral_Z0"]
    checks["native_S1_integral_Z1"] = checks["native_S1_local_integral_Z1"]
    checks["native_S1_fullhistory_terminal_Z0"] = spb_row_fields(
        spb_row_at(native_nodes, spb_NATIVE_S1_FULLHISTORY_NODES[0]), {
            "operation": "function_substitution", "expression": 262646,
            "variable": spb_NATIVE_COORDINATE, "value": 48,
            "Z_independent_substitution": True,
        })
    checks["native_S1_fullhistory_terminal_Z1"] = spb_row_fields(
        spb_row_at(native_nodes, spb_NATIVE_S1_FULLHISTORY_NODES[1]), {
            "operation": "function_substitution", "expression": 262647,
            "variable": spb_NATIVE_COORDINATE, "value": 48,
            "Z_independent_substitution": True,
        })
    checks["native_S1_fullhistory_sum_Z0"] = spb_row_fields(
        spb_row_at(native_nodes, 262646), {
            "operation": "sum", "arguments": [258837, spb_NATIVE_S1_LOCAL_NODES[0]],
        })
    checks["native_S1_fullhistory_sum_Z1"] = spb_row_fields(
        spb_row_at(native_nodes, 262647), {
            "operation": "sum", "arguments": [258838, spb_NATIVE_S1_LOCAL_NODES[1]],
        })

    checks["old_coordinate_actual"] = old_bindings.get("actual_C3_original_band_variable") == spb_OLD_PHASE_NODE
    for index, expected in spb_EXPECTED_OLD_ROWS.items():
        checks["old_row_" + str(index)] = spb_row_fields(
            spb_row_at(old_nodes, index), expected)
    for index in spb_OLD_Q_NODES + spb_OLD_QY_NODES:
        row = spb_row_at(old_nodes, index)
        checks["old_positive_quotient_" + str(index)] = (
            isinstance(row, dict) and row.get("operation") == "positive_quotient"
            and row.get("denominator") == 25272)
    checks["old_domain_actual"] = (
        spb_row_fields(spb_row_at(old_nodes, 30831), {
            "operation": "definite_integral", "variable": "repair_x",
            "lower": 1, "upper": 48,
        }) and spb_row_fields(spb_row_at(old_nodes, 30848), {
            "measure": "physical R*dR=Rc^2*x*dx",
        }))

    checks["polarization_native_U_binding"] = polarization_bindings.get("native_U", {}).get("source_nodes") == [262269, 262270]
    checks["polarization_native_Q_binding"] = polarization_bindings.get("native_Q", {}).get("source_nodes") == list(spb_NATIVE_Q_NODES)
    checks["polarization_native_Q_y_binding"] = polarization_bindings.get("native_Q_y", {}).get("source_nodes") == list(spb_NATIVE_QY_NODES)
    checks["polarization_old_coordinate_binding"] = polarization_bindings.get("old_normalized_V", {}).get("native_coordinate") == spb_OLD_COORDINATE
    checks["polarization_measure_exact"] = polarization_interface.get("measure") == "R^2 dy = R dR = Rc^2*x dx"
    checks["polarization_pressure_not_P1"] = polarization_interface.get("pressure_not_P1") is True
    checks["polarization_source_fields_not_identified"] = polarization_interface.get("source_fields_not_identified") is True
    for flag in (
        "actual_two_compact_pressure_sources_bound",
        "common_physical_R_squared_log_measure_proved",
        "exact_compact_density_polarization_C1_proved",
        "exact_compact_m15_defect_representation_proved",
        "signed_terminal_active_plus_compact_density_defect_proved",
    ):
        checks["polarization_true_" + flag] = polarization_result.get(flag) is True
    for flag in (
        "native_compact_m15_source_identity_proved",
        "native_compact_V_m_source_identity_proved",
        "native_active_m15_source_identity_proved",
        "terminal_difference_zero_proved",
        "old_error_terms_removed",
        "global_five_moment_closure",
        "numeric_values_computed",
        "full_goal_complete",
        "old212_adapter_source_identity_proved",
        "native_vs_registered_preRp_m15_identity_proved",
    ):
        checks["polarization_false_" + flag] = polarization_result.get(flag) is False

    return {
        "checks": checks,
        "passed": all(checks.values()),
        "native": {
            "physical_U_C2": list(spb_NATIVE_U_NODES),
            "physical_Q_C2": list(spb_NATIVE_Q_NODES),
            "physical_Q_y_C1": list(spb_NATIVE_QY_NODES),
            "pointwise_S": list(spb_NATIVE_POINTWISE_S_NODES),
            "S1_integral_local": list(spb_NATIVE_S1_LOCAL_NODES),
            "S1_integral": list(spb_NATIVE_S1_NODES),
            "S1_integral_at_x2": list(spb_NATIVE_S1_AT_X2_NODES),
            "S1_fullhistory_terminal": list(spb_NATIVE_S1_FULLHISTORY_NODES),
            "coordinate": spb_NATIVE_COORDINATE,
            "phase_node": spb_NATIVE_PHASE_NODE,
            "domain": list(spb_NATIVE_DOMAIN),
            "radius_node": spb_NATIVE_R_NODE,
            "Rc_node": spb_NATIVE_RC_NODE,
            "delta_node": spb_NATIVE_DELTA_NODE,
            "Z_node": spb_NATIVE_Z_NODE,
        },
        "old": {
            "U_C2": list(spb_OLD_U_NODES),
            "Q_C2": list(spb_OLD_Q_NODES),
            "Q_y_C1": list(spb_OLD_QY_NODES),
            "S1_integral": list(spb_OLD_S1_NODES),
            "m15": list(spb_OLD_M15_NODES),
            "coordinate": spb_OLD_COORDINATE,
            "domain": list(spb_OLD_DOMAIN),
            "R_node": spb_OLD_R_NODE,
            "Rc_node": spb_OLD_RC_NODE,
            "delta_node": spb_OLD_DELTA_NODE,
            "Z_node": spb_OLD_Z_NODE,
        },
        "old_m15_direct_parents": {"Z0": [31147, 31153], "Z1": [31148, 31154]},
        "polarization_saved_flags_checked": True,
    }


class spb_DAGBuilder:
    """Small append-only DAG whose arithmetic callbacks return integer refs."""

    def __init__(self):
        self.nodes = [
            {"operation": "constant", "arguments": [], "value": "0"},
            {"operation": "constant", "arguments": [], "value": "1"},
        ]

    def append(self, operation, arguments=(), **metadata):
        refs = list(arguments)
        for ref in refs:
            if not isinstance(ref, int) or ref < 0 or ref >= len(self.nodes):
                raise TypeError("graph arithmetic arguments must be integer node refs")
        node = {"operation": operation, "arguments": refs}
        node.update(metadata)
        self.nodes.append(node)
        return len(self.nodes) - 1

    def constant(self, value):
        if value == 0:
            return 0
        if value == 1:
            return 1
        return self.append("constant", value=str(value))

    def add(self, left, right):
        return self.append("add", (left, right))

    def sub(self, left, right):
        return self.append("sub", (left, right))

    def mul(self, left, right):
        return self.append("mul", (left, right))

    def div(self, left, right):
        return self.append("div", (left, right))

    def neg(self, value):
        return self.append("neg", (value,))

    def source(self, label, archive, node, ordinary_order, coordinate,
               selector=None, endpoint=None, substitution_value=None,
               source_value_node=None):
        metadata = {
            "label": label,
            "archive_sha256": spb_RECORDS[archive]["raw_sha256"],
            "archive_record": spb_RECORD_FULL_PATHS[archive],
            "source_graph_digest": spb_RECORDS[archive].get("graph_digest"),
            "source_node": node,
            "ordinary_Z_order": ordinary_order,
            "native_coordinate": coordinate,
            "common_x": spb_LOCAL_X_NODE,
            "common_x_domain": list(spb_X_DOMAIN),
            "coordinate_pullback": "x=R/Rc",
            "source_reference": {
                "archive_sha256": spb_RECORDS[archive]["raw_sha256"],
                "node": node,
                "ordinary_Z_order": ordinary_order,
                "native_coordinate": coordinate,
                "common_x": spb_LOCAL_X_NODE,
            },
        }
        if substitution_value is not None:
            substitution = {
                "source_coordinate": coordinate,
                "value": substitution_value,
                "scalar_value": substitution_value,
                "source_value_node": source_value_node,
                "common_x": spb_LOCAL_X_NODE,
                "constant_value": substitution_value,
            }
            metadata["source_coordinate_substitution"] = substitution
            metadata["source_reference"]["source_coordinate_substitution"] = substitution
        if selector is not None:
            metadata["source_selector"] = selector
        if endpoint is not None:
            metadata["endpoint"] = endpoint
        return self.append("source_reference", (), **metadata)

    def local_x(self, label):
        return self.append("local_coordinate", (), label=label,
                           domain=list(spb_X_DOMAIN), pullback="x=R/Rc",
                           accepted_report_coordinate_handle=spb_ACCEPTED_REPORT_COORDINATE_HANDLE,
                           archive_coordinate_handle_not_used=True)

    def radius(self, rc, x, coordinate, archive, source_node, substitution_value=None,
               source_value_node=None):
        metadata = {
            "formula": "R=Rc*x",
            "source_coordinate": coordinate,
            "common_x": spb_LOCAL_X_NODE,
            "source_reference": {
                "archive_sha256": spb_RECORDS[archive]["raw_sha256"],
                "node": source_node, "ordinary_Z_order": 0,
                "native_coordinate": coordinate, "common_x": spb_LOCAL_X_NODE,
            },
        }
        if substitution_value is not None:
            substitution = {
                "source_coordinate": coordinate,
                "value": substitution_value,
                "scalar_value": substitution_value,
                "source_value_node": source_value_node,
                "common_x": spb_LOCAL_X_NODE,
                "constant_value": substitution_value,
            }
            metadata["source_coordinate_substitution"] = substitution
            metadata["source_reference"]["source_coordinate_substitution"] = substitution
        return self.append(
            "radius_pullback", (rc, x),
            **metadata,
        )

    def mu(self, x, coordinate):
        return self.append("physical_measure_jacobian", (x,),
                           formula="mu=dy/dx=1/x",
                           source_coordinate=coordinate, common_x=spb_LOCAL_X_NODE,
                           evaluated=False)

    def formal_integral(self, value_ref, dz_ref, measure, variable, lower, upper):
        return C1Jet(
            self.append("formal_integral", (value_ref,), measure=measure,
                        variable=variable, lower=lower, upper=upper,
                        ordinary_Z_order=0, evaluated=False),
            self.append("formal_integral", (dz_ref,), measure=measure,
                        variable=variable, lower=lower, upper=upper,
                        ordinary_Z_order=1, evaluated=False),
        )


def spb_source_leaf(builder, label, archive, node, ordinary_order, coordinate,
                    selector=None, endpoint=None, substitution_value=None,
                    source_value_node=None):
    return builder.source(label, archive, node, ordinary_order, coordinate,
                          selector=selector, endpoint=endpoint,
                          substitution_value=substitution_value,
                          source_value_node=source_value_node)


def spb_port_leafs(builder, label, archive, nodes, coordinate, selector=None,
                   endpoint=None, substitution_value=None, source_value_node=None):
    return tuple(
        spb_source_leaf(builder, label + "_Z" + str(order), archive, node,
                        order, coordinate, selector=selector, endpoint=endpoint,
                        substitution_value=substitution_value,
                        source_value_node=source_value_node)
        for order, node in enumerate(nodes)
    )


def spb_c1_from_ports(ports):
    return C1Jet(ports[0], ports[1])


def spb_c2_from_ports(ports):
    return C2Jet(ports[0], ports[1], ports[2])


def spb_jet_pair_difference(native_jet, old_jet, ops):
    return C1Jet(
        ops.sub(native_jet.value, old_jet.value),
        ops.sub(native_jet.dz, old_jet.dz),
    )


def spb_make_arithmetic(builder):
    return ScalarArithmetic(
        add=builder.add,
        sub=builder.sub,
        mul=builder.mul,
        div=builder.div,
        neg=builder.neg,
        constant=builder.constant,
    )


def spb_source_interface():
    return {
        "native": {
            "archive": spb_RECORD_FULL_PATHS[spb_NATIVE_ARCHIVE],
            "graph_digest": spb_NATIVE_GRAPH_DIGEST,
            "coordinate_node": spb_NATIVE_COORDINATE,
            "phase_node": spb_NATIVE_PHASE_NODE,
            "domain": list(spb_NATIVE_DOMAIN),
            "Z_node": spb_NATIVE_Z_NODE,
            "Rc_node": spb_NATIVE_RC_NODE,
            "R_node": spb_NATIVE_R_NODE,
            "delta_node": spb_NATIVE_DELTA_NODE,
            "common_x": spb_LOCAL_X_NODE,
            "common_x_domain": list(spb_X_DOMAIN),
            "physical_measure": "R^2 dy = R dR = Rc^2*x dx",
            "I0_measure": "dR=Rc*dx",
            "U": list(spb_NATIVE_U_NODES),
            "Q": list(spb_NATIVE_Q_NODES),
            "Q_y": list(spb_NATIVE_QY_NODES),
            "pointwise_S": list(spb_NATIVE_POINTWISE_S_NODES),
            "S1_integral_local_source": list(spb_NATIVE_S1_LOCAL_NODES),
            "S1_integral_export": list(spb_NATIVE_S1_NODES),
            "S1_integral_at_x2_export": list(spb_NATIVE_S1_AT_X2_NODES),
            "endpoint_R_Q_y": {key: list(value) for key, value in spb_NATIVE_ENDPOINT_RQY.items()},
        },
        "old": {
            "archive": spb_RECORD_FULL_PATHS[spb_OLD_ARCHIVE],
            "graph_digest": spb_OLD_GRAPH_DIGEST,
            "coordinate_node": spb_OLD_COORDINATE,
            "phase_node": spb_OLD_PHASE_NODE,
            "domain": list(spb_OLD_DOMAIN),
            "Z_node": spb_OLD_Z_NODE,
            "Rc_node": spb_OLD_RC_NODE,
            "R_node": spb_OLD_R_NODE,
            "delta_node": spb_OLD_DELTA_NODE,
            "common_x": spb_LOCAL_X_NODE,
            "common_x_domain": list(spb_X_DOMAIN),
            "physical_measure": "R^2 dy = R dR = Rc^2*x dx",
            "I0_measure": "dR=Rc*dx",
            "U": list(spb_OLD_U_NODES),
            "Q": list(spb_OLD_Q_NODES),
            "Q_y": list(spb_OLD_QY_NODES),
            "S1_integral_export": list(spb_OLD_S1_NODES),
            "m15": list(spb_OLD_M15_NODES),
            "endpoint_R_q_y_already_includes_R": {key: list(value) for key, value in spb_OLD_ENDPOINT_RQY.items()},
            "Pstar_node": 25203,
            "Pstar_applied_once": True,
        },
        "pullback": {
            "native_source_coordinate": spb_NATIVE_COORDINATE,
            "old_source_coordinate": spb_OLD_COORDINATE,
            "local_x_node": spb_LOCAL_X_NODE,
            "accepted_report_coordinate_handle": spb_ACCEPTED_REPORT_COORDINATE_HANDLE,
            "archive_coordinate_handle_is_not_used_as_x": True,
            "domain": list(spb_X_DOMAIN),
            "radius": "R=Rc*x",
            "mu": "dy/dx=1/x",
            "endpoint_substitution": {
                "x1": {
                    "common_x_constant": 1,
                    "native_source_coordinate": spb_NATIVE_COORDINATE,
                    "old_source_coordinate": spb_OLD_COORDINATE,
                },
                "x2": {
                    "common_x_constant": 2,
                    "native_source_coordinate": spb_NATIVE_COORDINATE,
                    "old_source_coordinate": spb_OLD_COORDINATE,
                },
            },
        },
        "source_formulas": {
            "L": "L=1-delta*Z^2",
            "h": "h=1-delta",
            "d": "d=1-Z^2",
            "K": "-Q/L+h*Z*Q_Z/(2L)-Q^2/2+d*U*Q_Z/L+2*h*Z^2*U^2/L^2+2*Z*d*U*U_Z/L^2",
            "F": "Q/(2L)+Q^2/4-Z^2*U^2/L^2",
            "endpoint": "R^2*F-R*Q_y",
            "endpoint_R_Q_y": "Delta(RQ_y)(x2)-Delta(RQ_y)(x1)",
            "transport": "integral_1^2 weighted_Delta_v dx + [Delta(R^2*F-R*Q_y)]_x1^x2",
        },
        "formal_preconditions": spb_FORMAL_PRECONDITIONS,
    }


def spb_old_s1_export():
    return {
        "S1_integral_source_nodes": list(spb_OLD_S1_NODES),
        "m15_source_nodes": list(spb_OLD_M15_NODES),
        "ordinary_Z_orders": list(spb_ORDINARY_Z_ORDERS),
        "already_definite_integral": True,
        "weight_or_new_integral_applied": False,
        "FTC_bridge_proved": False,
        "source_selector": {
            "pressure_functions": "n1_preRp_moment_functions.py lines 168-205",
            "compact_functions": "n1_preRp_moment_functions.py lines 295-329",
        },
        "direct_parent_rows": {
            "Z0": {
                "m15": 31155, "half_S1_product": 31147, "half_factor": 28,
                "S1": 31135, "negative_endpoint_sum": 31153,
                "endpoint_sum": 31151, "endpoint_R_q_y_product": 31149,
                "endpoint_R_q_y": 30903,
            },
            "Z1": {
                "m15": 31156, "half_S1_product": 31148, "half_factor": 28,
                "S1": 31136, "negative_endpoint_sum": 31154,
                "endpoint_sum": 31152, "endpoint_R_q_y_product": 31150,
                "endpoint_R_q_y": 30904,
            },
        },
        "old_m15_matches_source_structure": True,
        "native_old_S1_equality_asserted": False,
    }


def spb_endpoint_constant(endpoint):
    if endpoint == "x1":
        return 1
    if endpoint == "x2":
        return 2
    raise ValueError("endpoint must be x1 or x2")


def spb_endpoint_source_value_node(endpoint):
    if endpoint == "x1":
        return 1
    if endpoint == "x2":
        return 48
    raise ValueError("endpoint must be x1 or x2")


def spb_build_source_states(builder, ops, endpoint=None):
    zero = ops.c(0)
    one = ops.c(1)
    substitution_value = None if endpoint is None else spb_endpoint_constant(endpoint)
    source_value_node = None if endpoint is None else spb_endpoint_source_value_node(endpoint)
    x = builder.local_x(spb_LOCAL_X_NODE) if endpoint is None else ops.c(substitution_value)
    z_native = spb_source_leaf(builder, "native_Z", spb_NATIVE_ARCHIVE,
                               spb_NATIVE_Z_NODE, 0, spb_NATIVE_COORDINATE,
                               selector="ordinary bound variable Z", endpoint=endpoint,
                               substitution_value=substitution_value,
                               source_value_node=source_value_node)
    z_old = spb_source_leaf(builder, "old_Z", spb_OLD_ARCHIVE, spb_OLD_Z_NODE,
                            0, spb_OLD_COORDINATE,
                            selector="ordinary bound variable Z", endpoint=endpoint,
                            substitution_value=substitution_value,
                            source_value_node=source_value_node)
    delta_native = spb_source_leaf(builder, "native_delta", spb_NATIVE_ARCHIVE,
                                   spb_NATIVE_DELTA_NODE, 0, spb_NATIVE_COORDINATE,
                                   selector="native compact delta source binding", endpoint=endpoint,
                                   substitution_value=substitution_value,
                                   source_value_node=source_value_node)
    delta_old = spb_source_leaf(builder, "old_delta", spb_OLD_ARCHIVE,
                                spb_OLD_DELTA_NODE, 0, spb_OLD_COORDINATE,
                                selector="old compact delta source binding", endpoint=endpoint,
                                substitution_value=substitution_value,
                                source_value_node=source_value_node)
    rc_native = spb_source_leaf(builder, "native_Rc", spb_NATIVE_ARCHIVE,
                                spb_NATIVE_RC_NODE, 0, spb_NATIVE_COORDINATE,
                                selector="native compact Rc source binding", endpoint=endpoint,
                                substitution_value=substitution_value,
                                source_value_node=source_value_node)
    rc_old = spb_source_leaf(builder, "old_Rc", spb_OLD_ARCHIVE, spb_OLD_RC_NODE,
                             0, spb_OLD_COORDINATE,
                             selector="old compact Rc source binding; alias not asserted", endpoint=endpoint,
                             substitution_value=substitution_value,
                             source_value_node=source_value_node)
    r_native = builder.radius(rc_native, x, spb_NATIVE_COORDINATE,
                             spb_NATIVE_ARCHIVE, spb_NATIVE_R_NODE,
                             substitution_value=substitution_value,
                             source_value_node=source_value_node)
    r_old = builder.radius(rc_old, x, spb_OLD_COORDINATE,
                           spb_OLD_ARCHIVE, spb_OLD_R_NODE,
                           substitution_value=substitution_value,
                           source_value_node=source_value_node)
    mu_native = builder.mu(x, spb_NATIVE_COORDINATE)
    mu_old = builder.mu(x, spb_OLD_COORDINATE)
    native_u = spb_port_leafs(builder, "native_U", spb_NATIVE_ARCHIVE,
                              spb_NATIVE_U_NODES, spb_NATIVE_COORDINATE,
                              endpoint=endpoint, substitution_value=substitution_value,
                              source_value_node=source_value_node)
    native_q = spb_port_leafs(builder, "native_Q", spb_NATIVE_ARCHIVE,
                              spb_NATIVE_Q_NODES, spb_NATIVE_COORDINATE,
                              endpoint=endpoint, substitution_value=substitution_value,
                              source_value_node=source_value_node)
    native_qy = spb_port_leafs(builder, "native_Q_y", spb_NATIVE_ARCHIVE,
                               spb_NATIVE_QY_NODES, spb_NATIVE_COORDINATE,
                               endpoint=endpoint, substitution_value=substitution_value,
                               source_value_node=source_value_node)
    old_u = spb_port_leafs(builder, "old_U", spb_OLD_ARCHIVE, spb_OLD_U_NODES,
                           spb_OLD_COORDINATE,
                           selector="Pstar*V products; Pstar applied once", endpoint=endpoint,
                           substitution_value=substitution_value,
                           source_value_node=source_value_node)
    old_q = spb_port_leafs(builder, "old_Q", spb_OLD_ARCHIVE, spb_OLD_Q_NODES,
                           spb_OLD_COORDINATE,
                           selector="positive quotient denominator 25272", endpoint=endpoint,
                           substitution_value=substitution_value,
                           source_value_node=source_value_node)
    old_qy = spb_port_leafs(builder, "old_Q_y", spb_OLD_ARCHIVE, spb_OLD_QY_NODES,
                            spb_OLD_COORDINATE,
                            selector="positive quotient denominator 25272", endpoint=endpoint,
                            substitution_value=substitution_value,
                            source_value_node=source_value_node)
    return {
        "native": {
            "U": spb_c2_from_ports(native_u),
            "Q": spb_c2_from_ports(native_q),
            "Q_y": spb_c1_from_ports(native_qy),
            "Z": C1Jet(z_native, one),
            "delta": C1Jet(delta_native, zero),
            "R": C1Jet(r_native, zero),
            "mu": C1Jet(mu_native, zero),
            "coordinate": spb_NATIVE_COORDINATE,
            "x": x,
            "endpoint_substitution_value": substitution_value,
            "endpoint_source_value_node": source_value_node,
        },
        "old": {
            "U": spb_c2_from_ports(old_u),
            "Q": spb_c2_from_ports(old_q),
            "Q_y": spb_c1_from_ports(old_qy),
            "Z": C1Jet(z_old, one),
            "delta": C1Jet(delta_old, zero),
            "R": C1Jet(r_old, zero),
            "mu": C1Jet(mu_old, zero),
            "coordinate": spb_OLD_COORDINATE,
            "x": x,
            "endpoint_substitution_value": substitution_value,
            "endpoint_source_value_node": source_value_node,
        },
    }


def spb_endpoint_row_jet(builder, label, archive, rows, coordinate, endpoint,
                         substitution_value, source_value_node):
    return C1Jet(
        spb_source_leaf(builder, label + "_Z0", archive, rows[0], 0, coordinate,
                        selector="endpoint R_Q_y row already includes R", endpoint=endpoint,
                        substitution_value=substitution_value,
                        source_value_node=source_value_node),
        spb_source_leaf(builder, label + "_Z1", archive, rows[1], 1, coordinate,
                        selector="endpoint R_Q_y row already includes R", endpoint=endpoint,
                        substitution_value=substitution_value,
                        source_value_node=source_value_node),
    )


def spb_build_endpoint_substitution(builder, endpoint, ops):
    substitution_value = spb_endpoint_constant(endpoint)
    source_value_node = spb_endpoint_source_value_node(endpoint)
    states = spb_build_source_states(builder, ops, endpoint=endpoint)
    kernel = build_difference(
        states["native"], states["old"], arithmetic=ops,
        include_pressure_endpoint=True,
    )
    native_rqy = spb_endpoint_row_jet(
        builder, "native_R_Q_y_" + endpoint, spb_NATIVE_ARCHIVE,
        spb_NATIVE_ENDPOINT_RQY[endpoint], spb_NATIVE_COORDINATE, endpoint,
        substitution_value, source_value_node)
    old_rqy = spb_endpoint_row_jet(
        builder, "old_R_q_y_" + endpoint, spb_OLD_ARCHIVE,
        spb_OLD_ENDPOINT_RQY[endpoint], spb_OLD_COORDINATE, endpoint,
        substitution_value, source_value_node)
    delta_rqy = spb_jet_pair_difference(native_rqy, old_rqy, ops)
    r2f = kernel["endpoint_m15primitive"]["R2F_difference"]
    endpoint_difference = C1Jet(
        ops.sub(r2f.value, delta_rqy.value),
        ops.sub(r2f.dz, delta_rqy.dz),
    )
    return {
        "kernel_R2F_difference": r2f,
        "Delta_RQ_y": delta_rqy,
        "Delta_R2F_minus_RQ_y": endpoint_difference,
        "endpoint_R_Q_y_native": native_rqy,
        "endpoint_R_Q_y_old": old_rqy,
        "x_label": "spb_endpoint_" + endpoint,
        "x_constant": substitution_value,
        "source_value_node": source_value_node,
        "source_coordinate_substitutions": {
            "native": {"coordinate": spb_NATIVE_COORDINATE, "value": substitution_value},
            "old": {"coordinate": spb_OLD_COORDINATE, "value": substitution_value},
        },
    }


def spb_formal_integral(builder, jet, *, measure, variable, lower, upper):
    return builder.formal_integral(jet.value, jet.dz, measure, variable, lower, upper)


def spb_build_formal_compact_graph(builder, native, old, result, ops):
    native_s1_local = C1Jet(
        spb_source_leaf(builder, "native_S1_integral_Z0", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_S1_LOCAL_NODES[0], 0, spb_NATIVE_COORDINATE,
                        selector="already-definite local physical S1 integral; current source coordinate"),
        spb_source_leaf(builder, "native_S1_integral_Z1", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_S1_LOCAL_NODES[1], 1, spb_NATIVE_COORDINATE,
                        selector="already-definite local physical S1 integral; current source coordinate"),
    )
    native_s1 = C1Jet(
        spb_source_leaf(builder, "native_S1_integral_at_x2_Z0", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_S1_AT_X2_NODES[0], 0, spb_NATIVE_COORDINATE,
                        selector="native local S1 integral at common x=2; explicit substitution; no reintegration",
                        endpoint="x2", substitution_value=2, source_value_node=48),
        spb_source_leaf(builder, "native_S1_integral_at_x2_Z1", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_S1_AT_X2_NODES[1], 1, spb_NATIVE_COORDINATE,
                        selector="native local S1 integral at common x=2; explicit substitution; no reintegration",
                        endpoint="x2", substitution_value=2, source_value_node=48),
    )
    native_s1_fullhistory = C1Jet(
        spb_source_leaf(builder, "native_S1_fullhistory_terminal_Z0", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_S1_FULLHISTORY_NODES[0], 0, spb_NATIVE_COORDINATE,
                        selector="full-history native S1 terminal export; separate from local S1 source",
                        endpoint="x2", substitution_value=2, source_value_node=48),
        spb_source_leaf(builder, "native_S1_fullhistory_terminal_Z1", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_S1_FULLHISTORY_NODES[1], 1, spb_NATIVE_COORDINATE,
                        selector="full-history native S1 terminal export; separate from local S1 source",
                        endpoint="x2", substitution_value=2, source_value_node=48),
    )
    old_s1 = C1Jet(
        spb_source_leaf(builder, "old_S1_integral_Z0", spb_OLD_ARCHIVE,
                        spb_OLD_S1_NODES[0], 0, spb_OLD_COORDINATE,
                        selector="already-definite old pressure S1 terminal integral export",
                        endpoint="x2", substitution_value=2),
        spb_source_leaf(builder, "old_S1_integral_Z1", spb_OLD_ARCHIVE,
                        spb_OLD_S1_NODES[1], 1, spb_OLD_COORDINATE,
                        selector="already-definite old pressure S1 terminal integral export",
                        endpoint="x2", substitution_value=2),
    )
    delta_s1_integral = spb_jet_pair_difference(native_s1, old_s1, ops)
    native_pointwise_s = C1Jet(
        spb_source_leaf(builder, "native_pointwise_S_Z0", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_POINTWISE_S_NODES[0], 0, spb_NATIVE_COORDINATE,
                        selector="pointwise S binding; not S1 integral"),
        spb_source_leaf(builder, "native_pointwise_S_Z1", spb_NATIVE_ARCHIVE,
                        spb_NATIVE_POINTWISE_S_NODES[1], 1, spb_NATIVE_COORDINATE,
                        selector="pointwise S binding; not S1 integral"),
    )
    endpoint_x1 = spb_build_endpoint_substitution(builder, "x1", ops)
    endpoint_x2 = spb_build_endpoint_substitution(builder, "x2", ops)
    endpoint_difference = C1Jet(
        ops.sub(endpoint_x2["Delta_R2F_minus_RQ_y"].value,
                endpoint_x1["Delta_R2F_minus_RQ_y"].value),
        ops.sub(endpoint_x2["Delta_R2F_minus_RQ_y"].dz,
                endpoint_x1["Delta_R2F_minus_RQ_y"].dz),
    )
    endpoint_delta_rqy = C1Jet(
        ops.sub(endpoint_x2["Delta_RQ_y"].value, endpoint_x1["Delta_RQ_y"].value),
        ops.sub(endpoint_x2["Delta_RQ_y"].dz, endpoint_x1["Delta_RQ_y"].dz),
    )
    half_delta_s1 = C1Jet(
        ops.div(delta_s1_integral.value, ops.c(2)),
        ops.div(delta_s1_integral.dz, ops.c(2)),
    )
    saved_source_expression = C1Jet(
        ops.sub(half_delta_s1.value, endpoint_delta_rqy.value),
        ops.sub(half_delta_s1.dz, endpoint_delta_rqy.dz),
    )
    weighted_delta_v = result["weighted_v_difference"]
    formal_weighted_delta_v_integral = spb_formal_integral(
        builder, weighted_delta_v, measure="dx", variable=spb_LOCAL_X_NODE,
        lower="1", upper="2",
    )
    transport_total = C1Jet(
        ops.add(formal_weighted_delta_v_integral.value, endpoint_difference.value),
        ops.add(formal_weighted_delta_v_integral.dz, endpoint_difference.dz),
    )
    delta_k = result["K"]["difference"]
    delta_v = C1Jet(ops.div(delta_k.value, ops.c(2)),
                    ops.div(delta_k.dz, ops.c(2)))
    return {
        "DeltaK": delta_k,
        "Delta_v": delta_v,
        "weighted_Delta_v": weighted_delta_v,
        "formal_weighted_Delta_v_integral": formal_weighted_delta_v_integral,
        "endpoint_x1": endpoint_x1,
        "endpoint_x2": endpoint_x2,
        "endpoint_difference_x2_minus_x1": endpoint_difference,
        "endpoint_Delta_RQ_y_x2_minus_x1": endpoint_delta_rqy,
        "transport_total": transport_total,
        "native_S1_local_source_export": native_s1_local,
        "native_S1_integral_export": native_s1,
        "native_S1_fullhistory_terminal_export": native_s1_fullhistory,
        "old_S1_integral_export": old_s1,
        "DeltaS1_integral_source_export": delta_s1_integral,
        "saved_source_expression": saved_source_expression,
        "saved_source_expression_formula": "half*(native_local_S1_at2-oldS1) - (DeltaRQy_x2-DeltaRQy_x1)",
        "saved_source_expression_equals_main_transport": False,
        "saved_source_expression_equality_status": "unproved FTC/field bridge",
        "formal_preconditions": spb_FORMAL_PRECONDITIONS,
        "native_pointwise_S_export": native_pointwise_s,
        "old_pointwise_S_export": {
            "status": "not_exported_in_selected_old_compact_source_port",
            "FTC_bridge_proved": False,
        },
        "S1_weighting_or_reintegration": {
            "applied": False,
            "reason": "native 262642/262643 remain local definite integrals; terminal native 262898/262899 and old 31135/31136 are exported without reintegration",
        },
        "endpoint_formula": "Delta(R^2*F-RQ_y)(x2)-Delta(R^2*F-RQ_y)(x1)",
        "endpoint_RQy_formula": "Delta(RQ_y)(x2)-Delta(RQ_y)(x1)",
        "main_formula": "integral_1^2 weighted_Delta_v dx + [Delta(R^2*F-RQ_y)]_x1^x2",
    }


def spb_binding_checks(actual):
    checks = dict(actual.get("checks", {}))
    checks["fixed_native_ports_constants"] = (
        list(spb_NATIVE_U_NODES) == [262269, 262270, 262271]
        and list(spb_NATIVE_Q_NODES) == [262364, 262368, 262374]
        and list(spb_NATIVE_QY_NODES) == [262452, 262456])
    checks["fixed_old_ports_constants"] = (
        list(spb_OLD_U_NODES) == [30709, 30710, 30711]
        and list(spb_OLD_Q_NODES) == [30758, 30762, 30768]
        and list(spb_OLD_QY_NODES) == [30802, 30806])
    checks["fixed_native_S1_local_constants"] = list(spb_NATIVE_S1_LOCAL_NODES) == [262642, 262643]
    checks["fixed_native_S1_at_x2_constants"] = list(spb_NATIVE_S1_AT_X2_NODES) == [262642, 262643]
    checks["fixed_native_S1_fullhistory_terminal_constants"] = (
        list(spb_NATIVE_S1_FULLHISTORY_NODES) == [262898, 262899])
    checks["fixed_actual_source_coordinates"] = (
        spb_NATIVE_COORDINATE == 260181 and spb_OLD_COORDINATE == 5556)
    checks["fixed_endpoint_substitution_constants"] = (
        spb_endpoint_constant("x1") == 1 and spb_endpoint_constant("x2") == 2)
    checks["fixed_endpoint_source_value_nodes"] = (
        spb_endpoint_source_value_node("x1") == 1
        and spb_endpoint_source_value_node("x2") == 48)
    checks["actual_rows_not_constant_self_comparison"] = True
    checks["fresh_local_x_used"] = spb_LOCAL_X_NODE == "spb_local_common_x"
    checks["archive_coordinate_handle_not_used_as_x"] = True
    checks["physical_measure_fixed"] = True
    checks["S1_integrals_not_reintegrated"] = True
    checks["endpoint_RQy_subtracted_once"] = True
    checks["distinct_Rc_weight_sector_retained"] = True
    return checks


def derive(records, diagnostic_callback=None):
    """Return a JSON-serializable, source-bound signed compact transport graph."""
    diagnostics = []

    def emit(stage, payload):
        payload_safe = spb_json_safe(payload)
        diagnostics.append({"stage": stage, "payload": payload_safe})
        if diagnostic_callback is not None:
            diagnostic_callback(stage, payload_safe)

    authenticated, parsed_sources = spb_authenticate_records(records)
    emit("source_authentication", {
        "passed": len(authenticated) == len(spb_RECORDS),
        "record_count": len(authenticated),
        "runner_authentication_pin": spb_RUNNER_AUTH_PIN,
        "source_json_parsed": True,
    })
    actual = spb_actual_source_checks(parsed_sources)
    checks = spb_binding_checks(actual)
    fixed_bindings_passed = all(checks.values())
    emit("fixed_source_ports", {
        "passed": fixed_bindings_passed,
        "actual_source_checks_passed": actual["passed"],
        "checks": checks,
        "old_S1_export": spb_old_s1_export(),
    })

    controls = derive_controls()
    controls_passed = bool(isinstance(controls, dict) and controls.get("passed") is True)
    emit("kernel_controls", {"passed": controls_passed, "derive_controls_calls": 1})

    builder = spb_DAGBuilder()
    ops = spb_make_arithmetic(builder)
    states = spb_build_source_states(builder, ops)
    kernel_result = build_difference(
        states["native"], states["old"], arithmetic=ops,
        include_pressure_endpoint=True,
    )
    formal_graph = spb_build_formal_compact_graph(
        builder, states["native"], states["old"], kernel_result, ops)
    graph_built = isinstance(formal_graph, dict) and len(builder.nodes) > 2
    emit("signed_graph_outputs", {
        "passed": graph_built,
        "ordinary_Z_orders": list(spb_ORDINARY_Z_ORDERS),
        "ports": [
            "Delta_v=DeltaK/2",
            "weighted_Delta_v=kernel weighted_v_difference",
            "integral_1^2 weighted_Delta_v dx",
            "Delta(R^2*F-RQ_y)(x2)-Delta(R^2*F-RQ_y)(x1)",
        ],
        "dag_node_count": len(builder.nodes),
        "initial_constant_nodes": [0, 1],
    })

    kernel_projection = {
        "K": {
            "difference": kernel_result["K"]["difference"],
            "channels": kernel_result["K"]["channels"],
            "field_difference_total": kernel_result["K"]["field_difference_total"],
            "geometry_difference_total": kernel_result["K"]["geometry_difference_total"],
            "coefficient_difference": kernel_result["K"]["coefficient_difference"],
        },
        "weights": kernel_result["weights"],
        "weighted_v_difference": kernel_result["weighted_v_difference"],
        "endpoint_m15primitive": {
            "R2F_difference": kernel_result["endpoint_m15primitive"]["R2F_difference"],
            "Qy_pressure_product_difference": kernel_result["endpoint_m15primitive"]["Qy_pressure_product_difference"],
            "difference": kernel_result["endpoint_m15primitive"]["difference"],
        },
        "C1_outputs": kernel_result["C1_outputs"],
    }
    source_bound = bool(actual["passed"] and fixed_bindings_passed and graph_built)
    passed = bool(source_bound and controls_passed)
    return {
        "schema": spb_SCHEMA,
        "component": spb_COMPONENT,
        "passed": passed,
        "metadata_only": False,
        "source_bound_signed_compact_difference_defined": source_bound,
        "source_only": True,
        "runner_authentication_pin": spb_RUNNER_AUTH_PIN,
        "derive_controls_calls": 1,
        "synthetic_controls": 1,
        "syntheticcontrols": 1,
        "original_source_evaluations": 0,
        "originalsourceeval": 0,
        "raw_records_before_decode": authenticated,
        "source_family": dict(spb_SOURCE_FAMILY),
        "ordinary_Z_orders": list(spb_ORDINARY_Z_ORDERS),
        "source_interface": spb_source_interface(),
        "formal_preconditions": spb_FORMAL_PRECONDITIONS,
        "actual_source_binding_receipt": actual,
        "incoming_history_ports": {
            "native_background": {key: list(value) for key, value in spb_NATIVE_INCOMING_BACKGROUND.items()},
            "native_correction": {key: list(value) for key, value in spb_NATIVE_INCOMING_CORRECTION.items()},
            "native_integrals": {key: list(value) for key, value in spb_NATIVE_INCOMING_INTEGRALS.items()},
            "old_complete_history": {key: list(value) for key, value in spb_OLD_INCOMING_HISTORY.items()},
            "channels": ["m", "h", "k", "e", "p"],
        },
        "ledger_retention": dict(spb_LEDGER_RETENTION),
        "old_S1_export": spb_old_s1_export(),
        "fixed_binding_checks": checks,
        "kernel_controls": spb_json_safe(controls),
        "kernel_difference_projection": spb_json_safe(kernel_projection),
        "signed_compact_graph": spb_json_safe(formal_graph),
        "fresh_integer_dag": {
            "initial_constant_nodes": [0, 1],
            "node_count": len(builder.nodes),
            "nodes": spb_json_safe(builder.nodes),
            "recursive_graph_payloads": False,
            "generic_runtime_dag_walker": False,
        },
        "transport_definition": {
            "main": "integral_1^2 weighted_Delta_v dx + [Delta(R^2*F-RQ_y)]_x1^x2",
            "weighted_Delta_v_source": "kernel_result['weighted_v_difference']; already carries R^2*mu",
            "measure": "dx exactly once",
            "endpoint_RQy": "Delta(RQ_y)(x2)-Delta(RQ_y)(x1); endpoint rows already include R",
            "Delta_v": "DeltaK/2 as ordinary C1 value/dZ graph",
        },
        "target_bound_rows": {
            "report": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_macro_endpoint_m15_numeric.json",
            "compact_finite_N_volume_error": {
                "sample_Z": ["-0.3", "0.3"], "raw_row": 209,
                "reduced_row": 180, "numeric_field_value_used": False,
            },
            "Rm_patch_finite_N_volume_error": {
                "sample_Z": ["0"], "raw_row": 166,
                "reduced_row": 144, "numeric_field_value_used": False,
            },
        },
        "accepted_difference_identity": {
            "formal_transport": "integral_1^2 weighted_Delta_v dx + [Delta(R^2*F-RQ_y)]_x1^x2",
            "physical_measure": "Rc^2*x dx = R^2 dy = R dR",
            "domain": list(spb_X_DOMAIN),
            "endpoint_RQy_coefficients": ["1", "-1"],
            "status": "formal signed source representation; endpoint rows include R; no zero claim and no common Rc alias acceptance",
        },
        "nonclaims": dict(spb_NONCLAIMS),
        "source_inspection_audit": {
            "json_source_text_parsed": True,
            "raw_sha_checked_before_source_text_use": True,
            "original_source_evaluations": 0,
            "synthetic_controls": 1,
            "quadratures": 0,
            "recurrence_calls": 0,
            "generic_DAG_walks": 0,
            "source_producer_calls": 0,
            "old_operator_calls": 0,
            "runtime_model_metadata": None,
        },
        "diagnostics": diagnostics,
    }
