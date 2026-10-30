"""Non-common P/F ordinary-``Z`` selector and moment jets.

This module is a caller-owned algebra boundary.  It accepts scalar values or
exact symbolic scalars and never reads a saved graph or follows a source edge.
Opaque source-reference mappings are rejected before arithmetic; a bare
integer has no provenance at this boundary and therefore cannot be classified
as either a caller scalar or a saved node id.  The rows in every input
jet are *ordinary derivatives*, not Taylor coefficients:

``(ab)_j = sum(k=0..j) binom(j,k) a_k b_(j-k)``.

The public builder keeps the four non-normalized rows ``Z0``--``Z3`` for
``c`` and ``m``.  P and F have independent incoming rows and independent
first-row weights.  It also emits both a direct difference and the conditional
correlated paths

``delta_c = w_P delta_C + delta_W C_F``

and

``delta_m = -(u_P delta_c + delta_u c_F)``.

The conditional paths are algebraic identities for the supplied caller data,
but source-function equality and source-owned difference consistency remain
explicitly unproved gates.  ``delta_P`` and ``delta_F`` are handled separately
in the normalized ``Q`` recovery; no derivative is silently replaced by zero.

The optional physical helpers accept independent ``T_P(Z)`` and ``T_F(Z)``
and retain the source-shaped products

``R_P^2 DeltaS + (R_P-R_F)(R_P+R_F) S_F``

and

``R_P DeltaQy + (R_P-R_F) QyF``.

The exact SymPy fixture is independent caller data only.  It is not run on
import and is not a numerical or source certificate.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import comb
from typing import Any


SCHEMA = "pf_noncommon_moment_jets_v1"
ORDERS = (0, 1, 2, 3)
Q_ORDERS = (0, 1, 2, 3)
PHYSICAL_ORDERS = (0, 1)


class PFNoncommonJetError(ValueError):
    """Base error for an invalid caller-owned jet or an impossible quotient."""


class MissingJetRowsError(PFNoncommonJetError):
    """Raised when an ordinary derivative row required by a formula is absent."""


class OpaqueNodeValueError(PFNoncommonJetError):
    """Raised when an opaque source reference is supplied where a scalar is needed."""


def _reject_opaque_scalar(value: Any, name: str) -> None:
    """Reject metadata/reference objects before any arithmetic is attempted.

    A plain integer is accepted as a caller scalar: this module has no source
    graph and therefore cannot infer whether an integer came from a graph.
    Saved reference dictionaries and reference-like objects are rejected
    explicitly, which prevents the common metadata-as-value error at this
    boundary but cannot classify bare integer provenance.
    """

    if isinstance(value, Mapping):
        keys = set(value)
        if keys.intersection(
            {
                "source_node_id",
                "node_id",
                "source_role",
                "namespace",
                "owner",
                "kind",
                "first_row_weight",
            }
        ):
            raise OpaqueNodeValueError(
                f"{name} is an opaque source reference, not a scalar value"
            )
        raise OpaqueNodeValueError(
            f"{name} is a mapping; source metadata cannot be used as a scalar"
        )
    if isinstance(value, (str, bytes, bytearray)):
        raise OpaqueNodeValueError(f"{name} is text/reference metadata, not a scalar")
    if value is None:
        raise PFNoncommonJetError(f"{name} is missing")
    for attribute in ("source_node_id", "node_id", "source_role"):
        if hasattr(value, attribute):
            raise OpaqueNodeValueError(
                f"{name} carries {attribute}, not a caller scalar value"
            )


def _jet(value: Any, name: str, minimum: int) -> tuple[Any, ...]:
    """Validate and preserve a caller ordinary-derivative jet."""

    if isinstance(value, (str, bytes, bytearray, Mapping)):
        _reject_opaque_scalar(value, name)
    try:
        values = tuple(value)
    except TypeError as exc:
        raise PFNoncommonJetError(
            f"{name} must be a sequence of ordinary derivative rows"
        ) from exc
    if len(values) < minimum:
        raise MissingJetRowsError(
            f"{name} has {len(values)} rows; at least {minimum} are required"
        )
    for index, item in enumerate(values):
        _reject_opaque_scalar(item, f"{name}[{index}]")
    return values


def _order(order: int, *, maximum: int, name: str = "order") -> int:
    if isinstance(order, bool) or not isinstance(order, int):
        raise PFNoncommonJetError(f"{name} must be an integer")
    if order < 0 or order > maximum:
        raise PFNoncommonJetError(f"{name} must lie in 0..{maximum}")
    return order


def _add(left: Sequence[Any], right: Sequence[Any], order: int) -> tuple[Any, ...]:
    return tuple(left[index] + right[index] for index in range(order + 1))


def _subtract(left: Sequence[Any], right: Sequence[Any], order: int) -> tuple[Any, ...]:
    return tuple(left[index] - right[index] for index in range(order + 1))


def _negate(value: Sequence[Any], order: int) -> tuple[Any, ...]:
    return tuple(-value[index] for index in range(order + 1))


def ordinary_product_jets(
    left: Sequence[Any],
    right: Sequence[Any],
    *,
    order: int = 3,
) -> tuple[Any, ...]:
    """Return ordinary Leibniz rows of a product through ``order``.

    No conversion to ``float`` is performed, so ``Fraction``, ``Decimal`` and
    exact SymPy values remain exact.  The helper is intentionally public: all
    coefficient, moment, Q and physical products use this same operation.
    """

    order = _order(order, maximum=4)
    left_rows = _jet(left, "left", order + 1)
    right_rows = _jet(right, "right", order + 1)
    return tuple(
        sum(
            comb(index, split)
            * left_rows[split]
            * right_rows[index - split]
            for split in range(index + 1)
        )
        for index in range(order + 1)
    )


# Short names are useful to callers porting the displayed product rule.
product_jets = ordinary_product_jets
ordinary_product_jet = ordinary_product_jets


def _difference(left: Sequence[Any], right: Sequence[Any], order: int) -> tuple[Any, ...]:
    return _subtract(left, right, order)


def _constant_jet(value: Any, order: int) -> tuple[Any, ...]:
    _reject_opaque_scalar(value, "constant")
    return (value,) + (value * 0,) * order


def _z_jet(z: Any, order: int) -> tuple[Any, ...]:
    _reject_opaque_scalar(z, "z")
    zero = z * 0
    if order == 0:
        return (z,)
    return (z, z * 0 + 1) + (zero,) * (order - 1)


def _inverse_jets(denominator: Sequence[Any], order: int, name: str) -> tuple[Any, ...]:
    """Return ordinary rows of ``1/denominator`` by the product recurrence."""

    rows = _jet(denominator, name, order + 1)
    first = rows[0]
    is_zero = getattr(first, "is_zero", None)
    if first == 0 or is_zero is True:
        raise PFNoncommonJetError(f"{name}[0] is zero; quotient is undefined")
    inverse = [1 / first]
    for index in range(1, order + 1):
        coupling = sum(
            comb(index, split) * rows[split] * inverse[index - split]
            for split in range(1, index + 1)
        )
        inverse.append(-coupling / first)
    return tuple(inverse)


def quotient_jets(
    numerator: Sequence[Any],
    denominator: Sequence[Any],
    *,
    order: int = 3,
) -> tuple[Any, ...]:
    """Return ordinary rows of ``numerator / denominator``."""

    order = _order(order, maximum=4)
    num = _jet(numerator, "numerator", order + 1)
    inv = _inverse_jets(denominator, order, "denominator")
    return ordinary_product_jets(num, inv, order=order)


def _named_rows(rows: Sequence[Any], count: int | None = None) -> dict[str, Any]:
    limit = len(rows) if count is None else min(count, len(rows))
    return {f"Z{index}": rows[index] for index in range(limit)}


def _coefficient_rows(
    w1: Sequence[Any],
    w2: Sequence[Any],
    C1: Sequence[Any],
    C2: Sequence[Any],
    order: int,
) -> tuple[Any, ...]:
    return _add(
        ordinary_product_jets(w1, C1, order=order),
        ordinary_product_jets(w2, C2, order=order),
        order,
    )


def _moment_rows(
    u: Sequence[Any],
    c: Sequence[Any],
    order: int,
) -> tuple[Any, ...]:
    return _negate(ordinary_product_jets(u, c, order=order), order)


def _delta_c_paths(
    w1_P: Sequence[Any],
    w2_P: Sequence[Any],
    w1_F: Sequence[Any],
    w2_F: Sequence[Any],
    C1_P: Sequence[Any],
    C2_P: Sequence[Any],
    C1_F: Sequence[Any],
    C2_F: Sequence[Any],
    c_P: Sequence[Any],
    c_F: Sequence[Any],
    order: int,
) -> dict[str, Any]:
    dC1 = _difference(C1_P, C1_F, order)
    dC2 = _difference(C2_P, C2_F, order)
    dW1 = _difference(w1_P, w1_F, order)
    dW2 = _difference(w2_P, w2_F, order)
    direct = _difference(c_P, c_F, order)
    wP_deltaC = _add(
        ordinary_product_jets(w1_P, dC1, order=order),
        ordinary_product_jets(w2_P, dC2, order=order),
        order,
    )
    deltaW_CF = _add(
        ordinary_product_jets(dW1, C1_F, order=order),
        ordinary_product_jets(dW2, C2_F, order=order),
        order,
    )
    conditional = _add(wP_deltaC, deltaW_CF, order)
    return {
        "direct": direct,
        "conditional": conditional,
        "conditional_terms": {
            "wP_deltaC": wP_deltaC,
            "deltaW_CF": deltaW_CF,
        },
        "delta_C": {"C1": dC1, "C2": dC2},
        "delta_W": {"w1": dW1, "w2": dW2},
        "formula": "delta_c=wP*deltaC+deltaW*CF",
        "difference_consistency_unproved": True,
        "ordinary_Z_rows": list(range(order + 1)),
    }


def _delta_m_paths(
    u_P: Sequence[Any],
    u_F: Sequence[Any],
    c_F: Sequence[Any],
    delta_c: Mapping[str, Any],
    m_P: Sequence[Any],
    m_F: Sequence[Any],
    order: int,
) -> dict[str, Any]:
    delta_u = _difference(u_P, u_F, order)
    direct = _difference(m_P, m_F, order)
    uP_delta_c = ordinary_product_jets(u_P, delta_c["conditional"], order=order)
    delta_u_cF = ordinary_product_jets(delta_u, c_F, order=order)
    conditional = _negate(_add(uP_delta_c, delta_u_cF, order), order)
    return {
        "direct": direct,
        "conditional": conditional,
        "conditional_terms": {
            "uP_delta_c": uP_delta_c,
            "delta_u_cF": delta_u_cF,
        },
        "delta_u": delta_u,
        "formula": "delta_m=-(uP*delta_c+delta_u*cF)",
        "difference_consistency_unproved": True,
        "ordinary_Z_rows": list(range(order + 1)),
    }


def _q_side_jets(
    *,
    z: Any,
    delta: Sequence[Any],
    m: Sequence[Any],
    order: int,
    side: str,
) -> dict[str, Any]:
    """Recover one side's Q rows from the non-normalized m rows."""

    order = _order(order, maximum=3, name="q_order")
    d = _jet(delta, f"delta_{side}", order + 1)
    # D_j contains m_(j+1), hence Q3 needs m4 explicitly.
    mm = _jet(m, f"m_{side}", order + 2)
    zjet = _z_jet(z, order)
    z2 = ordinary_product_jets(zjet, zjet, order=order)
    one = _constant_jet(1, order)
    one_minus_delta = _subtract(one, d, order)
    one_minus_z2 = _subtract(one, z2, order)
    zm = ordinary_product_jets(zjet, mm, order=order)
    m_z = tuple(mm[index + 1] for index in range(order + 1))
    D = _add(
        ordinary_product_jets(one_minus_delta, zm, order=order),
        ordinary_product_jets(one_minus_z2, m_z, order=order),
        order,
    )
    L = _subtract(one, ordinary_product_jets(d, z2, order=order), order)
    Q = _negate(quotient_jets(D, L, order=order), order)
    return {
        "delta": tuple(d[: order + 1]),
        "m": tuple(mm[: order + 2]),
        "D": D,
        "L": L,
        "Q": Q,
        "D_rows": _named_rows(D),
        "L_rows": _named_rows(L),
        "Q_rows": _named_rows(Q),
        "ordinary_Z_orders": list(range(order + 1)),
        "m_rows_required_through": order + 1,
        "source_values_selected": False,
    }


def correlated_delta_q_jets(
    *,
    z: Any,
    delta_P: Sequence[Any],
    delta_F: Sequence[Any],
    m_P: Sequence[Any],
    m_F: Sequence[Any],
    order: int = 2,
) -> dict[str, Any]:
    """Return direct and conditional P-minus-F Q jets.

    The conditional path uses complete ordinary jets in

    ``DeltaQ = -DeltaD/L_P + D_F*DeltaL/(L_P*L_F)``.

    It is valid with ``delta_P != delta_F``.  For order three the caller must
    provide ``m4`` on both sides; missing ``m4`` raises
    :class:`MissingJetRowsError` instead of dropping that derivative.
    """

    order = _order(order, maximum=3, name="q_order")
    P = _q_side_jets(z=z, delta=delta_P, m=m_P, order=order, side="P")
    F = _q_side_jets(z=z, delta=delta_F, m=m_F, order=order, side="F")
    ddelta = _difference(P["delta"], F["delta"], order)
    dm = _difference(P["m"], F["m"], order + 1)
    zjet = _z_jet(z, order)
    z2 = ordinary_product_jets(zjet, zjet, order=order)
    one = _constant_jet(1, order)
    one_minus_delta_P = _subtract(one, P["delta"], order)
    delta_D = _add(
        _subtract(
            ordinary_product_jets(
                ordinary_product_jets(one_minus_delta_P, zjet, order=order),
                dm,
                order=order,
            ),
            ordinary_product_jets(
                ordinary_product_jets(ddelta, zjet, order=order),
                F["m"],
                order=order,
            ),
            order,
        ),
        ordinary_product_jets(
            _subtract(one, z2, order),
            tuple(dm[index + 1] for index in range(order + 1)),
            order=order,
        ),
        order,
    )
    delta_L = _negate(
        ordinary_product_jets(ddelta, z2, order=order),
        order,
    )
    inverse_LP = _inverse_jets(P["L"], order, "L_P")
    inverse_LF = _inverse_jets(F["L"], order, "L_F")
    minus_delta_D_over_LP = _negate(
        ordinary_product_jets(delta_D, inverse_LP, order=order),
        order,
    )
    DF_delta_L_over_LP_LF = ordinary_product_jets(
        ordinary_product_jets(F["D"], delta_L, order=order),
        ordinary_product_jets(inverse_LP, inverse_LF, order=order),
        order=order,
    )
    conditional = _add(
        minus_delta_D_over_LP,
        DF_delta_L_over_LP_LF,
        order,
    )
    direct = _difference(P["Q"], F["Q"], order)
    return {
        "P": P,
        "F": F,
        "direct": direct,
        "conditional": conditional,
        "conditional_terms": {
            "minus_delta_D_over_LP": minus_delta_D_over_LP,
            "DF_delta_L_over_LP_LF": DF_delta_L_over_LP_LF,
        },
        "delta_D": delta_D,
        "delta_L": delta_L,
        "delta_m": dm,
        "delta_delta": ddelta,
        "formula": "delta_Q=-delta_D/L_P+D_F*delta_L/(L_P*L_F)",
        "difference_consistency_unproved": True,
        "conditional_gates": [
            "P_F_source_function_equality_unproved",
            "delta_P_delta_F_common_unproved",
            "ordinary_Z_transport_unproved",
            "L_P_nonzero_unproved",
            "L_F_nonzero_unproved",
        ],
        "ordinary_Z_orders": list(range(order + 1)),
        "source_values_selected": False,
    }


def _physical_input(
    value: Sequence[Any],
    name: str,
    order: int,
) -> tuple[Any, ...]:
    return _jet(value, name, order + 1)[: order + 1]


def _resolve_physical_scales(
    *,
    T_P: Sequence[Any] | None,
    T_F: Sequence[Any] | None,
    T: Sequence[Any] | None,
    order: int,
) -> tuple[tuple[Any, ...] | None, tuple[Any, ...] | None, list[str]]:
    """Resolve explicit side scales, without silently making them common."""

    if T is not None:
        if T_P is not None or T_F is not None:
            raise PFNoncommonJetError(
                "supply explicit T_P/T_F or common T shortcut, not both"
            )
        common = _physical_input(T, "T", order)
        return (
            common,
            common,
            ["T_P_T_F_common_conditional_shortcut_unproved"],
        )
    if T_P is None and T_F is None:
        return None, None, ["T_P_T_F_not_supplied"]
    if T_P is None or T_F is None:
        raise PFNoncommonJetError(
            "T_P and T_F are both required; a missing side scale is not zero"
        )
    return (
        _physical_input(T_P, "T_P", order),
        _physical_input(T_F, "T_F", order),
        [],
    )


def physical_difference_jets(
    *,
    T_P: Sequence[Any] | None = None,
    T_F: Sequence[Any] | None = None,
    T: Sequence[Any] | None = None,
    R_P: Sequence[Any],
    R_F: Sequence[Any],
    S_F: Sequence[Any],
    delta_S: Sequence[Any],
    Qy_F: Sequence[Any],
    delta_Qy: Sequence[Any],
    order: int = 1,
    endpoint_direction: int = -1,
) -> dict[str, Any]:
    """Return Z0--Z1 physical density and directed endpoint products.

    ``S_F`` and ``delta_S`` are caller-supplied physical scale jets.  Explicit
    ``T_P`` and ``T_F`` rows are retained and validated even on this direct
    path so both nonzero ``T_P,Z`` and ``T_F,Z`` remain visible in the result;
    :func:`physical_difference_from_q_jets` constructs the scale jets from
    those side scales and the Q rows when that source-shaped route is desired.
    """

    order = _order(order, maximum=1, name="physical_order")
    if type(endpoint_direction) is not int or endpoint_direction not in (-1, 1):
        raise PFNoncommonJetError("endpoint_direction must be +1 or -1")
    T_P_rows, T_F_rows, scale_gates = _resolve_physical_scales(
        T_P=T_P, T_F=T_F, T=T, order=order
    )
    rp = _physical_input(R_P, "R_P", order)
    rf = _physical_input(R_F, "R_F", order)
    sf = _physical_input(S_F, "S_F", order)
    ds = _physical_input(delta_S, "delta_S", order)
    qyf = _physical_input(Qy_F, "Qy_F", order)
    dqy = _physical_input(delta_Qy, "delta_Qy", order)
    sp = _add(sf, ds, order)
    rp2 = ordinary_product_jets(rp, rp, order=order)
    rf2 = ordinary_product_jets(rf, rf, order=order)
    direct_density = _subtract(
        ordinary_product_jets(rp2, sp, order=order),
        ordinary_product_jets(rf2, sf, order=order),
        order,
    )
    radius_difference = _subtract(rp, rf, order)
    radius_sum = _add(rp, rf, order)
    radius_product = ordinary_product_jets(radius_difference, radius_sum, order=order)
    radius_delta_term = ordinary_product_jets(rp2, ds, order=order)
    radius_SF_term = ordinary_product_jets(radius_product, sf, order=order)
    conditional_density = _add(radius_delta_term, radius_SF_term, order)
    qyp = _add(qyf, dqy, order)
    direct_endpoint = _subtract(
        ordinary_product_jets(rp, qyp, order=order),
        ordinary_product_jets(rf, qyf, order=order),
        order,
    )
    endpoint_radius_term = ordinary_product_jets(rp, dqy, order=order)
    endpoint_radius_difference_term = ordinary_product_jets(
        radius_difference,
        qyf,
        order=order,
    )
    endpoint_increment = _add(
        endpoint_radius_term,
        endpoint_radius_difference_term,
        order,
    )
    directed_endpoint = tuple(endpoint_direction * value for value in endpoint_increment)
    return {
        "T": T_P_rows if T_P_rows == T_F_rows else None,
        "T_P": T_P_rows,
        "T_F": T_F_rows,
        "R_P": rp,
        "R_F": rf,
        "S_F": sf,
        "delta_S": ds,
        "S_P": sp,
        "density_difference": {
            "direct": direct_density,
            "conditional": conditional_density,
            "conditional_terms": {
                "R_P_squared_delta_S": radius_delta_term,
                "radius_difference_times_radius_sum_times_S_F": radius_SF_term,
            },
            "formula": "R_P^2*delta_S+(R_P-R_F)*(R_P+R_F)*S_F",
            "difference_consistency_unproved": True,
        },
        "Qy_F": qyf,
        "delta_Qy": dqy,
        "Qy_P": qyp,
        "endpoint": {
            "direct": direct_endpoint,
            "increment": endpoint_increment,
            "directed": directed_endpoint,
            "conditional_terms": {
                "R_P_delta_Qy": endpoint_radius_term,
                "radius_difference_QyF": endpoint_radius_difference_term,
            },
            "formula": "R_P*delta_Qy+(R_P-R_F)*QyF",
            "direction": endpoint_direction,
            "direction_text": "minus (upper minus lower)"
            if endpoint_direction == -1
            else "plus (upper minus lower)",
            "difference_consistency_unproved": True,
        },
        "gates": [
            "F_radius_owner_unproved",
            "radius_difference_value_unproved",
            "physical_scale_source_identity_unproved",
            "T_Z_zero_not_assumed",
            *scale_gates,
        ],
        "ordinary_Z_orders": list(range(order + 1)),
        "source_values_selected": False,
    }


def _side_physical_scale_from_q(
    *,
    z: Any,
    T: Sequence[Any],
    delta: Sequence[Any],
    Q: Sequence[Any],
    order: int,
    side: str,
) -> dict[str, tuple[Any, ...]]:
    """Construct the source-shaped S and Qy jets for one P/F side."""

    if order != 1:
        raise PFNoncommonJetError("physical scale helper currently supports Z0--Z1")
    t = _physical_input(T, "T", order)
    d = _physical_input(delta, f"delta_{side}", order)
    q = _jet(Q, f"Q_{side}", order + 2)
    zjet = _z_jet(z, order)
    z2 = ordinary_product_jets(zjet, zjet, order=order)
    one = _constant_jet(1, order)
    L = _subtract(one, ordinary_product_jets(d, z2, order=order), order)
    q0 = (q[0], q[1])
    q1 = (q[1], q[2])
    one_minus_d = _subtract(one, d, order)
    linear_num = ordinary_product_jets(
        ordinary_product_jets(one_minus_d, zjet, order=order),
        q1,
        order=order,
    )
    linear = quotient_jets(
        linear_num,
        tuple(2 * value for value in L),
        order=order,
    )
    q0_square = ordinary_product_jets(q0, q0, order=order)
    quadratic = tuple(-value / 2 for value in q0_square)
    t2 = ordinary_product_jets(t, t, order=order)
    S = _add(
        ordinary_product_jets(t, linear, order=order),
        ordinary_product_jets(t2, quadratic, order=order),
        order,
    )
    Qy = _negate(ordinary_product_jets(t, q0, order=order), order)
    return {
        "slinear": linear,
        "squadratic": quadratic,
        "S": S,
        "Qy": Qy,
    }


def physical_difference_from_q_jets(
    *,
    z: Any,
    T_P: Sequence[Any] | None = None,
    T_F: Sequence[Any] | None = None,
    T: Sequence[Any] | None = None,
    delta_P: Sequence[Any],
    delta_F: Sequence[Any],
    Q_P: Sequence[Any],
    Q_F: Sequence[Any],
    R_P: Sequence[Any],
    R_F: Sequence[Any],
    endpoint_direction: int = -1,
) -> dict[str, Any]:
    """Construct physical Z0--Z1 products from two normalized Q sides.

    ``T_P`` and ``T_F`` are independent by default and both retain their Z1
    row.  The legacy ``T=...`` spelling is accepted only as an explicit common
    shortcut; the returned gate records that common-scale premise as unproved.
    """

    if type(endpoint_direction) is not int or endpoint_direction not in (-1, 1):
        raise PFNoncommonJetError("endpoint_direction must be +1 or -1")
    T_P_rows, T_F_rows, scale_gates = _resolve_physical_scales(
        T_P=T_P, T_F=T_F, T=T, order=1
    )
    if T_P_rows is None or T_F_rows is None:
        raise PFNoncommonJetError(
            "T_P and T_F are required for Q-derived physical products"
        )

    P = _side_physical_scale_from_q(
        z=z, T=T_P_rows, delta=delta_P, Q=Q_P, order=1, side="P"
    )
    F = _side_physical_scale_from_q(
        z=z, T=T_F_rows, delta=delta_F, Q=Q_F, order=1, side="F"
    )
    result = physical_difference_jets(
        T_P=T_P_rows,
        T_F=T_F_rows,
        R_P=R_P,
        R_F=R_F,
        S_F=F["S"],
        delta_S=_difference(P["S"], F["S"], 1),
        Qy_F=F["Qy"],
        delta_Qy=_difference(P["Qy"], F["Qy"], 1),
        order=1,
        endpoint_direction=endpoint_direction,
    )
    result["P_scale"] = P
    result["F_scale"] = F
    result["scale_difference"] = {
        "direct": _difference(P["S"], F["S"], 1),
        "ordinary_Z_orders": [0, 1],
    }
    result["gates"] = list(result["gates"]) + [
        "delta_P_delta_F_common_unproved",
        "T_Z_zero_not_assumed",
        *scale_gates,
    ]
    return result


def build_pf_noncommon_moment_jets(
    u_P: Sequence[Any],
    u_F: Sequence[Any],
    w1_P: Sequence[Any],
    w1_F: Sequence[Any],
    w2_P: Sequence[Any],
    w2_F: Sequence[Any],
    C1_P: Sequence[Any],
    C1_F: Sequence[Any],
    C2_P: Sequence[Any],
    C2_F: Sequence[Any],
    delta_P: Sequence[Any],
    delta_F: Sequence[Any],
    *,
    z: Any | None = None,
    Z: Any | None = None,
    q_order: int = 2,
    T_P: Sequence[Any] | None = None,
    T_F: Sequence[Any] | None = None,
    T: Sequence[Any] | None = None,
    R_P: Sequence[Any] | None = None,
    R_F: Sequence[Any] | None = None,
    endpoint_direction: int = -1,
) -> dict[str, Any]:
    """Build non-common P/F ``c``, ``m`` and conditional Q ordinary jets.

    The positional argument order mirrors the source roles: incoming rows,
    two independent first-row weight rows, then C1/C2 rows on each side, and
    finally the two ``delta(Z)`` jets.  All twelve required inputs need four
    rows.  ``q_order=3`` additionally requires five rows on every input that
    contributes to ``m4``; otherwise the Q3 request fails closed.  Physical
    products take independent ``T_P`` and ``T_F`` rows; ``T`` is retained only
    as an explicit common-scale shortcut with an unproved gate.
    """

    q_order = _order(q_order, maximum=3, name="q_order")
    if z is not None and Z is not None:
        raise PFNoncommonJetError("supply z or Z, not both")
    if z is None:
        z = Z
    if z is None:
        raise PFNoncommonJetError("a scalar z coordinate is required for Q recovery")
    input_values = {
        "u_P": u_P,
        "u_F": u_F,
        "w1_P": w1_P,
        "w1_F": w1_F,
        "w2_P": w2_P,
        "w2_F": w2_F,
        "C1_P": C1_P,
        "C1_F": C1_F,
        "C2_P": C2_P,
        "C2_F": C2_F,
        "delta_P": delta_P,
        "delta_F": delta_F,
    }
    rows = {name: _jet(value, name, 4) for name, value in input_values.items()}
    moment_input_rows = tuple(
        rows[name]
        for name in (
            "u_P", "u_F", "w1_P", "w1_F", "w2_P", "w2_F",
            "C1_P", "C1_F", "C2_P", "C2_F",
        )
    )
    requested_compute_order = (
        4 if q_order == 3 and all(len(value) >= 5 for value in moment_input_rows) else 3
    )
    c_P = _coefficient_rows(
        rows["w1_P"], rows["w2_P"], rows["C1_P"], rows["C2_P"], requested_compute_order
    )
    c_F = _coefficient_rows(
        rows["w1_F"], rows["w2_F"], rows["C1_F"], rows["C2_F"], requested_compute_order
    )
    m_P = _moment_rows(rows["u_P"], c_P, requested_compute_order)
    m_F = _moment_rows(rows["u_F"], c_F, requested_compute_order)
    c_delta = _delta_c_paths(
        rows["w1_P"], rows["w2_P"], rows["w1_F"], rows["w2_F"],
        rows["C1_P"], rows["C2_P"], rows["C1_F"], rows["C2_F"],
        c_P, c_F, requested_compute_order,
    )
    m_delta = _delta_m_paths(
        rows["u_P"], rows["u_F"], c_F, c_delta, m_P, m_F, requested_compute_order
    )
    if q_order == 3 and requested_compute_order < 4:
        raise MissingJetRowsError(
            "Q Z3 requires m4; provide five ordinary rows for every P/F "
            "incoming, weight and C jet"
        )
    q_delta = correlated_delta_q_jets(
        z=z,
        delta_P=rows["delta_P"],
        delta_F=rows["delta_F"],
        m_P=m_P,
        m_F=m_F,
        order=q_order,
    )
    gates = [
        "difference_consistency_unproved",
        "P_F_source_function_equality_unproved",
        "P_F_incoming_owner_unproved",
        "P_F_first_row_weight_owner_unproved",
        "delta_P_delta_F_common_unproved",
        "ordinary_Z_transport_unproved",
        "L_P_nonzero_unproved",
        "L_F_nonzero_unproved",
        "source_values_are_caller_owned_not_source_certificate",
    ]
    physical = None
    if any(value is not None for value in (T_P, T_F, T, R_P, R_F)):
        if (T_P is None and T_F is None and T is None) or R_P is None or R_F is None:
            raise PFNoncommonJetError(
                "T_P and T_F (or explicit common T shortcut), R_P and R_F "
                "are all required for physical products"
            )
        physical = physical_difference_from_q_jets(
            z=z,
            T_P=T_P,
            T_F=T_F,
            T=T,
            delta_P=rows["delta_P"],
            delta_F=rows["delta_F"],
            Q_P=q_delta["P"]["Q"],
            Q_F=q_delta["F"]["Q"],
            R_P=R_P,
            R_F=R_F,
            endpoint_direction=endpoint_direction,
        )
        gates.extend(physical["gates"])
    return {
        "schema": SCHEMA,
        "status": "conditional_caller_owned_ordinary_jet_algebra",
        "source_values_selected": False,
        "ordinary_Z_orders": list(ORDERS),
        "q_ordinary_Z_orders": list(range(q_order + 1)),
        "inputs": {
            "P": {
                "u": rows["u_P"],
                "w1": rows["w1_P"],
                "w2": rows["w2_P"],
                "C1": rows["C1_P"],
                "C2": rows["C2_P"],
                "delta": rows["delta_P"],
            },
            "F": {
                "u": rows["u_F"],
                "w1": rows["w1_F"],
                "w2": rows["w2_F"],
                "C1": rows["C1_F"],
                "C2": rows["C2_F"],
                "delta": rows["delta_F"],
            },
            "z": z,
        },
        "P": {
            "c": c_P,
            "c_rows": _named_rows(c_P, 4),
            "m": m_P,
            "m_rows": _named_rows(m_P, 4),
            "m4": m_P[4] if len(m_P) > 4 else None,
        },
        "F": {
            "c": c_F,
            "c_rows": _named_rows(c_F, 4),
            "m": m_F,
            "m_rows": _named_rows(m_F, 4),
            "m4": m_F[4] if len(m_F) > 4 else None,
        },
        "coefficient": {
            "P": {"c": c_P, "c_rows": _named_rows(c_P, 4)},
            "F": {"c": c_F, "c_rows": _named_rows(c_F, 4)},
            "difference": c_delta,
        },
        "moment": {
            "P": {"m": m_P, "m_rows": _named_rows(m_P, 4)},
            "F": {"m": m_F, "m_rows": _named_rows(m_F, 4)},
            "difference": m_delta,
        },
        "delta_c": c_delta,
        "delta_m": m_delta,
        "Q": q_delta,
        "physical": physical,
        "gates": sorted(set(gates)),
        "limitations": {
            "opaque_reference_mappings_rejected": True,
            "bare_integer_source_id_disambiguation": False,
            "P_F_function_equality_proved": False,
            "difference_consistency_proved": False,
            "F_radius_owner_proved": False,
            "physical_T_Z_zero_assumed": False,
            "T_P_equals_T_F_assumed": False,
            "delta_P_equals_delta_F_assumed": False,
            "L_P_nonzero_proved": False,
            "L_F_nonzero_proved": False,
            "numerical_or_error_certificate": False,
        },
    }


def build_noncommon_moment_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Convenience wrapper accepting explicit arguments or one mapping.

    Mapping aliases are consumed before dispatch.  Reference-like mapping
    leaves are rejected, while bare integer provenance remains unknowable.
    """

    if len(args) == 1 and isinstance(args[0], Mapping):
        payload = dict(args[0])
        aliases = {
            "uP": "u_P", "uF": "u_F", "w1P": "w1_P", "w1F": "w1_F",
            "w2P": "w2_P", "w2F": "w2_F", "C1P": "C1_P", "C1F": "C1_F",
            "C2P": "C2_P", "C2F": "C2_F", "deltaP": "delta_P", "deltaF": "delta_F",
        }
        for source, target in aliases.items():
            if target not in payload and source in payload:
                payload[target] = payload[source]
        for source in aliases:
            payload.pop(source, None)
        if "weights_P" in payload:
            payload.setdefault("w1_P", payload["weights_P"][0])
            payload.setdefault("w2_P", payload["weights_P"][1])
        if "weights_F" in payload:
            payload.setdefault("w1_F", payload["weights_F"][0])
            payload.setdefault("w2_F", payload["weights_F"][1])
        if "C_P" in payload:
            payload.setdefault("C1_P", payload["C_P"][0])
            payload.setdefault("C2_P", payload["C_P"][1])
        if "C_F" in payload:
            payload.setdefault("C1_F", payload["C_F"][0])
            payload.setdefault("C2_F", payload["C_F"][1])
        for key in ("weights_P", "weights_F", "C_P", "C_F"):
            payload.pop(key, None)
        return build_pf_noncommon_moment_jets(**payload)
    return build_pf_noncommon_moment_jets(*args, **kwargs)


compute_pf_noncommon_moment_jets = build_noncommon_moment_jets
pf_noncommon_moment_jets = build_noncommon_moment_jets


def verify_synthetic_fixture() -> dict[str, Any]:
    """Check the algebra with independent exact SymPy polynomial derivatives.

    The fixture intentionally uses distinct P/F polynomials and nonzero
    derivatives for every input family, including independent ``delta`` /
    ``T_P`` / ``T_F`` / radii.  It
    checks direct-vs-correlated differences, both endpoint signs, and the
    fail-closed Q3/m4 and missing-physical-row controls.  It returns exact
    ``sympy`` equality results and makes no numerical certification claim.
    """

    import sympy as sp

    z = sp.Symbol("z")
    z0 = sp.Rational(2, 11)

    def rows(expr: Any, count: int = 5) -> tuple[Any, ...]:
        return tuple(sp.diff(expr, z, index).subs(z, z0) for index in range(count))

    # Distinct degree-five polynomials supply nonzero P/F rows and m4.
    uP = 2 + 3 * z + 5 * z**2 + 7 * z**3 + 11 * z**4 + 13 * z**5
    uF = 3 + 2 * z + 7 * z**2 + 5 * z**3 + 13 * z**4 + 17 * z**5
    w1P = 1 + 2 * z + 3 * z**2 + 5 * z**3 + 7 * z**4 + 11 * z**5
    w1F = 2 + 3 * z + 5 * z**2 + 7 * z**3 + 11 * z**4 + 13 * z**5
    w2P = 3 + 5 * z + 7 * z**2 + 11 * z**3 + 13 * z**4 + 17 * z**5
    w2F = 5 + 7 * z + 11 * z**2 + 13 * z**3 + 17 * z**4 + 19 * z**5
    C1P = 7 + 11 * z + 13 * z**2 + 17 * z**3 + 19 * z**4 + 23 * z**5
    C1F = 11 + 13 * z + 17 * z**2 + 19 * z**3 + 23 * z**4 + 29 * z**5
    C2P = 13 + 17 * z + 19 * z**2 + 23 * z**3 + 29 * z**4 + 31 * z**5
    C2F = 17 + 19 * z + 23 * z**2 + 29 * z**3 + 31 * z**4 + 37 * z**5
    deltaP = sp.Rational(1, 7) + 2 * z + 3 * z**2 + 5 * z**3 + 7 * z**4
    deltaF = sp.Rational(2, 9) + 3 * z + 5 * z**2 + 7 * z**3 + 11 * z**4
    TP = sp.Rational(5, 4) + 2 * z + 3 * z**2 + 7 * z**3
    TF = sp.Rational(7, 6) + 3 * z + 5 * z**2 + 11 * z**3
    RP = 2 + 3 * z + 5 * z**2
    RF = 3 + 7 * z + 11 * z**2

    data = dict(
        u_P=rows(uP), u_F=rows(uF), w1_P=rows(w1P), w1_F=rows(w1F),
        w2_P=rows(w2P), w2_F=rows(w2F), C1_P=rows(C1P), C1_F=rows(C1F),
        C2_P=rows(C2P), C2_F=rows(C2F), delta_P=rows(deltaP), delta_F=rows(deltaF),
        z=z0, q_order=3, T_P=rows(TP), T_F=rows(TF),
        R_P=rows(RP, 2), R_F=rows(RF, 2),
    )
    result = build_pf_noncommon_moment_jets(**data)

    cP_expr = w1P * C1P + w2P * C2P
    cF_expr = w1F * C1F + w2F * C2F
    mP_expr = -uP * cP_expr
    mF_expr = -uF * cF_expr

    def q_expr(m_expr: Any, delta_expr: Any) -> Any:
        D = (1 - delta_expr) * z * m_expr + (1 - z**2) * sp.diff(m_expr, z)
        return -D / (1 - delta_expr * z**2)

    qP_expr = q_expr(mP_expr, deltaP)
    qF_expr = q_expr(mF_expr, deltaF)
    SP_expr = TP * (1 - deltaP) * z * sp.diff(qP_expr, z) / (2 * (1 - deltaP * z**2)) - TP**2 * qP_expr**2 / 2
    SF_expr = TF * (1 - deltaF) * z * sp.diff(qF_expr, z) / (2 * (1 - deltaF * z**2)) - TF**2 * qF_expr**2 / 2
    QyP_expr = -TP * qP_expr
    QyF_expr = -TF * qF_expr
    density_expr = RP**2 * SP_expr - RF**2 * SF_expr
    endpoint_expr = RP * QyP_expr - RF * QyF_expr

    checks: dict[str, bool] = {}

    def exact_rows(actual: Sequence[Any], expected: Any, count: int) -> bool:
        return all(
            sp.simplify(actual[index] - sp.diff(expected, z, index).subs(z, z0)) == 0
            for index in range(count)
        )

    checks["P_c_rows"] = exact_rows(result["P"]["c"], cP_expr, 4)
    checks["F_c_rows"] = exact_rows(result["F"]["c"], cF_expr, 4)
    checks["P_m_rows"] = exact_rows(result["P"]["m"], mP_expr, 4)
    checks["F_m_rows"] = exact_rows(result["F"]["m"], mF_expr, 4)
    checks["direct_delta_c"] = all(
        sp.simplify(
            result["delta_c"]["direct"][index]
            - sp.diff(cP_expr - cF_expr, z, index).subs(z, z0)
        ) == 0
        for index in range(4)
    )
    checks["conditional_delta_c"] = all(
        sp.simplify(
            result["delta_c"]["conditional"][index]
            - sp.diff(cP_expr - cF_expr, z, index).subs(z, z0)
        ) == 0
        for index in range(4)
    )
    checks["direct_delta_m"] = all(
        sp.simplify(
            result["delta_m"]["direct"][index]
            - sp.diff(mP_expr - mF_expr, z, index).subs(z, z0)
        ) == 0
        for index in range(4)
    )
    checks["conditional_delta_m"] = all(
        sp.simplify(
            result["delta_m"]["conditional"][index]
            - sp.diff(mP_expr - mF_expr, z, index).subs(z, z0)
        ) == 0
        for index in range(4)
    )
    checks["P_Q_rows"] = exact_rows(result["Q"]["P"]["Q"], qP_expr, 4)
    checks["F_Q_rows"] = exact_rows(result["Q"]["F"]["Q"], qF_expr, 4)
    checks["direct_delta_Q"] = all(
        sp.simplify(
            result["Q"]["direct"][index]
            - sp.diff(qP_expr - qF_expr, z, index).subs(z, z0)
        ) == 0
        for index in range(4)
    )
    checks["conditional_delta_Q"] = all(
        sp.simplify(
            result["Q"]["conditional"][index]
            - sp.diff(qP_expr - qF_expr, z, index).subs(z, z0)
        ) == 0
        for index in range(4)
    )
    checks["physical_density_direct"] = exact_rows(
        result["physical"]["density_difference"]["direct"], density_expr, 2
    )
    checks["physical_density_conditional"] = exact_rows(
        result["physical"]["density_difference"]["conditional"], density_expr, 2
    )
    checks["endpoint_direct"] = exact_rows(
        result["physical"]["endpoint"]["direct"], endpoint_expr, 2
    )
    checks["endpoint_minus_direction"] = all(
        result["physical"]["endpoint"]["directed"][index]
        == -result["physical"]["endpoint"]["increment"][index]
        for index in range(2)
    )
    plus = physical_difference_from_q_jets(
        z=z0, T_P=rows(TP), T_F=rows(TF),
        delta_P=rows(deltaP), delta_F=rows(deltaF),
        Q_P=result["Q"]["P"]["Q"], Q_F=result["Q"]["F"]["Q"],
        R_P=rows(RP, 2), R_F=rows(RF, 2), endpoint_direction=1,
    )
    checks["endpoint_plus_direction"] = all(
        plus["endpoint"]["directed"][index]
        == plus["endpoint"]["increment"][index]
        for index in range(2)
    )
    checks["independent_nonzero_T_rows"] = (
        TP.subs(z, z0) != TF.subs(z, z0)
        and sp.diff(TP, z).subs(z, z0) != 0
        and sp.diff(TF, z).subs(z, z0) != 0
    )
    checks["independent_T_rows_retained"] = (
        result["physical"]["T_P"][0] != result["physical"]["T_F"][0]
        and result["physical"]["T_P"][1] != 0
        and result["physical"]["T_F"][1] != 0
    )

    alias_payload = {
        "uP": data["u_P"], "uF": data["u_F"],
        "weights_P": (data["w1_P"], data["w2_P"]),
        "weights_F": (data["w1_F"], data["w2_F"]),
        "C_P": (data["C1_P"], data["C2_P"]),
        "C_F": (data["C1_F"], data["C2_F"]),
        "deltaP": data["delta_P"], "deltaF": data["delta_F"],
        "z": z0, "q_order": 3, "T_P": data["T_P"], "T_F": data["T_F"],
        "R_P": data["R_P"], "R_F": data["R_F"],
    }
    alias_result = build_noncommon_moment_jets(alias_payload)
    checks["mapping_aliases_consumed"] = all(
        alias_result["P"]["c"][index] == result["P"]["c"][index]
        for index in range(4)
    )

    missing_q3_rejected = False
    try:
        correlated_delta_q_jets(
            z=z0, delta_P=rows(deltaP, 4), delta_F=rows(deltaF, 4),
            m_P=rows(mP_expr, 4), m_F=rows(mF_expr, 4), order=3,
        )
    except MissingJetRowsError:
        missing_q3_rejected = True
    checks["missing_m4_rejected"] = missing_q3_rejected

    missing_physical_rejected = False
    try:
        physical_difference_jets(
            T_P=rows(TP, 1), T_F=rows(TF, 2),
            R_P=rows(RP, 2), R_F=rows(RF, 2),
            S_F=rows(SF_expr, 2), delta_S=rows(SP_expr - SF_expr, 2),
            Qy_F=rows(QyF_expr, 2), delta_Qy=rows(QyP_expr - QyF_expr, 2),
        )
    except MissingJetRowsError:
        missing_physical_rejected = True
    checks["missing_physical_derivative_rejected"] = missing_physical_rejected

    zero_L_rejected = False
    try:
        correlated_delta_q_jets(
            z=1, delta_P=(1, 0, 0), delta_F=(1, 0, 0),
            m_P=(1, 1, 1, 1), m_F=(2, 1, 1, 1), order=2,
        )
    except PFNoncommonJetError:
        zero_L_rejected = True
    checks["zero_L_rejected"] = zero_L_rejected

    opaque_dict_rejected = False
    try:
        ordinary_product_jets(
            ({"source_node_id": 17}, 0, 0, 0), (1, 0, 0, 0), order=3
        )
    except OpaqueNodeValueError:
        opaque_dict_rejected = True
    checks["opaque_mapping_rejected"] = opaque_dict_rejected

    bool_endpoint_rejected = False
    try:
        physical_difference_jets(
            T_P=(1, 1), T_F=(2, 1), R_P=(2, 1), R_F=(3, 1),
            S_F=(1, 1), delta_S=(2, 1), Qy_F=(1, 1), delta_Qy=(2, 1),
            endpoint_direction=True,
        )
    except PFNoncommonJetError:
        bool_endpoint_rejected = True
    checks["bool_endpoint_direction_rejected"] = bool_endpoint_rejected

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "scope": "independent exact SymPy synthetic polynomial rows",
        "numeric_approximation_used": False,
        "source_values_or_node_ids_used": False,
        "difference_consistency_gate_retained": True,
    }


__all__ = [
    "SCHEMA",
    "ORDERS",
    "Q_ORDERS",
    "PHYSICAL_ORDERS",
    "PFNoncommonJetError",
    "MissingJetRowsError",
    "OpaqueNodeValueError",
    "ordinary_product_jets",
    "ordinary_product_jet",
    "product_jets",
    "quotient_jets",
    "correlated_delta_q_jets",
    "physical_difference_jets",
    "physical_difference_from_q_jets",
    "build_pf_noncommon_moment_jets",
    "build_noncommon_moment_jets",
    "compute_pf_noncommon_moment_jets",
    "pf_noncommon_moment_jets",
    "verify_synthetic_fixture",
]
