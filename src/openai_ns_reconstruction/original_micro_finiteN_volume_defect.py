"""Pure append-only graph construction for the finite-``N`` volume defect.

The caller supplies the already admitted Q4 physical chart fields and the
prescribed-loop chart.  This module only ports their graph handles into the
signed pair builder and expands the six signed volume terms.  It deliberately
does not load a source chart, evaluate an integral, or evaluate a numerical
micro quantity.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import comb
from typing import Any

from .original_signed_graph_pairs import GraphPairs, Pair


_CERTIFICATE = "same positive physical L=1-delta*Z^2"
_DEFECT_NAMES = (
    "q",
    "q_Z",
    "q_squared",
    "U_q_Z",
    "U_squared",
    "U_U_Z",
)


def _graph_id(graph: GraphPairs, value: Any, *, name: str) -> int:
    """Validate one graph handle without allowing bool as an integer."""

    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer graph handle")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return graph._id(value, name=name)


def _rows(value: Any, *, name: str) -> tuple[Any, ...]:
    """Extract an explicit ordinary-Z row sequence from a field payload."""

    if isinstance(value, Mapping):
        for key in ("ordinary_Z", "ordinary_Z_rows", "rows", "values"):
            if key in value:
                value = value[key]
                break
        else:
            raise TypeError(f"{name} must contain an ordinary-Z row sequence")
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence of graph handles")
    return tuple(value)


def _field_rows(container: Mapping[str, Any], key: str, *, name: str) -> tuple[Any, ...]:
    if key not in container:
        raise KeyError(f"{name} is missing {key!r}")
    return _rows(container[key], name=f"{name}[{key!r}]")


def _checked_rows(
    graph: GraphPairs,
    rows: tuple[Any, ...],
    *,
    name: str,
    minimum: int,
) -> tuple[int, ...]:
    if len(rows) < minimum:
        raise ValueError(f"{name} requires at least {minimum} ordinary-Z rows")
    return tuple(_graph_id(graph, value, name=f"{name}[{index}]") for index, value in enumerate(rows))


def _same_length(jets: Sequence[Sequence[int]], *, name: str) -> int:
    if not jets:
        raise ValueError(f"{name} needs at least one jet")
    length = len(jets[0])
    if any(len(jet) != length for jet in jets[1:]):
        raise ValueError(f"{name} does not silently truncate unequal ordinary-Z rows")
    return length


def _jet_add(graph: GraphPairs, *jets: Sequence[int]) -> tuple[int, ...]:
    length = _same_length(jets, name="ordinary-Z addition")
    return tuple(graph.add(*(jet[index] for jet in jets)) for index in range(length))


def _jet_neg(graph: GraphPairs, jet: Sequence[int]) -> tuple[int, ...]:
    return tuple(graph.neg(value) for value in jet)


def _jet_sub(graph: GraphPairs, left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    _same_length((left, right), name="ordinary-Z subtraction")
    return _jet_add(graph, left, _jet_neg(graph, right))


def _jet_scale(graph: GraphPairs, jet: Sequence[int], scalar_id: int) -> tuple[int, ...]:
    scalar = _graph_id(graph, scalar_id, name="ordinary-Z scale scalar")
    return tuple(graph.mul(value, scalar) for value in jet)


def _jet_mul(
    graph: GraphPairs,
    left: Sequence[int],
    right: Sequence[int],
    order: int,
) -> tuple[int, ...]:
    """Multiply ordinary derivative rows with the binomial product rule."""

    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("ordinary-Z product order must be a nonnegative integer")
    if len(left) <= order or len(right) <= order:
        raise ValueError(f"ordinary-Z product requires rows through order {order}")
    result: list[int] = []
    for k in range(order + 1):
        terms = [
            graph.mul(
                graph.constant(comb(k, j)),
                left[j],
                right[k - j],
            )
            for j in range(k + 1)
        ]
        result.append(graph.add(*terms))
    return tuple(result)


def _jet_div(
    graph: GraphPairs,
    numerator: Sequence[int],
    denominator: Sequence[int],
    order: int,
    *,
    certificate: Any,
) -> tuple[int, ...]:
    """Divide ordinary derivative rows by the binomial quotient recurrence."""

    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("ordinary-Z quotient order must be a nonnegative integer")
    if len(numerator) <= order or len(denominator) <= order:
        raise ValueError(f"ordinary-Z quotient requires rows through order {order}")
    result: list[int] = [graph.div(numerator[0], denominator[0], certificate)]
    for k in range(1, order + 1):
        correction = [
            graph.mul(
                graph.constant(comb(k, j)),
                denominator[j],
                result[k - j],
            )
            for j in range(1, k + 1)
        ]
        residual = graph.add(numerator[k], graph.neg(graph.add(*correction)))
        result.append(graph.div(residual, denominator[0], certificate))
    return tuple(result)


def _ordinary_pair(jet: Sequence[int]) -> Pair:
    if len(jet) < 2:
        raise ValueError("a C1 pair requires ordinary-Z rows 0 and 1")
    return jet[0], jet[1]


def _symbolic_checks() -> dict[str, Any]:
    """Prove the six polynomial differences and their ordinary-Z1 rows."""

    try:
        import sympy as sp
    except ImportError:
        return {
            "available": False,
            "expansion": {name: False for name in _DEFECT_NAMES},
            "ordinary_Z1": {name: False for name in _DEFECT_NAMES},
            "passed": False,
        }

    z = sp.symbols("Z")
    U0 = sp.Function("U0")(z)
    Q0 = sp.Function("Q0")(z)
    dU = sp.Function("dU")(z)
    dQ = sp.Function("dQ")(z)
    U0z = sp.diff(U0, z)
    Q0z = sp.diff(Q0, z)
    dUz = sp.diff(dU, z)
    dQz = sp.diff(dQ, z)

    def terms(U: Any, Q: Any, Uz: Any, Qz: Any) -> dict[str, Any]:
        return {
            "q": Q,
            "q_Z": Qz,
            "q_squared": Q**2,
            "U_q_Z": U * Qz,
            "U_squared": U**2,
            "U_U_Z": U * Uz,
        }

    base = terms(U0, Q0, U0z, Q0z)
    full = terms(U0 + dU, Q0 + dQ, U0z + dUz, Q0z + dQz)
    expected = {
        "q": dQ,
        "q_Z": dQz,
        "q_squared": 2 * Q0 * dQ + dQ**2,
        "U_q_Z": U0 * dQz + Q0z * dU + dU * dQz,
        "U_squared": 2 * U0 * dU + dU**2,
        "U_U_Z": U0 * dUz + U0z * dU + dU * dUz,
    }

    def exact_zero(expression: Any) -> bool:
        return bool(sp.cancel(sp.together(expression)) == 0)

    expansion = {name: exact_zero(full[name] - base[name] - expected[name]) for name in _DEFECT_NAMES}
    ordinary_z1 = {
        name: exact_zero(sp.diff(full[name] - base[name] - expected[name], z))
        for name in _DEFECT_NAMES
    }
    return {
        "available": True,
        "expansion": expansion,
        "ordinary_Z1": ordinary_z1,
        "passed": all(expansion.values()) and all(ordinary_z1.values()),
    }


def build_defect(
    g: GraphPairs,
    physical: Mapping[str, Any],
    chart: Mapping[str, Any],
    Pstar_id: int = 33572,
    delta_id: int = 25220,
    Z_id: int = 20,
    N_id: int = 18,
) -> dict[str, Any]:
    """Build the exact six C1 signed finite-``N`` volume defect integrands.

    ``physical`` is the true Q4 physical field dictionary and ``chart`` is
    the prescribed-loop chart.  Every field is represented by graph handles;
    no source object is loaded and no numerical value is evaluated.
    """

    if not isinstance(g, GraphPairs):
        raise TypeError("g must be a GraphPairs instance")
    if not isinstance(physical, Mapping):
        raise TypeError("physical must be a mapping of true Q4 fields")
    if not isinstance(chart, Mapping):
        raise TypeError("chart must be a prescribed-loop chart mapping")

    Pstar = _graph_id(g, Pstar_id, name="Pstar_id")
    delta = _graph_id(g, delta_id, name="delta_id")
    Z = _graph_id(g, Z_id, name="Z_id")
    common_N = _graph_id(g, N_id, name="N_id")

    physical_U = _checked_rows(
        g,
        _field_rows(physical, "U", name="physical"),
        name="physical U",
        minimum=4,
    )
    physical_Q = _checked_rows(
        g,
        _field_rows(physical, "Q", name="physical"),
        name="physical Q",
        minimum=4,
    )

    finite_N = chart.get("finite_N")
    history = chart.get("history")
    if not isinstance(finite_N, Mapping):
        raise KeyError("chart is missing finite_N mapping")
    if not isinstance(history, Mapping):
        raise KeyError("chart is missing history mapping")

    delta_V = _checked_rows(
        g,
        _field_rows(finite_N, "delta_V", name="chart['finite_N']"),
        name="chart finite_N delta_V",
        minimum=3,
    )
    history_m = _checked_rows(
        g,
        _field_rows(history, "m", name="chart['history']"),
        name="chart history m",
        minimum=4,
    )

    zero = g.constant(0)
    one = g.constant(1)
    two = g.constant(2)

    # ``du`` uses the prescribed finite-N delta_V port, while ``dm`` retains
    # one extra ordinary row because dm_Z through Z2 consumes m_Z through Z3.
    du = _jet_scale(g, delta_V[:3], Pstar)
    dm = _jet_scale(g, history_m[:4], Pstar)
    zjet = (Z, one, zero, zero)
    deltajet = (delta, zero, zero, zero)
    onejet = (one, zero, zero, zero)
    z2 = _jet_mul(g, zjet, zjet, 3)
    Ljet = _jet_sub(g, onejet, _jet_mul(g, deltajet, z2, 3))
    Djet = _jet_sub(g, onejet, z2)
    hjet = _jet_sub(g, onejet, deltajet)
    dm_Z = dm[1:4]

    numerator = _jet_sub(
        g,
        _jet_scale(g, _jet_mul(g, zjet, du, 2), two),
        _jet_add(
            g,
            _jet_mul(g, _jet_mul(g, hjet, zjet, 2), dm, 2),
            _jet_mul(g, Djet, dm_Z, 2),
        ),
    )
    dq = _jet_div(g, numerator, Ljet, 2, certificate=_CERTIFICATE)

    # Subtract only the rows for which the finite-N correction is supplied;
    # source rows Z3 and above remain visible so a missing QZZ correction
    # cannot be mistaken for a closed estimate.
    U0 = _jet_sub(g, physical_U[:3], du)
    Q0 = _jet_sub(g, physical_Q[:3], dq)

    U0_pair = _ordinary_pair(U0)
    Q0_pair = _ordinary_pair(Q0)
    U0_Z_pair = (U0[1], U0[2])
    Q0_Z_pair = (Q0[1], Q0[2])
    du_pair = _ordinary_pair(du)
    dq_pair = _ordinary_pair(dq)
    du_Z_pair = (du[1], du[2])
    dq_Z_pair = (dq[1], dq[2])

    defect_integrands: dict[str, Pair] = {
        "q": dq_pair,
        "q_Z": dq_Z_pair,
        "q_squared": g.padd(
            g.pscale(g.pmul(Q0_pair, dq_pair), two),
            g.pmul(dq_pair, dq_pair),
        ),
        "U_q_Z": g.padd(
            g.pmul(U0_pair, dq_Z_pair),
            g.pmul(Q0_Z_pair, du_pair),
            g.pmul(du_pair, dq_Z_pair),
        ),
        "U_squared": g.padd(
            g.pscale(g.pmul(U0_pair, du_pair), two),
            g.pmul(du_pair, du_pair),
        ),
        "U_U_Z": g.padd(
            g.pmul(U0_pair, du_Z_pair),
            g.pmul(U0_Z_pair, du_pair),
            g.pmul(du_pair, du_Z_pair),
        ),
    }

    symbolic_checks = _symbolic_checks()
    if not symbolic_checks.get("passed", False):
        raise ValueError("exact raw finite-N volume defect symbolic checks failed")
    source_split_refs = {
        "baseline": {"U": U0, "Q": Q0},
        "correction": {"U": du, "Q": dq},
        "full_physical": {"U": physical_U, "Q": physical_Q},
        "Pstar": Pstar,
        "delta": delta,
        "Z": Z,
        "common_N": common_N,
        "delta_V_input": delta_V,
        "history_m_input": history_m,
        "du": du,
        "dm": dm,
        "dm_Z": dm_Z,
        "dq": dq,
    }
    qzz_gap = {
        "dq_ordinary_Z0_Z1_Z2_refs": dq,
        "dq_orders": (0, 1, 2),
        "mean_m_input_refs": dm,
        "mean_m_input_orders": (0, 1, 2, 3),
        "mean_m_Z3_input_ref": dm[3],
        "physical_Q_rows_retained": physical_Q,
        "physical_U_rows_retained": physical_U,
        "QZZ_gap_exposed": True,
        "ordinary_Z2_correction_is_available": True,
        "ordinary_Z3_correction_is_not_claimed": True,
        "no_Z3_baseline_claim": True,
        "noZ3baselineclaim": True,
    }
    flags = {
        "coefficient_profiles_depend_on_exact_common_N_through_phase_history": True,
        "exact_common_N_phase_history_dependency": True,
        "uniform_in_N_bound": False,
        "no_uniform_in_N_bound": True,
        "cap_as_value": False,
        "no_cap_as_value": True,
        "whole_micro_numeric": False,
        "no_whole_micro_numeric": True,
        "numeric_evaluation": False,
        "integrals_built": False,
        "source_chart_serialized": False,
        "old_actual_full_functions_retained": True,
        "no_Z3_baseline_claim": True,
        "noZ3baselineclaim": True,
        "V_equals_Vbackground_plus_delta_V_ast_guard_required": True,
    }

    return {
        "passed": bool(symbolic_checks["passed"]),
        "defect_integrands": defect_integrands,
        "integrands": defect_integrands,
        **defect_integrands,
        "du": du,
        "dm": dm,
        "dm_Z": dm_Z,
        "dq": dq,
        "U0": U0,
        "Q0": Q0,
        "full_physical": {"U": physical_U, "Q": physical_Q},
        "source_split_refs": source_split_refs,
        "source_split": source_split_refs,
        "qzz_gap": qzz_gap,
        "ordinary_Z0_Z1_Z2_dq_refs": dq,
        "mean_Z3_input_ref": dm[3],
        "source_provenance": {
            "du_source": "chart['finite_N']['delta_V']",
            "dm_source": "chart['history']['m']",
            "physical_U_source": "physical['U']",
            "physical_Q_source": "physical['Q']",
            "V_source_guard": "V=Vbackground+delta_V (AST guard supplied by parent)",
            "source_chart_serialized": False,
        },
        "symbolic_checks": symbolic_checks,
        "exact_symbolic_expansion_checks": symbolic_checks["expansion"],
        "ordinary_Z1_checks": symbolic_checks["ordinary_Z1"],
        "ordinary_Z1": symbolic_checks["ordinary_Z1"],
        "flags": flags,
        "coefficient_profiles_depend_on_exact_common_N_through_phase_history": True,
        "uniform_in_N_bound": False,
        "cap_as_value": False,
        "whole_micro_numeric": False,
        "common_N": common_N,
    }


__all__ = ["build_defect"]
