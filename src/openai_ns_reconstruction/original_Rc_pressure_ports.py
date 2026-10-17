"""Source-only pressure-port congruence for the quiet O3 / compact inlet.

The caller supplies the already saved graph tables and the two binding
objects.  This module deliberately imports no owner module and never builds,
replays, or numerically evaluates a source object.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import sympy as sp


class UnknownGraphOperation(ValueError):
    """Raised when the requested source slice is not an exact algebra DAG."""


def _binding(value: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    """Accept either a binding object or the corresponding archive wrapper."""
    if name in value and isinstance(value[name], Mapping):
        return value[name]
    return value


def _node_tables(nodes: Any) -> tuple[Sequence[Mapping[str, Any]], Sequence[Mapping[str, Any]]]:
    """Get outer and compact node tables without importing an archive loader."""
    def table(value: Any) -> Sequence[Mapping[str, Any]]:
        if isinstance(value, Mapping) and "exact_graph_nodes" in value:
            return value["exact_graph_nodes"]
        return value

    if isinstance(nodes, Mapping):
        if "outer" in nodes and "compact" in nodes:
            return table(nodes["outer"]), table(nodes["compact"])
        if "outer_nodes" in nodes and "compact_nodes" in nodes:
            return table(nodes["outer_nodes"]), table(nodes["compact_nodes"])
        if "true_outer_C6_C5_nodes" in nodes and "true_compact_C5_nodes" in nodes:
            return table(nodes["true_outer_C6_C5_nodes"]), table(nodes["true_compact_C5_nodes"])
    if isinstance(nodes, (tuple, list)) and len(nodes) == 2:
        return table(nodes[0]), table(nodes[1])
    raise TypeError("nodes must contain separate outer and compact exact_graph_nodes tables")


def _receipt_ok(receipt: Any, label: str) -> bool:
    """Require the caller's inherited source receipt and reject a failed one."""
    if not isinstance(receipt,Mapping) or receipt.get('passed') is not True or receipt.get('source_only') is not True:
        raise ValueError(f"{label} passed source-only receipt is required")
    for key in ("passed", "ok", "source_only"):
        if key in receipt and receipt[key] is False:
            raise ValueError(f"{label} receipt has {key}=false")
    return True


class _Normalizer:
    """Small exact SymPy evaluator for the admitted source operations."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]], replacements: Mapping[int, sp.Expr]):
        self.nodes = nodes
        self.replacements = dict(replacements)
        self._memo: dict[tuple[int, tuple[tuple[int, sp.Expr], ...]], sp.Expr] = {}

    def __call__(self, node: int, variables: Mapping[int, sp.Expr] | None = None) -> sp.Expr:
        variables = {} if variables is None else dict(variables)
        key = (int(node), tuple(sorted(variables.items(), key=lambda item: item[0])))
        if key in self._memo:
            return self._memo[key]
        value = self._eval(int(node), variables)
        self._memo[key] = value
        return value

    def _eval(self, node: int, variables: Mapping[int, sp.Expr]) -> sp.Expr:
        if node in variables:
            return variables[node]
        if node in self.replacements:
            return self.replacements[node]
        if not isinstance(node, int) or node < 0 or node >= len(self.nodes):
            raise UnknownGraphOperation(f"bad graph node {node}")
        record = self.nodes[node]
        operation = record.get("operation")
        if operation == "exact_rational":
            return sp.Rational(record["numerator"], record["denominator"])
        if operation == "sum":
            return sp.Add(*(self(arg, variables) for arg in record["arguments"]))
        if operation == "product":
            return sp.Mul(*(self(arg, variables) for arg in record["arguments"]))
        if operation == "negative":
            return -self(record["argument"], variables)
        if operation == "positive_quotient":
            denominator = self(record["denominator"], variables)
            if denominator == 0:
                raise ZeroDivisionError(f"zero positive_quotient denominator at node {node}")
            return sp.cancel(self(record["numerator"], variables) / denominator)
        if operation == "function_substitution":
            variable = int(record["variable"])
            value = self(record["value"], variables)
            child_variables = dict(variables)
            child_variables[variable] = value
            return self(record["expression"], child_variables)
        if operation == "bound_variable":
            # Z is the only free angular variable admitted by this proof.
            if record.get("name") == "Z" and node == 20:
                return _Z
            raise UnknownGraphOperation(
                f"unbound variable node {node} ({record.get('name')!r})"
            )
        raise UnknownGraphOperation(f"unsupported operation {operation!r} at node {node}")


_Z, _DELTA, _PSTAR, _RC = sp.symbols("Z delta Pstar Rc")
_MEAN = sp.symbols("M0:4")
_ZERO = sp.Integer(0)


def _side_replacements(
    *,
    binding: Mapping[str, Any],
    chart: Mapping[str, Any] | None,
    radius_node: int,
) -> dict[int, sp.Expr]:
    """Bind quiet-profile fields while preserving the raw pressure DAG."""
    fields = (chart or binding)["modified_fields"]
    radial = (chart or binding)["radial"]
    replacements: dict[int, sp.Expr] = {20: _Z, 25220: _DELTA, 33572: _PSTAR}
    for index, node in enumerate(fields["m"][:4]):
        replacements[int(node)] = _MEAN[index]
    for node in fields["V"]:
        replacements[int(node)] = _ZERO
    for node in radial["V_y"]:
        replacements[int(node)] = _ZERO
    replacements[int(radius_node)] = _RC
    return replacements


def _equal(left: sp.Expr, right: sp.Expr) -> bool:
    return sp.simplify(left - right) == 0


def _common_physical_ports() -> tuple[list[sp.Expr], list[sp.Expr]]:
    """Differentiate the physical recovery, using only the required M0..M3."""
    mean=sp.Function('M')(_Z)
    q=-_PSTAR*((1-_DELTA)*_Z*mean+(1-_Z**2)*sp.diff(mean,_Z))/(1-_DELTA*_Z**2)
    replacements={sp.diff(mean,_Z,k):_MEAN[k] for k in range(4)}
    rows=[sp.cancel(sp.diff(q,_Z,k).xreplace(replacements)) for k in range(3)]
    return rows,[-v for v in rows[:2]]


def _wrapper_rows(nodes,wrapped,raw,variable,value):
    return all(nodes[i].get('operation')=='function_substitution' and
        nodes[i].get('expression')==raw[k] and nodes[i].get('variable')==variable and
        nodes[i].get('value')==value and nodes[i].get('Z_independent_substitution') is True
        for k,i in enumerate(wrapped))


def _assert_all(checks: Mapping[str, bool]) -> bool:
    return all(bool(value) for value in checks.values())


def build_pressure_ports(
    nodes: Any,
    ob: Mapping[str, Any],
    cb: Mapping[str, Any],
    mean_receipt: Any,
    profile: Any,
) -> dict[str, Any]:
    """Prove literal source congruence of the physical O3 pressure ports.

    ``nodes`` must provide ``outer`` and ``compact`` exact graph tables.  The
    first three ``Q`` and first two ``Q_y`` rows are checked on each side,
    together with C0/C1 ``R_Q_y`` wrappers.  The remaining quiet-profile and
    mean-input facts are inherited only through the supplied receipts.
    """
    _receipt_ok(mean_receipt, "mean")
    _receipt_ok(profile, "profile")
    if not (len(mean_receipt.get('background_mean_identities',[]))==4 and
        len(mean_receipt.get('full_mean_identities',[]))==4 and
        all(v.get('passed') is True for v in mean_receipt['background_mean_identities']+mean_receipt['full_mean_identities']) and
        profile.get('full_V_and_V_y_exact_zero') is True and profile.get('baseline_V_and_V_y_exact_zero') is True):
        raise ValueError('Actual common mean and quiet velocity inputs required before normalization')
    outer_nodes, compact_nodes = _node_tables(nodes)
    outer = _binding(ob, "true_outer_C6_C5_bindings")
    compact = _binding(cb, "true_compact_C5_bindings")
    chart = outer["charts"]["O3_power"]
    endpoint = outer["endpoints"]["O3_power"]
    inlet = compact["inlet"]

    outer_radius = int(chart["radius"])
    compact_radius = int(compact["background"]["radius"])
    out_norm = _Normalizer(outer_nodes, _side_replacements(
        binding=outer, chart=chart, radius_node=outer_radius
    ))
    compact_norm = _Normalizer(compact_nodes, _side_replacements(
        binding=compact, chart=None, radius_node=compact_radius
    ))

    # These are the ten literal source rows: 3 Q + 2 Q_y on each side.
    out_x = {int(chart["coordinate"]): sp.Rational(2)}
    compact_x = {260181: sp.Rational(1)}
    out_q = [out_norm(node, out_x) for node in chart["physical"]["Q"][:3]]
    out_qy = [out_norm(node, out_x) for node in chart["physical"]["Q_y"][:2]]
    compact_q = [compact_norm(node, compact_x) for node in compact["physical"]["Q"][:3]]
    compact_qy = [compact_norm(node, compact_x) for node in compact["physical"]["Q_y"][:2]]
    out_u = [out_norm(node, out_x) for node in chart["physical"]["U"][:3]]
    compact_u = [compact_norm(node, compact_x) for node in compact["physical"]["U"][:3]]

    # M0 ... M3 are ordinary derivative rows. The scalar displayed
    # source identity is Q = -Pstar*((1-delta) Z M + (1-Z^2) M_Z)/L;
    # _common_physical_ports differentiates this formula independently.
    common_q, common_qy = _common_physical_ports()
    checks: dict[str, bool] = {
        "actual_source_bindings": all(key in chart["physical"] for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")),
        "outer_Q_common": all(_equal(value, common_q[index]) for index, value in enumerate(out_q)),
        "compact_Q_common": all(_equal(value, common_q[index]) for index, value in enumerate(compact_q)),
        "outer_Qy_common": all(_equal(value, common_qy[index]) for index, value in enumerate(out_qy)),
        "compact_Qy_common": all(_equal(value, common_qy[index]) for index, value in enumerate(compact_qy)),
        "outer_compact_Q_congruence": all(_equal(a, b) for a, b in zip(out_q, compact_q)),
        "outer_compact_Qy_congruence": all(_equal(a, b) for a, b in zip(out_qy, compact_qy)),
        "outer_U_zero": all(value == 0 for value in out_u),
        "compact_U_zero": all(value == 0 for value in compact_u),
    }

    # The stored outer physical endpoint must be the actual endpoint wrapper,
    # while the compact inlet has no stored physical field and is checked by a
    # virtual x -> 1 substitution of the raw physical source.
    endpoint_q = [out_norm(node) for node in endpoint["physical"]["Q"][:3]]
    endpoint_qy = [out_norm(node) for node in endpoint["physical"]["Q_y"][:2]]
    checks["outer_stored_physical_endpoint_wrapper"] = all(
        _equal(a, b) for a, b in zip(endpoint_q, out_q)
    ) and all(_equal(a, b) for a, b in zip(endpoint_qy, out_qy))

    out_rqy_raw = [out_norm(node, out_x) for node in chart["boundary"]["R_Q_y"][:2]]
    out_rqy_endpoint = [out_norm(node) for node in endpoint["boundary"]["R_Q_y"][:2]]
    compact_rqy_raw = [compact_norm(node, compact_x) for node in compact["boundary"]["R_Q_y"][:2]]
    compact_rqy_inlet = [compact_norm(node) for node in inlet["boundary"]["R_Q_y"][:2]]
    # The stored physical Q_y already carries Pstar.  To expose the source
    # factorization requested by the port contract, write it as Pstar times
    # the unscaled canonical Q_y, then multiply by R.
    canonical_unscaled_qy = [sp.cancel(value / _PSTAR) for value in common_qy]
    expected_rqy = [_PSTAR * canonical_unscaled_qy[index] * _RC for index in range(2)]
    checks.update({
        "outer_raw_R_Qy_is_P_Qy_R": all(_equal(a, b) for a, b in zip(out_rqy_raw, expected_rqy)),
        "outer_endpoint_R_Qy_C0_C1_wrapper": all(_equal(a, b) for a, b in zip(out_rqy_endpoint, expected_rqy)),
        "compact_raw_R_Qy_is_P_Qy_R": all(_equal(a, b) for a, b in zip(compact_rqy_raw, expected_rqy)),
        "compact_inlet_R_Qy_C0_C1_wrapper": all(_equal(a, b) for a, b in zip(compact_rqy_inlet, expected_rqy)),
    })

    L = 1 - _DELTA * _Z**2
    # These are F and its actual ordinary derivative, after U=U_Z=0.
    # Applying the value formula separately to Q_Z would not differentiate F.
    def F_rows(q):
        return [sp.cancel(q[0]/(2*L)+q[0]**2/4),
            sp.cancel(q[1]/(2*L)+q[0]*_DELTA*_Z/L**2+q[0]*q[1]/2)]
    out_f,compact_f=F_rows(out_q),F_rows(compact_q)
    checks['outer_exact_Q_Qy_substitution_wrappers']=all(_wrapper_rows(outer_nodes,
        endpoint['physical'][key][:count],chart['physical'][key][:count],chart['coordinate'],chart['domain'][1])
        for key,count in (('Q',3),('Q_y',2)))
    checks['outer_exact_RQy_substitution_wrappers']=_wrapper_rows(outer_nodes,
        endpoint['boundary']['R_Q_y'][:2],chart['boundary']['R_Q_y'][:2],chart['coordinate'],chart['domain'][1])
    checks['compact_exact_RQy_substitution_wrappers']=_wrapper_rows(compact_nodes,
        inlet['boundary']['R_Q_y'][:2],compact['boundary']['R_Q_y'][:2],260181,1)
    checks["F_source_compositional_congruence"] = all(
        _equal(a, b) for a, b in zip(out_f, compact_f)
    )

    return {
        "passed": _assert_all(checks),
        "source_only": True,
        "owners_replayed": False,
        "checks": checks,
        "rows": {"Q": 3, "Q_y": 2, "R_Q_y": 2, "F":2},
        "source_executed":False,"ordinary_derivatives_only":True,
        "actual_ports":dict(outer_Q=chart['physical']['Q'][:3],outer_Qy=chart['physical']['Q_y'][:2],
            compact_Q=compact['physical']['Q'][:3],compact_Qy=compact['physical']['Q_y'][:2],
            outer_RQy=endpoint['boundary']['R_Q_y'][:2],compact_RQy=inlet['boundary']['R_Q_y'][:2]),
        "outer": {"Q": [str(v) for v in out_q], "Q_y": [str(v) for v in out_qy], "F": [str(v) for v in out_f]},
        "compact": {"Q": [str(v) for v in compact_q], "Q_y": [str(v) for v in compact_qy], "F": [str(v) for v in compact_f]},
        "common": {"Q": [str(v) for v in common_q], "Q_y": [str(v) for v in common_qy]},
        "source_formula": "Q = -Pstar*((1-delta)*Z*M + (1-Z**2)*M_Z)/(1-delta*Z**2); Q_y = -Q",
        "canonical_unscaled_Q_y": [str(v) for v in canonical_unscaled_qy],
        "inherited_receipts": {"mean": True, "profile": True},
    }


__all__ = ["UnknownGraphOperation", "build_pressure_ports"]
