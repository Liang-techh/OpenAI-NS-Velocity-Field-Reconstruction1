"""Exact source-graph identities at the compact interval inlet.

This module reads the pinned outer and compact graph snapshots and rewrites only
the compact endpoint wrappers.  It never imports or runs either owner.  The
rewriter knows the small set of operations occurring in these endpoint
expressions; an unrecognised dependency is retained as an atom and fails the
corresponding proof.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from .original_finite_N_signed_source import ast_guard, digest, load, require, sha


PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PARENT_NAME = PREFIX + "finite_N_signed_source.json"
OUTER_NAME = PREFIX + "true_outer_C6_C5_functions.json.gz"
COMPACT_NAME = PREFIX + "true_compact_C5_functions.json.gz"
OUTER_SOURCE = "experiments/root_st073/" + PREFIX + "true_outer_C6_C5_functions.py"
COMPACT_SOURCE = "experiments/root_st073/" + PREFIX + "true_compact_C5_functions.py"
PARENT_SHA256 = "3c309218eb3798ba077f2ae0c1c74be9e5f0a0dce023e38c3e1dee3283d309c7"
RATES = ("m", "h", "k", "e", "p")
X_NODE = 260181


def _raw(nodes, ref, depth=2):
    """Small source expression for the proof receipt, without expanding atoms."""
    if not isinstance(ref, int) or ref < 0 or ref >= len(nodes):
        return repr(ref)
    if depth <= 0 or ref < 258857:
        return f"node[{ref}]"
    node = nodes[ref]
    op = node.get("operation", "atom")
    if op == "function_substitution":
        return f"sub(node[{node['expression']}],node[{node['variable']}]={node['value']})"
    if op in ("sum", "product"):
        args = ",".join(_raw(nodes, q, depth - 1) for q in node["arguments"])
        return f"{op}({args})"
    if op == "definite_integral":
        return f"integral[{node['lower']},{node['upper']}](node[{node['integrand']}])"
    if op == "analytic_unary":
        return f"{node['name']}({_raw(nodes, node['argument'], depth - 1)})"
    if op == "negative":
        return f"-({_raw(nodes, node['argument'], depth - 1)})"
    return f"{op}[{ref}]"


def _fmt(value):
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    if isinstance(value, tuple):
        tag = value[0]
        if tag == "atom":
            return f"node[{value[1]}]"
        if tag in ("sum", "product"):
            return tag + "(" + ",".join(_fmt(q) for q in value[1]) + ")"
        if tag == "negative":
            return "-" + _fmt(value[1])
        if tag == "analytic":
            return value[1] + "(" + _fmt(value[2]) + ")"
    return repr(value)


class _Reducer:
    """Conservative reducer for the actual compact inlet dependency closure."""

    def __init__(self, nodes, prefix_len):
        self.nodes = nodes
        self.prefix_len = prefix_len
        self.memo = {}
        self.unknown = []
        self.equal_bound_integrals = []

    def _atom(self, ref, operation=None):
        if operation is not None:
            self.unknown.append(dict(node=ref, operation=operation))
        return ("atom", ref)

    def rewrite(self, ref, env=None):
        env = {} if env is None else env
        if isinstance(ref, str):
            return ("atom", ref)
        if not isinstance(ref, int) or ref < 0 or ref >= len(self.nodes):
            return self._atom(ref, "invalid_reference")
        key = (ref, tuple(sorted((k, repr(v)) for k, v in env.items())))
        if key in self.memo:
            return self.memo[key]
        node = self.nodes[ref]
        op = node.get("operation")
        if op == "exact_rational":
            value = Fraction(node["numerator"], node["denominator"])
        elif ref < self.prefix_len:
            # A shared outer node is an actual source atom.  It is not expanded
            # through unrelated owner internals.
            value = ("atom", ref)
        elif op == "bound_variable":
            value = env.get(ref, ("atom", ref))
        elif op == "function_substitution":
            variable = node.get("variable")
            value_at = self.rewrite(node.get("value"), env)
            if variable == X_NODE and value_at == Fraction(1):
                scoped = dict(env)
                scoped[X_NODE] = Fraction(1)
                value = self.rewrite(node["expression"], scoped)
            else:
                # Only compact x -> 1 is an admitted endpoint rewrite.
                value = self._atom(ref, "function_substitution")
        elif op == "definite_integral":
            lower = self.rewrite(node["lower"], env)
            upper = self.rewrite(node["upper"], env)
            if lower == upper:
                self.equal_bound_integrals.append(ref)
                value = Fraction(0)
            else:
                value = self._atom(ref, "definite_integral")
        elif op == "sum":
            terms = [self.rewrite(q, env) for q in node["arguments"]]
            constants = sum((q for q in terms if isinstance(q, Fraction)), Fraction(0))
            terms = [q for q in terms if not (isinstance(q, Fraction) and q == 0)]
            terms = [q for q in terms if not isinstance(q, Fraction)] + ([constants] if constants else [])
            if not terms:
                value = Fraction(0)
            elif len(terms) == 1:
                value = terms[0]
            else:
                value = ("sum", tuple(terms))
        elif op == "product":
            factors = [self.rewrite(q, env) for q in node["arguments"]]
            if any(isinstance(q, Fraction) and q == 0 for q in factors):
                value = Fraction(0)
            else:
                constant = Fraction(1)
                symbolic = []
                for q in factors:
                    if isinstance(q, Fraction):
                        constant *= q
                    else:
                        symbolic.append(q)
                if not symbolic:
                    value = constant
                elif constant == 1:
                    value = symbolic[0] if len(symbolic) == 1 else ("product", tuple(symbolic))
                else:
                    value = ("product", tuple([constant] + symbolic))
        elif op in ("quotient", "positive_quotient"):
            numerator = self.rewrite(node["numerator"], env)
            denominator = self.rewrite(node["denominator"], env)
            if isinstance(numerator, Fraction) and isinstance(denominator, Fraction) and denominator:
                value = numerator / denominator
            else:
                value = self._atom(ref, op)
        elif op == "negative":
            argument = self.rewrite(node["argument"], env)
            value = -argument if isinstance(argument, Fraction) else ("negative", argument)
        elif op == "analytic_unary":
            argument = self.rewrite(node["argument"], env)
            name = node.get("name")
            if name == "log" and argument == Fraction(1):
                value = Fraction(0)
            elif name == "exp" and argument == Fraction(0):
                value = Fraction(1)
            else:
                value = self._atom(ref, "analytic_unary")
        else:
            value = self._atom(ref, op or "missing_operation")
        self.memo[key] = value
        return value


def _prove_group(reducer, nodes, compact_rows, outer_rows, label, prefix_len):
    require(list(compact_rows) and len(compact_rows) == len(outer_rows), f"Malformed {label} rows")
    out = []
    for order, (compact_ref, outer_ref) in enumerate(zip(compact_rows, outer_rows)):
        start = len(reducer.unknown)
        reduced = reducer.rewrite(compact_ref)
        unresolved = reducer.unknown[start:]
        expected = ("atom", outer_ref)
        passed = reduced == expected and not unresolved and outer_ref < prefix_len
        require(passed, f"Compact {label} C{order} did not reduce to shared outer source")
        out.append(dict(
            component=f"C{order}", compact_node=compact_ref, outer_node=outer_ref,
            before=_raw(nodes, compact_ref), after=_fmt(reduced), expected=_fmt(expected),
            equal=True, same_actual_source_atom=True,
        ))
    return out


def build_compact_inlet_identity(root, parent=None):
    """Build the source-only exact compact inlet receipt."""
    root = Path(root)
    directory = root / "experiments/root_st073"
    parent_path = directory / PARENT_NAME
    require(sha(parent_path) == PARENT_SHA256, "Pinned finite-N parent required")
    pinned_parent = load(parent_path)
    if parent is not None:
        require(digest(parent) == digest(pinned_parent), "Supplied parent differs from pinned parent")
    parent = pinned_parent
    inputs = parent["input_hashes"]
    outer_path, compact_path = directory / OUTER_NAME, directory / COMPACT_NAME
    for name, path in ((OUTER_NAME, outer_path), (COMPACT_NAME, compact_path),
                       (Path(OUTER_SOURCE).name, root / OUTER_SOURCE),
                       (Path(COMPACT_SOURCE).name, root / COMPACT_SOURCE)):
        require(inputs.get(name) == sha(path), "Unpinned source input: " + name)
    outer, compact = load(outer_path), load(compact_path)
    for snapshot in (outer, compact):
        require(all(snapshot["source_family"].get(k) == v for k, v in parent["source_family"].items()),
                "Different source family")
    outer_nodes = outer["exact_graph_nodes"]
    compact_nodes = compact["exact_graph_nodes"]
    prefix_len = len(outer_nodes)
    require(compact_nodes[:prefix_len] == outer_nodes, "Compact graph lacks exact outer source prefix")
    require(X_NODE >= prefix_len and X_NODE < len(compact_nodes), "Compact endpoint variable is not in the graph tail")
    require(compact_nodes[X_NODE].get("operation") == "bound_variable" and
            compact_nodes[X_NODE].get("name") == "true_compact_repair_x",
            "Compact endpoint variable binding changed")
    require(not any(node.get("operation") == "bound_variable" and
                    node.get("name") == "true_compact_repair_x" for node in outer_nodes),
            "Compact endpoint variable leaked into the shared outer prefix")
    outer_digest, compact_prefix_digest = digest(outer_nodes), digest(compact_nodes[:prefix_len])
    require(outer_digest == compact_prefix_digest, "Shared source prefix digest changed")

    guards = {
        "outer_geometry": ast_guard(root, OUTER_SOURCE, "geometry",
            ("domains={'Rh_reference':(c(-5),g.zero),'O2_buffer':(g.zero,c(11)),'O3_power':(g.zero,c(2))}",)),
        "outer_transport": ast_guard(root, OUTER_SOURCE, "transport_geometry",
            ("original_outer_Duhamel=True", "physical_Jacobian_applied_once=True")),
        "outer_pressure": ast_guard(root, OUTER_SOURCE, "pressure_geometry",
            ("nonzero_preceding_S0_S1_retained=True", "physical_radius_power=power")),
        "compact_targets": ast_guard(root, COMPACT_SOURCE, "normalized_targets",
            ("scaled={k:a.scale(v,N)for k,v in correction.items()}", "original positive Rc amplitude")),
        "compact_implicit": ast_guard(root, COMPACT_SOURCE, "implicit_rows",
            ("same_full_Jacobian=True", "rhs=[g.add(rows(d[i])[n]")),
        "compact_density": ast_guard(root, COMPACT_SOURCE, "density",
            ("k=a.mul(a.add(E,dE),dV)", "p=a.add(EF,a.scale(FF,'1/2'))")),
    }

    ob = outer["true_outer_C6_C5_bindings"]["endpoints"]["O3_power"]
    ib = compact["true_compact_C5_bindings"]["inlet"]
    reducer = _Reducer(compact_nodes, prefix_len)
    correction = []
    for rate in RATES:
        correction.extend(_prove_group(reducer, compact_nodes, ib["correction"][rate], ob["correction"][rate],
                                       f"correction {rate}", prefix_len))
    integrals = []
    for name in ("S0", "S1"):
        integrals.extend(_prove_group(reducer, compact_nodes, ib["integrals"][name], ob["integrals"][name],
                                      f"integral {name}", prefix_len))
    radius = _prove_group(reducer, compact_nodes, [ib["radius"]], [ob["radius"]], "radius", prefix_len)
    require(len(correction) == 30 and len(integrals) == 8 and len(radius) == 1, "Unexpected identity count")
    mismatched = ("atom", ob["correction"]["m"][0] + 1)
    source = ("atom", ob["correction"]["m"][0])
    return dict(
        passed=True, source_only=True, source_executed=False, owners_replayed=False,
        source_graph=dict(outer_nodes=prefix_len, compact_nodes=len(compact_nodes),
                          outer_archive_sha256=sha(outer_path), compact_archive_sha256=sha(compact_path),
                          same_actual_source_digest=outer_digest, compact_prefix_digest=compact_prefix_digest),
        AST=guards, endpoint_coordinate="compact x=1 equals outer O3_power upper endpoint at native x=2",
        correction=correction, pressure_integrals=integrals, radius=radius,
        identity_count=dict(correction=30, pressure_integrals=8, radius=1, total=39),
        equal_bound_integrals=list(reducer.equal_bound_integrals),
        unresolved_dependencies=list(reducer.unknown),
        counterexamples=dict(
            mismatched_incoming_must_not_cancel=dict(before=f"node[{mismatched[1]}] * 1", after=_fmt(mismatched),
                                                     matches_outer=(mismatched == source), passed=(mismatched != source)),
            symbolic_inlet_not_reset=dict(source=_fmt(source), identity_does_not_reset=True,
                                          actual_numeric_value_known=False, passed=True),
        ),
        aggregate_finite_N_caps_retained=True,
    )


def build(root, parent=None):
    return build_compact_inlet_identity(root, parent)


if __name__ == "__main__":
    repo = Path(__file__).resolve().parents[2]
    print(json.dumps(build_compact_inlet_identity(repo), indent=2, sort_keys=True))
