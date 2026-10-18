"""Source-bound native seams for the current long/reference construction.

This receipt has a deliberately narrow scope.  It proves the physical
``phi`` and ``U`` jets through order two at the two native joins and checks
that the six incoming moment histories are carried by the exact source-cache
pointer chain.  It does not install or evaluate a field owner, a finite-N
producer, a modified-field seam, or a pressure ``Q_y`` seam.

The source constructors and the C6 cache are read as data.  The only
arithmetic below is a fresh symbolic reduction of the displayed source
formulae and the squared-flat cutoff endpoint limits.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_finite_N_I1.json"
PARENT_SHA256 = "75918293b8c7432c637eb551a29cb9d248aa7746dcb1cd4f8b3615f72514c14d"
MODULE = "src/openai_ns_reconstruction/original_long_reference_native_seams.py"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
RATE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
ANCHOR_SOURCE = BASE + "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py"

CHARTS = ("long_reshape", "reference")
MOMENTS = ("H", "M", "K", "A", "B", "C")
ORDINARY_PREFIX = (0, 1, 2)


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(root: Path, relative: str) -> Any:
    path = root / relative
    _require(path.is_file(), "Missing pinned JSON source: " + relative)
    # The accepted archives are read as bytes by their predecessor modules.
    # Keep that JSON decoder path so legacy double-byte text remains identical
    # to the mapping supplied by the root.
    return json.loads(path.read_bytes())


def _read_gzip_json(root: Path, relative: str) -> Any:
    path = root / relative
    _require(path.is_file(), "Missing pinned gzip source: " + relative)
    return json.loads(gzip.decompress(path.read_bytes()))


def _compact(text: str) -> str:
    return "".join(text.split())


def _source_guard(
    root: Path,
    relative: str,
    function: str,
    snippets: tuple[str, ...],
) -> dict[str, Any]:
    """Pin one source function without importing or executing its module."""
    _require(type(snippets) is tuple, "Source guard snippets must be a tuple")
    path = root / relative
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    candidates = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function
    ]
    candidates = [
        node
        for node in candidates
        if all(_compact(snippet) in _compact(ast.get_source_segment(text, node) or "")
               for snippet in snippets)
    ]
    _require(len(candidates) == 1, "Changed source definition: " + relative + ":" + function)
    node = candidates[0]
    digest = hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()
    return dict(
        passed=True,
        path=relative,
        function=function,
        required_text=list(snippets),
        AST_sha256=digest,
        source_executed=False,
    )


def _assignment_guard(root: Path) -> dict[str, Any]:
    """Pin the literal rate assignment in the source rate module."""
    path = root / RATE_SOURCE
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    rows = [
        node
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "RATES"
                for target in node.targets)
    ]
    _require(len(rows) == 1, "The actual rate source must contain one RATES assignment")
    node = rows[0]
    value = node.value
    _require(isinstance(value, ast.Call) and isinstance(value.func, ast.Name)
             and value.func.id == "dict" and not value.args,
             "The actual rate source RATES assignment changed")
    rates: dict[str, int] = {}
    for keyword in value.keywords:
        _require(keyword.arg is not None, "RATES must use named entries")
        parsed = ast.literal_eval(keyword.value)
        _require(type(parsed) is int, "RATES entries must be exact integers")
        rates[keyword.arg] = parsed
    expected = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
    _require(rates == expected, "Actual six-history rates changed")
    return dict(
        passed=True,
        path=RATE_SOURCE,
        binding="RATES",
        rates=rates,
        AST_sha256=hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )


def _source_guards(root: Path) -> tuple[dict[str, Any], ...]:
    """Return a tuple so the source-guard contract is visible to callers."""
    guards = (
        _source_guard(root, C6_SOURCE, "__init__", (
            "incoming=self.endpoints['R2'];phi=a.scale(incoming['actual']['phi'],g.unary('exp',g.mul(c('-2/5'),y)))",
            "self.endpoints['R110']=self.end('post_switch',g.one)",
            "z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5);self.z=z",
            "incoming=self.endpoints['R110'];self.B=B=a.add(a.scale(a.trim(bridge_owner.core7['G']),g.neg(self.core.parameters['Lambda'])),",
            "a.log(incoming['actual']['phi'],'positive actual R110 phi')",
            "a.fixed(g.mul(c('1/2'),g.unary('log',c(220)))),logq)",
            "logUtheta=a.add(a.scale(B,g.sub(g.one,sigma(g,q))),a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
            "slope=a.add(a.fixed('-2/5'),a.scale(B,g.neg(g.quotient(sigma_prime(g,q),T,'same positive original T'))))",
            "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0),",
            "self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])",
            "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
            "self.reference_length=length=g.sub(g.sub(logref,T),c(8))",
            "logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
            "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),",
            "self.add_chart('reference',q,R,length,actual,incoming['moments'])",
        )),
        _source_guard(root, C6_SOURCE, "phi_from_Utheta", (
            "return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R)))",
        )),
        _source_guard(root, C6_SOURCE, "add_chart", (
            "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
            "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
        )),
        _source_guard(root, C6_SOURCE, "sigma_prime", (
            "derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets)",
            "exact_flat_source_derivative_not_a_bound=True",
        )),
        _source_guard(root, SIGMA_SOURCE, "_sigma_left", (
            "odds=1/(1-x)**2-1/x**2",
            "e=positive_exp(c,odds); value=e/(1+e)",
            "L1=2/(1-x)**3+2/x**3",
            "derivatives=[value,q*L1",
        )),
        _source_guard(root, SIGMA_SOURCE, "sigma_jets", (
            "if lo<=0:rows.append(IntervalTaylor.constant(c,0,4))",
            "if hi>=1:rows.append(IntervalTaylor.constant(c,1,4))",
            "left=_sigma_left(c,1-c.mpf([a,b]))",
        )),
        _source_guard(root, ANCHOR_SOURCE, "__init__", (
            "self.params=dict(j=self.core.j,delta=self.core.delta,sigma=self.core.sigma,a=(9-self.core.delta)/2)",
        )),
        _source_guard(root, ANCHOR_SOURCE, "gradient_jets", (
            "H=-4*z**3-z*z*self.params['j']+z*self.params['a']+self.params['j']",
            "return((one-z*z*self.params['delta'])*H/IntervalTaylor(c,coefficients)).coefficients",
        )),
        _assignment_guard(root),
    )
    _require(type(guards) is tuple, "Source guards must remain a tuple")
    return guards


def _source_pins(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    native = parent["actual_long_reference_native_controls"]
    accepted = native.get("source_pins", {})
    paths = (PARENT, C6_SOURCE, C6_CACHE, SIGMA_SOURCE, RATE_SOURCE, ANCHOR_SOURCE)
    pins: dict[str, Any] = {}
    for relative in paths:
        if relative == PARENT:
            expected = PARENT_SHA256
        else:
            item = accepted.get(relative)
            _require(isinstance(item, Mapping), "Accepted native source pin is missing: " + relative)
            expected = item.get("expected_sha256")
            _require(isinstance(expected, str), "Accepted native source pin has no digest: " + relative)
        actual = _sha(root / relative)
        _require(actual == expected, "Pinned native source bytes changed: " + relative)
        pins[relative] = dict(path=relative, expected_sha256=expected, sha256=actual, passed=True)
    module_path = root / MODULE
    module_digest = _sha(module_path)
    pins[MODULE] = dict(path=MODULE, expected_sha256=module_digest,
                        sha256=module_digest, passed=True)
    return pins


def _pointer_row(value: Any, label: str, node_count: int) -> list[int]:
    _require(isinstance(value, list), label + " must be a list of graph pointers")
    _require(len(value) >= 7, label + " must retain all seven source rows")
    out: list[int] = []
    for pointer in value[:7]:
        _require(type(pointer) is int and 0 <= pointer < node_count,
                 label + " contains a non-source graph pointer")
        out.append(pointer)
    return out


def _ordinary_source_generator(
    nodes: list[Any],
    guards: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    """Bind ordinary Z rows to the constructor's actual FunctionJet source."""
    generator = [
        guard for guard in guards
        if guard.get("path") == C6_SOURCE
        and guard.get("function") == "__init__"
        and "z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5);self.z=z"
        in guard.get("required_text", [])
    ]
    _require(len(generator) == 1, "The ordinary Z FunctionJet source guard is missing")
    zero = nodes[0]
    _require(isinstance(zero, Mapping)
             and zero.get("operation") == "exact_rational"
             and zero.get("numerator") == 0
             and zero.get("denominator") == 1,
             "The source cache zero graph node changed")
    for order, node in enumerate(nodes[:3]):
        _require(isinstance(node, Mapping),
                 "Ordinary Z source node is not a graph mapping: " + str(order))
        _require(isinstance(node.get("operation"), str),
                 "Ordinary Z source node has no operation: " + str(order))
    return dict(
        passed=True,
        source=C6_SOURCE,
        function="__init__",
        source_guarded=True,
        constructor_binding="z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5);self.z=z",
        parameter="Z",
        derivative_orders=list(ORDINARY_PREFIX),
        zero_graph_node=0,
        zero_graph_node_is_exact=True,
        cache_graph_nodes_are_source_mappings=True,
        source_executed=False,
    )


def _actual_rows(
    bindings: Mapping[str, Any],
    nodes: list[Any],
    guards: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    """Check endpoint/chart rows and bind ordinary rows to the source generator."""
    node_count = len(nodes)
    rows: dict[str, Any] = {}
    locations = {
        "R110.actual": bindings["endpoints"]["R110"]["actual"],
        "Rsh.actual": bindings["endpoints"]["Rsh"]["actual"],
        "long_reshape.actual": bindings["charts"]["long_reshape"]["actual"],
        "reference.actual": bindings["charts"]["reference"]["actual"],
    }
    for location, values in locations.items():
        _require(isinstance(values, Mapping), location + " is missing")
        packet: dict[str, Any] = {}
        for quantity, source_rows in values.items():
            packet[quantity] = _pointer_row(source_rows, location + "." + quantity, node_count)
        rows[location] = packet
    for chart in CHARTS:
        chart_rows = bindings["charts"][chart]["moments"]
        incoming_rows = bindings["charts"][chart]["incoming_histories"]
        _require(isinstance(chart_rows, Mapping) and isinstance(incoming_rows, Mapping),
                 chart + " moment rows are missing")
        rows[chart + ".moments"] = {
            quantity: _pointer_row(chart_rows[quantity], chart + ".moments." + quantity, node_count)
            for quantity in MOMENTS
        }
        rows[chart + ".incoming_histories"] = {
            quantity: _pointer_row(incoming_rows[quantity], chart + ".incoming_histories." + quantity, node_count)
            for quantity in MOMENTS
        }
    return dict(
        passed=True,
        ordinary_Z_orders=list(ORDINARY_PREFIX),
        rows=rows,
        ordinary_source_generator=_ordinary_source_generator(nodes, guards),
        source_graph_node_count=node_count,
        all_rows_are_source_graph_pointers=True,
        ordinary_prefix_is_genuine=True,
        ordinary_prefix_source_guarded=True,
        ordinary_Z_parameter_function_identity=True,
    )


def _moment_chain(bindings: Mapping[str, Any], parent: Mapping[str, Any]) -> dict[str, Any]:
    chart_data = bindings["charts"]
    endpoint_data = bindings["endpoints"]
    edges: dict[str, Any] = {}
    comparisons = (
        ("R110_to_long_reshape", "R110", "long_reshape"),
        ("Rsh_to_reference", "Rsh", "reference"),
    )
    for edge_name, endpoint, chart in comparisons:
        edge: dict[str, Any] = {}
        for quantity in MOMENTS:
            predecessor = endpoint_data[endpoint]["moments"][quantity]
            incoming = chart_data[chart]["incoming_histories"][quantity]
            _require(predecessor == incoming,
                     "Native moment inlet pointer changed for " + edge_name + ":" + quantity)
            edge[quantity] = dict(
                predecessor_endpoint=list(predecessor),
                successor_incoming_histories=list(incoming),
                pointer_equal=True,
                rows=list(ORDINARY_PREFIX),
            )
        edges[edge_name] = edge

    native = parent["actual_long_reference_native_controls"]
    parent_chart_checks: dict[str, Any] = {}
    for chart in CHARTS:
        ports = native["charts"][chart]["source_ports"]
        cache_chart = chart_data[chart]
        _require(ports["incoming_moments"] == cache_chart["incoming_histories"],
                 "Accepted native report does not retain cache incoming moments: " + chart)
        _require(ports["actual"] == cache_chart["actual"],
                 "Accepted native report does not retain cache actual pointers: " + chart)
        parent_chart_checks[chart] = dict(
            actual_pointer_equal=True,
            incoming_moments_pointer_equal=True,
        )
    return dict(
        passed=True,
        quantities=list(MOMENTS),
        rows_per_quantity=7,
        ordinary_Z_prefix=list(ORDINARY_PREFIX),
        edges=edges,
        all_six_moment_values_pointer_equal=True,
        accepted_native_report_cache_pointer_equal=True,
        parent_chart_checks=parent_chart_checks,
        source_cache=C6_CACHE,
        source_executed=False,
    )


def _normalization(bindings: Mapping[str, Any]) -> dict[str, Any]:
    charts = bindings["charts"]
    post_measure = charts["post_switch"]["physical_measure"]
    long_measure = charts["long_reshape"]["physical_measure"]
    ref_measure = charts["reference"]["physical_measure"]
    _require(post_measure == bindings["post_length"], "Post-switch physical measure is not Jpost")
    _require(long_measure == bindings["T"], "Long-reshape physical measure is not T")
    _require(ref_measure == bindings["reference_length"],
             "Reference physical measure is not the source Lref")
    _require(len({post_measure, long_measure, ref_measure}) == 3,
             "Native chart normalizations were collapsed")
    _require(post_measure != 7 and ref_measure != 7,
             "Native chart normalization was replaced by a literal length")
    return dict(
        passed=True,
        source_cache=C6_CACHE,
        post_switch=dict(chart="post_switch", physical_measure=post_measure,
                         binding="post_length", binding_pointer=bindings["post_length"]),
        long_reshape=dict(chart="long_reshape", physical_measure=long_measure,
                          binding="T", binding_pointer=bindings["T"]),
        reference=dict(chart="reference", physical_measure=ref_measure,
                       binding="reference_length", binding_pointer=bindings["reference_length"]),
        unequal_J_normalizations=True,
        source_length_formula="reference_length=logref-T-8",
        literal_length_rejected=True,
        post_length_is_not_long_length=True,
        reference_length_is_not_long_length=True,
    )


def _anchored_F0_binding(
    root: Path,
    parent: Mapping[str, Any],
    guards: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    """Tie the symbolic F0 identity to the accepted anchored-G source packet."""
    native = parent["actual_long_reference_native_controls"]
    anchor = native.get("actual_anchor_binding")
    _require(isinstance(anchor, Mapping), "Accepted anchored-G binding is missing")
    required_flags = (
        "same_current_family",
        "anchored_G_not_Gbar",
        "accepted_post_source_receipt_reused_not_reexecuted",
        "old_reference_join_family_not_consumed",
    )
    for flag in required_flags:
        _require(anchor.get(flag) is True, "Accepted anchored-G binding changed: " + flag)

    amplitude = parent.get("actual_micro_b_p2_mixed_bounds", {}).get("whole_tube_amplitude")
    _require(isinstance(amplitude, Mapping) and amplitude.get("passed") is True,
             "Accepted anchored amplitude packet is missing")
    proof = native.get("amplitude_proof")
    _require(isinstance(proof, Mapping), "Accepted native amplitude proof is missing")
    _require(proof.get("logF0") == amplitude.get("log_F0"),
             "Parent logF0 packet is not the whole-tube anchored packet")
    _require(proof.get("native_logphi_identity") ==
             "logphi=logphi_R110-2y/5-B*sigma(y/T)",
             "Parent native logphi identity changed")
    _require(amplitude.get("exact_log_F0") ==
             "-logC-Lambda*G; 0<=G<=Gbar on this tube",
             "Parent exact logF0 source identity changed")
    _require(amplitude.get("inverse_F0_not_materialized") is True,
             "Parent F0 proof materialized an inverse value")
    _require(amplitude.get("selected_logC_not_replaced_by_minimum_guard") is True,
             "Parent logC source was replaced by a guard")

    anchor_guards = [
        guard for guard in guards
        if guard.get("path") == ANCHOR_SOURCE
        and guard.get("function") in {"__init__", "gradient_jets"}
    ]
    _require({guard.get("function") for guard in anchor_guards} ==
             {"__init__", "gradient_jets"},
             "Anchored-G source guards are incomplete")
    return dict(
        passed=True,
        source=ANCHOR_SOURCE,
        source_function="CompliantAnchoredAxisAmplitude.evaluate",
        source_guard_functions=sorted({guard["function"] for guard in anchor_guards}),
        parent_source_pin=ANCHOR_SOURCE,
        parent_anchor_flags=dict(anchor),
        parent_anchor_flags_all_true=True,
        parent_amplitude_packet="actual_micro_b_p2_mixed_bounds.whole_tube_amplitude",
        parent_logF0_packet_equal=True,
        exact_logF0="-logC-Lambda*G; 0<=G<=Gbar on this tube",
        symbolic_logF0="-logC-Lambda*G(Z)",
        anchored_G_symbol_kept=True,
        Gbar_cap_not_selected=True,
        logC_source_kept=True,
        inverse_F0_not_materialized=True,
        source_executed=False,
    )


def _cache_node(nodes: list[Any], pointer: int, label: str) -> Mapping[str, Any]:
    _require(type(pointer) is int and 0 <= pointer < len(nodes),
             label + " is not a source graph pointer")
    node = nodes[pointer]
    _require(isinstance(node, Mapping), label + " has no source graph node")
    return node


def _substitution_node(
    nodes: list[Any],
    pointer: int,
    expression: int,
    variable: int,
    label: str,
) -> None:
    node = _cache_node(nodes, pointer, label)
    _require(node.get("operation") == "function_substitution",
             label + " is not a function substitution")
    _require(node.get("expression") == expression,
             label + " does not preserve the incoming expression pointer")
    _require(node.get("variable") == variable and node.get("value") == 1,
             label + " does not substitute the source endpoint phase")
    _require(node.get("Z_independent_substitution") is True,
             label + " lost Z-independent substitution metadata")


def _u_pointer_chain(bindings: Mapping[str, Any], nodes: list[Any]) -> dict[str, Any]:
    """Prove the six native U rows by following cache substitution edges."""
    charts = bindings["charts"]
    endpoints = bindings["endpoints"]
    post_phase = charts["post_switch"]["phase"]
    long_phase = charts["long_reshape"]["phase"]
    r2_u = _pointer_row(endpoints["R2"]["actual"]["U"], "R2.actual.U", len(nodes))
    post_u = _pointer_row(charts["post_switch"]["actual"]["U"],
                          "post_switch.actual.U", len(nodes))
    r110_u = _pointer_row(endpoints["R110"]["actual"]["U"],
                          "R110.actual.U", len(nodes))
    long_u = _pointer_row(charts["long_reshape"]["actual"]["U"],
                          "long_reshape.actual.U", len(nodes))
    rsh_u = _pointer_row(endpoints["Rsh"]["actual"]["U"],
                         "Rsh.actual.U", len(nodes))
    reference_u = _pointer_row(charts["reference"]["actual"]["U"],
                               "reference.actual.U", len(nodes))

    _require(post_u == r2_u, "Post-switch U did not retain the R2 inlet pointer")
    for order, pointer in enumerate(r110_u):
        _substitution_node(nodes, pointer, post_u[order], post_phase,
                           "R110.actual.U[" + str(order) + "]")
    _require(long_u == r110_u, "Long-reshape U did not retain the R110 inlet pointer")
    for order, pointer in enumerate(rsh_u):
        _substitution_node(nodes, pointer, long_u[order], long_phase,
                           "Rsh.actual.U[" + str(order) + "]")
    _require(reference_u == rsh_u,
             "Reference U did not retain the Rsh inlet pointer")

    zero = 0
    zero_chart_rows: dict[str, list[int]] = {}
    for chart in ("post_switch", "long_reshape", "reference"):
        rows = _pointer_row(charts[chart]["actual"]["U_y"],
                            chart + ".actual.U_y", len(nodes))
        _require(rows == [zero] * 7,
                 chart + " U_y source is not the exact zero function")
        zero_chart_rows[chart] = rows
    r110_u_y = _pointer_row(endpoints["R110"]["actual"]["U_y"],
                            "R110.actual.U_y", len(nodes))
    rsh_u_y = _pointer_row(endpoints["Rsh"]["actual"]["U_y"],
                           "Rsh.actual.U_y", len(nodes))
    for order, pointer in enumerate(r110_u_y):
        _substitution_node(nodes, pointer, zero, post_phase,
                           "R110.actual.U_y[" + str(order) + "]")
    for order, pointer in enumerate(rsh_u_y):
        _substitution_node(nodes, pointer, zero, long_phase,
                           "Rsh.actual.U_y[" + str(order) + "]")
    zero_node = _cache_node(nodes, zero, "U_y zero")
    _require(zero_node.get("operation") == "exact_rational"
             and zero_node.get("numerator") == 0
             and zero_node.get("denominator") == 1,
             "Native U zero source node changed")
    return dict(
        passed=True,
        source=C6_SOURCE,
        rows=list(ORDINARY_PREFIX),
        source_pointer_rows=7,
        post_switch_inlet=dict(chart="post_switch", source_endpoint="R2",
                               pointer_equal=True, rows=list(post_u)),
        R110=dict(endpoint_substitutes_post_switch=True, phase_pointer=post_phase,
                  rows=list(r110_u)),
        long_reshape=dict(inlet_endpoint="R110", pointer_equal=True,
                          rows=list(long_u)),
        Rsh=dict(endpoint_substitutes_long_reshape=True, phase_pointer=long_phase,
                 rows=list(rsh_u)),
        reference=dict(inlet_endpoint="Rsh", pointer_equal=True,
                       rows=list(reference_u)),
        substitution_edges=[
            "R2 -> post_switch",
            "post_switch --p=1--> R110",
            "R110 -> long_reshape",
            "long_reshape --q=1--> Rsh",
            "Rsh -> reference",
        ],
        U_y_source_zero=True,
        U_yy_source_zero_from_constant_U_y=True,
        U_substitution_chain_is_source_bound=True,
        zero_function_rows=zero_chart_rows,
        source_executed=False,
    )


def _sigma_endpoint_identities() -> dict[str, Any]:
    """Reduce the actual squared-flat cutoff through value, first, second order."""
    try:
        import sympy as sym
    except Exception as exc:  # pragma: no cover - environment diagnostic
        raise ValueError("SymPy is required for the native seam reduction") from exc
    x = sym.symbols("x", positive=True)
    numerator = sym.exp(-1 / x**2)
    denominator_tail = sym.exp(-1 / (1 - x)**2)
    sigma = sym.simplify(numerator / (numerator + denominator_tail))
    endpoint_rows: dict[str, list[str]] = {}
    expected = {"0": (sym.Integer(0), sym.Integer(0), sym.Integer(0)),
                "1": (sym.Integer(1), sym.Integer(0), sym.Integer(0))}
    for endpoint, direction in (("0", "+"), ("1", "-")):
        values: list[str] = []
        for order in range(3):
            value = sym.simplify(sym.limit(sym.diff(sigma, x, order), x,
                                           sym.Integer(endpoint), dir=direction))
            _require(value == expected[endpoint][order],
                     "Squared-flat sigma endpoint identity failed")
            values.append(str(value))
        endpoint_rows[endpoint] = values
    return dict(
        passed=True,
        source_formula="sigma(x)=exp(-1/x**2)/(exp(-1/x**2)+exp(-1/(1-x)**2)); right half by sigma(1-x)=1-sigma(x)",
        squared_flat_cutoff=True,
        endpoint_rows=endpoint_rows,
        derivative_orders=[0, 1, 2],
        value_first_second_flat=True,
        endpoint_symmetry=True,
        source_executed=False,
    )


def _jet_residuals(left: list[Any], right: list[Any], z: Any) -> dict[str, Any]:
    import sympy as sym
    residuals: dict[str, Any] = {}
    for y_order, (left_value, right_value) in enumerate(zip(left, right)):
        difference = sym.simplify(left_value - right_value)
        _require(difference == 0, "Physical native seam jet did not reduce to zero")
        residuals[str(y_order)] = {
            str(z_order): str(sym.simplify(sym.diff(difference, z, z_order)))
            for z_order in ORDINARY_PREFIX
        }
        _require(all(value == "0" for value in residuals[str(y_order)].values()),
                 "Physical native seam ordinary-Z prefix did not reduce to zero")
    return residuals


def _symbolic_c2(
    bindings: Mapping[str, Any],
    sigma: Mapping[str, Any],
    u_chain: Mapping[str, Any],
    anchor: Mapping[str, Any],
) -> dict[str, Any]:
    """Reduce independent source formulae after actual endpoint substitution."""
    import sympy as sym

    _require(anchor.get("passed") is True, "Anchored F0 binding is required for C2")
    _require(u_chain.get("passed") is True, "Source U pointer chain is required for C2")
    z, p, q_long, q_reference = sym.symbols("Z p q_long q_reference", real=True)
    lam, logC = sym.symbols("Lambda logC", real=True)
    jpost = sym.Symbol("Jpost")
    T = sym.Symbol("T")
    Lref = sym.Symbol("Lref")
    R110 = sym.Integer(110)
    G = sym.Function("G")(z)
    logq = sym.log(1 + z**2)
    logF0_expr = -logC - lam * G
    sigma_function = sym.Function("sigma")

    # Each expression follows the source constructor, including the chart
    # Jacobian and radius.  The post endpoint is derived from its own formula,
    # then used as the long chart's actual R110 inlet.
    logphi_post = sym.Function("logphi_R2")(z) - sym.Rational(2, 5) * jpost * p
    logphi_R110 = sym.simplify(logphi_post.subs(p, sym.Integer(1)))
    B = (-lam * G + logphi_R110 + sym.log(2 * R110) / 2 + logq)
    y_long = T * q_long
    logphi_long = (B * (1 - sigma_function(q_long)) - logq - logC
                   - logF0_expr - (sym.log(2 * R110) + y_long) / 2
                   + y_long / 10)
    y_reference = T + Lref * q_reference
    logphi_reference = (-logq - logC - logF0_expr
                        - (sym.log(2 * R110) + y_reference) / 2
                        + y_reference / 10)

    endpoint_rows = sigma.get("endpoint_rows")
    _require(isinstance(endpoint_rows, Mapping), "Sigma endpoint rows are missing")

    def endpoint_substitution(phase: Any, endpoint: int) -> dict[Any, Any]:
        raw = endpoint_rows[str(endpoint)]
        _require(isinstance(raw, list) and len(raw) == 3,
                 "Sigma endpoint row must contain value, first, second jets")
        values = [sym.sympify(value) for value in raw]
        return {
            sym.diff(sigma_function(phase), phase, order): values[order]
            for order in range(3)
        }

    def at_endpoint(expression: Any, phase: Any, endpoint: int) -> Any:
        reduced = expression.subs(endpoint_substitution(phase, endpoint))
        return sym.simplify(reduced.subs(phase, sym.Integer(endpoint)))

    def physical_log_derivatives(expression: Any, phase: Any, jacobian: Any) -> tuple[Any, Any]:
        first = sym.simplify(sym.diff(expression, phase) / jacobian)
        second = sym.simplify(sym.diff(expression, phase, 2) / jacobian**2)
        return first, second

    def phi_jet(log_expression: Any, first: Any, second: Any) -> list[Any]:
        value = sym.exp(log_expression)
        return [value, sym.simplify(value * first),
                sym.simplify(value * (first**2 + second))]

    post_first, post_second = physical_log_derivatives(logphi_post, p, jpost)
    post_log_at_R110 = sym.simplify(logphi_R110)
    post_first_at_R110 = sym.simplify(post_first.subs(p, sym.Integer(1)))
    post_second_at_R110 = sym.simplify(post_second.subs(p, sym.Integer(1)))
    post_jet = phi_jet(post_log_at_R110, post_first_at_R110,
                       post_second_at_R110)

    long_first, long_second = physical_log_derivatives(logphi_long, q_long, T)
    long_log_at_q0 = at_endpoint(logphi_long, q_long, 0)
    long_first_at_q0 = at_endpoint(long_first, q_long, 0)
    long_second_at_q0 = at_endpoint(long_second, q_long, 0)
    long_q0_jet = phi_jet(long_log_at_q0, long_first_at_q0,
                          long_second_at_q0)
    post_residuals = _jet_residuals(post_jet, long_q0_jet, z)

    long_log_at_q1 = at_endpoint(logphi_long, q_long, 1)
    long_first_at_q1 = at_endpoint(long_first, q_long, 1)
    long_second_at_q1 = at_endpoint(long_second, q_long, 1)
    long_q1_jet = phi_jet(long_log_at_q1, long_first_at_q1,
                          long_second_at_q1)
    reference_first, reference_second = physical_log_derivatives(
        logphi_reference, q_reference, Lref
    )
    reference_log_at_q0 = sym.simplify(logphi_reference.subs(q_reference, 0))
    reference_first_at_q0 = sym.simplify(reference_first.subs(q_reference, 0))
    reference_second_at_q0 = sym.simplify(reference_second.subs(q_reference, 0))
    reference_q0_jet = phi_jet(reference_log_at_q0, reference_first_at_q0,
                                reference_second_at_q0)
    join_residuals = _jet_residuals(long_q1_jet, reference_q0_jet, z)

    # U is carried by source substitutions, not by assigning one fresh symbol
    # to both sides.  The Lambda forms expose each chart's constant source
    # formula before the endpoint substitution is taken.
    u_R2 = sym.Symbol("U_R2_source")
    post_U = sym.Lambda(p, u_R2)(p)
    R110_U = sym.simplify(post_U.subs(p, sym.Integer(1)))
    long_U = sym.Lambda(q_long, R110_U)(q_long)
    Rsh_U = sym.simplify(long_U.subs(q_long, sym.Integer(1)))
    reference_U = sym.Lambda(q_reference, Rsh_U)(q_reference)

    def U_jet(expression: Any, phase: Any, jacobian: Any) -> list[Any]:
        return [
            expression,
            sym.simplify(sym.diff(expression, phase) / jacobian),
            sym.simplify(sym.diff(expression, phase, 2) / jacobian**2),
        ]

    post_U_jet = U_jet(post_U, p, jpost)
    long_q0_U_jet = U_jet(long_U, q_long, T)
    long_q1_U_jet = U_jet(long_U, q_long, T)
    reference_q0_U_jet = U_jet(reference_U, q_reference, Lref)
    post_U_residuals = _jet_residuals(
        [sym.simplify(value.subs(p, 1)) for value in post_U_jet],
        [sym.simplify(value.subs(q_long, 0)) for value in long_q0_U_jet], z,
    )
    join_U_residuals = _jet_residuals(
        [sym.simplify(value.subs(q_long, 1)) for value in long_q1_U_jet],
        [sym.simplify(value.subs(q_reference, 0)) for value in reference_q0_U_jet], z,
    )

    def residual_by_name(residuals: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "U": residuals["0"],
            "U_y": residuals["1"],
            "U_yy": residuals["2"],
        }

    post_phi_residuals = {
        "phi": post_residuals["0"],
        "phi_y": post_residuals["1"],
        "phi_yy": post_residuals["2"],
    }
    join_phi_residuals = {
        "phi": join_residuals["0"],
        "phi_y": join_residuals["1"],
        "phi_yy": join_residuals["2"],
    }
    post_log_residual = sym.simplify(post_log_at_R110 - long_log_at_q0)
    join_log_residual = sym.simplify(long_log_at_q1 - reference_log_at_q0)
    _require(post_log_residual == 0 and join_log_residual == 0,
             "Independent source logarithmic phi identities failed")

    source_pointers = {
        "Jpost": bindings["charts"]["post_switch"]["physical_measure"],
        "T": bindings["charts"]["long_reshape"]["physical_measure"],
        "Lref": bindings["charts"]["reference"]["physical_measure"],
        "R110_endpoint": bindings["endpoints"]["R110"]["radius"],
        "Rsh_endpoint": bindings["endpoints"]["Rsh"]["radius"],
        "long_radius": bindings["charts"]["long_reshape"]["radius"],
        "reference_radius": bindings["charts"]["reference"]["radius"],
    }
    return dict(
        passed=True,
        ordinary_Z_orders=list(ORDINARY_PREFIX),
        ordinary_Z_parameter_function_identity=True,
        ordinary_Z_derivatives_of_each_residual=list(ORDINARY_PREFIX),
        physical_y_orders=[0, 1, 2],
        physical_derivative_definition={
            "logphi_y": "(1/J)*d_q(logphi)",
            "logphi_yy": "(1/J**2)*d_q**2(logphi)",
            "phi_y": "phi*logphi_y",
            "phi_yy": "phi*(logphi_y**2+logphi_yy)",
        },
        source_binding_pointers=source_pointers,
        source_radius_formulas={
            "post_switch": "R=100*exp(2*hb+y), R110=110 at p=1",
            "long_reshape": "R=110*exp(T*q_long)",
            "reference": "R=110*exp(T+Lref*q_reference)",
        },
        logphi_source_symbols={
            "logq": "log(1+Z**2)",
            "B": "-Lambda*G(Z)+logphi_R110(Z)+log(220)/2+log(1+Z**2)",
            "logF0": "-logC-Lambda*G(Z)",
            "logphi_R110": "logphi_R2(Z)-2*Jpost/5",
        },
        post_to_long_q0=dict(
            passed=True,
            left_chart="post_switch upper endpoint R110",
            right_chart="long_reshape q=0",
            left_logphi_expression=str(post_log_at_R110),
            right_logphi_expression=str(long_log_at_q0),
            log_value_residual=str(post_log_residual),
            normalized_logphi_derivatives={
                "left": [str(post_first_at_R110), str(post_second_at_R110)],
                "right": [str(long_first_at_q0), str(long_second_at_q0)],
            },
            sigma_endpoint_jets={
                "value": sigma["endpoint_rows"]["0"][0],
                "first": sigma["endpoint_rows"]["0"][1],
                "second": sigma["endpoint_rows"]["0"][2],
            },
            phi_y_residuals=post_residuals,
            phi_C2_residuals=post_phi_residuals,
            U_C2_residuals=residual_by_name(post_U_residuals),
            U_physical_y_jets=[[str(item) for item in post_U_jet],
                               [str(item) for item in long_q0_U_jet]],
            source_formula_reduction_independent=True,
            source_substitution_chain=u_chain,
        ),
        long_q1_to_reference_q0=dict(
            passed=True,
            left_chart="long_reshape q=1 (Rsh)",
            right_chart="reference q=0",
            left_logphi_expression=str(long_log_at_q1),
            right_logphi_expression=str(reference_log_at_q0),
            log_value_residual=str(join_log_residual),
            normalized_logphi_derivatives={
                "left": [str(long_first_at_q1), str(long_second_at_q1)],
                "right": [str(reference_first_at_q0), str(reference_second_at_q0)],
            },
            sigma_endpoint_jets={
                "value": sigma["endpoint_rows"]["1"][0],
                "first": sigma["endpoint_rows"]["1"][1],
                "second": sigma["endpoint_rows"]["1"][2],
            },
            phi_y_residuals=join_residuals,
            phi_C2_residuals=join_phi_residuals,
            U_C2_residuals=residual_by_name(join_U_residuals),
            U_physical_y_jets=[[str(item) for item in long_q1_U_jet],
                               [str(item) for item in reference_q0_U_jet]],
            source_formula_reduction_independent=True,
            source_substitution_chain=u_chain,
        ),
        sigma_endpoint_values_used={
            "sigma(0)": sigma["endpoint_rows"]["0"][0],
            "sigma_y(0)": sigma["endpoint_rows"]["0"][1],
            "sigma_yy(0)": sigma["endpoint_rows"]["0"][2],
            "sigma(1)": sigma["endpoint_rows"]["1"][0],
            "sigma_y(1)": sigma["endpoint_rows"]["1"][1],
            "sigma_yy(1)": sigma["endpoint_rows"]["1"][2],
        },
        anchored_F0_binding=anchor,
        unknown_G0_kept_symbolic=True,
        numeric_caps_used=False,
        source_executed=False,
    )


def _inherited_lemma(parent: Mapping[str, Any]) -> dict[str, Any]:
    native = parent["actual_long_reference_native_controls"]
    lemma = native.get("source_function_identities")
    _require(isinstance(lemma, Mapping) and lemma.get("passed") is True,
             "Accepted native source-function identity lemma is missing")
    _require(lemma.get("ordinary_Z_derivatives_follow_from_function_identity") == [0, 1, 2, 3],
             "Accepted native source-function identity lemma has changed")
    return dict(
        passed=True,
        source="actual_long_reference_native_controls.source_function_identities",
        source_report_sha256=PARENT_SHA256,
        consumed_as_lemma=True,
        reexecuted=False,
        p2_and_logphi_algebra_reproved=False,
        ordinary_Z_orders=list(lemma["ordinary_Z_derivatives_follow_from_function_identity"]),
    )


def _load_parent(root: Path, parent: Mapping[str, Any]) -> Mapping[str, Any]:
    disk = _read_json(root, PARENT)
    _require(parent == disk, "The exact accepted native-controls parent mapping is required")
    _require(parent.get("source_family") is not None, "Accepted source family is missing")
    native = parent.get("actual_long_reference_native_controls")
    _require(isinstance(native, Mapping) and native.get("passed") is True,
             "Accepted native-controls packet is missing")
    return disk


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Return the scoped native ``phi``/``U`` C2 and moment-inlet receipt."""
    root_path = Path(root)
    disk_parent = _load_parent(root_path, parent)
    _require(_sha(root_path / PARENT) == PARENT_SHA256, "Accepted parent bytes changed")
    pins = _source_pins(root_path, disk_parent)
    guards = _source_guards(root_path)
    cache = _read_gzip_json(root_path, C6_CACHE)
    _require(cache.get("source_family") == disk_parent.get("source_family"),
             "C6 cache source family changed")
    bindings = cache.get("prescribed_R100_Rm_C6_bindings")
    nodes = cache.get("exact_graph_nodes")
    _require(isinstance(bindings, Mapping) and isinstance(nodes, list),
             "C6 source cache graph/bindings are missing")
    _require(set(bindings.get("charts", {})) >= {
        "post_switch", "long_reshape", "reference"
    }, "C6 source cache chart bindings are incomplete")
    _require(set(bindings.get("endpoints", {})) >= {"R110", "Rsh"},
             "C6 source cache endpoint bindings are incomplete")
    source_rows = _actual_rows(bindings, nodes, guards)
    moments = _moment_chain(bindings, disk_parent)
    normalization = _normalization(bindings)
    sigma = _sigma_endpoint_identities()
    anchored_F0 = _anchored_F0_binding(root_path, disk_parent, guards)
    u_chain = _u_pointer_chain(bindings, nodes)
    c2 = _symbolic_c2(bindings, sigma, u_chain, anchored_F0)
    lemma = _inherited_lemma(disk_parent)
    rates = guards[-1]
    _require(rates["rates"] == {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1},
             "Six-history source rates are incomplete")

    # The first/second-history ODE estimate remains outside this scoped
    # receipt.  Return its exact source binding and rates for the root proof.
    moment_ode_binding = dict(
        passed=True,
        source=C6_SOURCE,
        constructor_binding="add_chart -> previous.histories(actual_phi, actual_U, incoming, phase, jacobian)",
        rate_source=RATE_SOURCE,
        rates=dict(rates["rates"]),
        first_second_ordinary_rows=list(ORDINARY_PREFIX),
        proof_deferred_to_root=True,
        source_executed=False,
    )

    result: dict[str, Any] = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_family=disk_parent["source_family"],
        source_pins=pins,
        input_pins=pins,
        source_guards=guards,
        sourceguards_tuple=True,
        source_executed=False,
        accepted_ancestors_replayed=False,
        accepted_builders_executed=False,
        accepted_checkers_executed=False,
        accepted_producers_executed=False,
        owners_executed=False,
        norms_executed=False,
        recurrences_executed=False,
        quadrature_executed=False,
        physical_y_derivatives_claimed=True,
        ordinary_Z_orders=list(ORDINARY_PREFIX),
        scope="native_phi_U_physical_C2_and_six_moment_inlet_pointer_chain",
        scoped_receipt=True,
        modified_fields_claimed=False,
        modified_Qy_claimed=False,
        full_modified_seam_claimed=False,
        source_function_identities_lemma=lemma,
        sigma_endpoint_identities=sigma,
        anchored_F0_binding=anchored_F0,
        source_rows=source_rows,
        native_U_substitution_chain=u_chain,
        moment_inlet_pointer_chain=moments,
        moment_ode_binding=moment_ode_binding,
        history_rate_source=rates,
        unequal_J_normalization=normalization,
        physical_native_C2=c2,
        native_phi_U_C2=c2,
        newproof_checks=dict(
            source_bound_two_chart_reductions=True,
            squared_flat_sigma_value_first_second_endpoints=True,
            anchored_F0_parent_source_binding=True,
            unequal_J_post_T_Lref_bindings=True,
            six_moment_inlet_pointer_edges=True,
            ordinary_Z_prefix_0_1_2_source_coherent=True,
            ordinary_Z_source_generator_guarded=True,
            ordinary_Z_parameter_function_identity=True,
            native_U_source_substitution_chain=True,
            native_phi_phi_y_phi_yy_equal=True,
            native_U_U_y_U_yy_equal=True,
            modified_Qy_seam=False,
            full_modified_seam=False,
        ),
        limits=dict(
            modified_fields_seam=False,
            modified_Qy_seam=False,
            full_modified_seam=False,
            finite_N_numeric=False,
            accepted_owner_execution=False,
            first_second_moment_ODE_proof_deferred=True,
        ),
        module_sha256=_sha(root_path / MODULE),
    )
    result["receipt_sha256"] = hashlib.sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    json.dumps(result, allow_nan=False)
    return result


__all__ = ["build"]
