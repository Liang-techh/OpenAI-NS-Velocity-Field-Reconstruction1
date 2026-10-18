"""Source-bound periodic C5 loop seams for the long/reference charts.

This module reads the issued C5, C6, and Q4 source files and their frozen
graph archives as data.  It derives the exact integer phase drift of the
issued loop formula and records the two endpoint joins requested by the
long/reference source report.  It does not import or execute an archived
owner, checker, producer, norm, recurrence, quadrature, or numeric builder.

The returned certificate deliberately stops at source periodicity, graph
history binding, and the native C2/derivative source guards.  It does not
turn those guards into a modified ``U_y``/``Q_y`` closure or a full pressure
bound.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_finite_N_I1.json"
PARENT_SHA256 = "75918293b8c7432c637eb551a29cb9d248aa7746dcb1cd4f8b3615f72514c14d"
MODULE = "src/openai_ns_reconstruction/original_long_reference_loop_seams.py"

C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
C5_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
Q4_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"

PINNED_HASHES = {
    PARENT: PARENT_SHA256,
    C5_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    C5_CACHE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    C5_CHECK: "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    C6_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    C6_CACHE: "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    C6_CHECK: "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    Q4_SOURCE: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    Q4_CACHE: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    Q4_CHECK: "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f",
}

CORRECTION_RATES = ("m", "h", "k", "e", "p")
NATIVE_MOMENTS = ("H", "M", "K", "A", "B", "C")
SEAM_CHARTS = ("post_power", "long_reshape", "reference")

# Every snippet collection is a tuple.  The tuple shape is part of this
# module's source guard: a list here would make a future guard silently
# mutable.
SOURCE_GUARDS = (
    (C5_SOURCE, "loop_functions", None, (
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
        "numerator=m.add(cc,m.scale(rr,-1))",
        "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
        "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
        "t2=m.mul(t,t)",
        "ratio=a.add(one,a.mul(roots['t0'],roots['t0']),a.scale(a.mul(q,q),2))",
        "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
        "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
        "T2=FunctionJet(tuple(g.mul(c(math.factorial(k)),integrate(t2[k,0])) for k in range(6)))",
        "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
        "inverse=g.node('original_monotone_phase_inverse'",
        "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
        "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
        "chi=a.sub(a.fixed(phase),a.scale(psi,g.quotient(g.one,twopi,'positive 2*pi')))",
        "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))",
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
    )),
    (C5_SOURCE, "generic_roots", None, (
        "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)",
        "aa,bb=divide(C),divide(B)",
        "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
    )),
    (C5_SOURCE, "__init__", "CurrentOriginalPrescribedLoopC5Functions", (
        "self.phase_origin=g.mul(background.bridge.hb,self.sc,g.constant('1/2'))",
        "self.source_charts={",
        "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)",
        "phase=g.unary('fractional_part',g.mul(self.N,offset))",
        "history,parts,memory=transport(g,finite['five_signed_density_increments'],incoming,",
    )),
    (C5_SOURCE, "transport", None, (
        "for name,rate in RATES.items():",
        "predecessor_memory_never_reset=True",
        "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
    )),
    (C6_SOURCE, "__init__", "CurrentOriginalPrescribedR100RmC6Functions", (
        "self.post_length=g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))",
        "y=g.mul(self.post_length,p);R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))",
        "q=g.symbol('true_prescribed_long_reshape_fraction');y=g.mul(T,q);R=g.mul(c(110),g.unary('exp',y))",
        "self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])",
        "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
        "self.reference_length=length=g.sub(g.sub(logref,T),c(8))",
        "q=g.symbol('true_prescribed_reference_fraction');y=g.add(T,g.mul(length,q));R=g.mul(c(110),g.unary('exp',y))",
        "incoming=self.endpoints['Rsh']",
        "self.add_chart('reference',q,R,length,actual,incoming['moments'])",
    )),
    (C6_SOURCE, "add_chart", "CurrentOriginalPrescribedR100RmC6Functions", (
        "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
        "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
    )),
    (Q4_SOURCE, "derivative", "NativeDerivative", (
        "elif op=='definite_integral':",
        "elif op=='function_substitution':",
        "elif op=='original_monotone_phase_inverse':",
        "value=g.quotient(g.sub(self.derivative(target),at(held)),lam,'original positive phase lambda; exact implicit native derivative')",
    )),
    (Q4_SOURCE, "total_radial_functions", None, (
        "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))",
        "xy=a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
        "physical_y_equals_native_partial_divided_by_Jacobian=True",
    )),
    (Q4_SOURCE, "physical_pressure_fields", None, (
        "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
        "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)",
        "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
    )),
    (Q4_SOURCE, "__init__", "CurrentOriginalModifiedQ4PressureC3Functions", (
        "derivative=NativeDerivative(g,coordinate,chart['actual_phase']);self.derivatives[name]=derivative",
        "radial=total_radial_functions(g,chart,derivative,loop_owner.N)",
        "physical=physical_pressure_fields(g,chart['modified_fields'],radial,self.background.delta,self.background.Pstar)",
    )),
)


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path.name}")
    return value


def _gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path.name}")
    return value


def _compact(value: str) -> str:
    return "".join(value.split())


def _ast_target(tree: ast.AST, name: str, owner: str | None) -> ast.AST | None:
    if owner is None:
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
                return node
        return None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == owner:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name:
                    return child
            if name == owner:
                return node
    return None


def _source_guard(root: Path, path_name: str, name: str, owner: str | None,
                  snippets: tuple[str, ...]) -> dict[str, Any]:
    _require(isinstance(snippets, tuple), f"source snippets must be tuple: {name}")
    path = root / path_name
    raw = path.read_bytes()
    _require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
             f"strict LF single EOF guard failed: {path_name}")
    source = raw.decode("utf-8")
    tree = ast.parse(source, filename=str(path))
    node = _ast_target(tree, name, owner)
    _require(node is not None, f"AST target missing: {path_name}:{owner or ''}{name}")
    segment = ast.get_source_segment(source, node) or ""
    compact = _compact(segment)
    matched = all(_compact(token) in compact for token in snippets)
    _require(matched, f"source guard failed: {path_name}:{owner or ''}{name}")
    return {
        "path": path_name,
        "target": name,
        "owner": owner,
        "tokens": snippets,
        "token_count": len(snippets),
        "AST_sha256": hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        "passed": True,
        "source_executed": False,
    }


def _source_guards(root: Path) -> tuple[dict[str, Any], ...]:
    guards = tuple(_source_guard(root, path, name, owner, snippets)
                   for path, name, owner, snippets in SOURCE_GUARDS)
    _require(isinstance(guards, tuple), "source guards must remain a tuple")
    return guards


def _pin(root: Path, path_name: str, expected: str) -> dict[str, Any]:
    path = root / path_name
    _require(path.is_file(), f"pinned source missing: {path_name}")
    actual = _sha(path)
    _require(actual == expected, f"pinned source changed: {path_name}")
    return {"path": path_name, "expected_sha256": expected, "actual_sha256": actual,
            "passed": True, "source_executed": False}


def _load_parent(root: Path, parent: Mapping[str, Any] | str | Path | None) -> dict[str, Any]:
    path = root / PARENT
    _require(path.is_file() and _sha(path) == PARENT_SHA256, "pinned parent report changed")
    disk = _json(path)
    if parent is None:
        return disk
    if isinstance(parent, Mapping):
        _require(dict(parent) == disk, "exact accepted parent mapping is required")
        return dict(parent)
    parent_path = Path(parent)
    if not parent_path.is_absolute():
        parent_path = root / parent_path
    candidate = _json(parent_path)
    _require(candidate == disk, "parent path does not contain the pinned mapping")
    return disk


def _source_pins(root: Path) -> dict[str, dict[str, Any]]:
    pins = {path: _pin(root, path, digest) for path, digest in PINNED_HASHES.items()}
    module_path = root / MODULE
    _require(module_path.is_file(), "new seam module is missing")
    module_raw = module_path.read_bytes()
    _require(b"\r" not in module_raw and module_raw.endswith(b"\n") and not module_raw.endswith(b"\n\n"),
             "new seam module must use strict LF and one EOF newline")
    module_digest = _sha(module_path)
    pins[MODULE] = {"path": MODULE, "expected_sha256": module_digest,
                    "actual_sha256": module_digest, "passed": True, "source_executed": False}
    return pins


def _nodes(archive: Mapping[str, Any], key: str, length: int) -> tuple[list[dict[str, Any]], Mapping[str, Any]]:
    nodes = archive.get("exact_graph_nodes")
    bindings = archive.get(key)
    _require(isinstance(nodes, list) and len(nodes) == length, f"{key} exact graph length mismatch")
    _require(isinstance(bindings, Mapping), f"{key} bindings missing")
    return nodes, bindings


def _copy_json(value: Any) -> Any:
    # Archives contain only JSON values.  This keeps receipts detached from
    # the mutable input object while preserving tuple guards in the live map.
    return json.loads(json.dumps(value, separators=(",", ":")))


def _chart_graph_guard(c5: Mapping[str, Any], c6: Mapping[str, Any],
                       q4: Mapping[str, Any]) -> dict[str, Any]:
    c5_nodes, c5b = _nodes(c5, "prescribed_loop_C5_bindings", 121834)
    c6_nodes, c6b = _nodes(c6, "prescribed_R100_Rm_C6_bindings", 97295)
    q4_nodes, q4b = _nodes(q4, "modified_Q4_pressure_C3_bindings", 183700)
    _require(c5_nodes[:len(c6_nodes)] == c6_nodes, "C5 graph is not the exact C6 prefix")
    _require(q4_nodes[:len(c5_nodes)] == c5_nodes, "Q4 graph is not the exact C5 prefix")
    c5charts = c5b.get("charts"); c6charts = c6b.get("charts"); q4charts = q4b.get("charts")
    _require(isinstance(c5charts, Mapping) and isinstance(c6charts, Mapping), "C5/C6 chart bindings missing")
    _require(isinstance(q4charts, Mapping), "Q4 chart bindings missing")
    chart_map = {"post_power": "post_switch", "long_reshape": "long_reshape", "reference": "reference"}
    chart_receipts: dict[str, Any] = {}
    for c5name, c6name in chart_map.items():
        left = c5charts.get(c5name); right = c6charts.get(c6name)
        _require(isinstance(left, Mapping) and isinstance(right, Mapping), f"missing chart: {c5name}")
        ids = {key: left.get(key) for key in ("coordinate", "radius", "physical_measure", "global_offset", "actual_phase")}
        _require(left["coordinate"] == right["phase"] and left["radius"] == right["radius"]
                 and left["physical_measure"] == right["physical_measure"],
                 f"C5/C6 chart identity changed: {c5name}")
        q4row = q4charts.get(c5name)
        _require(isinstance(q4row, Mapping), f"Q4 chart missing: {c5name}")
        chart_receipts[c5name] = {
            "C5": ids,
            "C6": {key: right.get(key) for key in ("phase", "radius", "physical_measure")},
            "Q4_pressure_ports": {
                "radial": _copy_json(q4row.get("radial")),
                "physical": _copy_json(q4row.get("physical")),
                "pressure": _copy_json(q4row.get("pressure")),
            },
            "source_chart_binding": {"C5": c5name, "C6": c6name, "Q4": c5name},
        }
    return {"C5_nodes": len(c5_nodes), "C6_nodes": len(c6_nodes), "Q4_nodes": len(q4_nodes),
            "C5_prefix_equals_C6": True, "Q4_prefix_equals_C5": True,
            "charts": chart_receipts, "passed": True, "source_executed": False}


def _phase_graph_receipt(nodes: Sequence[Mapping[str, Any]], bindings: Mapping[str, Any]) -> dict[str, Any]:
    charts = bindings["charts"]
    N = bindings["N"]
    phase_origin = bindings["phase_origin"]
    _require(isinstance(N, int) and nodes[N].get("operation") == "exact_positive_integer_expression",
             "common C5 N graph handle is not the issued exact integer expression")
    _require(nodes[N].get("one_common_N_for_all_charts_and_control_denominators") is True,
             "C5 N is not the common source integer")
    offset_rows: dict[str, Any] = {}
    radius_denominator: int | None = None
    negative_origin: int | None = None
    for name in SEAM_CHARTS:
        chart = charts[name]
        offset = chart["global_offset"]; phase = chart["actual_phase"]; radius = chart["radius"]
        offrow = nodes[offset]
        _require(offrow.get("operation") == "sum" and len(offrow.get("arguments", ())) == 2,
                 f"global offset is not a two-term source sum: {name}")
        log_id, neg_id = offrow["arguments"]
        logrow = nodes[log_id]; negrow = nodes[neg_id]
        _require(logrow.get("operation") == "analytic_unary" and logrow.get("name") == "log",
                 f"global log offset missing: {name}")
        quotient = nodes[logrow["argument"]]
        _require(quotient.get("operation") == "positive_quotient" and quotient.get("numerator") == radius,
                 f"global offset radius binding changed: {name}")
        denom = quotient.get("denominator")
        _require(isinstance(denom, int), f"global offset denominator missing: {name}")
        if radius_denominator is None:
            radius_denominator = denom
        _require(denom == radius_denominator, "global offset does not retain one physical Ra")
        _require(negrow.get("operation") == "negative", f"global phase origin sign missing: {name}")
        origin_id = negrow.get("argument")
        if negative_origin is None:
            negative_origin = origin_id
        _require(origin_id == negative_origin and origin_id == phase_origin,
                 "global phase origin is not shared across source charts")
        phaserow = nodes[phase]
        _require(phaserow.get("operation") == "analytic_unary" and phaserow.get("name") == "fractional_part",
                 f"fractional phase binding missing: {name}")
        product = nodes[phaserow["argument"]]
        _require(product.get("operation") == "product" and product.get("arguments") == [N, offset],
                 f"common N fractional phase binding changed: {name}")
        offset_rows[name] = {"global_offset": offset, "radius": radius, "actual_phase": phase,
                             "radius_denominator_Ra": denom, "phase_origin_negative": neg_id,
                             "common_N": N}
    return {"common_N": N, "phase_origin": phase_origin,
            "physical_Ra_denominator": radius_denominator,
            "charts": offset_rows, "passed": True, "source_executed": False}


def _substitution(nodes: Sequence[Mapping[str, Any]], index: int, expression: int,
                  variable: int, value: int) -> bool:
    row = nodes[index]
    return (row.get("operation") == "function_substitution" and row.get("expression") == expression
            and row.get("variable") == variable and row.get("value") == value
            and row.get("Z_independent_substitution") is True)


def _history_receipt(c5: Mapping[str, Any], c6: Mapping[str, Any]) -> dict[str, Any]:
    _, c5b = _nodes(c5, "prescribed_loop_C5_bindings", 121834)
    c6nodes, c6b = _nodes(c6, "prescribed_R100_Rm_C6_bindings", 97295)
    c5charts = c5b["charts"]; c5ends = c5b["endpoints"]
    c6charts = c6b["charts"]; c6ends = c6b["endpoints"]
    chart_ports: dict[str, Any] = {}
    for chart_name in SEAM_CHARTS:
        chart = c5charts[chart_name]
        _require(isinstance(chart.get("incoming"), Mapping)
                 and isinstance(chart.get("history"), Mapping)
                 and isinstance(chart.get("memory"), Mapping),
                 f"C5 correction ports missing: {chart_name}")
        for rate in CORRECTION_RATES:
            incoming = chart["incoming"].get(rate); history = chart["history"].get(rate)
            memory = chart["memory"].get(rate)
            _require(isinstance(incoming, list) and len(incoming) == 6
                     and all(isinstance(index, int) for index in incoming),
                     f"C5 incoming correction pointers missing: {chart_name}:{rate}")
            _require(isinstance(history, list) and len(history) == 6
                     and all(isinstance(index, int) for index in history),
                     f"C5 history correction pointers missing: {chart_name}:{rate}")
            _require(isinstance(memory, int), f"C5 memory pointer missing: {chart_name}:{rate}")
        chart_ports[chart_name] = {
            "incoming": {rate: list(chart["incoming"][rate]) for rate in CORRECTION_RATES},
            "history": {rate: list(chart["history"][rate]) for rate in CORRECTION_RATES},
            "memory": {rate: chart["memory"][rate] for rate in CORRECTION_RATES},
        }
    correction: dict[str, Any] = {}
    for left, right in (("post_power", "long_reshape"), ("long_reshape", "reference")):
        rate_rows: dict[str, Any] = {}
        for rate in CORRECTION_RATES:
            endpoint = c5ends[left][rate]; incoming = c5charts[right]["incoming"][rate]
            _require(endpoint == incoming and len(endpoint) == 6,
                     f"C5 correction history pointer changed: {left}->{right}:{rate}")
            rate_rows[rate] = {"endpoint": list(endpoint), "incoming": list(incoming),
                               "equal": True, "rows": 6}
        correction[left + "_to_" + right] = {
            "endpoint_history": rate_rows,
            "incoming_history": {rate: list(c5charts[left]["incoming"][rate]) for rate in CORRECTION_RATES},
            "history": {rate: list(c5charts[left]["history"][rate]) for rate in CORRECTION_RATES},
            "memory": {rate: c5charts[left]["memory"][rate] for rate in CORRECTION_RATES},
            "all_five_rates_retained": True,
        }
    native: dict[str, Any] = {}
    for endpoint, right, label in (("R110", "long_reshape", "post_power_upper_to_long_lower"),
                                   ("Rsh", "reference", "long_upper_to_reference_lower")):
        rows: dict[str, Any] = {}
        for moment in NATIVE_MOMENTS:
            endpoint_rows = c6ends[endpoint]["moments"][moment]
            incoming_rows = c6charts[right]["incoming_histories"][moment]
            _require(endpoint_rows == incoming_rows and len(endpoint_rows) == 7,
                     f"C6 native history pointer changed: {label}:{moment}")
            rows[moment] = {"endpoint": list(endpoint_rows), "incoming": list(incoming_rows),
                            "equal": True, "rows": 7}
        native[label] = {"endpoint": endpoint, "incoming_chart": right, "moments": rows,
                         "all_six_native_moments_retained": True}
    native_chart_ports: dict[str, Any] = {}
    for chart_name in ("post_switch", "long_reshape", "reference"):
        chart = c6charts[chart_name]
        _require(isinstance(chart.get("incoming_histories"), Mapping)
                 and isinstance(chart.get("moments"), Mapping),
                 f"C6 native ports missing: {chart_name}")
        for moment in NATIVE_MOMENTS:
            incoming = chart["incoming_histories"].get(moment); moments = chart["moments"].get(moment)
            _require(isinstance(incoming, list) and len(incoming) == 7
                     and all(isinstance(index, int) for index in incoming),
                     f"C6 incoming native pointers missing: {chart_name}:{moment}")
            _require(isinstance(moments, list) and len(moments) == 7
                     and all(isinstance(index, int) for index in moments),
                     f"C6 native moment pointers missing: {chart_name}:{moment}")
        native_chart_ports[chart_name] = {
            "incoming_histories": {moment: list(chart["incoming_histories"][moment]) for moment in NATIVE_MOMENTS},
            "moments": {moment: list(chart["moments"][moment]) for moment in NATIVE_MOMENTS},
        }
    return {"correction_chart_ports": chart_ports, "correction": correction,
            "native_chart_ports": native_chart_ports, "native": native,
            "all_five_correction_incoming_history_memory_ports": True,
            "all_six_native_incoming_endpoint_moments": True,
            "passed": True, "source_executed": False}


def _sympy_endpoint_check() -> dict[str, Any]:
    """Reduce both C6 endpoint radius formulas independently of graph IDs."""

    try:
        import sympy as sp
    except ImportError:
        return {"available": False, "passed": False,
                "unproved": "SymPy is required for the exact endpoint reduction"}
    hb, T, logref = sp.symbols("hb T logref", real=True)
    post_upper = 100 * sp.exp(2 * hb + (sp.log(sp.Rational(11, 10)) - 2 * hb) * 1)
    long_lower = 110 * sp.exp(T * 0)
    long_upper = 110 * sp.exp(T * 1)
    reference_lower = 110 * sp.exp(T + (logref - T - 8) * 0)
    first_difference = sp.simplify(post_upper - long_lower)
    second_difference = sp.simplify(long_upper - reference_lower)
    _require(first_difference == 0 and second_difference == 0,
             "SymPy endpoint radius reductions are not exact")
    return {
        "available": True,
        "passed": True,
        "post_power_upper": str(sp.expand_log(post_upper, force=True)),
        "long_reshape_lower": str(sp.expand_log(long_lower, force=True)),
        "post_to_long_difference": str(first_difference),
        "long_reshape_upper": str(sp.expand_log(long_upper, force=True)),
        "reference_lower": str(sp.expand_log(reference_lower, force=True)),
        "long_to_reference_difference": str(second_difference),
        "source_formulas": {
            "post_power_upper": "100*exp(2*hb+(log(11/10)-2*hb)*1)",
            "long_reshape_lower": "110*exp(T*0)",
            "long_reshape_upper": "110*exp(T*1)",
            "reference_lower": "110*exp(T+(logref-T-8)*0)",
        },
        "source_tool": "SymPy exact symbols; no archived source execution",
    }


def _endpoint_radius_receipts(c5: Mapping[str, Any], c6: Mapping[str, Any],
                              phase_graph: Mapping[str, Any]) -> dict[str, Any]:
    c5nodes, c5b = _nodes(c5, "prescribed_loop_C5_bindings", 121834)
    c6nodes, c6b = _nodes(c6, "prescribed_R100_Rm_C6_bindings", 97295)
    c5charts = c5b["charts"]; c6charts = c6b["charts"]; c6ends = c6b["endpoints"]
    post = c6charts["post_switch"]; long = c6charts["long_reshape"]; ref = c6charts["reference"]
    r110 = c6ends["R110"]; rsh = c6ends["Rsh"]
    _require(_substitution(c6nodes, r110["radius"], post["radius"], post["phase"], 1),
             "C6 R110 endpoint substitution is not bound to post_switch")
    _require(_substitution(c6nodes, rsh["radius"], long["radius"], long["phase"], 1),
             "C6 Rsh endpoint substitution is not bound to long_reshape")
    # The actual graph rows are retained in the receipt, while equality is
    # proved from the source radius expressions below rather than by treating
    # differing graph nodes as numerically equal.
    phase_rows = phase_graph["charts"]
    common_origin = all(phase_rows[name]["phase_origin_negative"] == phase_rows["post_power"]["phase_origin_negative"]
                        for name in SEAM_CHARTS)
    _require(common_origin, "seam charts do not share the actual global phase origin")
    symbolic_endpoint = _sympy_endpoint_check()
    _require(symbolic_endpoint["passed"], "independent exact endpoint radius reduction failed")
    seams = (
        {
            "name": "post_power_upper_to_long_reshape_lower",
            "left_chart": "post_power", "left_fraction": 1,
            "right_chart": "long_reshape", "right_fraction": 0,
            "left_c6_endpoint": {"name": "R110", "radius": r110["radius"],
                                  "expression": post["radius"], "coordinate": post["phase"]},
            "right_c6_endpoint": {"name": "long_reshape_lower", "radius": long["radius"],
                                   "coordinate": long["phase"]},
            "source_radius_identity": (
                "100*exp(2*hb+(log(11/10)-2*hb)*1)=110=110*exp(T*0)"
            ),
            "physical_y_log_offset": "log(R/Ra)-phase_origin",
            "direct_endpoint_radius_equality": True,
            "integer_offset_difference_assumed": False,
            "global_offset_equality": True,
            "actual_phase_equality_under_common_N": True,
        },
        {
            "name": "long_reshape_upper_to_reference_lower",
            "left_chart": "long_reshape", "left_fraction": 1,
            "right_chart": "reference", "right_fraction": 0,
            "left_c6_endpoint": {"name": "Rsh", "radius": rsh["radius"],
                                  "expression": long["radius"], "coordinate": long["phase"]},
            "right_c6_endpoint": {"name": "reference_lower", "radius": ref["radius"],
                                   "coordinate": ref["phase"]},
            "source_radius_identity": (
                "110*exp(T*1)=110*exp(T+(logref-T-8)*0)"
            ),
            "physical_y_log_offset": "log(R/Ra)-phase_origin",
            "direct_endpoint_radius_equality": True,
            "integer_offset_difference_assumed": False,
            "global_offset_equality": True,
            "actual_phase_equality_under_common_N": True,
        },
    )
    _require(len(seams) == 2 and all(row["direct_endpoint_radius_equality"] for row in seams),
             "endpoint seam receipt incomplete")
    return {"seams": seams,
            "symbolic_endpoint_check": symbolic_endpoint,
            "same_global_phase_origin": common_origin,
            "physical_y_log_offset_retained": True,
            "global_offset_formula": "L(R)=log(R/Ra)-phase_origin",
            "fractional_phase_formula": "fractional_part(N*L(R))",
            "endpoint_phase_equality_reason": "direct equality of physical endpoint radii and shared phase_origin/Ra; no integer offset assumption",
            "integer_offset_difference_assumed": False,
            "passed": True, "source_executed": False}


def _root_and_loop_receipts(c5: Mapping[str, Any], c6: Mapping[str, Any]) -> dict[str, Any]:
    nodes, bindings = _nodes(c5, "prescribed_loop_C5_bindings", 121834)
    _, c6bindings = _nodes(c6, "prescribed_R100_Rm_C6_bindings", 97295)
    result: dict[str, Any] = {}
    c6_source_chart = {"post_power": "post_switch", "long_reshape": "long_reshape", "reference": "reference"}

    def integral_row(index: int, angle_id: int, label: str) -> Mapping[str, Any]:
        row = nodes[index]
        _require(row.get("operation") == "definite_integral", f"{label} is not an issued exact integral")
        _require(row.get("variable") == nodes[angle_id].get("name")
                 and row.get("lower") == 0 and row.get("upper") == angle_id,
                 f"{label} integral endpoints changed")
        _require(row.get("exact_function_integral") is True
                 and row.get("true_loop_C5_original_integrand") is True
                 and row.get("numerical_value_not_installed") is True
                 and row.get("angle_endpoints_independent_of_Z") is True,
                 f"{label} integral is not source-bound")
        return row

    def find_integral(index: int, angle_id: int, label: str, seen: set[int] | None = None) -> tuple[int, Mapping[str, Any]]:
        if seen is None:
            seen = set()
        _require(index not in seen, f"cyclic {label} graph")
        seen.add(index)
        row = nodes[index]
        if row.get("operation") == "definite_integral":
            return index, integral_row(index, angle_id, label)
        children: list[int] = []
        for key in ("arguments", "argument", "body", "expression"):
            value = row.get(key)
            if isinstance(value, int):
                children.append(value)
            elif isinstance(value, list):
                children.extend(value_item for value_item in value if isinstance(value_item, int))
        found = [find_integral(child, angle_id, label, seen) for child in children]
        _require(len(found) == 1, f"{label} must bind one exact integral")
        return found[0]

    for name in SEAM_CHARTS:
        chart = bindings["charts"][name]
        roots = chart["roots"]; loop = chart["loop"]
        for key in ("E", "a", "b", "p2", "t0", "Delta", "kappa"):
            _require(isinstance(roots.get(key), list) and len(roots[key]) == 6,
                     f"actual C5 root jet missing: {name}:{key}")
        for key in ("q", "inverse", "A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi",
                    "Phi", "phase_lambda_at_inverse", "T1_at_inverse", "direction_at_inverse"):
            _require(isinstance(loop.get(key), list) and len(loop[key]) == 6,
                     f"actual C5 loop jet missing: {name}:{key}")
        _require(isinstance(loop.get("angle"), int), f"actual C5 angle handle missing: {name}")
        angle_id = loop["angle"]
        for key in ("b", "t0"):
            _require(all(value == 0 for value in roots[key]),
                     f"actual C5 {key} root is not the issued exact zero specialization: {name}")
            _require(all(nodes[value].get("operation") == "exact_rational"
                         and nodes[value].get("numerator") == 0
                         and nodes[value].get("denominator") == 1 for value in roots[key]),
                     f"actual C5 {key} graph rows are not exact zero: {name}")
        for order, a_id in enumerate(roots["a"]):
            a_row = nodes[a_id]
            _require(a_row.get("operation") == "positive_quotient"
                     and a_row.get("denominator") == roots["E"][0]
                     and a_row.get("source_positive_certificate") == "same strictly positive actual prescribed swirl E",
                     f"actual C5 a=-? root positivity/division binding changed: {name}:{order}")
        c6chart = c6bindings["charts"][c6_source_chart[name]]
        c6_vy = c6chart["radial"]["V_y"]
        _require(isinstance(c6_vy, list) and len(c6_vy) == 7 and all(value == 0 for value in c6_vy),
                 f"native C6 V_y=0 source binding changed: {name}")
        t1_substitution = nodes[loop["T1_at_inverse"][0]]
        _require(t1_substitution.get("operation") == "substitute_original_inverse_angle",
                 f"actual C5 T1 inverse substitution missing: {name}")
        t1_id, t1_integral = find_integral(t1_substitution["body"], angle_id, f"{name}.T1")
        phi_row = nodes[loop["Phi"][0]]
        _require(phi_row.get("operation") == "product", f"actual C5 Phi product missing: {name}")
        p_id = next((value for value in phi_row["arguments"]
                     if isinstance(value, int) and nodes[value].get("operation") == "sum"
                     and angle_id in nodes[value].get("arguments", [])), None)
        _require(isinstance(p_id, int), f"actual C5 Phi=K*(angle+T2) binding missing: {name}")
        t2_id, t2_integral = find_integral(
            next(value for value in nodes[p_id]["arguments"] if value != angle_id),
            angle_id, f"{name}.T2")
        t1_integrand = t1_integral["integrand"]
        t2_integrand_row = nodes[t2_integral["integrand"]]
        _require(t2_integrand_row.get("operation") == "product"
                 and t2_integrand_row.get("arguments") == [t1_integrand, t1_integrand],
                 f"actual C5 T2 is not the square of actual T1 integrand: {name}")
        result[name] = {
            "roots": {key: list(roots[key]) for key in ("E", "a", "b", "p2", "t0", "Delta", "kappa")},
            "loop": {key: list(loop[key]) if isinstance(loop[key], list) else loop[key]
                     for key in ("q", "inverse", "A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi",
                                 "Phi", "phase_lambda_at_inverse", "T1_at_inverse", "direction_at_inverse", "angle")},
            "native_C6_V_y_graph_rows_are_exact_zero": list(c6_vy),
            "b_and_t0_graph_rows_are_exact_zero": True,
            "t0_source_formula": "t0=-b/a; graph specialization b=0 and issued a=positive_quotient(C,E)>0, hence t0=0",
            "T1_integral_graph": {"node": t1_id, "integrand": t1_integrand,
                                  "variable": t1_integral["variable"], "lower": t1_integral["lower"],
                                  "upper": t1_integral["upper"]},
            "T2_integral_graph": {"node": t2_id, "integrand": t2_integral["integrand"],
                                  "squared_T1_integrand": [t1_integrand, t1_integrand],
                                  "variable": t2_integral["variable"], "lower": t2_integral["lower"],
                                  "upper": t2_integral["upper"]},
            "actual_T1_T2_integral_definitions_bound": True,
        }
    return {"charts": result, "formula_handles_are_actual_C5_graph_nodes": True,
            "passed": True, "source_executed": False}


def _sympy_periodic_check() -> dict[str, Any]:
    """Check the rational kernel identities with exact symbolic arithmetic."""

    try:
        import sympy as sp
    except ImportError:
        return {"available": False, "passed": False,
                "unproved": "SymPy is required for the new exact rational check"}
    r, c, s = sp.symbols("r c s", real=True)
    z = c + sp.I * s
    zbar = c - sp.I * s
    denominator = 1 - 2 * r * c + r ** 2
    complex_real = sp.cancel((z / (1 - r * z) + zbar / (1 - r * zbar)) / 2)
    difference = sp.together(complex_real - (c - r) / denominator)
    numerator, _ = sp.fraction(difference)
    unit_circle_remainder = sp.groebner([c ** 2 + s ** 2 - 1], c, s, r, order="lex").reduce(numerator)[1]
    quotient = sp.factor(sp.cancel(numerator / (c ** 2 + s ** 2 - 1)))
    q, u, t0, a, b = sp.symbols("q u t0 a b", real=True)
    rho = u / sp.sqrt(1 + u ** 2)
    alpha = 2 * q / sp.sqrt(1 + u ** 2)
    alpha_identity = sp.simplify(alpha ** 2 / (1 - rho ** 2))
    K = 1 / (2 * sp.pi * (1 + t0 ** 2 + 2 * q ** 2))
    K_identity = sp.simplify(K * 2 * sp.pi * (1 + t0 ** 2 + 2 * q ** 2))
    M_identity = sp.simplify((-a * (2 * sp.pi * t0) / (2 * sp.pi) - b).subs(t0, -b / a))
    n = sp.symbols("n", integer=True)
    norm_raw = sp.summation(r ** (2 * (n - 1)) / 2, (n, 1, sp.oo))
    norm = sp.refine(norm_raw, sp.Q.positive(1 - r ** 2))
    passed = (unit_circle_remainder == 0 and sp.simplify(alpha_identity - 4 * q ** 2) == 0
              and sp.simplify(K_identity - 1) == 0 and sp.simplify(M_identity) == 0
              and sp.simplify(norm - 1 / (2 * (1 - r ** 2))) == 0)
    _require(passed, "exact SymPy C5 Fourier/rational identities failed")
    return {
        "available": True,
        "passed": True,
        "complex_real_part": "Re[z/(1-r*z)] with z=c+i*s",
        "rational_identity": "Re[z/(1-r*z)]-(c-r)/(1-2*r*c+r^2)",
        "unit_circle_relation": "c^2+s^2-1=0",
        "unit_circle_remainder": str(unit_circle_remainder),
        "rational_numerator_factor": str(quotient) + "*(c^2+s^2-1)",
        "alpha_squared_over_one_minus_r_squared": str(alpha_identity),
        "K_drift": str(K_identity),
        "M_drift_after_t0_minus_b_over_a": str(M_identity),
        "geometric_norm": "1/(2*(1-r**2))",
        "geometric_norm_sympy_normal_form": str(norm),
        "source_tool": "SymPy exact symbols; no archived source execution",
    }


def _periodic_formula_receipt() -> dict[str, Any]:
    symbolic = _sympy_periodic_check()
    return {
        "formula": {
            "t": "t0 + alpha*(cos(psi)-r)/(1-2*r*cos(psi)+r^2)",
            "r": "u/sqrt(1+u^2)",
            "alpha": "2*q/sqrt(1+u^2)",
            "K": "1/(2*pi*(1+t0^2+2*q^2))",
            "Phi": "K*(psi+T2)",
            "T1": "integral(t dpsi)",
            "T2": "integral(t^2 dpsi)",
            "chi": "phase-psi/(2*pi)",
            "M": "-a*T1hat/(2*pi)-b*phase",
            "A": "a*chi/2",
            "B_over_Pstar": "E*M/2",
            "t0": "-b/a",
        },
        "fourier_means": {
            "mean_(cos-r)/D": "0",
            "mean_((cos-r)/D)^2": "1/(2*(1-r^2))=(1+u^2)/2",
            "mean_t": "t0",
            "mean_t_squared": "t0^2+2*q^2",
            "mean_t_squared_expansion": "t0^2+2*t0*alpha*0+alpha^2/(2*(1-r^2))=t0^2+2*q^2",
        },
        "drifts": {
            "T1(psi+2*pi)-T1(psi)": "2*pi*t0",
            "T2(psi+2*pi)-T2(psi)": "2*pi*(t0^2+2*q^2)",
            "Phi(psi+2*pi)-Phi(psi)": "1",
            "psi(phase+1)-psi(phase)": "2*pi",
            "chi(phase+1)-chi(phase)": "0",
            "M(phase+1)-M(phase)": "-a*t0-b=0",
        },
        "geometric_complex_series": {
            "identity": "e^(i*psi)/(1-r*e^(i*psi))=sum(n>=1,r^(n-1)*e^(i*n*psi))",
            "real_part": "(cos(psi)-r)/(1-2*r*cos(psi)+r^2)",
            "hypothesis": "|r|<1",
            "orthogonality_mean_cos_n": "0 for n>=1",
            "orthogonality_mean_cos_n_cos_m": "delta_nm/2 for n,m>=1",
            "geometric_norm": "sum(n>=1,r^(2*n-2))/2=1/(2*(1-r^2))",
        },
        "sympy_exact_check": symbolic,
        "A_periodic_under_integer_phase_shift": True,
        "B_over_Pstar_periodic_under_integer_phase_shift": True,
        "A_periodic_under_integer_phase_shift_formal": True,
        "B_over_Pstar_periodic_under_integer_phase_shift_formal": True,
        "actual_inverse_periodicity_status": "conditional_on_source_positive_phase_controls",
        "actual_A_B_periodicity_status": "conditional_on_source_positive_phase_controls",
        "positive_phase_controls_proved_here": False,
        "active_cutoff_controls_proved_here": False,
        "source_derivative_positive_certificate": "K*(1+t^2)>0; same original a/q-flat construction",
        "positive_control_dependency": "eta,d_star,a,q and active flat-cutoff source controls are inherited, not reproved here",
        "positive_phase_lambda": "K*(1+t^2)>0",
        "implicit_function_theorem_smoothness": "conditional_on_source_positive_phase_controls",
        "implicit_function_theorem_smoothness_proved_here": False,
        "proof_kind": "exact Fourier/source formula; no numeric source evaluation",
        "passed": symbolic.get("passed") is True,
        "source_executed": False,
    }


def _q4_derivative_receipt(q4: Mapping[str, Any]) -> dict[str, Any]:
    _, bindings = _nodes(q4, "modified_Q4_pressure_C3_bindings", 183700)
    charts = bindings["charts"]
    ports: dict[str, Any] = {}
    for name in SEAM_CHARTS:
        row = charts[name]
        _require(isinstance(row.get("radial"), Mapping) and isinstance(row.get("physical"), Mapping),
                 f"Q4 native derivative ports missing: {name}")
        ports[name] = {"radial": _copy_json(row["radial"]), "physical": _copy_json(row["physical"]),
                       "pressure": _copy_json(row.get("pressure")),
                       "total_fast_phase_terms_included": row["radial"].get("total_fast_phase_terms_included") is True,
                       "physical_y_native_partial_divided_by_Jacobian": row["radial"].get(
                           "physical_y_equals_native_partial_divided_by_Jacobian") is True}
    return {
        "source_guard_bound": True,
        "actual_native_C2_source_ports": ports,
        "native_C2_source_guard_passed": True,
        "actual_derivative_identity_proved": False,
        "modified_Uy_Qy_closed": False,
        "natural_B_y_closed": False,
        "unproved": (
            "native physical-y derivative equality is source-bound but its full endpoint identity is open",
            "modified U_y/Q_y pressure closure is not proved by loop periodicity",
        ),
        "passed": True,
        "source_executed": False,
    }


def build(root: str | Path, parent: Mapping[str, Any] | str | Path | None = None) -> dict[str, Any]:
    """Build the source-only periodic seam receipt from a pinned parent map."""

    root_path = Path(root)
    parent_map = _load_parent(root_path, parent)
    pins = _source_pins(root_path)
    source_guards = _source_guards(root_path)
    c5 = _gzip_json(root_path / C5_CACHE)
    c6 = _gzip_json(root_path / C6_CACHE)
    q4 = _gzip_json(root_path / Q4_CACHE)
    graph = _chart_graph_guard(c5, c6, q4)
    _, c5bindings = _nodes(c5, "prescribed_loop_C5_bindings", 121834)
    phase_graph = _phase_graph_receipt(c5["exact_graph_nodes"], c5bindings)
    history = _history_receipt(c5, c6)
    endpoints = _endpoint_radius_receipts(c5, c6, phase_graph)
    roots = _root_and_loop_receipts(c5, c6)
    periodic = _periodic_formula_receipt()
    q4_derivative = _q4_derivative_receipt(q4)
    module_sha256 = _sha(root_path / MODULE)
    source_family = parent_map.get("source_family")
    checks = {
        "source_pins": True,
        "source_guards": True,
        "graph_prefix": graph["passed"],
        "actual_C5_loop_formula": periodic["passed"],
        "actual_C5_root_and_loop_handles": roots["passed"],
        "actual_C5_t0_equals_minus_b_over_a_graph_specialization": all(
            row["b_and_t0_graph_rows_are_exact_zero"]
            and row["actual_T1_T2_integral_definitions_bound"]
            for row in roots["charts"].values()),
        "global_phase_binding": phase_graph["passed"],
        "endpoint_source_radius_equalities": endpoints["passed"],
        "endpoint_radius_symbolic_reduction": endpoints["symbolic_endpoint_check"]["passed"],
        "native_history_binding": history["passed"],
        "correction_history_binding": history["passed"],
        "native_C2_source_guard": q4_derivative["native_C2_source_guard_passed"],
    }
    return {
        "passed": all(checks.values()),
        "source_only": True,
        "source_executed": False,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "module_sha256": module_sha256,
        "source_family": _copy_json(source_family),
        "source_pins": pins,
        "source_guards": source_guards,
        "sourceguards_tuple": True,
        "AST": {
            "C5_loop_formula": True,
            "C6_endpoint_and_chart_formula": True,
            "Q4_native_derivative_and_pressure_ports": True,
            "strict_lf_single_eof": True,
        },
        "checks": checks,
        "graph": graph,
        "phase_graph": phase_graph,
        "periodic_formula": periodic,
        "periodic_phase_lemma": periodic,
        "actual_C5_root_and_loop_handles": roots,
        "endpoint_identity_receipts": endpoints,
        "same_global_phase_endpoint_receipt": endpoints,
        "history_binding": history,
        "actual_C5_endpoint_history_binding": history,
        "native_C2_and_actual_derivative": q4_derivative,
        "native_C2": {
            "source_guard_passed": q4_derivative["native_C2_source_guard_passed"],
            "source_ports_bound": True,
            "endpoint_identity_proved": False,
        },
        "actual_derivative": {
            "source_guard_passed": q4_derivative["source_guard_bound"],
            "identity_proved": q4_derivative["actual_derivative_identity_proved"],
            "unproved": list(q4_derivative["unproved"]),
        },
        "unproved": list(q4_derivative["unproved"]) + [
            "actual inverse/A/B integer-shift periodicity is conditional on inherited positive phase and active-cutoff controls",
        ],
        "limits": {
            "modified_Uy_Qy_closed": False,
            "natural_B_y_closed": False,
            "full_pressure_global_closure": False,
            "global_five_moment_closure": False,
            "source_values_selected": False,
            "interval_or_cap_values_used": False,
            "G0_chosen": False,
            "N_materialized": False,
            "accepted_builder_checker_owner_producer_norm_recurrence_quadrature_executed": False,
            "ordinary_Z_orders": [0, 1, 2, 3, 4],
            "physical_y_derivatives_claimed": False,
            "periodic_history_seams_only": True,
            "actual_inverse_periodicity_conditional_on_inherited_controls": True,
        },
        "scope": {
            "seams": [row["name"] for row in endpoints["seams"]],
            "same_global_phase": True,
            "formal_C5_periodicity_under_explicit_controls": True,
            "actual_C5_periodicity": "conditional_on_source_positive_phase_controls",
            "actual_C5_history_binding": True,
            "native_C2_and_derivative_handled_separately": True,
            "runtime_model_metadata_exposed": False,
        },
        "source_contract": {
            "source_only": True,
            "common_N_graph_handle": phase_graph["common_N"],
            "N_materialized": False,
            "physical_y_log_offset_retained": endpoints["physical_y_log_offset_retained"],
            "integer_offset_difference_assumed": False,
            "all_five_correction_incoming_history_memory_ports": history[
                "all_five_correction_incoming_history_memory_ports"],
            "all_six_native_incoming_endpoint_moments": history[
                "all_six_native_incoming_endpoint_moments"],
            "modified_Uy_Qy_closed": False,
            "natural_B_y_closed": False,
            "full_pressure_global_closure": False,
        },
        "pure_checks": {
            "source_executed": False,
            "accepted_builder_checker_owner_producer_norm_recurrence_quadrature_executed": False,
            "source_values_selected": False,
            "interval_or_cap_values_used": False,
            "G0_chosen": False,
        },
    }


__all__ = ["build"]
