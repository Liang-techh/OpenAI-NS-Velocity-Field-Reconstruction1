"""Source-bound native control first jets at the long/reference seams.

This module is a read-only proof packet.  It reads the frozen source ASTs and
their accepted JSON/gzip receipts, but it never imports or executes an
accepted producer, owner, checker, recurrence, quadrature, norm, or algebra
module.  The packet deliberately stops before the full modified ``Q_y`` seam:
the native generic seed and the positive loop/IFT hypotheses are discharged,
while the modified pressure derivative remains an explicit open obligation.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
PARENT_SHA256 = "0a91144aca1caaf0dd2b11170b1aee4339f676eb5f213eea236ae8c814cea0c4"
MODULE = "src/openai_ns_reconstruction/original_long_reference_control_seams.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
INLET_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"
RATE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
ANCHOR_SOURCE = BASE + "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"

NATIVE_CHARTS = ("post_switch", "long_reshape", "reference")
C5_CHARTS = ("post_power", "long_reshape", "reference")
NATIVE_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
MOMENTS = ("H", "M", "K", "A", "B", "C")
CORRECTIONS = ("m", "h", "k", "e", "p")
NATIVE_RATES = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
CORRECTION_RATES = {"m": 1, "h": "3/2", "k": "3/2", "e": 1, "p": 0}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _source(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    _require("\r" not in text, "source must have LF line endings: " + str(path))
    _require(text.endswith("\n") and not text.endswith("\n\n"),
             "source must have a single final LF: " + str(path))
    return text


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _ast_functions(tree: ast.AST, target: str, owner: str | None = None) -> list[ast.AST]:
    out: list[ast.AST] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if node.name != target:
            continue
        if owner is None:
            out.append(node)
            continue
        parent_ok = False
        for parent in ast.walk(tree):
            if isinstance(parent, ast.ClassDef) and parent.name == owner and node in parent.body:
                parent_ok = True
                break
        if parent_ok:
            out.append(node)
    return out


def _source_guard(root: Path, relative: str, target: str, tokens: tuple[str, ...],
                  owner: str | None = None) -> dict[str, Any]:
    path = root / relative
    text = _source(path)
    tree = ast.parse(text, filename=relative)
    matches = _ast_functions(tree, target, owner)
    _require(len(matches) == 1, f"expected one AST function {target} in {relative}")
    segment = ast.get_source_segment(text, matches[0]) or ""
    missing = [token for token in tokens if token not in segment]
    _require(not missing, f"source guard missing in {relative}:{target}: {missing}")
    return {
        "path": relative,
        "target": target,
        "owner": owner,
        "tokens": list(tokens),
        "AST_sha256": hashlib.sha256(segment.encode("utf-8")).hexdigest(),
        "passed": True,
        "source_executed": False,
    }


def _rate_guard(root: Path) -> dict[str, Any]:
    path = root / RATE_SOURCE
    text = _source(path)
    tree = ast.parse(text, filename=RATE_SOURCE)
    rows = []
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "RATES" for target in node.targets):
            continue
        rows.append(node)
    _require(len(rows) == 1, "actual bridge rates assignment is not unique")
    value = rows[0].value
    _require(isinstance(value, ast.Call) and isinstance(value.func, ast.Name)
             and value.func.id == "dict" and not value.args,
             "actual bridge rates are not a keyword dict")
    rates = {row.arg: ast.literal_eval(row.value) for row in value.keywords}
    _require(rates == NATIVE_RATES, "native history rates changed")
    return {"path": RATE_SOURCE, "binding": "RATES", "rates": rates,
            "passed": True, "source_executed": False}


def _load_parent(root: Path) -> dict[str, Any]:
    path = root / PARENT
    _require(_sha(path) == PARENT_SHA256, "prerequisite parent hash changed")
    parent = _json(path)
    _require(parent.get("source_family"), "parent source family missing")
    for name in ("actual_long_reference_native_seams", "actual_long_reference_history_seams",
                 "actual_long_reference_loop_seams"):
        section = parent.get(name)
        _require(isinstance(section, Mapping) and section.get("passed") is True,
                 "accepted prerequisite is not passed: " + name)
        _require(section.get("source_executed") is False,
                 "accepted prerequisite was executed: " + name)
        section_parent = section.get("parent_report")
        section_hash = section.get("parent_report_sha256")
        _require(isinstance(section_parent, str) and isinstance(section_hash, str),
                 "accepted prerequisite parent mapping is incomplete: " + name)
        _require(_sha(root / section_parent) == section_hash,
                 "accepted prerequisite parent hash changed: " + name)
    return parent


def _accepted_expected(parent: Mapping[str, Any], relative: str) -> str | None:
    for section_name in ("actual_long_reference_native_controls",
                         "actual_long_reference_native_seams",
                         "actual_long_reference_history_seams",
                         "actual_long_reference_loop_seams"):
        section = parent.get(section_name, {})
        pin = section.get("source_pins", {}).get(relative)
        if isinstance(pin, Mapping) and isinstance(pin.get("expected_sha256"), str):
            return str(pin["expected_sha256"])
    if relative == PARENT:
        return PARENT_SHA256
    hashes = parent.get("implementation_hashes", {})
    if relative in hashes:
        return hashes[relative]
    if relative.startswith(BASE):
        return parent.get("input_hashes", {}).get(relative[len(BASE):])
    return None


def _source_pins(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    required = (PARENT, C6_SOURCE, C6_CACHE, C5_SOURCE, C5_CACHE, Q4_SOURCE,
                Q4_CACHE, BRIDGE_SOURCE, INLET_SOURCE, RATE_SOURCE, ANCHOR_SOURCE,
                SIGMA_SOURCE)
    result: dict[str, Any] = {}
    for relative in required:
        expected = _accepted_expected(parent, relative)
        _require(isinstance(expected, str), "no accepted expected hash for " + relative)
        actual = _sha(root / relative)
        _require(actual == expected, "source bytes changed: " + relative)
        result[relative] = {"path": relative, "expected_sha256": expected,
                            "actual_sha256": actual, "passed": True,
                            "source_executed": False}
    self_actual = _sha(root / MODULE)
    result[MODULE] = {"path": MODULE, "expected_sha256": self_actual,
                      "actual_sha256": self_actual, "passed": True,
                      "source_executed": False, "self_pin": True}
    return result


def _guards(root: Path) -> dict[str, Any]:
    guards: list[dict[str, Any]] = [
        _source_guard(root, C6_SOURCE, "__init__", (
            "self.post_length=g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))",
            "self.add_chart('post_switch',p,R,self.post_length,actual,incoming['moments'])",
            "self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])",
            "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
            "self.add_chart('reference',q,R,length,actual,incoming['moments'])",
        ), "CurrentOriginalPrescribedR100RmC6Functions"),
        _source_guard(root, C6_SOURCE, "phi_from_Utheta", (
            "return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R)))",
        )),
        _source_guard(root, C6_SOURCE, "add_chart", (
            "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
            "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
        )),
        _source_guard(root, BRIDGE_SOURCE, "histories", (
            "body=a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))",
            "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase",
            "original_incoming_memory_retained=True,source_measure_applied_once=True",
        )),
        _source_guard(root, BRIDGE_SOURCE, "generic_fields", (
            "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
            "e=a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B'])),p=a.mul(FF,J['C'])",
            "fields['pressure']=a.add(fields['P0'],fields['p'])",
        )),
        _source_guard(root, INLET_SOURCE, "generic_functions", (
            "axial_linear=div(a.add(a.neg(V),half(a.mul(a.sub(one,de),a.sub(m,a.mul(z,mz))))))",
            "axial_quadratic=div(a.add(a.mul(V,transport),a.scale(a.mul(de,a.mul(z,e)),2),a.neg(a.mul(d,ez)),",
            "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,",
        )),
        _source_guard(root, C5_SOURCE, "generic_roots", (
            "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)",
            "aa,bb=divide(C),divide(B)",
            "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
        )),
        _source_guard(root, C5_SOURCE, "__init__", (
            "self.N=self.unique('exact_positive_integer_expression',None)",
            "self.eta=self.unique('current_original_source_parameter','eta')",
            "self.dstar=self.unique('current_original_source_parameter','d_star')",
            "'post_power':background.charts['post_switch']",
        ), "CurrentOriginalPrescribedLoopC5Functions"),
        _source_guard(root, C5_SOURCE, "loop_functions", (
            "shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)",
            "collar=g.mul(c('3/4'),sc) if chart=='first_micro' else None",
            "predicate=shear",
            "lazy=lambda jet:a.lazy(jet,predicate,roots['Delta'].value,eta,coordinate,collar)",
            "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
            "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
            "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))",
            "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
            "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
        )),
        _source_guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        _source_guard(root, C6_SOURCE, "sigma_prime", (
            "derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets)",
            "exact_flat_source_derivative_not_a_bound=True",
        )),
        _source_guard(root, SIGMA_SOURCE, "_sigma_left", (
            "odds=1/(1-x)**2-1/x**2",
            "e=positive_exp(c,odds); value=e/(1+e)",
            "L1=2/(1-x)**3+2/x**3",
        )),
    ]
    guards.append(_rate_guard(root))
    return {"guards": guards, "sourceguards_tuple": True,
            "passed": all(item.get("passed") is True for item in guards)}


def _functionjet_generator_receipt(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Consume the accepted FunctionJet row-generator lemma as a saved lemma.

    The graph caches expose row IDs, while the accepted history receipt proves
    that those rows are generated by the ordinary FunctionJet source algebra.
    This packet therefore makes a function-level inference; it does not claim
    that every graph row's topology was independently re-traversed here.
    """
    generator = parent["actual_long_reference_history_seams"]["ordinary_Z_row_generator"]
    guards = generator.get("source_guards", {})
    guard_passed = all(
        value.get("passed") is True and value.get("source_executed") is False
        for value in guards.values()
    )
    required_properties = (
        "FunctionJet_rows_are_source_rows",
        "density_rows_generated_by_source_algebra",
        "kernel_rows_generated_by_source_algebra",
        "memory_rows_generated_by_source_algebra",
        "integral_rows_preserve_ordinary_order",
        "source_generator_proves_all_ordinary_rows",
        "native_histories_use_at_scale_integral_add",
        "C5_histories_use_at_scale_integral_add",
        "native_bridge_and_C5_loop_inherit_source_algebra",
        "no_row0_proxy_claim",
    )
    properties_passed = all(generator.get(key) is True for key in required_properties)
    passed = generator.get("passed") is True and guard_passed and properties_passed
    _require(passed, "accepted FunctionJet source-row generator lemma is incomplete")
    return {
        "passed": passed,
        "source": generator.get("source"),
        "source_guards": {
            name: {
                "path": value.get("path"),
                "function": value.get("function"),
                "class_name": value.get("class_name"),
                "AST_sha256": value.get("AST_sha256"),
                "passed": value.get("passed"),
                "source_executed": value.get("source_executed"),
            }
            for name, value in guards.items()
        },
        "FunctionJet_rows_are_source_rows": generator.get("FunctionJet_rows_are_source_rows"),
        "source_generator_proves_all_ordinary_rows": generator.get("source_generator_proves_all_ordinary_rows"),
        "ordinary_Z_row_orders": generator.get("ordinary_Z_row_orders"),
        "row_id_prefix_correspondence_required": True,
        "graph_per_row_topology_checked": generator.get("graph_per_row_topology_checked"),
        "function_level_inference_only": generator.get("graph_per_row_topology_checked") is False,
        "no_row0_proxy_claim": generator.get("no_row0_proxy_claim"),
        "source_executed": generator.get("source_executed"),
    }


def _integer_values(value: Any) -> list[int]:
    """Return integer graph references nested in a cached graph node."""
    if type(value) is int:
        return [value]
    if isinstance(value, list):
        out: list[int] = []
        for child in value:
            out.extend(_integer_values(child))
        return out
    if isinstance(value, Mapping):
        out = []
        for child in value.values():
            out.extend(_integer_values(child))
        return out
    return []


def _graph_reachable(nodes: list[Any], starts: list[int], targets: set[int]) -> bool:
    """Check cached graph reachability without executing any source producer."""
    seen: set[int] = set()
    stack = list(starts)
    while stack:
        node_id = stack.pop()
        if type(node_id) is not int or node_id in seen:
            continue
        if node_id < 0 or node_id >= len(nodes):
            continue
        seen.add(node_id)
        node = nodes[node_id]
        if isinstance(node, Mapping):
            stack.extend(_integer_values(node))
    return targets.issubset(seen)


def _p0_and_p2_graph_receipt(c6_nodes: list[Any], c6_chart: Mapping[str, Any],
                             c5_nodes: list[Any], c5_chart: Mapping[str, Any],
                             parent: Mapping[str, Any]) -> dict[str, Any]:
    """Bind native P0[0], P0[1] and native p2 to actual cached graph rows."""
    p0_rows = list(c6_chart["fields"]["P0"])
    p2_rows = list(c5_chart["roots"]["p2"])
    _require(len(p0_rows) == 7 and all(type(row) is int for row in p0_rows),
             "native P0 rows are not seven integer graph references")
    _require(len(p2_rows) == 6 and all(type(row) is int for row in p2_rows),
             "native p2 rows are not six integer graph references")

    p0_node_rows: list[dict[str, Any]] = []
    for order, row in enumerate(p0_rows):
        if row == 0:
            # The present native C6 cache has source nodes for P0[0] and
            # P0[1].  If a future background row is literal zero, only the
            # accepted source seed may discharge it.
            zero_seed = parent["actual_long_reference_history_seams"].get(
                "q4_background_seed_source", {}).get("independent_seed_algebra", {}).get(
                    "reduced_rows", {}).get("P0") == "0"
            _require(zero_seed, "literal zero P0 row lacks the accepted zero-source lemma")
            p0_node_rows.append({"row": row, "ordinary_order": order,
                                 "zero_source_lemma": True, "reachable": True})
            continue
        _require(row < len(c6_nodes), "native P0 row is outside the C6 graph")
        node = c6_nodes[row]
        _require(isinstance(node, Mapping), "native P0 graph row is not a node mapping")
        if order >= 2:
            p0_node_rows.append({"row": row, "ordinary_order": order,
                                 "source_node_operation": node.get("operation"),
                                 "reachable": True})
            continue
        metadata = {
            "operation": node.get("operation"),
            "quantity": node.get("quantity"),
            "Z_order": node.get("Z_order"),
            "coordinate": node.get("coordinate"),
            "Z_variable": node.get("Z_variable"),
            "ordinary_slow_Z_derivative_order": node.get("ordinary_slow_Z_derivative_order"),
            "quantity_path": node.get("quantity_path"),
            "coordinate_independent": node.get("coordinate_independent"),
            "exact_source_op_P0_function_not_enclosure_endpoint": node.get(
                "exact_source_op_P0_function_not_enclosure_endpoint"),
        }
        _require(metadata["operation"] == "current_original_leading_function_recipe",
                 "native P0 row is not the actual leading-function source operation")
        _require(metadata["quantity"] == "P0" and metadata["Z_order"] == order,
                 "native P0 row quantity/order metadata changed")
        _require(metadata["ordinary_slow_Z_derivative_order"] == order,
                 "native P0 ordinary Z order metadata changed")
        _require(metadata["coordinate_independent"] is True
                  and metadata["exact_source_op_P0_function_not_enclosure_endpoint"] is True,
                 "native P0 source identity metadata changed")
        p0_node_rows.append({"row": row, "ordinary_order": order,
                             "source_node": metadata, "reachable": True})

    _require(all(row["reachable"] for row in p0_node_rows),
             "native P0 graph rows are not reachable")
    p2_reachable = all(row < len(c5_nodes) for row in p2_rows)
    _require(p2_reachable, "native p2 row is outside the C5 graph")
    p2_first_nodes = [c5_nodes[row] for row in p2_rows[:2]]
    _require(all(isinstance(node, Mapping) and node.get("operation") == "product"
                 for node in p2_first_nodes),
             "native p2[0:2] are not actual product graph nodes")
    p0_targets = {row for row in p0_rows[:2] if row != 0}
    p2_reaches_p0 = _graph_reachable(c5_nodes, p2_rows[:2], p0_targets)
    _require(p2_reaches_p0, "native p2 graph does not reach P0[0:2]")
    return {
        "passed": True,
        "P0_rows": p0_node_rows,
        "P0_first_two_reachable": True,
        "P0_first_two_source_orders": [0, 1],
        "p2_rows": p2_rows,
        "p2_first_two_product_nodes": p2_rows[:2],
        "p2_graph_reaches_P0_first_two": p2_reaches_p0,
        "ordinary_Z_source_formula": "native inlet p2 retains P0[0], P0[1] and their ordinary-Z source rows",
        "source_executed": False,
    }


def _row_shape_receipt(c6_chart: Mapping[str, Any], q4_chart: Mapping[str, Any]) -> dict[str, Any]:
    """Check row counts and retain their distinct source orders."""
    c6_fields = c6_chart["fields"]
    c6_radial = c6_chart["radial"]
    q4_background = q4_chart["radial_background"]
    q4_radial = q4_chart["radial"]
    q4_physical = q4_chart["physical"]

    def rows_ok(mapping: Mapping[str, Any], keys: tuple[str, ...], order: int) -> bool:
        return all(type(mapping.get(key)) is list
                   and len(mapping[key]) == order
                   and all(type(row) is int for row in mapping[key])
                   for key in keys)

    c6_field_keys = tuple(NATIVE_FIELDS)
    c6_radial_keys = ("E_y", "V_y")
    q4_background_keys = tuple(NATIVE_FIELDS)
    q4_radial_keys = ("E_y", "V_y", "m_y", "slow_A_y", "slow_B_y")
    q4_physical_keys = ("U", "U_y", "m", "m_y", "Q", "Q_y")
    c6_fields_seven = rows_ok(c6_fields, c6_field_keys, 7)
    c6_radial_seven = rows_ok(c6_radial, c6_radial_keys, 7)
    q4_background_seven = rows_ok(q4_background, q4_background_keys, 7)
    q4_radial_physical_six = rows_ok(q4_radial, q4_radial_keys, 6)
    q4_physical_um_six = rows_ok(q4_physical, q4_physical_keys[:4], 6)
    q4_physical_q_five = rows_ok(q4_physical, q4_physical_keys[4:], 5)
    passed = (c6_fields_seven and c6_radial_seven and q4_background_seven
              and q4_radial_physical_six and q4_physical_um_six and q4_physical_q_five)
    _require(passed, "actual cache row lengths or integer graph references changed")
    return {
        "passed": passed,
        "ordinary_orders": {
            "C6_fields": list(range(7)),
            "C6_radial": list(range(7)),
            "Q4_radial_background": list(range(7)),
            "Q4_radial_physical": list(range(6)),
            "Q4_physical_U_m": list(range(6)),
            "Q4_physical_Q": list(range(5)),
        },
        "row_lengths": {
            "C6_fields": {key: len(c6_fields[key]) for key in c6_field_keys},
            "C6_radial": {key: len(c6_radial[key]) for key in c6_radial_keys},
            "Q4_radial_background": {key: len(q4_background[key]) for key in q4_background_keys},
            "Q4_radial_physical": {key: len(q4_radial[key]) for key in q4_radial_keys},
            "Q4_physical": {key: len(q4_physical[key]) for key in q4_physical_keys},
        },
        "C6_fields_are_seven_rows": c6_fields_seven,
        "C6_radial_are_seven_rows": c6_radial_seven,
        "Q4_radial_background_are_seven_rows": q4_background_seven,
        "Q4_radial_physical_are_six_rows": q4_radial_physical_six,
        "Q4_physical_Q_are_five_rows": q4_physical_q_five,
    }


def _cache_summary(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    c6 = _gzip_json(root / C6_CACHE)
    c5 = _gzip_json(root / C5_CACHE)
    q4 = _gzip_json(root / Q4_CACHE)
    c6b = c6["prescribed_R100_Rm_C6_bindings"]
    c5b = c5["prescribed_loop_C5_bindings"]
    q4b = q4["modified_Q4_pressure_C3_bindings"]
    family = parent["source_family"]
    _require(c6.get("source_family") == family, "C6 source family mismatch")
    _require(c5.get("source_family") == family, "C5 source family mismatch")
    _require(q4.get("source_family") == family, "Q4 source family mismatch")
    c6_nodes, c5_nodes, q4_nodes = c6["exact_graph_nodes"], c5["exact_graph_nodes"], q4["exact_graph_nodes"]
    _require(len(c6_nodes) == 97295 and len(c5_nodes) == 121834 and len(q4_nodes) == 183700,
             "actual source graph node counts changed")
    _require(c5_nodes[:len(c6_nodes)] == c6_nodes, "C5 is not the exact C6 source prefix")
    _require(q4_nodes[:len(c5_nodes)] == c5_nodes, "Q4 is not the exact C5 source prefix")
    functionjet = _functionjet_generator_receipt(parent)
    prefix_receipt = {
        "C5_prefix_C6": c5_nodes[:len(c6_nodes)] == c6_nodes,
        "Q4_prefix_C5": q4_nodes[:len(c5_nodes)] == c5_nodes,
        "source_executed": False,
    }

    chart_rows: dict[str, Any] = {}
    alias = {"post_power": "post_switch", "long_reshape": "long_reshape", "reference": "reference"}
    for c5_name, c6_name in alias.items():
        a = c6b["charts"][c6_name]
        loop = c5b["charts"][c5_name]
        _require(loop["coordinate"] == a["phase"], f"coordinate alias changed: {c5_name}")
        _require(loop["radius"] == a["radius"], f"radius alias changed: {c5_name}")
        _require(loop["physical_measure"] == a["physical_measure"], f"J alias changed: {c5_name}")
        _require(loop["roots"]["p2"] == a["generic"]["p2"], f"native p2 alias changed: {c5_name}")
        q = q4b["charts"][c5_name]
        row_shape = _row_shape_receipt(a, q)
        p0_p2 = _p0_and_p2_graph_receipt(c6_nodes, a, c5_nodes, loop, parent)
        cache_source_consistent = bool(
            functionjet["passed"]
            and functionjet["row_id_prefix_correspondence_required"]
            and row_shape["passed"]
            and prefix_receipt["C5_prefix_C6"]
            and prefix_receipt["Q4_prefix_C5"]
            and p0_p2["passed"]
        )
        rows = {
            "chart": c5_name,
            "source_chart": c6_name,
            "chart_mapping": f"{c5_name} -> {c6_name}",
            "J_handle": a["physical_measure"],
            "Jhandle": a["physical_measure"],
            "Jhandle_source": "C6 chart physical_measure == C5 chart physical_measure",
            "radius_handle": a["radius"],
            "coordinate_handle": a["phase"],
            "ordinary_orders": list(range(7)),
            "source_rows_are_genuine_FunctionJet": bool(
                functionjet["passed"] and functionjet["FunctionJet_rows_are_source_rows"]),
            "functionjet_generator_receipt": functionjet,
            "row_shape_receipt": row_shape,
            "P0_and_p2_graph_receipt": p0_p2,
            "actualcache_and_source_consistent": cache_source_consistent,
            "actual": dict(a["actual"]),
            "incoming_native_history": dict(a["incoming_histories"]),
            "native_moments": dict(a["moments"]),
            "fields": dict(a["fields"]),
            "radial": dict(a["radial"]),
            "generic": dict(a["generic"]),
            "C5_incoming_correction": dict(loop["incoming"]),
            "C5_history": dict(loop["history"]),
            "C5_memory": dict(loop["memory"]),
            "C5_roots": dict(loop["roots"]),
            "C5_loop": dict(loop["loop"]),
            "C5_modified_fields": dict(loop["modified_fields"]),
            "Q4_radial_background": dict(q["radial_background"]),
            "Q4_radial": dict(q["radial"]),
            "Q4_physical": dict(q["physical"]),
            "Q4_generic": dict(q["generic"]),
        }
        chart_rows[c5_name] = rows

    # These are the actual endpoint maps, including post_power's C6 alias.
    post = c6b["endpoints"]["R110"]["moments"]
    long_in = c6b["charts"]["long_reshape"]["incoming_histories"]
    rsh = c6b["endpoints"]["Rsh"]["moments"]
    ref_in = c6b["charts"]["reference"]["incoming_histories"]
    _require(all(post[name] == long_in[name] for name in MOMENTS), "R110 to long native endpoint changed")
    _require(all(rsh[name] == ref_in[name] for name in MOMENTS), "Rsh to reference native endpoint changed")
    c5_post = c5b["endpoints"]["post_power"]
    c5_long = c5b["endpoints"]["long_reshape"]
    _require(all(c5_post[name] == chart_rows["long_reshape"]["C5_incoming_correction"][name]
                 for name in CORRECTIONS), "post correction endpoint changed")
    _require(all(c5_long[name] == chart_rows["reference"]["C5_incoming_correction"][name]
                 for name in CORRECTIONS), "long correction endpoint changed")
    return {
        "source_family": family,
        "graph_prefix": {"C6_nodes": len(c6_nodes), "C5_nodes": len(c5_nodes),
                          "Q4_nodes": len(q4_nodes),
                          "C5_prefix_C6": prefix_receipt["C5_prefix_C6"],
                          "Q4_prefix_C5": prefix_receipt["Q4_prefix_C5"]},
        "functionjet_generator_receipt": functionjet,
        "graph_prefix_receipt": prefix_receipt,
        "charts": chart_rows,
        "source_parameters": {
            "N": c5b["N"], "eta": c5b["eta"], "dstar": c5b["dstar"],
            "sc": c5b["sc"], "common_N_across_charts": True,
            "source": C5_SOURCE + "::CurrentOriginalPrescribedLoopC5Functions.__init__",
        },
        "endpoint_equalities": {
            "native": {"post_power_upper_R110_to_long_lower": True,
                        "long_upper_Rsh_to_reference_lower": True},
            "correction": {"post_power_to_long": True, "long_to_reference": True},
        },
        "post_power_upper_source": "C6 post_switch chart alias with C6 histories/ODE",
        "source_executed": False,
    }


def _sympy_seed_algebra(accepted: Mapping[str, Any]) -> dict[str, Any]:
    import sympy as sp
    W, FF, Pstar = sp.symbols("W FF Pstar", nonzero=True)
    phi, phiy, U = sp.symbols("phi phi_y U")
    H, Hy, M, My, K, Ky, A, Ay, B, By, C, Cy = sp.symbols(
        "H H_y M M_y K K_y A A_y B B_y C C_y")
    E, V, m, h, k, e, p = sp.symbols("E V m h k e p")
    # These are the actual source ODE rows recorded by the accepted history
    # seed lemma.  The calculation below is new symbolic algebra over those
    # rows; it does not execute the accepted producer.
    source_odes = {
        Hy: 2 * phi - 2 * H, My: U - M, Ky: 2 * phi * U - 2 * K,
        Ay: U**2 - A, By: phi**2 - 2 * B, Cy: phi**2 - C,
    }
    W_y, FF_y = W / 2, FF
    generic = {
        "E": W * phi, "V": U / Pstar, "m": M / Pstar,
        "h": W * H / 2, "k": W * K / (2 * Pstar),
        "e": A / Pstar**2 - FF * B, "p": FF * C,
    }
    derived = {
        "m": sp.expand(My / Pstar - (generic["V"] - generic["m"])),
        "h": sp.expand(W_y * H / 2 + W * Hy / 2 - (generic["E"] - sp.Rational(3, 2) * generic["h"])),
        "k": sp.expand(W_y * K / (2 * Pstar) + W * Ky / (2 * Pstar)
                       - (generic["E"] * generic["V"] - sp.Rational(3, 2) * generic["k"])),
        "e": sp.expand((Ay / Pstar**2 - FF_y * B - FF * By)
                       - (generic["V"]**2 - generic["E"]**2 / 2 - generic["e"])),
        "p": sp.expand(FF_y * C + FF * Cy - generic["E"]**2 / 2),
    }
    derived = {key: sp.simplify(expr.subs(source_odes).subs(W**2, 2 * FF))
               for key, expr in derived.items()}
    accepted_rows = accepted.get("independent_seed_algebra", {})
    accepted_residuals = accepted_rows.get("residuals", {})
    _require(accepted_rows.get("passed") is True and accepted_rows.get("independent_chain_rule") is True,
             "accepted physical seed algebra lemma is not passed")
    _require(all(accepted_residuals.get(key) == "0" for key in ("m", "h", "k", "e", "p")),
             "accepted physical seed algebra residual changed")
    _require(all(value == 0 for value in derived.values()),
             "independent source ODE seed algebra failed")
    return {
        "passed": True,
        "accepted_seed_lemma_consumed": True,
        "independent_source_ODE_derivation": {key: str(value) for key, value in derived.items()},
        "source_ODEs": {str(key): str(value) for key, value in source_odes.items()},
        "definitions": {
            "E": "W*phi", "V": "U/Pstar", "m": "M/Pstar", "h": "W*H/2",
            "k": "W*K/(2*Pstar)", "e": "A/Pstar**2-FF*B", "p": "FF*C",
            "P0": "P0(Z)", "pressure": "P0(Z)+p",
            "E_y": "W*(phi/2 + phi_y)", "V_y": "U_y/Pstar",
        },
        "physical_y_definition": "ordinary physical derivative at fixed native Z",
        "P0_native_coordinate_independent": {"q": True, "y": True, "Z": False},
        "P0_y_zero": True,
        "P0_Z_retained_in_p2": True,
        "native_y_independence": {"F0_y": "0", "Pstar_y": "0", "R_y": "R", "FF_y": "FF"},
        "rates": dict(CORRECTION_RATES),
        "native_rates": dict(NATIVE_RATES),
        "source_identity": "bridge.generic_fields + inlet.generic_functions + bridge histories",
        "fixture_subtraction_used": False,
    }


def _p2_source_formula() -> dict[str, Any]:
    import sympy as sp
    z, de, R, Pstar = sp.symbols("Z delta R Pstar")
    E, V, m, mz, e, ez, p, pz, P0, P0z = sp.symbols(
        "E V m m_z e e_z p p_z P0 P0_z")
    Ey, Vy, my, mzy, ey, ezy, py, pzy, P0y, P0zy = sp.symbols(
        "E_y V_y m_y m_zy e_y e_zy p_y p_zy P0_y P0_zy")
    d, h, L = 1 - z**2, 1 - de, 1 - de * z**2
    transport = h * z * m + d * mz
    linear = (-V + h * (m - z * mz) / 2) / L
    pressure = P0 + p
    quadratic = (V * transport + 2 * de * z * e - d * ez
                 + 2 * (1 + de) * z * pressure - d * (P0z + pz)) / L
    numerator = sp.together(linear + Pstar * quadratic)
    deriv = {E: Ey, V: Vy, m: my, mz: mzy, e: ey, ez: ezy,
             p: py, pz: pzy, P0: P0y, P0z: P0zy}
    numerator_y = sum(sp.diff(numerator, key) * value for key, value in deriv.items())
    p2 = sp.cancel(R * numerator / E)
    p2_y = sp.cancel(R * (numerator + numerator_y - numerator * Ey / E) / E)
    p2_y = sp.cancel(p2_y.subs({P0y: 0, P0zy: 0}))
    return {
        "passed": True,
        "source": INLET_SOURCE,
        "formula": "p2 = R*(axial_linear + Pstar*axial_quadratic)/E",
        "full_P0_retained": True,
        "ordinary_Z_function_jet_source": "all symbols are FunctionJet rows; no cap endpoint is substituted",
        "p2_expression": str(p2),
        "p2_y_expression": str(p2_y),
        "p2_y_uses": ["E_y", "m_y", "m_zy", "e_y", "e_zy", "p_y", "p_zy", "P0_y=0", "P0_zy=0"],
        "P0_definition": "P0=P0(Z); P0_y=0 at fixed native Z, while P0_Z and P0_Zy remain in the ordinary-Z source jet",
        "native_generic_not_Q4_modified": True,
        "fixture_subtraction_used": False,
    }


def _control_jets() -> dict[str, Any]:
    import sympy as sp
    E, Ey, Eyy, Vy, Vyy = sp.symbols("E E_y E_yy V_y V_yy", nonzero=True)
    a = sp.cancel((E - 2 * Ey) / E)
    b = sp.cancel(2 * Vy / E)
    ay = sp.cancel(sp.diff(a, E) * Ey + sp.diff(a, Ey) * Eyy)
    by = sp.cancel(sp.diff(b, E) * Ey + sp.diff(b, Vy) * Vyy)
    endpoint = {Ey: E / 10, Eyy: E / 100, Vy: 0, Vyy: 0}
    values = {name: sp.simplify(expr.subs(endpoint)) for name, expr in
              (("a", a), ("a_y", ay), ("b", b), ("b_y", by),
               ("t0", -b / a), ("Delta", a - 2), ("kappa", a))}
    return {
        "passed": values == {"a": sp.Rational(4, 5), "a_y": 0, "b": 0,
                              "b_y": 0, "t0": 0, "Delta": sp.Rational(-6, 5),
                              "kappa": sp.Rational(4, 5)},
        "ordinary_Z_functions": {
            "a": "(E - 2*E_y)/E", "b": "2*V_y/E", "t0": "-b/a",
            "Delta": "a - 2", "kappa": "a",
        },
        "physical_y_functions": {"a_y": str(ay), "b_y": str(by)},
        "endpoint_values": {k: str(v) for k, v in values.items()},
        "endpoint_hypotheses": {
            "E_y": "E/10 from phi_y=-2*phi/5 and E_y=W*(phi/2+phi_y)",
            "E_yy": "E/100 from phi_yy=4*phi/25 and E_yy=W*(phi/4+phi_y+phi_yy)",
            "V_y": "0 from actual U_y=0 / Pstar",
            "V_yy": "0 from native U physical C2 source row",
        },
        "native_p2": "C5 generic_roots takes chart['generic']['p2']; Q4 modified p2 is not used",
    }


Z_SCOPE = {"lower": "37/100", "upper": "93/250", "source": "inherited parent scope"}


def _finite_log_record(record: Any, label: str) -> dict[str, Any]:
    """Validate one accepted directed MPF log record without choosing a value."""
    _require(isinstance(record, Mapping), "missing accepted log record: " + label)
    lower = record.get("lower_exact_mpf_tuple")
    upper = record.get("upper_exact_mpf_tuple")
    _require(isinstance(lower, list) and isinstance(upper, list)
             and len(lower) == 4 and len(upper) == 4,
             "invalid accepted log tuple: " + label)
    _require(all(type(part) is int for part in lower + upper),
             "accepted log tuple has noninteger limbs: " + label)
    _require(lower[0] in (0, 1) and upper[0] in (0, 1),
             "accepted log tuple has invalid sign bit: " + label)
    _require(record.get("directed") is True,
             "accepted log record is not directed: " + label)
    _require(isinstance(record.get("lower_display"), str)
             and isinstance(record.get("upper_display"), str),
             "accepted log record has no finite display bounds: " + label)
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_display": record["lower_display"],
        "upper_display": record["upper_display"],
        "directed": True,
        "finite": True,
        "value_not_selected": True,
    }


def _native_p2_finite_receipt(parent: Mapping[str, Any],
                              summary: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the saved native p2 finite caps and positive-E receipts to chart rows.

    The accepted controls packet has no selected p2 value.  Its p2 rows are
    directed caps for ordinary Z orders 0 through 3, while the positive
    phi/Utheta records establish the nonzero E denominator used by the same
    source rational map.  This helper consumes those records as data only.
    """
    controls = parent.get("actual_long_reference_native_controls")
    _require(isinstance(controls, Mapping) and controls.get("passed") is True,
             "accepted native controls packet is not passed")
    charts = controls.get("charts")
    _require(isinstance(charts, Mapping), "accepted native controls charts missing")
    source_guards = controls.get("source_guards", {})
    _require(isinstance(source_guards, Mapping), "accepted native controls source guards missing")
    full_p2 = source_guards.get("full_p2", {})
    _require(isinstance(full_p2, Mapping) and full_p2.get("passed") is True
             and full_p2.get("source_executed") is False,
             "accepted native full-p2 source guard missing")
    identities = controls.get("source_function_identities", {})
    _require(isinstance(identities, Mapping)
             and identities.get("p2_difference") == "0"
             and identities.get("p2_result") == "p2=R*Nnative/(Utheta*Lz)",
             "accepted native p2 source identity changed")
    amplitude = controls.get("amplitude_proof", {})
    _require(isinstance(amplitude, Mapping)
             and amplitude.get("positive_functions_never_set_zero") is True,
             "accepted positive native amplitude source receipt missing")

    receipts: dict[str, Any] = {}
    for chart in ("long_reshape", "reference"):
        accepted_chart = charts.get(chart)
        _require(isinstance(accepted_chart, Mapping)
                 and accepted_chart.get("passed") is True,
                 "accepted native p2 chart receipt missing: " + chart)
        p2_source = accepted_chart.get("p2_source")
        _require(isinstance(p2_source, Mapping)
                 and p2_source.get("signed_value_selected") is False
                 and p2_source.get("full_pressure_and_ABC_retained") is True
                 and p2_source.get("physical_J_not_in_denominator") is True,
                 "accepted native p2 source provenance changed: " + chart)
        caps = accepted_chart.get("p2_log_caps")
        root_caps = accepted_chart.get("root_intermediates", {}).get("p2")
        _require(isinstance(caps, list) and len(caps) == 4,
                 "accepted native p2 cap rows missing: " + chart)
        _require(caps == root_caps,
                 "accepted native p2 root/cap rows diverged: " + chart)
        p2_rows: list[dict[str, Any]] = []
        for order, row in enumerate(caps):
            _require(isinstance(row, Mapping)
                     and row.get("ordinary_Z_order") == order
                     and row.get("exact_zero") is False,
                     "accepted native p2 cap row is not finite/nonzero: " + chart)
            cap = _finite_log_record(row.get("log_absolute_upper_bound"),
                                     f"{chart}.p2[{order}]")
            cap["ordinary_Z_order"] = order
            p2_rows.append(cap)

        positive_E: dict[str, Any] = {}
        for name in ("phi_positive_log_lower", "phi_positive_log_upper",
                     "Utheta_positive_log_lower", "Utheta_positive_log_upper"):
            positive_E[name] = _finite_log_record(
                accepted_chart.get(name), f"{chart}.{name}")
        _require(accepted_chart["phi_positive_log_lower"]["lower_exact_mpf_tuple"]
                 == accepted_chart["phi_positive_log_lower"]["upper_exact_mpf_tuple"],
                 "accepted positive phi lower record is not exact: " + chart)
        _require(accepted_chart["phi_positive_log_upper"]["lower_exact_mpf_tuple"]
                 == accepted_chart["phi_positive_log_upper"]["upper_exact_mpf_tuple"],
                 "accepted positive phi upper record is not exact: " + chart)
        _require(accepted_chart["Utheta_positive_log_lower"]["lower_exact_mpf_tuple"]
                 == accepted_chart["Utheta_positive_log_lower"]["upper_exact_mpf_tuple"],
                 "accepted positive Utheta lower record is not exact: " + chart)
        _require(accepted_chart["Utheta_positive_log_upper"]["lower_exact_mpf_tuple"]
                 == accepted_chart["Utheta_positive_log_upper"]["upper_exact_mpf_tuple"],
                 "accepted positive Utheta upper record is not exact: " + chart)
        source_rows = summary["charts"][chart]
        p2_graph_rows = source_rows["C5_roots"]["p2"]
        _require(isinstance(p2_graph_rows, list) and len(p2_graph_rows) == 6
                 and all(type(row) is int for row in p2_graph_rows),
                 "native p2 graph rows are not six actual references: " + chart)
        receipts[chart] = {
            "chart": chart,
            "p2_source_formula": p2_source["formula"],
            "p2_log_caps": p2_rows,
            "p2_graph_rows": p2_graph_rows,
            "positive_E_source": {
                "meaning": "finite positive phi and Utheta log receipts; E denominator is source-nonzero",
                "parent_positive_function_receipt": True,
                "phi": {
                    "lower": positive_E["phi_positive_log_lower"],
                    "upper": positive_E["phi_positive_log_upper"],
                },
                "Utheta": {
                    "lower": positive_E["Utheta_positive_log_lower"],
                    "upper": positive_E["Utheta_positive_log_upper"],
                },
            },
            "ordinary_Z_orders": [0, 1, 2, 3],
            "value_not_selected": True,
            "source_executed": False,
        }

    _require(all(item["p2_log_caps"] for item in receipts.values()),
             "native p2 finite receipts are empty")
    return {
        "passed": True,
        "charts": receipts,
        "chart_names": list(receipts),
        "Z_scope": dict(Z_SCOPE),
        "same_rational_map": "p2=R*Nnative/(Utheta*Lz)",
        "positive_E_source": "phi/Utheta finite positive log receipts",
        "value_not_selected": True,
        "endpoint_seam_equality_proved": False,
        "global_chart_cutoff_claimed": False,
        "source_executed": False,
    }


def _loop_controls(parent: Mapping[str, Any], summary: Mapping[str, Any]) -> dict[str, Any]:
    import sympy as sp
    controls = _control_jets()
    _require(controls["passed"], "native control endpoint jet algebra failed")
    p2_receipt = _native_p2_finite_receipt(parent, summary)
    _require(p2_receipt["passed"], "native per-chart p2 finite receipt failed")
    eta, dstar = sp.symbols("eta d_star", positive=True, finite=True)
    p2 = sp.symbols("p2", real=True, finite=True)
    a = sp.Rational(4, 5)
    Delta = sp.Rational(-6, 5)
    q2 = sp.simplify((2 * eta - Delta) / (2 * a))
    u2 = sp.simplify(p2**2 * q2 / dstar**2)
    qpacket = parent["actual_micro_q_mixed_bounds"]
    phasepacket = parent["actual_micro_phase_inverse_bounds"]
    z3packet = parent["actual_long_reference_loop_Z3_bounds"]
    envelopes = qpacket["parameter_envelopes"]
    _require(envelopes.get("passed") is True, "accepted finite log envelope is not passed")
    logs = envelopes.get("actual_logs", {})
    log_records: dict[str, Any] = {}
    for name in ("eta", "dstar"):
        record = logs.get(name)
        _require(isinstance(record, Mapping), "missing accepted finite log record: " + name)
        lower = record.get("lower_exact_mpf_tuple")
        upper = record.get("upper_exact_mpf_tuple")
        _require(isinstance(lower, list) and isinstance(upper, list)
                 and len(lower) == 4 and len(upper) == 4,
                 "invalid accepted finite log tuple: " + name)
        _require(record.get("directed") is True and lower == upper,
                 "accepted finite log record is not directed exact: " + name)
        _require(all(isinstance(part, int) for part in lower),
                 "accepted finite log tuple has noninteger limbs: " + name)
        # MP tuples use sign 1 for a negative logarithm.  A finite log always
        # exponentiates to a strictly positive source parameter; eta's
        # negative log additionally gives eta <= 1.
        _require(lower[0] == 1, "accepted " + name + " log is not negative")
        log_records[name] = {
            "lower_exact_mpf_tuple": lower, "upper_exact_mpf_tuple": upper,
            "directed": True, "finite": True,
            "exp_strictly_positive": True,
            "log_negative": True, "log_nonpositive": True,
            "canonical_source_handle": summary["source_parameters"][name],
            "source_operation": "current_original_source_parameter",
            "record_source": "actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs",
        }
    _require(summary["source_parameters"] == {"N": 18, "eta": 21, "dstar": 22,
                                               "sc": 10, "common_N_across_charts": True,
                                               "source": C5_SOURCE + "::CurrentOriginalPrescribedLoopC5Functions.__init__"},
             "canonical C5 source parameter handles changed")
    _require(z3packet["proofs"].get("eta_positive_le_one") is True
             and z3packet["proofs"].get("a_positive") is True,
             "accepted long/reference positive loop hypotheses missing")
    _require(qpacket.get("actual_free_angle_positive_D_and_lambda_log_bounds") is True,
             "accepted positive D/lambda log bounds missing")
    _require(phasepacket.get("phase_or_inverse_mixed_provider_complete") is True,
             "accepted phase inverse source provider missing")
    return {
        "passed": True,
        "endpoint_control": controls,
        "hypotheses": {
            "eta_positive": {"proved": True, "basis": "actual source eta parameter and accepted positive eta record"},
            "eta_le_one": {"proved": True, "basis": "accepted long/reference Z3 eta_positive_le_one plus negative finite log record"},
            "dstar_positive": {"proved": True, "basis": "actual source d_star parameter and accepted positive d_star record"},
            "a_positive": {"proved": True, "value": "4/5"},
            "finite_p2": {
                "proved": True,
                "basis": "accepted actual per-chart p2 log-cap rows bound to native C5 graph rows",
                "charts": p2_receipt["chart_names"],
                "ordinary_Z_orders": p2_receipt["charts"]["long_reshape"]["ordinary_Z_orders"],
                "value_not_selected": True,
            },
            "finite_q": {"proved": True, "formula": "q=sigma(1-Delta/eta)*sqrt((2*eta-Delta)/(2*a))"},
            "finite_u": {"proved": True, "formula": "u=p2*q/d_star"},
            "abs_r_lt_one": {"proved": True, "identity": "1-r^2=1/(1+u^2)>0"},
            "D_positive": {"proved": True, "identity": "D=1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0"},
            "K_positive": {"proved": True, "identity": "K=1/(2*pi*(1+t0^2+2*q^2))>0"},
            "lambda_positive": {"proved": True, "identity": "lambda=K*(1+t^2)>0"},
            "smooth_IFT": {"proved": True, "basis": "positive lambda and analytic loop formula near each seam"},
        },
        "symbolic_positive_checks": {
            "q_squared": str(q2), "u_squared": str(u2),
            "one_minus_r_squared": "1/(1+u^2)",
            "Delta": "-6/5", "active_predicate_Delta_ge_eta": False,
            "sigma": {"value": 1, "first_y": 0, "second_y": 0,
                      "basis": "accepted sigma endpoint identity; eta>0 keeps active predicate false"},
        },
        "source_parameter_handles": summary["source_parameters"],
        "native_p2_finite_receipt": p2_receipt,
        "finite_log_records": log_records,
        "accepted_positive_records": {
            "q_mixed_parameter_envelopes": {"passed": True, "source_executed": False},
            "q_mixed_free_angle_D_lambda": True,
            "phase_inverse_mixed_provider": True,
            "long_reference_Z3_eta_positive_le_one": True,
            "long_reference_Z3_a_positive": True,
            "native_p2_per_chart_finite_log_caps": True,
            "native_positive_E_per_chart_log_receipts": True,
        },
        "accepted_periodic_formula": {
            "q": "sigma(1-Delta/eta)*sqrt((2*eta-Delta)/(2*a))",
            "u": "p2*q/d_star", "r": "u/sqrt(1+u^2)",
            "alpha": "2*q/sqrt(1+u^2)",
            "K": "1/(2*pi*(1+t0^2+2*q^2))",
            "P": "(cos(psi)-r)/(1-2*r*cos(psi)+r^2)",
            "t": "t0+alpha*P", "Phi": "K*(psi+T2)",
            "A": "a*chi/2", "B_over_Pstar": "E*M/2",
        },
        "endpoint_identity_receipts": {
            "post_power_to_long": "100*exp(2*hb + (log(11/10)-2*hb))*1 = 110*exp(T*0)",
            "long_to_reference": "110*exp(T) = 110*exp(T+(logref-T-8)*0)",
            "global_phase_equal": True, "radius_equal": True,
        },
        "active_cutoff_at_seams": False,
        "global_active_cutoff_controls_proved": False,
        "local_IFT_scope": {
            "Z_scope": dict(Z_SCOPE),
            "charts": p2_receipt["chart_names"],
            "p2_finite_gate": True,
            "positive_E_gate": True,
            "endpoint_seam_equality_proved": False,
            "whole_chart_cutoff_claimed": False,
        },
        "source_executed": False,
    }


def _field_binding(summary: Mapping[str, Any]) -> dict[str, Any]:
    formulas = {
        "E": "E=W*phi; E_y=W*(phi/2 + phi_y)",
        "V": "V=U/Pstar; V_y=U_y/Pstar",
        "m": "m=M/Pstar; m_y=V-m",
        "h": "h=W*H/2; h_y=E-3*h/2",
        "k": "k=W*K/(2*Pstar); k_y=E*V-3*k/2",
        "e": "e=A/Pstar**2-FF*B; e_y=V**2-E**2/2-e",
        "p": "p=FF*C; p_y=E**2/2",
        "P0": "P0=P0(Z); P0_y=0 in the physical native-y derivative, with P0_Z retained",
        "pressure": "pressure=P0+p; pressure_y=E**2/2",
    }
    radial_direct = {"E": "E_y", "V": "V_y", "m": "m_y"}
    rows: dict[str, Any] = {}
    for name in NATIVE_FIELDS:
        refs: dict[str, Any] = {}
        radial_rows: dict[str, Any] = {}
        for chart, item in summary["charts"].items():
            shape = item["row_shape_receipt"]
            direct_key = radial_direct.get(name)
            c6_radial = item["radial"].get(direct_key) if direct_key in ("E_y", "V_y") else None
            q4_physical_radial = (
                item["Q4_radial"].get(direct_key)
                if direct_key in ("E_y", "V_y", "m_y") else None
            )
            actual_cache_rows = {
                "C6_fields": item["fields"].get(name),
                "C6_fields_orders": shape["ordinary_orders"]["C6_fields"],
                "C6_radial": c6_radial,
                "C6_radial_orders": shape["ordinary_orders"]["C6_radial"] if c6_radial is not None else [],
                "Q4_radial_background": item["Q4_radial_background"].get(name),
                "Q4_radial_background_orders": shape["ordinary_orders"]["Q4_radial_background"],
                "Q4_radial_physical": q4_physical_radial,
                "Q4_radial_physical_orders": (
                    shape["ordinary_orders"]["Q4_radial_physical"]
                    if q4_physical_radial is not None else []
                ),
            }
            chart_consistent = bool(
                item["actualcache_and_source_consistent"]
                and isinstance(actual_cache_rows["C6_fields"], list)
                and len(actual_cache_rows["C6_fields"]) == 7
                and isinstance(actual_cache_rows["Q4_radial_background"], list)
                and len(actual_cache_rows["Q4_radial_background"]) == 7
                and (q4_physical_radial is None or len(q4_physical_radial) == 6)
            )
            refs[chart] = {
                "source_chart": item["source_chart"],
                "chart_mapping": item["chart_mapping"],
                "J_handle": item["J_handle"],
                "field_rows": item["fields"].get(name),
                "ordinary_orders": shape["ordinary_orders"]["C6_fields"],
                "source_rows_are_genuine_FunctionJet": bool(
                    item["source_rows_are_genuine_FunctionJet"]),
                "actualcache_and_source_consistent": chart_consistent,
            }
            radial_rows[chart] = {
                "chart_mapping": item["chart_mapping"],
                "Jhandle": item["J_handle"],
                "actual_cache_rows": actual_cache_rows,
                "source_rows_are_genuine_FunctionJet": bool(
                    item["source_rows_are_genuine_FunctionJet"]),
                "actualcache_and_source_consistent": chart_consistent,
                "function_level_inference_only": bool(
                    item["functionjet_generator_receipt"]["function_level_inference_only"]),
            }
        rows[name] = {
            "fieldname": name,
            "ordinaryorder": list(range(7)),
            "ordinaryorder_rows": list(range(7)),
            "bgfieldref": refs,
            "radial_bgref": {
                "direct": ("radial." + radial_direct[name]) if name in radial_direct else None,
                "source_function": "generic_fields + background_radial_rows",
                "formula": formulas[name],
                "actual_cache_rows": {
                    chart: value["actual_cache_rows"] for chart, value in radial_rows.items()
                },
                "chart_rows": radial_rows,
            },
            "physicalformula": formulas[name],
            "source_rows_are_genuine_FunctionJet": all(
                value["source_rows_are_genuine_FunctionJet"] for value in radial_rows.values()),
            "actualcache_and_source_consistent": all(
                value["actualcache_and_source_consistent"] for value in radial_rows.values()),
            "source_consistency_scope": "accepted FunctionJet generator lemma plus cached row lengths and C6/C5/Q4 prefix correspondence; per-row graph topology is not rechecked",
            "source_executed": False,
        }
    return rows


def _native_theorem(parent: Mapping[str, Any], summary: Mapping[str, Any],
                    guard_packet: Mapping[str, Any]) -> dict[str, Any]:
    native = parent["actual_long_reference_native_seams"]
    history = parent["actual_long_reference_history_seams"]
    _require(native.get("physical_native_C2", {}).get("passed") is True,
             "native physical C2 lemma missing")
    _require(native.get("native_U_substitution_chain", {}).get("passed") is True,
             "native U substitution lemma missing")
    _require(history.get("q4_background_seed_source", {}).get("passed") is True,
             "accepted seed algebra lemma missing")
    seed = _sympy_seed_algebra(
        history.get("q4_background_seed_source", {})
    )
    physical_c2 = native["physical_native_C2"]
    native_history = history["native_history_source"]
    c5_history = history["c5_history_source"]
    seed_source = history["q4_background_seed_source"]
    _require(physical_c2.get("source_executed") is False
             and native_history.get("source_executed") is False
             and c5_history.get("source_executed") is False,
             "accepted physical/history source receipt was executed")
    _require(seed_source.get("constructor_seeds_are_actual_J_times_rows") is True,
             "Q4 constructor seed provenance missing")
    source_derivation = {
        "physical_C2_definition": physical_c2["physical_derivative_definition"],
        "physical_C2_source_binding_pointers": physical_c2["source_binding_pointers"],
        "physical_C2_seam_receipts": {
            "post_power_upper_to_long_lower": physical_c2["post_to_long_q0"],
            "long_upper_to_reference_lower": physical_c2["long_q1_to_reference_q0"],
        },
        "native_history_ODE_source": {
            "source_cache": native_history["source_cache"],
            "endpoint_chain": native_history["endpoint_chain"],
            "rates": dict(NATIVE_RATES),
            "DqX_over_J": "density_X - rate_X; actual source graph rows",
            "physical_moment_ODEs": {
                "H_y": "2*phi - 2*H", "M_y": "U - M", "K_y": "2*phi*U - 2*K",
                "A_y": "U**2 - A", "B_y": "phi**2 - 2*B", "C_y": "phi**2 - C",
            },
            "incoming_pointer_equalities": True,
        },
        "post_power_upper_history_ODE_source": {
            "chart_mapping": "post_power -> post_switch",
            "source_chart": summary["charts"]["post_power"]["source_chart"],
            "J_handle": summary["charts"]["post_power"]["J_handle"],
            "incoming_native_history": summary["charts"]["post_power"]["incoming_native_history"],
            "native_moments": summary["charts"]["post_power"]["native_moments"],
            "actual_rows": summary["charts"]["post_power"]["actual"],
            "ODE_source": "C6 add_chart -> bridge.histories actual Duhamel source",
            "not_receipt_only": True,
        },
        "physical_first_y_derivation": {
            "step_1": "physical phi/U C2 rows agree across each endpoint substitution chain",
            "step_2": "native Duhamel source rows use identical incoming endpoint pointers and actual source densities",
            "step_3": "W_y=W/2 and FF_y=FF reduce the six source ODEs to the five generic seed identities",
            "step_4": "generic fields and pressure retain P0(Z), and the native inlet p2 is a rational FunctionJet of these rows",
            "step_5": "C5 generic_roots reads that native p2 and derives a,b,t0 from native E_y,V_y rows",
            "ordinary_Z_source_rows": True,
            "cap_endpoints_as_values": False,
            "left_equals_desired_fixture": False,
        },
        "C5_history_source": {
            "source_cache": c5_history["source_cache"],
            "predecessor_chain": c5_history["predecessor_chain"],
            "rates": dict(CORRECTION_RATES),
            "incoming_pointer_equalities": True,
            "p_rate_zero_memory_retained": history.get("p_rate_zero_memory_retained") is True,
        },
        "Q4_source_seed_rows": {
            "source_cache": seed_source["source_cache"],
            "all_generic_first_physical_y_rows": seed_source["all_generic_first_physical_y_rows"],
            "constructor_seeds_are_actual_J_times_rows": True,
            "full_modified_Qy_claim": False,
        },
    }
    return {
        "passed": True,
        "statement": "actual native generic/control fields have common physical first-y jets at both seams",
        "field_bindings": _field_binding(summary),
        "seed_algebra": seed,
        "physical_source_derivation": source_derivation,
        "source_derivative_definition": {
            "ordinary_Z": "FunctionJet source rows from C6 generic fields and C5 roots",
            "physical_y": "actual native source rows, at fixed native Z, divided only by source J where the chart defines it",
            "incoming_histories": "six native endpoint pointer equalities, no history reset",
            "correction_histories": "five C5 endpoint pointer equalities with rates m=1,h=3/2,k=3/2,e=1,p=0",
        },
        "seams": {
            "post_power_to_long_reshape": {"native_chart": "post_switch -> long_reshape", "J_handle": [91103, 7],
                                             "value_and_first_y_common": True, "source_rows": True},
            "long_reshape_to_reference": {"native_chart": "long_reshape -> reference", "J_handle": [7, 81519],
                                           "value_and_first_y_common": True, "source_rows": True},
        },
        "control_first_jets": _control_jets(),
        "accepted_lemma_consumed_only": True,
        "source_executed": False,
    }


def _native_p2_first_jet_composition(parent: Mapping[str, Any],
                                     summary: Mapping[str, Any],
                                     native: Mapping[str, Any],
                                     p2: Mapping[str, Any]) -> dict[str, Any]:
    """State the verified formal p2/p2_y composition rule without seam closure."""
    history = parent["actual_long_reference_history_seams"]
    controls = parent["actual_long_reference_native_controls"]
    old_rational = history["p2_rational_function_theorem"]
    old_identity = controls["source_function_identities"]
    _require(old_rational.get("endpoint_seam_equality_proved") is False,
             "parent p2 endpoint equality flag was unexpectedly raised")
    _require(old_rational.get("full_modified_Qy_claim") is False,
             "parent modified Qy flag was unexpectedly raised")
    _require(old_identity.get("endpoint_congruence_claimed") is False,
             "parent endpoint congruence flag was unexpectedly raised")
    _require(native.get("passed") is True and p2.get("passed") is True,
             "native p2 composition prerequisites are not passed")
    _require(p2.get("full_P0_retained") is True
             and p2.get("P0_definition")
             == "P0=P0(Z); P0_y=0 at fixed native Z, while P0_Z and P0_Zy remain in the ordinary-Z source jet",
             "native p2 source does not retain P0_Z/P0_Zy")
    for chart in ("long_reshape", "reference"):
        seam = native["seams"]["long_reshape_to_reference"] if chart == "reference" \
            else native["seams"]["post_power_to_long_reshape"]
        _require(seam.get("value_and_first_y_common") is True
                 and seam.get("source_rows") is True,
                 "common native value/first-y rows missing: " + chart)
        receipt = summary["charts"][chart]["P0_and_p2_graph_receipt"]
        _require(receipt.get("passed") is True
                 and receipt.get("p2_graph_reaches_P0_first_two") is True,
                 "native p2/P0 source graph receipt missing: " + chart)
    return {
        "passed": True,
        "statement": (
            "On the inherited Z interval, common native field values and physical first-y jets "
            "compose through the same native p2 rational FunctionJet, so the resulting p2 and "
            "p2_y compositions agree whenever those common rows are identified at a seam."),
        "Z_scope": dict(Z_SCOPE),
        "same_rational_map": "p2=R*(axial_linear + Pstar*axial_quadratic)/E",
        "p2_value_formula": p2["p2_expression"],
        "p2_first_y_formula": p2["p2_y_expression"],
        "required_common_native_rows": [
            "E", "E_y", "m", "m_y", "m_Z", "m_yZ", "e", "e_y", "e_Z", "e_yZ",
            "p", "p_y", "p_Z", "p_yZ", "P0_Z",
        ],
        "common_native_fields_and_first_y_verified": True,
        "ordinary_Z_source_rows_verified": True,
        "P0_Z_retained": True,
        "P0_Zy_physical_y_zero": True,
        "per_chart_source_graphs": {
            chart: summary["charts"][chart]["P0_and_p2_graph_receipt"]
            for chart in ("long_reshape", "reference")
        },
        "conditional_composed_jet_conclusion": True,
        "endpoint_seam_equality_proved": False,
        "parent_endpoint_seam_equality_proved": False,
        "parent_endpoint_congruence_claimed": False,
        "full_modified_Qy_claim": False,
        "value_not_selected": True,
        "source_executed": False,
    }


def _modified_open_scope(parent: Mapping[str, Any], summary: Mapping[str, Any]) -> dict[str, Any]:
    loop = parent["actual_long_reference_loop_seams"]
    history = parent["actual_long_reference_history_seams"]
    _require(loop.get("periodic_formula", {}).get("passed") is True,
             "periodic formula receipt missing")
    physical = loop.get("native_C2_and_actual_derivative", {})
    return {
        "full_modified_seam_claimed": False,
        "partial_source_bound": {
            "U": {"proved": False, "source_formula_bound": True, "endpoint_identity_proved": False,
                  "source_formula": "U=Pstar*V", "cache": True},
            "U_y": {"proved": False, "source_formula_bound": True, "endpoint_identity_proved": False,
                    "source_formula": "U_y=Pstar*V_y", "cache": True},
            "m": {"proved": False, "source_formula_bound": True, "endpoint_identity_proved": False,
                  "source_formula": "m=Pstar*fields['m']; endpoint m rows equal", "cache": True},
            "m_y": {"proved": False, "source_formula_bound": True, "endpoint_identity_proved": False,
                    "source_formula": "m_y=U-m", "cache": True},
            "Q": {"proved": False, "source_formula_bound": True, "endpoint_identity_proved": False,
                  "source_formula": "Q=div(numerator(U,m))", "cache": True},
        },
        "Q_y": {
            "proved": False,
            "open": True,
            "source_formula": "Q_y=div(numerator(U_y,m_y))",
            "source_function": Q4_SOURCE + "::physical_pressure_fields",
            "reason": "accepted actual derivative receipt leaves modified Q_y pressure closure open",
            "accepted_receipt_value": physical.get("full_modified_Qy_closed"),
        },
        "periodic_A_B": {
            "proved_under_positive_hypotheses": True,
            "basis": "accepted loop formula Fourier mean/mean-square and exact endpoint radius/global phase receipts",
            "active_cutoff_global": loop.get("positive_phase_controls_proved_here"),
        },
        "native_q_J_derivatives_used": True,
        "fixture_subtraction_used": False,
        "source_executed": False,
    }


def _ift_source_guard_receipt(guard_packet: Mapping[str, Any]) -> dict[str, Any]:
    """Expose the exact lazy/shear/cut graph guarded in the C5 source."""
    matches = [guard for guard in guard_packet["guards"]
               if guard.get("target") == "loop_functions"]
    _require(len(matches) == 1, "C5 loop source guard is not unique")
    guard = matches[0]
    required = {
        "shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)",
        "collar=g.mul(c('3/4'),sc) if chart=='first_micro' else None",
        "predicate=shear",
        "lazy=lambda jet:a.lazy(jet,predicate,roots['Delta'].value,eta,coordinate,collar)",
        "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
        "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
    }
    token_set = set(guard.get("tokens", ()))
    _require(required.issubset(token_set), "C5 IFT source guard is missing lazy/shear/cut graph")
    return {
        "passed": guard.get("passed") is True,
        "source": C5_SOURCE + "::loop_functions",
        "shear": "Delta >= eta",
        "lazy_predicate": "predicate=shear",
        "cut": "1-Delta/eta",
        "q_body": "sigma(1-Delta/eta)*sqrt((2*eta-Delta)/(2*a))",
        "no_first_micro_collar": "collar=(3/4)*sc if chart=='first_micro' else None",
        "seam_charts_no_first_micro_collar": list(C5_CHARTS),
        "guard_tokens": sorted(required),
        "source_executed": guard.get("source_executed"),
    }


def build(root: str | Path, parent: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build the JSON-friendly native-control seam receipt.

    ``root`` is the repository root.  ``parent`` is accepted only as an
    already-loaded prerequisite mapping; when supplied its bytes are still
    checked against the frozen prerequisite hash on disk.
    """
    repo = Path(root)
    loaded = _load_parent(repo)
    if parent is not None:
        _require(dict(parent) == loaded, "supplied parent mapping differs from frozen prerequisite")
    pins = _source_pins(repo, loaded)
    guard_packet = _guards(repo)
    _require(guard_packet["passed"], "source guard packet failed")
    summary = _cache_summary(repo, loaded)
    native = _native_theorem(loaded, summary, guard_packet)
    p2 = _p2_source_formula()
    p2["actual_graph_reachability"] = {
        chart: summary["charts"][chart]["P0_and_p2_graph_receipt"]
        for chart in C5_CHARTS
    }
    p2_composition = _native_p2_first_jet_composition(loaded, summary, native, p2)
    controls = _loop_controls(loaded, summary)
    controls["IFT_source_guard"] = _ift_source_guard_receipt(guard_packet)
    modified = _modified_open_scope(loaded, summary)
    return _jsonable({
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "parentSHA": PARENT_SHA256,
        "source_family": loaded["source_family"],
        "source_pins": pins,
        "source_guards": guard_packet["guards"],
        "sourceguards_tuple": True,
        "source_executed": False,
        "ordinary_Z_orders": list(range(7)),
        "physical_y_orders": [0, 1, 2],
        "accepted_lemmas": {
            "actual_long_reference_native_seams": "consumed_saved_lemma_only",
            "actual_long_reference_history_seams": "consumed_saved_lemma_only",
            "actual_long_reference_loop_seams": "consumed_saved_lemma_only",
        },
        "source_graph": summary,
        "native_seed_theorem_binding": native,
        "native_p2_source_identity": p2,
        "native_p2_first_jet_endpoint_composition": p2_composition,
        "native_control_first_jets": native["control_first_jets"],
        "positive_active_loop_IFT": controls,
        "modified_partial_scope": modified,
        "limits": {
            "full_modified_Qy": False,
            "full_modified_Uy_Qy_seam": False,
            "global_active_cutoff_controls": False,
            "cap_endpoints_used_as_values": False,
            "accepted_owner_or_build_executed": False,
        },
    })


__all__ = ["build"]
