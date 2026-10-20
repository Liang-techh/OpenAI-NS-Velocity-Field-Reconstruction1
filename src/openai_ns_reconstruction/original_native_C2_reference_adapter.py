"""Standalone constructive native C2 reference/source-function adapter.

The adapter consumes authenticated saved source records and an already pinned
graph as data.  It builds an unevaluated exact-expression DAG for the native
Rsh reference, restoration, and post-restoration charts.  No producer,
operator, recurrence, numerical quadrature, norm, or field evaluator is
imported or called here.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction


SCHEMA = "original_native_C2_reference_adapter_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
PRIOR_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_patch_C2_C6_source_projection.json"
PRIOR_SHA256 = "190d8b54bec0b9cf8e47d6db979fe68eaaf6594e28471f9da4ef86a8c22e88c9"
R110_REPORT = "lei_ren_part1_paper_compliant_current_original_R110_Rm_source_projection.json"
R110_SHA256 = "08ef44c6971d07a4c6583d481ef3a0ce1b37bf069995a732bdf2a9a7ccfa074e"

C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_SOURCE_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_ARCHIVE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
C6_DECODED_SHA256 = "18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a"
C6_CHECK = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"
C6_CHECK_SHA256 = "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6"

TRUE_C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py"
TRUE_C6_SOURCE_SHA256 = "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002"
TRUE_C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.json.gz"
TRUE_C6_ARCHIVE_SHA256 = "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"
TRUE_C6_DECODED_SHA256 = "b5671cc0c60821db9a5c591ea899c9616c3cc949cec4e9214cc57a31ae6cdc10"
TRUE_C6_CHECK = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions_check.json"
TRUE_C6_CHECK_SHA256 = "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7"

REF_SOURCE = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py"
REF_SOURCE_SHA256 = "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7"
REF_REPORT = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.json"
REF_REPORT_SHA256 = "47f21c94f6bb08ebef01288d7b506ef5135ecd2cf5e3fba9cfd15567f6cef440"
REF_CHECK = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions_check.json"
REF_CHECK_SHA256 = "94f5b78791bd571ae6759d2d53e6fdda47607b047138e9c4868e2611dff10b02"

SIGMA_SOURCE = "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
SIGMA_SOURCE_SHA256 = "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7"
SIGMA_REPORT = "lei_ren_part1_paper_compliant_flat_pulse_derivatives.json"
SIGMA_REPORT_SHA256 = "1a833c1f6f2a8d0b76fa1f888b24f9cf14288d8c3a4050493095b9c505eee3f3"
SIGMA_CHECK = "lei_ren_part1_paper_compliant_flat_pulse_derivatives_check.json"
SIGMA_CHECK_SHA256 = "572c50c5210ec51c9998780765bdaf3ca9331ead8f4bc8e968245c6a556f8468"

RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    PRIOR_REPORT: PRIOR_SHA256,
    R110_REPORT: R110_SHA256,
    C6_SOURCE: C6_SOURCE_SHA256,
    C6_ARCHIVE: C6_ARCHIVE_SHA256,
    C6_CHECK: C6_CHECK_SHA256,
    TRUE_C6_SOURCE: TRUE_C6_SOURCE_SHA256,
    TRUE_C6_ARCHIVE: TRUE_C6_ARCHIVE_SHA256,
    TRUE_C6_CHECK: TRUE_C6_CHECK_SHA256,
    REF_SOURCE: REF_SOURCE_SHA256,
    REF_REPORT: REF_REPORT_SHA256,
    REF_CHECK: REF_CHECK_SHA256,
    SIGMA_SOURCE: SIGMA_SOURCE_SHA256,
    SIGMA_REPORT: SIGMA_REPORT_SHA256,
    SIGMA_CHECK: SIGMA_CHECK_SHA256,
}
DECODED_RECORD_HASHES = {
    C6_ARCHIVE: C6_DECODED_SHA256,
    TRUE_C6_ARCHIVE: TRUE_C6_DECODED_SHA256,
}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")).hexdigest()


def _get(value, path, default=None):
    for part in path.split(".") if path else ():
        if not isinstance(value, Mapping) or part not in value:
            return default
        value = value[part]
    return value


def _text(records, name, gates):
    record = records.get(name)
    if not isinstance(record, Mapping):
        gates.append("missing_source_record:" + name)
        return ""
    if record.get("path") != "experiments/root_st073/" + name:
        gates.append("source_path_mismatch:" + name)
    if record.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name)
    value = record.get("source_text")
    if not isinstance(value, str):
        gates.append("source_text_missing:" + name)
        return ""
    if name in DECODED_RECORD_HASHES:
        if record.get("encoding") != "gzip" or record.get("decoded_sha256") != DECODED_RECORD_HASHES[name]:
            gates.append("encoded_archive_contract_mismatch:" + name)
            return ""
        try:
            raw = base64.b64decode(record.get("raw_bytes_base64", ""), validate=True)
            raw_hash = hashlib.sha256(raw).hexdigest()
            if raw_hash != RECORD_HASHES[name]:
                gates.append("archive_raw_hash_mismatch:" + name)
                return ""
            decoded = gzip.decompress(raw)
        except Exception as error:
            gates.append("archive_decode_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(decoded).hexdigest() != DECODED_RECORD_HASHES[name] or decoded != value.encode("utf-8"):
            gates.append("archive_decoded_bytes_mismatch:" + name)
            return ""
        return value
    if record.get("encoding") is not None:
        gates.append("unexpected_encoded_record:" + name)
    if hashlib.sha256(value.encode("utf-8")).hexdigest() != RECORD_HASHES[name]:
        gates.append("source_content_hash_mismatch:" + name)
        return ""
    return value


def _json(records, name, gates):
    text = _text(records, name, gates)
    if not text:
        return {}
    try:
        value = json.loads(text)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}
    return value if isinstance(value, Mapping) else {}


def _ast_source(records, name, gates):
    text = _text(records, name, gates)
    if not text:
        return None
    try:
        return ast.parse(text, filename=name)
    except SyntaxError as error:
        gates.append("source_ast_parse_error:" + name + ":" + str(error.lineno))
        return None


def _names(tree):
    if tree is None:
        return set()
    return {node.name for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def _source_has(text, *needles):
    return all(needle in text for needle in needles)


def _port(dag, archive, prefix, node, label, gates, *, role="graph"):
    ast_function = role == "source_AST_function"
    valid_prefix = isinstance(prefix, str) if ast_function else type(prefix) is int and prefix > 0
    valid_node = isinstance(node, str) if ast_function else type(node) is int and valid_prefix and 0 <= node < prefix
    if not valid_node:
        gates.append("invalid_source_port_node:" + label)
        node = "missing"
    if not valid_prefix:
        gates.append("invalid_source_port_prefix:" + label)
        prefix = "missing"
    return dag.atom("source_function_port", {
        "namespace": SCHEMA,
        "source_archive": archive,
        "archive_sha256": RECORD_HASHES.get(archive),
        "source_prefix": prefix,
        "source_node": node,
        "label": label,
        "role": role,
        "qualified": SCHEMA + "|graph=" + EXPECTED_NODES_DIGEST + "|archive=" + archive + "|prefix=" + str(prefix) + "|node=" + str(node),
    })


class _DAG:
    def __init__(self):
        self.nodes = []
        self.counter = 0

    def atom(self, op, value=None, **metadata):
        self.counter += 1
        ident = SCHEMA + ":d" + str(self.counter).zfill(6)
        node = {"id": ident, "namespace": SCHEMA, "op": op}
        if value is not None:
            node["value"] = value
        node.update(metadata)
        self.nodes.append(node)
        return ident

    def op(self, op, args=(), **metadata):
        self.counter += 1
        ident = SCHEMA + ":d" + str(self.counter).zfill(6)
        node = {"id": ident, "namespace": SCHEMA, "op": op, "args": list(args)}
        node.update(metadata)
        self.nodes.append(node)
        return ident


def _zero(dag):
    return dag.atom("exact_constant", "0")


def _one(dag):
    return dag.atom("exact_constant", "1")


def _const(dag, value):
    return dag.atom("exact_constant", value)


def _symbol(dag, name, *, domain=None):
    return dag.atom("exact_symbol", name, domain=domain)


def _add(dag, *args):
    return dag.op("add", args, exact=True)


def _sub(dag, a, b):
    return dag.op("sub", (a, b), exact=True)


def _mul(dag, a, b):
    return dag.op("mul", (a, b), exact=True)


def _div(dag, a, b):
    return dag.op("quotient", (a, b), exact=True)


def _neg(dag, a):
    return dag.op("neg", (a,), exact=True)


def _exp(dag, a):
    return dag.op("exp", (a,), exact=True)


def _pow(dag, a, exponent):
    exponent_node = _const(dag, exponent) if isinstance(exponent, str) else exponent
    return dag.op("power", (a, exponent_node), exact=True)


def _jet_const(dag, value):
    return [_const(dag, value), _zero(dag), _zero(dag)]


def _jet_add(dag, *jets):
    return [_add(dag, *(jet[i] for jet in jets)) for i in range(3)]


def _jet_sub(dag, a, b):
    return [_sub(dag, a[i], b[i]) for i in range(3)]


def _jet_scale(dag, a, scalar):
    return [_mul(dag, scalar, a[i]) for i in range(3)]


def _jet_mul(dag, a, b):
    return [
        _mul(dag, a[0], b[0]),
        _add(dag, _mul(dag, a[1], b[0]), _mul(dag, a[0], b[1])),
        _add(dag, _mul(dag, a[2], b[0]), _mul(dag, _mul(dag, _const(dag, "2"), a[1]), b[1]), _mul(dag, a[0], b[2])),
    ]


def _jet_inv(dag, a):
    a0, a1, a2 = a
    a0sq = _mul(dag, a0, a0)
    a0cu = _mul(dag, a0sq, a0)
    return [
        _div(dag, _one(dag), a0),
        _neg(dag, _div(dag, a1, a0sq)),
        _div(dag, _sub(dag, _mul(dag, _const(dag, "2"), _mul(dag, a1, a1)), _mul(dag, a0, a2)), a0cu),
    ]


def _jet_div(dag, a, b):
    return _jet_mul(dag, a, _jet_inv(dag, b))


def _jet_exp(dag, a):
    e0 = _exp(dag, a[0])
    return [e0, _mul(dag, e0, a[1]), _mul(dag, e0, _add(dag, a[2], _mul(dag, a[1], a[1])))]


def _jet_ports(dag, archive, prefix, values, label, gates):
    if not isinstance(values, list) or len(values) < 3:
        gates.append("source_port_jet_too_short:" + label)
        values = ["missing", "missing", "missing"]
    return [_port(dag, archive, prefix, values[i], label + "[" + str(i) + "]", gates) for i in range(3)]


def _jet_source(dag, source, label):
    """Lift a geometry source, whose Z independence is checked below."""
    return [source, _zero(dag), _zero(dag)]


def _chart_fields(dag, centered, zjet, pstar_jet, p0, swirl_normalized, vphysical):
    four_z = _jet_scale(dag, zjet, _const(dag, "4"))
    theta = _jet_add(dag, _jet_const(dag, "5/8"), centered["eh"])
    theta_z = _jet_add(dag, _jet_mul(dag, four_z, theta), centered["ek"])
    mean = _jet_add(dag, four_z, centered["em"])
    zz = _jet_mul(dag, zjet, zjet)
    axial = _jet_add(dag, _jet_scale(dag, zz, _const(dag, "16")), _jet_scale(dag, _jet_mul(dag, zjet, centered["em"]), _const(dag, "8")), centered["aa"])
    swirl = _jet_add(dag, _jet_const(dag, "5/6"), centered["eb"])
    pressure_shape = _jet_add(dag, _jet_const(dag, "5"), centered["ep"])
    e = _jet_sub(dag, _jet_div(dag, axial, _jet_mul(dag, pstar_jet, pstar_jet)), _jet_scale(dag, _jet_mul(dag, _jet_mul(dag, swirl_normalized, swirl_normalized), swirl), _const(dag, "1/2")))
    p = _jet_scale(dag, _jet_mul(dag, _jet_mul(dag, swirl_normalized, swirl_normalized), pressure_shape), _const(dag, "1/2"))
    return {
        "raw_shapes": {"mean": mean, "theta": theta, "theta_z": theta_z, "axial": axial, "swirl": swirl, "pressure_shape": pressure_shape},
        "normalized_fields": {
            "E": swirl_normalized,
            "V": _jet_div(dag, vphysical, pstar_jet),
            "m": _jet_div(dag, mean, pstar_jet),
            "h": _jet_mul(dag, swirl_normalized, theta),
            "k": _jet_div(dag, _jet_mul(dag, swirl_normalized, theta_z), pstar_jet),
            "e": e,
            "p": p,
            "P0": p0,
            "pressure": _jet_add(dag, p0, p),
        },
    }


def _integral(dag, rate, power, variable, sigma):
    s = _symbol(dag, "s", domain="[0,1]")
    sigma_s = dag.op("source_call", (sigma, s), exact=True, function="sigma_jets", jet_component=0, branch="lazy_flat_full_[0,1]")
    one_minus = _sub(dag, _one(dag), sigma_s)
    integrand = _mul(dag, _exp(dag, _neg(dag, _mul(dag, rate, _sub(dag, variable, s)))), _pow(dag, one_minus, power))
    return dag.op("integral", (integrand,), variable=s, bounds=[_zero(dag), variable], exact=True, unevaluated=True, quadrature=False, kernel="exp(-r*(t-s))*(1-sigma(s))^power", rate=rate, power=power)


def _first_three(values):
    return values[:3] if isinstance(values, list) and len(values) >= 3 else None


def _endpoint_functions(dag, chart, value):
    result = {}
    for group in ("centered", "raw_shapes", "normalized_fields"):
        result[group] = {
            name: [dag.op("function_substitution", (ref, chart["phase"], value),
                          ordinary_Z_order=i, radial_substitution_Z_independent=True)
                   for i, ref in enumerate(jet)] for name, jet in chart[group].items()
        }
    result["radius"] = [dag.op("function_substitution", (ref, chart["phase"], value),
                              ordinary_Z_order=i, radial_substitution_Z_independent=True)
                        for i, ref in enumerate(chart["radius"])]
    return result


def _definition_graph_checks(dag, charts, defects):
    seen = set()
    checks = {"unique_namespace_ids": True, "closed_topological_arguments": True,
              "rational_constants": True, "explicit_integral_binding": True,
              "all_field_and_defect_jets_have_three_rows": True}
    for node in dag.nodes:
        checks["unique_namespace_ids"] &= node["id"] not in seen and node["namespace"] == SCHEMA
        checks["closed_topological_arguments"] &= all(arg in seen for arg in node.get("args", []))
        if node["op"] == "exact_constant":
            try:
                Fraction(node["value"])
            except (ValueError, TypeError, ZeroDivisionError):
                checks["rational_constants"] = False
        if node["op"] == "integral":
            checks["explicit_integral_binding"] &= node.get("variable") in seen and all(
                endpoint in seen for endpoint in node.get("bounds", []))
        seen.add(node["id"])
    for chart in charts:
        for group in ("centered", "raw_shapes", "normalized_fields"):
            checks["all_field_and_defect_jets_have_three_rows"] &= all(
                isinstance(jet, list) and len(jet) == 3 and all(ref in seen for ref in jet)
                for jet in chart[group].values())
    checks["all_field_and_defect_jets_have_three_rows"] &= len(defects) == 5 and all(
        isinstance(jet, list) and len(jet) == 3 and all(ref in seen for ref in jet) for jet in defects)
    return checks


def derive(nodes, source_records, *, expected_nodes_digest):
    gates = []
    graph_digest = _digest(nodes)
    if len(nodes) != EXPECTED_NODE_COUNT:
        gates.append("node_count_mismatch")
    if graph_digest != EXPECTED_NODES_DIGEST or expected_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append("nodes_digest_mismatch")
    if not isinstance(source_records, Mapping):
        gates.append("source_records_not_mapping")
        source_records = {}

    texts = {name: _text(source_records, name, gates) for name in RECORD_HASHES}
    reports = {name: _json(source_records, name, gates) for name in (PARENT_REPORT, PRIOR_REPORT, R110_REPORT, REF_REPORT, SIGMA_REPORT)}
    r110 = reports[R110_REPORT]
    prior = reports[PRIOR_REPORT]
    c = r110.get("component_result", {}) if isinstance(r110, Mapping) else {}
    endpoint = _get(c, "endpoint_ports", {})
    params = _get(endpoint, "parameter_nodes", {})
    rsh = _get(endpoint, "Rsh", {})
    charts = _get(endpoint, "exact_chart_nodes", {})
    refchart = _get(charts, "reference", {})
    if not isinstance(endpoint, Mapping):
        gates.append("R110_endpoint_ports_missing")
        endpoint = {}
    if not isinstance(params, Mapping):
        gates.append("R110_parameter_nodes_missing")
        params = {}
    if not isinstance(rsh, Mapping):
        gates.append("R110_Rsh_ports_missing")
        rsh = {}
    if not isinstance(refchart, Mapping):
        gates.append("R110_reference_chart_ports_missing")
        refchart = {}

    try:
        native_archive = json.loads(texts[C6_ARCHIVE])
        native_prefix = len(native_archive.get("exact_graph_nodes", []))
    except Exception:
        native_archive, native_prefix = {}, 0
        gates.append("native_C6_archive_unreadable")
    try:
        true_archive = json.loads(texts[TRUE_C6_ARCHIVE])
        true_prefix = len(true_archive.get("exact_graph_nodes", []))
    except Exception:
        true_archive, true_prefix = {}, 0
        gates.append("true_C6_archive_unreadable")

    prior_component = prior.get("component_result", {}) if isinstance(prior, Mapping) else {}
    audit = prior_component.get("c6_archive_port_audit", {}) if isinstance(prior_component, Mapping) else {}
    audit_ok = audit.get("valid") is True and audit.get("accepted_C6_Rm_inlet_identity") is True
    r110_ok = all(c.get(key) is True for key in (
        "passed", "full_R110_Rm_background_projection_proved", "full_R110_Rm_own_history_projection_proved")) and (
        c.get("nodes_digest") == EXPECTED_NODES_DIGEST and c.get("ordinary_Z_orders") == [0, 1, 2]
        and c.get("numeric_values_computed") is False and c.get("old_operators_replayed") is False)
    if not audit_ok:
        gates.append("saved_C6_archive_port_audit_not_accepted")
    if not r110_ok:
        gates.append("saved_R110_own_history_not_accepted")

    source_ast = {
        "native_C6": _ast_source(source_records, C6_SOURCE, gates),
        "reference_restore": _ast_source(source_records, REF_SOURCE, gates),
        "flat_sigma": _ast_source(source_records, SIGMA_SOURCE, gates),
    }
    witness = {
        "native_C6_reference_measure_is_exact_G": _source_has(texts[C6_SOURCE], "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))", "self.reference_length=length=g.sub(g.sub(logref,T),c(8))", "physical_measure=jacobian"),
        "native_C6_actual_U_is_physical": _source_has(texts[C6_SOURCE], "actual=first_switch_fields", "actual", "self.Pstar"),
        "native_C6_Pstar_source_binding": _source_has(texts[C6_SOURCE], "self.Pstar=Pstar=bridge_owner.inlet.Pstar"),
        "native_C6_reference_restore_postrestore_sources": _source_has(texts[C6_SOURCE], "self.add_chart('reference'", "self.add_chart('restoration'", "self.add_chart('post_restore'"),
        "reference_centered_source_AST": _source_has(texts[REF_SOURCE], "mean_error=moments['mean']-z*4", "angular_error=moments['theta']-c.mpf(5)/8", "mixed_error=moments['theta_z']-(z*moments['theta'])*4", "axial_square=moments['axial']-(z*moments['mean'])*8+square(z)*16", "swirl_error=moments['swirl']-c.mpf(5)/6", "pressure_error=moments['pressure']-5"),
        "reference_formula_source_AST": _source_has(texts[REF_SOURCE], "reference_centered", "restore_centered", "sigma_jets", "actual_Rsh_centered_inlet"),
        "reference_terminal_source_AST": _source_has(texts[REF_SOURCE], "def postrestore", "restore_exit_to_Rm", "exact_reference_offset"),
        "flat_sigma_jets_source_AST": "sigma_jets" in _names(source_ast["flat_sigma"]),
        "flat_sigma_full_lazy_branch": _source_has(texts[SIGMA_SOURCE], "def sigma_jets", "0<=x<=1") or _source_has(texts[SIGMA_SOURCE], "def sigma_jets", "_sigma_left"),
    }
    source_ast["all_required"] = all(witness.values())
    if not source_ast["all_required"]:
        gates.append("source_AST_witness_incomplete")

    dag = _DAG()
    z = _symbol(dag, "Z", domain="[-1,1]")
    q = _symbol(dag, "q", domain="[0,1]")
    t = _symbol(dag, "t", domain="[0,1]")
    p = _symbol(dag, "p", domain="[0,1]")
    native_prefix_ok = native_prefix == 97295 and nodes[:native_prefix] == native_archive.get("exact_graph_nodes")
    true_prefix_ok = true_prefix == 196042 and nodes[:true_prefix] == true_archive.get("exact_graph_nodes")
    if not native_prefix_ok:
        gates.append("native_C6_prefix_mismatch")
    if not true_prefix_ok:
        gates.append("true_C6_prefix_mismatch")

    actual = _get(rsh, "actual", {})
    moments = _get(rsh, "moments", {})
    fields = _get(rsh, "fields", {})
    rsh_radius = _port(dag, C6_ARCHIVE, native_prefix, rsh.get("radius"), "Rsh.radius", gates)
    phi = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(actual, "phi"), "Rsh.actual.phi", gates)
    actual_u = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(actual, "U"), "Rsh.actual.U_physical", gates)
    H = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(moments, "H"), "Rsh.moments.H", gates)
    M = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(moments, "M"), "Rsh.moments.M", gates)
    K = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(moments, "K"), "Rsh.moments.K", gates)
    A = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(moments, "A"), "Rsh.moments.A", gates)
    B = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(moments, "B"), "Rsh.moments.B", gates)
    C = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(moments, "C"), "Rsh.moments.C", gates)
    swirl_input = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(fields, "E"), "Rsh.fields.E_native_swirl", gates)
    P0 = _jet_ports(dag, C6_ARCHIVE, native_prefix, _get(fields, "P0"), "Rsh.fields.P0", gates)
    zjet = [z, _one(dag), _zero(dag)]
    four_z = _jet_scale(dag, zjet, _const(dag, "4"))
    em = _jet_sub(dag, M, four_z)
    two_phi = _jet_scale(dag, phi, _const(dag, "2"))
    eh = _jet_sub(dag, _jet_div(dag, H, two_phi), _jet_const(dag, "5/8"))
    zH = _jet_mul(dag, _jet_const(dag, "4"), _jet_mul(dag, zjet, _jet_div(dag, H, two_phi)))
    ek = _jet_sub(dag, _jet_div(dag, K, two_phi), zH)
    aa = _jet_add(dag, _jet_sub(dag, A, _jet_scale(dag, _jet_mul(dag, zjet, M), _const(dag, "8"))), _jet_scale(dag, _jet_mul(dag, zjet, zjet), _const(dag, "16")))
    eb = _jet_sub(dag, _jet_div(dag, B, _jet_mul(dag, phi, phi)), _jet_const(dag, "5/6"))
    ep = _jet_sub(dag, _jet_div(dag, C, _jet_mul(dag, phi, phi)), _jet_const(dag, "5"))
    Eaxial = _jet_sub(dag, actual_u, four_z)

    pstar_node = params.get("Pstar")
    T_node = params.get("T")
    logref_node = params.get("logref")
    ref_measure_node = refchart.get("physical_measure", params.get("reference_length"))
    pstar = _port(dag, C6_ARCHIVE, native_prefix, pstar_node, "native_C6.parameter.Pstar_positive_Z_independent", gates, role="parameter")
    T = _port(dag, C6_ARCHIVE, native_prefix, T_node, "native_C6.parameter.T", gates, role="parameter")
    logref = _port(dag, C6_ARCHIVE, native_prefix, logref_node, "native_C6.parameter.logref_10_logC_plus_logP", gates, role="parameter")
    reference_jacobian = _port(dag, C6_ARCHIVE, native_prefix, ref_measure_node, "native_C6.reference.physicalJacobian_already_G", gates, role="reference_physical_measure")
    G = reference_jacobian
    dag.nodes[-1]["definition"] = "G = physicalJacobian/reference_length = 10*(logC+logP)-T-8; node 81519 is already this exact native reference length"
    Gjet = _jet_source(dag, G, "native_C6.reference.G")
    sigma = _port(dag, SIGMA_SOURCE, "AST:sigma_jets", "sigma_jets", "flat_pulse_derivatives.sigma_jets", gates, role="source_AST_function")
    dag.nodes[-1]["source_function_definition"] = {
        "function": "sigma_jets",
        "domain": "[0,1]",
        "analytic_branch": "original flat-pulse piecewise analytic source: left plateau, logistic middle, reflected middle, right plateau",
        "lazy_branch": "full_[0,1]",
        "interval_enclosure_callable_not_defining_values": True,
        "arbitrary_cutoff": False,
        "exact_definition": "sigma(x)=0 for x<=0; exp(q)/(1+exp(q)) for 0<x<=1/2, q=(1-x)^(-2)-x^(-2); 1-sigma(1-x) for 1/2<x<1; 1 for x>=1",
        "source_call_denotes_underlying_analytic_function": True,
    }
    E = Eaxial
    E_swirl = swirl_input
    pstar_jet = [pstar, _zero(dag), _zero(dag)]
    pstar_sq = _jet_mul(dag, pstar_jet, pstar_jet)
    exp_minus_three_fifths = _exp(dag, _const(dag, "-3/5"))
    exp_minus_three_fifths_jet = [exp_minus_three_fifths, _zero(dag), _zero(dag)]
    Am = _jet_div(dag, exp_minus_three_fifths_jet, _jet_add(dag, _jet_const(dag, "1"), _jet_mul(dag, zjet, zjet)))
    alpha = _jet_inv(dag, _jet_mul(dag, pstar_sq, _jet_mul(dag, Am, Am)))
    initial = {"em": em, "eh": eh, "ek": ek, "aa": aa, "eb": eb, "ep": ep}
    qjet = [q, _zero(dag), _zero(dag)]
    gjet = _jet_mul(dag, qjet, Gjet)
    def decay(rate, argument):
        return _exp(dag, _neg(dag, _mul(dag, _const(dag, rate), argument)))
    def jet_decay(rate, argument):
        return _jet_exp(dag, _jet_scale(dag, argument, _const(dag, "-" + rate)))
    d1 = jet_decay("1", gjet)
    d8 = jet_decay("8/5", gjet)
    d6 = jet_decay("6/5", gjet)
    d2 = jet_decay("1/5", gjet)
    ref_centered = {
        "em": _jet_add(dag, E, _jet_mul(dag, _jet_sub(dag, initial["em"], E), d1)),
        "eh": _jet_mul(dag, initial["eh"], d8),
        "ek": _jet_add(dag, _jet_mul(dag, initial["ek"], d8), _jet_div(dag, _jet_mul(dag, E, _jet_sub(dag, _jet_const(dag, "1"), d8)), _jet_const(dag, "8/5"))),
        "aa": _jet_add(dag, _jet_mul(dag, E, E), _jet_mul(dag, _jet_sub(dag, initial["aa"], _jet_mul(dag, E, E)), d1)),
        "eb": _jet_mul(dag, initial["eb"], d6),
        "ep": _jet_mul(dag, initial["ep"], d2),
    }
    vref = _jet_add(dag, four_z, E)
    one_minus_q = _jet_sub(dag, _jet_const(dag, "1"), qjet)
    swirl_exponent = _jet_scale(dag, _jet_sub(dag, _jet_const(dag, "-8"), _jet_mul(dag, one_minus_q, Gjet)), _const(dag, "1/10"))
    swref = _jet_mul(dag, swirl_input, _jet_exp(dag, _jet_scale(dag, gjet, _const(dag, "1/10"))))
    rsh_radius_jet = _jet_source(dag, rsh_radius, "native_C6.Rsh.radius")
    rref = _jet_mul(dag, rsh_radius_jet, _jet_exp(dag, gjet))
    reference = _chart_fields(dag, ref_centered, zjet, pstar_jet, P0, swref, vref)
    reference["centered"] = ref_centered
    reference["G"] = G
    reference["radius"] = rref
    reference["phase"] = q
    reference["physical_V"] = vref
    reference["normalized_swirl"] = swref
    reference["source_inputs"] = {"Eaxial_actual_U_minus_4Z": E, "native_swirl_E_port": swirl_input, "P0": P0}

    g_z = Gjet
    d1_z = jet_decay("1", g_z)
    d8_z = jet_decay("8/5", g_z)
    d6_z = jet_decay("6/5", g_z)
    d2_z = jet_decay("1/5", g_z)
    reference_at_q1 = {
        "em": _jet_add(dag, E, _jet_mul(dag, _jet_sub(dag, initial["em"], E), d1_z)),
        "eh": _jet_mul(dag, initial["eh"], d8_z),
        "ek": _jet_add(dag, _jet_mul(dag, initial["ek"], d8_z), _jet_div(dag, _jet_mul(dag, E, _jet_sub(dag, _jet_const(dag, "1"), d8_z)), _jet_const(dag, "8/5"))),
        "aa": _jet_add(dag, _jet_mul(dag, E, E), _jet_mul(dag, _jet_sub(dag, initial["aa"], _jet_mul(dag, E, E)), d1_z)),
        "eb": _jet_mul(dag, initial["eb"], d6_z),
        "ep": _jet_mul(dag, initial["ep"], d2_z),
    }
    reference["centered_at_q1"] = reference_at_q1

    def K_kernel(rate, power, variable=t):
        return _integral(dag, _const(dag, rate), power, variable, sigma)
    K1 = K_kernel("1", "1")
    K85 = K_kernel("8/5", "1")
    Ksq = K_kernel("1", "2")
    tjet = [t, _zero(dag), _zero(dag)]
    d1t = jet_decay("1", tjet)
    d8t = jet_decay("8/5", tjet)
    d6t = jet_decay("6/5", tjet)
    d2t = jet_decay("1/5", tjet)
    sigma_t = dag.op("source_call", (sigma, t), exact=True, function="sigma_jets", jet_component=0, branch="lazy_flat_full_[0,1]")
    restore_centered = {
        "em": _jet_add(dag, _jet_mul(dag, reference_at_q1["em"], d1t), _jet_scale(dag, E, K1)),
        "eh": _jet_mul(dag, reference_at_q1["eh"], d8t),
        "ek": _jet_add(dag, _jet_mul(dag, reference_at_q1["ek"], d8t), _jet_scale(dag, E, K85)),
        "aa": _jet_add(dag, _jet_mul(dag, reference_at_q1["aa"], d1t), _jet_scale(dag, _jet_mul(dag, E, E), Ksq)),
        "eb": _jet_mul(dag, reference_at_q1["eb"], d6t),
        "ep": _jet_mul(dag, reference_at_q1["ep"], d2t),
    }
    sigma_t_jet = [sigma_t, _zero(dag), _zero(dag)]
    vrestore = _jet_add(dag, four_z, _jet_mul(dag, E, _jet_sub(dag, _jet_const(dag, "1"), sigma_t_jet)))
    swrestore = _jet_div(dag, _jet_exp(dag, _jet_scale(dag, _jet_add(dag, _jet_const(dag, "-8"), tjet), _const(dag, "1/10"))), _jet_add(dag, _jet_const(dag, "1"), _jet_mul(dag, zjet, zjet)))
    rz_port = _port(dag, C6_ARCHIVE, native_prefix, _get(endpoint, "Rz.radius"), "native_C6.endpoint.Rz.radius", gates, role="endpoint_radius")
    rrestore = _jet_mul(dag, _jet_source(dag, rz_port, "native_C6.Rz.radius"), _jet_exp(dag, tjet))
    restoration = _chart_fields(dag, restore_centered, zjet, pstar_jet, P0, swrestore, vrestore)
    restoration.update({"centered": restore_centered, "radius": rrestore, "phase": t, "physical_V": vrestore, "normalized_swirl": swrestore, "kernels": {"mean": K1, "mixed": K85, "square": Ksq}, "sigma": sigma_t})

    rates = {"em": "1", "eh": "8/5", "ek": "8/5", "aa": "1", "eb": "6/5", "ep": "1/5"}
    t_one = _const(dag, "1")
    t_one_jet = [t_one, _zero(dag), _zero(dag)]
    K1_one = K_kernel("1", "1", t_one)
    K85_one = K_kernel("8/5", "1", t_one)
    Ksq_one = K_kernel("1", "2", t_one)
    restore_at_t1 = {
        "em": _jet_add(dag, _jet_mul(dag, reference_at_q1["em"], jet_decay("1", t_one_jet)), _jet_scale(dag, E, K1_one)),
        "eh": _jet_mul(dag, reference_at_q1["eh"], jet_decay("8/5", t_one_jet)),
        "ek": _jet_add(dag, _jet_mul(dag, reference_at_q1["ek"], jet_decay("8/5", t_one_jet)), _jet_scale(dag, E, K85_one)),
        "aa": _jet_add(dag, _jet_mul(dag, reference_at_q1["aa"], jet_decay("1", t_one_jet)), _jet_scale(dag, _jet_mul(dag, E, E), Ksq_one)),
        "eb": _jet_mul(dag, reference_at_q1["eb"], jet_decay("6/5", t_one_jet)),
        "ep": _jet_mul(dag, reference_at_q1["ep"], jet_decay("1/5", t_one_jet)),
    }
    post_pjet = [p, _zero(dag), _zero(dag)]
    post_centered = {name: _jet_mul(dag, restore_at_t1[name], jet_decay(rate, post_pjet)) for name, rate in rates.items()}
    p_one_jet = [_const(dag, "1"), _zero(dag), _zero(dag)]
    post_at_p1 = {name: _jet_mul(dag, restore_at_t1[name], jet_decay(rate, p_one_jet)) for name, rate in rates.items()}
    vpost = four_z
    swpost = _jet_div(dag, _jet_exp(dag, _jet_scale(dag, _jet_add(dag, _jet_const(dag, "-7"), post_pjet), _const(dag, "1/10"))), _jet_add(dag, _jet_const(dag, "1"), _jet_mul(dag, zjet, zjet)))
    rrestore_port = _port(dag, C6_ARCHIVE, native_prefix, _get(endpoint, "Rrestore.radius"), "native_C6.endpoint.Rrestore.radius", gates, role="endpoint_radius")
    rpost = _jet_mul(dag, _jet_source(dag, rrestore_port, "native_C6.Rrestore.radius"), _jet_exp(dag, post_pjet))
    postrestore = _chart_fields(dag, post_centered, zjet, pstar_jet, P0, swpost, vpost)
    postrestore.update({"centered": post_centered, "radius": rpost, "phase": p, "physical_V": vpost, "normalized_swirl": swpost})
    restoration["centered_at_t1"] = restore_at_t1
    postrestore["centered_at_p1"] = post_at_p1

    endpoint_substitutions = {
        name + suffix: _endpoint_functions(dag, chart, _const(dag, value))
        for name, chart, suffix, value in (
            ("reference", reference, "_q0", "0"), ("reference", reference, "_q1", "1"),
            ("restoration", restoration, "_t0", "0"), ("restoration", restoration, "_t1", "1"),
            ("postrestore", postrestore, "_p0", "0"), ("postrestore", postrestore, "_p1", "1"))
    }
    reference["endpoint_substitution_nodes"] = {"q=0": endpoint_substitutions["reference_q0"], "q=1": endpoint_substitutions["reference_q1"]}
    restoration["endpoint_substitution_nodes"] = {"t=0": endpoint_substitutions["restoration_t0"], "t=1": endpoint_substitutions["restoration_t1"]}
    postrestore["endpoint_substitution_nodes"] = {"p=0": endpoint_substitutions["postrestore_p0"], "p=1": endpoint_substitutions["postrestore_p1"]}

    rm_defect = [post_at_p1["em"], post_at_p1["ek"], post_at_p1["eh"], _jet_sub(dag, _jet_mul(dag, alpha, post_at_p1["aa"]), _jet_scale(dag, post_at_p1["eb"], _const(dag, "1/2"))), _jet_scale(dag, post_at_p1["ep"], _const(dag, "1/2"))]
    endpoint_ports = {
        "Rsh": {"archive": C6_ARCHIVE, "prefix": native_prefix, "radius": rsh.get("radius"), "actual": actual, "moments": moments, "fields": fields, "qualified": True},
        "parameter": {"Pstar": pstar_node, "T": T_node, "logref": logref_node, "reference_length": ref_measure_node, "archive": C6_ARCHIVE, "prefix": native_prefix},
        "Rz": {"source": "native_C6.endpoint.restoration", "qualified": True},
        "Rrestore": {"source": "native_C6.endpoint.post_restore", "qualified": True},
    }
    native_rsh_binding = _get(native_archive, "prescribed_R100_Rm_C6_bindings.endpoints.Rsh", {})
    native_rsh_binding_exact = isinstance(native_rsh_binding, Mapping) and native_rsh_binding == rsh
    if not native_rsh_binding_exact:
        gates.append("Rsh_native_archive_binding_mismatch")
    prior_qualified = prior_component.get("c6_archive_port_audit", {}).get("qualified_handles", {}) if isinstance(prior_component, Mapping) else {}
    same_binding = bool(audit_ok and r110_ok and native_rsh_binding_exact and native_prefix_ok and true_prefix_ok
                        and _get(c, "endpoint_ports.endpoint_substitution_identity.Rsh") is True
                        and isinstance(prior_qualified, Mapping))
    if not same_binding:
        gates.append("same_initial_Rsh_source_port_binding_unproved")
    source_premises = {
        "saved_R110_background_and_own_history": r110_ok,
        "saved_C6_port_audit": audit_ok,
        "same_initial_source_port_binding": same_binding,
        "native_prefix_is_saved_prefix": native_prefix_ok,
        "Rsh_matches_native_archive_endpoint": native_rsh_binding_exact,
        "true_C6_prefix_is_saved_prefix": true_prefix_ok,
        "source_AST_witness_all": source_ast["all_required"],
        "native_reference_node_is_already_G": witness["native_C6_reference_measure_is_exact_G"],
        "physical_U_distinct_from_swirl_E": True,
        "Pstar_positive_Z_independent_saved_parameter": pstar_node == 33572 and witness["native_C6_Pstar_source_binding"] and _get(c, "domain_evidence.pstar_exp_logP_recipe") is True,
        "geometry_Z_independent_native_parameters": T_node == 7 and logref_node == 81514 and ref_measure_node == 81519,
        "G_positive_common_domain_saved": _get(c, "domain_evidence.finite_positive_reference_length") is True,
        "positive_R110_phi_saved": _get(c, "domain_evidence.positive_R110_phi") is True,
        "full_phi_history_normalization_saved": _get(c, "domain_evidence.full_phi_history_normalization") is True,
        "phi_y_chart_qualification_saved": _get(c, "domain_evidence.phi_y_chart_qualification") is True,
        "six_history_rate_conjugation_saved": _get(c, "domain_evidence.six_history_rate_conjugation") is True,
        "common_Z_source_interval_saved": _get(c, "domain_evidence.common_Z_source_interval") is True,
        "pstar_exp_logP_recipe_saved": _get(c, "domain_evidence.pstar_exp_logP_recipe") is True,
        "lazy_flat_sigma_source": witness["flat_sigma_full_lazy_branch"],
        "saved_normalized_shape_definition": _get(c, "source_normalization_theorem.normalized_history_shapes") == {
            "theta": "H/(2phi)", "theta_z": "K/(2phi)", "pressure": "C/phi^2",
            "swirl": "B/phi^2", "mean": "M", "axial": "A"},
        "saved_centered_transport_source_theorem": all(_get(c, "ast_witnesses." + name) is True for name in (
            "six_history_rates_nonreset_and_conjugated", "smooth_sigma_and_triangular_fields", "c6_U_alias_and_mean_U_minus_mean_source")),
    }
    definitions = {
        "reference": reference,
        "restoration": restoration,
        "postrestore": postrestore,
        "rm_p1_defect_jet": rm_defect,
        "endpoint_substitutions": endpoint_substitutions,
        "ordinary_Z_jet_rules": {
            "orders": [0, 1, 2],
            "product_second": "a2*b0 + 2*a1*b1 + a0*b2",
            "reciprocal": "(1/a0, -a1/a0^2, (2*a1^2-a0*a2)/a0^3)",
            "exp": "(exp(a0), exp(a0)*a1, exp(a0)*(a2+a1^2))",
            "spatial_variable": "Z",
            "radial_variables": "q/t/p are native radial chart coordinates, independent of Z; they are not physical time",
            "temporal_recursion_claimed": False,
        },
    }
    graph_checks = _definition_graph_checks(dag, (reference, restoration, postrestore), rm_defect)
    for name, ok in graph_checks.items():
        if not ok:
            gates.append("definition_graph_check_failed:" + name)
    adapter_defined = bool(source_ast["all_required"] and native_prefix_ok and true_prefix_ok and all(graph_checks.values()) and not gates)
    identity_proved = bool(adapter_defined and all(source_premises.values()) and len(rm_defect) == 5)
    result = {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(adapter_defined and identity_proved and not gates),
        "nodes_digest": graph_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes),
        "native_C2_reference_adapter_defined": adapter_defined,
        "native_C2_reference_Rm_identity_proved": identity_proved,
        "legacy_feedback_reference_alias_proved": False,
        "definition_graph": {"namespace": SCHEMA, "active_graph_digest": EXPECTED_NODES_DIGEST, "nodes": dag.nodes, "node_count": len(dag.nodes), "outputs": {"reference": reference, "restoration": restoration, "postrestore": postrestore, "rm_p1_defect_jet": rm_defect}},
        "function_definitions": definitions,
        "qualified_endpoint_ports": endpoint_ports,
        "same_native_C6_Rsh_source_port_binding": same_binding,
        "source_ast_witness": witness,
        "definition_graph_checks": graph_checks,
        "source_prerequisite_audit": source_premises,
        "source_identity_argument": {
            "new_owner_only": True, "legacy_feedback_owner_identity_assumed": False,
            "initial_data": "The new owner consumes the exact same qualified Rsh phi,U,H,M,K,A,B,C,P0 ports.",
            "centered_ODE": {"em": "E_axial*a-em", "eh": "-8*eh/5", "ek": "E_axial*a-8*ek/5",
                             "aa": "E_axial^2*a^2-aa", "eb": "-6*eb/5", "ep": "-ep/5"},
            "axial_cutoff_a_by_chart": {"reference": "1", "restoration": "1-sigma(t)", "postrestore": "0"},
            "uniqueness": "The same initial functions and continuous linear centered equations on finite native intervals identify the full functions. Smooth parameter dependence gives ordinary Z0/1/2. All positivity and source-domain premises remain conditional as in the accepted theorem.",
            "numeric_equality_or_interval_overlap_used": False,
        },
        "six_rate_comparison": {"reference": {"em": "1", "eh": "8/5", "ek": "8/5", "aa": "1", "eb": "6/5", "ep": "1/5"}, "restoration": rates, "post_restore": rates, "source_formula": "same centered linear ODE normalization and uniqueness"},
        "radial_map_comparison": {"reference": "Rsh*exp(q*G)", "restoration": "Rz*exp(t)", "post_restore": "Rrestore*exp(p)", "offsets": {"reference": "-8", "restoration": "-8+t", "post_restore": "-7+p"}},
        "initial_source_port_binding": {"Rsh": endpoint_ports["Rsh"], "native_C6_same_construction": same_binding, "cross_archive_comparison": "qualified provenance relation only; no naked node comparison"},
        "sourceaudit1new": 1,
        "sourceaudit0ancestor": 0,
        "source_audit_1_new": 1,
        "source_audit_0_ancestor": 0,
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z_not_temporal_n": True,
        "source_only_true": True,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "patch_finite_N_source_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "time_order_recursion_claimed": False,
        "full_goal_complete": False,
        "runtime_modelmetadata": None,
        "runtime_model_metadata": None,
        "record_hashes": RECORD_HASHES,
        "decoded_record_hashes": DECODED_RECORD_HASHES,
        "source_inspection_audit": {
            "source_records_only": True,
            "source_AST_graph_saved_data_and_new_definitions_only": True,
            "nodes_hashed_only": False,
            "original_graph_nodes_inspected_not_evaluated": True,
            "new_function_definition_graph_constructed": True,
            "definition_graph_saved": True,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "producer_provider_builder_checker_replayed": False,
            "quadrature_or_norm_calls": 0,
            "numeric_field_evaluated": False,
            "numeric_values_computed": False,
            "old_operators_replayed": False,
            "legacy_feedback_owner_alias_assumed": False,
            "runtime_model_metadata": None,
        },
        "unresolved_gates": list(dict.fromkeys(gates)),
    }
    return result
