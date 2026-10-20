"""Source-only native-input C2 leading patch composition.

This consumer binds the accepted native C2 reference graph and the accepted
true Rm C6 source ports as qualified source atoms.  It constructs a new,
unevaluated expression DAG for the leading Rm patch.  No producer, solver,
operator, quadrature, norm, or project module is imported or executed.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction


SCHEMA = "original_native_C2_patch_composition_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
ADAPTER_REPORT = "lei_ren_part1_paper_compliant_current_original_native_C2_reference_adapter.json"
ADAPTER_SHA256 = "3fd6359a456345f15f3323c594558c9edb89169a21114837c160be3ef86e3db7"
PRIOR_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_patch_C2_C6_source_projection.json"
PRIOR_SHA256 = "190d8b54bec0b9cf8e47d6db979fe68eaaf6594e28471f9da4ef86a8c22e88c9"

C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py"
C6_SOURCE_SHA256 = "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002"
C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.json.gz"
C6_ARCHIVE_SHA256 = "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"
C6_ARCHIVE_DECODED_SHA256 = "b5671cc0c60821db9a5c591ea899c9616c3cc949cec4e9214cc57a31ae6cdc10"
C6_CHECK = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions_check.json"
C6_CHECK_SHA256 = "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7"

RECORD_HASHES = {
    ADAPTER_REPORT: ADAPTER_SHA256,
    PRIOR_REPORT: PRIOR_SHA256,
    C6_SOURCE: C6_SOURCE_SHA256,
    C6_ARCHIVE: C6_ARCHIVE_SHA256,
    C6_CHECK: C6_CHECK_SHA256,
}
DECODED_RECORD_HASHES = {C6_ARCHIVE: C6_ARCHIVE_DECODED_SHA256}

CENTERS = ("5/4", "3/2", "7/4")
RADIUS = "1/40"
RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
LEADING_HISTORY_NAMES = ("m", "h", "k", "e", "p")
RAW_INLET_NAMES = ("M", "H", "K", "A", "B", "C")


def _digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("utf-8")).hexdigest()


def _get(value, path, default=None):
    for part in path.split(".") if path else ():
        if isinstance(value, Mapping):
            if part not in value:
                return default
            value = value[part]
        elif isinstance(value, list) and part.isdigit():
            index = int(part)
            if index >= len(value):
                return default
            value = value[index]
        else:
            return default
    return value


def _record(records, name, gates):
    row = records.get(name) if isinstance(records, Mapping) else None
    if not isinstance(row, Mapping):
        gates.append("missing_source_record:" + name)
        return ""
    if row.get("path") != "experiments/root_st073/" + name:
        gates.append("source_path_mismatch:" + name)
    if row.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name)
    source = row.get("source_text")
    if not isinstance(source, str):
        gates.append("source_text_missing:" + name)
        return ""
    if name in DECODED_RECORD_HASHES:
        if row.get("encoding") != "gzip" or row.get("decoded_sha256") != DECODED_RECORD_HASHES[name]:
            gates.append("encoded_archive_contract_mismatch:" + name)
            return ""
        try:
            # Authenticate the raw archive before decompression or JSON parsing.
            raw = base64.b64decode(row.get("raw_bytes_base64", ""), validate=True)
        except Exception as error:
            gates.append("archive_raw_base64_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(raw).hexdigest() != RECORD_HASHES[name]:
            gates.append("archive_raw_hash_mismatch:" + name)
            return ""
        try:
            decoded = gzip.decompress(raw)
        except Exception as error:
            gates.append("archive_decode_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(decoded).hexdigest() != DECODED_RECORD_HASHES[name]:
            gates.append("archive_decoded_hash_mismatch:" + name)
        if decoded != source.encode("utf-8"):
            gates.append("archive_source_text_mismatch:" + name)
        return source
    if row.get("encoding") is not None or row.get("raw_bytes_base64") is not None:
        gates.append("unexpected_encoded_record:" + name)
    if hashlib.sha256(source.encode("utf-8")).hexdigest() != RECORD_HASHES[name]:
        gates.append("source_hash_mismatch:" + name)
        return ""
    return source


def _json(text, name, gates):
    if not text:
        return {}
    try:
        value = json.loads(text)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}
    return value if isinstance(value, Mapping) else {}


def _tree(text, name, gates):
    if not text:
        return None
    try:
        return ast.parse(text, filename=name)
    except SyntaxError as error:
        gates.append("source_ast_parse_error:" + name + ":" + str(error.lineno))
        return None


def _source_has(text, *needles):
    return all(needle in text for needle in needles)


def _function_names(tree):
    if tree is None:
        return set()
    return {node.name for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


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
    return dag.atom("exact_constant", str(value))


def _symbol(dag, name, domain=None):
    return dag.atom("exact_symbol", name, domain=domain)


def _add(dag, *args):
    return dag.op("add", args, exact=True)


def _sub(dag, a, b):
    return dag.op("sub", (a, b), exact=True)


def _mul(dag, *args):
    return dag.op("mul", args, exact=True)


def _div(dag, a, b):
    return dag.op("quotient", (a, b), exact=True)


def _neg(dag, a):
    return dag.op("neg", (a,), exact=True)


def _exp(dag, a):
    return dag.op("exp", (a,), exact=True)


def _log(dag, a):
    return dag.op("log", (a,), exact=True, positive_argument=True)


def _pow(dag, a, exponent):
    return dag.op("power", (a, _const(dag, exponent)), exact=True,
                  positive_base=True, exponent=str(exponent))


def _source_atom(dag, *, report_sha, namespace, node, label, source_archive,
                 role, path=None, **metadata):
    value = {
        "source_report_sha256": report_sha,
        "source_namespace": namespace,
        "source_node": node,
        "source_archive": source_archive,
        "archive_sha256": RECORD_HASHES.get(source_archive),
        "label": label,
        "role": role,
        "qualified": (SCHEMA + "|report=" + report_sha + "|namespace=" + namespace
                       + "|node=" + str(node) + "|label=" + label),
    }
    if path is not None:
        value["source_path"] = path
    value.update(metadata)
    return dag.atom("source_function_port", value)


def _saved_atom(dag, graph_ids, node, label, path, report_sha=ADAPTER_SHA256,
                namespace="original_native_C2_reference_adapter_v1"):
    valid = isinstance(node, str) and node in graph_ids
    if not valid:
        node = "missing"
    return _source_atom(dag, report_sha=report_sha, namespace=namespace, node=node,
                        label=label, source_archive=ADAPTER_REPORT,
                        role="saved_definition", path=path,
                        source_graph_digest=EXPECTED_NODES_DIGEST,
                        source_graph_node_valid=valid,
                        re_evaluated=False)


def _c6_atom(dag, archive_nodes, node, label, path, role="C6_source_port",
             bound_variable=None, archive=C6_ARCHIVE):
    valid = type(node) is int and 0 <= node < len(archive_nodes)
    if not valid:
        node = "missing"
    return _source_atom(dag, report_sha=PRIOR_SHA256,
                        namespace="original_Rm_patch_C2_C6_source_projection_v1",
                        node=node, label=label, source_archive=archive, role=role,
                        path=path, source_graph_prefix=196042,
                        source_node_valid=valid, bound_variable=bound_variable,
                        re_evaluated=False)


def _ast_function_atom(dag, label, function):
    return _source_atom(
        dag, report_sha=PRIOR_SHA256,
        namespace="original_Rm_patch_C2_C6_source_projection_v1",
        node=function, label=label, source_archive=C6_SOURCE,
        role="source_AST_function", path="function:" + function,
        re_evaluated=False)


def _jet_source(dag, source, label):
    return [source, _zero(dag), _zero(dag)]


def _jet_const(dag, value):
    return [_const(dag, value), _zero(dag), _zero(dag)]


def _jet_add(dag, *jets):
    return [_add(dag, *(jet[i] for jet in jets)) for i in range(3)]


def _jet_sub(dag, a, b):
    return [_sub(dag, a[i], b[i]) for i in range(3)]


def _jet_scale(dag, a, scalar):
    if isinstance(scalar, list):
        return _jet_mul(dag, a, scalar)
    return [_mul(dag, scalar, a[i]) for i in range(3)]


def _jet_mul(dag, a, b):
    return [
        _mul(dag, a[0], b[0]),
        _add(dag, _mul(dag, a[1], b[0]), _mul(dag, a[0], b[1])),
        _add(dag, _mul(dag, a[2], b[0]),
             _mul(dag, _const(dag, "2"), a[1], b[1]),
             _mul(dag, a[0], b[2])),
    ]


def _jet_inv(dag, a):
    a0, a1, a2 = a
    a0sq = _mul(dag, a0, a0)
    a0cu = _mul(dag, a0sq, a0)
    return [
        _div(dag, _one(dag), a0),
        _neg(dag, _div(dag, a1, a0sq)),
        _div(dag, _sub(dag, _mul(dag, _const(dag, "2"), a1, a1),
                       _mul(dag, a0, a2)), a0cu),
    ]


def _jet_div(dag, a, b):
    return _jet_mul(dag, a, _jet_inv(dag, b))


def _jet_exp(dag, a):
    e0 = _exp(dag, a[0])
    return [e0, _mul(dag, e0, a[1]),
            _mul(dag, e0, _add(dag, a[2], _mul(dag, a[1], a[1])))]


def _jet_sqrt(dag, a):
    s0 = dag.op("sqrt", (a[0],), exact=True, positive_argument=True)
    s1 = _div(dag, a[1], _mul(dag, _const(dag, "2"), s0))
    s2 = _div(dag, _sub(dag, a[2], _mul(dag, _const(dag, "2"), s1, s1)),
              _mul(dag, _const(dag, "2"), s0))
    return [s0, s1, s2]


def _jet_power(dag, a, exponent):
    """Ordinary Z jet for a positive base and exact rational exponent."""
    p = Fraction(str(exponent))
    p_minus = p - 1
    p_minus_two = p - 2
    a_p = _pow(dag, a[0], p)
    a_pm = _pow(dag, a[0], p_minus)
    a_pm2 = _pow(dag, a[0], p_minus_two)
    pnode = _const(dag, p)
    first = _mul(dag, pnode, a_pm, a[1])
    second = _add(
        dag,
        _mul(dag, pnode, _const(dag, p - 1), a_pm2, a[1], a[1]),
        _mul(dag, pnode, a_pm, a[2]),
    )
    return [a_p, first, second]


def _jet_external(dag, archive_nodes, values, label, path, *, role="C6_source_jet",
                  source_archive=C6_ARCHIVE, report_sha=PRIOR_SHA256,
                  namespace="original_Rm_patch_C2_C6_source_projection_v1"):
    if not isinstance(values, list) or len(values) < 3:
        values = ["missing", "missing", "missing"]
    return [
        _c6_atom(dag, archive_nodes, values[i], label + "[" + str(i) + "]",
                 path + "[" + str(i) + "]", role=role, archive=source_archive)
        for i in range(3)
    ]


def _source_call(dag, source, *arguments, **metadata):
    return dag.op("source_call", (source,) + tuple(arguments), exact=True, **metadata)


def _integral(dag, integrand, variable, lower, upper, **metadata):
    return dag.op("integral", (integrand,), variable=variable,
                  bounds=[lower, upper], exact=True, unevaluated=True,
                  quadrature=False, **metadata)


def _nonlinear(dag, weights, h, alpha):
    zero = _jet_const(dag, "0")
    cross = _jet_add(
        dag,
        _jet_scale(dag, _jet_mul(dag, h[0], h[2]), weights["fg"][0]),
        _jet_scale(dag, _jet_mul(dag, h[1], h[4]), weights["fg"][1]),
    )
    axial = _jet_add(
        dag,
        _jet_scale(dag, _jet_mul(dag, h[0], h[0]), weights["gg"][0]),
        _jet_scale(dag, _jet_mul(dag, h[1], h[1]), weights["gg"][1]),
    )
    swirl = _jet_add(dag, *[
        _jet_scale(dag, _jet_mul(dag, h[i + 2], h[i + 2]), weights["ff"][i])
        for i in range(3)])
    pressure = _jet_add(dag, *[
        _jet_scale(dag, _jet_mul(dag, h[i + 2], h[i + 2]), weights["ff_over_x"][i])
        for i in range(3)])
    half = _jet_const(dag, "1/2")
    return [zero, cross, zero,
            _jet_sub(dag, _jet_mul(dag, alpha, axial), _jet_scale(dag, swirl, half)),
            _jet_scale(dag, pressure, half)]


def _source_gamma(dag, archive_nodes, binding, x, normalization, index):
    argument = _add(dag, _mul(dag, _const(dag, "40"),
                               _sub(dag, x, _const(dag, CENTERS[index]))))
    beta_fn = _ast_function_atom(dag, "original_true_Rm.raw_beta", "raw_beta")
    beta = _source_call(dag, beta_fn, argument, function="raw_beta",
                        branch="original_lazy_flat_branch",
                        exact_definition="exp(-1/(1-u^2)) for |u|<1; exact 0 outside",
                        active_variable=argument, flat_value="0",
                        active_body_not_evaluated_on_flat_branch=True)
    beta_derivative_fn = _ast_function_atom(
        dag, "original_true_Rm.beta_first_two", "beta_first_two")
    d1 = _source_call(
        dag, beta_derivative_fn, argument, beta,
        function="beta_first_two", component="first",
        branch="original_lazy_flat_branch", flat_value="0",
        active_body_not_evaluated_on_flat_branch=True,
        exact_interior_formula="beta*(-2*u)/(1-u^2)^2")
    dsecond = _source_call(
        dag, beta_derivative_fn, argument, beta,
        function="beta_first_two", component="second",
        branch="original_lazy_flat_branch", flat_value="0",
        active_body_not_evaluated_on_flat_branch=True,
        exact_interior_formula=("beta*(-2/(1-u^2)^2-8*u^2/(1-u^2)^3"
                               "+4*u^2/(1-u^2)^4)"))
    one = _one(dag)
    scale = _div(dag, one, normalization)
    return (
        _mul(dag, _const(dag, "40"), scale, beta),
        _mul(dag, _const(dag, "1600"), scale, d1),
        _mul(dag, _const(dag, "64000"), scale, dsecond),
    )


def _weight_port(dag, archive_nodes, binding, key, label, x, *, table="partial_weights"):
    table_data = binding.get(table, {}) if isinstance(binding, Mapping) else {}
    node = _get(table_data, key)
    source_label = "Wpartial" if table == "partial_weights" else "W"
    path_prefix = "bindings." + table + "."
    port = _c6_atom(dag, archive_nodes, node, label, path_prefix + key,
                    role=source_label + "_integral",
                    bound_variable="true_Rm_patch_x" if table == "partial_weights" else None)
    if table == "partial_weights":
        return _source_call(dag, port, x,
                            source_variable="true_Rm_patch_x", source_variable_handle=185676,
                            source_coordinate_substitution=True,
                            shared_Z_coordinate="Z", source_Z_variable_handle=20)
    return port


def _partial_weights(dag, archive_nodes, binding, x):
    partial = binding.get("partial_weights", {})
    zero = _zero(dag)
    L = [
        [_weight_port(dag, archive_nodes, binding, "L.0.0", "Wpartial.L[0][0]", x),
         _weight_port(dag, archive_nodes, binding, "L.0.1", "Wpartial.L[0][1]", x), zero, zero, zero],
        [_weight_port(dag, archive_nodes, binding, "L.1.0", "Wpartial.L[1][0]", x),
         _weight_port(dag, archive_nodes, binding, "L.1.1", "Wpartial.L[1][1]", x), zero, zero, zero],
        [zero, zero,
         _weight_port(dag, archive_nodes, binding, "L.2.2", "Wpartial.L[2][2]", x),
         _weight_port(dag, archive_nodes, binding, "L.2.3", "Wpartial.L[2][3]", x),
         _weight_port(dag, archive_nodes, binding, "L.2.4", "Wpartial.L[2][4]", x)],
        [zero, zero,
         _weight_port(dag, archive_nodes, binding, "L.3.2", "Wpartial.L[3][2]", x),
         _weight_port(dag, archive_nodes, binding, "L.3.3", "Wpartial.L[3][3]", x),
         _weight_port(dag, archive_nodes, binding, "L.3.4", "Wpartial.L[3][4]", x)],
        [zero, zero,
         _weight_port(dag, archive_nodes, binding, "L.4.2", "Wpartial.L[4][2]", x),
         _weight_port(dag, archive_nodes, binding, "L.4.3", "Wpartial.L[4][3]", x),
         _weight_port(dag, archive_nodes, binding, "L.4.4", "Wpartial.L[4][4]", x)],
    ]
    def rows(name, count):
        return [_weight_port(dag, archive_nodes, binding, name + "." + str(i),
                             "Wpartial." + name + "[" + str(i) + "]", x)
                for i in range(count)]
    return {"L": L, "fg": rows("fg", 2), "gg": rows("gg", 2),
            "ff": rows("ff", 3), "ff_over_x": rows("ff_over_x", 3),
            "binding": partial, "zero": zero, "table": "partial_weights"}


def _full_weights(dag, archive_nodes, binding):
    zero = _zero(dag)

    def port(key, label):
        return _weight_port(dag, archive_nodes, binding, key, label, None, table="W")

    L = [
        [port("L.0.0", "W.L[0][0]"), port("L.0.1", "W.L[0][1]"), zero, zero, zero],
        [port("L.1.0", "W.L[1][0]"), port("L.1.1", "W.L[1][1]"), zero, zero, zero],
        [zero, zero, port("L.2.2", "W.L[2][2]"), port("L.2.3", "W.L[2][3]"), port("L.2.4", "W.L[2][4]")],
        [zero, zero, port("L.3.2", "W.L[3][2]"), port("L.3.3", "W.L[3][3]"), port("L.3.4", "W.L[3][4]")],
        [zero, zero, port("L.4.2", "W.L[4][2]"), port("L.4.3", "W.L[4][3]"), port("L.4.4", "W.L[4][4]")],
    ]

    def rows(name, count):
        return [port(name + "." + str(i), "W." + name + "[" + str(i) + "]")
                for i in range(count)]

    return {"L": L, "fg": rows("fg", 2), "gg": rows("gg", 2),
            "ff": rows("ff", 3), "ff_over_x": rows("ff_over_x", 3),
            "binding": binding.get("W", {}) if isinstance(binding, Mapping) else {},
            "zero": zero, "table": "W"}


def _map_residual(dag, weights, controls, defects, alpha):
    nonlinear = _nonlinear(dag, weights, controls, alpha)
    residual = []
    for i in range(5):
        linear = _jet_add(dag, *[
            _jet_scale(dag, controls[j], _jet_source(dag, weights["L"][i][j],
                                                      "L[%d][%d]" % (i, j)))
            for j in range(5) if weights["L"][i][j] != weights.get("zero")])
        residual.append(_jet_add(dag, defects[i], nonlinear[i], linear))
    return residual, nonlinear


def _patch_fields(dag, xjet, zjet, am, pstar, p0, defects, controls, alpha, weights, gamma):
    inv_pstar = _jet_inv(dag, pstar)
    changes = _nonlinear(dag, weights, controls, alpha)
    residual = []
    for i in range(5):
        linear = _jet_add(dag, *[
            _jet_scale(dag, controls[j], _jet_source(dag, weights["L"][i][j],
                                                      "Wpartial.L[%d][%d]" % (i, j)))
            for j in range(5) if weights["L"][i][j] != weights.get("zero")])
        residual.append(_jet_add(dag, defects[i], changes[i], linear))
    xpow_1_10 = _jet_power(dag, xjet, "1/10")
    xpow_m9_10 = _jet_power(dag, xjet, "-9/10")
    xpow_8_5 = _jet_power(dag, xjet, "8/5")
    xpow_6_5 = _jet_power(dag, xjet, "6/5")
    xpow_1_5 = _jet_power(dag, xjet, "1/5")
    xpow_m3_2 = _jet_power(dag, xjet, "-3/2")
    inv_x = _jet_inv(dag, xjet)
    H = _jet_add(dag, xpow_1_10, *[
        _jet_scale(dag, controls[i + 2], _jet_source(dag, gamma[i][0], "gamma[%d]" % i))
        for i in range(3)])
    axial = _jet_add(dag,
        _jet_scale(dag, controls[0], _jet_source(dag, gamma[0][0], "gamma[0]")),
        _jet_scale(dag, controls[1], _jet_source(dag, gamma[2][0], "gamma[2]")))
    Hx = _jet_add(dag, _jet_scale(dag, xpow_m9_10, _jet_const(dag, "1/10")), *[
        _jet_scale(dag, controls[i + 2], _jet_source(dag, gamma[i][1], "gamma_x[%d]" % i))
        for i in range(3)])
    axial_x = _jet_add(dag,
        _jet_scale(dag, controls[0], _jet_source(dag, gamma[0][1], "gamma_x[0]")),
        _jet_scale(dag, controls[1], _jet_source(dag, gamma[2][1], "gamma_x[2]")))
    four_z = _jet_scale(dag, zjet, _const(dag, "4"))
    mass = _jet_add(dag, four_z, _jet_mul(dag, residual[0], inv_x))
    theta = _jet_add(dag, residual[2], _jet_scale(dag, xpow_8_5, _const(dag, "5/8")))
    mixed = _jet_add(dag, _jet_scale(dag, _jet_mul(dag, zjet, theta), _const(dag, "4")), residual[1])
    energy = _jet_sub(dag, residual[3], _jet_scale(dag, xpow_6_5, _const(dag, "5/12")))
    pressure_moment = _jet_add(dag, residual[4], _jet_scale(dag, xpow_1_5, _const(dag, "5/2")))
    am2 = _jet_mul(dag, am, am)
    inv_pstar2 = _jet_mul(dag, inv_pstar, inv_pstar)
    first_energy = _jet_mul(
        dag, _jet_scale(dag, _jet_mul(dag, zjet, mass), _const(dag, "8")), inv_pstar2)
    second_energy = _jet_mul(
        dag, _jet_scale(dag, _jet_mul(dag, zjet, zjet), _const(dag, "-16")), inv_pstar2)
    raw_energy = _jet_add(
        dag, first_energy, second_energy,
        _jet_mul(dag, _jet_mul(dag, am2, energy), inv_x))
    fields = {
        "E": _jet_mul(dag, am, H),
        "V": _jet_scale(dag, _jet_add(dag, four_z, axial), inv_pstar),
        "m": _jet_scale(dag, mass, inv_pstar),
        "h": _jet_scale(dag, _jet_mul(dag, am, theta), xpow_m3_2),
        "k": _jet_scale(dag, _jet_mul(dag, am, mixed), _jet_mul(dag, xpow_m3_2, inv_pstar)),
        "e": raw_energy,
        "p": _jet_mul(dag, am2, pressure_moment),
        "P0": p0,
    }
    fields["pressure"] = _jet_add(dag, fields["P0"], fields["p"])
    physical = {
        "U": _jet_mul(dag, pstar, fields["V"]),
        "pressure": fields["pressure"],
        "P0": fields["P0"],
    }
    radial = {
        "E_y": _jet_scale(dag, _jet_mul(dag, am, Hx), xjet),
        "V_y": _jet_scale(dag, _jet_mul(dag, axial_x, inv_pstar), xjet),
    }
    primitives = {
        "actual_partial_defects": residual,
        "partial_changes": changes,
        "H": H,
        "axial_correction": axial,
        "physical_mass": mass,
        "theta": theta,
        "mixed": mixed,
        "centered_energy": energy,
        "pressure_moment": pressure_moment,
    }
    return fields, physical, radial, primitives


def _inlet_centered_histories(dag, archive_nodes, binding, zjet):
    inlet = binding.get("inlet", {})
    actual = inlet.get("actual", {})
    moments = inlet.get("moments", {})
    phi = _jet_external(dag, archive_nodes, actual.get("phi"), "Rm_inlet.phi",
                         "bindings.inlet.actual.phi")
    M = _jet_external(dag, archive_nodes, moments.get("M"), "Rm_inlet.M",
                      "bindings.inlet.moments.M")
    H = _jet_external(dag, archive_nodes, moments.get("H"), "Rm_inlet.H",
                      "bindings.inlet.moments.H")
    K = _jet_external(dag, archive_nodes, moments.get("K"), "Rm_inlet.K",
                      "bindings.inlet.moments.K")
    A = _jet_external(dag, archive_nodes, moments.get("A"), "Rm_inlet.A",
                      "bindings.inlet.moments.A")
    B = _jet_external(dag, archive_nodes, moments.get("B"), "Rm_inlet.B",
                      "bindings.inlet.moments.B")
    C = _jet_external(dag, archive_nodes, moments.get("C"), "Rm_inlet.C",
                      "bindings.inlet.moments.C")
    four_z = _jet_scale(dag, zjet, _const(dag, "4"))
    two_phi = _jet_scale(dag, phi, _const(dag, "2"))
    Hn = _jet_div(dag, H, two_phi)
    Kn = _jet_div(dag, K, two_phi)
    phi2 = _jet_mul(dag, phi, phi)
    return {
        "centered": {
            "em": _jet_sub(dag, M, four_z),
            "ek": _jet_sub(dag, Kn, _jet_scale(dag, _jet_mul(dag, zjet, Hn), _const(dag, "4"))),
            "eh": _jet_sub(dag, Hn, _jet_const(dag, "5/8")),
            "aa": _jet_add(dag, _jet_sub(dag, A, _jet_scale(dag, _jet_mul(dag, zjet, M), _const(dag, "8"))),
                            _jet_scale(dag, _jet_mul(dag, zjet, zjet), _const(dag, "16"))),
            "eb": _jet_sub(dag, _jet_div(dag, B, phi2), _jet_const(dag, "5/6")),
            "ep": _jet_sub(dag, _jet_div(dag, C, phi2), _jet_const(dag, "5")),
        },
        # These six signed inlet channels are retained as imported source data.
        # They are not propagated as six new history outputs in this stage.
        "retained_raw": {"M": M, "H": H, "K": K, "A": A, "B": B, "C": C, "phi": phi},
    }


def _adapter_jet(dag, adapter_ids, rows, label, path):
    """Bind a saved adapter endpoint jet without evaluating its definition."""
    if not isinstance(rows, list) or len(rows) < 3:
        rows = ["missing", "missing", "missing"]
    return [
        _saved_atom(dag, adapter_ids, rows[i], label + "[" + str(i) + "]",
                    path + "[" + str(i) + "]")
        for i in range(3)
    ]


def _jet_substitute(dag, jet, x, value, *, endpoint=None):
    """Create a source-preserving ordinary-Z jet with x replaced by value."""
    return [dag.op("function_substitution", (ref, x, value), exact=True,
                    ordinary_Z_order=i, radial_substitution_Z_independent=True,
                    endpoint=endpoint, source_function=True)
            for i, ref in enumerate(jet)]


def _own_histories(dag, xjet, fields, seeds):
    """Build the five normalized leading histories on the radial log scale.

    For each rate r, the source equation is x H_x + r H = source, hence
    H(x)=x^-r H(1)+integral_1^x (s/x)^r source(s) ds/s.  The source is
    the constructed patch field and is substituted x -> s before the
    integral.  No inlet generic transport rows are interpreted as histories.
    """
    x = xjet[0]
    s = _symbol(dag, "s", "[1,e]")
    one = _one(dag)
    source_fields = {
        "m": fields["V"],
        "h": fields["E"],
        "k": _jet_mul(dag, fields["E"], fields["V"]),
        "e": _jet_sub(dag, _jet_mul(dag, fields["V"], fields["V"]),
                       _jet_scale(dag, _jet_mul(dag, fields["E"], fields["E"]), _const(dag, "1/2"))),
        "p": _jet_scale(dag, _jet_mul(dag, fields["E"], fields["E"]), _const(dag, "1/2")),
    }
    histories = {}
    integral_nodes = {}
    for name in LEADING_HISTORY_NAMES:
        rate = RATES[name]
        source = _jet_substitute(dag, source_fields[name], x, s,
                                 endpoint="history_source_x_to_s")
        memory = _jet_power(dag, xjet, "-" + rate) if rate != "0" else _jet_const(dag, "1")
        propagated = _jet_mul(dag, seeds[name], memory)
        # The radial kernel is Z independent, but its zeroth row multiplies
        # every ordinary Z derivative of the source through Leibniz' rule.
        sjet = [s, _zero(dag), _zero(dag)]
        ratio = _jet_mul(dag, _jet_power(dag, sjet, rate), memory)
        ratio = _jet_mul(dag, ratio, _jet_inv(dag, sjet))
        integrand_jet = _jet_mul(dag, ratio, source)
        integral = []
        for order in range(3):
            integral.append(_integral(
                dag, integrand_jet[order], s, one, x,
                kernel="(s/x)^rate/s", rate=rate,
                exact_bound_interval="[1,e]", nonreset_incoming_history=True,
                source_function="constructed_patch_field", history=name,
                source_field=name, source_substitution="x->s",
                ordinary_Z_order=order))
        histories[name] = _jet_add(dag, propagated, integral)
        integral_nodes[name] = integral
    return histories, integral_nodes, source_fields


def _substitutions(dag, x, value, groups, label):
    def substitute(item):
        if isinstance(item, list) and len(item) == 3 and all(isinstance(ref, str) for ref in item):
            return [dag.op("function_substitution", (ref, x, value),
                           exact=True, ordinary_Z_order=i,
                           radial_substitution_Z_independent=True,
                           endpoint=label, source_function=True)
                    for i, ref in enumerate(item)]
        if isinstance(item, Mapping):
            return {name: substitute(child) for name, child in item.items()}
        if isinstance(item, list):
            return [substitute(child) for child in item]
        return item

    return {group_name: substitute(group) for group_name, group in groups.items()}


def _iter_jets(value):
    if isinstance(value, list) and len(value) == 3 and all(isinstance(ref, str) for ref in value):
        yield value
    elif isinstance(value, Mapping):
        for child in value.values():
            yield from _iter_jets(child)
    elif isinstance(value, list):
        for child in value:
            yield from _iter_jets(child)


def _graph_checks(dag, jets, histories, endpoint_substitutions):
    seen = set()
    checks = {
        "unique_namespace_ids": True,
        "closed_topological_arguments": True,
        "rational_constants": True,
        "integrals_unevaluated_and_bound": True,
        "qualified_source_atoms": True,
        "ordinary_Z_jets_have_three_rows": True,
        "no_numeric_or_solver_operation": True,
        "endpoint_substitutions_are_function_nodes": True,
    }
    forbidden = {"evaluate", "quadrature", "norm", "solve", "implicit_solve", "recurrence"}
    for node in dag.nodes:
        ident = node.get("id")
        checks["unique_namespace_ids"] &= ident not in seen and node.get("namespace") == SCHEMA
        checks["closed_topological_arguments"] &= all(arg in seen for arg in node.get("args", []))
        if node.get("op") == "exact_constant":
            try:
                Fraction(node.get("value"))
            except (ValueError, TypeError, ZeroDivisionError):
                checks["rational_constants"] = False
        if node.get("op") == "integral":
            bounds = node.get("bounds", [])
            checks["integrals_unevaluated_and_bound"] &= (
                node.get("unevaluated") is True and node.get("quadrature") is False
                and len(bounds) == 2 and all(bound in seen for bound in bounds)
                and node.get("variable") in seen)
        if node.get("op") == "source_function_port":
            value = node.get("value", {})
            checks["qualified_source_atoms"] &= (
                isinstance(value, Mapping) and isinstance(value.get("source_report_sha256"), str)
                and isinstance(value.get("source_namespace"), str)
                and "source_node" in value and value.get("source_node") != "missing"
                and value.get("re_evaluated") is False)
            if isinstance(value, Mapping) and value.get("role") != "source_AST_function":
                flags = [value[name] for name in ("source_graph_node_valid", "source_node_valid")
                         if name in value]
                checks["qualified_source_atoms"] &= bool(flags) and all(flag is True for flag in flags)
        checks["no_numeric_or_solver_operation"] &= node.get("op") not in forbidden
        if node.get("op") == "function_substitution":
            checks["endpoint_substitutions_are_function_nodes"] &= len(node.get("args", [])) == 3
        seen.add(ident)
    for collection in (jets, histories):
        for jet in _iter_jets(collection):
            checks["ordinary_Z_jets_have_three_rows"] &= isinstance(jet, list) and len(jet) == 3 and all(ref in seen for ref in jet)
    for endpoint in endpoint_substitutions.values():
        for jet in _iter_jets(endpoint):
            checks["endpoint_substitutions_are_function_nodes"] &= (
                isinstance(jet, list) and all(ref in seen for ref in jet))
    return checks


def _source_audit(adapter, prior, archive, c6_check, c6_tree, c6_text, gates):
    ac = adapter.get("component_result", {}) if isinstance(adapter, Mapping) else {}
    pc = prior.get("component_result", {}) if isinstance(prior, Mapping) else {}
    c6audit = pc.get("c6_archive_port_audit", {}) if isinstance(pc, Mapping) else {}
    branch = c6audit.get("branch", {}) if isinstance(c6audit, Mapping) else {}
    checks = c6audit.get("checks", {}) if isinstance(c6audit, Mapping) else {}
    binding = archive.get("true_Rm_patch_C6_bindings", {}) if isinstance(archive, Mapping) else {}
    graph = archive.get("exact_graph_nodes", []) if isinstance(archive, Mapping) else []
    p1_fields = _get(ac, "function_definitions.endpoint_substitutions.postrestore_p1.normalized_fields", {})
    control_ids = []
    for row in binding.get("controls", []) if isinstance(binding, Mapping) else []:
        if isinstance(row, list) and row:
            control_ids.append(row[0])
    source_witness = {
        "adapter_result_accepted": ac.get("passed") is True
            and ac.get("native_C2_reference_adapter_defined") is True
            and ac.get("native_C2_reference_Rm_identity_proved") is True
            and ac.get("legacy_feedback_reference_alias_proved") is False,
        "adapter_graph_saved": _get(ac, "definition_graph.namespace") == "original_native_C2_reference_adapter_v1"
            and _get(ac, "definition_graph.active_graph_digest") == EXPECTED_NODES_DIGEST
            and type(_get(ac, "definition_graph.node_count")) is int,
        "adapter_postrestore_p1_fields_saved": isinstance(p1_fields, Mapping)
            and all(isinstance(p1_fields.get(name), list) and len(p1_fields[name]) == 3
                    for name in LEADING_HISTORY_NAMES),
        "prior_C6_port_audit_accepted": c6audit.get("valid") is True
            and c6audit.get("accepted_C6_Rm_inlet_identity") is True,
        "prior_legacy_alias_remains_negative": pc.get("unresolved_gates") == ["feedback_reference_Rm_source_alias_unproved"],
        "c6_archive_prefix_exact": len(graph) == 196042 and c6audit.get("qualified_handles", {}).get("archive_graph_prefix_length") == 196042,
        "c6_branch_unique_common_box": checks.get("branch_common_box_unique") is True
            and checks.get("branch_transfer_source_identity") is True
            and checks.get("accepted_Rm_inlet_endpoint_identity") is True,
        "c6_system_exact": type(binding.get("system")) is int and binding.get("system") == 193470
            and isinstance(graph, list) and len(graph) > 193470
            and graph[193470].get("operation") == "original_Rm_unique_leading_implicit_system"
            and graph[193470].get("defining_map") == "L*h+Q(h,invAm2)+true_Rm_defects=0"
            and graph[193470].get("original_whole_Z_base_unique_branch_and_Jacobian_admitted") is True,
        "c6_controls_exact": len(control_ids) == 5 and all(
            type(node) is int and 0 <= node < len(graph)
            and graph[node].get("operation") == "original_Rm_unique_leading_implicit_control"
            and graph[node].get("system") == 193470
            and graph[node].get("component") == i
            and graph[node].get("exact_value_of_same_unique_original_common_box_solution") is True
            for i, node in enumerate(control_ids)),
        "c6_source_formula_corrected": _source_has(
            c6_text,
            "a.scale(a.mul(z,z),g.mul(g.constant(16),inv2))",
            "a.scale(a.mul(am2,energy),invx)",
            "fields=dict(E=a.mul(am,H),V=a.scale(a.add(a.scale(z,4),axial),inv)",
            "def inlet_defects",
            "def patch_fields",
            "def nonlinear",
            "def beta_first_two"),
        "c6_receipt_accepted": c6_check.get("all_passed") is True
            and c6_check.get("same_exact_original_unique_branch_theorem_bound") is True
            and c6_check.get("current_original_true_Rm_leading_implicit_control_C6_functions_installed") is True,
        "ordinary_Z_source_scope": ac.get("ordinary_Z_orders") == [0, 1, 2]
            and pc.get("ordinary_Z_orders") == [0, 1, 2]
            and ac.get("ordinary_Z_not_temporal_n") is True
            and pc.get("ordinary_Z_not_temporal_n") is True,
        "fixed_N0_exact_integer_premise": _get(binding, "branch_transfer.actual_leading_N0_equal_by_exact_integer") is True,
    }
    for name, ok in source_witness.items():
        if not ok:
            gates.append("source_witness_missing:" + name)
    return source_witness, binding, graph, control_ids, branch, checks


def derive(nodes, source_records, *, expected_nodes_digest):
    """Build an authenticated source-only native C2 leading patch DAG."""
    gates = []
    try:
        actual_digest = _digest(nodes)
    except Exception as error:
        actual_digest = ""
        gates.append("nodes_digest_error:" + type(error).__name__)
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append("expected_nodes_digest_mismatch")
    if actual_digest != EXPECTED_NODES_DIGEST:
        gates.append("nodes_digest_mismatch")
    if not isinstance(nodes, list) or len(nodes) != EXPECTED_NODE_COUNT:
        gates.append("node_count_mismatch")
    records = source_records if isinstance(source_records, Mapping) else {}
    if not isinstance(source_records, Mapping):
        gates.append("source_records_not_mapping")
    texts = {name: _record(records, name, gates) for name in RECORD_HASHES}
    adapter = _json(texts[ADAPTER_REPORT], ADAPTER_REPORT, gates)
    prior = _json(texts[PRIOR_REPORT], PRIOR_REPORT, gates)
    archive = _json(texts[C6_ARCHIVE], C6_ARCHIVE, gates)
    c6_check = _json(texts[C6_CHECK], C6_CHECK, gates)
    c6_tree = _tree(texts[C6_SOURCE], C6_SOURCE, gates)
    witness, binding, archive_nodes, control_ids, branch, c6_checks = _source_audit(
        adapter, prior, archive, c6_check, c6_tree, texts[C6_SOURCE], gates)
    adapter_component = adapter.get("component_result", {}) if isinstance(adapter, Mapping) else {}
    adapter_graph = adapter_component.get("definition_graph", {}) if isinstance(adapter_component, Mapping) else {}
    adapter_nodes = adapter_graph.get("nodes", []) if isinstance(adapter_graph, Mapping) else []
    adapter_ids = {row.get("id") for row in adapter_nodes if isinstance(row, Mapping)}
    dag = _DAG()
    z = _symbol(dag, "Z", "[-1,1]")
    x = _symbol(dag, "x", "[1,e]")
    zjet = [z, _one(dag), _zero(dag)]
    xjet = [x, _zero(dag), _zero(dag)]
    adapter_outputs = adapter_graph.get("outputs", {}) if isinstance(adapter_graph, Mapping) else {}
    adapter_defs = adapter_component.get("function_definitions", {}) if isinstance(adapter_component, Mapping) else {}
    endpoint_defs = adapter_defs.get("endpoint_substitutions", {}) if isinstance(adapter_defs, Mapping) else {}
    post_p1 = endpoint_defs.get("postrestore_p1", {}) if isinstance(endpoint_defs, Mapping) else {}
    post_fields = post_p1.get("normalized_fields", {}) if isinstance(post_p1, Mapping) else {}
    adapter_defect_rows = adapter_outputs.get("rm_p1_defect_jet", []) if isinstance(adapter_outputs, Mapping) else []
    defects = []
    for i, row in enumerate(adapter_defect_rows[:5]):
        values = row if isinstance(row, list) else []
        if len(values) < 3:
            values = ["missing", "missing", "missing"]
        defects.append([
            _saved_atom(dag, adapter_ids, values[order],
                        "adapter.Rm_defect[%d][%d]" % (i, order),
                        "definition_graph.outputs.rm_p1_defect_jet[%d][%d]" % (i, order))
            for order in range(3)
        ])
    while len(defects) < 5:
        defects.append(_jet_const(dag, "0"))
        gates.append("adapter_defect_jet_missing")
    pstar_port = _get(adapter_component, "qualified_endpoint_ports.parameter.Pstar")
    pstar_node = pstar_port if type(pstar_port) is int else _get(adapter_component, "qualified_endpoint_ports.parameter.Pstar.node", "missing")
    pstar = _jet_source(dag, _source_atom(
        dag, report_sha=ADAPTER_SHA256,
        namespace="original_native_C2_reference_adapter_v1",
        node=pstar_node, label="adapter.parameter.Pstar",
        source_archive=ADAPTER_REPORT, role="saved_parameter_port",
        path="qualified_endpoint_ports.parameter.Pstar", re_evaluated=False,
        source_node_valid=type(pstar_node) is int and pstar_node == 33572),
        "Pstar")
    am_numerator = _jet_exp(dag, _jet_const(dag, "-3/5"))
    am = _jet_div(dag, am_numerator,
                  _jet_add(dag, _jet_const(dag, "1"), _jet_mul(dag, zjet, zjet)))
    alpha = _jet_inv(dag, _jet_mul(dag, _jet_mul(dag, pstar, pstar), _jet_mul(dag, am, am)))
    inlet_data = _inlet_centered_histories(dag, archive_nodes, binding, zjet)
    history_seeds = {
        name: _adapter_jet(
            dag, adapter_ids, post_fields.get(name),
            "adapter.endpoint_substitutions.postrestore_p1.normalized_fields." + name,
            "function_definitions.endpoint_substitutions.postrestore_p1.normalized_fields." + name)
        for name in LEADING_HISTORY_NAMES
    }
    p0 = _adapter_jet(dag, adapter_ids, post_fields.get("P0"), "adapter.Rm.P0",
                      "function_definitions.endpoint_substitutions.postrestore_p1.normalized_fields.P0")
    rm_radius = _adapter_jet(dag, adapter_ids, post_p1.get("radius"), "adapter.Rm.radius",
                            "function_definitions.endpoint_substitutions.postrestore_p1.radius")
    geometry = {"Rm_radius": rm_radius, "radius": _jet_mul(dag, rm_radius, xjet),
                "physical_measure": _jet_inv(dag, xjet), "native_x": xjet}
    normalization = _c6_atom(dag, archive_nodes, binding.get("exact_beta_normalization"),
                             "original_Rm.beta_normalization", "bindings.exact_beta_normalization",
                             role="C6_exact_unevaluated_integral")
    gamma = [_source_gamma(dag, archive_nodes, binding, x, normalization, i) for i in range(3)]
    weights = _partial_weights(dag, archive_nodes, binding, x)
    control_rows = binding.get("controls", []) if isinstance(binding, Mapping) else []
    controls = []
    for i, row in enumerate(control_rows[:5]):
        controls.append(_jet_external(dag, archive_nodes,
                                       row[:3] if isinstance(row, list) else [],
                                       "Rm_control[%d]" % i,
                                       "bindings.controls[%d]" % i,
                                       role="C6_selected_unique_leading_control"))
    while len(controls) < 5:
        controls.append(_jet_const(dag, "0"))
        gates.append("c6_control_row_missing")
    full_weights = _full_weights(dag, archive_nodes, binding)
    map_residual, nonlinear = _map_residual(dag, full_weights, controls, defects, alpha)
    terminal_identity = [
        dag.op("source_identity", (map_residual[i][0],), exact=True,
                identity="L*h+Q(h,alpha)+D=0", component=i,
                source_control_node=control_ids[i] if i < len(control_ids) else "missing",
                endpoint_substitution="x=e",
                exact_full_weight_source_ports=True,
                conditional_common_box_branch=True, numerical_zero_not_claimed=True,
                solver_reentered=False)
        for i in range(5)
    ]
    fields, physical, radial, primitives = _patch_fields(
        dag, xjet, zjet, am, pstar, p0, defects, controls, alpha, weights, gamma)
    physical["Utheta"] = _jet_mul(dag, pstar, fields["E"])
    physical["Uz"] = physical["U"]
    own_histories, history_integrals, history_sources = _own_histories(
        dag, xjet, fields, history_seeds)
    all_jets = list(fields.values()) + list(physical.values()) + list(radial.values())
    all_jets += list(primitives.values()) + list(own_histories.values())
    all_jets += defects + controls + [am, alpha, pstar, p0] + list(geometry.values())
    endpoint_groups = {
        "fields": fields,
        "physical": physical,
        "radial": radial,
        "primitives": primitives,
        "own_histories": own_histories,
        "geometry": geometry,
    }
    x_one = _one(dag)
    x_e = _exp(dag, _one(dag))
    endpoint_substitutions = {
        "x=1": _substitutions(dag, x, x_one, endpoint_groups, "x=1"),
        "x=e": _substitutions(dag, x, x_e, endpoint_groups, "x=e"),
    }
    graph_checks = _graph_checks(
        dag, all_jets, own_histories, endpoint_substitutions)
    for name, ok in graph_checks.items():
        if not ok:
            gates.append("definition_graph_check_failed:" + name)
    adapter_accept = (witness["adapter_result_accepted"] and witness["adapter_graph_saved"]
                      and witness["adapter_postrestore_p1_fields_saved"])
    c6_accept = all(witness[name] for name in (
        "prior_C6_port_audit_accepted", "c6_archive_prefix_exact", "c6_branch_unique_common_box",
        "c6_system_exact", "c6_controls_exact", "c6_source_formula_corrected", "c6_receipt_accepted"))
    graph_ok = all(graph_checks.values())
    patch_formula_ok = all(
        isinstance(fields.get(name), list) and len(fields[name]) == 3
        for name in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")) and (
            "a.scale(a.mul(z,z),g.mul(g.constant(16),inv2))" in texts[C6_SOURCE]
            and "a.scale(a.mul(am2,energy),invx)" in texts[C6_SOURCE])
    history_ok = all(isinstance(own_histories.get(name), list) and len(own_histories[name]) == 3
                     for name in RATES) and all(history_integrals.values())
    same_control_map = bool(c6_accept and len(terminal_identity) == 5 and not gates
                            and witness["fixed_N0_exact_integer_premise"])
    leading_defined = bool(adapter_accept and c6_accept and patch_formula_ok and graph_ok
                           and len(fields) == 9 and len(primitives) == 9 and not gates)
    source_identity = bool(leading_defined and same_control_map and history_ok
                           and witness["ordinary_Z_source_scope"])
    own_history_identity = bool(leading_defined and history_ok
                               and witness["prior_C6_port_audit_accepted"]
                               and _get(adapter_component, "source_prerequisite_audit.saved_R110_background_and_own_history") is True)
    if not leading_defined:
        gates.append("native_C2_leading_patch_definitions_unproved")
    if not source_identity:
        gates.append("native_C2_C6_patch_source_identity_unproved")
    if not own_history_identity:
        gates.append("native_C2_patch_own_history_identity_unproved")
    if not same_control_map:
        gates.append("same_unique_leading_control_map_unproved")
    gates = list(dict.fromkeys(gates))
    retained = _get(adapter, "retained_construction_data", {})
    source_audit = {
        "source_records_only": True,
        "static_source_AST_graph_saved_data_and_new_definitions_only": True,
        "original_graph_nodes_inspected_not_evaluated": True,
        "old_graph_inspected_not_evaluated": True,
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
    }
    return {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(leading_defined and source_identity and own_history_identity and same_control_map and not gates),
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "native_C2_leading_patch_defined": leading_defined,
        "native_C2_C6_patch_source_identity_proved": source_identity,
        "native_C2_patch_own_history_identity_proved": own_history_identity,
        "same_unique_leading_control_map_proved": same_control_map,
        "legacy_feedback_reference_alias_proved": False,
        "definition_graph": {
            "namespace": SCHEMA,
            "active_graph_digest": EXPECTED_NODES_DIGEST,
            "nodes": dag.nodes,
            "node_count": len(dag.nodes),
            "outputs": {
                "normalized_fields": fields,
                "physical_fields": physical,
                "radial_fields": radial,
                "geometry": geometry,
                "primitives": primitives,
                "defects": defects,
                "controls": controls,
                "retained_inlet_data": inlet_data,
                "partial_weights": weights,
                "terminal_full_weights": full_weights,
                "leading_map_residual": map_residual,
                "leading_terminal_identity": terminal_identity,
                "own_leading_histories": own_histories,
                "own_history_source_fields": history_sources,
                "own_history_integrals": history_integrals,
                "endpoint_substitutions": endpoint_substitutions,
            },
        },
        "function_definitions": {
            "normalized_fields": fields,
            "physical_fields": physical,
            "radial_fields": radial,
            "geometry": geometry,
            "primitives": primitives,
            "retained_inlet_data": inlet_data,
            "own_leading_histories": own_histories,
            "own_history_source_fields": history_sources,
            "endpoint_substitutions": endpoint_substitutions,
            "ordinary_Z_jet_rules": {
                "orders": [0, 1, 2],
                "product_second": "a2*b0 + 2*a1*b1 + a0*b2",
                "reciprocal": "(1/a0, -a1/a0^2, (2*a1^2-a0*a2)/a0^3)",
                "exp": "(exp(a0), exp(a0)*a1, exp(a0)*(a2+a1^2))",
                "sqrt": "(sqrt(a0), a1/(2*sqrt(a0)), (a2-2*(a1/(2*sqrt(a0)))^2)/(2*sqrt(a0)))",
                "power": "(a0^p, p*a0^(p-1)*a1, p*(p-1)*a0^(p-2)*a1^2+p*a0^(p-1)*a2)",
                "spatial_variable": "Z",
                "radial_variable": "x in [1,e]",
                "radial_derivative": "D_y=x*D_x",
                "temporal_recursion_claimed": False,
            },
        },
        "qualified_source_ports": {
            "adapter_report_sha256": ADAPTER_SHA256,
            "adapter_namespace": "original_native_C2_reference_adapter_v1",
            "adapter_defect": "definition_graph.outputs.rm_p1_defect_jet",
            "C6_archive": C6_ARCHIVE,
            "C6_archive_sha256": C6_ARCHIVE_SHA256,
            "C6_prefix": 196042,
            "C6_system": 193470,
            "C6_control_rows": control_ids,
            "C6_inlet": _get(binding, "inlet", {}),
            "C6_partial_weight_rows": _get(binding, "partial_weights", {}),
            "source_control_rows": [
                {"component": i, "node": control_ids[i] if i < len(control_ids) else "missing",
                 "same_unique_common_box_root": i < len(control_ids)} for i in range(5)
            ],
        },
        "source_prerequisite_audit": witness,
        "source_identity_argument": {
            "new_owner_only": True,
            "legacy_feedback_owner_identity_assumed": False,
            "initial_data": "The new owner consumes the saved adapter Rm defect, normalized endpoint fields, and retained C6 Rm inlet ports.",
            "defect": "D=(em,ek,eh,alpha*aa-eb/2,ep/2) from the saved native adapter source definition",
            "canonical_amplitude": "Am=exp(-3/5)/(1+Z^2)",
            "alpha": "Pstar^-2*Am^-2",
            "leading_map": "L*h+Q(h,alpha)+D=0",
            "uniqueness": "Selected C6 controls are the exact original common-box roots; the branch and fixed N0 premises remain conditional source claims.",
            "numeric_equality_or_interval_overlap_used": False,
            "solver_reentered": False,
        },
        "full_profile_formulas": {
            "x_domain": "[1,e]",
            "mass": "4Z + residual[0]/x",
            "theta": "residual[2] + (5/8)*x^(8/5)",
            "mixed": "4Z*theta + residual[1]",
            "centered_energy": "residual[3] - (5/12)*x^(6/5)",
            "pressure_moment": "residual[4] + (5/2)*x^(1/5)",
            "raw_energy": "8Z*mass/Pstar^2 - 16Z^2/Pstar^2 + Am^2*centered_energy/x",
            "normalized_fields": "E=Am*H; V=(4Z+axial)/Pstar; m=mass/Pstar; h=Am*theta*x^(-3/2); k=Am*mixed*x^(-3/2)/Pstar; e=raw_energy; p=Am^2*pressure_moment; pressure=P0+p",
            "physical_velocity": "Uz=Pstar*V; Utheta=Pstar*E",
            "radial": "E_y=x*Am*H_x; V_y=x*axial_x/Pstar",
            "P0": "retained analytic native Rm inlet P0 source function",
        },
        "own_history_identity": {
            "names": list(RATES),
            "rates": dict(RATES),
            "incoming_seed": "saved adapter endpoint normalized fields m,h,k,e,p",
            "source_densities": {
                "m": "V",
                "h": "E",
                "k": "E*V",
                "e": "V^2-E^2/2",
                "p": "E^2/2",
            },
            "transport": "(s/x)^rate * constructed_patch_field(s) * ds/s",
            "bounds": "s in [1,x], x in [1,e]",
            "nonreset": True,
            "retained_raw_six_channels": ["M", "H", "K", "A", "B", "C"],
            "finite_N_section11_correction": False,
        },
        "beta_gamma_source": {
            "centers": list(CENTERS),
            "radius": RADIUS,
            "support_edges": ["49/40", "51/40", "59/40", "61/40", "69/40", "71/40"],
            "raw_beta": "exp(-1/(1-u^2)) for |u|<1; exact 0 outside",
            "gamma": "40*normalization^-1*beta(40*(x-center))",
            "gamma_x": "1600*normalization^-1*beta'(40*(x-center))",
            "gamma_xx": "64000*normalization^-1*beta''(40*(x-center))",
            "ordinary_x_derivative_binomial_factors": {"first": "40", "second": "40^2"},
            "integrals_unevaluated": True,
        },
        "branch_premises": {
            "actual_source_Z_domain": branch.get("actual_source_Z_domain", ["-1", "1"]),
            "common_box_unique": witness["c6_branch_unique_common_box"],
            "conditional_domain": True,
            "conditional_branch": True,
            "fixed_N0_exact_integer_premise": witness["fixed_N0_exact_integer_premise"],
            "numerical_margins_asserted": False,
            "interval_values_or_overlap_used": False,
        },
        "same_C6_source_map": {
            "map_system_node": 193470,
            "control_nodes": control_ids,
            "defining_map": "L*h+Q(h,invAm2)+true_Rm_defects=0",
            "same_unique_common_box_branch": same_control_map,
            "solver_rerun": False,
        },
        "definition_graph_checks": graph_checks,
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
        "P1_differs_from_m15": True,
        "physical_Q_differs_from_loop_q": True,
        "retained_construction_data": retained,
        "runtime_model_metadata": None,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "source_inspection_audit": source_audit,
        "unresolved_gates": gates,
    }


__all__ = ["SCHEMA", "EXPECTED_NODE_COUNT", "EXPECTED_NODES_DIGEST", "RECORD_HASHES", "derive"]
