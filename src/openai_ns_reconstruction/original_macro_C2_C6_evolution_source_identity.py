"""Static source identity for the original C2 to C6 macro evolution.

This consumer authenticates inherited source records and compares their ASTs.
It never imports a producer, evaluates a graph, materialises a numerical value,
or treats a Taylor enclosure as the function it encloses.  The positive result
is an exact source-function/ODE correspondence; finite-N, density, and Q_y/Z_1
seams deliberately remain outside this consumer.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction

SCHEMA = "original_macro_C2_C6_evolution_source_identity_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
POSITIVE_REPORT = "lei_ren_part1_paper_compliant_current_original_micro_evolution_source_projection.json"
POSITIVE_SHA256 = "9916b2aae49a624b5b505ba0c4f7b3cfbdf92f9a7a30d7852fb247efbd085dbd"
AFFINE_REPORT = "lei_ren_part1_paper_compliant_current_original_micro_affine_source_identity.json"
AFFINE_SHA256 = "3ba7c5829ebc1243ef63d04f092c17613cfda7e76913980e7e16b71c810fe37c"
C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz"
C6_ARCHIVE_SHA256 = "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92"
C6_DECODED_SHA256 = "f361ab39e7599f5ed770df630bccadc4056b3e7ad6c0e1de32c975963f070298"
FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    POSITIVE_REPORT: POSITIVE_SHA256,
    AFFINE_REPORT: AFFINE_SHA256,
    C6_ARCHIVE: C6_ARCHIVE_SHA256,
    "lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py": "74511ee00a655184021b0516478cffb7c8fdca830b5dc01ab44dfde2e9ad648b",
    "lei_ren_part1_paper_compliant_current_original_bridge_macro_functions_check.json": "076962a246a1e672673b132203b6b57ea3d694cec505df9593b44815bedf6f7a",
    "lei_ren_part1_paper_compliant_current_original_bridge_macro_moments.py": "ab0f25f7a9c9436f8fbaba9952d538a4c861695911d4e6596eebc219ac124b85",
    "lei_ren_part1_paper_compliant_current_original_bridge_macro_moments_check.json": "28434c4b09b11f9b0d78b5f85d26f3ed42670955194c8d793d346460f1965093",
    "lei_ren_part1_paper_compliant_comparison_point_integrals.py": "a7103469587709ff80b087b99fcf5490c8445afcd25c2794a4097b73963d98b7",
    "lei_ren_part1_paper_compliant_comparison_point_integrals_check.json": "3b2cd05bd8732acaf3b53b4f97422c55721a4ab9db9c760758b5cbdd0e31171b",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    "lei_ren_part1_paper_compliant_current_original_macro_finite_N.py": "13e60fc5968d1024552197ce06a3aa1b92b76ed073b7ab2179ceb545810ad317",
    "lei_ren_part1_paper_compliant_current_original_macro_finite_N_check.json": "0a44840b8fec9424445c512831847b119028352efe1cd9ee545ef9ecd39d7207",
    "lei_ren_part1_paper_compliant_macro_signed_integrals.py": "046b7982464eb00eb16c41c3a652c46a6ff82a1c7386fc3283dc5dd654377913",
    "lei_ren_part1_paper_compliant_current_original_micro_functions.py": "d4ef2385ff4627c91f3fd99d1c6819e5598db3e789093f02fe1dc965ec8be779",
    "lei_ren_part1_paper_compliant_current_original_micro_functions_check.json": "54d14670660ea41e90ffc25f536e65747da2a3008d8c34e790780cd2d9df8fac",
    "lei_ren_part1_paper_compliant_actual_bridge_integrals.py": "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3",
    "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json": "575fc6b5c1c78957bf9e0d178582bfcff26c167d4de037cd7aea41fc69f52145",
    "lei_ren_part1_paper_compliant_inner_bridge_profiles.py": "142945ac1dd8a5705c380e79da29ea6bd89374c18cfa50f0412f1ba70c520840",
    "lei_ren_part1_paper_compliant_inner_bridge_profiles_check.json": "ea360228f830af0dc78fdb38f43f8585cc68f6bf7b6e10b6bb3b4adce4af64ff",
    "lei_ren_part1_paper_compliant_frozen_comparison_field.py": "e38e99800be7ea5d6b142e9101e75e871c97e8dcf7a2b781927f729dcebd9bfb",
    "lei_ren_part1_paper_compliant_frozen_comparison_field_check.json": "68737b02622414f095ce8aecb6daf0dd5c222e0a4479a265cdb672249b00d163",
    "lei_ren_part1_paper_compliant_bridge_mixed_C4.py": "55569c8679beb43eaf4ba9681393a2604d34dfed1a83d66d53c6219cfa7e9187",
    "lei_ren_part1_paper_compliant_bridge_mixed_C4_check.json": "9842df322ae1d84ef79760b6bfe4a01a2d94877639a2be7130571bac19b801e8",
    "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py": "f3ad62f07b7631ccbfe67c787b4989c0d70db2ff49c96325ca3ab2be284ad87c",
    "lei_ren_part1_paper_compliant_anchored_axis_amplitude_check.json": "299b52671d3ff05b91968b1b78e9ab57804abb60815b32d6d469c4d3e47a00a9",
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py": "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json": "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    "lei_ren_part1_paper_compliant_core_transfer.json": "4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497",
    "lei_ren_part1_paper_compliant_K1_ledger.json": "38cd8e76c0d5f5519b510d9cb0502b8915446ffa1c6338d3a3a98c4aca2cdf7d",
    "lei_ren_part1_paper_compliant_core_uniform_bounds.py": "23d38d188d44e3169d9107bcf2eba201f6b866f4e408c4c83a84aa2f9d1420d7",
    "lei_ren_part1_paper_compliant_core_uniform_bounds.json": "7415cba5b9abc486b2078c450bd0f2e1b7798d1921f5863c00c04526620844e7",
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.json": "e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400",
    "lei_ren_part1_paper_compliant_current_original_Rm_native_transport_pressure_sources.json": "7c88052348a78ab7e26728f492296fb57858cead1794b9d8cb6e57ad873a7b9b",
}
DECODED_RECORD_HASHES = {C6_ARCHIVE: C6_DECODED_SHA256}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


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
        gates.append("plain_source_text_missing:" + name)
        return ""
    if name == C6_ARCHIVE:
        if record.get("encoding") != "gzip" or record.get("decoded_sha256") != C6_DECODED_SHA256:
            gates.append("encoded_archive_contract_mismatch:" + name)
            return ""
        try:
            raw = base64.b64decode(record.get("raw_bytes_base64", ""), validate=True)
            if hashlib.sha256(raw).hexdigest() != C6_ARCHIVE_SHA256:
                gates.append("archive_raw_hash_mismatch:" + name)
                return ""
            decoded = gzip.decompress(raw)
            if hashlib.sha256(decoded).hexdigest() != C6_DECODED_SHA256 or decoded != value.encode():
                gates.append("archive_decoded_bytes_mismatch:" + name)
                return ""
        except Exception as error:
            gates.append("archive_decode_error:" + type(error).__name__)
            return ""
        return value
    if record.get("encoding") is not None:
        gates.append("unexpected_encoded_record:" + name)
    if hashlib.sha256(value.encode()).hexdigest() != RECORD_HASHES[name]:
        gates.append("source_hash_mismatch:" + name)
    return value


def _tree(value, name, gates):
    try:
        return ast.parse(value, filename=name)
    except Exception as error:
        gates.append("ast_parse_error:" + name + ":" + type(error).__name__)
        return ast.Module(body=[], type_ignores=[])


def _json(value, name, gates):
    try:
        return json.loads(value, parse_int=lambda raw: int(raw), parse_float=float)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}


def _fn(tree, name, cls=None):
    if cls is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls:
                return next((child for child in node.body
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name), None)
        return None
    return next((node for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name), None)


def _norm(value):
    return "".join(value.split())


def _assign(fn, target, expression):
    if fn is None:
        return False
    want_target = _norm(target)
    want_value = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    for node in ast.walk(fn):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        if any(_norm(ast.unparse(item)) == want_target for item in targets):
            if ast.dump(node.value, include_attributes=False) == want_value:
                return True
    return False


def _has(fn, *needles):
    if fn is None:
        return False
    compact = _norm(ast.unparse(fn))
    return all(_norm(needle) in compact for needle in needles)


def _source_has(text, *needles):
    compact = _norm(text)
    return all(_norm(needle) in compact for needle in needles)


def _returns_mapping(fn, *keys):
    if fn is None:
        return False
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Dict):
            names = {key.value for key in node.value.keys if isinstance(key, ast.Constant)}
            if set(keys).issubset(names):
                return True
        if (isinstance(node, ast.Return) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name) and node.value.func.id == "dict"):
            names = {item.arg for item in node.value.keywords if item.arg is not None}
            if set(keys).issubset(names):
                return True
    return False


def _return_expression(fn, key, expression):
    if fn is None:
        return False
    expected = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    for node in ast.walk(fn):
        if not isinstance(node, ast.Return):
            continue
        if isinstance(node.value, ast.Dict):
            pairs = [(item.id if isinstance(item, ast.Name) else item.value if isinstance(item, ast.Constant) else None, value)
                for item, value in zip(node.value.keys, node.value.values)]
        elif isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == "dict":
            pairs = [(item.arg, item.value) for item in node.value.keywords]
        else:
            continue
        for name, value in pairs:
            if name == key and ast.dump(value, include_attributes=False) == expected:
                return True
    return False


def _receipt(receipts, name, required, forbidden, gates):
    value = receipts.get(name, {})
    good = value.get("all_passed") is True
    if not good:
        gates.append("receipt_not_accepted:" + name)
    for path in required:
        ok = _get(value, path) is True
        if not ok:
            gates.append("receipt_missing:" + name + ":" + path)
        good = good and ok
    for path in forbidden:
        if _get(value, path) is True:
            gates.append("receipt_open_seam_claimed:" + name + ":" + path)
            good = False
    return good


def _exact_tuple(value):
    if not isinstance(value, list) or len(value) != 4 or any(type(v) is not int for v in value):
        return None
    sign, mantissa, exponent, bits = value
    if sign not in (0, 1) or mantissa < 0 or bits != mantissa.bit_length() or abs(exponent) > 10000:
        return None
    result = Fraction(mantissa << exponent) if exponent >= 0 else Fraction(mantissa, 1 << (-exponent))
    return -result if sign else result


def _endpoint_ports(texts, nodes, gates):
    """Attach exact source graph functions, never an enclosure row, to C2."""
    archive = _json(texts[C6_ARCHIVE], C6_ARCHIVE, gates)
    if not isinstance(archive, Mapping):
        gates.append("C6_archive_mapping_missing")
        archive = {}
    saved = archive.get("exact_graph_nodes", [])
    bindings = archive.get("prescribed_bridge_C6_bindings", {})
    if not isinstance(bindings, Mapping):
        gates.append("C6_snapshot_mapping_missing")
        bindings = {}
    source_nodes = nodes if isinstance(nodes, list) else []
    checks = {"saved_C6_source_family": archive.get("source_family") == FAMILY,
        "same_73190_node_prefix": isinstance(saved, list) and len(saved) == 73190
            and isinstance(nodes, list) and nodes[:73190] == saved}
    def node(ref):
        return source_nodes[ref] if type(ref) is int and 0 <= ref < len(source_nodes) and isinstance(source_nodes[ref], Mapping) else {}
    def endpoints(endpoint, original, phase, order):
        if not isinstance(endpoint, list) or not isinstance(original, list) or len(endpoint) != order+1 or len(original) != order+1:
            return False
        for ref, expression in zip(endpoint, original):
            item = node(ref)
            value = node(item.get("value"))
            if (item.get("operation") != "function_substitution" or item.get("expression") != expression
                    or item.get("variable") != phase or item.get("Z_independent_substitution") is not True
                    or value.get("operation") != "exact_rational" or value.get("numerator") != 2
                    or value.get("denominator") != 1):
                return False
        return True
    actual = bindings.get("actual_micro_end", {})
    memories = bindings.get("actual_moment_end", {})
    comparison = bindings.get("comparison_micro_end", {})
    micro = bindings.get("micro", {})
    bar = bindings.get("comparison_micro", {})
    for name in ("phi", "U"):
        checks["actual_endpoint_"+name] = endpoints(actual.get(name), _get(micro, "actual."+name), micro.get("phase"), 6)
        checks["comparison_endpoint_"+name] = endpoints(comparison.get(name), bar.get(name), bar.get("phase"), 7)
    for name in ("H", "M", "K", "A", "B", "C"):
        checks["actual_memory_endpoint_"+name] = endpoints(memories.get(name), _get(micro, "moments."+name), micro.get("phase"), 6)
        checks["comparison_memory_endpoint_"+name] = endpoints(_get(comparison, "moments."+name), _get(bar, "moments."+name), bar.get("phase"), 7)
    hb = bindings.get("hb")
    checks["bridge_width_source_handle"] = (node(hb).get("operation") == "current_original_source_parameter"
        and node(hb).get("name") == "hbB" and node(hb).get("quantity_is_exact_original_parameter_not_range_endpoint") is True)
    switches = [i for i, item in enumerate(source_nodes) if isinstance(item, Mapping)
        and item.get("operation") == "current_original_source_parameter" and item.get("name") == "hbS"]
    checks["bridge_switch_source_handles_distinct"] = len(switches) == 1 and switches[0] != hb
    checks["macro_uses_bridge_width_handle"] = _get(bindings, "macro.chi") == hb
    passed = all(checks.values())
    if not passed:
        gates.extend("exact_endpoint_graph_port_missing:"+name for name, ok in checks.items() if not ok)
    return dict(schema="original_C2_macro_exact_inlet_function_ports_v1", graph_ports_valid=passed,
        checks=checks, graph_nodes_digest=EXPECTED_NODES_DIGEST, source_archive=C6_ARCHIVE,
        source_archive_sha256=C6_ARCHIVE_SHA256, source_prefix_nodes=73190,
        micro_projection_premise_sha256=POSITIVE_SHA256,
        C2_macro_initial_exact_fields={name: list(actual.get(name, [])) for name in ("phi", "U")},
        C2_macro_initial_exact_histories={name: list(memories.get(name, [])) for name in ("H", "M", "K", "A", "B", "C")},
        frozen_comparison_exact_functions=comparison, bridge_width_node=hb, switch_width_nodes=switches,
        bridge_switch_width_identity_proved=False,
        C2_attachment="The accepted exact C2 micro ODE projection identifies these symbolic s=2 endpoints with the C2 exact micro functions. New macro initial-function ports reference those functions; old second_exit rows are separate covers.",
        C2_exponential_route="q=phi_core/barphi2; phi_core*exp(ell_in+ell_macro)=phi_micro_end*exp(ell_macro). The micro prefix is used once.",
        enclosure_row_equality_claimed=False, numerical_values_computed=False,
        ordinary_Z_orders_admitted=[0, 1, 2], all_saved_C6_endpoint_rows_retained=True)


def _domain(texts, trees, gates):
    transfer = _json(texts["lei_ren_part1_paper_compliant_core_transfer.json"], "core_transfer", gates)
    ledger = _json(texts["lei_ren_part1_paper_compliant_K1_ledger.json"], "K1_ledger", gates)
    bounds = _json(texts["lei_ren_part1_paper_compliant_core_uniform_bounds.json"], "core_uniform_bounds", gates)
    lower = _exact_tuple(_get(transfer, "logLambda.lower_exact_mpf_tuple"))
    phi_lower = _exact_tuple(_get(bounds, "actual_Phi_lower.lower_exact_mpf_tuple"))
    double_h = _exact_tuple(_get(ledger, "decreasing_K_smallness_bounds.double_h.upper_exact_mpf_tuple"))
    run = _fn(trees["lei_ren_part1_paper_compliant_core_uniform_bounds.py"], "run")
    c6 = trees["lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"]
    c6_init = _fn(c6, "__init__", "CurrentOriginalPrescribedBridgeC6Functions")
    c6_actual = _fn(c6, "actual_fields")
    c6_hist = _fn(c6, "histories")
    ast_bound = all((_assign(run, "model_lower", "1-qmax/2+qmax*qmax/12-qmax**3/144"),
        _assign(run, "correction_value", "correction_norm*embedding(c,h,radius,0,0)"),
        _assign(run, "phi_lower", "model_lower-correction_value"),
        _has(run, "ifendpoints(phi_lower)[0]<=0", "full_real_axis_normalized_swirl_positive=True")))
    result = {
        "logLambda_lower_exact_gt_one": lower is not None and lower > 1,
        "double_h_upper_exact_le_one_hundredth": double_h is not None and double_h <= Fraction(1, 100),
        "positive_actual_Phi_lower_exact": phi_lower is not None and phi_lower > 0,
        "core_uniform_bounds_ast": ast_bound,
        "smallness_threshold_is_one_hundredth": _get(ledger, "smallness_thresholds.double_h") == ".01",
        "smallness_check_passed": _get(ledger, "smallness_checks.double_h") is True,
        "same_hb_epsilon_definition": _get(ledger, "parameter_definition.epsilon_b") == "h_b, the identical scalar",
        "positive_width_not_materialized": _get(ledger, "positive_width_not_materialized") is True,
        "axis_domain": _get(bounds, "real_axial_domain") == ["-1", "1"],
        "radial_domain": _get(bounds, "scaled_radial_domain") == ["0", "4.1"],
        "analytic_positive_field_not_sampled": _get(bounds, "full_real_axis_normalized_swirl_positive") is True and _get(bounds, "unique_analytic_field_bounded_not_sampled") is True and _get(bounds, "finite_coefficient_state_used") is False,
        "macro_length_exact_source": _has(c6_init, "self.Y=g.sub(g.unary('log',c(100)),g.unary('log',Ra))", "self.length=length=g.sub(self.Y,g.mul(c(2),hb))") and lower is not None and lower > 1 and double_h is not None and double_h <= Fraction(1, 100),
        "macro_radius_exact_source": _has(c6_init, "Rmacro=g.mul(Ra,g.unary('exp',g.add(g.mul(c(2),hb),g.mul(length,m))))", "self.macro=dict(phase=m,domain=[0,1],radius=Rmacro"),
        "frozen_barphi2_independent_of_macro_m": _has(c6_init, "barphi2,barU2=(self.comparison_micro_end[name]fornamein('phi','U'))", "equilibrium=dict(H=barphi2,M=barU2"),
        "macro_triangular_field_ode": _has(c6_actual, "phi=a.mul(initial['phi'],a.exp(logratio))", "a.sub(initial['U'],a.integral", "full_drive=a.add") and _has(c6_hist, "forname,rateinRATES.items()", "a.integral(body,t,g.zero,phase"),
        "macro_ode_uniqueness_premises": _has(c6_actual, "same positive actual comparison phi", "g.neg(jacobian)") and _has(c6_hist, "original_incoming_memory_retained=True", "source_measure_applied_once=True"),
    }
    result["proved"] = all(result.values())
    if not result["proved"]:
        gates.append("macro_domain_positivity_or_continuation_missing")
    result["derivation"] = {
        "interval": "m in [0,1], Z in [-1,1]; t=length*m=log(R/R0)",
        "length": "logLambda>1, Ra=4*exp(-logLambda), Y=log(25)+logLambda>1 and 2hb<=1/100 imply length=Y-2hb>99/100",
        "frozen_comparison": "barphi2=positive micro comparison phi(s=2) is independent of m. Each comparison moment is equilibrium+(inlet-equilibrium)*exp(-rate*length*m).",
        "coefficients": "The common source is analytic in Z. F0>0, barphi2>0 and L>=1/2 protect the comparison-direction denominators; finite sums of these functions and exp(-j*t) are continuous on the compact domain.",
        "radius": "0<R0=Ra*exp(2hb)<=R(m)<=100; the exact radius and length are not hupper covers",
        "fields": "The angular primitive is finite; phi=phi_micro_end*exp(ell_macro)>0. The axial equation is a forced integral with continuous coefficients and no division by actual phi.",
        "uniqueness": "The triangular linear field equations and then six linear history equations have unique full-interval solutions; analytic parameter dependence permits ordinary Z differentiation through order 2.",
        "bounds": "No macro numerical field, norm or enclosure radius is computed by this theorem. Old enclosure errors remain unchanged.",
    }
    return result


def _reports(texts, gates):
    parent = _json(texts[PARENT_REPORT], PARENT_REPORT, gates)
    positive = _json(texts[POSITIVE_REPORT], POSITIVE_REPORT, gates)
    affine = _json(texts[AFFINE_REPORT], AFFINE_REPORT, gates)
    if _get(parent, "source_family") != FAMILY:
        gates.append("parent_source_family_changed")
    if _get(positive, "source_family") != FAMILY:
        gates.append("positive_source_family_changed")
    if _get(positive, "parent_report_sha256") != PARENT_SHA256:
        gates.append("positive_parent_pin_changed")
    component = _get(positive, "component_result", {})
    required = {
        "passed": True, "source_only": True, "full_micro_background_evolution_projection_proved": True,
        "full_micro_own_history_projection_proved": True, "numeric_values_computed": False,
        "old_operators_replayed": False, "full_goal_complete": False,
    }
    for key, value in required.items():
        if component.get(key) is not value:
            gates.append("accepted_micro_component_changed:" + key)
    for key in ("same_exact_ODE_and_uniqueness", "same_six_history_rates_and_initials", "all_old_prefixes_are_enclosures", "F0_same_anchored_function"):
        if _get(component, "sublemmas." + key) is not True:
            gates.append("accepted_micro_sublemma_missing:" + key)
    for key in ("positivity_existence_uniqueness", "c2_geometry_and_chi", "complete_two_chart_control_memory", "formal_comparison_exact_global_functions", "complete_direction_operands"):
        if _get(component, "ast_witnesses." + key) is not True:
            gates.append("accepted_micro_source_premise_missing:" + key)
    for key in ("source_domain", "positive_denominators", "hb_positivity", "uniqueness"):
        if not isinstance(_get(component, "domain_evidence." + key), str) or not _get(component, "domain_evidence." + key):
            gates.append("accepted_micro_domain_premise_missing:" + key)
    scope = positive.get("scope", {})
    for key, value in {"ordinary_Z_orders": [0, 1, 2], "same_density_certificate": False,
        "finite_N_history_identity": False, "physical_Q_y_Z1": False,
        "global_five_moment_closure": False, "full_goal_complete": False}.items():
        if scope.get(key) != value:
            gates.append("accepted_micro_scope_changed:" + key)
    if _get(component, "nodes_digest") != EXPECTED_NODES_DIGEST or _get(scope, "nodes_digest") != EXPECTED_NODES_DIGEST:
        gates.append("accepted_micro_graph_changed")
    if _get(component, "accepted_source_stage_report") != AFFINE_REPORT or _get(component, "accepted_source_stage_sha256") != AFFINE_SHA256:
        gates.append("accepted_affine_pin_changed")
    affine_scope = _get(affine, "scope", {})
    if _get(affine, "nodes_digest") not in (None, EXPECTED_NODES_DIGEST):
        gates.append("affine_graph_changed")
    if affine_scope and affine_scope.get("same_density_certificate", affine_scope.get("offset_same_density_certificate", _get(affine_scope, "offset_flags.same_density_certificate"))) is not False:
        gates.append("affine_density_scope_changed")
    if _get(positive, "retained_construction_data.source_family") not in (None, FAMILY):
        gates.append("accepted_construction_family_changed")
    return parent, positive, component


def _receipts(texts, gates):
    receipts = {name: _json(text, name, gates) for name, text in texts.items() if name.endswith("_check.json")}
    specs = {
        "lei_ren_part1_paper_compliant_current_original_bridge_macro_functions_check.json": (
            ("original_frozen_macro_complete_nonlinear_FV_integral_functions_installed", "independent_complete_original_flow.passed", "genuine_native_source_functions.passed", "genuine_native_source_functions.micro_inlet_not_reset", "genuine_native_source_functions.full_bridge_moments_first_switch_and_whole_axis_remain_open", "guards.passed"), ()),
        "lei_ren_part1_paper_compliant_current_original_bridge_macro_moments_check.json": (
            ("original_frozen_macro_actual_six_moment_function_integrals_installed", "exact_signed_polynomial_kernels.passed", "independent_complete_nonlinear_actual_histories.passed", "genuine_native_source_histories.passed", "genuine_native_source_histories.B_source_phi_squared_with_rate2_preserved", "genuine_native_source_histories.micro_functions_whole_Z_switches_controls_and_recursion_remain_open", "guards.passed"), ()),
        "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json": (
            ("current_original_true_comparison_bridge_C7_functions_installed", "current_original_true_prescribed_micro_macro_C6_functions_installed", "current_original_true_prescribed_bridge_Q_inertial_p2_C5_functions_installed", "checks.actual_source.passed", "checks.actual_source.preserved_core_C6_before_actual_C7_comparison_source", "checks.actual_source.G_higher_rows_follow_original_anchored_FTC", "checks.actual_source.F0_is_Z_only_and_all_phase_bodies_have_no_free_core_rho", "checks.actual_source.original_P0_C6_preserved_and_executable_Z7_projection", "checks.original_direction.passed", "checks.original_direction.physical_swirl_F0_squared_dressing_and_external_scale_once", "checks.signed_transport.passed", "checks.signed_transport.nonunit_jacobian_chi_radius_Pstar_and_F0", "checks.signed_transport.full_axial_minus_sign_and_pressure_swirl_physical_scales_once", "checks.own_six_histories.passed", "checks.own_six_histories.all_nonzero_incoming_histories_and_sixth_derivatives_retained", "checks.own_six_histories.rates_H2_M1_K2_A1_B2_C1_and_B_forcing_not_halved", "checks.own_six_histories.nonunit_physical_jacobian_applied_once"), ("current_original_actual_modified_bridge_C6_functions_installed", "current_original_prescribed_bridge_C6_numeric_containment_installed")),
        "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json": (("full_finite_width_signed_bridge_integral_enclosures_available", "actual_own_six_moment_feedback_enclosures_available", "whole_real_axial_domain_R100_checked"), ("actual_bridge_mixed4_feedback_installed", "exact_production_point_parameters_selected")),
        "lei_ren_part1_paper_compliant_comparison_point_integrals_check.json": (("original_9_23_namespace_bound", "actual_atoms_and_formal_width_remainders_preserved"), ("original_signed_actual_bridge_integrals_resolved",)),
        "lei_ren_part1_paper_compliant_inner_bridge_profiles_check.json": (("angular_log_sign_bound_to_accepted_same_source_comparison_theorem", "source_bound_smoothed_comparison_axial6_enclosures_available", "actual_prescribed_shear_bridge_axial5_enclosures_available", "actual_bridge_radial_recovery_axial4_enclosures_available"), ("actual_bridge_radial_mixed4_certified",)),
        "lei_ren_part1_paper_compliant_frozen_comparison_field_check.json": (("functional_identities.passed", "actual_core_exit_five_primitive_axial5_enclosures_available", "frozen_comparison_all_five_primitives_available", "frozen_comparison_profile_mixed4_available"), ("actual_smooth_comparison_installed",)),
        "lei_ren_part1_paper_compliant_bridge_mixed_C4_check.json": (("phase1_actual_two_sided_physical_grids_identical", "phase2_actual_histories_retained", "actual_core_exit_R100_histories_and_original_P0_retained", "R100_physical_zeroth_radial_source_rows_identical", "exact_comparison_integral_and_six_own_moment_histories_bound", "core_bridge_and_R100_local_functional_joins_certified"), ("full_inner_interfaces_certified",)),
        "lei_ren_part1_paper_compliant_anchored_axis_amplitude_check.json": (("six_poles_and_original_unique_anchor_certified", "actual_anchored_G_resolved_with_directed_error", "positive_F0_retained_without_exp_logF0", "source_identities.passed", "source_identities.exact_rational_pole_decomposition", "source_identities.exact_simple_pole_residue_identity"), ("full_point_physical_field_evaluation",)),
        "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json": (("current_original_n1_true_core_fixed_point_functions_and_enclosure_bridge_installed", "source_functions_and_directed_enclosures_are_distinct", "checks.actual_oracle.exact_parameter_and_native_P0_identity.passed", "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.unique_H_root_in_original_bracket", "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.amplitude_uses_same_core", "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.amplitude_uses_same_nonlinear_rebuild", "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.anchor_leaf_uses_this_exact_source", "checks.actual_oracle.exact_parameter_and_native_P0_identity.interval_overlap_not_used_as_a_function_identity"), ("uniform_original_core_point_numerical_accuracy_claimed",)),
        "lei_ren_part1_paper_compliant_current_original_micro_functions_check.json": (("current_original_Ra_micro_hydrated_whole_source_functions_installed", "independent_original_micro_ODE.passed", "independent_original_micro_ODE.exact_core_memory_and_three_micro_control_seams_checked", "actual_live_source.passed", "actual_live_source.comparison_owner_hydrated_without_ancestor_constructor", "actual_live_source.finite_N_boundary_and_global_moment_closure_still_unsupplied", "original_source_bindings.original_micro_controls.first_chi.verified", "original_source_bindings.original_micro_controls.second_and_macro_chi.verified", "original_source_bindings.original_micro_controls.comparison_direction.verified", "original_source_bindings.original_micro_controls.actual_F_source_integral.verified", "original_source_bindings.original_micro_controls.actual_V_source_integral.verified", "original_source_bindings.recovered_norms_are_upper_bounds_not_source_values"), ("actual_R110_source_histories_closed",)),
    }
    status = {}
    for name, (required, forbidden) in specs.items():
        status[name] = _receipt(receipts, name, required, forbidden, gates)
    return status


def _witnesses(texts, trees, receipts, domain):
    macro = trees["lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py"]
    moments = trees["lei_ren_part1_paper_compliant_current_original_bridge_macro_moments.py"]
    comp = trees["lei_ren_part1_paper_compliant_comparison_point_integrals.py"]
    c6 = trees["lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"]
    micro = trees["lei_ren_part1_paper_compliant_current_original_micro_functions.py"]
    actual = trees["lei_ren_part1_paper_compliant_actual_bridge_integrals.py"]
    inner = trees["lei_ren_part1_paper_compliant_inner_bridge_profiles.py"]
    frozen = trees["lei_ren_part1_paper_compliant_frozen_comparison_field.py"]
    signed = trees["lei_ren_part1_paper_compliant_macro_signed_integrals.py"]
    finite = trees["lei_ren_part1_paper_compliant_current_original_macro_finite_N.py"]
    anchored = trees["lei_ren_part1_paper_compliant_anchored_axis_amplitude.py"]
    bounds = trees["lei_ren_part1_paper_compliant_core_uniform_bounds.py"]
    mf_init = _fn(macro, "__init__", "MacroFlow")
    mf_set = _fn(macro, "set_sources", "MacroFlow")
    mf_geo = _fn(macro, "geometry", "MacroFlow")
    mf_eval = _fn(macro, "evaluate", "MacroFlow")
    owner = _fn(macro, "owner", "OriginalBridgeMacroFunctions")
    mm_eval = _fn(moments, "evaluate", "ActualMacroMoments")
    mm_owner = _fn(moments, "owner", "OriginalBridgeMacroMoments")
    c6_init = _fn(c6, "__init__", "CurrentOriginalPrescribedBridgeC6Functions")
    c6_dirs = _fn(c6, "comparison_directions")
    c6_actual = _fn(c6, "actual_fields")
    c6_hist = _fn(c6, "histories")
    c2_dir = _fn(inner, "direction")
    c2_dress = _fn(frozen, "dress")
    micro_init = _fn(micro, "__init__", "MicroFunctions")
    micro_prefix = _fn(micro, "prefix", "MicroFunctions")
    micro_hist = _fn(micro, "histories", "MicroFunctions")
    micro_eval = _fn(micro, "evaluate", "MicroFunctions")
    actual_prepare = _fn(actual, "prepare", "CompliantActualBridgeIntegrals")
    actual_contrib = _fn(actual, "contributions", "CompliantActualBridgeIntegrals")
    actual_packet = _fn(actual, "packet", "CompliantActualBridgeIntegrals")
    actual_macro_moments = _fn(actual, "macro_moments", "CompliantActualBridgeIntegrals")
    comp_micro = _fn(comp, "micro", "CompliantComparisonPointIntegrals")
    comp_macro = _fn(comp, "macro", "CompliantComparisonPointIntegrals")
    width_packet = _fn(comp, "packet", "WidthPolynomial")
    frozen_dress = _fn(frozen, "dress")
    signed_f0 = _fn(signed, "F0_squared_log_enclosure")
    signed_target = _fn(signed, "_target_moments")
    signed_modes = _fn(signed, "_mode_rows")
    finite_init = _fn(finite, "__init__", "WholeMacroFlow")
    finite_geo = _fn(finite, "geometry", "WholeMacroFlow")
    anchored_eval = _fn(anchored, "evaluate", "CompliantAnchoredAxisAmplitude")
    out = {}
    out["c2_comparison_macro_zero_is_micro_two"] = all((_assign(owner, "entry", "comp['comparison_point_packets'][label]['macro'][0]"), _assign(actual_prepare, "end", "self.comparison.micro(Z,2,root)"), _has(comp_macro, "end=self.micro(Z,2,root)", "initial_smoothing_history_retained=True"), _has(c6_init, "self.comparison_micro_end", "a.at(barphi,s,c(2))", "a.at(barU,s,c(2))", "self.comparison_macro")))
    out["c2_second_exit_exact_function_ports"] = all((_has(owner, "frame=original['packets'][label]", "frame['core']['actual_phi_axial5']", "frame['core']['actual_raw_V_axial5']", "inlet=frame['second_exit']", "inlet_terms(inlet['signed_angular_integral_terms'])", "inlet_terms(inlet['signed_axial_integral_terms'][part],part)"), _has(mm_owner, "packet=self.fields.records['actual_bridge_integrals']['packets'][label]['second_exit']", "raw=packet['actual_own_six_moments_axial5']", "ActualMacroMoments(flow,inlet)"), _has(actual_prepare, "end=self.comparison.micro(Z,2,root)", "phi2=value_enclosure(end['phi'])", "moment2=named_moments"), _has(micro_eval, "chart=chart", "source_interval_functions_not_endpoint_hulls=True"), _returns_mapping(actual_packet, "actual_log_F_over_inlet_axial5", "actual_phi_axial5", "actual_raw_V_axial5", "actual_own_six_moments_axial5")))
    out["c2_second_exit_exact_function_ports"] = bool(out["c2_second_exit_exact_function_ports"]
        and _assign(actual_prepare, "data", "self.comparison.inlet(Z,root)")
        and _has(_fn(actual, "report", "CompliantActualBridgeIntegrals"), "second_exit=self.packet(Z,2,'second',root)")
        and _assign(c6_init, "self.actual_micro_end", "{name:b.at(actual[name],s,c(2)) for name in ('phi','U')}")
        and _assign(c6_init, "self.actual_moment_end", "{name:b.at(jet,s,c(2)) for name,jet in actualJ.items()}"))
    out["c2_c6_exact_parameter_ports"] = all((_assign(mf_init, "self.logs", "tuple((c.mpf(v) for v in (logh,logP2,logF02,logRa,0)))"), _assign(mf_init, "self.h", "self.factor((1,0,0,0,0))"), _assign(owner, "logh", "previous.read_interval(c,original['source_log_hb_enclosure'])"), _assign(owner, "logF02", "2*previous.read_interval(c,amplitude['logF0'])"), _assign(owner, "logP2", "2*previous.read_interval(c,core['logPstar'])"), _assign(owner, "logRa", "c.ln(4)-logLambda"), _has(c6_init, "self.hb=hb=inlet_owner.previous.lower.parameters['hb']", "self.Ra=Ra=inlet_owner.radius", "Pstar=inlet_owner.Pstar"), _has(micro_init, "self.Ra=flow.factor((0,0,0,1,0))")))
    out["comparison_alpha_rho_U0_and_common_core"] = all((_has(c6_init, "alpha=g.sub(g.one,sigma(g.sub(u,g.one)))", "rho=lambdaq:g.mul(c(4),g.unary('exp',g.mul(hb,q)))", "phi0=a.at(self.core7['Phi'],self.core.rho,c(4))", "U0=a.add(self.core7['U0'],a.scale(a.at(self.core7['Psi'],self.core.rho,c(4)),self.core.parameters['epsilon'])", "barphi=a.exp(logbar)", "barU=a.add(U0,a.integral"), _has(micro_eval, "chi=f.scalar(1-sig)+f.h*sig", "chi=f.h", "R=self.Ra*self.series.scalar_series(1,s)[0]", "physical_log_radius_over100=f.scalar(-f.Y)+f.h*s"), _has(c6_init, "barJ=histories(a,barphi,barU,self.initial7,s,hb,'comparison_micro')")))
    out["old_prefixes_are_enclosures_only"] = all((_source_has(texts["lei_ren_part1_paper_compliant_comparison_point_integrals.py"], "exact_source_representation='sum_{j=0}^2 hb^j a_j(Z)+hb^3 e(Z); sum_k |e_k|*weight^k<=remainder'", "width_not_materialized=True"), _source_has(texts["lei_ren_part1_paper_compliant_current_original_micro_functions_check.json"], "recovered_norms_are_upper_bounds_not_source_values"), _has(comp_macro, "moments[name]=decay*initial+(1-decay)*targets[name]", "exact_log_radius='y=2hb+fraction*(log(100/Ra)-2hb)'"), _has(actual_macro_moments, "signed_linear_volterra_terms", "quadratic_remainder"), _source_has(texts["lei_ren_part1_paper_compliant_actual_bridge_integrals.py"], "source_width_not_materialized=True", "cap_used_as_field_value=False")))
    out["c2_prefix_exponent_and_volterra_source"] = all((_has(micro_prefix, "angular.append(f.scale(d,-self.Rcover*mass*c.mpf('.5')))", "V=f.add(self.V0,*deltaV.values())", "full_exponential_error=error"), _has(micro_hist, "rhs=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),B=f.multiply(phi,phi),C=f.multiply(phi,phi))", "result[name]=f.add(f.scale(self.initial[name],decay),f.scale(rhs[name],mass))"), _has(micro_eval, "Vy=f.add(*[f.scale(f.multiply(quotient,row),-chi*f.factor((0,int(part=='pressure'),int(part=='swirl'),0,0)))forpart,rowindrives.items()])", "phi_y=f.scale(f.multiply(a,at['phi']),-c.mpf('.5'))", "V_y=Vy"), _has(actual_contrib, "angular.append(self.source_term(d['D_over_R']*(-Rmicro*mass/2),width_power))", "row=-(p['micros'][0]['quotient']*d['drive_'+name])*Rmicro**powerR[name]*mass")))
    out["no_double_inlet_factorization"] = all((_has(mf_eval, "ell=self.add(self.ell_in,*[self.scale(row,self.h*I[1,l]*-c.mpf('.5'))forl,rowinenumerate(self.d)])", "delta_phi=self.multiply(self.phi0,self.add(ell,self.error_rows(remainder)))", "product=self.multiply(self.q,self.drive[part][j])"), _has(owner, "phi0/phi.truncate(5)", "quotient='Phi_core_exit/barphi2; full micro log prefix enters exponential separately'"), _has(c6_actual, "phi=a.mul(initial['phi'],a.exp(logratio))", "ratio=a.div(phi,a.trim(barphi),'same positive actual comparison phi')"), _has(c6_init, "actualmacro=actual_fields(b,macroknown,self.actual_micro_end,m,length,hb,Rmacro,F0,Pstar,barphi2,'macro')")))
    out["c2_c6_direction_equations"] = all((_returns_mapping(c2_dir, "D_over_R", "drive_hydro", "drive_pressure", "drive_swirl"), _assign(c2_dir, "W", "1-z*m*(1-delta)-d*derivative(m)"), _assign(c2_dir, "angular", "h*(1-delta/2)-z*derivative(h)*((1-delta)/2)-d*derivative(k)+z*k*(2*delta-1)"), _return_expression(c2_dir, "drive_pressure", "(z*pressure*(2*(1+delta))-d*derivative(pressure))/(2*L)"), _returns_mapping(c6_dirs, "D_over_R", "drive_hydro", "drive_pressure", "drive_swirl"), _assign(c6_dirs, "F2", "b.mul(F0,F0)"), _assign(c6_dirs, "dressed_swirl", "divide(a.add(a.neg(a.scale(a.mul(de,a.mul(z,a.trim(B))),2)),a.mul(d,a.shift(B)),a.scale(a.mul(a.add(one,de),a.mul(z,a.trim(C))),2),a.neg(a.mul(d,a.shift(C)))))"), _has(c6_dirs, "swirl=a.div(dressed_swirl,a.trim(F2),'same strictly positive full F0 squared')")))
    out["full_physical_rhs_and_swirl_germ"] = all((_assign(c6_actual, "full_drive", "a.add(a.scale(known['drive_hydro'],radius),a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar)),a.mul(a.scale(a.mul(F0,F0),g.mul(radius,radius)),known['drive_swirl']))"), _assign(c6_actual, "Vy", "a.scale(a.mul(ratio,full_drive),g.neg(chi))"), _has(c6_actual, "a.sub(initial['U'],a.integral(a.scale(a.at(Vy,phase,v),g.neg(jacobian))"), _has(c2_dress, "jet*IntervalTaylor(c,[ratios[k]/math.factorial(k)for k in range(jet.order+1)])"), _has(frozen_dress, "jet*IntervalTaylor(c,[ratios[k]/math.factorial(k)for k in range(jet.order+1)])"), _has(owner, "F0_squared_ratios=ratios(2)", "normalized_swirl_rows_include_true_F0_squared_derivatives_once=True")))
    out["c2_mode_rows_match_c6_macro_modes"] = all((
        _return_expression(signed_target, "MTH", "phi"), _return_expression(signed_target, "MZ", "value"),
        _return_expression(signed_target, "MTHZ", "phi*value"),
        _return_expression(signed_target, "MZT", "{'axial':value*value,'swirl':phi*phi/2}"),
        _return_expression(signed_target, "MP", "phi*phi"),
        _assign(signed_modes, "target", "_target_moments(phi,value)"),
        _assign(signed_modes, "delta", "_sub_moments(inlet,target)"),
        _assign(signed_modes, "delta_one", "{MTH:_zero(phi),MZ:delta[MZ],MTHZ:_zero(phi),MZT:{'axial':delta[MZT]['axial'],'swirl':_zero(phi)},MP:delta[MP]}"),
        _assign(signed_modes, "delta_two", "{MTH:delta[MTH],MZ:_zero(phi),MTHZ:delta[MTHZ],MZT:{'axial':_zero(phi),'swirl':delta[MZT]['swirl']},MP:_zero(phi)}"),
        _assign(signed_modes, "base", "_mode_direction(bridge,z,inputs,phi,value,target)"),
        _assign(signed_modes, "one", "_mode_direction(bridge,z,inputs,phi,value,_add_moments(target,delta_one))"),
        _assign(signed_modes, "two", "_mode_direction(bridge,z,inputs,phi,value,_add_moments(target,delta_two))"),
        _return_expression(signed_modes, "D_over_R", "[base['D_over_R'],one['D_over_R']-base['D_over_R'],two['D_over_R']-base['D_over_R']]"),
        _return_expression(signed_modes, "drive_hydro", "[base['drive_hydro'],one['drive_hydro']-base['drive_hydro'],_zero(base['drive_hydro'])]"),
        _return_expression(signed_modes, "drive_pressure", "[base['drive_pressure'],_zero(base['drive_pressure']),_zero(base['drive_pressure'])]"),
        _return_expression(signed_modes, "drive_swirl", "[base['drive_swirl'],one['drive_swirl']-base['drive_swirl'],two['drive_swirl']-base['drive_swirl']]"),
        _assign(owner, "known", "modes._mode_rows(SimpleNamespace(ctx=c,delta=delta),Z,inputs,phi,V,bridge.named_moments(moment_rows))"),
        _assign(c6_init, "equilibrium", "dict(H=barphi2,M=barU2,K=a.mul(barphi2,barU2),A=a.mul(barU2,barU2),B=a.scale(a.mul(barphi2,barphi2),'1/2'),C=a.mul(barphi2,barphi2))"),
        _assign(c6_init, "macroknown", "comparison_directions(g,barphi2,barU2,macroJ,self.core7['F0'],self.P0,self.core.parameters['delta'])"),
        _has(c6_init, "macroJ[name]=a.add(a.scale(self.comparison_micro_end['moments'][name],decay),a.scale(equilibrium[name],g.sub(g.one,decay)))")))
    out["anchored_F0_same_function_and_ordinary_Z_germ"] = all((_assign(anchored_eval, "logF0", "-self.core.logC-self.core.Lambda*actual_G"), _has(anchored_eval, "G_gradient_taylor_coefficients=self.gradient_jets(z)", "exp_logF0_not_materialized=True", "positive_F0_retained_formally=True"), _has(owner, "def ratios(multiplier)", "rows=[c.mpf(1)]", "math.factorial(n)"), _source_has(texts["lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"], "same_implicit_native_P0=True", "common_P0_kept_independent_of_cumulative_pressure")))
    out["macro_modes_exp_minus_j_t_and_measure_conversion"] = all((_has(c6_init, "theta=g.unary('exp',g.neg(g.mul(length,m)))", "decay=g.unary('exp',g.mul(c(-rate),length,m))", "macroJ[name]=a.add", "actualmacro=actual_fields(b,macroknown,self.actual_micro_end,m,length,hb,Rmacro,F0,Pstar,barphi2,'macro')", "actualmacroJ=histories(b,actualmacro['phi'],actualmacro['U'],self.actual_moment_end,m,length,'actual_macro')"), _has(signed, "R*hydro", "R*Pstar^2*pressure", "R^2*F0^2*swirl"), _has(finite_geo, "R0=self.factor((0,0,0,1,0))*self.series.scalar_series(2,self.c.mpf(1))[0]", "S=self.ordinary_cover(self.length)*t")))
    out["six_history_rates_nonreset_and_derivatives"] = all((_assign(moments, "RATES", "dict(H=2,M=1,K=2,A=1,B=2,C=1)"), _has(mm_eval, "source_polys=dict(H=phi.scale(2),M=V,K=(phi*V).scale(2),A=V*V,B=phi*phi,C=phi*phi)", "inherited=f.scale(self.inlet[name],decay)", "derivatives[name]=f.add(source_end[name],f.scale(actual,-rate))"), _has(c6_hist, "decay=g.unary('exp',g.mul(g.constant(-rate),jacobian,phase))", "kernel=g.unary('exp',g.mul(g.constant(-rate),jacobian,g.sub(phase,t)))", "a.scale(initial[name],decay)", "actual_six_history_rate=rate", "original_incoming_memory_retained=True"), _has(micro_hist, "true_incoming_decay=decay", "true_positive_Volterra_mass=mass", "original_core_memory_retained=True"), _has(actual_macro_moments, "forname,rateinRATES.items()", "baseline[name]=start[name]*decay+targets[name]*(1-decay)", "feedback[name]")))
    out["frozen_comparison_positivity_and_ode_domain"] = bool(domain["proved"] and _has(c6_dirs, "same positive 2*(1-delta*Z^2)", "same positive 1-delta*Z^2", "same strictly positive full comparison F0*phi") and _has(c6_actual, "same positive actual comparison phi") and _source_has(texts["lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json"], "modified_finite_N_source_and_quantitative_C6_not_claimed"))
    out["finite_N_is_geometry_only"] = all((_has(finite_init, "self.original_flow=flow", "self.length=self.scalar(self.Y)-self.h*2"), _has(finite_geo, "R0=self.factor((0,0,0,1,0))*self.series.scalar_series(2,self.c.mpf(1))[0]", "R=self.scalar(100)ifterminalelseR0ifemptyelseself.radius(t)"), _source_has(texts["lei_ren_part1_paper_compliant_current_original_macro_finite_N_check.json"], "complete_local_micro_exit_Rm_operator_not_real_inlet_or_closure")))
    out["hupper_is_cover_separate_from_exact_geometry"] = all((
        _assign(mf_init, "self.hupper", "self.small_upper(self.h)"),
        _source_has(texts["lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py"], "An outward tail bound only; never a source value or selected width"),
        _assign(mf_geo, "S", "q*(self.Y-2*self.c.mpf([0,ep(self.hupper)[1]]))"),
        _assign(mf_geo, "R0", "self.factor((0,0,0,1,0),self.c.mpf([0,ep(2*self.hupper)[1]]))"),
        _assign(c6_init, "Rmacro", "g.mul(Ra,g.unary('exp',g.add(g.mul(c(2),hb),g.mul(length,m))))"),
        _has(c6_init, "self.macro=dict(phase=m,domain=[0,1],radius=Rmacro,physical_measure=length,chi=hb")))
    accepted_micro = _json(texts[POSITIVE_REPORT], POSITIVE_REPORT, [])
    out["no_fixture_as_actual_source"] = all(_get(accepted_micro, "component_result.fixture_limits." + name + ".actual_source_admission") is False
        for name in ("core_log_fixture", "varying_direction_fixture", "controls_fixture"))
    return out


def _result(gates, witnesses, receipt_status, nodes, actual_digest, reports, domain, ports):
    graph_ok = isinstance(nodes, list) and len(nodes) == EXPECTED_NODE_COUNT and actual_digest == EXPECTED_NODES_DIGEST
    background_keys = ("c2_comparison_macro_zero_is_micro_two", "c2_second_exit_exact_function_ports", "c2_c6_exact_parameter_ports", "comparison_alpha_rho_U0_and_common_core", "old_prefixes_are_enclosures_only", "c2_prefix_exponent_and_volterra_source", "no_double_inlet_factorization", "c2_c6_direction_equations", "c2_mode_rows_match_c6_macro_modes", "full_physical_rhs_and_swirl_germ", "anchored_F0_same_function_and_ordinary_Z_germ", "macro_modes_exp_minus_j_t_and_measure_conversion", "frozen_comparison_positivity_and_ode_domain", "finite_N_is_geometry_only", "no_fixture_as_actual_source")
    history_keys = ("six_history_rates_nonreset_and_derivatives",)
    receipts_ok = bool(receipt_status) and all(receipt_status.values())
    background = bool(graph_ok and receipts_ok and ports.get("graph_ports_valid") is True and not gates and all(witnesses.get(k) is True for k in background_keys))
    own = bool(background and witnesses.get("six_history_rates_nonreset_and_derivatives") is True)
    return {
        "schema": SCHEMA, "source_only": True, "passed": bool(background and own),
        "nodes_digest": actual_digest, "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "full_macro_background_evolution_projection_proved": background,
        "full_macro_own_history_projection_proved": own,
        "full_macro_evolution_identity_proved": bool(background and own),
        "full_macro_own_history_identity_proved": own, "full_macro_proved": bool(background and own),
        "ordinary_Z_orders": [0, 1, 2], "ordinary_orders": [0, 1, 2],
        "time_order_recursion_claimed": False, "timeorder_recursion_claimed": False,
        "same_density_certificate": False, "physical_Q_y_Z1": False,
        "finite_N_history_identity": False, "finite_N_sc_over_2": False,
        "global_five_moment_closure": False, "numeric_values_computed": False,
        "old_operators_replayed": False, "mathematical_operator_replayed": False,
        "full_goal_complete": False, "ast_witnesses": witnesses,
        "accepted_receipts": {k: True for k, v in receipt_status.items() if v},
        "receipt_status": receipt_status, "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_SHA256, "accepted_source_stage_report": POSITIVE_REPORT,
        "accepted_source_stage_sha256": POSITIVE_SHA256, "source_family": _get(reports[1], "source_family", {}),
        "projection_kind": "exact original source-function identity with C2 enclosure remainders retained; no finite-prefix value equality",
        "physical_swirl_identity": "F0_base^2*C2_drive_swirl_germ == [F0^2*C6_normalized_swirl]_germ; bare higher-Z swirl equality is not claimed",
        "parameter_ports": {"hbB": ports.get("bridge_width_node"), "hbS": ports.get("switch_width_nodes"), "bridge_switch_width_identity_proved": False, "same_Ra": witnesses.get("c2_c6_exact_parameter_ports") is True, "same_Pstar": witnesses.get("c2_c6_exact_parameter_ports") is True, "hb_symbols_kept_distinct": _get(ports, "checks.bridge_switch_source_handles_distinct") is True, "distinct_handles_do_not_prove_numerical_inequality": True, "differential_conversion": "d/dm=length*d/dt with t=log(R/R0)"},
        "exact_macro_inlet_function_ports": dict(ports, C2_source_projection_admitted=background),
        "endpoint_ports": {"comparison_macro_zero": "comparison micro source at s=2", "actual_second_exit": "exact C2 source-function port; serialized frame rows remain enclosures", "actual_micro_end": "C6 actual_fields/history evaluated at s=2 as a symbolic endpoint", "exact_function_definitions": {"ell_star": "ell_*=-1/2 integral_0^2 hb*chi(s)*R(s)*D_over_R(s) ds", "U_star": "U_*=V_core-integral_0^2 hb*chi(s)*(phi/barphi)*full_drive(s) ds", "J_star": "J_* = exp(-rate*hb*2) J_core + hb integral_0^2 exp(-rate*hb*(2-s))*rho_J(s) ds", "macro_factorization": "phi_core*exp(ell_in+ell_macro)=phi_*exp(ell_macro); q remains phi_core/barphi2 in the C2 route", "saved_rows": "actual frame second_exit rows are enclosure evidence only"}, "endpoint_rows_are_values": False},
        "domain_evidence": domain,
        "sublemmas": {"comparison_macro_zero_is_micro_two": witnesses.get("c2_comparison_macro_zero_is_micro_two", False), "actual_second_exit_function_port": witnesses.get("c2_second_exit_exact_function_ports", False), "same_analytic_comparison_core": witnesses.get("comparison_alpha_rho_U0_and_common_core", False), "same_full_physical_rhs": witnesses.get("full_physical_rhs_and_swirl_germ", False), "same_ordinary_F0_germ": witnesses.get("anchored_F0_same_function_and_ordinary_Z_germ", False), "same_six_history_rates_and_initials": own, "whole_interval_positivity_existence_uniqueness": witnesses.get("frozen_comparison_positivity_and_ode_domain", False), "all_old_prefixes_are_enclosures": witnesses.get("old_prefixes_are_enclosures_only", False)},
        "scope_flags": {"same_density_certificate": False, "physical_Q_y_Z1": False, "finite_N_sc_over_2": False, "global_five_moment_closure": False, "time_order_recursion": False, "full_goal_complete": False},
        "fixture_limits": {name: {"finite_fixture_only": True, "actual_source_admission": False} for name in ("core_log_fixture", "varying_direction_fixture", "controls_fixture")},
        "unresolved_gates": sorted(set(gates)) if gates or not (background and own) else [],
        "potential_gates": ["same_density_certificate", "physical_Q_y_Z1", "finite_N_history_identity", "global_five_moment_closure"],
        "record_hashes": dict(RECORD_HASHES), "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "source_inspection_audit": {"source_records_only": True, "static_source_AST_graph_and_saved_bound_inspection_only": True, "nodes_hashed_only": False, "source_graph_substitutions_checked": True, "new_source_consumer_calls": 1, "accepted_ancestor_operator_calls": 0, "producer_provider_builder_checker_replayed": False, "quadrature_or_norm_calls": 0, "numeric_field_evaluated": False},
        "worker_role": "luna_max_worker", "runtime_model_metadata": None,
    }


def derive(nodes, source_records, *, expected_nodes_digest):
    """Return the fail-closed source-function projection for the macro chart."""
    gates = []
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append("expected_nodes_digest_mismatch")
    try:
        actual_digest = _digest(nodes)
    except Exception as error:
        actual_digest = ""
        gates.append("nodes_digest_error:" + type(error).__name__)
    if actual_digest != EXPECTED_NODES_DIGEST:
        gates.append("nodes_digest_mismatch")
    if not isinstance(nodes, list) or len(nodes) != EXPECTED_NODE_COUNT:
        gates.append("node_count_mismatch")
    records = source_records if isinstance(source_records, Mapping) else {}
    texts = {name: _text(records, name, gates) for name in RECORD_HASHES}
    reports = _reports(texts, gates)
    receipt_status = _receipts(texts, gates)
    trees = {name: _tree(value, name, gates) for name, value in texts.items() if name.endswith(".py")}
    domain = _domain(texts, trees, gates)
    ports = _endpoint_ports(texts, nodes, gates)
    witnesses = _witnesses(texts, trees, receipt_status, domain)
    witnesses["exact_endpoint_function_graph_ports"] = ports["graph_ports_valid"]
    for name, ok in witnesses.items():
        if ok is not True:
            gates.append("ast_witness_missing:" + name)
    return _result(gates, witnesses, receipt_status, nodes, actual_digest, reports, domain, ports)
