"""Fresh Z-only first-switch background and C1 arithmetic.

This module reads accepted JSON and source text only.  It never imports or
executes an accepted producer, owner, checker, recurrence, quadrature, or
norm implementation.  The switch radius and all first-switch changes are
handled as directed Z-only majorants; no physical y-derivative claim is made.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath import make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json"
PARENT_SHA256 = "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c"
BACKGROUND = BASE + "lei_ren_part1_paper_compliant_current_original_macro_background_I1.json"
BACKGROUND_SHA256 = "5a77d3d8b9b294516cc0be05ed83667e2710d64ff8427c7e917b8aed3379aede"
R100_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
R100_SOURCE_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
R100_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.py"
R100_CHECK_SHA256 = "8ad5e2edc1daba05763b44dd031ac1671d74f30b1733fc80570873e2ca12089a"
R100_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
R100_CACHE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
FIRST_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.py"
FIRST_SOURCE_SHA256 = "a98753313fb1689f784c7a56cf09dc06df83e0b52cb738b5d4531c19f5402496"
FIRST_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.py"
FIRST_CHECK_SHA256 = "72675cf41c852d411535d191fb5677ed29736a9b78a8785f2f38b66edcde8cff"
FIRST_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.json"
FIRST_CACHE_SHA256 = "5ea141f4ee1befaf33d6d2d2b9ad871477b72ee3dbbfc6b1e55b82766278a53d"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
C6_SOURCE_SHA256 = "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605"
C6_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.py"
C6_CHECK_SHA256 = "c5e118a1d2691a878d535c5b0b39cf9e627994c0f1f6c8d36139b09d46fac836"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz"
C6_CACHE_SHA256 = "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
BRIDGE_SOURCE_SHA256 = "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3"
NATIVE_EXT = "src/openai_ns_reconstruction/original_macro_native_Z4_extension.py"
NATIVE_EXT_SHA256 = "ac0c29264fb791c52559035151a8cffc0653bb7189902d88b7e5de09575485ce"
R100_GEOMETRY_SOURCE_RECEIPT = (
    BASE + "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions_check.json"
)
R100_GEOMETRY_SOURCE_RECEIPT_SHA256 = "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158"
K1_LEDGER = BASE + "lei_ren_part1_paper_compliant_K1_ledger.json"
K1_LEDGER_SHA256 = "38cd8e76c0d5f5519b510d9cb0502b8915446ffa1c6338d3a3a98c4aca2cdf7d"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
SIGMA_SOURCE_SHA256 = "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7"
MODULE = "src/openai_ns_reconstruction/original_first_switch_background_I1_numeric.py"

ORDERS = (0, 1, 2, 3, 4)
C1_ORDERS = (0, 1)
RATES = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}


# This is intentionally a tuple.  Every source guard is AST-only and the
# expected digest binds the actual accepted function body, rather than this
# module's recipe text.
SOURCE_GUARDS = (
    (C6_SOURCE, "densities", (
        "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",
        "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi)",
    ), "2678e3d1926d5a81e43a9d040f8da8e9298b0693f9f6bc4c91b64803aa3ca4e1"),
    (C6_SOURCE, "histories", (
        "for name,rate in RATES.items():",
        "original_incoming_memory_retained=True,source_measure_applied_once=True",
    ), "f868b63668f662a5c7a07f44a09f47c7b1724b19fb474518a0f62fea072469c7"),
    (C6_SOURCE, "comparison_directions", (
        "F2=b.mul(F0,F0);B,C=(b.mul(F2,J[name]) for name in ('B','C'))",
        "D=a.div(a.add(a.neg(W),a.div(angular",
    ), "c849f0b83bc39effda3d3df88984cae566928e9b52c9a0455a96b30014d58ed6"),
    (R100_SOURCE, "first_switch_fields", (
        "g.mul(g.constant('-1/2'),hb,hb,radius)",
        "original_first_switch_full_current_radius_direction=True",
        "source_measure_applied_once=True",
    ), "ac140c0eaee3e5e3130d8bdd11fb5d27a3a25430233cd41d107dcd11ca23e6ac"),
    (R100_SOURCE, "full_drive", (
        "a.scale(known['drive_hydro'],radius)",
        "a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar))",
        "g.mul(radius,radius)",
    ), "e94f55c35f33828b20c98f65c5c5acbb20ab065d13d9ce0984e96247d4126bee"),
    (R100_SOURCE, "comparison_at_radius", (
        "eq=dict(H=phi,M=U,K=b.mul(phi,U)",
        "decay=g.unary('exp',g.mul(c(-rate),distance))",
        "J[name]=b.add(b.scale(initial[name],decay)",
    ), "1e4d5e9eadc05389a83bdcfad00a8ac4ae3c90c8a87159844bd4655759e7374b"),
    (R100_SOURCE, "add_chart", (
        "J=previous.histories(a,actual['phi'],actual['U'],incoming",
        "known_comparison_direction=known",
    ), "50194431065c07b29545a5f53c997f731f0ee79ee7e23e32641d789fae8ab2a9"),
    (R100_SOURCE, "parameter_geometry_guard", (
        "shared is sw.history.upstream",
        "switch_width_is_exact_source_reboxing=True",
        "Identical actual K1 width definition",
    ), "51bb38764a44ca9d7686b48cac5e172b6e29fdfc524e1c925bee1dcb0c708e01"),
    (R100_SOURCE, "__init__", (
        "self.hb=hb=self.native_parameter('hbS'",
        "R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))",
        "self.add_chart('first_switch'",
    ), "d875ec07da501f43644b9c89e2f0001976031e485c6637d3d018e6375cc796ef"),
    (FIRST_SOURCE, "fixed_comparison_denominator_binding", (
        "q=core_phi/barphi",
        "comparison_phi_constant_through_original_frozen_macro=True",
    ), "5e8a80582e5c6301bf36681f8a06335dfe60e87e060f860e1949b508fca7d9c4"),
)


def _require(value, message):
    if not value:
        raise ValueError(message)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _finite(value):
    try:
        return value._mpi_[0][0] in (0, 1) and value._mpi_[1][0] in (0, 1)
    except AttributeError:
        return False


def _record(value):
    _require(hasattr(value, "_mpi_"), "Directed interval required")
    return dict(
        lower_exact_mpf_tuple=list(value._mpi_[0]),
        upper_exact_mpf_tuple=list(value._mpi_[1]),
        lower_display=nstr(make_mpf(value._mpi_[0]), 24),
        upper_display=nstr(make_mpf(value._mpi_[1]), 24),
        directed=True,
    )


def _read_record(c, value):
    if isinstance(value, Mapping) and "value" in value and "lower_exact_mpf_tuple" not in value:
        value = value["value"]
    _require(isinstance(value, Mapping), "Directed exact record required")
    try:
        low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
        high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Malformed directed exact record") from exc
    out = c.mpf([low, high])
    _require(_finite(out), "Non-finite directed record")
    return out


def _lo(value):
    return make_mpf(value._mpi_[0])


def _hi(value):
    return make_mpf(value._mpi_[1])


def _point(c, value):
    if hasattr(value, "_mpi_"):
        low, high = make_mpf(value._mpi_[0]), make_mpf(value._mpi_[1])
        _require(low == high, "Point interval required")
        value = low
    return c.mpf([value, value])


def _unsigned_tuple(endpoint):
    exact = tuple(endpoint)
    if exact[0] == 0:
        return exact
    return (0, exact[1], exact[2], exact[3])


def _abs_upper(c, value):
    """Absolute upper endpoint, with no global-context abs or arithmetic."""
    low = c.mpf([make_mpf(_unsigned_tuple(value._mpi_[0]))] * 2)
    high = c.mpf([make_mpf(_unsigned_tuple(value._mpi_[1]))] * 2)
    bound = max(_hi(low), _hi(high))
    return c.mpf([bound, bound])


def _sym(c, value):
    upper = _abs_upper(c, value) if hasattr(value, "_mpi_") else _point(c, value)
    zero = c.mpf(0)
    neg = zero - upper
    return c.mpf([_lo(neg), _hi(upper)])


def _max_cap(c, values):
    values = list(values)
    _require(values, "Nonempty cap collection required")
    cap = c.mpf(0)
    for value in values:
        candidate = _abs_upper(c, value)
        if _hi(candidate) > _hi(cap):
            cap = candidate
    return cap


def _sum(c, values):
    out = c.mpf(0)
    for value in values:
        out = out + value
    return out


def _prod(c, left, right, order):
    return _sum(c, (c.mpf(math.comb(order, j)) * left[j] * right[order - j]
                    for j in range(order + 1)))


def _scale(c, scalar, value):
    return c.mpf(scalar) * value


def _neg(c, value):
    return c.mpf(0) - value


def _rows(c, value, count, name):
    _require(isinstance(value, (list, tuple)) and len(value) >= count,
             name + " must provide every requested ordinary row")
    found = {}
    for index, row in enumerate(value):
        if isinstance(row, Mapping) and "ordinary_Z_order" in row:
            order = row["ordinary_Z_order"]
            item = row.get("value", row)
        else:
            order, item = index, row
        _require(type(order) is int, name + " has a non-integer ordinary order")
        if order not in ORDERS:
            continue
        _require(order not in found, name + " repeats an ordinary order")
        found[order] = _read_record(c, item)
    _require(set(found) >= set(range(count)), name + " has missing ordinary rows")
    return [found[index] for index in range(count)]


def _reciprocal(c, values, count):
    _require(_lo(values[0]) > 0, "Positive reciprocal denominator required")
    one = _point(c, 1)
    result = [one / values[0]]
    for order in range(1, count):
        tail = _sum(c, (c.mpf(math.comb(order, j)) * values[j] * result[order - j]
                        for j in range(1, order + 1)))
        result.append(_neg(c, tail / values[0]))
    return result


def _contains(outer, inner):
    return _lo(outer) <= _lo(inner) and _hi(inner) <= _hi(outer)


def _assert_le(left, right, message):
    _require(_hi(left) <= _lo(right), message)


def _assert_lt(left, right, message):
    _require(_hi(left) < _lo(right), message)


def _log_cap(c, value):
    cap = _abs_upper(c, value)
    _require(_lo(cap) > 0, "Positive cap required before logarithm")
    return c.log(cap)


def _pack_rows(c, rows, kind=None):
    result = []
    for order, value in enumerate(rows):
        item = dict(ordinary_Z_order=order, value=_record(value))
        if kind is not None:
            item["kind"] = kind
        result.append(item)
    return result


def _jsonify(value):
    if hasattr(value, "_mpi_"):
        return _record(value)
    if isinstance(value, Mapping):
        return {key: _jsonify(child) for key, child in value.items()}
    if isinstance(value, tuple):
        return [_jsonify(child) for child in value]
    if isinstance(value, list):
        return [_jsonify(child) for child in value]
    return value


def _pin(root, path, expected):
    actual_path = root / path
    _require(actual_path.is_file(), "Pinned source is missing: " + path)
    actual = _sha(actual_path)
    _require(actual == expected, "Pinned source changed: " + path)
    return dict(path=path, expected_sha256=expected, passed=True, source_executed=False)


def _find_function(tree, name, path):
    found = [node for node in ast.walk(tree)
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
    _require(len(found) == 1, "Expected exactly one " + name + " in " + path)
    return found[0]


def _source_guard(root):
    output = []
    for path, function, snippets, expected_ast in SOURCE_GUARDS:
        text = (root / path).read_text(encoding="utf-8")
        tree = ast.parse(text, filename=path)
        node = _find_function(tree, function, path)
        segment = ast.get_source_segment(text, node) or ""
        compact = "".join(segment.split())
        actual_ast = hashlib.sha256(
            ast.dump(node, include_attributes=False).encode()
        ).hexdigest()
        _require(actual_ast == expected_ast, "Accepted AST changed: " + path + ":" + function)
        for snippet in snippets:
            _require("".join(snippet.split()) in compact,
                     "Accepted formula missing: " + path + ":" + function)
        output.append(dict(
            passed=True,
            path=path,
            function=function,
            required_text=list(snippets),
            AST_sha256=actual_ast,
            source_executed=False,
        ))
    return tuple(output)


def _load_parent(root, parent):
    path = root / PARENT
    _require(path.is_file(), "Pinned finite-N macro parent is missing")
    _require(_sha(path) == PARENT_SHA256, "Pinned finite-N macro parent changed")
    disk = json.loads(path.read_bytes())
    _require(parent == disk, "The exact accepted finite-N macro parent mapping is required")
    return disk


def _load_json(root, path):
    actual_path = root / path
    _require(actual_path.is_file(), "Pinned JSON source is missing: " + path)
    value = json.loads(actual_path.read_bytes())
    _require(isinstance(value, Mapping), "Pinned JSON source must be a mapping: " + path)
    return value


def _load_json_gz(root, path):
    actual_path = root / path
    _require(actual_path.is_file(), "Pinned gzip JSON source is missing: " + path)
    value = json.loads(gzip.open(actual_path, "rb").read())
    _require(isinstance(value, Mapping), "Pinned gzip JSON source must be a mapping: " + path)
    return value


def _exact_record_tuples(value):
    if isinstance(value, Mapping) and "value" in value and "lower_exact_mpf_tuple" not in value:
        value = value["value"]
    _require(isinstance(value, Mapping), "Exact record mapping required")
    return (tuple(value["lower_exact_mpf_tuple"]),
            tuple(value["upper_exact_mpf_tuple"]))


def _row_exact_tuple(row, name):
    _require(isinstance(row, Mapping), name + " row must be a mapping")
    order = row.get("ordinary_Z_order")
    _require(type(order) is int, name + " row ordinary order must be an integer")
    return order, _exact_record_tuples(row.get("value", row))


def _verify_extension_native_rows(extension):
    """Verify accepted macro prefixes and every new row against exact records."""
    native = extension.get("native_fields")
    accepted = extension.get("accepted_macro_source_rows")
    proofs = extension.get("new_field_row4_proofs")
    _require(isinstance(native, Mapping) and isinstance(accepted, Mapping),
             "Native extension row mappings are required")
    _require(isinstance(proofs, Mapping), "Native extension row-four proofs are required")
    fields = ("phi", "U", "M", "A", "B", "C")
    native_tuples = {}
    for name in fields:
        rows = native.get(name)
        _require(isinstance(rows, list) and len(rows) == 5,
                 "Native extension must provide five " + name + " rows")
        parsed = [_row_exact_tuple(row, "native " + name) for row in rows]
        _require([order for order, unused in parsed] == list(ORDERS),
                 "Native " + name + " rows are not ordered 0..4")
        native_tuples[name] = [record for unused, record in parsed]
    prefix = {}
    for name in ("phi", "U", "M"):
        rows = accepted.get(name)
        _require(isinstance(rows, list) and len(rows) == 4,
                 "Accepted macro source prefix must provide four " + name + " rows")
        parsed = [_row_exact_tuple(row, "accepted macro " + name) for row in rows]
        _require([order for order, unused in parsed] == list(range(4)),
                 "Accepted macro " + name + " rows are not ordered 0..3")
        values = [record for unused, record in parsed]
        _require(values == native_tuples[name][:4],
                 "Accepted macro " + name + " prefix does not equal native extension")
        prefix[name] = dict(orders=list(range(4)), exact_tuple_matches=True)
    row4 = {}
    for name in ("phi", "U", "M"):
        proof = proofs.get(name)
        _require(isinstance(proof, Mapping), "Missing row-four proof for " + name)
        _require(_exact_record_tuples(proof.get("value")) == native_tuples[name][4],
                 "Row-four proof does not equal native " + name)
        row4[name] = dict(value_exact_tuple_matches=True,
                          ordinary_Z_order=4)
    expected_density = {
        "A": ("U^2", "U^2", 1),
        "B": ("phi^2", "phi^2/2", 2),
        "C": ("phi^2", "phi^2", 1),
    }
    for name, (density, normalized, rate) in expected_density.items():
        proof = proofs.get(name)
        _require(isinstance(proof, Mapping), "Missing row-four proof for " + name)
        _require(proof.get("density") == density and
                 proof.get("normalized_density") == normalized and
                 proof.get("rate") == rate,
                 "Native " + name + " density/rate source receipt changed")
        rows = proof.get("positive_kernel_hull_rows")
        _require(isinstance(rows, list) and len(rows) == 5,
                 "Native " + name + " row-four hull must have five rows")
        hull = [_exact_record_tuples(row) for row in rows]
        _require(hull == native_tuples[name],
                 "Native " + name + " hull rows do not equal native fields")
        _require(proof.get("ordinary_product_binomial") is True and
                 proof.get("source_executed") is False,
                 "Native " + name + " arithmetic receipt changed")
        row4[name] = dict(
            positive_kernel_hull_exact_tuple_matches=True,
            ordinary_product_binomial=True,
            source_executed=False,
            ordinary_Z_orders=list(ORDERS),
        )
    return dict(
        passed=True,
        fields=list(fields),
        native_orders=list(ORDERS),
        accepted_macro_prefix=prefix,
        row4_proofs=row4,
        exact_tuple_comparison=True,
        source_executed=False,
    )


def _find_sigma_source_node(nodes):
    found = []
    for index, node in enumerate(nodes):
        if not isinstance(node, Mapping):
            continue
        if node.get("operation") != "exact_original_source_flat_sigma":
            continue
        definition = node.get("definition", "")
        if ("0 for x<=0" not in definition or "1 for x>=1" not in definition or
                "inside" not in definition):
            continue
        binding = node.get("source_binding")
        if not isinstance(binding, Mapping):
            continue
        required = ("module", "module_sha256", "function", "AST_sha256")
        if not all(key in binding for key in required):
            continue
        if (binding.get("module_sha256") != SIGMA_SOURCE_SHA256 or
                node.get("smooth_flat_endpoints") is not True or
                node.get("Z_independent_argument") is not True):
            continue
        found.append((index, node))
    _require(found, "Expected an exact original flat-sigma source node")
    index, node = found[0]
    _require(node["source_binding"]["module_sha256"] == SIGMA_SOURCE_SHA256,
             "Pinned sigma source digest changed")
    return index, node


def _verify_parent_drive(c, parent):
    drive = parent["actual_macro_background_I1_numeric"]["full_signed_axial_drive"]
    _require(drive.get("F0_squared_cancellation_exact") is True,
             "F0 squared cancellation receipt is required")
    _require(drive.get("no_inverse_F0") is True,
             "No-inverse-F0 receipt is required")
    units = "p0=P0/Pstar^2; Pstar^2 applied exactly once in pressure drive"
    _require(drive.get("pressure_source_units") == units,
             "Pressure source units receipt changed")
    pstar = _read_record(c, drive.get("Pstar_squared"))
    return dict(
        passed=True,
        F0_squared_cancellation_exact=True,
        no_inverse_F0=True,
        pressure_source_units=units,
        Pstar_squared_source=(
            "parent.actual_macro_background_I1_numeric.full_signed_axial_drive.Pstar_squared"
        ),
        Pstar_squared=_record(pstar),
        exact_log_F0=drive.get("exact_log_F0"),
        gaxis_source=drive.get("gaxis_source"),
        source_executed=False,
    )


def _build_geometry(root, parent, c, background, ea):
    r100_cache = _load_json_gz(root, R100_CACHE)
    first_cache = _load_json(root, FIRST_CACHE)
    bindings = r100_cache["prescribed_R100_Rm_C6_bindings"]
    certificate = bindings["geometry_certificate"]
    for key in (
        "shared_upstream_is_same_object",
        "bridge_logh_is_actual_upstream_object",
        "switch_width_is_exact_source_reboxing",
        "all_four_native_T_source_tuples_equal",
        "actual_native_T_and_logC_match_common_source",
        "bridge_hbB_switch_hbS_kept_distinct",
        "directed_T_positive",
        "directed_reference_gap_positive",
        "directed_post_length_positive",
        "bounds_are_not_defining_field_parameters",
    ):
        _require(certificate.get(key) is True,
                 "R100 geometry certificate missing " + key)
    graph = r100_cache["exact_graph_nodes"]
    hbB_handle = bindings["bridge_hb"]
    hbS_handle = bindings["hb"]
    _require(type(hbB_handle) is int and type(hbS_handle) is int and
             hbB_handle != hbS_handle, "Distinct hbB/hbS handles are required")
    hb_nodes = {index: graph[index] for index in (hbB_handle, hbS_handle)}
    _require(hb_nodes[hbB_handle].get("name") == "hbB" and
             hb_nodes[hbS_handle].get("name") == "hbS",
             "R100 hb handles do not identify hbB and hbS")
    hb_receipt = {}
    for handle, name in ((hbB_handle, "hbB"), (hbS_handle, "hbS")):
        node = hb_nodes[handle]
        _require(node.get("quantity_is_exact_original_parameter_not_range_endpoint") is True and
                 node.get("Z_independent") is True,
                 name + " source ownership receipt changed")
        _require(node.get("current_source_receipt") ==
                 "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions_check.json" and
                 node.get("current_source_receipt_sha256") == R100_GEOMETRY_SOURCE_RECEIPT_SHA256,
                 name + " current-source receipt changed")
        hb_receipt[name] = dict(
            graph_handle=handle,
            operation=node.get("operation"),
            name=name,
            definition=node.get("definition"),
            current_source_receipt=node.get("current_source_receipt"),
            current_source_receipt_sha256=node.get("current_source_receipt_sha256"),
            quantity_is_exact_original_parameter_not_range_endpoint=True,
            Z_independent=True,
        )
    ledger_hash = certificate.get("shared_K1_ledger_sha256")
    _require(ledger_hash == K1_LEDGER_SHA256,
             "R100 shared K1 ledger digest changed")
    _require(parent.get("input_hashes", {}).get(
        "lei_ren_part1_paper_compliant_K1_ledger.json") == ledger_hash,
             "Parent does not pin the shared K1 ledger")
    _require(r100_cache.get("input_hashes", {}).get(
        "lei_ren_part1_paper_compliant_K1_ledger.json") == ledger_hash,
             "R100 cache does not pin the shared K1 ledger")
    ledger_receipt = dict(
        path=K1_LEDGER,
        expected_sha256=ledger_hash,
        parent_input_hash_verified=True,
        R100_cache_input_hash_verified=True,
        shared_K1_ledger_sha256=ledger_hash,
        source_executed=False,
    )
    first_bindings = first_cache["source_control_bindings"]
    identities = first_bindings.get("identities", {})
    for key in ("first_angular_shear_is_unweighted_tinyD", "sig_source_binding",
                "complement_source_binding", "tinyD_source_binding"):
        _require(identities.get(key) is True,
                 "First-switch source identity missing " + key)
    _require(first_bindings.get("first_angular_weight") == "1" and
             first_bindings.get("first_axial_weight") == "1-sigma(s)",
             "First-switch angular/axial source weights changed")
    first_packet = first_cache["packets"]["0"][0]
    first_eval = first_packet["function_evaluation"]
    first_geometry = first_eval["geometry"]
    _require(first_packet.get("source_frame") == "0" and
             first_eval.get("true_sigma_function_used") is True,
             "First-switch source packet is not source-owned")
    _require(first_geometry.get("phase") == [0, 1] and
             first_geometry.get("exact_source") == "R=100*exp(hb*s); dy=hb*ds",
             "First-switch source geometry changed")
    _require(set(first_eval.get("ordinary_Z_orders", [])) >= set(range(6)),
             "First-switch source packet lacks its ordinary-Z rows")
    first_chart = bindings["charts"]["first_switch"]
    _require(first_chart.get("domain") == [0, 1] and
             first_chart.get("domain")[1] == 1 and
             first_chart.get("physical_measure") == hbS_handle,
             "R100 first-switch domain or endpoint changed")
    sigma_index, sigma_node = _find_sigma_source_node(graph)
    sigma_binding = sigma_node["source_binding"]
    sigma_receipt = dict(
        graph_node_index=sigma_index,
        operation=sigma_node["operation"],
        definition=sigma_node["definition"],
        source_binding=dict(sigma_binding),
        smooth_flat_endpoints=True,
        Z_independent=True,
        domain=[0, 1],
        range_enclosure="0<=sigma(s)<=1",
        endpoint_values={"s=0": 0, "s=1": 1},
        source_owned=True,
        source_executed=False,
    )
    geometric = background["geometric_measure"]
    loghbS = _read_record(c, geometric["loghb"])
    # The accepted geometry certificate proves this is one shared upstream
    # object; do not manufacture an independent second width record.
    loghb_upstream = loghbS
    threshold_hb = c.log(c.mpf(1) / c.mpf(1000))
    _assert_lt(loghbS, threshold_hb, "First-switch width must satisfy hbS < 0.001")
    r = c.power(c.mpf(2), -1024)
    log_r = c.log(r)
    _assert_lt(loghbS + c.log(c.mpf(4)), log_r,
               "Width is not small enough for normalized K enclosure")
    exp_half = c.exp(r / 2)
    _assert_lt(exp_half, c.mpf(1) + r,
               "Fresh normalized K exponential fixture does not close")
    _assert_lt(c.mpf(100) * c.exp(c.mpf("0.001")), c.mpf(101),
               "R <= 101 radius majorant failed")
    one = c.mpf(1)
    K_upper = one + r
    K = c.mpf([_lo(one), _hi(K_upper)])
    logW = c.log(c.mpf(10000)) + loghbS
    _require(ea["comparison_source_jet_proof"].get("radius_and_alpha_Z_independent") is True,
             "Accepted comparison alpha/radius Z-independence proof required")
    rho = ea["domain"].get("rho")
    _require(rho == ["4", "41/10"], "Accepted rho domain changed")
    sigma = ea["sigma_derivative_bound"]
    _require(sigma.get("passed") is True and
             sigma.get("original_sigma_not_its_sample_or_cap") is True,
             "Accepted original sigma derivative receipt required")
    return dict(
        loghbS=loghbS,
        loghb_upstream=loghb_upstream,
        loghb_shared_object=True,
        loghb_upstream_not_independent=True,
        loghb_source_relation=(
            "background.geometric_measure.loghb is the shared upstream width; "
            "R100 bridge hbB and switch hbS are distinct source objects"
        ),
        loghb_threshold=threshold_hb,
        hbS_source="R100_Rm_C6.native_parameter('hbS')",
        upstream_width_source="parameter_geometry_guard.shared.logh",
        source_nodes_distinct=True,
        source_node_aliasing_forbidden=True,
        hbB_hbS_receipt=dict(
            passed=True,
            hbB=hb_receipt["hbB"],
            hbS=hb_receipt["hbS"],
            handles_distinct=True,
            source_objects_distinct=True,
            shared_logh_object=True,
            switch_reboxing=certificate["switch_width_is_exact_source_reboxing"],
            shared_K1_ledger_sha256=ledger_hash,
        ),
        K1_ledger_receipt=ledger_receipt,
        R100_geometry_certificate=dict(certificate),
        R100_geometry_source_receipt=dict(
            path=R100_GEOMETRY_SOURCE_RECEIPT,
            expected_sha256=R100_GEOMETRY_SOURCE_RECEIPT_SHA256,
            shared_upstream_is_same_object=True,
            bridge_logh_is_actual_upstream_object=True,
            switch_width_is_exact_source_reboxing=True,
            source_executed=False,
        ),
        first_switch_cache_source_receipt=dict(
            path=FIRST_CACHE,
            expected_sha256=FIRST_CACHE_SHA256,
            phase=list(first_geometry["phase"]),
            endpoint_s=first_geometry["phase"][1],
            exact_source=first_geometry["exact_source"],
            source_frame=first_packet["source_frame"],
            true_sigma_function_used=True,
            source_executed=False,
        ),
        first_switch_endpoint_receipt=dict(
            chart="R100_Rm_C6.charts.first_switch",
            domain=list(first_chart["domain"]),
            endpoint_s=first_chart["domain"][1],
            physical_measure_graph_handle=first_chart["physical_measure"],
            source_executed=False,
        ),
        rho_domain=list(rho),
        alpha_Z_independent=True,
        sigma_range_receipt=dict(
            passed=True,
            original_sigma_not_its_sample_or_cap=True,
            source_owned=True,
            source_definition=sigma_receipt["definition"],
            source_binding=sigma_receipt["source_binding"],
            domain=[0, 1],
            range_enclosure="0<=sigma(s)<=1",
            endpoint_values={"s=0": 0, "s=1": 1},
            endpoint_flat=True,
            source_executed=False,
        ),
        sigma_source_receipt=sigma_receipt,
        first_switch_domain_receipt=dict(
            source_owned=True,
            domain=[0, 1],
            endpoint_s=1,
            endpoint_verified_from_R100_and_first_cache=True,
        ),
        R_upper=_point(c, 101),
        R_bound_formula="R=100*exp(hbS*s), 0<=s<=1",
        Wswitch_definition="10000*hbS",
        logWswitch=logW,
        physical_measure="R^2*hbS*ds",
        normalized_K=K,
        normalized_K_formula="integral_0^1 exp(2*hbS*s) ds",
        normalized_K_enclosure="[1,1+2^-1024]",
        normalized_K_proof=dict(
            log_width_plus_log4_lt_log_2minus1024=True,
            two_hbS_le_half_2minus1024=True,
            exp_half_width_lt_one_plus_2minus1024=True,
            K_is_enclosure_not_selected_value=True,
            hbS_not_materialized=True,
        ),
        radius_Z_independent=True,
        length_Z_independent=True,
        no_physical_hb_materialized=True,
        no_physical_R_materialized=True,
    )


def _fresh_exp_jet(c, radius):
    delta = [_sym(c, radius) for _ in ORDERS]
    result = [c.exp(delta[0])]
    for order in ORDERS[1:]:
        # Differentiate exp(delta): E_n = sum_{j=0}^{n-1}
        # C(n-1,j) E_j delta_(n-j).  This keeps ordinary derivatives
        # distinct from exponential Taylor coefficients.
        result.append(_sum(c, (c.mpf(math.comb(order - 1, j)) * result[j]
                               * delta[order - j] for j in range(order))))
    return delta, result


def _sector_and_drive(c, background, extension):
    bg_drive = background["full_signed_axial_drive"]
    proof = extension["comparison_fifth_radius_proof"]
    sector_source = {}
    for name, key in (("hydro", "hydro"), ("pressure", "pressure"),
                      ("dressed_swirl", "dressed_swirl")):
        source = _rows(c, bg_drive[key], 4, "accepted " + name + " sectors")
        fourth_key = {"hydro": "hydro4", "pressure": "pressure4",
                      "dressed_swirl": "dressedSwirl4"}[name]
        fourth = _read_record(c, proof[fourth_key])
        sector_source[name] = source + [fourth]
    pstar2 = _read_record(c, bg_drive["Pstar_squared"])
    full = []
    for order in ORDERS:
        value = _sum(c, (
            _scale(c, 101, _sym(c, sector_source["hydro"][order])),
            _scale(c, 101, pstar2 * _sym(c, sector_source["pressure"][order])),
            _scale(c, 10201, _sym(c, sector_source["dressed_swirl"][order])),
        ))
        full.append(_sym(c, value))
    P0source = _rows(c, bg_drive["P0"], 4, "accepted P0 rows")
    return dict(
        sectors={name: [_sym(c, value) for value in rows]
                 for name, rows in sector_source.items()},
        full_drive=full,
        Pstar_squared=pstar2,
        P0=P0source,
        source=dict(
            formula="fullDrive=101*hydro+101*Pstar^2*pressure+10201*dressed_swirl",
            R_upper=101,
            radius_scaling=dict(hydro=101, pressure=101, dressed_swirl=10201),
            F0_squared_cancellation_already_in_dressed_swirl=True,
            source_rows_0_to_3_are_symmetric_norm_caps=True,
            fourth_rows_from_accepted_macro_Z4_extension=True,
            signed_value_claim=False,
        ),
    )


def _comparison_rows(c, background, extension):
    macro = background["comparison"]
    proof = extension["comparison_fifth_radius_proof"]
    native_source_receipt = _verify_extension_native_rows(extension)
    native = extension["native_fields"]
    phi = _rows(c, proof["comparison_phi"], 5, "accepted frozen macro comparison Phi")
    U = _rows(c, macro["equilibrium_rows"]["M"], 5,
              "accepted frozen macro comparison U")
    history = {name: _rows(c, macro["history_rows"][name], 5,
                           "accepted frozen comparison " + name)
               for name in RATES}
    incoming = {name: _rows(c, native[name], 5,
                            "accepted macro-extension incoming " + name)
                for name in ("phi", "U", "M", "A", "B", "C")}
    directions = background["comparison_directions"]
    source = {name: _rows(c, directions[name], 4, "accepted direction " + name)
              for name in ("Z", "d", "L", "gaxis")}
    D = _rows(c, directions["D_over_R"], 4, "accepted D over R")
    D4 = _read_record(c, proof["D4cap"])
    D_with_four = D + [_sym(c, D4)]
    return dict(
        phi=phi,
        U=U,
        history=history,
        incoming=incoming,
        directions=source,
        D=D_with_four,
        D_norm=[_sym(c, value) for value in D_with_four],
        source=dict(
            frozen_macro_profile="comparison_at_radius: old comparison phi/U remain frozen",
            history_formula="decay*macro comparison endpoint+(1-decay)*same equilibrium",
            decay_Z_independent=True,
            rows_are_convex_hulls=True,
            D4_source="accepted macro native_Z4_extension.comparison_fifth_radius_proof.D4cap",
            input_boxes_unchanged=True,
            same_sector_formula_does_not_identify_fields=True,
            source_executed=False,
        ),
        native_source_receipt=native_source_receipt,
    )


def _first_switch_fields(c, extension, geometry, comparison, drive):
    incoming_phi = comparison["incoming"]["phi"]
    incoming_U = comparison["incoming"]["U"]
    incoming_histories = {name: comparison["incoming"][name]
                          for name in ("M", "A", "B", "C")}
    threshold = -c.mpf(1024) * c.log(c.mpf(2))
    loghb = _point(c, _hi(geometry["loghbS"]))
    Dlogs = []
    phi_log_rows = []
    for order in ORDERS:
        logD = _log_cap(c, comparison["D_norm"][order])
        bound = 2 * loghb + c.log(c.mpf(101)) + logD - c.log(c.mpf(2))
        _assert_lt(bound, threshold,
                   "First-switch angular row did not reach 2^-1024")
        Dlogs.append(bound)
        phi_log_rows.append(dict(
            ordinary_Z_order=order,
            log_absolute_upper_bound=bound,
            symmetric_radius=_point(c, c.power(c.mpf(2), -1024)),
            source_formula="-hbS^2/2*integral(R*D_over_R ds)",
        ))
    radius = c.power(c.mpf(2), -1024)
    delta_log, exp_delta = _fresh_exp_jet(c, _point(c, radius))
    phi = [_prod(c, incoming_phi, exp_delta, order) for order in ORDERS]

    ratio_inverse = _reciprocal(c, comparison["phi"], 5)
    ratio = [_prod(c, phi, ratio_inverse, order) for order in ORDERS]
    Ulog_rows = []
    U = []
    for order in ORDERS:
        product = _prod(c, ratio, drive["full_drive"], order)
        log_product = _log_cap(c, product)
        bound = 2 * loghb + log_product
        _assert_lt(bound, threshold,
                   "First-switch axial row did not reach 2^-1024")
        Ulog_rows.append(dict(
            ordinary_Z_order=order,
            log_absolute_upper_bound=bound,
            source_formula="-hbS*(1-sigma(s))*(phi/barphi)*fullDrive",
            sigma_range="0<=sigma<=1",
        ))
        U.append(incoming_U[order] + _sym(c, radius))

    densities = {
        "M": U,
        "A": [_prod(c, U, U, order) for order in ORDERS],
        "B": [_prod(c, phi, phi, order) for order in ORDERS],
        "C": [_prod(c, phi, phi, order) for order in ORDERS],
    }
    rates = {"M": 1, "A": 1, "B": 2, "C": 1}
    histories = {}
    history_proofs = {}
    for name in ("M", "A", "B", "C"):
        rows = []
        diff_logs = []
        rate = rates[name]
        for order in ORDERS:
            source_cap = _abs_upper(c, densities[name][order])
            incoming_cap = _abs_upper(c, incoming_histories[name][order])
            total = source_cap + c.mpf(rate) * incoming_cap
            _require(_lo(total) > 0, "Positive Duhamel source cap required")
            bound = loghb + c.log(total)
            _assert_lt(bound, threshold,
                       "First-switch history row did not reach 2^-1024")
            diff_logs.append(bound)
            rows.append(incoming_histories[name][order] + _sym(c, radius))
        histories[name] = rows
        history_proofs[name] = dict(
            rate=rate,
            density="U" if name == "M" else ("U^2" if name == "A" else "phi^2"),
            source_over_rate=("U" if name == "M" else
                              ("U^2" if name == "A" else
                               ("phi^2/2" if name == "B" else "phi^2"))),
            log_absolute_difference_rows=diff_logs,
            incoming_rows=incoming_histories[name],
            source_product_rows=densities[name],
            positive_kernel_bound="exp(-rate*hbS*(s-t))<=1",
            incoming_memory_retained=True,
            reset=False,
            source_executed=False,
        )
    return dict(
        phi=phi,
        U=U,
        histories=histories,
        densities=densities,
        ratio=ratio,
        delta_log=delta_log,
        exp_delta=exp_delta,
        phi_log_rows=phi_log_rows,
        U_log_rows=Ulog_rows,
        history_proofs=history_proofs,
        radius=radius,
        threshold=threshold,
    )


def _recover_Q(c, parent, fields, comparison):
    directions = comparison["directions"]
    Z, d, L = directions["Z"], directions["d"], directions["L"]
    delta = _read_record(c, parent["actual_core_signed_direction"]["parameters"]["delta"])
    h = _point(c, 1) - delta
    count = 3
    invL = _reciprocal(c, L, count)
    zU = [_prod(c, Z, fields["U"], n) for n in range(count)]
    zM = [_prod(c, Z, fields["histories"]["M"], n) for n in range(count)]
    mz = [fields["histories"]["M"][n + 1] for n in range(count)]
    hzm = [_scale(c, h, value) for value in zM]
    dmz = [_prod(c, d, mz, n) for n in range(count)]
    numerator = [_scale(c, 2, zU[n]) - hzm[n] - dmz[n] for n in range(count)]
    Q = [_prod(c, numerator, invL, n) for n in range(count)]
    return dict(
        Q=Q,
        Z=Z,
        d=d,
        L=L,
        inverse_L=invL,
        delta=delta,
        h=h,
        formula="Q_b=(2*Z*U-(1-delta)*Z*M-d*M_Z)/L",
        ordinary_product_binomial=True,
        source_executed=False,
    )


def _background_I1(c, parent, fields, qproof, geometry):
    """Compute the first-switch six integrals from the fresh Q/U rows."""
    U = fields["U"]
    Q = qproof["Q"]
    raw_integrands = dict(
        q=Q[:2],
        q_Z=Q[1:3],
        q_squared=[_prod(c, Q[:2], Q[:2], order) for order in C1_ORDERS],
        U_q_Z=[_prod(c, U[:2], Q[1:3], order) for order in C1_ORDERS],
        U_squared=[_prod(c, U[:2], U[:2], order) for order in C1_ORDERS],
        U_U_Z=[_prod(c, U[:2], U[1:3], order) for order in C1_ORDERS],
    )
    K = geometry["normalized_K"]
    integrals = {name: [K * value for value in rows]
                 for name, rows in raw_integrands.items()}
    Z = qproof["Z"][:2]
    d = qproof["d"][:2]
    L = qproof["L"][:2]
    delta = qproof["delta"]
    h = _point(c, 1) - delta
    invL = _reciprocal(c, L, 2)
    hZ = [_scale(c, h, Z[index]) for index in C1_ORDERS]
    invL2 = [_prod(c, invL, invL, order) for order in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, order) for order in C1_ORDERS]
    h_Z2 = [_scale(c, h, value) for value in Z2]
    Zd = [_prod(c, Z, d, order) for order in C1_ORDERS]
    coeff = {
        "minus_inv_L": [_neg(c, value) for value in invL],
        "hZ_over_2L": [_scale(c, c.mpf("0.5"), value)
                        for value in [_prod(c, hZ, invL, n) for n in C1_ORDERS]],
        "minus_half": [_point(c, "-0.5"), _point(c, 0)],
        "d_over_L": [_prod(c, d, invL, n) for n in C1_ORDERS],
        "two_h_Z2_over_L2": [_scale(c, 2, _prod(c, h_Z2, invL2, n))
                            for n in C1_ORDERS],
        "two_Zd_over_L2": [_scale(c, 2, _prod(c, Zd, invL2, n))
                           for n in C1_ORDERS],
    }
    D = []
    for order in C1_ORDERS:
        terms = (
            _prod(c, coeff["minus_inv_L"], integrals["q"], order),
            _prod(c, coeff["hZ_over_2L"], integrals["q_Z"], order),
            _prod(c, coeff["minus_half"], integrals["q_squared"], order),
            _prod(c, coeff["d_over_L"], integrals["U_q_Z"], order),
            _prod(c, coeff["two_h_Z2_over_L2"], integrals["U_squared"], order),
            _prod(c, coeff["two_Zd_over_L2"], integrals["U_U_Z"], order),
        )
        D.append(_sum(c, terms))
    return dict(
        normalized_by_Wmacro=False,
        normalized_by_Wswitch=True,
        Wmacro=None,
        positive_geometry_K=K,
        six_integrals={name: _pack_rows(c, rows) for name, rows in integrals.items()},
        D_rows=_pack_rows(c, D),
        coefficientsC1rows={name: _pack_rows(c, rows) for name, rows in coeff.items()},
        formula=("D=-I[q]/L+h*Z*I[q_Z]/(2L)-I[q_squared]/2+"
                 "d*I[U_q_Z]/L+2*h*Z^2*I[U_squared]/L^2+"
                 "2*Z*d*I[U_U_Z]/L^2"),
        coefficient_definitions=dict(h="1-delta", d="1-Z^2", L="1-delta*Z^2",
                                     hZ="(1-delta)*Z", h_Z2="h*Z^2",
                                     two_h_Z2_over_L2="2*h*Z^2/L^2", L2="L^2"),
        algebra_receipt=dict(
            passed=True,
            h_definition="h=1-delta",
            U_squared_coefficient="2*h*Z^2/L^2",
            source_match=True,
            avoids_h_squared_Z_squared=True,
        ),
        positive_geometry_multiplied_once=True,
        finite_N_correction_added=False,
        finite_N_correction_open=True,
        ordinary_Z_orders=list(C1_ORDERS),
        signed_C1_integrals_are_fresh_first_switch_enclosures=True,
        raw_first_switch_integrands={name: _pack_rows(c, rows)
                                     for name, rows in raw_integrands.items()},
        physical_weight_applied_once="R^2*hbS*ds / Wswitch = exp(2*hbS*s) ds",
        D_recomputed_with_fresh_ordinary_Leibniz=True,
        source_executed=False,
    )


def _source_pins(root):
    fixed = {
        PARENT: PARENT_SHA256,
        BACKGROUND: BACKGROUND_SHA256,
        R100_SOURCE: R100_SOURCE_SHA256,
        R100_CHECK: R100_CHECK_SHA256,
        R100_CACHE: R100_CACHE_SHA256,
        FIRST_SOURCE: FIRST_SOURCE_SHA256,
        FIRST_CHECK: FIRST_CHECK_SHA256,
        FIRST_CACHE: FIRST_CACHE_SHA256,
        C6_SOURCE: C6_SOURCE_SHA256,
        C6_CHECK: C6_CHECK_SHA256,
        C6_CACHE: C6_CACHE_SHA256,
        BRIDGE_SOURCE: BRIDGE_SOURCE_SHA256,
        NATIVE_EXT: NATIVE_EXT_SHA256,
        R100_GEOMETRY_SOURCE_RECEIPT: R100_GEOMETRY_SOURCE_RECEIPT_SHA256,
        K1_LEDGER: K1_LEDGER_SHA256,
        SIGMA_SOURCE: SIGMA_SOURCE_SHA256,
    }
    output = {path: _pin(root, path, digest) for path, digest in fixed.items()}
    output[MODULE] = _pin(root, MODULE, _sha(root / MODULE))
    return output


def build(root, parent):
    """Build the first-switch Z-only background packet from accepted data."""
    root = Path(root)
    disk_parent = _load_parent(root, parent)
    pins = _source_pins(root)
    guards = _source_guard(root)
    _require(isinstance(guards, tuple), "Source guards must remain a tuple")
    extension = disk_parent["actual_macro_loop_Z3_inputs"]["native_Z4_extension"]
    _require(extension.get("passed") is True, "Accepted macro Z4 extension is required")
    _require(extension.get("parentSHA") == BACKGROUND_SHA256,
             "Accepted macro Z4 extension parent hash mismatch")
    background = disk_parent["actual_macro_background_I1_numeric"]
    _require(background.get("passed") is True, "Accepted macro background packet is required")
    ea = disk_parent["actual_micro_Ea_mixed_bounds"]

    c = MPIntervalContext()
    c.prec = 1200
    drive_receipt = _verify_parent_drive(c, disk_parent)
    geometry = _build_geometry(root, disk_parent, c, background, ea)
    comparison = _comparison_rows(c, background, extension)
    drive = _sector_and_drive(c, background, extension)
    fields = _first_switch_fields(c, extension, geometry, comparison, drive)
    qproof = _recover_Q(c, disk_parent, fields, comparison)
    i1 = _background_I1(c, disk_parent, fields, qproof, geometry)

    native_fields = {
        "phi": _pack_rows(c, fields["phi"]),
        "U": _pack_rows(c, fields["U"]),
        "M": _pack_rows(c, fields["histories"]["M"]),
        "A": _pack_rows(c, fields["histories"]["A"]),
        "B": _pack_rows(c, fields["histories"]["B"]),
        "C": _pack_rows(c, fields["histories"]["C"]),
    }
    D_rows = _pack_rows(c, comparison["D"], kind="ordinary_D_enclosure")
    directions = {
        name: _pack_rows(c, values, kind="accepted_macro_direction")
        for name, values in comparison["directions"].items()
    }
    directions["D_over_R"] = D_rows
    directions["D_norm_caps"] = _pack_rows(c, comparison["D_norm"], kind="ordinary_norm_cap")
    sector_drive = {
        name: _pack_rows(c, rows, kind="ordinary_norm_cap")
        for name, rows in drive["sectors"].items()
    }
    sector_drive["full_drive"] = _pack_rows(c, drive["full_drive"], kind="ordinary_norm_cap")
    sector_drive["P0"] = _pack_rows(c, drive["P0"], kind="accepted_P0_enclosure")
    sector_drive["Pstar_squared"] = _record(drive["Pstar_squared"])

    comparison_packet = dict(
        frozen_phi=_pack_rows(c, comparison["phi"], kind="frozen_macro_comparison_phi"),
        frozen_U=_pack_rows(c, comparison["U"], kind="frozen_macro_comparison_U"),
        history_rows={name: _pack_rows(c, rows, kind="invariant_macro_history_hull")
                      for name, rows in comparison["history"].items()},
        decay_is_Z_independent=True,
        rows_are_convex_hulls=True,
        input_boxes_unchanged=True,
        source_binding=comparison["source"],
    )

    history_proofs = {
        name: _jsonify(proof) for name, proof in fields["history_proofs"].items()
    }
    row_proof = dict(
        phi=dict(
            formula="incoming_phi * exp(delta_logphi)",
            delta_log_rows=fields["phi_log_rows"],
            exp_increment_rows=_pack_rows(c, fields["exp_delta"], kind="fresh_exp_increment_jet"),
            incoming_macro_extension_rows=_pack_rows(
                c, comparison["incoming"]["phi"], kind="accepted_macro_endpoint_box"),
            value=native_fields["phi"],
            symmetric_error=_record(fields["radius"]),
            ordinary_product_binomial=True,
            source_executed=False,
        ),
        U=dict(
            formula="incoming_U +/- 2^-1024 after ratio*fullDrive transport test",
            log_rows=fields["U_log_rows"],
            incoming_macro_extension_rows=_pack_rows(
                c, comparison["incoming"]["U"], kind="accepted_macro_endpoint_box"),
            value=native_fields["U"],
            symmetric_error=_record(fields["radius"]),
            source_executed=False,
        ),
        M=dict(
            formula="incoming_M + Duhamel difference bound +/- 2^-1024",
            value=native_fields["M"],
            proof=history_proofs["M"],
            nonreset=True,
        ),
        A=history_proofs["A"],
        B=history_proofs["B"],
        C=history_proofs["C"],
    )

    ext_proof = extension["comparison_fifth_radius_proof"]
    f2_caps = extension["F2_derivative_log_caps0_4"]
    scalar_source = dict(
        D4cap=ext_proof["D4cap"],
        hydro4=ext_proof["hydro4"],
        pressure4=ext_proof["pressure4"],
        dressedSwirl4=ext_proof["dressedSwirl4"],
        fullDrive4=ext_proof["fullDrive4"],
        reused_accepted_macro_Z4_certificate=True,
        same_invariant_input_boxes=True,
        no_old_macro_arithmetic_executed=True,
    )
    background_I1 = i1.copy()
    background_I1["geometric_measure"] = dict(
        K=_record(geometry["normalized_K"]),
        loghbS=_record(geometry["loghbS"]),
        logWswitch=_record(geometry["logWswitch"]),
        R_upper=_record(geometry["R_upper"]),
        physical_measure=geometry["physical_measure"],
        normalized_by_Wswitch=True,
    )

    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parentSHA=PARENT_SHA256,
        ordinary_Z_orders=list(ORDERS),
        source_pins=pins,
        source_guards=guards,
        sourceguards_tuple=True,
        source_executed=False,
        physical_y_derivatives_claimed=False,
        global_five_moment_closure=False,
        finite_N_added=False,
        actual_phase_Z_independent=True,
        finite_N_phase_not_selected=True,
        geometry_measure=_jsonify(geometry),
        sourcegeometry=_jsonify(geometry),
        comparison=comparison_packet,
        comparison_invariance=dict(
            passed=True,
            frozen_phi_U=True,
            decay_Z_independent=True,
            convex_hull_rows=True,
            same_invariant_input_boxes=True,
            comparison_not_actual_field=True,
            history_rows_are_accepted_macro_hulls=True,
            source_guards_bind_C6_comparison_at_radius=True,
            native_extension_rows_exactly_verified=comparison["native_source_receipt"],
        ),
        comparison_directions=directions,
        full_signed_axial_drive=sector_drive,
        drive_source_proof=drive["source"],
        parent_axial_drive_source_receipt=drive_receipt,
        scalar_majorant_source=scalar_source,
        native_fields=native_fields,
        native_input_source=dict(
            path="actual_macro_loop_Z3_inputs.native_Z4_extension.native_fields",
            fields=("phi", "U", "M", "A", "B", "C"),
            ordinary_Z_orders=list(ORDERS),
            true_macro_endpoint_uniform_box=True,
            endpoint_rows_consumed_without_point_substitution=True,
            exact_tuple_rows_verified=comparison["native_source_receipt"],
        ),
        new_field_proofs=row_proof,
        Q_b=_pack_rows(c, qproof["Q"]),
        Q_b_proof=qproof,
        background_I1=background_I1,
        I1=background_I1,
        D=background_I1["D_rows"],
        F2_derivative_log_caps0_4=f2_caps,
        F2_log_caps0_4=f2_caps,
        source_contract=dict(
            passed=True,
            ordinary_Z_orders=list(ORDERS),
            first_switch_source="R100_Rm_C6.first_switch_fields",
            first_switch_radius="R=100*exp(hbS*s)",
            first_switch_measure="R^2*hbS*ds",
            normalized_by_Wswitch=True,
            comparison_frozen=True,
            incoming_memory_retained=True,
            signed_values_are_enclosures=True,
            physical_y_derivatives_claimed=False,
            global_five_moment_closure=False,
            finite_N_added=False,
            accepted_builders_checkers_producers_owners_norms_recurrence_quadrature_not_called=True,
            source_executed=False,
        ),
        limits=dict(
            physical_y_derivatives_claimed=False,
            global_five_moment_closure=False,
            finite_N_added=False,
            finite_N_phase_not_selected=True,
            actual_phase_Z_independent=True,
        ),
    )
    result = _jsonify(result)
    # Keep the runtime receipt a tuple while retaining ordinary JSON
    # serializability (json.dumps emits it as an array).
    result["source_guards"] = tuple(result["source_guards"])
    json.dumps(result)
    return result


def fixture_bound_check():
    """Fresh arithmetic-only checks for this module's Z-only propagation."""
    c = MPIntervalContext()
    c.prec = 180
    one = _point(c, 1)
    product = _prod(c, [one] * 5, [one] * 5, 4)
    _require(_hi(product) == 16, "Fresh ordinary binomial product fixture failed")
    radius = c.power(c.mpf(2), -32)
    delta, exponential = _fresh_exp_jet(c, _point(c, radius))
    _require(_lo(exponential[0]) > 0, "Fresh exponential increment fixture failed")
    _require(_hi(delta[4]) == radius, "Fresh symmetric increment fixture failed")
    K_upper = c.mpf(1) + radius
    K = c.mpf([_lo(c.mpf(1)), _hi(K_upper)])
    _require(_lo(K) == 1 and _hi(K) > 1, "Fresh normalized-K enclosure fixture failed")
    _require(isinstance(SOURCE_GUARDS, tuple), "Source guard recipe must be a tuple")
    return dict(passed=True, ordinary_binomial_order4=16,
                fresh_exp_increment=True, normalized_K_enclosure=True,
                sourceguards_tuple=True)


def self_check():
    return fixture_bound_check()
