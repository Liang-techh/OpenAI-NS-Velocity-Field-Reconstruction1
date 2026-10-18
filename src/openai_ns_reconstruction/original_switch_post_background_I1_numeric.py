"""Z-only second-switch and post-switch background arithmetic.

The packet is built from the accepted first-switch finite-N report and from
read-only source/cache receipts.  It does not import or execute any accepted
producer, owner, checker, norm, recurrence, quadrature, or exact-algebra
implementation.  All new field and integral rows are directed interval
enclosures; no physical y derivatives or phase values are selected.
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
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"
PARENT_SHA256 = "2cc8b4ff9670858564297bc1c27d5322c9af7d2cb85381939dc3e37a3eea346a"
FIRST_BG = "src/openai_ns_reconstruction/original_first_switch_background_I1_numeric.py"
FIRST_BG_SHA256 = "0c60e4bfdbd458c857363b119dcb61714a73d640230a3095c4c5f8583e74127d"
R100_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
R100_SOURCE_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
R100_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
R100_CACHE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
FIRST_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.json"
FIRST_CACHE_SHA256 = "5ea141f4ee1befaf33d6d2d2b9ad871477b72ee3dbbfc6b1e55b82766278a53d"
SECOND_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_second_switch_R110.py"
SECOND_SOURCE_SHA256 = "46df33f034c98ee885f6338141c9e8bbd7f29eb85999c229a17cda367e696023"
SECOND_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.py"
SECOND_CHECK_SHA256 = "a231bc6e70336cf850746a1472fdfb181b9c40ea9b97e1c1e53772e0d63e68b4"
SECOND_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_second_switch_R110.json"
SECOND_CACHE_SHA256 = "51b930434d4746c901eac1328eb0f1b8eae3ea2abb9fceca5cb8dc0b0ecc3751"
MODULE = "src/openai_ns_reconstruction/original_switch_post_background_I1_numeric.py"

ORDERS = (0, 1, 2, 3, 4)
C1_ORDERS = (0, 1)
Q_ORDERS = (0, 1, 2)
RATES = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}


# These are actual source AST bindings.  They are deliberately a tuple so
# callers can distinguish the guard recipe from any computed packet data.
SOURCE_GUARDS = (
    (R100_SOURCE, "second_switch_fields", (
        "a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))),",
        "g.mul(g.constant('-2/5'),hb,sig)",
        "original_second_switch_full_current_radius_direction_and_power=True",
    ), "af0fcc282e2569690f21f253e8b6a163390e63358db858c393fc0a447eb2115e"),
    (R100_SOURCE, "__init__", (
        "self.post_length=g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))",
        "R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))",
        "phi=a.scale(incoming['actual']['phi'],g.unary('exp',g.mul(c('-2/5'),y)))",
    ), "d875ec07da501f43644b9c89e2f0001976031e485c6637d3d018e6375cc796ef"),
    (SECOND_SOURCE, "source_bindings", (
        "original_second_V_zero_assignment=controls",
        "second_comparison_radius='100*exp(hb*(1+t)); fixed first-exit offset is retained'",
        "original_positive_post_power_assignments=transport",
    ), "d6a706023f3d13c1fef6421086bfbc3eddea82a060e46cd0e6e39d7c5c89d054"),
    (SECOND_SOURCE, "post", (
        "original_source_R2='100*exp(2hb)'",
        "original_post_power_a='4/5'",
        "actual_first_exit_V_preserved_exactly=True",
    ), "d788da7b1225cde10786af85af31a94394b2aa8b3b340840f99c57a04b690575"),
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
    low = c.mpf([make_mpf(_unsigned_tuple(value._mpi_[0]))] * 2)
    high = c.mpf([make_mpf(_unsigned_tuple(value._mpi_[1]))] * 2)
    bound = max(_hi(low), _hi(high))
    return c.mpf([bound, bound])


def _sym(c, value):
    upper = _abs_upper(c, value) if hasattr(value, "_mpi_") else _point(c, value)
    zero = c.mpf(0)
    neg = zero - upper
    return c.mpf([_lo(neg), _hi(upper)])


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


def _hull(c, left, right):
    return c.mpf([min(_lo(left), _lo(right)), max(_hi(left), _hi(right))])


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
        tail = _sum(c, (c.mpf(math.comb(order, j)) * values[j]
                        * result[order - j] for j in range(1, order + 1)))
        result.append(_neg(c, tail / values[0]))
    return result


def _log_cap(c, value):
    cap = _abs_upper(c, value)
    _require(_lo(cap) > 0, "Positive cap required before logarithm")
    return c.log(cap)


def _assert_lt(left, right, message):
    _require(_hi(left) < _lo(right), message)


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
        actual_ast = hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()
        _require(actual_ast == expected_ast, "Accepted AST changed: " + path + ":" + function)
        for snippet in snippets:
            _require("".join(snippet.split()) in compact,
                     "Accepted formula missing: " + path + ":" + function)
        output.append(dict(path=path, function=function, required_text=list(snippets),
                           AST_sha256=actual_ast, passed=True, source_executed=False))
    return tuple(output)


def _load_parent(root, parent):
    path = root / PARENT
    _require(path.is_file(), "Pinned first-switch finite-N parent is missing")
    _require(_sha(path) == PARENT_SHA256, "Pinned first-switch finite-N parent changed")
    disk = json.loads(path.read_bytes())
    _require(parent == disk, "The exact accepted first-switch finite-N parent mapping is required")
    return disk


def _load_json(root, path):
    actual_path = root / path
    _require(actual_path.is_file(), "Pinned JSON source is missing: " + path)
    value = json.loads(actual_path.read_bytes())
    _require(isinstance(value, Mapping), "Pinned JSON source must be a mapping: " + path)
    return value


def _load_json_gz(root, path):
    actual_path = root / path
    _require(actual_path.is_file(), "Pinned gzip source is missing: " + path)
    value = json.loads(gzip.open(actual_path, "rb").read())
    _require(isinstance(value, Mapping), "Pinned gzip source must be a mapping: " + path)
    return value


def _verify_first_background(root, parent):
    first = parent["actual_first_switch_background_I1_numeric"]
    _require(first.get("passed") is True, "Accepted first-switch background packet required")
    _require(first.get("parent_report_sha256") ==
             "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c",
             "First-switch background parent hash changed")
    pins = first.get("source_pins", {})
    _require(pins.get(FIRST_BG, {}).get("expected_sha256") == FIRST_BG_SHA256,
             "First-switch background implementation pin changed")
    _require(first.get("sourceguards_tuple") is True and
             first.get("source_executed") is False and
             first.get("actual_phase_Z_independent") is True,
             "First-switch background source contract changed")
    native = first.get("native_fields")
    _require(isinstance(native, Mapping), "First-switch native fields are required")
    for name in ("phi", "U", "M", "A", "B", "C"):
        _require(isinstance(native.get(name), list) and len(native[name]) == 5,
                 "First-switch native " + name + " rows are incomplete")
        _require([row.get("ordinary_Z_order") for row in native[name]] == list(ORDERS),
                 "First-switch native " + name + " rows are not ordered")
    geometry = first["geometry_measure"]
    for key in ("loghbS", "normalized_K", "hbB_hbS_receipt",
                "R100_geometry_source_receipt", "first_switch_cache_source_receipt",
                "sigma_source_receipt"):
        _require(key in geometry, "First-switch geometry receipt missing " + key)
    _require(geometry["hbB_hbS_receipt"].get("passed") is True and
             geometry["hbB_hbS_receipt"].get("handles_distinct") is True,
             "First-switch distinct hb source receipt changed")
    _require(geometry["R100_geometry_source_receipt"].get("expected_sha256") ==
             "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158",
             "First-switch R100 geometry source receipt changed")
    _require(geometry["first_switch_cache_source_receipt"].get("expected_sha256") ==
             FIRST_CACHE_SHA256,
             "First-switch cache source receipt changed")
    _require(geometry["sigma_source_receipt"].get("source_owned") is True and
             geometry["sigma_source_receipt"].get("domain") == [0, 1],
             "First-switch sigma source receipt changed")
    comparison = first.get("comparison")
    directions = first.get("comparison_directions")
    _require(isinstance(comparison, Mapping) and isinstance(directions, Mapping),
             "Frozen first-switch comparison packet required")
    _require(isinstance(directions.get("D_over_R"), list) and
             len(directions["D_over_R"]) == 5,
             "Frozen D over R rows 0..4 required")
    return first


def _verify_second_cache(root, parent):
    cache = _load_json(root, SECOND_CACHE)
    _require(parent.get("input_hashes", {}).get(
        "lei_ren_part1_paper_compliant_current_original_second_switch_R110.py") ==
             SECOND_SOURCE_SHA256,
             "Parent does not pin the accepted second-switch source")
    _require(parent.get("input_hashes", {}).get(
        "lei_ren_part1_paper_compliant_current_original_second_switch_R110.json") ==
             SECOND_CACHE_SHA256,
             "Parent does not pin the accepted second-switch cache")
    _require(cache.get("original_actual_second_switch_and_two_frame_R110_function_transport_installed") is True,
             "Accepted second-switch source cache is required")
    _require(cache.get("source_family"), "Second-switch source family is missing")
    bindings = cache.get("original_source_bindings")
    _require(isinstance(bindings, Mapping) and bindings.get("passed") is True,
             "Second-switch source bindings are required")
    _require(bindings.get("original_second_V_zero_assignment"),
             "Second-switch V=constant source binding is required")
    _require(bindings.get("original_positive_post_power_assignments"),
             "Post-power source binding is required")
    _require(bindings.get("second_comparison_radius") ==
             "100*exp(hb*(1+t)); fixed first-exit offset is retained",
             "Second-switch comparison radius source changed")
    second_packets = cache.get("second_switch_packets", {}).get("0")
    post_packets = cache.get("post_power_packets", {}).get("0")
    _require(isinstance(second_packets, list) and second_packets,
             "Second-switch source packet is missing")
    _require(isinstance(post_packets, list) and post_packets,
             "Post-switch source packet is missing")
    second = second_packets[0]
    second_eval = second.get("function_evaluation", {})
    second_geometry = second_eval.get("geometry", {})
    _require(second.get("source_frame") == "0" and
             second.get("source_family") == cache["source_family"] and
             second_geometry.get("phase") == [0, 1] and
             second_geometry.get("exact_source") == "R=100*exp(hb*(1+t)); dy=hb*dt",
             "Second-switch source geometry changed")
    _require(second_eval.get("ordinary_Z_orders") == list(range(6)) and
             second_eval.get("actual_first_exit_V_preserved_exactly") is True and
             second_eval.get("actual_moments_not_used_to_redefine_direction") is True,
             "Second-switch source packet row contract changed")
    post = post_packets[0]
    _require(post.get("source_frame") == "0" and
             post.get("source_family") == cache["source_family"] and
             post.get("actual_anchored_G_not_Gbar_used") is True and
             post.get("logCstar_cancelled_before_enclosure") is True and
             post.get("original_P0_retained_separately") is True,
             "Post-switch source ownership receipt changed")
    return dict(
        path=SECOND_CACHE,
        expected_sha256=SECOND_CACHE_SHA256,
        source_family=cache["source_family"],
        parent_input_hashes_verified=True,
        original_source_bindings=dict(bindings),
        second_switch=dict(
            source_frame=second["source_frame"],
            geometry=dict(second_geometry),
            ordinary_Z_orders=list(second_eval["ordinary_Z_orders"]),
            first_exit_preserved=True,
            direction_not_redefined=True,
        ),
        post_switch=dict(
            source_frame=post["source_frame"],
            anchored_G_not_Gbar=True,
            logCstar_cancelled_before_enclosure=True,
            original_P0_retained_separately=True,
        ),
        cache_values_not_consumed=True,
        source_executed=False,
    )


def _fresh_exp_jet(c, radius):
    delta = [_sym(c, radius) for unused in ORDERS]
    result = [c.exp(delta[0])]
    for order in ORDERS[1:]:
        result.append(_sum(c, (c.mpf(math.comb(order - 1, j)) * result[j]
                               * delta[order - j] for j in range(order))))
    return delta, result


def _directions(c, first):
    dirs = first["comparison_directions"]
    Z = _rows(c, dirs["Z"], 3, "inherited Z direction")
    d = _rows(c, dirs["d"], 3, "inherited d direction")
    Lz = _rows(c, dirs["L"], 3, "inherited L direction")
    D = _rows(c, dirs["D_over_R"], 5, "inherited D over R")
    delta = _read_record(c, first["Q_b_proof"]["delta"])
    h = _point(c, 1) - delta
    return dict(Z=Z, d=d, Lz=Lz, D=D, delta=delta, h=h)


def _second_fields(c, first, incoming, directions, geometry):
    """Fresh second-switch phi, U, and history boxes."""
    threshold = -c.mpf(1024) * c.log(c.mpf(2))
    loghb = geometry["loghbS"]
    R_upper = _point(c, 101)
    phi_in = incoming["phi"]
    log_rows = []
    for order in ORDERS:
        # sigma is bounded between zero and one by the source-owned receipt.
        logD = _log_cap(c, directions["D"][order])
        log_first = (c.mpf(2) * loghb + c.log(R_upper)
                     + logD - c.log(c.mpf(2)))
        log_second = loghb + c.log(c.mpf(2) / c.mpf(5))
        # x+y <= 2 max(x,y), so this stays entirely in the log domain.
        base = log_first if _hi(log_first) >= _hi(log_second) else log_second
        bound = base + c.log(c.mpf(2))
        _assert_lt(bound, threshold,
                   "Second-switch angular row did not reach 2^-1024")
        log_rows.append(dict(ordinary_Z_order=order,
                             log_absolute_upper_bound=bound,
                             source_formula=(
                                 "-hbS^2*R*D_over_R*(1-sigma)/2-2*hbS*sigma/5"),
                             sigma_range="0<=sigma<=1"))
    radius = c.power(c.mpf(2), -1024)
    delta_log, exp_delta = _fresh_exp_jet(c, _point(c, radius))
    phi = [_prod(c, phi_in, exp_delta, order) for order in ORDERS]
    U = list(incoming["U"])
    densities = {
        "M": U,
        "A": [_prod(c, U, U, order) for order in ORDERS],
        "B": [_prod(c, phi, phi, order) for order in ORDERS],
        "C": [_prod(c, phi, phi, order) for order in ORDERS],
    }
    histories = {}
    history_proofs = {}
    for name, rate in (("M", 1), ("A", 1), ("B", 2), ("C", 1)):
        rows = []
        logs = []
        for order in ORDERS:
            total = _abs_upper(c, densities[name][order]) + c.mpf(rate) * _abs_upper(c, incoming[name][order])
            bound = loghb + _log_cap(c, total)
            _assert_lt(bound, threshold,
                       "Second-switch history row did not reach 2^-1024")
            logs.append(bound)
            rows.append(incoming[name][order] + _sym(c, radius))
        histories[name] = rows
        history_proofs[name] = dict(
            rate=rate,
            density="U" if name == "M" else ("U^2" if name == "A" else "phi^2"),
            log_absolute_difference_rows=logs,
            incoming_rows=incoming[name],
            source_product_rows=densities[name],
            positive_kernel_bound="exp(-rate*hbS*t)<=1",
            incoming_memory_retained=True,
            reset=False,
            source_executed=False,
        )
    return dict(
        phi=phi,
        U=U,
        histories=histories,
        densities=densities,
        delta_log=delta_log,
        exp_delta=exp_delta,
        phi_log_rows=log_rows,
        history_proofs=history_proofs,
        radius=radius,
        threshold=threshold,
        source=dict(
            formula="R=100*exp(hbS*(1+t)); log(phi/phi_R1)=integral[-hbS^2*R*D*(1-sigma)/2-2*hbS*sigma/5]dt",
            U_preserved_exactly=True,
            sigma_source_owned=True,
            tiny_radius=_record(radius),
            source_executed=False,
        ),
    )


def _post_fields(c, incoming, geometry):
    """Post-switch positive transport with all incoming memory retained."""
    L0 = geometry["L0"]
    L = geometry["post_length"]
    phi_scale = c.mpf([c.exp(-c.mpf(2) * L0 / 5), c.mpf(1)])
    phi = [_scale(c, phi_scale, row) for row in incoming["phi"]]
    U = list(incoming["U"])
    # From the source log-width proof, hb<2^-1025.  This is an enclosure for
    # exp(-hb), not a selected or materialized physical width.
    exp_minus_hb = c.mpf([c.mpf(1) - geometry["r"], c.mpf(1)])
    exp_minus_hb_jet = [exp_minus_hb] + [_point(c, 0) for unused in range(4)]
    one_minus_exp_minus_hb_jet = [
        _point(c, 1) - exp_minus_hb
    ] + [_point(c, 0) for unused in range(4)]
    first_M = incoming["M_source"]
    # The same-source analytic endpoint is M_R2=U+(M_R1-U)*exp(-hbS).
    # Keeping the difference as an ordinary jet avoids replacing M_R1 by U.
    first_M_minus_U = [first_M[index] - U[index] for index in ORDERS]
    analytic_M = [U[order] + _prod(c, first_M_minus_U,
                                    exp_minus_hb_jet, order)
                  for order in ORDERS]
    # This is the independently retained-memory spelling of the same endpoint;
    # it remains a separate enclosure branch for the source receipt.
    memory_M = [_prod(c, first_M, exp_minus_hb_jet, order)
                + _prod(c, incoming["U"], one_minus_exp_minus_hb_jet, order)
                for order in ORDERS]
    incoming_M = [_hull(c, analytic_M[order], memory_M[order]) for order in ORDERS]
    e1 = c.exp(-L)
    e2 = c.exp(-c.mpf(2) * L)
    one = _point(c, 1)
    one_e1 = one - e1
    U2 = [_prod(c, U, U, order) for order in ORDERS]
    # The positive source over 0<=y<=L is phi_R2^2*exp(-4y/5).
    # `phi` above is the endpoint phi_R110 and must not be scaled twice.
    phi2 = [_prod(c, incoming["phi"], incoming["phi"], order)
            for order in ORDERS]
    phi2_scale = c.mpf([c.exp(-c.mpf(4) * L0 / 5), c.mpf(1)])
    # Exact positive convolution factors, with the source rate retained:
    # B'=phi^2-2B gives (5/6)(exp(-4L/5)-exp(-2L));
    # C'=phi^2-C gives 5(exp(-4L/5)-exp(-L)).
    source_decay = c.exp(-c.mpf(4) * L / 5)
    B_source_factor = c.mpf(5) / c.mpf(6) * (source_decay - e2)
    C_source_factor = c.mpf(5) * (source_decay - e1)
    _require(_lo(B_source_factor) > 0 and _lo(C_source_factor) > 0,
             "Post positive source convolution factors must be positive")
    phi2_uniform = [_scale(c, phi2_scale, row) for row in phi2]
    half_phi2_uniform = [_scale(c, c.mpf("0.5"), row)
                         for row in phi2_uniform]
    # These endpoint rows are retained for source review.  The native
    # histories below use full-segment convex hulls, since endpoint transport
    # alone is not a uniform enclosure for y in [0,L].
    endpoint_histories = {
        "M": list(incoming_M),
        "A": [_scale(c, e1, incoming["A"][order])
              + _scale(c, one_e1, U2[order]) for order in ORDERS],
        "B": [_scale(c, e2, incoming["B"][order])
              + _scale(c, B_source_factor, phi2[order])
              for order in ORDERS],
        "C": [_scale(c, e1, incoming["C"][order])
              + _scale(c, C_source_factor, phi2[order])
              for order in ORDERS],
    }
    # Positive-kernel transport keeps each ordinary row inside the hull of
    # its inlet and its nonnegative equilibrium/source density.  B has rate 2,
    # hence the equilibrium density is phi^2/2; C has rate 1.
    histories = {
        "M": [_hull(c, incoming_M[order], U[order]) for order in ORDERS],
        "A": [_hull(c, incoming["A"][order], U2[order])
              for order in ORDERS],
        "B": [_hull(c, incoming["B"][order], half_phi2_uniform[order])
              for order in ORDERS],
        "C": [_hull(c, incoming["C"][order], phi2_uniform[order])
              for order in ORDERS],
    }
    return dict(
        phi=phi,
        U=U,
        histories=histories,
        phi_scale=phi_scale,
        phi2_scale=phi2_scale,
        incoming_M=incoming_M,
        analytic_M=analytic_M,
        memory_M=memory_M,
        endpoint_histories=endpoint_histories,
        phi2_uniform=phi2_uniform,
        half_phi2_uniform=half_phi2_uniform,
        source_kernel_factors=dict(
            B=_record(B_source_factor),
            C=_record(C_source_factor),
        ),
        source=dict(
            formula="phi=phi_R2*exp(-2*y/5), U=U_R2; M=U+(M_R1-U)*exp(-hbS); M/A/B/C retain positive source transport",
            post_length_source="L=log(11/10)-2*hbS",
            phi_scale_enclosure="[exp(-2*log(11/10)/5),1]",
            phi2_scale_enclosure="[exp(-4*log(11/10)/5),1]",
            analytic_M_formula="M_R2=U+(M_R1-U)*exp(-hbS)",
            memory_M_formula="M_R1*exp(-hbS)+U*(1-exp(-hbS))",
            analytic_M_not_U_times_exp=True,
            incoming_M_same_source_endpoint=True,
            incoming_M_hull_includes_second_endpoint=True,
            B_source_convolution="phi_R2^2*(5/6)*(exp(-4*L/5)-exp(-2*L))",
            C_source_convolution="phi_R2^2*5*(exp(-4*L/5)-exp(-L))",
            B_rate=2,
            C_rate=1,
            positive_source_kernel_factors=True,
            native_histories_are_full_segment_hulls=True,
            native_M_hull="hull(incoming_M,U)",
            native_A_hull="hull(incoming_A,U^2)",
            native_B_hull="hull(incoming_B,phi^2_uniform/2)",
            native_C_hull="hull(incoming_C,phi^2_uniform)",
            endpoint_rows_retained_separately=True,
            scalar_width_jet_is_constant=True,
            reset=False,
            source_executed=False,
        ),
    )


def _background_B_rows(c, alpha, beta, U, directions, J, D, normalization):
    """Return the analytic normalized background ``B`` and ``S1=B+D`` rows.

    The constant-U alpha/beta representation is
    ``Q=alpha+beta*exp(-y)`` and ``Lz=1-delta*Z^2``.  Expanding the
    source algebra before interval boxing gives
    ``G=C0+C1*exp(-y)+C2*exp(-2*y)``.  The boundary calculation cancels the
    ``C2`` contribution symbolically, so the normalized row is exactly the
    requested ``2*C0*J0+C1*J1``.  This helper only uses the already enclosed
    normalized moments; it neither subtracts endpoints nor boxes an inverse
    tiny-width expression.
    """
    _require(normalization in ("Wswitch", "Wmacro"),
             "Background B normalization must be Wswitch or Wmacro")
    _require(len(alpha) >= 2 and len(beta) >= 2 and len(U) >= 2,
             "Background B requires alpha, beta, and U rows 0 and 1")
    _require(len(J) >= 2 and len(D) >= 2,
             "Background B requires J0, J1, and D rows 0 and 1")
    Z = directions["Z"]
    Lz = directions["Lz"]
    invL = _reciprocal(c, Lz, 2)
    invL2 = [_prod(c, invL, invL, order) for order in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, order) for order in C1_ORDERS]
    U2 = [_prod(c, U, U, order) for order in C1_ORDERS]
    alpha_over_L = [_prod(c, alpha, invL, order) for order in C1_ORDERS]
    alpha_squared_half = [
        _scale(c, c.mpf("0.5"), _prod(c, alpha, alpha, order))
        for order in C1_ORDERS
    ]
    Z2U2 = [_prod(c, Z2, U2, order) for order in C1_ORDERS]
    z2_u2_over_l2 = [
        _prod(c, Z2U2, invL2, order) for order in C1_ORDERS
    ]
    C0 = [
        _sum(c, (
            alpha_over_L[order],
            alpha_squared_half[order],
            _scale(c, -2, z2_u2_over_l2[order]),
        ))
        for order in C1_ORDERS
    ]
    beta_over_L = [_prod(c, beta, invL, order) for order in C1_ORDERS]
    alpha_beta = [_prod(c, alpha, beta, order) for order in C1_ORDERS]
    C1 = [_sum(c, (beta_over_L[order], alpha_beta[order]))
          for order in C1_ORDERS]
    C2 = [
        _scale(c, c.mpf("0.5"), _prod(c, beta, beta, order))
        for order in C1_ORDERS
    ]
    B = [
        _sum(c, (
            _scale(c, 2, _scale(c, J[0], C0[order])),
            _scale(c, J[1], C1[order]),
        ))
        for order in C1_ORDERS
    ]
    S1 = [_sum(c, (B[order], D[order])) for order in C1_ORDERS]
    W_formula = "10000*hbS" if normalization == "Wswitch" else "10000"
    receipt = dict(
        formula=(
            "Q=alpha+beta*exp(-y); Lz=1-delta*Z^2; "
            "G=Q/Lz+Q^2/2-2*Z^2*U^2/Lz^2; "
            "C0=alpha/Lz+alpha^2/2-2*Z^2*U^2/Lz^2; "
            "C1=beta/Lz+alpha*beta; C2=beta^2/2"
        ),
        boundary_formula="B=[R^2*G]_lo^hi",
        normalized_formula="B/W=2*C0*J0+C1*J1",
        S1_formula="S1=B+D",
        C0_rows=_pack_rows(c, C0, kind="background_C0_C1_jet"),
        C1_rows=_pack_rows(c, C1, kind="background_C1_C1_jet"),
        C2_rows=_pack_rows(c, C2, kind="background_C2_cancelled_jet"),
        C2_symbolically_cancelled=True,
        C2_cancelled_before_interval_boxing=True,
        C2_not_used_in_B_rows=True,
        no_endpoint_subtraction=True,
        no_inverse_tiny_hb_boxing=True,
        normalization=normalization,
        W_formula=W_formula,
        normalized_moments={
            "J0": _record(J[0]),
            "J1": _record(J[1]),
        },
        ordinary_Z_orders=list(C1_ORDERS),
        source_executed=False,
    )
    return dict(B=B, S1=S1, C0=C0, C1=C1, C2=C2, receipt=receipt)


def _q_and_i1(c, U, M, directions, J, decay_interval, normalization):
    """Compute pointwise Q and the six ordinary-C1 integrals from alpha+beta exp(-y).

    ``Q`` is a pointwise enclosure obtained from a Z-independent interval for
    ``exp(-y)``.  The ``q`` rows below are the integrated ``I[q]`` rows and
    retain the normalized moments ``J``.  These two objects must not be
    conflated: post-switch ``J0`` is about 0.105, while pointwise decay is
    bounded below by 10/11.
    """
    Z, d, Lz, h = directions["Z"], directions["d"], directions["Lz"], directions["h"]
    invL = _reciprocal(c, Lz, 3)
    delta = directions["delta"]
    one_plus_delta = _point(c, 1) + delta
    U3 = list(U[:3])
    U_Z = [U[1], U[2], U[3]]
    M3 = list(M[:3])
    M_minus_U = [M3[i] - U3[i] for i in range(3)]
    M_Z_minus_U_Z = [M[1] - U[1], M[2] - U[2], M[3] - U[3]]
    one_delta_Z = [one_plus_delta * Z[n] for n in Q_ORDERS]
    alpha_num = [_prod(c, one_delta_Z, U3, n)
                 - _prod(c, d, U_Z, n) for n in Q_ORDERS]
    alpha = [_prod(c, alpha_num, invL, n) for n in Q_ORDERS]
    hz = [h * Z[n] for n in Q_ORDERS]
    beta_num = [_neg(c, _prod(c, hz, M_minus_U, n))
                - _prod(c, d, M_Z_minus_U_Z, n) for n in Q_ORDERS]
    beta = [_prod(c, beta_num, invL, n) for n in Q_ORDERS]
    J0, J1, J2 = J
    q = [_scale(c, J0, alpha[n]) + _scale(c, J1, beta[n]) for n in Q_ORDERS]
    Q_uniform = [alpha[n] + beta[n] * decay_interval for n in Q_ORDERS]
    q_Z = [q[1], q[2]]
    q_squared = [(_scale(c, J0, _prod(c, alpha, alpha, n))
                  + _scale(c, c.mpf(2) * J1, _prod(c, alpha, beta, n))
                  + _scale(c, J2, _prod(c, beta, beta, n))) for n in C1_ORDERS]
    Uq_Z = [_prod(c, U[:2], q_Z, n) for n in C1_ORDERS]
    U_squared = [_scale(c, J0, _prod(c, U[:2], U[:2], n)) for n in C1_ORDERS]
    U_U_Z = [_scale(c, J0, _prod(c, U[:2], U_Z[:2], n)) for n in C1_ORDERS]
    invL2 = [_prod(c, invL, invL, n) for n in C1_ORDERS]
    hZ = [h * Z[n] for n in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, n) for n in C1_ORDERS]
    h_Z2 = [_scale(c, h, value) for value in Z2]
    Zd = [_prod(c, Z, d, n) for n in C1_ORDERS]
    coeff = {
        "minus_inv_L": [_neg(c, value) for value in invL[:2]],
        "hZ_over_2L": [_scale(c, c.mpf("0.5"), _prod(c, hZ, invL, n))
                        for n in C1_ORDERS],
        "minus_half": [_point(c, "-0.5"), _point(c, 0)],
        "d_over_L": [_prod(c, d, invL, n) for n in C1_ORDERS],
        "two_h_Z2_over_L2": [_scale(c, 2, _prod(c, h_Z2, invL2, n))
                              for n in C1_ORDERS],
        "two_Zd_over_L2": [_scale(c, 2, _prod(c, Zd, invL2, n))
                            for n in C1_ORDERS],
    }
    integrands = {
        "q": [_scale(c, J0, alpha[n]) + _scale(c, J1, beta[n]) for n in C1_ORDERS],
        "q_Z": q_Z,
        "q_squared": q_squared,
        "U_q_Z": Uq_Z,
        "U_squared": U_squared,
        "U_U_Z": U_U_Z,
    }
    D = []
    for n in C1_ORDERS:
        D.append(_sum(c, (
            _prod(c, coeff["minus_inv_L"], integrands["q"], n),
            _prod(c, coeff["hZ_over_2L"], integrands["q_Z"], n),
            _prod(c, coeff["minus_half"], integrands["q_squared"], n),
            _prod(c, coeff["d_over_L"], integrands["U_q_Z"], n),
            _prod(c, coeff["two_h_Z2_over_L2"], integrands["U_squared"], n),
            _prod(c, coeff["two_Zd_over_L2"], integrands["U_U_Z"], n),
        )))
    background = _background_B_rows(
        c, alpha, beta, U, directions, J, D, normalization,
    )
    return dict(
        Q=Q_uniform,
        Q_uniform=Q_uniform,
        integrated_q=q,
        decay_interval=decay_interval,
        alpha=alpha,
        beta=beta,
        integrals=integrands,
        D=D,
        coefficients=coeff,
        background_B=background["B"],
        background_S1=background["S1"],
        background_pressure_algebra=background["receipt"],
        Q_formula=("Q_pointwise=alpha+beta*decay_interval; "
                   "decay_interval is a Z-independent enclosure; "
                   "alpha=((1+delta)*Z*U-d*U_Z)/Lz; "
                   "beta=(-h*Z*(M_in-U)-d*(M_in_Z-U_Z))/Lz"),
        integrated_q_formula="I[q]: q=alpha*J0+beta*J1",
        formula=("Q_pointwise=alpha+beta*decay_interval; "
                 "decay_interval is a Z-independent enclosure; "
                 "alpha=((1+delta)*Z*U-d*U_Z)/Lz; "
                 "beta=(-h*Z*(M_in-U)-d*(M_in_Z-U_Z))/Lz"),
        D_formula=("D=-I[q]/Lz+h*Z*I[q_Z]/(2Lz)-I[q_squared]/2+"
                   "d*I[U_q_Z]/Lz+2*h*Z^2*I[U_squared]/Lz^2+"
                   "2*Z*d*I[U_U_Z]/Lz^2"),
        algebra_receipt=dict(
            h="1-delta",
            U_squared_coefficient="2*h*Z^2/Lz^2",
            canonical_key="two_h_Z2_over_L2",
            avoids_h_squared_Z_squared=True,
            coefficient_product_binomial=True,
        ),
        integral_vs_pointwise_receipt=dict(
            Q_rows_are_pointwise_decay_enclosures=True,
            integrated_q_rows_are_Iq_rows=True,
            J_not_applied_to_Q=True,
            decay_interval_Z_independent=True,
            Q_and_integrated_Iq_are_distinct=True,
        ),
        source_executed=False,
    )


def _weight_derivative_fixture(c, r, L0):
    """Execute the small derivative bounds used by the normalized weights.

    The source formulas are checked at the enclosing cap ``hbS <= r/100``;
    this is a fresh arithmetic certificate, not an execution of the accepted
    source implementation.
    """
    h_cap = r / c.mpf(100)
    e2 = c.exp(c.mpf(2) * h_cap)
    e4 = c.exp(c.mpf(4) * h_cap)
    second_cap = c.mpf(4) * e4
    _assert_lt(second_cap, c.mpf(100),
               "Second normalized-weight derivative fixture failed")
    post_d0 = c.mpf(2) * e4
    post_d1 = c.mpf("2.2") * e2 - c.mpf(4) * e4
    post_d2 = e4 * (c.mpf(4) * L0 - c.mpf(8) * h_cap - c.mpf(2))
    for value in (post_d0, post_d1, post_d2):
        _assert_lt(_abs_upper(c, value), c.mpf(100),
                   "Post normalized-weight derivative fixture failed")
    return dict(
        passed=True,
        h_cap=_record(h_cap),
        second_derivative_cap=_record(second_cap),
        post_derivatives=[_record(post_d0), _record(post_d1), _record(post_d2)],
        source_executed=False,
    )


def _geometry(c, first):
    source = first["geometry_measure"]
    loghb = _read_record(c, source["loghbS"])
    _require(source.get("loghb_shared_object") is True and
             source.get("loghb_upstream_not_independent") is True,
             "Shared first-switch log-width receipt required")
    _require(source.get("normalized_K_enclosure") == "[1,1+2^-1024]",
             "First-switch normalized-K source receipt changed")
    threshold = -c.mpf(1024) * c.log(c.mpf(2))
    r = c.power(c.mpf(2), -1024)
    _require(_lo(r) > 0, "Positive tiny radius required")
    log_r = c.log(r)
    _assert_lt(loghb + c.log(c.mpf(100)), log_r,
               "Source log-width does not close the tiny post/second radius")
    _assert_lt(loghb + c.log(c.mpf(2)), log_r,
               "Source log-width does not close the post-length enclosure")
    L0 = c.log(c.mpf(11) / c.mpf(10))
    weight_fixture = _weight_derivative_fixture(c, r, L0)
    # The actual hb is intentionally kept logarithmic.  The following boxes
    # use the proved hb<r/2 relation and never materialize exp(loghb).
    post_length = c.mpf([L0 - r, L0])
    _require(_lo(post_length) > 0, "Post-switch length must be positive")
    second_J = [c.mpf([c.mpf(1), c.mpf(1) + r]) for unused in range(3)]
    post_J = [
        c.mpf([c.mpf(21) / 200 - r, c.mpf(21) / 200]),
        c.mpf([c.mpf(1) / 10 - r, c.mpf(1) / 10 + r]),
        c.mpf([L0 - r, L0 + r]),
    ]
    _require(_lo(post_J[0]) > 0 and _lo(post_J[1]) > 0 and _lo(post_J[2]) > 0,
             "Post normalized moment lower endpoints must be positive")
    return dict(
        loghbS=loghb,
        r=r,
        threshold=threshold,
        L0=L0,
        post_length=post_length,
        second_J=second_J,
        post_J=post_J,
        weight_derivative_fixture=weight_fixture,
        executable_weight_derivative_fixture=True,
        second_radius_upper=_point(c, 101),
        post_radius_upper=_point(c, 110),
        source=dict(
            first_geometry_source="accepted first-switch geometry_measure",
            shared_loghb_receipt=True,
            hbS_kept_logarithmic=True,
            hbS_not_selected_or_materialized=True,
            second_formula="R=100*exp(hbS*(1+t)), t in [0,1]",
            post_formula="y=post_length*p, R=100*exp(2*hbS+y), post_length=log(11/10)-2*hbS",
            second_normalized_weight="Jk=exp(2*hbS)*integral_0^1 exp((2-k)*hbS*t)dt",
            post_normalized_weight="Jk=exp(4*hbS)*integral_0^post_length exp((2-k)*y)dy",
            second_J_enclosure="[1,1+2^-1024] for k=0,1,2",
            second_J_derivative_bound="|dJk/dhbS|<100 and Jk(hbS=0)=1",
            post_J_exact_formulas=dict(
                J0="(1.21-exp(4*hbS))/2",
                J1="1.1*exp(2*hbS)-exp(4*hbS)",
                J2="exp(4*hbS)*(log(11/10)-2*hbS)",
            ),
            post_J_derivative_bound="all three derivative absolute values <100",
            post_J_enclosure=(
                "J0=[21/200-r,21/200], J1=[1/10-r,1/10+r], "
                "J2=[log(11/10)-r,log(11/10)+r]"
            ),
            actual_weight_interval_not_selected=True,
            no_exp_of_huge_negative_loghb=True,
            source_executed=False,
        ),
    )


def _chart_geometry(c, geometry, J, normalization, R_upper, measure, J_name):
    return dict(
        normalization=normalization,
        normalized_moments=[_record(value) for value in J],
        mass=_record(J[0]),
        loghbS=_record(geometry["loghbS"]),
        R_upper=_record(R_upper),
        physical_measure=measure,
        physical_J_name=J_name,
        actual_weight_interval_not_selected=True,
        source_executed=False,
    )


def _incoming_fields(c, first):
    native = first["native_fields"]
    return {name: _rows(c, native[name], 5, "first-stage incoming " + name)
            for name in ("phi", "U", "M", "A", "B", "C")}


def _background_packet(c, qproof, geometry, normalization, J, source_kind):
    return dict(
        six_integrals={name: _pack_rows(c, rows)
                       for name, rows in qproof["integrals"].items()},
        D_rows=_pack_rows(c, qproof["D"], kind="ordinary_D_enclosure"),
        background_B_rows=_pack_rows(c, qproof["background_B"],
                                     kind="ordinary_background_B_over_W"),
        background_S1_rows=_pack_rows(c, qproof["background_S1"],
                                      kind="ordinary_background_S1_over_W"),
        background_pressure_algebra=qproof["background_pressure_algebra"],
        coefficientsC1rows={name: _pack_rows(c, rows)
                            for name, rows in qproof["coefficients"].items()},
        normalized_by_Wswitch=normalization == "Wswitch",
        normalized_by_Wmacro=normalization == "Wmacro",
        normalization=normalization,
        normalized_moments=[_record(value) for value in J],
        Q_rows_are_pointwise_decay_enclosures=True,
        Q_uniform=_pack_rows(c, qproof["Q_uniform"]),
        integrated_q=_pack_rows(c, qproof["integrated_q"]),
        decay_interval=_record(qproof["decay_interval"]),
        Q_formula=qproof["Q_formula"],
        integrated_q_formula=qproof["integrated_q_formula"],
        pointwise_Q_distinct_from_integrated_Iq=True,
        formula=qproof["formula"],
        D_formula=qproof["D_formula"],
        algebra_receipt=qproof["algebra_receipt"],
        source_kind=source_kind,
        source_executed=False,
    )


def _source_pins(root):
    fixed = {
        PARENT: PARENT_SHA256,
        FIRST_BG: FIRST_BG_SHA256,
        R100_SOURCE: R100_SOURCE_SHA256,
        R100_CACHE: R100_CACHE_SHA256,
        FIRST_CACHE: FIRST_CACHE_SHA256,
        SECOND_SOURCE: SECOND_SOURCE_SHA256,
        SECOND_CHECK: SECOND_CHECK_SHA256,
        SECOND_CACHE: SECOND_CACHE_SHA256,
    }
    output = {path: _pin(root, path, digest) for path, digest in fixed.items()}
    output[MODULE] = _pin(root, MODULE, _sha(root / MODULE))
    return output


def build(root, parent):
    """Build second-switch and post-switch Z-only background packets."""
    root = Path(root)
    disk_parent = _load_parent(root, parent)
    first = _verify_first_background(root, disk_parent)
    second_source_receipt = _verify_second_cache(root, disk_parent)
    pins = _source_pins(root)
    guards = _source_guard(root)
    _require(isinstance(guards, tuple), "Source guards must remain a tuple")
    c = MPIntervalContext()
    c.prec = 1200
    geometry = _geometry(c, first)
    dirs = _directions(c, first)
    incoming = _incoming_fields(c, first)
    second_fields = _second_fields(c, first, incoming, dirs, geometry)
    second_incoming = dict(incoming)
    second_incoming.update({name: second_fields["histories"][name]
                            for name in ("M", "A", "B", "C")})
    second_decay = c.mpf([c.mpf(1) - geometry["r"], c.mpf(1)])
    post_decay = c.mpf([c.mpf(10) / c.mpf(11), c.mpf(1)])
    second_q = _q_and_i1(c, second_fields["U"], incoming["M"], dirs,
                         geometry["second_J"], second_decay, "Wswitch")
    post_fields_input = dict(second_incoming)
    post_fields_input["phi"] = second_fields["phi"]
    post_fields_input["U"] = second_fields["U"]
    post_fields_input["M_source"] = incoming["M"]
    post_fields = _post_fields(c, post_fields_input, geometry)
    post_q = _q_and_i1(c, post_fields["U"], post_fields["incoming_M"], dirs,
                       geometry["post_J"], post_decay, "Wmacro")

    second_native = {
        "phi": second_fields["phi"],
        "U": second_fields["U"],
        "M": second_fields["histories"]["M"],
        "A": second_fields["histories"]["A"],
        "B": second_fields["histories"]["B"],
        "C": second_fields["histories"]["C"],
    }
    post_native = {
        "phi": post_fields["phi"],
        "U": post_fields["U"],
        "M": post_fields["histories"]["M"],
        "A": post_fields["histories"]["A"],
        "B": post_fields["histories"]["B"],
        "C": post_fields["histories"]["C"],
    }
    incoming_packet = {name: _pack_rows(c, rows, kind="first_stage_uniform_box")
                       for name, rows in incoming.items()}
    second_incoming_packet = {name: _pack_rows(c, rows, kind="second_stage_incoming_box")
                              for name, rows in second_incoming.items()
                              if name in ("phi", "U", "M", "A", "B", "C")}
    post_incoming_packet = {name: _pack_rows(c, rows, kind="post_stage_incoming_box")
                            for name, rows in post_fields_input.items()
                            if name in ("phi", "U", "M", "A", "B", "C")}
    second_chart = dict(
        passed=True,
        source_only_identity=False,
        new_evaluated_signed_enclosures=True,
        native_fields={name: _pack_rows(c, rows) for name, rows in second_native.items()},
        incoming_fields=incoming_packet,
        Q_b=_pack_rows(c, second_q["Q"]),
        background_I1=_background_packet(c, second_q, geometry, "Wswitch",
                                         geometry["second_J"], "second_switch"),
        background_B_rows=_pack_rows(
            c, second_q["background_B"], kind="ordinary_background_B_over_W"),
        background_S1_rows=_pack_rows(
            c, second_q["background_S1"], kind="ordinary_background_S1_over_W"),
        background_pressure_algebra=second_q["background_pressure_algebra"],
        geometry=_chart_geometry(c, geometry, geometry["second_J"], "Wswitch",
                                 geometry["second_radius_upper"],
                                 "R^2*hbS*dt / Wswitch", "hbS"),
        source_proof=dict(
            second_switch=second_fields["source"],
            history_proofs=second_fields["history_proofs"],
            Q_source=second_q["Q_formula"],
            integral_vs_pointwise=second_q["integral_vs_pointwise_receipt"],
            decay_interval=_record(second_q["decay_interval"]),
            integrated_q=_pack_rows(c, second_q["integrated_q"]),
            cache_source=second_source_receipt["second_switch"],
            native_rows_are_uniform_boxes=True,
            source_executed=False,
        ),
    )
    post_chart = dict(
        passed=True,
        source_only_identity=False,
        new_evaluated_signed_enclosures=True,
        native_fields={name: _pack_rows(c, rows) for name, rows in post_native.items()},
        incoming_fields=post_incoming_packet,
        Q_b=_pack_rows(c, post_q["Q"]),
        background_I1=_background_packet(c, post_q, geometry, "Wmacro",
                                         geometry["post_J"], "post_switch"),
        background_B_rows=_pack_rows(
            c, post_q["background_B"], kind="ordinary_background_B_over_W"),
        background_S1_rows=_pack_rows(
            c, post_q["background_S1"], kind="ordinary_background_S1_over_W"),
        background_pressure_algebra=post_q["background_pressure_algebra"],
        geometry=_chart_geometry(c, geometry, geometry["post_J"], "Wmacro",
                                 geometry["post_radius_upper"],
                                 "R^2*post_length*dp / Wmacro", "post_length"),
        source_proof=dict(
            post_switch=post_fields["source"],
            incoming_M_endpoint=dict(
                analytic_rows=_pack_rows(c, post_fields["analytic_M"]),
                same_source_memory_rows=_pack_rows(c, post_fields["memory_M"]),
                hulled_rows=_pack_rows(c, post_fields["incoming_M"]),
                includes_second_endpoint=True,
                source_executed=False,
            ),
            Q_source=post_q["Q_formula"],
            integral_vs_pointwise=post_q["integral_vs_pointwise_receipt"],
            decay_interval=_record(post_q["decay_interval"]),
            integrated_q=_pack_rows(c, post_q["integrated_q"]),
            cache_source=second_source_receipt["post_switch"],
            native_rows_are_positive_memory_hulls=True,
            no_reset=True,
            source_executed=False,
        ),
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
        actual_phase_Z_independent=True,
        finite_N_phase_not_selected=True,
        charts={"second_switch": second_chart, "post_switch": post_chart},
        native_fields={
            "second_switch": second_chart["native_fields"],
            "post_switch": post_chart["native_fields"],
        },
        incoming_fields=incoming_packet,
        Q_b={"second_switch": second_chart["Q_b"], "post_switch": post_chart["Q_b"]},
        comparison=first["comparison"],
        comparison_directions=first["comparison_directions"],
        F2_log_caps0_4=first.get("F2_log_caps0_4", first.get("F2_derivative_log_caps0_4")),
        geometry=_jsonify(geometry),
        sourceproof=dict(
            first_background_source_pins_verified=True,
            first_background_geometry_receipts_reused=True,
            second_switch_cache=second_source_receipt,
            second_switch_source_ast=guards,
            cache_values_not_consumed=True,
            source_only_identity_receipts=True,
            new_evaluated_signed_enclosures=True,
            source_executed=False,
        ),
        source_contract=dict(
            passed=True,
            ordinary_Z_orders=list(ORDERS),
            second_switch_source="R100_Rm_C6.second_switch_fields",
            post_switch_source="R100_Rm_C6.__init__ post_switch assignment",
            comparison_inherited_frozen=True,
            first_stage_uniform_boxes_retained=True,
            incoming_memory_retained=True,
            normalized_second="Wswitch",
            normalized_post="Wmacro",
            physical_y_derivatives_claimed=False,
            finite_N_phase_not_selected=True,
            source_executed=False,
            source_only_identity_receipts=True,
            new_evaluated_signed_enclosures=True,
        ),
        limits=dict(
            physical_y_derivatives_claimed=False,
            global_five_moment_closure=False,
            finite_N_phase_not_selected=True,
            source_executed=False,
        ),
    )
    result = _jsonify(result)
    result["source_guards"] = tuple(result["source_guards"])
    json.dumps(result)
    return result


def _post_transport_fixture(c):
    """Fresh scalar witness for post positive kernels and the M bridge."""
    ell = c.log(c.mpf(11) / c.mpf(10))
    e1 = c.exp(-ell)
    e2 = c.exp(-c.mpf(2) * ell)
    source_decay = c.exp(-c.mpf(4) * ell / 5)
    b_factor = c.mpf(5) / c.mpf(6) * (source_decay - e2)
    c_factor = c.mpf(5) * (source_decay - e1)
    expected_b = c.mpf([c.mpf("0.0834499753"), c.mpf("0.0834499754")])
    expected_c = c.mpf([c.mpf("0.0874767113"), c.mpf("0.0874767114")])
    b_from_zero = _point(c, 0) * e2 + b_factor * _point(c, 1)
    c_from_zero = _point(c, 0) * e1 + c_factor * _point(c, 1)
    _require(_lo(b_from_zero) >= _lo(expected_b) and
             _hi(b_from_zero) <= _hi(expected_b),
             "Fresh rate-2 B convolution fixture failed")
    _require(_lo(c_from_zero) >= _lo(expected_c) and
             _hi(c_from_zero) <= _hi(expected_c),
             "Fresh rate-1 C convolution fixture failed")
    wrong_b_factor = _point(c, 1) - e2
    _require(_hi(b_factor) < _lo(wrong_b_factor),
             "Fresh B fixture did not retain the rate-2 divisor")
    hb0 = _point(c, 0)
    exp_hb0 = c.exp(-hb0)
    M_R1 = _point(c, 4)
    U = _point(c, 2)
    corrected_M = U + (M_R1 - U) * exp_hb0
    memory_M = M_R1 * exp_hb0 + U * (_point(c, 1) - exp_hb0)
    wrong_M = U * exp_hb0
    _require(_lo(corrected_M) == 4 and _hi(corrected_M) == 4,
             "Fresh M bridge fixture did not retain M_R1")
    _require(_lo(memory_M) <= _lo(corrected_M) and
             _hi(memory_M) >= _hi(corrected_M),
             "Fresh M memory branch misses corrected endpoint")
    _require(_hi(wrong_M) < _lo(corrected_M),
             "Fresh M fixture did not reject U*exp(-hbS)")
    return dict(
        passed=True,
        rate2_B_factor=_record(b_factor),
        rate1_C_factor=_record(c_factor),
        B_expected_decimal_interval="[0.0834499753,0.0834499754]",
        C_expected_decimal_interval="[0.0874767113,0.0874767114]",
        B_rate2_divisor=True,
        M_corrected_formula="U+(M_R1-U)*exp(-hbS)",
        M_memory_branch_contains_corrected=True,
        M_wrong_U_times_exp_rejected=True,
        source_executed=False,
    )


def _background_B_fixture(c):
    """Fresh scalar witness for the C2 cancellation and S1 addition."""
    one = _point(c, 1)
    zero = _point(c, 0)
    two = _point(c, 2)
    directions = {
        "Z": [one, zero, zero],
        "Lz": [one, zero, zero],
    }
    alpha = [one, zero, zero]
    beta = [two, zero, zero]
    U = [one, zero, zero]
    J = [one, one, one]
    D = [zero, zero]
    result = _background_B_rows(
        c, alpha, beta, U, directions, J, D, "Wswitch")
    # C0=-1/2, C1=4, C2=2, so the prescribed 2*C0*J0+C1*J1 is 3;
    # including the cancelled C2*J2 term would incorrectly give 5.
    _require(_lo(result["B"][0]) == 3 and _hi(result["B"][0]) == 3,
             "Fresh background B C2-cancellation fixture failed")
    _require(_lo(result["S1"][0]) == 3 and _hi(result["S1"][0]) == 3,
             "Fresh background S1 fixture failed")
    receipt = result["receipt"]
    _require(receipt["C2_symbolically_cancelled"] is True and
             receipt["C2_cancelled_before_interval_boxing"] is True and
             receipt["no_endpoint_subtraction"] is True and
             receipt["no_inverse_tiny_hb_boxing"] is True,
             "Fresh background B algebra receipt failed")
    return dict(
        passed=True,
        B0=_record(result["B"][0]),
        S10=_record(result["S1"][0]),
        C2_not_added=True,
        source_executed=False,
    )


def fixture_bound_check():
    c = MPIntervalContext()
    c.prec = 180
    one = _point(c, 1)
    product = _prod(c, [one] * 5, [one] * 5, 4)
    _require(_hi(product) == 16, "Fresh ordinary binomial product fixture failed")
    radius = c.power(c.mpf(2), -32)
    delta, exponential = _fresh_exp_jet(c, _point(c, radius))
    _require(_lo(exponential[0]) > 0, "Fresh exponential increment fixture failed")
    _require(_hi(delta[4]) == radius, "Fresh symmetric increment fixture failed")
    weight_fixture = _weight_derivative_fixture(
        c, radius, c.log(c.mpf(11) / c.mpf(10)))
    _require(weight_fixture["passed"] is True,
             "Fresh normalized-weight derivative fixture failed")
    post_transport_fixture = _post_transport_fixture(c)
    _require(post_transport_fixture["passed"] is True,
             "Fresh post transport fixture failed")
    background_B_fixture = _background_B_fixture(c)
    _require(background_B_fixture["passed"] is True,
             "Fresh background B fixture failed")
    _require(isinstance(SOURCE_GUARDS, tuple), "Source guard recipe must be a tuple")
    return dict(passed=True, ordinary_binomial_order4=16,
                fresh_exp_increment=True,
                weight_derivative_fixture=weight_fixture,
                post_transport_fixture=post_transport_fixture,
                background_B_fixture=background_B_fixture,
                sourceguards_tuple=True)


def self_check():
    return fixture_bound_check()
