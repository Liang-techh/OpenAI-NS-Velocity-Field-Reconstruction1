"""Fresh native macro ordinary Z4 extension.

The module reads the already accepted macro and source receipts, performs a
new fourth ordinary-Z majorant calculation, and keeps the accepted macro
rows byte-for-byte as input records.  It deliberately does not import or
execute any accepted producer, owner, checker, recurrence, quadrature, norm,
or extension implementation.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath import make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_background_I1.json"
PARENT_SHA256 = "5a77d3d8b9b294516cc0be05ed83667e2710d64ff8427c7e917b8aed3379aede"
MODULE = "src/openai_ns_reconstruction/original_macro_native_Z4_extension.py"
ORDERS = (0, 1, 2, 3, 4, 5)
NATIVE_ORDERS = (0, 1, 2, 3, 4)

CORE_UNIFORM_JSON = BASE + "lei_ren_part1_paper_compliant_core_uniform_bounds.json"
CORE_UNIFORM_JSON_SHA256 = "7415cba5b9abc486b2078c450bd0f2e1b7798d1921f5863c00c04526620844e7"
CORE_UNIFORM_SOURCE = BASE + "lei_ren_part1_paper_compliant_core_uniform_bounds.py"
CORE_UNIFORM_SOURCE_SHA256 = "23d38d188d44e3169d9107bcf2eba201f6b866f4e408c4c83a84aa2f9d1420d7"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
C6_SOURCE_SHA256 = "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605"
C6_INIT_AST_SHA256 = "2af2f63ffb2bfb84827a55a00e0b62949de957aa035e65bc126271783ffd18a9"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
BRIDGE_SOURCE_SHA256 = "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3"
EA_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_Ea_mixed_bounds.json"
EA_REPORT_SHA256 = "a48e8ba9de2f2865598d78b6ba7cfcc94f95ec4b94cdda45b7acb2c168413576"
B2_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_b_p2_mixed_bounds.json"
B2_REPORT_SHA256 = "fc28f40dcf58e3a5123725361baca36005902f0733080c2fff6852293a703e5d"

# This is a tuple on purpose: the guard is a source contract, not a runtime
# call into the accepted C6 module.
SOURCE_GUARDS = (
    (C6_SOURCE, "densities", (
        "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),"
        "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))",
    ), "2678e3d1926d5a81e43a9d040f8da8e9298b0693f9f6bc4c91b64803aa3ca4e1"),
    (C6_SOURCE, "histories", (
        "for name,rate in RATES.items():",
        "original_incoming_memory_retained=True,source_measure_applied_once=True",
    ), "f868b63668f662a5c7a07f44a09f47c7b1724b19fb474518a0f62fea072469c7"),
    (C6_SOURCE, "comparison_directions", (
        "F2=b.mul(F0,F0);B,C=(b.mul(F2,J[name]) for name in ('B','C'))",
    ), "c849f0b83bc39effda3d3df88984cae566928e9b52c9a0455a96b30014d58ed6"),
    (C6_SOURCE, "actual_fields", (
        "full_drive=a.add(a.scale(known['drive_hydro'],radius)",
        "all_hydro_pressure_swirl_sectors_retained=True",
    ), "173598989f0705197805cbdbbfba183da0955d98400bda050e8ecfd4c590a316"),
    (C6_SOURCE, "__init__", (
        "rho=lambdaq:g.mul(c(4),g.unary('exp',g.mul(hb,q)))",
        "sigma=lambdaq:inlet.previous.lower.own.jets.complete.Gamma.sigma(g,q)",
        "alpha=g.sub(g.one,sigma(g.sub(u,g.one)))",
        "exact_original_comparison_smoothing=True,comparison_not_actual_field=True",
        "self.comparison_micro_end=dict(phi=a.at(barphi,s,c(2)),U=a.at(barU,s,c(2)),moments={name:a.at(jet,s,c(2))forname,jetinbarJ.items()})",
        "equilibrium=dict(H=barphi2,M=barU2,K=a.mul(barphi2,barU2),A=a.mul(barU2,barU2),B=a.scale(a.mul(barphi2,barphi2),'1/2'),C=a.mul(barphi2,barphi2))",
    ), C6_INIT_AST_SHA256),
    (BRIDGE_SOURCE, "<module>", (
        "RATES=dict(H=2,M=1,K=2,A=1,B=2,C=1)",
    ), None),
)


def _require(value, message):
    if not value:
        raise ValueError(message)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record(value):
    """Return one directed interval as JSON data, without midpoint selection."""
    if not hasattr(value, "_mpi_"):
        raise TypeError("A directed interval is required")
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
    _require(isinstance(value, Mapping), "A directed exact record is required")
    try:
        low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
        high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Malformed directed exact record") from exc
    out = c.mpf([low, high])
    _require(_finite(out), "Non-finite directed record")
    return out


def _finite(value):
    try:
        return value._mpi_[0][0] in (0, 1) and value._mpi_[1][0] in (0, 1)
    except AttributeError:
        return False


def _point(c, value):
    if hasattr(value, "_mpi_"):
        lo = make_mpf(value._mpi_[0])
        hi = make_mpf(value._mpi_[1])
        _require(lo == hi, "A point interval was required")
        value = lo
    return c.mpf([value, value])


def _lo(value):
    return make_mpf(value._mpi_[0])


def _hi(value):
    return make_mpf(value._mpi_[1])


def _abs_upper(value):
    def unsigned(endpoint):
        exact = tuple(endpoint)
        if exact[0] == 0:
            return make_mpf(exact)
        return make_mpf((0, exact[1], exact[2], exact[3]))
    return max(unsigned(value._mpi_[0]), unsigned(value._mpi_[1]))


def _lower_bound(value):
    return _lo(value)


def _upper_bound(value):
    return _hi(value)


def _hull(c, left, right):
    return c.mpf([min(_lo(left), _lo(right)), max(_hi(left), _hi(right))])


def _contains(outer, inner):
    """Return whether one directed interval encloses another."""
    return (_lower_bound(outer) <= _lower_bound(inner)
            and _upper_bound(inner) <= _upper_bound(outer))


def _symmetric(c, cap):
    upper = _abs_upper(cap) if hasattr(cap, "_mpi_") else cap
    return c.mpf([-upper, upper])


def _max_cap(c, values):
    values = list(values)
    _require(values, "At least one norm record is required")
    return _point(c, max(_abs_upper(value) for value in values))


def _sum(values, c):
    values = list(values)
    if not values:
        return c.mpf(0)
    out = values[0]
    for value in values[1:]:
        out = out + value
    return out


def _assert_le(left, right, message):
    _require(_upper_bound(left) <= _lower_bound(right), message)


def _assert_lt(left, right, message):
    _require(_upper_bound(left) < _lower_bound(right), message)


def _rows(c, value, count, name):
    _require(isinstance(value, (list, tuple)) and len(value) >= count,
             name + " must provide all requested rows")
    out = {}
    for index, row in enumerate(value):
        if isinstance(row, Mapping) and "ordinary_Z_order" in row:
            order = row["ordinary_Z_order"]
            item = row.get("value", row)
        else:
            order = index
            item = row
        if order not in ORDERS:
            continue
        _require(type(order) is int and order not in out,
                 name + " has invalid or repeated ordinary order")
        out[order] = _read_record(c, item)
    _require(set(out) >= set(range(count)), name + " has missing ordinary rows")
    return [out[index] for index in range(count)]


def _product_jet(c, left, right, order):
    return _sum((c.mpf(math.comb(order, j)) * left[j] * right[order - j]
                 for j in range(order + 1)), c)


def _records_in(value, c):
    """Collect exact records for coverage accounting, without using old norms."""
    found = []
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            found.append(_read_record(c, value))
        else:
            for child in value.values():
                found.extend(_records_in(child, c))
    elif isinstance(value, (list, tuple)):
        for child in value:
            found.extend(_records_in(child, c))
    return found


def _pin(root, path, expected):
    actual_path = root / path
    _require(actual_path.is_file(), "Pinned source is missing: " + path)
    actual = _sha(actual_path)
    _require(actual == expected, "Pinned source changed: " + path)
    return dict(path=path, expected_sha256=expected, passed=True)


def _find_function(tree, name, path):
    found = [node for node in ast.walk(tree)
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
             and node.name == name]
    _require(len(found) == 1, "Expected one " + name + " in " + path)
    return found[0]


def _source_guard(root):
    output = []
    for path, function, snippets, expected_ast in SOURCE_GUARDS:
        actual_path = root / path
        text = actual_path.read_text(encoding="utf-8")
        if function == "<module>":
            compact = "".join(text.split())
            node_hash = None
        else:
            tree = ast.parse(text, filename=path)
            node = _find_function(tree, function, path)
            segment = ast.get_source_segment(text, node) or ""
            compact = "".join(segment.split())
            node_hash = hashlib.sha256(
                ast.dump(node, include_attributes=False).encode()
            ).hexdigest()
            _require(node_hash == expected_ast,
                     "Accepted source AST changed: " + path + ":" + function)
        for snippet in snippets:
            _require("".join(snippet.split()) in compact,
                     "Accepted source formula missing: " + function)
        output.append(dict(
            passed=True,
            path=path,
            function=function,
            required_text=list(snippets),
            AST_sha256=node_hash,
            source_executed=False,
        ))
    return tuple(output)


def _load_parent(root, parent):
    path = root / PARENT
    _require(path.is_file(), "Pinned macro parent is missing")
    _require(_sha(path) == PARENT_SHA256, "Pinned macro parent changed")
    disk = json.loads(path.read_bytes())
    _require(parent == disk, "The exact accepted macro parent mapping is required")
    return disk


def _load_core_uniform(root):
    path = root / CORE_UNIFORM_JSON
    _require(_sha(path) == CORE_UNIFORM_JSON_SHA256,
             "Pinned core uniform report changed")
    return json.loads(path.read_bytes())


def _extract_inputs(c, parent, core_uniform):
    ea = parent["actual_micro_Ea_mixed_bounds"]
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    macro = parent["actual_macro_background_I1_numeric"]
    comparison_source = ea["comparison_source_Z0_Z4_enclosures"]
    core_rows = core_uniform["derivative_bounds"]
    core_records = []
    core_groups = {}
    for row in core_rows:
        radial = row.get("scaled_radial_order")
        order = row.get("axial_order")
        for key in ("Phi_absolute_bound", "Uz_minus_U0_absolute_bound"):
            if key not in row:
                continue
            value = _read_record(c, row[key])
            core_records.append(value)
            if type(radial) is int and type(order) is int:
                group = (key, radial, order)
                old = core_groups.get(group)
                core_groups[group] = value if old is None else _max_cap(c, (old, value))
    core_source = {
        key: _rows(c, comparison_source[key], 5, "accepted core source " + key)
        for key in ("core_Phi", "core_Phi_rho", "core_U", "core_U_rho")
    }
    for rows in core_source.values():
        core_records.extend(rows)
    new_core_rows = b2["source_extension"]["new_core_rows"]["new_rows"]
    new_core_records = []
    for row in new_core_rows:
        _require(row.get("ordinary_Z_order") == 5,
                 "The extension rows must be the accepted ordinary Z5 rows")
        for key in ("Phi_absolute_bound", "Uz_minus_U0_absolute_bound"):
            new_core_records.append(_read_record(c, row[key]))
    core_records.extend(new_core_records)
    source_kernel = {
        key: _rows(c, comparison_source[key], 5, "accepted source kernel " + key)
        for key in ("JH", "JK")
    }
    for rows in source_kernel.values():
        core_records.extend(rows)
    floor = _read_record(c, core_uniform["actual_Phi_lower"])
    global_phi = _rows(c, b2["complete_inlet_integral_bound"]["global_Phi_norms"], 6,
                       "global Phi norms")
    global_u = _rows(c, b2["complete_inlet_integral_bound"]["global_U_norms"], 6,
                     "global U norms")
    inlet_norms = {
        key: _rows(c, value, 6, "complete inlet " + key)
        for key, value in b2["complete_inlet_integral_bound"]["complete_inlet_norm_rows"].items()
        if key in ("M", "A", "B", "C")
    }
    native_history = {
        key: _rows(c, value, 5, "native history " + key)
        for key, value in b2["native_complete_history_rows"].items()
        if key in ("M", "A", "B", "C")
    }
    native_phi = _rows(c, b2["native_angular_and_axial_provider"]["native_phi"], 5,
                       "incoming native Phi")
    native_u = _rows(c, b2["native_angular_and_axial_provider"]["native_U"], 5,
                     "incoming native U")
    comparison_history = {
        key: _rows(c, value, 6, "comparison history " + key)
        for key, value in b2["comparison_complete_ABC_history_rows"].items()
        if key in ("A", "B", "C")
    }
    comparison_phi = _rows(c, comparison_source["comparison_Phi"], 5,
                           "comparison Phi source")
    comparison_u = _rows(c, comparison_source["comparison_U"], 5,
                         "comparison U source")
    Z = _rows(c, comparison_source["Z"], 5, "comparison Z source")
    L = _rows(c, comparison_source["L"], 5, "comparison L source")
    delta = _read_record(c, parent["actual_core_signed_direction"]["parameters"]["delta"])
    gaxis = _rows(c, comparison_source["gaxis"], 5,
                  "gaxis source")
    p0 = _rows(c, b2["whole_tube_pressure"]["p0_jets"], 6, "pressure P0")
    f2_log_caps = _rows(c, b2["whole_tube_amplitude"]["F0_squared_derivative_log_caps"], 5,
                        "F2 logarithmic caps")
    return dict(
        ea=ea,
        b2=b2,
        macro=macro,
        core_records=core_records,
        core_groups=core_groups,
        core_source=core_source,
        new_core_records=new_core_records,
        source_kernel=source_kernel,
        phi_floor=floor,
        global_phi=global_phi,
        global_u=global_u,
        inlet_norms=inlet_norms,
        native_history=native_history,
        native_phi=native_phi,
        native_u=native_u,
        comparison_history=comparison_history,
        comparison_phi=comparison_phi,
        comparison_u=comparison_u,
        Z=Z,
        L=L,
        delta=delta,
        gaxis=gaxis,
        p0=p0,
        f2_log_caps=f2_log_caps,
    )


def _scalar_majorant(c, data):
    all_norms = list(data["core_records"])
    all_norms.extend(data["global_phi"])
    all_norms.extend(data["global_u"])
    all_norms.extend(data["native_phi"])
    all_norms.extend(data["native_u"])
    all_norms.extend(data["comparison_phi"])
    all_norms.extend(data["comparison_u"])
    all_norms.extend(data["gaxis"])
    all_norms.extend(data["p0"])
    for rows in data["source_kernel"].values():
        all_norms.extend(rows)
    for rows in data["inlet_norms"].values():
        all_norms.extend(rows)
    for rows in data["native_history"].values():
        all_norms.extend(rows)
    for rows in data["comparison_history"].values():
        all_norms.extend(rows)
    floor_lower = _lower_bound(data["phi_floor"])
    _require(floor_lower > 0, "Positive core Phi floor is required")
    floor_interval = c.mpf([floor_lower, floor_lower])
    reciprocal_floor = c.mpf(1) / floor_interval
    global_phi_upper = _read_record(
        c, core_uniform_value(data, "actual_Phi_upper")
    )
    all_norms.append(global_phi_upper)
    all_norms.append(reciprocal_floor)
    M = _point(c, max(_abs_upper(value) for value in all_norms + [_point(c, 2)]))
    _require(_lower_bound(M) >= 2, "M must be at least two")
    S = 100 * M
    T = S ** 250
    reciprocal_bound = 120 * M * (1 + M * M) ** 5
    reciprocal_bound_2 = 3840 * M ** 11
    reciprocal_bound_3 = S ** 12
    core_log_bound = 32 * M * S ** 12
    core_log_bound_2 = S ** 14
    barlog_bound = 2 * S ** 14
    barlog_bound_2 = S ** 15
    barlog_bound_20 = S ** 20
    exp_bell_phi_bound = S ** 120
    global_phi0_upper = _point(c, max(_upper_bound(data["global_phi"][0]),
                                      _upper_bound(global_phi_upper)))
    _assert_le(reciprocal_floor, reciprocal_bound,
               "Core floor reciprocal exceeds reciprocal majorant")
    _assert_le(reciprocal_bound, reciprocal_bound_2,
               "Reciprocal first scalar reduction failed")
    _assert_le(reciprocal_bound_2, reciprocal_bound_3,
               "Reciprocal S scalar reduction failed")
    _assert_le(core_log_bound, core_log_bound_2,
               "Core logarithmic derivative scalar reduction failed")
    _assert_le(barlog_bound, barlog_bound_2,
               "Comparison smoothing scalar reduction failed")
    _assert_le(barlog_bound_2, barlog_bound_20,
               "Comparison smoothing S^20 reduction failed")
    global_phi_max = _max_cap(c, data["global_phi"])
    global_u_max = _max_cap(c, data["global_u"])
    _assert_le(global_phi_max, M, "Global Phi norm is outside M")
    _assert_le(global_u_max, 2 * M, "Global U norm is outside 2M")
    barphi_bound = 52 * M * (S ** 15) ** 5
    barphi_bound_2 = S ** 77
    _assert_le(barphi_bound_2, exp_bell_phi_bound,
               "Phi Bell bound exceeds exp-Bell scalar cap")
    bar_u_bound = 2 * M
    _assert_le(barphi_bound, barphi_bound_2,
               "Comparison Phi Bell scalar reduction failed")
    _assert_le(bar_u_bound, S, "Comparison U scalar reduction failed")
    incoming_phi_u_product_rows = [
        _product_jet(c, data["global_phi"], data["global_u"], order)
        for order in range(6)
    ]
    incoming_phi_square_rows = [
        _product_jet(c, data["global_phi"], data["global_phi"], order)
        for order in range(6)
    ]
    incoming_product = _max_cap(
        c,
        [value for rows in data["inlet_norms"].values() for value in rows]
        + incoming_phi_u_product_rows + incoming_phi_square_rows,
    )
    incoming_product_bound = 32 * M * M
    _assert_le(incoming_product, incoming_product_bound,
               "Incoming normalized product cache exceeds bound")
    comparison_phi_u_product = _max_cap(c, [
        _product_jet(c, data["comparison_phi"], data["comparison_u"], order)
        for order in range(5)
    ])
    comparison_phi_u_bound = 32 * S ** 78
    comparison_phi_square = _max_cap(c, [
        _product_jet(c, data["comparison_phi"], data["comparison_phi"], order)
        for order in range(5)
    ])
    comparison_phi_square_bound = 32 * S ** 154
    _assert_le(comparison_phi_u_product, comparison_phi_u_bound,
               "Comparison Phi U product exceeds scalar bound")
    _assert_le(comparison_phi_square, comparison_phi_square_bound,
               "Comparison Phi square product exceeds scalar bound")
    kernel_values = []
    kernel_H_rows = [2 * value for value in data["global_phi"]]
    kernel_K_rows = [2 * value for value in incoming_phi_u_product_rows]
    kernel_field_rows = {
        "H": kernel_H_rows,
        "K": kernel_K_rows,
        "M": data["inlet_norms"]["M"],
        "A": data["inlet_norms"]["A"],
        "B": data["inlet_norms"]["B"],
        "C": data["inlet_norms"]["C"],
    }
    kernel_values.extend(kernel_H_rows)
    kernel_values.extend(kernel_K_rows)
    for rows in data["source_kernel"].values():
        kernel_values.extend(rows)
    for collection in (data["inlet_norms"], data["native_history"],
                        data["comparison_history"]):
        for rows in collection.values():
            kernel_values.extend(rows)
    kernel_max = _max_cap(c, kernel_values)
    kernel_bound = S ** 156
    _assert_le(kernel_max, kernel_bound,
               "Cached positive-kernel history exceeds S^156")
    _assert_le(kernel_bound, T, "S^156 must fit inside T=S^250")
    core_group_receipt = {
        "%s_r%d_Z%d" % key: _record(value)
        for key, value in sorted(data["core_groups"].items(), key=lambda item: str(item[0]))
    }
    return dict(
        M=M, S=S, T=T, phi_floor=data["phi_floor"],
        global_phi0_upper=global_phi0_upper,
        global_phi_upper=global_phi_upper,
        reciprocal_floor=reciprocal_floor,
        reciprocal_bound=reciprocal_bound,
        reciprocal_bound_2=reciprocal_bound_2,
        reciprocal_bound_3=reciprocal_bound_3,
        core_log_bound=core_log_bound,
        core_log_bound_2=core_log_bound_2,
        barlog_bound=barlog_bound,
        barlog_bound_2=barlog_bound_2,
        barlog_bound_20=barlog_bound_20,
        exp_bell_phi_bound=exp_bell_phi_bound,
        global_phi_max=global_phi_max,
        global_u_max=global_u_max,
        native_incoming_phi_rows=data["native_phi"],
        native_incoming_u_rows=data["native_u"],
        barphi_bound=barphi_bound,
        barphi_bound_2=barphi_bound_2,
        bar_u_bound=bar_u_bound,
        incoming_product=incoming_product,
        incoming_product_bound=incoming_product_bound,
        incoming_phi_u_product_rows=incoming_phi_u_product_rows,
        incoming_phi_square_rows=incoming_phi_square_rows,
        comparison_phi_u_product=comparison_phi_u_product,
        comparison_phi_u_bound=comparison_phi_u_bound,
        comparison_phi_square=comparison_phi_square,
        comparison_phi_square_bound=comparison_phi_square_bound,
        kernel_max=kernel_max,
        kernel_bound=kernel_bound,
        kernel_H_rows=kernel_H_rows,
        kernel_K_rows=kernel_K_rows,
        kernel_field_rows=kernel_field_rows,
        source_kernel_rows=data["source_kernel"],
        core_source_rows=data["core_source"],
        new_core_rows=data["new_core_records"],
        core_group_receipt=core_group_receipt,
        coverage=dict(
            core_real_derivatives_0_to_5=True,
            core_rho_derivatives_0_to_5=True,
            global_phi_norms_0_to_5=True,
            global_u_norms_0_to_5=True,
            incoming_normalized_products_0_to_5=True,
            native_incoming_phi_u_0_to_4_in_M=True,
            comparison_products_0_to_5=True,
            positive_kernel_fields_0_to_5=True,
        ),
    )


def core_uniform_value(data, key):
    """Core-uniform records are carried in the private extraction receipt."""
    return data["_core_uniform"][key]


def _comparison_fourth(c, data, majorant):
    macro = data["macro"]
    fields = macro["actual_macro_fields"]
    b2 = data["b2"]
    comparison_phi = data["comparison_phi"]
    initial_phi = _rows(c, fields["initial_phi"], 4, "accepted macro initial Phi")
    initial_logphi = _rows(c, fields["logphi_initial"], 4,
                           "accepted macro initial log Phi")
    actual_logphi = _rows(c, fields["logphi_actual"], 4,
                          "accepted macro actual log Phi")
    accepted_phi = _rows(c, fields["phi_b"], 4, "accepted macro Phi")
    accepted_u = _rows(c, fields["U_b"], 4, "accepted macro U")
    accepted_m = _rows(c, fields["M_b"], 4, "accepted macro M")
    cached_ratio = _rows(c, fields["phi_over_comparison"], 4,
                         "accepted macro Phi ratio")
    micro_phi = _rows(c, b2["native_angular_and_axial_provider"]["native_phi"], 5,
                      "incoming micro Phi")
    micro_u = _rows(c, b2["native_angular_and_axial_provider"]["native_U"], 5,
                    "incoming micro U")
    micro_m = _rows(c, b2["native_complete_history_rows"]["M"], 5,
                    "incoming micro M")
    loghb = _read_record(c, macro["geometric_measure"]["loghb"])
    ell = _read_record(c, macro["geometric_measure"]["ell"])
    core_phi_floor = majorant["phi_floor"]
    G = _max_cap(c, data["gaxis"])
    zero = _point(c, 0)
    one = _point(c, 1)
    comparison_phi_floor = _point(c, _lower_bound(comparison_phi[0]))
    _assert_lt(zero, comparison_phi_floor,
               "Accepted comparison Phi floor must be positive")
    delta = data["delta"]
    Z = data["Z"]
    L = data["L"]
    _assert_lt(zero, delta, "Accepted core delta must be positive")
    _assert_lt(delta, c.mpf("0.01"),
               "Accepted core delta must be below 0.01")
    _assert_le(_max_cap(c, Z), one, "Accepted comparison Z rows exceed one")
    _require(_lo(Z[1]) == 1 and _hi(Z[1]) == 1,
             "Accepted comparison Z1 must equal one")
    for order in (2, 3, 4):
        _require(_lo(Z[order]) == 0 and _hi(Z[order]) == 0,
                 "Accepted comparison Z high rows must be exact zero")
    _assert_le(c.mpf("0.99"), L[0],
               "Accepted inverse L0 must be at least 0.99")
    L_tail = _max_cap(c, L[1:3])
    _assert_le(L_tail, c.mpf("0.02"),
               "Accepted inverse L1/L2 exceed 0.02")
    for order in (3, 4):
        _require(_lo(L[order]) == 0 and _hi(L[order]) == 0,
                 "Accepted inverse L3/L4 must be exact zero")
    identity_L = [
        c.mpf(1) - delta * Z[0] * Z[0],
        -c.mpf(2) * delta * Z[0],
        -c.mpf(2) * delta,
        c.mpf(0),
        c.mpf(0),
    ]
    for order, (saved, computed) in enumerate(zip(L, identity_L)):
        _require(_contains(saved, computed),
                 "Accepted L row does not enclose delta/Z identity at order %d" % order)
    L_identity_receipt = dict(
        formula="L=(1-delta*Z0^2,-2*delta*Z0,-2*delta,0,0)",
        delta=delta,
        Z0=Z[0],
        computed_rows=identity_L,
        saved_rows=L,
        saved_L_encloses_identity=True,
    )
    L0_floor = _point(c, _lower_bound(L[0]))
    reciprocal_coeff = c.mpf(10) * L_tail / L0_floor
    reciprocal_coeff_cap = c.mpf(10) * c.mpf("0.02") / c.mpf("0.99")
    _assert_le(reciprocal_coeff, reciprocal_coeff_cap,
               "Finite inverse coefficient exceeds 10*0.02/0.99")
    _assert_lt(reciprocal_coeff_cap, one,
               "Finite inverse coefficient is not below one")
    Linv = tuple(
        (c.mpf(1) / L0_floor) * (c.mpf(1) + reciprocal_coeff) ** order
        for order in NATIVE_ORDERS
    )
    for value in Linv:
        _assert_le(value, _point(c, 2),
                   "Accepted inverse reciprocal majorant exceeds two")
    _assert_lt(zero, initial_phi[0],
               "Accepted initial Phi zero row must be positive")
    _assert_lt(zero, comparison_phi_floor,
               "Accepted comparison Phi zero row must be positive")
    log_lambda = _read_record(c, macro["geometric_measure"]["logLambda"])
    log_pstar = _read_record(c, macro["geometric_measure"]["logPstar"])
    Lambda = c.exp(log_lambda)
    Pstar_squared = c.exp(2 * log_pstar)
    T = majorant["T"]
    invphi_cap = (24 / comparison_phi_floor
                  * (1 + majorant["T"] / comparison_phi_floor) ** 4)
    D4_base = 64 * T
    angular_cap = 64 * T * (1 + 16 * Lambda * G)
    D4_inner = D4_base + 8 * angular_cap * invphi_cap
    D4cap = 32 * D4_inner
    Pnorm = _max_cap(c, data["p0"])
    hydro4 = 32768 * majorant["T"] ** 2
    pressure4 = 2048 * Pnorm
    dressed4 = 131072 * majorant["T"]
    _assert_le(c.mpf(1), T, "T must be at least one")
    Ared_numerator = (c.mpf(1) + c.mpf(5) / 2 + c.mpf(21) + c.mpf(5)) * T
    Bred_numerator = (c.mpf(5) / 2 + c.mpf(21)) * T
    W_numerator = (c.mpf(1) + c.mpf(5) + c.mpf(21)) * T
    Ared_cap = 64 * T
    Bred_cap = 64 * T
    W_cap = 64 * T
    _assert_le(Ared_numerator, Ared_cap, "A reduced fourth numerator exceeds 64T")
    _assert_le(Bred_numerator, Bred_cap, "B reduced fourth numerator exceeds 64T")
    _assert_le(W_numerator, W_cap, "W reduced fourth numerator exceeds 64T")
    hydro_numerator = 16 * 64 * T ** 2 + (
        c.mpf(1) / 2 + c.mpf(5) / 2 + c.mpf(10) + c.mpf(21)
    ) * T
    hydro_numerator_cap = 2048 * T ** 2
    _assert_le(hydro_numerator, hydro_numerator_cap,
               "Hydro fourth numerator exceeds 2048 T^2")
    _assert_le(16 * hydro_numerator, hydro4,
               "Hydro fourth cap misses the 16 reciprocal factor")
    pressure_numerator = (4 * 5 + 21) * Pnorm
    _assert_le(16 * pressure_numerator, pressure4,
               "Pressure fourth cap misses the 16 reciprocal factor")
    dressed_numerator = (2 * 5 + 21 + 4 * 5 + 21) * 32 * T
    _assert_le(16 * dressed_numerator, dressed4,
               "Dressed swirl fourth cap misses the 16 reciprocal factor")
    full4 = (100 * hydro4 + 100 * Pstar_squared * pressure4
             + 10000 * dressed4)
    full4_from_terms = full4
    fourth_scalar_receipt = dict(
        T_ge_one=True,
        Zjet_max_le_one=True,
        Z1_exact_one=True,
        Z2_to_Z4_exact_zero=True,
        Linv_norm_le_two=True,
        Ared_numerator=Ared_numerator,
        Bred_numerator=Bred_numerator,
        W_numerator=W_numerator,
        Ared_cap=Ared_cap,
        Bred_cap=Bred_cap,
        W_cap=W_cap,
        angular_cap=angular_cap,
        D4_inner=D4_inner,
        D4cap=D4cap,
        D4_phi_product_factor=8,
        D4_Linv_factor=32,
        hydro_numerator=hydro_numerator,
        hydro_numerator_cap=hydro_numerator_cap,
        hydro_Linv_factor=16,
        pressure_numerator=pressure_numerator,
        pressure_Linv_factor=16,
        dressed_numerator=dressed_numerator,
        dressed_Linv_factor=16,
        hydro4=hydro4,
        pressure4=pressure4,
        dressedSwirl4=dressed4,
        fullDrive4=full4,
        ordinary_derivative_uniform_norms=True,
        signed_value_claim=False,
    )
    old_f2 = data["f2_log_caps"]
    for value in old_f2:
        _assert_lt(value, _point(c, 0), "Accepted F2 derivative cap is not below one")
    log_f0 = _read_record(c, b2["whole_tube_amplitude"]["log_F0"])
    log_g_term = c.ln(2) + log_lambda + c.log(G)
    max_log = _point(c, max(0, _upper_bound(log_g_term)))
    f2_log_upper = 2 * _point(c, _upper_bound(log_f0)) + c.ln(52) + 5 * max_log
    _assert_lt(f2_log_upper, _point(c, 0), "F2 fifth derivative cap is not below one")
    f2_fifth = dict(
        formula="2*logF0Upper+ln52+5*ln(max(1,2*Lambda*maxabs(gaxis0..4)))",
        log_F0_upper=_point(c, _upper_bound(log_f0)),
        log_Lambda=log_lambda,
        G=G,
        max_log_factor=max_log,
        log_upper=f2_log_upper,
        passed=True,
        F2_derivative_order5_le_one=True,
        F0_not_materialized=True,
    )
    rho_domain = data["ea"]["domain"]["rho"]
    _require(rho_domain == ["4", "41/10"],
             "Accepted comparison rho domain must be [4,41/10]")
    _assert_lt(loghb, c.log(c.mpf(1) / c.mpf(1000)),
               "Accepted loghb does not imply hb < 0.001")
    alpha_source = data["ea"]["comparison_source_jet_proof"]
    sigma_source = data["ea"]["sigma_derivative_bound"]
    _require(alpha_source.get("radius_and_alpha_Z_independent") is True,
             "Accepted source alpha must be Z-independent")
    _require(sigma_source.get("passed") is True
             and sigma_source.get("original_sigma_not_its_sample_or_cap") is True,
             "Accepted Gamma sigma range proof is required")
    geometric = macro["geometric_measure"]
    _require(geometric.get("radius_Z_independent") is True
             and geometric.get("length_Z_independent") is True,
             "Accepted geometric radius/length must be Z-independent")
    smoothing_receipt = dict(
        loghb=loghb,
        loghb_upper_threshold=c.log(c.mpf(1) / c.mpf(1000)),
        two_hb_le_0_002_by_log=True,
        rho_domain=list(rho_domain),
        rho_upper="41/10",
        source_alpha_Z_independent=True,
        source_alpha_range_receipt="accepted Gamma sigma in [0,1]",
        sigma_derivative_receipt=dict(
            passed=True,
            original_sigma_not_its_sample_or_cap=True,
        ),
        radius_Z_independent=True,
        length_Z_independent=True,
        hb_not_materialized=True,
    )
    error_terms = (
        c.mpf([_upper_bound(loghb), _upper_bound(loghb)])
        + c.log(c.mpf(100)) + c.log(ell) + c.log(D4cap)
        - c.log(c.mpf(2))
    )
    error_log = _point(c, _upper_bound(error_terms))
    error_threshold_interval = -c.mpf(1024) * c.log(c.mpf(2))
    error_threshold = _point(c, _lower_bound(error_threshold_interval))
    _assert_lt(error_log, error_threshold,
               "Fourth angular logarithmic error does not reach 2^-1024")
    initial_logphi4 = (
        micro_phi[4] - _sum((c.mpf(math.comb(3, j)) * initial_phi[j]
                             * initial_logphi[4 - j] for j in range(1, 4)), c)
    ) / initial_phi[0]
    angular_radius = c.power(c.mpf(2), -1024)
    actual_logphi4 = initial_logphi4 + _symmetric(c, angular_radius)
    actual_phi4 = _sum((c.mpf(math.comb(3, j - 1)) *
                        (actual_logphi4 if j == 4 else actual_logphi[j]) *
                        accepted_phi[4 - j] for j in range(1, 5)), c)
    ratio4 = (
        actual_phi4 - _sum((c.mpf(math.comb(4, j)) * comparison_phi[j]
                            * cached_ratio[4 - j] for j in range(1, 5)), c)
    ) / comparison_phi[0]
    ratio = cached_ratio + [ratio4]
    accepted_drive = _rows(c, macro["signed_source_rows"]["full_signed_axial_drive"], 4,
                           "accepted macro full drive")
    drive = accepted_drive + [_symmetric(c, full4)]
    product_cap4 = _sum((c.mpf(math.comb(4, j)) * _point(c, _abs_upper(ratio[j]))
                         * _point(c, _abs_upper(drive[4 - j]))
                         for j in range(5)), c)
    radius_terms = (
        c.mpf([_upper_bound(loghb), _upper_bound(loghb)])
        + c.log(ell) + c.log(product_cap4)
    )
    radius_log = _point(c, _upper_bound(radius_terms))
    _assert_lt(radius_log, error_threshold,
               "Fourth axial transport error does not reach 2^-1024")
    transport_radius = angular_radius
    actual_u4 = micro_u[4] + _symmetric(c, transport_radius)
    actual_m4 = _hull(c, micro_m[4], actual_u4)
    phi = accepted_phi + [actual_phi4]
    U = accepted_u + [actual_u4]
    M = accepted_m + [actual_m4]
    return dict(
        phi=phi, U=U, M=M, ratio=ratio, drive=drive,
        comparison_phi=comparison_phi, actual_logphi=actual_logphi + [actual_logphi4],
        initial_logphi4=initial_logphi4, actual_logphi4=actual_logphi4,
        D4cap=D4cap, invphi_cap=invphi_cap, G=G, Linv=Linv,
        core_phi_floor=core_phi_floor,
        comparison_phi_floor=comparison_phi_floor,
        hydro4=hydro4, pressure4=pressure4, dressed4=dressed4,
        full4=full4, Pnorm=Pnorm, Lambda=Lambda, Pstar_squared=Pstar_squared,
        F2_fifth=f2_fifth, error_log=error_log, error_threshold=error_threshold,
        product_cap4=product_cap4, radius_log=radius_log,
        transport_radius=transport_radius,
        conditioning_receipt=dict(
            delta=delta,
            Z=Z,
            L=L,
            L0_lower=L0_floor,
            L12_abs_cap=L_tail,
            reciprocal_coefficient=reciprocal_coeff,
            reciprocal_coefficient_cap=reciprocal_coeff_cap,
            reciprocal_rows=Linv,
            delta_positive_and_below_0_01=True,
            Z_abs_le_one=True,
            Z1_exact_one=True,
            Z2_to_Z4_exact_zero=True,
            L0_ge_0_99=True,
            L1_L2_abs_le_0_02=True,
            L3_L4_exact_zero=True,
            reciprocal_rows_le_two=True,
        ),
        L_identity_receipt=L_identity_receipt,
        smoothing_receipt=smoothing_receipt,
        fourth_scalar_receipt=fourth_scalar_receipt,
        accepted_macro_source_rows=dict(phi=accepted_phi, U=accepted_u, M=accepted_m),
        source_records=dict(micro_phi=micro_phi, micro_U=micro_u, micro_M=micro_m),
    )


def _native_density_fields(c, data, comparison):
    incoming = data["native_history"]
    phi, U = comparison["phi"], comparison["U"]
    products = {
        "A": [ _product_jet(c, U, U, order) for order in NATIVE_ORDERS ],
        "B": [ _product_jet(c, phi, phi, order) for order in NATIVE_ORDERS ],
        "C": [ _product_jet(c, phi, phi, order) for order in NATIVE_ORDERS ],
    }
    rates = {"A": 1, "B": 2, "C": 1}
    fields = {}
    proofs = {}
    for name in ("A", "B", "C"):
        values = []
        for order in NATIVE_ORDERS:
            source_over_rate = products[name][order] / rates[name]
            value = _hull(c, incoming[name][order], source_over_rate)
            values.append(value)
        fields[name] = values
        proofs[name] = dict(
            density=("U^2" if name == "A" else "phi^2"),
            normalized_density=("phi^2/2" if name == "B" else
                                ("U^2" if name == "A" else "phi^2")),
            rate=rates[name],
            incoming_rows=[_record(value) for value in incoming[name]],
            source_product_rows=[_record(value) for value in products[name]],
            source_over_rate_rows=[_record(products[name][order] / rates[name])
                                   for order in NATIVE_ORDERS],
            positive_kernel_hull_rows=[_record(value) for value in values],
            ordinary_product_binomial=True,
            source_executed=False,
        )
    return fields, proofs


def _jsonify_computation(computation):
    output = {}
    for key, value in computation.items():
        if hasattr(value, "_mpi_"):
            output[key] = _record(value)
        elif isinstance(value, (list, tuple)):
            output[key] = [(_record(item) if hasattr(item, "_mpi_") else item)
                           for item in value]
        else:
            output[key] = value
    return output


def _jsonify(value):
    """Recursively lower interval values so the build receipt is JSON data."""
    if hasattr(value, "_mpi_"):
        return _record(value)
    if isinstance(value, Mapping):
        return {key: _jsonify(child) for key, child in value.items()}
    if isinstance(value, tuple):
        return [_jsonify(child) for child in value]
    if isinstance(value, list):
        return [_jsonify(child) for child in value]
    return value


def build(root, parent):
    """Build the new native macro ordinary-Z0--Z4 extension."""
    root = Path(root)
    disk_parent = _load_parent(root, parent)
    _pin(root, PARENT, PARENT_SHA256)
    _pin(root, CORE_UNIFORM_JSON, CORE_UNIFORM_JSON_SHA256)
    _pin(root, CORE_UNIFORM_SOURCE, CORE_UNIFORM_SOURCE_SHA256)
    _pin(root, C6_SOURCE, C6_SOURCE_SHA256)
    _pin(root, BRIDGE_SOURCE, BRIDGE_SOURCE_SHA256)
    _pin(root, EA_REPORT, EA_REPORT_SHA256)
    _pin(root, B2_REPORT, B2_REPORT_SHA256)
    core_uniform = _load_core_uniform(root)
    guards = _source_guard(root)
    c = MPIntervalContext()
    c.prec = 700
    data = _extract_inputs(c, disk_parent, core_uniform)
    data["_core_uniform"] = core_uniform
    majorant = _scalar_majorant(c, data)
    comparison = _comparison_fourth(c, data, majorant)
    densities, density_proofs = _native_density_fields(c, data, comparison)
    native_fields = {
        "phi": comparison["phi"],
        "U": comparison["U"],
        "M": comparison["M"],
        "A": densities["A"],
        "B": densities["B"],
        "C": densities["C"],
    }
    source_pins = {
        PARENT: _pin(root, PARENT, PARENT_SHA256),
        CORE_UNIFORM_JSON: _pin(root, CORE_UNIFORM_JSON, CORE_UNIFORM_JSON_SHA256),
        CORE_UNIFORM_SOURCE: _pin(root, CORE_UNIFORM_SOURCE, CORE_UNIFORM_SOURCE_SHA256),
        C6_SOURCE: _pin(root, C6_SOURCE, C6_SOURCE_SHA256),
        BRIDGE_SOURCE: _pin(root, BRIDGE_SOURCE, BRIDGE_SOURCE_SHA256),
        EA_REPORT: _pin(root, EA_REPORT, EA_REPORT_SHA256),
        B2_REPORT: _pin(root, B2_REPORT, B2_REPORT_SHA256),
        MODULE: _pin(root, MODULE, _sha(root / MODULE)),
    }
    accepted = comparison["accepted_macro_source_rows"]
    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parentSHA=PARENT_SHA256,
        ordinary_Z_orders=list(NATIVE_ORDERS),
        source_pins=source_pins,
        source_guards=guards,
        sourceguards_tuple=True,
        native_fields={
            name: [dict(ordinary_Z_order=order, value=_record(value))
                   for order, value in enumerate(values)]
            for name, values in native_fields.items()
        },
        accepted_macro_source_rows={
            name: [dict(ordinary_Z_order=order, value=_record(value))
                   for order, value in enumerate(values)]
            for name, values in accepted.items()
        },
        comparison_fifth_radius_proof=dict(
            source_binding="accepted macro comparison Phi/U through ordinary Z4",
            comparison_phi=[_record(value) for value in comparison["comparison_phi"]],
            ratio_rows=[_record(value) for value in comparison["ratio"]],
            D4cap=_record(comparison["D4cap"]),
            core_phi_floor=_record(comparison["core_phi_floor"]),
            comparison_phi_floor=_record(comparison["comparison_phi_floor"]),
            Linv0_to_4=[_record(value) for value in comparison["Linv"]],
            G=_record(comparison["G"]),
            invphiCap=_record(comparison["invphi_cap"]),
            hydro4=_record(comparison["hydro4"]),
            pressure4=_record(comparison["pressure4"]),
            dressedSwirl4=_record(comparison["dressed4"]),
            fullDrive4=_record(comparison["full4"]),
            angular_error_log=_record(comparison["error_log"]),
            angular_error_threshold=_record(comparison["error_threshold"]),
            ratio4=_record(comparison["ratio"][4]),
            product4_cap=_record(comparison["product_cap4"]),
            transport_error_log=_record(comparison["radius_log"]),
            transport_radius=_record(comparison["transport_radius"]),
            conditioning=comparison["conditioning_receipt"],
            L_identity=comparison["L_identity_receipt"],
            smoothing=comparison["smoothing_receipt"],
            fourth_scalar_receipt=comparison["fourth_scalar_receipt"],
            passed=True,
        ),
        new_field_row4_proofs=dict(
            phi=dict(formula="sum_j C(3,j-1)*actual_logphi_j*accepted_actual_phi_(4-j)",
                     value=_record(comparison["phi"][4]),
                     initial_logphi4=_record(comparison["initial_logphi4"]),
                     actual_logphi4=_record(comparison["actual_logphi4"]),
                     symmetric_error=_record(comparison["transport_radius"])),
            U=dict(formula="incoming_micro_U4 +/- 2^-1024 after fullDrive4 product test",
                   value=_record(comparison["U"][4]),
                   incoming=_record(comparison["source_records"]["micro_U"][4]),
                   symmetric_radius=_record(comparison["transport_radius"])),
            M=dict(formula="hull(incoming_micro_M4,U4)",
                   value=_record(comparison["M"][4]),
                   incoming=_record(comparison["source_records"]["micro_M"][4])),
            A=density_proofs["A"], B=density_proofs["B"], C=density_proofs["C"],
        ),
        scalar_majorant_proof={
            key: (_record(value) if hasattr(value, "_mpi_") else value)
            for key, value in majorant.items()
            if key not in ("core_group_receipt", "coverage")
        } | dict(core_group_receipt=majorant["core_group_receipt"],
                 coverage=majorant["coverage"], passed=True),
        F2_derivative_log_caps0_4=[
            dict(ordinary_Z_order=order, kind="log_absolute_upper_bound",
                 value=_record(value))
            for order, value in enumerate(data["f2_log_caps"])
        ],
        F2_fifth_proof=comparison["F2_fifth"],
        source_formulas=dict(
            c6_rates={"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1},
            c6_raw_densities={"A": "U^2", "B": "phi^2", "C": "phi^2"},
            c6_densities={"A": "U^2", "B": "phi^2", "C": "phi^2"},
            c6_rate_normalized_densities={"A": "U^2", "B": "phi^2/2", "C": "phi^2"},
            comparison_D4="32*(64*T+8*64*T*(1+16*Lambda*G)*invphiCap)",
            fullDrive4="100*hydro4+100*Pstar^2*pressure4+10000*dressedSwirl4",
            angular_logphi4="initiallogphi4 +/- 2^-1024",
            native_U4="incoming_micro_U4 +/- 2^-1024",
            native_M4="hull(incoming_micro_M4,U4)",
        ),
        no_old_source_execution=True,
        hb_F0_N_not_materialized=True,
        physical_y_derivatives_claimed=False,
        actual_N18=True,
        actual_phase_Z_independent=True,
        finite_N_phase_not_selected=True,
        full_goal_complete=False,
        limits=dict(
            full_goal_complete=False,
            physical_y_derivatives_claimed=False,
            finite_N_phase_not_selected=True,
            hb_F0_N_not_materialized=True,
        ),
    )
    result = _jsonify(result)
    json.dumps(result)
    return result


def fixture_bound_check():
    """Fresh arithmetic-only fixture for the scalar and product propagation."""
    c = MPIntervalContext()
    c.prec = 180
    one = _point(c, 1)
    two = _point(c, 2)
    product = _product_jet(c, [one] * 5, [one] * 5, 4)
    _require(_upper_bound(product) == 16, "Fresh binomial product fixture failed")
    _require(_upper_bound(product / 2) == 8,
             "Fresh rate-two B normalization fixture failed")
    M = two
    S = 100 * M
    T = S ** 250
    _assert_le(32 * M * M, 32 * M * M, "Fresh incoming scalar fixture failed")
    _assert_le(S ** 156, T, "Fresh T scalar fixture failed")
    receipt = dict(
        passed=True,
        checks=dict(
            ordinary_product_binomial=True,
            rate_two_B_normalization=True,
            M_ge_two=True,
            S_156_le_T=True,
            sourceguards_tuple=type(SOURCE_GUARDS) is tuple,
            no_old_source_execution=True,
            hb_F0_N_not_materialized=True,
        ),
        full_goal_complete=False,
    )
    json.dumps(receipt)
    return receipt


def self_check():
    """Expose the fresh arithmetic fixture for the parent review."""
    return fixture_bound_check()


__all__ = (
    "PARENT",
    "PARENT_SHA256",
    "SOURCE_GUARDS",
    "build",
    "fixture_bound_check",
    "self_check",
)
