"""Current-N first-switch six-I1 and joined-volume enclosures.

The new chart arithmetic reads saved macro inlet errors and source-owned
first-switch background/loop packets. No accepted builder is executed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_mixed_jets import _logsum

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json"
PARENT_SHA256 = "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c"
MODULE = "src/openai_ns_reconstruction/original_first_switch_finite_N_I1_numeric.py"
C5 = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
Q4 = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
QLINEAR = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
ORDER = 3
I1_ORDER = 1
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")

def _load(root, path):
    return json.loads((root / path).read_bytes())


def _record_value(row, preferred):
    """Accept a direct record or the row wrappers used by saved packets."""
    if row is None:
        return None
    if isinstance(row, dict):
        if preferred in row:
            return row[preferred]
        if preferred == "log_absolute_upper_bound" and "value" in row:
            return row["value"]
    return row


def _ordered_rows(rows, label, length=ORDER + 1):
    require(isinstance(rows, list) and len(rows) >= length,
            "Incomplete " + label + " rows")
    rows = rows[:length]
    require([row.get("ordinary_Z_order") for row in rows
             if isinstance(row, dict)] == list(range(length)),
            "Unordered " + label + " rows")
    return rows


def _read_log_rows(c, rows, label, length=ORDER + 1):
    values = []
    for row in _ordered_rows(rows, label, length):
        if isinstance(row, dict) and row.get("exact_zero"):
            values.append(None)
            continue
        raw = _record_value(row, "log_absolute_upper_bound")
        require(raw is not None, "Missing nonzero logarithmic bound in " + label)
        values.append(read(c, raw))
    return values


def _read_signed_rows(c, rows, label, length):
    values = []
    for row in _ordered_rows(rows, label, length):
        raw = _record_value(row, "value")
        require(raw is not None, "Missing signed value in " + label)
        values.append(read(c, raw))
    return values


def _pack_logs(values):
    out = []
    for index, value in enumerate(values):
        row = {"ordinary_Z_order": index}
        if value is None:
            row.update(exact_zero=True, log_absolute_upper_bound=None)
        else:
            row["log_absolute_upper_bound"] = record(value)
        out.append(row)
    return out


def _pack_values(values):
    return [dict(ordinary_Z_order=index, value=record(value))
            for index, value in enumerate(values)]


def _absolute_log(c, value):
    upper = c.fabs(value).b
    if upper == 0:
        return None
    return c.ln(upper)


def _absolute_logs(c, values):
    return [_absolute_log(c, value) for value in values]


def _point(c, value):
    """Make a directed point from a known upper endpoint."""
    return c.mpf([value.b, value.b])


def _upper_max(c, values):
    values = [value for value in values if value is not None]
    require(values, "At least one absolute log-cap is required")
    return _point(c, max(value.b for value in values))


def _log_abs_scalar(c, scalar):
    upper = c.fabs(scalar).b
    if upper == 0:
        return None
    return c.ln(upper)


def _log_product(c, left, right):
    if left is None or right is None:
        return None
    return left + right


def _log_add(c, *values):
    values = [value for value in values if value is not None]
    return _logsum(c, values)


def _log_linear(c, terms):
    """Absolute log-majorant of a finite linear combination."""
    products = []
    for coefficient, value in terms:
        product = _log_product(c, _log_abs_scalar(c, coefficient), value)
        if product is not None:
            products.append(product)
    return _logsum(c, products)


def _log_scale(c, value, scalar):
    return _log_product(c, _log_abs_scalar(c, scalar), value)


def _log_jet_product(c, left, right, order=I1_ORDER):
    """Ordinary derivative product in log-cap form."""
    out = []
    for k in range(order + 1):
        terms = []
        for j in range(k + 1):
            term = _log_product(c, left[j], right[k - j])
            if term is not None:
                coefficient = math.comb(k, j)
                if coefficient > 1:
                    term += _log_abs_scalar(c, c.mpf(coefficient))
                terms.append(term)
        out.append(_logsum(c, terms))
    return out


def _ordinary_product_values(left, right, order=I1_ORDER):
    """Multiply two signed ordinary derivative rows by Leibniz' rule."""
    out = []
    for k in range(order + 1):
        total = None
        for j in range(k + 1):
            term = left[j] * right[k - j] * math.comb(k, j)
            total = term if total is None else total + term
        out.append(total)
    return out


def _symmetric(c, radius):
    return c.mpf([-radius.b, radius.b])


def _pin(root, parent, path):
    if path == PARENT:
        expected = PARENT_SHA256
    elif path.startswith(BASE):
        expected = parent.get("input_hashes", {}).get(Path(path).name)
    else:
        expected = parent.get("implementation_hashes", {}).get(path)
    require(expected is not None, "Parent has no hash for " + path)
    require(sha(root / path) == expected, "Changed saved source: " + path)
    return dict(path=path, expected_sha256=expected, passed=True)


def _loop_packet_rows(c, loop_packet):
    require(loop_packet.get("passed") is True, "Finite-N loop packet is not accepted")
    require(loop_packet.get("parent_report_sha256") == PARENT_SHA256,
            "Loop packet is tied to a different macro parent")
    require(loop_packet.get("ordinary_Z_orders") == [0, 1, 2, 3],
            "Four ordinary Z rows are required")
    require(loop_packet.get("physical_y_derivatives_claimed") is False,
            "This adapter accepts Z-only loop caps")
    scale = loop_packet["physical_scale_proof"]
    require(scale["passed"] is True and scale["B_quantity"] == "B/Pstar/frozen_E0"
            and scale["deltaU_formula"] ==
            "physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N"
            and scale["N_handle"] == "N18 graph handle"
            and all(scale[key] is True for key in (
                "actual_N18", "actual_N18_is_graph_handle", "actual_N_not_materialized",
                "phase_Z_independent", "source_bound_to_current_first_switch_packet"))
            and scale["physical_y_derivatives_claimed"] is False,
            "Missing current first-switch physical loop units")
    B = _read_log_rows(c, loop_packet["B_over_Pstar_over_frozen_E0_log_caps"],
                       "B_over_Pstar_over_frozen_E0", ORDER + 1)
    require(all(value is not None for value in B), "All B log caps are required")
    pstar_E_raw = _record_value(loop_packet.get("log_Pstar_E_upper"),
                                "log_absolute_upper_bound")
    logPstarE = read(c, pstar_E_raw)
    logN = read(c, _record_value(loop_packet.get("logN"),
                                 "log_absolute_upper_bound"))
    return B, logPstarE, logN


def _recover_delta_Q_logs(c, delta_U, delta_M, z, delta):
    """Recover absolute Q error logs through ordinary Z order two."""
    require(len(delta_U) >= 4 and len(delta_M) >= 4,
            "Q recovery requires physical U and mean rows Z0..Z3")
    h = c.mpf(1) - delta
    d = c.mpf(1) - z * z
    L = c.mpf(1) - delta * z * z
    require(L.a > 0, "Current macro L=1-delta*Z^2 is not positive")
    L1 = -2 * delta * z
    L2 = -2 * delta
    n0 = _log_linear(c, ((2 * z, delta_U[0]),
                         (-h * z, delta_M[0]),
                         (-d, delta_M[1])))
    n1 = _log_linear(c, ((c.mpf(2), delta_U[0]),
                         (2 * z, delta_U[1]),
                         (-h, delta_M[0]),
                         ((1 + delta) * z, delta_M[1]),
                         (-d, delta_M[2])))
    n2 = _log_linear(c, ((c.mpf(4), delta_U[1]),
                         (2 * z, delta_U[2]),
                         (2 * delta, delta_M[1]),
                         ((3 + delta) * z, delta_M[2]),
                         (-d, delta_M[3])))
    invL = -c.ln(L.a)
    q0 = n0 + invL
    q1 = _log_add(c, n1, _log_product(c, _log_abs_scalar(c, L1), q0)) + invL
    q2 = _log_add(c,
                  n2,
                  _log_scale(c, q1, 2 * L1),
                  _log_product(c, _log_abs_scalar(c, L2), q0)) + invL
    return [q0, q1, q2], dict(
        passed=True,
        ordinary_derivatives=True,
        no_top_row_padding=True,
        positive_L=record(L),
        L_Z=record(L1),
        L_ZZ=record(L2),
        numerator_log_caps=_pack_logs((n0, n1, n2)),
        formula=(
            "deltaQ=(2*Z*deltaU-(1-delta)*Z*deltaM-"
            "(1-Z^2)*deltaM_Z)/(1-delta*Z^2)"
        ),
    )


def _difference_logs(c, U, q, U_b, q_b):
    qz = q[1:3]
    qbz = q_b[1:3]
    Uz = U[1:3]
    Ubz = U_b[1:3]
    qb2 = _log_jet_product(c, q_b[:2], q[:2])
    q2 = _log_jet_product(c, q[:2], q[:2])
    q_squared = [_log_add(c, _log_scale(c, qb2[k], c.mpf(2)), q2[k])
                 for k in range(I1_ORDER + 1)]
    U_q_Z = []
    for k in range(I1_ORDER + 1):
        U_q_Z.append(_log_add(
            c,
            _log_jet_product(c, U_b[:2], qz)[k],
            _log_jet_product(c, U[:2], qbz)[k],
            _log_jet_product(c, U[:2], qz)[k],
        ))
    U_squared_b = _log_jet_product(c, U_b[:2], U[:2])
    U_squared_d = _log_jet_product(c, U[:2], U[:2])
    U_squared = [_log_add(c, _log_scale(c, U_squared_b[k], c.mpf(2)),
                          U_squared_d[k]) for k in range(I1_ORDER + 1)]
    UU_Z = []
    for k in range(I1_ORDER + 1):
        UU_Z.append(_log_add(
            c,
            _log_jet_product(c, U_b[:2], Uz)[k],
            _log_jet_product(c, U[:2], Ubz)[k],
            _log_jet_product(c, U[:2], Uz)[k],
        ))
    return dict(
        q=q[:2],
        q_Z=qz,
        q_squared=q_squared,
        U_q_Z=U_q_Z,
        U_squared=U_squared,
        U_U_Z=UU_Z,
    )


def _full_I1_and_D(c, background, difference, radius):
    signed = background["background_I1"]["six_integrals"]
    base = {name: _read_signed_rows(c, signed[name], "signed " + name, I1_ORDER + 1)
            for name in I1_NAMES}
    full_values = {name: [value + _symmetric(c, radius) for value in base[name]]
                   for name in I1_NAMES}
    coefficients = background["background_I1"]["coefficientsC1rows"]
    coeff = {name: _read_signed_rows(c, coefficients[name], "coefficient " + name,
                                     I1_ORDER + 1)
             for name in coefficients}
    squared_coefficient = coeff.get("two_h_Z2_over_L2", coeff.get("two_hZ2_over_L2"))
    require(squared_coefficient is not None, "Missing 2*h*Z^2/L^2 coefficient")
    coefficient_products = {
        "q": _ordinary_product_values(coeff["minus_inv_L"], full_values["q"]),
        "q_Z": _ordinary_product_values(coeff["hZ_over_2L"], full_values["q_Z"]),
        "q_squared": _ordinary_product_values(coeff["minus_half"],
                                               full_values["q_squared"]),
        "U_q_Z": _ordinary_product_values(coeff["d_over_L"], full_values["U_q_Z"]),
        "U_squared": _ordinary_product_values(squared_coefficient,
                                               full_values["U_squared"]),
        "U_U_Z": _ordinary_product_values(coeff["two_Zd_over_L2"],
                                           full_values["U_U_Z"]),
    }
    D_full = []
    for k in range(I1_ORDER + 1):
        terms = [
            coefficient_products["q"][k],
            coefficient_products["q_Z"][k],
            coefficient_products["q_squared"][k],
            coefficient_products["U_q_Z"][k],
            coefficient_products["U_squared"][k],
            coefficient_products["U_U_Z"][k],
        ]
        total = terms[0]
        for term in terms[1:]:
            total += term
        D_full.append(total)
    return (
        dict(
            passed=True,
            normalized_by_Wswitch=True,
            six_integrals={name: _pack_values(full_values[name]) for name in I1_NAMES},
            finite_N_absolute_error_log_bounds={
                name: _pack_logs(difference[name]) for name in I1_NAMES
            },
            signed_background_unchanged=True,
            symmetric_radius_added_once=record(radius),
            error_radius_is_not_a_selected_value=True,
        ),
        dict(
            passed=True,
            rows=_pack_values(D_full),
            D_rows=_pack_values(D_full),
            coefficientsC1rows=background["background_I1"]["coefficientsC1rows"],
            computed_from_full_I1_once=True,
            separate_D_error_not_added=True,
            normalized_by_Wswitch=True,
        ),
    )


def _guards(root):
    return dict(
        finite_N=source_guard(root, C5, "finite_N_functions", (
            "s=a.scale(loop['B_over_Pstar'],invN)", "V_N=a.add(V,s)",
        )),
        nonreset_mean=source_guard(root, C5, "transport", (
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "predecessor_memory_never_reset=True",
        )),
        physical_fields=source_guard(root, Q4, "physical_pressure_fields", (
            "U=a.scale(fields['V'],Pstar)", "m=a.scale(fields['m'],Pstar)",
        )),
        linear_Q=source_guard(root, QLINEAR, "new_controls", (
            "op=lambda u,m:(2*z*u-(1-de)*z*m-(1-z*z)*s.diff(m,z))/(1-de*z*z)",
        )),
    )


def _joined_volume(c, parent, switch_full, loghb, radius, threshold):
    # These are the complete signed first-switch integrals in Wswitch units.
    # The physical conversion to Wmacro is the exact unmaterialized hbS.
    switch = {name: _read_signed_rows(c, rows, "full first-switch " + name, 2)
              for name, rows in switch_full["six_integrals"].items()}
    whole_logs = {name: [None if _absolute_log(c, value) is None else
                        loghb + _absolute_log(c, value) for value in rows]
                  for name, rows in switch.items()}
    require(all(v is None or v.b < threshold.a
                for rows in whole_logs.values() for v in rows),
            "Complete first-switch volume does not fit the common-unit radius")
    macro = parent["actual_macro_finite_N_I1_numeric"]["I1_full"]
    pseudo = dict(background_I1=dict(
        six_integrals=macro["six_integrals"],
        coefficientsC1rows=parent["actual_macro_background_I1_numeric"]
                                 ["macro_I1"]["coefficientsC1rows"]))
    joined, D = _full_I1_and_D(c, pseudo, whole_logs, radius)
    joined["first_switch_full_volume_absolute_log_bounds"] = joined.pop(
        "finite_N_absolute_error_log_bounds")
    for item in (joined, D):
        item["normalized_by_Wswitch"] = False
        item["normalized_by_Wmacro"] = True
    joined.update(
        radius_purpose="complete first-switch signed volume in Wmacro units, once",
        source_relation="Ijoined/Wmacro=Imacro/Wmacro+exp(loghbS)*(Ifirst/Wswitch)",
        macro_finite_N_radius_not_added_again=True,
        first_switch_volume_not_selected_as_zero=True,
        first_switch_log_physical_relative_weight=record(loghb),
        signed_sum_is_interval_enclosure=True,
        signed_joined_center_value_claimed=False,
        only_volume_enclosed=True, outer_RQy_retained=True,
    )
    return dict(passed=True, I1=joined, D=D,
                full_first_switch_I1_in_own_units_retained=True,
                common_physical_weight_hbS_not_materialized=True,
                outer_RQy_retained=True,
                complete_pressure_enclosed=False)


def build(root, parent, background, loop_packet):
    root = Path(root)
    require(sha(root / PARENT) == PARENT_SHA256 and _load(root, PARENT) == parent,
            "Exact accepted macro finite-N parent required")
    require(background["passed"] is True
            and background["parent_report_sha256"] == PARENT_SHA256,
            "Source-owned first-switch background required")
    c = MPIntervalContext()
    c.prec = 1200
    B, logPE, logN = _loop_packet_rows(c, loop_packet)
    actual_N = read(c, parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]["N"])
    require(logN.a == actual_N.a and logN.b == actual_N.b,
            "First-switch loop must use the unchanged original huge N")
    geometry = background["geometry_measure"]
    require(geometry["hbS_source"] == "R100_Rm_C6.native_parameter('hbS')"
            and geometry["upstream_width_source"] == "parameter_geometry_guard.shared.logh"
            and geometry["source_nodes_distinct"] is True
            and geometry["source_node_aliasing_forbidden"] is True
            and geometry["physical_measure"] == "R^2*hbS*ds"
            and geometry["Wswitch_definition"] == "10000*hbS"
            and geometry["normalized_K_formula"] == "integral_0^1 exp(2*hbS*s) ds"
            and all(geometry["normalized_K_proof"][key] is True for key in (
                "log_width_plus_log4_lt_log_2minus1024", "two_hbS_le_half_2minus1024",
                "exp_half_width_lt_one_plus_2minus1024", "K_is_enclosure_not_selected_value",
                "hbS_not_materialized")), "Missing source-owned physical width/measure proof")
    require(background["native_input_source"]["path"] ==
            "actual_macro_loop_Z3_inputs.native_Z4_extension.native_fields"
            and background["native_input_source"]["true_macro_endpoint_uniform_box"] is True
            and background["Q_b_proof"]["formula"] ==
            "Q_b=(2*Z*U-(1-delta)*Z*M-d*M_Z)/L"
            and background["source_contract"]["incoming_memory_retained"] is True,
            "Missing native U/M/Q source binding")
    K = read(c, geometry["normalized_K"])
    loghb = read(c, geometry["loghbS"])
    radius = c.mpf(2) ** -1024
    require(K.a >= 1 and K.b <= (1 + radius).b
            and (loghb + c.ln(c.mpf(4))).b < (-1024*c.ln(c.mpf(2))).a,
            "First-switch normalized measure enclosure failed")
    delta = read(c, parent["actual_core_signed_direction"]["parameters"]["delta"])
    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    dU = [logPE + value - logN for value in B]
    incoming = _read_log_rows(c, parent["actual_macro_finite_N_I1_numeric"]
                             ["physical_error_log_bounds"]["delta_M"],
                             "macro physical delta M", 4)
    dM = [_upper_max(c, (incoming[k], dU[k])) for k in range(4)]
    dQ, Qproof = _recover_delta_Q_logs(c, dU, dM, z, delta)
    Ub = _absolute_logs(c, _read_signed_rows(c, background["native_fields"]["U"],
                                           "first-switch U", 4))
    Qb = _absolute_logs(c, _read_signed_rows(c, background["Q_b"], "first-switch Q", 3))
    density = _difference_logs(c, dU, dQ, Ub, Qb)
    defects = {name: [None if value is None else value + c.ln(K.b)
                       for value in rows] for name, rows in density.items()}
    threshold = -1024 * c.ln(c.mpf(2))
    require(all(value is None or value.b < threshold.a
                for rows in defects.values() for value in rows),
            "Current-N first-switch six-I1 defect exceeds 2^-1024")
    full, D = _full_I1_and_D(c, background, defects, radius)
    joined = _joined_volume(c, parent, full, loghb, radius, threshold)
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256)}
    for packet in (background, loop_packet):
        for key, value in packet["source_pins"].items():
            expected = value["expected_sha256"] if isinstance(value, dict) else value
            existing = pins.get(key, value)
            previous = existing["expected_sha256"] if isinstance(existing, dict) else existing
            require(previous == expected, "Conflicting new source pin: " + key)
            pins[key] = value
    for path in (C5, Q4, QLINEAR):
        pins[path] = _pin(root, parent, path)
    pins[MODULE] = dict(path=MODULE, expected_sha256=sha(root / MODULE), passed=True)
    return dict(
        passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        physical_y_derivatives_claimed=False,
        source_pins=pins, source_guards=_guards(root),
        source_family=parent["source_family"],
        source_contract=dict(
            passed=True, chart="first_switch", ordinary_Z_orders=[0, 1, 2, 3],
            N_handle=18, N_handle_is_not_numeric_N=True,
            physical_J_handle=9, physical_J_name="hbS",
            phase_Z_independent=True, physical_y_derivatives_claimed=False,
            physical_delta_U_local_only=True, no_incoming_U_added=True,
            incoming_physical_delta_M_from_macro=True,
            incoming_delta_M_path="parent.actual_macro_finite_N_I1_numeric.physical_error_log_bounds.delta_M",
            nonreset_mean_positive_kernel=True, physical_Pstar_applied_once=True,
            physical_R_not_a_mean_scale=True, Q_linear_source=Qproof,
            physical_measure="R^2*hbS ds, Wswitch=10000*hbS",
            finite_N_radius_added_once_to_background=True,
            signed_error_cap_not_selected_as_value=True,
        ),
        physical_error_log_bounds=dict(delta_U=_pack_logs(dU), delta_M=_pack_logs(dM),
            delta_Q=_pack_logs(dQ), logN=record(logN), log_Pstar_E_upper=record(logPE)),
        six_I1_difference_error_log_bounds={name: _pack_logs(rows) for name, rows in defects.items()},
        I1_full=full, D_full=D, joined_macro_first_switch_volume=joined,
        signed_I1_radius=record(radius), I1_error_threshold=record(threshold),
        all_normalized_I1_errors_below_2_power_minus1024=True,
        limits=dict(first_switch_full_finite_N_six_I1_C1_enclosed=True,
            outer_RQy_retained=True, outer_RQy_evaluated=False,
            complete_pressure_global_closure=False, global_five_moment_closure=False,
            physical_y_derivatives_claimed=False, old_212_and_signed_errors_unchanged=True),
        accepted_ancestors_replayed=False, full_goal_complete=False,
    )


__all__ = ["build"]
