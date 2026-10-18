"""Current-N second-switch and post-power signed volume enclosures.

This stage reads the accepted first-switch data and new chart-specific packets.
It never executes an accepted producer or arithmetic owner.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_mixed_jets import _logsum

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"
PARENT_SHA256 = "2cc8b4ff9670858564297bc1c27d5322c9af7d2cb85381939dc3e37a3eea346a"
MODULE = "src/openai_ns_reconstruction/original_switch_post_finite_N_I1_numeric.py"
C5 = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
Q4 = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
QLINEAR = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
ORDER = 3
I1_ORDER = 1
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
CHARTS = ("second_switch", "post_switch")

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
                "phase_Z_independent", "source_bound_to_current_switch_post_packet"))
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
            "my=a.sub(U,m)", "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        reduced_pressure=source_guard(root, Q4, "reduced_pressure_S", (
            "linear=div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy))",
            "quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))",
            "axial=div(a.mul(U,a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2))))",
        )),
        linear_Q=source_guard(root, QLINEAR, "new_controls", (
            "op=lambda u,m:(2*z*u-(1-de)*z*m-(1-z*z)*s.diff(m,z))/(1-de*z*z)",
        )),
    )


def _pressure_identity():
    """Derive B+D from the actual Q4 S, using m_y=U-m.

    This is an arbitrary-symbol identity, not a numerical endpoint fixture.
    It remains valid for finite-N U with nonzero physical U_y.
    """
    import sympy as sp
    z,de,U,Uz,Uy,m,mz,Q,Qz,Qy=sp.symbols("z de U Uz Uy m mz Q Qz Qy")
    L=1-de*z*z
    d=1-z*z
    h=1-de
    q_source=(2*z*U-h*z*m-d*mz)/L
    qy_source=(2*z*Uy-h*z*(U-m)-d*(Uz-mz))/L
    relation=(2*z*Uy+(1+de)*z*U-d*Uz)/L
    require(sp.factor(q_source+qy_source-relation)==0,
            "Q4 mean relation did not yield the Q+Q_y identity")
    S=(Q+h*z*Qz/2+Qy)/L+Q*(Q/2+Qy) \
      +U*(-2*z*Q+d*Qz-2*z*Qy)/L
    G=Q/L+Q*Q/2-2*z*z*U*U/(L*L)
    G_y=Qy/L+Q*Qy-4*z*z*U*Uy/(L*L)
    D=-Q/L+h*z*Qz/(2*L)-Q*Q/2+d*U*Qz/L \
      +2*h*z*z*U*U/(L*L)+2*z*d*U*Uz/(L*L)
    remainder=sp.factor((S-2*G-G_y-D).subs(Qy,relation-Q))
    require(remainder==0,"Actual Q4 S is not the claimed boundary-plus-D identity")
    return dict(passed=True,engine="sympy arbitrary-symbol rational identity",
        Q4_reduced_pressure_source_guarded=True,mean_relation="m_y=U-m",
        Q_plus_Qy="(2*Z*U_y+(1+delta)*Z*U-d*U_Z)/Lz",
        boundary_integrand="G=Q/Lz+Q^2/2-2*Z^2*U^2/Lz^2",
        identity="S=2*G+G_y+D_density; (R^2)_y=2*R^2",
        symbolic_remainder=str(remainder),finite_N_Uy_not_assumed_zero=True,
        ordinary_Z_identity=True,physical_y_numeric_bounds_claimed=False)


def _merge_pins(root, parent, *packets):
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256)}
    for packet in packets:
        require(packet["passed"] is True and
                packet["parent_report_sha256"] == PARENT_SHA256,
                "New packet must belong to the current first-switch parent")
        for key, value in packet["source_pins"].items():
            expected = value["expected_sha256"] if isinstance(value, dict) else value
            previous = pins.get(key, value)
            previous = previous["expected_sha256"] if isinstance(previous, dict) else previous
            require(previous == expected and sha(root/key) == expected,
                    "New source pin conflict or mutation: " + key)
            pins[key] = value
    for path in (C5, Q4, QLINEAR):
        pins[path] = _pin(root, parent, path)
    pins[MODULE] = dict(path=MODULE, expected_sha256=sha(root/MODULE), passed=True)
    return pins


def _normalization(packet, name):
    switch = name == "second_switch"
    packet["normalized_by_Wswitch"] = switch
    packet["normalized_by_Wmacro"] = not switch
    return packet


def _full_B_and_S1(c, name, background, densities, dQ, delta, z, D,
                   radius, threshold):
    """Enclose the boundary functional using its correlated analytic value.

    The background evaluates 2*C0*J0+C1*J1 before interval boxing.  Only the
    finite-N difference is bounded by two endpoint absolute values here.
    """
    L = 1 - delta*z*z
    invL = [1/L, 2*delta*z/(L*L)]
    two_Z2_L2 = [2*z*z/(L*L), 4*z/(L*L)+8*delta*z*z*z/(L*L*L)]
    linear = _log_jet_product(c, _absolute_logs(c, invL), dQ[:2])
    squared_U = _log_jet_product(c, _absolute_logs(c, two_Z2_L2),
                                 densities["U_squared"])
    deltaG = [_log_add(c, linear[k],
                       _log_scale(c, densities["q_squared"][k], c.mpf(1)/2),
                       squared_U[k]) for k in range(2)]
    geometry = background["geometry"]
    log_endpoint_factor = c.ln(c.mpf(2)) + 2*c.ln(read(c, geometry["R_upper"])) \
                          - c.ln(c.mpf(10000))
    if name == "second_switch":
        log_endpoint_factor -= read(c, geometry["loghbS"])
    errors = [None if value is None else value+log_endpoint_factor for value in deltaG]
    require(all(value is None or value.b < threshold.a for value in errors),
            "Current-N B boundary-functional error is not below 2^-1024: " + name)
    base = _read_signed_rows(c, background["background_B_rows"],
                            name+" analytic correlated B", 2)
    full = [value+_symmetric(c, radius) for value in base]
    D_rows = _read_signed_rows(c, D["rows"], name+" finite-N D", 2)
    B = _normalization(dict(passed=True, rows=_pack_values(full),
        finite_N_absolute_error_log_bounds=_pack_logs(errors),
        difference_G_log_bounds=_pack_logs(deltaG),
        endpoint_pair_log_factor=record(log_endpoint_factor),
        difference_formula=("deltaG=deltaQ/Lz+delta(Q^2)/2"
                            "-(2*Z^2/Lz^2)*delta(U^2); ordinary C1 products"),
        endpoint_error_formula="abs(deltaB)/W<=2*R_upper^2/W*sup(abs(deltaG))",
        second_W_has_actual_hbS_in_denominator=(name == "second_switch"),
        symmetric_radius_added_once=record(radius),
        error_radius_is_not_a_selected_value=True,
        correlated_background_not_endpoint_interval_subtraction=True,
        background_identity="B/W=2*C0*J0+C1*J1; C2 cancels before boxing"), name)
    S1 = _normalization(dict(passed=True,
        rows=_pack_values([full[k]+D_rows[k] for k in range(2)]),
        formula="S1_full=B_full+D_full", computed_from_B_and_D_once=True,
        separate_S1_error_not_added=True, physical_y_derivatives_claimed=False), name)
    return B, S1


def _chart(c, parent, name, background, loop, incoming, radius, threshold):
    require(background["passed"] is True, "New native chart did not pass: " + name)
    B, logPE, logN = _loop_packet_rows(c, loop)
    actual_N = read(c, parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]["N"])
    require(logN.a == actual_N.a and logN.b == actual_N.b,
            "Each chart must retain the unchanged original N")
    geometry = background["geometry"]
    normalization = "Wswitch" if name == "second_switch" else "Wmacro"
    require(geometry["normalization"] == normalization,
            "Wrong native chart physical normalization")
    mass = read(c, geometry["mass"])
    require(mass.a > 0 and geometry["mass"] == geometry["normalized_moments"][0],
            "Pointwise errors and signed integrals must use the same positive mass")
    delta = read(c, parent["actual_core_signed_direction"]["parameters"]["delta"])
    z = c.mpf([c.mpf(37)/100, c.mpf(93)/250])
    dU = [logPE + value - logN for value in B]
    dM = [_upper_max(c, (incoming[k], dU[k])) for k in range(4)]
    dQ, qproof = _recover_delta_Q_logs(c, dU, dM, z, delta)
    ub = _absolute_logs(c, _read_signed_rows(c, background["native_fields"]["U"],
                                           name + " physical U", 4))
    qb = _absolute_logs(c, _read_signed_rows(c, background["Q_b"], name + " physical Q", 3))
    densities = _difference_logs(c, dU, dQ, ub, qb)
    defects = {key: [None if value is None else value+c.ln(mass.b) for value in rows]
               for key, rows in densities.items()}
    require(all(value is None or value.b < threshold.a
                for rows in defects.values() for value in rows),
            "Current-N normalized volume error is not below 2^-1024: " + name)
    full, D = _full_I1_and_D(c, background, defects, radius)
    _normalization(full, name)
    _normalization(D, name)
    B_full, S1_full = _full_B_and_S1(c, name, background, densities, dQ,
                                    delta, z, D, radius, threshold)
    errors = dict(delta_U=_pack_logs(dU), delta_M=_pack_logs(dM), delta_Q=_pack_logs(dQ),
                  logN=record(logN), log_Pstar_E_upper=record(logPE))
    return dict(passed=True, chart=name, parent_report_sha256=PARENT_SHA256,
        ordinary_Z_orders=[0,1], physical_y_derivatives_claimed=False,
        I1_full=full, D_full=D, B_full=B_full, S1_full=S1_full,
        physical_error_log_bounds=errors,
        six_I1_difference_error_log_bounds={k:_pack_logs(v) for k,v in defects.items()},
        signed_I1_radius=record(radius), I1_error_threshold=record(threshold),
        normalized_mass=record(mass), geometry=geometry,
        source_contract=dict(passed=True, chart=name, common_N_handle=18,
            N_not_materialized=True, phase_Z_independent=True,
            physical_Pstar_applied_once=True, physical_R_not_a_mean_scale=True,
            physical_delta_U_local_only=True, incoming_delta_U_not_added=True,
            nonreset_positive_kernel_mean=True,
            incoming_mean_chart="first_switch" if name=="second_switch" else "second_switch",
            mean_bound="max(incoming deltaM,local deltaU); positive kernel mass<=1",
            Q_linear_source=qproof, normalization=normalization,
            source_chart_in_Q4="second_switch" if name=="second_switch" else "post_power",
            true_physical_mass_once=True, finite_N_radius_added_once=True,
            physical_y_derivatives_claimed=False),
        limits=dict(six_I1_finite_N_C1_enclosed=True, B_and_S1_finite_N_C1_enclosed=True,
            D_is_not_full_S1=True,
            outer_RQy_retained=True, outer_RQy_evaluated=False,
            B_endpoint_numeric=True, complete_pressure_global_closure=False,
            physical_y_derivatives_claimed=False, global_five_moment_closure=False))


def _joined_volume(c, parent, background, charts, radius, threshold):
    old = parent["actual_first_switch_finite_N_I1_numeric"]["joined_macro_first_switch_volume"]
    require(old["passed"] is True and old["complete_pressure_enclosed"] is False,
            "Inherited joined volume must remain a partial pressure term")
    old_values = {k:_read_signed_rows(c,v,"inherited macro/first "+k,2)
                  for k,v in old["I1"]["six_integrals"].items()}
    second = {k:_read_signed_rows(c,v,"second full "+k,2)
              for k,v in charts["second_switch"]["I1_full"]["six_integrals"].items()}
    post = {k:_read_signed_rows(c,v,"post full "+k,2)
            for k,v in charts["post_switch"]["I1_full"]["six_integrals"].items()}
    loghb = read(c, background["charts"]["second_switch"]["geometry"]["loghbS"])
    whole = {k:[None if _absolute_log(c,x) is None else loghb+_absolute_log(c,x)
                for x in rows] for k,rows in second.items()}
    require(all(x is None or x.b < threshold.a for rows in whole.values() for x in rows),
            "Whole second-switch physical volume does not fit the common-unit radius")
    combined = {k:[old_values[k][j]+_symmetric(c,radius)+post[k][j] for j in range(2)]
                for k in I1_NAMES}
    packed = {k:_pack_values(v) for k,v in combined.items()}
    coefficients = background["charts"]["post_switch"]["background_I1"]["coefficientsC1rows"]
    pseudo = dict(background_I1=dict(six_integrals=packed, coefficientsC1rows=coefficients))
    # A zero *additional* radius computes D from the already enclosed sum.
    # All finite-N and whole-second radii are already present in `combined`.
    unused, D = _full_I1_and_D(c, pseudo, {k:[None,None] for k in I1_NAMES}, c.mpf(0))
    D["normalized_by_Wswitch"] = False
    D["normalized_by_Wmacro"] = True
    return dict(passed=True, ordinary_Z_orders=[0,1], normalized_by_Wmacro=True,
        I1=dict(passed=True, six_integrals=packed, normalized_by_Wmacro=True,
            normalized_by_Wswitch=False, signed_sum_is_interval_enclosure=True,
            signed_center_value_claimed=False,
            whole_second_switch_absolute_log_bounds={k:_pack_logs(v) for k,v in whole.items()},
            whole_second_switch_radius_added_once=record(radius),
            physical_second_to_macro_log_weight=record(loghb),
            macro_and_first_radii_not_repeated=True, post_finite_N_radius_not_repeated=True,
            new_extra_joined_radius=False, error_cap_is_not_a_selected_value=True),
        D=D, coefficientsC1rows=coefficients,
        formula="Ijoined/Wmacro=Iparent/Wmacro+exp(loghbS)*Isecond/Wswitch+Ipost/Wmacro",
        parent_joined_volume_retained=True, second_own_signed_units_retained=True,
        post_physical_weight_has_no_extra_hbS=True,
        whole_micro_volume_not_part_of_this_numeric_packet=True,
        outer_RQy_retained=True, B_endpoint_numeric=False,
        complete_pressure_enclosed=False, global_five_moment_closure=False)


def build(root, parent, background, loop_packet):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and _load(root,PARENT)==parent,
            "Exact accepted first-switch parent required")
    pins=_merge_pins(root,parent,background,loop_packet)
    c=MPIntervalContext()
    c.prec=1200
    radius=c.mpf(2)**-1024
    threshold=-1024*c.ln(c.mpf(2))
    incoming=_read_log_rows(c,parent["actual_first_switch_finite_N_I1_numeric"]
        ["physical_error_log_bounds"]["delta_M"],"first-switch incoming deltaM",4)
    charts={}
    for name in CHARTS:
        charts[name]=_chart(c,parent,name,background["charts"][name],
                            loop_packet["charts"][name],incoming,radius,threshold)
        incoming=_read_log_rows(c,charts[name]["physical_error_log_bounds"]["delta_M"],
                                 name+" outgoing deltaM",4)
    joined=_joined_volume(c,parent,background,charts,radius,threshold)
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        source_pins=pins,source_guards=_guards(root), charts=charts,
        pressure_B_plus_D_identity=_pressure_identity(),
        joined_macro_switch_post_volume=joined,
        physical_y_derivatives_claimed=False,
        all_normalized_I1_errors_below_2_power_minus1024=True,
        source_contract=dict(passed=True,common_N_handle=18,phase_Z_independent=True,
            actual_source_chart_mapping={"second_switch":"second_switch","post_switch":"post_power"},
            first_to_second_to_post_mean_memory_retained=True,
            delta_U_has_no_incoming=True,physical_Pstar_applied_once=True,
            separate_Wswitch_and_Wmacro=True,physical_y_derivatives_claimed=False),
        limits=dict(outer_RQy_retained=True,outer_RQy_evaluated=False,
            chart_B_and_S1_numeric=True,joined_parent_B_numeric=False,
            D_is_not_full_S1=True,complete_pressure_global_closure=False,
            global_five_moment_closure=False,old_212_signed_errors_unchanged=True,
            physical_y_derivatives_claimed=False),
        accepted_ancestors_replayed=False,full_goal_complete=False)


__all__=["build"]
