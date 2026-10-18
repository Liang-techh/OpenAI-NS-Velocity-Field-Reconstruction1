"""Current-``N`` finite-error enclosures for long/reference ``I1`` rows.

This adapter consumes the accepted native constant-``U`` background and two
source-only current-``N`` loop packets.  It performs only directed logarithmic
error arithmetic.  The native signed rows stay the selected centres of their
interval enclosures; finite-``N`` rows are absolute error caps and are never
treated as signed values.  The two charts use their own ``W=R_upper**2``
normalization.
"""
from __future__ import annotations

import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_mixed_jets import _logsum


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_background_I1.json"
PARENT_SHA256 = "7fce21dc6e763f42b99ef7fcf7416decaf848a5124cda99741222dba9619eb1b"
MODULE = "src/openai_ns_reconstruction/original_long_reference_finite_N_I1_numeric.py"
C5 = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
Q4 = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
ORDER = 3
I1_ORDER = 1
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
CHARTS = ("long_reshape", "reference")


def _load(root, path):
    return json.loads((root / path).read_bytes())


def _record_value(row, preferred):
    """Accept a direct directed record or one of the saved row wrappers."""
    if row is None:
        return None
    if isinstance(row, Mapping):
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
             if isinstance(row, Mapping)] == list(range(length)),
            "Unordered " + label + " rows")
    return rows


def _read_log_rows(c, rows, label, length=ORDER + 1):
    values = []
    for row in _ordered_rows(rows, label, length):
        if isinstance(row, Mapping) and row.get("exact_zero"):
            values.append(None)
            continue
        raw = _record_value(row, "log_absolute_upper_bound")
        require(raw is not None, "Missing logarithmic bound in " + label)
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
    """Make a directed point at the known upper endpoint of ``value``."""
    if hasattr(value, "_mpi_"):
        value = value.b
    return c.mpf([value, value])


def _upper_max(c, values):
    values = [value for value in values if value is not None]
    require(values, "At least one absolute log cap is required")
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
    products = []
    for coefficient, value in terms:
        product = _log_product(c, _log_abs_scalar(c, coefficient), value)
        if product is not None:
            products.append(product)
    return _logsum(c, products)


def _log_scale(c, value, scalar):
    return _log_product(c, _log_abs_scalar(c, scalar), value)


def _log_jet_product(c, left, right, order=I1_ORDER):
    """Ordinary-Z derivative product in logarithmic cap form."""
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
    """Multiply two signed ordinary rows by the Leibniz rule."""
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


def _mapping(value, label):
    require(isinstance(value, Mapping), label + " mapping is required")
    return value


def _source_pins(root, parent, *packets):
    """Verify the parent and both source-only packets without executing them."""
    pins = {
        PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256,
                     passed=True, source_executed=False),
    }

    def add_packet(packet, label):
        require(packet.get("passed") is True, label + " did not pass")
        require(packet.get("parent_report_sha256") == PARENT_SHA256,
                label + " is tied to a different parent")
        raw = packet.get("source_pins")
        require(isinstance(raw, Mapping) and raw, label + ".source_pins required")
        for key, value in raw.items():
            require(isinstance(value, Mapping), label + " source pin must be a mapping")
            path = str(value.get("path", key))
            expected = value.get("expected_sha256") or value.get("sha256")
            require(value.get("passed") is True and isinstance(expected, str)
                    and len(expected) == 64, "Invalid source pin: " + path)
            local = Path(path)
            if not local.is_absolute():
                local = root / local
            require(local.is_file(), "Pinned source is missing: " + path)
            require(sha(local) == expected, "Pinned source changed: " + path)
            previous = pins.get(path)
            if previous is not None:
                old = previous.get("expected_sha256")
                require(old == expected, "Conflicting source pin: " + path)
            pins[path] = dict(path=path, expected_sha256=expected,
                              passed=True, source_executed=False)

    for packet, label in zip(packets, ("Native control packet", "Loop packet")):
        _mapping(packet, label)
        add_packet(packet, label)

    # Include the accepted native-background pins, but do not re-run its
    # pressure identity or any of its accepted producers.
    native_bg = _mapping(parent.get("actual_long_reference_background_I1_numeric"),
                         "Parent native-background report")
    for key, value in _mapping(native_bg.get("source_pins"),
                               "Parent native-background source_pins").items():
        require(isinstance(value, Mapping), "Parent source pin must be a mapping")
        path = str(value.get("path", key))
        expected = value.get("expected_sha256") or value.get("sha256")
        require(value.get("passed") is True and isinstance(expected, str)
                and len(expected) == 64, "Invalid parent source pin: " + path)
        local = Path(path)
        if not local.is_absolute():
            local = root / local
        require(local.is_file() and sha(local) == expected,
                "Changed parent source: " + path)
        previous = pins.get(path)
        if previous is not None:
            require(previous["expected_sha256"] == expected,
                    "Parent/source packet pin conflict: " + path)
        pins[path] = dict(path=path, expected_sha256=expected,
                          passed=True, source_executed=False)

    pins[MODULE] = dict(path=MODULE, expected_sha256=sha(root / MODULE),
                        passed=True, source_executed=False)
    return pins


def _source_guards(root):
    """Guard only the finite-N, full-mean transport, and physical-Q sources."""
    return dict(
        finite_N_C5_V_N=source_guard(root, C5, "finite_N_functions", (
            "s=a.scale(loop['B_over_Pstar'],invN)",
            "V_N=a.add(V,s)",
        )),
        full_mean_transport=source_guard(root, C5, "transport", (
            "body=a.scale(a.at(densities[name],coordinate,t),g.mul(jacobian,kernel))",
            "predecessor_memory_never_reset=True",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
        )),
        physical_Q4_Q=source_guard(root, Q4, "physical_pressure_fields", (
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
    )


def _loop_packet_rows(c, loop_packet, name):
    require(loop_packet.get("passed") is True,
            "Finite-N loop packet is not accepted")
    require(loop_packet.get("parent_report_sha256") == PARENT_SHA256,
            "Loop packet is tied to a different parent")
    require(loop_packet.get("ordinary_Z_orders") == [0, 1, 2, 3],
            "Four ordinary Z rows are required")
    require(loop_packet.get("physical_y_derivatives_claimed") is False,
            "Physical y derivatives must remain unclaimed")
    contract = _mapping(loop_packet.get("source_contract"),
                        "Loop packet source_contract")
    require(contract.get("finite_N_numeric") is False and
            contract.get("seams_proved") is False,
            "Loop packet finite-N/seam claims changed")
    scale = _mapping(loop_packet.get("charts", {}).get(name, {}).get(
        "physical_scale_proof"), name + ".physical_scale_proof")
    require(scale.get("passed") is True and
            scale.get("B_quantity") == "B/Pstar/frozen_E0" and
            scale.get("deltaU_formula") ==
            "physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N" and
            scale.get("current_N") == "N18 graph handle" and
            scale.get("actual_N18") is True and
            scale.get("actual_N18_is_graph_handle") is True and
            scale.get("actual_N_not_materialized") is True and
            scale.get("phase_Z_independent") is True and
            scale.get("source_bound_to_current_native_long_reference_controls") is True and
            scale.get("physical_y_derivatives_claimed") is False,
            "Missing current long/reference physical loop units")
    chart = _mapping(loop_packet.get("charts", {}).get(name),
                     name + ".loop chart")
    B = _read_log_rows(
        c, chart.get("B_over_Pstar_over_frozen_E0_log_caps4rows",
                     chart.get("B_over_Pstar_over_frozen_E0_log_caps")),
        name + ".B_over_Pstar_over_frozen_E0", ORDER + 1)
    require(all(value is not None for value in B),
            "All loop B log caps are required")
    logPE = read(c, _record_value(chart.get("log_Pstar_E_upper"),
                                  "log_absolute_upper_bound"))
    logN = read(c, _record_value(chart.get("logN"),
                                 "log_absolute_upper_bound"))
    return B, logPE, logN


def _recover_delta_Q_logs(c, delta_U, delta_M, z, delta):
    """Recover absolute ``Q`` error logs through ordinary Z order two."""
    require(len(delta_U) >= 4 and len(delta_M) >= 4,
            "Q recovery requires physical U and mean rows Z0..Z3")
    h = c.mpf(1) - delta
    d = c.mpf(1) - z * z
    L = c.mpf(1) - delta * z * z
    require(L.a > 0, "Current long/reference L=1-delta*Z^2 is not positive")
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
    q1 = _log_add(c, n1,
                  _log_product(c, _log_abs_scalar(c, L1), q0)) + invL
    q2 = _log_add(c, n2,
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
        formula=("deltaQ=(2*Z*deltaU-(1-delta)*Z*deltaM-"
                 "(1-Z^2)*deltaM_Z)/(1-delta*Z^2)"),
    )


def _difference_logs(c, U, q, U_b, q_b):
    """Bound the six ordinary C1 differences from the signed background rows."""
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
    return dict(q=q[:2], q_Z=qz, q_squared=q_squared,
                U_q_Z=U_q_Z, U_squared=U_squared, U_U_Z=UU_Z)


def _pick(mapping, *names):
    for name in names:
        if name in mapping:
            return mapping[name]
    raise KeyError(names[0])


def _full_I1_and_D(c, background_chart, difference, radius, normalization):
    arithmetic = _mapping(background_chart.get("arithmetic"),
                           "native arithmetic")
    signed = _mapping(arithmetic.get("six_integrals"),
                      "native six_integrals")
    base = {name: _read_signed_rows(c, signed[name], "signed " + name,
                                    I1_ORDER + 1) for name in I1_NAMES}
    full_values = {name: [value + _symmetric(c, radius) for value in base[name]]
                   for name in I1_NAMES}
    coefficients = _mapping(arithmetic.get("coefficients"),
                            "native coefficients")
    coeff = {name: _read_signed_rows(c, coefficients[name],
                                     "coefficient " + name,
                                     I1_ORDER + 1)
             for name in ("minus_inv_L", "hZ_over_2L", "minus_half",
                          "d_over_L", "two_Zd_over_L2")}
    coeff["two_h_Z2_over_L2"] = _read_signed_rows(
        c, _pick(coefficients, "two_h_Z2_over_L2", "two_hZ2_over_L2"),
        "coefficient two_h_Z2_over_L2", I1_ORDER + 1)
    coefficient_products = {
        "q": _ordinary_product_values(coeff["minus_inv_L"], full_values["q"]),
        "q_Z": _ordinary_product_values(coeff["hZ_over_2L"], full_values["q_Z"]),
        "q_squared": _ordinary_product_values(coeff["minus_half"],
                                               full_values["q_squared"]),
        "U_q_Z": _ordinary_product_values(coeff["d_over_L"],
                                           full_values["U_q_Z"]),
        "U_squared": _ordinary_product_values(coeff["two_h_Z2_over_L2"],
                                               full_values["U_squared"]),
        "U_U_Z": _ordinary_product_values(coeff["two_Zd_over_L2"],
                                           full_values["U_U_Z"]),
    }
    D_full = []
    for k in range(I1_ORDER + 1):
        terms = [coefficient_products[name][k] for name in I1_NAMES]
        total = terms[0]
        for term in terms[1:]:
            total += term
        D_full.append(total)
    return (
        dict(
            passed=True,
            normalized_by=normalization,
            six_integrals={name: _pack_values(full_values[name])
                           for name in I1_NAMES},
            finite_N_absolute_error_log_bounds={
                name: _pack_logs(difference[name]) for name in I1_NAMES
            },
            signed_background_unchanged=True,
            symmetric_radius_added_once=record(radius),
            error_radius_is_not_a_selected_value=True,
            Q_source="Q_pointwise",
            integrated_Q_used=False,
        ),
        dict(
            passed=True,
            rows=_pack_values(D_full),
            D_rows=_pack_values(D_full),
            coefficientsC1rows={name: _pack_values(values)
                               for name, values in coeff.items()},
            computed_from_full_I1_once=True,
            separate_D_error_not_added=True,
            normalized_by=normalization,
            no_extra_D_error=True,
        ),
    )


def _boundary_difference_logs(c, densities, dQ, delta, z):
    L = 1 - delta * z * z
    invL = [1 / L, 2 * delta * z / (L * L)]
    two_Z2_L2 = [2 * z * z / (L * L),
                 4 * z / (L * L) + 8 * delta * z * z * z / (L * L * L)]
    linear = _log_jet_product(c, _absolute_logs(c, invL), dQ[:2])
    squared_U = _log_jet_product(c, _absolute_logs(c, two_Z2_L2),
                                 densities["U_squared"])
    deltaG = [_log_add(c, linear[k],
                       _log_scale(c, densities["q_squared"][k], c.mpf(1) / 2),
                       squared_U[k]) for k in range(2)]
    return deltaG, [None if value is None else value + c.ln(2)
                    for value in deltaG]


def _radius_ready(values, log_radius):
    return all(value is None or value.b < log_radius.a for value in values)


def _chart(c, name, native_chart, native_packet, loop_packet, incoming, radius,
           log_radius, delta, z, threshold):
    require(native_chart.get("passed") is True,
            "Native background chart did not pass: " + name)
    require(native_chart.get("source_executed") is False,
            "Native background chart source was executed: " + name)
    require(native_chart.get("normalization") in ("Wlong", "Wref"),
            "Native chart normalization is not own W: " + name)
    expected_normalization = "Wlong" if name == "long_reshape" else "Wref"
    require(native_chart.get("normalization") == expected_normalization and
            native_chart.get("W_formula") == "W=R_upper^2; radius not materialized",
            "Wrong long/reference native normalization: " + name)
    arithmetic = _mapping(native_chart.get("arithmetic"),
                           name + ".arithmetic")
    require("Q_pointwise" in arithmetic and "integrated_q" in arithmetic,
            "Pointwise and integrated Q records are required: " + name)
    U_b = _absolute_logs(c, _read_signed_rows(
        c, native_chart.get("U_constant_rows"), name + " U_constant", 4))
    q_b = _absolute_logs(c, _read_signed_rows(
        c, arithmetic["Q_pointwise"], name + " Q_pointwise", 3))
    mass = read(c, native_chart["normalized_moments"][0])
    require(mass.a > 0 and mass.b <= c.mpf("0.5").b,
            "Native J0 mass must be positive and at most one half")
    B, logPE, logN = _loop_packet_rows(c, loop_packet, name)
    native_controls_N = read(c, native_packet["charts"][name]["logN"])
    # The native-control and loop records must agree exactly at the directed
    # endpoints for the common current-N graph handle.
    require(native_controls_N.a == logN.a and native_controls_N.b == logN.b,
            "Current N changed between native and loop packets")

    dU = [logPE + value - logN for value in B]
    dM = [_upper_max(c, (incoming[k], dU[k])) for k in range(4)]
    dQ, qproof = _recover_delta_Q_logs(c, dU, dM, z, delta)
    densities = _difference_logs(c, dU, dQ, U_b, q_b)
    mass_log = c.ln(mass.b)
    defects = {
        key: [None if value is None else value + mass_log for value in rows]
        for key, rows in densities.items()
    }
    require(all(value is None or value.b < threshold.a
                for rows in defects.values() for value in rows),
            "Current-N normalized six-I1 error is not below 2^-1024: " + name)
    deltaG, b_errors = _boundary_difference_logs(c, densities, dQ, delta, z)
    require(all(value is None or value.b < threshold.a for value in b_errors),
            "Current-N B boundary error is not below 2^-1024: " + name)
    all_actual_logs = [value for rows in defects.values() for value in rows]
    all_actual_logs += list(b_errors)
    require(_radius_ready(all_actual_logs, log_radius),
            "Finite-N defect exceeds the one-time 2^-1024 radius: " + name)

    background = dict(arithmetic=arithmetic)
    full, D = _full_I1_and_D(c, background, defects, radius,
                             expected_normalization)
    B_signed = _read_signed_rows(
        c, _pick(arithmetic, "B_over_W2", "B_over_W"),
        name + " signed B", 2)
    B_full_rows = [value + _symmetric(c, radius) for value in B_signed]
    D_rows = _read_signed_rows(c, D["rows"], name + " composed D", 2)
    B_full = dict(
        passed=True,
        rows=_pack_values(B_full_rows),
        finite_N_absolute_error_log_bounds=_pack_logs(b_errors),
        difference_G_log_bounds=_pack_logs(deltaG),
        endpoint_pair_log_factor=record(c.ln(2)),
        difference_formula=("deltaG=deltaQ/L+delta(Q^2)/2-"
                            "2*Z^2*delta(U^2)/L^2"),
        endpoint_error_formula="abs(deltaB)/W<=2*sup(abs(deltaG)); W=R_upper^2",
        own_W_formula="W=R_upper^2",
        symmetric_radius_added_once=record(radius),
        error_radius_is_not_a_selected_value=True,
        signed_background_unchanged=True,
        normalized_by=expected_normalization,
        Q_source="Q_pointwise",
        integrated_Q_used=False,
    )
    S1_rows = [B_full_rows[k] + D_rows[k] for k in range(2)]
    S1_full = dict(
        passed=True,
        rows=_pack_values(S1_rows),
        formula="S1_full=B_full+D_full",
        computed_from_B_and_D_once=True,
        separate_S1_error_not_added=True,
        no_extra_radius=True,
        normalized_by=expected_normalization,
        physical_y_derivatives_claimed=False,
    )
    return dict(
        passed=True,
        chart=name,
        parent_report_sha256=PARENT_SHA256,
        ordinary_Z_orders=[0, 1],
        physical_y_derivatives_claimed=False,
        I1_full=full,
        D_full=D,
        B_full=B_full,
        S1_full=S1_full,
        physical_error_log_bounds=dict(
            delta_U=_pack_logs(dU), delta_M=_pack_logs(dM),
            delta_Q=_pack_logs(dQ), logN=record(logN),
            log_Pstar_E_upper=record(logPE),
        ),
        six_I1_difference_error_log_bounds={
            key: _pack_logs(value) for key, value in defects.items()
        },
        B_difference_error_log_bounds=_pack_logs(b_errors),
        signed_I1_radius=record(radius),
        I1_error_threshold=record(threshold),
        normalized_mass=record(mass),
        geometry=dict(
            normalization=expected_normalization,
            log_R_upper=native_chart["log_R_upper"],
            log_W=native_chart["log_W"],
            W_formula="W=R_upper^2; radius not materialized",
            normalized_moments=native_chart["normalized_moments"],
            own_weight_used=True,
            no_extra_R_weight=True,
        ),
        source_contract=dict(
            passed=True,
            chart=name,
            common_N_handle=18,
            actual_current_N18=True,
            actual_current_N18_is_graph_handle=True,
            N_not_materialized=True,
            phase_Z_independent=True,
            physical_Pstar_applied_once=True,
            physical_R_not_a_mean_scale=True,
            physical_delta_U_local_only=True,
            incoming_delta_U_not_added=True,
            incoming_delta_M_retained=True,
            nonreset_positive_kernel_mean=True,
            mean_bound="max(incoming deltaM, local deltaU); positive kernel mass<=1",
            Q_source="native arithmetic.Q_pointwise",
            integrated_Q_used=False,
            own_W="R_upper^2",
            normalized_mass="J0 upper endpoint, <=1/2",
            true_physical_mass_once=True,
            finite_N_radius_added_once=True,
            signed_errors_are_caps_not_selected_values=True,
            physical_y_derivatives_claimed=False,
        ),
        limits=dict(
            six_I1_finite_N_C1_enclosed=True,
            B_and_S1_finite_N_C1_enclosed=True,
            D_is_not_full_S1=True,
            RQy_endpoints_uncomputed=True,
            outer_RQy_retained=True,
            outer_RQy_evaluated=False,
            seams_proved=False,
            B_functional_numeric=True,
            standalone_RQy_endpoint_numeric=False,
            complete_pressure_global_closure=False,
            full_pressure_global_closure=False,
            runtime_model_metadata_available=False,
            physical_y_derivatives_claimed=False,
            global_five_moment_closure=False,
        ),
    )


def build(root, parent, native_packet, loop_packet):
    """Build current-``N`` long-reshape/reference six-I1/B/D/S1 enclosures."""
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Current native-background parent bytes changed")
    disk_parent = _load(root, PARENT)
    require(parent == disk_parent, "Exact accepted native-background parent required")
    native_packet = _mapping(native_packet, "Native control packet")
    loop_packet = _mapping(loop_packet, "Loop packet")
    require(native_packet.get("parent_report_sha256") == PARENT_SHA256 and
            loop_packet.get("parent_report_sha256") == PARENT_SHA256,
            "New packets must belong to the native-background parent")
    require(native_packet.get("passed") is True and
            native_packet.get("physical_y_derivatives_claimed") is False,
            "Native control packet contract changed")
    native_contract = _mapping(native_packet.get("source_contract"),
                               "Native control source_contract")
    require(native_contract.get("finite_N_numeric") is False and
            native_contract.get("seams_proved") is False and
            native_contract.get("actual_native_long_reference_controls") is True,
            "Native controls must remain source-only")
    require(loop_packet.get("physical_y_derivatives_claimed") is False,
            "Loop packet physical y claim changed")

    pins = _source_pins(root, disk_parent, native_packet, loop_packet)
    guards = _source_guards(root)
    native_bg = disk_parent["actual_long_reference_background_I1_numeric"]
    lemma = _mapping(native_bg.get("inherited_lemma_reference"),
                     "Inherited pressure lemma")
    require(lemma.get("passed") is True and lemma.get("reexecuted") is False,
            "Only the inherited pressure identity may be reused")
    switch_post = _mapping(disk_parent.get("actual_switch_post_finite_N_I1_numeric"),
                          "Current post finite-N parent")
    require(switch_post.get("passed") is True,
            "Accepted post finite-N parent is required")
    post_chart = _mapping(switch_post.get("charts", {}).get("post_switch"),
                          "Accepted post chart")

    c = MPIntervalContext()
    c.prec = 1200
    radius = c.mpf(2) ** -1024
    log_radius = c.ln(radius)
    threshold = -1024 * c.ln(c.mpf(2))
    incoming = _read_log_rows(
        c, post_chart["physical_error_log_bounds"]["delta_M"],
        "post-switch incoming deltaM", 4)
    direction = _mapping(native_bg.get("directions"), "Native directions")
    delta = read(c, direction["delta"])
    z_rows = _read_signed_rows(c, direction["Z"], "source Z direction", 3)
    z = z_rows[0]
    charts = {}

    for name in CHARTS:
        charts[name] = _chart(
            c, name, native_bg["charts"][name], native_packet,
            loop_packet, incoming, radius, log_radius, delta, z, threshold)
        incoming = _read_log_rows(
            c, charts[name]["physical_error_log_bounds"]["delta_M"],
            name + " outgoing deltaM", 4)

    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_pins=pins,
        source_guards=guards,
        charts=charts,
        inherited_lemma_reference=lemma,
        pressure_identity_reexecuted=False,
        physical_y_derivatives_claimed=False,
        actual_current_N18=True,
        actual_current_N18_is_graph_handle=True,
        current_N_materialized=False,
        all_normalized_I1_errors_below_2_power_minus1024=True,
        source_contract=dict(
            passed=True,
            common_N_handle=18,
            actual_current_N18=True,
            actual_current_N18_is_graph_handle=True,
            N_materialized=False,
            actual_source_chart_mapping={
                "long_reshape": "long_reshape", "reference": "reference"
            },
            post_to_long_to_reference_mean_memory_retained=True,
            delta_U_has_no_incoming=True,
            physical_Pstar_applied_once=True,
            own_W_by_chart={"long_reshape": "R_upper^2", "reference": "R_upper^2"},
            normalized_mass_is_J0_upper_at_most_one_half=True,
            whole_Q_pointwise=True,
            integrated_Q_not_used=True,
            finite_N_radius_added_once=True,
            no_extra_D_error=True,
            no_extra_S1_radius=True,
            RQy_endpoints_uncomputed=True,
            seams_proved=False,
            full_pressure_global_closure=False,
            physical_y_derivatives_claimed=False,
            signed_errors_are_caps_not_selected_values=True,
        ),
        limits=dict(
            RQy_endpoints_uncomputed=True,
            outer_RQy_retained=True,
            outer_RQy_evaluated=False,
            seams_proved=False,
            chart_B_and_S1_numeric=True,
            D_is_not_full_S1=True,
            complete_pressure_global_closure=False,
            full_pressure_global_closure=False,
            runtime_model_metadata_available=False,
            global_five_moment_closure=False,
            accepted_pressure_identity_reexecuted=False,
            accepted_ancestors_replayed=False,
            physical_y_derivatives_claimed=False,
        ),
        accepted_ancestors_replayed=False,
        full_goal_complete=False,
    )
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    """Exercise only this module's pure logarithmic and Leibniz arithmetic."""
    c = MPIntervalContext()
    c.prec = 240
    zero = c.mpf(0)
    one = c.mpf(1)
    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    delta = c.mpf([0, c.mpf(1) / 1000])
    dU = [c.mpf(-20), c.mpf(-19), c.mpf(-18), c.mpf(-17)]
    dM = [c.mpf(-21), c.mpf(-20), c.mpf(-19), c.mpf(-18)]
    dQ, proof = _recover_delta_Q_logs(c, dU, dM, z, delta)
    require(proof["passed"] and len(dQ) == 3, "Fresh Q recovery fixture failed")
    U = [c.mpf(-20), c.mpf(-19), c.mpf(-18)]
    q = [c.mpf(-20), c.mpf(-19), c.mpf(-18)]
    U_b = [c.mpf(2), c.mpf(0), c.mpf(0)]
    q_b = [c.mpf(3), c.mpf(0), c.mpf(0)]
    densities = _difference_logs(c, U, q, U_b, q_b)
    require(set(densities) == set(I1_NAMES) and
            all(len(row) == 2 for row in densities.values()),
            "Fresh six-I1 difference fixture failed")
    deltaG, boundary = _boundary_difference_logs(c, densities, dQ, delta, z)
    require(len(deltaG) == len(boundary) == 2 and
            all(value is not None for value in boundary),
            "Fresh B boundary fixture failed")
    radius = c.mpf(2) ** -1024
    require(_radius_ready([c.ln(radius) - 1, None], c.ln(radius)),
            "Fresh one-time radius fixture failed")
    return dict(
        passed=True,
        Q_recovery=True,
        six_difference_rows=True,
        boundary_pair_formula=True,
        own_W_is_R_upper_squared=True,
        one_time_radius_gate=True,
        source_executed=False,
    )


def self_check():
    return fixture_bound_check()


__all__ = (
    "PARENT", "PARENT_SHA256", "MODULE", "CHARTS", "I1_NAMES",
    "build", "fixture_bound_check", "self_check",
)
