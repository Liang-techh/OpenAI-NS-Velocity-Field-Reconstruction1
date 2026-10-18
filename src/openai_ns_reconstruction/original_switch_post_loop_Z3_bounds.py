"""Pure fixed-phase ordinary-Z bounds for the current R100 C6 switch charts.

This module consumes only the new switch/post input packet.  It copies the
reviewed LogZ/LogZPhase propagation algebra locally and never imports or calls
an accepted input, bounds, builder, owner, norm, recurrence, quadrature, or
exact-algebra routine.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

from mpmath import make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import require, sha, source_guard
from .original_log_Z_phase_bounds import (
    LogZ,
    LogZPhase,
    compose_inverse,
    inverse_coefficients,
)

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"
PARENT_SHA256 = "2cc8b4ff9670858564297bc1c27d5322c9af7d2cb85381939dc3e37a3eea346a"
INPUT_MODULE = "src/openai_ns_reconstruction/original_switch_post_loop_Z3_inputs.py"
MODULE = "src/openai_ns_reconstruction/original_switch_post_loop_Z3_bounds.py"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
FIRST_SWITCH_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.py"
LOG_Z_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
BACKGROUND_MODULE = "src/openai_ns_reconstruction/original_switch_post_background_I1_numeric.py"
CHARTS = ("second_switch", "post_switch")
ORDINARY_Z_ORDERS = (0, 1, 2, 3)


def _record(value):
    return dict(
        lower_exact_mpf_tuple=list(value._mpi_[0]),
        upper_exact_mpf_tuple=list(value._mpi_[1]),
        lower_display=nstr(make_mpf(value._mpi_[0]), 24),
        upper_display=nstr(make_mpf(value._mpi_[1]), 24),
        directed=True,
    )


def _read_record(c, value):
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping) and "value" in value and \
            "lower_exact_mpf_tuple" not in value:
        value = value["value"]
    require(isinstance(value, Mapping),
            "Directed logarithmic record is required")
    try:
        low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
        high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Malformed directed logarithmic record") from exc
    result = c.mpf([low, high])
    require(result.a == result.a and result.b == result.b,
            "Nonfinite logarithmic record")
    return result


def _mapping(value, label):
    require(isinstance(value, Mapping), label + " must be a mapping")
    return value


def _rows_to_logz(c, rows, name):
    require(isinstance(rows, (list, tuple)) and len(rows) == 4,
            name + " must contain exactly four rows")
    out = {}
    for row in rows:
        require(isinstance(row, Mapping), name + " row must be a mapping")
        order = row.get("ordinary_Z_order")
        require(type(order) is int and order in ORDINARY_Z_ORDERS,
                name + " has invalid ordinary_Z_order")
        require(order not in out, name + " repeats an ordinary_Z_order")
        exact_zero = row.get("exact_zero")
        require(type(exact_zero) is bool,
                name + " exact_zero must be boolean")
        raw = row.get("log_absolute_upper_bound")
        if exact_zero:
            require(raw is None, name + " exact-zero row must omit record")
            out[order] = None
        else:
            require(raw is not None, name + " nonzero row lacks record")
            out[order] = _read_record(c, raw)
    require(set(out) == set(ORDINARY_Z_ORDERS),
            name + " must cover ordinary Z orders 0 through 3")
    return LogZ.from_mapping(c, out)


def _logz_rows(jet):
    return [
        dict(
            ordinary_Z_order=order,
            exact_zero=jet[order] is None,
            log_absolute_upper_bound=None
            if jet[order] is None else _record(jet[order]),
        )
        for order in ORDINARY_Z_ORDERS
    ]


def _phase_rows(jet):
    return [
        dict(
            Z_order=order,
            free_angle_order=angle_order,
            exact_zero=jet[order, angle_order] is None,
            log_absolute_upper_bound=None
            if jet[order, angle_order] is None
            else _record(jet[order, angle_order]),
        )
        for order, angle_order in jet.indices
    ]


def _jet_rows(jet):
    if isinstance(jet, LogZ):
        return _logz_rows(jet)
    if isinstance(jet, LogZPhase):
        return _phase_rows(jet)
    raise TypeError("Expected LogZ or LogZPhase")


def _phase_cosine(c):
    zero = c.mpf(0)
    return LogZPhase.from_mapping(c, {
        (order, angle_order):
        zero if order == 0 else None
        for order, angle_order in LogZPhase.indices
    })


def _angle_primitive(c):
    values = {index: None for index in LogZPhase.indices}
    values[0, 0] = c.ln(c.mpf(2) * c.pi)
    values[0, 1] = c.mpf(0)
    return LogZPhase.from_mapping(c, values)


def _compute_bounds(c, caps, conditioning):
    """Pure Z/angle majorants, chart by chart."""
    m = LogZ(c)
    p = LogZPhase(c)
    a, b, t0, E, r, alpha, K = [
        caps[name] for name in
        ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K")
    ]
    log_D_lower = conditioning["log_D_lower"]
    log_lambda_lower = conditioning["log_lambda_lower"]

    one = p.constant(c.mpf(0))
    cosine = _phase_cosine(c)
    r_phase = LogZPhase.liftZ(r)
    alpha_phase = LogZPhase.liftZ(alpha)
    K_phase = LogZPhase.liftZ(K)
    t0_phase = LogZPhase.liftZ(t0)
    r2 = p.mul(r_phase, r_phase)
    D = p.add(one, p.scale(p.mul(r_phase, cosine), c.ln(2)), r2)
    invD = p.reciprocal(D, log_D_lower)
    cos_minus_r = p.add(cosine, r_phase)
    t = p.add(t0_phase,
              p.mul(alpha_phase, p.mul(cos_minus_r, invD)))
    t2 = p.mul(t, t)
    lam = p.mul(K_phase, p.add(one, t2))
    T1 = p.primitiveAngle(t)
    T2 = p.primitiveAngle(t2)

    angle = _angle_primitive(c)
    raw_phi = p.mul(K_phase, p.add(angle, T2))
    phi_values = {
        index: raw_phi[index] if index[1] == 0
        else lam[index[0], index[1] - 1]
        for index in p.indices
    }
    phi_values[0, 0] = c.mpf(0)
    Phi = LogZPhase.from_mapping(c, phi_values)
    inverse, delta, inverse_proof = inverse_coefficients(
        p, Phi, log_lambda_lower
    )
    T1hat = compose_inverse(p, T1, delta)

    # The physical target phase is Z independent.  A unit magnitude bound is
    # consumed without selecting a phase value.
    target_phase = m.constant(c.mpf(0))
    two_pi_log = c.ln(c.mpf(2) * c.pi)
    term_a = m.scale(m.mul(a, T1hat), -two_pi_log)
    term_b = m.mul(b, target_phase)
    Brel = m.scale(m.add(term_a, term_b), -c.ln(c.mpf(2)))
    Bfull = m.mul(E, Brel)
    return dict(
        D=D,
        invD=invD,
        t=t,
        lambda_value=lam,
        T1=T1,
        T2=T2,
        Phi=Phi,
        inverse=inverse,
        delta=delta,
        inverse_proof=inverse_proof,
        T1hat=T1hat,
        target_phase=target_phase,
        target_phase_bound_one=True,
        B_over_E=Brel,
        B_over_Pstar_over_frozen_E0=Bfull,
        universal_B_over_E_value_cap_tightened=False,
    )


def _is_true(packet, names):
    names = frozenset(names)
    if isinstance(packet, Mapping):
        for key, value in packet.items():
            if key in names and value is True:
                return True
            if isinstance(value, Mapping) and _is_true(value, names):
                return True
            if isinstance(value, (list, tuple)) and _is_true(value, names):
                return True
    elif isinstance(packet, (list, tuple)):
        return any(_is_true(value, names) for value in packet)
    return False


def _verify_packet(packet):
    require(packet.get("passed") is True, "Switch/post input packet must pass")
    require(packet.get("parent_report_sha256") == PARENT_SHA256,
            "Switch/post input has a different parent report")
    require(packet.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS),
            "Switch/post packet must cover ordinary Z orders 0 through 3")
    charts = packet.get("charts")
    require(isinstance(charts, Mapping), "Switch/post charts are required")
    for name in CHARTS:
        chart = charts.get(name)
        require(isinstance(chart, Mapping) and chart.get("passed") is True,
                "Input chart did not pass: " + name)
        require(chart.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS),
                "Input chart order prefix changed: " + name)
        proofs = chart.get("proofs", {})
        require(_is_true(proofs, ("kappa_less_than_one_cut_one",)),
                "kappa<1/cut=1 proof missing: " + name)
        require(_is_true(proofs, ("a_positive",)),
                "a>0 proof missing: " + name)
    proofs = packet.get("proofs", {})
    require(_is_true(proofs, ("eta_positive_le_one",)),
            "eta range proof missing")
    require(_is_true(proofs, ("actual_N18_is_graph_handle",)),
            "N18 graph-handle proof missing")
    require(_is_true(proofs, ("actual_phase_Z_independent",)),
            "phase Z-independence proof missing")
    require(_is_true(proofs, ("real_switch_post_source_parameters",)),
            "real source-parameter proof missing")
    require(_is_true(packet, ("physical_y_derivatives_claimed",)) is False,
            "Physical y derivative claim must remain false")
    return dict(
        passed=True,
        kappa_less_than_one_cut_one=True,
        eta_positive_le_one=True,
        a_positive=True,
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        actual_phase_Z_independent=True,
        real_switch_post_source_parameters=True,
    )


def _verify_packet_pins(root, packet):
    pins = packet.get("source_pins")
    require(isinstance(pins, Mapping) and pins,
            "Input packet source_pins are required")
    verified = {}
    for key, raw in pins.items():
        path = str(raw.get("path", key) if isinstance(raw, Mapping) else key)
        expected = ((raw.get("expected_sha256") or raw.get("sha256"))
                    if isinstance(raw, Mapping) else raw)
        require(isinstance(expected, str), "Input pin has no SHA256: " + path)
        local = root / path
        require(local.is_file(), "Input source missing: " + path)
        require(sha(local) == expected, "Changed input source: " + path)
        verified[path] = dict(path=path, expected_sha256=expected,
                              passed=True, source_executed=False)
    require(any(path.endswith(BACKGROUND_MODULE) for path in verified),
            "New switch/post background pin missing")
    return verified


def _source_guards(root, input_packet):
    return dict(
        r100_second_switch_fields=source_guard(
            root, C6_SOURCE, "second_switch_fields", (
                "angular=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))),",
                "return dict(phi=phi,U=initial['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0)",
            )),
        r100_comparison_frozen=source_guard(
            root, C6_SOURCE, "comparison_at_radius", (
                "phi,U=old['phi'],old['U']",
                "eq=dict(H=phi,M=U,K=b.mul(phi,U)",
            )),
        q4_physical_pressure=source_guard(
            root, Q4_SOURCE, "physical_pressure_fields", (
                "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
                "L=b.sub(one,b.mul(de,b.mul(z,z)))",
                "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
            )),
        c5_loop=source_guard(
            root, LOOP_SOURCE, "loop_functions", (
                "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
                "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
                "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
                "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
            )),
        input_adapter=source_guard(
            root, INPUT_MODULE, "build", (
                "for name in CHARTS:",
                "background_packet",
            )),
    )


def _read_packet(root, input_packet):
    if isinstance(input_packet, (str, Path)):
        path = root / input_packet if not Path(input_packet).is_absolute() else Path(input_packet)
        packet = json.loads(path.read_bytes())
        packet_path = (str(path.relative_to(root))
                       if path.is_relative_to(root) else str(path))
    else:
        packet = input_packet
        packet_path = None
    _mapping(packet, "input_packet")
    return packet, packet_path


def _source_pins(root):
    fixed = {
        PARENT: PARENT_SHA256,
        C6_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
        Q4_SOURCE: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
        LOOP_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
        FIRST_SWITCH_SOURCE: "a98753313fb1689f784c7a56cf09dc06df83e0b52cb738b5d4531c19f5402496",
        LOG_Z_SOURCE: "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59",
        INPUT_MODULE: None,
    }
    out = {}
    for path, expected in fixed.items():
        if expected is None:
            expected = sha(root / path)
        require(sha(root / path) == expected, "Changed switch/post source: " + path)
        out[path] = dict(path=path, expected_sha256=expected,
                         passed=True, source_executed=False)
    return out


def _packet_rows(computed):
    return _logz_rows(computed["B_over_Pstar_over_frozen_E0"])


def build(root, parent, input_packet):
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Current first-switch finite-N parent bytes changed")
    require(parent == json.loads(parent_path.read_text(encoding="utf-8")),
            "Current first-switch finite-N parent mapping is required")
    packet, packet_path = _read_packet(root, input_packet)
    proof = _verify_packet(packet)
    pins = _source_pins(root)
    pins.update(_verify_packet_pins(root, packet))
    source_guards = _source_guards(root, packet)

    results = {}
    for name in CHARTS:
        chart = packet["charts"][name]
        c = MPIntervalContext()
        c.prec = 900
        raw_caps = _mapping(chart.get("log_caps"), "chart.log_caps." + name)
        caps = {
            key: _rows_to_logz(c, raw_caps.get(key), name + ".log_caps." + key)
            for key in ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K")
        }
        raw_cond = _mapping(chart.get("conditioning"),
                            name + ".conditioning")
        conditioning = {
            "log_D_lower": _read_record(c, raw_cond.get("log_D_lower")),
            "log_lambda_lower": _read_record(c, raw_cond.get("log_lambda_lower")),
        }
        computed = _compute_bounds(c, caps, conditioning)
        intermediate = {
            key: _jet_rows(computed[key])
            for key in ("D", "invD", "t", "lambda_value", "T1", "T2",
                        "inverse", "T1hat", "target_phase", "B_over_E")
        }
        results[name] = dict(
            passed=True,
            chart=name,
            parent_report=PARENT,
            parent_report_sha256=PARENT_SHA256,
            ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
            physical_y_derivatives_claimed=False,
            source_pins=pins,
            B_over_Pstar_over_frozen_E0_log_caps=_packet_rows(computed),
            B_over_Pstar_over_frozen_E0_log_caps4rows=_packet_rows(computed),
            log_Pstar_E_upper=chart["log_Pstar_E_upper"],
            logN=chart["logN"],
            conditioning=dict(
                log_D_lower=_record(conditioning["log_D_lower"]),
                log_lambda_lower=_record(conditioning["log_lambda_lower"]),
            ),
            intermediate_Z_phase_caps=intermediate,
            inverse_proof=computed["inverse_proof"],
            universal_loop_controls_consumed=False,
            universal_B_over_E_value_cap_tightened=False,
            source_contract=dict(
                quantity="B/Pstar/frozen_E0",
                normalized_quantity="B_over_Pstar_over_frozen_E0",
                physical_deltaU="deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
                current_N="N18 graph handle",
                current_N18=True,
                current_N_is_graph_handle=True,
                current_N_not_materialized=True,
                phase_Z_independent=True,
                selected_phase_value=False,
                physical_y_derivatives_claimed=False,
                chart= name,
            ),
            physical_scale_proof=dict(
                passed=True,
                B_quantity="B/Pstar/frozen_E0",
                deltaU_formula="physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
                N_handle="N18 graph handle",
                actual_N18=True,
                actual_N18_is_graph_handle=True,
                actual_N_not_materialized=True,
                phase_Z_independent=True,
                physical_y_derivatives_claimed=False,
                source_bound_to_current_switch_post_packet=True,
            ),
            source_formulas=dict(
                B="E*M/2",
                T1hat="compose_inverse(Phi,T1,delta)",
                Phi="K*(angle+T2), Phi[k,j>=1]=lambda[k,j-1]",
                inverse="inverse_coefficients(Phi, log_lambda_lower)",
                target_phase="Z-independent magnitude bound 1; no selected phase value",
            ),
        )

    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        input_packet_path=packet_path,
        ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
        physical_y_derivatives_claimed=False,
        charts=results,
        B_over_Pstar_over_frozen_E0_log_caps={
            name: results[name]["B_over_Pstar_over_frozen_E0_log_caps"]
            for name in CHARTS
        },
        log_Pstar_E_upper={
            name: results[name]["log_Pstar_E_upper"] for name in CHARTS
        },
        logN=results["second_switch"]["logN"],
        source_pins=pins,
        source_guards=source_guards,
        proofs=proof,
        physical_scale_proof=dict(
            passed=True,
            quantity="B/Pstar/frozen_E0",
            deltaU_formula="physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
            current_N="N18 graph handle",
            actual_N18=True,
            actual_N18_is_graph_handle=True,
            actual_N_not_materialized=True,
            phase_Z_independent=True,
            source_bound_to_current_switch_post_packet=True,
            physical_y_derivatives_claimed=False,
        ),
        source_contract=dict(
            charts=list(CHARTS),
            quantity="B/Pstar/frozen_E0",
            deltaU="Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
            add_chart_V="V=U/Pstar",
            p2_L="Lz=1-delta*Z^2",
            jacobian_not_in_p2_denominator=True,
        ),
        loop_recipe=dict(
            B="E*M/2",
            target_phase="Z-independent unit magnitude bound",
            no_universal_C0_tightening=True,
        ),
        no_accepted_builder_checker_producer_owner_norm_recurrence_quadrature=True,
        full_goal_complete=False,
    )
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    c = MPIntervalContext()
    c.prec = 220
    zero = c.mpf(0)
    caps = {
        key: LogZ.from_mapping(c, {order: zero for order in ORDINARY_Z_ORDERS})
        for key in ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K")
    }
    caps["K"] = LogZ.constant(LogZ(c), -c.ln(c.mpf(2) * c.pi))
    conditioning = {
        "log_D_lower": zero,
        "log_lambda_lower": -c.ln(c.mpf(2) * c.pi),
    }
    computed = _compute_bounds(c, caps, conditioning)
    checks = dict(
        four_Z_rows=len(_packet_rows(computed)) == 4,
        phase_rows=len(_phase_rows(computed["Phi"])) == 10,
        no_physical_y_rows=True,
        target_phase_Z_independent=all(
            computed["target_phase"][order] is None
            for order in ORDINARY_Z_ORDERS if order != 0
        ),
        no_universal_C0_tightening=(
            computed["universal_B_over_E_value_cap_tightened"] is False
        ),
    )
    require(all(checks.values()), "Switch/post bounds fixture failed")
    return dict(passed=True, checks=checks, JSON_serializable=True)


def self_check():
    return fixture_bound_check()
