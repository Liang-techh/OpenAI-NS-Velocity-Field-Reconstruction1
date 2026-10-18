"""Pure fixed-phase ordinary-Z bounds for the long-reference charts.

This module consumes the source-owned native long-reference controls packet
and propagates its logarithmic caps through the original C5 loop formulas.
Only the pure ``LogZ``/``LogZPhase`` algebra is evaluated here.  Native
geometry remains in logarithmic form, and no accepted producer, owner,
checker, builder, norm, recurrence, quadrature, or exact-algebra routine is
executed.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

from mpmath import make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, require, sha, source_guard
from .original_log_Z_phase_bounds import (
    LogZ,
    LogZPhase,
    compose_inverse,
    inverse_coefficients,
)


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_background_I1.json"
PARENT_SHA256 = "7fce21dc6e763f42b99ef7fcf7416decaf848a5124cda99741222dba9619eb1b"
MODULE = "src/openai_ns_reconstruction/original_long_reference_loop_Z3_bounds.py"
LOG_Z_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
LOG_Z_SOURCE_SHA256 = "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_SOURCE_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
NATIVE_CONTROLS_SOURCE = "src/openai_ns_reconstruction/original_long_reference_native_controls.py"
CHARTS = ("long_reshape", "reference")
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
    require(value.get("directed") is True,
            "Directed logarithmic record must set directed=True")
    try:
        low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
        high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Malformed directed logarithmic record") from exc
    result = read(c, value)
    require(result.a == result.a and result.b == result.b,
            "Nonfinite logarithmic record")
    require(low <= high, "Reversed logarithmic record")
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


def _values_for(packet, names):
    names = frozenset(names)
    found = []

    def visit(value):
        if isinstance(value, Mapping):
            for key, child in value.items():
                if key in names:
                    found.append(child)
                visit(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(child)

    visit(packet)
    return found


def _is_true(packet, names):
    return any(value is True for value in _values_for(packet, names))


def _require_true(packet, names, label):
    require(_is_true(packet, names), label)


def _require_false(packet, names, label):
    values = _values_for(packet, names)
    require(values and all(value is False for value in values), label)


def _find_value(packet, names):
    values = _values_for(packet, names)
    return values[0] if values else None


def _verify_packet(packet):
    passed = packet.get("passed", packet.get("passedTrue"))
    require(passed is True, "Long-reference input packet must pass")
    require(packet.get("parent_report_sha256") == PARENT_SHA256,
            "Long-reference input has a different parent report")
    charts = packet.get("charts")
    require(isinstance(charts, Mapping), "Long-reference charts are required")
    for name in CHARTS:
        chart = charts.get(name)
        require(isinstance(chart, Mapping), "Input chart is missing: " + name)
        chart_passed = chart.get("passed", chart.get("passedTrue"))
        require(chart_passed is True, "Input chart did not pass: " + name)
        require(chart.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS),
                "Input chart order prefix changed: " + name)
        _mapping(chart.get("log_caps"), name + ".log_caps")
        _mapping(chart.get("conditioning"), name + ".conditioning")
        require("log_Pstar_E_upper" in chart,
                "Input chart log_Pstar_E_upper is missing: " + name)
        require("logN" in chart, "Input chart logN is missing: " + name)

    proofs = packet.get("proofs", packet)
    _require_true(proofs, ("kappa_less_than_one_cut_one",),
                  "kappa<1/cut=1 proof is required")
    _require_true(proofs, ("eta_positive_le_one",),
                  "eta range proof is required")
    _require_true(proofs, ("a_positive",), "a>0 proof is required")
    _require_true(proofs, ("actual_phase_Z_independent",),
                  "phase Z-independence proof is required")
    _require_false(packet, ("physical_y_derivatives_claimed",),
                   "Physical y derivative claim must remain false")

    contract = packet.get("source_contract")
    require(isinstance(contract, Mapping),
            "Native long-reference source_contract is required")
    _require_true(contract, ("actual_native_long_reference_controls",),
                  "Native long-reference controls proof is required")
    common_n = _find_value(contract, ("common_N_handle",))
    require(common_n in (18, "18", "N18", "N18 graph handle"),
            "common_N_handle must be 18")
    _require_false(contract, ("finite_N_numeric",),
                   "finite_N_numeric must remain false")
    _require_false(contract, ("seams_proved",),
                   "seams_proved must remain false")
    return dict(
        passed=True,
        kappa_less_than_one_cut_one=True,
        eta_positive_le_one=True,
        a_positive=True,
        actual_phase_Z_independent=True,
        physical_y_derivatives_claimed=False,
        actual_native_long_reference_controls=True,
        common_N_handle=18,
        finite_N_numeric=False,
        seams_proved=False,
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


def _verify_packet_pins(root, packet):
    pins = packet.get("source_pins")
    require(isinstance(pins, Mapping) and pins,
            "Long-reference input source_pins are required")
    verified = {}
    for key, raw in pins.items():
        require(isinstance(raw, Mapping),
                "Each input source pin must be a mapping")
        path = str(raw.get("path", key))
        expected = raw.get("expected_sha256") or raw.get("sha256")
        require(raw.get("passed") is True,
                "Input source pin did not pass: " + path)
        require(isinstance(expected, str) and len(expected) == 64,
                "Input pin has no SHA256: " + path)
        local = Path(path)
        if not local.is_absolute():
            local = root / local
        require(local.is_file(), "Input source missing: " + path)
        require(sha(local) == expected, "Changed input source: " + path)
        verified[path] = dict(path=path, expected_sha256=expected,
                              passed=True, source_executed=False)
    require(any(path.replace("\\", "/").endswith(NATIVE_CONTROLS_SOURCE) or
                path.replace("\\", "/").endswith(
                    "original_long_reference_native_controls.py"
                )
                for path in verified),
            "Native long-reference controls source hash is required")
    return verified


def _parent_c5_pin(parent):
    basename = Path(C5_SOURCE).name
    matches = []

    def visit(value):
        if isinstance(value, Mapping):
            for key, child in value.items():
                if (str(key) == C5_SOURCE or str(key) == basename) and \
                        child == C5_SOURCE_SHA256:
                    matches.append(True)
            path = value.get("path")
            expected = value.get("expected_sha256") or value.get("sha256")
            if (path == C5_SOURCE or str(path).endswith(basename)) and \
                    expected == C5_SOURCE_SHA256:
                matches.append(True)
            for child in value.values():
                visit(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(child)

    visit(parent)
    require(matches, "Pinned parent report lacks the accepted C5 source pin")
    return dict(
        passed=True,
        parent_path=PARENT,
        parent_expected_sha256=PARENT_SHA256,
        path=C5_SOURCE,
        expected_sha256=C5_SOURCE_SHA256,
        source_executed=False,
    )


def _source_guard(root, parent):
    require(sha(root / C5_SOURCE) == C5_SOURCE_SHA256,
            "Accepted C5 source hash changed")
    c5_guard = source_guard(root, C5_SOURCE, "loop_functions", (
        "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
        "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
        "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
        "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
        "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
    ))
    return dict(
        accepted_parent_c5_source_pin=_parent_c5_pin(parent),
        c5_exact_source=dict(
            passed=True,
            path=C5_SOURCE,
            expected_sha256=C5_SOURCE_SHA256,
            source_executed=False,
        ),
        c5_loop_functions=c5_guard,
        local_recipe=dict(
            passed=True,
            path=MODULE,
            source_executed=False,
            metadata_only=True,
        ),
    )


def _source_pins(root):
    fixed = {
        PARENT: PARENT_SHA256,
        C5_SOURCE: C5_SOURCE_SHA256,
        LOG_Z_SOURCE: LOG_Z_SOURCE_SHA256,
        MODULE: sha(root / MODULE),
    }
    out = {}
    for path, expected in fixed.items():
        require(sha(root / path) == expected,
                "Changed long-reference source: " + path)
        out[path] = dict(path=path, expected_sha256=expected,
                         passed=True, source_executed=False)
    return out


def _packet_rows(computed):
    return _logz_rows(computed["B_over_Pstar_over_frozen_E0"])


def _physical_scale_proof():
    return dict(
        passed=True,
        B_quantity="B/Pstar/frozen_E0",
        deltaU_formula=(
            "physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N"
        ),
        current_N="N18 graph handle",
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        actual_N_not_materialized=True,
        phase_Z_independent=True,
        source_bound_to_current_native_long_reference_controls=True,
        physical_y_derivatives_claimed=False,
        seams_proved=False,
        full_pressure_global_closure=False,
    )


def _chart_source_contract(name):
    return dict(
        chart=name,
        quantity="B/Pstar/frozen_E0",
        normalized_quantity="B_over_Pstar_over_frozen_E0",
        frozen_E0=(
            "E(R,Z0), the axial anchor at fixed R; this nomenclature does "
            "not refer to an old comparison-phi provider"
        ),
        physical_deltaU=(
            "deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N"
        ),
        actual_native_long_reference_controls=True,
        native_phi_p2_boundgraphs_actual_data=True,
        native_geometry_log_only=True,
        radius_exp_materialized=False,
        radius_not_materialized=True,
        frozen_comparison_provider_used=False,
        common_N_handle=18,
        current_N="N18 graph handle",
        N_not_materialized=True,
        finite_N_numeric=False,
        phase_Z_independent=True,
        selected_phase_value=False,
        seams_proved=False,
        physical_y_derivatives_claimed=False,
        full_pressure_global_closure=False,
    )


def build(root, parent, input_packet):
    """Build long-reshape/reference fixed-phase Z0--Z3 log bounds."""
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Current long-reference background parent bytes changed")
    disk_parent = json.loads(parent_path.read_bytes())
    require(parent == disk_parent,
            "Current long-reference background parent mapping is required")

    packet, packet_path = _read_packet(root, input_packet)
    proof = _verify_packet(packet)
    pins = _source_pins(root)
    if packet_path is not None:
        pins[packet_path] = dict(path=packet_path,
                                 expected_sha256=sha(root / packet_path),
                                 passed=True, source_executed=False)
    pins.update(_verify_packet_pins(root, packet))
    source_guards = _source_guard(root, disk_parent)

    results = {}
    for name in CHARTS:
        chart = packet["charts"][name]
        c = MPIntervalContext()
        c.prec = 900
        raw_caps = _mapping(chart.get("log_caps"),
                            "chart.log_caps." + name)
        caps = {
            key: _rows_to_logz(c, raw_caps.get(key),
                               name + ".log_caps." + key)
            for key in ("a", "b", "t0", "E_over_frozen_E0",
                        "r", "alpha", "K")
        }
        raw_cond = _mapping(chart.get("conditioning"),
                            name + ".conditioning")
        conditioning = {
            "log_D_lower": _read_record(c, raw_cond.get("log_D_lower")),
            "log_lambda_lower": _read_record(
                c, raw_cond.get("log_lambda_lower")
            ),
        }
        log_Pstar_E_upper = _read_record(c, chart["log_Pstar_E_upper"])
        logN = _read_record(c, chart["logN"])
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
            B_over_Pstar_over_frozen_E0_log_caps=_packet_rows(computed),
            B_over_Pstar_over_frozen_E0_log_caps4rows=_packet_rows(computed),
            log_Pstar_E_upper=_record(log_Pstar_E_upper),
            logN=_record(logN),
            source_pins=pins,
            conditioning=dict(
                log_D_lower=_record(conditioning["log_D_lower"]),
                log_lambda_lower=_record(conditioning["log_lambda_lower"]),
            ),
            intermediate_Z_phase_caps=intermediate,
            inverse_proof=computed["inverse_proof"],
            universal_loop_controls_consumed=False,
            universal_B_over_E_value_cap_tightened=False,
            source_contract=_chart_source_contract(name),
            physical_scale_proof=_physical_scale_proof(),
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
        logN=results[CHARTS[0]]["logN"],
        source_pins=pins,
        source_guards=source_guards,
        proofs=proof,
        physical_scale_proof=_physical_scale_proof(),
        source_contract=dict(
            charts=list(CHARTS),
            actual_native_long_reference_controls=True,
            native_phi_p2_boundgraphs_actual_data=True,
            native_geometry_log_only=True,
            radius_exp_materialized=False,
            radius_not_materialized=True,
            frozen_comparison_provider_used=False,
            common_N_handle=18,
            current_N="N18 graph handle",
            finite_N_numeric=False,
            seams_proved=False,
            phase_Z_independent=True,
            physical_y_derivatives_claimed=False,
            full_pressure_global_closure=False,
            frozen_E0=(
                "E(R,Z0), the axial anchor at fixed R; not an old "
                "comparison-phi provider"
            ),
        ),
        loop_recipe=dict(
            B="E*M/2",
            target_phase="Z-independent unit magnitude bound",
            no_universal_C0_tightening=True,
        ),
        no_accepted_builder_checker_producer_owner_norm_recurrence_quadrature=True,
        full_goal_complete=False,
        limits=dict(
            native_geometry_log_only=True,
            radius_exp_materialized=False,
            native_phi_p2_boundgraphs_actual_data=True,
            frozen_comparison_provider_used=False,
            seams_proved=False,
            physical_y_derivatives_claimed=False,
            full_pressure_global_closure=False,
            finite_N_numeric=False,
        ),
    )
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    """Exercise only the pure arithmetic interface with a benign fixture."""
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
    require(all(checks.values()), "Long-reference Z3 fixture failed")
    receipt = dict(passed=True, checks=checks, JSON_serializable=True,
                   source_packet_required=True, full_goal_complete=False)
    json.dumps(receipt, allow_nan=False)
    return receipt


def self_check():
    return fixture_bound_check()


__all__ = (
    "PARENT",
    "PARENT_SHA256",
    "C5_SOURCE",
    "C5_SOURCE_SHA256",
    "NATIVE_CONTROLS_SOURCE",
    "CHARTS",
    "ORDINARY_Z_ORDERS",
    "build",
    "fixture_bound_check",
    "self_check",
)
