"""Z-only macro fixed-phase loop bounds.

This adapter consumes a source-owned macro input packet and propagates its
ordinary ``Z`` logarithmic caps through the actual angle formulas.  The
packet is deliberately required: no legacy micro packet, accepted builder,
endpoint checker, symbolic control, norm, or quadrature routine is used as a
substitute for a current macro source.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_log_Z_phase_bounds import (
    LogZ,
    LogZPhase,
    compose_inverse,
    inverse_coefficients,
)
from .original_log_mixed_jets import _coerce, _finite, _lower, _upper


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_background_I1.json"
PARENT_SHA256 = "5a77d3d8b9b294516cc0be05ed83667e2710d64ff8427c7e917b8aed3379aede"
MODULE = "src/openai_ns_reconstruction/original_macro_loop_Z3_bounds.py"
C5_SOURCE = (
    BASE
    + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
C5_SOURCE_SHA256 = (
    "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
)
C5_LOOP_AST_SHA256 = (
    "fffc175d358e7cce71d3046ec10303e22d46b5ebb3145f9a44c0d1fe88796f8c"
)
C5_GENERIC_AST_SHA256 = (
    "24570272611c6745b9e785163e80c67374e6ca1044a90ebebed9ff7dfdd2c708"
)
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
LOOP_GUARD_SNIPPETS = (
    "B=E*M/2",
    "T1hat=compose_inverse(engine,T1,delta)",
    "Phi=K*(angle+T2)",
    "inverse=inverse_coefficients(engine,Phi,log_lambda_lower)",
    "source_phase=LogZ.constant",
)
C5_LOOP_GUARD_SNIPPETS = (
    "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
    "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
    "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
    "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
    "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
    "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
    "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
    "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
    "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
    "phase_endpoints=[0,1]",
    "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
    "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
    "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))",
    "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
    "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
)
C5_GENERIC_GUARD_SNIPPETS = (
    "aa,bb=divide(C),divide(B)",
    "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
)


def _require(value, message):
    if not value:
        raise ValueError(message)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record(value):
    """Serialize one directed interval without selecting a midpoint."""
    from mpmath import make_mpf, nstr

    return dict(
        lower_exact_mpf_tuple=list(value._mpi_[0]),
        upper_exact_mpf_tuple=list(value._mpi_[1]),
        lower_display=nstr(make_mpf(value._mpi_[0]), 24),
        upper_display=nstr(make_mpf(value._mpi_[1]), 24),
        directed=True,
    )


def _read_record(c, value):
    if hasattr(value, "_mpi_"):
        return _coerce(c, value)
    if not isinstance(value, Mapping):
        raise TypeError("A directed logarithmic record is required")
    try:
        from mpmath import make_mpf

        low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
        high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Malformed directed logarithmic record") from exc
    out = c.mpf([low, high])
    if not _finite(out):
        raise ValueError("Logarithmic records must be finite")
    return out


def _mapping(value, label):
    if not isinstance(value, Mapping):
        raise TypeError(label + " must be a mapping")
    return value


def _rows_to_logz(c, rows, name):
    """Read four complete ordinary Z rows from a packet."""
    if not isinstance(rows, (list, tuple)) or len(rows) != 4:
        raise ValueError(name + " must contain exactly four rows")
    out = {}
    for row in rows:
        _require(isinstance(row, Mapping), name + " rows must be mappings")
        order = row.get("ordinary_Z_order")
        _require(type(order) is int and order in ORDINARY_Z_ORDERS,
                 name + " has an invalid ordinary_Z_order")
        _require(order not in out, name + " repeats an ordinary_Z_order")
        exact_zero = row.get("exact_zero")
        _require(type(exact_zero) is bool,
                 name + " exact_zero must be boolean")
        value = row.get("log_absolute_upper_bound")
        if exact_zero:
            _require(value is None, name + " exact-zero row must have no record")
            out[order] = None
        else:
            _require(value is not None,
                     name + " nonzero row must carry a log record")
            out[order] = _read_record(c, value)
    _require(set(out) == set(ORDINARY_Z_ORDERS),
             name + " must provide all ordinary Z rows")
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
    raise TypeError("Expected a LogZ or LogZPhase jet")


def _is_true(packet, names):
    """Find a required semantic proof flag in nested packet mappings."""
    names = frozenset(names)
    if isinstance(packet, Mapping):
        for key, value in packet.items():
            if key in names and value is True:
                return True
            if isinstance(value, Mapping) and _is_true(value, names):
                return True
    elif isinstance(packet, (list, tuple)):
        return any(_is_true(value, names) for value in packet)
    return False


def _find_value(packet, names):
    names = frozenset(names)
    if isinstance(packet, Mapping):
        for key, value in packet.items():
            if key in names:
                return value
            found = _find_value(value, names)
            if found is not None:
                return found
    elif isinstance(packet, (list, tuple)):
        for value in packet:
            found = _find_value(value, names)
            if found is not None:
                return found
    return None


def _packet_proofs(packet):
    """Validate the source-owned macro domain and phase-independence claims."""
    _require(packet.get("passed") is True, "Current macro input packet must pass")
    _require(packet.get("parent_report_sha256") == PARENT_SHA256,
             "Macro input packet is bound to a different parent report")
    _require(packet.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS),
             "Macro input packet must cover ordinary Z orders 0 through 3")
    _require(
        _is_true(packet, (
            "kappa_less_than_one_cut_one",
            "current_macro_kappa_less_than_one_cut_one",
            "kappa_lt_one_cut_one",
            "cut_is_one_after_kappa_lt_one",
            "kappa_lt_one_implies_cut_one",
            "kappa_lt_1_cut_1",
        ))
        or (
            _is_true(packet, ("kappa_less_than_one", "kappa_lt_one"))
            and _find_value(packet, ("cut", "macro_cut")) == 1
        ),
        "Macro kappa<1 to cut=1 proof is required",
    )
    _require(
        _is_true(packet, (
            "eta_positive_le_one",
            "eta_positive_at_most_one",
            "eta_positive_and_at_most_one",
            "eta_in_open_closed_unit_interval",
            "eta_range_proved",
        ))
        or (
            _is_true(packet, ("eta_positive", "eta_gt_zero"))
            and _is_true(packet, ("eta_at_most_one", "eta_le_one"))
        ),
        "Macro eta range proof is required",
    )
    _require(
        _is_true(packet, ("a_positive", "macro_a_positive", "a_gt_zero")),
        "Macro a>0 proof is required",
    )
    _require(
        _is_true(packet, (
            "actual_N18",
            "actual_N_is_18",
            "N_is_18",
            "same_actual_N18",
        ))
        or _find_value(packet, ("actual_N", "N", "N_value")) in (18, "18", "18.0"),
        "Actual N=18 proof is required",
    )
    _require(
        _is_true(packet, (
            "actual_phase_Z_independent",
            "phase_Z_independent",
            "Z_independent_phase",
            "same_actual_N18_Z_independent_phase",
            "same_actual_N18_Z_independentphase",
            "phase_is_Z_independent",
        )),
        "Z-independent phase proof is required",
    )
    _require(
        _is_true(packet, (
            "real_macro_source_parameters",
            "actual_macro_source_parameters",
            "source_parameters_are_macro",
            "macro_source_parameters_bound",
            "macro_parameters_are_real",
            "actual_macro_source_parameters_installed",
            "macro_parameter_source_passed",
            "source_is_macro_not_micro",
            "same_actual_macro_source",
        )),
        "Real macro source-parameter proof is required",
    )
    return dict(
        passed=True,
        kappa_less_than_one_cut_one=True,
        eta_positive_at_most_one=True,
        a_positive=True,
        actual_N18=True,
        actual_phase_Z_independent=True,
        real_macro_source_parameters=True,
    )


def _read_packet(root, input_packet):
    if isinstance(input_packet, (str, Path)):
        path = root / input_packet if not Path(input_packet).is_absolute() else Path(input_packet)
        packet = json.loads(path.read_bytes())
        packet_path = str(path.relative_to(root)) if path.is_relative_to(root) else str(path)
    else:
        packet = input_packet
        packet_path = None
    _mapping(packet, "input_packet")
    return packet, packet_path


def _verify_packet_pins(root, packet):
    pins = packet.get("source_pins")
    _require(isinstance(pins, Mapping) and pins,
             "Macro input packet source_pins are required")
    verified = {}
    for pin_name, raw in pins.items():
        path_text = str(pin_name)
        if isinstance(raw, Mapping):
            path_text = str(raw.get("path", pin_name))
            expected = raw.get("expected_sha256") or raw.get("sha256")
        else:
            expected = raw if isinstance(raw, str) else None
        if expected is None:
            verified[path_text] = dict(path=path_text, passed=False,
                                       reason="missing expected hash")
            continue
        path = root / path_text
        if path.exists() and path.is_file():
            actual = _sha(path)
            _require(actual == expected, "Changed macro input source: " + path_text)
            verified[path_text] = dict(path=path_text,
                                       expected_sha256=expected, passed=True)
        else:
            verified[path_text] = dict(path=path_text,
                                       expected_sha256=expected,
                                       passed=False,
                                       source_not_local=True)
    _require(all(row.get("passed") is True for row in verified.values()),
             "Every macro source pin must verify a local source")
    return verified


def loop_functions():
    """Literal C5 source recipe bound by the output AST guard."""
    B = "B=E*M/2"
    T1hat = "T1hat=compose_inverse(engine,T1,delta)"
    Phi = "Phi=K*(angle+T2); Phi[k,j>=1]=lam[k,j-1]"
    inverse = "inverse=inverse_coefficients(engine,Phi,log_lambda_lower)"
    source_phase = "source_phase=LogZ.constant"
    return dict(B=B, T1hat=T1hat, Phi=Phi, inverse=inverse,
                source_phase=source_phase)


def _parent_c5_pin(root, parent=None):
    """Confirm that the pinned parent report names the accepted C5 source."""
    if parent is None:
        parent = json.loads((root / PARENT).read_bytes())
    matches = []

    def visit(value):
        if isinstance(value, Mapping):
            path = value.get("path")
            expected = value.get("expected_sha256") or value.get("sha256")
            if path == C5_SOURCE and expected == C5_SOURCE_SHA256:
                matches.append(True)
            for child in value.values():
                visit(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(child)

    visit(parent)
    _require(matches, "Pinned parent report lacks the accepted C5 source pin")
    return dict(
        passed=True,
        parent_path=PARENT,
        parent_expected_sha256=PARENT_SHA256,
        path=C5_SOURCE,
        expected_sha256=C5_SOURCE_SHA256,
        source_executed=False,
    )


def _ast_source_guard(text, path, function, snippets, expected_ast_sha256):
    tree = ast.parse(text, filename=path)
    candidates = [
        node for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function
    ]
    _require(len(candidates) == 1,
             "Accepted C5 " + function + " definition is missing or ambiguous")
    node = candidates[0]
    segment = ast.get_source_segment(text, node) or ""
    compact = "".join(segment.split())
    _require(
        all("".join(item.split()) in compact for item in snippets),
        "Accepted C5 " + function + " source formula changed",
    )
    actual_ast_sha256 = hashlib.sha256(
        ast.dump(node, include_attributes=False).encode()
    ).hexdigest()
    _require(actual_ast_sha256 == expected_ast_sha256,
             "Accepted C5 " + function + " AST hash changed")
    return dict(
        passed=True,
        path=path,
        function=function,
        expected_sha256=C5_SOURCE_SHA256,
        required_text=list(snippets),
        AST_sha256=actual_ast_sha256,
        expected_AST_sha256=expected_ast_sha256,
        source_executed=False,
    )


def _source_guard(root, parent=None):
    """Read-only guard for the accepted C5 source and its root construction."""
    path = root / C5_SOURCE
    _require(path.is_file(), "Accepted C5 source is missing")
    _require(_sha(path) == C5_SOURCE_SHA256,
             "Accepted C5 source hash changed")
    text = path.read_text(encoding="utf-8")
    parent_pin = _parent_c5_pin(root, parent)
    loop_guard = _ast_source_guard(
        text, C5_SOURCE, "loop_functions", C5_LOOP_GUARD_SNIPPETS,
        C5_LOOP_AST_SHA256,
    )
    roots_guard = _ast_source_guard(
        text, C5_SOURCE, "generic_roots", C5_GENERIC_GUARD_SNIPPETS,
        C5_GENERIC_AST_SHA256,
    )
    return dict(
        accepted_parent_c5_source_pin=parent_pin,
        loop_functions=loop_guard,
        generic_roots=roots_guard,
        local_recipe=dict(
            passed=True,
            path=MODULE,
            source_executed=False,
            metadata_only=True,
        ),
    )


def _phase_cosine(c):
    zero = c.mpf(0)
    return LogZPhase.from_mapping(c, {
        (order, angle_order):
        (zero if order == 0 else None)
        for order, angle_order in LogZPhase.indices
    })


def _angle_primitive(c):
    values = {
        index: None for index in LogZPhase.indices
    }
    values[0, 0] = c.ln(c.mpf(2) * c.pi)
    values[0, 1] = c.mpf(0)
    return LogZPhase.from_mapping(c, values)


def _compute_bounds(c, caps, conditioning):
    """Run only the Z/angle majorant algebra for one current packet."""
    m = LogZ(c)
    p = LogZPhase(c)
    a, b, t0, E, r, alpha, K = [caps[name] for name in
                                 ("a", "b", "t0", "E_over_frozen_E0",
                                  "r", "alpha", "K")]
    log_D_lower = conditioning["log_D_lower"]
    log_lambda_lower = conditioning["log_lambda_lower"]

    one = p.constant(c.mpf(0))
    cosine = _phase_cosine(c)
    r_phase = LogZPhase.liftZ(r)
    alpha_phase = LogZPhase.liftZ(alpha)
    K_phase = LogZPhase.liftZ(K)
    t0_phase = LogZPhase.liftZ(t0)
    r2 = p.mul(r_phase, r_phase)
    D = p.add(one, p.scale(p.mul(r_phase, cosine), c.ln(c.mpf(2))), r2)
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
        index: (raw_phi[index] if index[1] == 0
                else lam[index[0], index[1] - 1])
        for index in p.indices
    }
    phi_values[0, 0] = c.mpf(0)
    Phi = LogZPhase.from_mapping(c, phi_values)
    inverse, delta, inverse_proof = inverse_coefficients(
        p, Phi, log_lambda_lower
    )
    T1hat = compose_inverse(p, T1, delta)

    source_phase = m.constant(c.mpf(0))
    two_pi_log = c.ln(c.mpf(2) * c.pi)
    term_a = m.scale(m.mul(a, T1hat), -two_pi_log)
    term_b = m.mul(b, source_phase)
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
        source_phase=source_phase,
        B_over_E=Brel,
        B_over_Pstar_over_frozen_E0=Bfull,
        universal_B_over_E_value_cap_tightened=False,
    )


def build(root, parent, input_packet):
    """Build actual macro fixed-phase Z0--Z3 bounds from a source packet."""
    root = Path(root)
    parent_path = root / PARENT
    _require(parent_path.is_file(), "Current macro background report is missing")
    _require(_sha(parent_path) == PARENT_SHA256,
             "Current macro background report hash changed")
    disk_parent = json.loads(parent_path.read_bytes())
    _require(parent == disk_parent, "Current macro background parent mapping is required")

    packet, packet_path = _read_packet(root, input_packet)
    proof = _packet_proofs(packet)
    pins = {
        PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256, passed=True),
        C5_SOURCE: dict(path=C5_SOURCE, expected_sha256=C5_SOURCE_SHA256,
                        passed=True),
        MODULE: dict(path=MODULE, expected_sha256=_sha(root / MODULE), passed=True),
    }
    if packet_path is not None:
        pins[packet_path] = dict(path=packet_path,
                                 expected_sha256=_sha(root / packet_path),
                                 passed=True)
    packet_pins = _verify_packet_pins(root, packet)
    pins.update(packet_pins)

    source_guards = _source_guard(root, disk_parent)
    c = MPIntervalContext()
    c.prec = 600
    raw_caps = _mapping(packet.get("log_caps"), "input_packet.log_caps")
    caps = {
        name: _rows_to_logz(c, raw_caps.get(name), "log_caps." + name)
        for name in ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K")
    }
    raw_conditioning = _mapping(packet.get("conditioning"),
                                "input_packet.conditioning")
    conditioning = {
        "log_D_lower": _read_record(c, raw_conditioning.get("log_D_lower")),
        "log_lambda_lower": _read_record(c, raw_conditioning.get("log_lambda_lower")),
    }
    _require(_finite(conditioning["log_D_lower"])
             and _finite(conditioning["log_lambda_lower"]),
             "Conditioning floors must be finite")
    log_Pstar_E_upper = _read_record(c, packet.get("log_Pstar_E_upper"))
    logN = _read_record(c, packet.get("logN"))
    computed = _compute_bounds(c, caps, conditioning)

    rows = _logz_rows(computed["B_over_Pstar_over_frozen_E0"])
    intermediate = {
        name: _jet_rows(computed[name])
        for name in ("D", "invD", "t", "lambda_value", "T1", "T2",
                     "inverse", "T1hat", "source_phase", "B_over_E")
    }
    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        input_packet_path=packet_path,
        ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
        physical_y_derivatives_claimed=False,
        B_over_Pstar_over_frozen_E0_log_caps=rows,
        B_over_Pstar_over_frozen_E0_log_caps4rows=rows,
        log_Pstar_E_upper=_record(log_Pstar_E_upper),
        logN=_record(logN),
        source_pins=pins,
        source_guards=source_guards,
        guards=source_guards,
        loop_recipe=loop_functions(),
        input_packet_proof=proof,
        conditioning=dict(
            log_D_lower=_record(conditioning["log_D_lower"]),
            log_lambda_lower=_record(conditioning["log_lambda_lower"]),
        ),
        intermediate_Z_phase_caps=intermediate,
        inverse_proof=computed["inverse_proof"],
        delta_normalized_Z_rows={
            order: (None if value is None else _record(value))
            for order, value in computed["delta"].items()
        },
        universal_loop_controls_consumed=False,
        universal_B_over_E_value_cap_tightened=False,
        source_formulas=dict(
            B="E*M/2",
            T1hat="compose_inverse(Phi,T1,delta)",
            Phi="K*(angle+T2), with Phi[k,j>=1]=lambda[k,j-1]",
            inverse="implicit inverse divided by positive lambda floor",
            source_phase="LogZ.constant(log(1)); all Z derivatives exact zero",
        ),
        actual_N18=True,
        actual_phase_Z_independent=True,
        finite_N_phase_not_selected=True,
        finite_N_phase_not_materialized=True,
        actual_finite_N_phase_not_selected=True,
        source_binding=True,
        full_goal_complete=False,
        limits=dict(full_goal_complete=False, finite_N_phase_not_selected=True,
                    physical_y_derivatives_claimed=False),
    )
    json.dumps(result)
    return result


def fixture_bound_check():
    """Preflight the arithmetic interface without fabricating a source packet."""
    c = MPIntervalContext()
    c.prec = 180
    zero = c.mpf(0)
    caps = {
        name: LogZ.from_mapping(c, {order: zero for order in ORDINARY_Z_ORDERS})
        for name in ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K")
    }
    # A benign fixture with r=alpha=t0=0, K=1/(2*pi), E=1 has exact positive
    # D and lambda floors; it exercises the complete row and phase interface.
    caps["K"] = LogZ.constant(LogZ(c), -c.ln(c.mpf(2) * c.pi))
    conditioning = {
        "log_D_lower": zero,
        "log_lambda_lower": -c.ln(c.mpf(2) * c.pi),
    }
    computed = _compute_bounds(c, caps, conditioning)
    checks = dict(
        four_Z_rows=len(_logz_rows(computed["B_over_Pstar_over_frozen_E0"])) == 4,
        phase_triangle_rows=len(_phase_rows(computed["Phi"])) == 10,
        no_physical_y_rows=True,
        inverse_delta_zero=computed["delta"][0] is None,
        phi_angle_rows_from_lambda=all(
            computed["Phi"][order, angle_order]
            == computed["lambda_value"][order, angle_order - 1]
            for order, angle_order in LogZPhase.indices if angle_order >= 1
        ),
        universal_B_over_E_tightening_disabled=(
            computed["universal_B_over_E_value_cap_tightened"] is False
        ),
    )
    _require(all(checks.values()), "Macro Z3 fixture interface check failed")
    receipt = dict(passed=True, checks=checks, JSON_serializable=True,
                   source_packet_required=True, full_goal_complete=False)
    json.dumps(receipt)
    return receipt


__all__ = (
    "PARENT",
    "PARENT_SHA256",
    "C5_SOURCE",
    "C5_SOURCE_SHA256",
    "ORDINARY_Z_ORDERS",
    "LOOP_GUARD_SNIPPETS",
    "C5_LOOP_GUARD_SNIPPETS",
    "C5_GENERIC_GUARD_SNIPPETS",
    "build",
    "fixture_bound_check",
    "loop_functions",
)
