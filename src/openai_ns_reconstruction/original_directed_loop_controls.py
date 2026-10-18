"""Pure fixed-phase C5 loop controls from the saved native packet.

This module is deliberately a small boundary between the accepted native
control packet and the next, still conditional, loop layer.  It reads only
the serialized ``raw_interval_tree`` supplied by the caller.  Every value is
an outward logarithmic majorant; no signed source value, midpoint, exponent
of a stored logarithm, old loop producer, or physical-y derivative is used.

The formulas are the prescribed C5 source formulas.  ``LogZPhase`` carries
ordinary Z/angle majorants through total order three, ``inverse_coefficients``
solves the fixed-target phase inverse, and ``compose_inverse`` evaluates the
three free-angle jets at that inverse.  The public rows are ordinary Z0--Z3
records, while the intermediate phase rows retain all ten canonical
(Z,angle) entries.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from fractions import Fraction

from mpmath import isfinite, make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext

from .original_log_Z_phase_bounds import (
    LogZ,
    LogZPhase,
    compose_inverse,
    inverse_coefficients,
)
from .original_log_mixed_jets import _coerce


BASE = "experiments/root_st073/"
MODULE = "src/openai_ns_reconstruction/original_directed_loop_controls.py"
C5_SOURCE = (
    BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
C5_SOURCE_SHA256 = (
    "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
)
LOG_Z_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
LOG_Z_SOURCE_SHA256 = (
    "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59"
)
NATIVE_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
NATIVE_REPORT_BYTES = 42177499
CHARTS = ("long_reshape", "reference", "restoration", "post_restore")
LONG_CHARTS = ("long_reshape", "reference")
RESTORATION_CHARTS = ("restoration", "post_restore")
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
PUBLIC_ROOT_NAMES = ("a", "b", "t0", "r", "alpha", "K")
PUBLIC_ROWS = (
    "A",
    "B_over_E",
    "B_over_Pstar_over_frozen_E0",
    "A_phi",
    "B_over_Pstar_over_frozen_E0_phi",
    "inverse_Psi",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _raw_endpoint(value):
    """Decode one exact MPF endpoint without rounding through a float."""
    return make_mpf(tuple(value))


def _read_record(c: MPIntervalContext, value, label: str = "interval"):
    """Read a directed interval and preserve both endpoints.

    ``_coerce`` is intentionally used at this boundary.  Its interval branch
    copies the lower and upper MPF endpoints into the supplied context, so the
    later reciprocal operation can negate the interval in directed arithmetic
    instead of negating a selected scalar endpoint.
    """
    if hasattr(value, "_mpi_"):
        out = _coerce(c, value)
    else:
        _require(isinstance(value, Mapping), label + " must be a directed record")
        _require(value.get("directed") is True,
                 label + " must set directed=True")
        try:
            lower = _raw_endpoint(value["lower_exact_mpf_tuple"])
            upper = _raw_endpoint(value["upper_exact_mpf_tuple"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Malformed " + label) from exc
        _require(lower <= upper, "Reversed " + label)
        out = _coerce(c, c.mpf([lower, upper]))
    _require(isfinite(out.a) and isfinite(out.b),
             label + " must have finite endpoints")
    _require(out.a <= out.b, "Reversed " + label)
    return out


def _record(value):
    """Serialize one interval without selecting a midpoint."""
    _require(hasattr(value, "_mpi_"), "An interval is required for serialization")
    lower = tuple(value._mpi_[0])
    upper = tuple(value._mpi_[1])
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_display": nstr(make_mpf(lower), 24),
        "upper_display": nstr(make_mpf(upper), 24),
        "directed": True,
    }


def _logz_rows(jet: LogZ) -> list[dict]:
    return [
        {
            "ordinary_Z_order": order,
            "exact_zero": jet[order] is None,
            "log_absolute_upper_bound": (
                None if jet[order] is None else _record(jet[order])
            ),
        }
        for order in ORDINARY_Z_ORDERS
    ]


def _phase_rows(jet: LogZPhase) -> list[dict]:
    return [
        {
            "Z_order": order,
            "free_angle_order": angle_order,
            "exact_zero": jet[order, angle_order] is None,
            "log_absolute_upper_bound": (
                None
                if jet[order, angle_order] is None
                else _record(jet[order, angle_order])
            ),
        }
        for order, angle_order in LogZPhase.indices
    ]


def _row_record(c: MPIntervalContext, row, label: str):
    """Decode a packed ordinary-Z row from the native packet."""
    _require(isinstance(row, Mapping), label + " row must be a mapping")
    order = row.get("ordinary_Z_order")
    _require(type(order) is int and order in ORDINARY_Z_ORDERS,
             label + " has an invalid ordinary_Z_order")
    exact_zero = row.get("exact_zero")
    _require(type(exact_zero) is bool,
             label + " exact_zero must be boolean")
    raw = row.get("log_absolute_upper_bound")
    if exact_zero:
        _require(raw is None, label + " exact-zero row carries a record")
        return order, None
    _require(raw is not None, label + " nonzero row lacks a record")
    return order, _read_record(c, raw, label + "[" + str(order) + "]")


def _rows_to_logz(c: MPIntervalContext, rows, label: str, *, count: int = 4):
    _require(isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)),
             label + " must be a packed row sequence")
    _require(len(rows) == count,
             label + " must contain exactly " + str(count) + " rows")
    values = {}
    for row in rows:
        order, value = _row_record(c, row, label)
        _require(order not in values, label + " repeats an ordinary_Z_order")
        values[order] = value
    _require(set(values) == set(ORDINARY_Z_ORDERS),
             label + " must cover ordinary Z orders 0 through 3")
    return LogZ.from_mapping(c, values)


def _native_section(native_report: Mapping) -> Mapping:
    """Resolve the v2 section after immutable parent admission."""
    _require(isinstance(native_report, Mapping),
             "Decoded native parent must be a mapping")
    section = native_report.get("actual_native_directed_controls_v2")
    _require(isinstance(section, Mapping),
             "actual_native_directed_controls_v2 section is required")
    _require(section.get("schema") == "original_native_directed_controls_v2_v1",
             "Native parent section schema changed")
    _require(section.get("passed") is True,
             "Native parent section did not pass")
    scope = section.get("scope")
    _require(isinstance(scope, Mapping), "Native parent section scope is required")
    _require(scope.get("full_corrected_native_controls_standalone") is True
             and scope.get("corrected_native_controls_verified") is True
             and scope.get("inherited_loop_revalidated") is False
             and scope.get("current_N_revalidated") is False
             and scope.get("inherited_log_caps_revalidated") is False,
             "Native parent section scope changed")
    tree = section.get("raw_interval_tree")
    _require(isinstance(tree, Mapping),
             "The source-pinned raw_interval_tree is required")
    _require(set(tree) == {"shared", "long_reference", "restoration"},
             "Unexpected native raw_interval_tree sections")
    return section


def _admit_native_report_bytes(native_report_bytes: bytes) -> Mapping:
    """Admit only the exact reviewed native parent before JSON decoding."""
    _require(type(native_report_bytes) is bytes,
             "compute requires native parent bytes")
    _require(len(native_report_bytes) == NATIVE_REPORT_BYTES,
             "Native parent byte length changed")
    digest = hashlib.sha256(native_report_bytes).hexdigest()
    _require(digest == NATIVE_REPORT_SHA256,
             "Native parent SHA256 changed")
    try:
        native_report = json.loads(native_report_bytes)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("Admitted native parent is not valid JSON") from exc
    return native_report


def _chart_source(section: Mapping, name: str) -> Mapping:
    tree = section["raw_interval_tree"]
    family = "long_reference" if name in LONG_CHARTS else "restoration"
    part = tree.get(family)
    _require(isinstance(part, Mapping), "Missing raw " + family + " interval tree")
    charts = part.get("charts")
    _require(isinstance(charts, Mapping), "Missing raw " + family + " charts")
    chart = charts.get(name)
    _require(isinstance(chart, Mapping), "Missing raw source chart: " + name)
    return chart


def _source_pin_receipt(section: Mapping) -> dict:
    """Return source identity without executing or loading a producer."""
    pins = section.get("source_pins")
    _require(isinstance(pins, Mapping), "Native source pins are required")
    found = {}
    for path, expected in (
        (C5_SOURCE, C5_SOURCE_SHA256),
        (LOG_Z_SOURCE, LOG_Z_SOURCE_SHA256),
    ):
        pin = pins.get(path)
        _require(isinstance(pin, Mapping), "Missing source pin: " + path)
        _require(pin.get("path") == path
                 and pin.get("expected_sha256") == expected
                 and pin.get("actual_sha256") == expected
                 and pin.get("passed") is True
                 and pin.get("source_executed") is False,
                 "Changed or executed source pin: " + path)
        found[path] = expected
    return found


def _phase_cosine(c: MPIntervalContext) -> LogZPhase:
    """Majorants for cos(psi): all angle derivatives have unit cap."""
    zero = c.mpf(0)
    return LogZPhase.from_mapping(
        c,
        {
            (order, angle_order): (zero if order == 0 else None)
            for order, angle_order in LogZPhase.indices
        },
    )


def _angle_primitive(c: MPIntervalContext) -> LogZPhase:
    values = {index: None for index in LogZPhase.indices}
    values[0, 0] = c.ln(c.mpf(2) * c.pi)
    values[0, 1] = c.mpf(0)
    return LogZPhase.from_mapping(c, values)


def _compute_chart(c: MPIntervalContext, caps: Mapping[str, LogZ],
                   E_relative: LogZ, conditioning: Mapping[str, object]):
    """Propagate one source chart through the prescribed C5 formulas."""
    m = LogZ(c)
    p = LogZPhase(c)
    a, b, t0, r, alpha, K = [caps[name] for name in PUBLIC_ROOT_NAMES]
    log_D_lower = conditioning["log_D_lower"]
    log_lambda_lower = conditioning["log_lambda_lower"]

    one_phase = p.constant(c.mpf(0))
    cosine = _phase_cosine(c)
    r_phase = LogZPhase.liftZ(r)
    alpha_phase = LogZPhase.liftZ(alpha)
    K_phase = LogZPhase.liftZ(K)
    t0_phase = LogZPhase.liftZ(t0)

    # D=1-2*r*cos(psi)+r^2.  Absolute majorants intentionally ignore
    # subtraction signs; positivity enters only through the saved lower floor.
    r2 = p.mul(r_phase, r_phase)
    D = p.add(
        one_phase,
        p.scale(p.mul(r_phase, cosine), c.ln(2)),
        r2,
    )
    invD = p.reciprocal(D, log_D_lower)
    cos_minus_r = p.add(cosine, r_phase)
    t = p.add(
        t0_phase,
        p.mul(alpha_phase, p.mul(cos_minus_r, invD)),
    )
    t2 = p.mul(t, t)
    lam = p.mul(K_phase, p.add(one_phase, t2))

    # T1 and T2 are the source angle primitives.  Phi's angle derivatives
    # come directly from lambda; its value remains the conservative raw bound
    # unless a separately reviewed endpoint-monotonicity theorem is consumed.
    T1 = p.primitiveAngle(t)
    T2 = p.primitiveAngle(t2)
    angle = _angle_primitive(c)
    raw_phi = p.mul(K_phase, p.add(angle, T2))
    phi_values = {
        index: (
            raw_phi[index]
            if index[1] == 0
            else lam[index[0], index[1] - 1]
        )
        for index in p.indices
    }
    # The value row is the uncentered raw Phi cap. The fixed-target inverse
    # recurrence cancels its constant term and never consumes Phi[0,0].
    Phi = LogZPhase.from_mapping(c, phi_values)

    inverse, delta, inverse_proof = inverse_coefficients(
        p, Phi, log_lambda_lower
    )
    T1hat = compose_inverse(p, T1, delta)
    that = compose_inverse(p, t, delta)
    lamhat = compose_inverse(p, lam, delta)

    # The issued target phase is Z independent and lies in [0,1].  This is a
    # magnitude cap, not a selected phase value, so all nonconstant rows are
    # exact zeros and the value row has log cap zero.
    target_phase = m.constant(c.mpf(0))
    two_pi_log = c.ln(c.mpf(2) * c.pi)
    half_log = c.ln(c.mpf(2))
    chi = m.add(target_phase, m.scale(inverse, -two_pi_log))
    A = m.scale(m.mul(a, chi), -half_log)

    # B_over_E=1/2[-a*T1hat/(2*pi)-b*phase], retaining the restoration b
    # and t0 rows even though signs are discarded by the majorant algebra.
    B_over_E = m.scale(
        m.add(
            m.scale(m.mul(a, T1hat), -two_pi_log),
            m.mul(b, target_phase),
        ),
        -half_log,
    )
    Bnorm = m.mul(E_relative, B_over_E)

    # Fast phase rows use the actual saved lambda lower floor; no 1/N factor
    # belongs in these controls.
    inv_lamhat = m.reciprocal(lamhat, log_lambda_lower)
    A_phi = m.scale(
        m.mul(a, m.add(m.constant(c.mpf(0)),
                       m.scale(inv_lamhat, -two_pi_log))),
        -half_log,
    )
    Bnorm_phi = m.mul(
        E_relative,
        m.scale(
            m.add(
                m.scale(m.mul(a, m.mul(that, inv_lamhat)), -two_pi_log),
                b,
            ),
            -half_log,
        ),
    )

    return {
        "D": D,
        "inverse_D": invD,
        "t": t,
        "lambda": lam,
        "T1": T1,
        "T2": T2,
        "Phi": Phi,
        "inverse": inverse,
        "delta": delta,
        "T1hat": T1hat,
        "that": that,
        "lamhat": lamhat,
        "target_phase": target_phase,
        "A": A,
        "B_over_E": B_over_E,
        "Bnorm": Bnorm,
        "A_phi": A_phi,
        "Bnorm_phi": Bnorm_phi,
        "inverse_proof": inverse_proof,
    }


def _chart_packet(c: MPIntervalContext, source: Mapping,
                  computed: Mapping[str, object], source_pin: Mapping[str, str]):
    """Serialize one chart while retaining source rows byte-for-byte."""
    root_caps_source = source["root_caps"]
    E_relative_source = source["E_relative"]
    conditioning_source = source["conditioning"]
    public_roots = {
        name: _logz_rows(computed["caps"][name])
        for name in PUBLIC_ROOT_NAMES
    }
    public_rows = {
        "A": _logz_rows(computed["A"]),
        "B_over_E": _logz_rows(computed["B_over_E"]),
        "B_over_Pstar_over_frozen_E0": _logz_rows(computed["Bnorm"]),
        "A_phi": _logz_rows(computed["A_phi"]),
        "B_over_Pstar_over_frozen_E0_phi": _logz_rows(computed["Bnorm_phi"]),
        "inverse_Psi": _logz_rows(computed["inverse"]),
    }
    # Only the free-angle objects are phase rows.  Composition at the
    # implicit inverse returns ordinary-Z jets, so T1hat/that/lamhat are
    # kept in the separate ordinary intermediate block below.
    phase_names = ("D", "inverse_D", "t", "lambda", "T1", "T2", "Phi")
    phase_caps = {name: _phase_rows(computed[name]) for name in phase_names}
    ordinary_intermediate = {
        name: _logz_rows(computed[name])
        for name in ("inverse", "T1hat", "that", "lamhat")
    }
    # The source E-relative packet has five rows; only the ordinary Z0..Z3
    # prefix participates in this component.  Preserve the complete input
    # under input_E_relative and expose the explicit truncated prefix too.
    E_relative4 = [copy.deepcopy(row) for row in E_relative_source[:4]]

    native_memory = source.get("incoming_memory")
    if native_memory is None:
        history = source.get("history")
        if isinstance(history, Mapping):
            native_memory = history.get("incoming")

    out = {
        "chart": source.get("chart"),
        "input_conditioning": copy.deepcopy(conditioning_source),
        "input_root_caps": copy.deepcopy(root_caps_source),
        "input_E_relative": copy.deepcopy(E_relative_source),
        "root_caps": public_roots,
        "E_relative": E_relative4,
        **public_rows,
        "phase_caps": phase_caps,
        "intermediate_Z_caps": ordinary_intermediate,
        "fixed_target_phase": True,
        "signed_value_selected": False,
        "phi_value_bound_one_used": False,
        "phi_value_raw_bound_retained": True,
        "psi_value_bound": "|Psi|<=2*pi; bound only, no value selected",
        "conditioning": copy.deepcopy(conditioning_source),
        "inverse_proof": copy.deepcopy(computed["inverse_proof"]),
        "formula_receipt": {
            "D": "1-2*r*cos(psi)+r^2",
            "t": "t0+alpha*(cos(psi)-r)/D",
            "lambda": "K*(1+t^2)",
            "T1": "angle primitive of t",
            "T2": "angle primitive of t^2",
            "Phi": "K*(psi+T2); angle derivatives from lambda",
            "inverse": "inverse_coefficients(Phi, log_lambda_lower)",
            "T1hat": "compose_inverse(T1, inverse_increment)",
            "that": "compose_inverse(t, inverse_increment)",
            "lamhat": "compose_inverse(lambda, inverse_increment)",
            "A": "a/2*(phase-Psi/(2*pi))",
            "B_over_E": "1/2*(-a*T1hat/(2*pi)-b*phase)",
            "Bnorm": "E_relative*B_over_E",
            "A_phi": "a/2*(1-1/(2*pi*lambda_hat))",
            "Bnorm_phi": (
                "E_relative/2*(-a*that/(2*pi*lambda_hat)-b)"
            ),
        },
        "source_hash": copy.deepcopy(source_pin),
        "native_memory_distinguished": True,
        "native_memory_input": copy.deepcopy(native_memory),
        "native_memory_fields": ["A", "B", "C"],
        "physical_y_derivatives_claimed": False,
    }
    # Keep the row names explicit for consumers that use the older naming
    # convention, while the six canonical public fields above remain the
    # admission interface.
    out["B_over_E_log_caps"] = copy.deepcopy(out["B_over_E"])
    out["B_over_Pstar_over_frozen_E0_log_caps"] = copy.deepcopy(
        out["B_over_Pstar_over_frozen_E0"]
    )
    out["B_over_Pstar_over_frozen_E0_phi_log_caps"] = copy.deepcopy(
        out["B_over_Pstar_over_frozen_E0_phi"]
    )
    out["inverse_Psi_log_caps"] = copy.deepcopy(out["inverse_Psi"])
    return out


def compute(c: MPIntervalContext, native_report_bytes: bytes) -> dict:
    """Compute from the exact reviewed native parent bytes only."""
    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    native_report = _admit_native_report_bytes(native_report_bytes)
    section = _native_section(native_report)
    source_pins = _source_pin_receipt(section)
    charts = {}
    for name in CHARTS:
        source = _chart_source(section, name)
        root_caps = source.get("root_caps")
        _require(isinstance(root_caps, Mapping),
                 "Missing native root_caps: " + name)
        caps = {
            key: _rows_to_logz(
                c, root_caps.get(key), name + ".root_caps." + key
            )
            for key in PUBLIC_ROOT_NAMES
        }
        E_relative_source = source.get("E_relative")
        _require(isinstance(E_relative_source, Sequence)
                 and not isinstance(E_relative_source, (str, bytes))
                 and len(E_relative_source) == 5,
                 name + ".E_relative must contain five rows")
        E_relative = _rows_to_logz(
            c, list(E_relative_source[:4]), name + ".E_relative[0:4]"
        )
        conditioning_source = source.get("conditioning")
        _require(isinstance(conditioning_source, Mapping),
                 "Missing native conditioning: " + name)
        conditioning = {
            "log_D_lower": _read_record(
                c, conditioning_source.get("log_D_lower"),
                name + ".conditioning.log_D_lower",
            ),
            "log_lambda_lower": _read_record(
                c, conditioning_source.get("log_lambda_lower"),
                name + ".conditioning.log_lambda_lower",
            ),
        }
        if name == "restoration":
            _require(caps["b"][0] is not None and caps["t0"][0] is not None,
                     "Restoration b and t0 value rows must remain nonzero")
        if name == "post_restore":
            _require(all(caps["b"][order] is None for order in ORDINARY_Z_ORDERS)
                     and all(caps["t0"][order] is None
                             for order in ORDINARY_Z_ORDERS),
                     "Post-restore b and t0 rows must be exact zero")
        computed = _compute_chart(c, caps, E_relative, conditioning)
        computed = dict(computed, caps=caps)
        charts[name] = _chart_packet(c, source, computed, source_pins)

    result = {
        "passed": True,
        "schema": "original_directed_loop_controls_v1",
        "sourcehash": C5_SOURCE_SHA256,
        "source_hashes": {
            C5_SOURCE: C5_SOURCE_SHA256,
            LOG_Z_SOURCE: LOG_Z_SOURCE_SHA256,
        },
        "native_parent_report_sha256": NATIVE_REPORT_SHA256,
        "native_parent_report_bytes": NATIVE_REPORT_BYTES,
        "charts": charts,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "phase_ordering": [list(index) for index in LogZPhase.indices],
        "checks": {
            "source_pinned_c5_loop_functions": True,
            "source_pinned_log_Z_phase_helper": True,
            "raw_interval_tree_consumed": True,
            "native_root_controls_recomputed": False,
            "native_conditioning_recomputed": False,
            "conditioning_records_copied_exactly": True,
            "E_relative_five_rows_consumed": True,
            "E_relative_truncated_to_first_four": True,
            "restoration_nonzero_b_and_t0_retained": True,
            "post_restore_b_and_t0_exact_zero": True,
            "majorant_signs_ignored": True,
            "stored_logs_exponentiated": False,
            "scalar_abs_or_endpoint_negation": False,
            "old_loop_compute_bounds_called": False,
            "old_loop_build_called": False,
            "old_fixture_called": False,
            "fixed_phase_Z0_through_Z3": True,
            "full_Z5_certified": False,
            "physical_y_derivatives_claimed": False,
            "exact_native_parent_bytes_admitted": True,
            "native_parent_hash_checked_before_json_load": True,
        },
        "scope": {
            "fixed_target_phase_Z0_through_Z3_verified": True,
            "legacy_loop_packets_recertified": False,
            "inherited_log_caps_revalidated": False,
            "full_C5_Z5_bounds": False,
            "physical_y_bounds": False,
            "current_N_revalidated": False,
            "full_pressure_revalidated": False,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
        },
        "execution_audit": {
            "accepted_ancestors_replayed": False,
            "native_builder_replayed": False,
            "old_loop_builders_called": False,
            "only_new_loop_inputs_propagated": True,
            "source_executed": False,
        },
    }
    # Keep the public result JSON serializable even when called directly from
    # a test with an in-memory source packet.
    json.dumps(result, allow_nan=False)
    return result


def _rational_log(c: MPIntervalContext, value) -> object:
    value = Fraction(value)
    if value == 0:
        return None
    numerator = value.numerator if value.numerator >= 0 else -value.numerator
    return c.ln(c.mpf(numerator) / value.denominator)


def _exact_phase_mapping(c: MPIntervalContext, expression, z, psi,
                         psi0: Fraction):
    import sympy as sp

    values = {}
    for order, angle_order in LogZPhase.indices:
        value = sp.diff(expression, z, order, psi, angle_order).subs(
            {z: 0, psi: sp.Rational(psi0.numerator, psi0.denominator)}
        )
        values[order, angle_order] = _rational_log(c, sp.Rational(value))
    return values


def _contains_log(c, cap, value) -> bool:
    value = Fraction(value)
    if value == 0:
        return cap is None
    if cap is None:
        return False
    expected = _rational_log(c, value)
    return cap.b >= expected.b


def _fixture_series(expression, z, psi, psi0, order=3):
    """Solve the fixed-target implicit series with independent SymPy algebra."""
    import sympy as sp

    target = expression.subs({z: 0, psi: psi0})
    coefficients = sp.symbols("c1:" + str(order + 1))
    series = psi0 + sum(coefficients[k - 1] * z**k
                        for k in range(1, order + 1))
    residual = sp.series(expression.subs(psi, series) - target,
                          z, 0, order + 1).removeO().expand()
    solved = {}
    for k, coefficient in enumerate(coefficients, start=1):
        equation = sp.expand(residual.subs(solved)).coeff(z, k)
        solved[coefficient] = sp.solve(equation, coefficient)[0]
    result = sp.expand(series.subs(solved))
    return result, target


def fresh_fixture() -> dict:
    """Run fresh rational implicit-phase and nonzero-b arithmetic controls.

    This fixture uses no saved native packet and is never called by
    :func:`compute`.  SymPy derives the phase inverse and composition rows
    independently before the logarithmic majorants are checked.
    """
    import sympy as sp

    c = MPIntervalContext()
    c.prec = 260
    z, psi = sp.symbols("z psi")
    psi0 = Fraction(1, 3)
    # The coefficient of psi is two at the base point, so this is a genuinely
    # nonunit implicit phase while all fixed-target derivatives remain rational.
    phase_expression = (
        (2 + z + z**2 / 2 + z**3 / 6) * psi
        + z * psi**2 / 10
        + z**2 * psi**3 / 20
    )
    phase = LogZPhase.from_mapping(
        c, _exact_phase_mapping(c, phase_expression, z, psi, psi0)
    )
    inverse, delta, proof = inverse_coefficients(phase_engine := LogZPhase(c),
                                                  phase, c.ln(c.mpf(2)))
    inverse_series, target_value = _fixture_series(
        phase_expression, z, psi, sp.Rational(psi0.numerator, psi0.denominator)
    )
    inverse_derivative_checks = {
        order: _contains_log(
            c,
            inverse[order],
            sp.diff(inverse_series, z, order).subs(z, 0),
        )
        for order in (1, 2, 3)
    }
    delta_checks = {
        0: delta[0] is None,
        **{
            order: _contains_log(
                c,
                delta[order],
                sp.diff(inverse_series, z, order).subs(z, 0)
                / math.factorial(order),
            )
            for order in (1, 2, 3)
        },
    }

    target_expression = psi**2 + z * psi + z**2 / 3 + z**3 / 7
    target = LogZPhase.from_mapping(
        c, _exact_phase_mapping(c, target_expression, z, psi, psi0)
    )
    composed = compose_inverse(phase_engine, target, delta)
    composed_series = sp.expand(target_expression.subs(psi, inverse_series))
    composition_checks = {
        order: _contains_log(
            c,
            composed[order],
            sp.diff(composed_series, z, order).subs(z, 0),
        )
        for order in ORDINARY_Z_ORDERS
    }

    # An algebra-isolation fixture exercises the restoration b term without
    # reading any source packet.  Here r=alpha=t0=0 exactly, a=1,
    # K=1/(2*pi), and E_relative=1.  This is intentionally independent
    # arithmetic: the source identity t0=-b/a is not asserted.  With
    # T1hat=that=0, both B_over_E and Bnorm_phi reduce to b/2.
    one = c.mpf(0)
    exact_zero = LogZ(c).zero()
    isolated_caps = {
        "a": LogZ.constant(LogZ(c), one),
        "b": LogZ.constant(LogZ(c), one),
        "t0": exact_zero,
        "r": exact_zero,
        "alpha": exact_zero,
        "K": LogZ.constant(LogZ(c), -c.ln(c.mpf(2) * c.pi)),
    }
    nonzero = _compute_chart(
        c,
        isolated_caps,
        LogZ.constant(LogZ(c), one),
        {"log_D_lower": one,
         "log_lambda_lower": -c.ln(c.mpf(2) * c.pi)},
    )
    expected_half_log = -c.ln(c.mpf(2))

    def _contains_interval(actual, expected):
        return (actual is not None and actual.a <= expected.a
                and actual.b >= expected.b)

    b_over_e_rows = nonzero["B_over_E"]
    bnorm_phi_rows = nonzero["Bnorm_phi"]
    nonzero_b = _contains_interval(b_over_e_rows[0], expected_half_log)
    nonzero_bphi = _contains_interval(bnorm_phi_rows[0], expected_half_log)
    b_over_e_derivatives_zero = all(
        b_over_e_rows[order] is None for order in (1, 2, 3)
    )
    bnorm_phi_derivatives_zero = all(
        bnorm_phi_rows[order] is None for order in (1, 2, 3)
    )
    checks = {
        "nonunit_phase_lambda": phase[0, 1].b > c.mpf(0),
        "inverse_derivatives_through_Z3": all(inverse_derivative_checks.values()),
        "normalized_inverse_through_Z3": all(delta_checks.values()),
        "composed_target_through_Z3": all(composition_checks.values()),
        "nonzero_b_term_retained": nonzero_b and nonzero_bphi,
        "isolated_B_over_E_equals_b_over_2": nonzero_b,
        "isolated_Bnorm_phi_equals_b_over_2": nonzero_bphi,
        "isolated_B_over_E_derivatives_zero": b_over_e_derivatives_zero,
        "isolated_Bnorm_phi_derivatives_zero": bnorm_phi_derivatives_zero,
        "inverse_proof_passed": proof["passed"] is True,
    }
    _require(all(checks.values()), "Fresh directed-loop fixture failed")
    result = {
        "passed": True,
        "checks": checks,
        "phase_expression": str(phase_expression),
        "fixed_phase_point": str(psi0),
        "implicit_target_value": str(target_value),
        "inverse_series": str(inverse_series),
        "inverse_derivative_checks": inverse_derivative_checks,
        "normalized_inverse_checks": delta_checks,
        "composition_checks": composition_checks,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "phase_rows": len(LogZPhase.indices),
        "independent_sympy_derivatives": True,
        "nonzero_b_fixture_scope": (
            "arithmetic isolation with r=alpha=t0=0 exactly, a=1, "
            "K=1/(2*pi), E_relative=1; t0=-b/a is not asserted"
        ),
        "nonzero_b_expected_log_half": _record(expected_half_log),
        "nonzero_b_rows": {
            "B_over_E": _record(b_over_e_rows[0]),
            "Bnorm_phi": _record(bnorm_phi_rows[0]),
        },
        "stored_logs_exponentiated": False,
        "old_loop_fixture_called": False,
    }
    json.dumps(result, allow_nan=False)
    return result


fixture_bound_check = fresh_fixture


__all__ = (
    "BASE",
    "MODULE",
    "C5_SOURCE",
    "C5_SOURCE_SHA256",
    "LOG_Z_SOURCE",
    "LOG_Z_SOURCE_SHA256",
    "CHARTS",
    "ORDINARY_Z_ORDERS",
    "compute",
    "fresh_fixture",
    "fixture_bound_check",
)
