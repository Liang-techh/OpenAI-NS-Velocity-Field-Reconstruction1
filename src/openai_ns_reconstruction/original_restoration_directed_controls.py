"""Corrected native controls for the restoration and post-restore charts.

The module is intentionally a consumer of source-bound data.  It accepts the
saved parent report, a shared corrected-input packet, and the new long/reference
packet.  It does not import or execute an accepted native producer, owner,
builder, checker, recurrence, quadrature, or self-check.  Signed source rows
are read again and converted with the tuple-safe ``logs`` primitive supplied by
``original_native_directed_common``.  The repaired local inverse-Utheta rows in
the parent report are read as numeric inputs; their producer is not replayed.

All returned rows are raw interval objects.  The coordinator serializes them
with the common ``pack`` helper after the complete interval arithmetic has
finished.  A logarithmic cap is never promoted to a signed function value.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any

from mpmath import make_mpf
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, require
from . import original_native_directed_common as common
from . import original_native_log_jet4 as jet4


MODULE = "src/openai_ns_reconstruction/original_restoration_directed_controls.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_control_directed_repair.json"
)
PARENT_REPORT_SHA256 = (
    "2ba16bd7c1bfd2fae90a5f13c69736d2764cccc049f8737e653fd10f9f8dbf3c"
)

CHARTS = ("restoration", "post_restore")
ORDERS = tuple(range(5))
NAMES = ("phi", "U", "M", "A", "B", "C")


def _path(root: Mapping[str, Any], dotted: str):
    """Read a dotted handle, allowing a literal key containing a colon."""
    value: Any = root
    for part in dotted.split("."):
        require(isinstance(value, Mapping) and part in value,
                "Missing source handle: " + dotted)
        value = value[part]
    return value


def _first(root: Mapping[str, Any], paths: Sequence[str], label: str):
    for dotted in paths:
        value: Any = root
        ok = True
        for part in dotted.split("."):
            if isinstance(value, Mapping) and part in value:
                value = value[part]
            else:
                ok = False
                break
        if ok:
            return value
    raise ValueError("Missing " + label)


def _interval(context: MPIntervalContext, value: Any):
    """Coerce a finite exact scalar/record/interval into ``context``."""
    if hasattr(value, "_mpi_"):
        rows = value._mpi_
        if len(rows) == 4 and all(type(part) is int for part in rows):
            rows = (rows, rows)
        require(len(rows) == 2, "Malformed interval input")
        return context.mpf([make_mpf(tuple(rows[0])), make_mpf(tuple(rows[1]))])
    if isinstance(value, Mapping):
        if "value" in value and "lower_exact_mpf_tuple" not in value:
            return _interval(context, value["value"])
        if "signed_interval" in value:
            return _interval(context, value["signed_interval"])
        if "lower_exact_mpf_tuple" in value:
            return read(context, value)
        raise ValueError("Malformed interval mapping")
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("Binary float interval inputs are not accepted")
    if isinstance(value, Fraction):
        return context.mpf(value.numerator) / value.denominator
    if type(value) is int or isinstance(value, str):
        return context.mpf(value)
    if isinstance(value, (list, tuple)) and len(value) == 2:
        left = _interval(context, value[0])
        right = _interval(context, value[1])
        return context.mpf([common.lower(left), common.upper(right)])
    raise TypeError("Unsupported interval input")


def _raw5(raw: Any):
    """Unwrap the row container used by worker packets."""
    if isinstance(raw, Mapping):
        for key in ("raw5", "raw_rows", "rows", "values", "native_caps", "log_caps"):
            if key in raw:
                return _raw5(raw[key])
    require(isinstance(raw, (list, tuple)) and len(raw) == 5,
            "Five ordinary-Z rows are required")
    return list(raw)


def _raw_rows(raw: Any, count: int):
    """Unwrap a worker row container with an explicitly requested length."""
    if isinstance(raw, Mapping):
        for key in ("raw5", "raw_rows", "rows", "values", "native_caps", "log_caps"):
            if key in raw:
                return _raw_rows(raw[key], count)
    require(isinstance(raw, (list, tuple)) and len(raw) == count,
            "Exactly " + str(count) + " ordinary-Z rows are required")
    return list(raw)


def _signed5(context: MPIntervalContext, raw: Any):
    """Decode five signed source rows without interpreting caps as values."""
    rows = _raw5(raw)
    out = []
    for order, row in enumerate(rows):
        if isinstance(row, Mapping):
            actual_order = row.get("ordinary_Z_order", row.get("Z_order", order))
            require(actual_order == order, "Signed source rows are misordered")
            if row.get("exact_zero") is True:
                out.append(context.mpf(0))
                continue
            if "value" in row:
                row = row["value"]
            elif "signed_interval" in row:
                row = row["signed_interval"]
            elif "raw" in row:
                row = row["raw"]
        out.append(_interval(context, row))
    return out


def _is_cap_row(row: Any):
    return isinstance(row, Mapping) and (
        "log_absolute_upper_bound" in row or "log_cap" in row
    )


def _cap5(context: MPIntervalContext, raw: Any):
    """Decode five already-repaired logarithmic rows."""
    rows = _raw5(raw)
    out = []
    for order, row in enumerate(rows):
        require(isinstance(row, Mapping), "Packed logarithmic row is required")
        require(row.get("ordinary_Z_order", row.get("Z_order", order)) == order,
                "Packed logarithmic rows are misordered")
        if row.get("exact_zero") is True:
            out.append(None)
            continue
        value = row.get("log_absolute_upper_bound", row.get("log_cap"))
        require(value is not None, "Nonzero logarithmic row is missing its cap")
        out.append(_interval(context, value))
    return out


def _raw_log_caps(context: MPIntervalContext, raw: Any, label: str):
    """Read explicit worker log-cap rows without applying another logarithm.

    The new long worker owns this boundary: ``native_caps[name]`` is a raw
    five-row sequence of logarithmic intervals or ``None`` for an exact zero.
    A signed ``value`` row is rejected rather than being guessed as a cap.
    """
    rows = _raw5(raw)
    require(len(rows) == 5, label + " must contain exactly five raw cap rows")
    out = []
    for order, row in enumerate(rows):
        if row is None:
            out.append(None)
            continue
        if isinstance(row, Mapping):
            actual_order = row.get("ordinary_Z_order", row.get("Z_order", order))
            require(actual_order == order, label + " rows are misordered")
            if row.get("exact_zero") is True:
                out.append(None)
                continue
            if "log_absolute_upper_bound" in row:
                row = row["log_absolute_upper_bound"]
            elif "log_cap" in row:
                row = row["log_cap"]
            elif "value" in row or "signed_interval" in row:
                raise ValueError(label + " contains a signed row, not a raw log cap")
        out.append(_interval(context, row))
    return out


def _shared(context: MPIntervalContext, shared: Mapping[str, Any], name: str):
    require(isinstance(shared, Mapping), "Corrected shared input mapping is required")
    require(name in shared, "Missing corrected shared input: " + name)
    return shared[name]


def _shared_scalar(context: MPIntervalContext, shared: Mapping[str, Any], name: str):
    value = _shared(context, shared, name)
    return _interval(context, value)


def _shared_rows(context: MPIntervalContext, shared: Mapping[str, Any], name: str,
                 count: int = 5, caps: bool = True):
    raw = _shared(context, shared, name)
    rows = _raw5(raw) if count == 5 else _raw_rows(raw, count)
    require(len(rows) == count, name + " must contain " + str(count) + " rows")
    if caps and count == 5 and all(_is_cap_row(row) for row in rows):
        return _cap5(context, rows)
    if count == 4 and all(_is_cap_row(row) for row in rows):
        return [
            None if row.get("exact_zero") is True else
            _interval(context, row.get("log_absolute_upper_bound", row.get("log_cap")))
            for row in rows
        ]
    return [_interval(context, row) if row is not None else None for row in rows]


def _parent_root(parent: Mapping[str, Any]):
    if "actual_native_control_directed_repair" in parent:
        return parent
    # A focused caller may pass the full section alongside its source inputs.
    require(isinstance(parent, Mapping), "Native directed parent mapping is required")
    return parent


def _source_rows(context: MPIntervalContext, parent: Mapping[str, Any], name: str):
    if name == "U_R110":
        raw = _first(parent, (
            "actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U",
            "actual_U_R110", "sources.U_R110",
        ), name)
    else:
        raw = _first(parent, (
            "actual_long_reference_background_I1_numeric.charts.reference.endpoint_M",
            "actual_M_Rz", "sources.M_Rz",
        ), name)
    return _signed5(context, raw)


def _inverse_utheta(context: MPIntervalContext, parent: Mapping[str, Any], chart: str):
    section = parent.get("actual_native_control_directed_repair")
    if not isinstance(section, Mapping):
        section = parent if isinstance(parent.get("charts"), Mapping) else None
    require(isinstance(section, Mapping), "Accepted directed native repair section is missing")
    if chart == "restoration":
        key = "restoration:restoration"
    else:
        key = "restoration:post_restore"
    charts = section.get("charts", {}) if isinstance(section, Mapping) else {}
    local = charts.get(key) if isinstance(charts, Mapping) else None
    require(isinstance(local, Mapping), "Accepted inverse-Utheta chart is missing")
    block = local.get("directed_inverse_Utheta_block")
    require(isinstance(block, Mapping), "Accepted inverse-Utheta block is missing")
    rows = block.get("directed_inverse_Utheta_rows")
    return _cap5(context, rows)


def _reference_caps(context: MPIntervalContext, long_packet: Mapping[str, Any]):
    """Read corrected reference A/B/C memories from the new worker packet."""
    require(isinstance(long_packet, Mapping), "New long/reference packet is required")
    planned = None
    for key in ("planned_charts", "plannedcharts", "charts"):
        candidate = long_packet.get(key)
        if isinstance(candidate, Mapping) and isinstance(candidate.get("reference"), Mapping):
            planned = candidate["reference"]
            break
    if planned is None and isinstance(long_packet.get("reference"), Mapping):
        planned = long_packet["reference"]
    require(isinstance(planned, Mapping), "New planned reference chart is missing")
    caps = planned.get("native_caps", planned.get("native_log_absolute_caps"))
    require(isinstance(caps, Mapping), "New reference native caps are missing")
    out = {}
    for name in ("A", "B", "C"):
        require(name in caps, "New reference cap is missing: " + name)
        out[name] = _raw_log_caps(context, caps[name],
                                  "New reference native cap " + name)
    return out


def _log_add_rows(context: MPIntervalContext, *rows):
    require(rows, "At least one logarithmic row jet is required")
    require(all(len(row) == 5 for row in rows), "Five logarithmic rows are required")
    return [common.logsum(context, *(row[k] for row in rows)) for k in ORDERS]


def _point_rows(context: MPIntervalContext, rows):
    return [common.point(context, common.upper(row)) if row is not None else None
            for row in rows]


def _physical_scales(context: MPIntervalContext, parent: Mapping[str, Any],
                    shared: Mapping[str, Any], chart: str):
    """Build physical logarithmic scales before extracting any endpoint."""
    log_p = _interval(context, shared["logP"] if "logP" in shared else
                      _first(parent, (
                          "actual_long_reference_background_I1_numeric.scalar_graph.logP",
                          "logP",
                      ), "logP"))
    log_ref = _interval(context, shared["logref"] if "logref" in shared else
                        _first(parent, (
                            "actual_long_reference_background_I1_numeric.scalar_graph.logref",
                            "logref",
                        ), "logref"))
    log_f0 = _shared_scalar(context, shared, "logF0")
    log_pstar = _shared_scalar(context, shared, "logPstar")
    q_values = _shared_rows(context, shared, "q_log_values5")
    require(q_values[0] is not None, "The signed log(1+Z^2) value row is missing")
    offsets = context.mpf(["-8", "-7"] if chart == "restoration" else ["-7", "-6"])
    log_utheta = log_p + offsets / context.mpf(10) - q_values[0]
    log_r = log_ref + context.ln(context.mpf(110)) + offsets
    log_phi = log_utheta - log_f0 - (
        context.ln(context.mpf(2)) + log_r
    ) / context.mpf(2)
    log_e = log_utheta - log_pstar
    log_w = context.mpf(2) * log_r
    return dict(
        logP=log_p,
        logref=log_ref,
        offsets=offsets,
        log_Utheta=log_utheta,
        log_Utheta_lower=common.point(context, common.lower(log_utheta)),
        log_Utheta_upper=common.point(context, common.upper(log_utheta)),
        log_R=log_r,
        log_R_lower=common.point(context, common.lower(log_r)),
        log_R_upper=common.point(context, common.upper(log_r)),
        log_phi=log_phi,
        log_phi_lower=common.point(context, common.lower(log_phi)),
        log_phi_upper=common.point(context, common.upper(log_phi)),
        log_E=log_e,
        log_E_lower=common.point(context, common.lower(log_e)),
        log_E_upper=common.point(context, common.upper(log_e)),
        log_W=log_w,
        log_W_upper=common.point(context, common.upper(log_w)),
    )


def _relative(context: MPIntervalContext, shared: Mapping[str, Any]):
    gradient = _shared_rows(context, shared, "gradient_logs4", count=4)
    q_derivative = _shared_rows(context, shared, "q_derivative_logs4", count=4)
    if "lambda_gradient_logs4" in shared:
        lambda_gradient = _shared_rows(context, shared, "lambda_gradient_logs4", count=4)
    else:
        log_lambda = _shared_scalar(context, shared, "logLambda")
        lambda_gradient = [None if value is None else value + log_lambda
                           for value in gradient]
    # Subtraction of two logarithmic derivative expressions is bounded by the
    # positive log sum of their absolute caps.  Bell majorants then remain
    # entirely in the new generic arithmetic path.
    phi_logderivs = [common.logsum(context, left, right)
                     for left, right in zip(lambda_gradient, q_derivative)]
    phi_relative = jet4.exp_relative(context, phi_logderivs)
    e_relative = jet4.exp_relative(context, q_derivative)
    return dict(
        gradient_logs4=gradient,
        lambda_gradient_logs4=lambda_gradient,
        q_derivative_logs4=q_derivative,
        phi_logderivs4=phi_logderivs,
        phi_relative=phi_relative,
        E_relative=e_relative,
    )


def _sigma_log(context: MPIntervalContext, parent: Mapping[str, Any],
               shared: Mapping[str, Any]):
    if "sigma_prime_upper" in shared:
        bound = _interval(context, shared["sigma_prime_upper"])
    else:
        bound = _interval(context, _first(parent, (
            "actual_long_reference_native_controls.slope_proof.sigma_prime_upper",
            "actual_restoration_native_controls.source_formula_inputs.sigma_prime_source.used_for_C6_restoration",
        ), "actual sigma-prime bound"))
    require(common.lower(bound) > 0, "The actual sigma-prime bound is not positive")
    return context.ln(bound), bound


def _native_chart(context: MPIntervalContext, name: str, z,
                  u_raw, m_raw, phi_logs, inverse_utheta,
                  incoming, sigma_log=None, incoming_m=None):
    four_z = [4 * z, context.mpf(4), context.mpf(0), context.mpf(0), context.mpf(0)]
    four_z_logs = common.logs(context, four_z)
    if name == "restoration":
        E_raw = [left - right for left, right in zip(u_raw, four_z)]
        K_raw = [left - right for left, right in zip(m_raw, u_raw)]
        E_logs = common.logs(context, E_raw)
        K_logs = common.logs(context, K_raw)
        U_logs = _log_add_rows(context, four_z_logs, E_logs)
        M_logs = _log_add_rows(context, four_z_logs, E_logs, K_logs)
        require(sigma_log is not None, "Restoration sigma-prime bound is required")
        Uy_logs = jet4.log_scale(context, E_logs, sigma_log)
        memory_source = "actual_M_Rz"
        raw_sources = {"U_R110": u_raw, "M_Rz": m_raw,
                       "four_Z": four_z, "E": E_raw, "K": K_raw}
    else:
        require(incoming_m is not None, "Post-restore incoming M memory is required")
        K_restore_logs = _log_add_rows(context, incoming_m, four_z_logs)
        U_logs = list(four_z_logs)
        M_logs = _log_add_rows(context, four_z_logs, K_restore_logs)
        Uy_logs = [None] * 5
        E_logs = [None] * 5
        K_logs = K_restore_logs
        E_raw = None
        K_raw = None
        memory_source = "restoration.M_caps"
        raw_sources = {"four_Z": four_z}
    density = {
        "A": jet4.log_product(context, U_logs, U_logs),
        "B": jet4.log_product(context, phi_logs, phi_logs),
        "C": jet4.log_product(context, phi_logs, phi_logs),
    }
    histories = {
        "A": common.history(context, incoming["A"], density["A"], 1),
        "B": common.history(context, incoming["B"], density["B"], 2),
        "C": common.history(context, incoming["C"], density["C"], 1),
    }
    native = {
        "phi": list(phi_logs), "U": U_logs, "M": M_logs,
        "A": histories["A"], "B": histories["B"], "C": histories["C"],
    }
    b_logs = jet4.log_scale(
        context,
        jet4.log_product(context, Uy_logs, inverse_utheta),
        context.ln(context.mpf(2)),
    )
    return native, dict(
        raw_sources=raw_sources,
        E_logs=E_logs,
        K_logs=K_logs,
        density=density,
        histories=histories,
        M_memory=M_logs,
        incoming_memory_source=memory_source,
        incoming_caps_are_bounds=True,
        cap_value_not_promoted=True,
        positive_Duhamel_rates={"A": 1, "B": 2, "C": 1},
        b_logs=b_logs,
        b_formula="b=2*U_y/Utheta",
        b_log_formula="log|b|=log(2)+log|U_y|-log(Utheta)",
        physical_U_y=True,
        source_executed=False,
    )


def _root_controls(context: MPIntervalContext, b_logs, p2, e_relative,
                   logeta, logdstar):
    # ``a`` is the exact source constant 4/5.  Its reciprocal is the exact
    # context value 5/4, so no legacy root recurrence is consumed.
    a_value = context.mpf(4) / context.mpf(5)
    inverse_a_value = context.mpf(5) / context.mpf(4)
    log_a = context.ln(a_value)
    log_inverse_a = context.ln(inverse_a_value)
    a_logs = [log_a, None, None, None, None]
    inverse_a_logs = [log_inverse_a, None, None, None, None]
    t0_logs = jet4.log_product(context, b_logs, inverse_a_logs)
    b_squared_logs = jet4.log_product(context, b_logs, b_logs)
    b_squared_over_a = jet4.log_product(context, b_squared_logs, inverse_a_logs)
    kappa_logs = [common.logsum(context, a_logs[k], b_squared_over_a[k])
                  for k in ORDERS]
    delta_logs = [common.logsum(
        context, a_logs[k], b_squared_over_a[k],
        context.ln(context.mpf(2)) if k == 0 else None,
    ) for k in ORDERS]
    if b_logs[0] is None:
        bmargin = dict(
            passed=True, exact_zero=True,
            log_abs_b_upper=None,
            log_abs_b_minus_log_0_1=None,
            bound="b=0",
        )
        kappa_upper_log = common.upper(log_a)
        kappa_upper_value = a_value
        rational_sharpening = "a=4/5 because post_restore U_y is exactly zero"
    else:
        margin_interval = b_logs[0] - context.ln(context.mpf("0.1"))
        margin_upper = common.upper(margin_interval)
        require(margin_upper < 0, "Actual restoration bound did not prove |b|<1/10")
        bmargin = dict(
            passed=True, exact_zero=False,
            log_abs_b_upper=common.point(context, common.upper(b_logs[0])),
            log_abs_b_minus_log_0_1=common.point(context, margin_upper),
            bound="|b|<1/10",
        )
        # From the actual threshold and a=4/5, b^2/a < 1/80 and hence
        # kappa < 13/16.  This is a directed rational consequence, not a
        # selected cap value.
        kappa_upper_log = common.upper(context.ln(context.mpf(13) / context.mpf(16)))
        kappa_upper_value = context.mpf(13) / context.mpf(16)
        rational_sharpening = "|b|<1/10 and a=4/5 imply kappa<13/16"
        kappa_logs[0] = common.point(context, kappa_upper_log)
    require(kappa_upper_log < 0, "Corrected kappa bound is not below one")
    # The active source branch has gamma=2*eta-Delta.  Delta is negative,
    # so its value is bounded by the positive absolute-log row already built
    # above.  Keep the value row as a cap; the derivative rows are the exact
    # constant-eta relation gamma_k=-kappa_k for k>=1.
    gamma_upper_log = common.logsum(
        context, context.ln(context.mpf(2)) + logeta, delta_logs[0]
    )
    gamma_logs = [gamma_upper_log] + list(kappa_logs[1:5])
    gamma_lower_value = context.mpf(2) - kappa_upper_value
    require(common.lower(gamma_lower_value) > 0,
            "Corrected gamma lower bound is not positive")
    q2_lower_log = (
        context.ln(gamma_lower_value)
        - context.ln(context.mpf(2) * a_value)
    )
    q_lower_log = q2_lower_log / context.mpf(2)

    # q=sqrt(q^2) through ordinary Z order 3.  The lower q-value is used
    # only in denominators; all cap rows remain logarithmic magnitudes.
    alg = common.LogZ(context)
    q2_logs = jet4.log_scale(
        context, gamma_logs,
        context.ln(context.mpf(5) / context.mpf(8)),
    )
    s = list(q2_logs[:4])
    add = lambda *terms: common.logsum(context, *terms)
    q_logs = [
        s[0] / context.mpf(2),
        None if s[1] is None else
        s[1] - context.ln(2) - q_lower_log,
        add(
            None if s[2] is None else
            s[2] - context.ln(2) - q_lower_log,
            None if s[1] is None else
            2 * s[1] - context.ln(4) - 3 * q_lower_log,
        ),
        add(
            None if s[3] is None else
            s[3] - context.ln(2) - q_lower_log,
            None if s[1] is None or s[2] is None else
            s[1] + s[2] + context.ln(3) - context.ln(4)
            - 3 * q_lower_log,
            None if s[1] is None else
            3 * s[1] + context.ln(3) - context.ln(8)
            - 5 * q_lower_log,
        ),
    ]
    p2_jet = alg.__class__(context, list(p2[:4]))
    q_jet = alg.__class__(context, q_logs)
    u_jet = alg.scale(
        alg.mul(p2_jet, q_jet),
        -common.point(context, common.lower(logdstar)),
    )
    u_logs = [u_jet[k] for k in range(4)]
    v_jet = alg.mul(u_jet, u_jet)
    v_logs = [v_jet[k] for k in range(4)]
    h_logs = [
        context.mpf(0),
        None if v_logs[1] is None else v_logs[1] - context.ln(2),
        add(
            None if v_logs[1] is None else
            2 * v_logs[1] + context.ln("0.75"),
            None if v_logs[2] is None else
            v_logs[2] - context.ln(2),
        ),
        add(
            None if v_logs[1] is None else
            3 * v_logs[1] + context.ln("1.875"),
            None if v_logs[1] is None or v_logs[2] is None else
            v_logs[1] + v_logs[2] + context.ln("2.25"),
            None if v_logs[3] is None else
            v_logs[3] - context.ln(2),
        ),
    ]
    h_jet = alg.__class__(context, h_logs)
    r_jet = alg.mul(u_jet, h_jet)
    r_logs = [r_jet[k] for k in range(4)]
    # |r|=|u|/sqrt(1+u^2)<1.  The value row is therefore the sharp log cap 0.
    r_logs[0] = context.mpf(0)
    alpha_jet = alg.scale(
        alg.mul(q_jet, h_jet), context.ln(2)
    )
    alpha_logs = [alpha_jet[k] for k in range(4)]
    t0_jet = alg.__class__(context, list(t0_logs[:4]))
    q2_jet = alg.__class__(context, list(q2_logs[:4]))
    denominator = alg.add(
        alg.constant(context.mpf(0)),
        alg.mul(t0_jet, t0_jet),
        alg.scale(q2_jet, context.ln(2)),
    )
    K_jet = alg.scale(
        alg.reciprocal(denominator, context.mpf(0)),
        -context.ln(2 * context.pi),
    )
    K_logs = [K_jet[k] for k in range(4)]
    # D=(cos(psi)-r)^2+sin(psi)^2 >= h^4/4.  Use an upper cap for |u|
    # before the directed lower endpoint is extracted from log(h).
    if u_logs[0] is None:
        h_log_upper = context.mpf(0)
    else:
        u0_upper = common.upper(u_logs[0])
        h_log_upper = common.logsum(
            context, context.mpf(0),
            2 * common.point(context, u0_upper),
        )
    h_lower_log = -context.mpf("0.5") * h_log_upper
    log_D_lower = (
        4 * common.point(context, common.lower(h_lower_log))
        - context.ln(4)
    )
    log_lambda_lower = (
        -context.ln(2 * context.pi)
        - common.point(context, common.upper(denominator[0]))
    )
    branchproof = dict(
        passed=True,
        kappa_upper_log=common.point(context, kappa_upper_log),
        kappa_upper_lt_one=True,
        Delta_upper_lt_zero=True,
        Delta_formula="Delta=a+b^2/a-2",
        eta_positive=True,
        eta_log_lower=logeta,
        cut_formula="cut=1-Delta/eta",
        cut_lower_gt_one=True,
        sigma_cut_branch="sigma(cut)=1 after the directed Delta<0 and eta>0 checks",
        derivatives_simplified_only_after_Delta_proof=True,
        rational_sharpening=rational_sharpening,
        gamma_formula="gamma=2*eta-Delta=2*eta+2-kappa",
        gamma_positive_lower=gamma_lower_value,
        gamma_ge_one=True,
        gamma_value_row_is_log_cap=True,
        gamma_derivatives_are_negative_kappa_rows=True,
        q_squared_formula="q^2=gamma/(2*a), a=4/5",
        q_lower_log=q_lower_log,
        D_lower_formula="D>=(1-|r|)^2>=h^4/4",
        lambda_lower_formula="lambda=K*(1+t^2)>=K",
    )
    root_caps = {
        "a": a_logs[:4], "b": b_logs[:4], "t0": t0_logs[:4],
        "p2": p2[:4], "Delta": delta_logs[:4], "kappa": kappa_logs[:4],
        "b_squared_over_a": b_squared_over_a[:4],
        "gamma": gamma_logs[:4], "q_squared": q2_logs[:4],
        "q": q_logs, "u": u_logs, "h": h_logs,
        "r": r_logs, "alpha": alpha_logs, "K": K_logs,
    }
    return dict(
        root_caps=root_caps,
        root_inputs=dict(
            a=a_logs, inverse_a=inverse_a_logs, b=b_logs, t0=t0_logs,
            p2=p2, Delta=delta_logs, kappa=kappa_logs,
            b_squared_over_a=b_squared_over_a,
            gamma=gamma_logs, q_squared=q2_logs, q=q_logs,
            u=u_logs, h=h_logs, r=r_logs, alpha=alpha_logs, K=K_logs,
            denominator=[denominator[k] for k in range(4)],
            log_D_lower=log_D_lower,
            log_lambda_lower=log_lambda_lower,
        ),
        conditioning=dict(
            log_D_lower=log_D_lower,
            log_lambda_lower=log_lambda_lower,
            positive_a_lower=a_value,
            positive_gamma_lower=gamma_lower_value,
            positive_q_squared_lower_log=q2_lower_log,
            gamma_ge_one=True,
            positive_dstar_source=True,
            eta_positive_source=True,
            eta_le_one=common.upper(logeta) <= 0,
            Delta_negative=True,
            denominator_positive_lower=context.mpf(1),
            positive_D_lower=True,
            positive_lambda_lower=True,
            D_positive_lower_formula="h^4/4",
            lambda_positive_lower_formula="K",
            no_cap_values_selected=True,
        ),
        bmargin=bmargin,
        branchproof=branchproof,
        E_relative=e_relative,
        a_value_interval=a_value,
        inverse_a_value_interval=inverse_a_value,
        eta_log=logeta,
        dstar_log=logdstar,
        ordinary_Z_orders=list(ORDERS),
        p2_ordinary_Z_orders=[0, 1, 2, 3],
    )


def _chart_output(context: MPIntervalContext, chart: str, parent,
                  shared, long_packet, relative, scales, incoming,
                  z, delta, F2, P0, inverseL, u_raw, m_raw,
                  sigma_log, sigma_bound):
    inverse_utheta = _inverse_utheta(context, parent, chart)
    if chart == "restoration":
        native, detail = _native_chart(
            context, chart, z, u_raw, m_raw,
            jet4.log_scale(context, relative["phi_relative"],
                           common.point(context, common.upper(scales["log_phi"]))),
            inverse_utheta, incoming, sigma_log=sigma_log,
        )
    else:
        native, detail = _native_chart(
            context, chart, z, u_raw, m_raw,
            jet4.log_scale(context, relative["phi_relative"],
                           common.point(context, common.upper(scales["log_phi"]))),
            inverse_utheta, incoming,
            incoming_m=incoming["M"],
        )
    p2, p2_detail = common.p2(
        context, native, F2, P0, z, delta,
        scales["log_R_upper"], inverse_utheta, inverseL,
    )
    roots = _root_controls(
        context, detail["b_logs"], p2, relative["E_relative"],
        _shared_scalar(context, shared, "logeta"),
        _shared_scalar(context, shared, "logdstar"),
    )
    phi_logs = native["phi"]
    chart_key = chart
    return dict(
        passed=True,
        chart=chart_key,
        coordinate="t" if chart == "restoration" else "q",
        domain=["0", "1"],
        ordinary_Z_orders=list(ORDERS),
        p2_ordinary_Z_orders=[0, 1, 2, 3],
        native_caps={name: common.pack(native[name]) for name in NAMES},
        raw5={name: common.pack(native[name]) for name in NAMES},
        phi_relative=common.pack(relative["phi_relative"]),
        E_relative=common.pack(relative["E_relative"]),
        inverse_Utheta=common.pack(inverse_utheta),
        p2=common.pack(p2),
        p2_detail=p2_detail,
        root_caps={name: common.pack(rows) for name, rows in roots["root_caps"].items()},
        root_inputs=roots["root_inputs"],
        conditioning=roots["conditioning"],
        bmargin=roots["bmargin"],
        branchproof=roots["branchproof"],
        scales=dict(
            log_Utheta=common.pack([scales["log_Utheta"]] + [None] * 4),
            log_Utheta_lower=scales["log_Utheta_lower"],
            log_Utheta_upper=scales["log_Utheta_upper"],
            log_R=scales["log_R"],
            log_R_lower=scales["log_R_lower"],
            log_R_upper=scales["log_R_upper"],
            log_phi=scales["log_phi"],
            log_phi_lower=scales["log_phi_lower"],
            log_phi_upper=scales["log_phi_upper"],
            log_E=scales["log_E"],
            log_E_lower=scales["log_E_lower"],
            log_E_upper=scales["log_E_upper"],
            log_W=scales["log_W"],
            log_W_upper=scales["log_W_upper"],
            offsets=scales["offsets"],
            physical_Utheta_formula="logP+offset/10-log(1+Z^2)",
            physical_radius_formula="logref+ln(110)+offset",
            context_interval_arithmetic_before_endpoint_extraction=True,
            radius_materialized=False,
        ),
        source_detail=detail,
        sigma_prime_upper=sigma_bound,
        incoming_memory={name: common.pack(incoming[name])
                         for name in ("A", "B", "C")},
        source_inputs=dict(
            actual_U_R110="actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U",
            actual_M_Rz="actual_long_reference_background_I1_numeric.charts.reference.endpoint_M",
            corrected_reference_memory="long_packet.planned_charts.reference.native_caps",
            post_incoming_M=("restoration.M_caps" if chart == "post_restore" else None),
            signed_value_selected=False,
            cap_value_not_promoted=True,
        ),
        source_execution=dict(
            source_executed=False,
            accepted_ancestors_replayed=False,
            accepted_builders_or_owners_replayed=False,
            accepted_checkers_replayed=False,
            accepted_recurrences_replayed=False,
            accepted_quadratures_replayed=False,
        ),
    )


def compute(context: MPIntervalContext, parent: Mapping[str, Any],
            shared: Mapping[str, Any], long_packet: Mapping[str, Any]):
    """Compute the corrected restoration/post-restore native cap packet."""
    require(isinstance(context, MPIntervalContext), "Interval context required")
    require(isinstance(parent, Mapping), "Native directed parent mapping required")
    root = _parent_root(parent)
    # The section is only read; accepted producers are never imported or run.
    local = root.get("actual_native_control_directed_repair")
    if isinstance(local, Mapping):
        require(local.get("source_executed") is not True,
                "Accepted directed parent claims source execution")
    z = _shared_scalar(context, shared, "z")
    delta = _shared_scalar(context, shared, "delta")
    require(common.lower(z) >= 0, "Corrected Z interval is outside the source domain")
    require(common.lower(delta) > 0 and common.upper(delta) < 1,
            "Corrected delta interval lost source positivity")
    relative = _relative(context, shared)
    F2 = _shared_rows(context, shared, "F2_logs5")
    P0 = _shared_rows(context, shared, "P0_logs5")
    inverseL = _shared_rows(context, shared, "inverse_L_logs5")
    incoming_reference = _reference_caps(context, long_packet)
    u_raw = _source_rows(context, root, "U_R110")
    m_raw = _source_rows(context, root, "M_Rz")
    sigma_log, sigma_bound = _sigma_log(context, root, shared)
    incoming = dict(incoming_reference)
    outputs = {}
    for chart in CHARTS:
        scales = _physical_scales(context, root, shared, chart)
        chart_incoming = dict(incoming)
        if chart == "post_restore":
            restoration = outputs["restoration"]
            chart_incoming["M"] = restoration["native_caps"]["M"]
            # The incoming M is retained as a cap only.  Decode it back to raw
            # intervals for the positive post-restore triangle bound.
            chart_incoming["M"] = _cap5(context, chart_incoming["M"])
        outputs[chart] = _chart_output(
            context, chart, root, shared, long_packet, relative, scales,
            chart_incoming, z, delta, F2, P0, inverseL, u_raw, m_raw,
            sigma_log, sigma_bound,
        )
        # Keep the actual raw interval histories for post_restore rather than
        # decoding serialized rows.  This is the only inter-chart memory edge.
        if chart == "restoration":
            restoration_native = outputs[chart]["source_detail"]["histories"]
            incoming = dict(restoration_native)
    return dict(
        passed=True,
        stage="original_restoration_post_restore_directed_controls",
        schema="original_restoration_directed_controls_v1",
        parent_report=PARENT_REPORT,
        parent_report_sha256=PARENT_REPORT_SHA256,
        charts=outputs,
        shared_inputs=dict(
            z=z, delta=delta,
            logLambda=_shared_scalar(context, shared, "logLambda"),
            logF0=_shared_scalar(context, shared, "logF0"),
            logPstar=_shared_scalar(context, shared, "logPstar"),
            logN=_shared_scalar(context, shared, "logN"),
            logeta=_shared_scalar(context, shared, "logeta"),
            logdstar=_shared_scalar(context, shared, "logdstar"),
            gradient_logs4=relative["gradient_logs4"],
            lambda_gradient_logs4=relative["lambda_gradient_logs4"],
            F2_logs5=F2, P0_logs5=P0,
            q_log_values5=_shared_rows(context, shared, "q_log_values5"),
            q_derivative_logs4=relative["q_derivative_logs4"],
            inverse_L_logs5=inverseL,
        ),
        source_inputs=dict(
            actual_U_R110=True,
            actual_M_Rz=True,
            corrected_reference_A_B_C_memory=True,
            repaired_inverse_Utheta_read_from_parent=True,
            source_executed=False,
            signed_values_selected=False,
            cap_values_promoted=False,
        ),
        checks=dict(
            passed=True,
            corrected_phi_U_M_A_B_C=True,
            corrected_p2=True,
            corrected_b_t0_Delta_kappa=True,
            b_threshold_rechecked=True,
            cut_branch_rechecked=True,
            physical_interval_scales_built_before_endpoint_extraction=True,
            accepted_local_inverse_Utheta_reused_without_rebuild=True,
            positive_kernel_histories=True,
            downstream_finite_N_recertified=False,
            downstream_full_N_certified=False,
            source_executed=False,
        ),
        limits=dict(
            downstream_finite_N_recertified=False,
            downstream_full_N_certified=False,
            loop_certified=False,
            Rm_certified=False,
            full_goal_complete=False,
        ),
        source_execution=dict(
            source_executed=False,
            accepted_ancestors_replayed=False,
            accepted_builders_or_owners_replayed=False,
            accepted_checkers_replayed=False,
            accepted_recurrences_replayed=False,
            accepted_quadratures_replayed=False,
        ),
    )


def fixture_bound_check():
    """Run only a fresh rational arithmetic fixture.

    This fixture is independent of the saved parent and does not invoke any
    accepted self-check or producer.
    """
    c = MPIntervalContext()
    c.prec = 256
    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    delta = c.mpf(1) / 1000
    q = [1 + z * z, 2 * z, c.mpf(2), c.mpf(0), c.mpf(0)]
    qlog = common.analytic_log_jet(c, q)
    qlogs = common.logs(c, qlog[1:])
    e = jet4.exp_relative(c, qlogs)
    require(len(e) == 5 and e[0] is not None, "Fresh Bell fixture failed")
    l = common.logs(c, [1 - delta * z * z, -2 * delta * z,
                        -2 * delta, c.mpf(0), c.mpf(0)])
    require(l[0] is not None and l[1] is not None, "Fresh L fixture failed")
    return dict(passed=True, fresh_rational=True, bell_rows=5,
                inverse_L_rows=5, historical_inputs_used=False,
                accepted_self_checks_called=False, source_executed=False)


__all__ = ["compute", "fixture_bound_check", "MODULE", "PARENT_REPORT",
           "PARENT_REPORT_SHA256"]
