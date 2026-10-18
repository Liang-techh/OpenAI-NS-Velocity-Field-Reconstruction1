"""Corrected ordinary-``Z`` native controls for the long/reference charts.

This module is the changed-input successor of
``original_long_reference_native_controls``.  It consumes the saved signed
source rows and the raw factors prepared by
``original_native_directed_common``.  The old native packet is used only for
source identity and branch hypotheses; none of its logarithmic caps is used
as a signed field value.

``compute`` deliberately returns interval objects.  The coordinator owns
serialization and the source/hash audit.  No accepted constructor, producer,
owner, checker, recurrence, quadrature, or self-check is imported here.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from mpmath.ctx_iv import MPIntervalContext

from . import original_native_directed_common as common
from .original_core_signed_direction import read, require


MODULE = "src/openai_ns_reconstruction/original_long_reference_directed_controls.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_control_directed_repair.json"
)
PARENT_REPORT_SHA256 = (
    "2ba16bd7c1bfd2fae90a5f13c69736d2764cccc049f8737e653fd10f9f8dbf3c"
)
COMMON_MODULE = "src/openai_ns_reconstruction/original_native_directed_common.py"
ORDERS = tuple(range(5))
CHARTS = ("long_reshape", "reference")
NATIVE_NAMES = ("phi", "U", "M", "A", "B", "C")


def _at(root: Mapping[str, Any], *keys: str) -> Any:
    value: Any = root
    for key in keys:
        require(isinstance(value, Mapping) and key in value,
                "Missing saved directed source handle: " + ".".join(keys))
        value = value[key]
    return value


def _interval(context: MPIntervalContext, value: Any):
    """Read one raw interval without selecting a midpoint.

    The coordinator normally passes raw intervals.  The small mapping branch
    keeps this API usable with an already packed in-memory row while retaining
    exact MPF endpoints.
    """
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value:
            return read(context, value)
        if "value" in value:
            return _interval(context, value["value"])
        if "log_absolute_upper_bound" in value:
            return _interval(context, value["log_absolute_upper_bound"])
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        return context.mpf(value)
    raise TypeError("Expected an interval or exact directed interval record")


def _finite_interval(context: MPIntervalContext, value: Any, label: str) -> None:
    lo, hi = common.lower(value), common.upper(value)
    require(not context.isinf(lo) and not context.isinf(hi)
            and not context.isnan(lo) and not context.isnan(hi),
            label + " must be a finite positive-source logarithm")


def _raw_rows(context: MPIntervalContext, rows: Any, *, label: str,
              signed: bool = False) -> list[Any]:
    require(isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)),
            label + " must be a five-row sequence")
    require(len(rows) == 5, label + " must contain five rows")
    out: list[Any] = []
    for order, row in enumerate(rows):
        if row is None:
            out.append(None)
            continue
        if isinstance(row, Mapping):
            source = row.get("value", row) if signed else row
            if not signed and "log_absolute_upper_bound" in row:
                source = row["log_absolute_upper_bound"]
            out.append(_interval(context, source))
        else:
            out.append(_interval(context, row))
    return out


def _raw_four(context: MPIntervalContext, rows: Any, *, label: str) -> list[Any]:
    require(isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)),
            label + " must be a four-row sequence")
    require(len(rows) == 4, label + " must contain four rows")
    return [_interval(context, row) if row is not None else None for row in rows]


def _signed_source(context: MPIntervalContext, rows: Any, label: str) -> list[Any]:
    """Read saved signed source rows through the common source-row adapter."""
    require(isinstance(rows, list), label + " must be a saved list")
    return common.signed_rows(context, rows)


def _log_source(context: MPIntervalContext, rows: Any, label: str) -> list[Any]:
    """Read raw logarithmic caps supplied by the new common input stage."""
    if (isinstance(rows, Sequence) and not isinstance(rows, (str, bytes))
            and all(isinstance(row, Mapping) and "ordinary_Z_order" in row
                    and "exact_zero" in row for row in rows)):
        # ``pack`` is the common provider's lossless transport form.  Decode
        # exact-zero rows before feeding the raw caps to the native formula.
        # During the raw stage the packed payload is still an interval; after
        # serialization it is an exact directed record.  Accept both forms.
        return [None if row["exact_zero"] else
                _interval(context, row["log_absolute_upper_bound"])
                for row in rows]
    return _raw_rows(context, rows, label=label)


def _max_positive_log(context: MPIntervalContext, value: Any):
    """Return an upper log bound for a positive interval without materialising it."""
    return context.ln(value)


def _log_cap_rows(context: MPIntervalContext, rows: Sequence[Any]) -> list[Any]:
    return list(rows[:5])


def _source_scope(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Check only non-execution/source-identity premises.

    ``passed`` flags from the historical native packet intentionally do not
    appear in the admission conditions.  A source packet marked as executed
    is rejected because this corrective chain is source-data-only.
    """
    sections = {
        "directed_repair": _at(parent, "actual_native_control_directed_repair"),
        "old_native": _at(parent, "actual_long_reference_native_controls"),
        "background": _at(parent, "actual_long_reference_background_I1_numeric"),
        "pressure": _at(parent, "actual_long_reference_pressure_source"),
        "post_switch": _at(parent, "actual_switch_post_background_I1_numeric"),
    }
    for name, section in sections.items():
        if parent.get("source_family") is not None and section.get("source_family") is not None:
            require(section.get("source_family") == parent.get("source_family"),
                    "Source-family handle changed: " + name)
        require(section.get("source_executed") is not True,
                "Historical source execution is not admissible: " + name)
        require(section.get("accepted_ancestors_replayed") is not True,
                "Historical ancestor replay is not admissible: " + name)
    return {
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "source_family": parent.get("source_family"),
        "historical_pass_flags_used_as_certificates": False,
        "saved_source_rows_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "accepted_builders_or_owners_replayed": False,
        "accepted_checkers_or_self_checks_replayed": False,
        "accepted_recurrences_or_quadratures_replayed": False,
        "sections_checked": list(sections),
    }


def _shared(context: MPIntervalContext, shared: Mapping[str, Any]) -> dict[str, Any]:
    required = (
        "z", "delta", "logLambda", "logF0", "logPstar", "logN", "logeta",
        "logdstar", "gradient_logs4", "F2_logs5", "P0_logs5",
        "q_log_values5", "q_derivative_logs4", "inverse_L_logs5",
    )
    for key in required:
        require(key in shared, "Missing corrected common input: " + key)
    out = dict(
        z=_interval(context, shared["z"]),
        delta=_interval(context, shared["delta"]),
        logLambda=_interval(context, shared["logLambda"]),
        logF0=_interval(context, shared["logF0"]),
        logPstar=_interval(context, shared["logPstar"]),
        logN=_interval(context, shared["logN"]),
        logeta=_interval(context, shared["logeta"]),
        logdstar=_interval(context, shared["logdstar"]),
        gradient_logs4=_raw_four(context, shared["gradient_logs4"],
                                 label="gradient_logs4"),
        F2_logs5=_log_source(context, shared["F2_logs5"], "F2_logs5"),
        P0_logs5=_log_source(context, shared["P0_logs5"], "P0_logs5"),
        q_log_values5=_raw_rows(context, shared["q_log_values5"],
                                label="q_log_values5"),
        q_derivative_logs4=_raw_four(context, shared["q_derivative_logs4"],
                                     label="q_derivative_logs4"),
        inverse_L_logs5=_log_source(context, shared["inverse_L_logs5"],
                                    "inverse_L_logs5"),
    )
    require(common.lower(out["z"]) >= common.lower(context.mpf("37/100"))
            and common.upper(out["z"]) <= common.upper(context.mpf("93/250")),
            "Corrected Z interval left the saved native window")
    require(common.upper(out["logeta"]) <= 0,
            "The canonical eta logarithm must establish eta <= 1")
    _finite_interval(context, out["logeta"], "The canonical eta logarithm")
    _finite_interval(context, out["logdstar"], "The canonical d_star logarithm")
    require(common.upper(out["logdstar"]) < 0,
            "The finite canonical d_star logarithm must establish 0 < d_star < 1")
    return out


def _inlet(context: MPIntervalContext, parent: Mapping[str, Any]) -> dict[str, list[Any]]:
    fields = _at(parent, "actual_switch_post_background_I1_numeric",
                 "charts", "post_switch", "native_fields")
    out = {name: _signed_source(context, fields[name], "inlet " + name)
           for name in NATIVE_NAMES}
    require(common.lower(out["phi"][0]) > 0,
            "Saved inlet phi must remain positive")
    require(common.lower(out["U"][0]) > 0,
            "Saved inlet U must remain positive")
    return out


def _log_phi(context: MPIntervalContext, inlet_phi: Sequence[Any]):
    # The closed-form helper is the only new log-jet calculation.  It returns
    # signed derivatives of log(phi); logs() is then applied to those changed
    # inputs to obtain magnitude caps.
    signed_log = common.analytic_log_jet(context, list(inlet_phi))
    caps = common.logs(context, signed_log)
    require(len(signed_log) == len(caps) == 5,
            "Corrected inlet phi log jet must contain five rows")
    return signed_log, caps


def _relative_factors(context: MPIntervalContext, shared: Mapping[str, Any],
                     phi_log_caps: Sequence[Any]):
    lambda_gradient = [
        None if value is None else value + shared["logLambda"]
        for value in shared["gradient_logs4"]
    ]
    qd = shared["q_derivative_logs4"]
    ell_long = [
        common.maxlog(context, (phi_log_caps[k],
                                common.logsum(context, lambda_gradient[k - 1],
                                              qd[k - 1])))
        for k in range(1, 5)
    ]
    ell_reference = [
        common.logsum(context, lambda_gradient[k - 1], qd[k - 1])
        for k in range(1, 5)
    ]
    e_long = [
        common.maxlog(context, (
            common.logsum(context, phi_log_caps[k], lambda_gradient[k - 1]),
            qd[k - 1],
        ))
        for k in range(1, 5)
    ]
    e_reference = list(qd)
    return dict(
        lambda_gradient=lambda_gradient,
        ell_long=ell_long,
        ell_reference=ell_reference,
        e_long=e_long,
        e_reference=e_reference,
        phi_relative_long=common.exp_relative(context, ell_long),
        phi_relative_reference=common.exp_relative(context, ell_reference),
        E_relative_long=common.exp_relative(context, e_long),
        E_relative_reference=common.exp_relative(context, e_reference),
    )


def _safe_B(context: MPIntervalContext, parent: Mapping[str, Any], shared: Mapping[str, Any],
            phi_log_caps: Sequence[Any], q_log_caps: Sequence[Any], factors: Mapping[str, Any]):
    """Rebuild B from the saved gradient bound and changed input logs."""
    amp = _at(parent, "actual_micro_b_p2_mixed_bounds", "whole_tube_amplitude")
    gbar = _interval(context, amp["G_modulus_upper"])
    gbar_log = common.logs(context, [gbar])[0]
    log_220_over_two = context.ln(context.mpf(220)) / 2
    log_220_over_two_cap = common.logs(context, [log_220_over_two])[0]
    B0 = common.logsum(
        context,
        shared["logLambda"] + gbar_log,
        phi_log_caps[0],
        log_220_over_two_cap,
        q_log_caps[0],
    )
    B = [B0]
    for k in range(1, 5):
        B.append(common.logsum(
            context,
            factors["lambda_gradient"][k - 1],
            phi_log_caps[k],
            shared["q_derivative_logs4"][k - 1],
        ))
    T = _interval(context, _at(parent, "actual_long_reference_background_I1_numeric",
                               "scalar_graph", "T"))
    margin = B0 - context.ln(T) + context.ln(200)
    require(common.upper(margin) < 0,
            "Corrected slope margin does not establish |B|/T <= 1/200")
    return B, dict(
        B0=B0,
        rows=B,
        margin=margin,
        ratio_bound="|B|/T<=1/200",
        T=T,
        G_modulus_upper=gbar,
        G_modulus_upper_is_source_bound=True,
        saved_gradient_rows_reused=True,
        old_gradient_recurrence_executed=False,
        safe_changed_input_logs=True,
    )


def _chart_phi_bounds(context: MPIntervalContext, chart: str,
                      signed_log_phi: Sequence[Any], shared: Mapping[str, Any]):
    T = _interval(context, shared["T"]) if "T" in shared else None
    if T is None:
        # T is a physical source scale, not a common arithmetic input.  Read it
        # only if the coordinator did not include it in shared.
        raise ValueError("The corrected chart needs shared physical scale T")
    logphi0_lower = common.point(context, common.lower(signed_log_phi[0]))
    logphi0_upper = common.point(context, common.upper(signed_log_phi[0]))
    lower = logphi0_lower - context.mpf("0.44") * T
    if chart == "long_reshape":
        return lower, logphi0_upper
    Lref = _interval(context, shared["reference_length"])
    return (lower - context.mpf("0.4") * Lref,
            logphi0_upper - context.mpf("0.36") * T)


def _inverse_utheta(context: MPIntervalContext, chart: str, logphi_bounds,
                    shared: Mapping[str, Any], logR):
    lower_phi, upper_phi = logphi_bounds
    # The source's monotone U_theta bounds are kept in their native chart
    # coordinates.  In particular, log R's saved upper endpoint is not a
    # lower bound for U_theta.  Compose each source expression first and box
    # its final lower endpoint before taking the reciprocal.
    logq0 = shared["logq0"]
    if chart == "long_reshape":
        log_utheta_lower = (
            common.point(context, common.lower(shared["logF0"]))
            + context.ln(220) / 2
            + common.point(context, common.lower(shared["logphi0"]))
        )
        log_utheta_upper = (
            shared["T"] / 10 - shared["logC"]
            - common.point(context, common.lower(logq0))
        )
    else:
        log_utheta_lower = (
            shared["T"] / 10 - shared["logC"]
            - common.point(context, common.upper(logq0))
        )
        log_utheta_upper = (
            shared["logP"] - context.mpf("0.8")
            - common.point(context, common.lower(logq0))
        )
    log_utheta_lower = common.point(
        context, common.lower(log_utheta_lower)
    )
    log_utheta_upper = common.point(
        context, common.upper(log_utheta_upper)
    )
    inverse = [
        None if row is None else row - common.point(context, common.lower(log_utheta_lower))
        for row in shared["E_relative"]
    ]
    return inverse, dict(
        formula="1/Utheta=E_relative/Utheta_lower; Utheta=F0*sqrt(2R)*phi",
        log_Utheta_lower=log_utheta_lower,
        log_Utheta_upper=log_utheta_upper,
        log_phi_lower=lower_phi,
        log_phi_upper=upper_phi,
        log_R_source=logR,
        saved_log_R_upper_not_used_as_Utheta_lower=True,
        chart=chart,
        physical_scale_not_materialized=True,
        scalar_endpoint_negation=False,
        signed_function_value_selected=False,
    )


def _root_controls(context: MPIntervalContext, chart: str, aa: list[Any],
                   a_lower, a_upper, p2_rows: Sequence[Any], shared: Mapping[str, Any]):
    """Build changed-input root and conditioning caps through ordinary Z0..Z3."""
    alg = common.LogZ(context)
    a = alg.__class__(context, aa[:4])
    inva = a.reciprocal(
        a,
        common.point(context, common.lower(context.ln(a_lower))),
    )

    # Delta=a-2 and gamma=2*eta-Delta=2*eta+2-a.  eta is retained in its
    # logarithmic source scale; exp(logeta) is never materialised.
    delta_lower = a_lower - 2
    delta_upper = a_upper - 2
    require(common.upper(delta_upper) < 0,
            "Corrected root branch must have Delta<0")
    require(common.upper(shared["logeta"]) <= 0,
            "Corrected root branch must have eta<=1")
    gamma_lower = 2 - a_upper
    require(common.lower(gamma_lower) > 0,
            "Corrected root branch lost positive gamma lower bound")
    gamma_upper_value = common.logsum(
        context,
        context.ln(2 - a_lower),
        context.ln(2) + shared["logeta"],
    )
    gamma = alg.__class__(context, [gamma_upper_value] + aa[1:4])
    q2 = alg.scale(alg.mul(gamma, inva), -context.ln(2))
    s = [q2[k] for k in range(4)]
    q2_lower_log = context.ln(gamma_lower) - context.ln(2) - context.ln(a_upper)
    q_lower_log = q2_lower_log / 2
    add = lambda *terms: common.logsum(context, *terms)
    q = [
        s[0] / 2,
        None if s[1] is None else s[1] - context.ln(2) - q_lower_log,
        add(
            None if s[2] is None else s[2] - context.ln(2) - q_lower_log,
            None if s[1] is None else 2 * s[1] - context.ln(4) - 3 * q_lower_log,
        ),
        add(
            None if s[3] is None else s[3] - context.ln(2) - q_lower_log,
            None if s[1] is None or s[2] is None else
            s[1] + s[2] + context.ln(3) - context.ln(4) - 3 * q_lower_log,
            None if s[1] is None else
            3 * s[1] + context.ln(3) - context.ln(8) - 5 * q_lower_log,
        ),
    ]

    p2 = alg.__class__(context, list(p2_rows[:4]))
    qjet = alg.__class__(context, q)
    ujet = alg.scale(alg.mul(p2, qjet),
                     -common.point(context, common.lower(shared["logdstar"])))
    u = [ujet[k] for k in range(4)]
    vjet = alg.mul(ujet, ujet)
    v = [vjet[k] for k in range(4)]
    h = [
        context.mpf(0),
        None if v[1] is None else v[1] - context.ln(2),
        add(
            None if v[1] is None else 2 * v[1] + context.ln("0.75"),
            None if v[2] is None else v[2] - context.ln(2),
        ),
        add(
            None if v[1] is None else 3 * v[1] + context.ln("1.875"),
            None if v[1] is None or v[2] is None else
            v[1] + v[2] + context.ln("2.25"),
            None if v[3] is None else v[3] - context.ln(2),
        ),
    ]
    hjet = alg.__class__(context, h)
    rjet = alg.mul(ujet, hjet)
    r = [rjet[k] for k in range(4)]
    # The zeroth root coordinate has the exact source identity
    # r = u / sqrt(1 + u^2), hence |r| < 1.  Keep that sharp bound while
    # retaining the full product jet for derivative orders 1 through 3.
    r[0] = context.mpf(0)
    alphajet = alg.scale(alg.mul(qjet, hjet), context.ln(2))
    alpha = [alphajet[k] for k in range(4)]
    denominator = alg.add(alg.constant(context.mpf(0)), alg.scale(q2, context.ln(2)))
    Kjet = alg.scale(
        alg.reciprocal(denominator, context.mpf(0)),
        -context.ln(2 * context.pi),
    )
    K = [Kjet[k] for k in range(4)]
    lh = -context.mpf("0.5") * add(context.mpf(0), 2 * u[0])
    log_D_lower = 4 * common.point(context, common.lower(lh)) - context.ln(4)
    log_lambda_lower = -context.ln(2 * context.pi) - common.point(
        context, common.upper(denominator[0])
    )
    delta_logs = common.logs(context, [context.mpf([delta_lower, delta_upper])])
    delta_rows = [delta_logs[0]] + list(aa[1:4])
    root_rows = {
        "a": list(aa[:4]),
        "b": [None] * 4,
        "t0": [None] * 4,
        "Delta": delta_rows,
        "kappa": list(aa[:4]),
        "gamma": [gamma[k] for k in range(4)],
        "q_squared": s,
        "q": q,
        "u": u,
        "h": h,
        "r": r,
        "alpha": alpha,
        "K": K,
        "p2": list(p2_rows[:4]),
    }
    conditioning = dict(
        log_D_lower=log_D_lower,
        log_lambda_lower=log_lambda_lower,
        positive_a_lower=a_lower,
        positive_gamma_lower=gamma_lower,
        positive_q_squared_lower_log=q2_lower_log,
        eta_log=shared["logeta"],
        eta_finite_positive_source=True,
        eta_positive_source=True,
        eta_le_one=common.upper(shared["logeta"]) <= 0,
        Delta_negative=True,
        active_sigma_branch="Delta<eta; sigma(1-Delta/eta)=1",
        denominator_formula="1+2*q^2",
    )
    formula = dict(
        a="(E-2*E_y)/E, bounded by corrected native slope rows",
        b="2*V_y/E=0 on the saved long/reference native branch",
        t0="-b/a=0",
        Delta="a-2",
        gamma="2*eta-Delta=2*eta+2-a",
        q="sigma(1-Delta/eta)*sqrt(gamma/(2*a)); active sigma=1",
        u="p2*q/d_star",
        r="u/sqrt(1+u^2)",
        alpha="2*q/sqrt(1+u^2)",
        K="1/(2*pi*(1+t0^2+2*q^2))",
    )
    return root_rows, conditioning, formula


def _sourceproof(parent: Mapping[str, Any], chart: str, rates: Mapping[str, int],
                 tail: Mapping[str, Any]) -> dict[str, Any]:
    control = _at(parent, "actual_long_reference_control_seams")
    native = _at(parent, "actual_long_reference_native_controls")
    return dict(
        parent_report=PARENT_REPORT,
        source_family=parent.get("source_family"),
        formula_source=(
            "original_long_reference_native_controls._native_p2 and "
            "prescribed_loop_C5_functions.generic_roots"
        ),
        chart=chart,
        actual_native_control_source_executed=False,
        historical_native_caps_used_as_signed_values=False,
        corrected_signed_inlet_rows=True,
        corrected_closed_form_phi_log_jet=True,
        corrected_gradient_rows_reused=True,
        full_P0_retained=True,
        F0_squared_BC_retained=True,
        physical_scales_retained=True,
        nonreset_memories=True,
        positive_history_rates=dict(rates),
        tail=tail,
        branch_hypotheses=dict(
            positive_eta=True,
            eta_le_one=True,
            positive_dstar=True,
            positive_a=True,
            positive_E=True,
            Delta_formula="a-2",
            active_sigma_formula="sigma(1-Delta/eta)=1",
            D_positive="(cos(psi)-r)^2+sin(psi)^2",
            K_positive="1/(2*pi*(1+t0^2+2*q^2))",
            lambda_positive="K*(1+t^2)",
            source_control_packet_present=isinstance(control, Mapping),
            historical_native_packet_is_formula_source=isinstance(native, Mapping),
        ),
        source_executed=False,
        accepted_ancestors_replayed=False,
    )


def _chart(context: MPIntervalContext, parent: Mapping[str, Any], shared: Mapping[str, Any],
           factors: Mapping[str, Any], signed_log_phi: Sequence[Any], phi_log_caps: Sequence[Any],
           safe_B: Sequence[Any], chart: str, incoming: Mapping[str, Sequence[Any]],
           inlet: Mapping[str, Sequence[Any]], tail: Mapping[str, Any]):
    long = chart == "long_reshape"
    ell = factors["ell_long"] if long else factors["ell_reference"]
    erel = factors["E_relative_long"] if long else factors["E_relative_reference"]
    phirel = factors["phi_relative_long"] if long else factors["phi_relative_reference"]
    logphi_bounds = _chart_phi_bounds(context, chart, signed_log_phi, shared)
    logphi_upper = common.upper(logphi_bounds[1])
    phi_caps = [
        None if value is None else value + common.point(context, logphi_upper)
        for value in phirel
    ]
    U_caps = common.logs(context, inlet["U"])
    M_rows = _signed_source(
        context,
        _at(parent, "actual_long_reference_background_I1_numeric", "charts", chart,
            "whole_segment_M_hull"),
        chart + " whole-segment M hull",
    )
    M_caps = common.logs(context, M_rows)
    # The signed inlet U rows are historical source data.  After their safe
    # logarithmic caps are formed, the new product helper supplies the
    # corrected U^2 density without replaying the accepted value product.
    density_A = common.log_product(context, U_caps, U_caps)
    density_B = common.log_product(context, phi_caps, phi_caps)
    density_C = list(density_B)
    rates = {"A": 1, "B": 2, "C": 1}
    native = dict(phi=phi_caps, U=U_caps, M=M_caps)
    for name, density in (("A", density_A), ("B", density_B), ("C", density_C)):
        native[name] = common.history(context, incoming[name], density, rates[name])
    logR = _interval(context, _at(parent, "actual_long_reference_background_I1_numeric",
                                   "charts", chart, "log_R_upper"))
    inverse_Utheta, inverse_detail = _inverse_utheta(
        context,
        chart,
        logphi_bounds,
        dict(shared, E_relative=erel),
        logR,
    )
    p2_rows, p2_detail = common.p2(
        context,
        native,
        shared["F2_logs5"],
        shared["P0_logs5"],
        shared["z"],
        shared["delta"],
        logR,
        inverse_Utheta,
        shared["inverse_L_logs5"],
    )
    if long:
        aa = [
            context.ln(context.mpf("0.88")),
            *[
                None if row is None else row + context.ln(16) -
                context.ln(_interval(
                    context,
                    _at(parent, "actual_long_reference_background_I1_numeric",
                        "scalar_graph", "T"),
                ))
                for row in safe_B[1:]
            ],
        ]
        a_lower, a_upper = context.mpf("0.72"), context.mpf("0.88")
    else:
        aa = [context.ln(context.mpf("0.8")), None, None, None, None]
        a_lower = a_upper = context.mpf("0.8")
    roots, conditioning, root_formula = _root_controls(
        context, chart, aa, a_lower, a_upper, p2_rows, shared,
    )
    scale = dict(
        chart=chart,
        logR=logR,
        logF0=shared["logF0"],
        logPstar=shared["logPstar"],
        logLambda=shared["logLambda"],
        logN=shared["logN"],
        logeta=shared["logeta"],
        logdstar=shared["logdstar"],
        T=shared["T"],
        reference_length=shared["reference_length"],
        logC=shared["logC"],
        logP=shared["logP"],
        physical_radius_used_once=True,
        physical_sqrt_R_over_two_retained=True,
        P0_and_F0_squared_BC_not_reset=True,
    )
    return dict(
        chart=chart,
        ordinary_Z_orders=list(ORDERS),
        root_ordinary_Z_orders=[0, 1, 2, 3],
        root_cap_order_note="Root intermediates are supplied only through Z3; native and p2 rows retain Z4.",
        native_caps=native,
        phi_relative=list(phirel),
        E_relative=list(erel),
        inverse_Utheta=list(inverse_Utheta),
        inverse_Utheta_detail=inverse_detail,
        p2=list(p2_rows),
        p2_detail=p2_detail,
        a=list(aa),
        root_caps=roots,
        conditioning=conditioning,
        scale=scale,
        sourceproof=_sourceproof(parent, chart, rates, tail),
        slope=dict(B0=safe_B[0], B=list(safe_B), margin=shared["slope_margin"]),
        native_sources=dict(
            signed_inlet_phi="actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.phi",
            signed_inlet_U="actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U",
            signed_whole_segment_M=(
                "actual_long_reference_background_I1_numeric.charts."
                + chart + ".whole_segment_M_hull"
            ),
        ),
        root_formula=root_formula,
        history=dict(
            incoming={name: list(incoming[name]) for name in ("A", "B", "C")},
            density_A=density_A,
            density_B=density_B,
            density_C=density_C,
            rates=rates,
            positive_kernel_maximum=True,
            endpoint_bound_not_selected_value=True,
        ),
        signed_function_value_selected=False,
        old_caps_reused_as_signed_inputs=False,
        finite_N_recertified=False,
        seam_proved=False,
    )


def compute(context: MPIntervalContext, parent: Mapping[str, Any],
            shared: Mapping[str, Any]) -> dict[str, Any]:
    """Return the corrected long/reference native cap chain as raw intervals."""
    require(isinstance(context, MPIntervalContext), "Interval context required")
    require(isinstance(parent, Mapping), "Parent source mapping required")
    require(isinstance(shared, Mapping), "Corrected common mapping required")
    scope = _source_scope(parent)
    data = _shared(context, shared)
    inlet = _inlet(context, parent)
    signed_log_phi, phi_log_caps = _log_phi(context, inlet["phi"])
    factors = _relative_factors(context, data, phi_log_caps)
    B, slope = _safe_B(context, parent, data, phi_log_caps,
                       common.logs(context, [data["q_log_values5"][0]]), factors)
    data["slope_margin"] = slope["margin"]
    data["T"] = slope["T"]
    data["logC"] = _interval(
        context,
        _at(parent, "actual_long_reference_background_I1_numeric", "scalar_graph",
            "logC"),
    )
    require(data["logC"] == _interval(
        context,
        _at(parent, "actual_micro_b_p2_mixed_bounds", "whole_tube_amplitude",
            "logC"),
    ), "Background and anchored amplitude logC handles differ")
    require(data["logN"] == _interval(
        context,
        _at(parent, "actual_switch_post_loop_Z3_inputs", "logN"),
    ), "Common N handle differs from the saved switch source")
    data["logP"] = _interval(
        context,
        _at(parent, "actual_long_reference_background_I1_numeric", "scalar_graph",
            "logP"),
    )
    data["logq0"] = data["q_log_values5"][0]
    data["logphi0"] = signed_log_phi[0]
    data["reference_length"] = _interval(
        context,
        _at(parent, "actual_long_reference_background_I1_numeric", "scalar_graph",
            "reference_length"),
    )
    incoming_signed = {
        name: common.logs(context, inlet[name]) for name in ("A", "B", "C")
    }
    charts: dict[str, Any] = {}
    tail = dict(
        strategy="new_safe_corrected_long_wholechart_endpoint_enclosure",
        formula="reference incoming A/B/C = corrected long wholechart native caps",
        historical_tail_formula="Aend=U^2+(A_in-U^2)*[0,2^-1024]",
        corrective_computation=False,
        corrected_wholechart_memory_propagation=True,
        wholechart_endpoint_enclosure=True,
        endpoint_contained=True,
        actual_endpoint_value_selected=False,
        saved_tail_data_preserved=True,
        direct_tail_recomputation=False,
        accepted_old_builder_replayed=False,
        source_handle="actual_switch_post_background_I1_numeric.incoming_fields.A",
    )
    charts["long_reshape"] = _chart(
        context, parent, data, factors, signed_log_phi, phi_log_caps, B,
        "long_reshape", incoming_signed, inlet, tail,
    )
    # Preserve nonreset histories into reference by passing the corrected long
    # whole-chart endpoint enclosures.  They contain the actual endpoint
    # memory, but remain conservative caps rather than selected values.
    reference_incoming = {
        "A": list(charts["long_reshape"]["native_caps"]["A"]),
        "B": list(charts["long_reshape"]["native_caps"]["B"]),
        "C": list(charts["long_reshape"]["native_caps"]["C"]),
    }
    charts["reference"] = _chart(
        context, parent, data, factors, signed_log_phi, phi_log_caps, B,
        "reference", reference_incoming, inlet, tail,
    )
    return dict(
        stage="original_long_reference_directed_controls",
        schema="original_long_reference_directed_controls_raw_v1",
        parent_report=PARENT_REPORT,
        parent_report_sha256=PARENT_REPORT_SHA256,
        common_module=COMMON_MODULE,
        source_scope=scope,
        ordinary_Z_orders=list(ORDERS),
        charts=charts,
        slope=slope,
        shared_inputs_consumed=[
            "z", "delta", "logLambda", "logF0", "logPstar", "logN", "logeta",
            "logdstar", "gradient_logs4", "F2_logs5", "P0_logs5",
            "q_log_values5", "q_derivative_logs4", "inverse_L_logs5",
        ],
        sourceproof=dict(
            source_family=parent.get("source_family"),
            source_executed=False,
            accepted_ancestors_replayed=False,
            old_native_module_formula_only=True,
            corrected_native_caps=True,
            corrected_roots_and_conditioning=True,
            all_source_branch_hypotheses_rechecked=True,
            eta_positive_source=True,
            eta_finite_positive_source=True,
            eta_le_one=common.upper(data["logeta"]) <= 0,
            dstar_positive_source=True,
            dstar_finite_positive_source=True,
            dstar_less_than_one=common.upper(data["logdstar"]) < 0,
            full_P0_retained=True,
            full_F0_squared_BC_retained=True,
            nonreset_histories=True,
            reference_memory_conservative_caps=True,
            reference_memory_actual_endpoint_not_selected=True,
            finite_N_recertified=False,
            modified_seams_proved=False,
        ),
        limits=dict(
            raw_interval_output=True,
            json_serialization_deferred=True,
            historical_caps_not_recertified=True,
            finite_N_numeric=False,
            Rm_numeric=False,
            global_closure=False,
            full_goal_complete=False,
        ),
    )


__all__ = (
    "MODULE",
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "COMMON_MODULE",
    "compute",
)
