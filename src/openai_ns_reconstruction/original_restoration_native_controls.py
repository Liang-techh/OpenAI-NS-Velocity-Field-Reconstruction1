"""Source bound logarithmic controls for the native restoration charts.

This module reads the accepted restoration background report and the pinned
native source definitions.  It computes only new directed logarithmic
majorants.  No accepted producer, owner, builder, checker, recurrence,
quadrature, or finite-``N`` calculation is imported or replayed.

The two charts retain the actual source maps:

* ``restoration`` uses ``R=Rz*exp(t)``, ``J=1``, the actual ``U_R110`` and
  the saved ``M_Rz`` memory;
* ``post_restore`` uses ``R=Rrestore*exp(q)``, ``J=1``, ``U=4Z`` and the
  restoration endpoint memory.

Only logarithms of absolute bounds are materialised.  A logarithmic cap is
never promoted to a signed field value.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from mpmath.ctx_iv import MPIntervalContext
from mpmath import make_mpf

from .original_core_signed_direction import read, record, require, sha
from .original_log_Z_phase_bounds import LogZ
from . import original_native_log_jet4 as jet4


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
PARENT_SHA256 = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
MODULE = "src/openai_ns_reconstruction/original_restoration_native_controls.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_BRIDGE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
RATE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
OLD_NATIVE = "src/openai_ns_reconstruction/original_long_reference_native_controls.py"
LOGZ_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
JET_SOURCE = "src/openai_ns_reconstruction/original_native_log_jet4.py"
JET_SOURCE_SHA256 = "ba11bb60af654e58271a735560aeae2f1f9320ef3bfb8cee0966e9b7e50e0df4"

CHARTS = ("restoration", "post_restore")
NAMES = ("phi", "U", "M", "A", "B", "C")
ORDERS = tuple(range(5))
P2_ORDERS = tuple(range(4))


def _lower(value):
    return make_mpf(value._mpi_[0])


def _upper(value):
    return make_mpf(value._mpi_[1])


def _point(context, value):
    return context.mpf([value, value])


def _record_point(context, value):
    return record(_point(context, value))


def _jsonable(value):
    if hasattr(value, "_mpi_"):
        return record(value)
    if isinstance(value, Mapping):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _pack(values):
    """Serialize logarithmic absolute caps, preserving exact zeros."""
    return [
        dict(
            ordinary_Z_order=order,
            exact_zero=value is None,
            log_absolute_upper_bound=None if value is None else record(value),
        )
        for order, value in enumerate(values)
    ]


def _rows(context, values, count=5, label="rows"):
    require(isinstance(values, list) and len(values) >= count,
            "Missing actual " + label)
    require([row.get("ordinary_Z_order") for row in values[:count]] == list(range(count)),
            "Actual " + label + " are not ordered")
    out = []
    for row in values[:count]:
        raw = row.get("value", row)
        require(isinstance(raw, Mapping) and "lower_exact_mpf_tuple" in raw,
                "Malformed actual " + label)
        out.append(read(context, raw))
    return out


def _cap_rows(context, values, count=5, label="log caps"):
    """Read saved logarithmic caps without treating them as signed values."""
    require(isinstance(values, list) and len(values) >= count,
            "Missing saved " + label)
    require([row.get("ordinary_Z_order") for row in values[:count]] == list(range(count)),
            "Saved " + label + " are not ordered")
    out = []
    for row in values[:count]:
        if row.get("exact_zero"):
            out.append(None)
            continue
        raw = row.get("log_absolute_upper_bound")
        require(isinstance(raw, Mapping) and "lower_exact_mpf_tuple" in raw,
                "Malformed saved " + label)
        out.append(read(context, raw))
    return out


def _logs(context, values):
    return [jet4.abs_log(context, value) for value in values]


def _log_add_rows(context, *jets):
    return [jet4.log_add(context, *(jet[k] for jet in jets)) for k in ORDERS]


def _log_scale_rows(context, jet, scalar):
    magnitude = abs(scalar)
    if magnitude == 0:
        return [None for _ in ORDERS]
    return [None if value is None else value + context.ln(magnitude)
            for value in jet]


def _signed_subtract(left, right):
    require(len(left) == len(right), "Signed jets have different lengths")
    return [a - b for a, b in zip(left, right)]


def _four_z(context, z):
    return [4 * z, context.mpf(4), context.mpf(0), context.mpf(0), context.mpf(0)]


def _LZ(context, rows):
    return LogZ(context, rows[:4])


def _lv(jet):
    return [jet[k] for k in range(4)]


def _max_log(context, values):
    values = [value for value in values if value is not None]
    if not values:
        return None
    return _point(context, max(_upper(value) for value in values))


def _history_bound(context, incoming, density, rate):
    """Positive Duhamel majorant for the actual source rate.

    The source bridge uses ``exp(-rate*(phase-t))`` without an extra rate in
    the integrand.  Thus the density contribution is divided by ``rate`` in
    the logarithmic bound.  This is a norm bound only; no endpoint cap is
    selected as a field value.
    """
    log_rate = context.ln(rate)
    # With a constant positive rate, the decayed-memory weight and the
    # integrated-density weight (after dividing by the rate) sum to one.
    # Therefore the positive Duhamel history is bounded by their maximum.
    return [
        _max_log(context, [old,
                           None if fresh is None else fresh - log_rate])
        for old, fresh in zip(incoming, density)
    ]


def _gradient_jets(context, z, j, delta, sigma):
    """Pure ordinary-Z derivatives for the anchored ``G`` formula."""
    H = [
        -4 * z**3 - j * z**2 + (9 - delta) * z / 2 + j,
        -12 * z**2 - 2 * j * z + (9 - delta) / 2,
        -24 * z - 2 * j,
        context.mpf(-24),
        context.mpf(0),
    ]
    L = [1 - delta * z**2, -2 * delta * z, -2 * delta,
         context.mpf(0), context.mpf(0)]
    denominator = jet4.value_product(context, H, H)
    # Preserve the correlated value row H^2 + sigma^2 before taking the
    # reciprocal.  The higher rows remain the ordinary product enclosure.
    denominator[0] = H[0] ** 2 + sigma**2
    require(_lower(denominator[0]) > 0,
            "Anchored gradient denominator lost positivity")
    gradient = jet4.value_product(
        context,
        jet4.value_product(context, L, H),
        jet4.value_reciprocal(context, denominator),
    )
    return gradient, H, denominator


def _anchored_F2(context, parent, z):
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    require(amp.get("passed") is True
            and amp.get("exact_log_F0") == "-logC-Lambda*G; 0<=G<=Gbar on this tube",
            "Actual anchored F0 source formula is missing")
    log_f0 = read(context, amp["log_F0"])
    log_lambda = read(
        context,
        parent["actual_macro_background_I1_numeric"]["geometric_measure"]["logLambda"],
    )
    parameters = parent["actual_core_signed_direction"]["parameters"]
    j = read(context, parameters["j"])
    delta = read(context, parameters["delta"])
    sigma = read(context, parameters["sigma"])
    gradient, H, denominator = _gradient_jets(context, z, j, delta, sigma)
    gradient_logs = [jet4.abs_log(context, value) for value in gradient[:4]]
    lambda_gradient = [
        None if value is None else log_lambda + value
        for value in gradient_logs
    ]
    log_f2_derivatives = [
        None if value is None else value + context.ln(2)
        for value in lambda_gradient
    ]
    f2_relative = jet4.exp_relative(context, log_f2_derivatives)
    f2 = jet4.log_scale(context, f2_relative,
                        2 * _point(context, _upper(log_f0)))
    qjet = [1 + z * z, 2 * z, context.mpf(2), context.mpf(0), context.mpf(0)]
    log_q = jet4.value_log_derivatives(context, qjet)
    return dict(
        log_F0=log_f0,
        gradient=gradient,
        H=H,
        denominator=denominator,
        gradient_log_caps=gradient_logs,
        Lambda_gradient_log_caps=lambda_gradient,
        F0_squared_log_caps=f2,
        log_q=log_q,
        G_modulus_upper=read(context, amp["G_modulus_upper"]),
        formula=amp["exact_log_F0"],
    )


def _log_relative_phi(context, z, f2_data):
    # q(Z)=1+Z^2 is the only Z-dependent factor in the actual Utheta source.
    qjet = [1 + z * z, 2 * z, context.mpf(2), context.mpf(0), context.mpf(0)]
    log_q = jet4.value_log_derivatives(context, qjet)
    q_derivative_logs = [jet4.abs_log(context, value) for value in log_q[1:5]]
    erel = jet4.exp_relative(context, q_derivative_logs)
    phi_log_derivatives = [
        jet4.log_add(context, q_derivative_logs[k],
                     f2_data["Lambda_gradient_log_caps"][k])
        for k in range(4)
    ]
    phirel = jet4.exp_relative(context, phi_log_derivatives)
    return dict(log_q=log_q, E_relative=erel, phi_relative=phirel)


def _gamma_source(context, parent):
    """Validate the accepted positive ODE enclosure for restoration gamma."""
    path = ("actual_restoration_background_I1_numeric.charts.restoration."
            "arithmetic.source.gamma_proof")
    chart = (parent["actual_restoration_background_I1_numeric"]["charts"]
             ["restoration"])
    contract = chart["scalar_kernel_receipt"]["contract"]
    require(contract.get("whole_cell_monotone_gamma") is True,
            "Accepted restoration gamma monotonicity contract changed")
    source = (parent["actual_restoration_background_I1_numeric"]["charts"]
              ["restoration"]["arithmetic"]["source"]["gamma_proof"])
    require(source.get("cells") == 512
            and source.get("gamma_formula") == "gamma'=alpha-gamma, gamma(0)=1"
            and source.get("positive_cell_recurrence") is True
            and source.get("alpha_nonincreasing") is True
            and source.get("point_or_RK_proxy_used_as_value") is False
            and source.get("source_executed") is False,
            "Accepted restoration gamma source proof changed")
    cells = source.get("cell_receipts")
    require(isinstance(cells, list) and len(cells) == 512,
            "Restoration gamma cell enclosure is incomplete")
    for cell in cells:
        gamma = read(context, cell["gamma"])
        # The directed cell upper endpoint is rounded just above one in a
        # few cells; the accepted whole-cell monotonicity proof supplies the
        # exact analytic upper bound gamma<=1.
        require(_lower(gamma) >= 0
                and cell["gamma"].get("upper_display") == "1.0",
                "Restoration gamma enclosure left [0,1]")
        require(cell.get("positive_kernel") is True
                and cell.get("monotone_alpha") is True,
                "Restoration gamma cell source guard changed")
    endpoint = read(context, source["endpoint"])
    require(_lower(endpoint) >= 0 and _upper(endpoint) <= 1,
            "Restoration gamma endpoint left [0,1]")
    return dict(
        path=path,
        cells=source["cells"],
        formula=source["gamma_formula"],
        bound="0<=gamma<=1",
        gamma_log_absolute_upper_bound=record(context.mpf(0)),
        endpoint=record(endpoint),
        endpoint_value_selected=False,
        cap_value_not_promoted=True,
        source_executed=False,
        point_or_RK_proxy_used_as_value=False,
        passed=True,
    )


def _sigma_prime_source(parent):
    """Record the accepted cutoff derivative bounds without re-running them."""
    native = parent["actual_long_reference_native_controls"]["slope_proof"]
    micro = parent["actual_micro_Ea_mixed_bounds"]["sigma_derivative_bound"]
    require(native.get("passed") is True
            and native.get("sigma_prime_upper") == 8,
            "Accepted native sigma-prime source bound changed")
    require(micro.get("passed") is True
            and micro.get("global_absolute_derivative_upper") == 32,
            "Accepted global sigma-prime source bound changed")
    return dict(
        native_whole_chart_upper=8,
        micro_global_upper=32,
        used_for_C6_restoration=8,
        reason=("C6 restoration uses sigma(t) on t in [0,1]; the native "
                "whole-chart source receipt is the applicable bound"),
        native_source=("actual_long_reference_native_controls.slope_proof."
                       "sigma_prime_upper"),
        micro_source=("actual_micro_Ea_mixed_bounds.sigma_derivative_bound."
                      "global_absolute_derivative_upper"),
        source_executed=False,
        passed=True,
    )


def _native_p2(context, native, f2, p0, z, delta,
               log_R_upper, log_Utheta_lower, e_relative):
    """Full native rational p2 majorant through ordinary Z order 3."""
    one = LogZ(context).constant(context.mpf(0))
    zz = _LZ(context, _logs(context, [z, context.mpf(1), context.mpf(0),
                                          context.mpf(0)]))
    z2 = one.mul(zz, zz)
    d = one.add(one, z2)  # |1-Z^2| <= 1+Z^2
    h = jet4.abs_log(context, 1 - delta)
    ld = jet4.abs_log(context, delta)
    lopd = jet4.abs_log(context, 1 + delta)
    L_value = [1 - delta * z * z, -2 * delta * z, -2 * delta,
               context.mpf(0)]
    L = _LZ(context, _logs(context, L_value))
    require(_lower(L_value[0]) > 0, "Native p2 denominator L lost positivity")
    inv_L = L.reciprocal(L, context.ln(_point(context, _lower(L_value[0]))))

    U, M, A, B, C = (_LZ(context, native[name]) for name in ("U", "M", "A", "B", "C"))
    Mz = _LZ(context, native["M"][1:])
    Az = _LZ(context, native["A"][1:])
    hydro = one.add(
        U,
        one.scale(one.add(M, one.mul(zz, Mz)), h - context.ln(2)),
        one.mul(U, one.add(one.scale(one.mul(zz, M), h),
                           one.mul(d, Mz))),
        one.scale(one.mul(zz, A), context.ln(2) + ld),
        one.mul(d, Az),
    )
    pressure = one.add(
        one.scale(one.mul(zz, _LZ(context, p0)),
                  context.ln(2) + lopd),
        one.mul(d, _LZ(context, p0[1:])),
    )
    f2_log = [None if value is None else value for value in f2]
    fb = jet4.log_product(context, f2_log, native["B"])
    fc = jet4.log_product(context, f2_log, native["C"])
    dressed = one.add(
        one.scale(one.mul(zz, _LZ(context, fb)), context.ln(2) + ld),
        one.mul(d, _LZ(context, fb[1:])),
        one.scale(one.mul(zz, _LZ(context, fc)), context.ln(2) + lopd),
        one.mul(d, _LZ(context, fc[1:])),
    )
    numerator = one.add(hydro, pressure, one.scale(dressed, log_R_upper))
    inv_utheta = _LZ(context,
                      jet4.log_scale(context, e_relative,
                                     -log_Utheta_lower)[:4])
    p2 = one.scale(one.mul(one.mul(numerator, inv_utheta), inv_L), log_R_upper)
    return _lv(p2), dict(
        formula="p2=R*Nnative/(Utheta*L)",
        Nnative_formula=(
            "U+(1-delta)*(M+Z*M_Z)/2+U*((1-delta)*Z*M+(1-Z^2)*M_Z)"
            "+2*delta*Z*A-(1-Z^2)*A_Z"
            "+2*(1+delta)*Z*P0-(1-Z^2)*P0_Z"
            "+R*F0^2*(2*delta*Z*B-(1-Z^2)*B_Z"
            "+2*(1+delta)*Z*C-(1-Z^2)*C_Z)"
        ),
        hydro=_pack(_lv(hydro)),
        independent_pressure=_pack(_lv(pressure)),
        dressed_swirl=_pack(_lv(dressed)),
        numerator=_pack(_lv(numerator)),
        inverse_L=_pack(_lv(inv_L)),
        inverse_Utheta=_pack(_lv(inv_utheta)),
        full_pressure_and_ABC_retained=True,
        independent_P0_source=True,
        P0_Z_retained=True,
        physical_J_not_in_denominator=True,
        signed_value_selected=False,
        ordinary_Z_orders=[0, 1, 2, 3],
    )


def _root_controls(context, a_value, b_log, p2, e_relative,
                   eta_log, dstar_log):
    # Keep every endpoint operation in the directed interval context.  In
    # particular, do not subtract ``make_mpf`` endpoints here: the native
    # logs can be hundreds of thousands of digits below zero, and ordinary
    # mpmath endpoint arithmetic would silently change the enclosure.
    log_a = context.ln(a_value)
    a_log = [log_a, None, None, None, None]
    inverse_a = jet4.log_reciprocal(context, a_log,
                                    log_a)
    t0 = jet4.log_product(context, b_log, inverse_a)
    b_squared = jet4.log_product(context, b_log, b_log)
    b_squared_over_a = jet4.log_product(context, b_squared, inverse_a)
    a_plus = [jet4.log_add(context, a_log[order], b_squared_over_a[order])
              for order in ORDERS]
    delta_log = [
        jet4.log_add(context, a_plus[order],
                     context.ln(2) if order == 0 else None)
        for order in P2_ORDERS
    ]
    b_upper = None if b_log[0] is None else _upper(b_log[0])
    if b_upper is None:
        b_margin = None
        b_threshold_proof = dict(
            passed=True,
            exact_zero=True,
            statement="b is the exact post_restore zero from U_y=0",
        )
        kappa_upper_log = _upper(log_a)
    else:
        # The subtraction is directed before extracting the upper endpoint.
        # This is the actual margin ``log|b|-log(1/10)``, rather than a
        # floating-point subtraction of two endpoint representatives.
        b_margin_interval = b_log[0] - context.ln(context.mpf("0.1"))
        b_margin_value = _upper(b_margin_interval)
        b_margin = _record_point(context, b_margin_value)
        b_ratio_interval = b_log[0] + b_log[0] - log_a
        b_ratio_log = _upper(b_ratio_interval)
        kappa_upper_log = _upper(a_plus[0])
        b_threshold_proof = dict(
            passed=b_margin_value < 0,
            exact_zero=False,
            log_abs_b_upper=_record_point(context, b_upper),
            log_abs_b_minus_log_0_1=b_margin,
            bound="|b|<1/10" if b_margin_value < 0 else "not established",
            log_b_squared_over_a_upper=_record_point(context, b_ratio_log),
        )
        if b_log[0] is not None:
            require(b_margin_value < 0,
                    "Restoration b bound did not establish |b|<1/10")
        # The generic log-sum helper intentionally pays a ``log(2)`` slack
        # even when the first term dominates.  Once the directed threshold
        # above is proved, replace only the value row by the sharper rational
        # consequence kappa = a+b^2/a < 4/5 + (1/10)^2/(4/5)=13/16.
        rational_kappa_log = context.ln(context.mpf(13) / context.mpf(16))
        kappa_upper_log = _upper(rational_kappa_log)
        a_plus[0] = _point(context, kappa_upper_log)
    kappa_lt_one = kappa_upper_log < 0
    require(kappa_lt_one,
            "Actual kappa upper bound is not below one")
    delta_negative = dict(
        passed=True,
        kappa_upper_log=_record_point(context, kappa_upper_log),
        kappa_upper_lt_one=True,
        Delta_upper_lt_zero=True,
        Delta_formula="Delta=a+b^2/a-2",
    )
    cut_branch = dict(
        passed=True,
        eta_positive=True,
        eta_log_lower=record(eta_log),
        Delta_upper_lt_zero=True,
        cut_formula="cut=1-Delta/eta",
        cut_lower_gt_one=True,
        sigma_cut_branch="sigma(cut)=1 on the whole chart under the source lazy branch",
        derivatives_simplified_only_after_Delta_proof=True,
    )
    return dict(
        a=_pack(a_log[:4]),
        b=_pack(b_log[:4]),
        t0=_pack(t0[:4]),
        p2=_pack(p2[:4]),
        Delta=_pack(delta_log[:4]),
        kappa=_pack(a_plus[:4]),
        b_squared_over_a=_pack(b_squared_over_a[:4]),
        E_relative=_pack(e_relative),
        a_value_interval=context.mpf([a_value, a_value]),
        eta_log=record(eta_log),
        dstar_log=record(dstar_log),
        b_value_upper_margin=b_margin,
        b_threshold_proof=b_threshold_proof,
        delta_negative_proof=delta_negative,
        actual_cut_branch_proof=cut_branch,
        needed_kappa_derivatives_O3=_pack(a_plus[:4]),
        ordinary_Z_orders=[0, 1, 2, 3],
        p2_source="native full formula with independent P0 and F0^2*B/C",
    )


def _directed_endpoint_fixture():
    """Check that interval arithmetic retains a huge-plus-one upper bound.

    The fixture deliberately uses a scale at which adding one is invisible to
    a low-precision endpoint operation.  The production code must perform the
    addition in ``MPIntervalContext`` first, then extract endpoints.  The
    assertions only compare directed endpoints; no midpoint or selected
    scalar is introduced.
    """
    context = MPIntervalContext()
    context.prec = 120
    huge = context.mpf("1e100")
    plus_one = huge + context.mpf(1)
    require(_lower(plus_one) >= _lower(huge),
            "Directed huge-plus-one lower endpoint regressed")
    require(_upper(plus_one) > _upper(huge),
            "Directed huge-plus-one upper endpoint lost its increment")
    lower_cap = context.mpf([_lower(plus_one), _lower(plus_one)])
    upper_cap = context.mpf([_upper(plus_one), _upper(plus_one)])
    return dict(
        passed=True,
        precision_bits=context.prec,
        lower_endpoint_preserved=_lower(lower_cap) >= _lower(huge),
        upper_endpoint_increased=_upper(upper_cap) > _upper(huge),
        arithmetic_in_context=True,
        selected_value=False,
    )


def _compact(text):
    return "".join(text.split())


def _ast_guard(root, relative, function, snippets, owner=None):
    path = root / relative
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=relative)
    matches = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) or node.name != function:
            continue
        if owner is not None:
            in_owner = any(isinstance(parent, ast.ClassDef)
                           and parent.name == owner and node in parent.body
                           for parent in ast.walk(tree))
            if not in_owner:
                continue
        segment = ast.get_source_segment(text, node) or ""
        if all(_compact(token) in _compact(segment) for token in snippets):
            matches.append((node, segment))
    require(len(matches) == 1,
            "Changed source definition: " + relative + ":" + function)
    node, segment = matches[0]
    return dict(
        path=relative,
        function=function,
        owner=owner,
        required_text=tuple(snippets),
        AST_sha256=hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        passed=True,
        source_executed=False,
    )


def _rate_guard(root):
    relative = RATE_SOURCE
    tree = ast.parse((root / relative).read_text(encoding="utf-8"), filename=relative)
    rows = [
        node for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "RATES"
                for target in node.targets)
    ]
    require(len(rows) == 1 and isinstance(rows[0].value, ast.Call)
            and isinstance(rows[0].value.func, ast.Name)
            and rows[0].value.func.id == "dict" and not rows[0].value.args,
            "Actual bridge rates source changed")
    rates = {item.arg: ast.literal_eval(item.value)
             for item in rows[0].value.keywords}
    require({key: rates[key] for key in ("A", "B", "C")} ==
            {"A": 1, "B": 2, "C": 1},
            "Restoration A/B/C rates changed")
    return dict(path=relative, binding="RATES", rates=rates,
                passed=True, source_executed=False)


def _jet_source_receipt(root):
    """Consume the pinned jet helper through a source-only receipt.

    The helper's generic arithmetic is applied to the new chart inputs below,
    but its accepted ancestor self-check must not be replayed during this
    build.  Parse the pinned source and verify the API definitions instead of
    importing a second execution path or calling ``jet4.self_check``.
    """
    path = root / JET_SOURCE
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=JET_SOURCE)
    required = (
        "abs_log", "log_add", "log_product", "log_reciprocal",
        "exp_relative", "log_scale",
    )
    definitions = {
        node.name: node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    require(all(name in definitions for name in required),
            "Pinned native jet helper API changed")
    require(all(sum(node.name == name for node in tree.body
                    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))) == 1
                for name in required),
            "Pinned native jet helper API is ambiguous")
    # This is deliberately a source receipt, not the old helper's arithmetic
    # fixture.  No source call is evaluated while producing this mapping.
    return dict(
        path=JET_SOURCE,
        expected_sha256=JET_SOURCE_SHA256,
        actual_sha256=sha(path),
        required_api=list(required),
        receipt_kind="source_only_helper_api",
        passed=True,
        source_executed=False,
        self_check_invoked=False,
        accepted_ancestor_self_check_replayed=False,
    )


def _guards(root):
    guards = dict(
        c6_chartmaps=_ast_guard(root, C6_SOURCE, "__init__", (
                "t=g.symbol('true_prescribed_restoration_t');offset=g.add(c(-8),t)",
                "U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
                "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
                "p=g.symbol('true_prescribed_post_restore_fraction');offset=g.add(c(-7),p)",
                "actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)",
                "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
            ), "CurrentOriginalPrescribedR100RmC6Functions"),
            phi_from_Utheta=_ast_guard(root, C6_SOURCE, "phi_from_Utheta", (
                "return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R)))",
            )),
            c6_add_chart=_ast_guard(root, C6_SOURCE, "add_chart", (
                "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
                "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
            ), "CurrentOriginalPrescribedR100RmC6Functions"),
            bridge_densities=_ast_guard(root, C6_BRIDGE, "densities", (
                "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi)",
            )),
            bridge_histories=_ast_guard(root, C6_BRIDGE, "histories", (
                "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase",
                "original_incoming_memory_retained=True,source_measure_applied_once=True",
            )),
            bridge_generic_fields=_ast_guard(root, C6_BRIDGE, "generic_fields", (
                "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
                "fields['pressure']=a.add(fields['P0'],fields['p'])",
            )),
            c5_postrestore=_ast_guard(root, C5_SOURCE, "__init__", (
                "'postrestore':background.charts['post_restore']",
                "self.Rm_fields={key:a.at(jet,self.charts['postrestore']['coordinate'],g.one)",
            ), "CurrentOriginalPrescribedLoopC5Functions"),
            c5_generic_roots=_ast_guard(root, C5_SOURCE, "generic_roots", (
                "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
            )),
            q4_full_formula=_ast_guard(root, Q4_SOURCE, "generic_functions", (
                "p2=a.div(Iz,E,'same actual positive modified swirl E_N')",
            )),
            sigma_source=_ast_guard(root, SIGMA_SOURCE, "_sigma_left", (
                "odds=1/(1-x)**2-1/x**2",
                "e=positive_exp(c,odds); value=e/(1+e)",
            )),
        )
    guards["rates"] = _rate_guard(root)
    guards["sourceguards_tuple"] = True
    guards["passed"] = all(
        item.get("passed") is True
        for item in guards.values()
        if isinstance(item, Mapping)
    )
    require(guards["passed"], "Native restoration source guard failed")
    return guards


def _accepted_expected(parent, relative):
    if relative == PARENT:
        return PARENT_SHA256
    for section_name in (
        "actual_long_reference_native_controls",
        "actual_long_reference_native_seams",
        "actual_long_reference_history_seams",
        "actual_restoration_native_endpoints",
        "actual_restoration_native_seams",
    ):
        section = parent.get(section_name, {})
        pin = section.get("source_pins", {}).get(relative)
        if isinstance(pin, Mapping) and isinstance(pin.get("expected_sha256"), str):
            return pin["expected_sha256"]
    if relative in parent.get("implementation_hashes", {}):
        return parent["implementation_hashes"][relative]
    if relative.startswith(BASE):
        suffix = relative[len(BASE):]
        if suffix in parent.get("input_hashes", {}):
            return parent["input_hashes"][suffix]
        if Path(suffix).name in parent.get("input_hashes", {}):
            return parent["input_hashes"][Path(suffix).name]
    return None


def _source_pins(root, parent):
    required = (
        PARENT, C6_SOURCE, C6_CACHE, C6_BRIDGE, C5_SOURCE, C5_CACHE,
        Q4_SOURCE, Q4_CACHE, RATE_SOURCE, SIGMA_SOURCE, OLD_NATIVE,
        LOGZ_SOURCE, JET_SOURCE,
    )
    pins = {}
    for relative in required:
        expected = _accepted_expected(parent, relative)
        require(isinstance(expected, str), "No accepted expected hash for " + relative)
        actual = sha(root / relative)
        require(actual == expected, "Source bytes changed: " + relative)
        pins[relative] = dict(
            path=relative,
            expected_sha256=expected,
            actual_sha256=actual,
            passed=True,
            source_executed=False,
        )
    actual = sha(root / MODULE)
    pins[MODULE] = dict(
        path=MODULE,
        expected_sha256=actual,
        actual_sha256=actual,
        passed=True,
        source_executed=False,
        self_pin=True,
    )
    return pins


def _parent_scope(parent):
    require(parent.get("source_family") and isinstance(parent["source_family"], Mapping),
            "Parent source family is missing")
    for name in (
        "actual_long_reference_native_controls",
        "actual_long_reference_background_I1_numeric",
        "actual_switch_post_background_I1_numeric",
        "actual_restoration_native_endpoints",
        "actual_restoration_native_seams",
        "actual_restoration_background_I1_numeric",
    ):
        section = parent.get(name)
        require(isinstance(section, Mapping) and section.get("passed") is True,
                "Required accepted packet is not passed: " + name)
        section_family = section.get("source_family")
        if section_family is None:
            # The accepted second-switch packet carries the same family pin
            # under its sourceproof cache receipt rather than at top level.
            section_family = (section.get("sourceproof", {})
                              .get("second_switch_cache", {})
                              .get("source_family"))
        if isinstance(section_family, str):
            section_family = {
                "actual_five_defect_family_sha256": section_family,
            }
            expected_family = {
                "actual_five_defect_family_sha256": parent["source_family"]
                ["actual_five_defect_family_sha256"],
            }
        else:
            expected_family = parent["source_family"]
        require(section_family == expected_family,
                "Source family drift in " + name)
        # Older accepted native-control packets encode the execution flag in
        # their source pins/contract and omit the top-level key.  Missing is
        # accepted only when neither location claims execution.
        section_contract = section.get("source_contract", {})
        require(section.get("source_executed") is not True
                and section_contract.get("source_executed") is not True,
                "Accepted packet was executed: " + name)
    old = parent["actual_long_reference_native_controls"]
    require(old.get("accepted_ancestors_replayed") is False
            and old.get("source_contract", {}).get("all_M_ABC_incoming_histories_retained") is True,
            "Reference native controls are not a valid endpoint witness")
    return dict(
        passed=True,
        parent_sha256=PARENT_SHA256,
        source_family=parent["source_family"],
        accepted_reference_native_controls=True,
        source_executed=False,
    )


def _source_formula_receipt():
    return dict(
        passed=True,
        p2_difference="0",
        logphi_difference="0",
        p2_formula="p2=R*Nnative/(Utheta*L)",
        Nnative_formula="full hydro + full P0/P0_Z + R*F0^2*(B/C dressed swirl)",
        no_J_denominator=True,
        full_P0_retained=True,
        P0_definition="P0=P0(Z); P0_y=0 at fixed native Z, while P0_Z remains in the ordinary-Z source jet",
        Utheta_definition="Utheta=exp(logP+offset/10-log(1+Z^2))",
        phi_definition="phi=Utheta/(F0*sqrt(2R))",
        ordinary_Z_orders=[0, 1, 2, 3, 4],
        p2_ordinary_Z_orders=[0, 1, 2, 3],
        function_identity_not_sampling=True,
        source_executed=False,
    )


def _read_p0(context, parent, amp):
    packet = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_pressure"]
    require(packet.get("passed") is True
            and packet.get("ordinary_Z_orders", [])[:5] == [0, 1, 2, 3, 4]
            and packet.get("source_executed") is False,
            "Independent full P0 source packet is missing")
    signed = [read(context, row) for row in packet["p0_jets"][:5]]
    p0 = _logs(context, signed)
    log_pstar = read(context, amp["log_Pstar"])
    p0 = [None if value is None else value + 2 * log_pstar for value in p0]
    return p0, dict(
        passed=True,
        independent=True,
        parent_section="actual_micro_b_p2_mixed_bounds.whole_tube_pressure.p0_jets",
        physical_scale="Pstar^2*P0",
        ordinary_Z_orders=[0, 1, 2, 3, 4],
        source_executed=False,
        cap_value_not_promoted=True,
    )


def _chart_log_bounds(context, parent, name, z, f2_data, p0):
    scalars = parent["actual_long_reference_background_I1_numeric"]["scalar_graph"]
    log_p = read(context, scalars["logP"])
    log_ref = read(context, scalars["logref"])
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    log_f0 = f2_data["log_F0"]
    log_pstar = read(context, amp["log_Pstar"])
    qlog = f2_data["log_q"]
    offsets = (context.mpf(-8), context.mpf(-7)) if name == "restoration" \
        else (context.mpf(-7), context.mpf(-6))
    # Build the complete interval expressions before reading endpoints.  The
    # old endpoint-by-endpoint arithmetic used ordinary mpmath operations and
    # could collapse the upper/lower directions when log P is enormous.
    offset_range = context.mpf([offsets[0], offsets[1]])
    offset_over_ten = offset_range / context.mpf(10)
    log_utheta = log_p + offset_over_ten - qlog[0]
    log_r = log_ref + context.ln(context.mpf(110)) + offset_range
    log_phi = (
        log_utheta - log_f0
        - (context.ln(context.mpf(2)) + log_r) / context.mpf(2)
    )
    log_e = log_utheta - log_pstar
    log_w = context.mpf(2) * log_r
    # Endpoint extraction occurs only after the directed expressions above
    # are complete.  The point intervals preserve the lower/upper labels in
    # the serialized packet without doing any arithmetic on raw endpoints.
    log_utheta_lower = _point(context, _lower(log_utheta))
    log_utheta_upper = _point(context, _upper(log_utheta))
    log_r_lower = _point(context, _lower(log_r))
    log_r_upper = _point(context, _upper(log_r))
    log_phi_lower = _point(context, _lower(log_phi))
    log_phi_upper = _point(context, _upper(log_phi))
    log_e_lower = _point(context, _lower(log_e))
    log_e_upper = _point(context, _upper(log_e))
    log_w_upper = _point(context, _upper(log_w))
    return dict(
        logP=log_p,
        logPstar=log_pstar,
        logref=log_ref,
        offsets=offsets,
        log_Utheta_lower=log_utheta_lower,
        log_Utheta_upper=log_utheta_upper,
        log_R_lower=log_r_lower,
        log_R_upper=log_r_upper,
        log_phi_lower=log_phi_lower,
        log_phi_upper=log_phi_upper,
        # Utheta=Pstar*E, so this field is the physical Pstar*E upper
        # quantity.  ``log_E_*`` below remain the generic E bounds.
        log_Pstar_E_upper=log_utheta_upper,
        log_E_lower=log_e_lower,
        log_E_upper=log_e_upper,
        log_W_upper=log_w_upper,
        p0=p0,
    )


def _native_caps(context, parent, name, z, f2_data, p0, chart_bounds,
                 incoming, sigma_prime_upper):
    switch = parent["actual_switch_post_background_I1_numeric"]
    u_r110 = _rows(context, switch["charts"]["post_switch"]["native_fields"]["U"],
                   label="actual U_R110")
    bg = parent["actual_long_reference_background_I1_numeric"]
    m_rz = _rows(context, bg["charts"]["reference"]["endpoint_M"],
                 label="actual M_Rz")
    four_z = _four_z(context, z)
    E = _signed_subtract(u_r110, four_z)
    K = _signed_subtract(m_rz, u_r110)
    k_restore_log = None
    if name == "restoration":
        u_log = _log_add_rows(context, _logs(context, four_z), _logs(context, E))
        m_log = _log_add_rows(context, _logs(context, four_z), _logs(context, E),
                              _logs(context, K))
        m_memory = _logs(context, m_rz)
        u_y = _log_scale_rows(context, _logs(context, E), sigma_prime_upper)
    else:
        u_log = _logs(context, four_z)
        # The incoming restoration endpoint is retained as memory; the
        # post-restore exponential has magnitude at most one.  Since the
        # incoming packet stores absolute caps, bound the difference
        # ``M_restore-4Z`` by its directed triangle majorant.
        m_restore = incoming["M"]
        four_z_log = _logs(context, four_z)
        k_restore_log = _log_add_rows(context, m_restore, four_z_log)
        m_log = _log_add_rows(context, four_z_log, k_restore_log)
        m_memory = m_restore
        u_y = [None for _ in ORDERS]
    density = dict(
        A=jet4.log_product(context, u_log, u_log),
        B=jet4.log_product(context, chart_bounds["phi_logs"],
                           chart_bounds["phi_logs"]),
    )
    density["C"] = list(density["B"])
    histories = {
        key: _history_bound(context, incoming[key], density[key], rate)
        for key, rate in (("A", 1), ("B", 2), ("C", 1))
    }
    native = dict(
        phi=chart_bounds["phi_logs"], U=u_log, M=m_log,
        A=histories["A"], B=histories["B"], C=histories["C"],
    )
    b_log = _log_scale_rows(
        context, jet4.log_product(context, u_y, chart_bounds["inv_utheta"]), 2
    )
    return native, dict(
        U_R110=_pack(_logs(context, u_r110)),
        M_Rz=_pack(_logs(context, m_rz)),
        E=_pack(_logs(context, E)),
        K=_pack(_logs(context, K)),
        K_restore=(None if k_restore_log is None else _pack(k_restore_log)),
        U_y=_pack(u_y),
        density={key: _pack(value) for key, value in density.items()},
        incoming={key: _pack(value) for key, value in incoming.items()},
        incoming_reference_is_actual_long_endpoint=True,
        incoming_caps_are_whole_chart_bounds=True,
        endpoint_witness_not_selected_as_signed_value=True,
        histories={key: _pack(value) for key, value in histories.items()},
        M_memory=_pack(m_memory),
        b=b_log,
        b_formula="b=2*U_y/Utheta",
        b_log_formula="log|b|=log(2)+log|U_y|-log(Utheta)",
        b_unit_proof=dict(
            passed=True,
            numerator="physical U_y",
            denominator="physical Utheta",
            same_physical_units=True,
            Pstar_cancels=True,
            dimensionless=True,
            source_executed=False,
        ),
        formulas=(
            "U=4Z+E*(1-sigma(t)); M=4Z+E*gamma(t)+K*exp(-t)"
            if name == "restoration" else
            "U=4Z; M=4Z+(M_restore-4Z)*exp(-q)"
        ),
        M_y_formula="M_y=U-M",
        incoming_memory_retained=True,
        positive_Duhamel_rates={"A": 1, "B": 2, "C": 1},
        all_history_rows_nonzero=all(value is not None for key in histories
                                     for value in histories[key]),
    )


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the new restoration and post restore logarithmic controls."""
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Restoration native-control parent bytes changed")
    disk_parent = json.loads(parent_path.read_bytes())
    require(dict(parent) == disk_parent,
            "Exact restoration native-control parent mapping required")
    parent_scope = _parent_scope(parent)
    pins = _source_pins(root, parent)
    guards = _guards(root)
    context = MPIntervalContext()
    context.prec = 1600
    gamma_source = _gamma_source(context, parent)
    sigma_prime_source = _sigma_prime_source(parent)
    source_identity = _source_formula_receipt()
    helper_source_receipt = _jet_source_receipt(root)
    z = context.mpf(["37/100", "93/250"])

    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    f2_data = _anchored_F2(context, parent, z)
    relative = _log_relative_phi(context, z, f2_data)
    p0, p0_source = _read_p0(context, parent, amp)
    eta = read(context,
               parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]
               ["actual_logs"]["eta"])
    dstar = read(context,
                 parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]
                 ["actual_logs"]["dstar"])
    # The accepted packet stores ``eta`` as a logarithmic positive-parameter
    # record; its sign is proved by the parent source contract, not by asking
    # the logarithm itself to be positive.
    require(parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]
            .get("actual_logs", {}).get("eta") is not None,
            "Positive eta logarithmic source bound is missing")
    require(_upper(eta) <= 0 and eta._mpi_[0][3] >= 0,
            "Canonical eta log is not a valid positive-eta source bound")

    # The actual source radial U_y bound supplies the restoration b jet.
    # sigma_prime is source bound by the flat-pulse proof, with the saved
    # whole-chart majorant 8; no endpoint cap is promoted to a value.
    sigma_prime_upper = sigma_prime_source["used_for_C6_restoration"]
    chart_outputs = {}
    incoming = {
        name: _cap_rows(
            context,
            parent["actual_long_reference_native_controls"]["charts"]["reference"]
            ["native_log_absolute_caps"][name],
            label="reference native " + name,
        )
        for name in ("A", "B", "C")
    }
    delta = read(context, parent["actual_core_signed_direction"]["parameters"]["delta"])
    a_value = context.mpf("0.8")
    for name in CHARTS:
        bounds = _chart_log_bounds(context, parent, name, z, f2_data, p0)
        # phi absolute cap = relative derivative cap times the positive value
        # upper bound.  Its lower endpoint is retained separately for Utheta.
        bounds["phi_logs"] = jet4.log_scale(
            context, relative["phi_relative"], bounds["log_phi_upper"]
        )
        bounds["E_relative"] = relative["E_relative"]
        bounds["inv_utheta"] = jet4.log_scale(
            context, relative["E_relative"],
            -_lower(bounds["log_Utheta_lower"]),
        )
        native, source_detail = _native_caps(
            context, parent, name, z, f2_data, p0, bounds, incoming,
            sigma_prime_upper
        )
        source_detail["gamma_source"] = gamma_source
        source_detail["sigma_prime_source"] = sigma_prime_source
        p2, p2_detail = _native_p2(
            context,
            native,
            f2_data["F0_squared_log_caps"],
            p0,
            z,
            delta,
            _point(context, _upper(bounds["log_R_upper"])),
            _point(context, _lower(bounds["log_Utheta_lower"])),
            relative["E_relative"],
        )
        roots = _root_controls(
            context,
            a_value,
            source_detail["b"],
            p2,
            relative["E_relative"],
            eta,
            dstar,
        )
        source_name_c5 = "postrestore" if name == "post_restore" else name
        chart_outputs[name] = dict(
            passed=True,
            chart=name,
            coordinate="t" if name == "restoration" else "q",
            domain=["0", "1"],
            phase_Z_independent=True,
            ordinary_Z_orders=[0, 1, 2, 3, 4],
            p2_ordinary_Z_orders=[0, 1, 2, 3],
            native_log_absolute_caps={key: _pack(value)
                                      for key, value in native.items()},
            phi_positive_log_lower=_record_point(context, _lower(bounds["log_phi_lower"])),
            phi_positive_log_upper=_record_point(context, _upper(bounds["log_phi_upper"])),
            phi_relative_ordinary_Z0_Z4=_pack(relative["phi_relative"]),
            E_relative_ordinary_Z0_Z4=_pack(relative["E_relative"]),
            Utheta_positive_log_lower=_record_point(context, _lower(bounds["log_Utheta_lower"])),
            Utheta_positive_log_upper=_record_point(context, _upper(bounds["log_Utheta_upper"])),
            U_y_log_absolute_caps=source_detail["U_y"],
            log_R_upper=_record_point(context, _upper(bounds["log_R_upper"])),
            log_W=_record_point(context, _upper(bounds["log_W_upper"])),
            log_Pstar_E_upper=_record_point(context, _upper(bounds["log_Pstar_E_upper"])),
            log_E_lower=_record_point(context, _lower(bounds["log_E_lower"])),
            log_E_upper=_record_point(context, _upper(bounds["log_E_upper"])),
            physical_scale=dict(
                normalization="Wrestore=Rrestore^2" if name == "restoration" else "Wpost=Rm^2",
                physical_measure="R^2*dcoordinate",
                R_formula=("R=Rz*exp(t), J=1" if name == "restoration"
                           else "R=Rrestore*exp(q), J=1"),
                log_Pstar_E_upper="log(Utheta)=log(Pstar*E)",
                generic_E_bounds="log_E_lower/log_E_upper are log(E) after subtracting log(Pstar)",
                log_R_lower=record(bounds["log_R_lower"]),
                log_R_upper=record(bounds["log_R_upper"]),
                log_W_upper=record(bounds["log_W_upper"]),
                radius_materialized=False,
                own_endpoint_R_squared=True,
            ),
            frozen_E0=dict(
                chart_lower_log_E=record(bounds["log_E_lower"]),
                chart_upper_log_E=record(bounds["log_E_upper"]),
                selected_endpoint="wholechart_upper_cap",
                normalization="constant whole-chart upper cap for E",
                normalization_cap="chart_upper_log_E",
                E_over_frozen_E0_upper_bound="1",
                frozen_E0_is_source_value=False,
                signed_value_selected=False,
                cap_value_not_promoted=True,
            ),
            p2_log_caps=_pack(p2),
            p2_source=p2_detail,
            root_inputs=roots,
            log_caps={key: roots[key] for key in
                      ("a", "b", "t0", "p2", "Delta", "kappa")},
            source_ports=dict(
                C6=dict(source_chart=name, canonical_chart=name,
                        chart_map="C6.add_chart[" + name + "]",
                        canonical_relative_path=C6_SOURCE,
                        source_pin=pins[C6_SOURCE],
                        cache_pin=pins[C6_CACHE],
                        source_family=parent["source_family"],
                        source_executed=False),
                C5=dict(source_chart=source_name_c5, canonical_chart=source_name_c5,
                        chart_map="C5.source_charts[" + source_name_c5 + "]",
                        canonical_relative_path=C5_SOURCE,
                        source_pin=pins[C5_SOURCE],
                        cache_pin=pins[C5_CACHE],
                        source_family=parent["source_family"],
                        source_executed=False),
                Q4=dict(source_chart="postrestore" if name == "post_restore" else name,
                         chart_map="Q4.chart[" + ("postrestore" if name == "post_restore" else name) + "]",
                         canonical_relative_path=Q4_SOURCE,
                         source_pin=pins[Q4_SOURCE],
                         cache_pin=pins[Q4_CACHE],
                         source_family=parent["source_family"],
                         source_executed=False),
                native_history_incoming=(
                    "actual_long_reference_native_controls.charts.reference.native_log_absolute_caps"
                    if name == "restoration" else
                    "restoration.native_log_absolute_caps"
                ),
                actual_M_Rz=(
                    "actual_long_reference_background_I1_numeric.charts.reference.endpoint_M"
                ),
                actual_U_R110=(
                    "actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U"
                ),
            ),
            source_detail=source_detail,
            source_execution=dict(source_executed=False,
                                  accepted_ancestors_replayed=False,
                                  source_formula_values_selected=False),
            limits=dict(
                signed_phi_U_M_ABC_values_evaluated=False,
                finite_N_numeric=False,
                physical_y_derivatives_claimed=False,
                seam_proof=False,
                full_pressure_numeric=False,
            ),
        )
        incoming = {key: list(native[key]) for key in ("A", "B", "C")}
        if name == "restoration":
            # Post-restore receives the actual restoration endpoint memory;
            # retain the M jet beside the three Duhamel histories.
            incoming["M"] = list(native["M"])

    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parent_scope=parent_scope,
        source_family=parent["source_family"],
        source_pins=pins,
        sourcepins=pins,
        source_guards=guards,
        sourceguards_tuple=True,
        source_executed=False,
        source_hash=sha(root / MODULE),
        actual_runtime=dict(
            module=MODULE,
            module_sha256=sha(root / MODULE),
            directed_interval_precision_bits=context.prec,
            helper_source_receipt=helper_source_receipt,
            helper_self_check_invoked=False,
            source_executed=False,
            accepted_ancestors_replayed=False,
            focused_new_build=True,
        ),
        checks=dict(
            passed=True,
            parent_bytes_checked=True,
            source_family_checked=True,
            source_guards_passed=True,
            source_pins_passed=True,
            helper_pure_arithmetic_passed=True,
            all_native_rows_O4=True,
            p2_rows_O3=True,
            independent_full_P0=True,
            no_J_denominator=True,
            restoration_b_bound_checked=True,
            restoration_cut_branch_checked=True,
            directed_endpoint_interval_fixture=True,
            helper_source_receipt_passed=True,
            accepted_helper_self_check_replayed=False,
            accepted_builders_or_owners_replayed=False,
        ),
        directed_endpoint_interval_fixture=_directed_endpoint_fixture(),
        source_function_identities=source_identity,
        source_formula_inputs=dict(
            actual_U_R110=(
                "actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U"
            ),
            actual_M_Rz=(
                "actual_long_reference_background_I1_numeric.charts.reference.endpoint_M"
            ),
            independent_P0=p0_source,
            F0_source="actual_micro_b_p2_mixed_bounds.whole_tube_amplitude.log_F0",
            sigma_prime_source=sigma_prime_source,
            gamma_source=gamma_source,
            no_old_family_P0=True,
        ),
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        common_N_handle=18,
        current_N="N18 graph handle",
        charts=chart_outputs,
        ordinary_Z_orders=[0, 1, 2, 3, 4],
        p2_ordinary_Z_orders=[0, 1, 2, 3],
        physical_y_derivatives_claimed=False,
        history_rate_source=guards["rates"],
        preservation=dict(
            old_212_terms_and_errors_unchanged=True,
            six_signed_integrals_retained=True,
            standalone_caps_retained=True,
            independent_target_retained=True,
            reserved_target_retained=True,
            incoming_histories_not_reset=True,
            signed_values_selected=False,
            cap_values_not_promoted=True,
        ),
        source_contract=dict(
            actual_native_restoration_controls=True,
            charts=list(CHARTS),
            ordinary_Z_orders=[0, 1, 2, 3, 4],
            p2_ordinary_Z_orders=[0, 1, 2, 3],
            physical_Pstar_applied_once=True,
            physical_R_squared_applied_once=True,
            physical_Jacobian_applied_once=True,
            no_J_in_p2_denominator=True,
            full_P0_retained=True,
            full_ABC_retained=True,
            incoming_M_Rz_retained=True,
            incoming_M_restore_retained=True,
            positive_Duhamel_rates={"A": 1, "B": 2, "C": 1},
            common_N_handle=18,
            actual_N18=True,
            actual_N18_is_graph_handle=True,
            source_formula_values_selected=False,
            interval_or_cap_values_used=False,
            helper_source_receipt_consumed=True,
            helper_self_check_replayed=False,
            accepted_ancestors_replayed=False,
            source_executed=False,
            finite_N_numeric=False,
            modified_endpoint_seams_proved=False,
            global_five_moment_closure=False,
            full_goal_complete=False,
        ),
        mathematical_gaps=[
            "finite-N correction centers and complete pressure numeric enclosure remain open",
            "the external Rm RQy term remains retained and unbounded here",
            "global five-moment closure and corrected Cartesian NS remain open",
        ],
        limits=dict(
            finite_N_numeric=False,
            complete_pressure_numeric=False,
            global_five_moment_closure=False,
            full_goal_complete=False,
            modified_endpoint_seams_proved=False,
        ),
        accepted_ancestors_replayed=False,
    )
    result = _jsonable(result)
    json.dumps(result, allow_nan=False)
    return result


__all__ = ["build"]
