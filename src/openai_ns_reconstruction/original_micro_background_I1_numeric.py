"""Signed background six-I1 source integrals on the two native micro charts.

This adapter evaluates only saved, source-bound intervals.  The common axis
source is ``4*Z+j``; the native U rows contain the saved actual physical
field, while the mean M rows are reconstructed from the positive Duhamel
equation ``M_y=U-M`` around that axis source.  No producer, owner, recurrence,
norm, or quadrature is run here.

The returned integral rows are normalized by ``W = Ra**2*h_b``.  The saved
``hb`` record is ``loghb``; physical ``h_b`` is ``exp(loghb)`` and
``log(W) = 2*log(Ra) + loghb``.  The exact astronomical physical scale is
retained as a logarithmic decomposition, so no unrepresentable epsilon or
radius is materialized.
"""
from __future__ import annotations

import ast
import hashlib
import json
import gzip
import math
from pathlib import Path

from .original_core_signed_direction import require, sha, read, record, source_guard
from .original_ordinary_jet_bounds import jet_add, jet_mul, jet_neg, jet_div, jet_reciprocal


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_macro_RQy_seam.json"
PARENT_SHA256 = "1ee3ff1234cad6bf79fb24a0a558be50a3c10058a7746c176ea5d4b54d260c3d"

B2_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_b_p2_mixed_bounds.json"
I1_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_micro_I1_defect_bounds.json"
UNIFORM = BASE + "lei_ren_part1_paper_compliant_core_uniform_bounds.json"
UNIFORM_SOURCE = BASE + "lei_ren_part1_paper_compliant_core_uniform_bounds.py"
LOWER_GRAPH = BASE + "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz"
LOWER_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py"
LOOP_GRAPH = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
CORRELATED_SOURCE = "src/openai_ns_reconstruction/original_n1_micro_correlated_m15.py"
JET_SOURCE = "src/openai_ns_reconstruction/original_ordinary_jet_bounds.py"
EA_SOURCE = "src/openai_ns_reconstruction/original_micro_Ea_mixed_bounds.py"
CORE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
CORE_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json"
INLET_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"
BRIDGE_RATES_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
BRIDGE_RATES_CHECK = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json"
MEAN_SOURCE = "src/openai_ns_reconstruction/original_micro_phase_mean_kernel.py"

NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")


def load(path):
    raw = Path(path).read_bytes()
    data = gzip.decompress(raw) if str(path).endswith(".gz") else raw
    return json.loads(data)


def _exact_tuple(c, value):
    """Turn one saved binary MPF tuple into an exact singleton interval."""
    from mpmath import make_mpf

    v = make_mpf(tuple(value))
    return c.mpf([v, v])


def _symmetric(c, radius):
    r = c.fabs(radius).b
    require(r >= 0, "A nonnegative interval radius is required")
    return c.mpf([-r, r])


def _pack(values):
    return [record(value) for value in values]


def _rows(values):
    return [dict(ordinary_Z_order=k, value=record(value)) for k, value in enumerate(values)]


def _pin(root, parent, path):
    """Pin a saved input through the latest report's exact hash table."""
    p = Path(path)
    group = "input_hashes" if path.startswith(BASE) else "implementation_hashes"
    key = p.name if path.startswith(BASE) else path
    expected = parent[group].get(key)
    require(expected is not None, "Latest report has no hash for " + path)
    require(sha(root / path) == expected, "Changed saved source: " + path)
    return dict(path=path, expected_sha256=expected)


def _bridge_rates_guard(root):
    """Read the pinned bridge RATES literal without importing or executing it."""
    path = root / BRIDGE_RATES_SOURCE
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text)
    values = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "RATES"
                for target in node.targets):
            if (isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Name)
                    and node.value.func.id == "dict"):
                values.append({keyword.arg: ast.literal_eval(keyword.value)
                               for keyword in node.value.keywords})
            else:
                values.append(ast.literal_eval(node.value))
    require(len(values) == 1 and values[0] == {"H": 2, "M": 1, "K": 2,
                                               "A": 1, "B": 2, "C": 1},
            "Pinned bridge RATES literal changed")
    return dict(passed=True, path=BRIDGE_RATES_SOURCE,
                rates={key: values[0][key] for key in ("H", "M", "K", "A", "B", "C")},
                source_executed=False)


def _jet_sub(c, left, right):
    return jet_add(c, left, jet_neg(c, right))


def _recover_Q(c, Z, U, M, delta, order=2):
    """Ordinary jet recovery from the exact physical Q4 source relation."""
    require(order == 2 and len(U) >= 3 and len(M) >= 4,
            "Q2 recovery requires U0..2 and M0..3")
    z = [Z, c.mpf(1), c.mpf(0)]
    one = [c.mpf(1), c.mpf(0), c.mpf(0)]
    z2 = jet_mul(c, z, z)
    d = _jet_sub(c, one, z2)
    h = c.mpf(1) - delta
    L = _jet_sub(c, one, [delta * value for value in z2])
    require(L[0].a > 0, "Physical Q denominator 1-delta*Z^2 must be positive")
    mm = list(M[:3])
    mz = list(M[1:4])
    numerator = jet_add(c,
        [2 * value for value in jet_mul(c, z, U[:3])],
        jet_neg(c, [h * value for value in jet_mul(c, z, mm)]))
    numerator = _jet_sub(c, numerator, jet_mul(c, d, mz))
    return jet_div(c, numerator, L)


def _source_contract(root, parent):
    """Check the saved graph identities used by this arithmetic adapter."""
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    native = b2["native_angular_and_axial_provider"]
    require(b2["passed"] and native["passed"] and native["ordinary_Z0_Z4"]
            and native["physical_y0_y1"]
            and native["original_core_U_inlet_and_all_corrections_retained"]
            and native["native_and_comparison_never_identified"]
            and native["physical_dy_hb_ds_once"],
            "Saved actual native U provider is required")
    source = parent["actual_micro_correlated_m15"]
    require(source["passed"] and source["physical_I1_measure"] == "R^2*hb*ds exactly once"
            and source["actual_full_Q4_six_I1_source_implemented"]
            and source["internal_RQy_and_reduced_B_seams_cancel_as_functions"],
            "Same-source physical six-I1 graph is required")
    i1 = parent["actual_micro_I1_defect_bounds"]
    require(i1["passed"] and i1["source_contract"]["original_physical_I1_measure"] ==
            "R^2*hb*ds exactly once" and i1["physical_R_squared_and_hb_applied_once"],
            "Saved finite-N I1 defect packet has the wrong physical measure")
    i1_contract = i1["source_contract"]
    wp = i1["weight_provenance"]
    require(wp["passed"] and wp["physical_measure_handle"] == 8
            and wp["physical_measure_is_original_hb_not_integer_eight"]
            and wp["loghb_exact_tuples_equal_accepted_source_basis"]
            and wp["same_actual_radius_handle"] == 68644
            and wp["native_chart_widths_at_most_one"]
            and wp["actual_second_upper_number"] == 2
            and wp["source_first_lower_is_sc_half"]
            and wp["positive_sc_upper_receipt_consumed"]
            and wp["Z_independent_measure_and_native_endpoints"]
            and wp["no_geometry_owner_or_hb_producer_executed"],
            "The accepted I1 weight provenance is required")
    bg = i1["background_native_norm_provider"]
    require(bg["passed"] and bg["inputs_are_saved_actual_native_bounds"]
            and bg["saved_U_M_already_physical"]
            and bg["inverse_Pstar_in_C6_restored_once_by_Q4"]
            and bg["no_extra_Pstar_on_saved_U_M"]
            and bg["native_mean_is_M_not_R_times_M"]
            and bg["native_not_identified_with_comparison"]
            and bg["new_background_Q_norm_computation"]
            and bg["accepted_U_M_producer_not_replayed"]
            and len(bg["physical_native_U"]) >= 3
            and len(bg["physical_native_mean"]) >= 4
            and len(bg["physical_native_Q"]) >= 3,
            "Saved native U/M physical units are required")
    qinc = parent["actual_micro_Q_increment_bounds"]
    rp = qinc["radius_provenance"]
    require(qinc["passed"] and rp["passed"]
            and rp["R_equals_epsilon_rho"]
            and rp["epsilon_Lambda_equals_one"]
            and rp["exact_logLambda"] == "4*logPstar+1000"
            and rp["actual_micro_radius_handle"] == 68644
            and rp["physical_R_applied_once"]
            and rp["physical_R_Z_independent"]
            and rp["same_loop_Q4_radius_graph_prefix"]
            and rp["actual_rho_source_AST_consumed"]
            and rp["core_recipe_receipt_consumed_not_replayed"]
            and not rp["source_producers_executed"],
            "The accepted Q-increment radius provenance is required")
    qinc_contract = qinc["source_contract"]
    require(qinc_contract["actual_physical_U_and_mean"]
            and qinc_contract["original_mean_transport_rate_one_AST_bound"]
            and qinc_contract["phase_Z_independent_and_total_y_phase_derivative_N"]
            and qinc_contract["same_source_Q_replacement_not_extra_ledger_contribution"]
            and qinc_contract["full_recovered_Q_equals_original_by_same_U_and_mean_function_identity"]
            and qinc_contract["old_source_function_and_all_error_packets_retained"]
            and not qinc_contract["accepted_recovery_controls_or_source_producers_replayed"],
            "The accepted physical U/mean Q-recovery receipt is required")
    mean = parent["actual_micro_mean_Z3_bounds"]
    mean_contract = mean["source_contract"]
    require(mean["passed"] and mean["new_history_controls"]["passed"]
            and mean["new_history_controls"]["physical_hb_applied_once"]
            and mean_contract["passed"]
            and mean_contract["accepted_source_proof_consumed"]
            and mean_contract["exact_transport_identity_consumed"]
            and mean_contract["ordinary_Z_seam_identity_consumed"]
            and mean_contract["second_memory_formula"] == "exp(-hb*(s-1))"
            and mean_contract["first_inlet_zero_only_at_source_owned_left_collar"]
            and mean_contract["second_history_not_reset"]
            and mean_contract["actual_phase_Z_independent"]
            and mean_contract["ordinary_Z_derivatives_have_no_factorials"]
            and mean_contract["source_phase_values_not_chosen"]
            and mean_contract["source_or_checker_not_executed"]
            and mean["physical_Pstar_applied_once_and_common_N_squared_retained"]
            and "hb" in mean["physical_scale_logs"]
            and mean["source_AST"]["phase"]["passed"]
            and mean["source_AST"]["moving"]["passed"]
            and mean["source_AST"]["joined"]["passed"],
            "The accepted nonreset phase-mean transport receipt is required")
    bridge_check = load(root / BRIDGE_RATES_CHECK)
    require(bridge_check["all_passed"]
            and bridge_check["implicit_source_sha256"] == parent["source_family"]["implicit_source_sha256"]
            and bridge_check["datum_enclosure_sha256"] == parent["source_family"]["datum_enclosure_sha256"],
            "The accepted bridge RATES receipt is required")
    rates = _bridge_rates_guard(root)
    q2 = parent["actual_correlated_q2_inlet"]
    sc = parent["actual_native_collar_transport"]["domain"]["sc_exact_mpf_tuple"]
    require(q2["shared_q2_inlet"]["exact_sc_tuple"] == sc
            and parent["actual_cutoff_complete_support_measure"]["source_sc_exact_mpf_tuple"] == sc,
            "The saved sc sources disagree")
    lower = load(root / LOWER_GRAPH)
    cert = lower["lower_collar_bindings"]["geometry_certificate"]
    require(cert["sc_exact_mpf_tuple"] == sc and cert["sc_and_hb_Z_independent"],
            "Original lower-collar sc graph is not the admitted source")
    require(i1_contract["same_native_background_on_both_micro_charts"]
            and i1_contract["source_graph_functions_not_cap_values"]
            and i1_contract["accepted_free_Z_internal_B_and_RQy_seam_consumed"]
            and set(i1_contract["six_signed_defect_handles"]["first_micro"]) == set(NAMES)
            and set(i1_contract["six_signed_defect_handles"]["second_micro"]) == set(NAMES),
            "The accepted joint six-source receipt is required")
    correlated = parent["actual_micro_correlated_m15"]
    require(correlated["same_source_error_function_shared_before_any_future_boxing"]
            and correlated["physical_I1_measure"] == "R^2*hb*ds exactly once"
            and correlated["viscous_endpoint_source_weights"]["joint"] ==
            {"lower": 1, "upper": -1},
            "The accepted joint source identity is required")
    direction = parent["actual_core_signed_direction"]
    require(direction["passed"] and direction["actual_core_only"]
            and direction["actual_core_Dbar_positive"]
            and direction["signed_actual_Phi_rho_negative"],
            "Current actual core direction receipt is required")
    core_check = load(root / CORE_CHECK)
    require(core_check["all_passed"]
            and core_check["source_family"] == parent["source_family"]
            and core_check["current_original_n1_true_core_fixed_point_functions_and_enclosure_bridge_installed"]
            and core_check["source_functions_and_directed_enclosures_are_distinct"]
            and core_check["checks"]["actual_oracle"]["same_checked_nonlinear_finite_packet_shared"],
            "Current core fixed-point receipt is required for the uniform correction rows")
    return dict(passed=True,
        actual_native_U_rows_are_signed_enclosures=True,
        source_axis="U0=4*Z+j; physical U=U0+epsilon*Psi; j retained in the baseline",
        mean_equation="M_y=U-M with positive exp(-h_b*(s-v)) Duhamel kernel",
        exact_Q="Q=(2*Z*U-(1-delta)*Z*M-(1-Z^2)*M_Z)/(1-delta*Z^2)",
        same_physical_measure="R^2*h_b*ds exactly once (source receipt: R^2*hb*ds exactly once)",
        sc_source_tuple=sc,sc_Z_independent=True,
        finite_N_defect_packet_same_W=True,
        actual_micro_I1_integrals_key="actual_micro_I1_defect_bounds",
        I1_weight_provenance=wp,
        native_background_unit_receipt=dict(
            saved_U_M_already_physical=True,
            no_extra_Pstar_on_saved_U_M=True,
            native_mean_is_M_not_R_times_M=True),
        Q_increment_radius_receipt=dict(
            actual_micro_radius_handle=rp["actual_micro_radius_handle"],
            exact_logLambda=rp["exact_logLambda"],
            physical_R_applied_once=True,
            physical_R_Z_independent=True),
        phase_mean_receipt=dict(
            second_memory_formula="exp(-h_b*(s-1))",
            second_memory_formula_source=mean_contract["second_memory_formula"],
            second_history_not_reset=True,
            new_history_controls_passed=mean["new_history_controls"]["passed"],
            physical_hb_weight_applied_once=mean["new_history_controls"]["physical_hb_applied_once"]),
        saved_hb_scale_receipt=dict(
            source_key="actual_micro_mean_Z3_bounds.physical_scale_logs.hb",
            saved_tuple_is_loghb=True,
            physical_hb="exp(loghb)",
            log_W="2*log_Ra+loghb"),
        bridge_rates_receipt=rates,
        joint_background_source_receipt=dict(
            same_native_background_on_both_micro_charts=True,
            source_graph_functions_not_cap_values=True,
            six_signed_defect_handles=i1_contract["six_signed_defect_handles"],
            joint_endpoint_weights=correlated["viscous_endpoint_source_weights"]["joint"]),
        current_core_direction_receipt=True,
        current_core_fixed_point_receipt=True,
        source_executed=False), dict(
        axis=source_guard(root, EA_SOURCE, "source_jets", (
            '4*z[0]+j+sym(cap("Uz_minus_U0_absolute_bound",0,0))',
            '4+sym(cap("Uz_minus_U0_absolute_bound",0,1))')),
        bridge_densities=source_guard(root, BRIDGE_SOURCE, "densities", (
            "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",)),
        bridge_histories=source_guard(root, BRIDGE_SOURCE, "histories", (
            "for name,rate in RATES.items()",
            "original_incoming_memory_retained=True", "source_measure_applied_once=True")),
        bridge_rates=rates,
        loop_transport=source_guard(root, LOOP_SOURCE, "transport", (
            "for name,rate in RATES.items():",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "predecessor_memory_never_reset=True")),
        mean_phase=source_guard(root, MEAN_SOURCE, "_phase_rows", (
            '"phase_primitive_is_direct_ordinary_Z_row": True',
            '"source_operation": "phase_integral_of_slow_B_y" if source_name == "Cy" else "phase_integral_of_B"')),
        mean_moving=source_guard(root, MEAN_SOURCE, "_moving_rows", (
            "exponent = graph.neg(graph.mul(hb, distance))",
            "integrand = graph.mul(kernel, primitive_at_t, hb)",
            '"physical_hb_weight_applied_once": True')),
        mean_joined=source_guard(root, MEAN_SOURCE, "build_micro_mean", (
            'formula="exp(-hb*(s-1))"',
            "_scale_rows(g, first_JC_at_one, second_memory)",
            "_scale_rows(g, first_JCy_at_one, second_memory)",
            "second_complete = _add_rows(g, second_known, second_unresolved)")),
        six_source=source_guard(root, CORRELATED_SOURCE, "build", (
            "source=dict(q=Q,q_Z=Q_Z,q_squared=mul(Q,Q),U_q_Z=mul(U,Q_Z),",
            'weight=fixed(g.mul(radius,radius,chart["physical_measure"]))')),
        reduced_D=source_guard(root, CORRELATED_SOURCE, "build", (
            'D=add(neg(divide(integrals["q"])),half(divide(mul(h,z,integrals["q_Z"])))',
            'g.pscale(divide(mul(z,d,integrals["U_U_Z"]),L2),two))')),
        physical_q=source_guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "Q=div(numerator(U,m))", "Qy=div(numerator(Uy,my))")),
        core_inlet_average=source_guard(root, INLET_SOURCE, "__init__", (
            "sample=g.mul(rho,t)",
            "A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi)",
            "self.moment_functions={key:a.integral(body,t,g.zero,g.one,actual_core_normalized_moment=key,")),
        uniform_transfer=source_guard(root, UNIFORM_SOURCE, "run", (
            "derivative_rows=[]", "Uz_minus_U0_absolute_bound=eps*psi_norm*factor",
            "correction_C2=sum((eps*psi_norm*embedding(c,h,radius,0,k)")),
        current_core_source=source_guard(root, CORE_SOURCE, "__init__", (
            "sigma=g.quotient(self.parameters['j'],c(500),'positive fixed scalar 500')",
            "logLambda=g.add(g.mul(c(4),logP),c(1000))",
            "Lambda,epsilon=g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda))",
            "U=a.add(U0,a.scale(Psi,epsilon));mean=a.add(U0,a.scale(M,epsilon)")))


def _saved_source_rows(c, root, parent):
    """Hydrate signed native U and source-bound core correction radii."""
    b2 = parent["actual_micro_b_p2_mixed_bounds"]
    native = b2["native_angular_and_axial_provider"]
    U = [read(c, row) for row in native["native_U"][:4]]
    Uy = [read(c, row) for row in native["native_U_minus_core_inlet_absolute"][:4]]
    uniform = load(root / UNIFORM)
    require(uniform["implicit_source_sha256"] == parent["source_family"]["implicit_source_sha256"]
            and uniform["datum_enclosure_sha256"] == parent["source_family"]["datum_enclosure_sha256"]
            and uniform["C2_convention"] ==
            "sum of separate suprema of derivatives order0,1,2 (no factorial discount)"
            and uniform["full_real_axis_Uz_and_average_C2_input_closeness_certified"]
            and uniform["unique_analytic_field_bounded_not_sampled"],
            "Core uniform source correction certificate is required")
    by_order = {(row["scaled_radial_order"], row["axial_order"]): row
                for row in uniform["derivative_bounds"]}
    require(all((0, k) in by_order for k in range(4)),
            "Current core derivative rows through ordinary Z3 are required")
    core_caps = [read(c, by_order[(0, k)]["Uz_minus_U0_absolute_bound"]) for k in range(4)]
    require(all(value.a >= 0 for value in core_caps),
            "Core source correction rows must be nonnegative bounds")
    return U, Uy, core_caps


def _signed_background(c, root, parent):
    """Build signed U/M/Q source jets from saved enclosures only."""
    U, native_y_error, core_caps = _saved_source_rows(c, root, parent)
    direction = parent["actual_core_signed_direction"]
    delta = read(c, direction["parameters"]["delta"])
    j = read(c, direction["parameters"]["j"])
    Z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    baseline = [4 * Z + j, c.mpf(4), c.mpf(0), c.mpf(0)]

    # For k <= 2, the normalized core mean inlet is a radial average of the
    # same U correction.  The positive Duhamel kernel has total mass <= 1,
    # therefore core correction + native y correction is a valid radius around
    # the exact 4Z+j source.  M3 is only needed as a radius for Q_ZZ and comes
    # from the accepted native physical M3 norm row, centered at zero.
    errors = [core_caps[k].b + native_y_error[k].b for k in range(3)]
    m3_row = parent["actual_micro_b_p2_mixed_bounds"]["native_complete_history_rows"]["M"][3]
    m3 = read(c, m3_row)
    require(m3_row["lower_exact_mpf_tuple"][0] == 1
            and m3_row["upper_exact_mpf_tuple"][0] == 0,
            "Accepted native M3 row must be a symmetric norm radius")
    errors.append(c.fabs(m3).b)
    M = [baseline[k] + _symmetric(c, errors[k]) if k < 3 else _symmetric(c, errors[k])
         for k in range(4)]
    Q = _recover_Q(c, Z, U, M, delta)
    integrands = {
        "q": Q[:2],
        "q_Z": Q[1:3],
        "q_squared": jet_mul(c, Q[:2], Q[:2]),
        "U_q_Z": jet_mul(c, U[:2], Q[1:3]),
        "U_squared": jet_mul(c, U[:2], U[:2]),
        "U_U_Z": jet_mul(c, U[:2], U[1:3]),
    }
    require(all(len(value) == 2 for value in integrands.values()),
            "Every six source must have ordinary Z0/Z1 rows")
    return dict(Z=Z, delta=delta, j=j, baseline=baseline, U=U, M=M, Q=Q,
                integrands=integrands, core_correction=core_caps,
                native_y_error=native_y_error, M_error_radius=errors,
                M3_norm_radius=m3)


def _physical_log_scale(c, parent, loghb):
    logP = read(c, parent["actual_core_signed_direction"]["parameters"]["logPstar"])
    # epsilon = exp(-(4 log Pstar + 1000)); Ra=4 epsilon.  Keep this as a
    # log, since exp(logRa) is far below the representable physical scale.
    logRa = c.ln(c.mpf(4)) - (4 * logP + 1000)
    logW = 2 * logRa + loghb
    return dict(logPstar=record(logP), loghb=record(loghb), log_Ra=record(logRa), log_W=record(logW),
                W_definition="Ra^2*h_b", radius_definition="R=Ra*exp(h_b*s)",
                epsilon_definition="exp(-(4*logPstar+1000))",
                saved_hb_is_loghb=True, physical_hb_definition="exp(loghb)",
                log_W_definition="2*log_Ra+loghb",
                no_astronomical_exp_materialized=True)


def _normalized_charts(c, source, sc, loghb):
    """Integrate a uniform signed source interval against exp(2*h_b*s)."""
    require(loghb.b < (-1024 * c.ln(c.mpf(2))).a,
            "Saved loghb must support the normalized geometry remainder")
    # If h_b <= 2^-1024 and 0<=s<=2, x=2*h_b*s <=2^-1022.  For x<=1/4,
    # exp(x)<=1/(1-x)<=1+2x, hence exp(2*h_b*s)<=1+2^-1021.
    tiny = c.mpf(2) ** (-1021)
    geometry = c.mpf(1) + c.mpf([0, tiny])
    lengths = {"first_micro": c.mpf(1) - sc / 2, "second_micro": c.mpf(1)}
    result = {}
    for chart, length in lengths.items():
        require(length.a > 0, "Native chart length must be positive")
        weight = length * geometry
        result[chart] = dict(
            normalized_integrals={name: _rows([value * weight for value in jet])
                                  for name, jet in source["integrands"].items()},
            chart_length=record(length), geometry_factor=record(geometry),
            normalized_measure="integral exp(2*h_b*s) ds; h_b=exp(loghb)",
            background_is_signed_enclosure=True,
            no_midpoint_or_cap_selected_as_function_value=True)
    return result, dict(passed=True, loghb=record(loghb), saved_hb_field="loghb",
        physical_hb="h_b=exp(loghb)",
        geometry_factor=record(geometry), first_length="1-sc/2", second_length="1",
        first_lower="sc/2", seam="s=1", upper="2",
        geometry_remainder_proof="loghb<=-1024*log(2) => h_b<=2^-1024; exp(2*h_b*s)<=1+2^-1021 for 0<=s<=2",
        physical_measure_reduced_once=True)


def _joint_background(c, charts):
    """Sum the two signed chart enclosures before adding the joint error once."""
    result = {}
    for name in NAMES:
        rows = []
        for k in range(2):
            first = read(c, charts["first_micro"]["normalized_integrals"][name][k]["value"])
            second = read(c, charts["second_micro"]["normalized_integrals"][name][k]["value"])
            rows.append(dict(ordinary_Z_order=k, value=record(first + second)))
        result[name] = rows
    return result


def _finite_N_radius(c, parent, log_scale):
    """Map the accepted same-source physical error to normalized W units."""
    packet = parent["actual_micro_I1_defect_bounds"]
    require(packet["physical_I1_error_C1_below_2_power_minus1024"]
            and packet["physical_R_squared_and_hb_applied_once"]
            and packet["source_contract"]["original_physical_I1_measure"] ==
            "R^2*hb*ds exactly once",
            "Accepted physical I1 error packet is required")
    loghb = read(c, parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]["hb"])
    logW = read(c, log_scale["log_W"])
    threshold = -1024 * c.ln(c.mpf(2))
    ratios = {}
    for name in NAMES:
        rows = packet["six_joint_I1_log_caps"][name]
        ratios[name] = []
        for row in rows:
            cap = read(c, row["log_absolute_upper_bound"])
            ratio = cap - logW
            require(ratio.b < threshold.a,
                    "Saved physical six-I1 error does not normalize below 2^-1024")
            ratios[name].append(record(ratio))
    radius = c.mpf(2) ** (-1024)
    return dict(passed=True, normalized_error_radius=record(radius),
                normalized_error_radius_definition="2^-1024 in W=Ra^2*h_b units; saved hb field is loghb",
                source_log_cap_minus_log_W=ratios,
                all_log_caps_below_threshold=True,
                threshold=record(threshold),
                same_source_function_error_not_a_selected_value=True,
                physical_error_packet_not_recomputed=True)


def _add_radius(c, packed_rows, radius):
    return [dict(ordinary_Z_order=row["ordinary_Z_order"],
                 value=record(read(c, row["value"]) + _symmetric(c, radius)))
            for row in packed_rows]


def _reduced_micro_volume_D(c, source, joint_rows):
    """Apply the accepted reduced-volume six-I1 coefficient identity."""
    Z = source["Z"]
    delta = source["delta"]
    zero = c.mpf(0)
    one = c.mpf(1)
    two = c.mpf(2)
    z = [Z, one]
    one_jet = [one, zero]
    h = _jet_sub(c, one_jet, [delta, zero])
    d = _jet_sub(c, one_jet, jet_mul(c, z, z))
    L = _jet_sub(c, one_jet, [delta * value for value in jet_mul(c, z, z)])
    require(L[0].a > 0, "Reduced-volume coefficient denominator must be positive")
    inv_L = jet_reciprocal(c, L)
    inv_L2 = jet_mul(c, inv_L, inv_L)
    half = [c.mpf("1/2"), zero]
    coefficients = {
        "minus_inv_L": jet_neg(c, inv_L),
        "hZ_over_2L": jet_mul(c, half, jet_mul(c, jet_mul(c, h, z), inv_L)),
        "minus_half": jet_neg(c, half),
        "d_over_L": jet_mul(c, d, inv_L),
        "two_hZ2_over_L2": jet_mul(c, [two, zero],
                                    jet_mul(c, jet_mul(c, h, jet_mul(c, z, z)), inv_L2)),
        "two_Zd_over_L2": jet_mul(c, [two, zero],
                                   jet_mul(c, jet_mul(c, z, d), inv_L2)),
    }
    integrals = {
        name: [read(c, row["value"]) for row in rows]
        for name, rows in joint_rows.items()
    }
    require(set(integrals) == set(NAMES)
            and all(len(values) == 2 for values in integrals.values()),
            "Joint full finite-N six-I1 rows are required")
    terms = [
        jet_mul(c, coefficients["minus_inv_L"], integrals["q"]),
        jet_mul(c, coefficients["hZ_over_2L"], integrals["q_Z"]),
        jet_mul(c, coefficients["minus_half"], integrals["q_squared"]),
        jet_mul(c, coefficients["d_over_L"], integrals["U_q_Z"]),
        jet_mul(c, coefficients["two_hZ2_over_L2"], integrals["U_squared"]),
        jet_mul(c, coefficients["two_Zd_over_L2"], integrals["U_U_Z"]),
    ]
    D = terms[0]
    for term in terms[1:]:
        D = jet_add(c, D, term)
    half_D = jet_mul(c, half, D)
    return dict(
        passed=True,
        ordinary_Z_orders=[0, 1],
        normalized_values=_rows(D),
        half_normalized_values=_rows(half_D),
        coefficientsC1rows={name: _rows(values) for name, values in coefficients.items()},
        same_W=True,
        formula=("D=-I[q]/L+hZ*I[q_Z]/(2L)-I[q_squared]/2+"
                 "d*I[U_q_Z]/L+2*hZ^2*I[U_squared]/L^2+"
                 "2*Z*d*I[U_U_Z]/L^2"),
        coefficient_definitions=dict(
            h="1-delta", d="1-Z^2", L="1-delta*Z^2",
            hZ="(1-delta)*Z", L2="L^2"),
        source_integrals="joint_I1.normalized_full_finite_N_integrals",
        finite_N_error_consumed_once=True,
        background_and_finite_rows_are_signed_enclosures=True,
        no_outer_boundary_claim=True)


def symbolic_controls():
    """Exact source and Duhamel identities used by the adapter."""
    import sympy as sp

    z, de, U, M, Mz = sp.symbols("z delta U M Mz")
    L = 1 - de * z ** 2
    q = (2 * z * U - (1 - de) * z * M - (1 - z ** 2) * Mz) / L
    require(sp.cancel(q * L - (2 * z * U - (1 - de) * z * M - (1 - z ** 2) * Mz)) == 0,
            "Q source identity failed")
    h, s, v = sp.symbols("h s v", positive=True)
    mass = sp.integrate(h * sp.exp(-h * (s - v)), (v, 0, s))
    require(sp.simplify(mass - (1 - sp.exp(-h * s))) == 0,
            "Positive Duhamel mass identity failed")
    return dict(passed=True, Q_source_identity=True,
                mean_source_identity="M_y=U-M",
                positive_Duhamel_mass="1-exp(-h_b*s), h_b=exp(loghb)",
                saved_hb_field="loghb",
                ordinary_factorial_normalization=False)


def build(root, parent):
    """Build the signed normalized background packet from saved intervals."""
    root = Path(root)
    disk_parent = load(root / PARENT)
    require(sha(root / PARENT) == PARENT_SHA256 and parent == disk_parent,
            "Latest accepted RQy seam report is required")
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256)}
    for path in (B2_REPORT, I1_REPORT, UNIFORM, UNIFORM_SOURCE, LOWER_GRAPH,
                 LOWER_SOURCE, LOOP_GRAPH, LOOP_SOURCE, BRIDGE_SOURCE,
                 Q4_SOURCE, CORRELATED_SOURCE, JET_SOURCE, EA_SOURCE,
                 CORE_SOURCE, CORE_CHECK, INLET_SOURCE, BRIDGE_RATES_SOURCE,
                 BRIDGE_RATES_CHECK, MEAN_SOURCE):
        pins[path] = _pin(root, parent, path)
    contract, guards = _source_contract(root, parent)

    from mpmath.ctx_iv import MPIntervalContext
    c = MPIntervalContext()
    c.prec = 600
    source = _signed_background(c, root, parent)
    loghb = read(c, parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]["hb"])
    sc = _exact_tuple(c, parent["actual_native_collar_transport"]["domain"]["sc_exact_mpf_tuple"])
    scale = _physical_log_scale(c, parent, loghb)
    charts, geometry = _normalized_charts(c, source, sc, loghb)
    finite = _finite_N_radius(c, parent, scale)
    joint_background = _joint_background(c, charts)
    radius = read(c, finite["normalized_error_radius"])
    joint_finite = {
        name: _add_radius(c, rows, radius)
        for name, rows in joint_background.items()
    }
    reduced_D = _reduced_micro_volume_D(c, source, joint_finite)

    return dict(
        passed=True,
        source_family=parent["source_family"],
        source_pins=pins,
        source_contract=contract,
        source_AST=guards,
        symbolic_controls=symbolic_controls(),
        domain=dict(s=["sc/2", "1", "2"], charts=["first_micro", "second_micro"],
                    Z=["37/100", "93/250"], ordinary_Z_orders=[0, 1],
                    physical_y="log(R/Ra)=h_b*s; h_b=exp(loghb)"),
        source_axis=dict(U0="4*Z+j", physical_U="U0+epsilon*Psi",
                         baseline=_rows(source["baseline"]),
                         j=record(source["j"]), ordinary_Z_orders=[0, 1, 2, 3],
                         source_rows_are_signed_enclosures=True,
                         axis_reference_not_selected_as_actual_function=True),
        signed_source_rows=dict(
            U=_rows(source["U"]), M=_rows(source["M"]), Q=_rows(source["Q"]),
            integrands={name: _rows(values) for name, values in source["integrands"].items()},
            core_correction=_rows(source["core_correction"]),
            native_y_error=_rows(source["native_y_error"]),
            M_error_radius=[record(_symmetric(c, value)) for value in source["M_error_radius"]],
            M3_norm_radius=record(source["M3_norm_radius"]),
            M3_center_is_zero=True,
            M0_M1_M2_are_axis_plus_symmetric_error=True,
            signed_function_values_are_enclosures_not_cap_centers=True),
        physical_scale=scale,
        geometry=geometry,
        charts=charts,
        joint_I1=dict(
            normalized_background_integrals=joint_background,
            normalized_full_finite_N_integrals=joint_finite,
            background_is_signed_enclosure=True,
            finite_N_error_radius_is_symmetric=True,
            finite_N_error_added_once=True,
            same_physical_measure="R^2*h_b*ds exactly once (source receipt: R^2*hb*ds exactly once)"),
        reduced_micro_volume_D=reduced_D,
        finite_N_error=finite,
        actual_signed_background_I1=True,
        actual_signed_finite_N_I1_enclosures=True,
        normalized_by_W=True,
        W_definition="W=Ra^2*h_b with h_b=exp(loghb); physical I1=W*normalized_integral",
        finite_N_error_added_once_as_normalized_radius=True,
        accepted_ancestors_replayed=False,
        owner_builder_fixture_recurrence_norm_quadrature_replayed=False)
