"""New directed value leaf for the original core parameter and axial gauge.

No ancestral provider, recurrence, root/pole solver or quadrature is called.
The actual functions remain exact recipes. Directed ranges are approximations
with errors, never selected function values. This is a local prerequisite;
the nonlinear core, finite-N O3 ports and physical norms are not evaluated.
"""
from __future__ import annotations

from typing import Any, Callable

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


COMPONENT = "native_O3_anchor_parameter_axial_jet"
SUCCESS_COUNTS = {
    "new_parameter_selection_family_calls": 1,
    "new_scaled_monotone_root_family_calls": 1,
    "new_explicit_axial_gauge_jet_family_calls": 1,
    "new_mathematical_proofs": 1,
    "numeric_evaluation_calls": 1,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "source_function_port_calls": 0,
    "N_materialization_calls": 0,
    "failed_family_calls": 0,
}
SUCCESS_FLAGS = {
    "exact_selected_logC_dyadic_bound": True,
    "exact_j_minimum_branch_selected": True,
    "same_source_unique_anchor_enclosed": True,
    "anchor_log_amplitude_DZ_DZZ_defined": True,
    "rational_gauge_DZ_DZZ_enclosed": True,
    "numeric_evaluation": True,
    "actual_parameter_admission": False,
    "actual_source_Z_jets_admitted": False,
    "full_nonlinear_core_values_evaluated": False,
    "finite_N_O3_values_evaluated": False,
    "physical_layer_completed": False,
    "global_axis_regularity_proved": False,
    "global_phase_seam_regularity_proved": False,
    "full_goal_complete": False,
}
CALLBACK_EVENT_COUNTS = {
    "parameter_selection_entered": 1,
    "parameter_selection": 1,
    "anchor_root_entered": 1,
    "anchor_root": 1,
    "axial_jets_entered": 1,
    "axial_jets": 1,
    "proof_bundle": 1,
}
COUNT_EVENT_MAP = {
    "new_parameter_selection_family_calls": ["parameter_selection_entered"],
    "new_scaled_monotone_root_family_calls": ["anchor_root_entered"],
    "new_explicit_axial_gauge_jet_family_calls": ["axial_jets_entered"],
    "numeric_evaluation_calls": ["parameter_selection_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
    "old_math_calls": [],
    "old_provider_calls": [],
    "source_function_port_calls": [],
    "N_materialization_calls": [],
    "failed_family_calls": [],
}


class KernelInputError(ValueError):
    pass


def _need(condition: Any, message: str) -> None:
    if not condition:
        raise KernelInputError(message)


def derive(capsule: dict[str, Any], emit: Callable | None = None) -> dict[str, Any]:
    """Enter once, using only authenticated, frozen literal source data."""
    _need(capsule.get("schema") == COMPONENT + "_source_capsule_v1", "capsule schema")
    _need(capsule.get("status", {}).get("mathematical_functions_called") == 0, "capsule extraction")
    _need(capsule.get("precision_dps") == 260 and capsule.get("bisection_steps") == 640, "frozen precision")
    bindings = capsule["source_bindings"]
    _need(bindings["ordinary_Z_order_is_not_sample_coordinate"] is True, "ordinary Z orders")
    _need(bindings["source_H"] == "-4*Z**3-j*Z**2+(9-delta)*Z/2+j", "source H")
    _need(bindings["source_g"] == "(1-delta*Z**2)*H/(H**2+(j/500)**2)", "source g")
    _need(bindings["source_logF0"] == "-logC-Lambda*integral_anchor^Z(g)", "source amplitude")
    _need(bindings["j_recipe"] == "min(epsilon0/4,min(1/(8*CA**2*CQ),1/(20000*CA*KN*CS))/100)/8", "j recipe")
    _need(bindings["exact_logP_recipe"] == "exp(40)+11", "logP recipe")
    _need(bindings["exact_delta_recipe"] == "exp(-4*logP-30)", "delta recipe")
    _need(bindings["exact_logLambda_recipe"] == "4*logP+1000", "logLambda recipe")
    _need(capsule["common_fixed_point_admission"]["all_passed"] is True, "core admission")
    _need(capsule["common_fixed_point_admission"]["current_core_actual_operator_Banach_fixed_point_admitted"] is True, "actual operator")
    _need(capsule["selected_family_admission"]["uniform_analytic_fixed_point_admission_extended"] is True, "selected family")
    _need(capsule["family_identity"] == capsule["common_fixed_point_admission"]["family_identity"], "source family")
    _need(capsule["selected_family_admission"]["implicit_source_sha256"] == capsule["family_identity"]["implicit_source_sha256"], "pressure source")
    _need(capsule["selected_family_admission"]["datum_enclosure_sha256"] == capsule["family_identity"]["datum_enclosure_sha256"], "pressure datum")
    transfer_family = capsule["core_transfer_family_identity"]
    _need(all(transfer_family[key] == capsule["family_identity"][key] for key in ("implicit_source_sha256", "datum_enclosure_sha256")), "transfer source and datum")
    _need(transfer_family["analytic_core_family_sha256"] == capsule["selected_family_admission"]["base_analytic_core_family_sha256"], "selected base analytic family")
    c = MPIntervalContext()
    c.dps = 260

    def lo(v):
        return mp.make_mpf(v._mpi_[0])

    def hi(v):
        return mp.make_mpf(v._mpi_[1])

    def read(row):
        lower, upper = row["lower_exact_mpf_tuple"], row["upper_exact_mpf_tuple"]
        _need(len(lower) == len(upper) == 4, "encoded MPF endpoints")
        return c.mpf([mp.make_mpf(tuple(lower)), mp.make_mpf(tuple(upper))])

    def record(v, meaning="directed enclosure of exact recipe; neither endpoint nor midpoint selected as a function value"):
        # Exact endpoint tuples carry the certificate; decimals are displays.
        with mp.workdps(280):
            width = c.mpf(hi(v)) - c.mpf(lo(v))
            return {
                "lower_exact_mpf_tuple": list(v._mpi_[0]),
                "upper_exact_mpf_tuple": list(v._mpi_[1]),
                "lower": mp.nstr(lo(v), 70),
                "upper": mp.nstr(hi(v), 70),
                "absolute_width_upper": mp.nstr(hi(width), 30),
                "meaning": meaning,
            }

    def contained(outer, inner):
        return lo(outer) <= lo(inner) <= hi(inner) <= hi(outer)

    callback = emit if emit is not None else lambda event, row: None
    callback("parameter_selection_entered", {"phase": "before new numerical parameter expressions", "numeric_evaluation": True})
    const = capsule["fixed_integer_constants"]
    _need(const == {"CA": 344, "CQ": 1081, "CS": 57749553, "KN": 18000}, "fixed source constants")
    CA, CQ, CS, KN = (const[k] for k in ("CA", "CQ", "CS", "KN"))
    e_source = c.mpf(1) / (20000 * CA * KN * CS)
    e_contraction = c.mpf(1) / (8 * CA**2 * CQ)
    eps0 = read(capsule["epsilon0_enclosure"])
    _need(hi(e_source) < lo(e_contraction), "minimum e source branch")
    _need(hi(e_source / 100) < lo(eps0 / 4), "minimum eta branch")
    j_denominator = 16000000 * CA * KN * CS
    j = c.mpf(1) / j_denominator
    _need(contained(read(capsule["required_j_enclosure"]), j), "exact j must lie inside source rounding certificate")
    selected = capsule["selected_logCstar"]
    _need(selected["lower_exact_mpf_tuple"] == selected["upper_exact_mpf_tuple"], "logC is not singleton")
    _need(capsule["logC_selection_policy"] == "bind the already selected singleton dyadic parameter, not a source field enclosure", "explicit selection policy")
    logC = read(selected)
    logP = c.exp(c.mpf(40)) + 11
    logLambda = 4 * logP + 1000
    delta = c.exp(-4 * logP - 30)
    Lambda = c.exp(logLambda)
    epsilon = c.exp(-logLambda)
    sigma = j / 500
    _need(0 < lo(j) and hi(j) < lo(c.mpf("1e-20")), "j domain")
    _need(0 < lo(delta) and hi(delta) < lo(c.mpf("1e-200")), "positive source delta")
    _need(contained(read(capsule["logP_enclosure"]), logP), "source logP enclosure")
    _need(contained(read(capsule["logLambda_enclosure"]), logLambda), "source logLambda enclosure")
    _need(contained(read(capsule["Lambda_enclosure"]), Lambda), "source Lambda enclosure")
    _need(contained(read(capsule["epsilon_enclosure"]), epsilon), "source epsilon enclosure")
    parameter_rows = {
        "j_exact": {"numerator": 1, "denominator": j_denominator, "definition": "same original minimum, source branch proved strictly selected"},
        "j_enclosure": record(j),
        "logC_exact_dyadic": {"sign": selected["lower_exact_mpf_tuple"][0], "mantissa": selected["lower_exact_mpf_tuple"][1], "binary_exponent": selected["lower_exact_mpf_tuple"][2], "bitcount": selected["lower_exact_mpf_tuple"][3]},
        "logC_enclosure": record(logC, "exact selected singleton dyadic construction parameter; not a field-value interval or a new endpoint choice"),
        "logP": record(logP), "logLambda": record(logLambda),
        "delta": record(delta), "Lambda": record(Lambda), "epsilon": record(epsilon), "sigma": record(sigma),
        "Lambda_epsilon_exact_identity": "exp(logLambda)*exp(-logLambda)=1 for the same exact logLambda; interval products are not used to prove this identity",
        "branch_inequalities": {"e_source": record(e_source), "e_contraction": record(e_contraction), "epsilon0_over_4": record(eps0 / 4), "e_source_over_100": record(e_source / 100)},
        "parameter_rounding_enclosures_are_not_selected_values": True,
        "global_parameter_admission": False,
    }
    callback("parameter_selection", parameter_rows)

    callback("anchor_root_entered", {"phase": "before new scaled monotone root expressions", "old_root_calls": 0})
    # H(j*s)/j = 1 +(9-delta)*s/2 -j^2*(s^2+4*s^3).
    # Fresh dyadic bisection in scaled s; no old root, pole or centre used.
    def scaled_H(s):
        return 1 + (9 - delta) * s / 2 - j**2 * (s**2 + 4 * s**3)

    left, right = c.mpf(-1) / 4, c.mpf(-1) / 8
    _need(hi(scaled_H(left)) < 0 < lo(scaled_H(right)), "scaled root initial signs")
    scaled_cell = c.mpf([lo(left), hi(right)])
    monotone = (9 - delta) / 2 - j**2 * (2 * scaled_cell + 12 * scaled_cell**2)
    _need(lo(monotone) > 4, "scaled root monotonicity")
    path = []
    for _ in range(640):
        middle = (left + right) / 2
        _need(middle._mpi_[0] == middle._mpi_[1], "bisection midpoint must be exact dyadic")
        test = scaled_H(middle)
        if hi(test) < 0:
            left = middle
            path.append("L")
        elif lo(test) > 0:
            right = middle
            path.append("R")
        else:
            raise KernelInputError("scaled root sign unresolved at frozen precision")
    sroot = c.mpf([lo(left), hi(right)])
    anchor = j * sroot
    _need(-lo(j) <= lo(anchor) < hi(anchor) < 0, "original anchor bracket")
    _need(contained(read(capsule["original_anchor_certificate"]), anchor), "same-root certificate containment")
    left_sign, right_sign = scaled_H(left), scaled_H(right)
    _need(hi(left_sign) < 0 < lo(right_sign), "final root endpoint signs")
    whole_bracket = c.mpf([-hi(j), 0])
    original_derivative = (9 - delta) / 2 - 2 * j * whole_bracket - 12 * whole_bracket**2
    _need(lo(original_derivative) > 4, "original root uniqueness on [-j,0]")
    anchor_row = {
        "exact_definition": "a is the unique real root H(a)=0 in [-j,0] for the exact j and delta recipes",
        "scaled_equation": "H(j*s)/j=1+(9-delta)*s/2-j^2*(s^2+4*s^3)",
        "scaled_root": record(sroot), "anchor": record(anchor),
        "initial_scaled_bracket": ["-1/4", "-1/8"],
        "bisection_steps": 640, "dyadic_decision_path": "".join(path),
        "final_left_sign": record(left_sign), "final_right_sign": record(right_sign),
        "scaled_derivative_lower_certificate": record(monotone),
        "original_derivative_lower_certificate": record(original_derivative),
        "original_root_interval_used_only_for_posthoc_containment": True,
        "old_root_or_nominal_center_used_as_input_to_solver": False,
        "root_interval_midpoint_selected": False,
    }
    callback("anchor_root", anchor_row)

    callback("axial_jets_entered", {"phase": "before new explicit ordinary axial jet expressions", "old_provider_calls": 0})
    def explicit_jet(z, root=False):
        H = -4 * z**3 - j * z**2 + (9 - delta) * z / 2 + j
        Hp = -12 * z**2 - 2 * j * z + (9 - delta) / 2
        Hpp = -24 * z - 2 * j
        L = 1 - delta * z**2
        Lp, Lpp = -2 * delta * z, -2 * delta
        if root:
            # Correlated H(a)=0 is an exact defining identity, not a box
            # subtraction or cancellation of any inherited error term.
            g = c.mpf(0)
            gp = L * Hp / sigma**2
            gpp = (2 * Lp * Hp + L * Hpp) / sigma**2
        else:
            num, den = L * H, H**2 + sigma**2
            np = Lp * H + L * Hp
            npp = Lpp * H + 2 * Lp * Hp + L * Hpp
            dp, dpp = 2 * H * Hp, 2 * Hp**2 + 2 * H * Hpp
            g = num / den
            gp = np / den - num * dp / den**2
            gpp = npp / den - num * dpp / den**2 - 2 * np * dp / den**2 + 2 * num * dp**2 / den**3
        ell1, ell2 = -Lambda * g, -Lambda * gp
        row = {
            "Z": record(z), "is_exact_implicit_anchor": root,
            "g": record(g), "ordinary_DZ_g": record(gp), "ordinary_DZZ_g": record(gpp),
            "ordinary_DZ_logF0": record(ell1), "ordinary_DZZ_logF0": record(ell2),
            "F0_Z_over_F0": record(ell1), "F0_ZZ_over_F0": record(ell2 + ell1**2),
            "absolute_F0_value_available": root,
            "absolute_F0_value_materialized": False,
            "absolute_F0_exact_definition": "exp(-logC)" if root else None,
            "G_value": {"exact_integer": 0} if root else None,
            "logF0_value": record(-logC, "exact negative of selected dyadic logC at the implicit root; follows from the anchored-integral identity G(a)=0") if root else None,
            "anchor_zero_and_value_identity_basis": "H(a)=0 by the unique-root definition; G(a)=0 by identical integral endpoints; g(a)=0 since sigma>0. These are exact implicit-root identities, not interval-centre evaluations." if root else None,
            "ordinary_Z_derivatives_at_fixed_common_source_parameters": True,
            "ordinary_Z_orders_are_not_temporal_coefficient_orders": True,
        }
        if root:
            _need(lo(gp) > 0 and hi(ell2) < 0, "anchor log-amplitude strict curvature")
            row["log10_abs_F0_ZZ_over_F0"] = record(c.ln(-ell2) / c.ln(10))
            row["strict_positive_F0_defined_by_exponential"] = True
            row["strict_negative_logF0_curvature"] = True
        return row

    jets = {"anchor": explicit_jet(anchor, True)}
    for label, numerator in (("Z_minus_quarter", -1), ("Z_zero", 0), ("Z_plus_quarter", 1)):
        jets[label] = explicit_jet(c.mpf(numerator) / 4)
    jet_packet = {
        "jets": jets,
        "ordinary_Z_orders_enclosed": [0, 1, 2],
        "sample_coordinates_are_in_open_original_axial_domain": True,
        "nonanchor_G_or_absolute_amplitude_not_evaluated": True,
        "finite_N_or_nonlinear_component_ports_not_evaluated": True,
        "exact_formula_packet": {
            "H": bindings["source_H"], "L": "1-delta*Z**2", "sigma": "j/500",
            "g": bindings["source_g"],
            "g_prime_at_anchor": "L(a)*Hprime(a)/sigma**2",
            "g_second_at_anchor": "(2*Lprime(a)*Hprime(a)+L(a)*Hsecond(a))/sigma**2",
            "logF0_at_anchor": "-logC", "DZ_logF0_at_anchor": "0",
            "DZZ_logF0_at_anchor": "-Lambda*L(a)*Hprime(a)/sigma**2",
            "F0ZZ_over_F0": "(DZ_logF0)**2+DZZ_logF0",
            "source_logF0": bindings["source_logF0"],
        },
    }
    callback("axial_jets", jet_packet)
    proof_rows = [
        {"identity": "original minimum selects exact rational j", "evidence": "strict directed inequalities and source-bound constructor recipes; source intervals only certify containment"},
        {"identity": "logC is the already selected exact singleton dyadic parameter", "evidence": "identical raw endpoint tuples, source selection policy and selected analytic family admission"},
        {"identity": "unique source anchor at exact common parameters", "evidence": "monotone scaled cubic, exact dyadic bisection, endpoint signs and original bracket derivative; no midpoint selected"},
        {"identity": "ordinary axial gauge/log-amplitude jets", "evidence": "explicit quotient/product calculus and correlated exact H(a)=0; independent new directed arithmetic"},
    ]
    callback("proof_bundle", {"proof_rows": proof_rows, "proof_count": 1, "full_goal_complete": False})
    return {
        "schema": COMPONENT + "_kernel_result_v1", "component": COMPONENT,
        "passed": True, "metadata_only": False, "source_only": True,
        "numeric_evaluation": True, "proof_count": 1,
        "actual_counts": dict(SUCCESS_COUNTS), "new_flags": dict(SUCCESS_FLAGS),
        "payload": {"parameters": parameter_rows, "anchor": anchor_row, "axial_gauge_jets": jet_packet, "proof_rows": proof_rows,
                    "family_identity": capsule["family_identity"],
                    "shape_summary": {"sample_points": 4, "ordinary_gauge_derivative_orders": [0, 1, 2], "anchor_log_amplitude_orders": [0, 1, 2], "proof_rows": 4},
                    "remaining": ["full same-source nonlinear core values and axial jets", "nonanchor anchored G values", "native finite-N O3 values and source jets", "signed m15 endpoint identities", "physical velocity/vorticity/norms and all global construction layers"]},
        "full_goal_complete": False, "physical_layer_completed": False,
        "runtime_model_metadata": None,
    }
