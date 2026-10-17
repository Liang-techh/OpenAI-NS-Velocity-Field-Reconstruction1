"""Bind the actual inner m15 point oracle to the independent pulse target.

Only saved admitted sources and fresh directed scalar arithmetic are used.
The shared pulse cancels as a function before the discrepancy is boxed.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
import mpmath as mp
import sympy as sp

BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "n1_signed_C_inner_m15.json"
PARENT_SHA256 = "b6c7aa8f10dcc3c1310949a030304da3bc2f54ce886439f959128475f210311c"
FULL = PREFIX + "n1_full_five_numeric_bridge.json.gz"
PULSE = PREFIX + "n1_pulse_numeric_bridge.json.gz"
FULL_CHECK = PREFIX + "n1_full_five_numeric_bridge_check.json"
PULSE_CHECK = PREFIX + "n1_pulse_numeric_bridge_check.json"
JOINT = PREFIX + "joint_reserved_m15_numeric.json"
FULL_SOURCE = PREFIX + "n1_full_five_numeric_bridge.py"
PULSE_SOURCE = PREFIX + "n1_pulse_numeric_bridge.py"
JOINT_SOURCE = "src/openai_ns_reconstruction/original_joint_reserved_m15.py"
NATIVE_SOURCE = "src/openai_ns_reconstruction/original_native_pulse_m15.py"
HELPERS = tuple("src/openai_ns_reconstruction/" + n + ".py" for n in (
    "original_finite_N_signed_source", "original_n1_selected_word_numeric"))
PRECISION_BITS = 1800
DYADIC_BITS = 4096


def exact_algebra():
    inner, lower, pre, pulse, unit = sp.symbols("inner lower pre pulse unit")
    i1, l1, r1, p1, u1 = sp.symbols("i1 l1 r1 p1 u1")
    full = inner + lower + pre + pulse
    full1 = i1 + l1 + r1 + p1
    checks = dict(
        shared_pulse_value=sp.expand(-full + pulse + inner + lower + pre) == 0,
        shared_pulse_ordinary_Z1=sp.expand(-full1 + p1 + i1 + l1 + r1) == 0,
        fixed_unit_value=sp.cancel((-full + pulse) / unit + (inner + lower + pre) / unit) == 0,
        fixed_unit_Z1=sp.cancel((-full1 + p1) / unit + (i1 + l1 + r1) / unit) == 0,
        quotient_unit_Z1_zero=sp.cancel(((i1 * unit - inner * u1) / unit**2).subs(u1, 0) - i1 / unit) == 0,
        reserved_is_target_plus_joint=sp.expand(-pulse - (inner + lower + pre) + full) == 0)
    if not all(checks.values()):
        raise ValueError("Reserved/pulse source algebra failed")
    return dict(passed=True, checks=checks)


def build(root, parent):
    # Pin every helper before trusting its saved-record or source-gating code.
    root = Path(root)
    sha0 = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if sha0(root / PARENT) != PARENT_SHA256:
        raise ValueError("Accepted actual inner point parent changed")
    pinned_parent = json.loads((root / PARENT).read_bytes())
    if pinned_parent != parent or parent["checks"]["passed"] is not True:
        raise ValueError("The exact accepted parent mapping is required")
    for path in HELPERS:
        if sha0(root / path) != pinned_parent["implementation_hashes"][path]:
            raise ValueError("Executable helper changed: " + path)
    from .original_finite_N_signed_source import ast_guard, digest, load, require, sha
    from .original_n1_selected_word_numeric import interval, endpoint, record

    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256, passed=True)}
    for path in (FULL, PULSE, FULL_CHECK, PULSE_CHECK, JOINT,
                 FULL_SOURCE, PULSE_SOURCE, JOINT_SOURCE, NATIVE_SOURCE, *HELPERS):
        expected = (parent["input_hashes"][Path(path).name] if path.startswith(BASE)
                    else parent["implementation_hashes"][path])
        require(sha(root / path) == expected, "Changed reserved point source: " + path)
        pins[path] = dict(path=path, expected_sha256=expected, passed=True)
    full, pulse, fc, pc, joint = [load(root / path) for path in (FULL, PULSE, FULL_CHECK, PULSE_CHECK, JOINT)]
    require(all(v["source_family"] == parent["source_family"] for v in (full, pulse, fc, pc, joint)),
            "One original source family required")
    for receipt, archive, gate in ((fc, FULL, "current_original_true_complete_five_defect_C1_numeric_bridge_installed"),
                                   (pc, PULSE, "current_original_true_pulse_pressure_C1_numeric_bridge_installed")):
        require(receipt["all_passed"] is True and receipt[gate] is True and
                receipt["input_hashes"][Path(archive).name] == pins[archive]["expected_sha256"],
                "Checked actual complete source enclosure required")
    require(fc["actual_physical_log_units_bound_to_source_recipes"] is True and
            fc["full_five_enclosures_are_unrepaired_defects_not_terminal_zeros"] is True and
            fc["exact_factored_ordinary_derivative_proof"]["passed"] is True and
            pc["exact_pressure_and_unit_proof"]["passed"] is True and
            pc["old_C1_double_normalization_error_corrected"] is True and
            pc["exact_pressure_and_unit_proof"]["all_remainders_and_finite_endpoints_are_enclosures_not_zero_functions"] is True,
            "Physical scale, ordinary derivative and finite-error admission required")
    src = full["source_identity"]
    require(all(src[k] is True for k in ("passed", "exact_five_signed_function_composition_replayed",
        "inner_lower_preRp_pulse_pressure_partitions_complete", "preRp_pulse_same_original_positive_reference_A_Rv_units",
        "physical_pressure_unit_logs_use_combined_inverse_mu_powers", "pressure_unit_Z1_exact_zero",
        "factored_source_Z1_already_includes_F0_relative_amplitude_derivatives", "m15_minus_already_in_original_signed_components")),
        "Admitted complete m15 source partition and fixed unit required")
    ps = pulse["source_identity"]
    require(all(ps[k] is True for k in ("passed", "exact_main_gap_end_pressure_C1_graph_and_endpoints_replayed",
        "all_three_charts_and_native_mean_histories_retained", "Pstar_applied_once_and_signed_m15_once",
        "all_nonzero_quadratic_finite_mu_and_weighted_viscous_errors_retained",
        "reference_is_positive_and_used_once", "energy_epsilon_already_in_reference_A_Rv_pressure_units",
        "corrected_C1_error_does_not_multiply_A_or_divide_reference_again")),
        "Complete pulse m15 with all nonzero remainders required")
    point, target_point = full["samples"]["point"], pulse["samples"]["point"]
    require(target_point["ordinary_Z0_Z1_enclosures_of_true_normalized_functions"] is True and
            target_point["complete_history_and_finite_mu_endpoints_retained"] is True and
            target_point["C1_energy_normalization_double_scaling_corrected"] is True and
            target_point["inner_lower_and_preRp_pressure_not_added"] is True and target_point["active_cells"] == [2],
            "Actual selected complete-pulse C1 source cell required")
    scaling = parent["core_pressure_source_identity"]["scaling"]
    require(scaling["epsilon_core"] == "exp(-logLambda)" and
            scaling["canonical"] == "physical_P0=Pstar^2*p" and
            scaling["forward"] == "P0_Z_taylor=epsilon_core*physical_P0_Taylor" and
            scaling["pressure_datum_epsilon_is_distinct"] is True and
            all(scaling["source_scaling_flags"][key] is True for key in (
                "Lambda_times_core_epsilon_exactly_one",
                "P0_step_equals_core_epsilon_times_Pstar_squared_times_normalized_pressure",
                "no_second_epsilon_applied_to_ell",
                "core_radial_epsilon_distinct_from_pressure_datum_parameter_epsilon")),
            "Exact core epsilon/Lambda and pressure-datum distinction required")
    actual = parent["actual_signed_C_inner_m15"]
    c = MPIntervalContext(); c.prec = PRECISION_BITS
    z, zp, zi = [interval(c, row) for row in (point["Z"], target_point["Z"], actual["selected_Z_interval"])]
    q = Fraction(371, 1000)
    require(z._mpi_ == zp._mpi_ and zi.a <= z.a <= z.b <= zi.b,
            "The common pulse point tube must be inside the actual inner tube")
    # Compare rational MPF endpoints exactly, without rounding decimal 0.371.
    def rational(bits):
        sign, mantissa, exponent, _ = bits
        value = Fraction((-1 if sign else 1) * mantissa)
        return value * 2**exponent if exponent >= 0 else value / 2**(-exponent)
    require(rational(z._mpi_[0]) <= q <= rational(z._mpi_[1]),
            "The common point tube must enclose exact Z=371/1000")
    require(actual["passed"] is True and actual["actual_C_inner_m15_value_enclosed"] is True and
            actual["C_inner_source_handles"] == [53161, 53162] and actual["complete_first_collar_added_once"] is True,
            "Same existing actual inner m15 function required")
    cached = joint["cached_component_source"]
    require(cached["passed"] is True and cached["report"] == Path(FULL).name and
            cached["report_sha256"] == pins[FULL]["expected_sha256"] and
            cached["source"] == Path(FULL_SOURCE).name and
            cached["source_sha256"] == pins[FULL_SOURCE]["expected_sha256"] and
            cached["source_errors_and_amplitude_derivatives_preserved"] is True and
            cached["unit_Z_derivative_exact_zero"] is True,
            "Joint identity must bind the identical admitted complete source")
    pair = joint["joint_reserved_minus_target_source"]
    expected_pairs = dict(inner=[53161, 53162], lower=[36592, 36593], preRp=[31221, 31220],
                          pulse=[35556, 35557], full=[53175, 53176])
    require(joint["checks"]["passed"] is True and pair["source_C1_pairs"] == expected_pairs and
            pair["source_identity_before_boxing"] is True and pair["independent_target_consumed_once"] is True and
            pair["exact_response"] == [0, 0, 1] and pair["selected_compact_x"] == "2",
            "Accepted exact reserved full-support joint source required")
    shared = pair["common_normalization"]
    require(shared["log_node"] == 53421 and shared["unit_node"] == 53427 and shared["square_node"] == 53430 and
            shared["U_Z_exact_zero"] is True and shared["ordinary_C1_is_scaled_physical_derivative"] is True,
            "Ordinary derivative in the shared fixed physical m15 unit required")
    nodes = pulse["normalization_graph_nodes"]; start = pulse["normalization_graph_span"][0]
    require(nodes[53421-start] == dict(operation="sum", arguments=[53409, 5, 53401, 31320, 53418, 53420, 31343]) and
            nodes[53427-start] == dict(operation="analytic_unary", name="exp", argument=53421) and
            nodes[53430-start] == dict(operation="product", arguments=[53427, 53427]),
            "Shared normalization graph changed")
    source_logs = parent["actual_correlated_leading_pressure_increment"]["first_four_source_log_bases"]
    require(actual["source_log_basis_names"] == ["hb", "F0_anchor", "Pstar_squared", "Lambda", "sigma_selected"],
            "Actual inner scale basis changed")
    unit_log = interval(c, point["actual_physical_pressure_unit_logs"]["m15"])
    local_log = interval(c, source_logs[2]) - 2 * interval(c, source_logs[3])
    ratio_log = local_log - unit_log
    parts = point["pressure_component_enclosures"]["m15"]
    for name in ("inner", "lower"):
        require(all(parts[name][k] is True for k in ("original_source_bases_are_fixed_within_Taylor_jets",
            "source_ordinary_Z1_already_contains_F0_amplitude_derivatives", "positive_unit_is_Z_independent",
            "nonzero_sector_not_selected_zero")), "Saved physical component derivative admission required")
    require(parts["pulse"] == target_point["normalized_true_pulse_C1_enclosures"]["m15"],
            "The independent complete pulse must be the identical source in the full partition")
    rows = {}
    dyadic = c.mpf(2) ** -DYADIC_BITS
    for k in (0, 1):
        saved = actual["ordinary_Z"][str(k)]["actual_C_inner_m15"]["normalized"]
        require(saved["normalizing_positive_factor"] == "Pstar^2*epsilon^2" and saved["strict_sign"] == "negative" and
                saved["upper_bounds_not_used_as_signed_values"] is True and
                saved["dominant_source_exponents"] == [0, 0, 1, -2, 0], "Signed local inner scale sector required")
        local = interval(c, saved["normalized_integral_enclosure"])
        require(local.b < 0, "Actual inner ordinary derivative must be strictly negative")
        log_abs = ratio_log + c.ln(abs(local))
        require(log_abs.b < -DYADIC_BITS * c.ln(2), "Inner value exceeds its admitted common-unit cover")
        # This is a one-sided enclosure of the negative function, never its value.
        common_inner = c.make_mpf(((-dyadic)._mpi_[0], (0, 0, 0, 0)))
        lower = interval(c, parts["lower"]["C1_enclosures"][k])
        pre = interval(c, parts["preRp"][k])
        pulse_value = interval(c, target_point["normalized_true_pulse_C1_enclosures"]["m15"][k])
        target = -pulse_value
        discrepancy = -(common_inner + lower + pre)
        reserved = target + discrepancy
        require(target.a > 0 and reserved.a > 0, "Independent target and reserved response signs unresolved")
        relative = abs(discrepancy) / target
        ctx = mp.mp.clone(); ctx.dps = 25
        display = lambda v: [ctx.nstr(ctx.make_mpf(t), 20) for t in v._mpi_]
        rows[str(k)] = dict(actual_inner_local_enclosure=record(local),
            inner_common_unit_log_absolute_enclosure=record(log_abs),
            inner_common_unit_enclosure=record(common_inner), inner_strict_sign="negative",
            inner_common_cover_dyadic="2^-" + str(DYADIC_BITS),
            lower_common_unit_enclosure=parts["lower"]["C1_enclosures"][k],
            preRp_common_unit_enclosure=parts["preRp"][k],
            complete_pulse_common_unit_enclosure=target_point["normalized_true_pulse_C1_enclosures"]["m15"][k],
            independent_negative_complete_pulse_target=record(target),
            reserved_full_support_response=record(reserved), correlated_reserved_minus_target=record(discrepancy),
            relative_correlated_discrepancy_enclosure=record(relative),
            display=dict(target=display(target), reserved=display(reserved), discrepancy=display(discrepancy), relative=display(relative)),
            cap_or_midpoint_selected_as_value=False, interval_overlap_used_as_closure=False)
    guards = dict(
        full_fixed_units=ast_guard(root, FULL_SOURCE, "sector_C1_in_fixed_pressure_unit", (
            "bound=log-unit_log+c.ln(c.mpf(magnitude))", "source_ordinary_Z1_already_contains_F0_amplitude_derivatives=True",
            "positive_unit_is_Z_independent=True", "nonzero_sector_not_selected_zero=True")),
        pulse_inlet=ast_guard(root, NATIVE_SOURCE, "inlet", (
            "normalized_true_pulse_C1_enclosures", "value=-self.J(self.iv,[self.native._decode(v) for v in rows])")),
        joint_source=ast_guard(root, JOINT_SOURCE, "negate_component_sum", (
            "for name in ('inner','lower','preRp')", "for k in (0,1)")))
    return dict(passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        parent_mapping_sha256=digest(parent), source_family=parent["source_family"], source_pins=pins,
        source_AST=guards, exact_algebra=exact_algebra(), ordinary_Z_orders=[0, 1], ordinary_Z=rows,
        precision_bits=PRECISION_BITS, selected_Z_interval=point["Z"], inner_enclosure_Z_interval=actual["selected_Z_interval"],
        common_point_tube_contained_in_inner_tube=True, exact_rational_Z="371/1000",
        shared_unit=dict(physical="reference*A*Rv^2", **shared),
        local_inner_unit="epsilon^2*Pstar^2", local_inner_unit_log=record(local_log),
        exact_core_pressure_scaling_checked=scaling, local_dominant_source_exponents_checked=[0, 0, 1, -2, 0],
        source_admission=dict(full_receipt=FULL_CHECK, pulse_receipt=PULSE_CHECK,
            checked_receipts_not_pre_admission_archive_flags=True, pulse_selected_active_cells=[2],
            pulse_is_actual_complete_finite_N_function_enclosure=True,
            pulse_nonzero_quadratic_C1_and_viscous_endpoint_errors_retained=True,
            cached_joint_full_source_hash_matches=True, positive_reference_and_energy_epsilon_applied_once=True),
        common_m15_unit_log=record(unit_log), local_to_common_unit_log_ratio=record(ratio_log),
        source_C1_pairs=expected_pairs, source_identity="reserved-(-pulse)=-(inner+lower+preRp)",
        independent_pulse_target_enclosure_admitted=True, reserved_response_enclosure_admitted=True,
        pulse_function_canceled_before_boxing=True, independent_target_consumed_once=True,
        inner_evaluation_replaced_once_in_new_point_view=True, lower_preRp_saved_errors_preserved=True,
        ordinary_Z1_is_physical_derivative_divided_by_fixed_unit=True,
        F0_connection_already_in_source_coefficients=True, inner_common_unit_cover_is_not_a_function_value=True,
        old_ledger_or_error_caps_changed=False, old_six_signed_contributions_retained=True,
        reserved_response_full_support_only=True, forward_Rp_discrepancy_enclosed=False,
        forward_six_source_function_identity=False,
        one_point_only=True, discrepancy_zero_identity_proved=False, independent_terminal_matching_closed=False,
        all_Z_functional_coverage=False, global_five_moment_closure=False, full_goal_complete=False,
        source_executed=False, accepted_ancestors_replayed=False, P1_is_distinct_from_m15=True)
