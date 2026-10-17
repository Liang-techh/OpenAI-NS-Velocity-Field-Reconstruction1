"""Compute directed finite-word enclosures at s=sc*q/4 without ancestor replay.

Signed source-sector coefficients are multiplied and added numerically.
The positive exponential source factors remain formal. Local source errors
are separately computed raw-unit radii, never signed increment values.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import ast_guard, digest, load, require, sha

BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_compliant_current_original_n1_"
PARENT = PREFIX + "cutoff_source_binding.json"
PARENT_SHA256 = "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856"
COUPLED = PREFIX + "coupled_first_step.json.gz"
DOMAIN = PREFIX + "inner_analytic_domain.json.gz"
VOLTERRA_SOURCE = PREFIX + "inner_volterra_functions.py"
DOMAIN_SOURCE = PREFIX + "inner_analytic_domain.py"
COUPLED_SOURCE = PREFIX + "coupled_first_step.py"
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
SELECTED_Q = (Fraction(0), Fraction(1, 2), Fraction(5, 8), Fraction(3, 4), Fraction(1))
PRECISION_BITS = 1800


def interval(c, record):
    lo, hi = tuple(record["lower_exact_mpf_tuple"]), tuple(record["upper_exact_mpf_tuple"])
    require(len(lo) == len(hi) == 4, "Exact four-entry MPF endpoint tuples required")
    v = c.make_mpf((lo, hi))
    require(v.a <= v.b, "Reversed interval")
    return v


def endpoint(c, value, side=1):
    bits = value._mpi_[side]
    return c.make_mpf((bits, bits))


def record(value):
    return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),
        upper_exact_mpf_tuple=list(value._mpi_[1]), directed=True)


def positive_bound(c, saved, side=1):
    v = endpoint(c, interval(c, saved), side)
    require(v > 0, "Strictly positive bound required")
    return v


def box(c, upper):
    u = endpoint(c, upper)
    require(u >= 0, "Nonnegative radius required")
    return c.make_mpf(((0, 0, 0, 0), u._mpi_[1]))


def add_sector(c, result, key, value, log_value):
    if key not in result:
        result[key] = (value, log_value)
        return
    old_value, old_log = result[key]
    lo = old_log._mpi_[0] if old_log.a <= log_value.a else log_value._mpi_[0]
    hi = old_log._mpi_[1] if old_log.b >= log_value.b else log_value._mpi_[1]
    result[key] = (old_value + value, c.make_mpf((lo, hi)))


def core_fields(c, sample, sections, order):
    result = {f: {} for f in FIELDS}
    for section in sections:
        for field in FIELDS:
            for term in sample[section][field]["terms"]:
                key = tuple(term["source_exponents"])
                require(len(key) == 4, "Original four-sector basis required")
                value = interval(c, term["ordinary_Z_derivatives"][order])
                log_value = interval(c, term["combined_source_log"])
                add_sector(c, result[field], key, value, log_value)
    return result


def linear_combination(c, *pairs):
    result = {}
    for coefficient, fields in pairs:
        for key, (value, log_value) in fields.items():
            add_sector(c, result, key, coefficient * value, log_value)
    return result


def transport(c, fields, y_upper):
    zero, one = c.mpf(0), c.mpf(1)
    y = box(c, y_upper)
    def decay(factor):
        lower = endpoint(c, one - factor * y_upper, 0)
        require(lower > 0, "Microscopic transport outside admitted positive range")
        return c.make_mpf((lower._mpi_[0], one._mpi_[1]))
    a, a3, ah = decay(one), decay(c.mpf(3) / 2), decay(one / 2)
    return dict(
        F1=linear_combination(c, (one, fields["F1"]), (y, fields["d_x_F1"])),
        Uz1=linear_combination(c, (one, fields["Uz1"]), (y, fields["d_x_Uz1"])),
        K1=linear_combination(c, (a, fields["K1"]), (-y, fields["d_x_Uz1"])),
        P1=linear_combination(c, (one, fields["P1"])),
        d_x_F1=linear_combination(c, (a3, fields["d_x_F1"])),
        d_x_Uz1=linear_combination(c, (ah, fields["d_x_Uz1"])),
    )


def export_fields(fields):
    return {f: [dict(source_exponents=list(key), coefficient_interval=record(value),
        combined_source_log_interval=record(log_value),
        meaning="coefficient_interval * exp(same exact original source log)")
        for key, (value, log_value) in sorted(sectors.items())]
        for f, sectors in fields.items()}


def local_radii(c, qsc, domain, tail, order):
    # The saved Ai was constructed as a uniform forcing cap times the full
    # phase length 1/2. Replacing that length by q*sc/4 multiplies it by q*sc/2.
    Ai = domain["actual_collar_integrated_source_upper"]
    G = {i: qsc * positive_bound(c, Ai[str(i)]) / 2 for i in (4, 5, 6)}
    # dx=2*expm1(hb*q*sc/8) <= q*sc/4, since hb<1/2, sc<2^-10.
    dx = qsc / 4
    M = positive_bound(c, tail["initial_base_common_tube_bound"])
    Cp = positive_bound(c, tail["off_kinematic_B0_norm_upper"])
    C1 = positive_bound(c, tail["B1_norm_upper"])
    loss = positive_bound(c, tail["total_available_Z_loss"], 0)
    H = Cp * M + C1 * M / loss
    r0 = [dx * G[5], dx * G[6], dx * G[6], G[4], G[5], G[6]]
    r1 = [dx**2 * H, dx**2 * H, dx**2 * H, dx * H, dx * H, dx * H]
    gap = loss / 2
    # g is bounded on the middle source tube. W0 is bounded on eta/4;
    # h0=C0 W0+B1 W0_Z is bounded on eta/8 by spending loss=eta/8.
    # Cauchy over a disk of radius gap=eta/16 controls its whole derivative,
    # including derivatives of C0/B1, without a separate product approximation.
    if order == 1:
        r0 = [v / gap for v in r0]
        r1 = [v / gap for v in r1]
    return dict(W0={f: box(c, v) for f, v in zip(FIELDS, r0)},
        W1={f: box(c, v) for f, v in zip(FIELDS, r1)},
        prefix={f: box(c, a + b) for f, a, b in zip(FIELDS, r0, r1)})


def pin(root, parent, path):
    expected = PARENT_SHA256 if path == PARENT else parent["input_hashes"][path[len(BASE):]]
    require(sha(root / path) == expected, "Changed selected-word source: " + path)
    return dict(path=path, expected_sha256=expected, passed=True)


def build(root, parent=None):
    root = Path(root)
    require(sha(root / PARENT) == PARENT_SHA256, "Accepted cutoff source report changed")
    saved = load(root / PARENT)
    if parent is None:
        parent = saved
    require(digest(parent) == digest(saved), "Different parent mapping")
    require(parent["checks"]["passed"] is True, "Passed source report required")
    raw = parent["actual_first_raw_n1_function"]
    remainder = parent["actual_first_raw_word_remainder"]
    require(raw["actual_first_raw_function_bound"] is True and
        raw["actual_first_raw_core_entry_bound"] is True and
        raw["actual_first_raw_all_word_convergence_bound"] is True and
        remainder["certified_raw_preconditioned_tail_majorant"] is True,
        "Actual raw same-source function and tail binding required")
    paths = (PARENT, COUPLED, DOMAIN, VOLTERRA_SOURCE, DOMAIN_SOURCE, COUPLED_SOURCE)
    pins = {path: pin(root, parent, path) for path in paths}
    archive, domain = load(root / COUPLED), load(root / DOMAIN)
    sample = archive["samples"]["core_exit"]
    geometry = parent["actual_native_collar_transport"]["domain"]
    require(archive["source_family"] == domain["source_family"] == parent["source_family"] and
        sample["source_family"] == parent["source_family"], "Different actual source family")
    require(sample["endpoint"]["chart"] == "core" and
        sample["endpoint"]["coordinate"]["lower_exact_mpf_tuple"] == [0, 1, 2, 1] and
        sample["computed_preconditioned_indices"] == [0, 1] and
        sample["ordinary_Z_prefix_derivative_order"] == 1,
        "Saved actual core-exit finite words required")
    require(geometry["passed"] is True and geometry["first_window_only"] is True and
        geometry["sc_positive_and_less_than_2_power_minus10"] is True and
        geometry["hb_positive_less_than_half_source_margin"]["lower_exact_mpf_tuple"][1] > 0 and
        "0<hb<1/2" in geometry["rho_bound_proof"],
        "Exact positive hb/sc source geometry required")
    guards = dict(
        green=ast_guard(root, VOLTERRA_SOURCE, "green_apply", (
            "complete_history_not_reset_at_core_seam=True",
            "distance=g.mul(hb,s)",
            "measure=g.mul(c('1/2'),hb,x)",
            "kernels=(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)")),
        uniform_forcing_measure=ast_guard(root, DOMAIN_SOURCE, "domain", (
            "integrated={key:row.shift(1)*(c.sqrt(rho_max)/4) for key,row in gb.items()}",
            "actual_collar_integrated_source_upper=integrated_upper",
            "Gg_Z_radius=eta/4,solution_Z_radius=eta/8")),
        feedback=ast_guard(root, COUPLED_SOURCE, "feedback", (
            "value=base[field] if label=='B0' else op.dz(base[field])",
            "rhs[str(row)]=rhs[str(row)]+entry*value")),
    )
    c = MPIntervalContext()
    c.prec = PRECISION_BITS
    bits = tuple(geometry["sc_exact_mpf_tuple"])
    sc = c.make_mpf((bits, bits))
    require(sc > 0 and sc < c.mpf(1) / 1024, "Exact microscopic sc required")
    tail = sample["analytic_preconditioned_tail"]
    loss = positive_bound(c, tail["total_available_Z_loss"], 0)
    source_eta = positive_bound(c, domain["common_tubes"]["outer_Z_radius"], 0)
    require(loss == source_eta / 8 and
        interval(c, sample["derivative_tail"]["Cauchy_gap"]) == loss / 2 and
        interval(c, domain["common_tubes"]["solution_Z_radius"]) == loss and
        tail["all_k_derivative_words_allowed"] is True,
        "Common holomorphic source/solution radii required")
    inlets = {order: dict(W0=core_fields(c, sample, ("initial_function",), order),
        W1=core_fields(c, sample, ("B0_feedback", "B1_derivative_feedback"), order))
        for order in (0, 1)}
    rows = []
    for q in SELECTED_Q:
        qsc = sc * q.numerator / q.denominator
        y_upper = endpoint(c, qsc / 8)
        s = qsc / 4
        orders = {}
        for order in (0, 1):
            w0 = transport(c, inlets[order]["W0"], y_upper)
            w1 = transport(c, inlets[order]["W1"], y_upper)
            prefix = {f: linear_combination(c, (c.mpf(1), w0[f]), (c.mpf(1), w1[f]))
                for f in FIELDS}
            errors = local_radii(c, qsc, domain, tail, order)
            orders[str(order)] = dict(W0=export_fields(w0), W1=export_fields(w1),
                prefix=export_fields(prefix),
                local_raw_error_radii={name: {f: record(v) for f, v in fields.items()}
                    for name, fields in errors.items()},
                enclosure="signed homogeneous sector enclosure + symmetric local raw error radius",
                allword_tail_included=False)
        rows.append(dict(q=str(q), s_exact_selected_enclosure=record(s),
            s_definition="sc*q/4", y_upper=record(y_upper),
            eta=("1" if q <= Fraction(1, 2) else "1/2" if q == Fraction(5, 8) else "0"),
            ordinary_Z=orders))
    return dict(passed=True, source_family=parent["source_family"],
        parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        parent_mapping_sha256=digest(parent), source_pins=pins, source_AST=guards,
        source_only=False, actual_selected_finite_word_enclosures_evaluated=True,
        actual_directed_sector_coefficients_and_local_radii_evaluated=True,
        physical_exponential_source_factors_materialized=False,
        ordinary_Z_orders=[0, 1], selected_coordinate_rows=rows,
        selected_Z_interval=sample["Z"], all_Z_functional_coverage=False,
        source_log_basis_names=sample["source_log_basis_names"],
        source_log_bases=sample["source_log_bases"],
        precision_bits=PRECISION_BITS,
        raw_local_error_proof=dict(
            Ai="saved uniform raw forcing cap times full phase length1/2",
            local_integrated_g="G_i=(q*sc/2)*Ai",
            dx_upper="dx=2*expm1(hb*q*sc/8)<=q*sc/4",
            local_measure="integral_0^s (hb*x(v)/2)dv=x(s)-2<=q*sc/4 by exact FTC",
            local_W0="[dx*G5,dx*G6,dx*G6,G4,G5,G6]",
            local_W1="[dx^2*H,dx^2*H,dx^2*H,dx*H,dx*H,dx*H]",
            H="Cp*Mbase+C1*Mbase/loss on Z radius eta/8",
            Z1="divide whole holomorphic local source bound by gap=eta/16",
            kernel_bounds="|K|,|u*log(X/u)|<=X-u<=dx; derivative and pressure kernels<=1",
            kinematic_transport="exact T(X,2), no seam reset or frozen core kernel"),
        core_quadrature_and_source_interval_errors_retained=True,
        old_actual_first_keep_point_values_used=False,
        old_core_polynomial_promoted_to_word_prefix=False,
        raw_allword_tail_kept_separate=True,
        useful_component_relative_accuracy_certified=False,
        signed_local_collar_integral_evaluated=False,
        signed_C_inner_m15_value_evaluated=False, numeric_m15_error_budget_instantiated=False,
        old_caps_or_six_signed_terms_removed=False, global_five_moment_closure=False,
        source_executed=False, accepted_ancestors_replayed=False)
