"""Complete-support cutoff measure, with the exact microscopic width factored.

This evaluates a positive integration weight, never a signed field value.
The core contributes exactly 4 in rho; the first collar contributes hb*sc*C.
"""
from __future__ import annotations

from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import ast_guard, digest, load, require, sha

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
PARENT_SHA256 = "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856"
FLAT = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
CUTOFF = BASE + "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py"


def build(root, parent):
    root = Path(root)
    require(sha(root / PARENT) == PARENT_SHA256, "Accepted cutoff source parent changed")
    require(digest(parent) == digest(load(root / PARENT)), "Different parent mapping")
    require(parent["checks"]["passed"] is True, "Passed source parent required")
    domain = parent["actual_native_collar_transport"]["domain"]
    require(domain["passed"] is True and
        domain["sc_positive_and_less_than_2_power_minus10"] is True and
        domain["ordinary_Z_independent_endpoints"] is True and
        parent["actual_n1_cutoff_history_identity"]["passed"] is True,
        "Actual Z-independent source geometry and cutoff required")
    width_margin = domain["hb_positive_less_than_half_source_margin"]["lower_exact_mpf_tuple"]
    require(width_margin[0] == 0 and width_margin[1] > 0 and
        "0<hb<1/2" in domain["rho_bound_proof"],
        "Pinned positive source width and hb<1/2 proof required")
    pins = {}
    for path in (PARENT, FLAT, CUTOFF):
        expected = (PARENT_SHA256 if path == PARENT else parent["input_hashes"][path[len(BASE):]])
        require(sha(root / path) == expected, "Cutoff measure source changed: " + path)
        pins[path] = dict(path=path, expected_sha256=expected, passed=True)
    guards = dict(
        cutoff=ast_guard(root, CUTOFF, "original_cutoff_rows", (
            "arg=g.sub(g.mul(c(4),q),c(2));sigma=own.jets.complete.Gamma.sigma(g,arg)",
            "value=g.sub(g.one,sigma)", "Z_independent_cutoff=True")),
        interior_sigma=ast_guard(root, FLAT, "_sigma_left", (
            "odds=1/(1-x)**2-1/x**2", "e=positive_exp(c,odds);value=e/(1+e)")),
        reflected_sigma=ast_guard(root, FLAT, "sigma_jets", (
            "left=_sigma_left(c,1-c.mpf([a,b]))", "1-left[0]")),
    )
    bits = tuple(domain["sc_exact_mpf_tuple"])
    require(bits[0] == 0 and bits[1] > 0 and bits[1].bit_length() == bits[3] and
        bits[2] + bits[3] <= -10, "Exact positive microscopic phase required")
    c = MPIntervalContext()
    c.dps = 100
    sc = c.make_mpf((bits, bits))
    lower = c.mpf(5) / 8
    upper = lower / (1 - c.mpf(3) * sc / 32)
    lo, hi = lower._mpi_[0], upper._mpi_[1]
    require(lo[0] == 0 and hi[0] == 0 and c.mpf(3) * sc / 32 < c.mpf(1) / 2,
        "Positive directed complete-support weight required")
    return dict(
        passed=True, source_family=parent["source_family"],
        parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        parent_mapping_sha256=digest(parent), source_pins=pins, source_AST=guards,
        selected_coordinate="s=sc*q/4", rho="4*exp(hb*sc*q/4)",
        eta="1-sigma(4*q-2)", integration_q_interval=["0", "1"], eta_support_q=["0", "3/4"],
        exact_selected_eta={"0": "1", "1/2": "1", "5/8": "1/2", "3/4": "0", "1": "0"},
        source_symmetry="sigma(t)+sigma(1-t)=1, including flat support endpoints",
        exact_unweighted_q_mass="1/2+(1/4)*(1/2)=5/8",
        first_measure="J_first=int_0^1 eta*d-rho=hb*sc*C",
        coefficient_lower="5/8",
        coefficient_upper_formula="(5/8)/(1-3*sc/32)",
        coefficient_bounds=dict(lower_exact_mpf_tuple=list(lo), upper_exact_mpf_tuple=list(hi),
            display=str(c.make_mpf((lo, hi))), directed_decimal_precision=100),
        inequality_proof=("0<=q<=3/4, 0<hb<1/2, so 0<=hb*sc*q/4<=3*sc/32; "
            "1<=exp(x)<=1/(1-x), and eta>=0"),
        complete_first_collar_support_included=True,
        core_measure="J_core=int_0^4 1*d-rho=4",
        complete_core_plus_first_measure="J_total=4+hb*sc*C",
        measure_is_Z_independent=True, source_sc_exact_mpf_tuple=list(bits),
        hb_kept_exact_formal_source_factor=True,
        scalar_weight_bound_is_not_signed_m15_value=True,
        numeric_m15_error_budget_instantiated=False,
        source_executed=False, accepted_ancestors_replayed=False,
        old_error_caps_removed=False, global_five_moment_closure=False,
    )
