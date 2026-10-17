"""Signed actual core direction from a new model sum and saved perturbation.

No accepted owner is executed.  The fixed-point error remains a symmetric
error around the directed model enclosure, never a selected function value.
"""
from __future__ import annotations
import ast
import hashlib
import json
from pathlib import Path

BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_"
PARENT = PREFIX + "compliant_current_original_n1_micro_phase_mean.json"
PARENT_SHA256 = "0005ff381670e7d4f453d43cedfcc1fcc42fd6f5ffee0b55f58c7af118b0d14c"
HELPER = "src/openai_ns_reconstruction/original_core_signed_model.py"
HELPER_SHA256 = "d5288dec6aa5a4a7fe1c7017973c432e53f6ada960183868b86e15dc0a4ba831"
UNIFORM = PREFIX + "compliant_core_uniform_bounds.json"
TRANSFER = PREFIX + "compliant_core_transfer.json"
TUBE = PREFIX + "shared_analytic_tube.json"
CORE_SOURCE = PREFIX + "compliant_current_original_n1_core_fixed_point_functions.py"
CORE_CHECK = PREFIX + "compliant_current_original_n1_core_fixed_point_functions_check.json"
BRIDGE_SOURCE = PREFIX + "compliant_current_original_prescribed_bridge_C6_functions.py"


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def record(v):
    from mpmath import make_mpf, nstr
    return dict(lower_exact_mpf_tuple=list(v._mpi_[0]), upper_exact_mpf_tuple=list(v._mpi_[1]),
                lower_display=nstr(make_mpf(v._mpi_[0]), 24), upper_display=nstr(make_mpf(v._mpi_[1]), 24),
                directed=True)


def read(c, row):
    from mpmath import make_mpf
    return c.mpf([make_mpf(tuple(row[k + "_exact_mpf_tuple"])) for k in ("lower", "upper")])


def source_guard(root, path, function, snippets):
    text = (root / path).read_text(encoding="utf-8")
    tree = ast.parse(text)
    candidates = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == function]
    # The core/bridge constructor is selected by its defining source text.
    candidates = [n for n in candidates if all("".join(s.split()) in
                  "".join(ast.get_source_segment(text, n).split()) for s in snippets)]
    require(len(candidates) == 1, "Changed source definition: " + path + ":" + function)
    return dict(passed=True, path=path, function=function, required_text=list(snippets),
                AST_sha256=hashlib.sha256(ast.dump(candidates[0], include_attributes=False).encode()).hexdigest(),
                source_executed=False)


def enclose_actual(c, model, uniform, rho):
    """Add the saved, source-bound derivative error once, with both signs."""
    actual = {}; rows = []
    bounds = {(r["scaled_radial_order"], r["axial_order"]): r for r in uniform["derivative_bounds"]}
    for key, v in model["derivatives"].items():
        i, k = key
        error = read(c, bounds[key]["Phi_minus_model_absolute_bound"])
        require(error.a >= 0, "Nonnegative actual-model norm bound required")
        # Only the upper endpoint is a bound.  The lower endpoint of the
        # bound record is not a lower bound on the size of the true error.
        symmetric = c.mpf([-error.b, error.b])
        value = v + symmetric
        actual[key] = value
        rows.append(dict(scaled_radial_order=i, ordinary_Z_order=k, model=record(v),
                         saved_absolute_error_bound=record(error), actual_enclosure=record(value),
                         error_center_selected_as_value=False, factorial_normalization=False))
    phi, pr = actual[0, 0], actual[1, 0]
    require(phi.a > 0 and pr.b < 0, "Actual core sign not established")
    A = -2 * rho * pr / phi
    AZ = -2 * rho * (actual[1, 1] / phi - pr * actual[0, 1] / phi**2)
    AZZ = -2 * rho * (actual[1, 2] / phi - 2 * actual[1, 1] * actual[0, 1] / phi**2
                      - pr * actual[0, 2] / phi**2 + 2 * pr * actual[0, 1]**2 / phi**3)
    require(A.a > 0, "Strict positive actual core Dbar required")
    return dict(passed=True, mixed_derivatives=rows, Phi=record(phi), Phi_rho=record(pr),
                Dbar_ordinary_Z=[record(v) for v in (A, AZ, AZZ)],
                Dbar_formula="-2*rho*Phi_rho/Phi",
                ordinary_Z_derivatives=True, derivative_of_error_cap_used=False,
                actual_fixed_point_error_added_once=True)


def direction_controls():
    import sympy as s
    r, z = s.symbols("rho Z", positive=True)
    p = s.Function("Phi")(r, z); pr = s.diff(p, r)
    A = -2 * r * pr / p
    AZ = -2 * r * (s.diff(pr, z) / p - pr * s.diff(p, z) / p**2)
    AZZ = -2 * r * (s.diff(pr, z, 2) / p - 2 * s.diff(pr, z) * s.diff(p, z) / p**2
                   - pr * s.diff(p, z, 2) / p**2 + 2 * pr * s.diff(p, z)**2 / p**3)
    checks = [s.cancel(v - s.diff(A, z, k)) == 0 for k, v in enumerate((A, AZ, AZZ))]
    require(all(checks), "Fixed-rho ordinary direction derivatives differ")
    return dict(passed=True, fixed_rho_ordinary_Z0_Z1_Z2=checks, rho_Z_exact_zero=True,
                no_factorial_normalization=True, error_cap_differentiated=False)


def build(root, parent):
    from mpmath.ctx_iv import MPIntervalContext
    root = Path(root)
    require(sha(root / PARENT) == PARENT_SHA256 and json.loads((root / PARENT).read_bytes()) == parent,
            "Exact accepted micro mean parent required")
    require(sha(root / HELPER) == HELPER_SHA256, "Reviewed signed model helper changed")
    from .original_core_signed_model import core_model, numeric_controls
    paths = (UNIFORM, TRANSFER, TUBE, CORE_SOURCE, CORE_CHECK, BRIDGE_SOURCE,
             PREFIX + "compliant_core_uniform_bounds.py", PREFIX + "compliant_core_physical_field.py")
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256),
            HELPER: dict(path=HELPER, expected_sha256=HELPER_SHA256)}
    for path in paths:
        expected = parent["input_hashes"][Path(path).name]
        require(sha(root / path) == expected, "Changed actual core input: " + path)
        pins[path] = dict(path=path, expected_sha256=expected)
    uniform, transfer, tube, receipt = [json.loads((root / path).read_bytes())
                                        for path in (UNIFORM, TRANSFER, TUBE, CORE_CHECK)]
    angular = parent["actual_native_core_angular_direction"]["angular_direction_identity"]
    identity = receipt["checks"]["actual_oracle"]["exact_parameter_and_native_P0_identity"]
    require(receipt["all_passed"] and receipt["source_family"] == parent["source_family"] and identity["passed"] and
            identity["actual_source_scalar_and_unique_anchor_identity"]["j_matches_admitted_source"] and
            identity["actual_source_scalar_and_unique_anchor_identity"]["sigma_is_same_j_over_500"] and
            identity["exact_graph_logP_delta_logLambda_epsilon_recipes_verified"] and
            identity["interval_overlap_not_used_as_a_function_identity"], "Consumed same-source scalar identity required")
    require(transfer["contraction_proved"] and transfer["analytic_core_exists_for_compliant_source"] and
            transfer["universal_tube_delta_range"] == ["0", "1e-200"] and
            all(uniform[k] == parent["source_family"][k] for k in ("implicit_source_sha256", "datum_enclosure_sha256")) and
            transfer["analytic_core_family_sha256"] == uniform["analytic_core_family_sha256"] == tube["analytic_core_family_sha256"] and
            uniform["unique_analytic_field_bounded_not_sampled"] and
            uniform["leading_model"] == "sum_n (-chi*r/2)^n/[n!(n+1)!], 0<=chi<=1" and
            uniform["scaled_radial_domain"] == ["0", "4.1"] and uniform["real_axial_domain"] == ["-1", "1"] and
            angular["symbolic_proof"]["Dbar_identity_exact"] and
            angular["nonvanishing_domain"]["Phi_receipt"]["path"] == UNIFORM,
            "Same actual analytic core and direction identity required")
    guards = dict(core_formula=source_guard(root, CORE_SOURCE, "__init__", (
        "U0=a.add(a.scale(Z,4),a.fixed(self.parameters['j']))",
        "H=a.add(a.scale(a.mul(a.add(one,a.neg(de)),Z),'1/2'),a.mul(d,U0))",
        "sigma=g.quotient(self.parameters['j'],c(500),'positive fixed scalar 500')",
        "chi=a.div(a.mul(H,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))")),
        actual_delta_recipe=source_guard(root, CORE_SOURCE, "parameter_identity", (
        "expected_logP=g.add(g.unary('exp',c(40)),c(11))",
        "expected_delta=g.unary('exp',g.neg(g.add(g.mul(c(4),p['logP']),c(30))))",
        "assert p['logP']==expected_logP and p['delta']==expected_delta")),
        uniform_error_convention=source_guard(root, PREFIX + "compliant_core_uniform_bounds.py", "run", (
        "Phi_minus_model_absolute_bound=correction_norm*factor", "scaled_radial_order=i,axial_order=k")),
        actual_micro_core_argument=source_guard(root, BRIDGE_SOURCE, "__init__", (
        "self.hb=hb=inlet_owner.previous.lower.parameters['hb']",
        "rho=lambda q:g.mul(c(4),g.unary('exp',g.mul(hb,q)))",
        "atcore=lambda jet,q:a.at(jet,self.core.rho,rho(q))")))
    c = MPIntervalContext(); c.prec = 600
    j = read(c, transfer["required_j"]); sigma = j / 500
    logP = c.exp(40) + 11; delta = c.exp(-4 * logP - 30)
    require(j.a > 0 and sigma.a > 0 and delta.a > 0 and delta.b < c.mpf("1e-200").a,
            "Admitted actual parameters required; tube delta is a cap, not the selected delta")
    cases = {}
    for name, rho, z in (("selected_core_exit", c.mpf(4), c.mpf(371) / 1000),
                          ("core_tube", c.mpf([4, c.mpf(41) / 10]), c.mpf([c.mpf(37) / 100, c.mpf(93) / 250]))):
        model = core_model(c, rho, z, j, delta, sigma, terms=64)
        exact_domain = dict(rho=["4", "4"], Z=["371/1000", "371/1000"]) if name == "selected_core_exit" else \
            dict(rho=["4", "41/10"], Z=["37/100", "93/250"])
        cases[name] = dict(exact_requested_domain=exact_domain, rho=record(rho), Z=record(z), model_parameters=model["parameters"],
                           series_certificate=model["certificate"], **enclose_actual(c, model, uniform, rho))
    selected = read(c, parent["actual_signed_lower_m15"]["selected_Z_interval"])
    require(selected.a >= c.mpf(37).a / 100 and selected.b <= c.mpf(93).b / 250,
            "Saved actual query must lie inside the new axial tube")
    controls = numeric_controls(); algebra = direction_controls()
    require(controls["passed"], "Independent model arithmetic controls required")
    basis = parent["actual_signed_lower_m15"]
    require(basis["source_log_basis_names"][0] == "hb", "Actual native width log required")
    loghb = read(c, basis["source_log_bases"][0])
    require(loghb.b < c.ln(c.mpf(1) / 1000).a, "Native width must be below 1/1000")
    # For 0<=s<=2: hb*s<1/500 and exp(x)<=1/(1-x).
    # No attempt is made to materialize hb=exp(loghb).
    radius_upper = c.mpf(2000) / 499
    require(radius_upper.b < (c.mpf(41) / 10).a, "Both micro core arguments must stay in the signed tube")
    geometry = dict(passed=True, native_s_interval=["0", "2"],
                    source_log_hb=record(loghb), hb_strict_upper="1/1000",
                    core_argument="rho=4*exp(hb*s)", core_argument_upper=record(radius_upper),
                    proof="0<=hb*s<1/500; exp(x)<=1/(1-x); rho<2000/499<41/10",
                    hb_materialized=False, applies_to_core_arguments_not_second_comparison_direction=True)
    return dict(passed=True, source_family=parent["source_family"], source_pins=pins, source_AST=guards,
                precision_bits=c.prec, model_terms=64, parameters=dict(j=record(j), sigma=record(sigma),
                delta=record(delta), logPstar=record(logP), delta_from_actual_recipe=True,
                shared_tube_delta_is_only_universal_cap=True), cases=cases, numeric_controls=controls,
                accepted_core_direction_identity_digest=digest(angular), micro_core_argument_domain=geometry,
                direction_symbolic_controls=algebra,
                signed_actual_Phi_rho_negative=True, actual_core_Dbar_positive=True,
                core_seed_chi_distinct_from_micro_control_chi=True,
                ordinary_Dbar_Z0_Z1_Z2=True, actual_core_only=True,
                comparison_second_micro_Dbar_positive=False, actual_C5_mixed_root_provider_complete=False,
                mean_Z3_numeric=False, global_five_moment_closure=False,
                old_signed_terms_and_all_errors_unchanged=True, accepted_ancestors_replayed=False,
                full_goal_complete=False)
