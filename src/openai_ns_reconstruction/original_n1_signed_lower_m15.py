"""Signed actual lower m15 from one inlet and paired pressure primitives.

The raw velocity cutoff is zero on q in [1,2].  Only the genuine leading
N1p source remains; no independently boxed Qpsi endpoints are subtracted.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "n1_reserved_point_m15.json"
PARENT_SHA256 = "b2e5014312f87eee67868cb20bf10c49f5b4823ef929a09c6ad8262ed65f15c6"
MEAN = PREFIX + "n1_leading_mean_epsilon_factor.json.gz"
MEAN_SOURCE = PREFIX + "n1_leading_mean_epsilon_factor.py"
DOMAIN = PREFIX + "n1_inner_analytic_domain.json.gz"
MAJOR = BASE + "lei_ren_part1_paper_compliant_core_transfer.json"
LOWER = PREFIX + "n1_lower_collar_functions.json.gz"
LOWER_SOURCE = PREFIX + "n1_lower_collar_functions.py"
CUT_SOURCE = PREFIX + "n1_inner_cutoff_functions.py"
FLAT_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
INNER_SOURCE = "src/openai_ns_reconstruction/original_n1_signed_C_inner_m15.py"
HELPERS = tuple("src/openai_ns_reconstruction/" + n + ".py" for n in (
    "original_finite_N_signed_source", "original_n1_selected_word_numeric",
    "original_n1_selected_pressure_primitive", "original_n1_correlated_pressure_increment",
    "original_n1_velocity_m15_increment"))
ZERO = (0, 0, 0, 0, 0, 0)
H = (1, 0, 0, 0, 0, 0)
PRECISION_BITS = 1800


def extend(poly, second=False):
    return {(k[:4] + (0, k[4]) if second else k + (0,)): v for k, v in poly.items()}


def normalize(c, signed, errors, bases, dominant, phase, record, endpoint, require):
    value = signed[dominant] / phase
    radius = c.mpf(0); checks = []
    require(value.a > 0 or value.b < 0, "Actual lower/prefix main coefficient sign unresolved")
    for kind, p in (("signed", signed), ("error", errors)):
        for key, coefficient in p.items():
            if kind == "signed" and key == dominant:
                continue
            bound = endpoint(c, abs(coefficient))
            if bound == 0:
                continue
            log_upper = c.ln(bound) - c.ln(phase)
            for power, base in zip((x-y for x,y in zip(key, dominant)), bases):
                if power:
                    log_upper += power * base
            require(log_upper.b < -1024*c.ln(2), "A lower/prefix source term exceeds the admitted tail")
            radius += c.mpf(2)**-1024
            checks.append(dict(kind=kind, exponents=list(key), log_absolute_normalized_upper=record(log_upper),
                               dyadic_upper="2^-1024"))
    box = c.make_mpf(((value-radius)._mpi_[0], (value+radius)._mpi_[1]))
    require(box.b < 0, "Actual signed lower/prefix enclosure must be negative")
    return dict(normalized_integral_enclosure=record(box), dominant_source_exponents=list(dominant),
        dominant_coefficient=record(value), strict_sign="negative", other_terms_absolute_upper=record(radius),
        relative_interval_width_upper=record((box.b-box.a)/min(abs(box.a),abs(box.b))),
        dyadic_term_bounds=checks, caps_not_selected_as_values=True)


def exact_algebra():
    h, a, rho1, rho2, K, V, S, sig1, sig2 = sp.symbols("h a rho1 rho2 K V S sig1 sig2")
    weight = (rho2**2-rho1**2)/(2*h)
    primitive = lambda r, sig: -K*(r*r-16)+(1-h)*V*r*r*sig+S*(r*r-16)/4
    difference = -2*h*K*weight+(1-h)*V*(rho2*rho2*sig2-rho1*rho1*sig1)+h*weight*S/2
    e, r = sp.symbols("epsilon rho", positive=True)
    b = sp.Function("b")(r); fullS = sp.Function("S")(r)
    Np = -fullS/2+2*sp.diff(b,r)+r*sp.diff(b,r,2)
    core, axial, swirl, Pa, Pb = sp.symbols("core axial swirl Pressure_a Pressure_2a")
    checks = dict(shared_inlet_primitive_cancels=sp.expand(primitive(rho2,sig2)-primitive(rho1,sig1)-difference)==0,
        physical_lower_FTC=sp.expand(-e**2*r*Np-sp.diff(-e**2*r*r*sp.diff(b,r),r)-e**2*r*fullS/2)==0,
        exact_weighted_geometry=sp.simplify(weight.subs({rho1:4*sp.exp(h*a),rho2:4*sp.exp(2*h*a)})
            -8*sp.exp(2*h*a)*(sp.exp(2*h*a)-1)/h)==0,
        physical_log_radial_viscosity=sp.expand(((r*sp.diff(b,r)+r*r*sp.diff(b,r,2))+r*sp.diff(b,r))/r
            -(2*sp.diff(b,r)+r*sp.diff(b,r,2)))==0,
        lower_cutoff_argument_minimum=4*sp.Integer(1)-2 >= 1,
        shared_q1_pressure_endpoint=sp.expand((core+axial+swirl+Pa)+(Pb-Pa)
            -(core+axial+swirl+Pb))==0)
    if not all(checks.values()):
        raise ValueError("Paired actual lower pressure algebra failed")
    return dict(passed=True, checks={k:bool(v) for k,v in checks.items()})


def build(root, parent):
    root = Path(root); sha0 = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if sha0(root/PARENT) != PARENT_SHA256 or json.loads((root/PARENT).read_bytes()) != parent:
        raise ValueError("Accepted reserved-point parent required")
    for path in HELPERS:
        if sha0(root/path) != parent["implementation_hashes"][path]:
            raise ValueError("Changed executable lower helper: "+path)
    from .original_finite_N_signed_source import ast_guard, digest, load, require, sha
    from .original_n1_selected_word_numeric import interval, endpoint, record
    from .original_n1_selected_pressure_primitive import add, scale, shift, product, export, inertial_norms, norm_polynomials, axis_jet
    from .original_n1_correlated_pressure_increment import model_jets, flat_atom, correlated_error
    from .original_n1_velocity_m15_increment import saved_polynomial
    require(parent["checks"]["passed"] is True, "Passed parent required")
    pins = {PARENT:dict(path=PARENT, expected_sha256=PARENT_SHA256, passed=True)}
    for path in (MEAN, MEAN_SOURCE, DOMAIN, MAJOR, LOWER, LOWER_SOURCE, CUT_SOURCE, FLAT_SOURCE, INNER_SOURCE, *HELPERS):
        expected = (parent["input_hashes"][Path(path).name] if path.startswith(BASE)
                    else parent["implementation_hashes"][path])
        require(sha(root/path) == expected, "Changed lower source: "+path)
        pins[path] = dict(path=path, expected_sha256=expected, passed=True)
    mean, domain, major, lower = [load(root/path) for path in (MEAN, DOMAIN, MAJOR, LOWER)]
    require(mean["source_family"] == domain["source_family"] == lower["source_family"] == parent["source_family"],
            "One actual source family required")
    native = parent["actual_native_collar_transport"]
    require(native["passed"] is True and native["native_defining_leading_function_equals_canonical_collar"] is True and
        native["actual_own_mean_and_pressure_equal_canonical_actual_histories"] is True and
        native["domain"]["canonical_q_interval"] == [0,2] and native["domain"]["alpha_identically_one"] is True and
        native["domain"]["ordinary_Z_independent_endpoints"] is True, "Actual whole q0..2 first-window transport required")
    hist = parent["actual_n1_cutoff_history_identity"]
    require(hist["passed"] is True and hist["cutoff_identity"]["leading_N1p_remains_after_eta_zero"] is True and
        hist["history_identity"]["leading_lower_m15_C0_C1_nodes"] == [36592,36593] and
        lower["lower_collar_bindings"]["functions"]["m15"] == [36592,36593], "Same genuine leading lower source required")
    require(hist["N1p_units"]["passed"] is True and
        hist["N1p_units"]["checks"]["normalization_applied_once"] is True and
        hist["N1p_units"]["checks"]["physical_Pstar_squared_canceled_once"] is True and
        hist["cutoff_identity"]["core_inlet_not_reset"] is True and
        parent["core_pressure_source_identity"]["scaling"]["epsilon_core"] == "exp(-logLambda)" and
        all(parent["core_pressure_source_identity"]["scaling"]["source_scaling_flags"].values()),
        "Physical N1p and distinct core epsilon/pressure units must be restored exactly once")
    binding = lower["lower_collar_bindings"]
    require(binding["pressure_bindings"]["direct_differential_form_not_endpoint_subtraction"] is True and
        binding["geometry_certificate"]["sc_and_hb_Z_independent"] is True,
        "Lower must use the same direct physical pressure source and Z-independent geometry")
    fields = binding["fields"]
    require(fields["domain"] == dict(lower=35824,upper=35825) and
        fields["all_histories_start_at_true_core_rho4"] is True and
        native["domain"]["rho_upper_rational"] == "4/(1-1/2048) < 41/10" and
        native["domain"]["sc_positive_and_less_than_2_power_minus10"] is True,
        "Shared q1 seam and lower domain with unchanged core history required")
    c = MPIntervalContext(); c.prec = PRECISION_BITS
    A, B = mean["samples"]["actual_first_inlet"], mean["samples"]["selected_actual_inner_exit"]
    saved = parent["actual_correlated_leading_pressure_increment"]
    require(A["Z"] == B["Z"] == saved["selected_Z_interval"] == parent["actual_signed_C_inner_m15"]["selected_Z_interval"],
            "Actual same-inlet axial interval required")
    a = interval(c,saved["selected_phase"]); t = 2*a
    sc = c.make_mpf((tuple(native["domain"]["sc_exact_mpf_tuple"]),)*2)
    require(a.a == a.b and 4*a == sc and 0 < a < t < c.mpf(1)/4096,
            "The exact lower interval must be sc/4..sc/2 within the admitted first window")
    proof = B["uniform_complex_bound_proof"]
    require(proof["actual_mean_integral_measure"] == "drho=hb*rho*ds; phase<=1/2" and
        proof["norm_bounds_not_selected_comparison_field"] is True and proof["all_original_pressure_atoms_retained"] is True and
        interval(c,native["domain"]["hb_positive_less_than_half_source_margin"]).a > 0,
        "Uniform whole first-collar source bounds and positive width margin required")
    require(interval(c,proof["original_outer_eta"])._mpi_ == interval(c,domain["common_tubes"]["outer_Z_radius"])._mpi_ and
        4/(1-c.mpf(1)/2048) < interval(c,domain["actual_inner_geometry"]["real_rho_upper_bound"]).a and
        interval(c,proof["Qpsi_and_Q_bound_tube_radius"]) == interval(c,proof["original_outer_eta"])/4,
        "Uniform radial/Z bounds must cover rho(2a) and the paired primitive")
    delta = interval(c,parent["actual_selected_leading_pressure_primitive"]["delta_enclosure"])
    norms = inertial_norms(c,B,domain); S1,S2 = norm_polynomials(c,norms,delta)
    require(norms["gap"] == interval(c,saved["error_proof"]["whole_error_Z1_gap"]), "Same complex Z Cauchy gap required")
    K,V = model_jets(c,A,delta); axis = axis_jet(c,A,delta)
    geometry = lambda s: (c.make_mpf(((4*s)._mpi_[0],(4*s/(1-s/2))._mpi_[1])),
                          c.make_mpf(((16*s)._mpi_[0],(16*s/(1-s))._mpi_[1])))
    cd1,cw1 = geometry(a); cd2,cw2 = geometry(t)
    delta_weight = c.make_mpf(((16*a)._mpi_[0],(16*a/(1-a)**2)._mpi_[1]))
    _,E1,gap,err1 = correlated_error(c,a,cd1,cw1,B,domain,major)
    _,E2,gap2,err2 = correlated_error(c,t,cd2,cw2,B,domain,major)
    require(gap == gap2 == norms["gap"], "Same whole-function derivative convention required")
    atom1,atom2 = flat_atom(c,a),flat_atom(c,t)
    require(atom1["ordinary_Z_derivative_is_zero"] is True and atom2["ordinary_Z_derivative_is_zero"] is True,
            "Both prescribed-shear flat atoms must be ordinary-Z independent")
    # Bind the accepted C_inner assembly to THIS same q1 leading pressure
    # primitive and these same two cut-velocity terms, before direct joining.
    assembly = {}
    for order in ("0","1"):
        q1 = saved["ordinary_Z"][order]["m15_leading_N1p_local"]
        first = parent["actual_complete_first_collar_m15_increment"]["ordinary_Z"][order]
        cp = first["complete_first_collar_m15"]
        signed_check = add(saved_polynomial(c,q1["signed_approximation"]),
            *(saved_polynomial(c,first[n]["signed_inlet_product_integral"])
              for n in ("cut_axial_velocity","cut_swirl_coupling")))
        error_check = add(saved_polynomial(c,q1["correlated_error_radius"]),
            *(saved_polynomial(c,first[n]["error_radius"])
              for n in ("cut_axial_velocity","cut_swirl_coupling")))
        ci = parent["actual_signed_C_inner_m15"]["ordinary_Z"][order]
        signed_inner = add(saved_polynomial(c,ci["actual_core_m15"]["signed_approximation"]),signed_check)
        error_inner = add(saved_polynomial(c,ci["actual_core_m15"]["error_radius"]),error_check)
        require(export(signed_check)==cp["signed_approximation"] and export(error_check,error=True)==cp["error_radius"] and
            export(signed_inner)==ci["actual_C_inner_m15"]["signed_approximation"] and
            export(error_inner,error=True)==ci["actual_C_inner_m15"]["error_radius"],
            "Accepted C_inner does not contain the identical q1 primitive plus cut velocity and core")
        assembly[order] = dict(passed=True,q1_pressure_digest=digest(q1),first_collar_digest=digest(first),
                              actual_core_and_inner_digest=digest(ci))
    bases = [interval(c,v) for v in saved["first_four_source_log_bases"]] + [
        interval(c,atom1["log_sigma_bounds"]),interval(c,atom2["log_sigma_bounds"])]
    scalar = lambda v: {ZERO:c.mpf(v)}
    one_minus_h = add(scalar(1),shift(scalar(-1),H))
    rho1sq = add(scalar(16),shift(scalar(2*cw1),H))
    rho2sq = add(scalar(16),shift(scalar(2*cw2),H))
    embed4 = lambda p: {k+(0,0):v for k,v in p.items()}
    rows = {}
    for k in (0,1):
        regular = shift(scale(extend(K[k]),-2*delta_weight),(1,0,0,-2,0,0))
        flat_start = shift(scale(product(product(one_minus_h,rho1sq),extend(V[k])),-1),(0,0,0,-2,1,0))
        flat_end = shift(product(product(one_minus_h,rho2sq),extend(V[k])),(0,0,0,-2,0,1))
        inertia = shift(scalar(delta_weight*axis[k]/2),(1,0,0,-2,0,0))
        factor = 1 if k == 0 else 1/gap
        pressure_error = scale(add(extend(E1),extend(E2,second=True)),factor)
        inertial_error = add(shift(scale(embed4(S1),delta_weight*factor/2),(1,0,0,-3,0,0)),
                             shift(scale(embed4(S2),delta_weight*factor/2),(1,0,0,-4,0,0)))
        signed = add(regular,flat_start,flat_end,inertia); error = add(pressure_error,inertial_error)
        # Join via the DIRECT 0..2a primitive. The identical q1 pressure
        # endpoint and its model error cancel as functions before boxing.
        # Core plus the two retained cut-velocity integrals are added once.
        p2reg = shift(scale(extend(K[k]),-2*cw2),(1,0,0,-2,0,0))
        p2axis = shift(scalar(cw2*axis[k]/2),(1,0,0,-2,0,0))
        p2signed = add(p2reg,flat_end,p2axis)
        p2error = add(scale(extend(E2,second=True),factor),
            shift(scale(embed4(S1),cw2*factor/2),(1,0,0,-3,0,0)),
            shift(scale(embed4(S2),cw2*factor/2),(1,0,0,-4,0,0)))
        core = parent["actual_signed_C_inner_m15"]["ordinary_Z"][str(k)]["actual_core_m15"]
        vel = parent["actual_complete_first_collar_m15_increment"]["ordinary_Z"][str(k)]
        velocity_signed = add(*(extend(saved_polynomial(c,vel[name]["signed_inlet_product_integral"]))
            for name in ("cut_axial_velocity","cut_swirl_coupling")))
        velocity_error = add(*(extend(saved_polynomial(c,vel[name]["error_radius"]))
            for name in ("cut_axial_velocity","cut_swirl_coupling")))
        joined = add(extend(saved_polynomial(c,core["signed_approximation"])),velocity_signed,p2signed)
        joined_error = add(extend(saved_polynomial(c,core["error_radius"])),velocity_error,p2error)
        norm_lower = normalize(c,signed,error,bases,(1,0,1,-2,0,0),a,record,endpoint,require)
        norm_joined = normalize(c,joined,joined_error,bases,(0,0,1,-2,0,0),c.mpf(1),record,endpoint,require)
        # Restrict the joined native-prefix enclosure to the existing narrow
        # pulse point tube, then consume the original target once.
        point = parent["actual_reserved_point_m15"]
        require(point["passed"] is True and point["local_inner_unit"] == "epsilon^2*Pstar^2" and
            point["shared_unit"]["physical"] == "reference*A*Rv^2" and point["shared_unit"]["U_Z_exact_zero"] is True,
            "Reserved comparison requires the identical fixed local and physical pressure units")
        zparent, zp = interval(c,A["Z"]),interval(c,point["selected_Z_interval"])
        require(zparent.a <= zp.a <= zp.b <= zparent.b, "Common pulse point must lie in the lower source tube")
        log_abs = interval(c,point["local_to_common_unit_log_ratio"]) + c.ln(abs(interval(c,norm_joined["normalized_integral_enclosure"])))
        require(log_abs.b < -4096*c.ln(2), "Joined native prefix exceeds its common-unit cover")
        prefix_box = c.mpf([-(c.mpf(2)**-4096),0])
        pr = point["ordinary_Z"][str(k)]
        pre = interval(c,pr["preRp_common_unit_enclosure"])
        target = interval(c,pr["independent_negative_complete_pulse_target"])
        diff = -(prefix_box+pre)
        require(target.a>0 and (target+diff).a>0, "Independent target/response sign unresolved")
        rows[str(k)] = dict(regular_viscous=export(regular),flat_start=export(flat_start),flat_end=export(flat_end),
            axis_inertial=export(inertia),paired_pressure_error_radius=export(pressure_error,error=True),
            lower_inertial_error_radius=export(inertial_error,error=True),
            actual_lower_m15=dict(signed_approximation=export(signed),error_radius=export(error,error=True),
                normalized=norm_lower),
            direct_zero_to_2a_leading_pressure=dict(signed_approximation=export(p2signed),error_radius=export(p2error,error=True)),
            retained_cut_velocity_signed=export(velocity_signed),retained_cut_velocity_error=export(velocity_error,error=True),
            actual_inner_plus_lower_m15=dict(signed_approximation=export(joined),error_radius=export(joined_error,error=True),
                normalized=norm_joined),
            point_reserved_comparison=dict(selected_Z_interval=point["selected_Z_interval"],
                joined_prefix_common_unit_log_absolute_enclosure=record(log_abs),
                joined_prefix_common_unit_enclosure=record(prefix_box),
                preRp_common_unit_enclosure=pr["preRp_common_unit_enclosure"],
                independent_negative_complete_pulse_target=pr["independent_negative_complete_pulse_target"],
                correlated_reserved_minus_target=record(diff),reserved_response=record(target+diff),
                relative_correlated_discrepancy=record(abs(diff)/target),
                common_unit="reference*A*Rv^2",prefix_cover_is_not_function_value=True,
                original_pulse_consumed_once=True,preRp_source_errors_unchanged=True))
    guards = dict(
        accepted_collar_assembly=ast_guard(root,HELPERS[-1],"build",(
            'P=pressure["ordinary_Z"][str(k)]["m15_leading_N1p_local"]',
            'ps=saved_polynomial(c,P["signed_approximation"])',
            'pe=saved_polynomial(c,P["correlated_error_radius"])',
            'signed=add(ps,Us[k],Fs[k]);errors=add(pe,Ue[k],Fe[k])')),
        accepted_core_and_collar=ast_guard(root,INNER_SOURCE,"build",(
            'collar=parent["actual_complete_first_collar_m15_increment"]["ordinary_Z"][str(k)]["complete_first_collar_m15"]',
            'full_signed=add(signed,saved_polynomial(c,collar["signed_approximation"]))',
            'full_error=add(errors,saved_polynomial(c,collar["error_radius"]))')),
        whole_mean_bounds=ast_guard(root,MEAN_SOURCE,"complex_mean_bounds",(
            "sig=domain.sigma_jets(c,c.mpf([0,'.5']))", "actual_mean_integral_measure='drho=hb*rho*ds; phase<=1/2'",
            "norm_bounds_not_selected_comparison_field=True,all_original_pressure_atoms_retained=True")),
        cutoff=ast_guard(root,CUT_SOURCE,"original_cutoff_rows",(
            "arg=g.sub(g.mul(c(4),q),c(2))", "value=g.sub(g.one,sigma)", "zero_fraction='3/4'")),
        flat_one=ast_guard(root,FLAT_SOURCE,"sigma_jets",(
            "if hi>=1:rows.append(IntervalTaylor.constant(c,1,4))", "if hi>0 and lo<=mp.mpf('.5'):",
            "if lo<1 and hi>=mp.mpf('.5'):")),
        preserved_pressure_history=ast_guard(root,CUT_SOURCE,"own_cutoff_fields",(
            "leading_N1p_not_cut=True,nonzero_core_pressure_inlet_retained=True", "actual_cutoff_own_mean=True")),
        lower_pressure=ast_guard(root,LOWER_SOURCE,"direct_pressure_functions",(
            "('m15',fields['weighted_physical_measure'],-1)","original_m15_minus_applied_once=True",
            "signed_integrand_not_a_cap=True")))
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,parent_mapping_sha256=digest(parent),
        source_family=parent["source_family"],source_pins=pins,source_AST=guards,exact_algebra=exact_algebra(),
        precision_bits=PRECISION_BITS,ordinary_Z_orders=[0,1],ordinary_Z=rows,selected_Z_interval=A["Z"],
        accepted_q1_pressure_core_and_velocity_assembly=assembly,
        accepted_C_inner_contains_identical_q1_pressure_primitive=True,
        lower_direct_differential_source_bound=True,lower_sc_and_hb_Z_independent=True,
        both_prescribed_flat_atoms_ordinary_Z_independent=True,reserved_comparison_fixed_units_checked=True,
        lower_phase_interval=[record(a),record(t)],lower_source_handles=[36592,36593],inner_source_handles=[53161,53162],
        source_log_basis_names=["hb","F0_anchor","Pstar_squared","Lambda","sigma_a","sigma_2a"],
        source_log_bases=[record(v) for v in bases],flat_start_atom=atom1,flat_end_atom=atom2,
        weighted_radial_difference_over_hb=record(delta_weight),complex_error_proofs=[err1,err2],
        lower_normalization="hb*epsilon^2*Pstar^2*(sc/4)",joined_normalization="epsilon^2*Pstar^2",
        lower_identity="-epsilon^2*integral_rho(a)^rho(2a) rho*N1p_leading drho",
        pressure_FTC="-epsilon^2*[rho^2*b_rho]_a^2a+epsilon^2/2*integral rho*S",
        exact_weighted_geometry="Cw_delta=8*exp(2*hb*a)*(exp(2*hb*a)-1)/hb",
        new_join_uses_direct_zero_to_2a_primitive=True,shared_q1_pressure_endpoint_canceled_as_function=True,
        independent_q1_pressure_model_error_box_not_added_to_new_join=True,
        joined_identity="core+cut_axial[0,a]+cut_swirl[0,a]+leading_pressure[0,2a] = C_inner+lower",
        raw_velocity_cutoff_identically_zero_on_lower=True,current_native_first_window_covers_q0_to_q2=True,
        same_core_inlet_source_used_for_both_lower_endpoints=True,independent_endpoint_Psi_uncertainty_subtraction_not_used=True,
        shared_Qaxis_Psi_and_mean_constants_canceled_before_boxing=True,
        actual_signed_lower_value_enclosed=True,actual_inner_plus_lower_value_enclosed=True,
        own_mean_and_pressure_histories_retained=True,distinct_nonzero_flat_endpoint_atoms_retained=True,
        old_six_signed_terms_and_all_error_caps_unchanged=True,old_finite_prefix_equals_actual_raw_function=False,
        old_six_numerical_source_join_closed=False,forward_Rp_discrepancy_enclosed=False,
        all_Z_functional_coverage=False,global_five_moment_closure=False,full_goal_complete=False,
        accepted_ancestors_replayed=False,source_executed=False,P1_is_distinct_from_m15=True,one_point_only=True)
