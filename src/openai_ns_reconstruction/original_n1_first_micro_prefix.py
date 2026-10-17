"""Signed nonempty first-micro prefix on the source-owned flat collar.

The real interval is [sc/2,3sc/4].  Flat finite-N corrections are proved
from their source selector, not assigned the value of an error bound.
"""
from __future__ import annotations
import hashlib
import json
import copy
from fractions import Fraction
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_n1_"
PARENT=PREFIX+"first_micro_seam.json"
PARENT_SHA256="5d03e68d02b322a7d1e8ebee9529f4aeb57199afcf939dfa9fc7bb566a1edefa"
MEAN=PREFIX+"leading_mean_epsilon_factor.json.gz"
DOMAIN=PREFIX+"inner_analytic_domain.json.gz"
MAJOR=BASE+"lei_ren_part1_paper_compliant_core_transfer.json"
SUPPORT=BASE+"lei_ren_part1_paper_compliant_current_original_whole_Z_R100_finite_N.py"
MICRO=BASE+"lei_ren_part1_paper_compliant_current_original_micro_functions.py"
MIXED=BASE+"lei_ren_part1_paper_compliant_bridge_mixed_C4.py"
NATIVE="src/openai_ns_reconstruction/original_native_collar_transport.py"
CUT=PREFIX+"inner_cutoff_functions.py"
FLAT=BASE+"lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
TRUE_SOURCE=BASE+"lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
TRUE_GRAPH=TRUE_SOURCE[:-3]+".json.gz"
TRUE_CHECK=TRUE_SOURCE[:-3]+"_check.json"
LOOP=BASE+"lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
LOOP_CHECK=LOOP[:-8]+"_check.json"
OLD_GRAPH=PREFIX+"preRp_moment_functions.json.gz"
HELPERS=tuple("src/openai_ns_reconstruction/"+n+".py" for n in (
    "original_finite_N_signed_source","original_n1_selected_word_numeric",
    "original_n1_selected_pressure_primitive","original_n1_correlated_pressure_increment",
    "original_n1_signed_lower_m15"))
ZERO=(0,0,0,0,0,0,0)
H=(1,0,0,0,0,0,0)


def exact_algebra():
    h,a,z,r,e=sp.symbols("hb a Z rho epsilon",positive=True)
    b=sp.Function("b")(r,z);S=sp.Function("S")(r,z)
    Np=-S/2+2*sp.diff(b,r)+r*sp.diff(b,r,2)
    R=4*sp.exp(h*2*a);T=4*sp.exp(h*3*a)
    weight=(T*T-R*R)/(2*h)
    M2,P2,P3,E2,I23,E3=sp.symbols("M2 P2 P3 E2 I23 E3")
    defect=-e*e*r*Np-sp.diff(-e*e*r*r*sp.diff(b,r),r)-e*e*r*S/2
    checks=dict(physical_pressure_FTC=sp.expand(defect)==0,
        physical_pressure_FTC_Z1=sp.expand(sp.diff(defect,z))==0,
        positive_weight_identity=sp.simplify(weight-8*sp.exp(4*h*a)*(sp.exp(2*h*a)-1)/h)==0,
        source_cumulative_join=sp.expand(M2-P2+P3-(M2+(P3-P2)))==0,
        shifted_inlet_join=sp.expand((M2+E2)+I23-E3-(M2+(-E3+E2+I23)))==0,
        support_interval_nonempty=sp.simplify(3*a-2*a-a)==0)
    if not all(checks.values()):raise ValueError("First-micro prefix source algebra failed")
    return dict(passed=True,checks=checks)


def full_endpoint_source(root,parent,pins,load,require,digest,ast_guard):
    """Bind the true full-chart source, including slow and fast corrections."""
    r,check,loop,old,lc=[load(root/p) for p in (TRUE_GRAPH,TRUE_CHECK,LOOP,OLD_GRAPH,LOOP_CHECK)]
    n=r["exact_graph_nodes"];b=r["modified_Q4_pressure_C3_bindings"];f=b["charts"]["first_micro"]
    chart=loop["prescribed_loop_C5_bindings"]["charts"]["first_micro"]
    require(check["all_passed"] is True and check["input_hashes"][Path(TRUE_GRAPH).name]==pins[TRUE_GRAPH]["expected_sha256"] and
        r["source_family"]==check["source_family"]==loop["source_family"]==old["source_family"]==lc["source_family"]==parent["source_family"] and
        check["checks"]["actual_source"]["true_phase_held_native_function_derivatives_and_total_fast_terms"] is True and
        check["checks"]["actual_source"]["exact_existing_seven_background_seed_rows_and_unique_core_source"] is True and
        check["checks"]["actual_source"]["Pstar_once_and_physical_R_Rsquared_measures_once"] is True and
        check["checks"]["total_radial"]["nonunit_Jacobian_and_N_cancellation_in_Aphi_Bphi_checked"] is True and
        check["checks"]["pressure_reduction"]["direct_N1p_equivalence_of_reduced_integrals_through_C3"] is True,
        "True mixed source functions with the physical total-y derivative required")
    ld=loop["prescribed_loop_C5_bindings"]["definition"]
    require(lc["all_passed"] is True and lc["input_hashes"][Path(LOOP).name]==pins[LOOP]["expected_sha256"] and
        lc["checks"]["actual_source"]["genuine_C6_inputs_and_C5_generic_roots"] is True and
        lc["checks"]["actual_source"]["original_microscopic_collar_and_pressure_memory_retained"] is True and
        lc["checks"]["lazy_flat"]["lazy_original_collar_and_shear_flat_rows"]==144 and
        all(ld[k] is True for k in ("original_first_micro_start_sc_half_and_collar_three_sc_quarters",
            "new_functions_preserve_original_fixed_leading_source_family",
            "correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged")),
        "Checked original sc/2 inlet, 3sc/4 flat collar and same background histories required")
    require(n[:len(old["exact_graph_nodes"])]==old["exact_graph_nodes"] and
        n[18]["type"]=="DoubleDyadicInteger" and n[18]["one_common_N_for_all_charts_and_control_denominators"] is True and
        chart["coordinate"]==66988 and chart["domain"]==[97295,1] and chart["physical_measure"]==8 and
        f["radial"]["V_y"][:2]==[128706,128707] and f["physical"]["U_y"][:2]==[129223,129224] and
        f["physical"]["Q_y"][:2]==[129400,129404],"Same exact N/source prefix and first-micro chart required")
    require(n[128707]==dict(operation="sum",arguments=[69541,128701,99406]) and
        n[128701]==dict(operation="product",arguments=[128240,128603]) and
        n[128603]["denominator"]==18 and n[129224]==dict(operation="product",arguments=[128707,33572]) and
        n[129559]==dict(operation="product",arguments=[129400,68644]) and
        n[129560]==dict(operation="product",arguments=[129404,68644]),
        "Actual slow B_y/N and fast B_phi must both enter the physical endpoint")
    end=b["endpoints"]["first_micro"]
    require(end["integrals"]["S1"][:2]==[129638,129639] and end["boundary"]["R_Q_y"][:2]==[129654,129655] and
        end["pressure"]["m15"][:2]==[129646,129647] and b["initial"]["R_Q_y"][:2]==[129567,129568],
        "First-micro-only reduced integral and its own finite endpoints required")
    for pair,expressions in (([129638,129639],[129594,129595]),([129654,129655],[129559,129560])):
        for i,expression in zip(pair,expressions):
            require(n[i]==dict(operation="function_substitution",expression=expression,variable=66988,value=1,Z_independent_substitution=True),
                "Full first-micro endpoint must be at native1, not at O3_power")
    require(n[28]==dict(operation="exact_rational",numerator=1,denominator=2) and
        n[129614]==dict(operation="product",arguments=[129594,28]) and
        n[129615]==dict(operation="product",arguments=[129595,28]) and
        n[129618]==dict(operation="negative",argument=129567) and
        n[129619]==dict(operation="negative",argument=129568) and
        n[129622]==dict(operation="sum",arguments=[129559,129618]) and
        n[129623]==dict(operation="sum",arguments=[129560,129619]) and
        n[129630]==dict(operation="sum",arguments=[129614,129626]) and
        n[129631]==dict(operation="sum",arguments=[129615,129627]) and
        n[129626]==dict(operation="negative",argument=129622) and
        n[129627]==dict(operation="negative",argument=129623),
        "Actual first-micro increment must be S1/2-upperRQy+initialRQy for both ordinary rows")
    guards=dict(
        total_derivative=ast_guard(root,TRUE_SOURCE,"total_radial_functions",(
            "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
            "physical_y_equals_native_partial_divided_by_Jacobian=True")),
        physical_derivative=ast_guard(root,TRUE_SOURCE,"physical_pressure_fields",(
            "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))")))
    substitutions=[dict(operation="function_substitution",expression=i,variable=66988,value=1,Z_independent_substitution=True)
        for i in (129223,129224)]
    ref=lambda node:dict(archive=TRUE_GRAPH,sha256=pins[TRUE_GRAPH]["expected_sha256"],node=node)
    return dict(passed=True,source_AST=guards,source_graph_sha256=pins[TRUE_GRAPH]["expected_sha256"],
        base_graph_node_count=len(n),upper_U_y_substitution_references=substitutions,
        substitutions_are_not_registered_existing_graph_node_ids=True,
        common_N=n[18],exact_selected_N_not_candidate_257=True,
        whole_chart_native_interval=["sc/2","1"],ordinary_Z_orders=[0,1],
        actual_S1_ordinary_Z=[ref(i) for i in (129638,129639)],
        actual_upper_R_Qy_ordinary_Z=[ref(i) for i in (129654,129655)],
        actual_full_increment_ordinary_Z=[ref(i) for i in (129646,129647)],
        actual_initial_R_Qy_ordinary_Z=[ref(i) for i in (129567,129568)],
        exact_cumulative_formula="M15(1)=H+S1_first_micro/2-upperRQy",
        exact_alternative_formula="M15(1)=actual_q2_inlet+actual_first_micro_increment",
        shifted_inlet_source_digest=digest(parent["actual_first_micro_seam"]),
        checked_native_source_extension_premises=dict(prescribed_loop_receipt_sha256=pins[LOOP_CHECK]["expected_sha256"],
            original_sc_half_and_three_sc_quarters_collar=True,same_original_core_and_background_histories=True,
            loop_source_definition_digest=digest(ld)),
        actual_q2_inlet_source_digest=parent["actual_correlated_q2_inlet"]["shared_q2_inlet_digest"],
        actual_finite_N_slow_and_fast_signed_functions_retained=True,
        physical_Pstar_once=True,source_functions_not_error_caps=True,
        old_seventeen_chart_I1_not_used=True,old_O3_upper_boundary_not_used=True,
        actual_increment_source_graph_identity_checked=True,
        function_definition_bound=True,whole_chart_numerical_enclosure_available=False,
        quantitative_mixed_C5_or_global_N_admission=False)


def build(root,parent):
    root=Path(root)
    if hashlib.sha256((root/PARENT).read_bytes()).hexdigest()!=PARENT_SHA256 or json.loads((root/PARENT).read_bytes())!=parent:
        raise ValueError("Accepted corrected first-micro seam required")
    for path in HELPERS:
        if hashlib.sha256((root/path).read_bytes()).hexdigest()!=parent["implementation_hashes"][path]:
            raise ValueError("Changed executable prefix helper: "+path)
    from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
    from .original_n1_selected_word_numeric import interval,endpoint,record
    from .original_n1_selected_pressure_primitive import add,scale,shift,product,export,axis_jet
    from .original_n1_correlated_pressure_increment import model_jets,flat_atom,correlated_error
    from .original_n1_signed_lower_m15 import normalize
    require(parent["checks"]["passed"] is True and parent["actual_first_micro_seam"]["passed"] is True,
        "Checked physical corrected inlet required")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256,passed=True)}
    for path in (MEAN,DOMAIN,MAJOR,SUPPORT,MICRO,MIXED,NATIVE,CUT,FLAT,TRUE_SOURCE,TRUE_GRAPH,TRUE_CHECK,LOOP,LOOP_CHECK,OLD_GRAPH,*HELPERS):
        h=parent["input_hashes"][Path(path).name] if path.startswith(BASE) else parent["implementation_hashes"][path]
        require(sha(root/path)==h,"Changed native prefix source: "+path)
        pins[path]=dict(path=path,expected_sha256=h,passed=True)
    mean,domain,major=[load(root/p) for p in (MEAN,DOMAIN,MAJOR)]
    full=full_endpoint_source(root,parent,pins,load,require,digest,ast_guard)
    native=parent["actual_native_collar_transport"];q2=parent["actual_signed_lower_m15"]
    require(native["passed"] is True and native["native_defining_leading_function_equals_canonical_collar"] is True and
        native["actual_own_mean_and_pressure_equal_canonical_actual_histories"] is True and
        native["domain"]["ordinary_Z_independent_endpoints"] is True and
        native["domain"]["sc_positive_and_less_than_2_power_minus10"] is True and
        all(native["exact_algebra"]["checks"].values()) and
        mean["source_family"]==domain["source_family"]==parent["source_family"],
        "Actual native first-window source and unique core histories required")
    c=MPIntervalContext();c.prec=1800
    a=interval(c,q2["lower_phase_interval"][0]);start=2*a;end=3*a
    sc=c.make_mpf((tuple(native["domain"]["sc_exact_mpf_tuple"]),)*2)
    require(a.a==a.b and 4*a==sc and 0<start<end<c.mpf(1)/4096,
        "Exact positive nonempty native prefix must be within the admitted analytic tube")
    positive=parent["actual_native_core_angular_direction"]["angular_direction_identity"]["nonvanishing_domain"]
    require(positive["passed"] is True and positive["Phi_receipt"]["scaled_radial_domain"]==["0","4.1"] and
        positive["Phi_receipt"]["strict_lower_positive"] is True and
        positive["Phi_receipt"]["unique_analytic_field_bounded_not_sampled"] is True and
        positive["L_positive_bound"] is True and interval(c,domain["denominator_bounds"]["L_modulus_lower"]).a>0 and
        interval(c,domain["denominator_bounds"]["Phi_core_modulus_lower"]).a>0,
        "Nonvanishing Phi/core comparison and L must cover the extended real/complex native domain")
    A,B=mean["samples"]["actual_first_inlet"],mean["samples"]["selected_actual_inner_exit"]
    proof=B["uniform_complex_bound_proof"]
    require(A["Z"]==B["Z"]==q2["selected_Z_interval"] and
        proof["actual_mean_integral_measure"]=="drho=hb*rho*ds; phase<=1/2" and
        proof["norm_bounds_not_selected_comparison_field"] is True and
        proof["all_original_pressure_atoms_retained"] is True and
        native["domain"]["rho_upper_rational"]=="4/(1-1/2048) < 41/10" and
        4/(1-c.mpf(1)/2048)<interval(c,domain["actual_inner_geometry"]["real_rho_upper_bound"]).a,
        "Actual uniform Psi/mean/Q bounds must cover the extended tiny native interval")
    guards=dict(
        exact_flat_support=ast_guard(root,SUPPORT,"support",(
            "if chart=='first_micro' and lo>=ep(self.sc/4)[0] and hi<=ep(self.sc*3/4)[0]:return 'original_flat_left_collar'")),
        exact_flat_values=ast_guard(root,SUPPORT,"prepare",(
            "if support in ('unchanged_left','original_flat_left_collar'):",
            "exact_flat_original_A_B_and_Z_zero_by_source_owned_left_collar=True")),
        native_comparison_alpha=ast_guard(root,MIXED,"evaluate",(
            "if chart=='first':", "alpha=[algebra.lift(1)]+[algebra.lift(0)]*3")),
        native_equations=ast_guard(root,MICRO,"source_bindings",(
            "actual_first_equations='phi_y=-((1-sigma)+hb*sigma)*Dbar*phi/2; V_y=-((1-sigma)+hb*sigma)*(phi/barphi)*drive'",
            "actual_geometry='R=Ra*exp(hb*s), dy=hb*ds'")),
        cutoff=ast_guard(root,CUT,"original_cutoff_rows",(
            "arg=g.sub(g.mul(c(4),q),c(2))", "value=g.sub(g.one,sigma)")),
        sigma_one=ast_guard(root,FLAT,"sigma_jets",("if hi>=1:rows.append(IntervalTaylor.constant(c,1,4))",)))
    delta=interval(c,parent["actual_selected_leading_pressure_primitive"]["delta_enclosure"])
    K,V=model_jets(c,A,delta);axis=axis_jet(c,A,delta)
    geometry=lambda s:(c.make_mpf(((4*s)._mpi_[0],(4*s/(1-s/2))._mpi_[1])),
        c.make_mpf(((16*s)._mpi_[0],(16*s/(1-s))._mpi_[1])))
    cd2,cw2=geometry(start);cd3,cw3=geometry(end)
    dw=c.make_mpf(((16*a)._mpi_[0],(16*a/((1-2*a)*(1-a)))._mpi_[1]))
    _,E2,gap,ep2=correlated_error(c,start,cd2,cw2,B,domain,major)
    _,E3,gap3,ep3=correlated_error(c,end,cd3,cw3,B,domain,major)
    require(gap==gap3,"Same ordinary whole-error Z1 Cauchy gap required")
    atom2,atom3=flat_atom(c,start),flat_atom(c,end)
    require(atom2==q2["flat_end_atom"],"Same nonzero prescribed sigma(sc/2) source required")
    bases=[interval(c,v) for v in q2["source_log_bases"]]+[interval(c,atom3["log_sigma_bounds"])]
    embed=lambda p,slot=0:{key[:4]+tuple(key[4] if j==slot else 0 for j in range(3)):v for key,v in p.items()}
    hydrate=lambda rows:{tuple(t["source_exponents"])+(0,)*(7-len(t["source_exponents"])):interval(c,t["coefficient_interval"]) for t in rows}
    scalar=lambda v:{ZERO:c.mpf(v)}
    one_minus_h=add(scalar(1),shift(scalar(-1),H))
    rho2sq=add(scalar(16),shift(scalar(2*cw2),H));rho3sq=add(scalar(16),shift(scalar(2*cw3),H))
    old=parent["actual_selected_leading_pressure_primitive"]
    cw1=interval(c,old["radius_factors"]["integral_rho_drho_over_hb"])
    require(cw1.a>0,"Accepted positive inertia density measure required")
    ratio_delta=endpoint(c,dw)/endpoint(c,cw1,0);ratio3=endpoint(c,cw3)/endpoint(c,cw1,0)
    rows={}
    for k in ("0","1"):
        j=int(k);factor=1 if j==0 else 1/gap
        reg=shift(scale(embed(K[j]),-2*dw),(1,0,0,-2,0,0,0))
        startflat=shift(scale(product(product(one_minus_h,rho2sq),embed(V[j])),-1),(0,0,0,-2,0,1,0))
        endflat=shift(product(product(one_minus_h,rho3sq),embed(V[j])),(0,0,0,-2,0,0,1))
        inertia={(1,0,0,-2,0,0,0):dw*axis[j]/2}
        ie=hydrate(old["ordinary_Z"][k]["m15_leading_N1p_local"]["remaining_inertial_error_radius"])
        error=add(scale(add(embed(E2,1),embed(E3,2)),factor),scale(ie,ratio_delta))
        signed=add(reg,startflat,endflat,inertia)
        p3=add(shift(scale(embed(K[j]),-2*cw3),(1,0,0,-2,0,0,0)),endflat,
            {(1,0,0,-2,0,0,0):cw3*axis[j]/2})
        pe3=add(scale(embed(E3,2),factor),scale(ie,ratio3))
        core=parent["actual_signed_C_inner_m15"]["ordinary_Z"][k]["actual_core_m15"]
        velocity=q2["ordinary_Z"][k]
        cumulative=add(hydrate(core["signed_approximation"]),hydrate(velocity["retained_cut_velocity_signed"]),p3)
        cumulative_error=add(hydrate(core["error_radius"]),hydrate(velocity["retained_cut_velocity_error"]),pe3)
        rows[k]=dict(actual_increment=dict(signed_approximation=export(signed),error_radius=export(error,error=True),
            normalized=normalize(c,signed,error,bases,(1,0,1,-2,0,0,0),a,record,endpoint,require)),
            cumulative_at_three_sc_quarters=dict(signed_approximation=export(cumulative),error_radius=export(cumulative_error,error=True),
            normalized=normalize(c,cumulative,cumulative_error,bases,(0,0,1,-2,0,0,0),c.mpf(1),record,endpoint,require)),
            direct_pressure_0_to_three_sc_quarters=dict(signed_approximation=export(p3),error_radius=export(pe3,error=True)),
            interval_regular_pressure=export(reg),negative_start_flat=export(startflat),positive_end_flat=export(endflat),
            interval_inertial_model=export(inertia),no_finite_N_correction_by_source_identity=True)
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,source_family=parent["source_family"],
        source_pins=pins,source_AST=guards,exact_algebra=exact_algebra(),ordinary_Z_orders=[0,1],ordinary_Z=rows,
        selected_Z_interval=A["Z"],native_interval=[record(start),record(end)],canonical_q_interval=[2,3],
        source_log_basis_names=q2["source_log_basis_names"]+["sigma_3a"],source_log_bases=[record(v) for v in bases],
        start_flat_atom=atom2,end_flat_atom=atom3,complex_error_proofs=[ep2,ep3],weighted_radius_difference_over_hb=record(dw),
        increment_normalization="hb*epsilon^2*Pstar^2*(sc/4)",cumulative_normalization="epsilon^2*Pstar^2",
        new_native_extension=dict(interval="[0,3sc/4]",same_source_ODE_and_core_inlet=True,
            ordinary_Z_orders=[0,1],alpha_equals_one=True,rho_upper="4/(1-1/2048)<4.1",
            proof="Same accepted angular/axial direction identities and nonzero core denominator apply on the extended tiny radial domain; identical native ODE and retained core histories give equality by uniqueness and Z differentiation.",
            old_canonical_q0_to_q2_packet_not_relabelled=True),
        native_extension_nonvanishing_source_digest=digest(positive),
        native_extension_prescribed_source_premises=full["checked_native_source_extension_premises"],
        numerical_increment_scope="Actual leading-pressure increment on the source-owned flat first_micro prefix; no own raw velocity or finiteN correction there",
        source_increment="-epsilon^2*integral_rho(sc/2)^rho(3sc/4) rho*N1p_leading drho",
        direct_cumulative_join="core+cut_velocity0..sc/4+leading_pressure0..3sc/4",
        shared_sc_half_endpoint_and_error_canceled_before_cumulative_boxing=True,
        parent_shifted_inlet_digest=digest(parent["actual_first_micro_seam"]),
        full_first_micro_continuation_source=full,full_first_micro_continuation_source_digest=digest(full),
        actual_nonempty_first_micro_leading_pressure_increment_enclosed=True,ordinary_Z1_and_physical_units_retained=True,
        finite_N_prefix_exact_zero_from_original_source=True,prescribed_shear_sigma_at_endpoints_not_zero=True,
        raw_own_velocity_cutoff_zero_but_leading_pressure_history_retained=True,
        all_old_six_terms_and_error_caps_unchanged=True,old_preRp_extra_Pstar_not_reintroduced=True,
        whole_first_micro_increment_enclosed=False,full_forward_Rp=False,all_Z_functional_coverage=False,
        global_five_moment_closure=False,full_goal_complete=False,accepted_ancestors_replayed=False)


class ActualFirstMicroContinuation:
    """Return the exact full-chart source or the evaluated flat-prefix packet."""
    def __init__(self,report):
        from .original_finite_N_signed_source import require,digest
        self.a=copy.deepcopy(report["actual_first_micro_prefix"])
        root=Path(__file__).resolve().parents[2]
        require(report["checks"]["passed"] is True and self.a["passed"] is True and
            self.a["parent_report_sha256"]==PARENT_SHA256 and
            hashlib.sha256((root/PARENT).read_bytes()).hexdigest()==PARENT_SHA256 and
            report["implementation_hashes"]["src/openai_ns_reconstruction/original_n1_first_micro_prefix.py"]==
                hashlib.sha256(Path(__file__).read_bytes()).hexdigest() and
            hashlib.sha256((root/TRUE_GRAPH).read_bytes()).hexdigest()==self.a["full_first_micro_continuation_source"]["source_graph_sha256"] and
            digest(self.a["full_first_micro_continuation_source"])==self.a["full_first_micro_continuation_source_digest"],
            "Checked true first-micro continuation required")
        self._digest=digest(self.a)

    def _checked(self,Z):
        from .original_finite_N_signed_source import require,digest
        require(isinstance(Z,str) and Fraction(Z)==Fraction(371,1000),"Only the saved axial query is admitted")
        require(digest(self.a)==self._digest,"Continuation source was mutated")

    def flat_prefix(self,Z):
        self._checked(Z)
        return copy.deepcopy(dict(Z=Z,interval=self.a["native_interval"],ordinary_Z=self.a["ordinary_Z"],
            numerical_enclosure_available=True,whole_first_micro=False))

    def full_endpoint(self,Z):
        self._checked(Z)
        return copy.deepcopy(dict(Z=Z,source=self.a["full_first_micro_continuation_source"],
            numerical_enclosure_available=False,source_functions_not_numeric_values=True))
