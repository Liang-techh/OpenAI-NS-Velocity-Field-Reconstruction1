"""Actual signed A-A0 primitive and a shared correlated q2 m15 inlet.

Six source functions share principal/remainder leaves.  Their exact source
weights cancel before any aggregate box is chosen.  Old six-row numerical
approximations and all their error caps are retained as separate evidence.
"""
from __future__ import annotations
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_"
PARENT=PREFIX+"n1_signed_lower_m15.json"
PARENT_SHA256="08a5c9cde3fbc33272de00f97deb88bd92832f7df3b55f51001c1632fc6e3e06"
MEAN=PREFIX+"n1_leading_mean_epsilon_factor.json.gz"
MEAN_CHECK=PREFIX+"n1_leading_mean_epsilon_factor_check.json"
DOMAIN=PREFIX+"n1_inner_analytic_domain.json.gz"
LOWER=PREFIX+"n1_lower_collar_functions.json.gz"
PRIMITIVE_SOURCE="src/openai_ns_reconstruction/original_collar_signed_primitive.py"
PARTITION_SOURCE="src/openai_ns_reconstruction/original_collar_pressure_partition.py"
FIRST_SOURCE="src/openai_ns_reconstruction/original_first_micro_m15.py"
PRE_SOURCE="src/openai_ns_reconstruction/original_preRp_signed_inner_m15.py"
HELPERS=tuple("src/openai_ns_reconstruction/"+n+".py" for n in (
    "original_finite_N_signed_source","original_n1_selected_word_numeric",
    "original_n1_selected_pressure_primitive","original_n1_correlated_pressure_increment",
    "original_n1_velocity_m15_increment","original_n1_signed_lower_m15"))
SIX_WEIGHTS=(dict(core_axial=1),dict(core_pressure=1),dict(cut_axial=1),
    dict(direct_pressure=1,cut_swirl=1,principal=-1,remainder=-1),dict(principal=1),dict(remainder=1))
SIX_NAMES=("core_axial_coupling","core_pressure_integral_and_finite_boundary","first_lower_axial_width",
    "first_lower_regular_pressure_width","first_lower_singular_pressure","first_lower_signed_A_minus_A0")
REDUCED_WEIGHTS=dict(core_axial=1,core_pressure=1,cut_axial=1,direct_pressure=1,cut_swirl=1)
ZERO=(0,0,0,0,0,0)


def reduce_weights(rows):
    result={}
    for row in rows:
        for name,weight in row.items():result[name]=result.get(name,0)+weight
    return {name:weight for name,weight in result.items() if weight}


def exact_algebra():
    h,r,z,L,R,p,pA=sp.symbols("h rho Z L ratio p pA")
    A=2*r*r*z*(1-h)*R*p/L; A0=32*z*pA/L
    bracket=(r*r-16)*R*p+16*(R-1)*p+16*(p-pA)
    leaves={name:sp.Symbol(name) for name in (*REDUCED_WEIGHTS,"principal","remainder")}
    expanded=sum(sum(w*leaves[n] for n,w in row.items()) for row in SIX_WEIGHTS)
    checks=dict(exact_signed_A_minus_A0_split=sp.cancel(A-A0+h*A0-(1-h)*2*z*bracket/L)==0,
        bracket_source_factorization=sp.expand(bracket-(r*r*R*p-16*pA))==0,
        exact_six_source_weights=sp.expand(expanded-sum(leaves[n] for n in REDUCED_WEIGHTS))==0,
        principal_cancels_before_boxing=sum(row.get("principal",0) for row in SIX_WEIGHTS)==0,
        remainder_and_its_error_cancel_before_boxing=sum(row.get("remainder",0) for row in SIX_WEIGHTS)==0,
        one_core_cut_velocity_pressure_source_each=reduce_weights(SIX_WEIGHTS)==REDUCED_WEIGHTS)
    if not all(checks.values()):raise ValueError("Correlated q2 source algebra failed")
    return dict(passed=True,checks=checks)


def build(root,parent):
    root=Path(root);sha0=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if sha0(root/PARENT)!=PARENT_SHA256 or json.loads((root/PARENT).read_bytes())!=parent:
        raise ValueError("Accepted actual lower/q2 parent required")
    for path in HELPERS:
        if sha0(root/path)!=parent["implementation_hashes"][path]:raise ValueError("Changed executable q2 helper: "+path)
    from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
    from .original_n1_selected_word_numeric import interval,endpoint,record
    from .original_n1_selected_pressure_primitive import add,scale,shift,export
    from .original_n1_correlated_pressure_increment import lift,model_jets
    from .original_n1_velocity_m15_increment import saved_polynomial
    from .original_n1_signed_lower_m15 import extend,normalize
    require(parent["checks"]["passed"] is True,"Passed q2 source parent required")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256,passed=True)}
    for path in (MEAN,MEAN_CHECK,DOMAIN,LOWER,PRIMITIVE_SOURCE,PARTITION_SOURCE,FIRST_SOURCE,PRE_SOURCE,*HELPERS):
        expected=(parent["input_hashes"][Path(path).name] if path.startswith(BASE) else parent["implementation_hashes"][path])
        require(sha(root/path)==expected,"Changed q2 inlet source: "+path)
        pins[path]=dict(path=path,expected_sha256=expected,passed=True)
    mean,mc,domain,graph=[load(root/p) for p in (MEAN,MEAN_CHECK,DOMAIN,LOWER)]
    require(mean["source_family"]==mc["source_family"]==domain["source_family"]==graph["source_family"]==parent["source_family"] and
        mc["all_passed"] is True and mc["checks"]["original_axis_and_identical_core_actual_first_inlet_factors"] is True and
        mc["input_hashes"][Path(MEAN).name]==pins[MEAN]["expected_sha256"],"Actual common core inlet and radial derivative source required")
    primitive=parent["canonical_collar_signed_primitive"];partition=parent["actual_collar_pressure_partition"]
    require(primitive["passed"] is True and all(primitive["radial_checks"].values()) and
        primitive["ordinary_C1"]["A0_formula"]=="32*Z*Psi_core_rho(4,Z)/L" and
        primitive["signed_primitive"]["physical_epsilon_squared_factor_count"]==1 and
        partition["passed"] is True and partition["full_pressure_partition_for_same_leading_fields"] is True and
        partition["physical_Pstar_squared_canceled_once"] is True and
        partition["signed_primitive_binding"]["packet_sha256"]==digest(primitive),
        "Exact physical signed radial-order-one A/A0 source partition required")
    native=parent["actual_native_collar_transport"];q2=parent["actual_signed_lower_m15"]
    require(native["passed"] is True and native["actual_quotient_equals_ratio"] is True and
        native["quotient_conventions"]["pointwise"]=="Phi_native(s,Z)/Phi_comparison(s,Z)=ratio(s,Z)" and
        native["quotient_conventions"]["extra_exponential_after_evolved_field"] is False and
        native["domain"]["canonical_q_interval"]==[0,2] and q2["passed"] is True and
        q2["new_join_uses_direct_zero_to_2a_primitive"] is True,"Actual native quotient and complete q2 joined source required")
    c=MPIntervalContext();c.prec=1800
    A=mean["samples"]["actual_first_inlet"];B=mean["samples"]["selected_actual_inner_exit"]
    require(A["Z"]==q2["selected_Z_interval"] and interval(c,A["geometry"]["phase"])==0,
            "Same actual core inlet Z and phase required")
    S=interval(c,q2["lower_phase_interval"][1]);sc=c.make_mpf((tuple(native["domain"]["sc_exact_mpf_tuple"]),)*2)
    require(S.a==S.b and 2*S==sc and S>0 and S<c.mpf(1)/4096,"Exact admitted S=sc/2 required")
    proof=B["uniform_complex_bound_proof"];ep=q2["complex_error_proofs"][1]
    eta=interval(c,proof["original_outer_eta"]);gap=interval(c,ep["whole_error_Z1_gap"])
    require(eta._mpi_==interval(c,domain["common_tubes"]["outer_Z_radius"])._mpi_ and gap==eta/16 and
        interval(c,proof["Psi_and_meanPsi_bound_tube_radius"])==eta/2 and
        interval(c,proof["Qpsi_and_Q_bound_tube_radius"])==eta/4 and
        proof["actual_mean_integral_measure"]=="drho=hb*rho*ds; phase<=1/2" and
        proof["norm_bounds_not_selected_comparison_field"] is True and
        interval(c,proof["actual_Phi_over_core_modulus_upper"]).a>=1 and
        interval(c,ep["ratio_exponent_cover"]).b<c.mpf(1)/2,
        "One uniform actual complex source tube and ratio factor-two cover required")
    Lmin=endpoint(c,interval(c,domain["denominator_bounds"]["L_modulus_lower"]),0)
    floor=interval(c,domain["denominator_bounds"]["Phi_core_modulus_lower"])
    require(Lmin>0 and floor.a>0,"Complex-modulus lower bounds required")
    B1=extend(lift(c,B["uniform_complex_factored_bounds"]["Psi"][1],upper=True))
    require(all(key[0]==0 and key[1]==0 for key in B1),"Core derivative cover cannot contain a width pole or varying amplitude base")
    B2=interval(c,ep["core_Psi_rho2_bound"]);Ac=interval(c,ep["core_log_derivative_upper"])
    rho=endpoint(c,interval(c,domain["actual_inner_geometry"]["real_rho_upper_bound"]))
    require(B2.a>=0 and Ac.a>0 and Ac._mpi_==interval(c,proof["actual_core_logarithmic_derivative_modulus_upper"])._mpi_,
            "Same actual core derivative variation bounds required")
    atom=q2["flat_end_atom"]
    require(atom["argument"]==q2["lower_phase_interval"][1] and atom["ordinary_Z_derivative_is_zero"] is True and
        atom["sigma_is_strictly_positive"] is True and atom["monotone_on_selected_real_interval_certified"] is True and
        atom["monotonicity_proof"]["strictly_positive"] is True and
        atom["monotonicity_proof"]["flat_left_endpoint"]=="sigma(0)=0 from pinned source branch",
        "Positive original d-sigma measure and nonzero retained sigma(S) required")
    Cd=c.make_mpf(((4*S)._mpi_[0],(4*S/(1-S/2))._mpi_[1]))
    Cw=c.make_mpf(((16*S)._mpi_[0],(16*S/(1-S))._mpi_[1]))
    scalar=lambda v:{ZERO:c.mpf(v)}
    variation=scale(add(scale(B1,4*Cw),scalar(16*Cd*B2),
        shift(scale(B1,32*S*rho*Ac),(0,0,0,0,0,1))),2*(1+eta)/Lmin)
    delta=interval(c,parent["actual_selected_leading_pressure_primitive"]["delta_enclosure"])
    _,V=model_jets(c,A,delta);bases=[interval(c,v) for v in q2["source_log_bases"]]
    rows={}
    for k in (0,1):
        principal=shift(scale(extend(V[k]),16),(0,0,0,-2,0,1))
        remainder=shift(scale(extend(V[k]),-16),(1,0,0,-2,0,1))
        remainder_error=shift(scale(variation,1 if k==0 else 1/gap),(1,0,0,-2,0,1))
        require(all(key[1]==0 for p in (principal,remainder,remainder_error) for key in p),
                "Ordinary F0 derivatives must remain in the saved coefficients")
        cert=normalize(c,remainder,remainder_error,bases,(1,0,1,-2,0,1),c.mpf(1),record,endpoint,require)
        core=parent["actual_signed_C_inner_m15"]["ordinary_Z"][str(k)]
        vel=parent["actual_complete_first_collar_m15_increment"]["ordinary_Z"][str(k)]
        lowerrow=q2["ordinary_Z"][str(k)]
        p2=lowerrow["direct_zero_to_2a_leading_pressure"]
        leaves=dict(
            core_axial=dict(signed_approximation=core["core_axial_axis_source_approximation"],error_radius=core["core_axial_error"]),
            core_pressure=dict(signed_approximation=export(add(extend(saved_polynomial(c,core["core_viscous_endpoint"])),
                extend(saved_polynomial(c,core["core_axis_inertia"])))),error_radius=export(add(
                extend(saved_polynomial(c,core["core_inertial_error"])),extend(saved_polynomial(c,core["core_swirl_coupling_unsigned_error"]))),error=True)),
            cut_axial=dict(signed_approximation=vel["cut_axial_velocity"]["signed_inlet_product_integral"],error_radius=vel["cut_axial_velocity"]["error_radius"]),
            cut_swirl=dict(signed_approximation=vel["cut_swirl_coupling"]["signed_inlet_product_integral"],error_radius=vel["cut_swirl_coupling"]["error_radius"]),
            direct_pressure=copy.deepcopy(p2),principal=dict(signed_approximation=export(principal),error_radius=[]),
            remainder=dict(signed_approximation=export(remainder),error_radius=export(remainder_error,error=True),normalized=cert))
        # Leaves with the historical four/five-factor export retain their own
        # basis size. They are never silently zipped into a six-factor basis.
        physical_sources=dict(
            core_axial="2epsilon*integral_0^4 Ulead*Uz1raw drho",
            core_pressure="-epsilon^2*integral_0^4 rho*(2Flead*F1raw+N1p) drho",
            cut_axial="2epsilon*integral_0^a hb*rho*Ulead*eta*Uz1raw ds",
            cut_swirl="-2epsilon^2*integral_0^a hb*rho^2*Flead*eta*F1raw ds",
            direct_pressure="-epsilon^2*integral_0^S hb*rho^2*N1p_leading ds",
            principal="epsilon^2*A0*sigma(S)",
            remainder="epsilon^2*integral_0^S(A-A0) d-sigma")
        for name,leaf in leaves.items():
            leaf["actual_function_binding"]=dict(name=name,physical_source=physical_sources[name],
                accepted_parent_sha256=PARENT_SHA256,ordinary_Z_order=k,
                actual_minus_model_error_is_same_function_in_all_references=True,
                radius_not_registered_as_signed_error_value=True)
            leaf["source_log_basis_count"]=len(leaf["signed_approximation"][0]["source_exponents"])
        rows[str(k)]=dict(source_leaves=leaves,actual_signed_A_minus_A0=leaves["remainder"],
            accepted_joined_q2_enclosure=lowerrow["actual_inner_plus_lower_m15"],
            aggregate_evaluated_after_source_weight_reduction=True,
            remainder_radius_subtracted_as_value=False)
    components=[dict(name=name,source_weights=weights,error_leaf_weights=weights,
        signed_function_not_legacy_numeric_error_slot=True) for name,weights in zip(SIX_NAMES,SIX_WEIGHTS)]
    phase=list(S._mpi_[0]);base_step=list(q2["lower_phase_interval"][0]["lower_exact_mpf_tuple"])
    packet=dict(packet_type="actual_correlated_q2_m15_inlet_v1",source_family=parent["source_family"],
        Z="371/1000",selected_Z_interval=q2["selected_Z_interval"],q="2",chart="selected_lower_collar",
        exact_sc_tuple=native["domain"]["sc_exact_mpf_tuple"],exact_phase_tuple=phase,q_scale_sc_quarter_tuple=base_step,source_coordinate="s=sc/2; rho=4exp(hb*sc/2)",
        ordinary_Z_orders=[0,1],components=components,reduced_source_weights=reduce_weights(SIX_WEIGHTS),
        source_leaves_digest=digest({k:v["source_leaves"] for k,v in rows.items()}),
        source_leaves_reference="actual_correlated_q2_inlet.ordinary_Z[k].source_leaves",
        component_four_exact_pressure_partition_digest=digest(partition),
        aggregate_source="accepted actual core+cut velocity0..a+leading pressure0..2a",
        aggregate_enclosures={k:v["accepted_joined_q2_enclosure"] for k,v in rows.items()},
        accepted_parent_report=PARENT,accepted_parent_sha256=PARENT_SHA256,
        normalization="epsilon^2*Pstar^2",ordinary_Z1_is_physical_derivative=True,
        old_six_numeric_rows_replaced_in_archives=False,legacy_six_entry_format=False,
        source_weights_reduce_before_boxing=True,do_not_sum_independent_component_error_radii=True)
    guards=dict(
        signed_source=ast_guard(root,PRIMITIVE_SOURCE,"build",(
            "pref, prod0, prod1 = two * rho_s**2 * (one - hb)",
            "integrand = (aa - bb) * spt", '"signed_A_minus_A0_remainder": True' )),
        first_micro=ast_guard(root,FIRST_SOURCE,"inlet",("self._inner_engine.evaluate(Z,'2')",)),
        preRp=ast_guard(root,PRE_SOURCE,"inner_lower",("self._engine.evaluate(Z,'2')",)))
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,parent_mapping_sha256=digest(parent),
        source_family=parent["source_family"],source_pins=pins,source_AST=guards,exact_algebra=exact_algebra(),
        ordinary_Z_orders=[0,1],ordinary_Z=rows,precision_bits=1800,selected_phase=record(S),
        selected_Z_interval=q2["selected_Z_interval"],source_log_basis_names=q2["source_log_basis_names"],
        source_log_bases=q2["source_log_bases"],sigma_S_atom=atom,
        signed_remainder_normalization="hb*epsilon^2*Pstar^2*sigma(sc/2)",
        signed_remainder_model="-hb*epsilon^2*A0*sigma(S)",
        variation_bound_per_hb=export(variation,error=True),
        variation_proof=dict(bracket="(rho^2-16)*ratio*pcore+16*(ratio-1)*pcore+16*(pcore-pA)",
            bound_per_hb="4*Cw*B1+32*S*rho_max*Ac*sigma(S)*B1+16*Cd*B2",
            ratio_modulus_upper=2,core_Psi_rho2_bound=record(B2),complex_L_modulus_lower=record(Lmin),
            complex_Z_gap=record(gap),positive_dsigma_integral_mass="sigma(S)-sigma(0)=sigma(S)",
            one_minus_hb_real_modulus_upper=1,whole_error_Z1_uses_Cauchy=True),
        actual_primitive_binding=dict(packet_digest=digest(primitive),
            signed_primitive_C0_C1=primitive["signed_primitive"]["signed_primitive_C0_C1"],
            principal_C0_C1=primitive["term4_A0_sigma"]["term4_A0_sigma_C0_C1"],
            upper_variable=primitive["signed_primitive"]["variable_scope"]["upper_s"],
            substitute_exact_sc_half=True,radial_Psi_order=1),
        shared_q2_inlet=packet,shared_q2_inlet_digest=digest(packet),
        signed_A_minus_A0_value_and_ordinary_Z1_enclosed=True,
        exact_six_source_weights_reduce_to_accepted_q2=True,shared_remainder_error_canceled_before_aggregate_boxing=True,
        original_nonzero_sigma_not_replaced_by_cap=True,old_six_and_all_error_caps_unchanged=True,
        old_finite_numeric_prefix_equals_exact_function=False,downstream_first_micro_increment_evaluated=False,
        forward_Rp_discrepancy_enclosed=False,all_Z_functional_coverage=False,global_five_moment_closure=False,
        full_goal_complete=False,one_point_only=True,accepted_ancestors_replayed=False,source_executed=False)


class ActualQ2Inlet:
    """Read one checked correlated inlet for both new downstream consumers."""
    def __init__(self,report):
        from .original_finite_N_signed_source import digest,require
        a=report["actual_correlated_q2_inlet"]
        require(a["passed"] is True and report["checks"]["passed"] is True and
            a["parent_report_sha256"]==PARENT_SHA256 and a["shared_q2_inlet_digest"]==digest(a["shared_q2_inlet"]),
            "Checked actual correlated q2 inlet required")
        packet=a["shared_q2_inlet"]
        root=Path(__file__).resolve().parents[2]
        require(hashlib.sha256((root/PARENT).read_bytes()).hexdigest()==PARENT_SHA256 and
            report["implementation_hashes"]["src/openai_ns_reconstruction/original_n1_correlated_q2_inlet.py"]==
            hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "Changed provider or accepted parent source")
        parent=json.loads((root/PARENT).read_bytes())
        require(packet["source_family"]==parent["source_family"] and
            packet["selected_Z_interval"]==parent["actual_signed_lower_m15"]["selected_Z_interval"] and
            packet["exact_sc_tuple"]==parent["actual_native_collar_transport"]["domain"]["sc_exact_mpf_tuple"] and
            packet["exact_phase_tuple"]==parent["actual_signed_lower_m15"]["lower_phase_interval"][1]["lower_exact_mpf_tuple"] and
            packet["q_scale_sc_quarter_tuple"]==parent["actual_signed_lower_m15"]["lower_phase_interval"][0]["lower_exact_mpf_tuple"] and
            packet["aggregate_enclosures"]=={k:parent["actual_signed_lower_m15"]["ordinary_Z"][k]["actual_inner_plus_lower_m15"]
                for k in ("0","1")}, "Actual q2 provider must retain the accepted aggregate source enclosures")
        require([v["name"] for v in packet["components"]]==list(SIX_NAMES) and
            reduce_weights([v["source_weights"] for v in packet["components"]])==REDUCED_WEIGHTS and
            packet["reduced_source_weights"]==REDUCED_WEIGHTS and packet["legacy_six_entry_format"] is False,
            "Exact shared-leaf reduction required; old six-row numerical format is distinct")
        self._leaves=copy.deepcopy({k:v["source_leaves"] for k,v in a["ordinary_Z"].items()})
        require(digest(self._leaves)==packet["source_leaves_digest"],"Different shared source/remainder/error leaf payload")
        self._packet=copy.deepcopy(packet);self._packet_digest=digest(self._packet)

    def _query(self,Z,consumer):
        if not isinstance(Z,str) or Fraction(Z)!=Fraction(371,1000):
            raise ValueError("Only the saved Z=371/1000 point is admitted")
        from .original_finite_N_signed_source import digest,require
        require(digest(self._packet)==self._packet_digest and digest(self._leaves)==self._packet["source_leaves_digest"],
                "Shared inlet or remainder/error payload was mutated")
        result=copy.deepcopy(self._packet);result["consumer"]=consumer;result["Z"]=Z
        result["source_leaves_by_ordinary_Z"]=copy.deepcopy(self._leaves)
        return result

    def first_micro_inlet(self,Z):return self._query(Z,"first_micro_actual_q2_inlet")

    def preRp_inner_lower(self,Z):return self._query(Z,"preRp_actual_inner_lower_inlet")
