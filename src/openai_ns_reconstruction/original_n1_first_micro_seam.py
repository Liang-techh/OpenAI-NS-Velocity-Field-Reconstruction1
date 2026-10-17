"""Contract the actual sc/2 viscous seam before first-micro boxing.

The shifted inlet is M15(sc/2)+Rminus*Q_y(sc/2), not a replacement
for the unresolved first-micro volume source or upper endpoint.
"""
from __future__ import annotations
import copy
import hashlib
import json
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_n1_"
PARENT=PREFIX+"correlated_q2_inlet.json"
PARENT_SHA256="420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13"
GRAPH=PREFIX+"preRp_moment_functions.json.gz"
CHECK=PREFIX+"preRp_moment_functions_check.json"
SOURCE=PREFIX+"preRp_moment_functions.py"
MEAN=PREFIX+"leading_mean_epsilon_factor.json.gz"
MEAN_SOURCE=PREFIX+"leading_mean_epsilon_factor.py"
ENDPOINT=PREFIX+"preRc_endpoint_C1.json.gz"
ENDPOINT_CHECK=PREFIX+"preRc_endpoint_C1_check.json"
SWITCH=BASE+"lei_ren_part1_paper_compliant_current_original_switch_finite_N.py"
GENERIC=BASE+"lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs.py"
RECOVERY=BASE+"lei_ren_part1_paper_compliant_current_original_long_reshape_finite_N.py"
MICRO=BASE+"lei_ren_part1_paper_compliant_current_original_micro_functions.py"
WHOLE=BASE+"lei_ren_part1_paper_compliant_current_original_whole_Z_micro_macro.py"
CORE_SOURCE="src/openai_ns_reconstruction/original_n1_signed_C_inner_m15.py"
LOWER_SOURCE="src/openai_ns_reconstruction/original_n1_signed_lower_m15.py"
HELPERS=tuple("src/openai_ns_reconstruction/"+n+".py" for n in (
    "original_finite_N_signed_source","original_n1_selected_word_numeric",
    "original_n1_selected_pressure_primitive","original_n1_velocity_m15_increment"))
ZERO=(0,0,0,0,0,0)


def exact_algebra():
    r,z,e,de=sp.symbols("rho Z epsilon delta",positive=True)
    p=sp.Function("Psi")(r,z);m=sp.Function("meanPsi")(r,z)
    u0=sp.Function("Uaxis")(z);L=1-de*z*z;d=1-z*z
    qa=((1+de)*z*u0-d*sp.diff(u0,z))/L
    b=(2*z*p-(1-de)*z*m-d*sp.diff(m,z))/L
    U=u0+e*p;Q=qa+e*b;Uy=e*r*sp.diff(p,r)
    qp=((1+de)*z*U+2*z*Uy-d*sp.diff(U,z))/L
    rules={sp.diff(m,r):(p-m)/r,sp.diff(m,r,z):(sp.diff(p,z)-sp.diff(m,z))/r}
    defect=(qp-Q-e*r*sp.diff(b,r)).subs(rules)
    Np=-sp.Function("S")(r,z)/2+2*sp.diff(b,r)+r*sp.diff(b,r,2)
    ftc=-e*e*r*Np-sp.diff(-e*e*r*r*sp.diff(b,r),r)-e*e*r*sp.Function("S")(r,z)/2
    C,A,W,I0,IS,B0,BS,J=sp.symbols("core_axial cut_axial cut_swirl core_inertia collar_inertia boundary0 boundaryS first_micro_S1")
    inlet=C+I0-B0+A+W+IS-BS+B0
    Pstar=sp.symbols("Pstar",positive=True)
    archived_qp=((1+de)*z*U+2*z*Pstar*Uy-d*sp.diff(U,z))/L
    checks=dict(actual_q_plus_qy_chain_rule=sp.cancel(defect)==0,
        actual_q_plus_qy_ordinary_Z1=sp.cancel(sp.diff(defect,z))==0,
        physical_leading_pressure_FTC=sp.cancel(ftc)==0,
        physical_leading_pressure_FTC_ordinary_Z1=sp.cancel(sp.diff(ftc,z))==0,
        shared_core_and_sc_half_endpoints=sp.expand(inlet+BS-(C+I0+A+W+IS))==0,
        archived_extra_Pstar_endpoint_defect=sp.cancel(e*r*(archived_qp-qp)-e*r*2*z*(Pstar-1)*Uy/L)==0)
    if not all(checks.values()):raise ValueError("Actual first-micro seam source algebra failed")
    return dict(passed=True,checks=checks)


def corrected_endpoint_graph(nodes):
    """Append a source graph patch; never rewrite the accepted old graph."""
    replacements={30365:30363,30366:30364};cache=dict(replacements);added=[]
    def visit(i):
        if i in cache:return cache[i]
        if i<30365:return i
        node=copy.deepcopy(nodes[i]);changed=False
        for key in ("arguments",):
            if key in node:
                old=node[key];new=[visit(v) for v in old];node[key]=new;changed|=old!=new
        for key in ("argument","numerator","denominator","expression"):
            if key in node:
                old=node[key];node[key]=visit(old);changed|=node[key]!=old
        if changed:
            cache[i]=len(nodes)+len(added);added.append(node)
        else:cache[i]=i
        return cache[i]
    pairs=dict(lower_R_Qy=[30442,30443],active_m15=[30707,30708],preRp_m15=[31221,31220])
    corrected={name:[visit(i) for i in pair] for name,pair in pairs.items()}
    return dict(base_graph_node_count=len(nodes),appended_nodes=added,original_pairs=pairs,corrected_pairs=corrected,
        changed_original_to_appended={str(k):v for k,v in cache.items() if k not in replacements and k!=v},
        exact_replacements={str(k):v for k,v in replacements.items()},old_graph_unchanged=True,
        corrected_Uy="actual_background_source.V_y; already physical, no additional Pstar",
        archived_minus_corrected_R_Qy="R*2Z*(Pstar-1)*V_y_raw/L",
        archived_minus_corrected_R_Qy_Z="R*2(Pstar-1)*[(V_y_raw+Z*V_y_raw_Z)/L+2delta*Z^2*V_y_raw/L^2]",
        Pstar_and_radius_Z_independent=True)


def graph_scale_defect_checks(nodes,patch):
    """Check actual old/new graph arithmetic against the physical defect."""
    z,de,R,P,Y,Yz,U,Uz,Uzz,q,qz=sp.symbols("Z delta R Pstar Vy VyZ U UZ UZZ q qZ")
    known={0:sp.Integer(0),1:sp.Integer(1),20:z,48:sp.Integer(2),220:sp.Integer(-1),
        25203:P,25220:de,25249:z*z,25246:2*z,30357:U,30358:Uz,30359:q,30360:qz,
        30361:R,30363:Y,30364:Yz,30372:Uzz}
    full=nodes+patch["appended_nodes"]
    admitted=set(range(30365,30415))|{int(v) for v in patch["changed_original_to_appended"]}|set(range(len(nodes),len(full)))
    def value(i):
        if i in known:return known[i]
        if i not in admitted:return sp.Symbol("unchanged_source_"+str(i))
        n=full[i];op=n["operation"]
        if op=="product":result=sp.prod(value(j) for j in n["arguments"])
        elif op=="sum":result=sum(value(j) for j in n["arguments"])
        elif op=="negative":result=-value(n["argument"])
        elif op=="positive_quotient":result=value(n["numerator"])/value(n["denominator"])
        else:raise ValueError("Unexpected changed endpoint graph operation: "+op)
        known[i]=result;return result
    L=1-de*z*z
    expected=[R*2*z*(P-1)*Y/L,R*2*(P-1)*((Y+z*Yz)/L+2*de*z*z*Y/L**2)]
    checks={name+"_Z"+str(k):sp.cancel(value(old)-value(new)-expected[k])==0
        for name,pair in patch["original_pairs"].items()
        for k,(old,new) in enumerate(zip(pair,patch["corrected_pairs"][name]))}
    if not all(checks.values()):raise ValueError("Actual graph repair changed more than the physical scale defect")
    return dict(passed=True,checks=checks,unchanged_integrals_are_independent_source_symbols=True)


def build(root,parent):
    root=Path(root)
    if hashlib.sha256((root/PARENT).read_bytes()).hexdigest()!=PARENT_SHA256 or json.loads((root/PARENT).read_bytes())!=parent:
        raise ValueError("Accepted actual correlated q2 parent required")
    for path in HELPERS:
        if hashlib.sha256((root/path).read_bytes()).hexdigest()!=parent["implementation_hashes"][path]:
            raise ValueError("Changed executable seam helper: "+path)
    from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
    from .original_n1_selected_word_numeric import interval,endpoint,record
    from .original_n1_selected_pressure_primitive import add,scale,export,axis_jet
    require(parent["checks"]["passed"] is True and parent["actual_correlated_q2_inlet"]["passed"] is True,
            "Checked actual q2 source required")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256,passed=True)}
    for path in (GRAPH,CHECK,SOURCE,MEAN,MEAN_SOURCE,ENDPOINT,ENDPOINT_CHECK,SWITCH,GENERIC,RECOVERY,MICRO,WHOLE,CORE_SOURCE,LOWER_SOURCE,*HELPERS):
        h=parent["input_hashes"][Path(path).name] if path.startswith(BASE) else parent["implementation_hashes"][path]
        require(sha(root/path)==h,"Changed first-micro seam source: "+path)
        pins[path]=dict(path=path,expected_sha256=h,passed=True)
    graph,check,mean,end,ec=[load(root/p) for p in (GRAPH,CHECK,MEAN,ENDPOINT,ENDPOINT_CHECK)]
    require(check["all_passed"] is True and graph["source_family"]==check["source_family"]==mean["source_family"]==parent["source_family"] and
        check["input_hashes"][Path(GRAPH).name]==pins[GRAPH]["expected_sha256"] and
        check["checks"]["actual_graph_relocation"]["true_lower_Uy_and_UyZ_provider_bound"] is True and
        check["checks"]["physical_pressure"]["both_finite_viscous_endpoints_retained"] is True and
        check["checks"]["physical_pressure"]["both_true_ordinary_Z_rows_checked_symbolically"] is True,
        "Checked physical pressure graph and actual lower V_y provider required")
    require(ec["all_passed"] is True and ec["input_hashes"][Path(ENDPOINT).name]==pins[ENDPOINT]["expected_sha256"] and
        ec["checks"]["lower_actual_V_y_rows_and_U_ZZ_checked"] is True and
        ec["checks"]["exact_four_Z_cell_cover"] is True and
        ec["checks"]["physical_common_Pstar_once"] is True and end["source_family"]==parent["source_family"],
        "Source-owned finite-N zero inlet and retained physical leading jets required")
    zero=end["samples"]["('0', '.5')"]["lower"]["proof"]
    require(zero["source_family"]==parent["source_family"] and zero["endpoint"]=="first_micro/sc_half singleton" and
        all(zero[k] is True for k in ("source_owned_collar_q_C0_Z_exact_zero","finite_N_initial_history_zero",
            "finite_N_predecessor_history_retained","finite_N_local_driver_exact_zero","actual_V_y_row_0_used",
            "actual_V_y_row_1_retained","actual_V_y_rows_0_through_3_factorial_converted_once","common_Pstar_applied_once")) and
        zero["incoming_history"]["all_four_incoming_N_scaled_mean_rows_exact_zero"] is True and
        all(v["exact_zero"] is True for v in zero["incoming_history"]["incoming_rows"]) and
        zero["flat_driver"]["source_owned_initial_condition"]["q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar"] is True and
        zero["flat_driver"]["source_owned_initial_condition"]["actual_background_incoming_histories_remain_nonzero"] is True,
        "Only finite-N corrections vanish; the actual leading lower field and histories remain")
    nodes=graph["exact_graph_nodes"];binding=graph["defining_bindings"]["active"]
    lower=binding["finite_endpoints"][0]["R_q_y"]
    require(lower==[30442,30443] and graph["signed_pieces"]["active"]["m15"]==[30707,30708] and
        nodes[30442]==dict(operation="product",arguments=[30361,30413]) and
        nodes[30443]==dict(operation="product",arguments=[30361,30414]),"Actual finite Rminus Q_y source ports changed")
    for k,i in enumerate((30363,30364)):
        node=nodes[i]
        require(node["operation"]=="current_original_physical_source_quantity" and node["quantity"]=="V_y" and
            node["source_family"]==parent["source_family"] and node["native_chart"]=="first_micro" and
            node["coordinate"]==29 and node["ordinary_slow_Z_derivative_order"]==k and
            node["source_projection_binding"]["module_sha256"]==pins[SOURCE]["expected_sha256"] and
            node["source_projection_binding"]["function"]=="lower_velocity_y_source" and
            node["defining_quantity_not_a_range_value"] is True and
            node["finite_N_lower_initial_history_and_local_driver_exact_zero"] is True,
            "Actual V_y must remain a defining source, not a selected cap or zero")
    require(nodes[25203]["name"]=="Pstar" and nodes[25203]["Z_independent"] is True and
        nodes[30365]==dict(operation="product",arguments=[25203,30363]) and
        nodes[30366]==dict(operation="product",arguments=[25203,30364]),
        "Archived endpoint's extra Pstar multiplication must be identified exactly")
    patch=corrected_endpoint_graph(nodes)
    graph_checks=graph_scale_defect_checks(nodes,patch)
    native=parent["actual_native_collar_transport"];q2=parent["actual_signed_lower_m15"]
    require(native["passed"] is True and native["native_defining_leading_function_equals_canonical_collar"] is True and
        native["actual_own_mean_and_pressure_equal_canonical_actual_histories"] is True and
        native["domain"]["native_first_s_interval"]=="[0,sc/2]" and native["domain"]["ordinary_Z_independent_endpoints"] is True and
        q2["current_native_first_window_covers_q0_to_q2"] is True and q2["own_mean_and_pressure_histories_retained"] is True,
        "Same native leading field and actual histories through the sc/2 endpoint required")
    guards=dict(
        pressure=ast_guard(root,SOURCE,"pressure_functions",("m15=a.sub(half(S1), boundary('R_q_y'))",)),
        lower_Vy=ast_guard(root,SOURCE,"lower_velocity_y_source",(
            "endpoint.c0._zero_first_incoming_proof(fields, tuple(ends))",
            "endpoint.c0._flat_driver_proof(fields, tuple(ends), 'first_micro')",
            "values = packet['actual_background_source']['V_y']")),
        lower_qy=ast_guard(root,SOURCE,"active_functions",(
            "Uy = source.C1Function(*(g.mul(S, v) for v in rows))",
            "qy = algebra.sub(qp, q)")),
        actual_mean=ast_guard(root,MEAN_SOURCE,"factored_leading_rows",(
            "w=[op.truncate(u,n)+eps(P) if i==0 else eps(P)",
            "qpsi=[(2*z*P-z*M*(1-physical.delta)-d*op.dz(M))/L")),
        core_FTC=ast_guard(root,CORE_SOURCE,"build",(
            'visc=shift(scale(lift(c,A["leading_Qpsi_rho_rows"][1],k),-16)',
            'inertia=shift(scalar(c,4*Saxis[k])')))
    guards.update(
        raw_physical_source=ast_guard(root,SWITCH,"recover_source",(
            "phi,V=source['fields']['phi'],source['fields']['V']",
            "velocity=dict(theta=[E,Ey],axial=[V,source['V_y']])")),
        normalized_generic_source=ast_guard(root,GENERIC,"recover_inputs",(
            "invS=f.factor((0,-.5,0,0,0));S=op.Pstar",
            "V=scale(raw['velocity']['axial'][0],invS);Vy=scale(raw['velocity']['axial'][1],invS)")),
        recovery_keeps_normalization=ast_guard(root,RECOVERY,"compile_recovery",(
            "elif name=='C':", "else:nodes.append(node)",
            "assert not skipping and len(replacements)==2")),
        actual_micro_derivative=ast_guard(root,MICRO,"evaluate",(
            "V_y=Vy", "fields=dict(phi=at['phi'],V=at['V'])")),
        actual_micro_owner=ast_guard(root,WHOLE,"owner",(
            "original_micro = micro.MicroFunctions(flow, series, decoder,",
            "source.bridge._encode(source_proof['fresh_actual_core'])")))
    c=MPIntervalContext();c.prec=1800
    S=interval(c,q2["lower_phase_interval"][1]);a=interval(c,q2["lower_phase_interval"][0])
    require(S.a==S.b and S==2*a and 0<S<1,"Exact sc/2 endpoint independent of Z required")
    require(zero["flat_driver"]["source_owned_initial_condition"]["selected_positive_s_c"]["lower_exact_mpf_tuple"]==
        native["domain"]["sc_exact_mpf_tuple"] and
        interval(c,q2["selected_Z_interval"]).a>0 and interval(c,q2["selected_Z_interval"]).b<c.mpf('.5'),
        "The actual saved Z tube and sc must lie in the proved flat-inlet source cell")
    p=parent["actual_selected_leading_pressure_primitive"]
    A,B=mean["samples"]["actual_first_inlet"],mean["samples"]["selected_actual_inner_exit"]
    require(A["Z"]==B["Z"]==q2["selected_Z_interval"] and
        B["uniform_complex_bound_proof"]["actual_mean_integral_measure"]=="drho=hb*rho*ds; phase<=1/2" and
        B["uniform_complex_bound_proof"]["all_original_pressure_atoms_retained"] is True,
        "Uniform accepted inertial density must cover the full new 0..sc/2 interval")
    cw1=interval(c,p["radius_factors"]["integral_rho_drho_over_hb"])
    require(cw1.a>0 and interval(c,p["selected_phase"])==a,"Positive accepted sc/4 inertia geometry required")
    cw2=c.make_mpf(((16*S)._mpi_[0],(16*S/(1-S))._mpi_[1]))
    # The accepted q1 inertial radius is Cw1 times the same uniform
    # complex source-density bound.  A positive upper ratio extends it
    # to 0..S without recomputing the accepted source norm producer.
    ratio=endpoint(c,cw2)/endpoint(c,cw1,0)
    axis=axis_jet(c,A,interval(c,p["delta_enclosure"]))
    extend=lambda rows:{tuple(t["source_exponents"])+(0,)*(6-len(t["source_exponents"])):interval(c,t["coefficient_interval"]) for t in rows}
    bases=[interval(c,v) for v in q2["source_log_bases"]]
    rows={}
    for k in ("0","1"):
        core=parent["actual_signed_C_inner_m15"]["ordinary_Z"][k]
        vel=parent["actual_complete_first_collar_m15_increment"]["ordinary_Z"][k]
        collar_axis={(1,0,0,-2,0,0):cw2*axis[int(k)]/2}
        collar_error=scale(extend(p["ordinary_Z"][k]["m15_leading_N1p_local"]["remaining_inertial_error_radius"]),ratio)
        signed=add(extend(core["core_axial_axis_source_approximation"]),extend(core["core_axis_inertia"]),collar_axis,
            *(extend(vel[n]["signed_inlet_product_integral"]) for n in ("cut_axial_velocity","cut_swirl_coupling")))
        errors=add(*(extend(core[n]) for n in ("core_inertial_error","core_axial_error","core_swirl_coupling_unsigned_error")),
            collar_error,*(extend(vel[n]["error_radius"]) for n in ("cut_axial_velocity","cut_swirl_coupling")))
        for poly in (signed,errors):require(all(len(key)==6 for key in poly),"Complete six-factor physical source basis required")
        # Bound the actual shifted inlet in the already fixed common unit.
        # Only an absolute enclosure is asserted: raw-core errors may
        # dominate the signed model, so their upper bounds are not values.
        logs=[];dominant=(0,0,1,-2,0,0)
        for kind,poly in (("signed_model",signed),("error_radius",errors)):
            for key,value in poly.items():
                bound=endpoint(c,abs(value))
                if bound==0:continue
                log=c.ln(bound)+interval(c,parent["actual_reserved_point_m15"]["local_to_common_unit_log_ratio"])
                for power,base in zip((x-y for x,y in zip(key,dominant)),bases):
                    if power:log+=power*base
                logs.append((kind,key,log))
        log_upper=max(endpoint(c,v) for _,_,v in logs)+c.ln(len(logs))
        rows[k]=dict(signed_approximation=export(signed),error_radius=export(errors,error=True),
            collar_axis_inertial_model=export(collar_axis),collar_inertial_error_radius=export(collar_error,error=True),
            core_viscous_and_collar_derivative_endpoint_not_boxed=True,
            absolute_common_unit_log_upper=record(log_upper),common_unit="reference*A*Rv^2",
            logarithmic_triangle_terms=[dict(kind=kind,source_exponents=list(key),log_absolute_upper=record(value)) for kind,key,value in logs],
            strict_sign_claimed=False,caps_not_selected_as_function_values=True)
    inlet=parent["actual_correlated_q2_inlet"]
    source_rows=dict(core_pressure=dict(E0=-1,Icore=1,core_swirl=1),
        direct_collar_pressure=dict(E0=1,ES=-1,Icollar=1),first_micro_lower_RQy=dict(ES=1))
    reduced={}
    for weights in source_rows.values():
        for name,weight in weights.items():reduced[name]=reduced.get(name,0)+weight
    reduced={name:weight for name,weight in reduced.items() if weight}
    require(reduced==dict(Icore=1,core_swirl=1,Icollar=1),"Only identical physical endpoints may cancel")
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,source_family=parent["source_family"],
        source_pins=pins,source_AST=guards,exact_algebra=exact_algebra(),ordinary_Z_orders=[0,1],ordinary_Z=rows,
        selected_Z_interval=q2["selected_Z_interval"],exact_phase=q2["lower_phase_interval"][1],
        source_log_basis_names=q2["source_log_basis_names"],source_log_bases=q2["source_log_bases"],
        accepted_q2_inlet_digest=inlet["shared_q2_inlet_digest"],
        finite_N_zero_inlet_proof_digest=digest(zero),native_leading_identity_digest=digest(native),
        corrected_physical_endpoint_graph=patch,corrected_physical_endpoint_graph_digest=digest(patch),
        actual_graph_scale_defect_checks=graph_checks,
        archived_preRp_endpoint_has_extra_Pstar=True,archived_preRp_m15_is_not_corrected_physical_source=True,
        shared_pressure_endpoint_source_reduction=dict(source_weight_rows=source_rows,reduced_weights=reduced,
            E0="epsilon^2*16*Qpsi_core_rho(4,Z)",ES="epsilon^2*rho(S)^2*Qpsi_actual_rho(S,Z)",
            Icore="epsilon^2/2*integral_0^4 rho*S_core drho",
            Icollar="epsilon^2/2*integral_4^rho(S) rho*S_leading drho",
            core_swirl="-2*epsilon^2*integral_0^4 rho*Flead*F1raw drho",
            endpoint_error_is_actual_function_minus_same_shared_model=True,
            no_endpoint_error_radius_subtracted_as_a_signed_value=True,
            original_q2_source_packet_digest=digest(inlet["shared_q2_inlet"])),
        new_cross_map=dict(actual_native_field_identity="U=Uaxis+epsilon*Psi; Q=Qaxis+epsilon*Qpsi",
            radial_chain_rule="d/dy=rho*d/drho; Uy=epsilon*rho*Psi_rho",
            physical_radius="R=epsilon*rho",mean_FTC="rho*meanPsi_rho=Psi-meanPsi",
            ordinary_Z1_mean_FTC="rho*meanPsi_rhoZ=Psi_Z-meanPsi_Z",
            smooth_native_identity_on_closed_first_window=True,finite_N_correction_lower_jet_zero=True,
            derivative_and_own_mean_source_not_replaced_by_comparison=True),
        lower_R_Qy_source_nodes=lower,lower_actual_Vy_source_nodes=[30363,30364],
        physical_seam_identity="corrected_Rminus*Q_yminus=epsilon^2*rho(S)^2*b_rho(S); b=Qpsi is epsilon-stripped, S=sc/2",
        shifted_inlet_identity="H=M15(sc/2)+Rminus*Q_yminus=core_axial+core_swirl+cut_axial+cut_swirl+core_inertia+collar_inertia",
        continuation_identity="M15(s)=H+S1_actual[sc/2,s]/2-R(s)*Q_y(s)",
        accepted_q1_inertia_geometry=record(cw1),new_q2_inertia_geometry=record(cw2),positive_inertia_radius_ratio_upper=record(ratio),
        same_error_function_at_shared_endpoint=True,endpoint_contraction_before_interval_boxing=True,
        lower_leading_Vy_kept_as_actual_source=True,leading_Qaxis_and_own_mean_not_reset=True,
        physical_F0_Z_connection_in_inherited_coefficients=True,no_second_F0_derivative_applied=True,
        old_six_numeric_terms_and_all_error_caps_unchanged=True,
        new_representation_only=True,one_saved_Z_tube_only=True,downstream_actual_S1_evaluated=False,
        first_micro_upper_R_Qy_evaluated=False,full_forward_Rp=False,all_Z_functional_closure=False,
        global_five_moment_closure=False,full_goal_complete=False,accepted_ancestors_replayed=False)
