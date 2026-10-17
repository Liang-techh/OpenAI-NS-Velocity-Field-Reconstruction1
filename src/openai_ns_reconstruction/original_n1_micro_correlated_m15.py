"""True first/second-micro pressure in a correlated six-I1 basis.

Accepted free-Z endpoint fingerprints are consumed, not replayed.  New
volume leaves contain U,Q and slow Z derivatives, with no total-y Q_y.
The physical outer R*Q_y is retained as an exact source function.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sympy as sp

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_"
PARENT=PREFIX+"n1_first_micro_prefix.json"
PARENT_SHA256="60f226cf67cd7f9ce07076459dd60ef981f98c0522a89e825b4fcb71c6b5d283"
TRUE_SOURCE=PREFIX+"modified_Q4_pressure_C3_functions.py"
TRUE_GRAPH=TRUE_SOURCE[:-3]+".json.gz"
TRUE_CHECK=TRUE_SOURCE[:-3]+"_check.json"
LOOP=PREFIX+"prescribed_loop_C5_functions.json.gz"
LOOP_SOURCE=PREFIX+"prescribed_loop_C5_functions.py"
SEAM=PREFIX+"micro_active_radial4_pressure_C3_functions.json.gz"
SEAM_CHECK=SEAM[:-8]+"_check.json"
IDENTITY=PREFIX+"micro_viscous_m15_numeric.json"
IDENTITY_SOURCE="src/openai_ns_reconstruction/original_micro_viscous_m15.py"
HELPER="src/openai_ns_reconstruction/original_signed_graph_pairs.py"
HELPER_SHA256="1e9699cd84a1c82dc7b6577f922f1da1362cbedcfef148e7cec1d146a99aaaff"
DEFECT="src/openai_ns_reconstruction/original_micro_finiteN_volume_defect.py"
DEFECT_SHA256="90e0d522c0591b4162653c0367efcc1c015497025972a10064894128beb31318"


def exact_algebra():
    y,z,de,R=sp.symbols("y Z delta R",real=True)
    U=sp.Function("U")(y,z);Q=sp.Function("Q")(y,z);m=sp.Function("m")(y,z)
    h=1-de;d=1-z*z;L=1-de*z*z
    qdef=(2*z*U-h*z*m-d*sp.diff(m,z))/L
    rules={sp.diff(m,y):U-m,sp.diff(m,y,z):sp.diff(U,z)-sp.diff(m,z)}
    qp=((1+de)*z*U+2*z*sp.diff(U,y)-d*sp.diff(U,z))/L
    recovery=(qdef+sp.diff(qdef,y)-qp).subs(rules)
    direct=(Q+h*z*sp.diff(Q,z)/2+sp.diff(Q,y))/L+Q*(Q/2+sp.diff(Q,y))
    direct+=U*(-2*z*Q+d*sp.diff(Q,z)-2*z*sp.diff(Q,y))/L
    primitive=Q/L+Q*Q/2-2*z*z*U*U/(L*L)
    density=-Q/L+h*z*sp.diff(Q,z)/(2*L)-Q*Q/2+d*U*sp.diff(Q,z)/L
    density+=2*h*z*z*U*U/(L*L)+2*z*d*U*sp.diff(U,z)/(L*L)
    defect=direct-(2*primitive+sp.diff(primitive,y)+density)
    factored=2*z*U*(Q+sp.diff(Q,y)-qp)/L
    E0,Es,E2,B0,Bs,B2,J1,J2=sp.symbols("E0 Es E2 B0 Bs B2 J1 J2")
    checks=dict(physical_Q_mean_recovery=sp.cancel(recovery)==0,
        physical_Q_mean_recovery_Z1=sp.cancel(sp.diff(recovery,z))==0,
        exact_R_squared_pressure_IBP=sp.cancel(defect+factored)==0,
        exact_R_squared_pressure_IBP_Z1=sp.cancel(sp.diff(defect+factored,z))==0,
        same_function_viscous_seam=sp.expand((Es-E0)+(E2-Es)-(E2-E0))==0,
        same_function_reduced_boundary_seam=sp.expand((Bs-B0+J1)+(B2-Bs+J2)-(B2-B0+J1+J2))==0)
    if not all(checks.values()):raise ValueError("True micro correlated pressure algebra failed")
    return dict(passed=True,checks=checks,
        direct_minus_reduced="-2*Z*U/L*(Q+Q_y-[(1+delta)ZU+2ZUy-dUZ]/L)",
        physical_radial_weight_derivative="d_y(R^2)=2R^2",
        no_Z_factor_on_U_QZ_volume_coefficient=True)


def seam_certificate(root,parent,graph,loop,pins,load,require,digest,ast_guard):
    active,check,identity=[load(root/p) for p in (SEAM,SEAM_CHECK,IDENTITY)]
    require(active["source_family"]==check["source_family"]==identity["source_family"]==parent["source_family"] and
        check["all_passed"] is True and check["current_original_n1_first_second_micro_full_mixed_source_join_installed"] is True and
        check["input_hashes"][Path(SEAM).name]==pins[SEAM]["expected_sha256"] and
        active["exact_graph_nodes"][:len(graph["exact_graph_nodes"])]==graph["exact_graph_nodes"],
        "Accepted same-source micro join must extend this exact Q4 graph")
    c=check["checks"]["actual_leading_micro_sources"]
    require(c["actual_first_complete_Duhamel_history_is_second_incoming"] is True and
        c["actual_first_complete_pressure_integrals_are_second_incoming"] is True and
        c["common_background_roots_phase_radius_and_periodic_loop_sources"] is True and
        c["leading_one_minus_delta_Q_gives_all_high_Q_rows_from_same_U_mean"] is True,
        "Original complete histories, phase, loop and leading Q convention required")
    checks=identity["checks"];trace=identity["samples"][0]["trace_identity"]
    require(checks["passed"] is True and checks["source_identity_rows_per_packet"]==14 and
        checks["direct_Qy_and_Z1_source_ids"] is True and checks["exact_cancelled_interface"]=="first_s1" and
        trace["passed"] is True and trace["signed_trace_centers_evaluated"] is False and
        trace["upper_caps_not_used_as_function_values"] is True and
        trace["exact_sc_tuple"]==parent["actual_native_collar_transport"]["domain"]["sc_exact_mpf_tuple"] and
        all(p["trace_identity"]==trace for p in identity["samples"]),
        "Use the accepted free-Z source fingerprints, never saved sample overlap")
    b=graph["modified_Q4_pressure_C3_bindings"];lb=loop["prescribed_loop_C5_bindings"]
    first,second=[lb["charts"][k] for k in ("first_micro","second_micro")]
    for key in ("coordinate","radius","physical_measure","actual_phase","roots","global_offset"):
        require(first[key]==second[key],"Actual chart seam source mismatch: "+key)
    require(first["domain"][1]==second["domain"][0]==1 and first["coordinate"]==trace["coordinate_node"] and
        first["radius"]==trace["radius_node"] and first["actual_phase"]==trace["phase_node"] and
        first["physical_measure"]==trace["physical_measure_node"] and lb["N"]==trace["selected_integer_node"]==18 and
        trace["seam_node"]==1 and trace["common_ordinary_Z_independent_radius"] is True and
        graph["exact_graph_nodes"][1]==dict(operation="exact_rational",numerator=1,denominator=1) and
        second["incoming"]==lb["endpoints"]["first_micro"] and
        b["charts"]["second_micro"]["incoming_integrals"]==b["endpoints"]["first_micro"]["integrals"],
        "Exact s=1 seam with complete Duhamel and pressure predecessor required")
    expected=[]
    for field in ("U","m","Q"):
        for radial in (0,1):
            key=field if radial==0 else field+"_y"
            for k in (0,1):
                expected.append((field,radial,k,b["charts"]["first_micro"]["physical"][key][k],
                    b["charts"]["second_micro"]["physical"][key][k]))
    for k in (0,1):expected.append(("R_Q_y",0,k,b["charts"]["first_micro"]["finite_boundary"]["R_Q_y"][k],
        b["charts"]["second_micro"]["finite_boundary"]["R_Q_y"][k]))
    rows=trace["source_identity_rows"]
    require(len(rows)==14 and [(p["field"],p["physical_y_order"],p["ordinary_Z_order"],p["left_node"],p["right_node"]) for p in rows]==expected and
        all(len(p["canonical_endpoint_sha256"])==64 for p in rows),"The accepted 14 source rows must be these exact Q4 field functions")
    guard=ast_guard(root,IDENTITY_SOURCE,"at",(
        "scope=scope if scope is not None else {self.b['coordinate']:self.one}",
        "elif op=='bound_variable':out=self.digest([op,self.name(row['name'])]),None"))
    return dict(passed=True,accepted_fingerprint_report_sha256=pins[IDENTITY]["expected_sha256"],
        source_receipt_sha256=pins[SEAM_CHECK]["expected_sha256"],source_identity_rows=rows,
        source_graph_prefix_and_all_14_field_ports_checked=True,free_Z_not_numeric_Z_substitution=True,
        ordinary_Z_orders=[0,1],physical_orders=[0,1],exact_s1_seam=True,
        complete_predecessor_history_retained=True,shared_phase_N_radius_and_Pstar=True,
        fingerprint_scope_AST=guard,accepted_fingerprint_computation_replayed=False,
        any_future_endpoint_error_cancellation_requires_one_shared_model=True)


def build(root,parent):
    root=Path(root)
    if hashlib.sha256((root/PARENT).read_bytes()).hexdigest()!=PARENT_SHA256 or json.loads((root/PARENT).read_bytes())!=parent:
        raise ValueError("Accepted first-micro numerical prefix and true source required")
    helper=root/HELPER
    if hashlib.sha256(helper.read_bytes()).hexdigest()!=HELPER_SHA256:raise ValueError("Changed new graph arithmetic helper")
    if hashlib.sha256((root/DEFECT).read_bytes()).hexdigest()!=DEFECT_SHA256:raise ValueError("Changed new signed volume defect helper")
    common="src/openai_ns_reconstruction/original_finite_N_signed_source.py"
    if hashlib.sha256((root/common).read_bytes()).hexdigest()!=parent["implementation_hashes"][common]:
        raise ValueError("Changed accepted source pin/AST helper")
    from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
    from .original_signed_graph_pairs import GraphPairs
    from .original_micro_finiteN_volume_defect import build_defect
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256,passed=True),
        HELPER:dict(path=HELPER,expected_sha256=HELPER_SHA256,passed=True),
        DEFECT:dict(path=DEFECT,expected_sha256=DEFECT_SHA256,passed=True)}
    for path in (TRUE_SOURCE,TRUE_GRAPH,TRUE_CHECK,LOOP,LOOP_SOURCE,SEAM,SEAM_CHECK,IDENTITY,IDENTITY_SOURCE,common):
        h=parent["input_hashes"][Path(path).name] if path.startswith(BASE) else parent["implementation_hashes"][path]
        require(sha(root/path)==h,"Changed true correlated micro source: "+path)
        pins[path]=dict(path=path,expected_sha256=h,passed=True)
    graph,check,loop=[load(root/p) for p in (TRUE_GRAPH,TRUE_CHECK,LOOP)]
    require(parent["checks"]["passed"] is True and graph["source_family"]==check["source_family"]==loop["source_family"]==parent["source_family"] and
        check["all_passed"] is True and check["input_hashes"][Path(TRUE_GRAPH).name]==pins[TRUE_GRAPH]["expected_sha256"] and
        check["checks"]["actual_source"]["true_phase_held_native_function_derivatives_and_total_fast_terms"] is True,
        "Actual finite-N mixed Q4 fields required")
    seam=seam_certificate(root,parent,graph,loop,pins,load,require,digest,ast_guard)
    nodes=graph["exact_graph_nodes"];b=graph["modified_Q4_pressure_C3_bindings"]
    require(nodes[:len(loop["exact_graph_nodes"])]==loop["exact_graph_nodes"],
        "Prescribed-loop field handles must belong to the identical Q4 graph prefix")
    lb=loop["prescribed_loop_C5_bindings"];g=GraphPairs(nodes)
    continuation=parent["actual_first_micro_prefix"]["full_first_micro_continuation_source"]
    require(continuation["passed"] is True and continuation["source_graph_sha256"]==pins[TRUE_GRAPH]["expected_sha256"] and
        continuation["common_N"]==nodes[18] and nodes[18]["type"]=="DoubleDyadicInteger" and
        nodes[18]["J"]==1358356628656378313 and nodes[18]["Z_independent"] is True and
        continuation["shifted_inlet_source_digest"]==digest(parent["actual_first_micro_seam"]) and
        [v["node"] for v in continuation["actual_initial_R_Qy_ordinary_Z"]]==
            b["initial"]["R_Q_y"][:2] and
        all(nodes[i]==dict(operation="exact_rational",numerator=0,denominator=1)
            for i in b["charts"]["first_micro"]["incoming_integrals"]["S1"]),
        "Bind the accepted physical shifted inlet to this same Q4 initial boundary")
    z=(20,g.constant(1));de=g.fixed(25220);one=g.fixed(g.constant(1));two=g.constant(2)
    mul=g.pmul;add=g.padd;neg=g.pneg;fixed=g.fixed
    d=add(one,neg(mul(z,z)));h=add(one,neg(de));L=add(one,neg(mul(de,z,z)));L2=mul(L,L)
    divide=lambda p,q=L:g.pdiv(p,q,"same actual positive physical L=1-delta*Z^2")
    half=lambda p:g.pscale(p,g.constant("1/2"))
    integrals={};local={};endpoint_fields={};defects={}
    names=("q","q_Z","q_squared","U_q_Z","U_squared","U_U_Z")
    for name in ("first_micro","second_micro"):
        f=b["charts"][name]["physical"];chart=lb["charts"][name]
        background_mean=[]
        for k in range(4):
            v=chart["modified_fields"]["V"][k];mv=chart["modified_fields"]["m"][k]
            require(nodes[v]==dict(operation="sum",arguments=[chart["roots"]["V"][k],chart["finite_N"]["delta_V"][k]]) and
                nodes[mv]["operation"]=="sum" and len(nodes[mv]["arguments"])==2 and
                nodes[mv]["arguments"][1]==chart["history"]["m"][k],
                "Actual V/m source must retain original background plus nonreset correction history")
            background_mean.append(nodes[mv]["arguments"][0])
        defects[name]=build_defect(g,f,chart)
        defects[name]["background_mean_ordinary_Z0_through_Z3_source_nodes"]=background_mean
        U,U_Z=tuple(f["U"][:2]),tuple(f["U"][1:3]);Q,Q_Z=tuple(f["Q"][:2]),tuple(f["Q"][1:3])
        source=dict(q=Q,q_Z=Q_Z,q_squared=mul(Q,Q),U_q_Z=mul(U,Q_Z),U_squared=mul(U,U),U_U_Z=mul(U,U_Z))
        radius=chart["radius"];weight=fixed(g.mul(radius,radius,chart["physical_measure"]))
        local[name]={}
        for key in names:
            value=g.integral(mul(weight,source[key]),chart["coordinate"],*chart["domain"],
                "actual_micro_correlated_"+name+"_I1_"+key,dict(original_chart=name,physical_radius_power=2,
                    physical_Jacobian_applied_once=True,no_total_y_derivative_in_volume=True))
            local[name][key]=value;integrals[key]=add(integrals.get(key,fixed(g.constant(0))),value)
        endpoint_fields[name]={}
        for side,end in zip(("lower","upper"),chart["domain"]):
            at=lambda pair:g.at(pair,chart["coordinate"],end)
            R2=fixed(g.mul(radius,radius))
            endpoint_fields[name][side]=dict(R_squared_q=at(mul(R2,Q)),
                R_squared_q_squared=at(mul(R2,mul(Q,Q))),R_squared_U_squared=at(mul(R2,mul(U,U))),
                R_Q_y=at(tuple(b["charts"][name]["finite_boundary"]["R_Q_y"][:2])))
    lower=endpoint_fields["first_micro"]["lower"];upper=endpoint_fields["second_micro"]["upper"]
    boundary=lambda key:add(upper[key],neg(lower[key]))
    B=add(divide(boundary("R_squared_q")),half(boundary("R_squared_q_squared")),
        neg(g.pscale(divide(mul(z,z,boundary("R_squared_U_squared")),L2),two)))
    D=add(neg(divide(integrals["q"])),half(divide(mul(h,z,integrals["q_Z"]))),
        neg(half(integrals["q_squared"])),divide(mul(d,integrals["U_q_Z"])),
        g.pscale(divide(mul(h,z,z,integrals["U_squared"]),L2),two),
        g.pscale(divide(mul(z,d,integrals["U_U_Z"]),L2),two))
    S1=add(B,D);increment=add(half(S1),neg(boundary("R_Q_y")))
    shifted_tail=add(half(S1),neg(upper["R_Q_y"]))
    guards=dict(
        finite_N_source=ast_guard(root,LOOP_SOURCE,"finite_N_functions",(
            "s=a.scale(loop['B_over_Pstar'],invN)","V_N=a.add(V,s)","delta_V=s")),
        nonreset_mean_history=ast_guard(root,LOOP_SOURCE,"transport",(
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "predecessor_memory_never_reset=True")),
        physical_mean=ast_guard(root,TRUE_SOURCE,"physical_pressure_fields",(
            "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)","my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))")),
        true_pressure_density=ast_guard(root,TRUE_SOURCE,"reduced_pressure_S",(
            "quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))",
            "axial=div(a.mul(U,a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2))))")))
    patch=dict(base_archive=TRUE_GRAPH,base_sha256=pins[TRUE_GRAPH]["expected_sha256"],base_graph_node_count=len(nodes),
        appended_nodes=g.appended,local_I1=local,joint_I1=integrals,endpoint_fields=endpoint_fields,
        finite_N_signed_integrand_defects=defects,
        reduced_boundary=B,reduced_volume=D,reduced_S1=S1,actual_micro_pair_increment=increment,
        shifted_inlet_tail=shifted_tail,output_formula="M15(s=2)=H+shifted_inlet_tail",
        actual_Q4_direct_cumulative_S1=b["endpoints"]["second_micro"]["integrals"]["S1"][:2],
        actual_Q4_direct_cumulative_m15=b["endpoints"]["second_micro"]["pressure"]["m15"][:2],
        exactly_one_representation_consumed=True,old_graph_unchanged=True)
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,source_family=parent["source_family"],
        source_pins=pins,source_AST=guards,exact_algebra=exact_algebra(),seam_certificate=seam,
        graph_patch=patch,graph_patch_digest=digest(patch),ordinary_Z_orders=[0,1],native_interval=["sc/2","2"],
        shifted_inlet_source_digest=digest(parent["actual_first_micro_seam"]),
        actual_full_Q4_six_I1_source_implemented=True,bulk_Qy_and_Uy_eliminated_by_exact_IBP=True,
        full_slow_and_fast_finiteN_functions_retained_at_outer_endpoint=True,
        common_N=nodes[18],six_I1_names=list(names),physical_I1_measure="R^2*hb*ds exactly once",
        internal_RQy_and_reduced_B_seams_cancel_as_functions=True,
        viscous_endpoint_source_weights=dict(first=dict(lower=1,seam=-1),second=dict(seam=1,upper=-1),joint=dict(lower=1,upper=-1)),
        same_source_error_function_shared_before_any_future_boxing=True,
        old_six_numeric_terms_and_all_error_caps_unchanged=True,accepted_seam_fingerprints_replayed=False,
        numerical_signed_joint_integral_evaluated=False,numerical_full_micro_pair_enclosure=False,
        all_Z_five_moment_closure=False,full_forward_Rp=False,full_goal_complete=False,accepted_ancestors_replayed=False)
