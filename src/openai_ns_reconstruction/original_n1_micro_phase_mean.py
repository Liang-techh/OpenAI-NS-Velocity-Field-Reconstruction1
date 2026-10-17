"""Actual micro mean with periodic phase primitives and nonreset histories.

The saved six-I1 source is consumed without rebuilding it.  All four
ordinary Z rows of the physical mean remain exact functions; the fixed-phase
slow derivative convolution is retained explicitly rather than boxed away.
"""
from __future__ import annotations
import ast
import hashlib
import json
from fractions import Fraction
from pathlib import Path

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_"
PARENT=PREFIX+"n1_micro_correlated_m15.json"
PARENT_SHA256="7868f1899398a54aaeba1b799b8b1015465718445444966648295e66ae6c97aa"
TRUE=PREFIX+"modified_Q4_pressure_C3_functions.json.gz"
TRUE_SOURCE=TRUE[:-8]+".py"
LOOP=PREFIX+"prescribed_loop_C5_functions.json.gz"
LOOP_SOURCE=LOOP[:-8]+".py"
LOOP_CHECK=LOOP[:-8]+"_check.json"
PARITY=BASE+"lei_ren_part1_paper_compliant_current_native_signed_averaging_check.json"
HELPER="src/openai_ns_reconstruction/original_micro_phase_mean_kernel.py"
HELPER_SHA256="3c22ebb63788ab1b8dd50dd8fa5923a67f579ec464ef00d6a9ec81eea7de1181"
ARITHMETIC="src/openai_ns_reconstruction/original_signed_graph_pairs.py"
C3=PREFIX+"C3_source_targets.json.gz"
C3_RANGE=PREFIX+"C3_target_ranges.json.gz"
FINGERPRINT="src/openai_ns_reconstruction/original_micro_viscous_m15.py"


def fixed_phase_guard(root,require,digest):
    source=(root/TRUE_SOURCE).read_text(encoding="utf-8")
    tree=ast.parse(source)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=="NativeDerivative")
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=="__init__")
    text="".join(ast.get_source_segment(source,fn).split())
    require("self.cache={phase.node:g.zero}" in text,"Fixed-phase derivative must annihilate the issued phase")
    return dict(passed=True,path=TRUE_SOURCE,class_name="NativeDerivative",method="__init__",
        phase_held_exactly_fixed=True,AST_digest=digest(ast.dump(fn,include_attributes=False)),source_executed=False)


def physical_prefactor(parent,N,require):
    """Directed logarithm of an actual coefficient, never a field/error value."""
    from mpmath.ctx_iv import MPIntervalContext
    from mpmath import make_mpf,nstr
    c=MPIntervalContext();c.prec=600
    basis=parent["actual_signed_lower_m15"]
    require(basis["source_log_basis_names"][2]=="Pstar_squared" and N["definition"]=="N=2^(2^J)",
        "Use the saved actual global Pstar and unchanged common integer")
    row=basis["source_log_bases"][2]
    logP2=c.mpf([make_mpf(tuple(row[k+"_exact_mpf_tuple"])) for k in ("lower","upper")])
    logN=c.mpf([make_mpf((sgn,man,e+N["J"],bits)) for sgn,man,e,bits in c.ln(2)._mpi_])
    record=lambda x:dict(lower_exact_mpf_tuple=list(x._mpi_[0]),upper_exact_mpf_tuple=list(x._mpi_[1]),directed=True)
    factor=logP2/2-2*logN
    require(factor.b<0,"Expected negative logarithm of actual Pstar/N^2")
    return dict(passed=True,log_Pstar_over_N_squared=record(factor),log_common_N=record(logN),
        log_Pstar_squared=record(logP2),upper_log_display=nstr(make_mpf(factor._mpi_[1]),20),
        exact_N_unchanged=True,coefficient_not_mean_value_or_error_radius=True,
        Pstar_over_N_squared_materialized=False,mean_smallness_does_not_follow_without_bracket_bound=True)


def recover_Q(g,U,m):
    """Physical Q0..2 from U0..2 and m0..3; ordinary derivatives."""
    z,de=20,25220;one=g.constant(1);two=g.constant(2)
    sub=lambda x,y:g.add(x,g.neg(y))
    h=sub(one,de);d=sub(one,g.mul(z,z));L=sub(one,g.mul(de,z,z))
    L1=g.neg(g.mul(two,de,z));L2=g.neg(g.mul(two,de))
    n0=sub(g.mul(two,z,U[0]),g.add(g.mul(h,z,m[0]),g.mul(d,m[1])))
    n1=g.add(g.mul(two,U[0]),g.mul(two,z,U[1]),g.neg(g.mul(h,m[0])),
        g.mul(g.add(one,de),z,m[1]),g.neg(g.mul(d,m[2])))
    n2=g.add(g.mul(g.constant(4),U[1]),g.mul(two,z,U[2]),g.mul(two,de,m[1]),
        g.mul(g.add(g.constant(3),de),z,m[2]),g.neg(g.mul(d,m[3])))
    cert="same actual positive L=1-delta*Z^2"
    q0=g.div(n0,L,cert);q1=g.div(sub(n1,g.mul(L1,q0)),L,cert)
    q2=g.div(sub(n2,g.add(g.mul(two,L1,q1),g.mul(L2,q0))),L,cert)
    return [q0,q1,q2]


def recovery_controls(require):
    import sympy as s
    z,de=s.symbols("Z delta");U=s.Function("U")(z);m=s.Function("m")(z)
    h,d,L=1-de,1-z*z,1-de*z*z
    n=2*z*U-h*z*m-d*s.diff(m,z)
    n1=2*U+2*z*s.diff(U,z)-h*m+(1+de)*z*s.diff(m,z)-d*s.diff(m,z,2)
    n2=4*s.diff(U,z)+2*z*s.diff(U,z,2)+2*de*s.diff(m,z)+(3+de)*z*s.diff(m,z,2)-d*s.diff(m,z,3)
    q0=n/L;q1=(n1-s.diff(L,z)*q0)/L
    q2=(n2-2*s.diff(L,z)*q1-s.diff(L,z,2)*q0)/L
    checks=[s.cancel(q-s.diff(n/L,z,k))==0 for k,q in enumerate((q0,q1,q2))]
    require(all(checks),"Physical ordinary Q recovery failed")
    return dict(passed=True,ordinary_Z0_Z1_Z2=checks,mean_Z3_retained=True,no_factorial_normalization=True)


def legacy_diagnostic(parent,legacy,ranges,loop,pins,require):
    """New arithmetic on saved caps; explicitly not a true-C5 source bound."""
    from mpmath.ctx_iv import MPIntervalContext
    from mpmath import make_mpf,nstr
    c=MPIntervalContext();c.prec=600
    read=lambda r:c.mpf([make_mpf(tuple(r[k+"_exact_mpf_tuple"])) for k in ("lower","upper")])
    record=lambda v:dict(lower_exact_mpf_tuple=list(v._mpi_[0]),upper_exact_mpf_tuple=list(v._mpi_[1]),directed=True)
    pref=physical_prefactor(parent,loop["exact_graph_nodes"][18],require)
    logN=read(pref["log_common_N"]);logP=read(pref["log_Pstar_squared"])/2
    loghb=read(parent["actual_signed_lower_m15"]["source_log_bases"][0])
    frames=[f for f in ranges["actual_four_Z_C3_target_transports"]
        if Fraction(f["exact_Z_cell"][0])<=Fraction(371,1000)<=Fraction(f["exact_Z_cell"][1])]
    require(len(frames)==1,"One saved axial cell containing the actual query required")
    out={}
    for name in ("first_micro","second_micro"):
        cells=next(x for x in frames[0]["actual_17_chart_C3_transports"] if x["chart"]==name)["cells"]
        logs=[read(x["actual_source_and_primitive_range"]["B"][3]["log_absolute_upper"]) for x in cells
            if not x["actual_source_and_primitive_range"]["B"][3]["exact_zero"]]
        require(logs,"Retain the genuine nonflat legacy cap")
        cover=c.mpf([max(x.a for x in logs).a,max(x.b for x in logs).b])
        direct=logP+2*loghb+cover-logN
        out[name]=dict(legacy_B_Z3_log_cover=record(cover),
            log_Pstar_hb_squared_B3_over_N=record(direct),
            upper_log_display=nstr(make_mpf(direct._mpi_[1]),20),
            radius_small=bool(direct.b<0),first_incoming_memory_in_second_not_reset=True,
            legacy_B_nodes=legacy["actual_C3_primitive_functions"][name]["B"],
            actual_C5_B_nodes=loop["prescribed_loop_C5_bindings"]["charts"][name]["loop"]["B_over_Pstar"][:4])
    require(out["first_micro"]["radius_small"] is False and out["second_micro"]["radius_small"] is True,
        "Legacy weighted cap diagnostic changed; inspect it rather than claiming a closure")
    return dict(passed=True,saved_Z_cell=frames[0]["exact_Z_cell"],ordinary_Z_order=3,charts=out,
        legacy_source_archive=C3,legacy_source_sha256=pins[C3]["expected_sha256"],
        cap_archive=C3_RANGE,cap_archive_sha256=pins[C3_RANGE]["expected_sha256"],
        numerical_bound_on_actual_C5_mean=False,same_function_cap_transfer_implemented=False,
        different_node_numbers_alone_do_not_prove_different_functions=True,
        source_family_equality_does_not_admit_cap_transfer=True,caps_used_as_signed_values=False,
        only_new_directed_arithmetic_on_saved_caps=True,old_producers_replayed=False)


def primitive_seam(nodes,tb,lb,require,digest):
    """Prove eight NEW B/slow-By traces; never rerun the accepted14 rows."""
    from .original_micro_viscous_m15 import EndpointFingerprint
    first,second=[lb["charts"][k] for k in ("first_micro","second_micro")]
    canonical=EndpointFingerprint(nodes,dict(coordinate=first["coordinate"],loop=lb["charts"]))
    phi=(canonical.digest(["bound_variable","new_free_issued_phase_phi"]),None)
    scope={first["coordinate"]:canonical.one,first["actual_phase"]:phi}
    rows=[]
    for field in ("B","slow_B_y"):
        if field=="B":left,right=first["loop"]["B_over_Pstar"],second["loop"]["B_over_Pstar"]
        else:left,right=[tb["charts"][k]["radial"]["slow_B_y"] for k in ("first_micro","second_micro")]
        for k in range(4):
            a,b=canonical.at(left[k],scope),canonical.at(right[k],scope)
            require(a[0]==b[0],"New free-Z/free-phase primitive seam differs: "+field+" Z"+str(k))
            rows.append(dict(field=field,ordinary_Z_order=k,left_node=left[k],right_node=right[k],
                canonical_source_sha256=a[0]))
    require(len(rows)==8 and canonical.collar,"Actual first-only collar must be discharged by its source guard")
    return dict(passed=True,new_source_rows=rows,rows=8,ordinary_Z_orders=[0,1,2,3],
        native_s_substituted_with_exact_one=True,issued_phase_kept_free=True,Z_kept_free=True,
        actual_phase_sample_not_evaluated=True,first_collar_predicates_discharged=len(canonical.collar),
        same_B_and_slow_By_imply_same_anchored_C_and_Cy=True,
        row_digest=digest(rows),accepted14_row_fingerprint_computation_replayed=False)


def build(root,parent):
    root=Path(root)
    if hashlib.sha256((root/PARENT).read_bytes()).hexdigest()!=PARENT_SHA256 or json.loads((root/PARENT).read_bytes())!=parent:
        raise ValueError("Accepted actual micro-pair graph required")
    if hashlib.sha256((root/HELPER).read_bytes()).hexdigest()!=HELPER_SHA256:
        raise ValueError("Changed new micro phase-mean kernel")
    common="src/openai_ns_reconstruction/original_finite_N_signed_source.py"
    for path in (ARITHMETIC,common,FINGERPRINT):
        if hashlib.sha256((root/path).read_bytes()).hexdigest()!=parent["implementation_hashes"][path]:
            raise ValueError("Changed accepted pure graph/pin helper")
    from .original_finite_N_signed_source import load,require,digest,sha,ast_guard
    from .original_signed_graph_pairs import GraphPairs
    from .original_micro_phase_mean_kernel import build_micro_mean
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256,passed=True),
        HELPER:dict(path=HELPER,expected_sha256=HELPER_SHA256,passed=True)}
    for path in (TRUE,TRUE_SOURCE,LOOP,LOOP_SOURCE,LOOP_CHECK,PARITY,C3,C3_RANGE,ARITHMETIC,common,FINGERPRINT):
        if path==PARITY:
            h=parent["actual_C5_signed_history_IBP"]["source_pins"]["zero_mean_receipt"]["expected_sha256"]
        else:h=parent["input_hashes"][Path(path).name] if path.startswith(BASE) else parent["implementation_hashes"][path]
        require(sha(root/path)==h,"Changed actual phase-mean source: "+path)
        pins[path]=dict(path=path,expected_sha256=h,passed=True)
    true,loop,receipt,parity=[load(root/p) for p in (TRUE,LOOP,LOOP_CHECK,PARITY)]
    a=parent["actual_micro_correlated_m15"];patch=a["graph_patch"]
    require(parent["checks"]["passed"] is True and digest(patch)==a["graph_patch_digest"] and
        patch["base_sha256"]==pins[TRUE]["expected_sha256"] and
        a["seam_certificate"]["passed"] is True and a["accepted_seam_fingerprints_replayed"] is False and
        true["exact_graph_nodes"][:len(loop["exact_graph_nodes"])]==loop["exact_graph_nodes"],
        "Consume the actual shared source graph and accepted seam, never rebuild it")
    require(all(r["source_family"]==parent["source_family"] for r in (true,loop,receipt,parity)) and
        receipt["all_passed"] is True and parity["all_passed"] is True,
        "Same actual source family and accepted loop/parity receipts required")
    theorem=parity["exact_source_averaging_theorem"]
    ids=theorem["original_primitive_parity_theorem"]["exact_identities"]
    require(theorem["passed"] is True and ids["B_reflection_zero_mean"] is True and
        ids["periodic_B_derivative"] is True and theorem["slow_y_not_total_y"] is True,
        "Original signed reflection identity and fixed-phase differentiation required")
    lb=loop["prescribed_loop_C5_bindings"];tb=true["modified_Q4_pressure_C3_bindings"]
    first,second=[lb["charts"][k] for k in ("first_micro","second_micro")]
    nodes=true["exact_graph_nodes"]+patch["appended_nodes"]
    from mpmath.libmp import mpf_cmp
    sc=tuple(parent["actual_native_collar_transport"]["domain"]["sc_exact_mpf_tuple"])
    require(mpf_cmp(sc,(0,0,0,0))>0 and mpf_cmp(sc,(0,1,-2,1))<0,
        "Original positive sc<1/4 makes the first inlet flat and the s=1 seam outside its collar")
    def scaled_sc(node,coefficient):
        row=nodes[node]
        if row["operation"]!="product" or len(row["arguments"])!=2 or lb["sc"] not in row["arguments"]:return False
        other=next(i for i in row["arguments"] if i!=lb["sc"]);v=nodes[other]
        return v["operation"]=="exact_rational" and Fraction(v["numerator"],v["denominator"])==coefficient
    require(scaled_sc(first["domain"][0],Fraction(1,2)) and
        scaled_sc(first["loop"]["original_collar_upper"],Fraction(3,4)) and Fraction(1,2)<Fraction(3,4),
        "The actual chart lower is sc/2 strictly within its source-owned 3sc/4 flat collar")
    require(patch["base_graph_node_count"]==len(true["exact_graph_nodes"]) and lb["N"]==18 and
        nodes[18]["Z_independent"] is True and
        first["roots"]==second["roots"] and first["actual_phase"]==second["actual_phase"] and
        first["coordinate"]==second["coordinate"] and first["radius"]==second["radius"] and
        first["physical_measure"]==second["physical_measure"] and
        first["domain"][1]==second["domain"][0]==1 and second["incoming"]==lb["endpoints"]["first_micro"] and
        all(c["total_Z_phase_derivative_exact_zero"] is True for c in (first,second)) and
        all(nodes[n]==dict(operation="exact_rational",numerator=0,denominator=1) for n in first["incoming"]["m"]),
        "One true phase, radius, Jacobian and nonreset correction history required")
    guards=dict(
        fixed_phase_derivative=fixed_phase_guard(root,require,digest),
        loop_formula=ast_guard(root,LOOP_SOURCE,"loop_functions",(
            "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))",
            "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
            "collar=g.mul(c('3/4'),sc) if chart=='first_micro' else None",
            "operator='<=',left=coordinate.node,right=collar.node")),
        root_relation=ast_guard(root,LOOP_SOURCE,"generic_roots",(
            "aa,bb=divide(C),divide(B)","t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))")),
        slow_source=ast_guard(root,TRUE_SOURCE,"total_radial_functions",(
            "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref)) for ref in rows(chart['loop'][key])))",
            "slow_B_y=slow['B_over_Pstar']")))
    radius=nodes[first["radius"]];phase=nodes[first["actual_phase"]]
    require(radius["operation"]=="product" and len(radius["arguments"])==2,
        "Original micro radius must be Ra*exp(hb*s)")
    Ra,ex=radius["arguments"]
    require(nodes[ex]["operation"]=="analytic_unary" and nodes[ex]["name"]=="exp" and
        nodes[nodes[ex]["argument"]]==dict(operation="product",arguments=[first["physical_measure"],first["coordinate"]]) and
        phase["operation"]=="analytic_unary" and phase["name"]=="fractional_part" and
        nodes[phase["argument"]]==dict(operation="product",arguments=[18,first["global_offset"]]),
        "Actual log-radius phase derivative is N, not N/hb or a chosen angle")
    offset=nodes[first["global_offset"]]
    logref,negative=offset["arguments"]
    require(offset["operation"]=="sum" and len(offset["arguments"])==2 and
        nodes[logref]["operation"]=="analytic_unary" and nodes[logref]["name"]=="log" and
        nodes[nodes[logref]["argument"]]["operation"]=="positive_quotient" and
        nodes[nodes[logref]["argument"]]["numerator"]==first["radius"] and
        nodes[nodes[logref]["argument"]]["denominator"]==Ra and
        nodes[negative]==dict(operation="negative",argument=lb["phase_origin"]),
        "Use the actual global log(R/Ra)-phase_origin")
    definition=lb["definition"]
    require(definition["original_first_micro_start_sc_half_and_collar_three_sc_quarters"] is True and
        definition["one_original_common_N_eta_dstar_and_global_log_R_rminus_phase"] is True and
        a["seam_certificate"]["exact_s1_seam"] is True,
        "Native first inlet is inside the source-owned flat collar; s=1 is outside it")
    require(definition["numeric_C5_containment_and_mixed_y_bounds_not_certified"] is True,
        "Do not silently transfer older C3 caps into genuine C5 mixed-source containment")
    seam=primitive_seam(nodes,tb,lb,require,digest)
    g=GraphPairs(nodes);kernel=build_micro_mean(g,tb,lb)
    require(kernel["passed"] is True and kernel["symbolic_controls"]["passed"] is True,
        "New source transport and nonreset seam algebra required")
    kernel.pop("graph")  # One authoritative suffix namespace, including the Q recovery below.
    from .original_micro_viscous_m15 import EndpointFingerprint
    lower_check=EndpointFingerprint(g.nodes,dict(coordinate=first["coordinate"],loop=lb["charts"]))
    zero_rows=[]
    for kind in ("JC","JCy"):
        for k,node in enumerate(kernel[kind]["second_micro"]):
            require(lower_check.at(node,{first["coordinate"]:lower_check.one})==lower_check.zero,
                "Second local convolution at its own lower endpoint must be an empty integral")
            zero_rows.append(dict(kind=kind,ordinary_Z_order=k,node=node))
    seam["new_second_local_convolutions_zero_at_s1"]=zero_rows
    scalars=kernel["graph_source_ports"]["scalars"]
    require(scalars["Pstar"]==33572 and scalars["N"]==lb["N"]==18 and
        scalars["hb"]==first["physical_measure"]==second["physical_measure"]==8 and
        nodes[33572]==dict(operation="analytic_unary",name="exp",argument=5) and
        nodes[8]["name"]=="hbB" and nodes[8]["Z_independent"] is True and
        a["common_N"]==nodes[18] and parent["actual_first_micro_prefix"]["full_first_micro_continuation_source"]["physical_Pstar_once"] is True,
        "Kernel fallbacks must equal the actual accepted physical Pstar, hb and common N")
    hb_receipt=nodes[8]["current_source_receipt"]
    require(parent["input_hashes"][hb_receipt]==nodes[8]["current_source_receipt_sha256"]==sha(root/BASE/hb_receipt),
        "Accepted positive hb source ledger receipt must match the selected parameter")
    for name in ("first_micro","second_micro"):
        for k in range(4):
            row=nodes[tb["charts"][name]["physical"]["U"][k]]
            require(row["operation"]=="product" and len(row["arguments"])==2 and
                set(row["arguments"])=={33572,lb["charts"][name]["modified_fields"]["V"][k]},
                "The selected Pstar must be the factor in the actual physical U rows")
    recovery={}
    for name in ("first_micro","second_micro"):
        U=tb["charts"][name]["physical"]["U"]
        bg=patch["finite_N_signed_integrand_defects"][name]["background_mean_ordinary_Z0_through_Z3_source_nodes"]
        mean=[g.add(g.mul(33572,bg[k]),kernel["complete_physical_delta_m"][name][k]) for k in range(4)]
        q=recover_Q(g,U,mean)
        up,uz=tuple(U[:2]),tuple(U[1:3]);qp,qz=tuple(q[:2]),tuple(q[1:3])
        integrands=dict(q=qp,q_Z=qz,q_squared=g.pmul(qp,qp),U_q_Z=g.pmul(up,qz),
            U_squared=g.pmul(up,up),U_U_Z=g.pmul(up,uz))
        recovery[name]=dict(physical_mean_ordinary_Z0_Z3=mean,physical_Q_ordinary_Z0_Z2=q,
            source_six_I1_C1_integrands=integrands,
            original_Q_source_nodes=tb["charts"][name]["physical"]["Q"][:3],
            same_source_replacement_not_extra_additive_contribution=True,
            numerical_QZZ_enclosure_available=False)
    graph=dict(base_archive=TRUE,base_sha256=pins[TRUE]["expected_sha256"],
        accepted_suffix_report=PARENT,accepted_suffix_sha256=PARENT_SHA256,
        base_graph_node_count=len(nodes),appended_nodes=g.appended,kernel=kernel,Q_recovery=recovery,
        old_Q4_and_624_node_source_suffix_unchanged=True)
    source_proof=dict(passed=True,ordinary_Z_orders=[0,1,2,3],
        original_t0_equals_minus_b_over_a=True,reflection_B_one_minus_phi_equals_minus_B_phi=True,
        C_phi_equals_B=True,C_0_equals_C_1_zero=True,C_y_is_fixed_phase_primitive_of_slow_B_y=True,
        source_owned_first_collar_makes_Cin_zero=True,only_original_rminus_mean_inlet_zero=True,
        common_roots_and_same_loop_recipe_after_inactive_first_collar_at_s1=True,
        first_second_C_and_Cy_seam_equal_as_free_Z_functions=True,
        ordinary_Z3_from_same_C5_functions_not_differentiated_caps=True,
        phase_total_y_derivative_equals_N=True,phase_integer_cuts_use_periodic_extension=True,
        actual_huge_N_phase_not_chosen_or_materialized=True,
        accepted_parity_receipt_consumed_not_replayed=True,
        accepted_endpoint_fingerprints_replayed=False)
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        source_family=parent["source_family"],source_pins=pins,source_AST=guards,source_proof=source_proof,
        new_free_phase_B_slow_By_seam=seam,
        graph_patch=graph,graph_patch_digest=digest(graph),ordinary_Z_orders=[0,1,2,3],
        exact_algebra=kernel["symbolic_controls"],actual_physical_mean_factor="Pstar/N^2",
        physical_Q_recovery_controls=recovery_controls(require),
        directed_physical_prefactor=physical_prefactor(parent,nodes[18],require),
        legacy_weighted_cap_diagnostic=legacy_diagnostic(parent,load(root/C3),load(root/C3_RANGE),loop,pins,require),
        actual_source_formula="delta_m=Pstar/N^2*(C-J_C-J_Cy), with nonreset joint convolutions",
        accepted_six_I1_source_digest=digest(a),same_source_internal_C_endpoint_canceled=True,
        complete_C_and_Cy_convolutions_retained=True,ordinary_mean_Z3_source_implemented=True,
        numeric_Cy_Z3_enclosure_available=False,whole_mean_Z3_numerical_enclosure_available=False,
        legacy_C3_caps_admitted_for_true_C5_source=False,numerical_QZZ_available=False,
        numerical_whole_micro_or_forward_Rp=False,global_five_moment_closure=False,
        old_six_and_all_error_caps_unchanged=True,accepted_ancestors_replayed=False,full_goal_complete=False)
