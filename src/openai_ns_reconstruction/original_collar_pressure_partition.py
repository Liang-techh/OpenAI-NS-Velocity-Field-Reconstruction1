"""Exact same-leading-source regular/singular collar pressure partition.

The existing forward regular frame omits only chi_rho.  This module derives
its missing physical pressure density and connects it to the accepted signed
A-A0 primitive.  The raw cut n1 fields remain separate input functions;
neither their finite inlet approximation nor a positive error cap is used as
an actual source value.  No inherited producer or quadrature is executed.
"""
from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp

from .original_finite_N_signed_source import load,sha,digest,require
from .original_signed_density_expression import SignedDensityExpressionGraph

PARENT='lei_ren_part1_paper_compliant_current_original_n1_forward_enclosure.json'
PARENT_SHA256='2906a9ca82c995ce8dd2419041835e5804550cf4af4d28292ad450d5919a7e54'
LOWER='lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz'
LOWER_SOURCE='lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py'
FORWARD='src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'


def _function(root,path,name):
    tree=ast.parse((root/path).read_text(encoding='utf-8'))
    nodes=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name]
    require(len(nodes)==1,'Unique actual source function required: '+name)
    node=nodes[0]
    return ast.unparse(node),dict(path=path,function=name,
        source_sha256=sha(root/path),
        AST_sha256=hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest(),
        source_executed=False)


def exact_algebra():
    """Fresh symbolic proof, with ordinary derivatives rather than Taylor jets."""
    r,h,e,P2,z,de,sigmas,chi,chir,A,Ar,ratio,pr,prr,phi=sp.symbols(
        'rho hb epsilon P2 Z delta sigma_s chi chi_rho A Ar ratio Psi_r Psi_rr Phi',nonzero=True)
    L=1-de*z*z
    # Differentiate ratio*chi*Psi_core,rho using ratio_rho=(chi-1)*A*ratio.
    psi_rr_full=ratio*((chi-1)*A*chi*pr+chir*pr+chi*prr)
    psi_rr_regular=ratio*(chi*(chi-1)*A*pr+chi*prr)
    phi_rr_regular=ratio*(chi*(chi-1)*A*(phi*A)+chi*(phi*(Ar+A*A)))
    checks={
        'Psi_rr_missing_only_ratio_chi_rho_Psi_core_r':sp.expand(
            psi_rr_full-psi_rr_regular-ratio*chir*pr)==0,
        'Phi_regular_row_equivalent':sp.expand(phi_rr_regular-ratio*phi*(chi*chi*A*A+chi*Ar))==0,
    }
    qrr_missing=2*z*e*ratio*pr*chir/L
    raw_missing=sp.cancel((r/e*qrr_missing).subs(chir,(h-1)*sigmas/(h*r)))
    expected_raw=2*z*(h-1)*ratio*pr*sigmas/(h*L)
    Aexact=2*r*r*z*(1-h)*ratio*pr/L
    checks['raw_N1p_missing_has_no_extra_rho']=sp.cancel(raw_missing-expected_raw)==0
    checks['weighted_physical_singular_density']=sp.cancel(-e*e*h*r*r*raw_missing-e*e*Aexact*sigmas)==0
    checks['normalized_P2_removed_once']=sp.cancel(-e*e*P2*h*r*r*(raw_missing/P2)-e*e*Aexact*sigmas)==0
    # Q,Q_rho,U and their slow derivatives are identical in full/regular
    # rows; the source is affine in Q_rhorho with coefficient rho/epsilon.
    q,qr,qrr,u,qz=sp.symbols('Q Q_r Q_rr U Q_Z')
    flow=(q+z*(1-de)*qz/2+r*qr)/L+q*(q/2+r*qr)+u*(-2*z*q+(1-z*z)*qz-2*z*r*qr)/L
    full=-flow/2+(2*qr+r*(qrr+qrr_missing))/e
    regular=-flow/2+(2*qr+r*qrr)/e
    checks['full_source_affine_difference']=sp.cancel(full-regular-r*qrr_missing/e)==0
    # Test ordinary Z1 on arbitrary same-source ratio/Psi_r functions, with
    # rho,hb,epsilon,delta and sigma_s held fixed as required by the geometry.
    R=sp.Function('ratio')(z);B=sp.Function('Psi_core_rho')(z)
    Aj=Aexact.subs({ratio:R,pr:B});Nj=expected_raw.subs({ratio:R,pr:B})
    checks['ordinary_Z1_physical_density']=sp.simplify(sp.diff(-e*e*h*r*r*Nj-e*e*Aj*sigmas,z))==0
    A0=32*z*B/L
    explicit_A0_Z=32*((B+z*sp.diff(B,z))/L+2*de*z*z*B/L**2)
    checks['A0_ordinary_Z1']=sp.simplify(sp.diff(A0,z)-explicit_A0_Z)==0
    # Native regular source includes the cut F1 coupling.  Its source field
    # is retained symbolically; this equality is scalar algebra, not a claim
    # that the two independent raw-cut producers have been identified.
    F0,fc,pc,Nreg=sp.symbols('F0 f_cut Phi_collar Nreg')
    S=e*e*F0*F0
    forward=-e*e*P2*h*r*r*(Nreg/P2+2*S*pc*fc/(e**3*P2))
    raw=-e*e*h*r*r*(Nreg+2*F0*F0*pc*fc/e)
    checks['cut_F_coupling_physical_scales']=sp.cancel(forward-raw)==0
    # Differentiate the signed split by FTC.  Its origin vanishes because
    # the original sigma(0)=0; no incoming pressure is reset here.
    t,x=sp.symbols('t s',real=True)
    af=sp.Function('A');a0=sp.Function('A0');sf=sp.Function('sigma')
    split=a0(z)*sf(x)+sp.Integral((af(t,z)-a0(z))*sp.diff(sf(t),t),(t,0,x))
    checks['signed_split_FTC_C0']=sp.simplify(sp.diff(split,x)-af(x,z)*sp.diff(sf(x),x))==0
    checks['signed_split_FTC_ordinary_Z1']=sp.simplify(sp.diff(sp.diff(split,x),z)-sp.diff(af(x,z),z)*sp.diff(sf(x),x))==0
    b1,b2,r1,r2=sp.symbols('boundary_at_q1 boundary_at_q2 remainder_at_q1 remainder_at_q2')
    checks['same_source_internal_endpoint_telescope']=sp.expand(
        (b1+r1)+(b2-b1+r2-r1)-(b2+r2))==0
    require(all(checks.values()),'An exact collar pressure identity failed')
    return dict(passed=True,checks=checks,
        chi_rho='(hb-1)*sigma_s/(hb*rho)',
        raw_N1p_singular=str(expected_raw),
        A=str(Aexact),
        raw_N1p='-Sflow/2+(2*Q_rho+rho*Q_rhorho)/epsilon',
        normalized_Np='raw_N1p/Pstar^2',
        weighted_pressure_density='-epsilon^2*hb*rho^2*raw_N1p',
        weighted_singular_density='epsilon^2*A*sigma_s',
        ordinary_Z_orders=[0,1],geometry_parameters_Z_independent=True,
        old_source_functions_executed=False)


def build(root,parent,leading):
    root=Path(root);base=root/'experiments/root_st073'
    require(sha(base/PARENT)==PARENT_SHA256 and digest(parent)==digest(load(base/PARENT)),
        'Actual accepted forward source parent required')
    require(leading.get('passed') is True and leading.get('parent_report_sha256')==PARENT_SHA256,
        'Same-source leading collar function receipt required')
    require(leading.get('checks') and all(v is True for v in leading['checks'].values()) and
        leading['source_family']==parent['source_family'] and
        leading['actual_leading_collar_function_identity']['new_smooth_object_equals_canonical_lower_collar'] is True and
        leading['actual_leading_collar_function_identity']['accepted_native_inlet_same_function'] is True and
        leading['actual_leading_collar_function_identity']['accepted_true_pressure_same_function'] is True,
        'Actual leading collar definition/graph identity required')
    require(leading['actual_leading_collar_function_identity']['old_forward_producer_equals_canonical'] is False and
        leading['actual_leading_collar_function_identity']['bridge_D_over_R_to_ratio_Psi_ODE_identity'] is False,
        'Native bridge scope changed; independently review its pressure attachment')
    require(sha(base/LOWER)==parent['input_hashes'][LOWER] and
        sha(base/LOWER_SOURCE)==parent['input_hashes'][LOWER_SOURCE] and
        sha(root/FORWARD)==parent['implementation_hashes'][FORWARD], 'Changed collar pressure source')
    fields_text,fields_pin=_function(root,'experiments/root_st073/'+LOWER_SOURCE,'collar_fields')
    pressure_text,pressure_pin=_function(root,'experiments/root_st073/'+LOWER_SOURCE,'direct_pressure_functions')
    frame_text,frame_pin=_function(root,FORWARD,'leading_regular_frame')
    regular_text,regular_pin=_function(root,FORWARD,'normalized_pressure_source')
    compact=lambda text:''.join(text.split())
    fields_text,pressure_text,frame_text,regular_text=map(compact,
        (fields_text,pressure_text,frame_text,regular_text))
    guards={
        'canonical_chi_rho':"g.mul(g.sub(hb,g.one),sigma_prime(s))" in fields_text,
        'canonical_own_mean_second_row':"a.add(psi[1],a.neg(a.scale(mean1,2)))" in fields_text,
        'canonical_pressure_Qyy_plus_Qy':"a.add(Qyy,Qy)" in pressure_text,
        'canonical_pressure_divides_physical_R':"g.quotient(g.one,fields['radius']" in pressure_text,
        'native_regular_chi':"chi=1-sigma+sigma*hb" in frame_text,
        'native_regular_second_Psi':"chi*(chi-1)*A*comparison['Psi'][1]+chi*comparison['Psi'][2]" in frame_text,
        'native_regular_pressure':"(qr*2+qrr*r)/(e*P2)" in regular_text,
    }
    require(all(guards.values()),'Actual pressure frame source formula changed')
    lower=load(base/LOWER);bindings=lower['lower_collar_bindings'];nodes=lower['exact_graph_nodes']
    require(bindings['geometry_certificate']['sc_and_hb_Z_independent'] is True and
        bindings['pressure_bindings']['N1p']==[36583,36584] and
        bindings['fields']['weighted_physical_measure']==36494,
        'Actual pressure/ordinary-Z/weighted measure bindings changed')
    primitive=parent['canonical_collar_signed_primitive']
    require(primitive['passed'] is True and primitive['source_binding_passed'] is True and
        primitive['local_nodes']==99 and primitive['new_exact_controls']['passed'] is True and
        primitive['full_singular_sum']['exact_split'] is True and
        primitive['signed_primitive']['physical_epsilon_squared_factor_count']==1 and
        primitive['graph']['source_graph_nodes_untouched'] is True,
        'Accepted signed A-A0 source primitive required, not its norm cap')
    require(all(primitive['radial_checks'].values()) and primitive['ordinary_C1']['exact_ordinary_C1_algebra'] is True,
        'Actual radial-order1 and ordinary-Z1 A/A0 sources required')
    algebra=exact_algebra()
    g=SignedDensityExpressionGraph()
    ref=lambda i:g.reference(LOWER,parent['input_hashes'][LOWER],i,source_nodes=nodes)
    eps,h,r,sigmap=(ref(i) for i in (35565,8,35851,35953))
    Arefs=[]
    for j,i in enumerate(primitive['A']['C0_C1']):
        require(0<=i<len(primitive['exact_local_expression_nodes']),'Actual signed A node missing')
        Arefs.append(g.insert(dict(operation='saved_local_function_reference',report=PARENT,
            report_sha256=PARENT_SHA256,packet='canonical_collar_signed_primitive',
            packet_sha256=digest(primitive),local_node=i,ordinary_Z_order=j,
            exact_function_not_numerical_value=True)))
    rows=[]
    for j,(i,Aref) in enumerate(zip(bindings['pressure_bindings']['N1p'],Arefs)):
        full=ref(i)
        nsing=-Aref*sigmap/(h*r*r)
        nreg=full-nsing
        density=-eps*eps*h*r*r*full
        singular=eps*eps*Aref*sigmap
        regular=density-singular
        rows.append(dict(ordinary_Z_order=j,raw_full_N1p=full.node,raw_regular_N1p=nreg.node,
            raw_singular_N1p=nsing.node,physical_full_pressure_density=density.node,
            physical_regular_pressure_density=regular.node,physical_singular_pressure_density=singular.node))
    require(bindings['fields']['domain']==dict(lower=35824,upper=35825) and
        bindings['parameters']['sc_quarter']==35824 and bindings['parameters']['sc_half']==35825,
        'Actual q1/q2 phase endpoints changed')
    def at_saved(node,value,order):
        return g.insert(dict(operation='saved_local_function_evaluation',report=PARENT,
            report_sha256=PARENT_SHA256,packet='canonical_collar_signed_primitive',
            packet_sha256=digest(primitive),local_expression_node=node,
            packet_variable_node=primitive['signed_primitive']['variable_scope']['upper_s'],
            value_node=value.node,ordinary_Z_order=order,
            exact_function_not_numerical_value=True))
    endpoint_rows=[]
    for j,(bn,rn) in enumerate(zip(primitive['term4_A0_sigma']['term4_A0_sigma_C0_C1'],
            primitive['signed_primitive']['signed_primitive_C0_C1'])):
        b1,b2=at_saved(bn,ref(35824),j),at_saved(bn,ref(35825),j)
        r1,r2=at_saved(rn,ref(35824),j),at_saved(rn,ref(35825),j)
        first=b1+r1;lower_increment=(b2-b1)+(r2-r1)
        endpoint_rows.append(dict(ordinary_Z_order=j,first_q0_q1=first.node,
            lower_q1_q2_increment=lower_increment.node,expanded_join=(first+lower_increment).node,
            reduced_q0_q2=(b2+r2).node,boundary_q1_common_node=b1.node,
            signed_remainder_q1_common_node=r1.node,
            common_source_internal_pair_cancels=True))
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        leading_collar_binding_sha256=digest(leading),
        old_forward_leading_function_identity=leading['actual_leading_collar_function_identity'][
            'full_forward_leading_function_identity'],
        native_D_over_R_bridge_equivalence_claimed=False,
        checks=dict(actual_source_AST=all(guards.values()),actual_ordinary_Z_geometry=True,
            exact_pressure_partition=algebra['passed'],accepted_signed_primitive_consumed=True),
        source_AST=dict(canonical_fields=fields_pin,canonical_pressure=pressure_pin,
            forward_regular_frame=frame_pin,forward_regular_pressure=regular_pin),
        actual_source_guards=guards,exact_algebra=algebra,
        source_function_graph=g.packet(),pressure_density_rows=rows,
        signed_singular_endpoint_rows=endpoint_rows,
        local_same_source_singular_q1_pair_contracted=True,
        saved_canonical_inner_plus_lower_or_ledger_contracted=False,
        signed_primitive_binding=dict(packet_sha256=digest(primitive),
            full_singular_sum_C0_C1=primitive['full_singular_sum']['full_singular_sum_C0_C1'],
            A0_sigma_C0_C1=primitive['term4_A0_sigma']['term4_A0_sigma_C0_C1'],
            signed_A_minus_A0_C0_C1=primitive['signed_primitive']['signed_primitive_C0_C1'],
            consumed_once=True,original_origin_sigma0_zero=True,incoming_core_pressure_history_retained=True,
            term5_norm_cap_selected_as_value=False),
        full_pressure_partition_for_same_leading_fields=True,
        physical_Pstar_squared_canceled_once=True,ordinary_Z0_Z1=True,
        existing_signed_A_minus_A0_function_reused=True,
        numerical_signed_A_minus_A0_value=False,
        first_chart_raw_cut_response_identity=False,
        complete_forward_regular_plus_cut_coupling_identity=False,
        six_signed_source_equality=False,six_terms_removed=False,all_old_error_caps_retained=True,
        global_moment_closure=False,accepted_ancestor_or_quadrature_replayed=False,
        limitations=['This identifies the pressure partition on the same actual leading collar fields',
            'The raw cut n1 response still needs a separate function identity; its finite inlet approximation is not a field value',
            'The existing signed A-A0 primitive is a function definition, not the saved positive term5 cap or a numerical signed value'])
