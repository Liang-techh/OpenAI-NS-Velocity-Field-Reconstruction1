"""Check actual signed C3 ownership, inverse calculus and response identity."""
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_complete_C3_defect_repair_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions_check as prior


def actual_source(owner):
    g=owner.g;a=current.Algebra(g,3)
    assert owner.inner_endpoint is owner.cutoff.endpoint
    assert owner.lower_pressure is owner.lower.pressure
    assert owner.preRp_pressure is owner.quiet.Rp_endpoint['pressure']
    assert owner.pulse_pressure is owner.pulse.pressure
    assert owner.linear is owner.full.linear
    assert all(len(current.rows(q))==4 for q in owner.moments.values())
    assert all(len(current.rows(q))==4 for q in owner.coefficients.values())
    assert all(len(current.rows(owner.repair_velocity[k]))==4 for k in ('Uz1','Utheta1','F1'))
    assert owner.moments==current.compose_full_moments(g,owner.cutoff.endpoint,
        owner.lower.pressure,owner.quiet.Rp_endpoint['pressure'],owner.pulse.pressure)
    for k in ('m11','m12','m14'):assert owner.moments[k]==a.trim(owner.cutoff.endpoint[k])
    inner=current.parent.linear.read(current.original.cutoff.RECEIPT)
    assert inner['all_passed'] and inner['checks']['own_history']['passed'] and inner['checks']['source']['passed']
    for q in owner.cutoff.endpoint.values():
        for ref in current.rows(q):
            n=g.nodes[ref.node]
            assert n['operation']=='function_substitution'
            assert n['variable']==owner.cutoff.inner.coordinates['first'].node
            assert n['value']==owner.cutoff.parameters['sc_quarter'].node
    linear=current.parent.linear.read(current.original.linear.RECEIPT)
    assert linear['determinant_identities']['axial_exact_difference_identity']
    assert linear['determinant_identities']['angular_exact_Vandermonde_identity']
    assert linear['actual_mu_positive_and_retained'] and linear['original_amplitude_positive_on_all_four_Z_cells']
    assert current.encode(owner.inverse_bindings['Astar'])==owner.linear.background['Astar_ordinary_Z0_Z3']
    assert owner.reserved_background==owner.full.reserved_background
    assert owner.beyond_Rv_support==owner.full.beyond_Rv_support
    assert all(all(ref==g.zero for ref in current.rows(q)) for q in owner.segments['beyond_Rv'].values())
    assert not owner.definition['quantitative_C3_bounds_transferred']
    assert not owner.definition['repaired_own_M_V_P_installed']
    return dict(passed=True,actual_retained_inner_C3_endpoint_rows=20,
        existing_inner_C3_receipt_reused_not_rebuilt=True,
        lower_preRp_pulse_same_graph_components_retained=True,
        original_Rin_endpoint_sc_hb_and_unchanged_N1p_history_retained=True,
        selected_operator_Astar_Z0_Z3_nonzero_mu_support_and_determinant_guards=True,
        full_physical_field_join_and_quantitative_C3_bounds_not_claimed=True)


def actual_signed_partition(owner):
    g=owner.g;bindings={};expected={};symbols={}
    for segment,fields in owner.segments.items():
        if segment=='beyond_Rv':continue
        for key,q in fields.items():
            symbols[segment,key]=[]
            for n,ref in enumerate(current.rows(q)):
                value=sy.Symbol(f'{segment}_{key}_{n}',real=True)
                assert ref.node not in bindings or bindings[ref.node]==value
                bindings[ref.node]=value;symbols[segment,key].append(value)
    cache={};count=0
    for key,q in owner.moments.items():
        for n,ref in enumerate(current.rows(q)):
            wanted=symbols['inner',key][n]
            if key in ('m13','m15'):wanted+=sum(symbols[s,key][n] for s in ('lower','preRp','pulse'))
            assert sy.expand(prior.symbolic(g,ref,bindings,cache)-wanted)==0,(key,n)
            count+=1
    return dict(passed=True,actual_saved_signed_partition_rows=count,
        all_exterior_pressure_components_counted_once=True,
        velocity_moments_unchanged_after_cutoff=True,m15_no_second_negation=True)


def actual_inverse_and_response(owner):
    """Arbitrary genuine functions, actual saved inverse and C3 graph rows."""
    g=owner.g;W=owner.linear.functions;b=owner.linear.background
    Z=sy.Symbol('Z',real=True);mu=sy.Symbol('original_mu',positive=True)
    Rc2=sy.Symbol('two_Rc',positive=True);amp=sy.Function('Astar')(Z)
    A=sy.Matrix(2,2,lambda i,j:sy.Symbol(f'axial_{i}_{j}',real=True))
    T=sy.Matrix(3,3,lambda i,j:sy.Symbol(f'angular_{i}_{j}',real=True))
    bindings={W['mu'].node:mu,W['log2Rc'].node:sy.log(Rc2),
        W['positive_axial_divided_determinant'].node:-A.det(),
        W['positive_angular_determinant'].node:T.det()}
    for n,node in enumerate(b['Astar_ordinary_Z0_Z3']):bindings[node]=sy.diff(amp,Z,n)
    for i,row in enumerate(W['axial_rows']):
        for j,q in enumerate(row):bindings[q.node]=A[i,j]
    for i,row in enumerate(W['angular_rows']):
        for j,q in enumerate(row):bindings[q.node]=T[i,j]
    for key,q in owner.moments.items():
        body=sy.Function(key)(Z)
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.diff(body,Z,n)
    cache={};calculus=[];response=[]
    for key,q in owner.coefficients.items():
        value=prior.symbolic(g,q.value,bindings,cache)
        for n,ref in enumerate(current.rows(q)[1:],1):
            actual=prior.symbolic(g,ref,bindings,cache)
            assert sy.cancel(sy.powsimp(actual-sy.diff(value,Z,n),force=True))==0,(key,n)
            calculus.append(f'{key}_Z{n}')
    for key,q in owner.linear_response_residuals.items():
        for n,ref in enumerate(current.rows(q)):
            value=prior.symbolic(g,ref,bindings,cache)
            assert sy.cancel(sy.powsimp(value,force=True))==0,(key,n)
            response.append(f'{key}_Z{n}')
    radius=sy.Symbol('R',positive=True);bindings[owner.repair_velocity['physical_radius'].node]=radius
    for j,q in enumerate(W['physical_bumps']):bindings[q.node]=sy.Symbol(f'physical_bump_{j}',real=True)
    cache={};repair=[]
    for key in ('Uz1','Utheta1','F1'):
        q=owner.repair_velocity[key];value=prior.symbolic(g,q.value,bindings,cache)
        for n,ref in enumerate(current.rows(q)[1:],1):
            assert sy.cancel(prior.symbolic(g,ref,bindings,cache)-sy.diff(value,Z,n))==0,(key,n)
            repair.append(f'{key}_Z{n}')
    for n in range(4):
        F=prior.symbolic(g,current.rows(owner.repair_velocity['F1'])[n],bindings,cache)
        theta=prior.symbolic(g,current.rows(owner.repair_velocity['Utheta1'])[n],bindings,cache)
        assert sy.cancel(F-theta/sy.sqrt(2*radius))==0,n
    return dict(passed=True,actual_coefficient_ordinary_derivative_rows=calculus,
        arbitrary_amplitude_and_defect_functions_not_samples=True,
        actual_physical_response_C3_zero_rows=response,
        accepted_exact_determinant_and_physical_bump_response_identities_consumed=True,
        original_finite_mu_divided_difference_and_Astar_third_derivatives_retained=True,
        actual_repair_ordinary_derivative_rows=repair,
        physical_F_equals_Utheta_over_sqrt_2R_C3=True,
        algebraic_response_only_not_terminal_own_field_identity=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['complete_C3_defect_repair_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_source',actual_source),('actual_signed_partition',actual_signed_partition),
                    ('actual_inverse_and_response',actual_inverse_and_response)):
        checks[name]=fn(owner);print('COMPLETE_C3_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('COMPLETE_SIGNED_C3_DEFECT_AND_RESERVED_REPAIR_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses checked source owners')
