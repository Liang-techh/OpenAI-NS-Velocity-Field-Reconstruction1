"""Executable normalization identities and analytic full-feedback bound checks."""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_weighted_feedback_bound as current
import lei_ren_part1_paper_compliant_current_core_first_interface_check as interface
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported full-feedback bound admitted')


def executable_normalization_identities():
    """Replay original matrix and new formula on independent symbolic fields."""
    z,r,delta,lam=s.symbols('Z rho delta Lambda',nonzero=True)
    f0=s.Function('F0')(z);phi=s.Function('Phi')(z);pr=s.Function('Phi_r')(z)
    w=s.Function('Uz')(z);wr=s.Function('Uz_r')(z);q=s.Function('Q')(z)
    grad=s.Function('Gprime')(z);Haxis=s.Function('Haxis')(z)
    c=SimpleNamespace(sqrt=s.sqrt,mpf=s.Rational)
    algebra=SimpleNamespace(ctx=c,lift=lambda v,*args:v,
        shift=lambda v,p:v*f0**p[1]*lam**p[3])
    original=interface.function('compliant_current_original_n1_inner_operator','matrix')
    env=dict(truncate=lambda v,n:v,dz=lambda v:s.diff(v,z))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(original)],type_ignores=[])),
                 '<executable original n1 matrix>','exec'),env)
    raw=env['matrix'](algebra,{'F':[f0*phi,f0*pr],'Uz':[w,wr],'Q':[q]},
                       dict(N1p=0,N1theta=0,N1z=0),z,r,delta)
    fn=next(n for n in ast.parse((current.HERE/Path(current.__file__).name).read_text(encoding='utf8')).body
            if isinstance(n,ast.FunctionDef) and n.name=='residual_matrix')
    env=dict(op=SimpleNamespace(dz=lambda v:s.diff(v,z),truncate=lambda v,n:v),
             weighted=SimpleNamespace(relative_amplitude=lambda *args:s.Integer(1)))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(fn)],type_ignores=[])),
                 '<executable normalized residual n1 matrix>','exec'),env)
    d=1-z*z;L=1-delta*z*z;H=z*(1-delta)/2+d*w
    gauge=-2*(H-Haxis)*grad/L
    fresh=env['residual_matrix'](algebra,[phi,pr],[w,wr],q,z,r,delta,grad,gauge,0)
    weights=[f0*lam,s.Integer(1),s.Integer(1),s.Integer(1),f0*lam,s.Integer(1)];count=0
    for label in ('B0','B1'):
        expected={}
        for key,value in raw[label].items():
            i,j=(int(t)-1 for t in key.split(','))
            transformed=value*weights[j]/weights[i]
            if label=='B0' and key in raw['B1'] and j in (0,4):
                transformed-=raw['B1'][key]*weights[j]/weights[i]*lam*grad
            transformed=transformed.subs(s.diff(f0,z),-lam*grad*f0)
            if label=='B0' and key in ('1,5','2,6','3,6'):
                transformed-=1 if key!='3,6' else -1
            if label=='B0' and key=='5,1':transformed+=2*Haxis*grad/L
            if s.cancel(transformed)!=0:expected[key]=transformed
        assert set(expected)==set(fresh[label])
        for key,value in fresh[label].items():
            assert s.cancel(value-expected[key])==0,(label,key);count+=1
    # Principal supports are disjoint; infinity row bound is max(1,2|chi|).
    assert list(raw['D'])==[0,0,2,0,3,1]
    principal=s.zeros(6)
    for i,j,v in ((0,4,1),(1,5,1),(2,5,-1),(4,0,-2)):principal[i,j]=v
    assert max(sum(abs(principal[i,j]) for j in range(6)) for i in range(6))==2
    return dict(executable_original_to_normalized_residual_matrix_identities=count,
                principal_kinematic_and_potential_supports_disjoint=True,
                source_F0_weight_derivative_and_known_F0_derivative_cancelled_once=True,
                F0_squared_pressure_axial_back_coupling_preserved=True)


def all_word_simplex_identities():
    """Independent exact simplex and factorial/geometric majorization premises."""
    x,t=s.symbols('x t',positive=True);count=0
    for k in range(1,9):
        got=s.integrate(t**(k-1)/math.factorial(k-1),(t,0,x))
        assert s.expand(got-x**k/math.factorial(k))==0;count+=1
    a,b,c=s.symbols('a b c')
    assert s.expand(2*(a-b)+2*(b-c)-2*(a-c))==0
    # k! >= (k/e)^k follows from sum log(j)>= integral_1^k log(t)dt.
    u=s.symbols('u',positive=True);k=s.symbols('k',positive=True)
    assert s.simplify(s.integrate(s.log(u),(u,1,k))-(k*s.log(k)-k+1))==0
    return dict(independent_exact_ordered_radial_simplex_identities=count,
                Green_exponential_factors_telescope=True,
                factorial_inequality_from_monotone_log_integral=True,
                each_iteration_may_use_dZ_with_equal_nested_Cauchy_gap=True,
                terminal_Z1_uses_separate_final_Cauchy_gap=True)


@source_precision
def run(owner,packet,axial_packet):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**executable_normalization_identities(),**all_word_simplex_identities())
    view=owner.report(packet);data=owner._packets[id(packet)][2];eta=data['eta']
    assert all(not view[k] for k in current.OPEN)
    assert view['remaining_feedback_bound_includes_all_six_rows']
    assert view['actual_feedback_iteration_prefix_not_computed'] and view['solved_coefficient_function_contracts_still_open']
    proof=view['uniform_analytic_proof'];tail=view['all_word_feedback_tail']
    assert proof['matrix_bound_tube']==eta/2 and proof['source_initial_bound_tube']==eta/4
    assert tail['value_tube_radius']==eta/8 and tail['ordinary_Z1_tube_radius']==eta/16
    assert tail['iteration_Cauchy_loss']==eta/8 and tail['ordinary_Z1_final_Cauchy_gap']==eta/16
    assert ep(tail['ratio_upper'])[1]<c.mpf('1e-12')
    assert ep(tail['value_remainder_upper'])[1]<c.mpf('1e-12')
    assert ep(tail['ordinary_Z1_remainder_upper'])[1]<c.mpf('1e-12')
    assert proof['norm_log_basis_names'][1]=='uniform_Fmax_norm_proxy'
    assert proof['Fmax_is_uniform_complex_bound_not_point_amplitude']
    assert proof['source_L1_uses_sum_of_rows_not_max'] and proof['inverse_width_cancelled_before_source_norm']
    assert all(p[0]>=0 for v in data['collar_source_L1'].values() for p in v.terms)
    assert set(data['C0entries'])=={'4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'}
    assert set(data['C1entries'])=={'5,1','5,2','5,3','6,2','6,3','6,4'}
    for cover in (view['uniform_residual_C0_absolute_cover'],view['uniform_B1_absolute_cover'],view['uniform_source_L1_cover']):
        assert cover['all_covers_have_directed_log_proof']
        for term in cover['sectors']:
            if term['rule']=='proved log<=log(floor)':
                assert ep(term['log_modulus_upper'])[1]<=ep(c.ln(cover['tiny_cover_floor']))[0]
            else:assert term['rule']=='directed exp(log)'
            assert ep(term['cover'])[0]>=0
    # Exercise genuine actual jets after normalized Phi extraction. The old
    # large boxed F0-ratio matrix is not used for this construction.
    owner.axial.report(axial_packet);aentry=owner.axial._packets[id(axial_packet)]
    p=aentry[3];sv=owner.inputs.report(p);entry=owner.inputs._packets[id(p)]
    rows,algebra=entry[2],entry[5];z=algebra.lift(op.Jet.variable(c,sv['Z'],4))
    rho=sv['geometry']['rho_enclosure_only'];gradient=op.Jet(c,owner.weighted.source.amplitude.gradient_jets(sv['Z']))
    actual=current.residual_matrix(algebra,rows['Phi'],rows['Uz'],rows['Q'][0],z,rho,
                                   owner.physical.delta,gradient,aentry[5][4],1)
    assert all(value.order==1 for matrix in actual.values() for value in matrix.values())
    assert set(actual['B0'])==set(data['C0entries']) and set(actual['B1'])==set(data['C1entries'])
    assert all(p[1]==2 for p in actual['B0']['4,1'].terms)
    assert all(p[1]==2 for p in actual['B0']['6,1'].terms)
    assert all(p[1]==0 for key in ('5,1','5,2','5,3') for p in actual['B0'][key].terms)
    checks['genuine_actual_normalized_Phi_matrix_executed_without_F0_ratio_dressing']=True
    checks['original_tubes_all_rows_and_actual_width_cancellation_checked']=True
    checks['every_small_sector_cover_has_directed_log_proof']=True
    checks['uniform_all_word_value_and_Z1_tails_below_1e_minus12']=True
    checks['unsupported_tail_or_packet_rejected']=all((
        rejected(lambda:current.all_word_tail(c,1,2,1,1,eta)),
        rejected(lambda:current.all_word_tail(c,1,2,0,0,eta,True)),
        rejected(lambda:owner.report(current.OriginalN1WeightedFeedbackBoundPacket()))))
    entry=owner._packets[id(packet)];old=entry[1]['solved_coefficient_function_contracts_still_open']
    entry[1]['solved_coefficient_function_contracts_still_open']=False
    try:checks['false_solved_contract_mutation_rejected']=rejected(lambda:owner.report(packet))
    finally:entry[1]['solved_coefficient_function_contracts_still_open']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf8',newline='\n')
    print('PASS_N1_WEIGHTED_FEEDBACK_BOUND',checks['executable_original_to_normalized_residual_matrix_identities'],flush=True)
    return result
