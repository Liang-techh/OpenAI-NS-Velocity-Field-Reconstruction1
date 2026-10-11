"""Executable original own-mean bindings and independent FTC/recovery checks."""
import ast
import copy
import json
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_leading_mean_epsilon_factor as current
import lei_ren_part1_paper_compliant_current_core_first_interface_check as interface
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported original own-mean factor admitted')


def executable_own_mean_bindings():
    specs=(
        ('compliant_core_physical_field','normalized_jets','average[index]','model/2+symmetric(c,avgerr)'),
        ('compliant_core_physical_field','normalized_jets','Vgrid[gridkey(i,k)]','base+self.epsilon*psi[gridkey(i,k)]'),
        ('compliant_core_physical_field','normalized_jets','Mgrid[gridkey(i,k)]','base+self.epsilon*average[gridkey(i,k)]'),
        ('compliant_actual_bridge_integrals','own_moment_feedback','targets','dict(H=phi0,M=V0,K=phi0*V0,A=square(V0),B=square(phi0)/2,C=square(phi0))'),
        ('compliant_actual_bridge_integrals','own_moment_feedback','baseline[name]','initial[name]*decay+targets[name]*(1-decay)'),
        ('compliant_current_original_n1_inner_operator','first_rows','mr','(w0-m0)/rho'),
        ('compliant_current_original_n1_inner_operator','first_rows','mrr','(wr-mr*2)/rho'))
    count=0
    for module,method,target,expression in specs:
        fn=interface.function(module,method);wanted=ast.dump(ast.parse(expression,mode='eval').body)
        wanted_target=ast.dump(ast.parse(target+' = None').body[0].targets[0])
        matches=[node for node in ast.walk(fn) if isinstance(node,ast.Assign)
                 and any(ast.dump(t)==wanted_target for t in node.targets) and ast.dump(node.value)==wanted]
        assert len(matches)==1,(module,method,target);count+=1
    assert current.domain.bridge.RATES['M']==1
    fn=interface.function('compliant_core_physical_field','normalized_jets')
    avg_assignments=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='avgerr' for t in n.targets)]
    assert {ast.unparse(n.value) for n in avg_assignments}=={'err / 2','err / (i + 1)'}
    # The inlet's own M atom is bound to the original scaled integral, not
    # a same-value comparison box. Its earlier accepted source hash remains.
    interface.binding('compliant_comparison_point_integrals','inlet','moments',
                      "{name:IntervalTaylor(c,values) for name,values in atoms['actual_core_atom_axial_coefficients'].items()}")
    return dict(executable_original_own_mean_source_bindings=count,
                original_mean_core_errors_and_actual_rate_one_verified=True,
                original_own_core_integral_inlet_bound=True)


def exact_ftc_recovery():
    r,z,eps,delta=s.symbols('rho Z epsilon delta',nonzero=True)
    u=s.Function('uaxis')(z);psi=s.Function('Psi')(r,z);inlet=s.Function('meanPsi4')(z)
    J=s.Function('I_Psi')(r,z);mean=(4*inlet+J)/r
    replacements={s.diff(J,r):psi,s.diff(J,r,2):s.diff(psi,r)}
    first=s.diff(mean,r).subs(replacements,simultaneous=True)
    second=s.diff(mean,r,2).subs(replacements,simultaneous=True)
    assert s.cancel(first-(psi-mean)/r)==0
    assert s.cancel(second-(s.diff(psi,r)-2*first)/r)==0
    assert s.cancel((r*mean).subs(r,4).subs(J.subs(r,4),0)-4*inlet)==0
    U=u+eps*psi;M=u+eps*mean;d=1-z*z;L=1-delta*z*z
    Q=(2*z*U-(1-delta)*z*M-d*s.diff(M,z))/L
    qa=((1+delta)*z*u-d*s.diff(u,z))/L
    qp=(2*z*psi-(1-delta)*z*mean-d*s.diff(mean,z))/L
    count=3
    for i in range(3):
        assert s.cancel(s.diff(Q,r,i)-s.diff(qa+eps*qp,r,i))==0;count+=1
    # Production factored-leading recovery is executed independently on
    # arbitrary smooth Psi and meanPsi functions for all three radial rows.
    fn=next(n for n in ast.parse((current.HERE/Path(current.__file__).name).read_text(encoding='utf8')).body
            if isinstance(n,ast.FunctionDef) and n.name=='factored_leading_rows')
    from types import SimpleNamespace
    algebra=SimpleNamespace(lift=lambda v,*args:v,shift=lambda v,p:v*eps**(-p[3]))
    env=dict(op=SimpleNamespace(dz=lambda v:s.diff(v,z),truncate=lambda v,n:v),ORDERS=current.ORDERS)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(fn)],type_ignores=[])),
                 '<live leading epsilon radial recovery>','exec'),env)
    rows=dict(axis=dict(z=z,u=u,d=d,L=L),Psi=[s.diff(psi,r,i) for i in range(3)],
              meanPsi=[s.diff(mean,r,i) for i in range(3)])
    got=env['factored_leading_rows'](algebra,rows,SimpleNamespace(delta=delta))
    for i in range(3):
        assert s.cancel(got['Uz'][i]-s.diff(U,r,i))==0
        assert s.cancel(got['mean'][i]-s.diff(M,r,i))==0
        assert s.cancel(got['Q'][i]-s.diff(Q,r,i))==0;count+=3
    return dict(independent_exact_own_FTC_and_leading_radial_recovery_identities=count,
                production_Q_decomposition_replayed_on_arbitrary_source_functions=True,
                leading_Q_axis_term_has_correct_one_plus_delta_coefficient=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**executable_own_mean_bindings(),**exact_ftc_recovery())
    axis=owner.report(packets['core_axis']);exit=owner.report(packets['selected_actual_inner_exit'])
    core=owner.report(packets['core_exit']);inlet=owner.report(packets['actual_first_inlet'])
    assert core['leading_Psi_rho_rows']==inlet['leading_Psi_rho_rows']
    assert core['leading_own_meanPsi_rho_rows']==inlet['leading_own_meanPsi_rho_rows']
    assert not axis['leading_Psi_rho_rows'][0]['terms'] and not axis['leading_own_meanPsi_rho_rows'][0]['terms']
    for packet in packets.values():
        v=owner.report(packet);entry=owner._packets[id(packet)][4]
        assert all(not v[k] for k in current.OPEN)
        assert v['actual_own_mean_history_not_comparison_mean'] and v['independent_raw_M_minus_uaxis_subtraction_not_used']
        assert v['component_weighted_n1_tail_still_open'] and v['n1_axial_pressure_source_integrals_still_open']
        assert [p.order for p in entry[1]['Psi']]==[5,4,3]
        assert [p.order for p in entry[1]['meanPsi']]==[5,4,3]
        assert [p.order for p in entry[2]['Q']]==[4,3,2]
        for label in ('Uz','mean','Q'):
            assert all(p[3]==-1 for row in entry[2][label][1:] for p in row.terms)
        for row in entry[2]['Qpsi']:
            assert all(p[3]==0 and p[1]==0 and p[2] in (0,1) for p in row.terms)
        assert all(p[0]>=0 for row in entry[4]['meanPsi'] for p in row.terms)
        assert all(p[0]>=-1 for row in entry[4]['Psi'] for p in row.terms)
        eta=owner.bound['common_tubes']['outer_Z_radius'];proof=v['uniform_complex_bound_proof']
        assert proof['original_outer_eta']==eta and proof['Psi_and_meanPsi_bound_tube_radius']==eta/2
        assert proof['Qpsi_and_Q_bound_tube_radius']==eta/4 and proof['point_increment_Taylor_Cauchy_radius']==eta/4
        assert all(f.order==2 for f in entry[3].values())
    assert len(exit['genuine_n1_forcing_from_factored_leading_rows'])==3
    checks['original_axis_and_identical_core_actual_first_inlet_factors']=True
    checks['actual_own_mean_and_Q_grade_five_four_three_and_four_three_two_available']=True
    checks['radial_leading_Uz_mean_Q_derivatives_retain_epsilon_pressure_sectors']=True
    checks['uniform_original_tube_and_actual_inverse_width_orders_checked']=True
    checks['genuine_three_n1_forcings_execute_from_exposed_leading_rows']=True
    checks['unsupported_endpoint_or_foreign_packet_rejected']=all((
        rejected(lambda:owner.evaluate('.371',c.mpf([0,1]))),
        rejected(lambda:owner.report(current.OriginalLeadingMeanEpsilonPacket()))))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])];old=entry[1]['component_weighted_n1_tail_still_open']
    entry[1]['component_weighted_n1_tail_still_open']=False
    try:checks['false_full_component_tail_scope_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['component_weighted_n1_tail_still_open']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf8',newline='\n')
    print('PASS_LEADING_OWN_MEAN_EPSILON_FACTOR',checks['independent_exact_own_FTC_and_leading_radial_recovery_identities'],flush=True)
    return result
