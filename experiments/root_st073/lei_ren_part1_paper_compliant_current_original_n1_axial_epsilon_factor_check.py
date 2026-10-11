"""Executable original-drive bindings and independent epsilon/gauge identities."""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_axial_epsilon_factor as current
import lei_ren_part1_paper_compliant_current_core_first_interface_check as interface
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported axial epsilon packet admitted')


def live_drive_binding():
    """Bind executable three-scale drive and replay its actual k0 control."""
    count=0
    specs=(
        ('compliant_bridge_mixed_C4','evaluate','vcore',"[sum((core['Uz'][i]*(S2[k][i]*rho**i) for i in range(1,k+1)),core['Uz'][0]*0) for k in range(1,4)]"),
        ('compliant_bridge_mixed_C4','evaluate','chi',"[algebra.lift(1-sigma[0])+h*sigma[0]]+[-algebra.lift(sigma[k])+h*sigma[k] for k in range(1,4)]"),
        ('compliant_bridge_mixed_C4','evaluate','drive',"[rate_rows(directions['hydro'],scale,k)*R+algebra.shift(rate_rows(directions['pressure'],scale,k)*R,(0,1,0,0))+algebra.shift(rate_rows(directions['swirl'],scale*2,k)*R**2,(0,0,1,0)) for k in range(4)]"),
        ('compliant_bridge_mixed_C4','evaluate','controls',"bridge_controls(algebra,scale,chi,Dbar,drive,phi/comparison['phi'].truncate(5),barlog)"),
        ('compliant_actual_bridge_integrals','contributions','scales',"dict(hydro=c.mpf(0),pressure=self.pressure_log,swirl=self.swirl_log)"),
        ('compliant_actual_bridge_integrals','contributions','powerR',"dict(hydro=1,pressure=1,swirl=2)"),
        ('compliant_actual_bridge_integrals','contributions','row',"-(p['micros'][0]['quotient']*d['drive_'+name])*Rmicro**powerR[name]*mass"),
        ('compliant_current_original_n1_inner_operator','first_rows','ratio','phi/pc[0]'),
        ('compliant_current_original_n1_inner_operator','first_rows','wr','ratio*chi*wc[1]'),
        ('compliant_core_physical_field','axis_inputs','slope','-g/(2*L)'),
        ('compliant_core_physical_field','normalized_jets','model',"(source['slope'][k]*math.factorial(k))*(rho if i==0 else 1 if i==1 else 0)"))
    for module,method,target,expression in specs:
        fn=interface.function(module,method);wanted=ast.dump(ast.parse(expression,mode='eval').body)
        matches=[node for node in ast.walk(fn) if isinstance(node,ast.Assign)
                 and any(ast.unparse(t)==target for t in node.targets) and ast.dump(node.value)==wanted]
        if len(matches)!=1:raise ValueError('Original drive source assignment differs: '+module+'.'+method+': '+target)
        count+=1
    fn=interface.function('compliant_bridge_mixed_C4','bridge_controls')
    wanted=ast.parse('value+=chi[i]*quotient_rows[j]*drive[ell]*weight').body[0]
    assert any(ast.dump(node)==ast.dump(wanted) for node in ast.walk(fn))
    wanted=ast.parse('V.append(-scale*value)').body[0]
    assert any(ast.dump(node)==ast.dump(wanted) for node in ast.walk(fn))
    # Source theorem: all three physical drives equal -V_R in the common core.
    proof=interface.integrated_core_equations_proof();assert proof['passed']
    inlet=interface.defining_boundary_bindings();assert inlet['passed']
    # Replay production coordinate control on the theorem, not a copied formula.
    hb,rho,eps,chi,ratio,psir=s.symbols('hb rho epsilon chi ratio Psi_rho')
    env=dict(math=math,product_rows=lambda a,b,k:sum(math.comb(k,j)*a[j]*b[k-j] for j in range(k+1)),
             exponential_derivatives=lambda rows:[s.Integer(1),s.Integer(0),s.Integer(0),s.Integer(0)])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(fn)],type_ignores=[])),
                 '<live first axial control>','exec'),env)
    controls=env['bridge_controls'](SimpleNamespace(lift=lambda value:value),hb,[chi,0,0,0],
                                    [0]*4,[-rho*eps*psir,0,0,0],ratio,[0]*4)
    assert s.expand(controls['Uz_positive_order_coordinate_derivatives'][0]-hb*rho*eps*chi*ratio*psir)==0
    return dict(executable_original_source_assignment_bindings=count,
                live_three_scale_stress_free_drive_theorem_replayed=True,
                live_original_first_coordinate_control_epsilon_rho_Jacobian_identity=True,
                accepted_original_first_inlet_function_identity_bound=True,
                source_metadata_strings_not_used_as_drive_proof=True)


def independent_gauge_identity():
    z,delta,lam,F0,P,H,sigma,Psi=s.symbols('Z delta Lambda F0 P H sigma Psi')
    d=1-z*z;L=1-delta*z*z;gradient=L*H/(H*H+sigma*sigma)
    eps=1/lam;Ha=H+d*eps*Psi
    full=-2*Ha*gradient/L+2*H*H/(H*H+sigma*sigma)
    residual=-2*eps*d*Psi*H/(H*H+sigma*sigma)
    assert s.cancel(full-residual)==0
    p0=s.Function('p0')(z);u=4*z+s.Symbol('j');Haxis=z*(1-delta)/2+d*u
    hydro=((1-2*z*u)*u*(1+delta)/2+4*Haxis)/(2*L)
    pressure=(d*s.diff(p0,z)-2*z*p0*(1+delta))/(2*L)
    old=-(-((1-2*z*u)*u*(1+delta)/2+4*Haxis+d*P*s.diff(p0,z)-2*z*P*p0*(1+delta)))/(2*L)
    assert s.cancel(hydro+P*pressure-old)==0
    return dict(independent_exact_residual_gauge_and_pressure_slope_identities=2)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**live_drive_binding(),**independent_gauge_identity())
    axis=owner.report(packets['core_axis']);exit=owner.report(packets['selected_actual_inner_exit'])
    inlet=owner.report(packets['actual_first_inlet']);core_exit=owner.report(packets['core_exit'])
    assert all(not axis[key]['terms'] for key in ('factored_Psi','factored_Uz_minus_uaxis','factored_Hactual_minus_Haxis','factored_remaining_swirl_gauge'))
    assert inlet['factored_Psi']==core_exit['factored_Psi']
    assert not inlet['factored_actual_collar_Psi_increment']['terms']
    assert len(exit['actual_collar_Psi_source_cells'])==4
    for packet in packets.values():
        v=owner.report(packet);assert all(not v[k] for k in current.OPEN)
        assert v['remaining_normalized_B0_B1_feedback_not_solved'] and v['full_weighted_residual_matrix_majorant_still_open']
        assert v['original_actual_field_not_replaced_by_comparison']
        for key in ('factored_Uz_minus_uaxis','factored_Hactual_minus_Haxis','factored_remaining_swirl_gauge'):
            assert all(row['source_exponents'][3]==-1 and row['source_exponents'][2] in (0,1)
                       and row['source_exponents'][1]==0 for row in v[key]['terms'])
        assert all(ep(row['modulus_upper'])[1]<c.mpf('1e-12') for row in v['uniform_gauge_absolute_budgets'].values())
        proof=v['uniform_complex_bound_proof'];eta=owner.bound['common_tubes']['outer_Z_radius']
        assert proof['original_outer_eta_preserved']==eta and proof['gauge_Z1_Z_tube_radius']==eta/4
        assert proof['source_pressure_scale_not_materialized']
    checks['actual_axis_zero_and_same_source_core_first_inlet_Psi_identity']=True
    checks['actual_four_Psi_collar_cells_preserve_pressure_and_epsilon_scales']=True
    checks['remaining_gauge_and_Z1_absolute_norms_below_1e_minus12_on_original_tubes']=True
    checks['unsupported_endpoint_partition_or_packet_rejected']=all((
        rejected(lambda:owner.evaluate('.371',c.mpf([0,1]))),
        rejected(lambda:owner.evaluate('.371','1','first',True)),
        rejected(lambda:owner.report(current.OriginalN1AxialEpsilonPacket()))))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])];old=entry[1]['full_weighted_residual_matrix_majorant_still_open']
    entry[1]['full_weighted_residual_matrix_majorant_still_open']=False
    try:checks['false_matrix_completion_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['full_weighted_residual_matrix_majorant_still_open']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_AXIAL_EPSILON_FACTOR',checks['executable_original_source_assignment_bindings'],flush=True)
    return result
