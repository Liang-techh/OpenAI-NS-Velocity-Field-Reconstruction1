"""Independent natural scaling, signed source identities and total budgets."""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_natural_component_enclosure as current
import lei_ren_part1_paper_compliant_current_core_first_interface_check as interface
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported natural component scope admitted')


def exact_scaling_and_primitives():
    z,x,r,eps,pi,delta=s.symbols('Z x rho epsilon Pi delta',nonzero=True)
    original=interface.function('compliant_current_original_n1_inner_operator','matrix')
    env=dict(truncate=lambda v,n:v,dz=lambda v:s.diff(v,z))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(original)],type_ignores=[])),
                 '<original forcing and principal rows>','exec'),env)
    alg=SimpleNamespace(ctx=SimpleNamespace(sqrt=s.sqrt,mpf=s.Rational),lift=lambda v,*args:v,
                        shift=lambda v,p:v*eps**(-p[3])*pi**p[2])
    Np,Nz,Ntheta=s.symbols('Np Nz Ntheta')
    raw=env['matrix'](alg,dict(F=[0,0],Uz=[0,0],Q=[0]),dict(N1p=Np,N1z=Nz,N1theta=Ntheta),z,r,delta)
    assert s.cancel(raw['g']['4']/(eps*pi)-2*s.sqrt(r)*Np/pi)==0
    assert s.cancel(raw['g']['6']/eps-(2*Nz-eps*4*z*r*Np/(1-delta*z*z)))==0
    T=[s.Integer(1),eps,eps,eps*pi,s.Integer(1),eps];count=2
    matrices={label:{key:s.Symbol(label+key.replace(',','_')) for key in raw[label]} for label in ('B0','B1')}
    got=current.conjugate_extra(alg,matrices)
    for label,rows in matrices.items():
        for key,v in rows.items():
            i,j=(int(k)-1 for k in key.split(','))
            assert s.cancel(got[label][key]-v*T[j]/T[i])==0;count+=1
    for i,j in ((0,4),(1,5),(2,5),(4,0)):assert s.cancel(T[j]/T[i]-1)==0
    assert all(s.diff(v,z)==0 and s.diff(v,x)==0 for v in T)
    q=s.Function('Qpsi')(r,z);primitive=q-q.subs(r,0)+r*s.diff(q,r)
    assert s.cancel(s.diff(primitive,r)-(2*s.diff(q,r)+r*s.diff(q,r,2)))==0
    assert s.cancel(primitive.subs(r,0))==0
    C=s.Function('C')(z);U=C*x*x/4;K=-C*x*x/8;D=C*x/2
    assert s.simplify(s.diff(U,x)-D)==0
    assert s.simplify(s.diff(D,x)+D/x-C)==0
    assert s.simplify(s.diff(K,x)+2*K/x+D)==0
    assert s.simplify(s.integrate(U,(x,0,x))-C*x**3/12)==0
    return dict(executable_original_forcing_and_extra_diagonal_scaling_identities=count,
                extra_T_constant_in_Z_and_x=True,principal_kinematic_and_swirl_potential_unchanged=True,
                exact_viscous_pressure_FTC_primitive_verified=True,
                exact_axis_forced_regular_axial_U_K_dxU_equations_verified=True)


def rational(c,v):
    v=s.cancel(v)
    return c.mpf(int(s.numer(v)))/int(s.denom(v))


@source_precision
def independent_seed_fixtures(c):
    z,r,dt,j=s.symbols('z rho delta j')
    u=4*z+j;d=1-z*z;L=1-dt*z*z
    psi=(2+3*z+z*z)*r+(1-z)*r*r
    mu=s.cancel(s.integrate(psi,(r,0,r))/r)
    qp=s.cancel((2*z*psi-(1-dt)*z*mu-d*s.diff(mu,z))/L)
    T=lambda weight,f:(weight*z*f+d*s.diff(f,z)-2*z*r*s.diff(f,r))/L
    C=s.cancel(-2*T(-2,T(-1-dt,u)));U=C*r/4;K=-C*r/8;M=C*r/8
    expected=dict(Uz1=U,K1=K,P1=qp-qp.subs(r,0)+r*s.diff(qp,r),d_x_Uz1=C*s.sqrt(r)/2,
                  mean1=M,Q1=((1-dt)*z*U-(1+dt)*z*K-d*s.diff(M,z))/L)
    expected['V1']=r*expected['Q1'];count=0
    for zz,rr,dd,jj in (('0.25','1','0.2','0.1'),('-0.4','4','0.15','-0.2')):
        sub={z:s.Rational(zz),r:s.Rational(rr),dt:s.Rational(dd),j:s.Rational(jj)}
        algebra=op.N1LogAlgebra(c,(c.mpf(0),)*4,[])
        jet=lambda f,n:op.Jet(c,[rational(c,s.diff(f,z,k).subs(sub))/math.factorial(k) for k in range(n+1)])
        rows=dict(axis=dict(u=jet(u,5)))
        fields=dict(Qpsi=[algebra.lift(jet(s.diff(qp,r,i),4-i)) for i in range(3)])
        values=current.seeds(algebra,rows,fields,c.mpf(zz),c.mpf(rr),c.mpf(dd))
        for key,v in values.items():
            for k in range(v.order+1):
                got=sum((row[k] for row in v.terms.values()),c.mpf(0))
                want=rational(c,s.diff(expected[key],z,k).subs(sub))/math.factorial(k)
                assert ep(got)[0]<=ep(want)[1] and ep(got)[1]>=ep(want)[0],(key,k);count+=1
        # Check the triangle operator on independently signed radial/Z jets.
        f=(1-z+2*z*z)*r+(2+z)*r*r
        exact=-T(-2,T(-1-dt,f))
        vals=[abs(rational(c,v.subs(sub))) for v in (f,s.diff(f,r),s.diff(f,r,2),s.diff(f,z),s.diff(f,z,2),s.diff(f,r,z))]
        bound=current.axial_majorant(*vals,abs(c.mpf(zz)),abs(1-c.mpf(zz)**2),1-c.mpf(dd)*c.mpf(zz)**2,c.mpf(rr),abs(c.mpf(dd)))
        assert ep(bound)[1]>=ep(abs(rational(c,exact.subs(sub))))[1]
    return dict(independent_signed_axial_pressure_mean_Q_V_ordinary_Z_fixture_rows=count,
                axial_majorant_covers_independently_signed_second_operator_fixtures=True)


@source_precision
def run(owner,domain_packet,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**exact_scaling_and_primitives(),**independent_seed_fixtures(c))
    v=owner.report(domain_packet);data=owner._data;eta=data['eta'];tail=v['all_word_feedback_tail'];proof=v['uniform_proof']
    assert all(not v[k] for k in current.OPEN)
    assert proof['original_outer_eta']==owner.bound['common_tubes']['outer_Z_radius']
    assert proof['matrix_bound_tube']==eta/4 and proof['source_M0_tube']==eta/8
    assert tail['iteration_Cauchy_loss']==eta/16 and tail['value_tube_radius']==eta/16
    assert tail['ordinary_Z_final_Cauchy_gap']==eta/32 and tail['ordinary_Z1_Z2_tube_radius']==eta/32
    assert ep(tail['ordinary_Z2_remainder_upper'])[1]<c.mpf('1e-12')
    assert set(data['C0entries'])=={'4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'}
    assert set(data['C1entries'])=={'5,1','5,2','5,3','6,2','6,3','6,4'}
    assert all(p[0]>=0 for value in data['collar_L1'].values() for p in value.terms)
    for cover in (data['C0cover'],data['C1cover'],data['source_cover'],data['pressure_rest'],data['axial_rest']):
        assert cover['all_covers_have_directed_log_proof']
        for row in cover['sectors']:
            if row['rule']=='proved log<=log(floor)':assert ep(row['log_modulus_upper'])[1]<=ep(c.ln(cover['tiny_cover_floor']))[0]
            else:assert row['rule']=='directed exp(log)'
    axis=owner.report(packets['core_axis'])
    assert axis['exact_zero_axis_data']
    assert all(not row['terms'] for row in axis['normalized_signed_seed'].values())
    assert all(not row['terms'] for row in axis['normalized_total_error'].values())
    assert axis['passing_normalized_width_rows']==0
    for name,packet in packets.items():
        view=owner.report(packet);values,errors,physical,physical_errors=owner._packets[id(packet)][3]
        assert all(not view[k] for k in current.OPEN)
        assert view['full_coefficient_function_and_equation_contracts_still_open'] and view['n1_own_five_moment_repair_and_n2_still_open']
        if name!='core_axis':
            assert view['passing_normalized_width_rows']==14 and view['passing_physical_width_rows']==14
            for label,rows in view['normalized_component_budgets'].items():
                for key,row in rows.items():assert row['remaining_B0_B1_feedback_included'] and row['all_nonseed_sources_included']
        for label,value in values.items():
            assert value.order==(1 if label in ('Q1','V1') else 2)
            assert op.export(physical[label])==op.export(value.algebra.shift(value,current.SCALES[label]))
            assert op.export(physical_errors[label])==op.export(value.algebra.shift(errors[label],current.SCALES[label]))
    inlet=owner.report(packets['actual_first_inlet']);exit=owner.report(packets['selected_actual_inner_exit'])
    assert inlet['normalized_signed_seed']==owner.report(packets['core_exit'])['normalized_signed_seed']
    # Execute the exact newly scaled residual on genuine own-mean/Q fields.
    mp=owner._packets[id(packets['selected_actual_inner_exit'])][2];me=owner.mean._packets[id(mp)]
    mv=owner.mean.report(mp);rows,fields=me[4][1:3];algebra=fields['Uz'][0].algebra
    ap=owner.mean.axial.evaluate('.371','1');ae=owner.mean.axial._packets[id(ap)]
    z=algebra.lift(op.Jet.variable(c,mv['Z'],5));rho=mv['geometry']['rho_enclosure_only']
    grad=op.Jet(c,owner.weighted.source.amplitude.gradient_jets(mv['Z']))
    actual=current.conjugate_extra(algebra,current.feedback.residual_matrix(algebra,rows['Phi'],fields['Uz'],fields['Q'][0],
        z,rho,owner.physical.delta,grad,ae[5][4],1))
    assert all(p[1:3]==(2,-1) and p[3] in (1,2) for p in actual['B0']['4,1'].terms)
    assert all(p[1:]==(2,-1,1) for p,row in actual['B0']['4,1'].terms.items() if ep(row[0])!=(0,0))
    assert all(p[1:3]==(2,0) and p[3] in (0,1) for p in actual['B0']['6,1'].terms)
    assert all(p[1:]==(2,0,0) for p,row in actual['B0']['6,1'].terms.items() if ep(row[0])!=(0,0))
    checks['genuine_actual_own_mean_and_Q_natural_residual_executed']=True
    checks['original_eta_and_all_six_row_tail_through_ordinary_Z2_checked']=True
    checks['every_norm_floor_has_original_directed_log_proof']=True
    checks['actual_inlet_and_corrected_exit_fourteen_normalized_and_physical_width_rows_pass']=True
    checks['physical_epsilon_pressure_scales_and_exact_zero_axis_preserved']=True
    checks['no_false_whole_Z_or_repaired_n1_scope_promoted']=True
    checks['foreign_packet_rejected']=rejected(lambda:owner.report(current.OriginalN1NaturalComponentPacket()))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])];old=entry[1]['full_coefficient_function_and_equation_contracts_still_open']
    entry[1]['full_coefficient_function_and_equation_contracts_still_open']=False
    try:checks['false_full_function_scope_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['full_coefficient_function_and_equation_contracts_still_open']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf8',newline='\n')
    print('PASS_N1_NATURAL_COMPONENT_ENCLOSURE',checks['independent_signed_axial_pressure_mean_Q_V_ordinary_Z_fixture_rows'],flush=True)
    return result
