"""Directed source replay and independent order-one operator/axis checks."""
import copy
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace

import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_core_hierarchy_source as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError):return True
    return False


def restore(c,value):return c.mpf(current.ends(value))


def independent_axial_second(c,grid,rho,z,delta,a):
    """Independent fully collected coefficients of Z_(a-1+delta) Z_a."""
    L=1-delta*z*z;d=1-z*z;b=a-1+delta
    # Derive the coefficient rows directly from the paper operators, rather
    # than differentiating or calling the producer's intermediate values.
    weights={(0,0):a*b*z*z/L**2+a*d/L**2+2*a*delta*z*z*d/L**3,
             (0,1):(a+b-2)*z*d/L**2+2*delta*z*d*d/L**3,
             (0,2):d*d/L**2,
             (1,0):(-2*(a+b-2)*z*z*rho-2*d*rho)/L**2-4*delta*z*z*d*rho/L**3,
             (1,1):-4*z*d*rho/L**2,
             (2,0):4*z*z*rho*rho/L**2}
    return sum((coefficient*grid[current.core_module.gridkey(i,k)] for (i,k),coefficient in weights.items()),c.mpf(0))


def independent_pressure_source(c,grid,rho,z,delta,Lambda):
    Q=current.core_module.Q;V=current.core_module.V;k=current.core_module.gridkey
    q,qr,qrr,qz=(grid[Q][k(i,j)] for i,j in ((0,0),(1,0),(2,0),(0,1)))
    uz=grid[V][k(0,0)];L=1-delta*z*z
    # Directly collected Omega0/R, with each R derivative converted.
    constant=(q+(1-delta)*z*qz/2+rho*qr+uz*(-2*z*q+(1-z*z)*qz-2*z*rho*qr))/L
    return -constant/2-q*q/4-q*rho*qr/2+Lambda*(2*qr+rho*qrr)


def exact_operator_fixture():
    rho,z,delta,a=s.symbols('rho Z delta a',real=True)
    L=1-delta*z*z;d=1-z*z
    f=1+rho/7+rho*rho*z/11+z*z/13+rho*z**3/17
    Za=lambda weight,fn:(weight*z*fn+d*s.diff(fn,z)-2*z*rho*s.diff(fn,rho))/L
    direct=Za(a-1+delta,Za(a,f))
    c=current.original.macro.MPIntervalContext();c.dps=800
    rational=lambda v:c.mpf(int(v.p))/int(v.q)
    count=0
    for rr in (s.Rational(0),s.Rational(1),s.Rational(41,10)):
        for zz in (s.Rational(0),s.Rational(7,19)):
            subs={rho:rr,z:zz,delta:s.Rational(1,17),a:s.Rational(-35,17)}
            grid={current.core_module.gridkey(i,k):rational(s.diff(f,rho,i,z,k).subs(subs))
                  for i in range(3) for k in range(3-i)}
            expected=rational(s.cancel(direct.subs(subs)))
            produced=current.axial_second(c,grid,rational(rr),rational(zz),rational(subs[delta]),rational(subs[a]))
            independent=independent_axial_second(c,grid,rational(rr),rational(zz),rational(subs[delta]),rational(subs[a]))
            assert current.contains(produced,expected) and current.contains(independent,expected)
            count+=1
    R=s.symbols('R',nonnegative=True);Q=1+R/7+R*R*z/11+z*z/13;Uz=2+z/3+R*z/5
    radial=R*Q
    Ta=lambda weight,fn:(-weight*fn/2+(1-delta)*z*s.diff(fn,z)/2+R*s.diff(fn,R))/L
    Zr=lambda weight,fn:(weight*z*fn+d*s.diff(fn,z)-2*z*R*s.diff(fn,R))/L
    omega=Ta(0,radial)+radial*(s.diff(radial,R)-radial/(2*R))+Uz*Zr(0,radial)-2*R*s.diff(radial,R,2)
    regular=s.cancel(-omega/(2*R))
    assert not regular.has(s.zoo) and s.denom(regular).subs(R,0)!=0
    for rr in (s.Rational(0),s.Rational(1),s.Rational(41,10)):
        zz=s.Rational(7,19);dt=s.Rational(1,17);lam=s.Integer(3)
        q=Q.subs(R,rho/lam);u=Uz.subs(R,rho/lam);sub={rho:rr,z:zz,delta:dt}
        grids={label:{current.core_module.gridkey(i,k):rational(s.diff(fn,rho,i,z,k).subs(sub))
               for i in range(3) for k in range(3-i)} for label,fn in ((current.core_module.Q,q),(current.core_module.V,u))}
        expected=rational(regular.subs({R:rr/lam,z:zz,delta:dt}))
        core=SimpleNamespace(delta=rational(dt),Lambda=rational(lam))
        profile=dict(rho=rational(rr),Z=rational(zz),ordinary_mixed_profile_grids=grids)
        assert current.contains(current.regular_radial_forcing(core,profile),expected)
        assert current.contains(independent_pressure_source(c,grids,rational(rr),rational(zz),core.delta,core.Lambda),expected)
        count+=1
    return dict(independent_exact_polynomial_operator_and_regular_axis_fixtures=count,
                ordered_second_axial_operator_not_square=True,regular_Omega0_over_R_identity=True)


@source_precision
def run(owner,values):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['source_family']==owner.family and raw['input_hashes']==owner.hashes
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    assert set(values)==set(current.SAMPLES) and all(not raw[k] for k in current.OPEN)
    fixture=exact_operator_fixture()
    replay=current.original.OriginalWholeZBridgeSource(800);core=current.OriginalCoreContext(replay);c=replay.c
    assert (replay.source,replay.family,replay.datum)==(owner.source.source,owner.source.family,owner.source.datum)
    function_rows=pressure_parts=forcing_rows=axis_rows=0
    evidence={}
    for name,value in values.items():
        entry=owner._fields[id(value)];packet=entry[2];view=owner.report(value)
        assert current.first.encode(current.original.serialized(view))==raw['actual_core_function_and_forcing_samples'][name]
        # Replay exact requested coordinates at the new precision. Reusing a
        # rounded 500-digit enclosure of 4.1 as an uncertain 800-digit input
        # would include points outside the exact endpoint.
        rho,z=(c.mpf(query) for query in current.SAMPLES[name])
        profile=core.profiles(rho,z);grid=profile['ordinary_mixed_profile_grids']
        amplitude=replay.amplitude.evaluate(z);logF=amplitude['logF0']
        for field,key,scale in (('F_0',current.core_module.F,logF),('Uz_0',current.core_module.V,c.mpf(0)),
                                ('V_0_over_R',current.core_module.Q,c.mpf(0))):
            for i in range(5):
                for k in range(5-i):
                    index=current.core_module.gridkey(i,k);row=packet['actual_leading_function_derivatives'][field][index]
                    assert row['radial_derivative']==i and row['ordinary_Z_derivative']==k
                    assert current.contains(row['coefficient'],grid[key][index]),(name,field,index)
                    assert current.contains(row['log_scale'],scale+i*core.logLambda)
                    function_rows+=1
        for i in range(5):
            for k in range(5-i):
                index=current.core_module.gridkey(i,k);row=packet['actual_leading_function_derivatives']['P_0'][index]
                for part,key,scale in (('axis_datum',current.core_module.PD,2*core.logP+i*core.logLambda),
                                      ('centrifugal_increment',current.core_module.PI,2*logF+(i-1)*core.logLambda)):
                    assert current.contains(row[part]['coefficient'],grid[key][index])
                    assert current.contains(row[part]['log_scale'],scale);pressure_parts+=1
                if i>0:assert current.ends(row['axis_datum']['coefficient'])==(0,0)
        expected=dict(N1theta=-independent_axial_second(c,grid[current.core_module.F],rho,z,core.delta,-2-core.delta),
                      N1z=-independent_axial_second(c,grid[current.core_module.V],rho,z,core.delta,-1-core.delta),
                      N1p=independent_pressure_source(c,grid,rho,z,core.delta,core.Lambda))
        for key,bound in expected.items():
            row=packet['genuine_order_one_known_forcing'][key]
            assert current.contains(row['coefficient'],bound),(name,key)
            assert current.contains(row['log_scale'],logF if key=='N1theta' else c.mpf(0));forcing_rows+=1
        if current.ends(rho)==(0,0):
            for k in range(5):
                index=current.core_module.gridkey(0,k)
                assert current.ends(grid[current.core_module.PI][index])==(0,0)
                # The source axis law is Uz0(0,Z)=4Z+j, not positive-order axis data.
                uz=(4*z+core.j if k==0 else c.mpf(4) if k==1 else c.mpf(0))
                assert current.contains(packet['actual_leading_function_derivatives']['Uz_0'][index]['coefficient'],uz)
                axis_rows+=1
        assert packet['radial_forcing_has_no_division_by_R'] and all(not packet[k] for k in current.OPEN)
        evidence[name]=dict(real_domain=packet['real_common_domain'],forcing_signs={key:
            1 if current.ends(row['coefficient'])[0]>0 else -1 if current.ends(row['coefficient'])[1]<0 else None
            for key,row in packet['genuine_order_one_known_forcing'].items()})
    first=next(iter(values.values()));entry=owner._fields[id(first)]
    invalid=dict(copied_packet=rejected(lambda:owner.report(copy.copy(first))),
                 outside_radial_domain=rejected(lambda:owner.evaluate('4.2','0')),
                 outside_Z_domain=rejected(lambda:owner.evaluate('0','1.1')))
    row=entry[2]['genuine_order_one_known_forcing']['N1z'];saved=row['coefficient']
    try:
        row['coefficient']=owner.ctx.mpf(0)
        invalid['changed_known_forcing']=rejected(lambda:owner.report(first))
    finally:row['coefficient']=saved
    saved=owner.core.delta
    try:
        owner.core.delta=owner.ctx.mpf(0)
        invalid['changed_original_delta']=rejected(lambda:owner.report(first))
    finally:owner.core.delta=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    receipt=dict(all_passed=True,source_family=owner.family,**fixture,
        independently_replayed_leading_function_derivative_rows=function_rows,
        independently_replayed_pressure_parts=pressure_parts,
        independently_checked_genuine_order_one_forcing_rows=forcing_rows,
        independently_checked_original_axis_rows=axis_rows,
        source_replay_digits=800,real_source_evidence=evidence,invalid_inputs_rejected=invalid,
        inherited_analytic_fixed_point_bounds_not_new_point_solution=True,
        common_holomorphic_extension_not_promoted=True,input_hashes=hashes,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.first.encode(current.original.serialized(receipt)),indent=2)+'\n',
                                          encoding='utf8',newline='\n')
    print('PASS_ORIGINAL_CORE_HIERARCHY_SOURCE',function_rows,pressure_parts,forcing_rows,axis_rows,flush=True)
    return receipt
