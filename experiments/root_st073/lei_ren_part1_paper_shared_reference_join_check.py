"""Independent moment-to-stress identities and log-interface checks.

Derives (9.32)/(9.35) from the actual five primitives and (9.13), rather
than testing copies of the scalar ODEs. Fixtures check the callable log
interfaces; they do not establish the whole-axis cone.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_shared_reference_join_bounds import (
    shape_log_jet,reference_log_jet,restore_axial_jet)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def moment_stress_identities():
    R,Z=s.symbols('R Z',positive=True)
    delta,a,zeta,Vy=s.symbols('delta a zeta Vy')
    u=s.Function('u')(R,Z);V=s.Function('V')(R,Z);P=s.Function('P')(R,Z)
    mt,mz,mtz,mzt=[s.Function(n)(R,Z) for n in ('Mtheta','Mz','Mtheta_z','Mztheta')]
    d=1-Z**2;L=1-delta*Z**2
    Aop=lambda q:(1-delta)*Z*q+d*s.diff(q,Z)
    Pop=lambda q:2*(1+delta)*Z*q-d*s.diff(q,Z)
    W=1-Aop(mz/R);Hv=(1-delta)*Z/2+d*V
    # Direct inertial stress numerator, independently normalized from9.13.
    angular=((1-delta/2)*mt-(1-delta)*Z*s.diff(mt,Z)/2
             -d*s.diff(mtz,Z)+(2*delta-1)*Z*mtz)
    Q=(-R+Aop(mz))/R+angular/(s.sqrt(2)*R**s.Rational(3,2)*u)
    N=((-R+Aop(mz))*V+(1-delta)*(mz-Z*s.diff(mz,Z))/2
       +2*delta*Z*mzt-d*s.diff(mzt,Z)+R*Pop(P))/R
    # The five exact axis-primitive RHS, with pressure=P0+Mp.
    rhs={mt:s.sqrt(2*R)*u,mz:V,mtz:s.sqrt(2*R)*u*V,
         mzt:V**2-u**2/2,P:u**2/(2*R)}
    updates={}
    for q,f in rhs.items():
        updates[s.diff(q,R)]=f
        updates[s.diff(s.diff(q,Z),R)]=s.diff(f,Z)
    updates[s.diff(u,R)]=(1-a)*u/(2*R)
    updates[s.diff(V,R)]=Vy/R
    # First restore primitives, then log derivative u_Z=zeta*u.
    q_res=R*s.diff(Q,R)+(2-a/2)*Q+W*(1-a/2)+delta*(1-2*Z*V)/2+Hv*zeta
    n_res=R*s.diff(N,R)+N-Z*u**2-Pop(P)+W*Vy+(1+delta)*(1-2*Z*V)*V/2+Hv*s.diff(V,Z)
    reduce=lambda q:s.factor(s.together(q.xreplace(updates).subs(s.diff(u,Z),zeta*u)))
    qzero=reduce(q_res);nzero=reduce(n_res)
    if qzero!=0 or nzero!=0:
        raise ArithmeticError('Five actual primitives do not imply the stress ODEs: '+str((qzero,nzero)))
    t,j,sigma,Z0=s.symbols('t j sigma Z0',real=True)
    H0=(1-delta)*Z/2+(1-Z**2)*(4*Z+j)
    factored=(4*Z+j)*(1-Z**2)+(1-delta)*Z/2
    g=L*H0/(H0**2+sigma**2)
    primitive=s.Integral(g.subs(Z,t),(t,Z0,Z))
    G_derivative=s.factor(s.diff(primitive,Z)-g)
    if s.expand(H0-factored)!=0 or G_derivative!=0 or primitive.subs(Z,Z0).doit()!=0:
        raise ArithmeticError('Implicit real-integral G branch identity failed')
    return dict(angular_9_32_residual=str(qzero),axial_9_35_residual=str(nzero),
                H0_sign_factorization_residual='0',G_integral_derivative_residual=str(G_derivative),
                G_integral_at_anchor='0',pole_log_primitive_used=False,
                pressure_radial_primitive='P_R=u^2/(2R), same axis datum',
                V_radial_derivative_arbitrary=True,
                independent_five_primitive_to_stress_derivation_passed=True)


def log_interfaces():
    c=MPIntervalContext();c.dps=100
    with mp.workdps(140):
        z=IntervalTaylor.variable(c,c.mpf('.2'),2)
        B=IntervalTaylor(c,[c.mpf('-3'),c.mpf('.7'),c.mpf('-.2')])
        A=c.mpf(10);logC=c.mpf(400);logP=c.mpf(1)
        start=shape_log_jet(c,0,A,logC,logP,z,B)
        end=shape_log_jet(c,1,A,logC,logP,z,B)
        expected=reference_log_jet(c,400*A-10*(logC+logP),z)
        def overlap(x,y):
            xl,xh=endpoints(x);yl,yh=endpoints(y)
            return max(xl,yl)<=min(xh,yh)
        checks=dict(shape_start_slope=endpoints(start['angular_shear_a'])==endpoints(c.mpf('.8')),
                    shape_end_slope=endpoints(end['angular_shear_a'])==endpoints(c.mpf('.8')),
                    reference_at_shape_end=all(overlap(x,y) for x,y in
                        zip(end['log_Utheta_over_Pstar_coefficients'],expected.coefficients)))
        v1=z*4+IntervalTaylor(c,[c.mpf('.01'),c.mpf('-.02'),c.mpf('.03')])
        vr0=restore_axial_jet(c,0,z,v1);vr1=restore_axial_jet(c,1,z,v1)
        checks['restore_start_inherits_v1']=all(overlap(x,y) for x,y in zip(vr0['Uz_coefficients'],v1.coefficients))
        checks['restore_end_is_4Z']=all(overlap(x,y) for x,y in zip(vr1['Uz_coefficients'],(z*4).coefficients))
        checks['restore_endpoint_radial_derivatives_zero']=all(endpoints(x)==(mp.mpf(0),mp.mpf(0))
             for v in (vr0,vr1) for x in v['Uz_y_coefficients'])
        # Test exact offset use on a base so large its rounded subtraction
        # cannot resolve a one-unit difference. The API never needs the base.
        base=c.mpf('1e200')
        if endpoints((base-7)-(base-8))[0]>1 or endpoints((base-7)-(base-8))[1]<1:
            raise ArithmeticError('Directed arithmetic lost its containment')
        minus8=reference_log_jet(c,-8,z);minus7=reference_log_jet(c,-7,z)
        difference=minus7[0]-minus8[0]
        checks['one_unit_reference_offset_preserved']=overlap(difference,c.mpf('.1')) and endpoints(difference)[1]-endpoints(difference)[0]<mp.mpf('1e-90')
        if not all(checks.values()):
            raise ArithmeticError('Callable log-interface check failed: '+str(checks))
        return checks


def run():
    name=PREFIX+'shared_reference_join_bounds.json'
    receipt=json.loads((HERE/name).read_bytes())
    for n,digest in receipt['input_hashes'].items():
        if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
            raise ValueError('Reference join source changed: '+n)
    result=dict(moment_stress_identities=moment_stress_identities(),log_interface_checks=log_interfaces(),
                reference_join_family_sha256=receipt['reference_join_family_sha256'],
                symbolic_derivation_not_sampling=True,fixtures_are_not_global_cone_proof=True,
                full_velocity_or_heat_exterior_built=False,temporal_recursion=False,
                input_hashes={**receipt['input_hashes'],name:hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
                    Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                symbolic_runtime_version=s.__version__)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Independent actual five-moment -> angular/axial stress identities:2 PASS',flush=True)
    print('Callable log-interface checks:',len(result['log_interface_checks']),'PASS',flush=True)
    return result


if __name__=='__main__':
    run()
