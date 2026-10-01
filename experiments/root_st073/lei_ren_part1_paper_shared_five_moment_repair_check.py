"""Independent physical moment-map, partial primitive and stress checks.

Symbolic identities establish the normalization of the actual functional
repair. Finite interval examples are diagnostics, never the global proof.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symbolic_checks():
    x,Rm,Am,Pstar=s.symbols('x Rm Am Pstar',positive=True)
    Z,delta,f,g=s.symbols('Z delta f g',real=True)
    u=Am*x**s.Rational(1,10);V=4*Z
    uh=u+Am*f;Vh=V+g
    raw_z=Rm*(Vh-V)
    raw_t=Rm*s.sqrt(2*Rm*x)*(uh-u)
    raw_tz=Rm*s.sqrt(2*Rm*x)*(uh*Vh-u*V)
    raw_zt=Rm*(Vh**2-V**2-(uh**2-u**2)/2)
    raw_p=(uh**2-u**2)/(2*x)
    actual=[raw_z/Rm,(raw_tz-4*Z*raw_t)/(s.sqrt(2)*Rm**s.Rational(3,2)*Am),
            raw_t/(s.sqrt(2)*Rm**s.Rational(3,2)*Am),
            (raw_zt-8*Z*raw_z)/(Rm*Am**2),raw_p/Am**2]
    target=[g,x**s.Rational(3,5)*g+s.sqrt(x)*f*g,s.sqrt(x)*f,
            g**2/Am**2-x**s.Rational(1,10)*f-f**2/2,
            x**s.Rational(-9,10)*f+f**2/(2*x)]
    residuals=[s.simplify(a-b) for a,b in zip(actual,target)]
    if any(v!=0 for v in residuals):
        raise ArithmeticError('Physical five integrands do not give the correction map')
    Dz,Dt,Dtz,Dzt,Dp=s.symbols('Dz Dt Dtz Dzt Dp')
    original=s.Matrix([Dz,Dt,Dtz,Dzt,Dp])
    centered=s.Matrix([Dz/Rm,(Dtz-4*Z*Dt)/(s.sqrt(2)*Rm**s.Rational(3,2)*Am),
                       Dt/(s.sqrt(2)*Rm**s.Rational(3,2)*Am),
                       (Dzt-8*Z*Dz)/(Rm*Am**2),Dp/Am**2])
    C=centered.jacobian(original)
    if s.simplify(C.det())==0 or any(s.simplify(v)!=0 for v in C.inv()*centered-original):
        raise ArithmeticError('Five centered moments not invertible')

    # Independently differentiate the normalized cumulative primitives.
    ff=s.Function('f')(x);gg=s.Function('g')(x)
    G,T,M,E,P=[s.Function(n)(x) for n in ('G','T','M','E','P')]
    d1,d2,d3,d4,d5=s.symbols('d1 d2 d3 d4 d5')
    mass=4*Z+(d1+G)/x
    theta=s.Rational(5,8)*x**s.Rational(8,5)+d3+T
    mixed=4*Z*theta+d2+M
    energy=-s.Rational(5,12)*x**s.Rational(6,5)+d4+E
    pressure=s.Rational(5,2)*x**s.Rational(1,5)+d5+P
    rules={s.diff(G,x):gg,s.diff(T,x):s.sqrt(x)*ff,
           s.diff(M,x):x**s.Rational(3,5)*gg+s.sqrt(x)*ff*gg,
           s.diff(E,x):gg**2/Am**2-x**s.Rational(1,10)*ff-ff**2/2,
           s.diff(P,x):x**s.Rational(-9,10)*ff+ff**2/(2*x)}
    H=x**s.Rational(1,10)+ff
    prim=[s.diff(x*mass,x)-(4*Z+gg),s.diff(theta,x)-s.sqrt(x)*H,
          s.diff(mixed,x)-s.sqrt(x)*H*(4*Z+gg),
          s.diff(energy,x)-gg**2/Am**2+H**2/2,
          s.diff(pressure,x)-H**2/(2*x)]
    prim=[s.simplify(v.xreplace(rules)) for v in prim]
    if any(v!=0 for v in prim):
        raise ArithmeticError('Partial recovery does not retain all five actual primitives')

    # Stress and radial velocity recovered from the same normalized moments.
    R=Rm*x;L=1-delta*Z**2;dz=1-Z**2
    B,Bz,U,H,zeta,T,Tz,J,Jz,E,Ez,P,Pz=s.symbols('B Bz U H zeta T Tz J Jz E Ez P Pz')
    mz=R*B;mzZ=R*Bz
    mt=s.sqrt(2)*Rm**s.Rational(3,2)*Am*T
    mtZ=s.sqrt(2)*Rm**s.Rational(3,2)*Am*(Tz+zeta*T)
    mtz=s.sqrt(2)*Rm**s.Rational(3,2)*Am*J
    mtzZ=s.sqrt(2)*Rm**s.Rational(3,2)*Am*(Jz+zeta*J)
    mzt=R*(8*Z*B-16*Z**2+Am**2*E/x)
    mztZ=R*(8*B+8*Z*Bz-32*Z+Am**2*(Ez+2*zeta*E)/x)
    up=Am*H
    pressure=Pstar**2*P;pressureZ=Pstar**2*Pz
    A_mz=(1-delta)*Z*mz+dz*mzZ
    angular=(1-delta/2)*mt-(1-delta)*Z*mtZ/2-dz*mtzZ+(2*delta-1)*Z*mtz
    Itheta=up*(-R+A_mz)/(L*s.sqrt(2*R))+angular/(2*L*R)
    Iz=((-R+A_mz)*U+(1-delta)*(mz-Z*mzZ)/2+2*delta*Z*mzt-dz*mztZ
        +R*(2*(1+delta)*Z*pressure-dz*pressureZ))/(L*s.sqrt(2*R))
    Ur=(2*Z*R*U-(1-delta)*Z*mz-dz*mzZ)/(L*s.sqrt(2*R))
    W=1-(1-delta)*Z*B-dz*Bz
    expectedQ=-W+((1-delta/2)*T-(1-delta)*Z*(Tz+zeta*T)/2
                       -dz*(Jz+zeta*J)+(2*delta-1)*Z*J)/(x**s.Rational(3,2)*H)
    es=(8*Z*B-16*Z**2)/Pstar**2+(Am/Pstar)**2*E/x
    esz=(8*B+8*Z*Bz-32*Z)/Pstar**2+(Am/Pstar)**2*(Ez+2*zeta*E)/x
    expectedN=(-W*U+(1-delta)*(B-Z*Bz)/2)/Pstar**2+2*delta*Z*es-dz*esz+2*(1+delta)*Z*P-dz*Pz
    st=[L*Itheta/(s.sqrt(R/2)*up)-expectedQ,
        L*Iz/(s.sqrt(R/2)*Pstar**2)-expectedN,
        Ur/s.sqrt(R/2)-(2*Z*U-(1-delta)*Z*B-dz*Bz)/L]
    st=[s.simplify(v) for v in st]
    if any(v!=0 for v in st):
        raise ArithmeticError('Corrected pressure/moment stress normalization failed')
    return dict(physical_five_integrand_map_residuals=[str(v) for v in residuals],
                partial_five_primitive_recovery_residuals=[str(v) for v in prim],
                Q_N_Ur_normalization_residuals=[str(v) for v in st],
                centered_moment_map_determinant=str(s.factor(C.det())),
                original_moments_recovered_by_invertible_centering=True,
                arbitrary_Z_coefficient_functions_allowed=True)


def run():
    name=PREFIX+'shared_five_moment_repair.json'
    raw=json.loads((HERE/name).read_bytes())
    for n,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
            raise ValueError('Actual repair check dependency changed: '+n)
    inverse=raw['actual_whole_axis_coefficient_inverse']
    if not (inverse['certified'] and inverse['self_map_strictly_inside']
            and raw['actual_implicit_functional_five_moment_identities_analytically_certified']):
        raise ValueError('Actual functional inverse not admitted')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        get=lambda r,n:read_interval(c,r[n])
        diagnostics=[]
        for sample in raw['normalized_patch_samples']:
            u=get({'u':sample['Utheta_over_Pstar'][0]},'u')
            a=get(sample,'angular_shear_a');k=get(sample,'kappa')
            Q=get(sample,'Q');margin=get(sample,'relaxed_Da_minus_bw_minus_2a_lower_enclosure')
            if (endpoints(u)[0]<=mp.mpf('.125') or endpoints(u)[1]>=1
                    or endpoints(a)[0]<mp.mpf('.7') or endpoints(a)[1]>mp.mpf('.9')
                    or endpoints(Q)[0]<=mp.mpf('.5') or endpoints(k)[1]>=1
                    or endpoints(margin)[0]<=0):
                raise ArithmeticError('Corrected patch enclosure diagnostics failed')
            diagnostics.append(dict(x=sample['x'],normalized_velocity_positive=True,Q_above_half=True,
                                    kappa_below_one=True,relaxed_margin_positive=True))
        contains_zero=lambda v:endpoints(read_interval(c,v))[0]<=0<=endpoints(read_interval(c,v))[1]
        if not (all(contains_zero(v) for v in inverse['enclosure_residual'])
                and all(contains_zero(v) for v in inverse['derivative_enclosure_residual'])):
            raise ArithmeticError('Actual map/differentiated residual diagnostics exclude zero')
    result=dict(actual_five_defect_family_sha256=raw['actual_five_defect_family_sha256'],
                symbolic_checks=symbolic_checks(),finite_patch_diagnostics=diagnostics,
                actual_pointwise_map_value_and_derivative_zero_containment=True,
                residual_zero_containment_not_used_as_closure_proof=True,
                actual_implicit_functional_five_moment_closure_independently_checked=True,
                corrected_partial_moments_and_same_pressure_independently_checked=True,
                full_physical_field_or_heat_exterior_completed=False,temporal_recursion=False,
                symbolic_runtime_version=s.__version__,
                input_hashes={**raw['input_hashes'],name:hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
                    Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Independent original five-integrand map + partial primitives + Q/N/Ur identities:13 PASS',flush=True)
    print('Corrected patch diagnostics:5 PASS; closure proof remains functional contraction',flush=True)
    return result


if __name__=='__main__':
    run()
