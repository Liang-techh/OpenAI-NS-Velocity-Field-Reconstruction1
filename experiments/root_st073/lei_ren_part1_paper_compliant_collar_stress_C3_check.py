"""Focused actual-history, full future moment and signed collar stress check."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_collar_stress_C3 import (
    CompliantCollarStressC3, collar_defect_rows, collar_stress_rows,
    collar_moment_stress_identities, collar_stress_source_bridge,collar_Gamma_endpoint_binding)
from lei_ren_part1_paper_compliant_collar_pressure_C4 import CompliantCollarPressureC4
from lei_ren_part1_paper_compliant_heat_terminal_history_bridge import terminal_history_bridge
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_full_moment_fixture():
    """Differentiate full defining integrals in Z, then apply paper formulas.

    Gamma tails are evaluated to infinity via compact changes of variable.
    The angular Gamma tail uses its full-function ODE/FTC antiderivative.
    This fixture has moderate parameters, not the extreme project scales.
    """
    with mp.workdps(42):
        c=MPIntervalContext(); c.dps=100
        a,S,eps,Z=map(mp.mpf,('.15','.004','.002','.4')); delta=2*a; k=1-a; p=1+delta; b=(1-delta)/2
        heat=CompliantCollarGammaC4.__new__(CompliantCollarGammaC4)
        heat.ctx=c; heat.a=c.mpf(a); heat.Scap=c.mpf(S); heat.S=c.mpf([0,S]); heat.eps=c.mpf(eps)
        heat.k=c.mpf(k); heat.delta=c.mpf(delta); heat.prate=c.mpf(p)
        heat.cells=64; heat.shape_cache={}; heat.tail_cache={}; heat.gamma_cache={}
        def H(x,n=0):
            if not x:return (-1)**n*mp.rf(a,n)*mp.rf(1+a,n)
            return (-1)**n*mp.rf(a,n)*mp.rf(1+a,n)*x**(-1-a-n)*mp.hyperu(1+a+n,2,1/x)
        def sigma(v):
            if v<=0:return mp.mpf(0)
            if v>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/v**2-1/(1-v)**2))
        def Kjet(v,z,n):
            sig=sigma(v); phi=mp.exp(-4/(3-v)**2) if v<3 else mp.mpf(0)
            ss=S*mp.exp(-v); x=2*(1-z*z)*ss
            value=mp.factorial(n)*sum(H(x,n-j)*(-4*z*ss)**(n-2*j)*(-2*ss)**j/
                                     (mp.factorial(n-2*j)*mp.factorial(j)) for j in range(n//2+1))
            return sig*(1-eps*phi)*value+((1-sig)*(1-eps) if n==0 else 0)
        def squarejet(v,z,n):
            return sum(mp.binomial(n,j)*Kjet(v,z,j)*Kjet(v,z,n-j) for j in range(n+1))
        def gamma_A(z):
            x=2*(1-z*z)*S*mp.exp(-3)
            return (H(x)+x*((1+a)*H(x)+x*H(x,1)))/k
        def gamma_square(z,q,rate,n):
            # q=x^(1/rate) removes the original q^(rate-1) endpoint weight.
            v=3-mp.log(q)/rate if q else mp.inf
            if not q:return mp.mpf(1) if n==0 else mp.mpf(0)
            return squarejet(v,z,n)/rate
        def encloses(value,reference,label):
            lo,hi=endpoints(value)
            if not lo<=reference<=hi:raise ArithmeticError('Independent collar stress/moment not enclosed: '+label)
        future_checks=stress_checks=nonzero=0
        for t in map(mp.mpf,('.5','2')):
            print('Independent full collar moments at offset '+str(t),flush=True)
            breaks=[t]+([mp.mpf(1)] if t<1 else [])+[mp.mpf(3)]
            A=[]; E=[]; P=[]
            for n in range(3):
                Aj=mp.quad(lambda v:mp.exp(k*(v-t))*((1-Kjet(v,Z,0)) if n==0 else -Kjet(v,Z,n)),breaks)
                Aj+=mp.exp(k*(3-t))*(mp.diff(gamma_A,Z,n)-(1/k if n==0 else 0))
                A.append(Aj+(1/k if n==0 else 0))
                Ej=mp.quad(lambda v:mp.exp(-delta*v)*squarejet(v,Z,n),breaks)
                Ej+=mp.exp(-3*delta)*mp.quad(lambda q:gamma_square(Z,q,delta,n),[0,1])
                E.append(mp.exp(delta*t)*Ej)
                Pj=mp.quad(lambda v:mp.exp(-p*v)*squarejet(v,Z,n)/2,breaks)
                Pj+=mp.exp(-3*p)*mp.quad(lambda q:gamma_square(Z,q,p,n),[0,1])/2
                P.append(mp.exp(p*t)*Pj)
            shape=heat.shape(c.mpf(Z),c.mpf(t)); tails=heat.collar_tails(c.mpf(Z),c.mpf(t))
            defects=collar_defect_rows(heat,shape,tails,c.mpf(t))
            rows=collar_stress_rows(heat,shape,defects,c.mpf(Z),c.mpf(t))
            for label,jet,values,baseline in (
                    ('angular',defects['angular_defect_rows'][0],A,1/k),
                    ('energy',defects['energy_defect_rows'][0],E,1/delta),
                    ('pressure',defects['pressure_defect_rows'][0],P,1/(2*p))):
                for n in range(3):
                    encloses(jet[n]*math.factorial(n),values[n]-(baseline if n==0 else 0),label+str(n))
                    future_checks+=1
            R=mp.exp(t)/S; B=mp.mpf('1.3')*mp.exp(-(mp.mpf('.5')+a)*t)
            qt=mp.sqrt(R/2)*B; qz=mp.sqrt(R/2)*B*B
            L=1-delta*Z*Z; d=1-Z*Z; Kj=Kjet(t,Z,0); Kz=Kjet(t,Z,1)
            # Direct original unnormalized Mtheta, Mztheta and absolute P.
            M=mp.sqrt(2)*R**mp.mpf('1.5')*B*A[0]; Mz=mp.sqrt(2)*R**mp.mpf('1.5')*B*A[1]
            Me=R*B*B*E[0]/2; Mez=R*B*B*E[1]/2; Pr=-B*B*P[0]; Prz=-B*B*P[1]
            U=B*Kj; Uy=B*(mp.diff(lambda v:Kjet(v,Z,0),t)-(mp.mpf('.5')+a)*Kj)
            Ttheta=(k*M-b*Z*Mz-R*mp.sqrt(2*R)*U)/(2*L*R)+mp.sqrt(2*R)/R*Uy-U/mp.sqrt(2*R)
            Tz=(2*delta*Z*Me-d*Mez+R*(2*p*Z*Pr-d*Prz))/(L*mp.sqrt(2*R))
            theta_numerator=k*A[0]-b*Z*A[1]-Kj
            theta_numerator_z=(k-b)*A[1]-b*Z*A[2]-Kz
            shear_z=2/R*(mp.diff(lambda v:Kjet(v,Z,1),t)-(1+a)*Kz)
            theta_z=theta_numerator_z/L+theta_numerator*(2*delta*Z)/L**2+shear_z
            axial_numerator=delta*Z*E[0]-d*E[1]/2-2*p*Z*P[0]+d*P[1]
            axial_numerator_z=delta*E[0]+(delta+1)*Z*E[1]-d*E[2]/2-2*p*P[0]-(2*p+2)*Z*P[1]+d*P[2]
            axial_z=axial_numerator_z/L+axial_numerator*(2*delta*Z)/L**2
            for label,reference,n in (('theta',Ttheta/qt,0),('axial',Tz/qz,0),
                                       ('theta',theta_z,1),('axial',axial_z,1)):
                encloses(rows[label][0][n]*math.factorial(n),reference,label+'_Z'+str(n)); stress_checks+=1
                if n==0 and abs(reference)>mp.mpf('1e-20'):nonzero+=1
            # Independently use paper (3.13) to differentiate the direct
            # inertial stress, rather than the producer's defect recurrences.
            Ky=mp.diff(lambda v:Kjet(v,Z,0),t)
            Kyy=mp.diff(lambda v:Kjet(v,Z,0),t,2)
            Kyyy=mp.diff(lambda v:Kjet(v,Z,0),t,3)
            Kzy=mp.diff(lambda v:Kjet(v,Z,1),t)
            It=theta_numerator/L; Iz=axial_numerator/L
            Nt=-(Ky+b*Z*Kz)/L; Nty=-(Kyy+b*Z*Kzy)/L
            Cs=2/R*(Ky-(1+a)*Kj)
            Csy=2/R*(Kyy-(2+a)*Ky+(1+a)*Kj)
            Csyy=2/R*(Kyyy-(3+a)*Kyy+(3+2*a)*Ky-(1+a)*Kj)
            Nz=(d*P[1]-2*Z*(p*P[0]-Kj*Kj/2))/L
            Nzy=(d*(p*P[1]-Kj*Kz)-2*Z*(p*p*P[0]-p*Kj*Kj/2-Kj*Ky))/L
            radial=(('theta',1,-It+Nt+Csy-a*Cs),
                    ('theta',2,It-(1+a)*Nt+Nty+Csyy-2*a*Csy+a*a*Cs),
                    ('axial',1,-Iz/2+Nz),
                    ('axial',2,Iz/4+Nzy-(1+delta)*Nz))
            for label,j,reference in radial:
                encloses(rows[label][j][0],reference,label+'_y'+str(j)); stress_checks+=1
        if nonzero!=4:raise ArithmeticError('Fixture did not exercise both nonzero collar stresses')
        return dict(full_future_moment_Z_derivative_checks=future_checks,
                    direct_original_moment_to_stress_checks=stress_checks,
                    nonzero_stress_components_exercised=nonzero,
                    finite_moderate_fixture_only=True, passed=True)


def run():
    name=PREFIX+'collar_stress_C3.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
    for source,digest in hashes.items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Collar stress dependency changed: '+source)
    proof=collar_moment_stress_identities()
    if proof!=record['collar_moment_stress_identities']:raise ValueError('Collar stress exact proof changed')
    pressure=CompliantCollarPressureC4(); history=terminal_history_bridge()
    bridge=collar_stress_source_bridge(pressure,history)
    if bridge!=record['collar_stress_source_bridge']:raise ValueError('Collar stress actual history changed')
    endpoint_proof=collar_Gamma_endpoint_binding()
    if endpoint_proof!=record['actual_pre_override_Gamma_endpoint_binding']:raise ValueError('Actual pre-override Gamma endpoint source changed')
    c=MPIntervalContext(); c.dps=270; finite=zeros=0
    points=record['samples']+[record['whole_Z_collar'],record['whole_Z_terminal']]
    for point in points:
        for flag in ('actual_full_collar_moments_same_source_verified','actual_original_collar_similarity_stress_recovered',
                     'original_collar_and_Gamma_stress_mixed3_join_verified','full_sigma_phi_Gamma_future_integral_used',
                     'exact_unit_baseline_cancellation_before_enclosure','original_angular_energy_pressure_histories_retained',
                     'collar_shear_strength_kappa_gt2_certified','collar_angular_shear_negative_certified'):
            if not point[flag]:raise ValueError('Collar source/stress gate missing: '+flag)
        for flag in ('collar_cone_certified','global_admissible_stress_lift_constructed','whole_outer_cone_certified',
                     'physical_energy_integral_certified','temporal_recursion'):
            if point[flag]:raise ValueError('Collar similarity-stress scope overclaimed: '+flag)
        if not point['actual_physical_map_and_remainder_transfer_pending']:raise ValueError('Physical remainder not yet built')
        if endpoints(read_interval(c,point['collar_shear_strength_kappa_minus2_Taylor']['coefficients'][0]))[0]<=0:
            raise ArithmeticError('Actual collar kappa>2 margin not positive')
        for label,grid in point['collar_similarity_stress_mixed3_factored'].items():
            if len(grid)!=10:raise ValueError('Stress mixed3 coverage missing')
            for key,value in grid.items():
                lo,hi=endpoints(read_interval(c,value))
                if not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite stress coefficient')
                finite+=1
                if point['Gamma_endpoint_stress_exact_zero_from_same_moments']:
                    if lo!=0 or hi!=0:raise ArithmeticError('Same-source Gamma stress join nonzero')
                    zeros+=1
        for jet in point['collar_meridional_moments_and_velocities'].values():
            for value in jet['coefficients']:
                if endpoints(read_interval(c,value))!=(0,0):raise ArithmeticError('Original zero meridional history lost')
    fixture=independent_full_moment_fixture()
    hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
                implicit_source_sha256=record['implicit_source_sha256'],all_passed=True,
                actual_original_collar_similarity_stress_recovered=True,
                original_collar_and_Gamma_stress_mixed3_join_verified=True,
                actual_finite_signed_stress_mixed_bounds_checked=finite,
                exact_Gamma_join_zero_derivatives_checked=zeros,
                whole_collar_shear_strength_kappa_gt2_certified=True,
                exact_original_moment_stress_identities=proof,
                actual_collar_stress_source_bridge=bridge,independent_full_moment_fixture=fixture,
                actual_pre_override_Gamma_endpoint_binding=endpoint_proof,
                physical_remainder_transfer_pending=True,collar_cone_certified=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual collar stress/history mixed3 and independent full-moment fixture PASS',flush=True)
    return result


if __name__=='__main__':run()
