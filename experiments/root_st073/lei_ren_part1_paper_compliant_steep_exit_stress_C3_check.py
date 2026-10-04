"""Steep-exit source/FTC, signed full moments and functional join check."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from scipy.integrate import quad
from lei_ren_part1_paper_compliant_steep_exit_stress_C3 import (
    CompliantSteepExitStressC3,exit_signed_kernels,exit_shape,exit_defect_rows,
    exit_pressure_rows,exit_transport_identities,exit_source_bridge,quotient_rows)
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def moderate_original_sigma_fixture():
    """Independent original-sigma quadrature and unnormalized stress jets.

    A convergent analytic future follows the constant waiting segment. Its
    complete integrals are known exactly. Original sigma quadrature supplies
    the finite exit interval only. Local primitives of scalar K Taylor jets
    differentiate the original full moment/stress equations through order
    four; they do not call the production defect or stress recurrences.
    This numerical fixture does not certify actual Gamma history or global NS.
    """
    with mp.workdps(65):
        a,eps,S,wait,Z=map(mp.mpf,('.15','.002','.004','1.3','.4'))
        delta=2*a; k=1-a; p=1+delta; bh=mp.mpf('.5')+a; b=(1-delta)/2; K0=1-eps
        c=MPIntervalContext(); c.dps=110
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),eps=c.mpf(eps),S=c.mpf(S),delta=c.mpf(delta),
                             k=c.mpf(k),prate=c.mpf(p),pressure_scale=c.mpf('1.7'),
                             steep=SimpleNamespace(wait=c.mpf(wait)))
        def g(z):return mp.mpf('.003')*(1+z*z)
        def future(label,z):
            if label=='A':return 1/k+eps/(2-k)-g(z)/(2-k)**2
            rate=delta if label=='E' else p; v=g(z)
            full=1/rate-2*eps/(rate+2)+2*v/(rate+2)**2+eps**2/(rate+4)-2*eps*v/(rate+4)**2+2*v*v/(rate+4)**3
            return full*(1 if label=='E' else mp.mpf('.5'))
        base={}
        for label,key,unit in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
            coefficients=[mp.diff(lambda z:future(label,z),Z,n)/math.factorial(n) for n in range(6)]
            coefficients[0]-=unit
            base[key+'_defect_rows']=[IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-55'),v+mp.mpf('1e-55')]) for v in coefficients])]
        def sigma(v):
            if v<=0:return mp.mpf(0)
            if v>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/v**2-1/(1-v)**2))
        def sigma_float(v):
            if v<=0:return 0.
            if v>=1:return 1.
            odds=1/(1-v)**2-1/v**2
            if odds>=0:return 1/(1+math.exp(-odds))
            ee=math.exp(odds); return ee/(1+ee)
        def f_float(v):
            return quad(lambda x:1-sigma_float(x),v,1,epsabs=2e-13,epsrel=2e-13)[0]
        def K_float(v):return float(K0)*math.exp(float(k)*f_float(v))
        errors=[]; integral_checks=stress_checks=pressure_checks=moment_checks=quotient_checks=0
        tolerance=mp.mpf('2e-7')
        def encloses(row,value,label):
            lo,hi=endpoints(row); error=max(lo-value,value-hi,mp.mpf(0)); errors.append(error)
            if error>tolerance:raise ArithmeticError('Steep exit independent fixture mismatch: '+label+' '+mp.nstr(error,12))
        for t0 in map(mp.mpf,('.25','.7')):
            tfloat=float(t0); length=1-t0; q0=-wait-1+t0
            kernels=exit_signed_kernels(c,c.mpf(t0),c.mpf(k),c.mpf(delta),c.mpf(eps))
            f0=mp.mpf(f_float(tfloat)); encloses(kernels['f'],f0,'remaining sigma primitive'); integral_checks+=1
            for name,rate,power in (('angular',-k,1),('energy',delta,2),('pressure',p,2)):
                value=quad(lambda v:math.exp(-float(rate)*(v-tfloat))*(K_float(v)**power-1),tfloat,1,
                           epsabs=2e-12,epsrel=2e-12)[0]
                encloses(kernels[name],mp.mpf(value),name+' signed integral'); integral_checks+=1
            # Direct positive energy/pressure integrals, independently of the
            # signed-defect cancellation used by the implementation.
            IA=mp.mpf(quad(lambda v:math.exp(float(k)*(v-tfloat))*(K_float(v)-1),tfloat,1,
                          epsabs=2e-12,epsrel=2e-12)[0])
            IE=mp.mpf(quad(lambda v:math.exp(-float(delta)*(v-tfloat))*K_float(v)**2,tfloat,1,
                          epsabs=2e-12,epsrel=2e-12)[0])
            IP=mp.mpf(quad(lambda v:math.exp(-float(p)*(v-tfloat))*K_float(v)**2/2,tfloat,1,
                          epsabs=2e-12,epsrel=2e-12)[0])
            def reference(label,z):
                if label=='A':
                    return 1/k-IA+eps*(mp.exp(k*(length+wait))-mp.exp(k*length))/k+mp.exp(k*(length+wait))*(future('A',z)-1/k)
                rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
                finite=IE if label=='E' else IP
                return finite+factor*K0*K0*mp.exp(-rate*length)*(-mp.expm1(-rate*wait))/rate+mp.exp(-rate*(length+wait))*future(label,z)
            # Taylor polynomial of the actual remaining sigma primitive.
            # Its first through fourth derivatives come from sigma itself,
            # rather than the production shape/ordinary-row helper.
            fcoeff=[f0,sigma(t0)-1]+[mp.diff(sigma,t0,n-1)/math.factorial(n) for n in range(2,5)]
            def localK(h):return K0*mp.exp(k*sum(fcoeff[n]*h**n for n in range(5)))
            densityA=mp.taylor(lambda h:mp.exp(k*h)*localK(h),0,4)
            densityE=mp.taylor(lambda h:mp.exp(-delta*h)*localK(h)**2,0,4)
            densityP=mp.taylor(lambda h:mp.exp(-p*h)*localK(h)**2/2,0,4)
            def primitive(coeff,h):return sum(v*h**(n+1)/(n+1) for n,v in enumerate(coeff))
            def moment(qq,z,label):
                h=qq-q0
                if label=='A':return mp.exp(-k*h)*(reference(label,z)+primitive(densityA,h))
                rate=delta if label=='E' else p; density=densityE if label=='E' else densityP
                return mp.exp(rate*h)*(reference(label,z)-primitive(density,h))
            def fullstress(qq,z,label):
                R=mp.exp(qq)/S; B=mp.mpf('1.3')*mp.exp(-bh*qq)
                A=moment(qq,z,'A'); E=moment(qq,z,'E'); Pr=moment(qq,z,'P'); K=localK(qq-q0)
                Mt=mp.sqrt(2)*R**mp.mpf('1.5')*B*A
                Mtz=mp.sqrt(2)*R**mp.mpf('1.5')*B*mp.diff(lambda zz:moment(qq,zz,'A'),z)
                Me=R*B*B*E/2; Mez=R*B*B/2*mp.diff(lambda zz:moment(qq,zz,'E'),z)
                pressure=-B*B*Pr; pressurez=-B*B*mp.diff(lambda zz:moment(qq,zz,'P'),z)
                L=1-delta*z*z; d=1-z*z; U=B*K
                if label=='theta':
                    inertial=(k*Mt-b*z*Mtz-R*mp.sqrt(2*R)*U)/(2*L*R)
                    Uy=mp.diff(lambda yy:mp.mpf('1.3')*mp.exp(-bh*yy)*localK(yy-q0),qq)
                    shear=mp.sqrt(2*R)/R*Uy-U/mp.sqrt(2*R)
                    return inertial+shear
                return (2*delta*z*Me-d*Mez+R*(2*p*z*pressure-d*pressurez))/(L*mp.sqrt(2*R))
            shape=exit_shape(heat,c.mpf(t0),kernels); defects=exit_defect_rows(heat,base,c.mpf(t0),kernels,shape)
            rows=collar_stress_rows(heat,shape,defects,c.mpf(Z),c.mpf(q0))
            pressure=exit_pressure_rows(heat,defects,shape,c.mpf(q0))
            Arows=[IntervalTaylor.constant(c,1/k,5)+defects['angular_defect_rows'][0]]+defects['angular_defect_rows'][1:]
            Erows=[IntervalTaylor.constant(c,1/delta,5)+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
            xr=quotient_rows(Arows,shape['K_rows']); er=quotient_rows(Erows,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
            for label,factor in (('theta',mp.sqrt(mp.exp(q0)/S/2)*mp.mpf('1.3')*mp.exp(-bh*q0)),
                                 ('axial',mp.sqrt(mp.exp(q0)/S/2)*mp.mpf('1.3')**2*mp.exp(-2*bh*q0))):
                for j in range(4):
                    for n in range(4-j):
                        value=mp.diff(lambda zz:mp.diff(lambda qq:fullstress(qq,zz,label),q0,j),Z,n)/factor
                        encloses(rows[label][j][n]*math.factorial(n),value,label+str((j,n))); stress_checks+=1
            for j in range(5):
                for n in range(5-j):
                    value=mp.diff(lambda zz:mp.diff(lambda qq:-mp.mpf('1.7')*mp.exp(-p*qq)*moment(qq,zz,'P'),q0,j),Z,n)
                    encloses(pressure[j][n]*math.factorial(n),value,'absolute pressure'+str((j,n))); pressure_checks+=1
                for label,key,unit in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
                    value=mp.diff(lambda qq:moment(qq,Z,label),q0,j)-(unit if j==0 else 0)
                    encloses(defects[key+'_defect_rows'][j][0],value,key+str(j)); moment_checks+=1
                for label,grid in (('A',xr),('E',er)):
                    value=mp.diff(lambda qq:moment(qq,Z,label)/(localK(qq-q0) if label=='A' else 2*localK(qq-q0)**2),q0,j)
                    encloses(grid[j][0],value,'original normalized quotient'+label+str(j)); quotient_checks+=1
        return dict(original_sigma_signed_integral_checks=integral_checks,original_unnormalized_stress_mixed3_checks=stress_checks,
                    absolute_pressure_mixed4_checks=pressure_checks,complete_future_moment_derivative_checks=moment_checks,
                    original_X_and_half_energy_quotient_derivative_checks=quotient_checks,
                    max_unenclosed_numeric_error=mp.nstr(max(errors),12),numeric_comparison_tolerance=str(tolerance),
                    finite_original_sigma_quadrature_is_numeric_diagnostic=True,
                    analytic_convergent_future_fixture_only=True,actual_Gamma_terminal_history_proved_by_fixture=False,
                    actual_source_history_consumed_separately=True,passed=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'steep_exit_stress_C3.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Steep exit source changed: '+path)
        provider=CompliantSteepExitStressC3(); c=provider.ctx
        if exit_transport_identities()!=record['transport_identities']:raise ValueError('Exit FTC proof changed')
        if exit_source_bridge(provider.waiting_source)!=record['steep_exit_source_bridge']:raise ValueError('Exit source binding changed')
        finite=pressure_rows=zeros=flat=0
        for point in record['samples']+[record['whole_steep_exit']]:
            for grid in point['steep_exit_similarity_stress_mixed3_factored'].values():
                if len(grid)!=10:raise ValueError('Exit stress mixed3 incomplete')
                for value in grid.values():
                    if any(not mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite exit stress')
                    finite+=1
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                if len(jet['coefficients'])!=6:raise ValueError('Exit pressure axial5 incomplete')
                for n in range(5-j):
                    if any(not mp.isfinite(v) for v in endpoints(read_interval(c,jet['coefficients'][n]))):raise ArithmeticError('Nonfinite exit pressure')
                    pressure_rows+=1
            for jet in point['steep_exit_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(0,0):raise ArithmeticError('Exit meridional history lost')
                    zeros+=1
            if endpoints(read_interval(c,point['steep_exit_shear_strength_kappa_minus2']))[0]<=0:raise ArithmeticError('Exit source shear margin lost')
        kernels=provider.kernels(c.mpf(1)); shape=exit_shape(provider.heat,c.mpf(1),kernels)
        for jet in shape['K_rows'][1:]:
            for value in jet.coefficients:
                if endpoints(value)!=(0,0):raise ArithmeticError('Actual exit K endpoint not flat')
                flat+=1
        for flag in ('steep_exit_cone_certified','steep_exit_physical_axial_viscosity_exact_zero',
                     'global_admissible_stress_lift_constructed','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Steep exit scope overclaimed: '+flag)
        fixture=moderate_original_sigma_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
                    all_passed=True,actual_original_steep_exit_similarity_stress_recovered=True,
                    actual_steep_exit_absolute_pressure_same_source_mixed4_available=True,
                    steep_exit_waiting_stress_mixed3_join_verified=True,steep_exit_waiting_pressure_mixed4_join_verified=True,
                    actual_finite_signed_stress_rows=finite,actual_finite_pressure_rows=pressure_rows,
                    exact_meridional_history_zero_rows=zeros,actual_flat_endpoint_K_derivative_checks=flat,
                    independent_original_sigma_full_future_fixture=fixture,steep_exit_cone_certified=False,
                    steep_exit_physical_axial_viscosity_exact_zero=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS actual steep-exit full-future stress/pressure, exact waiting joins and independent original-sigma fixture',flush=True)
        return result


if __name__=='__main__':run()
