"""Original steep-entry source, full-future stress and functional joins."""
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from scipy.integrate import quad
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import (
    CompliantSteepEntryStressC3,entry_signed_kernels,entry_shape,entry_angular_rows,
    entry_defect_rows,entry_stress_rows,entry_pressure_rows,quotient_rows)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def moderate_original_sigma_fixture():
    """Independent finite quadrature, complete analytic future and native stress.

    After t=1 choose K(v)=1+(KS-1+g(Z)*v^2)*exp(-rho*v).
    Its complete angular-defect/energy/pressure integrals are exact. Local
    primitives of independently differentiated sigma/K densities provide
    ordinary moment and unnormalized stress derivatives. Production defect,
    quotient and stress recurrences never define the reference functions.
    This is separate from actual Gamma-history admission and global NS error.
    """
    with mp.workdps(65):
        a,mu,KS,rho,S,wait,Z=map(mp.mpf,('.15','.2','1.4','1.05','.004','1.3','.4'))
        delta=2*a; k=1-a; p=1+delta; r=1-mu; bh=mp.mpf('.5')+a; b=(1-delta)/2
        Ts=4*mp.log(2/delta); C=mp.mpf('3e-8'); D=KS-1
        c=MPIntervalContext(); c.dps=110
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),mu=c.mpf(mu),delta=c.mpf(delta),k=c.mpf(k),
            prate=c.mpf(p),S=c.mpf(S),pressure_scale=c.mpf(C),steep=SimpleNamespace(wait=c.mpf(wait),Ts=c.mpf(Ts)))
        def g(z):return mp.mpf('.003')*(1+z*z)
        def future(label,z):
            v=g(z)
            if label=='A':return 1/k-D/(rho-k)-2*v/(rho-k)**3
            rate=delta if label=='E' else p
            full=1/rate+2*D/(rate+rho)+4*v/(rate+rho)**3+D*D/(rate+2*rho)
            full+=4*D*v/(rate+2*rho)**3+24*v*v/(rate+2*rho)**5
            return full*(1 if label=='E' else mp.mpf('.5'))
        terminal={'KS':c.mpf(KS)}
        def jet_at(fn):
            vals=[mp.diff(fn,Z,n)/math.factorial(n) for n in range(6)]
            return IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-55'),v+mp.mpf('1e-55')]) for v in vals])
        for label,key,unit in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
            terminal[key+'_defect_rows']=[jet_at(lambda z:future(label,z)-unit)]
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
        @lru_cache(maxsize=None)
        def f_float(v):return quad(sigma_float,v,1,epsabs=2e-13,epsrel=2e-13)[0]
        def K_float(v):return float(KS)*math.exp(float(mu-a)*(1-v)+float(r)*f_float(v))
        errors=[]; integral_checks=stress_checks=pressure_checks=moment_checks=quotient_checks=0
        tolerance=mp.mpf('2e-7')
        def encloses(row,value,label):
            lo,hi=endpoints(row); error=max(lo-value,value-hi,mp.mpf(0)); errors.append(error)
            if error>tolerance:raise ArithmeticError('Independent original entry fixture mismatch: '+label+' '+mp.nstr(error,12))
        for t0 in map(mp.mpf,('.25','.7')):
            tf=float(t0); ell=1-t0; q0=-wait-Ts-2+t0
            kernels=entry_signed_kernels(c,c.mpf(t0),c.mpf(KS),c.mpf(mu),c.mpf(delta))
            f0=mp.mpf(f_float(tf)); encloses(kernels['f'],f0,'remaining original sigma'); integral_checks+=1
            for key,rate,power in (('angular',-k,1),('energy',delta,2),('pressure',p,2)):
                value=quad(lambda v:math.exp(-float(rate)*(v-tf))*(K_float(v)**power-1),tf,1,
                           epsabs=2e-12,epsrel=2e-12)[0]
                encloses(kernels[key],mp.mpf(value),key+' signed remaining integral'); integral_checks+=1
            IA=mp.mpf(quad(lambda v:math.exp(float(k)*(v-tf))*K_float(v),tf,1,epsabs=2e-12,epsrel=2e-12)[0])
            IE=mp.mpf(quad(lambda v:math.exp(-float(delta)*(v-tf))*K_float(v)**2,tf,1,epsabs=2e-12,epsrel=2e-12)[0])
            IP=mp.mpf(quad(lambda v:math.exp(-float(p)*(v-tf))*K_float(v)**2/2,tf,1,epsabs=2e-12,epsrel=2e-12)[0])
            def reference(label,z):
                if label=='A':return mp.exp(k*ell)*future(label,z)-IA
                rate=delta if label=='E' else p
                return mp.exp(-rate*ell)*future(label,z)+(IE if label=='E' else IP)
            fcoeff=[f0,-sigma(t0)]+[-mp.diff(sigma,t0,n-1)/math.factorial(n) for n in range(2,5)]
            def localK(h):return KS*mp.exp((mu-a)*(ell-h)+r*sum(fcoeff[n]*h**n for n in range(5)))
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
            X0=jet_at(lambda z:reference('A',z)/localK(0))
            shape=entry_shape(heat,terminal,c.mpf(t0),kernels)
            xr=entry_angular_rows(heat,X0,shape)
            defects=entry_defect_rows(heat,terminal,c.mpf(t0),kernels,shape,X0)
            rows=entry_stress_rows(heat,shape,defects,xr,c.mpf(Z),c.mpf(q0))
            pressure=entry_pressure_rows(heat,defects,shape,c.mpf(q0))
            Erows=[IntervalTaylor.constant(c,1/delta,5)+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
            er=quotient_rows(Erows,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
            for label,factor in (('theta',mp.sqrt(mp.exp(q0)/S/2)*mp.mpf('1.3')*mp.exp(-bh*q0)),
                                 ('axial',mp.sqrt(mp.exp(q0)/S/2)*mp.mpf('1.3')**2*mp.exp(-2*bh*q0))):
                for j in range(4):
                    for n in range(4-j):
                        value=mp.diff(lambda zz:mp.diff(lambda qq:fullstress(qq,zz,label),q0,j),Z,n)/factor
                        encloses(rows[label][j][n]*math.factorial(n),value,label+str((j,n))); stress_checks+=1
            for j in range(5):
                for n in range(5-j):
                    value=mp.diff(lambda zz:mp.diff(lambda qq:-C*mp.exp(-p*qq)*moment(qq,zz,'P'),q0,j),Z,n)
                    encloses(pressure[j][n]*math.factorial(n),value,'absolute pressure'+str((j,n))); pressure_checks+=1
                for label,key,unit in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
                    value=mp.diff(lambda qq:moment(qq,Z,label),q0,j)-(unit if j==0 else 0)
                    encloses(defects[key+'_defect_rows'][j][0],value,key+str(j)); moment_checks+=1
                for label,grid in (('A',xr),('E',er)):
                    value=mp.diff(lambda qq:moment(qq,Z,label)/(localK(qq-q0) if label=='A' else 2*localK(qq-q0)**2),q0,j)
                    encloses(grid[j][0],value,'native quotient'+label+str(j)); quotient_checks+=1
        return dict(original_sigma_signed_integral_checks=integral_checks,original_unnormalized_stress_mixed3_checks=stress_checks,
            absolute_pressure_mixed4_checks=pressure_checks,complete_future_moment_derivative_checks=moment_checks,
            original_X_and_half_energy_quotient_derivative_checks=quotient_checks,
            max_unenclosed_numeric_error=mp.nstr(max(errors),12),numeric_comparison_tolerance=str(tolerance),
            numerical_finite_quadrature_diagnostic_not_actual_Gamma_proof=True,
            complete_convergent_analytic_future_fixture_only=True,actual_source_history_consumed_separately=True,passed=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'steep_entry_stress_C3.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Entry source changed: '+path)
        provider=CompliantSteepEntryStressC3(); c=provider.ctx
        if provider.proof!=record['transport_identities'] or provider.bridge!=record['steep_entry_source_bridge']:
            raise ValueError('Original entry source/FTC proof changed')
        finite=pressure_rows=zeros=0
        for point in record['samples']+[record['whole_steep_entry']]:
            for grid in point['steep_entry_similarity_stress_mixed3_factored'].values():
                if len(grid)!=10:raise ValueError('Entry stress mixed3 incomplete')
                for value in grid.values():
                    if any(not mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite entry stress')
                    finite+=1
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                if len(jet['coefficients'])!=6:raise ValueError('Entry pressure axial5 incomplete')
                for n in range(5-j):
                    if any(not mp.isfinite(v) for v in endpoints(read_interval(c,jet['coefficients'][n]))):raise ArithmeticError('Nonfinite entry pressure')
                    pressure_rows+=1
            for jet in point['steep_entry_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(0,0):raise ArithmeticError('Entry meridional history lost')
                    zeros+=1
            if endpoints(read_interval(c,point['steep_entry_shear_strength_kappa_minus2']))[0]<=0:
                raise ArithmeticError('Entry source shear margin lost')
        for flag in ('steep_entry_physical_decomposition_constructed','steep_entry_cone_certified',
            'upstream_angular_stress_companion_constructed','global_admissible_stress_lift_constructed',
            'physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Entry scope overclaimed: '+flag)
        fixture=moderate_original_sigma_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            all_passed=True,actual_original_steep_entry_similarity_stress_recovered=True,
            actual_steep_entry_absolute_pressure_same_source_mixed4_available=True,
            steep_entry_power_stress_mixed3_join_verified=True,steep_entry_power_pressure_mixed4_join_verified=True,
            native_angular_entry_left_field_pressure_join_verified=True,actual_finite_signed_stress_rows=finite,
            actual_finite_pressure_rows=pressure_rows,exact_meridional_history_zero_rows=zeros,
            actual_entry_power_formula_join_identity_count=len(provider.bridge['actual_entry_power_formula_join']['identities']),
            independent_original_sigma_complete_future_fixture=fixture,
            steep_entry_physical_decomposition_constructed=False,steep_entry_cone_certified=False,
            upstream_angular_stress_companion_constructed=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original steep-entry same full-future stress/pressure, functional power joins and independent sigma fixture',flush=True)
        return result


if __name__=='__main__':run()
