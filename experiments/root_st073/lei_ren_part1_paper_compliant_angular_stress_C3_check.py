"""Original angular full-moment recovery, source joins and independent beta fixture."""
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_angular_stress_C3 import (
    CompliantAngularStressC3,angular_future_changes,angular_shape,angular_X_rows,
    angular_defect_rows,angular_stress_rows,angular_pressure_rows)
from lei_ren_part1_paper_compliant_outer_angular_repair import bump_weights
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_original_beta_fixture():
    """Direct cumulative/complete-future integrals and unnormalized stresses.

    Both original compact bumps have nonconstant axial polynomial coefficients.
    Beyond s=0 the same baseline exponential gives convergent quadratic tails.
    Reference functions use independent mpmath quadrature and differentiation,
    not production moment/stress recurrences. This is not global NS validation.
    """
    with mp.workdps(65):
        a,mu,KR,S,C,v0,z0=map(mp.mpf,('.1','.3','2','.004','.08','-1.02','.2'))
        delta=2*a; k=1-a; p=1+delta; r=1-mu; d=a-mu; po=1+2*mu
        bh=mp.mpf('.5')+a; b=(1-delta)/2; qR=mp.mpf('-2.3'); ell=mp.mpf('.15')
        def raw(u):return mp.exp(-1/(1-u*u)) if abs(u)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1])
        def beta(x):return raw(x/ell)/(ell*norm)
        def beta1(x):
            u=x/ell
            return beta(x)*(-2*u)/(ell*(1-u*u)**2) if abs(u)<1 else mp.mpf(0)
        def coeff(z):return (mp.mpf('.025')*(1+mp.mpf('.2')*z+mp.mpf('.05')*z*z),
                             -mp.mpf('.018')*(1-mp.mpf('.1')*z+mp.mpf('.03')*z*z))
        def coeff_z(z):return (mp.mpf('.025')*(mp.mpf('.2')+mp.mpf('.1')*z),
                               -mp.mpf('.018')*(-mp.mpf('.1')+mp.mpf('.06')*z))
        def weight(local,rate,power,future):
            return weight_integral(max(-ell,min(ell,local)),rate,power,future,170)
        @lru_cache(maxsize=512)
        def weight_integral(stop,rate,power,future,precision_digits):
            lo,hi=(stop,ell) if future else (-ell,stop)
            if lo==hi:return mp.mpf(0)
            with mp.workdps(precision_digits):
                return mp.quad(lambda x:mp.exp(rate*x)*beta(x)**power,[lo,(lo+hi)/2,hi])
        def reference(v,z):
            cs=coeff(z); cz=coeff_z(z); G=KR*mp.exp(d*v)
            F=mp.mpf(1); Fv=mp.mpf(0); N=1/r; Nz=mp.mpf(0)
            Ec=mp.mpf(0); Ez=mp.mpf(0); Pc=mp.mpf(0); Pz=mp.mpf(0)
            for dj,dz,center in zip(cs,cz,(-3,-1)):
                local=v-center; F+=dj*beta(local); Fv+=dj*beta1(local)
                wa=weight(local,r,1,False)*mp.exp(-r*(v-center))
                N+=dj*wa; Nz+=dz*wa
                ie=weight(local,-2*mu,1,True); iff=weight(local,-2*mu,2,True)
                ib=weight(local,-po,1,True); idd=weight(local,-po,2,True)
                Ec+=(2*dj*ie+dj*dj*iff)*mp.exp(-2*mu*center)
                Ez+=(2*dz*ie+2*dj*dz*iff)*mp.exp(-2*mu*center)
                Pc+=(dj*ib+dj*dj*idd/2)*mp.exp(-po*center)
                Pz+=(dz*ib+dj*dz*idd)*mp.exp(-po*center)
            K=G*F; Kv=G*(d*F+Fv); A=G*N; Az=G*Nz
            E=G*G/(2*mu)+KR*KR*mp.exp(delta*v)*Ec
            E_z=KR*KR*mp.exp(delta*v)*Ez
            Pr=G*G/(2*po)+KR*KR*mp.exp(p*v)*Pc
            Pr_z=KR*KR*mp.exp(p*v)*Pz
            L=1-delta*z*z; q=qR+v
            Ct=(k*A-b*z*Az-K)/L+2*S*mp.exp(-q)*(Kv-(1+a)*K)
            Cz=(delta*z*E-(1-z*z)*E_z/2-2*p*z*Pr+(1-z*z)*Pr_z)/L
            B=mp.sqrt(C)*mp.exp(-bh*q); qt=mp.sqrt(mp.exp(q)/(2*S))*B
            return dict(X=N/F,theta=qt*Ct,axial=qt*B*Cz,pressure=-C*mp.exp(-p*q)*Pr)
        c=MPIntervalContext(); c.dps=110
        h=SimpleNamespace(ctx=c,a=c.mpf(a),mu=c.mpf(mu),delta=c.mpf(delta),k=c.mpf(k),
                          prate=c.mpf(p),S=c.mpf(S),pressure_scale=c.mpf(C))
        z=IntervalTaylor.variable(c,c.mpf(z0),5); one=IntervalTaylor.constant(c,1,5)
        cs=[(one+z*c.mpf('.2')+z*z*c.mpf('.05'))*c.mpf('.025'),
            (one-z*c.mpf('.1')+z*z*c.mpf('.03'))*c.mpf('-.018')]
        normbox=c.mpf([norm-mp.mpf('1e-55'),norm+mp.mpf('1e-55')])
        weights=bump_weights(c,h.mu,normbox,cells=512)
        outer=SimpleNamespace(ctx=c,mu=h.mu,prate=1+2*h.mu,normalization=normbox,weights=weights,cells=512)
        H=[one*0 for _ in range(5)]
        for dj,center in zip(cs,(-3,-1)):
            jets=beta_jets(c,c.mpf(v0-center)/c.mpf('.15'))
            for j in range(5):H[j]+=dj*(jets[j]*math.factorial(j)/(c.mpf('.15')**(j+1)*normbox))
        packet=dict(swirl_factor_one_plus_h_Taylor=H[0]+one,actual_angular_bump_y_derivatives=H)
        terminal=dict(KR=c.mpf(KR),energy_defect_rows=[one*(c.mpf(KR)**2/(2*h.mu)-1/h.delta)],
                      pressure_defect_rows=[one*(c.mpf(KR)**2/(2*(1+2*h.mu))-1/(2*h.prate))])
        X=IntervalTaylor(c,[c.mpf([mp.diff(lambda zz:reference(v0,zz)['X'],z0,n)/math.factorial(n)-mp.mpf('1e-50'),
                                  mp.diff(lambda zz:reference(v0,zz)['X'],z0,n)/math.factorial(n)+mp.mpf('1e-50')]) for n in range(6)])
        shape=angular_shape(h,terminal,c.mpf(v0),packet)
        xr=angular_X_rows(h,shape,X); kernels=angular_future_changes(outer,dict(coeff=cs),c.mpf(v0))
        defects=angular_defect_rows(h,terminal,c.mpf(v0),shape,xr[0],kernels)
        q=c.mpf(qR+v0); stress=angular_stress_rows(h,shape,defects,c.mpf(z0),q)
        pressure=angular_pressure_rows(h,defects,shape,q)
        pivot=reference(v0,z0); B=mp.sqrt(C)*mp.exp(-bh*(qR+v0)); qt=mp.sqrt(mp.exp(qR+v0)/(2*S))*B
        counts=dict(stress_mixed3=0,pressure_mixed4=0); maxmiss=mp.mpf(0); tolerance=mp.mpf('1e-35')
        step=mp.mpf('1e-26'); max_step_change=mp.mpf(0)
        def independent_derivative(label,j,n,factor=1):
            nonlocal max_step_change
            fn=lambda vv,zz:reference(vv,zz)[label]
            value=mp.diff(fn,(v0,z0),(j,n),h=step)/factor
            doubled=mp.diff(fn,(v0,z0),(j,n),h=step*2)/factor
            change=abs(value-doubled); max_step_change=max(max_step_change,change)
            if change>tolerance:raise ArithmeticError('Angular reference step stability failed: '+str((label,j,n,change)))
            return value
        def compare(value,box,label):
            nonlocal maxmiss
            lo,hi=endpoints(box); miss=max(lo-value,value-hi,mp.mpf(0)); maxmiss=max(maxmiss,miss)
            if miss>tolerance:
                print('Angular fixture mismatch',label,'value',mp.nstr(value,18),'box',(mp.nstr(lo,18),mp.nstr(hi,18)),'miss',mp.nstr(miss,18),flush=True)
                raise ArithmeticError('Independent angular original beta fixture failed: '+label)
        for label,factor in (('theta',qt),('axial',qt*B)):
            for j in range(4):
                for n in range(4-j):
                    value=independent_derivative(label,j,n,factor)
                    compare(value,stress[label][j][n]*math.factorial(n),label+str((j,n))); counts['stress_mixed3']+=1
        for j in range(5):
            for n in range(5-j):
                value=independent_derivative('pressure',j,n)
                compare(value,pressure[j][n]*math.factorial(n),'pressure'+str((j,n))); counts['pressure_mixed4']+=1
        return dict(all_passed=True,fixture='Original compact beta, two nonconstant axial coefficients, convergent quadratic future',
                    comparison_counts=counts,maximum_positive_enclosure_miss=str(maxmiss),comparison_tolerance=str(tolerance),
                    reference_quadrature_digits=170,reference_derivative_step=str(step),
                    maximum_doubled_step_change=str(max_step_change),reference_step_stability_comparisons=35,
                    actual_source_history_admission_is_separate=True,fixture_is_global_NS_validation=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'angular_stress_C3.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for filename,digest in hashes.items():
            if hashlib.sha256((HERE/filename).read_bytes()).hexdigest()!=digest:raise ValueError('Angular stress source changed: '+filename)
        fresh=CompliantAngularStressC3().report()
        if encode(pack(fresh))!=record:raise ValueError('Current original angular source/moments/joins differ')
        c=MPIntervalContext(); c.dps=300; signed=0; pressure=0; zeros=0
        for point in record['samples']+[record['whole_angular']]+record['support_crossings']:
            for label,grid in point['angular_similarity_stress_mixed3_factored'].items():
                for value in grid.values():
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite angular stress row')
                    signed+=1
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                for value in jet['coefficients'][:5-j]:
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite angular pressure row')
                # The checked mixed4 grid uses only j+n<=4.
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                pressure+=len(jet['coefficients'][:5-j])
            for jet in point['angular_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Actual angular meridional zero lost')
                    zeros+=1
        for flag in ('actual_original_angular_similarity_stress_recovered','actual_angular_absolute_pressure_same_source_mixed4_available',
                     'angular_entry_stress_mixed3_and_pressure_mixed4_join_verified','actual_angular_K_Z_dependence_retained',
                     'original_selected_angular_coefficient_functions_retained'):
            if not record[flag]:raise ValueError('Angular stress admission missing: '+flag)
        for flag in ('angular_physical_decomposition_constructed','angular_cone_certified','preceding_power_stress_companion_constructed',
                     'global_admissible_stress_lift_constructed','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Angular stress scope overclaimed: '+flag)
        print('Current angular source, signed rows and arbitrary-data right joins admitted; independent beta fixture running',flush=True)
        fixture=independent_original_beta_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=fresh['actual_five_defect_family_sha256'],
                    implicit_source_sha256=fresh['implicit_source_sha256'],
                    finite_signed_stress_mixed_rows=signed,finite_pressure_mixed_rows=pressure,exact_meridional_zero_coefficients=zeros,
                    actual_original_angular_similarity_stress_recovered=True,actual_angular_absolute_pressure_same_source_mixed4_available=True,
                    angular_entry_stress_mixed3_and_pressure_mixed4_join_verified=True,
                    actual_angular_K_Z_dependence_retained=True,original_selected_angular_coefficient_functions_retained=True,
                    actual_full_native_history_and_source_AST_bridge_verified=True,
                    actual_angular_entry_formula_join_identities=len(fresh['source_bridge']['actual_angular_entry_formula_join']['identities']),
                    independent_original_beta_fixture=fixture,
                    angular_physical_decomposition_constructed=False,angular_cone_certified=False,
                    preceding_power_stress_companion_constructed=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original angular full moments, axial-dependent stress, absolute pressure and right join; physical/cone pending',flush=True)
        return result


if __name__=='__main__':run()
