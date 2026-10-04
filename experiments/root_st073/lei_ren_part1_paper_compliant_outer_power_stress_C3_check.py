"""Original outer-power full moments, ordinary-q source join and bounded fixture."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_outer_power_stress_C3 import (
    CompliantOuterPowerStressC3,outer_power_shape,outer_power_defect_rows,
    outer_power_stress_rows,outer_power_pressure_rows)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_piecewise_exponential_fixture():
    """Independent past/future integrals with axial-dependent boundary moments.

    K is independent of Z on the tested power chart. Earlier and later profiles
    depend on Z and leave nonzero A_Z/E_Z/P_Z here. Moderate parameters test the
    transfer and stress/pressure operators, not original global NS accuracy.
    """
    with mp.workdps(85):
        a,mu,KR,S,C,y0,z0,L,h=map(mp.mpf,('.1','.3','2','.004','.08','-1.1','.2','2','.45'))
        delta=2*a; k=1-a; p=1+delta; d=a-mu; r=1-mu
        bh=mp.mpf('.5')+a; b=(1-delta)/2; qR=mp.mpf('-2.3'); past_rate=mp.mpf('.6')
        Kleft=KR*mp.exp(-d*L)
        def g(z):return mp.mpf('.08')*(1+mp.mpf('.2')*z+mp.mpf('.05')*z*z)
        def ga(z):return mp.mpf('.03')*(1-mp.mpf('.1')*z+mp.mpf('.04')*z*z)
        def right_future(rate,z):
            t=rate+2*h
            return KR**2*(1/t+2*g(z)/(t+1)**2+2*g(z)**2/(t+2)**3)
        def left_past(z):
            return Kleft*(1/(k+past_rate)-ga(z)/(k+past_rate+1)**2)
        def moments(y,z):
            Ky=KR*mp.exp(d*y)
            A=left_past(z)*mp.exp(-k*(y+L))+Ky*(-mp.expm1(-r*(y+L)))/r
            E=mp.exp(delta*y)*right_future(delta,z)+Ky**2*(-mp.expm1(2*mu*y))/(2*mu)
            Pr=mp.exp(p*y)*right_future(p,z)/2+Ky**2*(-mp.expm1((1+2*mu)*y))/(2*(1+2*mu))
            return A,E,Pr
        def actual_K(v,z):
            if v<-L:
                w=v+L
                return Kleft*mp.exp(past_rate*w)*(1+ga(z)*w*mp.exp(w))
            if v<=0:return KR*mp.exp(d*v)
            return KR*mp.exp(-h*v)*(1+g(z)*v*mp.exp(-v))
        integral_count=0; integral_miss=mp.mpf(0)
        for yy in (-L,y0,mp.mpf(0)):
            A,E,Pr=moments(yy,z0)
            directA=mp.quad(lambda v:actual_K(v,z0)*mp.exp(k*(v-yy)),[-mp.inf,-L,yy]) if yy>-L else mp.quad(lambda v:actual_K(v,z0)*mp.exp(k*(v-yy)),[-mp.inf,yy])
            directE=mp.quad(lambda v:actual_K(v,z0)**2*mp.exp(-delta*(v-yy)),[yy,0,mp.inf]) if yy<0 else mp.quad(lambda v:actual_K(v,z0)**2*mp.exp(-delta*v),[0,mp.inf])
            directP=mp.quad(lambda v:actual_K(v,z0)**2*mp.exp(-p*(v-yy))/2,[yy,0,mp.inf]) if yy<0 else mp.quad(lambda v:actual_K(v,z0)**2*mp.exp(-p*v)/2,[0,mp.inf])
            for value,direct in zip((A,E,Pr),(directA,directE,directP)):
                miss=abs(value-direct);integral_miss=max(integral_miss,miss);integral_count+=1
                if miss>mp.mpf('1e-65'):raise ArithmeticError('Independent whole-history integral mismatch')
        def reference(y,z):
            Ky=KR*mp.exp(d*y);A,E,Pr=moments(y,z)
            Az=mp.exp(-k*(y+L))*mp.diff(left_past,z)
            Ez=mp.exp(delta*y)*mp.diff(lambda zz:right_future(delta,zz),z)
            Pz=mp.exp(p*y)*mp.diff(lambda zz:right_future(p,zz)/2,z)
            den=1-delta*z*z;q=qR+y
            Ct=(k*A-b*z*Az-Ky)/den+2*S*mp.exp(-q)*(d-(1+a))*Ky
            Cz=(delta*z*E-(1-z*z)*Ez/2-2*p*z*Pr+(1-z*z)*Pz)/den
            B=mp.sqrt(C)*mp.exp(-bh*q);qt=mp.sqrt(mp.exp(q)/(2*S))*B
            return dict(theta=qt*Ct,axial=qt*B*Cz,pressure=-C*mp.exp(-p*q)*Pr)
        c=MPIntervalContext();c.dps=110
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),mu=c.mpf(mu),delta=c.mpf(delta),k=c.mpf(k),
            prate=c.mpf(p),S=c.mpf(S),pressure_scale=c.mpf(C))
        z=IntervalTaylor.variable(c,c.mpf(z0),5);one=IntervalTaylor.constant(c,1,5)
        gj=(one+z*c.mpf('.2')+z*z*c.mpf('.05'))*c.mpf('.08')
        gaj=(one-z*c.mpf('.1')+z*z*c.mpf('.04'))*c.mpf('.03')
        def future_jet(rate):
            t=rate+2*c.mpf(h)
            return (one/t+gj*(2/(t+1)**2)+gj*gj*(2/(t+2)**3))*c.mpf(KR)**2
        terminal=dict(Kright=c.mpf(KR),energy_defect_rows=[future_jet(heat.delta)-one/heat.delta],
            pressure_defect_rows=[future_jet(heat.prate)/2-one/(2*heat.prate)])
        Aleft=(one/(heat.k+c.mpf(past_rate))-gaj/(heat.k+c.mpf(past_rate)+1)**2)*c.mpf(Kleft)
        Ky=c.mpf(KR)*c.exp((heat.a-heat.mu)*c.mpf(y0))
        X=(Aleft*c.exp(-heat.k*c.mpf(y0+L))+one*(Ky*(-c.expm1(-c.mpf(r)*c.mpf(y0+L)))/c.mpf(r)))/Ky
        shape=outer_power_shape(heat,terminal,c.mpf(y0))
        defects=outer_power_defect_rows(heat,terminal,c.mpf(y0),shape,X)
        stress=outer_power_stress_rows(heat,shape,defects,c.mpf(z0),c.mpf(qR+y0))
        pressure=outer_power_pressure_rows(heat,defects,shape,c.mpf(qR+y0))
        B=mp.sqrt(C)*mp.exp(-bh*(qR+y0));qt=mp.sqrt(mp.exp(qR+y0)/(2*S))*B
        counts=dict(stress_mixed3=0,pressure_mixed4=0);maxmiss=mp.mpf(0);tolerance=mp.mpf('1e-55')
        def compare(value,box,label):
            nonlocal maxmiss
            lo,hi=endpoints(box);miss=max(lo-value,value-hi,mp.mpf(0));maxmiss=max(maxmiss,miss)
            if miss>tolerance:raise ArithmeticError('Independent power stress/pressure mismatch: '+label)
        for label,factor in (('theta',qt),('axial',qt*B)):
            for j in range(4):
                for n in range(4-j):
                    value=mp.diff(lambda yy,zz:reference(yy,zz)[label],(y0,z0),(j,n))/factor
                    compare(value,stress[label][j][n]*math.factorial(n),label+str((j,n)))
                    counts['stress_mixed3']+=1
        for j in range(5):
            for n in range(5-j):
                value=mp.diff(lambda yy,zz:reference(yy,zz)['pressure'],(y0,z0),(j,n))
                compare(value,pressure[j][n]*math.factorial(n),'pressure'+str((j,n)))
                counts['pressure_mixed4']+=1
        axial=[mp.diff(lambda zz:moments(y0,zz)[j],z0) for j in range(3)]
        if not all(v!=0 for v in axial):raise ArithmeticError('Nonzero moment axial fixture dependence lost')
        return dict(all_passed=True,fixture='Piecewise exponential K: axial-dependent past and complete future, Z-independent current K',
            independent_whole_history_integrals=integral_count,maximum_integral_error=str(integral_miss),
            comparison_counts=counts,maximum_positive_enclosure_miss=str(maxmiss),comparison_tolerance=str(tolerance),
            nonzero_A_Z_E_Z_P_Z_verified=True,actual_source_history_admission_is_separate=True,
            fixture_is_global_NS_validation=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'outer_power_stress_C3.json';record=json.loads((HERE/name).read_bytes());hashes=dict(record['input_hashes'])
        for filename,digest in hashes.items():
            if hashlib.sha256((HERE/filename).read_bytes()).hexdigest()!=digest:raise ValueError('Outer power source changed: '+filename)
        fresh=CompliantOuterPowerStressC3().report()
        if encode(pack(fresh))!=record:raise ValueError('Current original outer-power source/moments/join differ')
        c=MPIntervalContext();c.dps=300;signed=0;pressure=0;zeros=0
        for point in record['samples']+[record['whole_outer_power']]:
            for grid in point['outer_power_similarity_stress_mixed3_factored'].values():
                for value in grid.values():
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite outer-power stress')
                    signed+=1
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                for value in jet['coefficients'][:5-j]:
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite outer-power pressure')
                    pressure+=1
            for jet in point['outer_power_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Outer-power meridional zero lost')
                    zeros+=1
            if not point['ordinary_q_derivatives_not_unscaled_phase_derivatives']:raise ValueError('Ordinary q derivative metadata lost')
        for flag in ('actual_original_outer_power_similarity_stress_recovered','actual_outer_power_absolute_pressure_same_source_mixed4_available',
            'outer_power_angular_stress_mixed3_and_pressure_mixed4_join_verified','original_selected_angular_full_future_retained',
            'actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained'):
            if not record[flag]:raise ValueError('Outer power admission missing: '+flag)
        for flag in ('outer_power_physical_decomposition_constructed','outer_power_cone_certified','left_flatten_stress_companion_constructed',
            'global_admissible_stress_lift_constructed','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Outer power stress scope overclaimed: '+flag)
        print('Current outer-power source, signed rows and arbitrary-data right join admitted; bounded fixture running',flush=True)
        fixture=independent_piecewise_exponential_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=fresh['actual_five_defect_family_sha256'],
            implicit_source_sha256=fresh['implicit_source_sha256'],finite_signed_stress_mixed_rows=signed,
            finite_pressure_mixed_rows=pressure,exact_meridional_zero_coefficients=zeros,
            actual_original_outer_power_similarity_stress_recovered=True,actual_outer_power_absolute_pressure_same_source_mixed4_available=True,
            outer_power_angular_stress_mixed3_and_pressure_mixed4_join_verified=True,original_selected_angular_full_future_retained=True,
            actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained=True,
            actual_full_native_history_and_source_AST_bridge_verified=True,
            actual_outer_power_angular_formula_join_identities=len(fresh['source_bridge']['actual_outer_power_angular_formula_join']['identities']),
            independent_piecewise_exponential_fixture=fixture,outer_power_physical_decomposition_constructed=False,
            outer_power_cone_certified=False,left_flatten_stress_companion_constructed=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original outer-power full moments, stress, absolute pressure and angular join; physical/cone/flatten pending',flush=True)
        return result


if __name__=='__main__':run()
