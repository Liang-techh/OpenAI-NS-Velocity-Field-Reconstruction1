"""Focused original waiting history/FTC, signed derivative and pressure check."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_waiting_stress_C3 import (
    CompliantWaitingStressC3,waiting_defect_rows,waiting_stress_rows,
    waiting_pressure_rows,waiting_identities,waiting_source_bridge)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def moderate_future_fixture():
    """Independent full convergent future-integral transport fixture.

    A moderate analytic future K tests the backward integral algebra and
    derivative factors. It is not the actual Gamma source or a proof of its
    terminal history: those are consumed separately with current hashes.
    """
    with mp.workdps(65):
        a,eps,S,Z=map(mp.mpf,('.15','.002','.004','.4'))
        delta=2*a; k=1-a; p=1+delta; b=(1-delta)/2; K0=1-eps
        def g(z):return mp.mpf('.003')*(1+z*z)
        def future(v,z):return 1-eps*mp.exp(-2*v)+g(z)*v*mp.exp(-2*v)
        def terminal(label,z):
            if label=='A':
                return 1/k+mp.quad(lambda v:mp.exp((k-2)*v)*(eps-g(z)*v),[0,1,mp.inf])
            rate=delta if label=='E' else p
            return mp.quad(lambda v:mp.exp(-rate*v)*future(v,z)**2,[0,1,mp.inf])*(1 if label=='E' else mp.mpf('.5'))
        def exact_terminal(label,z):
            v=g(z)
            if label=='A':return 1/k+eps/(2-k)-v/(2-k)**2
            rate=delta if label=='E' else p
            value=1/rate-2*eps/(rate+2)+2*v/(rate+2)**2+eps**2/(rate+4)-2*eps*v/(rate+4)**2+2*v*v/(rate+4)**3
            return value*(1 if label=='E' else mp.mpf('.5'))
        quadrature=[]
        for label in ('A','E','P'):
            error=terminal(label,Z)-exact_terminal(label,Z)
            if abs(error)>mp.mpf('1e-52'):raise ArithmeticError('Complete future quadrature failed')
            quadrature.append(mp.nstr(error,10))
        c=MPIntervalContext(); c.dps=110
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),eps=c.mpf(eps),S=c.mpf(S),delta=c.mpf(delta),
                             k=c.mpf(k),prate=c.mpf(p),pressure_scale=c.mpf('1.7'))
        base={}
        for label,key,baseline in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
            values=[mp.diff(lambda z:exact_terminal(label,z),Z,n)/math.factorial(n) for n in range(6)]
            values[0]-=baseline
            # Controlled fixture input enclosure covers numerical terminal
            # evaluation; actual source endpoints are not selected this way.
            jet=IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-55'),v+mp.mpf('1e-55')]) for v in values])
            base[key+'_defect_rows']=[jet]
        def moment(q,z,label):
            if label=='A':return mp.exp(-k*q)*(exact_terminal('A',z)-K0*(1-mp.exp(k*q))/k)
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            return mp.exp(rate*q)*(exact_terminal(label,z)+factor*K0*K0*(mp.exp(-rate*q)-1)/rate)
        def fullstress(q,z,label):
            R=mp.exp(q)/S; B=mp.mpf('1.3')*mp.exp(-(mp.mpf('.5')+a)*q)
            A=moment(q,z,'A'); E=moment(q,z,'E'); P=moment(q,z,'P')
            Az=mp.diff(lambda zz:moment(q,zz,'A'),z)
            Ez=mp.diff(lambda zz:moment(q,zz,'E'),z); Pz=mp.diff(lambda zz:moment(q,zz,'P'),z)
            L=1-delta*z*z; d=1-z*z
            if label=='theta':return mp.sqrt(R/2)*B*((k*A-b*z*Az-K0)/L-2/R*(1+a)*K0)
            return mp.sqrt(R/2)*B*B*(delta*z*E-d*Ez/2-2*p*z*P+d*Pz)/L
        errors=[]; stress_checks=pressure_checks=moment_checks=0
        def encloses(row,value,label):
            lo,hi=endpoints(row); tol=mp.mpf('1e-48')
            error=max(lo-value,value-hi,mp.mpf(0))
            if error>tol:raise ArithmeticError('Waiting independent derivative mismatch: '+label+' '+str(error))
            errors.append(error)
        for q in map(mp.mpf,('-2','-.7','0')):
            rows=waiting_stress_rows(heat,base,c.mpf(Z),c.mpf(q))
            defects=waiting_defect_rows(heat,base,c.mpf(q)); pressure=waiting_pressure_rows(heat,base,c.mpf(q))
            B=mp.mpf('1.3')*mp.exp(-(mp.mpf('.5')+a)*q); R=mp.exp(q)/S
            for label,factor in (('theta',mp.sqrt(R/2)*B),('axial',mp.sqrt(R/2)*B*B)):
                for j in range(4):
                    for n in range(4-j):
                        value=mp.diff(lambda zz:mp.diff(lambda qq:fullstress(qq,zz,label),q,j),Z,n)/factor
                        encloses(rows[label][j][n]*math.factorial(n),value,label+str((j,n))); stress_checks+=1
            for label,key,baseline in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
                for j in range(3):
                    value=mp.diff(lambda qq:moment(qq,Z,label),q,j)-(baseline if j==0 else 0)
                    encloses(defects[key+'_defect_rows'][j][0],value,key+str(j)); moment_checks+=1
            for j in range(5):
                for n in range(5-j):
                    value=mp.diff(lambda zz:mp.diff(lambda qq:-mp.mpf('1.7')*mp.exp(-p*qq)*moment(qq,zz,'P'),q,j),Z,n)
                    encloses(pressure[j][n]*math.factorial(n),value,'pressure'+str((j,n))); pressure_checks+=1
        return dict(full_future_quadrature_errors=quadrature,original_stress_mixed3_derivative_checks=stress_checks,
                    absolute_pressure_mixed4_derivative_checks=pressure_checks,backward_moment_FTC_checks=moment_checks,
                    max_unenclosed_numeric_error=mp.nstr(max(errors),12),numeric_comparison_tolerance='1e-48',
                    analytic_moderate_future_fixture_only=True,actual_Gamma_terminal_history_proved_by_fixture=False,
                    actual_source_history_consumed_separately=True,passed=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'waiting_stress_C3.json'; record=json.loads((HERE/name).read_bytes())
        hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Waiting dependency changed: '+source)
        provider=CompliantWaitingStressC3()
        if waiting_identities()!=record['waiting_identities']:raise ValueError('Waiting exact proof changed')
        if waiting_source_bridge(provider.collar_source)!=record['waiting_source_bridge']:raise ValueError('Waiting source bridge changed')
        c=provider.ctx; finite=pressure_rows=zeros=joins=0
        for point in record['samples']+[record['whole_waiting']]:
            for label,grid in point['waiting_similarity_stress_mixed3_factored'].items():
                if len(grid)!=10:raise ValueError('Waiting stress mixed3 incomplete')
                for value in grid.values():
                    if any(not mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Waiting stress nonfinite')
                    finite+=1
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                if len(jet['coefficients'])!=6:raise ValueError('Waiting pressure axial jets incomplete')
                for n in range(5-j):
                    if any(not mp.isfinite(v) for v in endpoints(read_interval(c,jet['coefficients'][n]))):raise ArithmeticError('Waiting pressure nonfinite')
                    pressure_rows+=1
            for jet in point['waiting_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(0,0):raise ArithmeticError('Waiting meridional source history not zero')
                    zeros+=1
            if endpoints(read_interval(c,point['waiting_shear_strength_kappa_minus2']))[0]<=0:raise ArithmeticError('Waiting shear margin lost')
        # Source flat K jets, plus arbitrary-endpoint FTC identities above,
        # give the functional stress join. Overlap is only a diagnostic.
        for Z in ('0','.5',[-1,1]):
            shape=provider.heat.shape(c.mpf(Z),0)
            if endpoints(shape['K_rows'][0][0]-(1-provider.heat.eps))[0]>0 or endpoints(shape['K_rows'][0][0]-(1-provider.heat.eps))[1]<0:
                raise ArithmeticError('Waiting/collar K endpoint differs')
            for jet in shape['K_rows'][1:]:
                for value in jet.coefficients:
                    if endpoints(value)!=(0,0):raise ArithmeticError('Actual collar K endpoint not flat')
                    joins+=1
        for flag in ('waiting_cone_certified','global_admissible_stress_lift_constructed','physical_energy_integral_certified',
                     'temporal_recursion','steep_exit_waiting_stress_mixed3_join_verified'):
            if record[flag]:raise ValueError('Waiting scope overclaimed: '+flag)
        fixture=moderate_future_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
                    all_passed=True,actual_original_waiting_similarity_stress_recovered=True,
                    actual_waiting_absolute_pressure_same_source_mixed4_available=True,
                    waiting_collar_stress_mixed3_join_verified=True,
                    actual_finite_signed_stress_rows=finite,actual_finite_pressure_rows=pressure_rows,
                    exact_meridional_history_zero_rows=zeros,actual_flat_endpoint_K_derivative_checks=joins,
                    independent_full_future_transport_fixture=fixture,
                    steep_exit_waiting_stress_mixed3_join_verified=False,waiting_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,
                    temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS actual waiting full-future FTC/stress/pressure source and independent transport fixture',flush=True)
        return result


if __name__=='__main__':run()
