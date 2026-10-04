"""Steep power: full-source FTC, functional exit join and native fixture."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_steep_power_stress_C3 import (
    CompliantSteepPowerStressC3,power_shape,power_defect_rows,power_stress_rows,
    power_pressure_rows,power_transport_identities,power_source_bridge,DOMAIN,quotient_rows,power_angular_rows)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def moderate_full_future_fixture():
    """Original finite power integrals plus an independent infinite future.

    Complete future moments come from KQ*exp(-rho*v)*(1+g(Z)*v^2).
    Native original unnormalized stresses and pressure are differentiated
    independently; production defect/stress recurrences are not their oracle.
    Actual sigma/Gamma history and joins are consumed separately by the source
    bridge. This moderate fixture is not a global NS or cone certificate.
    """
    with mp.workdps(80):
        a,eps,S,wait,Z,C=map(mp.mpf,('.15','.002','.004','1.3','.4','1.7'))
        delta=2*a; k=1-a; p=1+delta; bh=mp.mpf('.5')+a; b=(1-delta)/2
        Ts=4*mp.log(2/delta); qQ=-wait-1; KQ=(1-eps)*mp.exp(k/2); rho=k+mp.mpf('.07')
        c=MPIntervalContext(); c.dps=120
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),eps=c.mpf(eps),S=c.mpf(S),delta=c.mpf(delta),
                             k=c.mpf(k),prate=c.mpf(p),pressure_scale=c.mpf(C),
                             steep=SimpleNamespace(wait=c.mpf(wait),Ts=c.mpf(Ts)))
        def g(z):return mp.mpf('.003')*(1+z*z)
        def future(label,z):
            gg=g(z)
            if label=='A':
                rate=rho-k
                return KQ*(1/rate+2*gg/rate**3)
            rate=(delta if label=='E' else p)+2*rho
            return KQ*KQ*(1/rate+4*gg/rate**3+24*gg*gg/rate**5)*(1 if label=='E' else mp.mpf('.5'))
        terminal=dict(KQ=c.mpf(KQ))
        xscoeff=[mp.diff(lambda z:future('A',z)/KQ-Ts,Z,n)/math.factorial(n) for n in range(6)]
        terminal['original_XS']=IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-70'),v+mp.mpf('1e-70')]) for v in xscoeff])
        for label,key,unit in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
            coefficients=[mp.diff(lambda z:future(label,z),Z,n)/math.factorial(n) for n in range(6)]
            coefficients[0]-=unit
            terminal[key+'_defect_rows']=[IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-70'),v+mp.mpf('1e-70')]) for v in coefficients])]
        def original_K(qq):return KQ*mp.exp(-k*(qq-qQ))
        def moment(label,qq,z):
            x=qq-qQ; K=original_K(qq)
            if label=='A':return mp.exp(-k*x)*(future('A',z)+KQ*x)
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            decay_rate=rate+2*k
            # Original positive finite integral, evaluated in closed form.
            finite=K*K*(-mp.expm1(decay_rate*x))/decay_rate*factor
            return finite+mp.exp(rate*x)*future(label,z)
        def native_stress(qq,z,label):
            R=mp.exp(qq)/S; B=mp.sqrt(C)*mp.exp(-bh*qq); K=original_K(qq)
            L=1-delta*z*z; d=1-z*z
            A=moment('A',qq,z); E=moment('E',qq,z); Pr=moment('P',qq,z)
            Mt=mp.sqrt(2)*R**mp.mpf('1.5')*B*A
            Mtz=mp.sqrt(2)*R**mp.mpf('1.5')*B*mp.diff(lambda zz:moment('A',qq,zz),z)
            Me=R*B*B*E/2; Mez=R*B*B/2*mp.diff(lambda zz:moment('E',qq,zz),z)
            pressure=-B*B*Pr; pressurez=-B*B*mp.diff(lambda zz:moment('P',qq,zz),z)
            if label=='theta':
                inertial=(k*Mt-b*z*Mtz-R*mp.sqrt(2*R)*B*K)/(2*L*R)
                derivative=mp.diff(lambda vv:mp.sqrt(C)*mp.exp(-bh*vv)*original_K(vv),qq)
                return inertial+mp.sqrt(2*R)/R*derivative-B*K/mp.sqrt(2*R)
            return (2*delta*z*Me-d*Mez+R*(2*p*z*pressure-d*pressurez))/(L*mp.sqrt(2*R))
        errors=[]; counts=dict(finite_original_integrals=0,full_moment_derivatives=0,
                              normalized_quotient_derivatives=0,native_original_stress_mixed3=0,
                              native_absolute_pressure_mixed4=0)
        tolerance=mp.mpf('1e-60')
        def encloses(row,value,label):
            lo,hi=endpoints(row); error=max(lo-value,value-hi,mp.mpf(0)); errors.append(error)
            if error>tolerance:raise ArithmeticError('Steep power independent fixture mismatch: '+label+' '+mp.nstr(error,12))
        for phase in map(mp.mpf,('.2','.8')):
            qq=qQ-Ts*(1-phase); ell=qQ-qq
            shape=power_shape(heat,terminal,c.mpf(phase)); defects=power_defect_rows(heat,terminal,c.mpf(phase),shape)
            rows=power_stress_rows(heat,terminal,c.mpf(phase),shape,c.mpf(Z))
            pressure_rows=power_pressure_rows(heat,terminal,c.mpf(phase),shape)
            one=IntervalTaylor.constant(c,1,5)
            ordinary=dict(A=[one/heat.k+defects['angular_defect_rows'][0]]+defects['angular_defect_rows'][1:],
                          E=[one/heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:],
                          P=[one/(2*heat.prate)+defects['pressure_defect_rows'][0]]+defects['pressure_defect_rows'][1:])
            angular=quotient_rows(ordinary['A'],shape['K_rows'])
            energy=quotient_rows(ordinary['E'],[v*2 for v in product_rows(shape['K_rows'],shape['K_rows'])])
            # Independently integrate each actual finite moment density.
            for label in ('A','E','P'):
                if label=='A':
                    finite=mp.quad(lambda v:mp.exp(k*(v-qq))*original_K(v),[qq,qQ])
                    direct=mp.exp(k*ell)*future(label,Z)-finite
                else:
                    rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
                    finite=mp.quad(lambda v:mp.exp(-rate*(v-qq))*original_K(v)**2*factor,[qq,qQ])
                    direct=finite+mp.exp(-rate*ell)*future(label,Z)
                encloses(ordinary[label][0][0],direct,label+' full finite integral'); counts['finite_original_integrals']+=1
                for j in range(5):
                    encloses(ordinary[label][j][0],mp.diff(lambda qv:moment(label,qv,Z),qq,j),label+' y'+str(j))
                    counts['full_moment_derivatives']+=1
            for label,values in (('angular',angular),('energy',energy)):
                if label=='angular':reference=lambda qv:moment('A',qv,Z)/original_K(qv)
                else:reference=lambda qv:moment('E',qv,Z)/(2*original_K(qv)**2)
                for j in range(5):
                    encloses(values[j][0],mp.diff(reference,qq,j),label+' quotient y'+str(j))
                    counts['normalized_quotient_derivatives']+=1
            R=mp.exp(qq)/S; B=mp.sqrt(C)*mp.exp(-bh*qq)
            for label in ('theta','axial'):
                factor=mp.sqrt(R/2)*B*(1 if label=='theta' else B)
                for j in range(4):
                    for n in range(4-j):
                        native=mp.diff(lambda qv,zv:native_stress(qv,zv,label),(qq,Z),(j,n))/factor
                        encloses(rows[label][j][n]*math.factorial(n),native,label+' y'+str(j)+' Z'+str(n))
                        counts['native_original_stress_mixed3']+=1
            for j in range(5):
                for n in range(5-j):
                    native=mp.diff(lambda qv,zv:-C*mp.exp(-p*qv)*moment('P',qv,zv),(qq,Z),(j,n))
                    encloses(pressure_rows[j][n]*math.factorial(n),native,'pressure y'+str(j)+' Z'+str(n))
                    counts['native_absolute_pressure_mixed4']+=1
        return dict(passed=True,checks=counts,tolerance=mp.nstr(tolerance,12),
                    maximum_interval_violation=mp.nstr(max(errors),16),
                    complete_convergent_future_used=True,actual_sigma_Gamma_history_proved_by_source_bridge_separately=True,
                    original_unnormalized_stress_differentiated_natively=True,global_NS_or_actual_cone_certified=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'steep_power_stress_C3.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power source changed: '+source)
        provider=CompliantSteepPowerStressC3(); c=provider.ctx
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(provider.family,provider.source):
            raise ValueError('Power common family/source differs')
        if provider.proof!=record['transport_identities']:raise ValueError('Power exact FTC/mode proof changed')
        if provider.bridge!=record['steep_power_source_bridge']:raise ValueError('Power original history/join binding changed')
        finite=pressure=zero=0
        for point in record['samples']+[record['whole_steep_power']]:
            for label,rows in point['steep_power_similarity_stress_mixed3_factored'].items():
                if set(rows)!={'y'+str(j)+'_Z'+str(n) for j in range(4) for n in range(4-j)}:
                    raise ValueError('Power mixed3 stress grid incomplete: '+label)
                for value in rows.values():
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Power nonfinite stress')
                    finite+=1
            if len(point['pressure_y_derivative_axial5_Taylor'])!=5:raise ValueError('Power pressure rows incomplete')
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                if len(jet['coefficients'])!=6:raise ValueError('Power pressure axial5 incomplete')
                for n in range(5-j):
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,jet['coefficients'][n]))):raise ArithmeticError('Power nonfinite pressure')
                    pressure+=1
            for jet in point['steep_power_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(0,0):raise ArithmeticError('Power meridional source zero lost')
                    zero+=1
            if endpoints(read_interval(c,point['steep_power_shear_strength_kappa_minus2']))!=(2,2):
                raise ArithmeticError('Actual pure power shear strength lost')
            angular=point['angular_y_derivative_Taylor']
            for j in range(1,5):
                for n,value in enumerate(angular[j]['coefficients']):
                    expected=(1,1) if j==1 and n==0 else (0,0)
                    if endpoints(read_interval(c,value))!=expected:raise ArithmeticError('Exact source power angular slope/flat jets lost')
        if record['domain']!=DOMAIN:raise ValueError('Whole original power domain changed')
        for flag in ('actual_original_steep_power_similarity_stress_recovered',
                     'actual_steep_power_absolute_pressure_same_source_mixed4_available',
                     'steep_power_exit_stress_mixed3_join_verified','steep_power_exit_pressure_mixed4_join_verified',
                     'actual_same_full_angular_future_recovered_with_source_XS_identity'):
            if not record[flag]:raise ValueError('Power recovery flag missing: '+flag)
        for flag in ('steep_power_cone_certified','steep_power_physical_stress_remainder_identity_verified',
                     'global_admissible_stress_lift_constructed','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Power scope overclaimed: '+flag)
        fixture=moderate_full_future_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
                    actual_original_steep_power_similarity_stress_recovered=True,
                    actual_steep_power_absolute_pressure_same_source_mixed4_available=True,
                    steep_power_exit_stress_mixed3_join_verified=True,steep_power_exit_pressure_mixed4_join_verified=True,
                    actual_same_full_angular_future_recovered_with_source_XS_identity=True,
                    domain=DOMAIN,actual_finite_signed_stress_rows=finite,actual_finite_pressure_rows=pressure,
                    exact_meridional_source_zero_rows=zero,independent_complete_future_power_fixture=fixture,
                    steep_power_cone_certified=False,steep_power_physical_stress_remainder_identity_verified=False,
                    global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,
                    temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS actual complete Ts steep power stress/pressure and original functional exit join',flush=True)
        return result


if __name__=='__main__':run()
