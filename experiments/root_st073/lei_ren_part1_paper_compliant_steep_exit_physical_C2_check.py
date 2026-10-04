"""Original steep-exit physical map, completed tensor and regional error."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_steep_exit_physical_C2 import (
    steep_exit_physical_identities,steep_exit_adapter_binding,steep_exit_divergence_grids,
    exit_pressure_fluctuation)
from lei_ren_part1_paper_compliant_collar_physical_C2 import collar_velocity_bracket
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def original_sigma_positive_pressure_fixture():
    """Direct original-sigma quadrature for the new positive source kernel."""
    from scipy.integrate import quad
    c=MPIntervalContext(); c.dps=110
    def sigma(v):
        if v<=0:return 0.
        if v>=1:return 1.
        odds=1/(1-v)**2-1/v**2
        if odds>=0:return 1/(1+math.exp(-odds))
        value=math.exp(odds); return value/(1+value)
    def primitive(v):return quad(lambda q:1-sigma(q),v,1,epsabs=2e-13,epsrel=2e-13)[0]
    k=.85; eps=.002; p=1.3; K0=1-eps; checks=0
    for t in (.3,.7):
        record=exit_pressure_fluctuation(c,c.mpf(str(t)),c.mpf(str(k)),c.mpf(str(eps)),c.mpf(str(p)))
        value=quad(lambda v:math.exp(-p*(v-t))*K0*K0*math.expm1(2*k*primitive(v)),t,1,
                   epsabs=2e-12,epsrel=2e-12)[0]
        for label,expected in (('integral',value),('current_delta_square',K0*K0*math.expm1(2*k*primitive(t)))):
            lo,hi=endpoints(record[label])
            if max(lo-expected,expected-hi,mp.mpf(0))>mp.mpf('1e-10'):raise ArithmeticError('Original-sigma positive pressure kernel mismatch')
            checks+=1
    terminal=exit_pressure_fluctuation(c,c.mpf(1),c.mpf(str(k)),c.mpf(str(eps)),c.mpf(str(p)))
    if any(endpoints(terminal[key])!=(0,0) for key in ('integral','current_delta_square','f')):
        raise ArithmeticError('Original positive pressure kernel endpoint not zero')
    return dict(original_sigma_direct_pressure_fluctuation_checks=checks,exact_endpoint_zero_checks=3,
                numeric_diagnostic_tolerance='1e-10',formal_integral_split_proved_separately=True,passed=True)


def independent_native_cartesian_fixture():
    """Complete future and variable K, differentiated in native Cartesian space.

    A moderate analytic polynomial exit K is used to differentiate every
    tensor entry, including r*Tz_z, without repeatedly differentiating sigma
    quadrature. The actual original sigma and Gamma histories are admitted
    separately by the source bridge and current stress receipts. This fixture
    tests the physical general-K transfer and the new divergence cancellation.
    """
    with mp.workdps(60):
        a,eps,S,B0,wait,corr=map(mp.mpf,('.15','.002','.004','1.3','1.3','.06'))
        delta=2*a; k=1-a; p=1+delta; bh=mp.mpf('.5')+a; b=(1-delta)/2; K0=1-eps; beta=-2-delta
        def future(label,Z):
            g=mp.mpf('.003')*(1+Z*Z)
            if label=='A':return 1/k+eps/(2-k)-g/(2-k)**2
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            return factor*(1/rate-2*eps/(rate+2)+2*g/(rate+2)**2+eps**2/(rate+4)
                           -2*eps*g/(rate+4)**2+2*g*g/(rate+4)**3)
        def left(label,Z):
            if label=='A':return mp.exp(k*wait)*future(label,Z)-K0*mp.expm1(k*wait)/k
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            return mp.exp(-rate*wait)*future(label,Z)+factor*K0*K0*(-mp.expm1(-rate*wait))/rate
        def weighted_power(rate,h,n):
            # Exact integral_0^h exp(rate*v)*v^n dv, from its finite primitive.
            series=sum((-1)**j*mp.factorial(n)/mp.factorial(n-j)*h**(n-j)/rate**(j+1) for j in range(n+1))
            return mp.exp(rate*h)*series-(-1)**n*mp.factorial(n)/rate**(n+1)
        def K(q):return K0+corr*(-wait-q)**5
        def moment(q,Z,label):
            h=-wait-q
            if label=='A':
                integral=K0*weighted_power(-k,h,0)+corr*weighted_power(-k,h,5)
                return mp.exp(k*h)*(left(label,Z)-integral)
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            integral=K0*K0*weighted_power(rate,h,0)+2*K0*corr*weighted_power(rate,h,5)+corr*corr*weighted_power(rate,h,10)
            return mp.exp(-rate*h)*(left(label,Z)+factor*integral)
        def coordinates(r,z,tau,nu):
            zs=z/mp.sqrt(nu)
            lam=mp.findroot(lambda value:value*value-zs*zs*value**(2*delta)-tau,mp.sqrt(tau)+abs(zs))
            Z=zs/lam**(1-delta); R=r*r/(2*nu*lam*lam)
            return lam,mp.log(R*S),Z,R
        def swirl(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu)
            return mp.sqrt(nu)*lam**(-1-delta)*B0*mp.exp(-bh*q)*K(q)
        def pressure(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu)
            return -nu*lam**(-2-2*delta)*B0*B0*mp.exp(-p*q)*moment(q,Z,'P')
        def stress(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu); B=B0*mp.exp(-bh*q); L=1-delta*Z*Z; d=1-Z*Z
            A=moment(q,Z,'A'); E=moment(q,Z,'E'); P=moment(q,Z,'P')
            Az=mp.diff(lambda value:moment(q,value,'A'),Z)
            Ez=mp.diff(lambda value:moment(q,value,'E'),Z); Pz=mp.diff(lambda value:moment(q,value,'P'),Z)
            Ct=(k*A-b*Z*Az-K(q))/L+2/R*(mp.diff(K,q)-(1+a)*K(q))
            Cz=(delta*Z*E-d*Ez/2-2*p*Z*P+d*Pz)/L
            return nu*lam**beta*mp.sqrt(R/2)*B*Ct,nu*lam**beta*mp.sqrt(R/2)*B*B*Cz
        errors=[]; divergences=[]; remainder_errors=[]; div_errors=[]; nonzero=0
        tolerance=mp.mpf('1e-36')
        for phase,nu in ((mp.mpf('.3'),mp.mpf('.01')),(mp.mpf('.75'),mp.mpf('.7'))):
            print('Independent native Cartesian variable-K exit phase '+str(phase)+' nu '+str(nu),flush=True)
            q=-wait-1+phase; Z=mp.mpf('.4'); tau=mp.mpf('.7'); lam=mp.sqrt(tau/(1-Z*Z)); angle=mp.mpf('.6')
            R=mp.exp(q)/S; r=mp.sqrt(nu)*lam*mp.sqrt(2*R); z=mp.sqrt(nu)*lam**(1-delta)*Z
            point=[r*mp.cos(angle),r*mp.sin(angle),z]
            def vel(x,y,zz,ta=tau):
                rr=mp.sqrt(x*x+y*y); value=swirl(rr,zz,ta,nu)
                return [-y/rr*value,x/rr*value,mp.mpf(0)]
            def pres(x,y,zz):return pressure(mp.sqrt(x*x+y*y),zz,tau,nu)
            def tensor(x,y,zz):
                rr=mp.sqrt(x*x+y*y); cs=x/rr; sn=y/rr; tt,tz=stress(rr,zz,tau,nu)
                diagonal=rr*mp.diff(lambda value:stress(rr,value,tau,nu)[1],zz)
                return [[-2*cs*sn*tt+sn*sn*diagonal,(cs*cs-sn*sn)*tt-cs*sn*diagonal,cs*tz],
                        [(cs*cs-sn*sn)*tt-cs*sn*diagonal,2*cs*sn*tt+cs*cs*diagonal,sn*tz],
                        [cs*tz,sn*tz,mp.mpf(0)]]
            def partial(fn,j,n=1):
                def evaluate(value):
                    args=list(point); args[j]=value; return fn(*args)
                return mp.diff(evaluate,point[j],n)
            u=vel(*point); Etheta=-nu*mp.diff(lambda zz:swirl(r,zz,tau,nu),z,2)
            if abs(Etheta)<mp.mpf('1e-16'):raise ArithmeticError('Fixture does not exercise nonzero axial viscosity')
            nonzero+=1; Ecart=[-mp.sin(angle)*Etheta,mp.cos(angle)*Etheta,mp.mpf(0)]
            for component in range(3):
                dt=-mp.diff(lambda value:vel(*point,value)[component],tau)
                conv=sum(u[j]*partial(lambda x,y,zz:vel(x,y,zz)[component],j) for j in range(3))
                lap=sum(partial(lambda x,y,zz:vel(x,y,zz)[component],j,2) for j in range(3))
                div=sum(partial(lambda x,y,zz:tensor(x,y,zz)[component][j],j) for j in range(3))
                error=dt+conv+partial(pres,component)-nu*lap+div-Ecart[component]
                if abs(error)>tolerance:raise ArithmeticError('Exit completed Cartesian decomposition failed: '+mp.nstr(error,12))
                errors.append(mp.nstr(error,12))
            divu=sum(partial(lambda x,y,zz:vel(x,y,zz)[j],j) for j in range(3))
            if abs(divu)>tolerance:raise ArithmeticError('Native physical divergence failed')
            divergences.append(mp.nstr(divu,12))
            # Independent full native moments above vs the actual new reduced
            # divergence code. Analytic fixture supplies its exact Iplus.
            c=MPIntervalContext(); c.dps=110
            heat=SimpleNamespace(ctx=c,a=c.mpf(a),eps=c.mpf(eps),S=c.mpf(S),delta=c.mpf(delta),k=c.mpf(k),prate=c.mpf(p))
            base={'Z':c.mpf(Z)}
            values=[mp.diff(lambda zz:future('P',zz),Z,n)/math.factorial(n) for n in range(6)]; values[0]-=1/(2*p)
            base['pressure_defect_rows']=[IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-52'),v+mp.mpf('1e-52')]) for v in values])]
            Krows=[IntervalTaylor.constant(c,mp.diff(K,q,j),5) for j in range(5)]
            Qrows=[IntervalTaylor.constant(c,mp.diff(lambda qq:K(qq)**2-1,q,j),5) for j in range(5)]
            shape=dict(K_rows=Krows,K_squared_defect_rows=Qrows); h=-wait-q
            positive=2*K0*corr*h**5+corr*corr*h**10
            Iplus=mp.exp(-p*h)*(2*K0*corr*weighted_power(p,h,5)+corr*corr*weighted_power(p,h,10))
            fluctuation=dict(current_delta_square=c.mpf(positive),integral=c.mpf(Iplus))
            grids=steep_exit_divergence_grids(heat,base,shape,c.mpf(q),fluctuation)
            B=B0*mp.exp(-bh*q); Qt=mp.sqrt(R/2)*B; Qz=mp.sqrt(R/2)*B*B
            def native_div(rr,zz,label):
                component=0 if label=='theta' else 1
                return mp.diff(lambda value:stress(value,zz,tau,nu)[component],rr)+(2 if component==0 else 1)*stress(rr,zz,tau,nu)[component]/rr
            def native_error(rr,zz):return -nu*mp.diff(lambda value:swirl(rr,value,tau,nu),zz,2)
            for i in range(3):
                for j in range(3-i):
                    for label,factor in (('theta',Qt),('axial',Qz)):
                        value=mp.diff(lambda zz:mp.diff(lambda rr:native_div(rr,zz,label),r,i),z,j)
                        coefficient=physical_bracket(c,grids[label],i,j,c.mpf(Z),c.mpf(delta),c.mpf(beta-1))
                        scale=nu**(mp.mpf('.5')-mp.mpf(i+j)/2)*lam**(beta-1-i+j*(delta-1))*factor*(2/R)**(mp.mpf(i+1)/2)
                        lo,hi=endpoints(coefficient*c.mpf(scale)); error=max(lo-value,value-hi,mp.mpf(0))
                        if error>tolerance:raise ArithmeticError('Exit reduced divergence mixed2 mismatch: '+str((label,i,j,error)))
                        div_errors.append(mp.nstr(error,12))
                    value=mp.diff(lambda zz:mp.diff(lambda rr:native_error(rr,zz),r,i),z,j)
                    coefficient=-collar_velocity_bracket(c,Krows,i,j+2,c.mpf(Z),c.mpf(delta))
                    scale=nu**(mp.mpf('.5')-mp.mpf(i+j)/2)*lam**(-1-delta-i+(j+2)*(delta-1))*B*(2/R)**(mp.mpf(i)/2)
                    lo,hi=endpoints(coefficient*c.mpf(scale)); error=max(lo-value,value-hi,mp.mpf(0))
                    if error>tolerance:raise ArithmeticError('Exit axial viscosity mixed2 mismatch: '+str((i,j,error)))
                    remainder_errors.append(mp.nstr(error,12))
        return dict(Cartesian_decomposition_errors=errors,physical_divergence_errors=divergences,
                    full_completed_tensor_native_derivatives_exercised=True,nonzero_axial_viscosity_remainders_exercised=nonzero,
                    reduced_divergence_mixed2_checks=len(div_errors),reduced_divergence_mixed2_errors=div_errors,
                    nonzero_remainder_mixed2_checks=len(remainder_errors),remainder_mixed2_errors=remainder_errors,
                    numeric_tolerance=str(tolerance),viscosities=['.01','.7'],native_implicit_coordinates_differentiated=True,
                    polynomial_nonconstant_K_and_analytic_full_future_fixture_only=True,
                    actual_sigma_and_Gamma_history_consumed_separately=True,passed=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'steep_exit_physical_C2.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Exit physical source changed: '+path)
        if record['steep_exit_physical_identities']!=steep_exit_physical_identities():raise ValueError('Exit physical proof changed')
        if record['steep_exit_adapter_binding']!=steep_exit_adapter_binding():raise ValueError('Exit adapter changed')
        c=MPIntervalContext(); c.dps=300; finite=zeros=terminal_zeros=0
        def rowcheck(row,require_zero=False):
            nonlocal finite,zeros,terminal_zeros
            lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
            if not all(mp.isfinite(v) for v in (lo,hi)):raise ArithmeticError('Exit signed physical row nonfinite')
            if row['exact_zero']:
                if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Invalid exit zero row')
                zeros+=1; terminal_zeros+=require_zero
            else:
                if require_zero:raise ArithmeticError('Exit endpoint remainder not exactly zero')
                if any(not mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):raise ArithmeticError('Exit physical log bound nonfinite')
                finite+=1
        for point in record['samples']+[record['whole_steep_exit']]:
            for flag in ('actual_regional_physical_steep_exit_stress_remainder_identity_verified',
                         'steep_exit_homogeneous_A_E_and_stationary_pressure_divergence_cancelled_before_enclosure',
                         'steep_exit_waiting_physical_stress_mixed3_and_remainder_mixed2_join_verified'):
                if not point[flag]:raise ValueError('Exit physical gate missing: '+flag)
            for group,count in (('physical_cylindrical_stress_mixed3',10),('physical_cylindrical_stress_divergence_mixed2',6)):
                for grid in point[group].values():
                    if len(grid)!=count:raise ValueError('Exit physical mixed coverage incomplete')
                    for row in grid.values():rowcheck(row)
            terminal=point['steep_exit_terminal_remainder_and_axial_time_velocity_jets_exact_zero']
            for row in point['completed_theta_theta_stress_mixed2'].values():rowcheck(row)
            for row in point['physical_angular_axial_viscosity_remainder_mixed2'].values():rowcheck(row,terminal)
            if terminal and not point['physical_axial_and_time_velocity_rows_all_exact_zero']:raise ValueError('Exit flat endpoint velocity zero jets missing')
            if terminal:
                if endpoints(read_interval(c,point['steep_exit_pressure_fluctuation']['integral']))!=(0,0):raise ArithmeticError('Stationary pressure not removed at waiting join')
            for label in ('exact_physical_radial_momentum_residual','exact_completed_tensor_radial_divergence',
                          'exact_physical_radial_remainder','exact_physical_axial_remainder','exact_physical_divergence'):
                if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Exit structural zero lost')
                zeros+=1
            for flag in ('steep_exit_cone_certified','steep_exit_regional_remainder_exact_zero','whole_outer_cone_certified',
                         'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
                         'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion'):
                if point[flag]:raise ValueError('Exit physical scope overclaimed: '+flag)
        positive_fixture=original_sigma_positive_pressure_fixture()
        fixture=independent_native_cartesian_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],implicit_source_sha256=record['implicit_source_sha256'],
                    all_passed=True,actual_regional_physical_steep_exit_stress_remainder_identity_verified=True,
                    steep_exit_waiting_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                    actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
                    actual_terminal_remainder_mixed2_exact_zero_rows=terminal_zeros,independent_native_Cartesian_fixture=fixture,
                    independent_original_sigma_positive_pressure_fixture=positive_fixture,
                    steep_exit_regional_remainder_exact_zero=False,steep_exit_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
                    physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original steep-exit physical completed stress, nonzero remainder and native Cartesian decomposition',flush=True)
        return result


if __name__=='__main__':run()
