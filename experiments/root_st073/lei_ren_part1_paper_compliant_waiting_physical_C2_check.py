"""Waiting physical map, completed tensor, exact remainder and divergence check."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_waiting_physical_C2 import (
    waiting_physical_identities,waiting_adapter_binding,waiting_divergence_grids)
from lei_ren_part1_paper_compliant_waiting_stress_C3 import waiting_stress_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def moderate_cartesian_fixture():
    """Independent native implicit map and full completed tensor derivatives.

    Uses the convergent analytic future from the waiting FTC fixture. The
    actual Gamma terminal history is separately consumed, not replaced.
    All physical tensor entries, including r*d_z(Tz), are differentiated.
    """
    with mp.workdps(60):
        a,eps,S,B0=map(mp.mpf,('.15','.002','.004','1.3'))
        delta=2*a; k=1-a; bh=mp.mpf('.5')+a; p=1+delta; b=(1-delta)/2; K0=1-eps
        def terminal(label,Z,n=0):
            g=mp.mpf('.003')*(1+Z*Z); gz=mp.mpf('.006')*Z
            if label=='A':return (1/k+eps/(2-k)-g/(2-k)**2) if n==0 else -gz/(2-k)**2
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            if n==0:return factor*(1/rate-2*eps/(rate+2)+2*g/(rate+2)**2+eps**2/(rate+4)-2*eps*g/(rate+4)**2+2*g*g/(rate+4)**3)
            return factor*(2*gz/(rate+2)**2-2*eps*gz/(rate+4)**2+4*g*gz/(rate+4)**3)
        def moment(q,Z,label,n=0):
            if label=='A':return mp.exp(-k*q)*(terminal(label,Z,n)-(K0*(1-mp.exp(k*q))/k if n==0 else 0))
            rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
            return mp.exp(rate*q)*(terminal(label,Z,n)+(factor*K0*K0*(mp.exp(-rate*q)-1)/rate if n==0 else 0))
        def coordinates(r,z,tau,nu):
            zs=z/mp.sqrt(nu)
            lam=mp.findroot(lambda l:l*l-zs*zs*l**(2*delta)-tau,mp.sqrt(tau)+abs(zs))
            Z=zs/lam**(1-delta); R=r*r/(2*nu*lam*lam)
            return lam,mp.log(R*S),Z,R
        def swirl(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu)
            return mp.sqrt(nu)*lam**(-1-delta)*B0*mp.exp(-bh*q)*K0
        def pressure(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu)
            return -nu*lam**(-2-2*delta)*B0*B0*mp.exp(-p*q)*moment(q,Z,'P')
        def stress(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu); B=B0*mp.exp(-bh*q)
            Qt=mp.sqrt(R/2)*B; Qz=mp.sqrt(R/2)*B*B; L=1-delta*Z*Z; d=1-Z*Z
            A=moment(q,Z,'A'); Az=moment(q,Z,'A',1)
            E=moment(q,Z,'E'); Ez=moment(q,Z,'E',1); P=moment(q,Z,'P'); Pz=moment(q,Z,'P',1)
            Ct=(k*A-b*Z*Az-K0)/L-2/R*(1+a)*K0
            Cz=(delta*Z*E-d*Ez/2-2*p*Z*P+d*Pz)/L
            return nu*lam**(-2-delta)*Qt*Ct,nu*lam**(-2-delta)*Qz*Cz
        errors=[]; deriv_errors=[]; remainders=[]; div_errors=[]; nonzero_stress=0
        for q,nu in ((mp.mpf('-2'),mp.mpf('.01')),(mp.mpf('-.7'),mp.mpf('.7'))):
            print('Independent native Cartesian waiting fixture q '+str(q)+' nu '+str(nu),flush=True)
            Z=mp.mpf('.4'); tau=mp.mpf('.7'); lam=mp.sqrt(tau/(1-Z*Z)); angle=mp.mpf('.6')
            R=mp.exp(q)/S; r=mp.sqrt(nu)*lam*mp.sqrt(2*R); z=mp.sqrt(nu)*lam**(1-delta)*Z
            point=[r*mp.cos(angle),r*mp.sin(angle),z]
            def vel(x,y,zz,ta=tau):
                rr=mp.sqrt(x*x+y*y); u=swirl(rr,zz,ta,nu)
                return [-y/rr*u,x/rr*u,mp.mpf(0)]
            def pres(x,y,zz):return pressure(mp.sqrt(x*x+y*y),zz,tau,nu)
            def tensor(x,y,zz):
                rr=mp.sqrt(x*x+y*y); cs=x/rr; sn=y/rr
                tt,tz=stress(rr,zz,tau,nu)
                diagonal=rr*mp.diff(lambda w:stress(rr,w,tau,nu)[1],zz)
                return [[-2*cs*sn*tt+sn*sn*diagonal,(cs*cs-sn*sn)*tt-cs*sn*diagonal,cs*tz],
                        [(cs*cs-sn*sn)*tt-cs*sn*diagonal,2*cs*sn*tt+cs*cs*diagonal,sn*tz],
                        [cs*tz,sn*tz,mp.mpf(0)]]
            def partial(fn,j,n=1):
                def evaluate(value):
                    args=list(point); args[j]=value; return fn(*args)
                return mp.diff(evaluate,point[j],n)
            u=vel(*point)
            for i in range(3):
                dt=-mp.diff(lambda ta:vel(*point,ta)[i],tau)
                conv=sum(u[j]*partial(lambda x,y,zz:vel(x,y,zz)[i],j) for j in range(3))
                lap=sum(partial(lambda x,y,zz:vel(x,y,zz)[i],j,2) for j in range(3))
                momentum=dt+conv+partial(pres,i)-nu*lap
                divT=sum(partial(lambda x,y,zz:tensor(x,y,zz)[i][j],j) for j in range(3))
                error=momentum+divT
                if abs(error)>mp.mpf('1e-43'):raise ArithmeticError('Waiting Cartesian decomposition failed: '+str(error))
                errors.append(mp.nstr(error,12))
            e=-nu*mp.diff(lambda zz:swirl(r,zz,tau,nu),z,2)
            if abs(e)>mp.mpf('1e-43'):raise ArithmeticError('Waiting axial viscosity should vanish')
            remainders.append(mp.nstr(e,12))
            divu=sum(partial(lambda x,y,zz:vel(x,y,zz)[j],j) for j in range(3))
            if abs(divu)>mp.mpf('1e-43'):raise ArithmeticError('Waiting physical divergence failed')
            div_errors.append(mp.nstr(divu,12))
            nonzero_stress+=sum(abs(v)>mp.mpf('1e-20') for v in stress(r,z,tau,nu))
            # Check the actual NEW reduced divergence implementation and its
            # mixed2 pullback against native differentiation of full stress.
            c=MPIntervalContext(); c.dps=100
            heat=SimpleNamespace(ctx=c,a=c.mpf(a),delta=c.mpf(delta),k=c.mpf(k),prate=c.mpf(p),eps=c.mpf(eps),S=c.mpf(S))
            base={}
            for label,key,baseline in (('A','angular',1/k),('E','energy',1/delta),('P','pressure',1/(2*p))):
                values=[mp.diff(lambda zz:terminal(label,zz),Z,n)/math.factorial(n) for n in range(6)]
                values[0]-=baseline
                base[key+'_defect_rows']=[IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-50'),v+mp.mpf('1e-50')]) for v in values])]
            fields=waiting_stress_rows(heat,base,c.mpf(Z),c.mpf(q)); grids=waiting_divergence_grids(heat,fields)
            B=B0*mp.exp(-bh*q); Qt=mp.sqrt(R/2)*B; Qz=mp.sqrt(R/2)*B*B; beta=-2-delta
            def native_div(rr,zz,label):
                index=0 if label=='theta' else 1
                return mp.diff(lambda radial:stress(radial,zz,tau,nu)[index],rr)+(2 if index==0 else 1)*stress(rr,zz,tau,nu)[index]/rr
            for label,factor in (('theta',Qt),('axial',Qz)):
                for i in range(3):
                    for j in range(3-i):
                        value=mp.diff(lambda zz:mp.diff(lambda rr:native_div(rr,zz,label),r,i),z,j)
                        coefficient=physical_bracket(c,grids[label],i,j,c.mpf(Z),c.mpf(delta),c.mpf(beta-1))
                        scale=nu**(mp.mpf('.5')-mp.mpf(i+j)/2)*lam**(beta-1-i+j*(delta-1))*factor*(2/R)**(mp.mpf(i+1)/2)
                        lo,hi=endpoints(coefficient*c.mpf(scale))
                        error=max(lo-value,value-hi,mp.mpf(0))
                        if error>mp.mpf('1e-42'):raise ArithmeticError('Waiting reduced divergence mixed2 failed: '+str((label,i,j,error)))
                        deriv_errors.append(mp.nstr(error,12))
        if nonzero_stress!=4:raise ArithmeticError('Both nonzero waiting stresses were not exercised')
        return dict(Cartesian_decomposition_errors=errors,full_completed_tensor_native_derivatives_exercised=True,
                    exact_zero_remainder_numeric_errors=remainders,physical_divergence_errors=div_errors,
                    reduced_divergence_mixed2_checks=len(deriv_errors),reduced_divergence_mixed2_errors=deriv_errors,
                    nonzero_stress_components_exercised=int(nonzero_stress),viscosities=['.01','.7'],
                    native_implicit_coordinates_differentiated=True,analytic_moderate_future_fixture_only=True,
                    actual_source_history_consumed_separately=True,passed=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'waiting_physical_C2.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Waiting physical dependency changed: '+source)
        if record['waiting_physical_identities']!=waiting_physical_identities():raise ValueError('Waiting physical proof changed')
        if record['waiting_adapter_binding']!=waiting_adapter_binding():raise ValueError('Waiting adapter source changed')
        c=MPIntervalContext(); c.dps=300; finite=zeros=0
        def rowcheck(row,require_zero=False):
            nonlocal finite,zeros
            lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
            if not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Waiting physical signed row nonfinite')
            if row['exact_zero']:
                if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Invalid exact waiting zero')
                zeros+=1
            else:
                if require_zero:raise ArithmeticError('Waiting physical remainder is not exact zero')
                if any(not mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):raise ArithmeticError('Waiting physical source log nonfinite')
                finite+=1
        for point in record['samples']+[record['whole_waiting']]:
            for flag in ('actual_regional_physical_waiting_stress_remainder_identity_verified',
                         'waiting_physical_axial_viscosity_exact_zero','waiting_regional_remainder_exact_zero',
                         'physical_axial_and_time_velocity_rows_all_exact_zero',
                         'waiting_inertial_and_energy_divergence_cancelled_before_enclosure'):
                if not point[flag]:raise ValueError('Waiting physical gate missing: '+flag)
            for group in ('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_divergence_mixed2'):
                for label,grid in point[group].items():
                    if len(grid)!=(10 if group.endswith('mixed3') else 6):raise ValueError('Waiting mixed derivative coverage incomplete')
                    for row in grid.values():rowcheck(row)
            for row in point['completed_theta_theta_stress_mixed2'].values():rowcheck(row)
            for row in point['physical_angular_axial_viscosity_remainder_mixed2'].values():rowcheck(row,True)
            for label in ('exact_physical_radial_momentum_residual','exact_completed_tensor_radial_divergence',
                          'exact_physical_radial_remainder','exact_physical_axial_remainder','exact_physical_divergence'):
                if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Waiting exact structural zero lost')
                zeros+=1
            for flag in ('waiting_cone_certified','whole_outer_cone_certified','global_admissible_stress_lift_constructed',
                         'independently_bounded_global_flat_remainder','physical_energy_integral_certified',
                         'full_background_NS_validation','temporal_recursion','steep_exit_waiting_stress_mixed3_join_verified'):
                if point[flag]:raise ValueError('Waiting physical scope overclaimed: '+flag)
        fixture=moderate_cartesian_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
                    implicit_source_sha256=record['implicit_source_sha256'],all_passed=True,
                    actual_regional_physical_waiting_stress_remainder_identity_verified=True,
                    waiting_regional_remainder_exact_zero=True,actual_finite_signed_physical_rows=finite,
                    exact_physical_zero_rows=zeros,independent_native_Cartesian_fixture=fixture,
                    waiting_collar_stress_mixed3_join_verified=True,steep_exit_waiting_stress_mixed3_join_verified=False,
                    waiting_cone_certified=False,global_admissible_stress_lift_constructed=False,
                    independently_bounded_global_flat_remainder=False,physical_energy_integral_certified=False,
                    temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS actual waiting physical stress, exact-zero remainder and independent Cartesian decomposition',flush=True)
        return result


if __name__=='__main__':run()
