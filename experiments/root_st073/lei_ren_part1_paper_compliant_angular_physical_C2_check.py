"""Original angular repair: axial-dependent physical transfer and Cartesian oracle."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_angular_physical_C2 import (
    CompliantAngularStressC3,angular_physical_identities,angular_adapter_binding,angular_divergence_grids)
from lei_ren_part1_paper_compliant_collar_physical_C2 import collar_velocity_bracket
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_native_cartesian_fixture():
    """Complete future and variable K, differentiated in native Cartesian space.

    A moderate K(q,Z)=exp((a-mu)*q)*(1+.05*Z+.03*Z^2) has analytic
    cumulative angular and convergent full quadratic future moments. Every
    Cartesian tensor entry, including r*Tz_z, is differentiated independently.
    K_Z and K_ZZ are nonzero. Original beta/Gamma histories are admitted
    separately by the current source receipts; this tests the general-K map
    and the new axial-dependent divergence cancellation.
    """
    with mp.workdps(60):
        a,mu,S,B0=map(mp.mpf,('.1','.35','.004','1.3'))
        delta=2*a; k=1-a; p=1+delta; bh=mp.mpf('.5')+a; b=(1-delta)/2; rate=a-mu; beta=-2-delta
        def K(q,Z):return mp.exp(rate*q)*(1+mp.mpf('.05')*Z+mp.mpf('.03')*Z*Z)
        def moment(q,Z,label):
            if label=='A':return K(q,Z)/(k+rate)
            if label=='E':return K(q,Z)**2/(delta-2*rate)
            return K(q,Z)**2/(2*(p-2*rate))
        def coordinates(r,z,tau,nu):
            zs=z/mp.sqrt(nu)
            lam=mp.findroot(lambda value:value*value-zs*zs*value**(2*delta)-tau,mp.sqrt(tau)+abs(zs))
            Z=zs/lam**(1-delta); R=r*r/(2*nu*lam*lam)
            return lam,mp.log(R*S),Z,R
        def swirl(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu)
            return mp.sqrt(nu)*lam**(-1-delta)*B0*mp.exp(-bh*q)*K(q,Z)
        def pressure(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu)
            return -nu*lam**(-2-2*delta)*B0*B0*mp.exp(-p*q)*moment(q,Z,'P')
        def stress(r,z,tau,nu):
            lam,q,Z,R=coordinates(r,z,tau,nu); B=B0*mp.exp(-bh*q); L=1-delta*Z*Z; d=1-Z*Z
            A=moment(q,Z,'A'); E=moment(q,Z,'E'); P=moment(q,Z,'P')
            Az=mp.diff(lambda value:moment(q,value,'A'),Z)
            Ez=mp.diff(lambda value:moment(q,value,'E'),Z); Pz=mp.diff(lambda value:moment(q,value,'P'),Z)
            Ct=(k*A-b*Z*Az-K(q,Z))/L+2/R*(mp.diff(lambda value:K(value,Z),q)-(1+a)*K(q,Z))
            Cz=(delta*Z*E-d*Ez/2-2*p*Z*P+d*Pz)/L
            return nu*lam**beta*mp.sqrt(R/2)*B*Ct,nu*lam**beta*mp.sqrt(R/2)*B*B*Cz
        errors=[]; divergences=[]; remainder_errors=[]; div_errors=[]; nonzero=0
        tolerance=mp.mpf('1e-36')
        for phase,nu in ((mp.mpf('.3'),mp.mpf('.01')),(mp.mpf('.75'),mp.mpf('.7'))):
            print('Independent native Cartesian axially dependent angular K sample '+str(phase)+' nu '+str(nu),flush=True)
            q=mp.mpf("-1.5")+phase; Z=mp.mpf('.4'); tau=mp.mpf('.7'); lam=mp.sqrt(tau/(1-Z*Z)); angle=mp.mpf('.6')
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
                if abs(error)>tolerance:raise ArithmeticError('Angular completed Cartesian decomposition failed: '+mp.nstr(error,12))
                errors.append(mp.nstr(error,12))
            divu=sum(partial(lambda x,y,zz:vel(x,y,zz)[j],j) for j in range(3))
            if abs(divu)>tolerance:raise ArithmeticError('Native physical divergence failed')
            divergences.append(mp.nstr(divu,12))
            # Independent full native moments above vs the actual new reduced
            # divergence code. Analytic fixture supplies its exact Iplus.
            c=MPIntervalContext(); c.dps=110
            heat=SimpleNamespace(ctx=c,a=c.mpf(a),S=c.mpf(S),delta=c.mpf(delta),k=c.mpf(k),prate=c.mpf(p))
            values=[mp.diff(lambda zz:moment(q,zz,'P'),Z,n)/math.factorial(n) for n in range(6)]; values[0]-=1/(2*p)
            Pd=IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-52'),v+mp.mpf('1e-52')]) for v in values])
            def sourcejet(fn):
                values=[mp.diff(fn,Z,n)/math.factorial(n) for n in range(6)]
                return IntervalTaylor(c,[c.mpf([v-mp.mpf('1e-52'),v+mp.mpf('1e-52')]) for v in values])
            Krows=[sourcejet(lambda zz:mp.diff(lambda qq:K(qq,zz),q,j)) for j in range(5)]
            Qrows=[sourcejet(lambda zz:mp.diff(lambda qq:K(qq,zz)**2-1,q,j)) for j in range(5)]
            if abs(mp.diff(lambda zz:K(q,zz),Z))<mp.mpf('1e-10') or abs(mp.diff(lambda zz:K(q,zz),Z,2))<mp.mpf('1e-10'):
                raise ArithmeticError('Independent angular fixture must exercise K_Z and K_ZZ')
            shape=dict(K_rows=Krows,K_squared_defect_rows=Qrows)
            grids=angular_divergence_grids(heat,dict(pressure_defect_rows=Pd),shape,c.mpf(q),c.mpf(Z))
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
                        if error>tolerance:raise ArithmeticError('Angular reduced divergence mixed2 mismatch: '+str((label,i,j,error)))
                        div_errors.append(mp.nstr(error,12))
                    value=mp.diff(lambda zz:mp.diff(lambda rr:native_error(rr,zz),r,i),z,j)
                    coefficient=-collar_velocity_bracket(c,Krows,i,j+2,c.mpf(Z),c.mpf(delta))
                    scale=nu**(mp.mpf('.5')-mp.mpf(i+j)/2)*lam**(-1-delta-i+(j+2)*(delta-1))*B*(2/R)**(mp.mpf(i)/2)
                    lo,hi=endpoints(coefficient*c.mpf(scale)); error=max(lo-value,value-hi,mp.mpf(0))
                    if error>tolerance:raise ArithmeticError('Angular axial viscosity mixed2 mismatch: '+str((i,j,error)))
                    remainder_errors.append(mp.nstr(error,12))
        return dict(Cartesian_decomposition_errors=errors,physical_divergence_errors=divergences,
                    full_completed_tensor_native_derivatives_exercised=True,nonzero_axial_viscosity_remainders_exercised=nonzero,
                    reduced_divergence_mixed2_checks=len(div_errors),reduced_divergence_mixed2_errors=div_errors,
                    nonzero_remainder_mixed2_checks=len(remainder_errors),remainder_mixed2_errors=remainder_errors,
                    numeric_tolerance=str(tolerance),viscosities=['.01','.7'],native_implicit_coordinates_differentiated=True,
                    axially_nonconstant_exponential_K_and_analytic_full_future_fixture_only=True,
                    actual_original_beta_and_Gamma_history_consumed_separately=True,passed=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'angular_physical_C2.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Angular physical source changed: '+path)
        provider=CompliantAngularStressC3()
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(provider.family,provider.source):
            raise ValueError('Angular physical source family changed')
        if record['angular_physical_identities']!=angular_physical_identities():raise ValueError('Angular physical proof changed')
        if record['angular_adapter_binding']!=angular_adapter_binding(provider):raise ValueError('Angular adapter changed')
        c=MPIntervalContext(); c.dps=300; finite=zeros=0
        def rowcheck(row):
            nonlocal finite,zeros
            lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
            if not all(mp.isfinite(v) for v in (lo,hi)):raise ArithmeticError('Angular signed physical row nonfinite')
            if row['exact_zero']:
                if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Invalid angular zero row')
                zeros+=1
            else:
                if any(not mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):raise ArithmeticError('Angular physical log bound nonfinite')
                finite+=1
        for point in record['samples']+[record['whole_angular']]+record['support_crossings']:
            for flag in ('actual_K_Z_and_K_ZZ_retained','actual_regional_physical_angular_stress_remainder_identity_verified',
                         'angular_homogeneous_A_E_and_pressure_baselines_cancelled_before_enclosure',
                         'angular_entry_physical_stress_mixed3_and_remainder_mixed2_join_verified'):
                if not point[flag]:raise ValueError('Angular physical gate missing: '+flag)
            for group,count in (('physical_cylindrical_stress_mixed3',10),('physical_cylindrical_stress_divergence_mixed2',6)):
                for grid in point[group].values():
                    if len(grid)!=count:raise ValueError('Angular physical mixed coverage incomplete')
                    for row in grid.values():rowcheck(row)
            for group in ('completed_theta_theta_stress_mixed2','physical_angular_axial_viscosity_remainder_mixed2'):
                if len(point[group])!=6:raise ValueError('Angular physical operator mixed2 coverage incomplete')
                for row in point[group].values():rowcheck(row)
            for label in ('exact_physical_radial_momentum_residual','exact_completed_tensor_radial_divergence',
                          'exact_physical_radial_remainder','exact_physical_axial_remainder','exact_physical_divergence'):
                if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Angular structural zero lost')
                zeros+=1
            for flag in ('angular_cone_certified','angular_regional_remainder_exact_zero','whole_outer_cone_certified',
                         'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
                         'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion','preceding_power_physical_stress_join_verified'):
                if point[flag]:raise ValueError('Angular physical scope overclaimed: '+flag)
        print('Current angular physical source/joins/whole signed rows admitted; independent axial-dependent Cartesian fixture running',flush=True)
        fixture=independent_native_cartesian_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],implicit_source_sha256=record['implicit_source_sha256'],
                    all_passed=True,actual_regional_physical_angular_stress_remainder_identity_verified=True,
                    actual_K_Z_and_K_ZZ_retained=True,entry_Z_independent_reduction_used=False,
                    angular_entry_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                    actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,independent_native_Cartesian_fixture=fixture,
                    angular_regional_remainder_exact_zero=False,angular_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
                    preceding_power_physical_stress_join_verified=False,physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original angular completed physical stress, axial-dependent divergence, nonzero remainder and native Cartesian decomposition',flush=True)
        return result


if __name__=='__main__':run()
