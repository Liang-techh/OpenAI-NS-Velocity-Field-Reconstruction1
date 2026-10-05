"""Whole main/exit physical source checks and independent Cartesian oracle."""
import hashlib
import json
import math
from pathlib import Path
from types import MethodType

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_main_exit_physical_C2 import (
    CompliantPulseMainExitPhysicalC2,source_and_join_binding,compiled_main_exit_lift,
    FALSE_FLAGS,PREFIX,DOMAIN,source_precision,SourceAST)
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import CompliantPulseMainExitSimilarityC4
import lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 as similarity_source
from lei_ren_part1_paper_compliant_pulse_main_exit_fixture_integrals import run as load_fixture_integrals
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent


def independent_main_Cartesian_oracle():
    """Original startup integrals and exact main plateau, with nonzero Uz.

    Smooth axial fixtures test all operators, rather than selected production
    parameters or the actual-source cone. Only the two scalar integral inputs
    in an exact replay of the current producer are replaced by point fixtures.
    """
    with mp.workdps(90):
        c=MPIntervalContext();c.dps=130
        v,z=s.symbols('y Z',real=True)
        mu=s.Rational(1,5);delta=s.Rational(1,2);bp=s.Rational(1,2)+mu
        B0,D0,Rp,Xp=map(s.Rational,('.45','.07','2.8','.3'))
        C=1/(1+z*z);xx=mu*v;p=1+2*mu
        ap=s.Rational(105,100)+3*z/100-z*z/50
        incoming=(s.Rational(1,10)*(z+z**3),s.Rational(7,100)*(z+z**3))
        energy_in=s.Rational(15,1000)+(z*z+2*z**4+z**6)/500
        J0=s.Rational(1,5)+z*z/20;P0=-7*C*C/20+z/50
        cache=load_fixture_integrals();vals=cache['values']
        constant=lambda value:s.Float(value,110)
        startup=[constant(value) for value in vals['entrance_linear']]
        start_energy=constant(vals['entrance_energy']);Kpulse=constant(vals['Kpulse'])
        xleft=s.Rational(1,50)
        primitive=lambda a:-s.exp(-2*a)*((a-s.Rational(1,100))**2/2
            +(a-s.Rational(1,100))/2+s.Rational(1,4))
        partial_energy=start_energy+primitive(xx)-primitive(xleft)
        kernels=[]
        for i in (1,2):
            rate=s.Rational(1,2)-i*mu
            head=startup[i-1]/mu-s.exp(rate*xleft/mu)*(
                (xleft-s.Rational(1,100))/rate-mu/rate**2)
            kernels.append((xx-s.Rational(1,100))/rate-mu/rate**2+s.exp(-rate*v)*head)
        m,n=[ap*kernels[i-1]+s.exp(-(s.Rational(1,2)-i*mu)*v)*incoming[i-1] for i in (1,2)]
        Bh=ap*(xx-s.Rational(1,100))
        energy=s.exp(2*xx)*(energy_in+ap*ap*partial_energy/mu-(1-s.exp(-2*xx))/(4*mu))
        future=s.exp(26)/mu*(ap*ap*Kpulse-(1-s.exp(-26))/4+mu*energy_in)+D0*D0*J0
        B=B0*s.exp(-bp*v);R=Rp*s.exp(v);H=s.exp(-(1-mu)*v)
        Ut=B*C;Uz=Ut*Bh
        Mt=s.sqrt(2)*R**s.Rational(3,2)*Ut*(1/(1-mu)+(Xp-1/(1-mu))*H)
        Mz=R*Ut*m;Mtz=s.sqrt(2)*R**s.Rational(3,2)*Ut*Ut*n;Mzt=R*Ut*Ut*energy
        P=B*B*(-C*C/(2*p)+s.exp(-p*(13-xx)/mu)*(P0+C*C/(2*p)))
        d,L=1-z*z,1-delta*z*z
        transport=-R+(1-delta)*z*Mz+d*s.diff(Mz,z)
        It=Ut*transport/(L*s.sqrt(2*R))+((1-delta/2)*Mt
            -(1-delta)*z*s.diff(Mt,z)/2-d*s.diff(Mtz,z)+(2*delta-1)*z*Mtz)/(2*L*R)
        Iz=(transport*Uz+(1-delta)*(Mz-z*s.diff(Mz,z))/2+2*delta*z*Mzt
            -d*s.diff(Mzt,z)+R*(2*(1+delta)*z*P-d*s.diff(P,z)))/(L*s.sqrt(2*R))
        Tt=It+(2*s.diff(Ut,v)-Ut)/s.sqrt(2*R)
        Tz=Iz+s.sqrt(2*R)*s.diff(Uz,v)/R
        Ur=(2*z*R*Uz-(1-delta)*z*Mz-d*s.diff(Mz,z))/(L*s.sqrt(2*R))
        direct={name:s.lambdify((v,z),expr,'mpmath',cse=True) for name,expr in
            (('Ur',Ur),('Ut',Ut),('Uz',Uz),('P',P),('Tt',Tt),('Tz',Tz))}
        xi0,z0=mp.mpf(1),mp.mpf('.31');v0=xi0/mp.mpf('.2')
        def jet(expr):
            return IntervalTaylor(c,[c.mpf(s.lambdify(z,s.diff(expr,z,k),'mpmath')(z0)/math.factorial(k))
                for k in range(6)])
        sim=object.__new__(CompliantPulseMainExitSimilarityC4)
        sim.ctx=c;sim.mu=c.mpf('.2');sim.delta=c.mpf('.5');sim.Xp=c.mpf('.3');sim.cells=128
        sim.ap=jet(ap);sim.incoming=[jet(expr) for expr in incoming]
        sim.future=jet(future);sim.J0=jet(J0);sim.P0=jet(P0)
        sim.Pin=c.mpf('.12');sim.logP=c.ln(c.mpf('.45'));sim.logU=c.mpf(0)
        sim.logRp=c.ln(c.mpf('2.8'));sim.finite=c.ln(c.mpf('.07'))+1/sim.mu
        actual_integrals=[c.mpf(s.lambdify(v,expr,'mpmath')(v0)) for expr in kernels]
        remaining=c.mpf(s.lambdify(v,Kpulse-partial_energy,'mpmath')(v0))
        def exact_kernel(ctx,mm,x,rate,cells):
            i=0 if endpoints(rate)[0]>mp.mpf('.2') else 1
            return dict(enclosure=actual_integrals[i])
        asts=SourceAST();env=dict(vars(similarity_source))
        env.update(partial_linear_kernel=exact_kernel,partial_future_energy=lambda *args:remaining)
        replay=asts.replay('pulse_main_exit_similarity_C4','main_exit',env)
        sim.main_exit=MethodType(replay,sim)
        fixture=object.__new__(CompliantPulseMainExitPhysicalC2)
        fixture.ctx=c;fixture.similarity=sim;fixture.lift,_=compiled_main_exit_lift()
        checks = {};worst = mp.mpf(0);tol = mp.mpf('1e-55')
        def verify(name,actual,value):
            nonlocal worst
            lo,hi = endpoints(actual);miss = max(lo-value,value-hi,mp.mpf(0));worst=max(worst,miss)
            if not mp.isfinite(value) or miss>tol*max(1,abs(value)):
                raise ArithmeticError('Full Cartesian oracle mismatch '+name+': '+mp.nstr(miss,12))
            checks[name] = True
        theta = mp.mpf('.37');cs,sn=mp.cos(theta),mp.sin(theta)
        tau = mp.mpf('.8');t0 = 1-tau;lam0=mp.sqrt(tau/(1-z0*z0))
        nonzero = {}
        for nu in (mp.mpf('.01'),mp.mpf('.7')):
            point = fixture.main_exit(c.mpf(z0),c.mpf(xi0),log_tau=mp.nstr(mp.log(tau),95),theta=mp.nstr(theta,95),viscosity=mp.nstr(nu,95))
            r0 = mp.sqrt(2*nu*mp.mpf(str(Rp))*mp.exp(v0))*lam0
            axial0 = mp.sqrt(nu)*z0*mp.sqrt(lam0)
            def coordinates(x,y,axial,t=t0):
                rr = mp.sqrt(x*x+y*y);zt = axial/mp.sqrt(nu);ta = 1-t
                lam = (zt*zt+mp.sqrt(zt**4+4*ta))/2
                radial = rr*rr/(2*nu*lam*lam);Z = zt/mp.sqrt(lam)
                return rr,lam,mp.log(radial/mp.mpf(str(Rp))),Z
            def fields(x,y,axial,t=t0):
                rr,lam,sv,Z = coordinates(x,y,axial,t);cosine,sine=x/rr,y/rr
                ur=mp.sqrt(nu)/lam*direct['Ur'](sv,Z)
                ut=mp.sqrt(nu)*lam**(-mp.mpf('1.5'))*direct['Ut'](sv,Z)
                uz=mp.sqrt(nu)*lam**(-mp.mpf('1.5'))*direct['Uz'](sv,Z)
                pressure=nu*lam**(-3)*direct['P'](sv,Z)
                return (cosine*ur-sine*ut,sine*ur+cosine*ut,uz,pressure)
            def physical_pair(x,y,axial):
                rr,lam,sv,Z = coordinates(x,y,axial)
                return nu*lam**(-mp.mpf('2.5'))*direct['Tt'](sv,Z),nu*lam**(-mp.mpf('2.5'))*direct['Tz'](sv,Z)
            def tensor(x,y,axial):
                rr=mp.sqrt(x*x+y*y);cosine,sine=x/rr,y/rr
                tt,tz=physical_pair(x,y,axial)
                diagonal=rr*mp.diff(lambda zphys:physical_pair(x,y,zphys)[1],axial)
                return ((-2*cosine*sine*tt+sine*sine*diagonal,(cosine*cosine-sine*sine)*tt-cosine*sine*diagonal,cosine*tz),
                    ((cosine*cosine-sine*sine)*tt-cosine*sine*diagonal,2*cosine*sine*tt+cosine*cosine*diagonal,sine*tz),
                    (cosine*tz,sine*tz,mp.mpf(0)))
            def divtensor(x,y,axial,component):
                values=(x,y,axial)
                result=mp.mpf(0)
                for axis in range(3):
                    def axisfn(value):
                        coords=list(values);coords[axis]=value
                        return tensor(*coords)[component][axis]
                    result+=mp.diff(axisfn,values[axis])
                return result
            def residual(x,y,axial,component):
                coords=(x,y,axial);base=fields(*coords)
                value=mp.diff(lambda time:fields(*coords,t=time)[component],t0)
                for axis in range(3):
                    def velocity_axis(position):
                        args=list(coords);args[axis]=position
                        return fields(*args)[component]
                    def pressure_axis(position):
                        args=list(coords);args[axis]=position
                        return fields(*args)[3]
                    value+=base[axis]*mp.diff(velocity_axis,coords[axis])
                    value-=nu*mp.diff(velocity_axis,coords[axis],2)
                    if axis==component:value+=mp.diff(pressure_axis,coords[axis])
                return value
            def cylindrical_div(rr,axial,label):
                x,y=rr*cs,rr*sn
                dx,dy,dz=(divtensor(x,y,axial,k) for k in range(3))
                return dict(radial=cs*dx+sn*dy,theta=-sn*dx+cs*dy,axial=dz)[label]
            def cylindrical_error(rr,axial,label):
                x,y=rr*cs,rr*sn
                ex,ey,ez=(residual(x,y,axial,k)+divtensor(x,y,axial,k) for k in range(3))
                return dict(radial=cs*ex+sn*ey,theta=-sn*ex+cs*ey,axial=ez)[label]
            def materialize(row):
                coefficient=row['signed_coefficient']
                return coefficient*c.exp(sum(row['actual_source_log_parts'].values(),c.mpf(0))
                    +row['radial_log_prefactor']+row['physical_lambda_exponent']*c.ln(c.mpf(mp.nstr(lam0,95)))
                    +row['physical_viscosity_exponent']*c.ln(c.mpf(mp.nstr(nu,95))))
            prefix='nu'+mp.nstr(nu,3)+'/'
            for label,component in (('theta',0),('axial',1)):
                for i in range(4):
                    for j in range(4-i):
                        key='r'+str(i)+'_z'+str(j)
                        actual=sum((materialize(grid[key]) for grid in point['physical_cylindrical_stress_mixed3'][label].values()),c.mpf(0))
                        expected=mp.diff(lambda rr,az:physical_pair(rr*cs,rr*sn,az)[component],(r0,axial0),(i,j))
                        verify(prefix+'stress/'+label+'/'+key,actual,expected)
            for i in range(3):
                for j in range(3-i):
                    key='r'+str(i)+'_z'+str(j)
                    actual=sum((materialize(grid[key]) for grid in point['completed_theta_theta_stress_mixed2'].values()),c.mpf(0))
                    expected=mp.diff(lambda rr,az:tensor(rr*cs,rr*sn,az)[0][0]*sn*sn-2*tensor(rr*cs,rr*sn,az)[0][1]*cs*sn
                        +tensor(rr*cs,rr*sn,az)[1][1]*cs*cs,(r0,axial0),(i,j))
                    verify(prefix+'diagonal/'+key,actual,expected)
                    for label in ('theta','axial'):
                        actual=sum((materialize(grid[key]) for grid in point['physical_cylindrical_stress_divergence_mixed2'][label].values()),c.mpf(0))
                        verify(prefix+'Cartesian_div/'+label+'/'+key,actual,mp.diff(lambda rr,az:cylindrical_div(rr,az,label),(r0,axial0),(i,j)))
                    for label in ('radial','theta','axial'):
                        actual=sum((materialize(grid[key]) for grid in point['physical_three_component_remainder_mixed2'][label].values()),c.mpf(0))
                        expected=mp.diff(lambda rr,az:cylindrical_error(rr,az,label),(r0,axial0),(i,j))
                        verify(prefix+'Cartesian_NS_error/'+label+'/'+key,actual,expected)
                        if i==j==0:
                            if abs(expected)<mp.mpf('1e-20'):raise ArithmeticError('Fixture lost nonzero main error '+label)
                            nonzero[prefix+label]=mp.nstr(expected,20)
            # Fixed-x incompressibility, including meridional moving basis.
            x,y=r0*cs,r0*sn
            divu=mp.diff(lambda xx:fields(xx,y,axial0)[0],x)+mp.diff(lambda yy:fields(x,yy,axial0)[1],y)+mp.diff(lambda az:fields(x,y,az)[2],axial0)
            verify(prefix+'Cartesian_incompressibility',c.mpf(0),divu)
            print('Independent main Cartesian oracle checked viscosity '+str(nu),flush=True)
        return dict(all_passed=True,checks=len(checks),identities=checks,tolerance=str(tol),
            maximum_positive_enclosure_miss=mp.nstr(worst,30),
            independent_Cartesian_time_convection_pressure_full_laplacian_and_tensor_divergence=True,
            mixed_stress_order=3,mixed_divergence_diagonal_remainder_order=2,
            fixture_errors_with_all_three_nonzero_components=nonzero,
            original_startup_integrals_and_exact_plateau_primitives_used=True,
            full_nonzero_axial_shear_and_incoming_radial_cross_products_retained=True,
            exact_point_fixture_integrals_replayed_into_current_source_AST=True,
            fixture_parameters_are_not_actual_source_values=True)

@source_precision
def run():
    companion=CompliantPulseMainExitPhysicalC2()
    name=PREFIX+'pulse_main_exit_physical_C2.json';raw=(HERE/name).read_bytes()
    record=json.loads(raw);actual=json.loads(json.dumps(encode(pack(companion.report()))))
    if actual!=record:raise ValueError('Current main/exit physical report differs from source')
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('Main/exit physical source changed: '+path)
    if record['domain']!=DOMAIN:raise ValueError('Original main/exit domain required')
    if record['source_and_physical_join_binding']!=source_and_join_binding(companion.records):
        raise ValueError('Current main/exit physical source join differs')
    c=MPIntervalContext();c.dps=300;finite=zeros=0
    def rowcheck(row):
        nonlocal finite,zeros
        lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
        if not all(mp.isfinite(value) for value in (lo,hi)):
            raise ArithmeticError('Nonfinite physical signed coefficient')
        if row['exact_zero']:
            if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Incorrect structural zero')
            zeros+=1
        else:
            if not all(mp.isfinite(value) for value in endpoints(read_interval(c,row['log_absolute_upper']))):
                raise ArithmeticError('Nonfinite physical log bound')
            finite+=1
    for pointkey in ('whole_original_main','whole_original_exit','common_main_exit','original_exit_gap'):
        point=record[pointkey]
        for flag in ('actual_pulse_main_exit_physical_decomposition_constructed',
            'main_exit_and_exit_gap_completed_physical_interfaces_verified',
            'all_local_incoming_nonlinear_history_cross_products_retained',
            'frozen_source_logs_not_differentiated_again','regional_error_not_claimed_flat'):
            if not point[flag]:raise ValueError('Missing physical scope '+flag)
        for group,order in (('physical_cylindrical_stress_mixed3',3),
            ('physical_cylindrical_stress_divergence_mixed2',2),
            ('physical_three_component_remainder_mixed2',2)):
            required={'r'+str(i)+'_z'+str(j) for i in range(order+1) for j in range(order+1-i)}
            for sectors in point[group].values():
                for grid in sectors.values():
                    if set(grid)!=required:raise ValueError('Incomplete physical mixed rows')
                    for row in grid.values():rowcheck(row)
        required={'r'+str(i)+'_z'+str(j) for i in range(3) for j in range(3-i)}
        for grid in point['completed_theta_theta_stress_mixed2'].values():
            if set(grid)!=required:raise ValueError('Incomplete diagonal mixed2')
            for row in grid.values():rowcheck(row)
        for label in ('exact_completed_tensor_radial_divergence','exact_physical_divergence'):
            if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Lost exact structural zero')
        for flag in FALSE_FLAGS:
            if point[flag] or record[flag]:raise ValueError('Physical scope overclaimed '+flag)
        if point['source_caps_used_as_defining_field_values']:raise ValueError('Cap chosen as field')
    end=record['original_exit_gap']['physical_three_component_remainder_mixed2']
    if any(not row['exact_zero'] for grid in end['axial'].values() for row in grid.values()):
        raise ArithmeticError('Flat xi11 axial remainder is not zero')
    if not any(not row['exact_zero'] for grid in end['radial'].values() for row in grid.values()):
        raise ArithmeticError('Nonzero xi11 radial history was removed')
    print('Current main/exit sources, complete rows and functional joins admitted',flush=True)
    oracle=independent_main_Cartesian_oracle()
    hashes=dict(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for stem in ('pulse_main_exit_fixture_integrals.py','pulse_main_exit_fixture_integrals.json'):
        path=PREFIX+stem;hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=companion.family,
        implicit_source_sha256=companion.source,domain=DOMAIN,input_hashes=hashes,
        actual_pulse_main_exit_physical_decomposition_constructed=True,
        main_exit_and_exit_gap_completed_physical_interfaces_verified=True,
        actual_main_exit_full_meridional_history_and_three_component_remainder_preserved=True,
        actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
        current_whole_main_exit_physical_report_recomputed=True,
        admitted_full_physical_operator_identities=record['source_and_physical_join_binding']['admitted_full_physical_operator_identities'],
        actual_source_and_join_identities=len(record['source_and_physical_join_binding']['identities']),
        independent_main_Cartesian_oracle=oracle,source_caps_used_as_defining_field_values=False,
        **{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS original main/exit full physical tensor, all remainder components and completed gap join; cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
