"""Whole inactive-gap physical source checks and independent Cartesian oracle."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_gap_physical_C2 import (
    CompliantPulseGapPhysicalC2,source_and_join_binding,FALSE_FLAGS,PREFIX,DOMAIN,source_precision)
from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import CompliantPulseGapSimilarityC4
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent


def independent_gap_Cartesian_oracle():
    """Independent physical NS oracle with complete beta history and Uz=0.

    Cartesian time, material derivative, pressure, Laplacian and completed
    tensor divergence use the original implicit lambda relation. Moderate
    fixture values test the adapter and units; they are not actual sources.
    """
    with mp.workdps(90):
        c = MPIntervalContext();c.dps = 100
        v,z,q = s.symbols('s Z q',real=True)
        mu = s.Rational(13,100);delta = s.Rational(1,2);bp = s.Rational(1,2)+mu
        B0,D,R0,Xp = map(s.Rational,('0.45','0.07','2.8','0.3'))
        H0=s.exp(-13*(1-mu)/mu)
        C=1/(1+z*z)
        mm=mp.mpf('0.13');ell=mp.mpf('.15')
        normal=mp.quad(lambda q:mp.exp(-1/(1-q*q)),[-1,0,1])
        def beta(t):
            q=t/ell
            return mp.exp(-1/(1-q*q))/(ell*normal) if abs(q)<1 else mp.mpf(0)
        controls=(s.Rational(1,10)*(1+z),-s.Rational(7,100)*(1-z*z))
        weights=[s.Float(mp.nstr(mp.quad(
            lambda q,rate=mp.mpf('.5')-i*mm:mp.exp(rate*q)*beta(q),
            [-ell,0,ell]),100),100) for i in (1,2)]
        gram=s.Float(mp.nstr(mp.quad(lambda q:mp.exp(-2*mm*q)*beta(q)**2,[-ell,0,ell]),100),100)
        M1,M2=[-sum(controls[j]*s.exp((s.Rational(1,2)-i*mu)*center)*weights[i-1]
            for j,center in enumerate((-3,-1))) for i in (1,2)]
        J0=sum(controls[j]**2*s.exp(-2*mu*center)*gram for j,center in enumerate((-3,-1)))
        Bh=s.Integer(0)
        m=M1*s.exp(-(s.Rational(1,2)-mu)*v)
        n=M2*s.exp(-(s.Rational(1,2)-2*mu)*v)
        J=J0*s.exp(2*mu*v)
        fu = s.Rational(4,5)+z*z/10+z**3/50
        e0 = fu*s.exp(2*mu*v)+(1-s.exp(2*mu*v))/(4*mu)
        pp = 1+2*mu
        P0 = -2*C*C/5-3*z*z/100
        Pc = P0*s.exp(pp*v)+C*C*(s.exp(pp*v)-1)/(2*pp)
        R = R0*s.exp(v);B = B0*s.exp(-bp*v);H = H0*s.exp(-(1-mu)*v)
        Ut = B*C;Uz = B*C*D*Bh
        Mt = s.sqrt(2)*R**s.Rational(3,2)*Ut*(1/(1-mu)+(Xp-1/(1-mu))*H)
        Mz = R*Ut*D*m
        Mtz = s.sqrt(2)*R**s.Rational(3,2)*Ut*Ut*D*n
        Mzt = R*Ut*Ut*(e0-D*D*J);P = B*B*Pc
        d,L = 1-z*z,1-delta*z*z
        transport = -R+(1-delta)*z*Mz+d*s.diff(Mz,z)
        It = Ut*transport/(L*s.sqrt(2*R))+((1-delta/2)*Mt-(1-delta)*z*s.diff(Mt,z)/2-d*s.diff(Mtz,z)+(2*delta-1)*z*Mtz)/(2*L*R)
        Iz = (transport*Uz+(1-delta)*(Mz-z*s.diff(Mz,z))/2+2*delta*z*Mzt-d*s.diff(Mzt,z)+R*(2*(1+delta)*z*P-d*s.diff(P,z)))/(L*s.sqrt(2*R))
        Tt = It+(2*s.diff(Ut,v)-Ut)/s.sqrt(2*R)
        Tz = Iz+s.sqrt(2*R)*s.diff(Uz,v)/R
        Ur = (2*z*R*Uz-(1-delta)*z*Mz-d*s.diff(Mz,z))/(L*s.sqrt(2*R))
        direct = {name:s.lambdify((v,z),expr,'mpmath',cse=True) for name,expr in
            (('Ur',Ur),('Ut',Ut),('Uz',Uz),('P',P),('Tt',Tt),('Tz',Tz))}
        # Exercise the actual new gap packet and its physical adapter.
        distance=mp.mpf('1.2');v0,z0=-distance/mm,mp.mpf('.31')
        def jet(expr):
            return IntervalTaylor(c,[c.mpf(s.lambdify(z,s.diff(expr,z,k),'mpmath')(z0)/math.factorial(k))
                for k in range(6)])
        sim=object.__new__(CompliantPulseGapSimilarityC4)
        sim.ctx=c;sim.mu=c.mpf('.13');sim.delta=c.mpf('.5');sim.Xp=c.mpf('.3')
        sim.M=[jet(M1),jet(M2)];sim.future=jet(fu);sim.J0=jet(J0);sim.P0=jet(P0)
        sim.Pin=c.mpf('.12');sim.U=c.mpf(1)
        sim.logP=c.ln(c.mpf('.45'))+13/(2*sim.mu)+13
        sim.logU=c.mpf(0);sim.logRp=c.ln(c.mpf('2.8'))-13/sim.mu
        sim.finite=c.ln(c.mpf('.07'))+1/sim.mu
        fixture=object.__new__(CompliantPulseGapPhysicalC2)
        fixture.ctx=c;fixture.similarity=sim
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
            point = fixture.gap(c.mpf(z0),c.mpf(distance),log_tau=mp.nstr(mp.log(tau),95),theta=mp.nstr(theta,95),viscosity=mp.nstr(nu,95))
            r0 = mp.sqrt(2*nu*mp.mpf(str(R0))*mp.exp(v0))*lam0
            axial0 = mp.sqrt(nu)*z0*mp.sqrt(lam0)
            def coordinates(x,y,axial,t=t0):
                rr = mp.sqrt(x*x+y*y);zt = axial/mp.sqrt(nu);ta = 1-t
                lam = (zt*zt+mp.sqrt(zt**4+4*ta))/2
                radial = rr*rr/(2*nu*lam*lam);Z = zt/mp.sqrt(lam)
                return rr,lam,mp.log(radial/mp.mpf(str(R0))),Z
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
                            if label!='axial' and abs(expected)<mp.mpf('1e-20'):raise ArithmeticError('Fixture lost nonzero gap error '+label)
                            if label=='axial' and abs(expected)>tol:raise ArithmeticError('Gap axial remainder is not zero')
                            nonzero[prefix+label]=mp.nstr(expected,20)
            # Fixed-x incompressibility, including meridional moving basis.
            x,y=r0*cs,r0*sn
            divu=mp.diff(lambda xx:fields(xx,y,axial0)[0],x)+mp.diff(lambda yy:fields(x,yy,axial0)[1],y)+mp.diff(lambda az:fields(x,y,az)[2],axial0)
            verify(prefix+'Cartesian_incompressibility',c.mpf(0),divu)
            print('Independent gap Cartesian oracle checked viscosity '+str(nu),flush=True)
        return dict(all_passed=True,checks=len(checks),identities=checks,tolerance=str(tol),
            maximum_positive_enclosure_miss=mp.nstr(worst,30),
            independent_Cartesian_time_convection_pressure_full_laplacian_and_tensor_divergence=True,
            mixed_stress_order=3,mixed_divergence_diagonal_remainder_order=2,
            fixture_errors_with_nonzero_radial_and_theta=nonzero,
            complete_original_beta_integrals_used=True,structural_axial_velocity_and_remainder_zero=True,
            fixture_parameters_are_not_actual_source_values=True)

@source_precision
def run():
    companion=CompliantPulseGapPhysicalC2()
    name=PREFIX+'pulse_gap_physical_C2.json';raw=(HERE/name).read_bytes()
    record=json.loads(raw);actual=json.loads(json.dumps(encode(pack(companion.report()))))
    if actual!=record:raise ValueError('Current whole-gap physical report differs from source')
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('Gap physical source changed: '+path)
    if record['domain']!=DOMAIN:raise ValueError('Whole original gap domain required')
    if record['source_and_physical_join_binding']!=source_and_join_binding(companion.records):
        raise ValueError('Current gap physical source join differs')
    c=MPIntervalContext();c.dps=300;finite=zeros=0
    def rowcheck(row):
        nonlocal finite,zeros
        lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
        if not all(mp.isfinite(value) for value in (lo,hi)):
            raise ArithmeticError('Nonfinite physical gap signed coefficient')
        if row['exact_zero']:
            if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Incorrect structural gap zero')
            zeros+=1
        else:
            if not all(mp.isfinite(value) for value in endpoints(read_interval(c,row['log_absolute_upper']))):
                raise ArithmeticError('Nonfinite exact-source gap physical log bound')
            finite+=1
    points=record['samples']+[record['right_end_source'],record['whole_original_gap']]
    for point in points:
        for flag in ('actual_pulse_gap_physical_decomposition_constructed',
            'gap_end_completed_physical_interface_verified','effective_radial_history_scale_is_D1',
            'full_nonzero_meridional_velocity_and_radial_remainder_retained','regional_error_not_claimed_flat',
            'gap_right_radial_remainder_not_assumed_zero'):
            if not point[flag]:raise ValueError('Missing gap physical scope '+flag)
        for group,count in (('physical_cylindrical_stress_mixed3',10),
            ('physical_cylindrical_stress_divergence_mixed2',6),
            ('physical_three_component_remainder_mixed2',6)):
            for sectors in point[group].values():
                for grid in sectors.values():
                    if len(grid)!=count:raise ValueError('Incomplete physical gap mixed rows')
                    for row in grid.values():rowcheck(row)
        for grid in point['completed_theta_theta_stress_mixed2'].values():
            if len(grid)!=6:raise ValueError('Incomplete gap diagonal mixed2')
            for row in grid.values():rowcheck(row)
        for label in ('exact_completed_tensor_radial_divergence','exact_physical_divergence'):
            if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Lost exact gap structural zero')
        if any(not row['exact_zero'] for grid in point['physical_three_component_remainder_mixed2']['axial'].values()
            for row in grid.values()):raise ArithmeticError('Uz=0 gap axial remainder is not zero')
        if not any(not row['exact_zero'] for grid in point['physical_three_component_remainder_mixed2']['radial'].values()
            for row in grid.values()):raise ArithmeticError('Nonzero gap radial history removed')
        for flag in FALSE_FLAGS:
            if point[flag] or record[flag]:raise ValueError('Gap physical scope overclaimed '+flag)
        if point['source_caps_used_as_defining_field_values']:raise ValueError('Cap chosen as gap field')
    if record['right_end_source']['s']!=-4:
        raise ArithmeticError('Formal coupled gap endpoint is not exact end s=-4')
    print('Current physical gap sources, full rows and source-functional join admitted',flush=True)
    oracle=independent_gap_Cartesian_oracle()
    hashes=dict(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=companion.family,
        implicit_source_sha256=companion.source,domain=DOMAIN,input_hashes=hashes,
        actual_pulse_gap_physical_decomposition_constructed=True,
        gap_end_completed_physical_interface_verified=True,
        actual_whole_gap_nonzero_radial_history_and_remainder_preserved=True,
        actual_gap_axial_velocity_and_remainder_structural_zero=True,
        actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
        current_whole_gap_physical_report_recomputed=True,
        admitted_full_physical_operator_identities=record['source_and_physical_join_binding']['admitted_full_physical_operator_identities'],
        actual_source_and_join_identities=len(record['source_and_physical_join_binding']['identities']),
        independent_gap_Cartesian_oracle=oracle,source_caps_used_as_defining_field_values=False,
        **{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS whole inactive-gap full physical tensor, radial remainder and completed end join; cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
