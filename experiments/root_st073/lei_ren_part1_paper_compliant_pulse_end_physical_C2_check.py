"""Independent full-meridional Cartesian NS oracle for pulse-end physical lift."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import (
    full_physical_identities, source_and_join_binding,
    pulse_velocity_rows, lift_physical_packet, FALSE_FLAGS, PREFIX)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import pulse_coefficients
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE = Path(__file__).parent


def independent_Cartesian_oracle():
    """Smooth fixture with nonzero Ur/Uz, all moments and signed memory.

    The original physical lambda relation has a closed positive root at
    delta=1/2. Cartesian time, convection, pressure, full Laplacian and
    completed tensor divergence are differentiated independently by mp.diff.
    Fixture parameters do not replace actual source parameters.
    """
    with mp.workdps(90):
        c = MPIntervalContext();c.dps = 100
        v,z,q = s.symbols('s Z q',real=True)
        mu = s.Rational(13,100);delta = s.Rational(1,2);bp = s.Rational(1,2)+mu
        B0,D,H0,R0,Xp = map(s.Rational,('0.45','0.07','0.08','2.8','0.3'))
        C = 1/(1+z*z);g = 1+v/5+v*v/20
        coefficient = s.Rational(1,10)+3*z/100+z*z/50
        Bh = coefficient*g
        def backward(rate,power=1):
            primitive = s.integrate(s.exp(rate*q)*g.subs(v,q)**power,q)
            return coefficient**power*s.exp(-rate*v)*(primitive.subs(q,0)-primitive.subs(q,v))
        m = -backward(s.Rational(1,2)-mu)
        n = -backward(s.Rational(1,2)-2*mu)
        J = backward(-2*mu,2)
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
        # Build the actual production coefficient packets with independent source jets.
        v0,z0 = mp.mpf('-1.2'),mp.mpf('.31')
        def source_rows(expr):
            rows = []
            for k in range(5):
                rowexpr = s.diff(expr,v,k)
                vals = []
                for l in range(6):
                    fn = s.lambdify((v,z),s.diff(rowexpr,z,l),'mpmath',cse=True)
                    vals.append(c.mpf(mp.nstr(fn(v0,z0)/math.factorial(l),95)))
                rows.append(IntervalTaylor(c,vals))
            return rows
        bhrows,mrows,nrows,erows,jrows,prows = [source_rows(expr) for expr in (Bh,m,n,e0,J,Pc)]
        zz = IntervalTaylor.variable(c,mp.nstr(z0,95),5)
        cc = (1+zz*zz).reciprocal()
        parts = pulse_coefficients(c.mpf('.5'),c.mpf('.13'),zz,cc,c.mpf('.3'),bhrows,mrows,nrows,erows,jrows,prows)
        logB = dict(fixture_common_reference=c.ln(c.mpf('.45'))-c.mpf('.63')*c.mpf(mp.nstr(v0,95)))
        lr = c.ln(c.mpf('2.8'))+c.mpf(mp.nstr(v0,95));ld=c.ln(c.mpf('.07'))
        lh = c.ln(c.mpf('.08'))-c.mpf('.87')*c.mpf(mp.nstr(v0,95))
        sectors = {}
        for label,values in parts.items():
            sectors[label] = {}
            for name,part in values.items():
                rp,bpower,dp,hp = part['mode']
                logs = dict(logR=rp*lr,logB=bpower*sum(logB.values()),logD=dp*ld,logH=hp*lh,normalization=-c.ln(2)/2)
                sectors[label][name] = dict(exact_source_log_parts=logs,
                    full_stress_mixed3_coefficient_enclosures={'s'+str(k)+'_Z'+str(l):row[l]*math.factorial(l)
                        for k,row in enumerate(part['full_derivative_rows']) for l in range(4-k)})
        packet = dict(Z=c.mpf(mp.nstr(z0,95)),s=c.mpf(mp.nstr(v0,95)),exact_logR=lr,
            exact_pulse_reference_logB_parts=logB,exact_logD=ld,exact_logH=lh,
            full_meridional_stress_log_sectors=sectors)
        velocity = pulse_velocity_rows(c,c.mpf('.5'),c.mpf('.13'),zz,cc,bhrows,mrows)
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
            point = lift_physical_packet(c,packet,c.mpf('.5'),velocity,mp.nstr(mp.log(tau),95),mp.nstr(theta,95),mp.nstr(nu,95))
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
                            if abs(expected)<mp.mpf('1e-20'):raise ArithmeticError('Fixture lost nonzero error '+label)
                            nonzero[prefix+label]=mp.nstr(expected,20)
            # Fixed-x incompressibility, including meridional moving basis.
            x,y=r0*cs,r0*sn
            divu=mp.diff(lambda xx:fields(xx,y,axial0)[0],x)+mp.diff(lambda yy:fields(x,yy,axial0)[1],y)+mp.diff(lambda az:fields(x,y,az)[2],axial0)
            verify(prefix+'Cartesian_incompressibility',c.mpf(0),divu)
            print('Independent full Cartesian oracle checked viscosity '+str(nu),flush=True)
        return dict(all_passed=True,checks=len(checks),identities=checks,tolerance=str(tol),
            maximum_positive_enclosure_miss=mp.nstr(worst,30),
            independent_Cartesian_time_convection_pressure_full_laplacian_and_tensor_divergence=True,
            mixed_stress_order=3,mixed_divergence_diagonal_remainder_order=2,
            fixture_nonzero_three_component_errors=nonzero,
            fixture_parameters_are_not_actual_source_values=True)


def run():
    with mp.workdps(300):
        name=PREFIX+'pulse_end_physical_C2.json'
        record=json.loads((HERE/name).read_bytes());hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Physical source changed: '+path)
        if record['full_physical_operator_proof']!=full_physical_identities():raise ValueError('Physical source proof changed')
        records={stem:json.loads((HERE/(PREFIX+stem+'.json')).read_bytes()) for stem in ('pulse_end_flatten_join','pulse_end_flatten_join_check')}
        if record['source_and_physical_join_binding']!=source_and_join_binding(records):raise ValueError('Physical join source changed')
        c=MPIntervalContext();c.dps=300;finite=zeros=0
        def rowcheck(row):
            nonlocal finite,zeros
            lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
            if not all(mp.isfinite(value) for value in (lo,hi)):raise ArithmeticError('Nonfinite signed physical coefficient')
            if row['exact_zero']:
                if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Incorrect structural zero')
                zeros+=1
            else:
                if not all(mp.isfinite(value) for value in endpoints(read_interval(c,row['log_absolute_upper']))):raise ArithmeticError('Nonfinite exact-source physical log bound')
                finite+=1
        for point in record['samples']+[record['whole_original_end']]:
            for flag in ('actual_pulse_end_physical_decomposition_constructed','pulse_end_flatten_completed_physical_interface_verified',
                'full_nonzero_meridional_velocity_and_radial_remainder_retained','regional_error_not_claimed_flat'):
                if not point[flag]:raise ValueError('Missing physical gate '+flag)
            for group,count in (('physical_cylindrical_stress_mixed3',10),('physical_cylindrical_stress_divergence_mixed2',6),
                ('physical_three_component_remainder_mixed2',6)):
                for sectors in point[group].values():
                    for grid in sectors.values():
                        if len(grid)!=count:raise ValueError('Incomplete physical mixed rows')
                        for row in grid.values():rowcheck(row)
            for grid in point['completed_theta_theta_stress_mixed2'].values():
                if len(grid)!=6:raise ValueError('Incomplete diagonal mixed2')
                for row in grid.values():rowcheck(row)
            for label in ('exact_completed_tensor_radial_divergence','exact_physical_divergence'):
                if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Lost exact structural zero')
            for flag in FALSE_FLAGS:
                if point[flag]:raise ValueError('Physical scope overclaimed '+flag)
            if point['source_caps_used_as_defining_field_values']:raise ValueError('Source cap selected as field')
        endpoint=record['samples'][-1]
        for label in ('radial','axial'):
            if any(not row['exact_zero'] for grid in endpoint['physical_three_component_remainder_mixed2'][label].values() for row in grid.values()):
                raise ArithmeticError('Actual empty-support endpoint meridional remainder not zero')
        print('Current physical source hashes, three-component sectors and exact endpoint admitted',flush=True)
        oracle=independent_Cartesian_oracle()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
            implicit_source_sha256=record['implicit_source_sha256'],
            actual_pulse_end_physical_decomposition_constructed=True,pulse_end_flatten_completed_physical_interface_verified=True,
            actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
            actual_full_physical_operator_identities=len(record['full_physical_operator_proof']['identities']),
            actual_source_and_join_identities=len(record['source_and_physical_join_binding']['identities']),
            independent_full_meridional_Cartesian_oracle=oracle,source_caps_used_as_defining_field_values=False,
            **{flag:False for flag in FALSE_FLAGS},input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS full meridional pulse-end physical tensor, three-component error and functional flatten physical join; cone/flat/global pending',flush=True)
        return result


if __name__=='__main__':
    run()
