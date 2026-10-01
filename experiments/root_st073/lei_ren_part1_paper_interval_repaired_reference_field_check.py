"""Independent physical moment quadrature for the C1 annulus adapter.

Moderate-scale synthetic controls exercise nonlinear terms and axial
derivatives. This fixture is not a certificate for the production controls.
"""
import hashlib
import json
from pathlib import Path
from fractions import Fraction
import mpmath as mp
from lei_ren_part1_paper_interval_partial_five_bump_map import IntervalPartialFiveBumpMap
from lei_ren_part1_paper_interval_repaired_reference_field import evaluate_reference
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress

HERE=Path(__file__).parent

def run():
    backend=IntervalPartialFiveBumpMap();c=backend.ctx
    with mp.workdps(90):
        z0=mp.mpf('.5');rm=mp.mpf(10);x=mp.mpf('1.5');delta=mp.mpf('.01')
        hv=list(map(mp.mpf,('.001','-.002','.003','-.004','.005')))
        hd=list(map(mp.mpf,('.0002','-.0003','.0004','-.0005','.0006')))
        dv=list(map(mp.mpf,('.01','-.02','.03','-.04','.05')))
        dd=list(map(mp.mpf,('.002','-.003','.004','-.005','.006')))
        def jet(v,d):return IntervalTaylor(c,[c.mpf(v),c.mpf(d)])
        z=jet(z0,1);am=(1+z*z).reciprocal()*2;p0=z*z+3
        actual=evaluate_reference(c,z,c.mpf(delta),c.mpf(rm),am,p0,
            [jet(v,d) for v,d in zip(hv,hd)],[jet(v,d) for v,d in zip(dv,dd)],Fraction(3,2),backend)
        def raw(t):return mp.exp(-1/(1-t*t)) if abs(t)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1]);radius=mp.mpf(1)/40
        centers=list(map(mp.mpf,('1.25','1.5','1.75')))
        def bumps(y):return [raw((y-b)/radius)/(radius*norm) for b in centers]
        A=2/(1+z0*z0);AZ=-4*z0/(1+z0*z0)**2
        def integrand(y,key,order):
            beta=bumps(y);f=sum(hv[i+2]*beta[i] for i in range(3));fz=sum(hd[i+2]*beta[i] for i in range(3))
            g=hv[0]*beta[0]+hv[1]*beta[2];gz=hd[0]*beta[0]+hd[1]*beta[2]
            base=A*y**mp.mpf('.1');baseZ=AZ*y**mp.mpf('.1')
            u=base+A*f;uZ=baseZ+AZ*f+A*fz;v=4*z0+g;vZ=4+gz
            root=mp.sqrt(2*rm*y)
            vals=dict(z=g,theta=root*A*f,theta_z=root*(u*v-base*4*z0),
                z_theta=v*v-(4*z0)**2-(u*u-base*base)/2,p=(u*u-base*base)/(2*rm*y))
            ders=dict(z=gz,theta=root*(AZ*f+A*fz),theta_z=root*(uZ*v+u*vZ-baseZ*4*z0-base*4),
                z_theta=2*v*vZ-32*z0-(u*uZ-base*baseZ),p=(u*uZ-base*baseZ)/(rm*y))
            return rm*(vals if order==0 else ders)[key]
        def uncorrected(Z):
            a=2/(1+Z*Z);D=[v+(Z-z0)*d for v,d in zip(dv,dd)];R=rm*x;pow=x**mp.mpf('.1');scale=mp.sqrt(2)*rm**mp.mpf('1.5')*a
            return dict(z=4*Z*R+rm*D[0],theta=mp.mpf(5)/8*mp.sqrt(2)*R**mp.mpf('1.5')*a*pow+scale*D[2],
                theta_z=4*Z*(mp.mpf(5)/8*mp.sqrt(2)*R**mp.mpf('1.5')*a*pow)+scale*(D[1]+4*Z*D[2]),
                z_theta=16*Z*Z*R-mp.mpf(5)/12*R*a*a*pow*pow+rm*a*a*D[3]+8*Z*rm*D[0],p=mp.mpf('2.5')*a*a*pow*pow+a*a*D[4])
        expected={};count=0
        for key in ('z','theta','theta_z','z_theta','p'):
            vals=[]
            for order in (0,1):
                value=uncorrected(z0)[key] if order==0 else mp.diff(lambda Z:uncorrected(Z)[key],z0)
                for center in centers:
                    left=center-radius;right=min(x,center+radius)
                    if right>left:value+=mp.quad(lambda y:integrand(y,key,order),[left,(left+right)/2,right])
                lo,hi=endpoints(actual['physical_moments'][key][order])
                if not lo<=value<=hi:raise AssertionError(('physical integral excluded',key,order,mp.nstr(value,15)))
                vals.append(value);count+=1
            expected[key]=vals
        beta=bumps(x);f=sum(hv[i+2]*beta[i] for i in range(3));fz=sum(hd[i+2]*beta[i] for i in range(3))
        g=hv[0]*beta[0]+hv[1]*beta[2];gz=hd[0]*beta[0]+hd[1]*beta[2]
        u=A*(x**mp.mpf('.1')+f);uZ=AZ*(x**mp.mpf('.1')+f)+A*fz
        v=4*z0+g;vZ=4+gz;R=rm*x
        # x is a bump centre: all radial bump derivatives vanish exactly.
        uy=A*mp.mpf('.1')*x**mp.mpf('.1');vy=mp.mpf(0)
        P=3+z0*z0+expected['p'][0];PZ=2*z0+expected['p'][1]
        stress=evaluate_mp_stress(mp.log(R),z0,delta,Utheta=u,Uz=v,Utheta_y=uy,Utheta_Z=uZ,Uz_y=vy,Uz_Z=vZ,
            moments={k:v[0] for k,v in expected.items()},moments_Z={k:v[1] for k,v in expected.items()},P=P,P_Z=PZ,precision=90)
        for key in ('I_theta','I_z','S_theta','S_z','T_theta','T_z'):
            lo,hi=endpoints(actual['stress'][key])
            if not lo<=stress[key]<=hi:raise AssertionError(('stress excluded',key))
        for order,value in enumerate((P,PZ)):
            lo,hi=endpoints(actual['P'][order])
            if not lo<=value<=hi:raise AssertionError(('pressure excluded',order))
        lo,hi=endpoints(actual['Ur_value'])
        if not lo<=stress['U_r']<=hi:raise AssertionError('radial velocity excluded')
    result=dict(independent_physical_moment_coefficients_contained=count,independent_stress_components_contained=6,
        same_pressure_value_derivative_contained=True,radial_velocity_contained=True,synthetic_fixture_only=True,
        production_cone_or_whole_axis_certified=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_repaired_reference_field.py','lei_ren_part1_paper_interval_partial_five_bump_map.py','lei_ren_part1_paper_mp_stress.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Independent physical field fixture passed:',result,flush=True)
    return result

if __name__=='__main__':run()
