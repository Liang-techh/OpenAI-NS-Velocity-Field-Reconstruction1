"""Independent frozen primitive ODEs, radial recovery and direction units."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_frozen_comparison_field import (
    frozen_algebra,profile_grids,IntervalTaylor,MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_frozen_comparison_field.json'


def functional_identities():
    R,r,Z=s.symbols('R r Z',positive=True)
    delta=s.symbols('delta',real=True)
    f,v,p0=[s.Function(n)(Z) for n in ('f','v','p0')]
    inlet=[s.Function(n)(Z) for n in ('mt','mz','mtz','mzt','mp')]
    mt,mz,mtz,mzt,mp0=inlet
    moments=[mt+f*(R**2-r**2),mz+v*(R-r),
             mtz+f*v*(R**2-r**2),mzt+v*v*(R-r)-f*f*(R**2-r**2)/2,
             mp0+f*f*(R-r)]
    rhs=[2*R*f,v,2*R*f*v,v*v-R*f*f,f*f]
    for q,target,base in zip(moments,rhs,inlet):
        if s.simplify(s.diff(q,R)-target)!=0 or s.simplify(q.subs(R,r)-base)!=0:
            raise ArithmeticError('Original frozen primitive RHS/exit identity failed')
    mean=moments[1]/R
    if s.simplify(mean+R*s.diff(mean,R)-v)!=0:raise ArithmeticError('Frozen mean does not recover V')
    d=1-Z**2;L=1-delta*Z**2
    Q=(2*Z*v-(1-delta)*Z*mean-d*s.diff(mean,Z))/L
    div=L*(Q+R*s.diff(Q,R))-(1+delta)*Z*v+d*s.diff(v,Z)
    if s.simplify(div)!=0:raise ArithmeticError('Recovered frozen physical divergence failed')
    rho=s.symbols('rho',nonnegative=True)
    if s.integrate(rho,(rho,0,4))!=8 or s.integrate(1,(rho,0,4))!=4:
        raise ArithmeticError('Core primitive positive integration weights changed')
    return dict(five_primitive_RHS_and_exit_identities=10,mean_recovery=1,
                structural_divergence=1,positive_radial_weight_masses=2,passed=True)


def independent_fixture():
    c=MPIntervalContext();c.dps=100;tol=mp.mpf('1e-60')
    rho,z=s.symbols('rho z',real=True);rr=s.Rational(4,5);eps=s.Rational(1,5)
    phi=1-rho/25+z/10+rho*z/100
    vel=s.Rational(3,10)+2*z+7*rho/100+rho*z*z/50
    f0=s.exp(17*z/100+3*z*z/100)
    expressions=dict(phi=phi.subs(rho,4),v=vel.subs(rho,4),
        mean=s.integrate(vel,(rho,0,4))/4,H=s.integrate(rho*phi,(rho,0,4))/8,
        K=s.integrate(rho*phi*vel,(rho,0,4))/8,A=s.integrate(vel**2,(rho,0,4))/4,
        B=s.integrate(rho*phi**2,(rho,0,4))/16,C=s.integrate(phi**2,(rho,0,4))/4,
        p0=-2-z*z/10)
    delta=mp.mpf('.03');Pstar=mp.mpf('1.4');r=mp.mpf('.8');epsilon=mp.mpf('.2')
    fn={n:s.lambdify(z,v,'mpmath') for n,v in expressions.items()};amp=s.lambdify(z,f0,'mpmath')
    R,ZZ=s.symbols('R ZZ',positive=True)
    fe=f0*expressions['phi'];ve=expressions['v']
    # Independently integrate the actual polynomial core, retaining all
    # inherited moments, and add the original frozen radial integrands.
    coremom=[2*eps**2*f0*s.integrate(rho*phi,(rho,0,4)),
        eps*s.integrate(vel,(rho,0,4)),
        2*eps**2*f0*s.integrate(rho*phi*vel,(rho,0,4)),
        eps*s.integrate(vel**2,(rho,0,4))-eps**2*f0**2*s.integrate(rho*phi**2,(rho,0,4)),
        eps*f0**2*s.integrate(phi**2,(rho,0,4))]
    raw=[coremom[0]+fe*(R**2-rr**2),coremom[1]+ve*(R-rr),
        coremom[2]+fe*ve*(R**2-rr**2),
        coremom[3]+ve**2*(R-rr)-fe**2*(R**2-rr**2)/2,
        coremom[4]+fe**2*(R-rr)]
    d=1-z*z;L=1-s.Rational(3,100)*z*z
    Aop=lambda q:(1-s.Rational(3,100))*z*q+d*s.diff(q,z)
    Pop=lambda q:2*(1+s.Rational(3,100))*z*q-d*s.diff(q,z)
    u=s.sqrt(2*R)*fe
    angular=(1-s.Rational(3,200))*raw[0]-(1-s.Rational(3,100))*z*s.diff(raw[0],z)/2
    angular-=d*s.diff(raw[2],z)-(2*s.Rational(3,100)-1)*z*raw[2]
    Itheta=fe*(-R+Aop(raw[1]))/L+angular/(2*L*R)
    pressure=s.Rational(49,25)*expressions['p0']+raw[4]
    numerator=(-R+Aop(raw[1]))*ve+(1-s.Rational(3,100))*(raw[1]-z*s.diff(raw[1],z))/2
    numerator+=2*s.Rational(3,100)*z*raw[3]-d*s.diff(raw[3],z)+R*Pop(pressure)
    targets={MTH:raw[0]/(f0*R**2),MZ:raw[1]/R,MTHZ:raw[2]/(f0*R**2),
        MZT:raw[3]/R,MP:raw[4]/(R*f0**2),'D_over_R':Itheta/(fe*R),'axial_drive':numerator/(2*L)}
    targetfn={n:[s.lambdify((R,z),s.diff(expr,z,k),'mpmath') for k in range(5)] for n,expr in targets.items()}
    physical_profiles=dict(Ur=(2*z*R*ve-(1-s.Rational(3,100))*z*raw[1]-d*s.diff(raw[1],z))/(L*s.sqrt(2*R)),
        Utheta=u,Uz=ve,P0=s.Rational(49,25)*expressions['p0'],PI=raw[4])
    mixedfn={}
    for name,expression in physical_profiles.items():
        for i in range(5):
            if i:expression=R*s.diff(expression,R)
            for k in range(5-i):mixedfn[(name,i,k)]=s.lambdify((R,z),s.diff(expression,z,k),'mpmath')
    counts=dict(moment_derivatives=0,direction_derivatives=0,mixed_profile_derivatives=0)
    for zz in (mp.mpf(0),mp.mpf('.4')):
        inputs={n:IntervalTaylor(c,[c.mpf([v-tol,v+tol])/math.factorial(k)
                for k in range(6) for v in (mp.diff(f,zz,k),)]) for n,f in fn.items()}
        inputs['F0_ratios']=[c.mpf(mp.diff(amp,zz,k)/amp(zz)) for k in range(6)]
        inputs['F0_squared_ratios']=[c.mpf(mp.diff(lambda v:amp(v)**2,zz,k)/amp(zz)**2) for k in range(6)]
        for radius in (r,mp.mpf(1),mp.mpf(4)):
            data=frozen_algebra(c,c.mpf(r/radius),c.mpf(zz),c.mpf(delta),inputs)
            def compare(interval,target):
                lo,hi=endpoints(interval)
                if not lo-tol*100<=target<=hi+tol*100:
                    raise ArithmeticError('Independent frozen comparison failed: '+str((zz,radius,target,lo,hi)))
            moments=data['moments']
            for name in (MTH,MZ,MTHZ,MP):
                for k in range(5):
                    compare(moments[name][k]*math.factorial(k),targetfn[name][k](radius,zz));counts['moment_derivatives']+=1
            # Dress split moments only when reconstructing their physical
            # derivative; testing derivative(A-R F0^2 B) requires F0^2 jets.
            fs=IntervalTaylor(c,[inputs['F0_squared_ratios'][k]/math.factorial(k) for k in range(6)])*amp(zz)**2
            joint=moments[MZT]['axial']-moments[MZT]['swirl']*fs*radius
            for k in range(5):
                compare(joint[k]*math.factorial(k),targetfn[MZT][k](radius,zz));counts['moment_derivatives']+=1
                compare(data['frozen_D_over_R'][k]*math.factorial(k),targetfn['D_over_R'][k](radius,zz))
                drive=data['axial_direction_drive_parts']
                value=(drive['hydro'][k]+drive['pressure'][k]*Pstar**2)*radius+drive['swirl'][k]*(radius**2*amp(zz)**2)
                compare(value*math.factorial(k),targetfn['axial_drive'][k](radius,zz));counts['direction_derivatives']+=2
            grids=profile_grids(c,data,c.mpf(r/radius))
            prefactors=dict(Ur=mp.sqrt(radius/2),Utheta=mp.sqrt(2*radius)*amp(zz),Uz=1,P0=Pstar**2,PI=radius*amp(zz)**2)
            for (name,i,k),target in mixedfn.items():
                compare(grids[name]['y'+str(i)+'_Z'+str(k)]*prefactors[name],target(radius,zz))
                counts['mixed_profile_derivatives']+=1
    return dict(**counts,finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def run():
    with mp.workdps(280):
        receipt=json.loads((HERE/NAME).read_bytes())
        for name,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Frozen comparison source changed: '+name)
        counts=0
        packets=[receipt[n] for n in ('whole_domain','exit_domain','terminal_domain')]+receipt['samples']
        for packet in packets:
            if not packet['actual_unsmoothed_frozen_comparison'] or packet['actual_prescribed_shear_bridge']:
                raise ValueError('Frozen field was incorrectly promoted to actual bridge')
            for grid in packet['mixed_profile_derivatives_total_order_le4'].values():
                if len(grid)!=15:raise ValueError('Frozen total-order-four grid incomplete')
                for value in grid.values():
                    lo,hi=endpoints(read_interval(mp.iv,value))
                    if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:
                        raise ArithmeticError('Nonfinite actual frozen profile enclosure')
                    counts+=1
        result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
            implicit_source_sha256=receipt['implicit_source_sha256'],
            functional_identities=functional_identities(),independent_fixture=independent_fixture(),
            actual_frozen_mixed_profile_bounds_checked=counts,
            actual_core_exit_five_primitive_axial5_enclosures_available=True,
            frozen_comparison_all_five_primitives_available=True,
            frozen_comparison_profile_mixed4_available=True,
            frozen_inertial_direction_axial4_enclosures_available=True,
            actual_smooth_comparison_installed=False,actual_prescribed_shear_bridge_installed=False,
            core_inner_annulus_interfaces_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
            all_passed=True,input_hashes={**receipt['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Five source primitives/exit/divergence and independent frozen direction units PASS',flush=True)
    return result


if __name__=='__main__':run()
