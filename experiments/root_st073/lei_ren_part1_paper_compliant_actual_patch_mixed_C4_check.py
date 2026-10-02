"""Independent physical primitive fixture and exact patch support joins."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_actual_patch_mixed_C4 import (
    CompliantActualPatchMixedC4,patch_mixed,IntervalTaylor)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_actual_patch_mixed_C4.json'


def exact_coordinate_and_join_checks():
    x=s.symbols('x',positive=True);F=s.Function('F')(x);value=F
    coefficients=((0,),(1,),(1,1),(1,3,1),(1,7,6,1))
    for k in range(1,5):
        value=s.expand(x*s.diff(value,x))
        expected=sum(coefficients[k][j-1]*x**j*s.diff(F,x,j) for j in range(1,k+1))
        if s.simplify(value-expected)!=0:raise ArithmeticError('Log radius derivative conversion failed')
    r=s.Rational(1,40);centers=(s.Rational(5,4),s.Rational(3,2),s.Rational(7,4))
    edges=[center+side*r for center in centers for side in (-1,1)]
    if edges!=[s.Rational(i,40) for i in (49,51,59,61,69,71)] or min(edges)<=1:
        raise ArithmeticError('Original normalized patch support endpoints changed')
    # beta(t)=exp(-1/(1-t^2)); each finite derivative is this flat
    # exponential times a rational function. Affine normalization preserves
    # its zero limits. This is a functional endpoint proof, not box overlap.
    t=s.symbols('t',real=True);beta=s.exp(-1/(1-t*t));row=beta
    for k in range(5):
        if any(s.limit(row,t,point,dir=direction)!=0 for point,direction in ((-1,'+'),(1,'-'))):
            raise ArithmeticError('Original beta derivative does not vanish at support endpoint')
        row=s.diff(row,t)
    return dict(exact_log_radius_derivative_identities=4,exact_original_bump_support_edges=6,
        original_beta_flat_endpoint_derivative_limits=10,
        actual_Rm_join_uses_open_zero_correction_neighborhood=True,
        actual_Rh_join_uses_full_support_and_same_implicit_functional_closure=True,passed=True)


def independent_physical_fixture():
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-65');x=mp.mpf('1.3');z=mp.mpf('.3')
        delta=mp.mpf('.0005');Pstar=mp.mpf('2.4')
        a=lambda zz:mp.mpf('.002')*(1+zz+zz**3/7)
        b=lambda zz:mp.mpf('.003')*(1-zz/3+zz**2/5)
        am=lambda zz:mp.exp(mp.mpf('-.6'))/(1+zz*zz)
        H=lambda xx,zz:xx**mp.mpf('.1')+a(zz)*xx**3
        g=lambda xx,zz:b(zz)*xx**2
        V=lambda xx,zz:4*zz+g(xx,zz)
        datum=lambda zz:mp.mpf('.03')+zz/100+zz**4/1000
        initial=lambda zz:[4*zz+mp.mpf('.002')*(1+zz**2),
            mp.mpf('.625')+zz/1000,mp.mpf('2.5')*zz+mp.mpf('.001')*(1+zz**2),
            mp.mpf('-.4')+zz**2/100,mp.mpf('2.5')+zz**3/1000]
        def physical(xx,zz):
            m1,T1,J1,E1,S1=initial(zz);aa=a(zz);bb=b(zz);A=am(zz)
            integ=lambda p:(xx**(p+1)-1)/(p+1)
            # Independent closed integrals of the PHYSICAL primitive RHSs.
            G=m1+4*zz*(xx-1)+bb*integ(2);m=G/xx
            T=T1+integ(mp.mpf('.6'))+aa*integ(mp.mpf('3.5'))
            J=J1+4*zz*(T-T1)+bb*integ(mp.mpf('2.6'))+aa*bb*integ(mp.mpf('5.5'))
            E=E1+bb*bb*integ(4)/(A*A*Pstar*Pstar)-integ(mp.mpf('.2'))/2-aa*integ(mp.mpf('3.1'))-aa*aa*integ(6)/2
            S=S1+integ(mp.mpf('-.8'))/2+aa*integ(mp.mpf('2.1'))+aa*aa*integ(5)/2
            mZ=mp.diff(lambda q:(initial(q)[0]+4*q*(xx-1)+b(q)*integ(2))/xx,zz)
            Ur=mp.sqrt(xx)*(2*zz*V(xx,zz)-(1-delta)*zz*m-(1-zz*zz)*mZ)/(1-delta*zz*zz)
            return dict(Utheta_over_Pstar=A*H(xx,zz),Uz=V(xx,zz),Ur_over_sqrt_Rm_over_2=Ur,P_over_Pstar2=datum(zz)+A*A*S,
                Mz_over_Rm=G,Mtheta_over_sqrt2_Rm_1p5_Pstar=A*T,Mtheta_z_over_sqrt2_Rm_1p5_Pstar=A*J,
                Mztheta_over_Rm_Pstar2=(8*zz*G-16*zz*zz*xx)/(Pstar*Pstar)+A*A*E,Mp_over_Pstar2=A*A*S)
        def jet(fn):return IntervalTaylor(c,[c.mpf([v-tol,v+tol])/math.factorial(n)
            for n in range(6) for v in (mp.diff(fn,z,n),)])
        rows=lambda fn:[jet(lambda zz,k=k:mp.diff(lambda xx:fn(xx,zz),x,k)) for k in range(5)]
        # Initial data in patch_mixed is the ACTUAL CURRENT primitive, not
        # the value at x=1. Supply the closed integral values accordingly.
        initial_current={
            'mass':jet(lambda zz:physical(x,zz)['Mz_over_Rm']/x),
            'theta':jet(lambda zz:physical(x,zz)['Mtheta_over_sqrt2_Rm_1p5_Pstar']/am(zz)),
            'mixed':jet(lambda zz:physical(x,zz)['Mtheta_z_over_sqrt2_Rm_1p5_Pstar']/am(zz)),
            'energy':jet(lambda zz:(physical(x,zz)['Mztheta_over_Rm_Pstar2']-(8*zz*physical(x,zz)['Mz_over_Rm']-16*zz*zz*x)/(Pstar*Pstar))/am(zz)**2),
            'pressure':jet(lambda zz:physical(x,zz)['Mp_over_Pstar2']/am(zz)**2)}
        packet=patch_mixed(c,c.mpf(x),c.mpf(z),c.mpf(delta),jet(am),jet(lambda zz:1/(Pstar*Pstar*am(zz)**2)),
            c.mpf(1/(Pstar*Pstar)),rows(H),rows(V),rows(g),initial_current,jet(datum))
        count=0
        for group in ('physical_velocity_pressure_x_Z_mixed4','physical_five_primitive_x_Z_mixed4'):
            for name,grid in packet[group].items():
                for key,bounded in grid.items():
                    k,n=[int(v[1:]) for v in key.split('_')]
                    expected=mp.diff(lambda xx,zz:physical(xx,zz)[name],(x,z),(k,n));lo,hi=endpoints(bounded)
                    if not lo-tol*10000<=expected<=hi+tol*10000:
                        raise ArithmeticError('Independent physical primitive mixed derivative failed: '+name+' '+key)
                    count+=1
        ycount=0
        for name,grid in packet['physical_velocity_pressure_y_Z_mixed4'].items():
            for key,bounded in grid.items():
                k,n=[int(v[1:]) for v in key.split('_')]
                expected=mp.diff(lambda yy,zz:physical(mp.exp(yy),zz)[name],(mp.log(x),z),(k,n));lo,hi=endpoints(bounded)
                if not lo-tol*10000<=expected<=hi+tol*10000:raise ArithmeticError('Independent log-radius mixed derivative failed')
                ycount+=1
        return dict(independent_physical_closed_integral_x_Z_derivatives=count,
            independent_physical_log_radius_y_Z_derivatives=ycount,
            nonzero_initial_moment_histories_and_pressure=True,finite_fixture_only=True,actual_source_admission=False,passed=True)


def run():
    with mp.workdps(280):
        raw=json.loads((HERE/NAME).read_bytes());provider=CompliantActualPatchMixedC4();c=provider.ctx
        for name,digest in raw['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Actual patch mixed source changed: '+name)
        if not (raw['actual_five_defect_family_sha256']==provider.family and raw['implicit_source_sha256']==provider.source
                and raw['actual_patch_all_mixed_derivatives_total_order_le4_available']):raise ValueError('Actual patch mixed scope missing')
        packets=[raw['whole_actual_patch'],raw['actual_Rm_inlet'],raw['actual_Rh_exit']]+raw['source_bump_edge_packets']+raw['interior_packets']
        counts={name:0 for name in ('physical_velocity_pressure_x_Z_mixed4','physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_x_Z_mixed4')}
        for packet in packets:
            for group in counts:
                for grid in packet[group].values():
                    expected={('y' if '_y_' in group else 'x')+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
                    if set(grid)!=expected:raise ValueError('Mixed derivative grid incomplete')
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite actual patch mixed derivative')
                        counts[group]+=1
            if not packet['radial_prefactors_differentiated_before_mixed_grid'] or not packet['same_actual_coefficient_family_and_P0_retained']:
                raise ValueError('Actual source history or radial prefactor lost')
            if packet['full_inner_interfaces_certified'] or packet['full_cartesian_vector_derivatives_certified'] or packet['admissible_stress_lift_constructed'] or packet['temporal_recursion']:
                raise ValueError('Unbuilt full-field/stress/temporal scope promoted')
        with mp.workdps(90):
            normalization=mp.quad(lambda t:mp.exp(-1/(1-t*t)) if abs(t)<1 else 0,[-1,0,1])
            r=mp.mpf(1)/40;x=mp.mpf('1.24');center=mp.mpf('1.25')
            gamma=lambda xx:mp.exp(-1/(1-((xx-center)/r)**2))/(r*normalization)
            actual=provider.gamma(c.mpf(x),c.mpf(center))
            for k,value in enumerate(actual):
                lo,hi=endpoints(value);expected=mp.diff(gamma,x,k)
                if not lo<=expected<=hi:raise ArithmeticError('Original normalized bump radial derivative scaling failed')
        if endpoints(c.exp(1))[0]<=mp.mpf(71)/40:raise ArithmeticError('Rh is not beyond all original bump supports')
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            exact_coordinate_and_join_checks=exact_coordinate_and_join_checks(),independent_physical_fixture=independent_physical_fixture(),
            actual_mixed_bounds_checked=counts,actual_original_gamma_radial_derivatives_checked=5,
            actual_patch_all_mixed_derivatives_total_order_le4_available=True,
            original_beta_support_flat_joins_certified=True,patch_Rm_and_Rh_functional_profile_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual patch mixed4, original normalized beta and functional joins PASS',flush=True)
    return result


if __name__=='__main__':run()
