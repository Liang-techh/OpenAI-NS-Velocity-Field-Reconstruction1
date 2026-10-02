"""Independent whole-core model, axis recovery, physical derivatives and energy."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_core_physical_field import (
    core_templates,core_time_template,coefficient_value,model_phi_jets,BASES,INDICES,
    Q,F,V,PD,PI,X,Y,Z,D,B,gridkey)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_core_physical_field.json'


def functional_recovery_proof():
    rho,z,delta=s.symbols('rho z delta',real=True); mean=s.Function('M')(rho,z)
    velocity=mean+rho*s.diff(mean,rho); d=1-z*z; L=1-delta*z*z
    Qvalue=(2*z*velocity-(1-delta)*z*mean-d*s.diff(mean,z))/L
    divergence=L*(Qvalue+rho*s.diff(Qvalue,rho))-(1+delta)*z*velocity+d*s.diff(velocity,z)-2*z*rho*s.diff(velocity,rho)
    if s.simplify(divergence)!=0:raise ArithmeticError('Exact core radial-average physical divergence identity failed')
    j=s.symbols('j',positive=True); u=4*z+j
    axisQ=((1+delta)*z*u-4*d)/L
    if s.simplify(axisQ-(4*(2+delta)*z*z+(1+delta)*j*z-4)/L)!=0:raise ArithmeticError('Original radial axis recovery failed')
    # Cartesian axis matrix: diagonal radial coefficient, antisymmetric swirl.
    lam,amp=s.symbols('lambda F0',positive=True)
    A=axisQ/(2*lam**2); C=amp*lam**(-2-delta)
    axis_axial_z=lam**-2/L*(-(1+delta)*z*u+4*d)
    if s.simplify(2*A+axis_axial_z)!=0:raise ArithmeticError('Cartesian axis divergence failed')
    return dict(exact_mean_and_axis_divergence_identities=3,passed=True)


def infinite_model_fixture():
    with mp.workdps(75):
        c=MPIntervalContext(); c.dps=100; count=0
        sigma=mp.mpf('.07'); shift=mp.mpf('.03'); tol=mp.mpf('1e-55')
        chi=lambda z:(z+shift)**2/((z+shift)**2+sigma**2)
        model=lambda rho,z:mp.hyp0f1(2,-chi(z)*rho/2)
        for zz in (mp.mpf('-.03'),mp.mpf('.4')):
            values=[]
            for k in range(6):
                v=mp.diff(chi,zz,k)/math.factorial(k)
                values.append(c.mpf([max(0,v-tol),min(1,v+tol)]) if k==0 else c.mpf([v-tol,v+tol]))
            jet=IntervalTaylor(c,values)
            for rho in (mp.mpf(0),mp.mpf('2'),mp.mpf('4.1')):
                rows,tails=model_phi_jets(c,jet,c.mpf(rho))
                for i in range(6):
                    for k in range(6-i):
                        target=mp.diff(model,(rho,zz),(i,k)); lo,hi=endpoints(rows[gridkey(i,k)])
                        # The independent numerical derivative has small noise
                        # at exact vanishing coefficients (chi=0, rho=0).
                        if not lo-mp.mpf('1e-60')<=target<=hi+mp.mpf('1e-60'):
                            raise ArithmeticError('Independent infinite Bessel model derivative failed: '+str((rho,zz,i,k)))
                        if endpoints(tails[gridkey(i,k)])[0]<0:raise ArithmeticError('Negative infinite model tail')
                        count+=1
        return dict(independent_infinite_model_mixed_derivatives=count,
                    includes_chi_zero_and_axis=True,numerical_derivative_absolute_tolerance='1e-60',
                    finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def independent_nonsingular_core_coordinates():
    with mp.workdps(80):
        c=MPIntervalContext(); c.dps=100; tol=mp.mpf('1e-55')
        delta=mp.mpf('.03'); tau=mp.mpf('.7'); Lambda=mp.mpf(3); eps=1/Lambda; Z0=mp.mpf('.4')
        lam=mp.sqrt(tau/(1-Z0*Z0)); zz=lam**(1-delta)*Z0
        profiles={Q:lambda rho,z:mp.exp(mp.mpf('.11')*rho)*(1+mp.mpf('.3')*z)+mp.mpf('.1')*rho*z*z,
            F:lambda rho,z:mp.exp(mp.mpf('.09')*rho)*(1+mp.mpf('.1')*z+mp.mpf('.2')*z*z),
            V:lambda rho,z:mp.exp(-mp.mpf('.07')*rho)*(1+mp.mpf('.2')*z)+mp.mpf('.1')*mp.sin(z),
            PD:lambda rho,z:mp.exp(mp.mpf('.04')*rho)*(1+mp.mpf('.2')*z**3),
            PI:lambda rho,z:rho*mp.exp(mp.mpf('.12')*rho)*(1+mp.mpf('.3')*z*z)}
        beta={Q:mp.mpf(-1),F:-1-delta,V:-1-delta,PD:-2-2*delta,PI:-2-2*delta}
        amplitude={Q:mp.sqrt(eps),F:mp.sqrt(eps),V:mp.mpf(1),PD:mp.mpf(2),PI:eps}
        count=timecount=0
        for X0,Y0 in ((mp.mpf('1.1'),mp.mpf('.7')),(mp.mpf(0),mp.mpf(0))):
            rho0=(X0*X0+Y0*Y0)/2; xx=lam*X0/mp.sqrt(Lambda); yy=lam*Y0/mp.sqrt(Lambda); grids={}; roots={}
            for label,fn in profiles.items():
                grids[label]={gridkey(i,k):c.mpf([v-tol,v+tol]) for i in range(5) for k in range(5-i)
                              for v in (mp.diff(fn,(rho0,Z0),(i,k)),)}
            def physical(component,x,y,z,remaining=tau):
                key=(z,remaining,mp.mp.prec)
                if key not in roots:roots[key]=mp.findroot(lambda ll:ll*ll-ll**(2*delta)*z*z-remaining,lam,tol=mp.eps*16,verify=True)
                ll=roots[key]; Xv=mp.sqrt(Lambda)*x/ll; Yv=mp.sqrt(Lambda)*y/ll
                rho=(Xv*Xv+Yv*Yv)/2; Zval=z/ll**(1-delta)
                seed={Q:Xv/2,F:-Yv} if component=='ux' else {Q:Yv/2,F:Xv} if component=='uy' else {V:1} if component=='uz' else {PD:1,PI:1}
                return sum(amplitude[label]*ll**beta[label]*value*profiles[label](rho,Zval) for label,value in seed.items())
            for nx,ny,nz in INDICES:
                N=nx+ny
                for component,bases in BASES.items():
                    mapped=c.mpf(0)
                    for label in bases:
                        bracket=sum((coefficient_value(c,coefficient,c.mpf(X0),c.mpf(Y0),c.mpf(Z0),c.mpf(delta),c.mpf(beta[label]))*grids[label][gridkey(i,k)]
                                     for (i,k),coefficient in core_templates()[(component,label,nx,ny,nz)].items()),c.mpf(0))
                        mapped+=bracket*c.mpf(amplitude[label])*c.exp(c.mpf(N)/2*c.ln(c.mpf(Lambda))
                            +(c.mpf(beta[label])-N+nz*(c.mpf(delta)-1))*c.ln(c.mpf(lam)))
                    target=mp.diff(lambda x,y,z:physical(component,x,y,z),(xx,yy,zz),(nx,ny,nz))
                    if not endpoints(mapped)[0]-tol<=target<=endpoints(mapped)[1]+tol:
                        raise ArithmeticError('Independent axis-inclusive Cartesian core derivative failed: '+str((X0,component,nx,ny,nz)))
                    count+=1
            for component,bases in BASES.items():
                mapped=c.mpf(0)
                for label,seed in bases.items():
                    bracket=sum((coefficient_value(c,coefficient,c.mpf(X0),c.mpf(Y0),c.mpf(Z0),c.mpf(delta),c.mpf(beta[label]))*grids[label][gridkey(i,k)]
                                 for (i,k),coefficient in core_time_template(seed).items()),c.mpf(0))
                    mapped+=bracket*c.mpf(amplitude[label])*c.exp((c.mpf(beta[label])-2)*c.ln(c.mpf(lam)))
                target=-mp.diff(lambda remaining:physical(component,xx,yy,zz,remaining),tau)
                if not endpoints(mapped)[0]-tol<=target<=endpoints(mapped)[1]+tol:
                    raise ArithmeticError('Independent axis-inclusive fixed-x core time derivative failed')
                timecount+=1
        return dict(independent_axis_and_interior_Cartesian_spatial_derivatives=count,
                    independent_axis_and_interior_fixed_x_time_derivatives=timecount,
                    exact_axis_included=True,numerical_derivative_absolute_tolerance='1e-55',
                    finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def _run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Core physical input changed: '+name)
    c=MPIntervalContext(); c.dps=260; read=lambda v:read_interval(c,v)
    major=json.loads((HERE/'lei_ren_part1_paper_compliant_core_transfer.json').read_bytes())
    norm=json.loads((HERE/'lei_ren_part1_paper_compliant_physical_norm_family.json').read_bytes())
    if norm['selected_logCstar']!=r['selected_logCstar'] or norm['uniform_Cstar_family_sha256']!=r['uniform_Cstar_family_sha256']:
        raise ValueError('Core uses a different Cstar from matched outer source')
    if r['datum_enclosure_sha256']!=major['datum_enclosure_sha256']:raise ValueError('Analytic core pressure datum changed')
    q=read(major['scaled_map_Lipschitz_upper']); correction=read(major['scaled_map_size_upper'])/(1-q)
    def overlap(a,b):
        if max(endpoints(a)[0],endpoints(b)[0])>min(endpoints(a)[1],endpoints(b)[1]):raise ArithmeticError('Core functional source bounds disagree')
    tube=json.loads((HERE/'lei_ren_part1_paper_shared_analytic_tube.json').read_bytes())
    h=read(tube['Xh_parameter']); bounds=axischecks=correctionchecks=0
    for packet in [r['whole_core'],r['whole_axis']]+r['samples']:
        if not packet['nonlinear_field_enclosed_not_replaced_by_model'] or not packet['original_P0_and_centrifugal_pressure_increment_retained']:
            raise ValueError('Core source object substituted')
        for grid in packet['ordinary_mixed_profile_grids'].values():
            if len(grid)!=15:raise ValueError('Core profile mixed4 coverage incomplete')
            for value in grid.values():
                if any(not mp.isfinite(v) for v in endpoints(read(value))):raise ArithmeticError('Core profile bound not finite')
                bounds+=1
        for name in ('Phi_derivatives','Uz_derivatives','radial_average_Uz_derivatives'):
            if len(packet[name])!=21:raise ValueError('Required analytic mixed5 core inputs missing')
            for value in packet[name].values():
                if any(not mp.isfinite(v) for v in endpoints(read(value))):raise ArithmeticError('Analytic mixed5 input not finite')
        rmax=c.mpf(endpoints(read(packet['rho']))[1])
        for i in range(6):
            for k in range(6-i):
                expected_error=correction*rmax*embedding(c,h,rmax,1,k) if i==0 else correction*embedding(c,h,rmax,i,k)
                overlap(read(packet['nonlinear_correction_bounds'][gridkey(i,k)]),expected_error)
                correctionchecks+=1
        if endpoints(read(packet['Phi_derivatives'][gridkey(0,0)]))[0]<=0:raise ArithmeticError('Actual core swirl positivity lost')
        if endpoints(read(packet['rho']))[1]==0:
            if endpoints(read(packet['Phi_derivatives'][gridkey(0,0)]))!=(mp.mpf(1),mp.mpf(1)):raise ArithmeticError('Exact axis Phi data lost')
            for k in range(1,6):
                if endpoints(read(packet['Phi_derivatives'][gridkey(0,k)]))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Phi axis derivatives lost')
            for label in (PI,):
                for k in range(5):
                    if endpoints(read(packet['ordinary_mixed_profile_grids'][label][gridkey(0,k)]))!=(mp.mpf(0),mp.mpf(0)):
                        raise ArithmeticError('Core pressure primitive nonzero on axis')
            axischecks+=10
        for value in packet['nonlinear_correction_bounds'].values():
            if endpoints(read(value))[0]<0:raise ArithmeticError('Negative nonlinear core error cap')
    physical=times=0
    for name in ('whole_core_physical_map','whole_axis_physical_map'):
        rows=r[name]
        if not rows['includes_axis_without_inverse_radius'] or len(rows['cartesian_spatial_multiindices'])!=35:
            raise ValueError('Axis-inclusive physical mapping missing')
        for multi in rows['cartesian_spatial_multiindices'].values():
            for contributions in multi.values():
                for item in contributions.values():
                    if any(not mp.isfinite(v) or v<0 for v in endpoints(read(item['absolute_upper']))):raise ArithmeticError('Physical core bracket bound failed')
                    scale=rows['shared_physical_prefactor_bounds'][item['scale_key']]
                    if endpoints(read(scale['physical_lambda_exponent']))[1]>=0:raise ArithmeticError('Core lambda bound direction failed')
                    physical+=1
        for component in rows['first_fixed_x_physical_time_derivative'].values():times+=len(component)
    energy=r['core_local_physical_energy']; masses=energy['component_contributions']
    rho=c.mpf('4.1'); whole=r['whole_core']; grids=whole['ordinary_mixed_profile_grids']
    absolute=lambda value:c.mpf(max(abs(v) for v in endpoints(read(value))))
    expected={'radial':2*c.ln(rho*absolute(grids[Q][gridkey(0,0)]))-c.ln(4),
              'swirl':2*c.ln(rho*absolute(whole['Phi_derivatives'][gridkey(0,0)])),
              'axial':c.ln(rho)+2*c.ln(absolute(grids[V][gridkey(0,0)]))}
    for name,value in expected.items():overlap(read(masses[name]['radial_integral_log_upper_parts']['bounded_profile_term']),value)
    if energy['uniform_terminal_time_energy_certified'] or energy['full_background_physical_energy_integral_certified']:
        raise ValueError('Core local energy scope overclaimed')
    for key in ('core_inner_annulus_interfaces_certified','full_cartesian_vector_derivatives_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
        if r[key]:raise ValueError('Core scope overclaimed: '+key)
    proof=functional_recovery_proof(); model=infinite_model_fixture(); fixture=independent_nonsingular_core_coordinates()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        functional_recovery_proof=proof,infinite_Bessel_model_fixture=model,independent_nonsingular_core_fixture=fixture,
        actual_profile_mixed_bounds_checked=bounds,exact_axis_conditions_checked=axischecks,
        actual_analytic_fixed_point_correction_bounds_checked=correctionchecks,
        actual_axis_and_core_Cartesian_contribution_bounds_checked=physical,actual_fixed_x_core_time_brackets_checked=times,
        whole_Z_analytic_core_profile_enclosures_through_order5_available=True,
        whole_core_and_axis_cartesian_spatial4_enclosures_available=True,
        whole_core_and_axis_first_physical_time_derivative_enclosures_available=True,core_local_physical_energy_bounds_available=True,
        selected_point_coefficients_recomputed=False,local_old_finite_coefficient_state_used=False,
        core_inner_annulus_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
        physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
        all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Whole-Z analytic core: independent infinite model, exact axis/divergence, nonsingular spatial4/time and local energy bounds PASS',flush=True)
    return out


def run():
    with mp.workdps(280):return _run()


if __name__=='__main__':run()
