"""Independent physical Jacobian/energy integral and domain obstruction checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_physical_energy import physical_energy_weights,volume_jacobian
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_physical_energy.json'


def symbolic_volume_and_mass():
    R,Z,tau=s.symbols('R Z tau',positive=True); delta=s.symbols('delta',positive=True)
    d=1-Z**2; L=1-delta*Z**2; lam=s.sqrt(tau/d)
    r=lam*s.sqrt(2*R); z=lam**(1-delta)*Z
    jac=s.diff(r,R)*s.diff(z,Z)-s.diff(r,Z)*s.diff(z,R)
    if s.simplify(s.powsimp(r*jac/lam**(3-delta))-L/d)!=0:
        raise ArithmeticError('Physical r dr dz Jacobian failed')
    # Independent scale algebra cancels both inverse-mu and waiting factors.
    lr,lp,lu,mu,Ts,Lrel,wait,logone,lt=s.symbols('lr lp lu mu Ts Lrel wait logone lt',real=True)
    bh=(1+delta)/2
    logtail=lr+13/mu+102+Lrel+Ts+wait
    logEv=lp+lu-13/(2*mu)-13
    logc=logEv+lt-bh*wait-logone+bh*logtail
    original=2*logc-delta*(logtail+3)-s.log(delta)
    recovered=lr+2*lp+2*lu-26+2*lt+102+Lrel+Ts-delta*wait-2*logone-3*delta-s.log(delta)
    if s.simplify(original-recovered)!=0:raise ArithmeticError('Gamma mass factor cancellation failed')
    if s.simplify(logtail-wait+2*logEv-(lr+2*lp+2*lu-26+102+Lrel+Ts))!=0:
        raise ArithmeticError('Post-Rv inverse-mu cancellation failed')
    # Exact physical heat velocity loses its Z dependence after substitution.
    rp,ll,dd=s.symbols('r lambda d',positive=True)
    profile_radius=rp**2/(2*ll**2)
    factor=ll**(-1-delta)*profile_radius**(-(1+delta)/2)
    if s.simplify(s.powsimp(factor)-(s.sqrt(2)/rp)**(1+delta))!=0:
        raise ArithmeticError('Physical heat amplitude scaling failed')
    if s.simplify(2*dd/profile_radius-4*(ll**2*dd)/rp**2)!=0:
        raise ArithmeticError('Physical Gamma argument failed')
    return dict(symbolic_coordinate_mass_and_heat_checks=5,passed=True)


def independent_physical_volume_integral():
    """A nonzero three-component field integrated in physical r,z coordinates."""
    with mp.workdps(60):
        delta=mp.mpf('.03'); tau=mp.mpf('.7'); za=mp.mpf('-.4'); zb=mp.mpf('.5')
        ra=mp.mpf('.8'); rb=mp.mpf('2.3'); c=MPIntervalContext(); c.dps=80; count=0
        def lam_of_Z(Z):return mp.sqrt(tau/(1-Z*Z))
        def physical_z(Z):return lam_of_Z(Z)**(1-delta)*Z
        def direct_slice(z):
            # Independent physical implicit root; no production volume weight.
            lam=mp.findroot(lambda l:l*l-l**(2*delta)*z*z-tau,mp.sqrt(tau+z*z),tol=mp.eps*16,verify=True)
            Z=z/lam**(1-delta); rlo=lam*mp.sqrt(2*ra); rhi=lam*mp.sqrt(2*rb)
            ur_square_coeff=(1+Z)**2/(2*lam**4)
            ut_square_coeff=(1-mp.mpf('.2')*Z)**2/(2*lam**(4+2*delta))
            uz_square=lam**(-2-2*delta)*(1+Z*Z)**2
            return mp.pi*((ur_square_coeff+ut_square_coeff)*(rhi**4-rlo**4)/4+uz_square*(rhi*rhi-rlo*rlo)/2)
        direct=mp.quad(direct_slice,[physical_z(za),0,physical_z(zb)])
        def transformed(Z):
            Ir=(1+Z)**2*(rb*rb-ra*ra)/2
            It=(1-mp.mpf('.2')*Z)**2*(rb*rb-ra*ra)/2
            Iz=(1+Z*Z)**2*(rb-ra)
            weights=physical_energy_weights(c,c.mpf(Z),c.mpf(delta),c.mpf(mp.log(tau)))
            # Point midpoint only for the independent moderate-parameter quadrature.
            return mp.pi*(sum(endpoints(weights['radial']))/2*Ir+sum(endpoints(weights['angular_axial']))/2*(It+Iz))
        mapped=mp.quad(transformed,[za,0,zb])
        if abs(direct-mapped)>mp.mpf('1e-48')*abs(direct):raise ArithmeticError('Independent physical kinetic integral disagrees')
        for Z in (za,mp.mpf(0),zb):
            lam=lam_of_Z(Z); L=1-delta*Z*Z; d=1-Z*Z
            jac=volume_jacobian(c,c.mpf(Z),c.mpf(delta),c.mpf(mp.log(tau)))
            target=2*mp.pi*lam**(3-delta)*L/d
            if not endpoints(jac)[0]-mp.mpf('1e-50')<=target<=endpoints(jac)[1]+mp.mpf('1e-50'):
                raise ArithmeticError('Physical volume point fixture failed')
            count+=1
        return dict(independent_physical_volume_points=count,independent_three_component_physical_energy_integrals=1,
                    relative_integral_error=mp.nstr(abs(direct-mapped)/abs(direct),12),
                    finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def original_gamma_fixture():
    with mp.workdps(40):
        a=mp.mpf('.025'); rb=mp.mpf(8); floor=1-a*(1+a)*2/rb; count=0
        for Z in (-1,mp.mpf('-.4'),0,mp.mpf('.7'),1):
            for R in (rb,2*rb,10*rb):
                xi=2*(1-Z*Z)/R
                H=mp.mpf(1) if xi==0 else xi**(-1-a)*mp.hyperu(1+a,2,1/xi)
                if not floor<=H<=1+mp.mpf('1e-35'):raise ArithmeticError('Positive Gamma floor fixture failed')
                count+=1
        return dict(independent_positive_Gamma_floor_points=count,finite_parameter_fixture_only=True,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Energy source changed: '+name)
    if hashlib.sha256(Path(r['paper_text_path']).read_bytes()).hexdigest()!=r['paper_text_sha256']:raise ValueError('Authoritative energy paper text changed')
    c=MPIntervalContext(); c.dps=260; read=lambda v:read_interval(c,v)
    delta=read(r['selected_delta']); p=(1-3*delta)/2; alpha=p+1
    if endpoints(delta)[0]<=0 or endpoints(delta)[1]>=mp.mpf(1)/200:raise ArithmeticError('Actual energy delta domain failed')
    source=json.loads((HERE/'lei_ren_part1_paper_compliant_flatten_mixed_C4.json').read_bytes())
    sourceJ=source['whole_Z_inlet']['remaining_swirl_energy_in_Rv_Ev0_squared']['coefficients'][0]
    if sourceJ!=r['exact_complete_post_Rv_radial_integral']['normalized_in_Rv_Ev0_squared']:raise ValueError('Physical energy no longer uses complete original radial future integral')
    J=read(sourceJ); mass=r['exact_complete_post_Rv_radial_integral']['exact_positive_mass_log_parts']
    if endpoints(read(mass['inverse_mu_term']))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Exact inverse-mu cancellation lost')
    count=0
    def overlap(actual,expected):
        if max(endpoints(actual)[0],endpoints(expected)[0])>min(endpoints(actual)[1],endpoints(expected)[1]):raise ArithmeticError('Energy source-bound inequality scale failed')
    for zstar,rows in r['similarity_sector_kinetic_energy_bounds'].items():
        z=c.mpf(zstar); dmin=1-z*z
        for logtau,row in rows.items():
            if row['shared_radial_mass_log_parts']!=mass:raise ValueError('Sector energy amplitude reset')
            upper=c.ln(2*c.pi*z*c.mpf(endpoints(J)[1]))+p*c.mpf(logtau)-alpha*c.ln(dmin)
            lower=c.ln(2*c.pi*z*(1-delta)*c.mpf(endpoints(J)[0]))+p*c.mpf(logtau)
            overlap(read(row['additional_energy_log_upper']),upper); overlap(read(row['additional_energy_log_lower']),lower)
            if not row['kinetic_energy_finite_on_this_domain']:raise ValueError('Accepted bounded sector energy lost')
            count+=2
    # At |z|<=1 and tau<=1, lambda>=sqrt(3) would imply tau>1.
    if endpoints(3-c.exp(delta*c.ln(3)))[0]<=1:raise ArithmeticError('Fixed physical strip lambda upper barrier failed')
    for logtau,row in r['fixed_physical_axial_strip_postpulse_energy_bounds'].items():
        overlap(read(row['additional_energy_log_upper']),c.ln(2*c.pi*c.mpf(endpoints(J)[1]))+alpha*c.ln(3)-c.mpf(logtau))
        if row['uniform_bound_as_tau_to_zero_certified']:raise ValueError('Terminal-time energy overclaimed')
        count+=1
    gamma=r['exact_Gamma_radial_integral']; floor=read(gamma['uniform_H_lower'])
    steep=json.loads((HERE/'lei_ren_part1_paper_compliant_steep_waiting_C4.json').read_bytes())
    following=json.loads((HERE/'lei_ren_part1_paper_compliant_power_angular_C4.json').read_bytes())
    angular=json.loads((HERE/'lei_ren_part1_paper_compliant_outer_angular_repair.json').read_bytes())
    thetaT_log=sum((read(value) for value in steep['whole_Z']['waiting']['inlet']['exact_relative_velocity_log_parts'].values()),c.mpf(0))
    expected_Gamma_finite=2*thetaT_log+102+read(following['postflatten_length'])+read(steep['original_steep_power_length'])
    expected_Gamma_finite-=delta*read(steep['original_refined_waiting_root'])+2*read(steep['waiting_log_one_minus_epsilon'])+3*delta+c.ln(delta)
    overlap(read(gamma['exact_positive_mass_log_parts']['finite_Gamma_tail_term']),expected_Gamma_finite)
    expected_floor=1-delta/2*(1+delta/2)*2*read(angular['strong_inverse_radius_positive_cap'])*c.exp(-3)
    overlap(floor,expected_floor)
    count+=2
    if endpoints(floor)[0]<=0 or endpoints(floor)[1]>1:raise ArithmeticError('Gamma whole-Z lower factor failed')
    overlap(read(gamma['relative_radial_integral_lower']),floor**2)
    for part in gamma['exact_positive_mass_log_parts'].values():
        if any(not mp.isfinite(v) for v in endpoints(read(part))):raise ArithmeticError('Formal Gamma energy mass no longer finite')
    if endpoints(alpha)[0]<=1 or not r['whole_space_obstruction']['unlocalized_whole_space_kinetic_energy_is_infinite']:
        raise ArithmeticError('Original Gamma global-energy obstruction hidden')
    away=r['fixed_axial_strip_postpulse_time_integrated_energy']
    overlap(read(away['additional_spacetime_energy_log_upper']),c.ln(2*c.pi*c.mpf(endpoints(J)[1]))+alpha*c.ln(3)+c.ln(99))
    if away['time_integrability_through_tau_zero_certified']:raise ValueError('Terminal spacetime energy overclaimed')
    count+=1
    for key in ('full_background_physical_energy_integral_certified','global_physical_energy_integral_certified',
                'physical_energy_integral_certified','core_axis_interfaces_certified','full_outer_C4_certified','whole_outer_cone_certified','temporal_recursion'):
        if r[key]:raise ValueError('Energy domain scope overclaimed: '+key)
    proof=symbolic_volume_and_mass(); fixture=independent_physical_volume_integral(); gamma_fixture=original_gamma_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        symbolic_checks=proof,independent_physical_integral_fixture=fixture,independent_Gamma_fixture=gamma_fixture,
        actual_local_energy_bound_checks=count,physical_volume_and_kinetic_energy_functional_restored=True,
        complete_postpulse_local_physical_energy_bounds_available=True,
        original_unlocalized_whole_space_kinetic_energy_is_infinite=True,
        full_background_physical_energy_integral_certified=False,global_physical_energy_integral_certified=False,
        physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
        all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Physical energy: independent volume/kinetic integral, original Gamma positivity, exact scale cancellation and local/global domain distinction PASS',flush=True)
    return out


if __name__=='__main__':run()
