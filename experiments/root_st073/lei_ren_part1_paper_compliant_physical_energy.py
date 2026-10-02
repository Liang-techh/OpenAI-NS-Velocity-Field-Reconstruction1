"""Physical volume/energy map and source-bound postpulse energy domains.

Radial profile moments are not physical volume energy. The exact paper heat
tail proves finite radial energy but infinite whole-space kinetic energy for
the unlocalized source. Compact physical domains and bounded similarity
sectors are reported separately, without inserting an unauthorized cutoff.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_cartesian_field import CompliantCartesianField
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def physical_energy_weights(c,Z,delta,log_tau):
    """Multipliers of Ir and Iperp in E=pi*integral(wr*Ir+wp*Iperp)dZ."""
    Z=c.mpf(Z); d=1-Z**2; L=1-delta*Z**2
    if endpoints(d)[0]<=0:raise ValueError('Physical energy weights require |Z|<1')
    pr=(1-delta)/2; pp=(1-3*delta)/2
    return dict(radial=c.exp(pr*log_tau)*L/d**(pr+1),
                angular_axial=c.exp(pp*log_tau)*L/d**(pp+1))


def volume_jacobian(c,Z,delta,log_tau):
    Z=c.mpf(Z); d=1-Z**2; L=1-delta*Z**2
    if endpoints(d)[0]<=0:raise ValueError('Physical volume requires |Z|<1')
    log_lambda=(log_tau-c.ln(d))/2
    return 2*c.pi*c.exp((3-delta)*log_lambda)*L/d


class CompliantPhysicalEnergy:
    def __init__(self):
        self.field=CompliantCartesianField(); self.ctx=c=self.field.ctx
        self.delta=self.field.delta; self.mu=self.field.mu
        self.family=self.field.family; self.source=self.field.source; self.hashes=dict(self.field.hashes)
        checkname=PREFIX+'compliant_cartesian_field_check.json'; check=json.loads((HERE/checkname).read_bytes())
        if not check['all_passed'] or not check['accepted_outer_chart_cartesian_spatial_C4_mapped']:
            raise ValueError('Accepted actual outer Cartesian map required')
        for name,digest in check['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Physical energy dependency changed: '+name)
        self.hashes.update(check['input_hashes']); self.hashes[checkname]=hashlib.sha256((HERE/checkname).read_bytes()).hexdigest()
        if endpoints(self.delta)[0]<=0:
            raise ArithmeticError('Positive original delta required')
        if endpoints(self.delta)[1]>=mp.mpf(1)/200:raise ArithmeticError('Original small-delta domain required')
        inlet=self.field.records['flatten_mixed_C4']['whole_Z_inlet']
        self.J=read_interval(c,inlet['remaining_swirl_energy_in_Rv_Ev0_squared']['coefficients'][0])
        if endpoints(self.J)[0]<=0:raise ArithmeticError('Positive complete post-Rv radial swirl integral required')
        angularname=PREFIX+'compliant_outer_angular_repair.json'
        angular=json.loads((HERE/angularname).read_bytes())
        if angular['implicit_source_sha256']!=self.source or angular['actual_five_defect_family_sha256']!=self.family:
            raise ValueError('Gamma inverse-radius source mismatch')
        self.Scap=read_interval(c,angular['strong_inverse_radius_positive_cap'])
        self.hashes[angularname]=hashlib.sha256((HERE/angularname).read_bytes()).hexdigest()
        self.logone=read_interval(c,self.field.records['steep_waiting_C4']['waiting_log_one_minus_epsilon'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def postpulse_mass_parts(self):
        """Rv*Ev0^2: cancel +13/mu and -13/mu BEFORE interval arithmetic."""
        p=self.field.pulse; ev=self.field.Evparts
        return dict(logCstar_term=p.logRp_parts['logCstar'],
                    logPstar_term=2*p.logP+p.logRp_parts['logPstar'],
                    finite_radius_term=p.logRp_parts['finite_outer_offset'],
                    inlet_amplitude_term=ev['inlet_log']+ev['finite_offset'],
                    inverse_mu_term=self.ctx.mpf(0))

    def gamma_radial_mass_parts(self):
        """log(c_infty^2 Rb^-delta/delta), keeping exact small scales.

        Rtail=Rv exp(102+Lrel+Ts+wait), Rb=Rtail exp(3),
        theta_base=thetaT exp(-bh*wait-log(1-epsilon)).
        The wait factors combine to -delta*wait; inverse-mu terms cancel.
        """
        c=self.ctx; bp=c.mpf('.5')+self.mu; k=1-self.delta/2
        logthetaT=-100*bp-bp*self.field.Lrel-c.ln(2)-bp-(1-self.mu)/2
        logthetaT-=c.mpf('1.5')*self.field.Ts+c.mpf('1.5')-k/2
        return dict(**self.postpulse_mass_parts(),
            finite_Gamma_tail_term=2*logthetaT+102+self.field.Lrel+self.field.Ts
                -self.delta*self.field.wait-2*self.logone-3*self.delta-c.ln(self.delta))

    def report(self):
        c=self.ctx; p=(1-3*self.delta)/2; alpha=p+1
        mass=self.postpulse_mass_parts(); upperJ=c.mpf(endpoints(self.J)[1]); lowerJ=c.mpf(endpoints(self.J)[0])
        sectors={}
        for zstar in ('.5','.9','.99'):
            z=c.mpf(zstar); dmin=1-z*z
            sectors[zstar]={}
            for logtau in ('-1','-10','-100'):
                lo=c.ln(c.pi*2*z*(1-self.delta))+p*c.mpf(logtau)+c.ln(lowerJ)
                hi=c.ln(c.pi*2*z)+p*c.mpf(logtau)-alpha*c.ln(dmin)+c.ln(upperJ)
                sectors[zstar][logtau]=dict(shared_radial_mass_log_parts=mass,
                    additional_energy_log_lower=lo,additional_energy_log_upper=hi,
                    physical_axial_halfwidth_log=(1-self.delta)*(c.mpf(logtau)-c.ln(dmin))/2+c.ln(z),
                    domain='R>=Rv, |Z|<='+zstar+', all theta, fixed tau=exp('+logtau+')',
                    postpulse_velocity_has_Ur_Uz_exactly_zero=True,
                    kinetic_energy_finite_on_this_domain=True)
        cylinder={}
        for logtau in ('-1','-10','-100'):
            # |z|<=1, tau<=1, delta<1/200 imply lambda^2<3 and d>=tau/3.
            # alpha=p+1 gives tau^p*(tau/3)^-alpha=3^alpha/tau.
            cylinder[logtau]=dict(shared_radial_mass_log_parts=mass,
                additional_energy_log_upper=c.ln(2*c.pi*upperJ)+alpha*c.ln(3)-c.mpf(logtau),
                lower_bound_one_minus_Z_squared_log=c.mpf(logtau)-c.ln(3),
                domain='R>=Rv intersect physical |z|<=1, all physical radii and theta, fixed positive tau',
                finite_at_each_positive_tau=True,uniform_bound_as_tau_to_zero_certified=False)
        awaytime=dict(shared_radial_mass_log_parts=mass,
            additional_spacetime_energy_log_upper=c.ln(2*c.pi*upperJ)+alpha*c.ln(3)+c.ln(99),
            physical_domain='R>=Rv intersect |z|<=1; exp(-100)<=tau<=exp(-1)',
            proof='integral_tau_min^tau_max tau^-1 d_tau=log(tau_max/tau_min)=99',
            time_integral_finite_on_this_domain=True,time_integrability_through_tau_zero_certified=False)
        a=self.delta/2; Hfloor=1-a*(1+a)*2*self.Scap*c.exp(-3)
        if endpoints(Hfloor)[0]<=0:raise ArithmeticError('Uniform positive Gamma floor failed')
        if endpoints(alpha)[0]<=1:raise ArithmeticError('Whole-space endpoint divergence threshold not met')
        paper=HERE.parents[1]/'work_paper_cache'/'2609.35406v2.txt'
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            selected_delta=self.delta,selected_mu=self.mu,
            physical_map=dict(tau='1-t_phys',lambda_relation='sqrt(tau/(1-Z^2))',
                radial='r=lambda*sqrt(2R)',axial='z=lambda^(1-delta)*Z',
                volume='dV=2*pi*lambda^(3-delta)*(1-delta*Z^2)/(1-Z^2)*dR*dZ',
                kinetic_energy='E=pi*integral L/d*[lambda^(1-delta)*Ur^2+lambda^(1-3delta)*(Utheta^2+Uz^2)]dR dZ',
                radial_weight_exponent=(3-self.delta)/2,angular_axial_weight_exponent=alpha,
                source='paper 2.2-2.9, pp11-14'),
            exact_complete_post_Rv_radial_integral=dict(normalized_in_Rv_Ev0_squared=self.J,
                exact_positive_mass_log_parts=mass,inverse_mu_cancellation_is_exact=True,
                source='accepted original flatten inlet complete future swirl integral; Ur=Uz=0 after Rv'),
            similarity_sector_kinetic_energy_bounds=sectors,
            fixed_physical_axial_strip_postpulse_energy_bounds=cylinder,
            fixed_axial_strip_postpulse_time_integrated_energy=awaytime,
            exact_Gamma_radial_integral=dict(exact_positive_mass_log_parts=self.gamma_radial_mass_parts(),
                uniform_H_lower=Hfloor,relative_radial_integral_lower=Hfloor**2,relative_radial_integral_upper=c.mpf(1),
                formula='c_infty^2*Rb^-delta/delta times a factor in[Hdelta(2/Rb)^2,1]',
                source='paper5.4-5.10; full positive Gamma expectation retained, Rb=Rtail*exp(3)'),
            whole_space_obstruction=dict(radial_integral_has_uniform_strictly_positive_lower_bound=True,
                axial_endpoint_weight_exponent=alpha,divergence_threshold='alpha>=1; here delta<1/200 implies alpha>1',
                Z_endpoint_domain='Z=+-1 means physical z=+-infinity at each fixed tau>0',
                unlocalized_whole_space_kinetic_energy_is_infinite=True,
                no_axial_or_physical_space_cutoff_inserted=True,
                full_space_finite_energy_variant_requires_new_construction=True),
            paper_domain_notes=dict(global_finite_kinetic_energy_is_not_a_claim_of_Theorem_1_1=True,
                profile_Mp_and_Mztheta_are_not_physical_kinetic_energy=True,
                source='Theorem1.1 and1.2, pp6-7; scale cutoffs13.1-13.2 and16.7 retain leading Gamma tail;17.28-17.29 restrict endpoint comparisons'),
            physical_volume_and_kinetic_energy_functional_restored=True,
            complete_postpulse_local_physical_energy_bounds_available=True,
            full_background_physical_energy_integral_certified=False,
            global_physical_energy_integral_certified=False,
            physical_energy_integral_certified=False,core_axis_interfaces_certified=False,
            full_outer_C4_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
            next_dependency='Restore core/axis physical field and local energy including radial/axial components; specify finite physical domain and terminal-time integrability; stress/flat remainder and actual coefficient recursion',
            paper_text_path=str(paper),paper_text_sha256=hashlib.sha256(paper.read_bytes()).hexdigest(),input_hashes=self.hashes)


def run():
    with mp.workdps(270):result=CompliantPhysicalEnergy().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Physical energy: exact volume weights, complete postpulse local bounds and original whole-space Gamma obstruction generated',flush=True)
    return result


if __name__=='__main__':run()
