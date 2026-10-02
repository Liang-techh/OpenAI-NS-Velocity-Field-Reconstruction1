"""Resolve the admitted axis-amplitude peak without adding tiny root offsets.

The exact shared H(a)=0 is imposed before interval enclosure. Finite scaled
K=Lambda*G is integrated as a polynomial numerator plus a directed positive-
denominator remainder. This is a local F0 envelope, not temporal recursion or
a measured whole-vortex aspect ratio. The absolute F0 remains unmaterialized.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_anchored_axis_amplitude import CompliantAnchoredAxisAmplitude
from lei_ren_part1_paper_compliant_core_physical_field import symmetric,relative_amplitude_jet
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


class CompliantRootCenteredPeak:
    def __init__(self):
        self.amplitude=CompliantAnchoredAxisAmplitude();self.core=self.amplitude.core;self.ctx=c=self.core.ctx
        self.hashes=dict(self.amplitude.hashes)
        name=PREFIX+'anchored_axis_amplitude_check.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or receipt['implicit_source_sha256']!=self.core.source:
            raise ValueError('Accepted current-source anchored amplitude required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Anchored amplitude prerequisite changed: '+path)
        self.hashes.update(receipt['input_hashes'])
        for path in (name,Path(__file__).name,'lei_ren_part1_paper_compliant_core_uniform_bounds.py',
                     'lei_ren_part1_paper_compliant_core_physical_field.py','lei_ren_part1_paper_interval_taylor.py'):
            self.hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
        self.anchor=self.amplitude.roots['anchor_interval'];self.delta=self.core.delta;self.sigma=self.core.sigma
        self.anchor_rows=None
        self.initialize_chart()

    def initialize_chart(self):
        c=self.ctx
        # epsilon is the SAME exact reciprocal of Lambda. Source intervals
        # are propagated conservatively, without selecting their midpoints.
        self.epsilon=1/self.core.Lambda
        self.b=self.sigma*c.sqrt(self.epsilon)
        self.log_b=c.ln(self.sigma)-self.core.logLambda/2
        self.hp=(9-self.delta)/2-12*self.anchor**2-2*self.core.j*self.anchor
        self.q1=-12*self.anchor-self.core.j;self.La=1-self.delta*self.anchor**2
        self.domain=c.mpf(4);self.offset_upper=self.b*self.domain
        self.q_change=abs(self.q1)*self.offset_upper+4*self.offset_upper**2
        self.L_change=abs(self.delta)*(2*abs(self.anchor)*self.offset_upper+self.offset_upper**2)
        self.qmax=abs(self.hp)+self.q_change;self.Lmax=abs(self.La)+self.L_change
        if endpoints(self.hp-self.q_change)[0]<=0 or endpoints(self.La-self.L_change)[0]<=0:
            raise ArithmeticError('Root-centered numerator is not positive')
        if endpoints(abs(self.anchor)+self.offset_upper)[1]>=1 or endpoints(self.b)[0]<=0:
            raise ArithmeticError('Positive width in original finite axial domain required')
        # Exact L(a+b*s)*q(b*s), before dividing by 1+epsilon*s^2*q^2.
        a,b,d,A,L,h=self.anchor,self.b,self.delta,self.q1,self.La,self.hp
        self.numerator=[L*h,b*(L*A-2*d*a*h),b**2*(-4*L-2*d*a*A-d*h),
                        b**3*(8*d*a-d*A),4*d*b**4]
        self.denominator_remainder_coefficient=self.epsilon*self.Lmax*self.qmax**3/4
        self.gaussian_error_coefficient=(self.L_change*self.qmax+abs(self.La)*self.q_change
                                        +self.epsilon*self.domain**2*self.Lmax*self.qmax**3)
        self.gaussian_uniform_error=self.domain**2*self.gaussian_error_coefficient/2

    def derivative_jets(self,xi,order=4):
        c=self.ctx;s=IntervalTaylor.variable(c,xi,order);one=IntervalTaylor.constant(c,1,order)
        q=one*self.hp+s*(self.b*self.q1)-s*s*(4*self.b**2)
        L=one*self.La-s*(2*self.delta*self.anchor*self.b)-s*s*(self.delta*self.b**2)
        denominator=one+s*s*q*q*self.epsilon
        coefficients=list(denominator.coefficients)
        coefficients[0]=1+self.epsilon*c.mpf(xi)**2*q[0]**2
        return (s*L*q/IntervalTaylor(c,coefficients)).coefficients

    def evaluate(self,xi):
        c=self.ctx;s=c.mpf(xi)
        if not all(mp.isfinite(v) for v in endpoints(s)) or endpoints(abs(s))[1]>4:
            raise ValueError('Finite scaled root offset xi in [-4,4] required')
        polynomial=sum((coefficient*s**(k+2)/(k+2) for k,coefficient in enumerate(self.numerator)),c.mpf(0))
        remainder=self.denominator_remainder_coefficient*abs(s)**4
        # L and q are positive on the entire real chart. Omitting the
        # positive denominator always increases the anchored integral on
        # either side (the negative path reverses orientation as well).
        raw=polynomial-c.mpf([0,endpoints(remainder)[1]])
        K=c.mpf([max(mp.mpf(0),endpoints(raw)[0]),endpoints(raw)[1]])
        if endpoints(s)==(mp.mpf(0),mp.mpf(0)):K=c.mpf(0)
        normalized=c.exp(-K);gradient=self.derivative_jets(s)
        relative=relative_amplitude_jet(c,[-value for value in gradient])
        return dict(xi=s,coordinate_definition='Z=a+b*xi, shared H(a)=0, b=sigma/sqrt(Lambda)',
            root_interval=self.anchor,offset_as_separate_positive_width_times_xi=dict(width=self.b,xi=s),
            actual_scaled_G=K,integrated_numerator_polynomial=polynomial,
            denominator_remainder_abs_upper=remainder,
            gaussian_leading=self.La*self.hp*s**2/2,
            gaussian_error_abs_upper=self.gaussian_error_coefficient*s**2/2,
            normalized_F0_relative_to_shared_anchor=normalized,
            log_F0_relative_to_shared_anchor=-K,
            absolute_log_F0_as_separate_terms=[-self.core.logC,-K],
            K_derivative_orders=list(range(1,6)),
            K_derivatives=[math.factorial(k)*value for k,value in enumerate(gradient)],
            normalized_F0_derivatives=[normalized*value for value in relative],
            exact_anchor_zero=endpoints(s)==(mp.mpf(0),mp.mpf(0)),
            root_offset_not_added_to_anchor=True,selected_logCstar_preserved=True,
            relative_amplitude_is_not_absolute_F0=True,absolute_F0_not_materialized=True,
            actual_source_uniform_enclosure=True,measured_blowup_dynamics=False,temporal_recursion=False)

    def level_width(self,log_drop):
        """Uniform certified left/right locations of F0/F0(a)=exp(-drop)."""
        c=self.ctx;drop=c.mpf(log_drop);E=self.gaussian_uniform_error
        if endpoints(drop-E)[0]<=0:raise ValueError('Positive level outside the enclosure error required')
        p=self.numerator[0]
        # A directed global local-chart comparison encloses BOTH root sides.
        lower=c.sqrt(2*(drop-E)/c.mpf(endpoints(p)[1]))
        upper=c.sqrt(2*(drop+E)/c.mpf(endpoints(p)[0]))
        # Small explicit enlargement creates strict source-uniform signs;
        # it is reported as bracket tolerance, not absorbed into source data.
        tolerance=c.mpf('1e-170')
        low=c.mpf(endpoints(lower-tolerance)[0]);high=c.mpf(endpoints(upper+tolerance)[1])
        if endpoints(low)[0]<=0 or endpoints(high)[1]>4:raise ValueError('Level outside local chart')
        sides={}
        for side in (-1,1):
            low_K=self.evaluate(side*low)['actual_scaled_G'];high_K=self.evaluate(side*high)['actual_scaled_G']
            if not (endpoints(low_K)[1]<endpoints(drop)[0] and endpoints(high_K)[0]>endpoints(drop)[1]):
                raise ArithmeticError('Strict source-uniform peak level signs failed')
            sides['left' if side<0 else 'right']=dict(sign=side,
                absolute_xi_root_interval=c.mpf([endpoints(low)[0],endpoints(high)[1]]),
                inner_endpoint_K=low_K,outer_endpoint_K=high_K,
                monotone_from_anchor_proved_by_positive_L_and_q=True,
                unique_level_for_each_admitted_source=True)
        xi_width=sides['left']['absolute_xi_root_interval']+sides['right']['absolute_xi_root_interval']
        # Integral mean of dz/dZ without subtracting near-coincident z values.
        zcover=self.anchor+symmetric(c,self.offset_upper)
        D=1-zcover**2;L=1-self.delta*zcover**2
        J=L/c.power(D,(3-self.delta)/2)
        anchor_D=1-self.anchor**2
        log_axial_prefactor=self.log_b+c.ln(xi_width)+c.ln(J)
        log_radial_core_radius_prefactor=(c.ln(c.mpf('8.2'))-self.core.logLambda-c.ln(anchor_D))/2
        # Cancel sqrt(Lambda) symbolically, before subtracting huge logs.
        log_aspect_prefactor=c.ln(self.sigma)+c.ln(xi_width)+c.ln(J)+(c.ln(anchor_D)-c.ln(c.mpf('8.2')))/2
        return dict(log_drop=drop,normalized_level=c.exp(-drop),sides=sides,
            xi_full_width=xi_width,bracket_enlargement_each_endpoint=tolerance,
            physical_axial_integral_mean_factor=J,log_axial_width_prefactor=log_axial_prefactor,
            log_radial_core_cutoff_radius_prefactor=log_radial_core_radius_prefactor,
            log_axial_width_over_radial_core_cutoff_radius_prefactor=log_aspect_prefactor,
            axial_width_definition='Full axial F0 fractional-level width under original physical z map',
            radial_reference_definition='Core domain cutoff rho=4.1 at shared anchor, not measured vortex radial width',
            physical_width_time_exponent=dict(base_half=c.mpf('.5'),offset=-self.delta/2),
            local_F0_peak_levels_certified=True,whole_vortex_aspect_ratio_measured=False)

    def time_geometry(self,width,log_tau):
        c=self.ctx;t=c.mpf(log_tau)
        if not all(mp.isfinite(v) for v in endpoints(t)):raise ValueError('Finite logarithmic time required')
        common=t/2;elongation=-self.delta*t/2
        return dict(log_tau=t,physical_axial_width_log_terms=[common,elongation,width['log_axial_width_prefactor']],
            physical_radial_core_cutoff_radius_log_terms=[common,width['log_radial_core_cutoff_radius_prefactor']],
            axial_F0_width_over_core_cutoff_radius_log_terms=[width['log_axial_width_over_radial_core_cutoff_radius_prefactor'],elongation],
            relative_elongation_log_factor=elongation,
            split_logs_keep_tiny_delta_time_effect=True,coordinate_scaling_observation_not_fitted_dynamics=True,
            temporal_recursion=False)

    def swirl_value(self,xi,rho,log_tau='-1'):
        """Same nonlinear Phi, transported from the shared root with Xh bounds."""
        c=self.ctx;r=c.mpf(rho);peak=self.evaluate(xi);s=peak['xi'];logtau=c.mpf(log_tau)
        if endpoints(r)[0]<=0 or endpoints(r)[1]>endpoints(c.mpf('4.1'))[1]:
            raise ValueError('Positive core rho<=4.1 required; axis swirl is structurally zero')
        if not all(mp.isfinite(v) for v in endpoints(logtau)):raise ValueError('Finite logarithmic time required')
        if self.anchor_rows is None:self.anchor_rows=self.amplitude.rebuild.rebuild(self.anchor)
        anchor_profile=self.amplitude.rebuild.profile(self.anchor_rows,r)
        anchor_phi=anchor_profile['source_profile_enclosures']['Phi']
        if endpoints(anchor_phi)[0]<=0:raise ArithmeticError('Shared-root actual Phi positivity unresolved')
        # Uniform analytic derivative of the SAME admitted fixed point,
        # not an extrapolation of finite radial rows or a model replacement.
        derivative_bound=self.amplitude.rebuild.phi_norm*embedding(c,self.core.h,r,0,1)
        phi_change=derivative_bound*abs(self.b*s)
        phi=anchor_phi+symmetric(c,phi_change)
        phi_ratio_change=phi_change/c.mpf(endpoints(anchor_phi)[0])
        if endpoints(phi_ratio_change)[1]>=1:raise ArithmeticError('Nontrivial Phi transport unresolved')
        da=1-self.anchor**2;shift=-2*self.anchor*self.b*s-self.b**2*s**2
        dcover=1-(self.anchor+symmetric(c,self.offset_upper))**2
        if endpoints(dcover/da)[0]<mp.mpf('.5') or endpoints(self.delta)[0]<0 or endpoints(self.delta)[1]>1:
            raise ArithmeticError('Physical lambda ratio comparison outside admitted range')
        # For q=(1+delta)/2 in[.5,1] and t>=.5,
        # |t^q-1|<=2|t-1|. Keep the tiny increment factored.
        geometry_ratio_change=2*abs(shift)/da
        phi_ratio=1+symmetric(c,phi_ratio_change);geometry_ratio=1+symmetric(c,geometry_ratio_change)
        normalized=peak['normalized_F0_relative_to_shared_anchor']*phi_ratio*geometry_ratio
        loglambda=(logtau-c.ln(dcover))/2
        return dict(xi=s,rho=r,anchor_Phi=anchor_phi,transported_Phi_enclosure=phi,
            actual_Phi_axial_derivative_bound=derivative_bound,Phi_transport_abs_error_upper=phi_change,
            Phi_relative_to_shared_anchor_error_upper=phi_ratio_change,
            physical_lambda_ratio_error_upper=geometry_ratio_change,
            physical_utheta_relative_to_shared_anchor=normalized,
            physical_utheta_coefficient=phi*peak['normalized_F0_relative_to_shared_anchor'],
            physical_utheta_positive_scale_log_terms=[-self.core.logC,(c.ln(2*r)-self.core.logLambda)/2,
                                                      -loglambda,-self.delta*loglambda],
            same_nonlinear_Phi_from_fresh_radial_rows=True,shared_root_transport_uses_full_analytic_norm=True,
            absolute_utheta_not_materialized=True,normalization_compares_same_rho_and_tau=True,
            full_swirl_peak_location_measured=False,temporal_recursion=False)

    def report(self):
        samples={text:self.evaluate(text) for text in ('-4','-2','-1','-.5','0','.5','1','2','4')}
        levels={name:self.level_width(drop) for name,drop in (
            ('half_maximum',self.ctx.ln(2)),('one_over_e',self.ctx.mpf(1)),
            ('one_millionth',self.ctx.ln(10**6)))}
        geometry=[self.time_geometry(levels['half_maximum'],text) for text in ('-1','-10','-10000')]
        # A finite slow-time point preserves the original positive delta;
        # it is not a substituted larger delta or a numerical fit.
        geometry.append(self.time_geometry(levels['half_maximum'],-2/self.delta))
        swirl={xi:[self.swirl_value(xi,rho) for rho in ('2','4.1')] for xi in ('-1','0','1')}
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,root_interval=self.anchor,
            original_selected_logCstar=self.core.logC,original_delta=self.delta,
            scaled_chart=dict(domain=[-4,4],positive_b=self.b,log_b=self.log_b,
                exact_shared_anchor_factorization='H(a+x)=x*(Hprime(a)+(-12a-j)*x-4*x^2)',
                epsilon=self.epsilon,Hprime_at_anchor=self.hp,q_linear=self.q1,L_at_anchor=self.La,
                numerator_polynomial_coefficients=self.numerator,
                whole_chart_positive_q_lower=self.hp-self.q_change,
                whole_chart_positive_L_lower=self.La-self.L_change,
                denominator_remainder_coefficient=self.denominator_remainder_coefficient,
                gaussian_uniform_error_abs_upper=self.gaussian_uniform_error),
            peak_packets=samples,level_widths=levels,physical_geometry_packets=geometry,
            root_centered_physical_swirl_packets=swirl,
            actual_scaled_G_peak_resolved=True,normalized_axis_amplitude_peak_resolved=True,
            original_width_not_numerically_capped=True,absolute_F0_not_materialized=True,
            tiny_root_offset_and_log_amplitude_correction_kept_separate=True,
            same_nonlinear_Phi_weighted_local_swirl_values_available=True,
            full_swirl_peak_including_Phi_measured=False,whole_vortex_aspect_ratio_measured=False,
            full_point_physical_field_evaluation=False,all_annular_source_values_resolved=False,
            measured_blowup_dynamics=False,physical_energy_integral_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(400):result=CompliantRootCenteredPeak().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Original root-centered scaled G/F0 peak and local level widths generated',flush=True)
    return result


if __name__=='__main__':run()
