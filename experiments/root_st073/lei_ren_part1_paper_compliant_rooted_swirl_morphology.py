"""Actual local physical swirl extrema on fixed-physical-radius slices.

Uses the same nonlinear Phi, original G and physical mapping. The second
centered chart Z=a+b^2*y resolves the peak shift and log-height gain without
adding microscopic offsets. Radial monotonicity determines the next needed
annular values. No cutoff radius is relabeled as measured vortex width.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_rooted_core_field import CompliantRootedCoreField
from lei_ren_part1_paper_compliant_core_physical_field import symmetric
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


class CompliantRootedSwirlMorphology:
    def __init__(self):
        self.field=CompliantRootedCoreField();self.peak=self.field.peak;self.core=self.field.core;self.ctx=self.field.ctx
        self.hashes=dict(self.field.hashes);name=PREFIX+'rooted_core_field_check.json'
        receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or receipt['implicit_source_sha256']!=self.core.source:
            raise ValueError('Accepted same-source rooted vector required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Rooted vector changed: '+path)
        self.hashes.update(receipt['input_hashes'])
        for path in (name,Path(__file__).name):self.hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
        self.cache={}

    def slice_bounds(self,rho_anchor):
        c=self.ctx;r=c.mpf(rho_anchor)
        if endpoints(r)[0]<=0 or endpoints(r)[1]>4:raise ValueError('Positive anchor rho<=4 with core-domain buffer required')
        key=r._mpi_
        if key in self.cache:return self.cache[key]
        a=self.peak.anchor;b=self.peak.b;xi=c.mpf([-4,4]);z=a+symmetric(c,4*b)
        Da=1-a*a;D=1-z*z
        if endpoints(D/Da)[0]<mp.mpf('.5') or endpoints(self.core.delta)[0]<0 or endpoints(self.core.delta)[1]>1:
            raise ArithmeticError('Physical log-weight comparison outside admitted positive domain')
        # Fixed physical r at each tau: rho(Z)=rho_anchor*D(Z)/D(a).
        rho=r*D/Da
        if endpoints(rho)[0]<=0 or endpoints(rho)[1]>endpoints(c.mpf('4.1'))[1]:
            raise ArithmeticError('Physical-radius slice exits admitted core radius')
        jets={}
        for i,k in ((0,0),(1,0),(0,1),(0,2),(1,1),(2,0)):
            jets[(i,k)]=self.field.rooted_profile(xi,rho,i,k)['source_jets']['Phi']
        P,Pr,Pz,Pzz,Prz,Prr=(jets[key] for key in ((0,0),(1,0),(0,1),(0,2),(1,1),(2,0)))
        if endpoints(P)[0]<=0:raise ArithmeticError('Nonzero same-source Phi required')
        E=Pz-2*z*rho*Pr/D
        B=E/P-(2+self.core.delta)*z/D
        along_second=Pzz-4*z*rho*Prz/D+4*z*z*rho*rho*Prr/D**2-2*rho*Pr/D
        Bz=along_second/P-(E/P)**2-(2+self.core.delta)*(1+z*z)/D**2
        curvature=self.peak.derivative_jets(xi,1)[1]
        log_curvature=-curvature+b*b*Bz
        if endpoints(B)[0]<=0 or endpoints(curvature)[0]<=0 or endpoints(log_curvature)[1]>=0:
            raise ArithmeticError('Actual slice peak positivity/strict concavity unresolved')
        shift_D=2*abs(a)*4*b+(4*b)**2
        shift_rho=r*shift_D/Da
        anchor_phi=self.field.rooted_profile(0,r)['source_jets']['Phi']
        anchor_pr=self.field.rooted_profile(0,r,1,0)['source_jets']['Phi']
        anchor_pz=self.field.rooted_profile(0,r,0,1)['source_jets']['Phi']
        anchor_B=(anchor_pz-2*a*r*anchor_pr/Da)/anchor_phi-(2+self.core.delta)*a/Da
        if endpoints(anchor_B)[0]<=0:raise ArithmeticError('Actual positive anchor log-slope unresolved')
        change=abs(Pz)*4*b+abs(Pr)*shift_rho
        relative=change/c.mpf(endpoints(anchor_phi)[0])
        if endpoints(relative)[1]>=mp.mpf('.5'):raise ArithmeticError('Shared Phi log-ratio comparison unresolved')
        # |log(Phi/Phi(a))|<=2*relative; geometry exponent in[1,1.5].
        log_weight_error=2*relative+3*shift_D/Da
        result=dict(rho_anchor=r,fixed_physical_radius_definition='r=lambda(a)*sqrt(2*rho_anchor/Lambda) at requested tau',
            whole_slice_Z_cover=z,whole_slice_rho_cover=rho,Phi_cover=P,
            same_Phi_mixed_jets={str(key):value for key,value in jets.items()},log_weight_B=B,along_fixed_radius_B_Z=Bz,
            K_second_derivative=curvature,actual_log_swirl_xi_curvature=log_curvature,
            anchor_log_weight_B=anchor_B,actual_negative_log_swirl_xi_curvature=-log_curvature,
            Phi_anchor=anchor_phi,whole_slice_log_weight_error=log_weight_error,
            strictly_concave_log_swirl_on_rooted_slice=True)
        self.cache[key]=result;return result

    def peak_location(self,rho_anchor):
        c=self.ctx;bounds=self.slice_bounds(rho_anchor);B=bounds['log_weight_B'];h=bounds['K_second_derivative']
        lower=c.mpf(endpoints(B)[0])/c.mpf(endpoints(h)[1]);upper=c.mpf(endpoints(B)[1])/c.mpf(endpoints(h)[0])
        tolerance=c.mpf('1e-100');low=c.mpf(endpoints(lower-tolerance)[0]);high=c.mpf(endpoints(upper+tolerance)[1])
        if endpoints(low)[0]<=0:raise ArithmeticError('Positive true peak shift not resolved')
        if endpoints(self.peak.b*high)[1]>=4:raise ArithmeticError('Peak bracket exits local root chart')
        def scaled_slope(y):
            # xi=b*y, hence K_xi/b=y*L(a+b^2 y)*q(b^2 y)/denominator.
            x=self.peak.b*self.peak.b*y;z=self.peak.anchor+x
            L=1-self.core.delta*z*z;q=self.peak.hp+self.peak.q1*x-4*x*x
            coefficient=L*q/(1+self.peak.epsilon*self.peak.b**2*y*y*q*q)
            return B-y*coefficient
        inner=scaled_slope(low);outer=scaled_slope(high)
        if endpoints(inner)[0]<=0 or endpoints(outer)[1]>=0:raise ArithmeticError('Actual peak shift bracket signs not strict')
        y=c.mpf([endpoints(low)[0],endpoints(high)[1]])
        # H'(0)=b*B0 and -kappa_hi<=H''<=-kappa_lo on the slice.
        # Integrating these actual log-curvature bounds gives two parabolas;
        # the lower trial xi=b*B0_lo/kappa_hi must remain in this chart.
        B0=bounds['anchor_log_weight_B'];kappa=bounds['actual_negative_log_swirl_xi_curvature']
        trial=self.peak.b*c.mpf(endpoints(B0)[0])/c.mpf(endpoints(kappa)[1])
        if endpoints(trial)[0]<=0 or endpoints(trial)[1]>=4:
            raise ArithmeticError('Lower log-height comparison trial exits root chart')
        gain=c.mpf([endpoints(c.mpf(endpoints(B0)[0])**2/(2*c.mpf(endpoints(kappa)[1])))[0],
                    endpoints(c.mpf(endpoints(B0)[1])**2/(2*c.mpf(endpoints(kappa)[0])))[1]])
        return dict(rho_anchor=bounds['rho_anchor'],peak_y_interval=y,
            exact_peak_coordinate='Z_peak=a+b^2*y_peak; xi_peak=b*y_peak, b=sigma/sqrt(Lambda)',
            xi_peak_as_separate_factors=[self.peak.b,y],Z_peak_offset_as_separate_factors=[self.peak.b,self.peak.b,y],
            inner_scaled_log_slope=inner,outer_scaled_log_slope=outer,
            bracket_enlargement_each_endpoint=tolerance,
            normalized_log_peak_gain_over_b_squared=gain,
            actual_log_peak_gain_as_separate_factors=[self.peak.b,self.peak.b,gain],
            height_gain_definition='log[utheta(peak,r,tau)/utheta(anchor,r,tau)]=b^2*positive_gain',
            height_gain_uses_actual_log_curvature_and_anchor_slope=True,
            lower_height_comparison_trial_xi=trial,
            unique_actual_local_physical_swirl_peak=True,peak_is_strictly_to_right_of_F0_anchor=True,
            Phi_and_physical_lambda_derivatives_included=True,rooted_slice_bounds=bounds,
            global_swirl_peak_across_annuli_measured=False)

    def level_width(self,rho_anchor,log_drop):
        c=self.ctx;bounds=self.slice_bounds(rho_anchor);location=self.peak_location(rho_anchor)
        drop=c.mpf(log_drop);p=self.peak.numerator[0]
        gain=self.peak.b**2*location['normalized_log_peak_gain_over_b_squared']
        # Compare actual Phi/geometry/true-peak-normalized log profile with
        # the original source quadratic, retaining every directed correction.
        error=self.peak.gaussian_uniform_error+bounds['whole_slice_log_weight_error']+gain
        if endpoints(drop-error)[0]<=0:raise ValueError('Level is not separated from uniform comparison error')
        tol=c.mpf('1e-170')
        low=c.mpf(endpoints(c.sqrt(2*(drop-error)/c.mpf(endpoints(p)[1]))-tol)[0])
        high=c.mpf(endpoints(c.sqrt(2*(drop+error)/c.mpf(endpoints(p)[0]))+tol)[1])
        if endpoints(low)[0]<=0 or endpoints(high)[1]>4:raise ValueError('Actual level outside rooted chart')
        sides={}
        for sign,label in ((-1,'left'),(1,'right')):
            inner=-self.peak.evaluate(sign*low)['actual_scaled_G']+symmetric(c,bounds['whole_slice_log_weight_error'])-gain
            outer=-self.peak.evaluate(sign*high)['actual_scaled_G']+symmetric(c,bounds['whole_slice_log_weight_error'])-gain
            if endpoints(inner)[0]<=-endpoints(drop)[0] or endpoints(outer)[1]>=-endpoints(drop)[1]:
                raise ArithmeticError('Actual true-peak-normalized fractional level signs failed')
            sides[label]=dict(sign=sign,absolute_xi_level_interval=c.mpf([endpoints(low)[0],endpoints(high)[1]]),
                inner_log_ratio_to_true_peak=inner,outer_log_ratio_to_true_peak=outer,
                unique_by_strict_log_concavity=True)
        width=sides['left']['absolute_xi_level_interval']+sides['right']['absolute_xi_level_interval']
        z=bounds['whole_slice_Z_cover'];D=1-z*z;J=(1-self.core.delta*z*z)/c.power(D,(3-self.core.delta)/2)
        return dict(rho_anchor=bounds['rho_anchor'],log_drop=drop,normalized_to_actual_physical_swirl_peak=True,
            sides=sides,full_xi_width=width,physical_axial_width_prefactor_log=self.peak.log_b+c.ln(width)+c.ln(J),
            physical_axial_time_log_terms='logtau/2 and -delta*logtau/2 kept separately',
            radius_is_fixed_during_each_axial_scan=True,whole_vortex_aspect_ratio_measured=False)

    def radial_monotonicity(self):
        c=self.ctx;xi=c.mpf([-4,4]);r=c.mpf([0,'4.1'])
        P=self.field.rooted_profile(xi,r)['source_jets']['Phi']
        Pr=self.field.rooted_profile(xi,r,1,0)['source_jets']['Phi']
        coefficient=P+2*r*Pr
        if endpoints(coefficient)[0]<=0:raise ArithmeticError('Actual radial angular velocity derivative sign unresolved')
        return dict(xi_domain=[-4,4],rho_domain=[0,'4.1'],positive_radial_swirl_derivative_coefficient=coefficient,
            physical_identity='partial_r utheta=F0*lambda^(-2-delta)*(Phi+2rho Phi_rho)',
            actual_swirl_strictly_increases_radially_through_rooted_core=True,
            no_interior_radial_swirl_maximum_in_rooted_core=True,
            core_cutoff_cannot_be_reported_as_measured_radial_peak=True,
            global_radial_swirl_peak_requires_original_annular_values=True,
            vorticity_core_width_not_inferred_from_swirl_monotonicity=True)

    def report(self):
        locations={r:self.peak_location(r) for r in ('1','2','4')}
        widths={r:{name:self.level_width(r,drop) for name,drop in (
            ('half_maximum',self.ctx.ln(2)),('one_over_e',self.ctx.mpf(1)),('one_millionth',self.ctx.ln(10**6)))} for r in locations}
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,actual_physical_swirl_peak_locations=locations,
            true_peak_normalized_axial_widths=widths,radial_monotonicity=self.radial_monotonicity(),
            actual_Phi_weighted_local_physical_swirl_peak_certified=True,
            actual_local_axial_fractional_widths_available=True,
            complete_vortex_radial_width_measured=False,whole_vortex_aspect_ratio_measured=False,
            all_annular_source_values_resolved=False,full_point_physical_field_evaluation=False,
            measured_blowup_dynamics=False,physical_energy_integral_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(400):result=CompliantRootedSwirlMorphology().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual Phi-weighted physical swirl extrema, axial levels and radial monotonicity generated',flush=True)
    return result


if __name__=='__main__':run()
