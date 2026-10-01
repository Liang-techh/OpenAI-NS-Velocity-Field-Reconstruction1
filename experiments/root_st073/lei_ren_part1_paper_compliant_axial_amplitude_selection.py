"""Actual Section 7.34 positive axial pulse amplitude and end coefficients.

The full corrected future swirl energy supplies the fixed right side.
End coefficients are affine in a with a formal super-small exponential.
The positive quadratic root and its C1 derivative are directed enclosures
of the actual smooth functions, not midpoint or nominal a=1 choices.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_future_swirl_energy import CompliantFutureSwirlEnergy
from lei_ren_part1_paper_compliant_outer_pulse_map import SharedOuterPulseMap
from lei_ren_part1_paper_compliant_outer_angular_repair import magnitude
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


class CompliantAxialAmplitude:
    def __init__(self):
        self.future=CompliantFutureSwirlEnergy();self.pulse=SharedOuterPulseMap()
        self.ctx=c=self.future.ctx;self.mu=self.future.mu;self.K=self.pulse.Kpulse
        self.hashes=dict(self.future.hashes);self.hashes.update(self.pulse.hashes)
        for stem in ('compliant_future_swirl_energy','compliant_future_swirl_energy_check',
                     'compliant_outer_pulse_map','compliant_outer_pulse_map_check'):
            name=PREFIX+stem+'.json';record=json.loads((HERE/name).read_bytes())
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Amplitude prerequisite changed: '+source)
                self.hashes[source]=digest
            if record['actual_five_defect_family_sha256']!=self.future.angular.initial.family:
                raise ValueError('Amplitude family mismatch')
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            self.base=(1-c.exp(-26))/4
            self.log_end_scale=self.pulse.logscale
            self.log_energy_cap=2*self.future.params.log_mu-1000
            self.energy_cap=c.exp(self.log_energy_cap)
            self.energy_log_margins=[]
            for Kj in self.pulse.basis['end_energy_weights']:
                logexact=self.future.params.log_mu+c.ln(Kj)+2*self.log_end_scale
                margin=self.log_energy_cap-logexact
                if endpoints(margin)[0]<=0:raise ArithmeticError('Formal end-energy finite cap failed')
                self.energy_log_margins.append(margin)
            # nu_j=mu*Kj*exp(2log_end_scale) is exact and positive.
            # A [0,cap] enclosure is used without exponentiating log_end_scale.
            self.nu=[c.mpf([0,endpoints(self.energy_cap)[1]]) for _ in range(2)]
            Pi=self.pulse.rows['scaled_full_rows'];self.v=self.linear_inverse(-self.mu*Pi[0],-(Pi[1]-Pi[0]))
            # Bind row-normalized incoming functions to their actual jets.
            # Qi=mi*exp(-13lambda_i/mu-common_logpref), lambda_i=.5-i*mu.
            # One fixed finite factor cap is admitted over the whole Z domain;
            # multiplying mi and mi_Z by it encloses the same true source.
            trial_whole=self.pulse.coefficients([-1,1],1)
            self.incoming_factor_log_caps=[];self.incoming_factor_caps=[]
            self.incoming_factor_log_definitions=[];self.incoming_factor_margins=[]
            for row,jet in enumerate(trial_whole['actual_incoming_moment_Taylor_enclosures'],1):
                norm=sum((c.mpf(magnitude(v)) for v in jet.coefficients),c.mpf(0))
                ln_norm=endpoints(c.ln(norm))[1]
                cut=c.mpf(-1100)-c.mpf(max(mp.mpf(0),ln_norm))
                exact_log_factor=-13/(2*self.mu)+13*row-self.pulse.rows['common_logpref']
                margin=cut-exact_log_factor
                if endpoints(margin)[0]<=0 or endpoints(c.ln(norm)+cut)[1]>=-1000:
                    raise ArithmeticError('Actual incoming row factor/C1 cap not proved')
                self.incoming_factor_log_caps.append(cut)
                self.incoming_factor_caps.append(c.mpf([0,endpoints(c.exp(cut))[1]]))
                self.incoming_factor_log_definitions.append(exact_log_factor)
                self.incoming_factor_margins.append(margin)
            whole=self.select([-1,1]);a=whole['selected_ap_Taylor']
            if not (mp.mpf('.9')<endpoints(a[0])[0]<=endpoints(a[0])[1]<mp.mpf('1.2')):
                raise ArithmeticError('Actual selected amplitude outside paper interval')
            self.whole=whole
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def linear_inverse(self,r1,r2):
        A=self.pulse.basis['first_row'];D=self.pulse.basis['exact_divided_difference_row']
        det=self.pulse.basis['divided_difference_determinant']
        return [(r1*D[1]-r2*A[1])/det,(r2*A[0]-r1*D[0])/det]

    def affine(self,Z):
        """Actual c_j=exp(log_end_scale)*(u_j(Z)+v_j*a)."""
        c=self.ctx;trial=self.pulse.coefficients(Z,1)
        incoming=[IntervalTaylor(c,jet.coefficients)*factor for jet,factor in zip(
            trial['actual_incoming_moment_Taylor_enclosures'],self.incoming_factor_caps)]
        u=self.linear_inverse(incoming[0]*(-self.mu),-(incoming[1]-incoming[0]))
        return u,trial,incoming

    def select(self,Z):
        c=self.ctx;Z=c.mpf(Z);u,trial,row_incoming=self.affine(Z)
        energy=self.future.future(Z);future=energy['Section7_34_weighted_future_Taylor']
        incoming=IntervalTaylor(c,trial['actual_normalized_incoming_energy'])*self.mu
        target=future-incoming+self.base
        if target.order!=1:raise ArithmeticError('Amplitude target C1 lost')
        A2=self.K;A1=u[0]*0;A0=-target
        for nu,uj,vj in zip(self.nu,u,self.v):
            A2+=nu*vj*vj
            A1+=uj*(2*nu*vj)
            A0+=uj*uj*nu
        disc=A1[0]**2-4*A2*A0[0]
        if endpoints(A2)[0]<=0 or endpoints(A0[0])[1]>=0 or endpoints(disc)[0]<=0:
            raise ArithmeticError('Positive quadratic amplitude root not proved')
        # Completed-square form avoids repeating the uncertain K twice in
        # a sqrt(K)/K enclosure; no small subtractive cancellation occurs.
        halflinear=A1[0]/(2*A2)
        ap=-halflinear+c.sqrt(halflinear**2-A0[0]/A2)
        denominator=2*A2*ap+A1[0]
        if endpoints(denominator)[0]<=0:raise ArithmeticError('Amplitude derivative inverse failed')
        apZ=-(A1[1]*ap+A0[1])/denominator
        selected=IntervalTaylor(c,[ap,apZ])
        scaled=[uj+selected*vj for uj,vj in zip(u,self.v)]
        if not endpoints(scaled[0][0])[1]<0<endpoints(scaled[1][0])[0]:
            raise ArithmeticError('Actual end coefficient signs unresolved')
        low=c.mpf('.9');high=c.mpf('1.2')
        low_residual=A2*low**2+A1[0]*low+A0[0]
        high_residual=A2*high**2+A1[0]*high+A0[0]
        if endpoints(low_residual)[1]>=0 or endpoints(high_residual)[0]<=0:
            raise ArithmeticError('Whole-Z endpoint root bracket failed')
        end_energy=scaled[0]*0
        for nu,coeff in zip(self.nu,scaled):end_energy+=coeff*coeff*nu
        return dict(Z=Z,selected_ap_Taylor=selected,
            selected_scaled_end_coefficient_Taylor=scaled,
            end_coefficient_definition='cj=exp(log_end_scale)*(uj(Z)+vj*ap(Z)); scale is formal nonzero, not materialized',
            common_log_end_coefficient_scale=self.log_end_scale,
            actual_affine_incoming_Taylor_enclosures=u,affine_pulse_slopes=self.v,
            actual_incoming_moment_Taylor_enclosures=trial['actual_incoming_moment_Taylor_enclosures'],
            actual_row_normalized_incoming_Taylor_enclosures=row_incoming,
            actual_row_incoming_definition='Qi=mi(Z)*exp(-13lambda_i/mu-common_logpref); Qi_Z=mi_Z times the SAME fixed factor',
            actual_weighted_incoming_energy_Taylor=incoming,actual_weighted_future_energy_Taylor=future,
            energy_target_Taylor=target,quadratic_coefficients=dict(A2=A2,A1=A1,A0=A0),
            positive_discriminant=disc,positive_root_derivative_denominator=denominator,
            low_endpoint_energy_residual=low_residual,high_endpoint_energy_residual=high_residual,
            end_energy_positive_enclosure=end_energy,
            exact_end_energy_weights='nu_j=mu*Kj*exp(2log_end_scale)>0; finite cap is only an enclosure',
            exact_selected_energy_equation='Kpulse*ap^2+sum_j nu_j*(uj+vj*ap)^2=(1-exp(-26))/4-mu*incoming_energy+weighted_full_future_energy',
            actual_ap_selected=True,actual_end_linear_and_energy_moment_closure=True,
            unique_small_positive_smooth_branch=True,all_incoming_and_end_tails_retained=True,
            full_outer_five_moment_field_built=False,whole_outer_cone_certified=False,temporal_recursion=False)

    def report(self):
        with mp.workdps(210):
            return dict(actual_five_defect_family_sha256=self.future.angular.initial.family,
                implicit_source_sha256=self.future.angular.initial.datum.source_sha,
                datum_enclosure_sha256=self.future.angular.initial.datum.datum_sha,
                fixed_pulse_K=self.K,common_log_end_scale=self.log_end_scale,
                finite_log_end_energy_cap=self.log_energy_cap,positive_end_energy_cap=self.energy_cap,
                exact_end_energy_log_margin_checks=self.energy_log_margins,
                actual_incoming_row_factor_log_definitions=self.incoming_factor_log_definitions,
                finite_incoming_row_factor_log_caps=self.incoming_factor_log_caps,
                actual_incoming_row_factor_positive_caps=self.incoming_factor_caps,
                actual_incoming_row_factor_log_margins=self.incoming_factor_margins,
                actual_incoming_jets_and_first_derivatives_explicitly_scaled=True,
                whole_Z_C1_selected_amplitude=self.whole,samples=[self.select(z) for z in ('-1','0','.5','1')],
                actual_ap_selected=True,full_future_corrected_swirl_energy_used=True,
                selected_actual_c1_c2_functions_defined=True,actual_O4_partial_pulse_field_installed=False,
                full_outer_five_moment_match=False,whole_outer_cone_certified=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    amplitude=CompliantAxialAmplitude();result=amplitude.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    a=result['whole_Z_C1_selected_amplitude']['selected_ap_Taylor'][0]
    print('Actual whole-Z positive axial amplitude selected:',[mp.nstr(v,14) for v in endpoints(a)],'; formal end tails retained; full pulse/outer field pending',flush=True)
    return result


if __name__=='__main__':run()
