"""Selected-source exact Gamma heat field and controlled moment deficits.

H is the true positive Gamma integral, not an infinite Taylor series or H=1.
Finite Taylor inequalities enclose its deficit while the selected inverse
radius is kept as a formal nonzero log scale. This component does not claim
that incoming moments or the Section7 correction thresholds are satisfied.
"""
# Recomputed for the distinct compliant pressure/moment family.
# Formula origin: lei_ren_part1_paper_shared_exact_heat_component.py; old .01 receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_outer_angular_candidate import SharedOuterAngularCandidate
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


class SharedExactHeatComponent:
    def __init__(self,cells=256):
        self.angular=SharedOuterAngularCandidate(cells);self.ctx=c=self.angular.ctx
        self.params=self.angular.params;self.delta=self.angular.delta;self.a=self.delta/2
        self.epsilon=self.angular.epsilon;self.hashes=dict(self.angular.hashes)
        for n in ('compliant_outer_angular_candidate','compliant_outer_angular_candidate_check'):
            name=PREFIX+n+'.json';record=json.loads((HERE/name).read_bytes())
            if record['actual_five_defect_family_sha256']!=self.angular.initial.family:raise ValueError('Heat source family mismatch')
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Heat component input changed: '+source)
                self.hashes[source]=digest
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            if not 0<endpoints(self.delta)[0]<=endpoints(self.delta)[1]<=mp.mpf('.5'):
                raise ValueError('Exact heat requires0<delta<=1/2')
            norms=self.angular.initial.repair.records['compliant_physical_norm_family']
            self.logC=read_interval(c,norms['selected_logCstar'])
            self.logRref=c.ln(110)+10*(self.logC+self.params.logPstar)
            self.tail_finite=self.params.yd+1+self.params.Tw+100-30*self.params.log_mu+2+self.params.Ts+self.angular.waiting
            self.logradius_terms=dict(selected_reference=self.logRref,pulse_term=13/self.params.mu,finite_offset=self.tail_finite)
            if endpoints(-self.logRref-13/self.params.mu-self.tail_finite)[1]>=-1000:
                raise ArithmeticError('Formal selected inverse-radius cap failed')
            self.Scap=c.exp(-1000)
            self.logS_terms={n:-v for n,v in self.logradius_terms.items()}
            self.heat_definition='H_delta(xi)=Gamma(1+a)^-1 integral_0^infinity exp(-v)*v^a*(1+xi*v)^-a dv, a=delta/2'
            self.logtail_relative=self.angular.waiting_field('0',1)['relative_log_profile'][0]
            self.logone=self.angular.waiting_logone
            self.preheat_atoms=self.preheat_backward_atoms()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def preheat_backward_atoms(self,cells=768):
        """Small epsilon corrections are stored separately from1/rate."""
        c=self.ctx;atoms={n:dict(W=c.mpf(0),W_squared=c.mpf(0)) for n in ('energy','pressure')}
        for i in range(cells):
            left=c.mpf(3)*i/cells;right=c.mpf(3)*(i+1)/cells
            t=c.mpf([endpoints(left)[0],endpoints(right)[1]])
            sig=stable_sigma(c,t)[0];edge=(3-t)/2
            phi=c.exp(-1/edge**2) if endpoints(edge)[0]>0 else c.mpf([0,endpoints(c.exp(-1/edge**2))[1]])
            W=1-sig+sig*phi
            W=c.mpf([max(mp.mpf(0),endpoints(W)[0]),min(mp.mpf(1),endpoints(W)[1])])
            for name,rate in (('energy',self.delta),('pressure',1+self.delta)):
                weight=(right-left)*c.mpf([endpoints(c.exp(-rate*right))[0],endpoints(c.exp(-rate*left))[1]])
                atoms[name]['W']+=weight*W;atoms[name]['W_squared']+=weight*W**2
        return atoms

    def deficit(self,Z,offset):
        """Exact 1-H = (a/Rtail)*Dhat, including first axial derivative."""
        c=self.ctx;Z=c.mpf(Z);t=c.mpf(offset)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(t)[0]<0:
            raise ValueError('Z in[-1,1], log(R/Rtail)>=0 required')
        decay=c.exp(-t);dz=1-Z**2;xi_cap=2*self.Scap*decay
        # H''<=a*(1+a)^2*(2+a) from the positive Gamma integral.
        B=(1+self.a)**2*(2+self.a)
        F=c.mpf([endpoints(1+self.a-B*xi_cap/2)[0],endpoints(1+self.a)[1]])
        G=c.mpf([endpoints(1+self.a-B*xi_cap)[0],endpoints(1+self.a)[1]])
        value=2*dz*decay*F;derivative=-4*Z*decay*G
        if endpoints(F)[0]<=0:raise ArithmeticError('Exact heat deficit lower factor lost')
        jet=IntervalTaylor(c,[value,derivative])
        magnitude=self.a*self.Scap*max(endpoints(value)[1],mp.mpf(0))
        H=c.mpf([max(mp.mpf(0),endpoints(1-magnitude)[0]),1])
        return dict(exact_heat_definition=self.heat_definition,
            xi_definition='xi=2*(1-Z^2)*exp(-offset)/Rtail',xi_positive_cap=xi_cap,
            deficit_scale=dict(log_a=c.ln(self.a),log_inverse_radius_terms=self.logS_terms),
            deficit_scaled_Taylor=jet,meaning='1-H_delta(xi)=(a/Rtail)*deficit_scaled_Taylor',
            H_value_enclosure=H,H_not_replaced_by_one=True,finite_point_Gamma_quadrature_evaluated=False,
            finite_remainder_inequality='a(1+a)xi-a(1+a)^2(2+a)xi^2/2 <=1-H<=a(1+a)xi; (1+a)-(1+a)^2(2+a)xi <=-Hprime/a<=1+a',
            infinite_Taylor_series_used=False)

    def profile(self,Z,offset):
        c=self.ctx;t=c.mpf(offset);h=self.deficit(Z,t);D=h['deficit_scaled_Taylor']
        sig=stable_sigma(c,t)[0]
        if endpoints(t)[0]>=3:phi=c.mpf(0)
        elif endpoints(t)[1]<=3:
            edge=(3-t)/2
            if endpoints(edge)[0]<=0:phi=c.mpf([0,endpoints(c.exp(-1/edge**2))[1]])
            else:phi=c.exp(-1/edge**2)
        else:raise ValueError('Subdivide across the collar endpoint t=3')
        W=1-sig+sig*phi;pre=1-self.epsilon*W;C=1-self.epsilon*phi
        Dactual=D*(sig*C)
        return dict(offset=t,Z=c.mpf(Z),heat=h,
            exact_bracket_definition='K=(1-sigma(t))*(1-epsilon)+sigma(t)*H_delta(xi)*(1-epsilon*f((3-t)/2)) for0<=t<=3; K=H_delta(xi) for t>=3',
            preheat_bracket=pre,actual_bracket_deficit_scaled=Dactual,
            velocity_definition='Utheta=c_infinity*R^(-(1+delta)/2)*[preheat_bracket-(a/Rtail)*actual_bracket_deficit_scaled]',
            heat_amplitude_definition='c_infinity=Utheta(Rtail,0)*Rtail^((1+delta)/2)/(1-epsilon)',
            heat_amplitude_log_terms=dict(formal_Utheta_Rv0_origin=self.angular.packet(self.angular.waiting_field('0',1))['formal_log_amplitude_origin'],
                logPstar=self.params.logPstar,log_tail_relative=self.logtail_relative,
                minus_log_one_minus_epsilon=-self.logone,radial_power=(1+self.delta)/2,
                logRtail_terms=self.logradius_terms),
            pressure_definition='Pheat(R,Z)=-integral_R^infinity Utheta_heat(rho,Z)^2/(2rho) drho; collar pressure uses the actual bracket and backward integral',
            Ur_exterior=0,Uz_exterior=0,exact_exterior_profile_specified=True,
            exterior_matches_incoming_moments=False,exterior_global_stress_zero_not_yet_certified=True,
            actual_heat_factor_or_deficit_not_zeroed=True,temporal_recursion=False)

    def future_defects(self,Z,cells=768):
        """Positive heat-vs-preheat defects in the Rtail normalizations.

        Angular heat difference = S*Theta_hat, pressure = a*S*Pressure_hat,
        and swirl energy = a*S*Energy_hat. S=1/Rtail is formal nonzero.
        These are enclosing exact positive integrals, not fitted moments.
        """
        c=self.ctx;Z=c.mpf(Z);dz=1-Z**2
        angular=c.mpf(0);pressure=c.mpf(0);energy=c.mpf(0)
        angular_Z=c.mpf(0);pressure_Z=c.mpf(0);energy_Z=c.mpf(0)
        absZ=max(abs(v) for v in endpoints(Z));B=(1+self.a)**2*(2+self.a)
        for i in range(cells):
            t=c.mpf([mp.mpf(3)*i/cells,mp.mpf(3)*(i+1)/cells]);dt=c.mpf(3)/cells
            sig=stable_sigma(c,t)[0];edge=(3-t)/2
            if endpoints(edge)[0]<=0:phi=c.mpf([0,endpoints(c.exp(-1/edge**2))[1]])
            else:phi=c.exp(-1/edge**2)
            C=1-self.epsilon*phi;W=1-sig+sig*phi;K0=1-self.epsilon*W
            h=self.deficit(Z,t);D=h['deficit_scaled_Taylor'];xi_cap=h['xi_positive_cap']
            # Kpre^2-Kactual^2 = a*S*Dhat*sig*C*(2Kpre-a*S*Dhat*sig*C).
            dcap=self.a*self.Scap*max(mp.mpf(0),endpoints(D[0])[1])
            square_factor=c.mpf([endpoints(2*K0-dcap*sig*C)[0],endpoints(2*K0)[1]])
            angular+=dt*c.exp((1-self.a)*t)*self.a*sig*C*D[0]
            angular_Z+=dt*c.exp((1-self.a)*t)*self.a*sig*C*D[1]
            pressure+=dt*c.exp(-(1+self.delta)*t)*sig*C*D[0]*square_factor/2
            energy+=dt*c.exp(-self.delta*t)*sig*C*D[0]*square_factor
            # Exact derivative of the square difference. K0,C have noZ
            # dependence and the same heat source supplies D_Z.
            derivative_factor=c.mpf([endpoints(2*K0-2*dcap*sig*C)[0],endpoints(2*K0)[1]])
            pressure_Z+=dt*c.exp(-(1+self.delta)*t)*sig*C*D[1]*derivative_factor/2
            energy_Z+=dt*c.exp(-self.delta*t)*sig*C*D[1]*derivative_factor
        # Entire exterior: use a finite remainder, integrate all power
        # weights exactly, and keep the positive leading terms.
        lead=2*dz*(1+self.a)*c.exp(-3*self.a)
        angular_error=2*self.a*(1+self.a)*(2+self.a)*dz**2*self.Scap*c.exp(-3*(1+self.a))
        angular+=c.mpf([max(mp.mpf(0),endpoints(lead-angular_error)[0]),endpoints(lead)[1]])
        angular_Z+=-4*Z*(1+self.a)*c.exp(-3*self.a)+c.mpf([-1,1])*8*self.a*(1+self.a)*(2+self.a)*absZ*self.Scap
        def square_tail(rate):
            leading=4*dz*(1+self.a)*c.exp(-3*(rate+1))/(rate+1)
            # 2D-D^2 differs from2a(1+a)xi by at most
            # [a*B+a^2(1+a)^2]*xi^2; divide by a*S.
            error=4*(B+self.a*(1+self.a)**2)*dz**2*self.Scap*c.exp(-3*(rate+2))/(rate+2)
            value=c.mpf([max(mp.mpf(0),endpoints(leading-error)[0]),endpoints(leading)[1]])
            derivative=-8*Z*(1+self.a)*c.exp(-3*(rate+1))/(rate+1)
            errZ=16*(B+self.a*(1+self.a)**2)*absZ*self.Scap*c.exp(-3*(rate+2))/(rate+2)
            return value,derivative+c.mpf([-endpoints(errZ)[1],endpoints(errZ)[1]])
        pv,pz=square_tail(1+self.delta);ev,ez=square_tail(self.delta)
        pressure+=pv/2;pressure_Z+=pz/2;energy+=ev;energy_Z+=ez
        return dict(Z=Z,angular_heat_difference_scaled=IntervalTaylor(c,[angular,angular_Z]),
            pressure_heat_difference_scaled=IntervalTaylor(c,[pressure,pressure_Z]),
            swirl_energy_heat_difference_scaled=IntervalTaylor(c,[energy,energy_Z]),
            definitions=dict(angular='int_Rtail^infinity sqrt(2R)*(Upre-Uactual)dR =sqrt(2)*c_infinity*Rtail^(1-a)*(1/Rtail)*angular_scaled',
                pressure='int_Rtail^infinity (Upre^2-Uactual^2)/(2R)dR =c_infinity^2*Rtail^(-1-delta)*(a/Rtail)*pressure_scaled',
                energy='int_Rtail^infinity (Upre^2-Uactual^2)dR =c_infinity^2*Rtail^(-delta)*(a/Rtail)*energy_scaled'),
            exact_backward_target_formulas=dict(
                angular='Mtheta_target/(sqrt(2)*c_infinity*Rtail^(1-a))=1/(1-a)+epsilon*Jcollar+(1/Rtail)*angular_scaled',
                pressure_integral='int_Rtail^infinity Uactual^2/(2R)dR /(c_infinity^2*Rtail^(-1-delta))=1/(2*(1+delta))-epsilon*IW_pressure+epsilon^2*IW2_pressure/2-(a/Rtail)*pressure_scaled',
                energy_integral='int_Rtail^infinity Uactual^2 dR /(c_infinity^2*Rtail^(-delta))=1/delta-2*epsilon*IW_energy+epsilon^2*IW2_energy-(a/Rtail)*energy_scaled'),
            backward_target_baselines=dict(angular=1/(1-self.a),pressure=1/(2*(1+self.delta)),energy=1/self.delta),
            preheat_epsilon_correction_atoms=self.preheat_atoms,
            angular_preheat_collar_J=self.angular.collarJ,
            finite_swirl_energy_tail_normalization_because_delta_positive=True,
            log_inverse_radius_terms=self.logS_terms,log_a=c.ln(self.a),
            exact_Gamma_heat_source_retained=True,infinite_exterior_integrals_included=True,
            incoming_angular_pressure_matching_completed=False)

    def report(self):
        c=self.ctx
        with mp.workdps(210):
            defects=[self.future_defects(Z) for Z in ('-1','0','.5','1')]
            current_ceps=self.epsilon/self.delta
            from fractions import Fraction
            exact_ceps=Fraction(self.angular.initial.datum.definition['c_epsilon'])
            paper_gate=0<exact_ceps<=Fraction(1,1000)
            return dict(actual_five_defect_family_sha256=self.angular.initial.family,
                implicit_source_sha256=self.angular.initial.datum.source_sha,datum_enclosure_sha256=self.angular.initial.datum.datum_sha,
                exact_heat_definition=self.heat_definition,delta=self.delta,formal_selected_Rtail_log_terms=self.logradius_terms,
                positive_inverse_radius_cap=self.Scap,profiles=[self.profile('.5',t) for t in (0,1,3,4)],
                future_heat_defects=defects,
                exact_selected_heat_component_and_collar_specified=True,
                directed_positive_heat_deficit_and_three_future_defects_callable=True,
                finite_remainder_used_without_convergent_series=True,
                exterior_swirl_radial_energy_integrable_because_delta_positive=True,
                full_physical_kinetic_energy_certified=False,
                paper_Section7_c_epsilon_requirement='c_epsilon<=.001',current_c_epsilon=current_ceps,
                paper_Section7_c_epsilon_gate=paper_gate,
                c_epsilon_gate_proof='exact source rational epsilon/delta=1/1000, not an independently divided interval box',
                all_Section7_hypotheses_certified=False,
                existing_source_not_silently_changed=True,
                angular_pressure_corrections_installed=False,actual_ap_selected=False,
                heat_exterior_matched_to_incoming_five_moments=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    field=SharedExactHeatComponent();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual Gamma heat component/collar and three positive future-defect enclosures generated; global matching pending',flush=True)
    print('Section7 c_epsilon<=.001 gate:',result['paper_Section7_c_epsilon_gate'],'source unchanged',flush=True)
    return result


if __name__=='__main__':run()
