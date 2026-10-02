"""Original 100..110 switches and actual reshape inlet, axial enclosures.

The microscopic radii remain 100*exp(hb*s), with the admitted positive hb.
Enclosure caps are never substituted into the source equations. This module
continues the ACTUAL bridge histories, not the comparison's moments.
Radial mixed derivatives and whole-field stress/temporal claims remain open.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    CompliantInnerBridgeProfiles, IntervalTaylor, cumulative_enclosures,
    logarithm, derivative, square, dress, coefficient_lists, symmetric,
    MTH, MZ, MTHZ, MZT, MP)
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def as_initial(moments):
    return dict(H=moments[MTH],mean=moments[MZ],K=moments[MTHZ],
                A=moments[MZT]['axial'],B=moments[MZT]['swirl'],C=moments[MP])


def power_transport(c,theta,initial,phi2,v2):
    """Exact positive-integral weights for a=4/5,b=0 after source R2.

    Shapes are normalized by original F0, not F2. Initial histories remain.
    theta=R2/R; intersect differences with their known nonnegative ranges
    so directed rounding near theta=1 cannot invent negative integral mass.
    """
    if endpoints(theta)[0]<=0 or endpoints(theta)[1]>1:
        raise ValueError('Post-switch transport requires 0<theta<=1')
    def difference(first,second):
        value=theta**first-theta**second
        return c.mpf([max(mp.mpf(0),endpoints(value)[0]),endpoints(value)[1]])
    angular=difference(c.mpf(2)/5,2)*c.mpf(5)/4
    swirl=difference(c.mpf(4)/5,2)*c.mpf(5)/6
    pressure=difference(c.mpf(4)/5,1)*5
    return {MTH:initial['H']*theta**2+phi2*angular,
        MZ:initial['mean']*theta+v2*(1-theta),
        MTHZ:initial['K']*theta**2+(phi2*v2)*angular,
        MZT:dict(axial=initial['A']*theta+square(v2)*(1-theta),
                 swirl=initial['B']*theta**2+square(phi2)*swirl),
        MP:initial['C']*theta+square(phi2)*pressure}


class CompliantInnerSwitchProfiles:
    def __init__(self):
        self.bridge=CompliantInnerBridgeProfiles();self.core=self.bridge.core
        self.ctx=c=self.bridge.ctx;self.family=self.bridge.family;self.source=self.bridge.source
        self.logh=self.bridge.logh;self.cap=self.bridge.cap;self.hashes=dict(self.bridge.hashes)
        name=PREFIX+'inner_bridge_profiles_check.json'
        self.check=json.loads((HERE/name).read_bytes())
        for path,digest in self.check['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Switch prerequisite changed: '+path)
        if not (self.check['all_passed'] and self.check['actual_five_defect_family_sha256']==self.family
                and self.check['implicit_source_sha256']==self.source
                and self.check['actual_prescribed_shear_bridge_axial5_enclosures_available']
                and self.check['angular_log_sign_bound_to_accepted_same_source_comparison_theorem']):
            raise ValueError('Actual same-source bridge receipt required')
        ledger=self.bridge.records['K1_ledger']
        self.a_upper=read_interval(c,ledger['decreasing_K_smallness_bounds']['initial_switch_a'])
        if endpoints(self.a_upper)[1]>=endpoints(c.mpf(4)/5)[0]:
            raise ValueError('Original switch angular shear bound not admitted')
        if endpoints(self.logh)[1]>=endpoints(c.ln(c.ln(c.mpf(110)/100)/2))[0]:
            raise ValueError('Formal source R2 is not below 110')
        self.hashes.update(self.check['input_hashes'])
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        self.cache={}

    def inputs(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=Z._mpi_
        if key in self.cache:return self.cache[key]
        incoming=self.bridge.actual(Z,self.bridge.r/100)
        jet=lambda values:IntervalTaylor(c,values)
        phi=jet(incoming['F_actual_over_F0_axial5_coefficients'])
        v=jet(incoming['Uz_actual_axial5_coefficients'])
        stored=incoming['actual_moment_shape_axial5_coefficients']
        moments={k:({n:jet(q) for n,q in row.items()} if k==MZT else jet(row)) for k,row in stored.items()}
        source=self.bridge.inputs(Z)
        comparison=self.bridge.comparison(Z,c.mpf([endpoints(self.bridge.r/110)[0],endpoints(self.bridge.r/100)[1]]))
        # JD has total positive weight 1 + 1/2. Dbar=R*(Dbar/R), R<=110.
        h2=[self.bridge.width_product(comparison['direction']['D_over_R'][k]*c.mpf('82.5'),self.logh)
            for k in range(6)]
        det=self.bridge.width_product(c.mpf('.4'))
        # This bounds every PREFIX of both switches. Sigma contributes only
        # a Z-independent term; positive Dbar makes the zeroth log <=0.
        ell=IntervalTaylor(c,[c.mpf([-endpoints(h2[0]+det)[1],0])]
                          +[symmetric(c,value) for value in h2[1:]])
        pc=phi*ell.exp()
        quotient=pc/comparison['phi'].truncate(5)
        drive={name:quotient*value for name,value in comparison['direction'].items() if name.startswith('drive_')}
        inc=[]
        for k in range(6):
            # First switch integral weight <=1; b=0 from phase1 onward.
            terms=[self.bridge.width_product(drive['drive_hydro'][k]*110,self.logh),
                   self.bridge.width_product(drive['drive_pressure'][k]*110,self.logh+2*self.core.logP),
                   self.bridge.width_product(drive['drive_swirl'][k]*12100,self.logh-2*self.core.logC)]
            inc.append(symmetric(c,sum(terms,c.mpf(0))))
        vc=v+IntervalTaylor(c,inc)
        result=dict(incoming=incoming,phi=phi,v=v,moments=moments,source=source,
                    comparison=comparison,h2_bounds=h2,short_log_cover=ell,
                    phi_cover=pc,v_cover=vc,velocity_increment=inc)
        self.cache[key]=result;return result

    def packet(self,Z,phi,v,moments,**metadata):
        c=self.ctx;inputs=self.inputs(Z);source=inputs['source']
        z=IntervalTaylor.variable(c,c.mpf(Z),5);d=1-square(z);L=1-square(z)*self.core.delta
        m=moments[MZ];Q=(2*z*v-(z*m)*(1-self.core.delta)-d*derivative(m))/L
        return dict(Z=c.mpf(Z),F_actual_over_F0_axial5_coefficients=list(phi.coefficients),
            F_actual_true_axial5_divided_by_F0=list(dress(phi,source['F0_ratios']).coefficients),
            Uz_actual_axial5_coefficients=list(v.coefficients),
            actual_moment_shape_axial5_coefficients=coefficient_lists(moments),
            actual_Q_axial4_coefficients=list(Q.coefficients),
            pressure_axis_axial5_coefficients=list(source['p0'].truncate(5).coefficients),
            pressure_increment_true_axial5_divided_by_R_F0_squared=list(dress(moments[MP],source['F0_squared_ratios']).coefficients),
            actual_velocity_prefactors=dict(Ur='sqrt(R/2)*Q',Utheta='sqrt(2R)*F0*phi',Uz='V'),
            original_P0_retained=True,all_actual_moments_inherited=True,
            formal_positive_hb_used_not_numerical_cap=True,comparison_moments_substituted=False,
            switch_radial_mixed4_certified=False,**metadata)

    def phase(self,Z,phase):
        """Actual source R=100*exp(hb*phase), 0<=phase<=2, axial5 cover.

        An interval phase encloses its whole microscopic source collar.
        Separate endpoints express exact source histories/identities.
        """
        c=self.ctx;phase=c.mpf(phase);slo,shi=endpoints(phase)
        if slo<0 or shi>2:raise ValueError('Original switch phase in[0,2] required')
        inp=self.inputs(Z);start=slo==shi==0;end=slo==shi==2
        if start:
            phi,v,moments=inp['phi'],inp['v'],inp['moments'];theta=c.mpf(1)
            ell=IntervalTaylor.constant(c,0,5)
        else:
            # hb and source R are not evaluated. These are enclosure-only
            # theta bounds with an explicit exact source expression.
            theta=c.mpf([endpoints(c.exp(-phase*self.cap))[0],1])
            ell=inp['short_log_cover']
            if end:
                hterm=self.bridge.width_product(c.mpf('.2'))
                ell=IntervalTaylor(c,[c.mpf([-endpoints(inp['h2_bounds'][0]+hterm)[1],0])]
                    +[symmetric(c,value) for value in inp['h2_bounds'][1:]])
            phi=inp['phi']*ell.exp();v=inp['v_cover']
            moments=cumulative_enclosures(c,theta,as_initial(inp['moments']),inp['phi_cover'],inp['v_cover'])
        return self.packet(Z,phi,v,moments,phase=phase,theta_100_over_R_enclosure=theta,
            exact_source_radius='100*exp(hb*phase)',
            log_F_actual_over_F100_coefficients=list(ell.coefficients),
            exact_R100_history=start,Uz_constant_from_phase1_onward=slo>=1,
            exact_Uz_source='V100 - hb^2*integral_0^min(phase,1)(1-sigma(t))*(phi_actual/barphi)*drive(100exp(hb*t),Z)dt',
            angular_source='-.5*hb^2*[integral_0^min(phase,1)Dbar(100exp(hb*t),Z)dt + integral_0^max(phase-1,0)(1-sigma(t))*Dbar(100exp(hb*(1+t)),Z)dt] -.4*hb*integral_0^max(phase-1,0)sigma(t)dt',
            exact_R2_log_identity='log(F2/F100)=-.2*hb-.5*hb^2*JD' if end else None,
            current_radius_Dbar_Ebar_retained=True)

    def post(self,Z,R):
        """Actual a=4/5,b=0 interval from formal R2 through fixed R<=110."""
        c=self.ctx;R=c.mpf(R)
        if endpoints(R)[0]<=100 or endpoints(R)[1]>110:raise ValueError('Fixed post-switch R in(100,110] required')
        y=c.ln(R/100)
        # Compare loghb directly to log(y/2), not to its numerical cap.
        if endpoints(self.logh)[1]>=endpoints(c.ln(y/2))[0]:
            raise ValueError('Fixed radius is not certified after formal R2; use phase()')
        inp=self.inputs(Z);p2=self.phase(Z,2)
        jet=lambda row:IntervalTaylor(c,row)
        phi2=jet(p2['F_actual_over_F0_axial5_coefficients'])
        v2=jet(p2['Uz_actual_axial5_coefficients'])
        stored=p2['actual_moment_shape_axial5_coefficients']
        moments2={k:({n:jet(q) for n,q in row.items()} if k==MZT else jet(row)) for k,row in stored.items()}
        # Enclosure of formal theta=(100/R)*exp(2hb), keeping tiny offsets.
        theta=intersection(c,(c.mpf([1,endpoints(c.exp(2*self.cap))[1]]))*(100/R),c.mpf([0,1]))
        moments=power_transport(c,theta,as_initial(moments2),phi2,v2)
        plus=self.bridge.width_product(c.mpf('.6'))
        correction=IntervalTaylor(c,[c.mpf([-endpoints(inp['h2_bounds'][0])[1],endpoints(plus)[1]])]
            +[symmetric(c,value) for value in inp['h2_bounds'][1:]])
        phi=(inp['phi']*correction.exp())*c.exp(-y*c.mpf(2)/5)
        return self.packet(Z,phi,v2,moments,R=R,theta_R2_over_R_enclosure=theta,
            exact_source_radius='R (fixed); inherited source R2=100*exp(2hb)',
            log_correction_to_F100_power_coefficients=list(correction.coefficients),
            exact_post_log_identity='log(F/F100)=-(2/5)*log(R/100)+.6*hb-.5*hb^2*JD',
            exact_Uz_source='V(R,Z)=V(100exp(hb),Z) for every R>=100exp(hb)',
            Uz_constant_from_phase1_onward=True,source_R2_not_rounded_to_100=True,
            exact_power_moment_transport_from_actual_R2=True)

    def inlet(self,Z):
        c=self.ctx;result=self.post(Z,110)
        phi=IntervalTaylor(c,result['F_actual_over_F0_axial5_coefficients'])
        z=IntervalTaylor.variable(c,c.mpf(Z),5);source=self.core.axis_inputs(c.mpf(Z))
        gradient=source['L']*source['H']/(square(source['H'])+self.core.sigma**2)
        # G is anchored at the admitted unique H0 root. Its value is kept
        # within the admitted positive bound; its derivatives are exact jets.
        G=IntervalTaylor(c,[c.mpf([0,endpoints(self.core.Gbar)[1]])]
                         +[gradient[k-1]/k for k in range(1,6)])
        B=G*(-self.core.Lambda)+logarithm(phi)+logarithm(1+square(z))+c.ln(c.mpf(220))/2
        result.update(actual_R110_log_shape_axial5_coefficients=list(B.coefficients),
            log_shape_source='B=log(Cstar*u110*(1+Z^2))=-Lambda*G+log(phi110)+.5log(220)+log(1+Z^2)',
            logCstar_cancelled_before_enclosure=True,actual_R110_inlet_axial5_available=True,
            long_reshape_with_this_inlet_installed=False)
        return result

    def report(self):
        whole=self.phase([-1,1],[0,2]);start=self.phase([-1,1],0)
        first=self.phase([-1,1],1);second=self.phase([-1,1],2);inlet=self.inlet([-1,1])
        samples=[self.phase(z,s) for z in ('0','.5') for s in ('.5','1.5')]
        samples += [self.post(z,105) for z in ('0','.5')]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            admitted_inner_parameter_family_sha256=self.bridge.records['K1_ledger']['admitted_inner_parameter_family_sha256'],
            original_prescription='(4.37)/(9.30): first a=hb*Dbar,b=-hb*Ebar*(1-sigma(s)); second a=hb*Dbar*(1-sigma(s-1))+.8*sigma(s-1),b=0; then a=.8,b=0 to110',
            source_log_hb_enclosure=self.logh,source_hb_positive_not_materialized=True,cap_is_not_width=True,
            sigma_symmetry='sigma(t)+sigma(1-t)=1, integral_0^1 sigma=1/2',
            JD_source='integral_0^1 Dbar(100exp(hb*t),Z)dt + integral_1^2(1-sigma(t-1))*Dbar(100exp(hb*t),Z)dt',
            JD_positive_weight='1+1/2=3/2',
            first_switch=first,second_switch=second,whole_short_switches=whole,
            actual_R100_inlet=start,actual_R110_inlet=inlet,samples=samples,
            angular_sign_and_shear_source=dict(certificate=PREFIX+'global_exit_certificate.json',
                ledger=PREFIX+'K1_ledger.json',same_source_Dbar_positive=True,
                initial_hb_Dbar_upper=self.a_upper,second_a_is_convex_interpolation=True),
            source_log_product_cap_proofs=self.bridge.cap_proofs,
            actual_100_110_switch_axial5_enclosures_available=True,
            actual_100_110_radial_recovery_axial4_available=True,
            actual_R110_inlet_axial5_available=True,
            switch_radial_mixed4_certified=False,bridge_radial_mixed4_certified=False,
            newly_recomputed_point_coefficients=False,long_reshape_with_this_inlet_installed=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantInnerSwitchProfiles().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Original 100..110 switch axial5 enclosures and ACTUAL 110 reshape inlet generated',flush=True)
    return result


if __name__=='__main__':run()
