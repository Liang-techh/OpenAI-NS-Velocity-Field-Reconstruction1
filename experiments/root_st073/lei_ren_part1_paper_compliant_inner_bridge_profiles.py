"""Smoothed comparison and actual prescribed-shear bridge axial enclosures.

Source hb=cstar*K^-100 stays exact and formal. All bounds use its admitted
log enclosure before capping width-times-pressure/gradient products. These
are enclosures of the original integrations, not midpoint coefficients.
Radial mixed4 and the switches at100..110 are not certified here.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_frozen_comparison_field import (
    CompliantFrozenComparisonField,axial_jet,derivative,square,dress,MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def norm(c,value):return c.mpf(max(abs(v) for v in endpoints(value)))


def symmetric(c,bound):
    upper=endpoints(bound)[1]
    return c.mpf([-upper,upper])


def logarithm(jet):
    quotient=derivative(jet)/jet
    return IntervalTaylor(jet.ctx,[jet.ctx.ln(jet[0])]+[quotient[k-1]/k for k in range(1,jet.order+1)])


def cumulative_enclosures(c,theta,initial,phi_cover,v_cover):
    """Own primitives from exact positive radial weights, not frozen resets."""
    return {MTH:initial['H']*theta**2+phi_cover*(1-theta**2),
        MZ:initial['mean']*theta+v_cover*(1-theta),
        MTHZ:initial['K']*theta**2+(phi_cover*v_cover)*(1-theta**2),
        MZT:dict(axial=initial['A']*theta+square(v_cover)*(1-theta),
                 swirl=initial['B']*theta**2+square(phi_cover)*((1-theta**2)/2)),
        MP:initial['C']*theta+square(phi_cover)*(1-theta)}


def direction(c,Z,delta,phi,v,moments,pressure,ratios,ratios2):
    """(9.13), including true amplitude derivatives, with axial drive split."""
    z=IntervalTaylor.variable(c,Z,phi.order);d=1-square(z);L=1-square(z)*delta
    m=moments[MZ];a=moments[MZT]['axial']
    f=dress(phi,ratios);h=dress(moments[MTH],ratios);k=dress(moments[MTHZ],ratios)
    b=dress(moments[MZT]['swirl'],ratios2);p=dress(moments[MP],ratios2)
    W=1-(z*m)*(1-delta)-d*derivative(m)
    angular=h*(1-delta/2)-(z*derivative(h))*((1-delta)/2)-d*derivative(k)+(z*k)*(2*delta-1)
    return dict(D_over_R=(-W+angular/(2*f))/L,
        drive_hydro=(-W*v+(m-z*derivative(m))*((1-delta)/2)+(z*a)*(2*delta)-d*derivative(a))/(2*L),
        drive_pressure=((z*pressure)*(2*(1+delta))-d*derivative(pressure))/(2*L),
        drive_swirl=(-(z*b)*(2*delta)+d*derivative(b)+(z*p)*(2*(1+delta))-d*derivative(p))/(2*L))


def coefficient_lists(value):
    if isinstance(value,IntervalTaylor):return list(value.coefficients)
    return {key:coefficient_lists(row) for key,row in value.items()}


class CompliantInnerBridgeProfiles:
    def __init__(self):
        self.frozen=CompliantFrozenComparisonField();self.core=self.frozen.core
        self.ctx=c=self.core.ctx;self.family=self.core.family;self.source=self.core.source
        self.delta=self.core.delta;self.r=self.frozen.r;self.hashes=dict(self.frozen.hashes)
        self.records={}
        for part in ('frozen_comparison_field_check','K1_ledger','global_exit_certificate'):
            name=PREFIX+part+'.json';record=json.loads((HERE/name).read_bytes())
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Bridge source prerequisite changed: '+path)
            self.records[part]=record;self.hashes.update(record['input_hashes'])
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check=self.records['frozen_comparison_field_check'];ledger=self.records['K1_ledger'];gate=self.records['global_exit_certificate']
        if not (check['all_passed'] and check['actual_five_defect_family_sha256']==self.family
                and check['implicit_source_sha256']==self.source
                and ledger['uniform_Cstar_family_sha256']==self.core.records['physical_norm_family']['uniform_Cstar_family_sha256']
                and ledger['base_analytic_core_family_sha256']==self.core.records['core_transfer']['analytic_core_family_sha256']
                and gate['admitted_inner_parameter_family_sha256']==ledger['admitted_inner_parameter_family_sha256']
                and gate['actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified']
                and ledger['positive_width_not_materialized'] and ledger['h_b_equals_epsilon_b_by_definition']):
            raise ValueError('Same admitted width, core and exact bridge prescription required')
        self.logh=read_interval(c,ledger['shared_positive_width_log_enclosure'])
        self.cap=c.mpf('1e-180');self.cap_proofs=[];self.cache={}
        if endpoints(self.logh)[1]>=endpoints(c.ln(self.cap))[0]:
            raise ValueError('Formal source width does not admit enclosure-only cap')
        if not endpoints(2*self.cap)[1]<endpoints(c.mpf('.01'))[0]:
            raise ValueError('Core continuation width must cover smoothing')
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def width_product(self,value,extra_log=0):
        c=self.ctx;upper=norm(c,value)
        if endpoints(upper)[1]==0:return c.mpf(0)
        logbound=self.logh+c.mpf(extra_log)+c.ln(upper)
        if endpoints(logbound)[1]>endpoints(c.ln(self.cap))[0]:
            raise ArithmeticError('Width-weighted source bound needs a larger analytic enclosure')
        self.cap_proofs.append(dict(source_log_upper=logbound,cap=self.cap,
            input_absolute_upper=upper,additional_source_log=c.mpf(extra_log),
            exact_source='hb times positive source factor; cap does not define hb',passed=True))
        return c.mpf([0,endpoints(self.cap)[1]])

    def augment_core(self,packet,name):
        c=self.ctx;jet=axial_jet(c,packet[name])
        major=self.core.records['core_transfer'];linear=self.core.records['shared_linear_resolvent']
        if name=='Phi':bound=read_interval(c,linear['Phi_model_Xh_norm_upper'])+self.core.correction
        else:bound=(read_interval(c,major['fresh_Psi_model_Xh_norm_upper'])+self.core.correction)*self.core.epsilon
        sixth=bound*embedding(c,self.core.h,c.mpf('4.1'),0,6)
        return IntervalTaylor(c,list(jet.coefficients)+[symmetric(c,sixth)/math.factorial(6)])

    def inputs(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=Z._mpi_
        if key in self.cache:return self.cache[key]
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Z in[-1,1] required')
        at=self.core.normalized_jets(4,Z);radial=self.core.normalized_jets([0,4],Z)
        cont=self.core.normalized_jets([4,'4.1'],Z)
        phi=self.augment_core(at,'Phi');v=self.augment_core(at,'Uz');mean=self.augment_core(at,'Mz_over_R')
        pc=self.augment_core(radial,'Phi');vc=self.augment_core(radial,'Uz')
        lc=logarithm(self.augment_core(cont,'Phi'))
        # IBP: alpha*g(y)+int(-alpha')*g is a positive, Z-independent
        # average with total weight1. It encloses all six axial derivatives.
        barphi=lc.exp()
        barphi=IntervalTaylor(c,[intersection(c,barphi[0],c.mpf([endpoints(self.core.phi_floor)[0],endpoints(self.core.phi_ceiling)[1]]))]
            +list(barphi.coefficients[1:]))
        barv=self.augment_core(cont,'Uz')
        source=at['source'];H=source['H'];den=square(H)+self.core.sigma**2
        gradient=source['L']*H/den
        # Bell recurrence at order6; no materialization of F0 or Cstar.
        def ratios(multiplier):
            out=[c.mpf(1)]
            for n in range(1,7):
                out.append(sum((-multiplier*self.core.Lambda*gradient[j]*out[n-1-j] for j in range(n)),c.mpf(0))/n)
            return [value*math.factorial(k) for k,value in enumerate(out)]
        datum=self.core.datum.normalized_jets(endpoints(Z),6)
        pressure=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in datum['normalized_pressure_coefficients']])
        result=dict(phi=phi,v=v,mean=mean,H=pc,K=pc*vc,A=square(vc),B=square(pc)/2,C=square(pc),
            p0=pressure,F0_ratios=ratios(1),F0_squared_ratios=ratios(2),
            comparison_phi_cover=barphi,comparison_v_cover=barv,
            comparison_log_phi_cover=lc)
        self.cache[key]=result;return result

    def comparison(self,Z,theta):
        c=self.ctx;inputs=self.inputs(Z)
        phi=inputs['comparison_phi_cover'];v=inputs['comparison_v_cover']
        if endpoints(theta)==(mp.mpf(1),mp.mpf(1)):phi=inputs['phi'];v=inputs['v']
        moments=cumulative_enclosures(c,theta,inputs,inputs['comparison_phi_cover'],inputs['comparison_v_cover'])
        coefficients=direction(c,Z,self.delta,phi,v,moments,inputs['p0'],inputs['F0_ratios'],inputs['F0_squared_ratios'])
        return dict(phi=phi,v=v,moments=moments,direction=coefficients)

    def actual(self,Z,theta):
        c=self.ctx;Z=c.mpf(Z);theta=c.mpf(theta)
        floor=endpoints(self.r/100)[0]
        if endpoints(theta)[0]<floor or endpoints(theta)[1]>1:raise ValueError('Ra<=R<=100 required')
        inputs=self.inputs(Z);start=endpoints(theta)==(mp.mpf(1),mp.mpf(1))
        whole_theta=c.mpf([floor,1]);comparison=self.comparison(Z,whole_theta)
        direction_bounds=comparison['direction'];ymax=c.ln(100/self.r)
        weights=1+ymax
        ell=[]
        for k in range(6):
            bound=c.mpf(0) if start else self.width_product(direction_bounds['D_over_R'][k]*weights*50)
            ell.append(c.mpf([-endpoints(bound)[1],0]) if k==0 else symmetric(c,bound))
        relative=IntervalTaylor(c,ell).exp()
        actualphi=inputs['phi'].truncate(5)*relative
        quotient=actualphi/comparison['phi'].truncate(5)
        drive={name:quotient*value for name,value in direction_bounds.items() if name.startswith('drive_')}
        velocity=[]
        for k in range(6):
            bounds=[c.mpf(0)]*3 if start else [
                self.width_product(drive['drive_hydro'][k]*weights*100),
                self.width_product(drive['drive_pressure'][k]*weights*100,2*self.core.logP),
                self.width_product(drive['drive_swirl'][k]*weights*10000,-2*self.core.logC)]
            velocity.append(inputs['v'][k]+symmetric(c,sum(bounds,c.mpf(0))))
        actualv=IntervalTaylor(c,velocity)
        # Same uniform-in-radius field bounds enclose the actual integrated
        # histories at every radius. No comparison moments are substituted.
        initial={name:value.truncate(5) if isinstance(value,IntervalTaylor) else value for name,value in inputs.items()}
        moments=cumulative_enclosures(c,theta,initial,actualphi,actualv)
        z=IntervalTaylor.variable(c,Z,5);d=1-square(z);L=1-square(z)*self.delta;m=moments[MZ]
        Q=(2*z*actualv-(z*m)*(1-self.delta)-d*derivative(m))/L
        pressure=dress(moments[MP],inputs['F0_squared_ratios'])
        return dict(Z=Z,theta_r_over_R=theta,log_remaining_radial_length=ymax,
            log_F_actual_over_exit_F_coefficients=list(ell),
            F_actual_over_F0_axial5_coefficients=list(actualphi.coefficients),
            F_actual_true_axial5_divided_by_F0=list(dress(actualphi,inputs['F0_ratios']).coefficients),
            Uz_actual_axial5_coefficients=list(actualv.coefficients),
            actual_moment_shape_axial5_coefficients=coefficient_lists(moments),
            actual_Q_axial4_coefficients=list(Q.coefficients),
            pressure_axis_axial5_coefficients=list(inputs['p0'].truncate(5).coefficients),
            pressure_increment_true_axial5_divided_by_R_F0_squared=list(pressure.coefficients),
            actual_velocity_prefactors=dict(Ur='sqrt(R/2)*Q',Utheta='sqrt(2R)*F0*phi_actual',Uz='V_actual'),
            core_exit_exact=start,original_P0_retained=True,all_actual_moments_inherited=True,
            formal_positive_hb_used_not_numerical_cap=True,
            actual_bridge_radial_mixed4_certified=False,
            source_integral_definitions=dict(F='F=f*exp(-.5*integral_0^y chi*Dbar dt)',
                V='V=v-integral_0^y chi*(phi_actual/phi_bar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl) dt',
                chi='1-(1-hb)*sigma(y/hb); hb=cstar*K^-100',
                comparison='alpha=1-sigma((y-hb)/hb); logFbar andVbar integrate alpha*core radial derivatives; own five moments'))

    def report(self):
        c=self.ctx;theta=c.mpf([endpoints(self.r/100)[0],1]);whole=self.actual([-1,1],theta)
        comparison=self.comparison(c.mpf([-1,1]),theta)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            admitted_inner_parameter_family_sha256=self.records['K1_ledger']['admitted_inner_parameter_family_sha256'],
            exact_width_source='hb=epsilon_b=cstar*K^-100 with exact physical norm sumK(9.16)',
            source_log_hb_enclosure=self.logh,positive_width_materialized=False,cap_is_not_width=True,
            whole_bridge=whole,exit=self.actual([-1,1],1),terminal=self.actual([-1,1],self.r/100),
            samples=[self.actual(z,self.r/R) for z in ('0','.5') for R in ('1','100')],
            comparison_phi_axial6_coefficients=list(comparison['phi'].coefficients),
            comparison_v_axial6_coefficients=list(comparison['v'].coefficients),
            comparison_own_moment_axial6_coefficients=coefficient_lists(comparison['moments']),
            comparison_direction_axial5_coefficients=coefficient_lists(comparison['direction']),
            comparison_IBP_weights='alpha*g(y)+integral(-alpha_prime)*g; nonnegative, sum1, independent of Z',
            integrated_chi_bound='integral_0^y chi dt<=hb*(1+y); early length<=hb, laterchi=hb',
            source_log_product_cap_proofs=self.cap_proofs,
            source_bound_smoothed_comparison_axial6_enclosures_available=True,
            actual_prescribed_shear_bridge_axial5_enclosures_available=True,
            actual_bridge_radial_recovery_axial4_enclosures_available=True,
            actual_bridge_radial_mixed4_certified=False,bridge_100_110_switches_installed=False,
            newly_recomputed_point_coefficients=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,
            next_dependency='Actual bridge phase/radial mixed derivatives with formal inverse-hb prefactors, core functional joins, original100..110 switches and110 inlet -> reshape/restoration/moment patch',
            input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantInnerBridgeProfiles().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Same-source smoothed comparison axial6 and actual prescribed-shear bridge axial5 enclosures generated',flush=True)
    return result


if __name__=='__main__':run()
