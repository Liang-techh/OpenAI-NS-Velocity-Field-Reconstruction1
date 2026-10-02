"""Current-source anchored axis primitive and logarithmic core values.

Six certified complex poles evaluate the exact rational G primitive. Root
centers are provisional; uniform directed Rouche disks enclose roots for all
admitted parameter values. Original anchor, Cstar, pressure and fresh radial
rows remain shared. exp(logF0) is deliberately never materialized.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_core_coefficient_rebuild import CompliantCoreCoefficientRebuild
from lei_ren_part1_paper_candidate_exact_amplitude import (
    _H,_Hprime,_L,_point_complex,_disk_complex,_encode)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def contains_zero(value):
    lo,hi=endpoints(value)
    return lo<=0<=hi


class CompliantAnchoredAxisAmplitude:
    def __init__(self):
        self.rebuild=CompliantCoreCoefficientRebuild();self.core=self.rebuild.core;self.ctx=c=self.core.ctx
        self.hashes=dict(self.rebuild.hashes)
        name=PREFIX+'core_coefficient_rebuild_check.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or receipt['implicit_source_sha256']!=self.core.source:
            raise ValueError('Accepted fresh compliant core rows required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Fresh core input changed: '+path)
        self.hashes.update(receipt['input_hashes'])
        for path in (name,Path(__file__).name,'lei_ren_part1_paper_candidate_exact_amplitude.py'):
            self.hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
        self.params=dict(j=self.core.j,delta=self.core.delta,sigma=self.core.sigma,a=(9-self.core.delta)/2)
        self.radius=c.mpf('1e-180');self.rows={};self.cache={}
        with mp.workdps(400):self.roots=self.root_data()

    def certify_root(self,center,sign):
        c=self.ctx;r=self.radius;point=_point_complex(c,center)
        target=c.mpc(0,sign*self.params['sigma']);disk=_disk_complex(c,center,endpoints(r)[1])
        residual=_H(c,point,self.params,target);derivative=_Hprime(c,point,self.params)
        # Exact cubic Taylor remainder about the provisional center.
        lhs=abs(residual)+(12*abs(point)+abs(self.params['j']))*r**2+4*r**3
        rhs=abs(derivative)*r
        if endpoints(lhs)[1]>=endpoints(rhs)[0]:raise ArithmeticError('Uniform source Rouche disk failed')
        return dict(target_sign=sign,nominal_center=point,disk=disk,radius=r,
            residual_abs_upper=abs(residual),linear_derivative_abs_lower=abs(derivative),
            cubic_remainder_bound=(12*abs(point)+abs(self.params['j']))*r**2+4*r**3,
            rouche_lhs=lhs,rouche_rhs=rhs,strict_uniform_rouche=True,
            provisional_center_is_not_selected_source_root=True)

    def root_data(self):
        c=self.ctx
        # Centers only guide root disks. The subsequent certificate uses the
        # ORIGINAL interval parameters, not these provisional scalar values.
        nominal=lambda value:sum(endpoints(value))/2
        j=nominal(self.params['j']);delta=nominal(self.params['delta']);sigma=j/500
        a=(9-delta)/2
        nominal_roots=lambda target:mp.polyroots([-4,-j,a,j-target],maxsteps=2000,error=False)
        z0=min(nominal_roots(0),key=abs)
        if mp.im(z0)!=0:raise ArithmeticError('Nominal anchor is not real')
        anchor=self.certify_root(mp.mpc(z0),0);anchor_interval=anchor['disk'].real
        jhi=endpoints(self.params['j'])[1];bracket=c.mpf([-jhi,0])
        left=_H(c,c.mpc(c.mpf(-jhi),0),self.params).real
        right=_H(c,c.mpc(0,0),self.params).real
        derivative=_Hprime(c,c.mpc(bracket,0),self.params).real
        if not (endpoints(left)[1]<0<endpoints(right)[0] and endpoints(derivative)[0]>0):
            raise ArithmeticError('Original unique real anchor bracket failed')
        if not (-jhi<=endpoints(anchor_interval)[0]<=endpoints(anchor_interval)[1]<=0):
            raise ArithmeticError('Certified anchor disk outside original bracket')
        x=c.mpf([-1,1])-anchor_interval
        # H(a+x)=x*[H'(a)+(-12a-j)x-4x^2], using H(a)=0.
        # A positive factor proves the sign of H across the full real domain,
        # and hence the anchored primitive is nonnegative on both sides.
        sign_factor=_Hprime(c,c.mpc(anchor_interval,0),self.params).real+(-12*anchor_interval-self.params['j'])*x-4*x**2
        if endpoints(sign_factor)[0]<=0:raise ArithmeticError('Whole real domain H/root sign factor is not positive')
        L_lower=1-self.params['delta']
        if endpoints(self.params['delta'])[0]<0 or endpoints(L_lower)[0]<=0:
            raise ArithmeticError('L is not positive across the original real domain')
        poles=[]
        for center in nominal_roots(mp.mpc(0,sigma)):
            for sign,root in ((1,mp.mpc(center)),(-1,mp.mpc(center).conjugate())):
                packet=self.certify_root(root,sign)
                lo,hi=endpoints(packet['disk'].imag)
                if not (lo>0 or hi<0):raise ArithmeticError('Certified pole disk meets real path')
                packet['imaginary_sign_certified']=True;poles.append(packet)
        all_disks=[anchor]+poles;separations=[]
        for i,left in enumerate(all_disks):
            for right in all_disks[i+1:]:
                margin=abs(left['nominal_center']-right['nominal_center'])-2*self.radius
                if endpoints(margin)[0]<=0:raise ArithmeticError('Certified root disks are not disjoint')
                separations.append(margin)
        return dict(anchor=anchor,anchor_interval=anchor_interval,anchor_bracket=bracket,
            anchor_Hprime_lower=derivative,unique_real_anchor_certified=True,
            full_real_domain_H_sign_factor=sign_factor,G_nonnegative_on_real_domain=True,
            full_real_domain_L_lower=L_lower,
            root_real_by_conjugation_and_disk_uniqueness=True,poles=poles,
            pairwise_separation_margins=separations,all_six_poles_certified=True,
            uniform_interval_parameter_root_certificates=True,
            nominal_parameters_used_only_for_provisional_centers=True)

    def gradient_jets(self,Z,order=5):
        c=self.ctx;z=IntervalTaylor.variable(c,Z,order);one=IntervalTaylor.constant(c,1,order)
        H=-4*z**3-z*z*self.params['j']+z*self.params['a']+self.params['j']
        denominator=H*H+self.params['sigma']**2
        coefficients=list(denominator.coefficients);coefficients[0]=H[0]**2+self.params['sigma']**2
        return ((one-z*z*self.params['delta'])*H/IntervalTaylor(c,coefficients)).coefficients

    def evaluate(self,Z):
        c=self.ctx;anchored=isinstance(Z,str) and Z=='anchor'
        z=self.roots['anchor_interval'] if anchored else c.mpf(Z)
        if endpoints(z)[0]<-1 or endpoints(z)[1]>1:raise ValueError('Original real axial domain [-1,1] required')
        key=(anchored,z._mpi_)
        if key in self.cache:return self.cache[key]
        anchor=self.roots['anchor_interval'];G=c.mpc(0,0);derivative=c.mpc(0,0);terms=[]
        for pole in self.roots['poles']:
            r=pole['disk'];residue=_L(c,r,self.params)/(2*_Hprime(c,r,self.params))
            numerator=c.mpc(z,0)-r;denominator=c.mpc(anchor,0)-r
            # Each argument remains in one open imaginary half-plane over
            # the complete real path. Principal logs individually therefore
            # have a continuous difference with zero value at the anchor.
            for argument in (numerator,denominator):
                lo,hi=endpoints(argument.imag)
                if not (lo>0 or hi<0):raise ArithmeticError('Log primitive lost half-plane branch certificate')
            logdifference=c.log(numerator)-c.log(denominator)
            G+=residue*logdifference;derivative+=residue/numerator
            terms.append(dict(target_sign=pole['target_sign'],pole=r,residue=residue,
                log_difference=logdifference,individual_logs_stay_in_fixed_halfplane=True))
        direct=_L(c,c.mpc(z,0),self.params)*_H(c,c.mpc(z,0),self.params)/(
            _H(c,c.mpc(z,0),self.params)**2+self.params['sigma']**2)
        difference=derivative-direct
        if not contains_zero(G.imag) or not contains_zero(difference.real) or not contains_zero(difference.imag):
            raise ArithmeticError('Anchored rational primitive consistency check failed')
        actual_G=c.mpf(0) if anchored else c.mpf([max(mp.mpf(0),endpoints(G.real)[0]),endpoints(G.real)[1]])
        logF0=-self.core.logC-self.core.Lambda*actual_G
        result=dict(Z=z,is_shared_anchor=anchored,G=actual_G,complex_primitive_diagnostic=G,
            logF0=logF0,log_F0_relative_to_anchor=-self.core.Lambda*actual_G,
            selected_logCstar=self.core.logC,original_anchored_primitive='integral_a^Z L H/(H^2+sigma^2), H(a)=0 in[-j,0]',
            signed_pole_log_terms=terms,Gprime_from_poles=derivative,Gprime_direct=direct,
            Gprime_difference=difference,G_gradient_taylor_coefficients=self.gradient_jets(z),
            anchor_exact_zero_used=anchored,exact_source_anchor_retained=True,
            root_parameter_uncertainty_propagated=True,log_branches_certified=True,
            exp_logF0_not_materialized=True,positive_F0_retained_formally=True,
            full_point_physical_field_evaluation=False,temporal_recursion=False)
        self.cache[key]=result;return result

    def core_value(self,Z,rho,log_tau='-1',theta='0'):
        c=self.ctx;amplitude=self.evaluate(Z);z=amplitude['Z'];r=c.mpf(rho)
        if not (-1<endpoints(z)[0]<=endpoints(z)[1]<1):raise ValueError('Finite physical core point requires |Z|<1')
        key=z._mpi_
        if key not in self.rows:self.rows[key]=self.rebuild.rebuild(z)
        profile=self.rebuild.values(self.rows[key],r)
        logtau=c.mpf(log_tau)
        if not all(mp.isfinite(v) for v in endpoints(logtau)):raise ValueError('Finite log(tau) required')
        loglambda=(logtau-c.ln(1-z**2))/2;delta=self.core.delta
        component=lambda coefficient,scale:dict(coefficient_enclosure=coefficient,positive_scale_log=scale,
            value_definition='coefficient_enclosure * exp(positive_scale_log)',source_value_enclosure_not_point_selection=True)
        axis=endpoints(r)==(mp.mpf(0),mp.mpf(0))
        if axis:
            radial=swirl=component(c.mpf(0),c.mpf(0));log_radius=None
        else:
            if endpoints(profile['Phi'])[0]<=0:raise ArithmeticError('Fresh core normalized swirl positivity unresolved')
            logR=c.ln(r)-self.core.logLambda
            radial=component(profile['radial_recovery_Q'],(logR-c.ln(2))/2-loglambda)
            swirl=component(profile['Phi'],amplitude['logF0']+(logR+c.ln(2))/2-(1+delta)*loglambda)
            log_radius=loglambda+(logR+c.ln(2))/2
        axial=component(profile['Uz'],-(1+delta)*loglambda)
        pressure=component(profile['P_scaled'],self.core.logLambda-2*(1+delta)*loglambda)
        cosine=c.cos(c.mpf(theta));sine=c.sin(c.mpf(theta))
        return dict(Z=z,rho=r,requested_log_tau=logtau,log_lambda=loglambda,physical_log_r=log_radius,
            physical_z=component(z,(1-delta)*loglambda),axis_exact_zero=axis,
            cylindrical_components=dict(ur=radial,utheta=swirl,uz=axial,p=pressure),
            cartesian_velocity_terms=dict(ux=[component(cosine*radial['coefficient_enclosure'],radial['positive_scale_log']),
                                              component(-sine*swirl['coefficient_enclosure'],swirl['positive_scale_log'])],
                uy=[component(sine*radial['coefficient_enclosure'],radial['positive_scale_log']),
                    component(cosine*swirl['coefficient_enclosure'],swirl['positive_scale_log'])],uz=[axial]),
            amplitude_packet=amplitude,core_profile_packet=profile,
            output_kind='Physical core source-value enclosures in signed logarithmic factors; no float underflow substitution',
            fresh_core_and_actual_anchored_amplitude_combined=True,
            full_point_physical_field_evaluation=False,all_annular_source_values_resolved=False,
            measured_blowup_dynamics=False,admissible_stress_lift_constructed=False,temporal_recursion=False)

    def report(self):
        samples={Z:self.evaluate(Z) for Z in ('-1','-.5','-.3','0','.3','.5','1','anchor')}
        values={Z:[self.core_value(Z,r) for r in ('0','2','4.1')] for Z in ('-.5','0','.3','.5')}
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,root_certificates=self.roots,
            anchored_amplitude_packets=samples,physical_core_value_packets=values,
            same_selected_Cstar_and_current_source=True,actual_anchored_G_resolved_with_directed_error=True,
            positive_F0_retained_in_logarithmic_form=True,original_candidate_fixed_parameters_not_used=True,
            full_point_physical_field_evaluation=False,all_annular_source_values_resolved=False,
            measured_blowup_dynamics=False,physical_energy_integral_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(400):result=CompliantAnchoredAxisAmplitude().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Current anchored G and logarithmic physical core source values generated',flush=True)
    return result


if __name__=='__main__':run()
