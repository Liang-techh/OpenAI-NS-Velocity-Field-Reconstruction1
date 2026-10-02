"""Entire post-flatten power and actual two-support angular mixed C4.

Complete future energy is factored BEFORE cancelling q^2. A backward
positive-tail representation keeps the original angular repair, epsilon
atoms and entire Gamma energy. It does not subtract the long power interval
from an unrelated whole-Z future-energy enclosure or reset any history.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import CompliantFlattenMixedC4,flatten_mixed
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import FlatPulseDerivatives
from lei_ren_part1_paper_compliant_outer_angular_repair import bump_weights
from lei_ren_part1_paper_compliant_corrected_outer_field import future_bump_weights
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_compliant_angular_high_jets import log_taylor
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def quotient_log_rates(rows):
    """True y derivatives0..3 of F_y/F from F derivatives0..4."""
    if len(rows)!=5 or any(v.order!=5 for v in rows):raise ValueError('F axial5 at y orders0..4 required')
    if endpoints(rows[0][0])[0]<=0:raise ValueError('Positive actual swirl factor required')
    rates=[]
    for n in range(4):
        value=rows[n+1]
        for j in range(1,n+1):value-=rows[j]*rates[n-j]*math.comb(n,j)
        rates.append(value/rows[0])
    return rates


class CompliantPowerAngularC4:
    def __init__(self,cells=128):
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive directed integration cell count required')
        self.flatten=CompliantFlattenMixedC4(cells); self.ctx=c=self.flatten.ctx
        self.mu=self.flatten.mu; self.delta=self.flatten.delta; self.rate=1-self.mu
        self.bp=c.mpf('.5')+self.mu; self.prate=1+2*self.mu
        self.family=self.flatten.family; self.source=self.flatten.source
        self.cells=cells; self.hashes=dict(self.flatten.hashes); self.cache={}
        name=PREFIX+'compliant_flatten_mixed_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['full_pulse_C4_installed'] or not receipt['entire_original_100_unit_flatten_mixed_C4_available']:
            raise ValueError('Admitted local pulse/entire flatten prerequisites required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power/angular source changed: '+source)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        ename=PREFIX+'compliant_future_swirl_energy.json'; energy=json.loads((HERE/ename).read_bytes())
        aname=PREFIX+'compliant_outer_angular_repair.json'; angular=json.loads((HERE/aname).read_bytes())
        if any(r['actual_five_defect_family_sha256']!=self.family or r['implicit_source_sha256']!=self.source for r in (energy,angular)):
            raise ValueError('Actual angular/complete future source mismatch')
        self.hashes[ename]=hashlib.sha256((HERE/ename).read_bytes()).hexdigest()
        self.hashes[aname]=hashlib.sha256((HERE/aname).read_bytes()).hexdigest()
        self.Lrel=read_interval(c,energy['postflatten_length'])
        if endpoints(self.Lrel)[0]<=4:raise ValueError('Original power/last4 chart split not available')
        pref={k:read_interval(c,v) for k,v in energy['relative_energy_prefactors'].items()}
        kernels=energy['postrel_energy_kernels']; atoms=energy['whole_Z_C1']['separate_preheat_tail_atoms']
        self.tail_multiplier=read_interval(c,energy['tail_multiplier'])
        self.post_scalar=read_interval(c,kernels['steep_in'])+pref['Ns']*read_interval(c,energy['steep_power_energy_integral'])
        self.post_scalar+=pref['Nq']*read_interval(c,kernels['steep_out'])+pref['Nt']*read_interval(c,energy['waiting_energy_integral'])
        self.post_scalar+=self.tail_multiplier*sum((read_interval(c,v) for v in atoms.values()),c.mpf(0))
        self.S_cap=read_interval(c,angular['strong_inverse_radius_positive_cap'])
        self.heatcap=self.tail_multiplier*(self.delta/2)*self.S_cap
        self.weights={k:read_interval(c,v) for k,v in angular['bump_weights'].items()}
        self.flat=FlatPulseDerivatives(c); self.normalization=self.flat.normalization
        self.fifth=self.flatten.fifth
        self.hashes.update(self.flat.hashes)
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def data(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Z in[-1,1] required')
        if key in self.cache:return self.cache[key]
        source=self.fifth['whole_Z']['angular']
        for sample in self.fifth['samples']:
            if endpoints(read_interval(c,sample['angular']['Z']))==key:source=sample['angular']; break
        jet=lambda v:IntervalTaylor(c,[read_interval(c,x) for x in v['coefficients']])
        coeff=[jet(v) for v in source['physical_coefficient_Taylor']]
        heat=jet(source['scaled_Gamma_future_defect_Taylor']['energy'])
        post=IntervalTaylor.constant(c,self.post_scalar,5)-heat*c.mpf([0,endpoints(self.heatcap)[1]])
        full_change=post*0
        for dj,center in zip(coeff,(-3,-1)):
            full_change+=(dj*(2*self.weights['E'])+dj*dj*self.weights['F'])*c.exp(-2*self.mu*center)
        f=self.flatten.flatten(Z,100)
        if endpoints(post[0])[0]<=0 or endpoints((post+full_change)[0])[0]<=0:
            raise ArithmeticError('Complete post-angular future energy positivity lost')
        value=dict(coeff=coeff,post=post,full_angular_change=full_change,
            flatten_exit_X=f['angular_Taylor'],flatten_exit_pressure=f['pressure'],
            flatten_exit_energy=f['energy_Taylor'],complete_future_energy_input=self.flatten.future_energy(Z))
        self.cache[key]=value; return value

    def _packet(self,Z,theta,X,energy,pressure,rates,coordinate,extra):
        c=self.ctx; one=IntervalTaylor.constant(c,1,5)
        mixed=flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)
        if endpoints(energy[0])[0]<=0:raise ArithmeticError('Backward complete future energy no longer positive')
        mixed.update(Z=c.mpf(Z),coordinate=coordinate,theta_over_Ev0_Taylor=theta,
            angular_Taylor=X,energy_Taylor=energy,pressure_over_Pstar_squared_Taylor=pressure,
            log_theta_y_rates_Taylor=rates,
            actual_terminal_zero_linear_and_radial_histories_inherited=True,
            original_pressure_datum_and_positive_complete_future_energy_preserved=True,
            correlated_q_squared_factors_cancelled_before_enclosure=True,
            exact_positive_Ev0Pstar2_log_parts=self.flatten.logEv2_parts,
            positive_Ev0Pstar2_cap=self.flatten.Ev2,
            velocity_units='U/Ev0; exact formal Ev0/Pstar retained separately',
            pressure_units='P/Pstar^2; enclosure of original P0+Mp',
            full_pulse_C4_installed=True,full_outer_C4_certified=False,
            physical_energy_integral_certified=False,whole_outer_cone_certified=False,
            temporal_recursion=False,**extra)
        return mixed

    def power(self,Z,phase):
        c=self.ctx; phase=c.mpf(phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Power phase in[0,1] required')
        data=self.data(Z); y=(self.Lrel-4)*phase
        D=4+(self.Lrel-4)*(1-phase)
        energy=(data['post']+data['full_angular_change'])*(c.exp(-2*self.mu*D)/2)
        energy+=decay_integral(c,2*self.mu,D)/2
        one=IntervalTaylor.constant(c,1,5); X=one/self.rate+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)
        theta=one*(c.exp(-100*self.bp-self.bp*y)/2)
        deltaP=decay_integral(c,self.prate,y)*(self.flatten.Ev2*c.exp(-100*self.prate)/8)
        pressure=data['flatten_exit_pressure']['P_over_Pstar_squared']+deltaP
        return self._packet(Z,theta,X,energy,pressure,[one*0 for _ in range(4)],
            dict(origin='Rf',phase=phase,offset=y,remaining_to_Rrel=D),
            dict(entire_following_power_high_mixed_derivatives_available=True,
                energy_representation='[D(2mu,Lrel-y)+exp(-2mu*(Lrel-y))*(Post+full_angular_change)]/2',
                exact_relative_velocity_log_parts=dict(flatten_exit=-100*self.bp-c.ln(2),power_offset=-self.bp*y),
                no_subtraction_of_long_power_from_unrelated_future_energy_boxes=True))

    def angular(self,Z,s):
        c=self.ctx; s=c.mpf(s)
        if endpoints(s)[0]<-4 or endpoints(s)[1]>0:raise ValueError('Angular chart s in[-4,0] required')
        data=self.data(Z); zero=data['post']*0; h=[zero for _ in range(5)]
        pastA=zero; pastP=zero; futureE=zero
        for dj,center in zip(data['coeff'],(-3,-1)):
            local=s-center; beta=self.flat.beta(local)
            for k in range(5):h[k]+=dj*(beta[k]*math.factorial(k))
            lo,hi=endpoints(local); ell=c.mpf('.15')
            if hi<=-endpoints(ell)[1]:past={k:c.mpf(0) for k in self.weights}; future=self.weights
            elif lo>=endpoints(ell)[1]:past=self.weights; future={k:c.mpf(0) for k in self.weights}
            else:
                past=bump_weights(c,self.mu,self.normalization,local,cells=self.cells)
                future=future_bump_weights(c,self.mu,self.normalization,local,cells=self.cells)
            pastA+=dj*(c.exp(self.rate*center)*past['A'])
            pastP+=(dj*past['B']+dj*dj*past['D']/2)*c.exp(-self.prate*center)
            futureE+=(dj*(2*future['E'])+dj*dj*future['F'])*c.exp(-2*self.mu*center)
        F=[h[0]+1]+h[1:]; rates=quotient_log_rates(F)
        y=self.Lrel+s; eq=IntervalTaylor.constant(c,1/self.rate,5)
        Xbase=eq+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)
        X=(Xbase+pastA*c.exp(-self.rate*s))/F[0]
        energy=((data['post']+futureE)*c.exp(2*self.mu*s)+decay_integral(c,2*self.mu,-s))/(F[0]*F[0]*2)
        theta=F[0]*(c.exp(-100*self.bp-self.bp*self.Lrel-self.bp*s)/2)
        baseline=decay_integral(c,self.prate,y)*(self.flatten.Ev2*c.exp(-100*self.prate)/8)
        pressure=data['flatten_exit_pressure']['P_over_Pstar_squared']+baseline
        pressure+=pastP*(self.flatten.Ev2*c.exp(-self.prate*(100+self.Lrel))/4)
        return self._packet(Z,theta,X,energy,pressure,rates,dict(origin='Rrel',offset=s),
            dict(actual_angular_source_C5_coefficients=data['coeff'],
                actual_angular_bump_y_derivatives=h,
                swirl_factor_one_plus_h_Taylor=F[0],actual_past_angular_change=pastA,
                actual_past_pressure_change=pastP,actual_future_energy_change=futureE,
                complete_postangular_future_energy_Taylor=data['post'],
                both_actual_angular_supports_high_mixed_derivatives_available=True,
                disjoint_beta_support_cross_energy_exact_zero=True,
                angular_coefficients_not_axial_pulse_coefficients=True,
                energy_representation='[exp(2mu*s)*(Post+future_angular_change)+D(2mu,-s)]/[2*(1+h)^2]',
                exact_relative_velocity_log_parts=dict(flatten_exit=-100*self.bp-c.ln(2),power_length=-self.bp*self.Lrel,
                    angular_offset=-self.bp*s,actual_repair_log_Taylor=log_taylor(F[0]))))

    def report(self):
        with mp.workdps(270):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                selected_mu=self.mu,selected_delta=self.delta,postflatten_length=self.Lrel,
                scalar_postangular_complete_future_energy=self.post_scalar,
                exact_positive_heat_factor_cap=self.heatcap,
                power_samples=[self.power(z,p) for z in ('-1','0','.5','1') for p in ('0','.5','1')],
                angular_samples=[self.angular(z,s) for z in ('-1','0','.5','1') for s in ('-4','-3','-2','-1','0')],
                whole_Z_power_box=self.power([-1,1],[0,1]),whole_Z_power_inlet=self.power([-1,1],0),
                whole_Z_power_terminal=self.power([-1,1],1),
                whole_Z_angular_box=self.angular([-1,1],[-4,0]),
                whole_Z_angular_inlet=self.angular([-1,1],-4),whole_Z_angular_terminal=self.angular([-1,1],0),
                angular_support_crossings=[self.angular([-1,1],[center+side*.15-.01,center+side*.15+.01])
                    for center in (-3,-1) for side in (-1,1)],
                entire_following_power_high_mixed_derivatives_available=True,
                both_actual_angular_supports_high_mixed_derivatives_available=True,
                original_signed_angular_changes_and_infinite_heat_energy_retained=True,
                flatten_power_and_power_angular_joins_certified=False,
                full_pulse_C4_installed=True,full_outer_C4_certified=False,
                physical_energy_integral_certified=False,whole_outer_cone_certified=False,
                temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantPowerAngularC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Entire O6 power and both actual angular supports: correlated positive future energy and mixed order<=4 generated',flush=True)
    return result


if __name__=='__main__':run()
