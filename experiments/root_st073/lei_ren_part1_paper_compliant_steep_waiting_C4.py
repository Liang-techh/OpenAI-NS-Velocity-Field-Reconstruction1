"""Original O7 steep entry/power/exit and waiting, mixed spatial order<=4.

The same angular-terminal histories and exact epsilon/Gamma future source
are transported throughout. Backward energy is normalized BEFORE enclosure
so the enormous steep-power attenuation never has to be divided back out.
This is leading-profile spatial recovery, not temporal coefficient recursion.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_power_angular_C4 import CompliantPowerAngularC4
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_axial_pulse_field import sigma_enclosure,decay_integral
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def transition_kernels(c,t,mu,delta,kind,cells=128):
    """Directed original J, forward angular/pressure and backward energy.

    Every cell has length t/cells or (1-t)/cells, preserving its positivity
    even for an entire interval of endpoints. No difference of independent
    endpoint boxes represents a length or an exponential weight.
    """
    t=c.mpf(t)
    if kind not in ('in','out') or endpoints(t)[0]<0 or endpoints(t)[1]>1:
        raise ValueError('Original unit transition kind/domain required')
    if not isinstance(cells,int) or cells<1:raise ValueError('Positive cell count required')
    rate=1-mu if kind=='in' else 1-delta/2
    J=c.mpf(0); angular=c.mpf(0); pressure=c.mpf(0)
    def cell(a,b):
        return c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(1),endpoints(b)[1])])
    if endpoints(t)[1]>0:
        ds=t/cells
        for i in range(cells):
            a=t*i/cells; b=t*(i+1)/cells; v=cell(a,b)
            nextJ=J+ds*sigma_enclosure(c,v)
            jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
            if kind=='in':
                angular+=ds*c.exp(rate*(v-jc))
                pressure+=ds*c.exp(-(1+2*mu)*v-2*rate*jc)/2
            else:
                angular+=ds*c.exp(rate*jc)
                pressure+=ds*c.exp(-3*v+2*rate*jc)/2
            J=nextJ
    if endpoints(t)==(mp.mpf(1),mp.mpf(1)):J=c.mpf('.5')
    endpointJ=J; remaining=c.mpf(0); length=1-t
    if endpoints(length)[1]>0:
        ds=length/cells
        for i in range(cells):
            a=t+length*i/cells; b=t+length*(i+1)/cells; v=cell(a,b)
            nextJ=J+ds*sigma_enclosure(c,v)
            jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
            if kind=='in':remaining+=ds*c.exp(-2*mu*v-2*rate*jc)
            else:remaining+=ds*c.exp(-2*v+2*rate*jc)
            J=nextJ
    return dict(J=endpointJ,angular=angular,pressure=pressure,remaining_energy=remaining,
        exact_J_definition='integral_0^t original_sigma(v)dv; J(1)=1/2',
        positive_correlated_cell_lengths=True)


class CompliantSteepWaitingC4:
    def __init__(self,cells=128):
        self.outer=CompliantPowerAngularC4(cells); self.ctx=c=self.outer.ctx
        self.mu=self.outer.mu; self.delta=self.outer.delta; self.rate=1-self.mu
        self.k=1-self.delta/2; self.bp=c.mpf('.5')+self.mu; self.bh=c.mpf('.5')+self.delta/2
        self.family=self.outer.family; self.source=self.outer.source; self.cells=cells
        self.hashes=dict(self.outer.hashes); self.cache={}; self.kernel_cache={}
        name=PREFIX+'compliant_power_angular_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['flatten_power_and_power_angular_joins_certified']:
            raise ValueError('Accepted same-source angular-terminal provider required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Steep prerequisite changed: '+source)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        name=PREFIX+'compliant_outer_angular_candidate.json'; waiting=json.loads((HERE/name).read_bytes())
        if waiting['actual_five_defect_family_sha256']!=self.family or waiting['implicit_source_sha256']!=self.source:
            raise ValueError('Waiting root family mismatch')
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.wait=read_interval(c,waiting['refined_same_source_waiting_root'])
        self.logone=read_interval(c,waiting['waiting_log_one_minus_epsilon'])
        # This is the original source Ts=4(log2-logdelta). Epsilon is not
        # read from the old parameter class's legacy .01 relation.
        self.Ts=4*(c.ln(2)-c.ln(self.delta)); self.epsilon=c.mpf('.001')*self.delta
        energy=json.loads((HERE/(PREFIX+'compliant_future_swirl_energy.json')).read_bytes())
        self.atoms={k:read_interval(c,v) for k,v in energy['whole_Z_C1']['separate_preheat_tail_atoms'].items()}
        self.preheat_scalar=sum(self.atoms.values(),c.mpf(0))
        self.tail_normalization=c.exp(-2*self.logone)
        self.S=read_interval(c,json.loads((HERE/(PREFIX+'compliant_outer_angular_repair.json')).read_bytes())['strong_inverse_radius_positive_cap'])
        self.infull=self.kernels(1,'in'); self.outfull=self.kernels(1,'out')
        self.inenergy=self.kernels(0,'in')['remaining_energy']; self.outenergy=self.kernels(0,'out')['remaining_energy']
        self.thetaR=c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2
        self.thetaS=self.thetaR*c.exp(-self.bp-self.rate/2)
        self.thetaQ=self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)
        self.thetaT=self.thetaQ*c.exp(-c.mpf('1.5')+self.k/2)
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def kernels(self,t,kind):
        key=(kind,tuple(endpoints(self.ctx.mpf(t))))
        if key not in self.kernel_cache:self.kernel_cache[key]=transition_kernels(self.ctx,t,self.mu,self.delta,kind,self.cells)
        return self.kernel_cache[key]

    def data(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if key in self.cache:return self.cache[key]
        terminal=self.outer.angular(Z,0)
        # The original infinite Gamma deficit is an actual C5 function.
        # S is positive in its definition; [0,S_cap] is only its enclosure.
        source=self.outer.fifth['whole_Z']['angular']
        for sample in self.outer.fifth['samples']:
            if endpoints(read_interval(c,sample['angular']['Z']))==key:source=sample['angular']; break
        heat=IntervalTaylor(c,[read_interval(c,v) for v in source['scaled_Gamma_future_defect_Taylor']['energy']['coefficients']])
        preheat_tail=IntervalTaylor.constant(c,self.preheat_scalar,5)-heat*c.mpf([0,endpoints(self.delta*self.S/2)[1]])
        H=preheat_tail*self.tail_normalization
        waiting_future=H*c.exp(-self.delta*self.wait)+decay_integral(c,self.delta,self.wait)
        after_power=waiting_future*c.exp(-1-self.delta/2)+self.outenergy
        if endpoints(H[0])[0]<=0 or endpoints(after_power[0])[0]<=c.mpf('.5'):
            raise ArithmeticError('Original heat tail / uniform power energy floor not proved')
        after_entry=(after_power*c.exp(-2*self.Ts)+decay_integral(c,2,self.Ts))*c.exp(-1-self.mu)
        XR=terminal['angular_Taylor']; PR=terminal['pressure_over_Pstar_squared_Taylor']
        XS=(XR+self.infull['angular'])*c.exp(-self.rate/2)
        XQ=XS+self.Ts; XT=(XQ+self.outfull['angular'])*c.exp(-self.k/2)
        PS=PR+self.infull['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)
        PQ=PS+decay_integral(c,3,self.Ts)*(self.outer.flatten.Ev2*self.thetaS**2/2)
        PT=PQ+self.outfull['pressure']*(self.outer.flatten.Ev2*self.thetaQ**2)
        item=dict(H=H,waiting_future=waiting_future,after_power=after_power,after_entry=after_entry,
            XR=XR,XS=XS,XQ=XQ,XT=XT,PR=PR,PS=PS,PQ=PQ,PT=PT,
            original_angular_terminal=terminal)
        self.cache[key]=item; return item

    def packet(self,Z,theta,X,energy,pressure,rates,stage,coordinate,logs,extra=None):
        c=self.ctx; one=IntervalTaylor.constant(c,1,5)
        result=self.outer._packet(Z,theta,X,energy,pressure,[one*r for r in rates],coordinate,
            dict(stage=stage,entire_steep_waiting_high_mixed_derivatives_available=True,
                same_actual_angular_terminal_histories_retained=True,
                exact_Gamma_and_both_epsilon_atoms_retained=True,
                exact_relative_velocity_log_parts=logs,
                no_forward_subtraction_of_unrelated_long_future_energy=True))
        if extra:result.update(extra)
        return result

    def steep_in(self,Z,t):
        c=self.ctx; t=c.mpf(t); data=self.data(Z); kernels=self.kernels(t,'in'); J=kernels['J']
        one=IntervalTaylor.constant(c,1,5)
        theta=one*self.thetaR*c.exp(-self.bp*t-self.rate*J)
        X=(data['XR']+kernels['angular'])*c.exp(-self.rate*(t-J))
        energy=(data['after_entry']+kernels['remaining_energy'])*c.exp(2*self.mu*t+2*self.rate*J)/2
        pressure=data['PR']+kernels['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)
        sig=sigma_jets(c,t); rates=[-self.rate*sig[j]*math.factorial(j) for j in range(4)]
        return self.packet(Z,theta,X,energy,pressure,rates,'O7 steep entry',dict(origin='Rrel',offset=t),
            dict(flatten_exit=-100*self.bp-c.ln(2),power_length=-self.bp*self.outer.Lrel,
                transition_offset=-self.bp*t,transition_J=-self.rate*J),dict(original_transition_kernels=kernels))

    def steep_power(self,Z,phase):
        c=self.ctx; phase=c.mpf(phase); self._phase(phase)
        data=self.data(Z); t=self.Ts*phase; left=self.Ts*(1-phase); one=IntervalTaylor.constant(c,1,5)
        theta=one*self.thetaS*c.exp(-c.mpf('1.5')*t); X=data['XS']+t
        # Algebraically [I2(left)+exp(-2left)*after_power]/2, with
        # after_power>1/2 proved above. This gives a uniform 1/4 floor.
        energy=(data['after_power']-c.mpf('.5'))*c.exp(-2*left)/2+c.mpf('.25')
        pressure=data['PS']+decay_integral(c,3,t)*(self.outer.flatten.Ev2*self.thetaS**2/2)
        rates=[-self.rate,c.mpf(0),c.mpf(0),c.mpf(0)]
        return self.packet(Z,theta,X,energy,pressure,rates,'O7 steep power',dict(origin='Rs',phase=phase,offset=t,remaining=left),
            dict(flatten_exit=-100*self.bp-c.ln(2),power_length=-self.bp*self.outer.Lrel,
                entry_offset=-self.bp-self.rate/2,steep_power_offset=-c.mpf('1.5')*t),
            dict(uniform_normalized_swirl_energy_floor=c.mpf('.25')))

    def steep_out(self,Z,t):
        c=self.ctx; t=c.mpf(t); data=self.data(Z); kernels=self.kernels(t,'out'); J=kernels['J']
        one=IntervalTaylor.constant(c,1,5)
        theta=one*self.thetaQ*c.exp(-c.mpf('1.5')*t+self.k*J)
        X=(data['XQ']+kernels['angular'])*c.exp(-self.k*J)
        energy=(data['waiting_future']*c.exp(-1-self.delta/2)+kernels['remaining_energy'])*c.exp(2*t-2*self.k*J)/2
        pressure=data['PQ']+kernels['pressure']*(self.outer.flatten.Ev2*self.thetaQ**2)
        sig=sigma_jets(c,t)
        rates=[-self.rate+self.k*sig[0]]+[self.k*sig[j]*math.factorial(j) for j in range(1,4)]
        return self.packet(Z,theta,X,energy,pressure,rates,'O7 steep exit',dict(origin='Rq',offset=t),
            dict(flatten_exit=-100*self.bp-c.ln(2),power_length=-self.bp*self.outer.Lrel,
                entry_offset=-self.bp-self.rate/2,steep_power=-c.mpf('1.5')*self.Ts,
                exit_offset=-c.mpf('1.5')*t,exit_J=self.k*J),dict(original_transition_kernels=kernels))

    @staticmethod
    def _phase(phase):
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Phase in[0,1] required')

    def waiting(self,Z,phase):
        c=self.ctx; phase=c.mpf(phase); self._phase(phase)
        data=self.data(Z); t=self.wait*phase; left=self.wait*(1-phase); one=IntervalTaylor.constant(c,1,5)
        theta=one*self.thetaT*c.exp(-self.bh*t)
        X=(data['XT']-1/self.k)*c.exp(-self.k*t)+1/self.k
        energy=(data['H']*c.exp(-self.delta*left)+decay_integral(c,self.delta,left))/2
        pressure=data['PT']+decay_integral(c,1+self.delta,t)*(self.outer.flatten.Ev2*self.thetaT**2/2)
        rates=[self.mu-self.delta/2,c.mpf(0),c.mpf(0),c.mpf(0)]
        return self.packet(Z,theta,X,energy,pressure,rates,'O7 waiting',dict(origin='Rt',phase=phase,offset=t,remaining=left),
            dict(flatten_exit=-100*self.bp-c.ln(2),power_length=-self.bp*self.outer.Lrel,
                entry_offset=-self.bp-self.rate/2,steep_power=-c.mpf('1.5')*self.Ts,
                exit_offset=-c.mpf('1.5')+self.k/2,waiting_offset=-self.bh*t),
            dict(refined_original_waiting_root=self.wait,waiting_log_one_minus_epsilon=self.logone,
                exact_collar_future_energy_in_waiting_exit_units=data['H']))

    def report(self):
        with mp.workdps(270):
            methods={'entry':self.steep_in,'power':self.steep_power,'exit':self.steep_out,'waiting':self.waiting}
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                selected_mu=self.mu,selected_delta=self.delta,actual_epsilon=self.epsilon,
                epsilon_relation='epsilon=.001*delta',original_steep_power_length=self.Ts,
                original_refined_waiting_root=self.wait,waiting_log_one_minus_epsilon=self.logone,
                preheat_tail_atoms=self.atoms,exact_Gamma_factor_positive_cap=self.delta*self.S/2,
                samples={name:[method(z,p) for z in ('-1','0','.5','1') for p in ('0','.5','1')]
                    for name,method in methods.items()},
                whole_Z={name:dict(inlet=method([-1,1],0),entire=method([-1,1],[0,1]),terminal=method([-1,1],1))
                    for name,method in methods.items()},
                entire_steep_waiting_high_mixed_derivatives_available=True,
                angular_steep_and_internal_joins_certified=False,
                full_pulse_C4_installed=True,full_outer_C4_certified=False,
                collar_and_Gamma_high_mixed_derivatives_available=False,
                physical_energy_integral_certified=False,whole_outer_cone_certified=False,
                temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepWaitingC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual O7 steep entry/power/exit and waiting: same histories, positive backward energy and mixed<=4 generated',flush=True)
    return result


if __name__=='__main__':run()
