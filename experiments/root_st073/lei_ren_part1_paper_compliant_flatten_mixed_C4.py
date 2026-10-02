"""Original 100-unit O.5 flatten, axial C5 and mixed profile order<=4.

Uses the same complete future energy and canonical original pressure.
Flat original sigma supplies y derivatives through4. The exact formal Ev0
scale stays positive; numerical caps only enclose it. There is no temporal
coefficient recursion or physical energy/stress-cone claim in this layer.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_power_inlet_C4 import CompliantPowerInletC4
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_angular_high_jets import log_taylor
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def flatten_mixed(c,mu,rho,sigma_derivatives,theta,X,energy,pressure,Ev2):
    """Ordinary y derivatives from variable-rate logarithmic/primitive ODEs.

    sigma_derivatives are true d_y^j sigma(y/100), not Taylor coefficients.
    Every input axial jet is order5. y rows are truncated only at output.
    """
    if len(sigma_derivatives)!=5 or any(v.order!=5 for v in (rho,theta,X,energy,pressure)):
        raise ValueError('Axial5 and original sigma y derivatives0..4 required')
    one=IntervalTaylor.constant(c,1,5); bp=c.mpf('.5')+mu
    rates=[rho*sigma_derivatives[j+1] for j in range(4)]
    lograte=list(rates); lograte[0]=lograte[0]-bp
    angularrate=list(rates); angularrate[0]=angularrate[0]+1-mu
    energyrate=[v*(-2) for v in rates]; energyrate[0]=energyrate[0]+2*mu
    theta_rows=[theta]; Xrows=[X]; erows=[energy]; prows=[pressure]
    for k in range(4):
        th=theta*0; xx=one if k==0 else one*0; ee=-one/2 if k==0 else one*0
        square=theta*0
        for j in range(k+1):
            th+=lograte[j]*theta_rows[k-j]*math.comb(k,j)
            xx-=angularrate[j]*Xrows[k-j]*math.comb(k,j)
            ee+=energyrate[j]*erows[k-j]*math.comb(k,j)
            square+=theta_rows[j]*theta_rows[k-j]*math.comb(k,j)
        theta_rows.append(th); Xrows.append(xx); erows.append(ee); prows.append(square*(Ev2/2))
    zero=theta*0
    rows={UZ:[zero for _ in range(5)],UT:theta_rows,UR:[zero for _ in range(5)],P:prows}
    grid={label:{'y'+str(k)+'_Z'+str(n):jet[n]*math.factorial(n)
        for k,jet in enumerate(values) for n in range(5-k)} for label,values in rows.items()}
    return dict(physical_mixed_derivatives_total_order_le4=grid,
        physical_velocity_and_pressure_y_derivative_Taylor={label:[v.truncate(4-k) for k,v in enumerate(values)] for label,values in rows.items()},
        angular_y_derivative_Taylor=[v.truncate(5-k) for k,v in enumerate(Xrows)],
        energy_y_derivative_Taylor=[v.truncate(5-k) for k,v in enumerate(erows)],
        log_theta_y_rate_Taylor=lograte,angular_ODE_rate_Taylor=angularrate,energy_ODE_rate_Taylor=energyrate)


class CompliantFlattenMixedC4:
    def __init__(self,cells=128):
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive directed integration cell count required')
        self.inlet=CompliantPowerInletC4(); self.ctx=c=self.inlet.ctx; self.cells=cells
        self.mu=self.inlet.mu; self.delta=self.inlet.delta; self.rate=1-self.mu
        self.bp=c.mpf('.5')+self.mu; self.prate=1+2*self.mu
        self.family=self.inlet.family; self.source=self.inlet.source; self.hashes=dict(self.inlet.hashes)
        name=PREFIX+'compliant_power_inlet_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['two_sided_O3_pulse_join_certified']:raise ValueError('O3/pulse inlet join required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten source changed: '+source)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.fifth=json.loads((HERE/(PREFIX+'compliant_fifth_axial_jets.json')).read_bytes())
        pulse=json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())
        self.Xv=read_interval(c,pulse['whole_Z_terminal']['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][0])
        self.U=self.inlet.constants['U']
        self.logEv2_parts=dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))
        logcap=2*c.ln(self.mu)-1000; exactlog=sum(self.logEv2_parts.values(),c.mpf(0))
        self.log_pressure_decay_parts=dict(inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))
        self.log_decay=sum(self.log_pressure_decay_parts.values(),c.mpf(0))
        if endpoints(logcap-exactlog)[0]<=0 or endpoints(logcap-self.log_decay)[0]<=0:
            raise ArithmeticError('Exact positive Ev0/pressure-decay scale cap not proved')
        self.Ev2=c.mpf([0,endpoints(c.exp(logcap))[1]])
        self.pressure_decay=c.mpf([0,endpoints(c.exp(logcap))[1]])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def future_energy(self,Z):
        c=self.ctx; Z=c.mpf(Z)
        # Each receipt coefficient is a bound on derivatives of the same
        # actual complete future integral. Global boxes are not evaluated
        # as a fitted polynomial when the caller requests another Z.
        packet=self.fifth['whole_Z']['energy']
        for sample in self.fifth['samples']:
            if endpoints(read_interval(c,sample['energy']['Z']))==endpoints(Z):packet=sample['energy']; break
        return IntervalTaylor(c,[read_interval(c,v)/2 for v in packet['complete_future_energy_Taylor']['coefficients']])

    def flatten(self,Z,t):
        c=self.ctx; Z=c.mpf(Z); t=c.mpf(t)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(t)[0]<0 or endpoints(t)[1]>100:
            raise ValueError('Original flatten domain Z in[-1,1], t in[0,100] required')
        data=self.inlet.incoming(Z); q=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0])
        rho=log_taylor(q)-c.ln(2); sig=sigma_jets(c,t/100)
        sj=[sig[k]*math.factorial(k)/100**k for k in range(5)]
        F=(rho*sj[0]).exp(); ev=self.future_energy(Z)
        Xint=q*0; Eint=q*0; Pint=q*0
        if endpoints(t)[1]!=0:
            for i in range(self.cells):
                a=t*i/self.cells; b=t*(i+1)/self.cells
                v=c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(100),endpoints(b)[1])])
                Fc=(rho*sigma_jets(c,v/100)[0]).exp()
                length=t/self.cells
                # Exact exponential weights; interval Fc encloses every
                # source derivative on the cell, including support endpoints.
                Xint+=Fc*(c.exp(-self.rate*(t*(self.cells-i-1)/self.cells))*decay_integral(c,self.rate,length))
                Eint+=Fc*Fc*(c.exp(-2*self.mu*a)*decay_integral(c,2*self.mu,length))
                Pint+=Fc*Fc/(q*q)*(c.exp(-self.prate*a)*decay_integral(c,self.prate,length)/2)
        X=(Xint+self.Xv*c.exp(-self.rate*t))/F
        energy=(ev-Eint/2)*c.exp(2*self.mu*t)/(F*F)
        Mp_v=data['Mp']+data['u']*data['u']*((1-self.pressure_decay)/(2*self.prate))
        Mp=Mp_v+Pint*self.Ev2; pressure=Mp+data['P0']
        theta=F/q*c.exp(-self.bp*t)
        if endpoints(t)==(mp.mpf(100),mp.mpf(100)):
            # Original sigma is exactly1, so F/q=1/2 identically in Z.
            theta=IntervalTaylor.constant(c,c.exp(-100*self.bp)/2,5)
        mixed=flatten_mixed(c,self.mu,rho,sj,theta,X,energy,pressure,self.Ev2)
        if endpoints(energy[0])[0]<=0:raise ArithmeticError('Complete positive future energy lost through flatten')
        mixed.update(Z=Z,t=t,F_Taylor=F,sigma_y_derivatives=sj,
            theta_over_Ev0_Taylor=theta,angular_Taylor=X,energy_Taylor=energy,
            selected_terminal_energy_Taylor=ev,
            pressure=dict(P0_over_Pstar_squared=data['P0'],Mp_over_Pstar_squared=Mp,P_over_Pstar_squared=pressure),
            remaining_swirl_energy_in_Rv_Ev0_squared=(ev-Eint/2)*2/(q*q),
            future_energy_input_is_derivative_enclosure_not_a_polynomial_fit=True,
            exact_positive_Ev0Pstar2_log_parts=self.logEv2_parts,positive_Ev0Pstar2_cap=self.Ev2,
            exact_positive_terminal_pressure_decay_log_parts=self.log_pressure_decay_parts,
            positive_terminal_pressure_decay_cap=self.pressure_decay,
            terminal_pressure_is_an_enclosure_of_the_original_exact_source=True,
            exact_terminal_Mp_definition='Pin/q^2+(U/q)^2*(1-exp(-13/mu-26))/(2*(1+2mu)); no numerical cap replaces the defining function',
            actual_terminal_zero_linear_and_radial_histories_inherited=True,
            original_pressure_datum_and_complete_positive_future_energy_preserved=True,
            independently_evaluated_O5_flatten_high_mixed_derivatives=True,
            component_units=dict(velocity='U/Ev0; exact formal Ev0/Pstar kept separately',pressure='P/Pstar^2'),
            full_pulse_C4_installed=False,full_outer_C4_certified=False,
            physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False)
        return mixed

    def report(self):
        with mp.workdps(270):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                selected_mu=self.mu,selected_delta=self.delta,Ev0_squared_scale_log_parts=self.logEv2_parts,
                terminal_pressure_decay_log_parts=self.log_pressure_decay_parts,
                samples=[self.flatten(z,t) for z in ('-1','0','.5','1') for t in ('0','50','100')],
                whole_Z_inlet=self.flatten([-1,1],0),whole_Z_flatten_box=self.flatten([-1,1],[0,100]),
                whole_Z_exit=self.flatten([-1,1],100),
                independently_evaluated_O5_flatten_high_mixed_derivatives=True,
                entire_original_100_unit_flatten_mixed_C4_available=True,
                following_power_high_mixed_derivatives_available=False,
                two_sided_pulse_O5_join_certified=False,full_pulse_C4_installed=False,full_outer_C4_certified=False,
                physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
                input_hashes=self.hashes)


def run():
    result=CompliantFlattenMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Original O5 entire flatten: canonical swirl, angular/positive energy and original pressure mixed derivatives through4 generated',flush=True)
    return result


if __name__=='__main__':run()
