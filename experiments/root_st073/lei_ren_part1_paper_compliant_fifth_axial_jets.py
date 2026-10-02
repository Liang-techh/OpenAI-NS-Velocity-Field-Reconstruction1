"""Separate fifth-order extension of the SAME admitted leading functions.

All orders zero through four are inherited unchanged. Only the fifth
coefficient is recovered from the exact equations and true integral
derivatives. No C4 admission guard is bypassed or old receipt rewritten.
These are leading-profile derivatives, not temporal coefficient recursion.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_axial_high_jets import CompliantAxialHighJets
from lei_ren_part1_paper_compliant_angular_high_jets import heat_future_jets, log_taylor
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_compliant_outer_buffer import decay_integral
from lei_ren_part1_paper_compliant_outer_angular_repair import magnitude
from lei_ren_part1_paper_compliant_future_energy_high_jets import copy_jet
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def append_fifth(c, prior, coefficient):
    if prior.order != 4: raise ValueError('An admitted order-four prefix is required')
    return IntervalTaylor(c,list(copy_jet(c,prior).coefficients)+[coefficient])


def angular_fifth(c, prior, b1, b2, p, q, nonlinear):
    x,y=prior; J21=1+2*nonlinear*x[0]; J22=q*(1+2*nonlinear*y[0])
    det=p*J22-J21
    if endpoints(det)[0] <= 0 <= endpoints(det)[1]: raise ArithmeticError('Fifth-order angular inverse lost')
    known=nonlinear*sum((x[i]*x[5-i]+q*y[i]*y[5-i] for i in range(1,5)),c.mpf(0))
    value=b2-known
    return [append_fifth(c,x,(J22*b1-value)/det),append_fifth(c,y,(p*value-J21*b1)/det)],det


def selected_fifth(c, prior, A2, A1, A0, denominator):
    if prior.order != 4 or A1.order != 5 or A0.order != 5:
        raise ValueError('Actual selected C4 prefix and C5 quadratic data required')
    if endpoints(denominator)[0] <= 0: raise ArithmeticError('Fifth-order positive root inverse lost')
    cross=A2*sum((prior[i]*prior[5-i] for i in range(1,5)),c.mpf(0))
    mixed=sum((A1[i]*prior[5-i] for i in range(1,6)),c.mpf(0))
    return append_fifth(c,prior,-(cross+mixed+A0[5])/denominator)


def preheat_fifth(c,Z,r):
    Z=c.mpf(Z); rate=c.mpf(endpoints(r.rate))
    q=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]); logq=log_taylor(q)
    Xf=q.reciprocal()*(2*c.mpf(endpoints(r.angular.Xv))*c.exp(-100*rate))
    for i in range(r.cells):
        a=c.mpf(100)*i/r.cells; b=c.mpf(100)*(i+1)/r.cells
        sig=stable_sigma(c,c.mpf([endpoints(a)[0],endpoints(b)[1]])/100)[0]; power=1-sig
        weight=(c.exp(-rate*(100-b))-c.exp(-rate*(100-a)))/rate
        Xf+=(logq*(-power)+power*c.ln(2)).exp()*weight
    return -Xf[5]


class CompliantFifthAxialJets:
    def __init__(self):
        self.fourth=CompliantAxialHighJets(); self.ctx=self.fourth.ctx
        self.angular4=self.fourth.energy.angular
        self.hashes=dict(self.fourth.hashes); self.angular_cache={}; self.energy_cache={}; self.cache={}
        for stem,gate in (('compliant_angular_high_jets_check','angular_coefficient_C4_available'),
                          ('compliant_future_energy_high_jets_check','complete_future_corrected_energy_C4_available'),
                          ('compliant_axial_high_jets_check','actual_selected_ap_c1_c2_C4_available')):
            name=PREFIX+stem+'.json'; r=json.loads((HERE/name).read_bytes())
            if not r.get(gate) or r['actual_five_defect_family_sha256'] != self.fourth.base.future.angular.initial.family:
                raise ValueError('Fifth-order prerequisite/family mismatch: '+name)
            for source,digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Fifth-order source changed: '+source)
                self.hashes[source]=digest
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def angular(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if key in self.angular_cache: return self.angular_cache[key]
        with mp.workdps(210):
            r=self.angular4.repair; old=self.angular4.coefficients(endpoints(Z))
            box=lambda v:c.mpf(endpoints(v))
            pre=append_fifth(c,old['preheat_difference_scaled_Taylor'],preheat_fifth(c,Z,r))
            rawheat=heat_future_jets(c,Z,box(r.delta)/2,box(r.heat.epsilon),box(r.strong_S_cap),5,self.angular4.cells)
            heat={name:append_fifth(c,old['scaled_Gamma_future_defect_Taylor'][name],rawheat[name][5])
                  for name in ('theta','pressure','energy')}
            rh5=heat['theta'][5]*c.mpf([0,endpoints(r.theta_heat_over_scale_cap)[1]])
            sh5=heat['pressure'][5]*c.mpf([0,endpoints(r.pressure_heat_over_scale_cap)[1]])
            rhs=append_fifth(c,old['r_scaled_Taylor'],pre[5]+rh5)
            sh=append_fifth(c,old['s_scaled_Taylor'],sh5)
            prior=[copy_jet(c,j) for j in old['scaled_coefficient_Taylor']]
            coeffs,det=angular_fifth(c,prior,rhs[5]*(c.exp(box(r.rate))/box(r.weights['A'])),
                sh[5]*(c.exp(-3*box(r.prate))/box(r.weights['B'])),box(r.p),box(r.q),box(r.k)*box(r.scale))
            physical=[j*box(r.scale) for j in coeffs]
            bound=sum((math.factorial(5)*c.mpf(magnitude(j[5])) for j in physical),c.mpf(0))
            if endpoints(bound-box(r.mu)/100)[1] >= 0: raise ArithmeticError('Angular fifth derivative smallness failed')
            out=dict(Z=Z,scaled_coefficient_Taylor=coeffs,physical_coefficient_Taylor=physical,
                preheat_difference_scaled_Taylor=pre,scaled_Gamma_future_defect_Taylor=heat,
                r_scaled_Taylor=rhs,s_scaled_Taylor=sh,derivative_jacobian_determinant=det,
                fifth_physical_derivative_sum_bound=bound,angular_coefficient_C5_available=True,
                ordinary_Taylor_order=5,admitted_C4_prefix_preserved=True,
                exact_positive_S_not_replaced_by_cap=True,entire_Gamma_tail_included=True)
            self.angular_cache[key]=out; return out

    def future(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if key in self.energy_cache: return self.energy_cache[key]
        with mp.workdps(210):
            f=self.fourth.energy.base; high=self.angular(Z); box=lambda v:c.mpf(endpoints(v))
            q=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]); logq=log_taylor(q)-c.ln(2)
            flatten=q*0
            for i in range(f.cells):
                left=c.mpf(100)*i/f.cells; right=c.mpf(100)*(i+1)/f.cells
                t=c.mpf([endpoints(left)[0],endpoints(right)[1]])
                sig=stable_sigma(c,t/100)[0]
                flatten+=(logq*(2*sig)).exp()*((right-left)*c.exp(-2*box(f.mu)*t))
            Nf=q*q*(c.exp(-200*box(f.mu))/4)
            Nrel=q*q*(c.exp(-2*box(f.mu)*(100+box(f.Lrel)))/4)
            power=Nf*decay_integral(c,2*box(f.mu),box(f.Lrel))
            angular_change=q*0; w=f.repair.weights
            for dj,center in zip(high['physical_coefficient_Taylor'],(-3,-1)):
                angular_change+=(dj*(2*box(w['E']))+dj*dj*box(w['F']))*c.exp(-2*box(f.mu)*center)
            angular_change=angular_change*Nrel
            atoms=f.heat.preheat_atoms['energy']; eps=box(f.heat.epsilon)
            postrel=Nrel*(box(f.kernels['steep_in'])+box(f.Ns)*box(f.steep_power)
                +box(f.Nq)*box(f.kernels['steep_out'])+box(f.Nt)*box(f.waiting_energy)
                +box(f.tail_multiplier)*(1/box(f.delta)-2*eps*box(atoms['W'])+eps**2*box(atoms['W_squared'])))
            heatcap=box(f.tail_multiplier)*(box(f.delta)/2)*box(f.repair.strong_S_cap)
            heatdifference=Nrel*(high['scaled_Gamma_future_defect_Taylor']['energy']*c.mpf([0,endpoints(heatcap)[1]]))
            total=flatten+power+angular_change+postrel-heatdifference
            old=self.fourth.energy.future(endpoints(Z))
            admitted=append_fifth(c,old['complete_future_energy_Taylor'],total[5])
            out=dict(Z=Z,complete_future_energy_Taylor=admitted,
                Section7_34_weighted_future_Taylor=admitted*box(f.weighted_factor),
                fifth_pieces=dict(flatten=flatten[5],power_buffer=power[5],angular_bump_energy_change=angular_change[5],
                    steep_waiting_epsilon_and_infinite_tail=postrel[5],Gamma_deficit=heatdifference[5]),
                units=old['units'],ordinary_Taylor_order=5,complete_future_corrected_energy_C5_available=True,
                admitted_C4_prefix_preserved=True,entire_infinite_heat_tail_included=True)
            self.energy_cache[key]=out; return out

    def select(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if key in self.cache: return self.cache[key]
        with mp.workdps(210):
            old=self.fourth.select(Z); b=self.fourth.base
            incoming=old['incoming']; k=self.fourth.constants
            moments=[append_fifth(c,j,c.mpf(0)) for j in incoming['moment_Taylor']]
            energy=append_fifth(c,incoming['energy_Taylor'],6*Z*k['C_E'])
            rows=[j*factor for j,factor in zip(moments,b.incoming_factor_caps)]
            u=b.linear_inverse(rows[0]*(-b.mu),-(rows[1]-rows[0]))
            future=self.future(Z)['Section7_34_weighted_future_Taylor']
            target=future-energy*b.mu+b.base
            quad=old['quadratic_coefficients']; A2=quad['A2']
            A1=append_fifth(c,quad['A1'],sum((2*nu*uj[5]*vj for nu,uj,vj in zip(b.nu,u,b.v)),c.mpf(0)))
            A05=-target[5]+sum((nu*(uj*uj)[5] for nu,uj in zip(b.nu,u)),c.mpf(0))
            A0=append_fifth(c,quad['A0'],A05)
            ap=selected_fifth(c,old['selected_ap_Taylor'],A2,A1,A0,old['positive_root_derivative_denominator'])
            controls=[append_fifth(c,prior,uj[5]+vj*ap[5])
                      for prior,uj,vj in zip(old['selected_scaled_end_coefficient_Taylor'],u,b.v)]
            result=dict(old)
            result.update(selected_ap_Taylor=ap,selected_scaled_end_coefficient_Taylor=controls,
                incoming=dict(moment_Taylor=moments,energy_Taylor=energy,ordinary_Taylor_order=5,
                    exact_shapes=incoming['exact_shapes'],Z_independent_constant_definitions=k),
                actual_affine_incoming_Taylor=u,actual_row_normalized_incoming_Taylor=rows,
                energy_target_Taylor=target,actual_weighted_incoming_energy_Taylor=energy*b.mu,
                actual_weighted_future_energy_Taylor=future,quadratic_coefficients=dict(A2=A2,A1=A1,A0=A0),
                actual_selected_ap_c1_c2_C5_available=True,ordinary_Taylor_order=5,
                admitted_C4_prefix_preserved=True,sixth_derivative_Taylor_remainder_available=False)
            self.cache[key]=result; return result

    def report(self):
        with mp.workdps(210):
            return dict(samples=[dict(angular=self.angular(z),energy=self.future(z),selected=self.select(z))
                        for z in ('-1','0','.5','1')],
                whole_Z=dict(angular=self.angular([-1,1]),energy=self.future([-1,1]),selected=self.select([-1,1])),
                actual_five_defect_family_sha256=self.fourth.base.future.angular.initial.family,
                implicit_source_sha256=self.fourth.base.future.angular.initial.datum.source_sha,
                angular_coefficient_C5_available=True,complete_future_corrected_energy_C5_available=True,
                actual_selected_ap_c1_c2_C5_available=True,ordinary_Taylor_order=5,
                admitted_C4_prefix_preserved=True,sixth_derivative_Taylor_remainder_available=False,
                full_pulse_C4_installed=False,full_outer_C4_certified=False,whole_outer_cone_certified=False,
                temporal_recursion=False,input_hashes=self.hashes)


def run():
    r=CompliantFifthAxialJets().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(r)),indent=2)+'\n',encoding='utf-8')
    print('Separate SAME-source fifth angular/complete energy/selected ap/end coefficients generated; pulse radial C4 next',flush=True)
    return r


if __name__=='__main__': run()
