"""Complete corrected swirl energy from Rv through the infinite Gamma tail.

All pieces use Rv*Utheta(Rv,Z)^2 units. The actual angular corrections,
epsilon atoms and formal nonzero heat deficit are retained. This supplies
Section 7.34's fixed future energy, not a full physical energy certificate.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_outer_angular_repair import CompliantAngularRepair, magnitude
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_compliant_outer_buffer import decay_integral
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_functional_defects import logjet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def steep_energy_kernels(c,mu,delta,cells=512):
    """Closed cells enclose both unit transitions' actual energy kernels."""
    J=c.mpf(0);inside=c.mpf(0);outside=c.mpf(0)
    for i in range(cells):
        a=c.mpf(i)/cells;b=c.mpf(i+1)/cells
        t=c.mpf([endpoints(a)[0],endpoints(b)[1]])
        sig=stable_sigma(c,t)[0];nextJ=J+(b-a)*sig
        jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
        inside+=(b-a)*c.exp(-2*mu*t-2*(1-mu)*jc)
        outside+=(b-a)*c.exp(-2*t+2*(1-delta/2)*jc)
        J=nextJ
    return dict(steep_in=inside,steep_out=outside,endpoint_J_exact=c.mpf('.5'),cells=cells)


class CompliantFutureSwirlEnergy:
    def __init__(self,cells=512):
        self.repair=CompliantAngularRepair();self.ctx=c=self.repair.ctx
        self.angular=self.repair.angular;self.heat=self.repair.heat;self.params=self.repair.params
        self.mu=self.repair.mu;self.delta=self.repair.delta;self.cells=cells
        self.hashes=dict(self.repair.hashes)
        for stem in ('compliant_outer_angular_repair','compliant_outer_angular_repair_check'):
            name=PREFIX+stem+'.json';record=json.loads((HERE/name).read_bytes())
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Future energy repair dependency changed: '+source)
                self.hashes[source]=digest
            if record['implicit_source_sha256']!=self.angular.initial.datum.source_sha:
                raise ValueError('Future energy source mismatch')
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            self.Lrel=-30*self.params.log_mu
            self.kernels=steep_energy_kernels(c,self.mu,self.delta,cells)
            self.Ns=c.exp(-1-self.mu)
            self.Nq=self.Ns*c.exp(-2*self.params.Ts)
            self.Nt=self.Nq*c.exp(-1-self.delta/2)
            self.Ntail=self.Nt*c.exp(-self.delta*self.angular.waiting)
            self.steep_power=decay_integral(c,c.mpf(2),self.params.Ts)
            self.waiting_energy=decay_integral(c,self.delta,self.angular.waiting)
            self.tail_multiplier=self.Ntail*c.exp(-2*self.angular.waiting_logone)
            self.pulse_energy_attenuation=c.exp(-26)
            self.weighted_factor=self.mu*self.pulse_energy_attenuation/2
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def future(self,Z):
        c=self.ctx;Z=c.mpf(Z)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Z in [-1,1] required')
        q=IntervalTaylor(c,[1+Z**2,2*Z]);logq=logjet(q)-c.ln(2)
        flatten=q*0
        for i in range(self.cells):
            left=c.mpf(100)*i/self.cells;right=c.mpf(100)*(i+1)/self.cells
            t=c.mpf([endpoints(left)[0],endpoints(right)[1]])
            sig=stable_sigma(c,t/100)[0]
            flatten+=(logq*(2*sig)).exp()*((right-left)*c.exp(-2*self.mu*t))
        Nf=q*q*(c.exp(-200*self.mu)/4)
        Nrel=q*q*(c.exp(-2*self.mu*(100+self.Lrel))/4)
        power=Nf*decay_integral(c,2*self.mu,self.Lrel)
        d=self.repair.coefficients(Z)['physical_coefficient_Taylor']
        w=self.repair.weights;angular_change=q*0
        for dj,center in zip(d,(-3,-1)):
            angular_change+=(dj*(2*w['E'])+dj*dj*w['F'])*c.exp(-2*self.mu*center)
        angular_change=angular_change*Nrel
        atoms=self.heat.preheat_atoms['energy'];eps=self.heat.epsilon
        baseline=1/self.delta
        epsilon_atom=-2*eps*atoms['W'];epsilon_squared_atom=eps**2*atoms['W_squared']
        preheat_tail=baseline+epsilon_atom+epsilon_squared_atom
        postrel=Nrel*(self.kernels['steep_in']+self.Ns*self.steep_power
            +self.Nq*self.kernels['steep_out']+self.Nt*self.waiting_energy
            +self.tail_multiplier*preheat_tail)
        heat=self.heat.future_defects(Z,cells=256)['swirl_energy_heat_difference_scaled']
        heat_cap=self.tail_multiplier*(self.delta/2)*self.repair.strong_S_cap
        heat_jet=IntervalTaylor(c,[c.mpf([0,max(mp.mpf(0),endpoints(heat[0]*heat_cap)[1])]),
            c.mpf([-magnitude(heat[1]*heat_cap),magnitude(heat[1]*heat_cap)])])
        if endpoints(Z)==(mp.mpf(0),mp.mpf(0)):heat_jet=IntervalTaylor(c,[heat_jet[0],0])
        if endpoints(Z) in ((mp.mpf(-1),mp.mpf(-1)),(mp.mpf(1),mp.mpf(1))):
            heat_jet=IntervalTaylor(c,[0,heat_jet[1]])
        heat_difference=Nrel*heat_jet
        total=flatten+power+angular_change+postrel-heat_difference
        if total.order!=1 or endpoints(total[0])[0]<=0:
            raise ArithmeticError('Complete future energy positivity/C1 lost')
        weighted=total*self.weighted_factor
        return dict(Z=Z,complete_future_energy_Taylor=total,
            units='integral_Rv^infinity Utheta_corrected^2 dR /(Rv*Utheta(Rv,Z)^2)',
            Section7_34_weighted_future_Taylor=weighted,
            Section7_34_units='mu/(2*Rp*Utheta(Rp,Z)^2) integral_Rv^infinity Utheta_corrected^2 dR',
            pieces=dict(flatten=flatten,power_buffer=power,angular_bump_energy_change=angular_change,
                steep_waiting_preheat_collar_and_infinite_tail=postrel,
                positive_Gamma_heat_energy_deficit_enclosure=heat_difference),
            exact_heat_energy_deficit_definition='Nrel*Ntail/(1-epsilon)^2 * (a*S) * Energy_hat(Z), S=1/Rtail; strong_S_cap is only an upper bound',
            separate_preheat_tail_atoms=dict(baseline=baseline,epsilon_atom=epsilon_atom,epsilon_squared_atom=epsilon_squared_atom),
            Nf=Nf,Nrel=Nrel,complete_future_corrected_swirl_energy_defined=True,
            infinite_heat_tail_included=True,signed_angular_energy_change_retained=True,
            actual_heat_factor_and_formal_S_not_zeroed=True,actual_ap_selected=False,
            full_physical_kinetic_energy_certified=False,temporal_recursion=False)

    def report(self):
        with mp.workdps(210):
            return dict(actual_five_defect_family_sha256=self.angular.initial.family,
                implicit_source_sha256=self.angular.initial.datum.source_sha,
                datum_enclosure_sha256=self.angular.initial.datum.datum_sha,
                samples=[self.future(z) for z in ('-1','0','.5','1')],whole_Z_C1=self.future([-1,1]),
                postrel_energy_kernels=self.kernels,relative_energy_prefactors=dict(Ns=self.Ns,Nq=self.Nq,Nt=self.Nt,Ntail=self.Ntail),
                postflatten_length=self.Lrel,steep_power_energy_integral=self.steep_power,
                waiting_energy_integral=self.waiting_energy,tail_multiplier=self.tail_multiplier,
                pulse_energy_attenuation=self.pulse_energy_attenuation,Section7_34_weighted_factor=self.weighted_factor,
                complete_future_corrected_swirl_energy_available=True,
                actual_ap_selected=False,full_outer_five_moment_match=False,whole_outer_cone_certified=False,
                full_physical_kinetic_energy_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    energy=CompliantFutureSwirlEnergy();result=energy.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Complete corrected future swirl energy: flatten, power, actual angular bumps, steep/waiting, epsilon collar and infinite Gamma tail C1 enclosures generated',flush=True)
    return result


if __name__=='__main__':run()
