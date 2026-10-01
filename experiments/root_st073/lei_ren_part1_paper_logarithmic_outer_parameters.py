"""Directed outer parameters and segmented radii, without Decimal expansion.

Huge pulse lengths remain an exponential term, so subtracting neighboring
checkpoints never discards the 100-unit flatten or 3-unit collar. Waiting
is enclosed from the continuous source equation, with explicit integral
bounds; no float bump weights or nominal ODE factory are instantiated.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


class LogarithmicOuterParameters:
    stages=('reference','d','w','p','v','f','rel','s','q','t','tail','b')

    def __init__(self,Md,precision=120):
        self.ctx=c=MPIntervalContext();c.dps=precision
        self.precision=precision;self.Md=str(Md)
        with mp.workdps(precision+40):
            self.md=c.mpf(Md)
            if endpoints(self.md)[0]<2:raise ValueError('Md>=2 required')
            self.logPstar=c.exp(self.md)+11
            self.log_mu=c.ln(c.mpf('.001'))-4*self.logPstar
            limit=-200*c.ln(10);choice=-4*self.logPstar-30
            if endpoints(choice)[1]<endpoints(limit)[0]:self.log_delta=choice
            elif endpoints(limit)[1]<endpoints(choice)[0]:self.log_delta=limit
            else:raise ValueError('ambiguous min branch requires a tighter enclosure')
            self.log_epsilon=c.ln(c.mpf('.01'))+self.log_delta
            self.mu=c.exp(self.log_mu);self.delta=c.exp(self.log_delta)
            self.epsilon=c.exp(self.log_epsilon)
            self.Tw=-60*self.log_mu;self.Ts=4*(c.ln(2)-self.log_delta)
            self.yd=c.exp(self.md)+11
            self.waiting=self.waiting_root_enclosure()
            self.edges=[self.yd,c.mpf(1),self.Tw,None,c.mpf(100),
                -30*self.log_mu,c.mpf(1),self.Ts,c.mpf(1),self.waiting['root_interval'],c.mpf(3)]

    def waiting_root_enclosure(self):
        c=self.ctx
        a=1-self.mu;k=1-self.delta/2;ed=c.exp(a/2);er=c.exp(k/2);eq=1/k
        if endpoints(self.mu)[1]>mp.mpf('.001'):raise ValueError('mu upper bound required')
        # Before y_rel, the Z=0 raw angular moment solves X'+rX=1,
        # X(0)=5/8, r>=1-mu-32 log(2)/100>1/2. Hence 0<X<=2.
        rate_floor=1-self.mu-32*c.ln(2)/100
        if endpoints(rate_floor)[0]<=c.mpf('.5'):raise ValueError('incoming moment bound not justified')
        incoming=c.mpf([0,2])
        Lrestore=c.mpf([1,endpoints(er)[1]])
        Ldrop=c.mpf([1,endpoints(ed)[1]])
        K=c.mpf([0,endpoints((c.exp(3*k)-1)/k)[1]])
        constant=ed*(er*eq-Lrestore-self.Ts)-Ldrop
        positive=incoming-constant
        # Never form the very small exponential coefficient before taking
        # its logarithm. This is the same equation as waiting_length.py.
        log_coefficient=a/2+k/2+self.log_epsilon-c.ln(1-self.epsilon)+c.ln(eq+K)
        if endpoints(positive)[0]<=0:raise ValueError('positive waiting target required')
        root=(c.ln(positive)-log_coefficient)/k
        if endpoints(root)[0]<0:raise ValueError('nonnegative waiting root not enclosed')
        return dict(root_interval=root,incoming_X_interval=incoming,
            incoming_ODE_rate_lower=rate_floor,restore_integral=Lrestore,drop_integral=Ldrop,
            collar_deficit_integral_K=K,constant_term=constant,
            log_exponential_coefficient=log_coefficient,growth_rate=k,
            continuous_waiting_equation_has_unique_positive_root=True,
            root_encloses_all_continuous_source_integrals_under_recorded_bounds=True,
            explicit_hypotheses=['mu<=.001','Tf=100','0<=sigma_prime<=32',
                '0<=flat_edge<=1','uncorrected incoming X at Z=0 starts at 5/8'],
            exact_root_value_selected=False,angular_moment_repair_completed=False,
            heat_moment_repair_completed=False)

    def distance(self,start,end):
        """Return an affine distance; never collapse the huge pulse term."""
        c=self.ctx;i=self.stages.index(start);j=self.stages.index(end)
        sign=1 if j>=i else -1;lo=min(i,j);hi=max(i,j)
        finite=c.mpf(0);terms=[]
        for n in range(lo,hi):
            if self.edges[n] is None:
                terms.append(dict(coefficient=13*sign,log_scale=-self.log_mu,
                    scale_symbol='inverse_mu'))
            else:finite+=self.edges[n]*sign
        return dict(start=start,end=end,finite_offset=finite,exponential_terms=terms)

    def local_distance(self,start,end):
        packet=self.distance(start,end)
        if packet['exponential_terms']:
            raise ValueError('cross-pulse distances require the segmented representation')
        return packet['finite_offset']

    def report(self):
        c=self.ctx
        with mp.workdps(self.precision+40):
            hierarchy_margin=c.ln(c.mpf('.001'))+self.log_mu-self.log_delta
            if endpoints(hierarchy_margin)[0]<=0:raise ValueError('delta<c_delta*mu failed')
            for start,end,length in [('d','w',1),('v','f',100),('tail','b',3)]:
                actual=self.local_distance(start,end)
                if endpoints(actual)!=(mp.mpf(length),mp.mpf(length)):
                    raise ValueError('exact short stage was lost')
            return dict(Md=self.Md,logPstar=self.logPstar,log_mu=self.log_mu,
                log_delta=self.log_delta,log_epsilon=self.log_epsilon,
                actual_epsilon_relation='epsilon=.01*delta',Tw=self.Tw,Ts=self.Ts,
                waiting_root=self.waiting,delta_hierarchy_log_margin=hierarchy_margin,
                coordinates={name:self.distance('reference',name) for name in self.stages},
                pulse_length_times_mu=13,short_stage_lengths_preserved_exactly=True,
                transcendental_parameter_rounding_enclosed=True,
                precision_independent_of_inverse_mu_decimal_digit_count=True,
                source_fourteen_pressure_atoms_generated=False,new_core_generated=False,
                background_admissibility_certified=False,temporal_recursion=False)


def run():
    rows=[]
    with mp.workdps(180):
        old=json.loads((HERE/'lei_ren_part1_paper_large_Md_screen.json').read_bytes())
        for md in ('4','5','6','40','128'):
            p=LogarithmicOuterParameters(md);row=p.report()
            if md in ('4','5','6'):
                trial=next(t for t in old['trials'] if t['Md']==md)
                value=mp.mpf(trial['new_schedule_inputs']['waiting_length'])
                lo,hi=endpoints(row['waiting_root']['root_interval'])
                if not lo<=value<=hi:raise ValueError('existing nominal root outside analytic enclosure')
                row['previous_nominal_waiting_root_contained']=True
            rows.append(row)
            print('Logarithmic Md',md,'waiting enclosed; short stages exact',flush=True)
    names=(Path(__file__).name,'lei_ren_part1_paper_waiting_length.py',
        'lei_ren_part1_paper_large_Md_screen.json','lei_ren_part1_paper_outer.py',
        'lei_ren_part1_paper_interval_long_reshape_field.py')
    result=dict(parameter_families=rows,input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':run()
