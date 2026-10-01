"""Same-source angular candidate O.5-O.7 and continuous preheat waiting.

The angular moment X=Mtheta/(sqrt(2)*R^1.5*Utheta) is transported from
the actual Rp inlet. Pulse axial coefficients do not change this primitive.
Super-large amplitude/radius origins remain separated formal log terms.
This is angular candidate data, not completed five-moment/heat matching.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_shared_outer_buffer import SharedOuterBuffer
from lei_ren_part1_paper_shared_five_moment_repair import pack
from lei_ren_part1_paper_shared_outer_initial import stable_sigma
from lei_ren_part1_paper_interval_outer_slope_field import fraction_box
from lei_ren_part1_paper_interval_functional_defects import logjet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def unit_kernels(c,t,rate,kind,cells=256):
    """Directed J and integrating-factor primitive in a unit transition."""
    t=c.mpf(t)
    if endpoints(t)[0]<0 or endpoints(t)[1]>1:raise ValueError('Unit-transition t in [0,1] required')
    J=c.mpf(0);I=c.mpf(0)
    for i in range(cells):
        a=t*i/cells;b=t*(i+1)/cells;da=b-a
        sa=stable_sigma(c,a)[0];sb=stable_sigma(c,b)[0]
        nextJ=J+da*c.mpf([endpoints(sa)[0],endpoints(sb)[1]])
        jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
        if kind=='in':I+=(c.exp(rate*b)-c.exp(rate*a))/rate*c.exp(-rate*jc)
        elif kind=='out':I+=da*c.exp(rate*jc)
        else:raise ValueError('kind in or out required')
        J=nextJ
    if endpoints(t)==(mp.mpf(1),mp.mpf(1)):J=c.mpf('.5')
    return J,I


def collar_preheat_integral(c,k,cells=768):
    """True J in7.10: integral_0^3 exp(k t)[1-sigma+sigma*f]dt."""
    result=c.mpf(0)
    def f(t):
        y=(3-t)/2
        if endpoints(y)[1]<=0:return c.mpf(0)
        return c.exp(-1/y**2)
    for i in range(cells):
        a=c.mpf(3)*i/cells;b=c.mpf(3)*(i+1)/cells
        sa=stable_sigma(c,a)[0];sb=stable_sigma(c,b)[0]
        sig=c.mpf([endpoints(sa)[0],endpoints(sb)[1]])
        fv=c.mpf([endpoints(f(b))[0],endpoints(f(a))[1]])
        bracket=1-sig+sig*fv
        bracket=c.mpf([max(mp.mpf(0),endpoints(bracket)[0]),min(mp.mpf(1),endpoints(bracket)[1])])
        result+=(c.exp(k*b)-c.exp(k*a))/k*bracket
    return result


class SharedOuterAngularCandidate:
    def __init__(self,cells=256):
        self.buffer=SharedOuterBuffer();self.initial=self.buffer.initial;self.ctx=c=self.buffer.ctx
        self.params=self.buffer.params;self.mu=self.params.mu;self.delta=self.initial.delta;self.epsilon=self.params.epsilon
        self.rate=1-self.mu;self.restore_rate=1-self.delta/2;self.cells=cells;self.hashes=dict(self.buffer.hashes)
        for n in ('shared_outer_buffer','shared_outer_buffer_check','shared_outer_pulse_map','shared_outer_pulse_map_check'):
            name=PREFIX+n+'.json';r=json.loads((HERE/name).read_bytes())
            if r.get('actual_five_defect_family_sha256')!=self.initial.family:raise ValueError('Angular/pulse family mismatch')
            for source,digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Angular candidate source changed: '+source)
                self.hashes[source]=digest
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            inlet=self.buffer.power('0',1,cells)
            self.loguRp0=c.ln(inlet['Utheta_over_Pstar'][0])
            Xp=inlet['Mtheta_over_sqrt2_R_3half_Pstar'][0]/inlet['Utheta_over_Pstar'][0]
            eq=1/self.rate;diff=Xp-eq;absdiff=max(abs(v) for v in endpoints(diff))
            if absdiff and endpoints(c.ln(c.mpf(absdiff))-13*self.rate/self.mu)[1]>=-1000:
                raise ArithmeticError('Actual pulse angular-memory cap failed')
            self.pulse_memory_cap=c.exp(-1000)
            cap=endpoints(self.pulse_memory_cap)[1]
            self.Xv=eq+c.mpf([-cap,cap])
            self.Xp=Xp
            self.collarJ=collar_preheat_integral(c,self.restore_rate)
            terminal=self.before_waiting('0')
            Xt=terminal['X'][0];k=self.restore_rate;eq=1/k
            if endpoints(Xt-eq)[0]<=0:raise ArithmeticError('Waiting source numerator not positive')
            logone=c.mpf([endpoints(-self.epsilon/(1-self.epsilon))[0],endpoints(-self.epsilon)[1]])
            self.waiting=(c.ln(Xt-eq)+logone-self.params.log_epsilon-c.ln(eq+self.collarJ))/k
            old=self.params.waiting['root_interval']
            if not endpoints(old)[0]<=endpoints(self.waiting)[0]<=endpoints(self.waiting)[1]<=endpoints(old)[1]:
                raise ArithmeticError('Actual refined waiting root outside same-source prior enclosure')
            if endpoints(self.waiting)[0]<=0:raise ArithmeticError('Selected implicit waiting length not positive')
            self.waiting_logone=logone
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def qlog(self,Z):
        c=self.ctx;z,qi=self.initial.coordinate(Z)
        q=IntervalTaylor(c,[1+z[0]**2,2*z[0]])
        return z,logjet(q)

    def flatten(self,Z,phase):
        c=self.ctx;phase=fraction_box(c,phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Flatten phase in[0,1] required')
        t=100*phase;z,lq=self.qlog(Z);ratio=lq-c.ln(2)
        sig=stable_sigma(c,phase)[0];F=(ratio*sig).exp()
        integral=z*0
        for i in range(self.cells):
            a=t*i/self.cells;b=t*(i+1)/self.cells
            sa=stable_sigma(c,a/100)[0];sb=stable_sigma(c,b/100)[0]
            ss=c.mpf([endpoints(sa)[0],endpoints(sb)[1]])
            factor=(ratio*ss).exp()
            integral+=factor*((c.exp(self.rate*b)-c.exp(self.rate*a))/self.rate)
        X=(integral+self.Xv)*c.exp(-self.rate*t)/F
        logfactor=ratio*sig-lq+(-c.mpf('.5')-self.mu)*t
        return dict(X=X,relative_log_profile=logfactor,stage='O.5 angular flatten candidate',
                    coordinate=dict(origin='v',local_offset=t,phase=phase),flat_sigma=sig)

    def reference_buffer(self,Z,phase=1):
        c=self.ctx;phase=fraction_box(c,phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Post-flatten phase in[0,1] required')
        inlet=self.flatten(Z,1);t=-30*self.params.log_mu*phase;eq=1/self.rate
        departure=(inlet['X']-eq)*c.exp(-self.rate*t)
        X=departure+eq;logfactor=inlet['relative_log_profile']+(-c.mpf('.5')-self.mu)*t
        return dict(X=X,relative_log_profile=logfactor,stage='O.6 post-flatten angular candidate',
            coordinate=dict(origin='f',local_offset=t,phase=phase),X_departure=departure)

    def steep_in(self,Z,offset=1):
        c=self.ctx;t=fraction_box(c,offset);inlet=self.reference_buffer(Z)
        J,I=unit_kernels(c,t,self.rate,'in',self.cells)
        X=(inlet['X']+I)*c.exp(-self.rate*(t-J))
        logfactor=inlet['relative_log_profile']+(-c.mpf('.5')-self.mu)*t-self.rate*J
        return dict(X=X,relative_log_profile=logfactor,stage='O.7 steep transition in',
                    coordinate=dict(origin='rel',local_offset=t),J=J,integrating_factor_primitive=I)

    def steep_power(self,Z,phase=1):
        c=self.ctx;phase=fraction_box(c,phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Steep-power phase in[0,1] required')
        inlet=self.steep_in(Z);t=self.params.Ts*phase
        return dict(X=inlet['X']+t,relative_log_profile=inlet['relative_log_profile']-c.mpf('1.5')*t,
                    stage='O.7 steep power',coordinate=dict(origin='s',local_offset=t,phase=phase))

    def before_waiting(self,Z,offset=1):
        c=self.ctx;t=fraction_box(c,offset);inlet=self.steep_power(Z)
        J,I=unit_kernels(c,t,self.restore_rate,'out',self.cells)
        X=(inlet['X']+I)*c.exp(-self.restore_rate*J)
        logfactor=inlet['relative_log_profile']-c.mpf('1.5')*t+self.restore_rate*J
        return dict(X=X,relative_log_profile=logfactor,stage='O.7 steep transition out',
                    coordinate=dict(origin='q',local_offset=t),J=J,integrating_factor_primitive=I)

    def waiting_field(self,Z,phase):
        c=self.ctx;phase=fraction_box(c,phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Waiting phase in[0,1] required')
        inlet=self.before_waiting(Z);t=self.waiting*phase;eq=1/self.restore_rate
        departure=(inlet['X']-eq)*c.exp(-self.restore_rate*t)
        return dict(X=departure+eq,X_departure=departure,
                    relative_log_profile=inlet['relative_log_profile']+(-c.mpf('.5')-self.delta/2)*t,
                    stage='O.7 waiting candidate',coordinate=dict(origin='t',local_offset=t,phase=phase))

    def packet(self,record):
        c=self.ctx
        return dict(**record,formal_log_amplitude_origin=dict(
            log_Utheta_Rp0_over_Pstar=self.loguRp0,inverse_mu_term=dict(numerator='-13/2',denominator='mu'),finite_offset=-13),
            meaning='log(Utheta/Pstar)=formal origin plus relative log profile; sum of huge origins is not materialized',
            X_units='Mtheta/(sqrt(2)*R^1.5*Utheta)',
            angular_moment_inherited_through_axial_pulse=True,
            pulse_memory_tail_not_replaced_by_zero=True,
            angular_pressure_corrections_installed=False,actual_ap_selected=False,
            exact_heat_exterior_installed=False,whole_outer_cone_certified=False,temporal_recursion=False)

    def report(self):
        c=self.ctx
        with mp.workdps(210):
            examples=[self.flatten('.5',p) for p in ('0','.5','1')]
            examples += [self.reference_buffer('.5'),self.steep_in('.5'),self.steep_power('.5'),self.before_waiting('.5')]
            examples += [self.waiting_field('.5',p) for p in ('0','1')]
            return dict(actual_five_defect_family_sha256=self.initial.family,
                implicit_source_sha256=self.initial.datum.source_sha,datum_enclosure_sha256=self.initial.datum.datum_sha,
                angular_Xp_at_pulse_inlet=self.Xp,angular_Xv_enclosure=self.Xv,
                retained_pulse_angular_memory_positive_cap=self.pulse_memory_cap,
                directed_preheat_collar_J=self.collarJ,refined_same_source_waiting_root=self.waiting,
                prior_same_source_waiting_root=self.params.waiting['root_interval'],
                waiting_log_one_minus_epsilon=self.waiting_logone,waiting_rate=self.restore_rate,
                samples=[self.packet(r) for r in examples],
                actual_angular_candidate_through_waiting_callable=True,
                actual_continuous_preheat_waiting_root_more_tightly_enclosed=True,
                waiting_root_definition_not_changed=True,
                angular_pressure_corrections_installed=False,actual_ap_selected=False,
                complete_five_moment_outer_built=False,exact_heat_exterior_installed=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    field=SharedOuterAngularCandidate();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Same-source angular flatten/steep/waiting candidate callable; actual preheat waiting root refined',flush=True)
    return result


if __name__=='__main__':run()
