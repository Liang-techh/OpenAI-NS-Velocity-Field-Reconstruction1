"""Conditioned C0 Section11 phase inverse and primitives on native boxes.

The two monotone angle coordinates keep narrow Poisson peaks executable.
Positive rho/s and tiny q remain log-factored. An interval returned by the
inverse encloses every source root in its input box; no midpoint defines a
field. Phase is an explicit free candidate parameter here, not N*log(R/r_-).
"""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_correlated_shear_q as current

prior=current.prior;packets=current.packets;native=current.native
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;ZERO=(0,0)
NAME=PREFIX+'current_native_conditioned_phase.json'
RECEIPT=PREFIX+'current_native_conditioned_phase_check.json'
GATE='current_native_conditioned_C0_phase_inverse_and_primitives_executed'


def clipped(c,value,lower,upper):
    lo,hi=ep(value);lo=max(lo,ep(c.mpf(lower))[0]);hi=min(hi,ep(c.mpf(upper))[1])
    if lo>hi:raise ArithmeticError('Analytic range and directed enclosure disagree')
    return c.mpf((lo,hi))


def bounded_value(value):
    """Materialize bounded values/tails only; never exponentiate a huge log."""
    if value.zero:return value.ctx.mpf(0)
    if ep(value.scale.evaluate())[1]<-1000:
        return value.coefficient*value.bounded_exp(value.scale.evaluate())
    return value.finite_interval()


def absolute(value):
    lo,hi=ep(value.coefficient)
    coeff=value.ctx.mpf((0 if lo<=0<=hi else min(abs(lo),abs(hi)),max(abs(lo),abs(hi))))
    return prior.ScaledEnclosure(value.scale,coeff,value.ledger)


def transformed_identity():
    q,t0,s,r,E,psi=sy.symbols('q t0 s r E psi',nonzero=True)
    W1=(E-psi)/(2*r)
    W2=((2/s-3)*E+psi+2*r*sy.sin(E)/s)/(4*r*r)
    original=t0*t0*psi+4*t0*q*sy.sqrt(s)*W1+4*q*q*s*W2
    conditioned=t0*t0*psi+2*t0*q*sy.sqrt(s)/r*(E-psi)+q*q/r**2*((2-3*s)*E+s*psi+2*r*sy.sin(E))
    if sy.simplify(original-conditioned)!=0:raise ArithmeticError('Original W2 transformation failed')
    return dict(passed=True,original_T2_equivalence=True,s_factor_retained=True,
        source='current_generic_shear_loop._w_integrals/integrals; paper11.10/11.13',
        monotonicity='Phi_psi=(1+t^2)/(2*pi*(1+t0^2+2*q^2))>0; Mobius E_psi=s/D>0')


class ConditionedPhase:
    def __init__(self,query,dstar_log):
        self.source=query;self.q=query['q'];self.roots=query['roots'];self.c=c=self.q.ctx
        self.point=mp.mp.clone();self.point.dps=c.dps+20
        self.a=self.roots['a'][ZERO];self.t0=self.roots['t0'][ZERO];self.E=self.roots['E'][ZERO]
        self.scalar=self.q.scalar;self.dstar=prior.ScaledEnclosure(prior.FormalScale(self.q.scale.bases,offset=dstar_log),1,self.q.ledger)
        self.u=(self.roots['p2'][ZERO]*self.q).positive_divide(self.dstar,dstar_log)
        self.nu=self.scalar(1)+current.square(self.t0)+current.square(self.q)*2
        self.flat=self.q.zero;self.geometry=None
        if self.flat:self.geometry='flat';return
        absu=absolute(self.u);lo,hi=ep(absu.coefficient)
        logupper=ep(absu.scale.evaluate()+c.ln(c.mpf(hi)))[1] if hi else -mp.inf
        loglower=ep(absu.scale.evaluate()+c.ln(c.mpf(lo)))[0] if lo else -mp.inf
        if logupper<=ep(c.ln(c.mpf('.25')))[0]:
            self.geometry='small_r_series';uf=bounded_value(self.u);h=c.sqrt(1+uf**2)
            self.r=uf/h;self.s=1/(1+uf**2);self.hinv=self.scalar(1/h)
            self.rho=self.scalar(1/(h*(h+abs(uf))))
        elif lo>0 and loglower>=ep(c.ln(c.mpf('.125')))[1]:
            self.geometry='signed_Mobius';self.sign=1 if ep(self.u.coefficient)[0]>0 else -1
            # Rebase to a bounded coefficient. This endpoint chooses arithmetic
            # coordinates only; it is not a source value or source selection.
            magnitude=hi;L=absu.scale+prior.FormalScale(absu.scale.bases,offset=c.ln(c.mpf(magnitude)))
            coeff=absu.coefficient/c.mpf(magnitude)
            z=self.q.bounded_exp((-L-L).evaluate());H=c.sqrt(coeff**2+z)
            if ep(H)[0]<=0:raise ArithmeticError('Large-u positive coefficient lost')
            self.hinv=prior.ScaledEnclosure(-L,1/H,self.q.ledger)
            self.rho=prior.ScaledEnclosure(-L-L,1/(H*(H+coeff)),self.q.ledger)
            self.s_source=prior.ScaledEnclosure(-L-L,1/(coeff**2+z),self.q.ledger)
            rho=clipped(c,bounded_value(self.rho),0,1)
            ar=clipped(c,1-rho, c.mpf('.125')/c.sqrt(1+c.mpf('.125')**2),1)
            self.r=ar*self.sign;self.s=clipped(c,bounded_value(self.s_source),0,1)
        else:
            self.geometry='requires_signed_source_refinement';return
        self.hinv_finite=clipped(c,bounded_value(self.hinv),0,1)
        self.n0=self.normalized(self.scalar(1),1,True)
        self.nt=self.normalized(current.square(self.t0),1,True)
        self.nq=self.normalized(current.square(self.q),c.mpf('.5'),True)
        self.ntq=self.normalized(self.t0*self.q,1/(2*c.sqrt(2)),False)

    def normalized(self,value,cap,positive):
        ratio=value.positive_divide(self.nu,0)
        try:finite=bounded_value(ratio)
        except ArithmeticError:finite=self.c.mpf((0 if positive else -ep(cap)[1],ep(cap)[1]))
        return clipped(self.c,finite,0 if positive else -cap,cap)

    def geometry_record(self):
        result=dict(branch=self.geometry,original_u=self.u.record(),q=self.q.record(),
            materialized_huge_exponentials=False,flat=self.flat)
        if self.geometry not in ('flat','requires_signed_source_refinement'):
            result.update(r_interval=self.r,rho=self.rho.record(),rho_positive_source_retained=True,
                one_minus_r_squared=self.s,s_computed_without_subtracting_r_squared=True)
        return result

    def angle_endpoint(self,x,inverse=False):
        c=self.c
        if x in (0,c.mpf('.5'),1):return c.mpf(x)
        reflect=x>c.mpf('.5');argument=c.mpf(x)
        if reflect:argument=1-argument
        rho=clipped(c,bounded_value(self.rho),0,1)
        signed=self.sign*(-1 if inverse else 1)
        plus,minus=(2-rho,rho) if signed>0 else (rho,2-rho)
        angle=clipped(c,c.atan2(plus*c.sin(c.pi*argument),minus*c.cos(c.pi*argument))/c.pi,0,1)
        return 1-angle if reflect else angle

    def angle_map(self,x,inverse=False):
        lo,hi=ep(self.c.mpf(x))
        return self.c.mpf((ep(self.angle_endpoint(lo,inverse))[0],ep(self.angle_endpoint(hi,inverse))[1]))

    def small_integrals(self,x):
        c=self.c;psi=2*c.pi*x;r=self.r;s=self.s
        R=max(abs(v) for v in ep(r));M=48
        if not R:return c.sin(psi),psi/2+c.sin(2*psi)/4
        W1=c.mpf(0);W2=psi/(2*s)
        for k in range(1,M+1):
            sine=c.sin(k*psi)/k
            W1+=r**(k-1)*sine
            coeff=r**k/s
            if k>=2:coeff+=(k-1)*r**(k-2)/2
            W2+=coeff*sine
        rc=c.mpf(R);smin=ep(s)[0]
        tail1=ep(rc**M/((M+1)*(1-rc)))[1]
        tail2=ep(rc**(M+1)/(c.mpf(smin)*(M+1)*(1-rc))+rc**(M-1)/(2*(1-rc)))[1]
        return W1+c.mpf((-tail1,tail1)),W2+c.mpf((-tail2,tail2))

    def angles(self,x,chart):
        if self.geometry=='signed_Mobius':
            return (self.angle_map(x,True),x) if chart=='E' else (x,self.angle_map(x))
        return x,None

    def phase(self,x,chart='psi'):
        c=self.c;x=c.mpf(x);lo,hi=ep(x)
        if lo==hi and lo in (0,mp.mpf('.5'),1):return x
        if self.flat:return x
        if self.geometry=='requires_signed_source_refinement':return c.mpf((0,1))
        psi,E=self.angles(x,chart)
        if self.geometry=='small_r_series':
            W1,W2=self.small_integrals(psi)
            value=(self.n0+self.nt)*psi+2*self.ntq*self.hinv_finite*W1/c.pi+2*self.nq*self.s*W2/c.pi
        else:
            value=(self.n0+self.nt)*psi+2*self.ntq*self.hinv_finite/self.r*(E-psi)
            value+=self.nq/self.r**2*((2-3*self.s)*E+self.s*psi+self.r/c.pi*c.sin(2*c.pi*E))
        return clipped(c,value,0,1)

    def inverse_bracket(self,phase,chart,bits):
        c=self.c;target=c.mpf(phase);pl,ph=ep(target)
        lo,hi=self.point.mpf(0),self.point.mpf(1)
        for _ in range(bits):
            mid=(lo+hi)/2
            if ep(self.phase(c.mpf(mid),chart))[1]<pl:lo=mid
            else:hi=mid
        lower=lo;lo,hi=self.point.mpf(0),self.point.mpf(1)
        for _ in range(bits):
            mid=(lo+hi)/2
            if ep(self.phase(c.mpf(mid),chart))[0]>ph:hi=mid
            else:lo=mid
        bracket=c.mpf((lower,hi));image=self.phase(bracket,chart)
        return dict(chart=chart,coordinate_interval=bracket,phase_image=image,
            phase_width=c.mpf(ep(image)[1])-c.mpf(ep(image)[0]),bracket_proof='directed endpoint inequalities and exact strict source monotonicity')

    def primitives(self,coordinate,chart='psi'):
        c=self.c;x=c.mpf(coordinate);lo,hi=ep(x)
        if self.geometry=='requires_signed_source_refinement':raise ValueError('Refine original signed source before evaluating primitives')
        if self.flat or (lo==hi and lo in (0,mp.mpf('.5'),1)):
            return dict(A=self.scalar(0),B_over_Pstar=self.scalar(0))
        psi,E=self.angles(x,chart);q=self.q;t0=self.t0;a=self.a
        if self.geometry=='small_r_series':
            W1,W2=self.small_integrals(psi)
            numerator=t0*q*self.hinv*(4*W1)+current.square(q)*(4*self.s*W2-4*c.pi*psi)
            A=(a*numerator).positive_divide(self.nu,0)*(1/(4*c.pi))
            B=self.E*(t0*A-a*q*self.hinv*(W1/(2*c.pi)))
        else:
            difference=E-psi
            numerator=t0*q*self.hinv*(2*difference/self.r)
            numerator+=current.square(q)*(((2-3*self.s)*difference+self.r/c.pi*c.sin(2*c.pi*E))/self.r**2)
            A=(a*numerator).positive_divide(self.nu,0)*c.mpf('.5')
            B=self.E*(t0*A-a*q*self.hinv*(difference/(2*self.r)))
        return dict(A=A,B_over_Pstar=B)

    def evaluate(self,phase,bits=40):
        c=self.c;phase=c.mpf(phase);lo,hi=ep(phase)
        if type(bits) is not int or not 4<=bits<self.point.prec-8:raise ValueError('Bisection bits must fit the source arithmetic context')
        if not 0<=lo<=hi<=1:raise ValueError('One explicit closed-period phase box required; no interval floor')
        if self.geometry=='requires_signed_source_refinement':
            return dict(status=self.geometry,phase=phase,geometry=self.geometry_record(),inverse_installed_on_this_box=False)
        if lo==hi and lo in (0,mp.mpf('.5'),1):
            best=dict(chart='psi',coordinate_interval=phase,phase_image=phase,phase_width=c.mpf(0),
                bracket_proof='exact periodic/half-period symmetry')
        elif self.flat:
            best=dict(chart='psi',coordinate_interval=phase,phase_image=phase,phase_width=c.mpf(hi)-c.mpf(lo),bracket_proof='exact flat Phi=psi/(2pi)')
        else:
            candidates=[self.inverse_bracket(phase,chart,bits) for chart in (('psi','E') if self.geometry=='signed_Mobius' else ('psi',))]
            best=min(candidates,key=lambda x:ep(x['phase_width'])[1])
        psi,E=self.angles(best['coordinate_interval'],best['chart']);values=self.primitives(best['coordinate_interval'],best['chart'])
        return dict(status='enclosed',phase=phase,selected_inverse=best,psi_fraction_interval=psi,
            Mobius_fraction_interval=E,geometry=self.geometry_record(),
            primitives={k:v.record() for k,v in values.items()},inverse_installed_on_this_box=True,
            free_phase_parameter_not_spatial_phase=True,original_common_N_and_radius_phase_bound=False,
            slow_derivatives_or_moment_integrals_installed=False,source_midpoint_used_as_field_value=False)


class NativeConditionedPhase:
    def __init__(self,owner):
        if type(owner) is not current.NativeCorrelatedShearQ:raise ValueError('Same live original q backend required')
        admitted=json.loads((HERE/current.RECEIPT).read_bytes())
        if not admitted['all_passed'] or admitted['source_family']!=owner.family:raise ValueError('Accepted same q source family required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.service.bind_hashes(admitted['input_hashes']);self.service.bind_hashes({current.RECEIPT:sha(current.RECEIPT),Path(__file__).name:sha(Path(__file__).name),PREFIX+'current_generic_shear_loop.py':sha(PREFIX+'current_generic_shear_loop.py')})
        self.dstar_log=packets.interval(self.ctx,owner.owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
    def query(self,chart,Z,coordinate):
        source=self.owner.query(chart,Z,coordinate)
        return source,ConditionedPhase(source,self.dstar_log)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeConditionedPhase(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))
        records={}
        for chart,coordinate in (('inner_reference','.1337'),('O2_slope','.1337'),('O2_buffer','5.337'),('O3_slope_mu','.537'),('O3_power',None)):
            if coordinate is None:
                old=json.loads((HERE/current.NAME).read_bytes())['actual_correlated_shear_and_q_records'][chart]['source_provenance']['coordinate_box']
                coordinate=packets.interval(owner.ctx,old)
            source,loop=owner.query(chart,('.5','.5'),coordinate)
            phases={phase:loop.evaluate(phase) for phase in ('0','.137','.337','.5','.663','.863','1')}
            records[chart]=dict(source=source['record'],conditioned_geometry=loop.geometry_record(),candidate_phase_queries=phases)
            print('Live conditioned phase:',chart,loop.geometry,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_source_box_count=len(records),
        actual_C0_phase_inverse_and_primitive_records=records,transformed_original_antiderivative_identity=transformed_identity(),
        native_source_midpoints_or_caps_used_as_field_values=False,spatial_phase_binding_installed=False,
        actual_changed_defect_integral_functions_installed=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='C0 two-coordinate inverse and original A/B enclosures at explicit free periodic phases on five original source boxes. Not whole-chart/radius-phase/slow-jets/densities/integrals/repair/common-N/recursive/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
