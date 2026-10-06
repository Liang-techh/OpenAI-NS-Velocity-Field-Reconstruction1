"""Current full collar stress mixed4, including the actual two constants.

Extends the original sigma/phi/Gamma shape to radial5. The same infinite
future moments and actual current pressure/angular histories are retained.
"""
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_heat_pressure_stress import (
    CurrentHeatPressureStress,HERE,PREFIX,sha,pack,encode,endpoints,mixed,
    constant_stress_rows,SCOPES,IntervalTaylor,binding)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import (
    phi_jets,PHI_POLYNOMIALS,gamma_deficit_mixed,product_rows)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import (
    sigma_jets,flat_tail_exp,polynomial_add,polynomial_mul,symmetric)
from lei_ren_part1_paper_compliant_collar_stress_C3 import (
    axial_derivative,shifted_rows,collar_defect_rows)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_collar_stress_mixed_C4.json'
RECEIPT=PREFIX+'current_collar_stress_mixed_C4_check.json'
GATE='current_collar_actual_stress_mixed4_recovered'
SHAPE_GATE='current_full_collar_K_radial5_certified'
JOIN_GATE='current_collar_Gamma_actual_stress_mixed4_join_certified'
SIGMA5_CONSTANT=44704
VIEWS={'whole_collar':([-1,1],[0,3]),'inlet':([-1,1],0),
    'sigma_crossing':([-1,1],[0,1]),'phi_crossing':([-1,1],['2.99',3]),
    'Gamma_join':([-1,1],3),'fresh_collar':('.407','.41')}


def phi_fifth_polynomial():
    p=PHI_POLYNOMIALS[-1];n=4
    derivative=[j*p[j] for j in range(1,len(p))] or [0]
    return polynomial_add(polynomial_add([v*(-8) for v in p],polynomial_mul([0,0,3*n],p)),
        polynomial_mul([0,0,0,-1],derivative))


PHI5_POLYNOMIAL=phi_fifth_polynomial()


def sigma_fifth_formula(value,x):
    """Exact fifth derivative of the original logistic odds function."""
    q=value*(1-value);a=1-2*value
    L1=2/(1-x)**3+2/x**3;L2=6/(1-x)**4-6/x**4
    L3=24/(1-x)**5+24/x**5;L4=120/(1-x)**6-120/x**6
    L5=720/(1-x)**7+720/x**7
    return q*(L5+a*(5*L1*L4+10*L2*L3)
        +(1-6*q)*(10*L1**2*L3+15*L1*L2**2)
        +a*(1-12*q)*10*L1**3*L2+(1-30*q+120*q**2)*L1**5)


def sigma_fifth_left(c,x):
    lo,hi=endpoints(x)
    if hi<=0:return c.mpf(0)
    # q<=sigma<=exp(4-x^-2); Bell5 and |L_j|<=2(j+1)!x^(-j-2)
    # give 44704*exp(4-x^-2)*x^-15 on the entire left half.
    if lo<=0 or endpoints(1/(1-x)**2-1/x**2)[1]<-1000:
        width=c.mpf(hi);peak=c.sqrt(c.mpf(2)/15)
        rho=width if endpoints(width)[1]<endpoints(peak)[0] else peak
        bound=flat_tail_exp(c,c.ln(c.mpf(SIGMA5_CONSTANT))+4-1/rho**2-15*c.ln(rho),width)
        return symmetric(c,bound)
    return sigma_fifth_formula(sigma_jets(c,x)[0],x)


def sigma_C5(c,x):
    x=c.mpf(x);lo,hi=endpoints(x);values=[]
    if lo<=0 or hi>=1:values.append(c.mpf(0))
    if hi>0 and lo<=mp.mpf('.5'):
        values.append(sigma_fifth_left(c,c.mpf([max(mp.mpf(0),lo),min(mp.mpf('.5'),hi)])))
    if lo<1 and hi>=mp.mpf('.5'):
        # sigma(1-x)=1-sigma(x): the fifth derivative has positive parity.
        values.append(sigma_fifth_left(c,1-c.mpf([max(mp.mpf('.5'),lo),min(mp.mpf(1),hi)])))
    fifth=c.mpf([min(endpoints(v)[0] for v in values),max(endpoints(v)[1] for v in values)])
    prior=sigma_jets(c,x)
    return IntervalTaylor(c,list(prior.coefficients)+[fifth/math.factorial(5)])


def phi_C5(c,t):
    t=c.mpf(t);lo,hi=endpoints(t);prior=phi_jets(c,t)
    if lo>=3:fifth=c.mpf(0)
    else:
        dist=c.mpf(3)-c.mpf(lo);width=c.mpf(endpoints(dist)[1])
        flat=endpoints(-4/width**2)[1]<-1000
        if hi>=3 or flat:
            coefficient=sum(abs(v)*3**j for j,v in enumerate(PHI5_POLYNOMIAL))
            peak=c.sqrt(c.mpf(8)/15)
            rho=width if endpoints(width)[1]<endpoints(peak)[0] else peak
            bound=flat_tail_exp(c,c.ln(c.mpf(coefficient))-4/rho**2-15*c.ln(rho),width)
            fifth=symmetric(c,bound)
        else:
            distance=IntervalTaylor.variable(c,3-t,5)
            distance=IntervalTaylor(c,[distance[0],-1,0,0,0,0])
            fifth=((distance**(-2)*(-4)).exp())[5]*math.factorial(5)
    return IntervalTaylor(c,list(prior.coefficients)+[fifth/math.factorial(5)])


def shape_radial5(heat,Z,t):
    c=heat.ctx;Z=c.mpf(Z);t=c.mpf(t)
    prior=heat.shape(Z,t);sig=sigma_C5(c,t);phi=phi_C5(c,t)
    zero=IntervalTaylor.constant(c,0,5);one=zero+1
    sr=[one*(sig[j]*math.factorial(j)) for j in range(6)]
    fr=[one*(phi[j]*math.factorial(j)) for j in range(6)]
    unity=[one]+[zero]*5
    W=[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]
    W[0]=prior['W_rows'][0]
    C=[unity[j]-fr[j]*heat.eps for j in range(6)]
    pre=[unity[j]-W[j]*heat.eps for j in range(6)]
    dr=gamma_deficit_mixed(c,Z,heat.a,heat.Scap,t,5)
    D=product_rows(product_rows(sr,C),dr)
    K=[pre[j]-D[j]*(heat.a*heat.S) for j in range(6)]
    # Existing bounds remain literally unchanged. All six rows enclose the
    # same original product/chain derivatives by the source identities.
    for name,rows in (('K_rows',K),('D_rows',D),('W_rows',W),('pre_rows',pre)):
        rows[:5]=prior[name]
    return dict(K_rows=K,D_rows=D,W_rows=W,pre_rows=pre,sigma=sig,phi=phi,
        original_first_five_shape_rows_retained=True,full_positive_Gamma_used=True)


def collar_stress_mixed4(heat,shape,defects,Z,t):
    c=heat.ctx;z=IntervalTaylor.variable(c,Z,5);d=1-z*z;L=1-z*z*heat.delta
    b=(1-heat.delta)/2;A=defects['angular_defect_rows'];E=defects['energy_defect_rows']
    P=defects['pressure_defect_rows'];Km=defects['K_defect_rows'];K=shape['K_rows']
    inertial=[(A[j]*heat.k-z*axial_derivative(A[j])*b-Km[j])/L for j in range(5)]
    axial=[(z*E[j]*heat.delta-d*axial_derivative(E[j])/2
        -z*P[j]*(2*heat.prate)+d*axial_derivative(P[j]))/L for j in range(5)]
    viscous=[K[j+1]-K[j]*(1+heat.a) for j in range(5)]
    shear=[v*(2*heat.S*c.exp(-t)) for v in shifted_rows(viscous,-1,4)]
    theta=[inertial[j]+shear[j] for j in range(5)]
    return dict(theta=shifted_rows(theta,-heat.a,4),
        axial=shifted_rows(axial,-c.mpf('.5')-heat.delta,4))


class CurrentCollarStressMixedC4:
    @source_precision
    def __init__(self,companion=None,require_checked=True):
        self.companion=companion if companion is not None else CurrentHeatPressureStress()
        if not self.companion.acceptance_loaded:raise ValueError('Checked actual current heat companion required')
        self.heat=self.companion.heat;self.ctx=self.heat.ctx
        self.family=self.companion.family;self.source=self.companion.source;self.datum_sha=self.companion.datum_sha
        self.hashes=dict(self.companion.hashes);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            record=accepted(RECEIPT,self.family,self.source,GATE);_verify_hashes(record)
            if not record[SHAPE_GATE] or not record[JOIN_GATE] or any(record[k] for k in SCOPES) or record['datum_enclosure_sha256']!=self.datum_sha:
                raise ValueError('Current collar mixed4 source/scope differs')
            self.hashes.update(record['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def evaluate(self,Z,t):
        c=self.ctx;Z=c.mpf(Z);t=c.mpf(t)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(t)[0]<0 or endpoints(t)[1]>3:
            raise ValueError('Original collar Z in[-1,1], t in[0,3] required')
        key=(Z._mpi_,t._mpi_)
        if key in self.cache:return self.cache[key]
        before=self.companion.evaluate('heat_collar',Z,t)
        shape=shape_radial5(self.heat,Z,t);tails=self.heat.collar_tails(Z,t)
        defects=collar_defect_rows(self.heat,shape,tails,t)
        canonical=collar_stress_mixed4(self.heat,shape,defects,Z,t)
        constants=before['actual_current_terminal_constants']
        extra=constant_stress_rows(self.heat,Z,t,constants['angular_defect'],constants['pressure_infinity'],4)
        stress=dict(theta_Qtheta=mixed([canonical['theta'][j]+extra['theta'][j] for j in range(5)],4),
            axial_Qz=mixed(canonical['axial'],4),axial_pressure_Qpressure=mixed(extra['axial_pressure_constant'],4))
        result=dict(Z=Z,coordinate=t,actual_stress_factored_mixed4=stress,
            current_full_collar_shape_radial5=shape,
            actual_current_terminal_constants=constants,
            current_absolute_pressure_mixed4=before['stable_same_forward_pressure_mixed4'],
            positive_source_stress_factors=before['positive_source_stress_factors'],
            same_current_common_heat_source_graph=self.companion.graph,
            same_original_full_future_integrals=True,angular_and_pressure_constants_retained=True,
            **dict.fromkeys((GATE,SHAPE_GATE,JOIN_GATE),self.acceptance_loaded),**dict.fromkeys(SCOPES,False))
        self.cache[key]=result;return result


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCollarStressMixedC4(require_checked=False)
    views={name:field.evaluate(Z,t) for name,(Z,t) in VIEWS.items()}
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,actual_current_collar_stress_mixed4_views=views,
        shape_radial_order=5,stress_total_mixed_order=4,
        actual_current_constants_not_reset=True,
        **dict.fromkeys((GATE,SHAPE_GATE,JOIN_GATE),False),**dict.fromkeys(SCOPES,False),input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(encode(pack(result)),indent=2)+'\n').encode('utf8'))
    print('Built full current collar radial5 and actual stress mixed4, constants retained',flush=True)
    return result


if __name__=='__main__':run()
