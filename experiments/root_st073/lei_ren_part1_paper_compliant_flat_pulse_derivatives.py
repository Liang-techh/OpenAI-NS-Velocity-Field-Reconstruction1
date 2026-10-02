"""Original O.4 flat shapes: directed ordinary derivatives through order four.

Endpoint-crossing boxes use analytic flat-tail majorants, never a reciprocal
of an interval containing zero. Tiny exponentials are enclosed by positive
caps; the exact sigma/beta definitions are not replaced by those caps.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_axial_pulse_field import gp as original_gp
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
SIGMA_CONSTANTS=(1,4,28,256,3104)


def polynomial_add(a,b):
    return [(a[k] if k<len(a) else 0)+(b[k] if k<len(b) else 0) for k in range(max(len(a),len(b)))]


def polynomial_mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    while len(out)>1 and out[-1]==0:out.pop()
    return out


def beta_polynomials(order=4):
    rows=[[1]]; w=[1,0,-1]
    for n in range(order):
        derivative=[k*rows[-1][k] for k in range(1,len(rows[-1]))] or [0]
        rows.append(polynomial_add(polynomial_mul(polynomial_mul(w,w),derivative),
            polynomial_mul([0,4*n-2,0,-4*n],rows[-1])))
    return rows


BETA_POLYNOMIALS=beta_polynomials()
BETA_CONSTANTS=tuple(sum(abs(v) for v in row) for row in BETA_POLYNOMIALS)


def upper_point(c,value):return c.mpf(endpoints(value)[1])


def symmetric(c,bound):
    upper=endpoints(bound)[1]
    return c.mpf([-upper,upper])


def positive_exp(c,logvalue):
    """Finite outward enclosure, with exp(-1000) only an upper cap."""
    lo,hi=endpoints(logvalue)
    cap=c.exp(-1000)
    if hi<=-1000:return c.mpf([0,endpoints(cap)[1]])
    lower=c.mpf(0) if lo<=-1000 else c.exp(c.mpf(lo))
    return c.mpf([endpoints(lower)[0],endpoints(c.exp(c.mpf(hi)))[1]])


def flat_tail_exp(c,logvalue,distance):
    """A vanishing positive cap, admitted by its exact log inequality.

    For sufficiently small W the analytic exponential is <=W*exp(-1000).
    Materializing that cap is safe even when exp(logvalue) is not. This
    preserves a vanishing numerical upper envelope as W tends to zero.
    """
    distance=upper_point(c,distance)
    cap_log=-1000+c.ln(distance)
    if endpoints(logvalue)[1]<=endpoints(cap_log)[0]:
        cap=distance*c.exp(-1000)
        return c.mpf([0,endpoints(cap)[1]])
    return positive_exp(c,logvalue)


def sigma_tail_bound(c,n,W):
    """Supremum on 0<m<=W of Cn exp(4-m^-2)m^-3n."""
    W=upper_point(c,c.mpf(W))
    if endpoints(W)[1]<=0:return c.mpf(0)
    W=c.mpf(min(endpoints(W)[1],mp.mpf('.5')))
    if n==0:return flat_tail_exp(c,4-1/W**2,W)
    peak=c.sqrt(c.mpf(2)/(3*n))
    rho=W if endpoints(W)[1]<endpoints(peak)[0] else peak
    return flat_tail_exp(c,c.ln(SIGMA_CONSTANTS[n])+4-1/rho**2-3*n*c.ln(rho),W)


def beta_tail_bound(c,n,W):
    """Supremum on 0<w<=W; beta peak is w=1/(2n), not sqrt."""
    W=upper_point(c,c.mpf(W))
    if endpoints(W)[1]<=0:return c.mpf(0)
    W=c.mpf(min(endpoints(W)[1],mp.mpf(1)))
    if n==0:return flat_tail_exp(c,-1/W,W)
    peak=c.mpf(1)/(2*n)
    eta=W if endpoints(W)[1]<endpoints(peak)[0] else peak
    return flat_tail_exp(c,c.ln(BETA_CONSTANTS[n])-1/eta-2*n*c.ln(eta),W)


def hull_jets(c,rows):
    return IntervalTaylor(c,[c.mpf([min(endpoints(row[n])[0] for row in rows),
        max(endpoints(row[n])[1] for row in rows)]) for n in range(5)])


def _sigma_left(c,x):
    """The left half [0,.5], including its flat endpoint."""
    lo,hi=endpoints(x)
    if hi<=0:return IntervalTaylor.constant(c,0,4)
    if lo<=0:
        value=c.mpf([0,min(mp.mpf('.5'),endpoints(sigma_tail_bound(c,0,hi))[1])])
        return IntervalTaylor(c,[value]+[symmetric(c,sigma_tail_bound(c,n,hi))/math.factorial(n) for n in range(1,5)])
    odds=1/(1-x)**2-1/x**2
    if endpoints(odds)[1]<-1000:
        return IntervalTaylor(c,[c.mpf([0,endpoints(sigma_tail_bound(c,0,hi))[1]])]+
            [symmetric(c,sigma_tail_bound(c,n,hi))/math.factorial(n) for n in range(1,5)])
    e=positive_exp(c,odds); value=e/(1+e)
    # Monotonic endpoint evaluations tighten C0 without assuming interval
    # independence of the numerator and denominator.
    values=[]
    for v in (lo,hi):
        xx=c.mpf(v); ee=positive_exp(c,1/(1-xx)**2-1/xx**2)
        values.append(ee/(1+ee))
    value=c.mpf([max(mp.mpf(0),endpoints(values[0])[0]),min(mp.mpf('.5'),endpoints(values[1])[1])])
    L1=2/(1-x)**3+2/x**3
    L2=6/(1-x)**4-6/x**4
    L3=24/(1-x)**5+24/x**5
    L4=120/(1-x)**6-120/x**6
    q=value*(1-value); a=1-2*value
    derivatives=[value,q*L1,q*(L2+a*L1**2),
        q*(L3+3*a*L1*L2+(1-6*q)*L1**3),
        q*(L4+a*(4*L1*L3+3*L2**2)+6*(1-6*q)*L1**2*L2+a*(1-12*q)*L1**4)]
    if lo==hi==mp.mpf('.5'):derivatives[0]=c.mpf('.5'); derivatives[2]=derivatives[4]=c.mpf(0)
    return IntervalTaylor(c,[v/math.factorial(n) for n,v in enumerate(derivatives)])


def sigma_jets(c,x):
    """Ordinary Taylor coefficients at every point of a real interval."""
    x=c.mpf(x); lo,hi=endpoints(x); rows=[]
    if lo<=0:rows.append(IntervalTaylor.constant(c,0,4))
    if hi>=1:rows.append(IntervalTaylor.constant(c,1,4))
    if hi>0 and lo<=mp.mpf('.5'):
        rows.append(_sigma_left(c,c.mpf([max(mp.mpf(0),lo),min(mp.mpf('.5'),hi)])))
    if lo<1 and hi>=mp.mpf('.5'):
        a=max(mp.mpf('.5'),lo); b=min(mp.mpf(1),hi)
        left=_sigma_left(c,1-c.mpf([a,b]))
        rows.append(IntervalTaylor(c,[1-left[0]]+[(-1)**(n+1)*left[n] for n in range(1,5)]))
    return hull_jets(c,rows)


def beta_jets(c,r):
    """Original exp(-1/(1-r²)) compact bump, including support crossings."""
    r=c.mpf(r); lo,hi=endpoints(r)
    if hi<=-1 or lo>=1:return IntervalTaylor.constant(c,0,4)
    near=mp.mpf(0) if lo<=0<=hi else min(abs(lo),abs(hi))
    far=max(abs(lo),abs(hi)); wmax=1-c.mpf(near)**2
    upper=beta_tail_bound(c,0,wmax)
    lower=c.mpf(0) if far>=1 else positive_exp(c,-1/(1-c.mpf(far)**2))
    value=c.mpf([endpoints(lower)[0],endpoints(upper)[1]])
    if lo<=-1 or hi>=1 or endpoints(-1/wmax)[1]<-1000:
        return IntervalTaylor(c,[value]+[symmetric(c,beta_tail_bound(c,n,wmax))/math.factorial(n) for n in range(1,5)])
    w=1-r**2; derivatives=[value]
    for n,row in enumerate(BETA_POLYNOMIALS[1:],1):
        polynomial=c.mpf(0)
        for coefficient in reversed(row):polynomial=polynomial*r+coefficient
        derivatives.append(value*polynomial/w**(2*n))
    if lo==hi==0:derivatives[1]=derivatives[3]=c.mpf(0)
    return IntervalTaylor(c,[v/math.factorial(n) for n,v in enumerate(derivatives)])


def _gp_piece(c,xi):
    """Primitive(sigma(50 xi))*sigma(11-xi), with exact flat endpoints."""
    xi=c.mpf(xi); lo,hi=endpoints(xi)
    if hi<=0 or lo>=11:return IntervalTaylor.constant(c,0,4)
    entrance=sigma_jets(c,50*xi)
    primitive=original_gp(c,xi)['value'] if hi<=10 else None
    # Original gp value equals the primitive before the exit transition;
    # beyond xi=.02 the primitive is exactly xi-.01 by cutoff symmetry.
    if lo>=mp.mpf('.02'):primitive=xi-c.mpf('.01')
    elif hi>10:raise ValueError('Split a box spanning the entrance and exit transitions')
    coefficients=[primitive]
    for n in range(1,5):coefficients.append(entrance[n-1]*50**(n-1)/n)
    P=IntervalTaylor(c,coefficients); exit_shape=sigma_jets(c,11-xi)
    exit_shape=IntervalTaylor(c,[(-1)**n*exit_shape[n] for n in range(5)])
    return P*exit_shape


def gp_jets(c,xi):
    """Split entrance/middle/exit boxes before applying the original product."""
    xi=c.mpf(xi); lo,hi=endpoints(xi)
    if lo==hi:return _gp_piece(c,xi)
    cuts=[lo]+[mp.mpf(v) for v in ('0','.02','10','11') if lo<mp.mpf(v)<hi]+[hi]
    return hull_jets(c,[_gp_piece(c,c.mpf([a,b])) for a,b in zip(cuts,cuts[1:])])


class FlatPulseDerivatives:
    def __init__(self,ctx=None):
        self.ctx=ctx or MPIntervalContext()
        if ctx is None:self.ctx.dps=180
        name=PREFIX+'compliant_pulse_radial_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['pulse_radial_velocity_axial_C4_installed']:
            raise ValueError('Actual fifth-order pulse prerequisite required')
        self.hashes=dict(receipt['input_hashes'])
        for source,digest in self.hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Flat derivative source changed: '+source)
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        pulse_name=PREFIX+'compliant_axial_pulse_field.json'; pulse=json.loads((HERE/pulse_name).read_bytes())
        self.hashes[pulse_name]=hashlib.sha256((HERE/pulse_name).read_bytes()).hexdigest()
        self.normalization=read_interval(self.ctx,pulse['bump_normalization'])
        if endpoints(self.normalization)[0]<=0:raise ValueError('Positive original beta normalization required')
        self.family=receipt['actual_five_defect_family_sha256']; self.source=receipt['implicit_source_sha256']
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def beta(self,s):
        c=self.ctx; ell=c.mpf('.15'); raw=beta_jets(c,c.mpf(s)/ell)
        return IntervalTaylor(c,[raw[n]/(ell**(n+1)*self.normalization) for n in range(5)])

    def report(self):
        c=self.ctx
        with mp.workdps(210):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                sigma_derivative_constants=SIGMA_CONSTANTS,beta_derivative_polynomials=BETA_POLYNOMIALS,
                beta_derivative_constants=BETA_CONSTANTS,bump_normalization=self.normalization,
                sigma_samples=[dict(x=c.mpf(x),jets=sigma_jets(c,x)) for x in
                    ('0','.01','.25','.5','.9','1',['-.001','.001'],['.999','1.001'])],
                beta_samples=[dict(r=c.mpf(x),jets=beta_jets(c,x)) for x in
                    ('-1','-.99','-.5','0','.75','1',['-1.001','-.999'],['.999','1.001'])],
                gp_samples=[dict(xi=c.mpf(x),jets=gp_jets(c,x)) for x in
                    ('0','.01','.02','1','10','10.5','11',['-.0001','.0001'],['10.999','11.001'])],
                normalized_beta_support_crossings=[dict(s=c.mpf(x),jets=self.beta(x)) for x in
                    (['-.151','-.149'],['.149','.151'])],
                sigma_global_derivative_bounds=[sigma_tail_bound(c,n,'.5') for n in range(5)],
                beta_global_derivative_bounds=[beta_tail_bound(c,n,1) for n in range(5)],
                endpoint_tail_bounds=[dict(W=c.mpf(W),sigma=[sigma_tail_bound(c,n,W) for n in range(5)],
                    beta=[beta_tail_bound(c,n,W) for n in range(5)]) for W in ('.1','.03','.01')],
                retained_radial_shape_order=4,original_flat_shapes_retained=True,
                support_crossing_derivatives_available=True,
                quantitative_flat_tail_definition=dict(sigma='Cn exp(4-1/m²)m^(-3n), m=min(x,1-x)',
                    sigma_order_zero_scope='endpoint deviation: sigma near0, 1-sigma near1; global |sigma|<=1',
                    beta='An exp(-1/w)w^(-2n), w=1-r²',
                    analytic_endpoint_deviations_and_derivatives_through4_tend_to_zero=True,
                    numerical_tail_upper_bounds_tend_to_zero=True,
                    positive_tail_cap='W*exp(-1000), used only after proving logbound<=-1000+log(W)',
                    numerical_cap_is_only_an_enclosure=True),
                full_pulse_C4_installed=False,full_outer_C4_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=FlatPulseDerivatives().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Original sigma/gp/beta radial derivatives through4 and quantitative flat support tails generated',flush=True)
    return result


if __name__=='__main__':run()
