"""Independent original-shape derivatives and analytic endpoint envelopes."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import (
    FlatPulseDerivatives, sigma_jets, beta_jets, gp_jets, sigma_tail_bound,
    beta_tail_bound, SIGMA_CONSTANTS, BETA_POLYNOMIALS, BETA_CONSTANTS)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_flat_pulse_derivatives.json'


def symbolic_envelopes():
    r,w,n=s.symbols('r w n',real=True); x=s.symbols('x',positive=True)
    beta=s.exp(-1/(1-r*r)); checks={}
    for k,row in enumerate(BETA_POLYNOMIALS):
        polynomial=sum(s.Integer(v)*r**j for j,v in enumerate(row))
        actual=s.diff(beta,r,k)/beta*(1-r*r)**(2*k)
        if s.simplify(actual-polynomial)!=0:raise ArithmeticError('Original beta derivative polynomial mismatch')
        if sum(abs(v) for v in row)!=BETA_CONSTANTS[k]:raise ArithmeticError('Beta coefficient norm mismatch')
        checks['beta_direct_derivative_'+str(k)]=True
    # Logarithmic envelope derivatives establish the location of the global
    # maximum and monotone decay as the support endpoint is approached.
    if s.simplify(s.diff(-1/w-2*n*s.log(w),w)-(1-2*n*w)/w**2)!=0:
        raise ArithmeticError('Beta tail maximum is incorrect')
    if s.simplify(s.diff(4-1/x**2-3*n*s.log(x),x)-(2-3*n*x**2)/x**3)!=0:
        raise ArithmeticError('Sigma tail maximum is incorrect')
    checks['beta_peak_1_over_2n']=True; checks['sigma_peak_sqrt_2_over_3n']=True
    # q<=exp(4-m^-2); |Lj|<=2(j+1)!m^(-j-2), m<=.5.
    # The following sums dominate every term after reducing all powers to
    # m^(-3n). |1-2sigma|<=1, |1-6q|<=1, |1-12q|<=2.
    constants=(4,12+16,48+3*4*12+64,
        240+4*4*48+3*12**2+6*4**2*12+2*4**4)
    if constants!=SIGMA_CONSTANTS[1:]:raise ArithmeticError('Sigma envelope term constants mismatch')
    checks['sigma_termwise_constants']=True
    for power in range(13):
        if s.limit(s.exp(-1/x**2)*x**(-power),x,0,dir='+')!=0:
            raise ArithmeticError('Sigma flat limit failed')
    for power in range(9):
        if s.limit(s.exp(-1/x)*x**(-power),x,0,dir='+')!=0:
            raise ArithmeticError('Beta flat limit failed')
    checks['all_required_endpoint_flat_limits']=True
    return checks


def direct_fixture():
    with mp.workdps(85):
        c=MPIntervalContext(); c.dps=110; counts={}
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            a=mp.exp(-1/x**2); b=mp.exp(-1/(1-x)**2)
            return a/(a+b)
        def beta(x):return mp.exp(-1/(1-x*x)) if abs(x)<1 else mp.mpf(0)
        def primitive(x):
            if x<=0:return mp.mpf(0)
            if x>=mp.mpf('.02'):return x-mp.mpf('.01')
            return mp.quad(lambda v:sigma(50*v),[0,x/2,x])
        def gp(x):return primitive(x)*sigma(11-x)
        def target_derivative(label,function,x,k):
            if label!='gp':return mp.diff(function,x,k)
            if x<=0 or x>=11:return mp.mpf(0)
            # Independently integrate P, then use the fundamental theorem
            # and Leibniz rule. Differentiating adaptive quadrature inside
            # a high-order finite difference wastes precision and runtime.
            P=primitive(x); total=mp.mpf(0)
            for j in range(k+1):
                pj=P if j==0 else 50**(j-1)*mp.diff(sigma,50*x,j-1)
                total+=math.comb(k,j)*pj*(-1)**(k-j)*mp.diff(sigma,11-x,k-j)
            return total
        for label,provider,function,values in (
            ('sigma',sigma_jets,sigma,('.025','.2','.5','.8','.975')),
            ('beta',beta_jets,beta,('-.98','-.3','0','.4','.98')),
            ('gp',gp_jets,gp,('.01','.02','1','10','10.4','10.95'))):
            count=0
            for value in values:
                x=mp.mpf(value); jet=provider(c,c.mpf(value))
                for k in range(5):
                    target=target_derivative(label,function,x,k)/math.factorial(k); lo,hi=endpoints(jet[k])
                    # Independent quadrature/differentiation has finite
                    # precision. Its allowance is kept separate from the
                    # directed interval admission of the actual source.
                    allowance=mp.mpf('1e-65')*max(mp.mpf(1),abs(target))
                    if not lo-allowance<=target<=hi+allowance:
                        raise ArithmeticError('Direct original derivative failed: '+label+' '+str(x)+' '+str(k))
                    count+=1
            counts[label+'_independent_point_derivatives']=count
        # Explicit crossing boxes exercise exterior zeros and interior
        # tails on both sides; a broad box also exercises all chart splits.
        for label,provider,function,boxes in (
            ('sigma',sigma_jets,sigma,(('-.001','.001'),('.999','1.001'),('0','1'))),
            ('beta',beta_jets,beta,(('-1.001','-.999'),('.999','1.001'),('-1','1'))),
            ('gp',gp_jets,gp,(('-.0001','.0001'),('10.999','11.001'),('0','11')))):
            count=0
            for a,b in boxes:
                a,b=mp.mpf(a),mp.mpf(b); jet=provider(c,c.mpf([a,b]))
                for x in (a,a+(b-a)/4,(a+b)/2,b-(b-a)/4,b):
                    for k in range(5):
                        target=target_derivative(label,function,x,k)/math.factorial(k); lo,hi=endpoints(jet[k])
                        allowance=mp.mpf('1e-65')*max(mp.mpf(1),abs(target))
                        if not lo-allowance<=target<=hi+allowance:raise ArithmeticError('Support-crossing derivative failed: '+label)
                        count+=1
            counts[label+'_independent_box_derivatives']=count
        flat=FlatPulseDerivatives(c)
        for x in ('-.14','0','.14'):
            actual=flat.beta(x); raw=beta_jets(c,c.mpf(x)/c.mpf('.15'))
            for k in range(5):
                expected=raw[k]/(c.mpf('.15')**(k+1)*flat.normalization)
                if endpoints(actual[k])!=endpoints(expected):raise ArithmeticError('Normalized end beta derivative scale lost')
        # Sub-precision entrance values must use finite positive caps.
        tiny=c.exp(-10000)
        for provider in (sigma_jets,gp_jets):
            jet=provider(c,tiny)
            if endpoints(jet[0])[1]<=0:raise ArithmeticError('Exact tiny source was replaced by zero')
            for v in jet.coefficients:
                if any(not mp.isfinite(t) for t in endpoints(v)):raise ArithmeticError('Tiny source derivative unbounded')
        for k in range(5):
            for bound in (sigma_tail_bound(c,k,'.01'),beta_tail_bound(c,k,'.01')):
                if endpoints(bound)[1]<=0 or not mp.isfinite(endpoints(bound)[1]):raise ArithmeticError('Flat cap invalid')
            # The admitted numerical cap must also vanish with distance,
            # rather than stop at a fixed exp(-1000) floor.
            for bound_function in (sigma_tail_bound,beta_tail_bound):
                large=bound_function(c,k,c.exp(-10000))
                small=bound_function(c,k,c.exp(-20000))
                if not 0<endpoints(small)[1]<endpoints(large)[1]:raise ArithmeticError('Flat tail cap does not decrease with distance')
                W=c.exp(-20000); allowed=W*c.exp(-1000)
                if endpoints(small)[1]>endpoints(allowed)[1]:raise ArithmeticError('Log-admitted vanishing cap not retained')
        return dict(checks=counts,normalization_derivative_checks=15,tiny_source_positive_caps_checked=True,
            distance_dependent_positive_tail_caps_checked=True,
            independent_fixture_allowance='1e-65 relative/absolute; not an actual-source interval widening',
            actual_Md40_source_admission=False,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Flat derivative source changed: '+name)
    proof=symbolic_envelopes()
    print('Original beta symbolic derivatives and flat-envelope extrema PASS',flush=True)
    fixture=direct_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        analytic_original_beta_and_flat_envelope_checks=proof,independent_original_shape_fixture=fixture,
        original_radial_shape_derivatives_C4_available=True,quantitative_flat_support_majorants_available=True,
        support_crossing_derivatives_available=True,original_flat_shapes_retained=True,all_passed=True,
        full_pulse_C4_installed=False,full_outer_C4_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Original flat pulse derivatives: direct derivatives, support crossings, normalized scales and analytic flat envelopes PASS',flush=True)
    return result


if __name__=='__main__':run()
