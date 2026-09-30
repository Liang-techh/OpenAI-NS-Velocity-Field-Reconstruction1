"""Actual row-normalized pulse integrals in (7.31), for small mu.

Near the flat endpoint xi=11 the dominant width is O(mu**(2/3)).
A saddle-centered change of variable evaluates the integral without losing
the narrow peak or materializing exp(-lambda*2/mu). Omitted positive pieces
have explicit upper bounds; finite quadrature error is reported separately.
"""
import json
import math
from pathlib import Path
import mpmath as mp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from lei_ren_part1_paper_outer import PaperOuterSchedule


def normalized_pulse_integral(mu, row, *, order=128, band=20., precision=120):
    if row not in (1,2):
        raise ValueError('row must be 1 or 2')
    with mp.workdps(precision):
        mu=mp.mpf(str(mu)); lam=mp.mpf('.5')-row*mu
        if not 0<mu<=mp.mpf('1e-6'):
            raise ValueError('This centered evaluator requires 0<mu<=1e-6')
        k=lam/mu; u0=mp.root(2/k,3); L=1/(u0*u0)
        relative_width=1/mp.sqrt(6*L)
        if band*relative_width>=mp.mpf('.25'):
            raise ValueError('Increase centering accuracy for this mu/band')
        width=float(relative_width); center=float(u0)
        if center==0 or width==0 or center*center==0:
            raise ValueError('Centered local coordinates underflow; arbitrary precision local quadrature is required')
        nodes,weights=leggauss(order)
        total=0.
        for node,weight in zip(nodes,weights):
            x=band*float(node); v=1+x*width; u=center*v
            phase_difference=x*x/6*(2*v+1)/(v*v)
            smooth=1/(1-u)**2
            tiny=-1/u**2+smooth
            denominator=1+math.exp(tiny) if tiny>-745 else 1.
            total+=band*float(weight)*(10.99-u)*math.exp(smooth-phase_difference)/denominator
        if total<=0:
            raise ArithmeticError('Centered quadrature lost pulse mass')
        # du/dx=u0^2/sqrt(6), dxi=-du; original integral contains 1/mu.
        log_prefactor=2*mp.log(u0)-mp.log(6)/2-mp.log(mu)
        log_integral=-2*k-3*L+log_prefactor+mp.log(total)
        # Positive omitted pieces: band complement in u in (0,.5),
        # u in [.5,1], and xi<=10. No tails are silently set to zero.
        phases=[]
        for sign in (-1,1):
            v=1+sign*band*relative_width
            phases.append(mp.mpf(str(band))**2/6*(2*v+1)/(v*v))
        log_bounds=[mp.log(mp.mpf('5.5'))+4-mp.log(mu)-2*k-3*L-min(phases),
                    mp.log(mp.mpf('5.5'))-mp.log(mu)-mp.mpf('2.5')*k,
                    mp.log(11/lam)-3*k]
        # Report each bound separately; logsumexp is valid at arbitrary exponent.
        largest=max(log_bounds)
        log_bound=largest+mp.log(sum(mp.exp(v-largest) for v in log_bounds))
        fmt=lambda x:mp.nstr(x,precision)
        return {'row':row,'mu':fmt(mu),'lambda':fmt(lam),
                'log_normalized_pulse_integral':fmt(log_integral),
                'log_relative_omitted_positive_bound':fmt(log_bound-log_integral),
                'log_omitted_piece_bounds':[fmt(v) for v in log_bounds],
                'center_u':fmt(u0),'centered_integral':total,
                'quadrature_order':order,'band':band,'precision':precision,
                'scope':'Actual centered pulse integral with omitted positive bounds; quadrature uncertainty is separate.'}


def run():
    schedule=PaperOuterSchedule(logPstar=14,logRref=10,delta='1e-32',Md='.5',
                               c_mu='.001',c_delta='.001',c_epsilon='.01')
    rows=[]
    for index in (1,2):
        coarse=normalized_pulse_integral(schedule.mu,index,order=64)
        fine=normalized_pulse_integral(schedule.mu,index,order=128)
        fine['centered_quadrature_relative_difference']=abs(
            fine['centered_integral']-coarse['centered_integral'])/fine['centered_integral']
        if fine['centered_quadrature_relative_difference']>1e-6:
            raise ArithmeticError('Centered pulse quadrature failed refinement')
        rows.append(fine)
    # Independent original-u adaptive quadrature at a representable test mu.
    # Its exponent is computed directly, not by the centered phase identity.
    check=normalized_pulse_integral('1e-6',1,order=192)
    test_mu=1e-6; k=(.5-test_mu)/test_mu; u0=(2/k)**(1/3)
    peak=-k*u0-1/u0**2+1/(1-u0)**2
    def original_integrand(u):
        if u<=0 or u>=.5:return 0.
        exponent=-k*u-1/u**2+1/(1-u)**2-peak
        return (10.99-u)*math.exp(exponent) if exponent>-745 else 0.
    integral=0.; error=0.
    for left,right in ((0.,u0),(u0,2*u0),(2*u0,.5)):
        value,err=quad(original_integrand,left,right,epsabs=1e-13,epsrel=1e-11,limit=150)
        integral+=value; error+=err
    with mp.workdps(120):
        direct_log=-2*mp.mpf(str(k))+mp.mpf(str(peak))+mp.log(mp.mpf(str(integral))/mp.mpf('1e-6'))
        relative_difference=float(abs(mp.expm1(direct_log-mp.mpf(check['log_normalized_pulse_integral']))))
    if relative_difference>1e-8:
        raise ArithmeticError('Independent original-u integral disagrees')
    report={'source':'https://arxiv.org/html/2609.35406v1','equations':'7.28,7.31',
            'schedule':schedule.metadata(),'normalized_pulse_rows':rows,
            'independent_original_u_check':{'mu':test_mu,
                'relative_difference':relative_difference,'quadrature_error':error},
            'full_axial_closure':False,
            'scope':'Pulse RHS rows only; incoming moments, affine correction and actual energy target are separate.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'refinement':[v['centered_quadrature_relative_difference'] for v in rows],
                      'omitted_relative_log_bounds':[v['log_relative_omitted_positive_bound'] for v in rows]}))
    return report


if __name__=='__main__':run()
