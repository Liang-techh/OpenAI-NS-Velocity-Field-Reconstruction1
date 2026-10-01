"""Directed global real-axis jet enclosures for the stored Section 8 data.

The G integral is bounded analytically by its path length and
|L H/(H^2+sigma^2)| <= sup|L|/(2 sigma). No quadrature is used.
This is conservative; original parameter derivations are excluded.
"""
import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,constant
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def uniform_axis_jets(iv,*,radius,j,Lambda,logC,delta,length):
    a=iv.mpf(radius);j=iv.mpf(j);lam=iv.mpf(Lambda);lc=iv.mpf(logC);dt=iv.mpf(delta)
    if not 0<endpoints(a)[0]<=endpoints(a)[1]<=1:raise ValueError('Require radius in (0,1]')
    if not 0<endpoints(j)[0]<=endpoints(j)[1]<endpoints(a)[0]:raise ValueError('Require 0<j<radius')
    if endpoints(lam)[0]<=0 or endpoints(dt)[0]<=0 or endpoints(dt)[1]>=1:raise ValueError('Positive Lambda and 0<delta<1 required')
    order=length-1;z=IntervalTaylor.variable(iv,iv.mpf([-endpoints(a)[1],endpoints(a)[1]]),order)
    one=constant(iv,1,order);u=4*z+j;L=one-(z*z)*dt;H=z*((1-dt)/2)+(one-z*z)*u
    # Root bracket: H(-j)<0<H(0), H' positive on [-j,0].
    left=-j;Hleft=(1-dt)*left/2+(1-left**2)*(4*left+j)
    derivative_lower=(9-dt)/2-12*j*j-2*j*j
    if endpoints(Hleft)[1]>=0 or endpoints(derivative_lower)[0]<=0:raise ValueError('Axis anchor root bracket not proved')
    sigma=j/500
    path_upper=a+j;Lupper=1+dt*a*a
    Gupper=path_upper*Lupper/(2*sigma)
    gmax=endpoints(Gupper)[1]
    f0=iv.exp(-lc-lam*iv.mpf([-gmax,gmax]))
    denominator=H*H+sigma*sigma
    # Preserve H^2>=0 at each real center. Interval multiplication of two
    # copies of a crossing-zero interval loses this exact square dependency.
    coefficients=list(denominator.coefficients)
    coefficients[0]=H.coefficients[0]**2+sigma*sigma
    denominator=IntervalTaylor(iv,coefficients)
    gradient=(L*H)/denominator
    f=[f0]
    for n in range(length-1):
        f.append(-lam*sum((gradient.coefficients[k]*f[n-k] for k in range(n+1)),iv.mpf(0))/(n+1))
    return dict(F0=f,U0=list(u.coefficients),G_absolute_upper=Gupper,
        F0_interval=f0,gradient_coefficients=list(gradient.coefficients),
        axis_anchor_root_bracket=(left,iv.mpf(0)),root_bracket_proved=True,
        uniform_axis_jets_enclosed=True,axis_integral_quadrature_used=False,
        original_parameter_errors_enclosed=False,nonlinear_core_remainder_enclosed=False)
