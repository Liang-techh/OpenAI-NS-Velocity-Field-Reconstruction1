"""Directed complex-axis pole exclusion and conditional X_h radial tails.

Implements the coefficient norm in paper (8.30). The nonlinear X_h norm is
an explicit missing input; a conditional tail factor is not a core proof.
"""
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def tail_factor(ctx,*,degree,radial_order,axial_order,radius,h):
    """Bound mixed derivative of the tail divided by its X_h norm.

For n>N, derivative terms have ratio at most
    radius/20 * (1+(i+k)/(N+2-i)).
The norm includes derivative_m/m!, binomial(n+m,m), and square weights.
All radial derivatives here use scaled R=Lambda R, not x=R/Ra.
"""
    N=int(degree);i=int(radial_order);k=int(axial_order)
    if min(N,i,k)<0 or i>N+1:raise ValueError('Unsupported tail index')
    a=ctx.mpf(radius);h=ctx.mpf(h)
    if endpoints(a)[0]<=0 or endpoints(h)[0]<=0:raise ValueError('Positive radius and h required')
    n=N+1
    rho=a/20*(1+ctx.mpf(i+k)/(N+2-i))
    if endpoints(rho)[1]>=1:raise ValueError('Geometric tail majorant does not contract')
    first=(ctx.mpf(math.factorial(k))*h**(-k)*ctx.mpf(math.factorial(n)//math.factorial(n-i))
           *a**(n-i)*ctx.mpf(math.comb(n+k,k))
           /(ctx.mpf(20)**n*(n+1)**2*(k+1)**2))
    return dict(first_term=first,ratio_upper=rho,tail_per_Xh_norm=first/(1-rho))


def run():
    ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        j=ctx.mpf('1e-14');dt=ctx.mpf('1e-200');sigma=j/500
        eta=sigma/40;h=eta/8;a=1+eta
        # H0=(9-delta)Z/2+j-4Z^3-jZ^2. On every disk
        # |Z-x|<=eta, x real in[-1,1], the path derivative bound follows.
        H_derivative=(9+dt)/2+2*j*a+12*a*a
        displacement=eta*H_derivative
        pole_factor_lower=sigma-displacement
        if endpoints(H_derivative)[1]>=20 or endpoints(pole_factor_lower)[0]<=0:
            raise AssertionError('Common complex tube exclusion failed')
        Llower=1-dt*a*a
        Qlower=1-eta*(2+eta)
        if min(endpoints(Llower)[0],endpoints(Qlower)[0])<=0:
            raise AssertionError('L or pressure q pole not excluded')
        # H0 real on the real axis, so distance to +/-i sigma >=sigma;
        # variation controls both complex denominator factors uniformly.
        chi_upper=1+sigma*sigma/pole_factor_lower**2
        beta_upper=(3+dt/2+j*a)/Llower
        Hupper=(9+dt)*a/2+j+4*a**3+j*a*a
        Lupper=1+dt*a*a
        complex_gradient_upper=Lupper*Hupper/pole_factor_lower**2
        # Follow the real path from the root in[-j,0], then a segment of
        # length<=eta inside its disk. The capsule is simply connected and
        # pole-free, so the analytic primitive is the same axis G.
        real_G_upper=(1+j)*(1+dt)/(2*sigma)
        complex_G_upper=real_G_upper+eta*complex_gradient_upper
        lam=ctx.mpf('1e36');logC=ctx.mpf('5e151')
        logC_required=2*ctx.log(lam)+lam*complex_G_upper
        Cstar_margin=logC-logC_required
        if endpoints(Cstar_margin)[0]<=0:raise AssertionError('Axis Cstar guard not proved')
        axis_F_log_modulus_upper=-logC+lam*complex_G_upper
        rows=[]
        for N in (18,20):
            for total in range(4):
                for i in range(total+1):
                    k=total-i
                    bound=tail_factor(ctx,degree=N,radial_order=i,axial_order=k,radius='4.1',h=h)
                    rows.append(dict(degree=N,scaled_radial_order=i,axial_order=k,**bound))
        report=dict(j=j,delta=dt,sigma=sigma,complex_tube_radius=eta,Xh_parameter=h,
             H0_derivative_modulus_upper=H_derivative,H0_variation_upper=displacement,
             denominator_factor_modulus_lower=pole_factor_lower,L_modulus_lower=Llower,
             pressure_q_modulus_lower=Qlower,chi_modulus_upper=chi_upper,beta_modulus_upper=beta_upper,
             H0_modulus_upper=Hupper,complex_gradient_modulus_upper=complex_gradient_upper,
             G_modulus_upper=complex_G_upper,A_Omega_upper=complex_G_upper,
             logC_required_upper=logC_required,Cstar_log_margin=Cstar_margin,
             axis_F_log_modulus_upper=axis_F_log_modulus_upper,
             Cstar_complex_axis_guard_certified=True,
             common_complex_axis_poles_excluded=True,
             tube_definition='union of disks |Z-x|<=eta over real x in[-1,1]',
             norm_definition='paper(8.30):20^n h^m (n+1)^2(m+1)^2 |d_Z^m f_n|/[m! binomial(n+m,m)]',
             tail_factors=rows,
             conditional_tail_inequality='|d_scaledR^i d_Z^k tail_N| <= certified_Xh_norm * tail_per_Xh_norm',
             actual_nonlinear_Xh_norm_certified=False,analytic_P0_modulus_certified=False,
             nonlinear_contraction_certified=False,infinite_core_remainder_enclosed=False,
             original_parameter_derivation_certified=False,temporal_recursion=False,
             next_dependency='Bound the accepted analytic pressure and normalized coupled nonlinear fixed point in X_h; do not supply an arbitrary norm constant')
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('axis poles excluded; eta',mp.nstr(endpoints(eta)[1],20),'h',mp.nstr(endpoints(h)[1],20),flush=True)
        print('conditional mixed C3 tail factors',len(rows),'actual nonlinear norm remains open',flush=True)
        return report


if __name__=='__main__':run()
