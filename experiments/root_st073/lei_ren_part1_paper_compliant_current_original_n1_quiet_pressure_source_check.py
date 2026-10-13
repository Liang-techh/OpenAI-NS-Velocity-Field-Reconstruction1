"""Independent source/correlation and normalized integral checks."""
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source as current


def independent_kernel_reference():
    """Direct original sigma integrand, independently normalized for quadrature."""
    with mp.workdps(90):
        Y=mp.exp(40)
        def D(u):
            x=-mp.log1p(-u/Y)/40
            return 1/x**2-1/(1-x)**2
        root=mp.findroot(lambda u:mp.diff(lambda w:w+D(w),u),
                         (2*(40*Y)**2)**(mp.mpf(1)/3))
        reference=root+D(root)
        scale=mp.sqrt(root)
        def integrand(t):
            u=root+scale*t
            # This is exactly exp(-u)*sigma(x), with a shared normalizing log.
            return mp.exp(reference-u-D(u))/(1+mp.exp(-D(u)))
        integral=scale*mp.quad(integrand,[-20,-5,0,5,20])
        A=mp.log(integral)-reference
        F=A+mp.log1p(mp.exp(1-Y-A))
        return +A,+F


def _derivative_and_convexity_proof():
    u,Y=s.symbols('u Y',positive=True)
    x=-s.log(1-u/Y)/40
    D=1/x**2-1/(1-x)**2
    Dx=-2/x**3-2/(1-x)**3
    Dxx=6/x**4-6/(1-x)**4
    first=1+Dx/(40*(Y-u))
    second=(Dxx/1600+Dx/40)/(Y-u)**2
    assert s.simplify(s.diff(u+D,u)-first)==0
    assert s.simplify(s.diff(u+D,u,2)-second)==0
    # Uniform lower bound for Dxx+40Dx when 0<x<=1/100.
    lower=(6-s.Rational(80,100))*100**4
    negative=6*s.Rational(100,99)**4+80*s.Rational(100,99)**3
    assert lower>negative
    return dict(actual_H1_H2_derivative_identities=True,
                convexity_full_0_to_x_one_percent_proved_by_positive_rational_bound=True)


@current.source_precision
def run(owner=None, view=None):
    began=time.monotonic()
    if owner is None or view is None:
        owner,view=current.run()
    owner.assert_source()
    checks=_derivative_and_convexity_proof()
    # A correlated Rc/Pstar/Tw elimination from the actual source expression.
    kappa=owner.primitives.kappa
    symbols={str(q):q for q in kappa.free_symbols}
    Tw,logP=symbols['Tw'],symbols['logP']
    Rc=s.Symbol('Rc',positive=True)
    F=s.cancel(kappa*s.exp(Tw+logP+12)/4)
    physical_K=s.cancel(Rc*s.exp(Tw-2)*s.exp(logP)*kappa)
    assert s.simplify(physical_K-4*Rc*s.exp(-14)*F)==0
    # Derive the two source-scaled coefficients from the physical primitive.
    L=s.Symbol('L',positive=True); z=s.Symbol('Z',real=True)
    packet=owner.primitives.primitive(4*Rc*s.exp(-14)*F,2*Rc,2*Rc*s.exp(L))
    assert s.simplify(packet['m13']/Rc-2*s.exp(-28)*F**2*(1-s.exp(-L)))==0
    assert s.simplify(packet['m15']/Rc**2+4*s.exp(-28)*F**2*L)==0
    assert s.diff(physical_K,z)==0
    A,F_reference=independent_kernel_reference()
    c=current.MPIntervalContext();c.dps=140
    logA=current.read_interval(c,view['kernel']['log_kernel'])
    logF=current.read_interval(c,view['kernel']['log_F'])
    assert current.ep(logA)[0]<A<current.ep(logA)[1]
    assert current.ep(logF)[0]<F_reference<current.ep(logF)[1]
    actualL=current.read_interval(c,view['original_quiet_length'])
    log13=current.read_interval(c,view['m13']['log_coefficient'])
    log15=current.read_interval(c,view['m15']['log_absolute_coefficient'])
    with mp.workdps(90):
        l=(current.ep(actualL)[0]+current.ep(actualL)[1])/2
        ref13=mp.log(2)-28+2*F_reference+mp.log(-mp.expm1(-l))
        ref15=mp.log(4)-28+2*F_reference+mp.log(l)
    assert current.ep(log13)[0]<ref13<current.ep(log13)[1]
    assert current.ep(log15)[0]<ref15<current.ep(log15)[1]
    assert view['m13']['sign']=='positive' and view['m15']['sign']=='negative'
    assert view['m13']['radial_scale_power']==1 and view['m15']['radial_scale_power']==2
    assert view['normalized_quiet_coefficient_relative_width_passed']
    assert current.ep(current.read_interval(c,view['normalized_quiet_coefficient_relative_width']))[1]<mp.mpf('.001')
    assert view['source_scaled_Z1_exact_zero']
    assert all(current.ep(current.read_interval(c,q))[0]>0 for q in view['kernel']['normalized_positive_tails'])
    assert all(not view[k] for k in current.OPEN)
    checks.update(actual_Pstar_Rp_Tw_correlation_cancelled_before_evaluation=True,
        actual_source_physical_K_and_Rc_Rc_squared_pressure_scales_verified=True,
        production_kernel_encloses_independent_original_sigma_quadrature=True,
        production_signed_pressure_coefficients_enclose_independent_physical_primitives=True,
        both_normalized_signed_quiet_coefficients_relative_width_below_0_001=True,
        original_positive_mean_memory_sigma_denominator_and_integral_tails_retained=True,
        original_absolute_Rc_graph_node_retained_without_radius_materialization=True,
        quiet_contribution_only_full_tails_and_recursion_not_promoted=True)
    hashes=dict(owner.hashes)
    hashes[current.NAME]=current.sha(current.NAME)
    hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,separators=(',',':'))+'\n').encode())
    print('PASS_QUIET_ACTUAL_SOURCE_PRESSURE',view['kernel']['cells'],flush=True)
    return receipt


if __name__=='__main__':run()
