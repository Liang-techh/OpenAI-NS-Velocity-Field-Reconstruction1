"""Focused mathematical checks for actual-collar bounds and tail admission."""
import copy
import json
import math
from pathlib import Path
import time

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:fn()
    except (ValueError,TypeError):return True
    return False


def exact_collar_identities():
    """Directly differentiate actual ODEs, without the positive recurrence."""
    rho=s.symbols('rho',positive=True)
    A=s.Function('A')(rho);chi=s.Function('chi')(rho)
    phi=s.Function('phi')(rho);r=s.Function('ratio')(rho);w=s.Function('w')(rho)
    core_w=s.Function('core_w')(rho)
    rules={s.diff(phi,rho):chi*A*phi,s.diff(r,rho):(chi-1)*A*r,
           s.diff(w,rho):chi*r*s.diff(core_w,rho)}
    phi_second=s.diff(rules[s.diff(phi,rho)],rho).subs(rules)
    w_second=s.diff(rules[s.diff(w,rho)],rho).subs(rules)
    assert s.expand(phi_second-phi*(s.diff(chi,rho)*A+chi**2*A**2+chi*s.diff(A,rho)))==0
    assert s.expand(w_second-r*(s.diff(chi,rho)*s.diff(core_w,rho)
                  +chi*(chi-1)*A*s.diff(core_w,rho)+chi*s.diff(core_w,rho,2)))==0
    hb,eps,rho0=s.symbols('hb epsilon_b rho0',positive=True)
    sig=s.Function('sigma');phase=s.log(rho/rho0)/hb
    cutoff=1-(1-eps)*sig(phase)
    t=s.symbols('t')
    sigma_prime=s.Subs(s.diff(sig(t),t),t,phase)
    assert s.simplify(s.diff(cutoff,rho)+(1-eps)*sigma_prime/(hb*rho))==0
    # The source measure cancels an inverse width, never its numerical cap.
    x=s.symbols('x',positive=True);a,b=s.symbols('a b',nonnegative=True)
    assert s.expand((a/hb+b)*(x*hb/2)-(a+b*hb)*x/2)==0
    return 4


def tail_majorant_checks():
    """Independent high-precision scalar summation, including alpha crossings."""
    rows=0
    with mp.workdps(120):
        for alpha in map(mp.mpf,('0','1e-100','.001','.02','1','4')):
            # Check the derived two-step majorant across the max(1,alpha*m)
            # switch, independent of the interval implementation.
            for k in range(1,160):
                m=(k+1)//2
                a=max(mp.mpf(1),alpha*m)**m
                b=max(mp.mpf(1),alpha*(m+1))**(m+1)
                assert b/a<=mp.e*(1+alpha*(k+3)/2)
                rows+=1
        c=MPIntervalContext();c.dps=100
        for alpha,N in (('0',64),('.02',64),('1',128),('4',512)):
            result=current.neumann_tail(c,c.mpf(3),c.mpf(2),c.mpf(1),c.mpf(alpha),c.mpf(1),N)
            def direct(k):
                m=(k+1)//2
                return 3*mp.mpf(4)**k/mp.factorial(k)*max(mp.mpf(1),mp.mpf(alpha)*m)**m
            total=mp.fsum(direct(k) for k in range(N+1,N+1200))
            assert total<=current.ep(result['omitted_sum_bound'])[1]
            rows+=1
    return rows


@source_precision
def run(owner,packet):
    began=time.monotonic();view=owner.report(packet);c=owner.c
    identities=exact_collar_identities();tail_rows=tail_majorant_checks()
    scalar=lambda v:current.WidthBound.scalar(c,c.mpf(v))
    w=current.WidthBound(c,{-1:c.mpf(3),0:c.mpf(5)})
    controls=dict(uncancelled_inverse_width_not_materialized=rejected(lambda:w.numeric_upper(owner.source.cap)),
        shifted_measure_cancels_inverse_width=set(w.shift(1).rows)=={0,1},
        wrong_packet_identity=rejected(lambda:owner.report(copy.copy(packet))),
        wrong_owner_type=rejected(lambda:current.CurrentOriginalN1InnerAnalyticDomain(copy.copy(owner))),
        negative_bound=rejected(lambda:current.WidthBound(c,{0:c.mpf(-1)})),
        unsupported_tail_degree=rejected(lambda:current.neumann_tail(c,c.mpf(1),c.mpf(2),c.mpf(1),c.mpf(1),c.mpf(1),0)))
    # Check true original nozero/tube guard and every retained inverse width.
    geom=view['actual_inner_geometry'];tube=view['common_tubes'];den=view['denominator_bounds']
    assert current.ep(den['Phi_core_modulus_lower'])[0]>0
    assert current.ep(den['L_modulus_lower'])[0]>0
    assert current.ep(tube['solution_Z_radius'])[1]<=current.ep(tube['Gg_Z_radius']/2)[1]
    for group in view['actual_collar_regular_forcing_bounds'].values():assert min(map(int,group))>=-1
    for group in view['actual_collar_source_integral_bounds'].values():assert min(map(int,group))>=0
    assert geom['exact_positive_width_preserved'] and geom['actual_bridge_not_replaced_by_comparison']
    assert all(not view[key] for key in current.OPEN)
    # The actual matrix has no inverse width; its accepted core norm must
    # fit inside the new bound at one representative inner point.
    core_packet=owner.n1.system('2','.371');core=owner.n1.report(core_packet)
    Amax=view['original_axis_amplitude']['complex_modulus_upper']
    norms=[]
    for which,norm_key in (('B0','regular_B0_infinity_norm_upper'),('B1','regular_B1_infinity_norm_upper')):
        rows=[c.mpf(0) for _ in range(6)]
        for ij,sectors in core['regular_system'][which].items():
            i=int(ij.split(',')[0])-1
            rows[i]+=sum((abs(v)*Amax**int(p) for p,v in sectors.items()),c.mpf(0))
        assert max(current.ep(v)[1] for v in rows)<=current.ep(view[norm_key])[1]
        norms.append(which)
    # Snapshot guards must reject width, context and family replacement.
    saved=owner.source.logh
    try:
        owner.source.logh=saved+1
        controls['changed_exact_width']=rejected(lambda:owner.report(packet))
    finally:owner.source.logh=saved
    saved=owner.family
    try:
        owner.family={}
        controls['changed_family']=rejected(lambda:owner.report(packet))
    finally:owner.family=saved
    saved=owner.c
    try:
        owner.c=MPIntervalContext()
        controls['changed_context']=rejected(lambda:owner.report(packet))
    finally:owner.c=saved
    assert all(controls.values())
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,
        exact_differentiated_ODE_and_width_measure_identities=identities,
        independent_tail_ratio_and_sum_rows=tail_rows,
        accepted_core_matrix_norm_containment=norms,
        source_width_is_not_a_materialized_parameter=True,
        common_tube_nozero_guards_passed=True,
        integrated_source_has_no_negative_width_power=True,
        radial_bound_orders=[0,4],bounds_do_not_supply_computed_iteration_terms=True,
        invalid_inputs_rejected=controls,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.parent.parent.first.encode(
        current.parent.parent.original.serialized(receipt)),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_ACTUAL_COLLAR_N1_ANALYTIC_BOUND_AND_TAIL',identities,tail_rows,flush=True)
    return receipt
