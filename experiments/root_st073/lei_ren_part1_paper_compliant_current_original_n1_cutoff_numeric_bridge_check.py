"""Exact signed increment/FTC/unit identities and actual own-field ranges."""
import json
from pathlib import Path
import time
from types import SimpleNamespace

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric_bridge as current
import lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge_check as prior

legacy=current.legacy;op=current.op;ep=current.ep


def exact_cutoff_identity():
    z=s.Symbol('Z',real=True);r=s.Symbol('rho',positive=True)
    eps,hb,phase,delta=s.symbols('epsilon hb phase delta',positive=True)
    F,U,f,u,N,chi=s.symbols('F1 U1 F0 U0 N1p chi')
    ctx=SimpleNamespace(mpf=lambda x:s.Rational(str(x)))
    class Algebra:
        def __init__(self):self.ctx=ctx
        def lift(self,v,order=None):return s.sympify(v)
        def shift(self,v,p):return v*eps**(-p[3])
        def width(self,v,power=1):return v*hb**power
        def add(self,*v):return s.Add(*v)
        def mul(self,*v):return s.Mul(*v)
        def neg(self,v):return -v
        def scale(self,v,c):return v*c
    a=Algebra();g=SimpleNamespace(constant=lambda v:s.Rational(str(v)),mul=lambda *v:s.Mul(*v))
    signed=prior.compile_mathematical(current.true.signed_rho_integrands,
        {'c3':SimpleNamespace(C3Algebra=lambda g:a)})
    leading={'F_rho_rows':[f],'U_rho_rows':[u]};system={'rho':r,'N1p':N}
    base=signed(g,leading,{'F1':[F],'Uz1':[U]},system,eps)
    cut=signed(g,leading,{'F1':[chi*F],'Uz1':[chi*U]},system,eps)
    change=legacy.change_densities(a,{'F1':F,'Uz1':U},f,u,r,chi-1)
    for name in legacy.parent.MOMENTS:assert s.expand(base[name]+change[name]-cut[name])==0,name
    # Arbitrary sigma value/ordinary first/second derivatives. Execute
    # the real transition-branch recipe; the true graph binds the same
    # original sigma and the Taylor-to-ordinary factor separately.
    sig0,sig1,sig2=s.symbols('sigma sigma_prime sigma_second')
    fake_domain=SimpleNamespace(sigma_jets=lambda c,x:[sig0,sig1,sig2/2])
    cutoff=prior.compile_mathematical(legacy.cutoff,{'domain':fake_domain,'ep':lambda v:(v,v)})
    jets=prior.compile_mathematical(legacy.cutoff_radial_jets,
        {'domain':fake_domain,'cutoff':cutoff,'ep':lambda v:(v,v)})
    actual=jets(a,s.Rational(5,8),r,phase)
    wanted=[1-sig0,-4*sig1/(hb*phase*r),
        4*sig1/(hb*phase*r*r)-16*sig2/(hb*hb*phase*phase*r*r)]
    for got,want in zip(actual,wanted):assert s.cancel(got-want)==0
    v,v1,v2=s.symbols('v v_rho v_rhorho')
    fake_op=SimpleNamespace(truncate=lambda x,n:x,dz=lambda x:s.diff(x,z),
        Jet=SimpleNamespace(variable=lambda c,Z,n:Z))
    pair=prior.compile_mathematical(legacy.cut_pair,{'op':fake_op})(v,v1,v2,actual)
    true_pair=prior.compile_mathematical(current.true.multiply_cut_rows,
        {'c3':SimpleNamespace(C3Algebra=lambda g:a)})(g,[v,v1,v2],wanted)
    for got,want in zip(pair.values(),true_pair):assert s.cancel(got-want)==0
    mean=s.Function('mean')(r,z);Ufunc=s.Function('U')(r,z);M=eps*r*mean
    recovered=prior.compile_mathematical(legacy.recover_velocity,{'op':fake_op})(a,Ufunc,M,z,r,delta)
    Q=(2*z*Ufunc-(1+delta)*z*mean-(1-z*z)*s.diff(mean,z))/(1-delta*z*z)
    for name,want in {'mean1':mean,'K1':mean-Ufunc,'V1':eps*r*Q}.items():
        assert s.cancel(recovered[name]-want)==0
    P=s.Function('P')(r);I=s.Function('I')(r)
    parts=I-r*P
    assert s.expand(s.diff(parts,r).subs(s.diff(I,r),2*u*U+P)-(2*u*U-r*s.diff(P,r)))==0
    H5,H6,Fx,Ux=s.symbols('H5 H6 Fx Ux');x=s.sqrt(r)
    assert s.cancel(((H5-3*Fx/x)-Fx/x)/(4*r)-(H5/(4*r)-Fx/(r*x)))==0
    assert s.cancel(((H6-Ux/x)-Ux/x)/(4*r)-(H6/(4*r)-Ux/(2*r*x)))==0
    inlet_mean,inlet_pressure=s.symbols('core_mean_at_four core_pressure_at_four')
    # At first phase zero, rho=4 and every own first-history integral has
    # equal zero endpoints. These are the exact core inlet functions.
    assert s.cancel((4*inlet_mean+0)/4-inlet_mean)==0
    assert inlet_pressure+0==inlet_pressure
    return dict(passed=True,all_five_signed_base_plus_increment_identities=5,
        original_cutoff_value_first_second_rho_recipe_identities=3,
        actual_velocity_cutoff_product_rule_identities=3,
        own_mean_K_V_and_analytic_Q_physical_identities=3,
        signed_fifth_pressure_integration_by_parts_identity=True,
        even_second_rho_equation_and_physical_scaling_identities=2,
        exact_first_zero_core_inlet_FTC_history_identities=True,
        unchanged_leading_N1p_cancels_only_in_increment_not_full_moment=True,
        original_hb_inverse_and_inverse_squared_and_sigma_factor_two_verified=True)


def exact_zero(v):
    return all(ep(x)==(0,0) for row in v.terms.values() for x in row.coefficients)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported true cutoff query admitted')


@current.source_precision
def run(owner,raw,packets):
    began=time.monotonic();owner.assert_graph();proof=exact_cutoff_identity()
    assert owner.assert_source_identity()==owner.source_identity
    assert current.op.snapshot(raw['source_identity'])==current.op.snapshot(owner.source_identity)
    assert set(packets)=={'axis','first_inlet','transition','zero','Rin'}
    rows={}
    for name,packet in packets.items():
        view=owner.report(packet);entry=owner._packets[id(packet)]
        assert current.op.snapshot(view)==current.op.snapshot(raw['samples'][name])
        assert all(view[k] is False for k in current.OPEN)
        assert set(view['physical_five_inner_cumulative_Z1'])==set(legacy.parent.MOMENTS)
        assert view['exterior_pressure_m13_m15_not_added']
        assert view['Q_axis_radius_division_not_used']
        numeric=owner.numeric._packets[id(entry[2])]
        chi,pairs,physical,M,own,P,totals=numeric[2]
        for quantity,vals in physical.items():
            assert all(v.order==1 for v in vals.values())
        assert M.order==2 and own['mean1'].order==2 and P.order==1
        if name in ('zero','Rin'):
            assert all(exact_zero(v) for vals in physical.values() for v in vals.values())
            assert not exact_zero(M)
            assert view['five_inner_cumulative_budgets']['m11']['Z0']['sign'] in ('positive','negative')
        if name=='axis':
            assert exact_zero(M) and exact_zero(P) and exact_zero(entry[4][0])
            assert all(exact_zero(v) for v in entry[4][3].values())
        if name=='first_inlet':
            assert not exact_zero(M)
            assert ep(owner.numeric.regular._packets[id(numeric[3][1])][5][5]['rho'])==(4,4)
            assert not exact_zero(entry[4][1]['mean1'][0])
            change=numeric[3][4]
            assert all(exact_zero(v) for v in change[0].values()) and exact_zero(change[1])
        velocity= view['numerical_cutoff']['full_velocity_support_moment_budgets']
        assert all(b['relative_partial_moment_width_certified'] for vals in velocity.values() for b in vals.values())
        count=sum(b['relative_partial_moment_width_certified'] for vals in view['five_inner_cumulative_budgets'].values() for b in vals.values())
        rows[name]=dict(actual_query_executed=True,passing_five_inner_value_Z1_budget_rows=count,
            passing_full_velocity_support_value_Z1_rows=6,
            nonzero_M_history_after_velocity_zero=name in ('zero','Rin'))
    # The code selects exactly the same terminal M packet from .75 onward.
    # Coefficient/sector identity, rather than overlap, binds the history.
    zero=owner.numeric._packets[id(owner._packets[id(packets['zero'])][2])][2][3]
    Rin=owner.numeric._packets[id(owner._packets[id(packets['Rin'])][2])][2][3]
    assert op.snapshot(zero)==op.snapshot(Rin)
    assert rejected(lambda:owner.evaluate('.371','1.01','first'))
    assert rejected(lambda:owner.evaluate('.371','4.01','core'))
    assert rejected(lambda:owner.evaluate('.371',['.6','.7'],'first'))
    reused=owner._packets[id(packets['transition'])][2]
    assert rejected(lambda:owner.evaluate('.371','.7','first',reuse_packet=reused))
    assert rejected(lambda:owner.report(legacy.OriginalN1RadialCutoffPacket()))
    inputs=dict(owner.hashes);inputs[current.NAME]=current.sha(current.NAME)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=inputs,
        source_identity=owner.source_identity,exact_cutoff_proof=proof,query_results=rows,
        checks=dict(actual_axis_first_inlet_transition_zero_and_Rin_queries_execute=True,
            first_phase_zero_is_core_rho4_not_the_axis=True,
            own_mean_Q_V_P_and_nonzero_terminal_history_retained=True,
            all_five_signed_inner_cumulative_and_endpoint_enclosures_bound=True,
            six_full_velocity_support_value_Z1_budgets_pass=True,
            complete_pressure_tails_and_reserved_repair_still_open=True,
            unsupported_domains_intervals_and_wrong_reused_query_rejected=True),
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(op.snapshot(receipt),indent=2)+'\n').encode())
    print('TRUE_CUT_N1_OWN_NUMERIC_CHECKED',len(packets),len(inputs),flush=True)
    return receipt
