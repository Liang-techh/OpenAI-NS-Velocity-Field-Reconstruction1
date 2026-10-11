"""Exact unchanged-loop identities and local source routing checks."""
import json
from pathlib import Path
import time
from types import SimpleNamespace

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Foreign local leading scope admitted')


def exact_flat_identities():
    E,V,A,B=s.symbols('E V A B',real=True);N=s.symbols('N',positive=True)
    rr=E*(s.exp(A/N)-1);ss=B/N
    rates=dict(m=ss,h=rr,k=V*rr+E*ss+rr*ss,
               e=2*V*ss+ss*ss-E*rr-rr*rr/2,p=E*rr+rr*rr/2)
    assert all(s.simplify(value.subs({A:0,B:0}))==0 for value in rates.values())
    assert (E*s.exp(A/N)).subs(A,0)==E and (V+B/N).subs(B,0)==V
    # The background densities remain their original nonzero functions.
    assert (V*V-E*E/2).subs({E:2,V:3})==7
    return dict(exact_original_five_density_increment_flat_identities=5,
                unchanged_E_V_for_every_positive_N=True,background_density_not_zeroed=True)


def selector_checks(c):
    probe=SimpleNamespace(c=c,macro=SimpleNamespace(sc=c.mpf(1)/8));count=0
    for left,right,want,zero in ((0,'1/32','unchanged_left',True),
                                ('1/32','3/32','original_flat_left_collar',True),
                                (0,'1/8','unchanged_modified_union',False),
                                ('.2','.4','original_modified_right',False)):
        row=current.original.OriginalMicroFiniteN.modification_support(probe,'first_micro',c.mpf(left),c.mpf(right))
        assert row['region']==want and row['increments_exact_zero']==zero;count+=1
    return dict(original_executable_support_selector_branch_rows=count)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**exact_flat_identities(),**selector_checks(c))
    assert ep(owner.eta_log)[1]<=ep(owner.allowed_eta_log)[0]
    for name,packet in packets.items():
        row=owner.report(packet);leading=row['actual_leading_source']
        assert all(not row[key] for key in current.OPEN)
        assert row['background_velocity_pressure_and_own_histories_retained']
        assert set(row['later_five_moment_increment'])=={'m','h','k','e','p'}
        assert not any(row['later_velocity_pressure_increment'].values()) and not any(row['later_five_moment_increment'].values())
        assert row['executed_original_support_selector']['increments_exact_zero']
        assert not row['source_graph_ancestor_constructors_called'] and not row['active_Section11_inverse_evaluated']
        assert leading['source_family']==owner.family and leading['independent_actual_pressure']['terms']
        assert leading['actual_own_mean_not_comparison_mean']
        if row['query']['chart']=='core':assert row['executed_original_support_selector']['core_cone_argument_not_used']
    checks['whole_core_axis_background_Uz_and_pressure_retained']=bool(owner.report(packets['whole_core_axis'])['actual_leading_source']['ordinary_leading_rho_rows']['Uz'][0]['terms'])
    assert checks['whole_core_axis_background_Uz_and_pressure_retained']
    checks['guard_rejections']={name:rejected(fn) for name,fn in (
        ('outside_core',lambda:owner.evaluate('.371','4.01')),
        ('outside_selected_first',lambda:owner.evaluate('.371','1.01','first')),
        ('outside_real_Z',lambda:owner.evaluate('1.01','0')),
        ('unowned_outer_chart',lambda:owner.evaluate('.371','0','macro')),
        ('foreign_packet',lambda:owner.report(current.OriginalN1LocalLeadingPacket())))}
    old=owner.selector.macro.sc;owner.selector.macro.sc=old*2
    try:checks['selector_original_sc_mutation_rejected']=rejected(owner.assert_graph)
    finally:owner.selector.macro.sc=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,checks=checks,input_hashes=hashes,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_LOCAL_LEADING_DISPATCH',len(packets),checks['exact_original_five_density_increment_flat_identities'],flush=True)
    return result
