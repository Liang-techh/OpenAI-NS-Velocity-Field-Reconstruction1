"""Focused exact common-inner geometry and delegated field checks."""
import json
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Invalid common-inner scope admitted')


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c;view=owner.geometry()
    h,sc=s.symbols('hb sc',positive=True);a,p=s.symbols('logC logP',positive=True)
    keep=h*sc/8;inner=h*sc/4;left=h*sc/2
    assert s.simplify(inner-keep-h*sc/8)==0 and s.simplify(left-inner-h*sc/4)==0
    rm=10*a+5*(2*p)+s.log(110)-6;ra=s.log(4)-4*p-1000
    assert s.expand(rm-ra-inner-(10*a+14*p+s.log(s.Rational(110,4))+994-inner)).equals(0)
    assert 0<ep(owner.phase_keep)[0]<ep(owner.phase_in)[0]<ep(owner.phase_minus)[0]<ep(c.mpf('.5'))[0]
    assert view['selected_sc_exact_mpf_tuple']==list(ep(owner.sc)[0]._mpf_)
    assert ep(view['log_Rm_over_Rin_strict_lower'])[0]>0
    for packet in packets.values():
        row=owner.report(packet);field=row['canonical_n1_field']
        assert field[current.coupled.GATE]
        assert all(not row[key] and not field[key] for key in current.OPEN)
        assert row['query']['actual_source_coordinate']==owner.phase_in
        assert field['ordinary_Z_prefix_derivative_order']==1 and field['initial_Z_derivative_order']==2
        assert field['radial_recovery_uses_temporal_q1_equals_1_plus_delta']
        assert field['canonical_inner_series_values_include_remaining_uniform_absolute_tail']
        assert set(field['n1_radial_recovery_series_enclosures'])=={'mean1','Q1','V1'}
    guards={name:rejected(fn) for name,fn in (
        ('outside_corrected_fraction',lambda:owner.evaluate('.371','1.01','first')),
        ('negative_corrected_fraction',lambda:owner.evaluate('.371','-.01','first')),
        ('outside_actual_core',lambda:owner.evaluate('.371','4.01')),
        ('foreign_chart',lambda:owner.evaluate('.371','0','macro')),
        ('foreign_packet',lambda:owner.report(current.OriginalN1CommonInnerPacket())))}
    original=owner.sc;owner.sc=original*2
    try:guards['selected_original_sc_mutation']=rejected(owner.assert_graph)
    finally:owner.sc=original
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,
        checks=dict(exact_positive_common_inner_gap_identities=True,
            original_Rm_second_basis_logP_factor_checked=True,original_sc_exact_binary_endpoint_retained=True,
            accepted_coupled_field_computed_on_corrected_common_inner_exit=True,
            completed_leading_dispatch_not_promoted=True,guards_rejected=guards),
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_CORRECTED_COMMON_INNER_SELECTION',len(packets),len(guards),flush=True)
    return result
