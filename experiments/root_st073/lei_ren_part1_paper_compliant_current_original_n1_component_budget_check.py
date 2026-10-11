"""Independent finite-sector normalization and radial tail bounds."""
import itertools
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_component_budget as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def contains(got,want):return ep(got)[0]<=ep(want)[0] and ep(got)[1]>=ep(want)[1]


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported component budget admitted')


def sector_fixture(c):
    bases=(3,5,7,11);algebra=op.N1LogAlgebra(c,tuple(c.ln(x) for x in bases),[])
    a=algebra.shift(algebra.lift(op.Jet(c,[c.mpf([9,11])/10,2])),(0,1,0,0))
    b=algebra.shift(algebra.lift(op.Jet(c,[-c.mpf([1,2])/10,-c.mpf(1)/3])),(0,0,0,-1))
    d=algebra.shift(algebra.lift(op.Jet(c,[c.mpf(1)/2,0])),(1,0,1,0));value=a+b+d;count=0
    precise=type(c)();precise.dps=c.dps*2
    for k in (0,1):
        row=current.component_budget(value,k,c.mpf(1)/100)
        p=row['reference_source_exponents'];M=precise.mpf(row['directed_reference_coefficient_norm'])
        unit=M
        for base,power in zip(bases,p):unit*=precise.mpf(base)**power
        direct=5*precise.mpf([9,11])/10-precise.mpf([1,2])/110+precise.mpf(21)/2 if k==0 else precise.mpf(10)-precise.mpf(1)/33
        assert contains(row['signed_prefix_normalized_enclosure'],direct/unit)
        assert contains(row['full_series_normalized_enclosure'],(direct+precise.mpf([-1,1])/100)/unit);count+=2
    # The global scalar context is intentionally below the owner context.
    # A negative non-dyadic endpoint must not round through scalar abs.
    negative=algebra.lift(op.Jet(c,[-c.mpf(1)/3]))
    row=current.component_budget(negative,0,c.mpf(0))
    assert contains(row['signed_prefix_normalized_enclosure'],
        -precise.mpf(1)/3/precise.mpf(row['directed_reference_coefficient_norm']));count+=1
    huge=op.N1LogAlgebra(c,(c.mpf(0),-c.mpf(10)**1000,c.mpf(0),c.mpf(0)),[])
    value=huge.shift(huge.lift(c.mpf([98,102])/100,1),(0,1,0,0))+huge.shift(huge.lift(1,1),(0,3,0,0))
    row=current.component_budget(value,0,c.mpf('1e-900'))
    assert row['reference_source_exponents']==[0,1,0,0]
    assert row['remaining_tail_normalized_norm_upper'] is None
    assert row['remaining_tail_bound_rule']=='positive_relative_scale_retained_logarithmically'
    assert not row['relative_width_certified'] and not row['source_amplitude_materialized']
    assert any('directed_small_log_cover' in part['endpoint_bound_rules'] for part in row['source_sector_contributions'])
    assert rejected(lambda:current.component_budget(value,2,c.mpf(1)))
    zero=current.component_budget(huge.lift(0,1),0,c.mpf(0));assert zero['total_exact_zero']
    return dict(independent_rational_sector_sum_and_tail_containment_rows=count,
        directed_negative_endpoint_magnitude_kept_at_owner_precision=True,
        huge_equal_source_logs_cancel_by_exponent_difference=True,
        unweighted_tail_swamping_F0_sector_detected_without_exponentiation=True,
        unsupported_derivative_rejected_and_exact_zero_trace_retained=True)


def radial_fixture(c):
    count=0
    for z in ('-.4','0','.371'):
        z=c.mpf(z);delta=c.mpf(1)/5;E=c.mpf(1)/100;E1=c.mpf(1)/20
        bound=((abs(1-delta)+abs(1+delta))*abs(z)*E+2*abs(1-z*z)*E1)/(1-delta*z*z)
        for signs in itertools.product((-1,1),repeat=4):
            u,k,uz,kz=(c.mpf(sign)*norm for sign,norm in zip(signs,(E,E,E1,E1)))
            q=((1-delta)*z*u-(1+delta)*z*k-(1-z*z)*(uz+kz))/(1-delta*z*z)
            assert ep(abs(q))[1]<=ep(bound)[1];count+=1
    return dict(independent_real_Q1_tail_extreme_error_rows=count)


@source_precision
def run(owner,packet):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**sector_fixture(c),**radial_fixture(c));view=owner.report(packet)
    assert all(not view[key] for key in current.OPEN)
    rows=[row for field in view['component_budgets'].values() for row in field.values()]
    assert len(rows)==view['total_component_derivative_rows']==16
    assert sum(row['relative_width_certified'] for row in rows)==view['passing_relative_total_width_rows']
    assert view['source_uncertainty_and_integration_error_still_combined'] and view['sector_preserving_analytic_tail_still_open']
    assert all(not row['source_amplitude_materialized'] for row in rows)
    assert view['component_budgets']['V1']['Z0']['remaining_tail_source_exponents']==[0,0,0,-1]
    checks['physical_tail_Lambda_factor_retained']=True
    checks['foreign_common_inner_packet_rejected']=rejected(lambda:owner.evaluate(current.selection.OriginalN1CommonInnerPacket()))
    old=owner._packets[id(packet)][1]['passing_relative_total_width_rows']
    owner._packets[id(packet)][1]['passing_relative_total_width_rows']=old+1
    try:checks['live_budget_mutation_rejected']=rejected(lambda:owner.report(packet))
    finally:owner._packets[id(packet)][1]['passing_relative_total_width_rows']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        passing_relative_total_width_rows=view['passing_relative_total_width_rows'],total_component_derivative_rows=16,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_COMPONENT_BUDGET',checks['independent_rational_sector_sum_and_tail_containment_rows'],checks['independent_real_Q1_tail_extreme_error_rows'],flush=True)
    return result
