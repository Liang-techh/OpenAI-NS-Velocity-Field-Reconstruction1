"""Independent homogeneous transport fixture and actual collar source checks."""
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_correlated_inner_source as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported correlated source admitted')


def homogeneous_fixture(c):
    """Exact chi=0 solution and chi-variation, independent of the gain formula."""
    algebra=op.N1LogAlgebra(c,(c.mpf(0),)*4,[]);count=0
    for a,b,az,bz in ((1,2,3,-1),(-3,1,-2,4)):
        vals={'F1':algebra.lift(op.Jet(c,[a,az])),'d_x_F1':algebra.lift(op.Jet(c,[b,bz]))}
        errors={field:algebra.lift(0,1) for field in current.FIELDS}
        norm=current.vector_norm(vals,errors)
        for xvalue in ('2','2.001','2.02'):
            x=c.mpf(xvalue);cz=c.mpf(1)/3
            got,_=current.homogeneous_transfer_error(c,norm,cz,x)
            if xvalue=='2':assert not got.terms;continue
            bound=got.terms[(0,0,0,0)]
            y=c.mpf(a)+b*(1-4/x**2);dy=8*b/x**3
            # Differentiate the exact chi=0 Cauchy ODE in chi.
            q=-(a+b)*(x*x-4)/4+4*b*c.ln(x/2)+4*(a-b)*(c.mpf(1)/4-1/x**2)
            dq=-(a+b)*x/2+4*b/x+8*(a-b)/x**3
            yz=c.mpf(az)+bz*(1-4/x**2)+cz*q;dyz=8*bz/x**3+cz*dq
            for k,changes in ((0,(y-a,dy-b)),(1,(yz-az,dyz-bz))):
                for value in changes:
                    assert ep(abs(value))[1]<=ep(bound[k])[1];count+=1
    return dict(independent_exact_homogeneous_and_parameter_variation_bound_rows=count,
                transfer_inlet_gain_exact_zero=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c;checks=homogeneous_fixture(c)
    old=json.loads(gzip.decompress((current.HERE/current.weighted.NAME).read_bytes()))
    old_exit=old['samples']['selected_actual_inner_exit']
    assert old['source_family']==op.snapshot(owner.family)
    assert all(row['prefix']['prefix_width_over_abs_lower'] is None for field in old_exit['normalized_kernel_sector_budgets'].values() for row in field.values())
    checks['previous_accepted_whole_cell_swirl_prefix_four_widths_uncertified']=True
    for name,packet in packets.items():
        view=owner.report(packet);entry=owner._packets[id(packet)]
        assert all(not view[k] for k in current.OPEN)
        assert view['initial_source_integral_only'] and view['remaining_normalized_B0_B1_feedback_not_solved']
        assert view['original_actual_collar_not_comparison_or_Bessel_model']
        assert view['passing_normalized_initial_width_rows']==4
        for label in ('normalized_initial_integral_budgets','physical_initial_integral_budgets'):
            assert all(row['relative_initial_integral_width_certified'] and not row['remaining_B0_B1_feedback_included']
                       for field in view[label].values() for row in field.values())
        for value in (*entry[6].values(),*entry[7].values()):
            assert all(p[0]>=0 and p[1]==0 and p[2]==0 and p[3]<=0 for p in value.terms)
        if name=='actual_first_inlet':
            assert not view['actual_collar_source_cells']
            assert not view['normalized_homogeneous_core_to_exit_error']['terms']
        else:
            assert len(view['actual_collar_source_cells'])==4
            assert view['query']['actual_source_coordinate']==owner.local.inner.sc/4
            assert ep(view['homogeneous_transfer_proof']['length_upper'])[1]>0
            for cell in view['actual_collar_source_cells']:
                assert cell['actual_normalized_Phi_used_before_F0_dressing'] and cell['inverse_hb_cancelled_before_enclosure']
    checks['actual_original_corrected_exit_four_initial_width_targets_pass']=True
    checks['physical_ordinary_Z1_error_transfer_and_actual_four_collar_cells_retained']=True
    checks['original_fixed_phase_geometry_and_all_inverse_hb_cancelled']=True
    checks['core_source_transport_and_collar_tail_separate_from_coupled_feedback']=True
    checks['unsupported_fraction_partition_or_foreign_packet_rejected']=all((
        rejected(lambda:owner.evaluate('.371','1.1')),
        rejected(lambda:owner.evaluate('.371',c.mpf([0,1]))),
        rejected(lambda:owner.evaluate('.371','1',True)),
        rejected(lambda:owner.report(current.OriginalN1CorrelatedInnerSourcePacket()))))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])];oldflag=entry[1]['initial_source_integral_only']
    entry[1]['initial_source_integral_only']=False
    try:checks['false_coupled_scope_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['initial_source_integral_only']=oldflag
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_CORRELATED_INNER_SOURCE',checks['independent_exact_homogeneous_and_parameter_variation_bound_rows'],flush=True)
    return result
