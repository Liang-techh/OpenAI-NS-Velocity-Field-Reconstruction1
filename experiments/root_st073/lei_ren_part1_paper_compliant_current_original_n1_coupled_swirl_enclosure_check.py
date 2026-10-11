"""Full-tail consumption, physical derivative gauge and actual width checks."""
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_coupled_swirl_enclosure as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported coupled swirl enclosure admitted')


def physical_error_fixture(c):
    """Independent error Jacobian W=S What, S_Z=-Lambda*Gprime*S."""
    algebra=op.N1LogAlgebra(c,(c.mpf(0),c.ln(5),c.mpf(0),c.ln(7)),[])
    gradient=op.Jet(c,[c.mpf(1)/3,c.mpf(0)])
    e0,e1=c.mpf(2)/11,c.mpf(3)/13
    error=algebra.lift(op.Jet(c,[e0,e1]))
    factor=current.weighted.absolute_jet(current.weighted.relative_amplitude(algebra,gradient,1,1))
    physical=algebra.shift(factor*error,(0,1,0,1))
    expected=(35*e0,35*(e1+7*e0/3));count=0
    for k,want in enumerate(expected):
        got=sum((v[k]*5**p[1]*7**p[3] for p,v in physical.terms.items()),c.mpf(0))
        assert ep(got)[0]<=ep(want)[0] and ep(got)[1]>=ep(want)[1];count+=1
    return dict(independent_physical_error_Jacobian_rows=count,
                ordinary_Z_error_contains_S_Z_times_value_error_once=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=physical_error_fixture(c);tail=owner.bound['all_word_feedback_tail']
    assert owner.bound[current.feedback.GATE]
    assert tail['retained_indices']==[0,0] and tail['covers_all_derivative_words']
    for packet in packets.values():
        v=owner.report(packet);assert all(not v[k] for k in current.OPEN)
        assert v['full_six_row_remaining_feedback_uniform_bound_included'] and v['unknown_feedback_not_set_to_zero']
        assert v['ordinary_Z1_tail_tube_radius']==tail['ordinary_Z1_tube_radius']
        assert v['only_F1_and_d_x_F1_components_enclosed_here'] and v['solved_coefficient_function_contracts_still_open']
        assert v['axial_pressure_mean_Q_accuracy_still_open'] and v['n1_own_repair_and_n2_still_open']
        assert v['physical_F0_Lambda_Z_derivative_gauge_added_once']
        assert v['passing_normalized_coupled_width_rows']==4 and v['passing_physical_coupled_width_rows']==4
        state=owner._packets[id(packet)][3];error=state[5]
        assert set(error.terms)=={(0,0,0,0)}
        assert error.terms[(0,0,0,0)][0]==tail['value_remainder_upper']
        assert error.terms[(0,0,0,0)][1]==tail['ordinary_Z1_remainder_upper']
        assert all(ep(e)[0]>0 for e in error.terms[(0,0,0,0)].coefficients)
        for field in current.FIELDS:
            for k in (0,1):
                row=v['normalized_coupled_swirl_budgets'][field]['Z'+str(k)]
                physical=v['physical_coupled_swirl_budgets'][field]['Z'+str(k)]
                assert row['remaining_B0_B1_feedback_included'] and physical['remaining_B0_B1_feedback_included']
                assert row['relative_coupled_swirl_width_certified'] and physical['relative_coupled_swirl_width_certified']
                assert ep(row['total_width_over_abs_lower'])[1]<=c.mpf('.001')
                assert ep(physical['total_width_over_abs_lower'])[1]<=c.mpf('.001')
                assert any(t['source_exponents']==[0,0,0,0] for t in row['error_source_sector_contributions'])
    checks['actual_first_inlet_and_corrected_exit_four_coupled_swirl_width_rows_pass']=True
    checks['full_all_word_value_Z1_tail_consumed_without_setting_feedback_to_zero']=True
    checks['physical_and_normalized_budgets_include_remaining_feedback']=True
    checks['axial_pressure_repair_full_solution_scope_remains_open']=True
    checks['unsupported_endpoint_or_foreign_packet_rejected']=all((
        rejected(lambda:owner.evaluate('.371',c.mpf([0,1]))),
        rejected(lambda:owner.report(current.OriginalN1CoupledSwirlEnclosurePacket()))))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])]
    old=entry[1]['axial_pressure_mean_Q_accuracy_still_open']
    entry[1]['axial_pressure_mean_Q_accuracy_still_open']=False
    try:checks['false_complete_velocity_scope_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['axial_pressure_mean_Q_accuracy_still_open']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf8',newline='\n')
    print('PASS_N1_COUPLED_SWIRL_ENCLOSURE',flush=True)
    return result
