"""Actual current inlet/selected/flatten object and defining-equation checks."""
import json
import copy
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_Rp_selected_pulse as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import CompliantFlattenMixedC4


def overlap(a,b,label):
    al,ah=current.inlet.endpoints(a);bl,bh=current.inlet.endpoints(b)
    if max(al,bl)>min(ah,bh):raise ArithmeticError('Current consumer diagnostic disjoint: '+label)


def rows(left,right,label):
    assert left.order==right.order==5
    for n in range(6):overlap(left[n],right[n],label+'/'+str(n))
    return 6


def reject_stale_owner_links(owner):
    rejected=[]
    for label in ('old_pulse_fifth','old_high_constants','old_flatten_pulse'):
        candidate=copy.copy(owner);candidate.pulse=copy.copy(owner.pulse)
        candidate.flatten=copy.copy(owner.flatten)
        candidate.dispatch=candidate
        if label=='old_pulse_fifth':candidate.pulse.fifth=owner.seed.fifth
        elif label=='old_high_constants':
            candidate.pulse.high=copy.copy(owner.pulse.high)
            candidate.pulse.high.constants=owner.seed.axial4.constants
        else:candidate.flatten.pulse=owner.seed.pulse
        try:candidate.assert_graph()
        except ValueError:rejected.append(label)
        else:raise AssertionError('Stale native owner link accepted: '+label)
    return rejected


@source_precision
def run(owner=None):
    began=time.monotonic();raw=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=owner if owner is not None else current.CurrentOriginalRpSelectedPulse(require_checked=False)
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.PENDING)
    assert raw['source_family']==owner.family_record and all(owner.assert_graph().values())
    assert raw['actual_same_object_graph']==owner.assert_graph()
    rejected=reject_stale_owner_links(owner)
    admitted=json.loads((current.HERE/current.selected.RECEIPT).read_bytes())
    assert owner.current_selection_source_proof==admitted['current_future_selected_source_proof']
    assert owner.current_selection_source_proof['passed']
    original=current.selected
    methods=dict(actual_C4_incoming=owner.axial4.incoming.__func__ is original.CompliantAxialHighJets.incoming,
        actual_C4_selected=owner.axial4.select.__func__ is original.CompliantAxialHighJets.select,
        actual_C5_selected=owner.fifth.select.__func__ is original.CompliantFifthAxialJets.select,
        actual_C5_future=owner.fifth.future.__func__ is original.CompliantFifthAxialJets.future,
        actual_pulse_inlet_consumer=owner.pulse.data.__func__ is CompliantPulseRadialC4.data,
        actual_mixed_pulse=owner.pulse._high_packet.__func__ is original.CompliantPulseMixedC4._high_packet,
        actual_flatten_equations=owner.flatten.flatten.__func__ is CompliantFlattenMixedC4.flatten)
    assert all(methods.values())
    assert owner.flatten_binding['passed']
    assert owner.flatten_binding['retained_canonical_identity_counts']==dict(functional=155,inlet_datum=34,endpoint=24)
    assert not owner.flatten_binding['interval_overlap_used_as_join_proof']
    assert not owner.flatten_binding['source_caps_used_as_defining_field_values']
    assert type(owner.pulse.flat) is current.FlatPulseDerivatives
    assert owner.pulse.flat.beta.__func__ is current.FlatPulseDerivatives.beta
    beta=owner.pulse.flat.beta(0)
    assert beta.ctx is owner.ctx and beta.order==4 and current.inlet.endpoints(beta[0])[0]>0
    assert owner.fifth.angular4 is owner.energy4.angular is owner.angular4
    assert owner.angular4.cache is not owner.seed.fifth.angular4.cache
    c=owner.ctx;Z=c.mpf('.371')
    selected,inlet,u,incoming_energy,incoming_moments=owner.pulse.data(Z)
    assert selected['incoming']['Z_independent_constant_definitions'] is owner.constants
    assert selected['incoming']['ordinary_Taylor_order']==5
    assert u.ctx is incoming_energy.ctx is incoming_moments[0].ctx is incoming_moments[1].ctx is c
    assert all(jet.ctx is c for jet in selected['selected_scaled_end_coefficient_Taylor'])
    active=owner.evaluate('pulse_end',Z,-3)['source_packet']
    assert active['Uz_over_Utheta'].ctx is c
    ap=selected['selected_ap_Taylor'];quad=selected['quadratic_coefficients']
    residual=ap*ap*quad['A2']+ap*quad['A1']+quad['A0']
    for n in range(6):
        lo,hi=current.inlet.endpoints(residual[n]);assert lo<=0<=hi,('selected quadratic',n)
    assert current.inlet.endpoints(selected['positive_root_derivative_denominator'])[0]>0
    controls=selected['selected_scaled_end_coefficient_Taylor']
    assert current.inlet.endpoints(controls[0][0])[1]<0<current.inlet.endpoints(controls[1][0])[0]
    entrance=owner.evaluate('pulse_entrance',Z,0)['source_packet']
    expected=owner.inlet.current_power(Z,0)['canonical_incoming_Taylor']
    count=0
    for key,common in (
        ('Mz_over_R_Utheta','m1'),('Mtheta_z_over_sqrt2_R_3half_Utheta_squared','m2'),
        ('Mtheta_over_sqrt2_R_3half_Utheta','X'),('Mztheta_over_R_Utheta_squared','energy')):
        count+=rows(entrance[key],expected[common],'Rp/'+common)
    for key,common in (('P0_over_Pstar_squared','P0'),('Mp_over_Pstar_squared','Mp'),('P_over_Pstar_squared','pressure')):
        count+=rows(entrance['pressure'][key],expected[common],'Rp/'+common)
    assert count==42
    assert owner.flatten.inlet.datum is owner.datum is owner.future.angular.initial.datum
    terminal=owner.evaluate('pulse_end',Z,0)['source_packet']
    flat=owner.evaluate('flatten',Z,0)['source_packet']
    joined=0
    # Flatten publishes angular/energy under its own canonical names. Its
    # linear moments vanish by the retained terminal theorem, rather than
    # by fabricated fields added to the returned packet.
    for key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
        assert terminal[key].order==5
        for n in range(6):
            lo,hi=current.inlet.endpoints(terminal[key][n]);assert lo<=0<=hi
            joined+=1
    for key,target in (('Mtheta_over_sqrt2_R_3half_Utheta','angular_Taylor'),
                       ('Mztheta_over_R_Utheta_squared','energy_Taylor')):
        joined+=rows(terminal[key],flat[target],'Rv/'+key)
    for key in ('P0_over_Pstar_squared','Mp_over_Pstar_squared','P_over_Pstar_squared'):
        joined+=rows(terminal['pressure'][key],flat['pressure'][key],'Rv/'+key)
    assert joined==42
    future=owner.fifth.future(Z)['complete_future_energy_Taylor']/2
    rows(terminal['Mztheta_over_R_Utheta_squared'],future,'terminal complete future/2')
    assert all(owner.assert_graph().values())
    for chart,(method,domain) in current.EXPECTED.items():
        assert owner.provider(chart) is owner.pulse
        assert owner.chain_routes[chart]['method']==method and owner.chain_routes[chart]['domain']==domain
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        actual_rebound_current_object_graph=owner.assert_graph(),unchanged_actual_native_methods=methods,
        stale_native_owner_mutations_rejected=rejected,
        active_beta_center_and_end_bump_use_current_context=True,
        current_angular_C4_cache_is_separate=True,shared_future_and_angle_are_cache_free=True,
        original_selected_defining_equations_and_AST_proof=owner.current_selection_source_proof,
        actual_same_selected_pulse_flatten_binding=owner.flatten_binding,
        actual_fresh_Z='.371',actual_Rp_entrance_Taylor_diagnostic_rows=count,
        actual_Rv_flatten_join_Taylor_diagnostic_rows=joined,
        terminal_linear_zero_diagnostics_from_retained_exact_theorem=12,
        actual_selected_C5_quadratic_residual_contains_zero=True,
        positive_selected_root_and_end_control_signs_retained=True,
        terminal_energy_uses_same_current_complete_future_half=True,
        six_actual_routes_share_one_current_pulse_owner=True,
        current_selected_P0_object_used_not_replaced_by_saved_datum=True,
        source_function_ownership_not_inferred_from_interval_overlap=True,
        scope='Restricted actual selected pulse and flatten; full chart mixed contracts/postpulse physical assembly open',
        **dict.fromkeys(current.PENDING,False),input_hashes={**owner.hashes,
            current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.inlet.encode(current.inlet.pack(result)),indent=2)+'\n',
        encoding='utf8',newline='\n')
    print('Actual current selected pulse/flatten: input ownership, selected C5 and Rp/Rv joins passed',flush=True)
    return result


if __name__=='__main__':run()
