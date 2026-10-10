"""Focused eight current raw postpulse function/scale seam checks."""
import copy
import gzip
import json
from pathlib import Path
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_postpulse_mixed_seams as current
from lei_ren_part1_paper_compliant_current_original_Rp_pulse_mixed_seams_check import interpretation
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def exact_scales_and_diagnostics(owner,views):
    reader,mu=interpretation(owner)
    L,T,W=s.symbols('same_original_Lrel same_original_Ts same_true_original_waiting',positive=True)
    for key,value in (('Lrel',L),('Ts',T),('waiting',W)):
        reader.bindings[owner.radius.functions[key].node]=value
    J={'outer_power':L-4,'steep_power':T,'waiting':W}
    scales=native=count=0;differences={}
    ends=current.mixed.pulse.radius.post.selected.inlet.endpoints
    for name,view in views.items():
        left,right=view['left'],view['right'];differences[name]={}
        assert s.simplify(reader.at(left['geometry']['logR'])-reader.at(right['geometry']['logR']))==0,name
        JL=reader.at(left['geometry']['native_to_log_radius_jacobian'])
        JR=reader.at(right['geometry']['native_to_log_radius_jacobian'])
        assert s.simplify(JL-J.get(left['chart'],1))==0 and s.simplify(JR-J.get(right['chart'],1))==0
        for key,grid in left['log_radius_mixed_rows'].items():
            differences[name][key]={}
            for label,a in grid.items():
                b=right['log_radius_mixed_rows'][key][label]
                assert a.derivative==b.derivative and a.powers==b.powers and a.source_units==b.source_units
                alog=sum(reader.at(v.node) for _,v in a.log_scale_parts)
                blog=sum(reader.at(v.node) for _,v in b.log_scale_parts)
                assert s.simplify(alog-blog)==0,(name,key,label,'exact common raw scale')
                scales+=1;k,j=a.derivative
                for side,base,jacobian in ((left,a,JL),(right,b,JR)):
                    row=side['native_coordinate_mixed_rows'][key]['n%d_Z%d'%(k,j)]
                    assert row.coefficients is base.coefficients and row.powers[-1]==k
                    total=sum(reader.at(v.node) for _,v in row.log_scale_parts)
                    expected=sum(reader.at(v.node) for _,v in base.log_scale_parts)+k*s.log(jacobian)
                    assert s.simplify(total-expected)==0,(name,key,label,'exact native Jacobian power')
                    native+=1
                al,ah=ends(a.coefficients[0]);bl,bh=ends(b.coefficients[0])
                assert max(al,bl)<=min(ah,bh),(name,key,label,'directed current source inconsistency')
                difference=a.coefficients[0]-b.coefficients[0];lo,hi=ends(difference)
                assert lo<=0<=hi,(name,key,label)
                differences[name][key][label]=dict(signed_common_scale_difference_enclosure=difference,
                    directed_difference_width_bound=owner.ctx.mpf(hi)-owner.ctx.mpf(lo))
                count+=1
    assert (scales,native,count)==(1200,2400,1200)
    return dict(passed=True,exact_lazy_absolute_radius_seams=8,
        exact_common_raw_scale_identities=scales,exact_native_Jacobian_power_identities=native,
        actual_common_scale_overlap_diagnostics=count,
        native_phase_rows_compared_only_after_exact_coordinate_conversion=True),differences


@source_precision
def run(before=None,observed_owner=None,observed_views=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpPostpulseMixedSeams(before,require_checked=False)
    assert not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record and candidate['actual_source_graph']==owner.assert_graph()
    assert candidate['original_power_steep_function_theorems']==owner.canonical
    assert candidate['actual_current_flatten_full_future_function_theorem']==owner.flatten_theorem
    assert candidate['same_actual_C4_C5_complete_future_prefix_binding']==current.mixed.pulse.raw.packed(owner.prefix)
    assert candidate['same_actual_current_angular_steep_forward_history_binding']==current.mixed.pulse.raw.packed(owner.steep_binding)
    assert candidate['current_original_function_transfer']==current.mixed.pulse.raw.packed(owner.function_transfer)
    assert candidate['accepted_original_function_certificate_flags']==owner.canonical_flags
    assert candidate['fourteen_interface_source_chain_registry']==owner.interface_registry()
    assert len(owner.interface_registry())==14
    if observed_owner is not None:
        assert type(observed_owner) is type(owner) and observed_owner.before is before
        assert observed_owner.transport is owner.transport and observed_owner.graph is owner.graph
        assert observed_owner.hashes==owner.hashes and all(observed_owner.assert_graph().values())
        assert set(observed_views)==set(current.SEAMS)
        views=observed_views
    else:views={name:owner.evaluate(name,candidate['fresh_Z']) for name in current.SEAMS}
    for name,view in views.items():
        assert current.mixed.pulse.raw.packed(current.report(view))==candidate['actual_eight_postpulse_mixed_seam_views'][name]
        assert not any(view[k] for k in current.GATES+current.OPEN)
        owner.assemble(name,candidate['fresh_Z'],view['left'],view['right'])
    exact,differences=exact_scales_and_diagnostics(owner,views)
    rejected=[]
    try:owner.evaluate('unknown','.521')
    except ValueError:rejected.append('unknown_seam')
    else:raise AssertionError('Unknown current postpulse seam accepted')
    try:owner.evaluate_interface('unknown','.521')
    except ValueError:rejected.append('unknown_chain_interface')
    else:raise AssertionError('Unknown chain interface accepted')
    bad=copy.copy(owner);bad.prefix=dict(owner.prefix)
    bad.prefix['first_five_complete_future_rows_are_actual_C4_prefix']=False
    try:bad.source_function_transfer()
    except ValueError:rejected.append('changed_current_complete_future_prefix')
    else:raise AssertionError('Changed actual C4/C5 future prefix accepted')
    bad=copy.copy(owner);bad.steep_binding=dict(owner.steep_binding)
    bad.steep_binding['current_absolute_pressure_and_nonzero_angular_histories_retained']=False
    try:bad.source_function_transfer()
    except ValueError:rejected.append('reset_original_forward_pressure_or_angular_memory')
    else:raise AssertionError('Reset original forward history accepted')
    bad=copy.copy(owner);bad.canonical_flags=copy.deepcopy(owner.canonical_flags)
    bad.canonical_flags['power']['flatten_power_and_power_angular_joins_certified']=False
    try:bad.assert_graph()
    except ValueError:rejected.append('missing_original_power_function_certificate')
    else:raise AssertionError('Missing original function certificate accepted')
    wrong=dict(views['flatten_power']['left']);wrong['geometry']=dict(wrong['geometry'])
    wrong['geometry']['exact_native_coordinate']={'numerator':99,'denominator':1}
    try:owner.assemble('flatten_power','.521',wrong,views['flatten_power']['right'])
    except ValueError:rejected.append('wrong_flatten_boundary')
    else:raise AssertionError('Wrong current flatten boundary accepted')
    wrong=dict(views['power_angular']['left']);wrong['log_radius_mixed_rows']=dict(wrong['log_radius_mixed_rows'])
    wrong['log_radius_mixed_rows']['Ur']=dict(wrong['log_radius_mixed_rows']['Ur'])
    del wrong['log_radius_mixed_rows']['Ur']['y4_Z0']
    try:owner.assemble('power_angular','.521',wrong,views['power_angular']['right'])
    except ValueError:rejected.append('missing_actual_order4_radial_velocity_row')
    else:raise AssertionError('Missing actual Ur row accepted')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        current_six_nonheat_primitive_mixed4_function_identities=len(owner.function_transfer['current_six_nonheat_primitive_mixed4_function_identities']),
        accepted_two_closed_heat_primitive_mixed4_function_identities=120,
        accepted_fourteen_interface_source_chain_registry=owner.interface_registry(),
        same_current_complete_future_and_nonzero_forward_histories_identified=True,
        original_raw_cumulative_density_equations=owner.density_theorem,
        exact_geometry_scale_and_native_transfer=exact,signed_common_scale_difference_diagnostics=differences,
        same_current_ten_raw_outputs_and_independent_P0_preserved=True,
        interval_overlap_is_only_a_consistency_diagnostic=True,
        actual_same_source_typed_observations_not_receipt_scalar_rows_used=True,
        rejected_sources_and_domains=rejected,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.mixed.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_POSTPULSE_MIXED_SEAMS all eight current original raw function joins',flush=True)
    return result


if __name__=='__main__':run()
