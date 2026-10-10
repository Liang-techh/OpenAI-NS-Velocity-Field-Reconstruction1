"""Focused original raw heat seam function, scale and source checks."""
import copy
import gzip
import json
from pathlib import Path
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_closed_heat_mixed_seams as current
from lei_ren_part1_paper_compliant_current_original_Rp_pulse_mixed_seams_check import interpretation
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def exact_scales_and_diagnostics(owner,views):
    reader,mu=interpretation(owner)
    W=s.Symbol('same_true_original_waiting',positive=True)
    reader.bindings[owner.radius.functions['waiting'].node]=W
    scales=native=count=0;differences={}
    ends=current.mixed.pulse.radius.post.selected.inlet.endpoints
    for name,view in views.items():
        left,right=view['left'],view['right'];differences[name]={}
        assert s.simplify(reader.at(left['geometry']['logR'])-reader.at(right['geometry']['logR']))==0,name
        JL=reader.at(left['geometry']['native_to_log_radius_jacobian'])
        JR=reader.at(right['geometry']['native_to_log_radius_jacobian'])
        assert JL==(W if name=='waiting_collar' else 1) and JR==1
        for key,grid in left['log_radius_mixed_rows'].items():
            differences[name][key]={}
            for label,a in grid.items():
                b=right['log_radius_mixed_rows'][key][label]
                assert a.derivative==b.derivative and a.powers==b.powers and a.source_units==b.source_units
                alog=sum(reader.at(v.node) for _,v in a.log_scale_parts)
                blog=sum(reader.at(v.node) for _,v in b.log_scale_parts)
                assert s.simplify(alog-blog)==0,(name,key,label,'original common scale')
                scales+=1;k,j=a.derivative
                for side,base,J in ((left,a,JL),(right,b,JR)):
                    row=side['native_coordinate_mixed_rows'][key]['n%d_Z%d'%(k,j)]
                    assert row.coefficients is base.coefficients and row.powers[-1]==k
                    total=sum(reader.at(v.node) for _,v in row.log_scale_parts)
                    expected=sum(reader.at(v.node) for _,v in base.log_scale_parts)+k*s.log(J)
                    assert s.simplify(total-expected)==0,(name,key,label,'true native Jacobian power')
                    native+=1
                al,ah=ends(a.coefficients[0]);bl,bh=ends(b.coefficients[0])
                assert max(al,bl)<=min(ah,bh),(name,key,label,'directed source inconsistency')
                difference=a.coefficients[0]-b.coefficients[0];lo,hi=ends(difference)
                assert lo<=0<=hi,(name,key,label)
                differences[name][key][label]=dict(signed_common_scale_difference_enclosure=difference,
                    directed_difference_width_bound=owner.ctx.mpf(hi)-owner.ctx.mpf(lo))
                count+=1
    assert (scales,native,count)==(300,600,300)
    return dict(passed=True,exact_lazy_radius_seams=2,exact_original_common_scale_identities=scales,
        exact_native_Jacobian_power_identities=native,common_scale_overlap_diagnostics=count,
        ordinary_logR_rows_compared_before_native_phase_rows=True),differences


@source_precision
def run(before=None,observed_owner=None,observed_views=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpClosedHeatMixedSeams(before,require_checked=False)
    assert not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record and candidate['actual_source_graph']==owner.assert_graph()
    assert candidate['original_heat_function_theorem']==owner.canonical
    assert candidate['same_current_complete_Gamma_future_binding']==owner.shared_future
    assert candidate['actual_closed_output_AST_bindings']==owner.output_bindings
    assert candidate['current_original_function_transfer']==owner.function_transfer
    assert candidate['accepted_full_Gamma_function_certificate_flags']==owner.canonical_flags
    if observed_owner is not None:
        assert type(observed_owner) is type(owner) and observed_owner.before is before
        assert observed_owner.transport is owner.transport and observed_owner.graph is owner.graph
        assert observed_owner.hashes==owner.hashes and all(observed_owner.assert_graph().values())
        assert set(observed_views)==set(current.SEAMS)
        views=observed_views
    else:views={name:owner.evaluate(name,candidate['fresh_Z']) for name in current.SEAMS}
    for name,view in views.items():
        assert current.mixed.pulse.raw.packed(current.report(view))==candidate['actual_two_closed_heat_mixed_seam_views'][name]
        assert not any(view[k] for k in current.GATES+current.OPEN)
        owner.assemble(name,candidate['fresh_Z'],view['left'],view['right'])
    exact,differences=exact_scales_and_diagnostics(owner,views)
    rejected=[]
    try:owner.evaluate('unknown','.521')
    except ValueError:rejected.append('unknown_seam')
    else:raise AssertionError('Unknown seam accepted')
    bad=copy.copy(owner);bad.canonical_flags=dict(owner.canonical_flags)
    bad.canonical_flags['full_infinite_Gamma_source_and_formal_nonzero_S_retained']=False
    try:bad.assert_graph()
    except ValueError:rejected.append('missing_full_Gamma_function_certificate')
    else:raise AssertionError('Missing full Gamma function certificate accepted')
    bad=copy.copy(owner);bad.shared_future=dict(owner.shared_future)
    bad.shared_future['both_energy_epsilon_atoms_retained']=False
    try:bad.source_function_transfer()
    except ValueError:rejected.append('changed_complete_energy_normalization')
    else:raise AssertionError('Changed current full energy normalization accepted')
    wrong=dict(views['waiting_collar']['left']);wrong['geometry']=dict(wrong['geometry'])
    wrong['geometry']['exact_native_coordinate']={'numerator':0,'denominator':1}
    try:owner.assemble('waiting_collar','.521',wrong,views['waiting_collar']['right'])
    except ValueError:rejected.append('wrong_waiting_boundary')
    else:raise AssertionError('Wrong waiting boundary accepted')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        current_original_primitive_mixed4_function_identities=len(owner.function_transfer['canonical_current_primitive_mixed4_identities']),
        canonical_whole_Z_base_axial5_function_identities_consumed=True,
        current_original_Dtheta_Cp_closure_and_complete_Gamma_energy_functions_consumed=True,
        same_current_ten_raw_outputs_and_independent_P0_preserved=True,
        exact_geometry_scale_and_native_transfer=exact,signed_common_scale_difference_diagnostics=differences,
        interval_overlap_only_used_as_consistency_diagnostic=True,
        actual_same_source_typed_observations_not_receipt_scalar_rows_used=True,
        rejected_sources_and_domains=rejected,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.mixed.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_CLOSED_HEAT_MIXED_SEAMS two current original raw function joins',flush=True)
    return result


if __name__=='__main__':run()
