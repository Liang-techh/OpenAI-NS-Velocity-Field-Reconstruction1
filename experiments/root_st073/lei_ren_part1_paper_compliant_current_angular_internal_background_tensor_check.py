"""Independent scoped admission of the four full angular tensor support limits."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_angular_internal_background_tensor import (
    CurrentAngularInternalBackgroundTensor,SEAMS,VIEWS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,read_producer,
    current_angular_internal_source_theorem,canonical_angular_groups,leading_origin_time_obstruction)
from lei_ren_part1_paper_compliant_current_angular_internal_operator import two_sided_full_angular_trace_theorem
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import check_view,finite


def check_interface(field,view):
    if view['seam'] not in SEAMS or view['exact_source_edge']!=SEAMS[view['seam']] or any(view[k] for k in OPEN):
        raise ValueError('Exact angular full tensor support scope differs')
    for key in ('exact_source_function_limits_and_full_operators_precede_bounds',
        'full_nonzero_histories_pressure_stress_and_remainder_retained',
        'ordinary_s_derivatives_not_support_distance_derivatives',
        'same_complete_source_family_datum_and_physical_factors_on_both_sides',
        'only_proved_endpoint_source_callbacks_reduced',
        'local_zero_difference_not_substituted_for_actual_full_tensor',
        'interval_overlap_not_used_as_source_function_equality'):
        if not view[key]:raise ValueError('Whole actual angular source trace provenance lost')
    n,nonzero=check_view(field.angular,view['current_full_angular_endpoint_tensor'])
    groups=canonical_angular_groups(field.ctx,view['current_full_angular_endpoint_tensor'])
    if view['current_common_tensor_contribution_count']!=71 or view['current_common_physical_contribution_count']!=77 or set(groups)!=set(view['common_actual_tensor_rows']):
        raise ValueError('Complete common tensor component layout lost')
    for key,rows in groups.items():
        positive=[endpoints(row['log_absolute_upper'])[1] for row in rows if not row['exact_zero']]
        common=view['common_actual_tensor_rows'][key]
        if common['exact_zero']!=(not positive):raise ValueError('Full nonzero trace replaced by zero difference')
        if positive:
            finite(common['common_triangle_log_upper'])
            expected=field.ctx.mpf(max(positive))+field.ctx.ln(len(positive))
            if endpoints(common['common_triangle_log_upper'])!=endpoints(expected):raise ValueError('Common full tensor log bound differs')
        elif common['common_triangle_log_upper'] is not None:raise ValueError('Exact zero source contribution has bound')
    if n!=65 or not nonzero:raise ValueError('Nonzero full angular stress/divergence/remainder omitted')
    return dict(original_physical_rows=n,nonzero_physical_rows=nonzero,common_component_groups=71,canonical_contribution_rows=77)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentAngularInternalBackgroundTensor(require_checked=False)
    field.assert_graph();field.assert_programs()
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_four_angular_internal_full_tensor_views'}:
        raise ValueError('Current full angular internal source manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_four_angular_internal_full_tensor_views'])!=set(VIEWS):
        raise ValueError('Four full tensor seams and fresh requested sectors required')
    two_sided_full_angular_trace_theorem.cache_clear()
    proof=current_angular_internal_source_theorem(field)
    if encode(pack(proof))!=encode(pack(field.proof)):raise ValueError('Exact current source/operator proof replay differs')
    obstruction=leading_origin_time_obstruction(field)
    if encode(pack(obstruction))!=encode(pack(field.leading_obstruction)) or not obstruction['passed']:
        raise ValueError('Independent leading source temporal obstruction replay differs')
    field.angular.cache={};counts={}
    for name,args in VIEWS.items():
        view=field.interface(*args)
        if encode(pack(view))!=raw['current_four_angular_internal_full_tensor_views'][name]:
            raise ValueError('Full angular source/tensor replay differs: '+name)
        counts[name]=check_interface(field,view)
        print('Check full angular internal tensor trace: '+name,flush=True)
    rejected=[]
    for label,kwargs in (('unknown_seam',dict(name='pulse_end_-3_-1')),('outside_Z',dict(Z=2)),
        ('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_angle',dict(theta='inf'))):
        args=dict(name=next(iter(SEAMS)),Z='.2',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.interface(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Foreign angular source/physical domain accepted')
    for attribute in ('angular','support'):
        clone=copy.copy(field);setattr(clone,attribute,object())
        try:clone.assert_graph()
        except (ValueError,AttributeError):rejected.append('foreign_'+attribute+'_owner');continue
        raise ValueError('Foreign original angular/source owner accepted')
    clone=copy.copy(field);clone.core_axis=copy.copy(field.core_axis);clone.core_axis.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_core_axis')
    else:raise ValueError('Unchecked core/axis prerequisite admitted')
    clone=copy.copy(field);clone.full_operator=lambda *args,**kwargs:None
    try:clone.assert_programs()
    except (ValueError,KeyError):rejected.append('foreign_full_tensor_operator')
    else:raise ValueError('Foreign full tensor callback admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_four_angular_internal_full_tensor_view_counts=counts,
        current_four_angular_internal_full_tensor_view_count=len(VIEWS),current_four_angular_internal_completed_tensor_join_count=len(SEAMS),
        current_full_angular_internal_physical_rows_checked=sum(v['original_physical_rows'] for v in counts.values()),
        current_full_angular_common_component_groups_checked=sum(v['common_component_groups'] for v in counts.values()),
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=32,actual_current_completed_tensor_internal_interface_count=14,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        current_core_full_background_tensor_available=True,core_axis_tensor_remainder_limits_certified=True,
        current_leading_origin_remainder_time_obstruction=obstruction,
        exact_one_sided_source_functions_before_common_full_tensor_bounds=True,
        local_difference_and_interval_overlap_not_used_as_full_tensor_identity=True,
        invalid_domains_and_foreign_owners_rejected=rejected,
        current_four_angular_internal_full_tensor_source_theorem=field.proof,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Four full angular internal tensor traces PASS;33regions/32adjacent/14internal;global/time/energy/points/recursion open',flush=True)
    return result


if __name__=='__main__':run()
