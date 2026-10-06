"""Admit current full main/exit tensors and their two source-function joins."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_main_exit_background_tensor import (
    CurrentPulseMainExitBackgroundTensor,VIEWS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    if any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','current_selected_complete_history_graph_retained',
            'pressure_remaining_reduction_before_enclosure','all_local_incoming_cross_and_square_radial_terms_retained',
            'positive_omitted_forward_tails_retained','source_bounds_not_resolved_physical_point_values')):
        raise ValueError('Current full main/exit source scope lost')
    full=view['source_full_main_exit_history_rows']
    for key in ('Bh','ml','mi','nl','ni','e0','J','pressure_baseline_rows','pressure_memory_rows'):
        if len(full[key])!=5 or any(row.order!=5 for row in full[key]):raise ValueError('Full axial5 ordinary history required')
        for row in full[key]:
            for value in row.coefficients:finite(value)
    for row in [view['current_selected_ap'],view['current_complete_C5_future'],
            view['current_signed_absolute_Rv_pressure'],view['current_incoming_energy']]+view['current_incoming_moments']+view['current_selected_controls']:
        if row.order!=5:raise ValueError('Current full C5 selection/history/pressure source lost')
        for value in row.coefficients:finite(value)
    for kernel in view['original_partial_linear_kernel_bounds']:
        if not kernel['omitted_history_not_deleted'] or endpoints(kernel['positive_omitted_tail_bound'])[1]<=0:
            raise ValueError('Original positive omitted forward history lost')
        finite(kernel['enclosure']);finite(kernel['positive_omitted_tail_bound'])
    source=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(source['theta'])!=6 or len(source['axial'])!=9:raise ValueError('Full main/exit incoming stress splitter incomplete')
    for parts in source.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=508 or len(canonical_tensor_groups(view))!=71:raise ValueError('Full actual main/exit tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Source zero carries nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Current exact physical source factors lost')
    radial=view['physical_three_component_remainder_mixed2']['radial']
    if not any(not row['exact_zero'] for sectors in radial.values() for row in sectors.values()):raise ValueError('Full radial remainder collapsed')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentPulseMainExitBackgroundTensor(require_checked=False)
    if not all(field.proof['current_fixed_gp_provider_bindings'].values()):raise ValueError('Current GP provider admission lost')
    for key in ('original_generic_full_main_exit_stress_theorem','original_generic_full_main_exit_physical_and_exit_gap_theorem'):
        if not field.proof[key]['exact_reviewed_legacy_statement_allowlist_enforced']:raise ValueError('Exact legacy receipt allowlist lost')
    field.assert_graph();omitted={'current_actual_main_exit_tensor_views','current_actual_two_main_exit_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current main/exit defining source/manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_main_exit_tensor_views'])!=set(VIEWS):raise ValueError('Producer main/exit acceptance/views exceed scope')
    counts={}
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_main_exit_tensor_views'][name]:raise ValueError('Actual main/exit source replay differs: '+name)
        counts[name]=check_view(value)
        if name=='exit_gap':
            for row in value['source_full_main_exit_history_rows']['Bh']:
                if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in row.coefficients):raise ValueError('Original xi11 flat shear lost')
            if endpoints(value['original_partial_future_energy_bound'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Exact xi11 empty remaining energy lost')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_two_main_exit_tensor_interfaces'][name]:raise ValueError('Actual main/exit tensor function trace replay differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('outside_main',('pulse_main','.2',11)),('outside_exit',('pulse_exit','.2',9)),
            ('outside_Z',('pulse_main',2,4)),('foreign_chart',('pulse_end','.2',4)),('nonfinite_coordinate',('pulse_main','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual main/exit domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('pulse_main','.2',4,**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual main/exit physical sector admitted')
    clone=copy.copy(field);clone.pulse=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_or_old_native_pulse')
    else:raise ValueError('Foreign actual main/exit pulse owner admitted')
    clone=copy.copy(field);clone.gap_tensor=copy.copy(field.gap_tensor);clone.gap_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_gap_tensor')
    else:raise ValueError('Unchecked current gap tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_main_exit_tensor_view_count=len(VIEWS),current_actual_main_exit_view_counts=counts,
        current_actual_main_exit_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_main_exit_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_two_main_exit_tensor_interface_rows=joins,current_actual_two_fresh_main_exit_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=13,
        actual_current_completed_tensor_internal_interface_count=4,current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
        current_source_function_tensor_trace_theorem=field.proof,current_full_C5_future_selected_incoming_histories_retained=True,
        all_local_incoming_cross_and_square_radial_terms_retained=True,original_positive_omitted_kernel_tails_retained=True,
        gzip_contains_complete_unpruned_producer=True,invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='Actual current main/exit completed physical tensors and two tensor function joins.14 regions/13 adjacent/4 internal tensor traces; entrance/core/axis/angular internal/global cone/lift/temporal remainder/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual main/exit tensor and two joins PASS;14 regions/13 adjacent/4 internal;global/time open',flush=True)
    return result


if __name__=='__main__':run()
