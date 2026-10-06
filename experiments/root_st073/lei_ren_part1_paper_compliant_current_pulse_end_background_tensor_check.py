"""Admit the current full pulse-end tensor and five completed traces."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor import (
    CurrentPulseEndBackgroundTensor,VIEWS,EDGES,GATES,OPEN,NAME,RECEIPT,HERE,TENSOR_KEYS,
    sha,pack,encode,source_precision,endpoints)
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes


def physical_rows(view):
    def leaves(value):
        if isinstance(value,dict) and 'signed_coefficient' in value:return [value]
        return sum((leaves(v) for v in (value.values() if isinstance(value,dict) else value)),[])
    return sum((leaves(view[key]) for key in TENSOR_KEYS),[])


def check_view(view):
    if any(view[k] for k in OPEN) or not view['actual_full_stress_not_local_difference'] or not view['current_selected_complete_history_graph_retained'] or not view['pressure_remaining_reduction_before_enclosure']:
        raise ValueError('Actual current full end history/decomposition scope lost')
    for key in ('source_formal_Bhat_rows','source_formal_Mz_rows','source_formal_Mtheta_z_rows',
            'source_unperturbed_complete_energy_rows','source_selected_backward_energy_loss_rows','signed_pressure_relative_pure_swirl_reference_rows'):
        if len(view[key])!=5 or any(v.order!=5 for v in view[key]):raise ValueError('Full axial5 ordinary rows required')
        for row in view[key]:
            for value in row.coefficients:finite(value)
    for label,parts in view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors'].items():
        if len(parts)!=(4 if label=='theta' else 6):raise ValueError('Original full meridional stress sector lost')
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Original stress mixed3 incomplete')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=325:raise ValueError('Full sector tensor/divergence/remainder layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Source zero has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact source factors lost')
    return len(rows),nonzero


def check_common(view,count):
    if any(view[k] for k in OPEN) or len(view['common_actual_tensor_rows'])!=count or view['current_common_tensor_contribution_count']!=count:raise ValueError('Actual common tensor trace layout/scope differs')
    nonzero=0
    for row in view['common_actual_tensor_rows'].values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Zero common trace has bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
    return dict(rows=count,nonzero=nonzero)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseEndBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_pulse_end_tensor_views','current_actual_end_flatten_tensor_interface','current_actual_four_pulse_end_support_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current end tensor manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_pulse_end_tensor_views'])!=set(VIEWS):raise ValueError('Producer acceptance/views scope differs')
    count=nonzero=0;view_counts={}
    for name,args in VIEWS.items():
        view=field.end(*args)
        if encode(pack(view))!=raw['current_actual_pulse_end_tensor_views'][name]:raise ValueError('Current end tensor replay differs: '+name)
        n,m=check_view(view);count+=n;nonzero+=m;view_counts[name]=dict(rows=n,nonzero=m)
    center_views=('first_center','second_center','fresh_first_support','fresh_second_support')
    if any(not view_counts[k]['nonzero'] for k in center_views):raise ValueError('Actual selected support tensor collapsed')
    join=field.interface()
    if encode(pack(join))!=raw['current_actual_end_flatten_tensor_interface']:raise ValueError('Current end-flatten common bounds differ')
    join_count=check_common(join,65)
    fresh=field.interface(Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2');fresh_count=check_common(fresh,65)
    edges={};fresh_edges={}
    for edge in EDGES:
        value=field.support_interface(edge)
        if encode(pack(value))!=raw['current_actual_four_pulse_end_support_tensor_interfaces'][edge['exact_edge']]:raise ValueError('Exact current support tensor trace replay differs')
        if not value['both_one_sided_actual_tensor_traces_identified'] or not value['actual_nonzero_shared_boundary_histories_preserved'] or not value['exact_empty_full_support_evaluation_before_enclosure']:raise ValueError('Actual edge history/source lost')
        edges[edge['exact_edge']]=check_common(value,325);check_view(value['current_actual_full_tensor_and_meridional_remainder_endpoint'])
        value=field.support_interface(edge,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2')
        fresh_edges[edge['exact_edge']]=dict(counts=check_common(value,325),bounds=value['common_actual_tensor_rows'])
    rejected=[]
    for label,kwargs in (('outside_Z',dict(Z=2)),('outside_end',dict(offset=1)),('nonfinite_coordinate',dict(offset='-inf')),
            ('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('foreign_edge',dict(edge={})),
            ('wrong_exact_edge_coordinate',dict(edge=EDGES[0],offset=-3))):
        args=dict(Z='.2',offset='-2.937',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.end(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid current end tensor request admitted')
    clone=copy.copy(field);clone.pulse=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_or_old_native_dispatcher_pulse')
    else:raise ValueError('Foreign/old native pulse admitted')
    clone=copy.copy(field);clone.flatten_power=copy.copy(field.flatten_power);clone.flatten_power.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_flatten_power_tensor')
    else:raise ValueError('Unchecked current tensor chain admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_pulse_end_tensor_view_count=len(VIEWS),current_actual_pulse_end_physical_rows_checked=count,
        current_actual_pulse_end_nonzero_physical_contributions_checked=nonzero,current_actual_pulse_end_view_counts=view_counts,
        current_actual_end_flatten_tensor_interface_rows=join_count,current_actual_fresh_end_flatten_tensor_interface_rows=fresh_count,
        current_actual_four_pulse_end_support_tensor_interface_rows=edges,current_actual_fresh_four_support_tensor_interface_rows=fresh_edges,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=9,
        actual_current_completed_tensor_internal_interface_count=4,current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
        source_theorems={k:raw[k] for k in raw if k.endswith('_theorem')},
        current_selected_controls_and_full_five_moment_graph_retained=True,radial_material_transport_and_three_component_remainder_retained=True,
        exact_current_reduced_pressure_no_tiny_cap_division=True,actual_boundary_tensor_not_local_difference=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='Actual current meridional end tensor, end-flatten attachment and4internal support tensor traces.10actualregions/9adjacent/4internal tensor traces. Remaining pulse/core/axis/angular internal/global cone/lift/NS/independent temporal remainder/energy/points/n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual end tensor/end-flatten/four support traces PASS;10regions/9adjacent/4internal;global/time open',flush=True)
    return result


if __name__=='__main__':run()
