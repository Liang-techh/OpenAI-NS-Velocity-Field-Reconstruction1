"""Scoped one-graph routing admission; global physical location stays open."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_tensor_registry import (
    CurrentTensorRegistry,NAME,RECEIPT,GATES,OPEN,HERE,sha,pack,encode,endpoints,
    source_precision,_verify_hashes,defining_method)

COORDINATES=dict(core_positive_radius='.137',bridge_first='.417',bridge_second='1.417',bridge_macro='.417',
    switch_first='.417',switch_second='1.417',switch_power='.417',reshape='.417',
    inner_reference='.417',axial_restore='.417',restore_buffer='-6.417',actual_patch='1.37',
    Rh_reference='-2.6',O2_slope='.417',O2_axial='.417',O2_buffer='5.37',O3_slope_mu='.417',
    O3_power='.417',pulse_entrance='.007',pulse_main='.537',pulse_exit='10.417',
    pulse_gap='11.417',pulse_gap_end='.417',pulse_end='-2.6',flatten='41.7',outer_power='.417',
    outer_angular='-2.6',steep_entry='.417',steep_power='.417',steep_exit='.417',waiting='.417',
    heat_collar='1.7',heat_exterior='4.3')

def bounds(value):
    return endpoints(value) if hasattr(value,'_mpi_') else (mp.mpf(value),mp.mpf(value))

def check_full_view(view,expected_groups=71):
    groups=view['canonical_signed_component_groups']
    if len(groups)!=expected_groups or not view['original_full_source_sectors_and_factors_retained'] or any(view[k] for k in OPEN):
        raise ValueError('Complete original tensor layout/scope required')
    rows=nonzero=0
    for label,parts in groups.items():
        if not parts:raise ValueError('Signed source sector omitted: '+label)
        for row in parts:
            for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent',
                        'physical_viscosity_log_prefactor','radial_log_prefactor'):
                if not all(mp.isfinite(v) for v in bounds(row[key])):raise ValueError('Nonfinite exact signed factor')
            if not row['source_factors_combined_before_enclosure'] or not row['positive_source_factors_not_materialized']:
                raise ValueError('Original factored physical source convention lost')
            if any(not all(mp.isfinite(v) for v in bounds(value)) for value in row['actual_source_log_parts'].values()):
                raise ValueError('Nonfinite original amplitude/radius factor')
            if row['exact_zero']:
                if bounds(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)) or row['log_absolute_upper'] is not None:
                    raise ValueError('Declared zero full source row is nonzero')
            else:
                if not all(mp.isfinite(v) for v in bounds(row['log_absolute_upper'])):raise ValueError('Missing full contribution bound')
                nonzero+=1
            rows+=1
    return dict(groups=len(groups),signed_contributions=rows,nonzero_signed_contributions=nonzero,layout=view['layout'])

@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentTensorRegistry(require_checked=False)
    field.assert_graph()
    if encode(pack(field.manifest()))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Native source registry manifest differs')
    if set(COORDINATES)!=set(field.registry):raise ValueError('All 33 actual native regions required')
    edges=[(row['left_region'],row['right_region']) for row in field.adjacent.values()]
    vertices=set(sum(([a,b] for a,b in edges),[]))
    if len(edges)!=32 or vertices!=set(field.registry) or len(set(edges))!=32:
        raise ValueError('Regional adjacency ledger is incomplete or duplicate')
    if any(edges[j][1]!=edges[j+1][0] for j in range(31)):raise ValueError('Named radial adjacency chain has a gap')
    kinds={name:sum(row['region']==name for row in field.internal.values()) for name in ('actual_patch','pulse_end','outer_angular')}
    if kinds!=dict(actual_patch=6,pulse_end=4,outer_angular=4):raise ValueError('Support trace inventory or axis counting differs')
    for key,binding in field.bindings.items():
        alias,method=key.split('.')
        if defining_method(field.owners[alias],method)!=binding:raise ValueError('Original method AST/signature replay differs')
    counts={}
    for region,coordinate in COORDINATES.items():
        counts[region]=check_full_view(field.native(region,'.537',coordinate,'-2.6','.41','.8'))
        print('Check current tensor native route: '+region,flush=True)
    axis=field.native('core_positive_radius','.537',0,'-2.6','.41','.8')
    if not axis['original_complete_view']['exact_axis']:raise ValueError('rho0 not routed through exact nonsingular axis')
    counts['core_axis_extension']=check_full_view(axis,210)
    counts['core_Cartesian_disk']=check_full_view(field.core_cartesian('.317','.537','-.419','-2.7','.2'),210)
    tail=field.exterior(Z='.731',log_tau='-2.9',theta='.37',viscosity='.2')
    counts['unbounded_Gamma']=check_full_view(tail)
    if counts['unbounded_Gamma']['nonzero_signed_contributions']:
        raise ValueError('Exact full Gamma zero tensor/NS source lost')
    # The admitted endpoint recipes, rather than rounded native coordinates,
    # own support limits and reciprocal/named radial joins.
    traces={}
    for name in ('gap_coordinate','end_flatten','angular_entry','patch_support_49',
                 'pulse_end_-63/20','outer_angular_-3_-1'):
        value=field.trace(name,Z='.419',log_tau='-2.7',theta='.37',viscosity='.2')
        if not value['common_actual_tensor_rows'] or any(value[k] for k in OPEN if k in value):raise ValueError('Admitted full source trace scope lost')
        traces[name]=len(value['common_actual_tensor_rows'])
    rejected=[]
    for label,args in (('unknown_region',('core_axis','.2',0)),('outside_core',('core_positive_radius','.2',5)),
                       ('negative_core_radius',('core_positive_radius','.2','-.1'))):
        try:field.native(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid native source domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_nu',dict(viscosity=0)),('nonfinite_axis_angle',dict(theta='inf'))):
        try:field.native('core_positive_radius','.2',0,**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid axis physical input admitted')
    clone=copy.copy(field);clone.owners=dict(field.owners);clone.owners['heat']=object()
    try:clone.assert_graph()
    except (ValueError,AttributeError):rejected.append('foreign_heat_owner')
    else:raise ValueError('Foreign tensor graph admitted')
    clone=copy.copy(field);clone.routes=dict(field.routes);clone.routes['core_positive_radius']=clone.routes['heat_exterior']
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_native_route')
    else:raise ValueError('Native owner route change admitted')
    try:field.trace('core_axis')
    except ValueError:rejected.append('axis_not_counted_as_support_trace')
    else:raise ValueError('Axis relabelled as a counted support boundary')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_native_tensor_route_count=33,current_native_tensor_route_full_view_counts=counts,
        current_original_method_AST_bindings_replayed=len(field.bindings),current_actual_fresh_trace_route_counts=traces,
        actual_current_tensor_regions_available=list(field.registry),actual_current_completed_tensor_adjacent_interface_count=32,
        actual_current_completed_tensor_internal_interface_count=14,internal_count_composition=kinds,
        core_axis_is_separate_analytic_extension_not_a_counted_internal_support=True,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        exact_original_ordinary_derivative_and_signed_source_factors_retained=True,
        return_is_native_source_enclosure_not_resolved_physical_point=True,invalid_domains_and_foreign_routes_rejected=rejected,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current 33 native tensor routes/32 adjacent/14 support registry PASS; physical locator/global cover/cone/recursion open',flush=True)
    return result

if __name__=='__main__':run()
