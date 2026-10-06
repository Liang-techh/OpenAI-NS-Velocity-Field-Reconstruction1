"""Accept full33 source routing while keeping global NS/recursion open."""
import gzip
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch import (
    CurrentModifiedVelocityPressureDispatch,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,
    INHERITED_GATES,GATES,OPEN,native_domains,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch_operator import exact_modified_velocity_dispatch_theorem
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces_operator import ORIGINAL_LABELS,ORIGINAL_POWERS
from lei_ren_part1_paper_compliant_cartesian_field import INDICES

def numeric_bounds(value):
    if isinstance(value,dict):
        if 'lower' in value:return mp.mpf(value['lower']),mp.mpf(value['upper'])
        if 'value' in value:return mp.mpf(value['value']),mp.mpf(value['value'])
    return endpoints(value)

def check_layout(view):
    """Check every stored physical contribution; core keeps its native type."""
    expected={'x%d_y%d_z%d'%index for index in INDICES}
    core=view['physical_layout']=='core_native_Cartesian'
    spatial=view['requested_core_spatial_log_bounds'] if core else view['physical_spatial_cartesian_mixed4']
    time=view['requested_core_time_log_bounds'] if core else view['first_fixed_x_physical_time_derivative']
    if set(spatial)!=expected or set(time)!={'ux','uy','uz','p'}:
        raise ValueError('Full spatial4 and fixed-x time1 velocity/pressure layout required')
    count=0
    for components in list(spatial.values())+[time]:
        if set(components)!={'ux','uy','uz','p'}:raise ValueError('Velocity/absolute pressure component omitted')
        for parts in components.values():
            for row in parts.values():
                count+=1
                bound=row if core else row['log_absolute_upper']
                if bound is not None and not all(mp.isfinite(v) for v in numeric_bounds(bound)):
                    raise ValueError('Nonfinite actual physical derivative contribution')
                if not core:
                    if row['exact_zero']!=(bound is None) or (row['exact_zero'] and row['terms']):
                        raise ValueError('Signed nonzero source contribution erased')
                    if not row['positive_source_exponentials_not_materialized'] or numeric_bounds(row['physical_lambda_exponent'])[1]>=0:
                        raise ValueError('Factored units or whole-Z lambda upper bound changed')
    if core:
        if not view['core_axis_nonsingular_native_map']['includes_axis_without_inverse_radius']:
            raise ValueError('Native nonsingular core map required')
    elif view.get('actual_original_geometry_chart')=='switch_second':
        # Original uncapped switch2 publication leaves the30 identically
        # zero radial derivatives of uz as empty native contributions.
        for i,j,b in INDICES:
            key='x%d_y%d_z%d'%(i,j,b)
            if bool(spatial[key]['uz'])!=(i+j==0):
                raise ValueError('Original switch2 axial radial-zero layout changed')
        if count!=186:raise ValueError('Actual switch2 native zero layout differs')
    elif count!=216:raise ValueError('Full216 spatial/time contribution groups required')
    if view.get('full_point_physical_field_evaluation',False) or view.get('resolved_point_coefficient_values_available',False):
        raise ValueError('Function enclosures cannot claim resolved point values')
    return count

def check_query(view):
    if view['resolved_physical_point_values_available'] or any(view[k] for k in OPEN):
        raise ValueError('Native source query exceeded scoped admission')
    if not view['pieces_are_alternative_source_charts_not_summed_fields'] or not view['source_pieces']:
        raise ValueError('Directed source pieces must remain alternative charts')
    return [check_layout(piece['complete_velocity_pressure_source']) for piece in view['source_pieces']]

def check_original_power_packet(view):
    packet=view['actual_original_power_mixed4_packet'];source=view['actual_original_power_source']['current_original_pre_source']
    if packet['signed_Pstar_sector_inventory']!={label:[power] for label,power in ORIGINAL_POWERS.items()}:
        raise ValueError('Original signed velocity/pressure units not preserved')
    for label,grid in packet['grids'].items():
        for key,terms in grid.items():
            actual=source['physical_velocity_pressure_y_Z_mixed4'][ORIGINAL_LABELS[label]][key]
            if len(terms)!=1 or endpoints(terms[0][1])!=endpoints(actual):
                raise ValueError('Original ordinary mixed4 rows rescaled or replaced')
    return 60

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedVelocityPressureDispatch(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[key] for key in INHERITED_GATES+GATES+OPEN):
        raise ValueError('Actual full33 producer source/scope changed')
    proof=exact_modified_velocity_dispatch_theorem(field)
    if encode(pack(proof))!=raw['exact_actual_modified_velocity_pressure_dispatch_theorem']:
        raise ValueError('Actual canonical source/coordinate/function proof changed')
    saved=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    extra={'axis','unbounded_exterior','straddling_q2','exact_power_endpoint','exact_gap_start','exact_gap_end'}
    if set(saved)!=set(native_domains(field))|extra or len(field.registry.routes)!=33:
        raise ValueError('All33 actual native views plus six coverage views required')
    counts={name:check_query(view) for name,view in saved.items()}
    for name in field.registry.routes:
        if saved[name]['region']!=name:raise ValueError('Native source region relabeled')
    for name in ('O3_power','pulse_gap_end','straddling_q2'):
        pieces=saved[name]['source_pieces']
        if len(pieces)!=2:raise ValueError('Actual local/original or gap source cover requires two pieces')
        if numeric_bounds(pieces[0]['source_coordinate_box'])[1]<numeric_bounds(pieces[1]['source_coordinate_box'])[0]:
            raise ValueError('Directed source coordinate coverage has a gap')
    if numeric_bounds(saved['O3_power']['source_pieces'][0]['same_actual_q_offset_enclosure'])[1]>=2:
        raise ValueError('Local phase piece exceeds same accepted compact domain')
    if not saved['axis']['source_pieces'][0]['complete_velocity_pressure_source']['core_axis_nonsingular_native_map']['axis_only']:
        raise ValueError('Exact axis route did not use native nonsingular axis source')
    exterior=saved['unbounded_exterior']['source_pieces'][0]['complete_velocity_pressure_source']
    if not exterior['unbounded_coordinate_is_whole_domain_box_not_point_at_infinity'] or not exterior['same_checked_full_Gamma_exterior_source_and_absolute_pressure_retained']:
        raise ValueError('Whole unbounded heat query lost its actual Gamma source')
    switch=field.geometry.dispatch.evaluate('switch_second',(-1,1),(1,2))['source_packet']
    radial_zero_rows=[row for row in switch['final_factored_physical_row_ledgers']
      if row['physical_row'].startswith('Uz/s') and int(row['physical_row'].split('/s')[1].split('_')[0])>0]
    if len(radial_zero_rows)!=10 or any(row['terms'] for row in radial_zero_rows):
        raise ValueError('Omitted switch2 physical terms are not actual native zero radial rows')
    fresh={};fresh_views={}
    queries=(('entrance',lambda:field.native('pulse_entrance','.537','.013337','-2.337','.337')),
      ('gap_lower',lambda:field.gap_phase('.537','0.000013337','-2.337','.337')),
      ('gap_upper',lambda:field.gap_phase('.537','.337','-2.337','.337')),
      ('gap_true_end',lambda:field.gap_phase('.537',1,'-2.337','.337')),
      ('original_power',lambda:field.power_offset('.537','2.413337','-2.337','.337')),
      ('power_true_end',lambda:field.power_phase('.537',1,'-2.337','.337')),
      ('axis',lambda:field.axis('.537','-2.337','.337')))
    original_rows=0
    for name,query in queries:
        view=query();fresh_views[name]=view;fresh[name]=check_query(view)
        for piece in view['source_pieces']:
            inner=piece['complete_velocity_pressure_source']
            if 'actual_original_power_mixed4_packet' in inner:original_rows+=check_original_power_packet(inner)
        print('Checked fresh native velocity/pressure route: '+name,flush=True)
    # Compare entrance bridge with the actual canonical callable, and
    # original continuation with its actual pre-pulse source (not a cone).
    entrance=fresh_views['entrance']['source_pieces'][0]['complete_velocity_pressure_source']
    direct=field.geometry.evaluate('pulse_entrance','.537',field.ctx.mpf('.013337')/field.geometry.pulse.mu,'-2.337','.337')
    for key in ('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','source_logR_enclosure'):
        if encode(pack(entrance[key]))!=encode(pack(direct[key])):raise ValueError('Entrance canonical source bridge differs')
    offset=fresh_views['original_power']['source_pieces'][0]['complete_velocity_pressure_source']
    source=field.pre.power(field.ctx.mpf('.537'),field.ctx.mpf('2.413337')/field.Tw)
    if encode(pack(offset['actual_original_power_source']['current_original_pre_source']))!=encode(pack(source)):
        raise ValueError('Continuation lost actual original pre/P0/histories')
    for action in (lambda:field.native('foreign',0,0),lambda:field.gap_phase(0,'1.001'),
      lambda:field.power_phase(0,'-0.001'),lambda:field.native('heat_exterior',0,mp.inf)):
        try:action()
        except ValueError:pass
        else:raise ValueError('Invalid native or point-at-infinity query admitted')
    result=dict(actual_five_defect_family_sha256=field.registry.family,implicit_source_sha256=field.registry.source,
      datum_enclosure_sha256=field.registry.datum_sha,finite_integer_N=field.N,
      modified_velocity_pressure_dispatch_definition_sha256=field.modified_velocity_pressure_dispatch_definition_sha256,
      new_exact_source_and_coordinate_routing_identities=len(proof['identities']),
      complete_native_velocity_pressure_chart_count=33,actually_changed_native_route_count=3,
      complete_whole_Z_time_velocity_pressure_view_counts=counts,fresh_actual_query_contribution_counts=fresh,
      independent_original_power_packet_ordinary_rows_checked=original_rows,
      consumed_twelve_affected_velocity_pressure_interface_identities=len(field.interfaces.theorem['identities']),
      true_correlated_power_endpoint_and_both_reciprocal_gap_endpoints_available=True,
      core_axis_and_whole_unbounded_heat_box_keep_distinct_native_layouts=True,
      original_switch_second_thirty_zero_axial_radial_contributions_preserved=True,
      physical_spatial_multiindices=35,
      global_velocity_interfaces_locator_energy_common_N_cones_recursion_resolved_NS_remain_open=True,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(INHERITED_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current full33 native velocity/absolute-pressure dispatcher PASS; actual recursion/resolved NS open',flush=True)
    return result

if __name__=='__main__':run()
