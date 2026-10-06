"""Check twelve actual two-germ mixed4 functions and their physical traces."""
import gzip
import json
import math
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces import (
    CurrentModifiedVelocityPressureInterfaces,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,
    INHERITED_GATES,GATES,OPEN,SEAMS,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces_operator import (
    physical_contribution_groups,exact_modified_velocity_pressure_interfaces_theorem,ORIGINAL_LABELS,ORIGINAL_POWERS)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_operator import LABELS,UT,UZ,UR,P
from lei_ren_part1_paper_compliant_cartesian_field import INDICES

def check_side(view):
    expected={'x%d_y%d_z%d'%index for index in INDICES}
    if set(view['physical_spatial_cartesian_mixed4'])!=expected:
        raise ValueError('Full spatial4 Cartesian multiindex inventory required')
    if view['resolved_point_coefficient_values_available']:
        raise ValueError('Boundary source enclosures cannot claim resolved point coefficients')
    packet=view['actual_boundary_mixed4_packet'];source=view['actual_boundary_source']
    original=view['actual_boundary_source_kind']=='original_native_germ'
    inventory=({label:[power] for label,power in ORIGINAL_POWERS.items()} if original
      else {UT:[1],UZ:[1],UR:[0,1],P:[2]})
    if packet['signed_Pstar_sector_inventory']!=inventory:
        raise ValueError('Actual original/own signed source units changed')
    keys={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}
    if set(packet['grids'])!=set(inventory) or any(set(grid)!=keys for grid in packet['grids'].values()):
        raise ValueError('Four-label full ordinary mixed4 source inventory required')
    # Check every public adapter row against the actual source consumed
    # by this side. Original rows are ordinary already; own Z jets get n!.
    count=0
    for label,grid in packet['grids'].items():
        for j in range(5):
            for n in range(5-j):
                terms=grid['y%d_Z%d'%(j,n)]
                if original:
                    values=[source['original_native_velocity_pressure_source']['physical_velocity_pressure_y_Z_mixed4'][ORIGINAL_LABELS[label]]['y%d_Z%d'%(j,n)]]
                elif label==P:
                    values=[source['modified_absolute_pressure_over_Pstar2_ordinary_logR_rows'][j][n]*math.factorial(n)]
                else:
                    component=next(name for name,key in LABELS.items() if key==label)
                    values=[part['ordinary_logR_rows'][j][n]*math.factorial(n)
                      for part in source['modified_cylindrical_velocity_source_log_sectors'][component]]
                if len(terms)!=len(values):raise ValueError('Signed source sector omitted')
                for term,value in zip(terms,values):
                    if endpoints(term[1])!=endpoints(value):raise ValueError('Actual ordinary boundary source adapter differs')
                    count+=1
    if count!=(60 if original else 75):raise ValueError('Actual signed mixed4 input count differs')
    if view['named_source_edge']!='buffer_transition':
        if any(endpoints(term[1])!=(0,0) for terms in packet['grids'][UZ].values() for term in terms):
            raise ValueError('Actual flat boundary axial mixed4 rows must be exactly zero')
    for groups in list(view['physical_spatial_cartesian_mixed4'].values())+[view['first_fixed_x_physical_time_derivative']]:
        if set(groups)!={'ux','uy','uz','p'}:raise ValueError('Four physical component traces required')
        for component,parts in groups.items():
            if set(parts)!=({'ux':{UR,UT},'uy':{UR,UT},'uz':{UZ},'p':{P}}[component]):
                raise ValueError('Actual moving-basis physical source labels changed')
            for row in parts.values():
                if endpoints(row['physical_lambda_exponent'])[1]>=0:raise ValueError('Whole-Z lambda bound requires negative powers')
                if not row['positive_source_exponentials_not_materialized'] or not row['original_source_factors_combined_before_enclosure']:
                    raise ValueError('Positive large source units capped or materialized')
                if row['exact_zero']:
                    if row['terms'] or row['log_absolute_upper'] is not None:raise ValueError('Nonzero source erased')
                elif not all(mp.isfinite(v) for v in endpoints(row['log_absolute_upper'])):
                    raise ValueError('Nonfinite physical contribution bound')
                for term in row['terms']:
                    if not all(mp.isfinite(v) for v in endpoints(term['signed_coefficient'])+endpoints(term['log_absolute_upper'])):
                        raise ValueError('Nonfinite signed physical source term')
    return count

def check_interface(view):
    if any(view[key] for key in OPEN) or view['resolved_physical_point_values_available']:
        raise ValueError('Twelve interfaces cannot admit full33 dispatch, energy or recursion')
    if not view['source_function_and_FTC_identities_precede_physical_bounds'] or not view['interval_overlap_or_midpoint_not_used_as_function_identity']:
        raise ValueError('Function equality must precede common interval bounds')
    left=view['left_complete_velocity_pressure_source'];right=view['right_complete_velocity_pressure_source']
    count=check_side(left)+check_side(right)
    a=physical_contribution_groups(left);b=physical_contribution_groups(right)
    common=view['common_physical_velocity_pressure_rows']
    if set(common)!=set(a) or set(a)!=set(b) or len(common)!=216:
        raise ValueError('All 216 spatial4/time1 physical contribution groups required')
    for key,row in common.items():
        uppers=[endpoints(part['log_absolute_upper'])[1] for part in (a[key],b[key]) if not part['exact_zero']]
        if row['exact_zero']!=(not uppers) or not row['both_actual_side_bounds_retained']:
            raise ValueError('Common trace incorrectly cleared a signed source')
        if uppers:
            if endpoints(row['log_absolute_upper'])[0]<max(uppers):raise ValueError('Common trace does not enclose both actual sides')
        elif row['log_absolute_upper'] is not None:raise ValueError('Exact zero trace must have no finite logarithm')
    return count

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedVelocityPressureInterfaces(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[key] for key in INHERITED_GATES+GATES+OPEN):
        raise ValueError('Actual twelve-interface producer scope/source changed')
    proof=exact_modified_velocity_pressure_interfaces_theorem(field)
    if encode(pack(proof))!=raw['exact_actual_modified_velocity_pressure_interface_theorem']:
        raise ValueError('Actual source-program mixed4 comparison differs')
    coupled=proof['actual_coupled_velocity_pressure_packet_source_joins']
    if len(coupled['identities'])!=720 or not all(coupled['identities'].values()) or set(coupled['actual_source_germ_input_cases'])!=set(SEAMS):
        raise ValueError('All twelve actual two-packet source programs must agree in 60 mixed4 rows each')
    saved=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if set(saved)!=set(SEAMS):raise ValueError('Twelve whole-Z/time interface views required')
    counts={}
    for name in SEAMS:
        view=field.interface(name)
        if encode(pack(view))!=saved[name]:raise ValueError('Actual complete two-germ view changed: '+name)
        counts[name]=check_interface(view)
        print('Checked actual complete physical velocity/pressure interface: '+name,flush=True)
    fresh={}
    for name in ('buffer_inlet','transition_modulation_end','repair0_out','repair_exit'):
        view=field.interface(name,'.537','-2.337','.337')
        fresh[name]=check_interface(view)
    for name in ('transition_modulation_end','repair0_out'):
        case=coupled['actual_source_germ_input_cases'][name]
        if not case['nonzero_cumulative_moment_and_pressure_kept']:
            raise ValueError('Flat local perturbation must retain actual cumulative histories')
    try:field.interface('heat_exterior')
    except ValueError:pass
    else:raise ValueError('Affected-interface API accepted a foreign native region')
    result=dict(actual_five_defect_family_sha256=field.velocity.heat.registry.family,
      implicit_source_sha256=field.velocity.heat.registry.source,datum_enclosure_sha256=field.velocity.heat.registry.datum_sha,
      finite_integer_N=field.N,modified_velocity_pressure_interfaces_definition_sha256=field.modified_velocity_pressure_interfaces_definition_sha256,
      new_source_and_two_germ_operator_identities=len(proof['identities']),
      actual_coupled_packet_mixed4_source_identities=len(coupled['identities']),
      consumed_modified_tensor_source_interface_identities=proof['consumed_checked_modified_source_interface_identities'],
      consumed_full_velocity_pressure_source_identities=proof['consumed_checked_velocity_pressure_source_identities'],
      complete_whole_Z_time_physical_interface_signed_source_rows=counts,
      fresh_actual_boundary_signed_source_rows=fresh,
      common_physical_contribution_groups_per_interface=216,physical_spatial_multiindices=35,
      whole_Z_common_functions_proved_before_bounds=True,
      exact_original_geometry_absolute_P0_and_nonzero_cumulative_memory_retained=True,
      twelve_affected_velocity_pressure_interfaces_do_not_admit_global_energy_cones_or_recursion=True,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(INHERITED_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Twelve actual velocity/absolute-pressure spatial4/time1 physical interfaces PASS; full33 velocity dispatch/energy/cones/recursion open',flush=True)
    return result

if __name__=='__main__':run()
