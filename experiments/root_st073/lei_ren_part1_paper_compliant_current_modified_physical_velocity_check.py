"""Check physical candidate composition, actual lambda and constant nu."""
import gzip
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_current_modified_physical_velocity import (
    CurrentModifiedPhysicalVelocity,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,
    INHERITED_GATES,GATES,OPEN,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_physical_velocity_operator import exact_modified_physical_velocity_theorem
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch_check import check_query
from lei_ren_part1_paper_compliant_cartesian_field import INDICES

def number(c,value):
    if isinstance(value,dict):
        if 'lower_exact_mpf_tuple' in value:return read_interval(c,value)
        if 'exact_mpf_tuple' in value:return c.mpf(mp.make_mpf(tuple(value['exact_mpf_tuple'])))
    return c.mpf(value)

def check_physical_query(field,query):
    if query['resolved_physical_point_values_available'] or any(query[key] for key in OPEN):
        raise ValueError('Physical source enclosures exceeded actual scope')
    location=query['location'];c=field.ctx
    candidates=location['candidates'];views=query['current_candidate_velocity_pressure_views']
    if not candidates or len(candidates)!=len(views) or not query['candidate_union_and_native_pieces_are_alternative_charts_not_summed_fields']:
        raise ValueError('Every directed physical candidate must remain in the union')
    q=number(c,location['actual_log_lambda']);lt=number(c,location['requested_log_tau'])
    nu=number(c,location['physical_viscosity']);lnnu=c.ln(nu)
    if endpoints(nu)[0]<=0 or not all(mp.isfinite(v) for v in endpoints(q)):
        raise ValueError('Positive fixed nu and finite actual lambda bracket required')
    count=0
    for candidate,view in zip(candidates,views):
        if encode(pack(candidate))!=encode(pack(view['actual_physical_candidate'])):
            raise ValueError('Physical candidate dropped, modified or reordered')
        native=view['current_native_velocity_pressure_query'];check_query(native)
        if len(native['source_pieces'])!=len(view['actual_physical_velocity_pressure_pieces']):
            raise ValueError('Native source partitions dropped in physical mapping')
        for piece,physical in zip(native['source_pieces'],view['actual_physical_velocity_pressure_pieces']):
            if encode(pack(piece))!=encode(pack(physical['native_source_piece'])):
                raise ValueError('Actual signed source/function payload changed')
            actual=physical['actual_physical_velocity_pressure_bounds'];source=piece['complete_velocity_pressure_source']
            core=source['physical_layout']=='core_native_Cartesian'
            spatial=actual['actual_lambda_viscosity_spatial4_bounds'];time=actual['actual_lambda_viscosity_fixed_x_time1_bounds']
            if set(spatial)!={'x%d_y%d_z%d'%index for index in INDICES}:
                raise ValueError('All35 actual physical spatial derivatives required')
            for i,j,k in INDICES:
                index='x%d_y%d_z%d'%(i,j,k);order=i+j+k
                original=(source['core_axis_nonsingular_native_map']['cartesian_spatial_multiindices'][index] if core
                  else source['physical_spatial_cartesian_mixed4'][index])
                for component,parts in original.items():
                    expected_power=c.mpf((2 if component=='p' else 1)-order)/2
                    if set(spatial[index][component])!=set(parts):raise ValueError('Actual source contribution layout changed')
                    for label,row in parts.items():
                        result=spatial[index][component][label];count+=1
                        power=number(c,result['constant_viscosity_power'])
                        if endpoints(power)!=endpoints(expected_power):raise ValueError('Constant-viscosity spatial chain rule differs')
                        if core:
                            scale=source['core_axis_nonsingular_native_map']['shared_physical_prefactor_bounds'][row['scale_key']]
                            norm=endpoints(number(c,row['absolute_upper']))[1];gamma=number(c,scale['physical_lambda_exponent'])
                            expected=None if not norm else c.ln(c.mpf(norm))+number(c,scale['logLambda_term'])+number(c,scale['amplitude_log_upper'])+gamma*q+power*lnnu
                        else:
                            gamma=number(c,row['physical_lambda_exponent'])
                            expected=None if row['exact_zero'] else number(c,row['log_absolute_upper'])+gamma*(q-lt/2)+power*lnnu
                        check_bound(c,result,expected)
            original=(source['core_axis_nonsingular_native_map']['first_fixed_x_physical_time_derivative'] if core
              else source['first_fixed_x_physical_time_derivative'])
            for component,parts in original.items():
                power=c.mpf(2 if component=='p' else 1)/2
                if set(time[component])!=set(parts):raise ValueError('Fixed-x time derivative contribution omitted')
                for label,row in parts.items():
                    result=time[component][label];count+=1;gamma=number(c,row['physical_lambda_exponent'])
                    if endpoints(number(c,result['constant_viscosity_power']))!=endpoints(power):raise ValueError('Constant-viscosity fixed-x time chain rule differs')
                    if core:
                        norm=endpoints(number(c,row['absolute_upper']))[1]
                        expected=None if not norm else c.ln(c.mpf(norm))+number(c,row['logLambda_term'])+number(c,row['amplitude_log_upper'])+gamma*q+power*lnnu
                    else:expected=None if row['exact_zero'] else number(c,row['log_absolute_upper'])+gamma*(q-lt/2)+power*lnnu
                    check_bound(c,result,expected)
    return count

def check_bound(c,result,expected):
    if result['exact_zero']!=(expected is None):raise ValueError('Actual physical mapping erased nonzero signed source')
    if expected is None:
        if result['log_absolute_upper'] is not None:raise ValueError('Exact zero physical row has a logarithm')
    else:
        bound=number(c,result['log_absolute_upper'])
        if endpoints(bound)!=endpoints(expected) or not all(mp.isfinite(v) for v in endpoints(bound)):
            raise ValueError('Actual lambda/nu physical contribution bound differs')

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedPhysicalVelocity(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[key] for key in INHERITED_GATES+GATES+OPEN):
        raise ValueError('Actual physical velocity producer/source/scope changed')
    proof=exact_modified_physical_velocity_theorem(field)
    if encode(pack(proof))!=raw['exact_actual_modified_physical_velocity_theorem']:
        raise ValueError('Actual source/dilation/coordinate proof changed')
    saved=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if set(saved)!={'physical_axis','physical_cartesian','source_transition','source_repair','source_heat'}:
        raise ValueError('Five Cartesian/axis/modified/Gamma actual queries required')
    counts={name:check_physical_query(field,query) for name,query in saved.items()}
    heat=saved['source_heat']['current_candidate_velocity_pressure_views']
    if not any(view['actual_physical_candidate']['region']=='heat_exterior' for view in heat):raise ValueError('Actual finite Gamma velocity query missing')
    for view in heat:
        if view['actual_physical_candidate']['region']=='heat_exterior':
            coordinate=number(field.ctx,view['actual_physical_candidate']['native_coordinate_enclosure'])
            if not all(mp.isfinite(v) for v in endpoints(coordinate)) or endpoints(coordinate)[0]<3:
                raise ValueError('Finite physical heat coordinate changed into infinite domain or zero tensor')
            native=view['current_native_velocity_pressure_query']
            if native['requested_coordinate_kind']!='original_native_coordinate':raise ValueError('Actual finite heat source route lost')
    fresh=field.log_radius('-inf','.137','-1.337','.537','2.3')
    fresh_count=check_physical_query(field,fresh)
    if not any(view['current_native_velocity_pressure_query']['source_pieces'][0]['complete_velocity_pressure_source']['core_axis_nonsingular_native_map']['axis_only']
      for view in fresh['current_candidate_velocity_pressure_views']):raise ValueError('Fresh physical axis did not use native analytic axis')
    for action in (lambda:field.cartesian(0,0,0,1),lambda:field.cartesian(mp.inf,0,0,0),lambda:field.log_radius(0,0,-1,0,0)):
        try:action()
        except ValueError:pass
        else:raise ValueError('Invalid physical point/time/viscosity admitted')
    result=dict(actual_five_defect_family_sha256=field.registry.family,implicit_source_sha256=field.registry.source,
      datum_enclosure_sha256=field.registry.datum_sha,finite_integer_N=field.N,
      modified_physical_velocity_definition_sha256=field.modified_physical_velocity_definition_sha256,
      new_physical_source_coordinate_and_viscosity_identities=len(proof['identities']),
      actual_physical_query_contribution_counts=counts,fresh_log_radius_axis_contributions=fresh_count,
      all33_native_coordinate_routes_bound_to_actual_global_physical_cover=True,
      exact_zero_native_layout_and_nonzero_heat_velocity_retained=True,
      physical_viscosity_not_restricted_to_one=True,actual_lambda_not_replaced_by_sqrt_tau=True,
      resolved_points_global_interfaces_energy_common_N_cones_recursion_full_NS_remain_open=True,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(INHERITED_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Physical Cartesian/source-correlated velocity/pressure, actual lambda and positive constant nu PASS; actual recursion/resolved NS open',flush=True)
    return result

if __name__=='__main__':run()
