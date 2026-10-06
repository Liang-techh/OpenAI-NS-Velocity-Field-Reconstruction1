"""Check own mixed4 source/physical outputs without admitting global joins."""
import gzip
import json
import math
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure import (
    CurrentModifiedVelocityPressure,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,INHERITED_GATES,GATES,OPEN,
    QUERIES,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_operator import (
    exact_modified_velocity_pressure_theorem,modified_velocity_pressure_packet,LABELS,UT,UZ,UR,P)
from lei_ren_part1_paper_compliant_cartesian_field import INDICES
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

def check_physical(field,view):
    if len(view['physical_spatial_cartesian_mixed4'])!=35:raise ValueError('Incomplete actual Cartesian spatial4 inventory')
    expected={'x%d_y%d_z%d'%index for index in INDICES}
    if set(view['physical_spatial_cartesian_mixed4'])!=expected:raise ValueError('Wrong spatial multiindices')
    if any(view[k] for k in OPEN) or view['resolved_point_coefficient_values_available']:
        raise ValueError('Local modified packet cannot admit physical interfaces or resolved global field')
    packet=view['actual_own_velocity_pressure_mixed4_packet']
    if packet['signed_Pstar_sector_inventory']!={UT:[1],UZ:[1],UR:[0,1],P:[2]}:
        raise ValueError('Signed own velocity/absolute pressure unit sector lost')
    keys={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}
    for grid in packet['grids'].values():
        if set(grid)!=keys:raise ValueError('Own source mixed4 input incomplete')
    count=0
    for groups in list(view['physical_spatial_cartesian_mixed4'].values())+[view['first_fixed_x_physical_time_derivative']]:
        if set(groups)!={'ux','uy','uz','p'}:raise ValueError('Physical field labels incomplete')
        for component,parts in groups.items():
            if set(parts)!=({'ux':{UR,UT},'uy':{UR,UT},'uz':{UZ},'p':{P}}[component]):
                raise ValueError('Actual moving Cartesian basis labels changed')
            for row in parts.values():
                count+=1
                if endpoints(row['physical_lambda_exponent'])[1]>=0:raise ValueError('Invalid lambda lower-bound use')
                if not row['positive_source_exponentials_not_materialized'] or not row['original_source_factors_combined_before_enclosure']:
                    raise ValueError('Actual large units were capped or materialized')
                if row['exact_zero']:
                    if row['terms'] or row['log_absolute_upper'] is not None:raise ValueError('Nonzero source silently cleared')
                elif not all(mp.isfinite(v) for v in endpoints(row['log_absolute_upper'])):
                    raise ValueError('Nonfinite actual physical contribution bound')
                for term in row['terms']:
                    if not all(mp.isfinite(v) for v in endpoints(term['signed_coefficient'])+endpoints(term['log_absolute_upper'])):
                        raise ValueError('Nonfinite own signed source coefficient')
    if count!=216:raise ValueError('Actual physical spatial/time contribution groups incomplete')
    return count

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedVelocityPressure(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[k] for k in INHERITED_GATES+GATES+OPEN):
        raise ValueError('Own velocity/pressure manifest or source scope changed')
    if encode(pack(exact_modified_velocity_pressure_theorem(field)))!=raw['actual_modified_velocity_pressure_source_theorem']:
        raise ValueError('Own mixed4 source/absolute pressure/original physical consumer proof changed')
    saved=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()));counts={}
    for key,region,q in QUERIES:
        view=field.physical(region,(-1,1),q,('-3','-1'),None)
        if encode(pack(view))!=saved[key]:raise ValueError('Actual full-domain modified source physical view changed: '+key)
        counts[key]=check_physical(field,view)
    # A fresh query consumes the actual same-root partial source; verify
    # all four labels/75 signed input rows against those defining jets.
    fresh=field.physical('quiet_O3_power','.537','1.217337','-2.337','.337')
    check_physical(field,fresh);view=fresh['actual_own_history_pressure_source']
    packet=fresh['actual_own_velocity_pressure_mixed4_packet'];ordinary=0
    for name,label in LABELS.items():
        for j in range(5):
            for n in range(5-j):
                for term,part in zip(packet['grids'][label]['y%d_Z%d'%(j,n)],view['modified_cylindrical_velocity_source_log_sectors'][name]):
                    if endpoints(term[1])!=endpoints(part['ordinary_logR_rows'][j][n]*math.factorial(n)):
                        raise ValueError('Fresh own ordinary velocity row differs')
                    ordinary+=1
    for j in range(5):
        for n in range(5-j):
            value=packet['grids'][P]['y%d_Z%d'%(j,n)][0][1]
            if endpoints(value)!=endpoints(view['modified_absolute_pressure_over_Pstar2_ordinary_logR_rows'][j][n]*math.factorial(n)):
                raise ValueError('Fresh actual absolute P0+Cp mixed row differs')
            ordinary+=1
    if ordinary!=75:raise ValueError('Incomplete signed four-label mixed4 source')
    broken=dict(view,modified_cylindrical_velocity_source_log_sectors={
      **view['modified_cylindrical_velocity_source_log_sectors'],
      'radial':view['modified_cylindrical_velocity_source_log_sectors']['radial'][:1]})
    try:modified_velocity_pressure_packet(field,broken)
    except ValueError:pass
    else:raise ValueError('Adapter accepted a missing own radial Pstar1 source')
    try:field.physical('heat_exterior',0,4)
    except ValueError:pass
    else:raise ValueError('Local own packet incorrectly claims downstream/full33 coverage')
    result=dict(actual_five_defect_family_sha256=field.heat.registry.family,
      implicit_source_sha256=field.heat.registry.source,datum_enclosure_sha256=field.heat.registry.datum_sha,
      finite_integer_N=field.N,modified_velocity_pressure_definition_sha256=field.modified_velocity_pressure_definition_sha256,
      new_source_mixed4_and_physical_operator_identities=len(field.theorem['identities']),
      absolute_pressure_source_identities=len(field.theorem['actual_absolute_pressure_mixed4_source_theorem']['identities']),
      consumed_original_variable_velocity_identities=len(field.theorem['consumed_original_variable_velocity_mixed4_theorem']['identities']),
      complete_local_physical_query_groups=counts,fresh_actual_signed_mixed4_source_rows=ordinary,
      source_pressure_uses_Dy_to_order_j_minus_one_of_actual_Unew_square=True,
      own_spatial4_and_fixed_x_time1_bounds_do_not_admit_affected_interfaces=True,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(INHERITED_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Own modified four-label mixed4 pressure/velocity and Cartesian spatial4/time1 PASS; quantified interfaces/energy/cones/recursion/global NS open',flush=True)
    return result

if __name__=='__main__':run()
