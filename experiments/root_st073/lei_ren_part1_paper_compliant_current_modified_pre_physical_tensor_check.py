"""Accept local completed physical source, not cone/energy/global NS admission."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor import (
    CurrentModifiedPrePhysicalTensor,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,QUERIES,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_operator import (
    exact_modified_pre_physical_theorem,modified_physical_packet,compiled_modified_pre_physical_lift)
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import (
    raw_pre_velocity_rows,raw_pre_stress_rows,raw_pre_remainder_sectors,copy_jet)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

GROUPS=('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_divergence_mixed2',
  'completed_theta_theta_stress_mixed2','physical_three_component_remainder_mixed2',
  'physical_completed_stress_tensor_cartesian','physical_completed_stress_divergence_cartesian',
  'physical_remainder_cartesian','physical_momentum_residual_decomposition_cartesian')

def nonzero_rows(rows):
    return sorted(json.dumps(encode(pack(row)),sort_keys=True) for row in rows if not row['exact_zero'])
def compare_exit_physical_groups(actual,reference):
    # The extra split Pstar sectors are exact zero at exit. Canonicalize
    # only by removing those zero rows, preserving every signed other row.
    for name in GROUPS[:4]:
        if name=='completed_theta_theta_stress_mixed2':
            grids_actual=[actual[name]];grids_reference=[reference[name]]
        else:
            grids_actual=[actual[name][label] for label in reference[name]]
            grids_reference=[reference[name][label] for label in reference[name]]
        for a,b in zip(grids_actual,grids_reference):
            keys=set(next(iter(b.values())))
            if not all(set(grid)==keys for grid in a.values()):raise ValueError('Physical mixed row inventory differs at exit')
            for key in keys:
                if nonzero_rows([grid[key] for grid in a.values()])!=nonzero_rows([grid[key] for grid in b.values()]):
                    raise ValueError('Actual physical completed source differs from original at exit: '+name+'.'+key)
    for name in GROUPS[4:]:
        if set(actual[name])!=set(reference[name]):raise ValueError('Cartesian source component inventory differs')
        for component in reference[name]:
            if nonzero_rows(actual[name][component])!=nonzero_rows(reference[name][component]):
                raise ValueError('Cartesian actual/source exit mismatch: '+name+'.'+component)

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedPrePhysicalTensor(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[key] for key in GATES+OPEN):raise ValueError('Modified physical source manifest/scope changed')
    if encode(pack(exact_modified_pre_physical_theorem(field)))!=raw['exact_actual_modified_physical_source_theorem']:
        raise ValueError('Original operator/geometry or own M physical divergence source changed')
    _,adaptation=compiled_modified_pre_physical_lift()
    if encode(pack(adaptation))!=raw['source_remainder_only_original_physical_lift_adaptation']:
        raise ValueError('More than original remainder input changed in physical lift')
    views={key:field.tensor(region,(-1,1),q,lt,theta,nu) for key,region,q,lt,theta,nu in QUERIES}
    if encode(pack(views))!=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes())):
        raise ValueError('Whole/full-Z modified completed physical views differ')
    for view in views.values():
        if any(view[key] for key in OPEN) or not view['local_physical_source_bounds_not_resolved_physical_point_values']:
            raise ValueError('Local completed source cannot inherit global cone/NS/energy scope')
        if endpoints(view['exact_completed_tensor_radial_divergence'])!=(0,0) or endpoints(view['exact_physical_divergence'])!=(0,0):
            raise ValueError('Actual tensor/own-moment velocity divergence identity omitted')
        packet=view['current_actual_modified_physical_packet'];signed=view['current_actual_modified_signed_source']
        if packet['parent_modified_stress_definition_sha256']!=field.source.modified_stress_definition_sha256:
            raise ValueError('Physical source uses foreign signed columns/remainder')
        if encode(pack(packet['actual_modified_source_remainder_sectors']))!=encode(pack(signed['full_signed_source_remainder2_sectors'])):
            raise ValueError('Own full signed remainder not supplied to physical mapper')
        if [len(view['physical_cylindrical_stress_mixed3'][key]) for key in ('theta','axial')]!=[7,8]:
            raise ValueError('Physical source stress sectors omitted')
        if len(view['completed_theta_theta_stress_mixed2'])!=8 or [len(view['physical_three_component_remainder_mixed2'][key]) for key in ('radial','theta','axial')]!=[9,1,1]:
            raise ValueError('Complete physical diagonal/remainder sector inventory changed')
        for parts in view['physical_cylindrical_stress_mixed3'].values():
            if not all(len(grid)==10 for grid in parts.values()):raise ValueError('Physical stress mixed3 source rows missing')
        for group in ('physical_cylindrical_stress_divergence_mixed2','physical_three_component_remainder_mixed2'):
            for parts in view[group].values():
                if not all(len(grid)==6 for grid in parts.values()):raise ValueError('Physical divergence/remainder mixed2 rows missing')
        if not all(len(grid)==6 for grid in view['completed_theta_theta_stress_mixed2'].values()):raise ValueError('Completed diagonal mixed2 rows missing')
        for component,count in (('xx',15),('xy',15),('xz',8),('yy',15),('yz',8),('zz',1)):
            if len(view['physical_completed_stress_tensor_cartesian'][component])!=count:raise ValueError('Completed Cartesian sector inventory changed')
    # Compare the physical q=2 source to independently recovered original
    # raw histories/velocity/pressure and their original coefficient programs.
    exit_view=views['exact_exit'];signed=exit_view['current_actual_modified_signed_source'];source=signed['actual_source']
    c=field.ctx;z=IntervalTaylor.variable(c,c.mpf((-1,1)),5)
    velocity=raw_pre_velocity_rows(c,source['current_original_pre_source'])
    histories={key:[copy_jet(c,row) for row in rows] for key,rows in source['current_original_pre_source']['actual_normalized_primitive_y_derivative_axial5'].items()}
    original=dict(signed,
      full_signed_paper_theta_axial_stress3_sectors=raw_pre_stress_rows(c,field.delta,z,velocity['theta'],velocity['axial'],histories,source['original_absolute_pressure_over_Pstar2_ordinary_logR_rows']),
      full_signed_source_remainder2_sectors=raw_pre_remainder_sectors(c,field.delta,z,velocity))
    original_packet,_=modified_physical_packet(field,original)
    reference=field.lift(c,original_packet,field.delta,None,c.mpf('-1'),None,c.mpf(1))
    compare_exit_physical_groups(exit_view,reference)
    fresh=field.tensor('quiet_O3_power','.537','1.223337','-2.337','.337','.8')
    if encode(pack(fresh['current_actual_modified_signed_source']))!=encode(pack(field.source.stress('quiet_O3_power','.537','1.223337'))):
        raise ValueError('Fresh physical source does not consume the same repaired field')
    result=dict(actual_five_defect_family_sha256=field.source.source.histories.family,
      implicit_source_sha256=field.source.source.histories.source,datum_enclosure_sha256=field.source.source.histories.datum_sha,
      parent_modified_source_definition_sha256=field.source.source.modified_source_definition_sha256,
      parent_modified_stress_definition_sha256=field.source.modified_stress_definition_sha256,
      modified_physical_tensor_definition_sha256=field.modified_physical_tensor_definition_sha256,finite_integer_N=field.N,
      new_exact_physical_source_packet_radius_and_own_divergence_identities=len(field.theorem['identities']),
      consumed_original_full_physical_operator_identities=len(field.theorem['consumed_original_arbitrary_source_full_physical_operator']['identities']),
      actual_physical_stress_sector_counts=dict(theta=7,axial=8),completed_diagonal_sector_count=8,
      actual_physical_remainder_sector_counts=dict(radial=9,theta=1,axial=1),
      only_actual_remainder_input_changes_in_original_physical_lift=True,
      completed_diagonal_radial_divergence_and_own_velocity_divergence_source_certified=True,
      q2_full_Z_completed_physical_source_equals_independent_original_source=True,
      signed_source_sectors_and_true_radius_Pstar_factors_retained=True,
      local_source_bounds_are_not_resolved_physical_points_or_complete_33_chart_dispatch=True,
      common_cone_frequency_global_admissibility_energy_and_corrected_NS_remain_open=True,
      current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual local modified completed physical tensor/divergence/remainder PASS; modified global cone/energy/NS open',flush=True)
    return result

if __name__=='__main__':run()
