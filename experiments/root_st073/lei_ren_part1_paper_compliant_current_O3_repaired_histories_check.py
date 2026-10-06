"""Accept the source-connected local repaired field and five-moment exit."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_O3_repaired_histories import (
    CurrentO3RepairedHistories,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_repaired_histories_operator import actual_repaired_history_theorem

def zeros(rows):return all(all(endpoints(v)==(0,0) for v in row.coefficients) for row in rows)

@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3RepairedHistories(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Actual repaired source manifest/scope differs')
    if encode(pack(actual_repaired_history_theorem()))!=raw['exact_repaired_partial_history_and_exit_source_theorem']:
        raise ValueError('Actual callable partial field/implicit exit source bridge differs')
    views={key:field.history('quiet_O3_power',(-1,1),q) for key,q in (
      ('whole_repair_band',(1,2)),('first_overlap_interior','1.201337'),('middle_swirl_interior','1.501337'),
      ('last_overlap_interior','1.801337'),('exact_exit',2))}
    if encode(pack(views))!=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes())):
        raise ValueError('Actual full-Z repaired source views differ')
    for view in views.values():
        if any(view[k] for k in OPEN) or not view['same_exact_control_vector_for_all_queries']:
            raise ValueError('Changed-source field cannot inherit a whole modified cone')
        if view['same_unique_implicit_repair_definition_sha256']!=field.repair.repair_definition_sha256:
            raise ValueError('Partial field and exit use different implicit controls')
        for name in ('m','k'):
            if [part['Pstar_power'] for part in view['modified_five_histories_in_original_normalized_units_Pstar_sectors'][name]]!=[0,1]:
                raise ValueError('Repaired M/K extra Pstar sector omitted')
    before=field.history('quiet_O3_power','.537',1);original=field.histories.history('quiet_O3_power','.537',1)
    for key in ('exact_modified_history_increment_rows','modified_absolute_pressure_over_Pstar2_ordinary_logR_rows',
                'modified_cylindrical_velocity_source_log_sectors'):
        if encode(pack(before[key]))!=encode(pack(original[key])):raise ValueError('Independent repair changes incoming source below support')
    exit_view=views['exact_exit'];delta=exit_view['exact_modified_history_increment_rows']
    if not exit_view['exact_functional_exit_reduction_applied'] or not all(zeros(rows) for rows in delta.values()):
        raise ValueError('Five actual source history increments do not close at the exit')
    if not zeros(exit_view['pressure_correction_from_own_cumulative_Cp_ordinary_logR_rows']):
        raise ValueError('Same-axis pressure defect did not close at the exit')
    velocity=exit_view['modified_cylindrical_velocity_source_log_sectors']
    if not zeros(velocity['axial'][0]['ordinary_logR_rows']) or not zeros(velocity['radial'][1]['ordinary_logR_rows']):
        raise ValueError('Own source axial/radial exit does not match original')
    if endpoints(exit_view['actual_positive_modulation_kinetic_history_retained'])[0]<=0:
        raise ValueError('Signed energy closure incorrectly erased nonzero kinetic source')
    W=exit_view['current_same_source_partial_bump_weights'];repair=field.repair;matrix=repair.matrix
    c=field.ctx;ci=repair.weights['centers'];mu=field.mu
    expected_full=dict(mass=[c.mpf(1)]*3,D=repair.weights['divided_axial_rows'],
      I=[repair.weights['H']['I']*c.exp(x/2) for x in ci],
      S=[repair.weights['H']['S']*c.exp((-c.mpf('.5')-mu)*x) for x in ci],
      Cp=[repair.weights['H']['Cp']*c.exp((-c.mpf('1.5')-mu)*x) for x in ci],
      cross=matrix['cross_weights'],energy=matrix['energy_weights'],pressure=matrix['pressure_weights'])
    if encode(pack(W))!=encode(pack(expected_full)):
        raise ValueError('Full partial weights differ from the very same defining implicit repair matrix')
    crossed=field.history('quiet_O3_power','.7',('1.99','2'))
    if crossed['exact_functional_exit_reduction_applied']:
        raise ValueError('Interval crossing the exit was incorrectly assigned zero throughout')
    for region in ('O2_buffer','O3_slope_mu'):
        q=11 if region=='O2_buffer' else 0
        new=field.history(region,'.537',q);old=field.histories.history(region,'.537',q)
        if encode(pack(new['modified_cylindrical_velocity_source_log_sectors']))!=encode(pack(old['modified_cylindrical_velocity_source_log_sectors'])):
            raise ValueError('Independent later repair changed the modulation source')
    result=dict(actual_five_defect_family_sha256=field.histories.family,
      implicit_source_sha256=field.histories.source,datum_enclosure_sha256=field.histories.datum_sha,
      parent_repair_definition_sha256=repair.repair_definition_sha256,
      modified_source_definition_sha256=field.modified_source_definition_sha256,
      finite_integer_N=field.N,partial_range_subdivisions=field.partial_cells,
      exact_source_partial_FTC_unit_jet_and_exit_identities=len(field.theorem['identities']),
      same_exact_implicit_controls_reproduce_full_partial_matrix_and_close_all_five=True,
      actual_same_axis_pressure_and_own_radial_exit_source_rows_match_original=True,
      full_Z_functional_exit_not_sample_fitting_or_interval_midpoint_substitution=True,
      nonzero_kinetic_terms_kept_while_signed_total_energy_defect_closes=True,
      before_support_and_original_O2_O3_source_preserved=True,
      exit_crossing_box_does_not_use_exact_zero_branch=True,
      modified_complete_tensor_and_heat_compatibility_not_yet_admitted=True,
      common_cone_frequency_not_admitted_by_repair_only_N=True,
      current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Own repaired partial five histories/pressure/radial and functional exact exit PASS; modified tensor/cone open',flush=True)
    return result

if __name__=='__main__':run()
