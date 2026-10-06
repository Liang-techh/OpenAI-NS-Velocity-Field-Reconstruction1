"""Accept actual signed local columns/remainder, not a completed cone tensor."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_modified_pre_stress import (
    CurrentModifiedPreStress,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,QUERIES,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_pre_stress_operator import exact_modified_pre_stress_theorem
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import (
    raw_pre_velocity_rows,raw_pre_stress_rows,raw_pre_remainder_sectors,copy_jet)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

def zeros(rows):return all(all(endpoints(value)==(0,0) for value in row.coefficients) for row in rows)
def same_retained_rows(left,right,order):
    if len(left)!=len(right):return False
    for j,(a,b) in enumerate(zip(left,right)):
        common=min(a.order,b.order)
        if common<order-j:raise ValueError('Required mixed stress/remainder source order is missing')
        if encode(pack(a.truncate(common)))!=encode(pack(b.truncate(common))):return False
    return True
def exact_exit_groups(actual,original,rows_key):
    for label,parts in original.items():
        for name,part in parts.items():
            chosen={key.rsplit('_P',1)[1]:value for key,value in actual[label].items() if key.rsplit('_P',1)[0]==name}
            power=str(part['mode'][1])
            if zeros(part[rows_key]):
                if not all(zeros(value[rows_key]) for value in chosen.values()):
                    raise ValueError('Exact exit zero original stress/remainder source changed: '+label+'.'+name)
            else:
                if power not in chosen or not same_retained_rows(chosen[power][rows_key],part[rows_key],3 if rows_key=='full_derivative_rows' else 2):
                    raise ValueError('Actual signed source does not equal original at exit: '+label+'.'+name)
                if any(not zeros(value[rows_key]) for key,value in chosen.items() if key!=power):
                    raise ValueError('Nonzero foreign Pstar source sector at exact exit')

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedPreStress(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Actual modified stress manifest/scope differs')
    if encode(pack(exact_modified_pre_stress_theorem()))!=raw['exact_signed_modified_stress_remainder_source_theorem']:
        raise ValueError('Original variable source paper/remainder bridge changed')
    views={key:field.stress(region,(-1,1),q) for key,region,q in QUERIES}
    if encode(pack(views))!=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes())):
        raise ValueError('Actual full-Z signed modified stress/remainder views differ')
    for view in views.values():
        if any(view[k] for k in OPEN):raise ValueError('Local modified source cannot inherit completed tensor/cone')
        if view['modified_source_definition_sha256']!=field.source.modified_source_definition_sha256:
            raise ValueError('Stress and remainder use a foreign source')
        columns=view['full_signed_paper_theta_axial_stress3_sectors'];remainder=view['full_signed_source_remainder2_sectors']
        if [len(columns[k]) for k in ('theta','axial')]!=[7,8] or [len(remainder[k]) for k in ('radial','theta','axial')]!=[9,1,1]:
            raise ValueError('Full variable source term inventory differs')
        for parts in columns.values():
            if not all(len(part['full_derivative_rows'])==4 for part in parts.values()):raise ValueError('Ordinary stress3 rows missing')
        for parts in remainder.values():
            if not all(len(part['rows'])==3 for part in parts.values()):raise ValueError('Ordinary remainder2 rows missing')
        if [remainder['radial']['nonlinear_transport_P'+str(power)]['mode'][1] for power in (0,1,2)]!=[0,1,2]:
            raise ValueError('Radial nonlinear incoming/cross/square sectors missing')
    exit_view=views['exact_exit'];source=exit_view['actual_source'];c=field.ctx;z=IntervalTaylor.variable(c,c.mpf((-1,1)),5)
    velocity=raw_pre_velocity_rows(c,source['current_original_pre_source'])
    histories={key:[copy_jet(c,row) for row in rows] for key,rows in source['current_original_pre_source']['actual_normalized_primitive_y_derivative_axial5'].items()}
    original_columns=raw_pre_stress_rows(c,field.delta,z,velocity['theta'],velocity['axial'],histories,source['original_absolute_pressure_over_Pstar2_ordinary_logR_rows'])
    original_remainder=raw_pre_remainder_sectors(c,field.delta,z,velocity)
    exact_exit_groups(exit_view['full_signed_paper_theta_axial_stress3_sectors'],original_columns,'full_derivative_rows')
    exact_exit_groups(exit_view['full_signed_source_remainder2_sectors'],original_remainder,'rows')
    fresh=field.stress('quiet_O3_power','.537','1.223337')
    if encode(pack(fresh['actual_source']))!=encode(pack(field.source.history('quiet_O3_power','.537','1.223337'))):
        raise ValueError('Fresh signed columns do not consume the actual repaired field')
    result=dict(actual_five_defect_family_sha256=field.source.histories.family,
      implicit_source_sha256=field.source.histories.source,datum_enclosure_sha256=field.source.histories.datum_sha,
      parent_modified_source_definition_sha256=field.source.modified_source_definition_sha256,
      modified_stress_definition_sha256=field.modified_stress_definition_sha256,finite_integer_N=field.N,
      new_exact_signed_source_stress_remainder_identities=len(field.theorem['identities']),
      consumed_original_variable_full_paper_stress_identities=len(field.theorem['consumed_arbitrary_variable_original_full_paper_stress_theorem']['identities']),
      actual_signed_stress_sector_counts=dict(theta=7,axial=8),actual_signed_remainder_sector_counts=dict(radial=9,theta=1,axial=1),
      full_actual_source_and_original_paper_programs_bound=True,
      full_Z_q2_signed_columns_and_remainder_equal_original_source=True,
      exact_exit_comparison_uses_all_common_retained_axial_rows_without_discarding_higher_source_rows=True,
      all_Pstar0_Pstar1_Pstar2_radial_cross_and_square_terms_retained=True,
      physical_Uz_Pstar_not_normalized_Vhat_used=True,
      inherited_original_cones_not_admission_of_modified_source=True,
      own_complete_energy_tensor_divergence_and_physical_joins_remain_open=True,
      current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual modified signed paper stress15/remainder11 sectors PASS; own energy/complete tensor/cone open',flush=True)
    return result

if __name__=='__main__':run()
