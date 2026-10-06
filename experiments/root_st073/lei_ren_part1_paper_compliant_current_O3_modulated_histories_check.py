"""Independent scope/units and continuous source-history acceptance checks."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_O3_modulated_histories import (
    CurrentO3ModulatedHistories,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import SHAPES

def zeros(rows):
    return all(all(endpoints(v)==(0,0) for v in row.coefficients) for row in rows)

@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3ModulatedHistories(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Modified source producer/scope differs')
    if not field.theorem['passed'] or not all(field.theorem['identities'].values()):raise ValueError('Exact cumulative source theorem required')
    if not field.bounds['all_finite_integer_N_at_least1'] or field.bounds['finite_uniform_N_admissible_cone_choice']:
        raise ValueError('Uniform error constants are not an admitted frequency')
    views={key:field.history(region,(-1,1),coordinate) for key,region,coordinate in (
      ('buffer_whole_support','O2_buffer',(9,11)),
      ('transition_whole','O3_slope_mu',(0,1)),
      ('quiet_power_repair_preview','quiet_O3_power',(1,2)))}
    published=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if encode(pack(views))!=published:raise ValueError('Actual changed full source views differ')
    c=field.ctx;K=field.bounds['untransported_abs_integral_times_N_upper']
    for name,value in K.items():
        if endpoints(value)[0]<=0:raise ValueError('Source error constant must be positive: '+name)
    for view in views.values():
        if any(view[k] for k in OPEN) or not view['terminal_five_moment_restore_not_performed']:
            raise ValueError('Unrepaired modified source cannot admit cone or exterior matching')
        histories=view['modified_five_histories_in_original_normalized_units_Pstar_sectors']
        if set(histories)!=set(SHAPES):raise ValueError('All five own modified histories required')
        for name in ('m','k'):
            if [part['Pstar_power'] for part in histories[name]]!=[0,1]:
                raise ValueError('Extra actual Pstar in M/K omitted')
        velocity=view['modified_cylindrical_velocity_source_log_sectors']
        if velocity['axial'][0]['Pstar_power']!=1 or [part['Pstar_power'] for part in velocity['radial']]!=[0,1]:
            raise ValueError('Actual axial/radial Pstar sectors lost')
        t=view['actual_shared_logR_offset']
        for name,(_,_,rate) in SHAPES.items():
            transported=K[name]/field.N*(c.exp(-c.mpf(str(rate))*t) if rate else 1)
            actual=view['signed_actual_cumulative_scalar_enclosures'][name]
            if max(abs(v) for v in endpoints(actual))>endpoints(transported)[1]:
                raise ValueError('Signed range integral exceeds source uniform bound: '+name)
        if endpoints(t)[0]>=0 and endpoints(view['signed_actual_cumulative_scalar_enclosures']['e_kinetic'])[0]<=0:
            raise ValueError('Exact nonzero kinetic mass from the unit buffer plateau omitted')
    left=field.history('O2_buffer','.537',11);right=field.history('O3_slope_mu','.537',0)
    for key in ('exact_modified_history_increment_rows','modified_absolute_pressure_over_Pstar2_ordinary_logR_rows',
                'modified_cylindrical_velocity_source_log_sectors'):
        if encode(pack(left[key]))!=encode(pack(right[key])):raise ValueError('Common modified O2/O3 source seam differs: '+key)
    before=field.history('O2_buffer','.537',9)
    if not all(zeros(rows) for rows in before['exact_modified_history_increment_rows'].values()):
        raise ValueError('Increment is not exactly zero before support')
    if not zeros(before['modified_cylindrical_velocity_source_log_sectors']['radial'][1]['ordinary_logR_rows']):
        raise ValueError('Original radial source altered below support')
    after=field.history('quiet_O3_power','.537','1.337')
    if not zeros(after['modified_cylindrical_velocity_source_log_sectors']['axial'][0]['ordinary_logR_rows']):
        raise ValueError('Modified axial profile did not return to its original quiet value')
    if not zeros(after['pressure_correction_from_own_cumulative_Cp_ordinary_logR_rows'][1:]):
        raise ValueError('Pressure defect must persist as an X-constant after support')
    if endpoints(after['signed_actual_cumulative_scalar_enclosures']['e_kinetic'])[0]<=0:
        raise ValueError('Transported positive kinetic source was reset')
    fresh=field.history('O3_slope_mu','.7','.1337000000001337')
    if endpoints(fresh['pressure_correction_from_own_cumulative_Cp_ordinary_logR_rows'][1][0])[1]>=0:
        raise ValueError('Changed centrifugal pressure derivative was copied from the old source')
    for N,cells in ((0,8),(True,8),(1.5,8),(11,3),(11,512)):
        try:CurrentO3ModulatedHistories(candidate=field.candidate,N=N,cells=cells,require_checked=False)
        except ValueError:pass
        else:raise ValueError('Invalid frequency/range source accepted')
    for region,z,v in (('foreign',0,0),('O3_slope_mu',1.1,0),('quiet_O3_power',0,3)):
        try:field.history(region,z,v)
        except ValueError:pass
        else:raise ValueError('Foreign source domain accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
      datum_enclosure_sha256=field.datum_sha,modified_source_definition_sha256=field.modified_source_definition_sha256,
      finite_integer_N=field.N,
      signed_range_subdivisions_per_fixed_piece=field.cells,
      exact_source_history_transport_pressure_radial_and_uniform_axial_identities=len(field.theorem['identities']),
      same_own_integrated_moments_and_common_axis_pressure_radial_recovery=True,
      extra_Pstar_M_K_and_axial_radial_source_modes_retained=True,
      whole_Z_shapes_and_continuous_axial_O1_over_N_bounds_certified=True,
      positive_kinetic_buffer_mass_retained_for_all_finite_integer_N=True,
      signed_cell_ranges_not_phase_samples_or_cycle_resolution=True,
      unrepaired_outgoing_five_defects_and_pressure_are_not_reset=True,
      original_source_cone_counts_do_not_admit_modified_source=True,
      current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Modified own five histories, same-axis pressure, radial recovery and uniform bounds PASS; repair/cone open',flush=True)
    return result

if __name__=='__main__':run()
