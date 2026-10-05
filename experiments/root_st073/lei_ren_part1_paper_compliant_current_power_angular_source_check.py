"""Focused current flatten -> O6 power/angular source acceptance.

Retain the existing algorithm theorem; newly check source acquisition,
complete future decomposition and original whole-domain source packets.
"""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_power_angular_source import (
    CurrentPowerAngularSourceAssembly, NEW_CHARTS, CHARTS, GATE, THEOREM,
    UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_power_angular_source.json'


def run():
    raw=json.loads((HERE/NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Current power/angular dependency changed: '+name)
    field=CurrentPowerAngularSourceAssembly(require_checked=False)
    omitted={'whole_current_power_angular_source_maps','current_angular_support_crossings','current_flatten_exit'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:
        raise ValueError('Fresh current power/angular source manifest differs')
    binding=raw['current_power_angular_source_bindings']
    prefix=binding['current_prefix_and_future_decomposition']
    radius=binding['original_radius_functional_join']
    parameter=binding['current_native_parameter_source_bridge']
    if (not binding['passed'] or not all(binding['current_source_object_graph'].values())
            or not all(binding['actual_C4_prefix_parameter_object_graph'].values())
            or not parameter['passed'] or not parameter['same_exact_parameter_definition_and_source_hashes']
            or not all(parameter['native_current_parameter_object_and_copy_graph'].values())
            or not prefix['passed'] or not all(prefix['actual_complete_future_factor_identities'].values())
            or not prefix['first_five_scaled_angular_and_Gamma_rows_copied_before_appending_fifth']
            or not prefix['physical_prefix_uses_same_exact_repair_scale']
            or not prefix['first_five_complete_future_rows_are_actual_C4_prefix']
            or not radius['passed'] or not radius['flatten100_equals_power0']
            or not radius['power1_equals_angular_minus4']
            or not binding['complete_post_future_decomposition_retained']
            or binding['source_caps_used_as_defining_field_values']
            or binding['interval_overlap_used_as_join_proof']):
        raise ValueError('Current source/prefix/decomposition functional binding omitted')
    if (set(raw['ordered_current_chart_registry'])!=set(CHARTS)
            or raw['current_downstream_chart_owner_count']!=23
            or raw[GATE] or raw[UNIFORM] or raw['full_pulse_C4_installed']
            or raw['current_power_angular_physical_owner_installed']
            or raw['all_profile_source_charts_callable'] or any(raw[k] for k in SCOPES+OPEN)):
        raise ValueError('Current power/angular source scope widened')
    if any(field.provider(chart) is not field.outer for chart in NEW_CHARTS):
        raise ValueError('Both new charts must use the same current provider')
    c=field.outer.ctx
    def interval(row):
        value=read_interval(c,row);lo,hi=endpoints(value)
        if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:
            raise ValueError('Nonfinite or unordered current power/angular bounds')
        return value
    def exact_zero(row):
        if endpoints(interval(row))!=(mp.mpf(0),mp.mpf(0)):
            raise ValueError('Current terminal zero meridional source lost')
    indices={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
    maps=raw['whole_current_power_angular_source_maps']
    expected={'power_inlet':('outer_power','phase',0),'power_domain':('outer_power','phase',[0,1]),
        'power_exit':('outer_power','phase',1),'angular_inlet':('outer_angular','offset',-4),
        'angular_domain':('outer_angular','offset',[-4,0]),'angular_exit':('outer_angular','offset',0)}
    if set(maps)!=set(expected) or len(raw['current_angular_support_crossings'])!=4:
        raise ValueError('Whole original power/angular coverage omitted')
    rows_checked={};zero_rows=0
    for name,packet in list(maps.items())+[(f'angular_support_crossing_{i}',p)
            for i,p in enumerate(raw['current_angular_support_crossings'])]:
        if (packet[GATE] or packet[UNIFORM] or packet['full_pulse_C4_installed']
                or packet['current_power_angular_physical_owner_installed']
                or not packet['current_power_angular_source_functional_joins_proved']
                or packet['datum_enclosure_sha256']!=field.datum_sha
                or any(packet[k] for k in SCOPES+OPEN)):
            raise ValueError('Current packet source/scope differs')
        source=packet['source_packet']
        if name in expected:
            chart,key,target=expected[name]
            wanted=target if isinstance(target,list) else [target,target]
            if packet['chart']!=chart or endpoints(interval(source['coordinate'][key]))!=tuple(map(mp.mpf,wanted)):
                raise ValueError('Original power/angular coordinate coverage differs')
        for flag in ('current_power_angular_source_inputs_used','entire_infinite_heat_tail_and_epsilon_atoms_retained',
                'actual_terminal_zero_linear_and_radial_histories_inherited',
                'original_pressure_datum_and_positive_complete_future_energy_preserved',
                'correlated_q_squared_factors_cancelled_before_enclosure'):
            if not source[flag]:raise ValueError('Actual complete future/history source omitted: '+flag)
        for flag in ('full_pulse_C4_installed','full_outer_C4_certified','physical_energy_integral_certified',
                'whole_outer_cone_certified','temporal_recursion'):
            if source[flag]:raise ValueError('Current packet scope promoted: '+flag)
        for key in ('energy_Taylor','angular_Taylor','theta_over_Ev0_Taylor','pressure_over_Pstar_squared_Taylor'):
            coefficients=source[key]['coefficients']
            if len(coefficients)!=6:raise ValueError('Current axial5 source row omitted')
            for row in coefficients:interval(row)
        if endpoints(interval(source['energy_Taylor']['coefficients'][0]))[0]<=0:
            raise ValueError('Complete current future energy positivity lost')
        histories=source['current_native_Rv_five_histories']
        for key,rows in histories.items():
            if len(rows['coefficients'])!=6:raise ValueError('Actual native Rv history axial5 omitted')
            for row in rows['coefficients']:
                interval(row)
                if key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):exact_zero(row)
        for key in ('Mtheta_over_sqrt2_R_3half_Utheta','Mztheta_over_R_Utheta_squared'):
            if endpoints(interval(histories[key]['coefficients'][0]))[0]<=0:
                raise ValueError('Nonzero actual Rv angular/energy history lost')
        grid=source['physical_mixed_derivatives_total_order_le4']
        if set(grid)!={UZ,UT,UR,P}:raise ValueError('Current velocity/pressure component omitted')
        count=0
        for label,rows in grid.items():
            if set(rows)!=indices:raise ValueError('Current ordinary mixed4 index omitted')
            for row in rows.values():
                interval(row);count+=1
                if label in (UZ,UR):exact_zero(row);zero_rows+=1
        rows_checked[name]=count
    # Reuse the already proved canonical algorithms. Exercise only the newly
    # changed live source acquisition at a fresh axial argument.
    fresh=field.outer.data('0.271')
    if (not fresh['live_angular_C5_with_admitted_C4_prefix']
            or not fresh['post_scalar_from_actual_prefix_future_object']
            or fresh['post'].order!=5 or fresh['complete_future_energy_input'].order!=5
            or endpoints(fresh['post'][0])[0]<=0):
        raise ValueError('Fresh unsaved current source acquisition failed')
    for chart,coordinate in (('outer_power',0),('outer_angular',-4)):
        packet=field.evaluate(chart,'0.271',coordinate)['source_packet']
        if endpoints(packet['energy_Taylor'][0])[0]<=0:raise ValueError('Fresh positive energy lost')
    for chart,Z,value in (('outer_power',[-1.001,1],0),('outer_power',0,-1),
            ('outer_power',0,2),('outer_angular',0,-5),('outer_angular',0,1)):
        try:field.evaluate(chart,Z,value)
        except ValueError:pass
        else:raise ValueError('Original current power/angular domain widened')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_downstream_chart_owners_checked=23,
        newly_checked_current_source_charts=list(NEW_CHARTS),current_power_angular_mixed4_rows_checked=rows_checked,
        total_new_power_angular_source_rows_checked=sum(rows_checked.values()),
        exact_zero_meridional_source_rows_checked=zero_rows,
        actual_prefix_future_factor_identities_checked=len(prefix['actual_complete_future_factor_identities']),
        retained_canonical_functional_identity_count=len(field.functional_proof),retained_canonical_endpoint_theorem=THEOREM,
        current_C4_prefix_source_and_parameter_bridge_checked=True,
        current_complete_post_future_decomposition_checked=True,current_radius_functional_joins_checked=True,
        fresh_unsaved_Z_current_source_acquisition_checked=True,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,
        **{GATE:True,UNIFORM:False},full_pulse_C4_installed=False,
        current_power_angular_physical_owner_installed=False,
        **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current flatten -> power/angular: twenty-three profile owners, C4/C5 source bridge and 600 new mixed4 rows',flush=True)
    return result


if __name__=='__main__':run()
