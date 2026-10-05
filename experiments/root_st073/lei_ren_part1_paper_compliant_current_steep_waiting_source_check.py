"""Focused current angular -> original O7 source admission."""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_steep_waiting_source import (
    CurrentSteepWaitingSourceAssembly,NEW_CHARTS,CHARTS,GATE,PROVED,THEOREM,
    UNIFORM,SCOPES,OPEN,HERE,PREFIX,sha)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_steep_waiting_source.json'


def run():
    raw=json.loads((HERE/NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Current steep/waiting dependency changed: '+name)
    field=CurrentSteepWaitingSourceAssembly(require_checked=False)
    omitted={'whole_current_steep_waiting_source_maps','current_angular_terminal'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:
        raise ValueError('Fresh current steep/waiting source manifest differs')
    binding=raw['current_steep_waiting_source_bindings']
    kernels=binding['exact_current_kernel_and_parameter_source_bindings']
    if (not binding['passed'] or not all(binding['current_defining_object_graph'].values())
            or not all(binding['current_complete_future_factor_identities'].values())
            or not binding['current_waiting_root_and_heat_atoms_not_read_from_saved_packets']
            or not binding['current_absolute_pressure_and_nonzero_angular_histories_retained']
            or not kernels['passed'] or not kernels['same_exact_original_logistic_sigma_function']
            or not kernels['independent_cell_enclosures_do_not_define_different_integrals']
            or not all(kernels['actual_shared_kernel_functional_identities'].values())
            or not binding['actual_current_angular_s0_future_support_bindings']['passed']
            or binding['source_caps_used_as_defining_field_values']
            or binding['interval_overlap_used_as_join_proof']):
        raise ValueError('Actual current O7 source ownership or full future decomposition omitted')
    if (set(raw['ordered_current_chart_registry'])!=set(CHARTS)
            or raw['current_downstream_chart_owner_count']!=27
            or raw[GATE] or not raw[PROVED] or raw[UNIFORM] or raw['full_pulse_C4_installed']
            or raw['current_steep_waiting_physical_owner_installed']
            or any(raw[k] for k in SCOPES+OPEN)):
        raise ValueError('Current O7 source scope widened')
    c=field.steep.ctx
    def interval(row):
        value=read_interval(c,row);lo,hi=endpoints(value)
        if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:raise ValueError('Nonfinite or unordered O7 source row')
        return value
    def exact_zero(row):
        if endpoints(interval(row))!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Inherited actual meridional zero lost')
    indices={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
    maps=raw['whole_current_steep_waiting_source_maps'];counts={};zeros=0
    if set(maps)!=set(NEW_CHARTS):raise ValueError('Original complete O7 source chart coverage omitted')
    for chart,views in maps.items():
        if set(views)!={'inlet','whole_domain','exit'} or field.provider(chart) is not field.steep:
            raise ValueError('Original O7 whole-domain coverage or common provider differs')
        for name,packet in views.items():
            if (packet['chart']!=chart or packet[GATE] or not packet[PROVED] or packet[UNIFORM]
                    or packet['full_pulse_C4_installed'] or packet['current_steep_waiting_physical_owner_installed']
                    or packet['datum_enclosure_sha256']!=field.datum_sha or any(packet[k] for k in SCOPES+OPEN)):
                raise ValueError('Current O7 packet source or unfinished scope differs')
            source=packet['source_packet']
            for flag in ('current_steep_waiting_source_inputs_used','original_infinite_heat_tail_and_epsilon_atoms_retained',
                    'same_actual_angular_terminal_histories_retained','exact_Gamma_and_both_epsilon_atoms_retained',
                    'no_forward_subtraction_of_unrelated_long_future_energy'):
                if not source[flag]:raise ValueError('Original complete O7 input omitted: '+flag)
            for flag in ('full_pulse_C4_installed','full_outer_C4_certified','physical_energy_integral_certified',
                    'whole_outer_cone_certified','temporal_recursion'):
                if source[flag]:raise ValueError('O7 source scope promoted: '+flag)
            coord='offset' if chart in ('steep_entry','steep_exit') else 'phase'
            target={'inlet':(0,0),'whole_domain':(0,1),'exit':(1,1)}[name]
            if endpoints(interval(source['coordinate'][coord]))!=tuple(map(mp.mpf,target)):
                raise ValueError('Original O7 coordinate domain differs')
            for key in ('theta_over_Ev0_Taylor','angular_Taylor','energy_Taylor','pressure_over_Pstar_squared_Taylor'):
                rows=source[key]['coefficients']
                if len(rows)!=6:raise ValueError('Whole-Z O7 axial5 source omitted')
                for row in rows:interval(row)
            floor=endpoints(interval(source['energy_Taylor']['coefficients'][0]))[0]
            if floor<=0 or (chart=='steep_power' and floor<mp.mpf('.25')):
                raise ValueError('Positive complete future / steep-power floor lost')
            if chart in ('steep_entry','steep_exit') and name!='whole_domain':
                target=mp.mpf(0) if name=='inlet' else mp.mpf('.5')
                if endpoints(interval(source['original_transition_kernels']['J']))!=(target,target):
                    raise ValueError('Original exact J transition endpoint lost')
            grid=source['physical_mixed_derivatives_total_order_le4'];count=0
            if set(grid)!={UZ,UT,UR,P}:raise ValueError('O7 velocity/pressure component omitted')
            for label,rows in grid.items():
                if set(rows)!=indices:raise ValueError('Whole O7 mixed4 index omitted')
                for row in rows.values():
                    interval(row);count+=1
                    if label in (UZ,UR):exact_zero(row);zeros+=1
            counts[chart+'_'+name]=count
    # Only the changed live acquisition is exercised anew. Original directed
    # kernel and variable-rate fixtures are retained in the source-hashed
    # canonical theorem instead of rerun.
    fresh=field.steep.data('0.271')
    if not fresh['current_live_Gamma_C5_and_admitted_C4_prefix_used'] or endpoints(fresh['H'][0])[0]<=0:
        raise ValueError('Fresh unsaved current Gamma/future acquisition failed')
    for chart in NEW_CHARTS:
        point=field.evaluate(chart,'0.271',0)['source_packet']
        if endpoints(point['energy_Taylor'][0])[0]<=0:raise ValueError('Fresh positive O7 source lost')
        for value in (-1,2):
            try:field.evaluate(chart,0,value)
            except ValueError:pass
            else:raise ValueError('Original O7 domain widened')
    try:field.evaluate('waiting',[-1.001,1],0)
    except ValueError:pass
    else:raise ValueError('Original axial domain widened')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_downstream_chart_owners_checked=27,
        newly_checked_current_source_charts=list(NEW_CHARTS),current_steep_waiting_mixed4_rows_checked=counts,
        total_new_steep_waiting_source_rows_checked=sum(counts.values()),exact_zero_meridional_rows_checked=zeros,
        current_complete_future_factor_identities_checked=len(binding['current_complete_future_factor_identities']),
        current_shared_exact_kernel_functional_identities_checked=len(kernels['actual_shared_kernel_functional_identities']),
        retained_canonical_functional_identity_count=len(field.functional_proof),retained_canonical_endpoint_theorem=THEOREM,
        fresh_unsaved_Z_current_source_acquisition_checked=True,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,
        **{GATE:True,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
        current_steep_waiting_physical_owner_installed=False,
        **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current27 source owners: current angular -> steep/waiting, complete heat future and 720 new mixed4 rows',flush=True)
    return result


if __name__=='__main__':run()
