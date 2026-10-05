"""Focused current-object pulse-terminal / O5 flatten source acceptance."""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_pulse_flatten_source import (
    CurrentPulseFlattenSourceAssembly, CHARTS, GATE, PROVED, THEOREM,
    UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_pulse_flatten_source.json'


def run():
    raw=json.loads((HERE/NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Current flatten dependency changed: '+name)
    field=CurrentPulseFlattenSourceAssembly(require_checked=False)
    manifest=encode(pack(field.manifest()))
    omitted={'current_native_terminal_packet','current_flatten_source_packets'}
    if manifest!={k:v for k,v in raw.items() if k not in omitted}:
        raise ValueError('Fresh current-object defining manifest differs')
    binding=raw['current_pulse_terminal_flatten_source_bindings']
    if (not binding['passed'] or not all(binding['current_defining_object_graph'].values())
            or not all(binding['exact_terminal_function_and_unit_identities'].values())
            or not binding['current_future_callable_used_for_every_Z']
            or not binding['no_saved_terminal_Xv_or_energy_lookup']
            or not binding['actual_Rv_five_history_adapter_available']
            or not binding['incoming_power_data_are_Rp_reference_not_terminal_histories']
            or binding['source_caps_used_as_defining_field_values']
            or binding['interval_overlap_used_as_join_proof']):
        raise ValueError('Current functional endpoint binding omitted')
    if (set(raw['ordered_current_chart_registry'])!=set(CHARTS)|{'flatten'}
            or raw['current_downstream_chart_owner_count']!=21
            or raw[GATE] or not raw[PROVED] or raw[UNIFORM]
            or raw['current_flatten_physical_owner_installed']
            or raw['full_pulse_C4_installed'] or any(raw[k] for k in SCOPES+OPEN)):
        raise ValueError('Current flatten source scope differs')
    c=field.flatten.ctx
    def interval(row):
        value=read_interval(c,row);lo,hi=endpoints(value)
        if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:
            raise ValueError('Nonfinite or unordered flatten source bounds')
        return value
    def exact_zero(row):
        if endpoints(interval(row))!=(mp.mpf(0),mp.mpf(0)):
            raise ValueError('Empty-support terminal moment was not exactly zero')
    endpoint=raw['current_native_terminal_packet']['source_packet']
    for flag in ('terminal_linear_zero_by_empty_future_supports','positive_energy_target_at_Rv'):
        if not endpoint[flag]:raise ValueError('Actual native endpoint source not used')
    for key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
        rows=endpoint[key]['coefficients']
        if len(rows)!=6:raise ValueError('Terminal axial5 history omitted')
        for row in rows:exact_zero(row)
    if endpoints(interval(endpoint['Mztheta_over_R_Utheta_squared']['coefficients'][0]))[0]<=0:
        raise ValueError('Positive complete future terminal energy lost')
    if endpoints(interval(endpoint['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][0]))[0]<=0:
        raise ValueError('Nonzero terminal angular memory lost')
    for key in ('Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared'):
        rows=endpoint['pressure'][key]['coefficients']
        if len(rows)!=6:raise ValueError('Absolute pressure axial5 history omitted')
        for row in rows:interval(row)
    indices={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
    packets=raw['current_flatten_source_packets'];counts={}
    if set(packets)!={'inlet','whole_domain','exit'}:raise ValueError('Original flatten coverage omitted')
    for name,packet in packets.items():
        if (packet[GATE] or not packet[PROVED] or packet[UNIFORM]
                or not packet['current_native_terminal_provider_used']
                or not packet['current_future_energy_callable_used']
                or packet['datum_enclosure_sha256']!=field.datum_sha
                or any(packet[k] for k in SCOPES+OPEN)):
            raise ValueError('Flatten packet source or scope differs')
        source=packet['source_packet'];grid=source['physical_mixed_derivatives_total_order_le4']
        histories=source['current_terminal_histories']
        if (not source['incoming_power_data_are_Rp_reference_not_terminal_histories']
                or not source['current_terminal_zero_linear_histories_from_native_supports']):
            raise ValueError('Rp reference data was relabeled as Rv terminal histories')
        for key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
            if len(histories[key]['coefficients'])!=6:raise ValueError('Actual terminal linear axial5 rows omitted')
            for row in histories[key]['coefficients']:exact_zero(row)
        for key in ('Mtheta_over_sqrt2_R_3half_Utheta','Mztheta_over_R_Utheta_squared',
                'Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared'):
            if len(histories[key]['coefficients'])!=6:raise ValueError('Actual terminal nonzero history omitted')
            for row in histories[key]['coefficients']:interval(row)
        if set(grid)!={UZ,UT,UR,P}:raise ValueError('Flatten velocity/pressure component omitted')
        count=0
        for label,rows in grid.items():
            if set(rows)!=indices:raise ValueError('Flatten ordinary mixed4 row omitted')
            for row in rows.values():
                interval(row)
                if label in (UZ,UR):exact_zero(row)
                count+=1
        if not source['original_pressure_datum_and_complete_positive_future_energy_preserved']:
            raise ValueError('Common energy/datum history lost')
        if endpoints(interval(source['energy_Taylor']['coefficients'][0]))[0]<=0:
            raise ValueError('Whole flatten positive future energy lost')
        counts[name]=count
    # A new Z request exercises the current future callable, not a saved sample.
    fresh=field.evaluate('flatten','0.271',0)['source_packet']
    terminal=field.pulse.end('0.271',0)
    for key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
        if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in terminal[key].coefficients):
            raise ValueError('Fresh current endpoint empty supports differ')
    if endpoints(fresh['energy_Taylor'][0])[0]<=0:
        raise ValueError('Fresh current future integral lost positivity')
    for Z,t in ((['-1.001',1],0),(0,-1),(0,101)):
        try:field.evaluate('flatten',Z,t)
        except ValueError:pass
        else:raise ValueError('Original flatten domain widened')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_downstream_chart_owners_checked=21,
        current_flatten_mixed4_rows_checked=counts,total_new_flatten_source_rows_checked=sum(counts.values()),
        retained_canonical_endpoint_theorem=THEOREM,retained_canonical_identity_counts=binding['retained_canonical_identity_counts'],
        fresh_unsaved_Z_current_future_callable_checked=True,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,
        **{GATE:True,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
        current_flatten_physical_owner_installed=False,
        **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current native terminal -> flatten: twenty-one profile owners, retained exact histories and 180 mixed4 source rows',flush=True)
    return result


if __name__=='__main__':run()
