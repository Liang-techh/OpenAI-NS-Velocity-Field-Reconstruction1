"""Focused admission of live current collar/full Gamma sources."""
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_heat_source import (
    CurrentHeatSourceAssembly, CHARTS, NEW_CHARTS, VIEWS, GATE, PROVED,
    THEOREM, GAMMA_THEOREM, UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_global_physical_assembly import UZ,UT,UR,P
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_heat_source.json'


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentHeatSourceAssembly(require_checked=False)
        manifest={k:v for k,v in raw.items() if k not in ('whole_current_heat_source_maps','current_waiting_terminal')}
        if encode(pack(field.manifest()))!=manifest:raise ValueError('Fresh current heat source manifest differs')
        binding=raw['current_heat_source_bindings'];shared=binding['shared_exact_Gamma_future_binding']
        if (not all(binding['current_defining_object_graph'].values()) or not binding['passed']
                or not shared['passed'] or not all(shared['exact_full_Gamma_future_identities'].values())
                or binding['source_caps_used_as_defining_field_values'] or binding['interval_overlap_used_as_join_proof']
                or raw[GATE] or not raw[PROVED] or raw[UNIFORM] or raw['full_pulse_C4_installed']
                or any(raw[k] for k in SCOPES+OPEN)
                or tuple(raw['ordered_current_chart_registry'])!=CHARTS):
            raise ValueError('Current heat defining source or scope differs')
        c=MPIntervalContext();c.dps=240
        def interval(row):
            value=read_interval(c,row);lo,hi=endpoints(value)
            if not lo<=hi or any(not mp.isfinite(v) for v in (lo,hi)):raise ValueError('Nonfinite source derivative bound')
            return value
        def exact_zero(row):
            if endpoints(interval(row))!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Inherited exact zero source lost')
        groups=raw['whole_current_heat_source_maps'];indices={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
        if set(groups)!=set(NEW_CHARTS):raise ValueError('Both current heat source owners required')
        counts={};zeros=0
        for chart,views in groups.items():
            if set(views)!=set(VIEWS[chart]):raise ValueError('Full collar/infinite exterior coverage omitted')
            for name,wrapped in views.items():
                source=wrapped['source_packet']
                if (wrapped['source_provider']!=field.registry[chart]['provider'] or wrapped[GATE]
                        or not wrapped[PROVED] or wrapped['datum_enclosure_sha256']!=field.datum_sha
                        or any(wrapped[k] for k in SCOPES+OPEN)):
                    raise ValueError('Current heat packet owner/scope differs')
                for flag in ('current_live_waiting_heat_source_used','current_nonzero_angular_and_absolute_pressure_histories_retained',
                        'original_flat_collar_and_exact_Gamma_source_retained','same_waiting_terminal_angular_pressure_histories_retained',
                        'full_infinite_Gamma_tail_included','actual_Gamma_not_defined_by_finite_S_series'):
                    if not source[flag]:raise ValueError('Live current full heat source omitted: '+flag)
                if source['full_pulse_C4_installed'] or any(source[k] for k in ('full_outer_C4_certified',
                        'physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion')):
                    raise ValueError('Current heat source scope promoted')
                target=VIEWS[chart][name];target=target if isinstance(target,list) else [target,target]
                lo,hi=endpoints(read_interval(c,source['coordinate']['offset']))
                if (lo,hi)!=endpoints(field.heat.ctx.mpf(target)):raise ValueError('Original heat coordinate coverage differs')
                if chart=='heat_exterior' and not source['unbounded_exterior_covered']:
                    raise ValueError('Full exact exterior source omitted')
                for key in ('theta_over_Ev0_Taylor','angular_Taylor','energy_Taylor','pressure_over_Pstar_squared_Taylor'):
                    rows=source[key]['coefficients']
                    if len(rows)!=6:raise ValueError('Canonical whole-Z axial5 source omitted')
                    for row in rows:interval(row)
                for jet in source['heat_bracket_y_derivatives']:
                    for row in jet['coefficients']:interval(row)
                if (endpoints(interval(source['energy_Taylor']['coefficients'][0]))[0]<=0
                        or endpoints(interval(source['heat_bracket_y_derivatives'][0]['coefficients'][0]))[0]<=0):
                    raise ValueError('Full future energy / heat bracket positivity lost')
                grid=source['physical_mixed_derivatives_total_order_le4'];count=0
                if set(grid)!={UZ,UT,UR,P}:raise ValueError('Velocity/pressure source component omitted')
                for label,rows in grid.items():
                    if set(rows)!=indices:raise ValueError('Heat mixed4 source index omitted')
                    for row in rows.values():
                        interval(row);count+=1
                        if label in (UZ,UR):exact_zero(row);zeros+=1
                counts[chart+'_'+name]=count
        inlet=groups['heat_collar']['inlet']['source_packet'];terminal=raw['current_waiting_terminal']['source_packet']
        if inlet['pressure_over_Pstar_squared_Taylor']!=terminal['pressure_over_Pstar_squared_Taylor']:
            raise ValueError('Actual live waiting absolute pressure not inherited exactly')
        for jet in inlet['heat_bracket_y_derivatives'][1:]:
            for row in jet['coefficients']:exact_zero(row)
        for view in ('exit','phi_crossing'):
            source=groups['heat_collar'][view]['source_packet']
            if view=='exit':
                for row in source['original_collar_phi_y_Taylor']['coefficients']:exact_zero(row)
        # One new unsaved live Z acquisition. Unchanged Gamma/flat/rate
        # fixtures are consumed by source hash, not repeatedly rerun.
        fresh=field.heat.data('0.271')
        if fresh['waiting_terminal']['pressure_over_Pstar_squared_Taylor'].coefficients!=fresh['Ptail'].coefficients:
            raise ValueError('Fresh unsaved absolute pressure history differs')
        for chart,value in (('heat_collar','1.2'),('heat_exterior','4.2')):
            packet=field.evaluate(chart,'0.271',value)['source_packet']
            if endpoints(packet['energy_Taylor'][0])[0]<=0:raise ValueError('Fresh positive heat source lost')
        for chart,value in (('heat_collar',-1),('heat_collar','3.1'),('heat_exterior','2.9')):
            try:field.evaluate(chart,0,value)
            except ValueError:pass
            else:raise ValueError('Original heat chart domain widened')
        try:field.evaluate('heat_collar',[-1.001,1],0)
        except ValueError:pass
        else:raise ValueError('Original axial domain widened')
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_downstream_chart_owners_checked=29,
            newly_checked_current_source_charts=list(NEW_CHARTS),current_heat_mixed4_rows_checked=counts,
            total_new_heat_source_rows_checked=sum(counts.values()),exact_zero_meridional_rows_checked=zeros,
            current_shared_exact_Gamma_future_identities_checked=len(shared['exact_full_Gamma_future_identities']),
            retained_canonical_functional_identity_count=len(field.functional_proof),retained_canonical_endpoint_theorem=THEOREM,
            retained_full_Gamma_enclosure_theorem=GAMMA_THEOREM,
            entire_unbounded_exterior_and_flat_phi_crossing_checked=True,
            fresh_unsaved_Z_current_source_acquisition_checked=True,
            source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,
            **{GATE:True,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
            current_heat_physical_owner_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current29 source owners: waiting/collar/full Gamma, 360 new mixed4 rows',flush=True)
    return result


if __name__=='__main__':run()
