"""Focused current core/axis graph and spatial4/time1 source admission."""
import ast
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_core_physical_assembly import (
    CurrentCorePhysicalAssembly,CHARTS,VIEWS,GATE,COMPOSED,SOURCE_GATE,
    THEOREM,UNIFORM,SCOPES,OPEN,HERE,PREFIX,sha)
from lei_ren_part1_paper_compliant_core_physical_field import Q,F,V,PD,PI,gridkey,BASES
from lei_ren_part1_paper_compliant_global_physical_assembly import INDICES,COMPONENTS
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_core_physical_assembly.json'


def actual_current_core_call_bindings():
    tree=ast.parse((HERE/(PREFIX+'current_core_physical_assembly.py')).read_text(encoding='utf8'))
    expected=('self.core.profiles(rho,z)',
        'BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)')
    for expression in expected:
        target=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(isinstance(n,ast.Call) and ast.dump(n)==target for n in ast.walk(tree))!=1:
            raise ValueError('Actual current core/physical call differs: '+expression)
    return dict(actual_nested_current_core_profiles_call_bound=True,
        actual_original_physical_evaluate_with_axis_mode_call_bound=True,passed=True)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentCorePhysicalAssembly(require_checked=False)
        if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='whole_current_core_physical_maps'}:
            raise ValueError('Fresh current30 core manifest differs')
        if (tuple(raw['current_physical_chart_owners'])!=CHARTS or not all(raw['current_provider_graph_identity'].values())
                or not raw['original_core_operator_bindings']['passed'] or raw[GATE] or raw[COMPOSED] or raw[SOURCE_GATE]
                or raw['current_core_axis_physical_owner_installed'] or not raw['current_twenty_nine_downstream_physical_source_ownership_certified']
                or raw['core_inner_annulus_interfaces_certified'] or raw[UNIFORM] or raw['full_pulse_C4_installed']
                or any(raw[k] for k in SCOPES+OPEN)):
            raise ValueError('Current core graph/admission scope differs')
        c=MPIntervalContext();c.dps=240
        def interval(row):
            value=read_interval(c,row);lo,hi=endpoints(value)
            if not lo<=hi or not all(mp.isfinite(v) for v in (lo,hi)):raise ValueError('Nonfinite current core bound')
            return value
        def exact_zero(row):
            if endpoints(interval(row))!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Exact analytic axis zero lost')
        views=raw['whole_current_core_physical_maps'];counts={};sourcecounts={};zero_count=0
        indices={'x'+str(i)+'_y'+str(j)+'_z'+str(k) for i,j,k in INDICES}
        profileindices={gridkey(i,k) for i in range(5) for k in range(5-i)}
        if set(views)!=set(VIEWS):raise ValueError('Current core/axis/exit coverage omitted')
        for name,wrapped in views.items():
            source=wrapped['source']['source_packet'];physical=wrapped['physical'];rho,axis=VIEWS[name]
            rho=rho if isinstance(rho,list) else [rho,rho]
            if (endpoints(interval(source['rho']))!=tuple(mp.mpf(v) for v in rho)
                    or endpoints(interval(source['Z']))!=(-1,1)
                    or not source['actual_selected_Cstar_source_retained']
                    or not source['nonlinear_field_enclosed_not_replaced_by_model']
                    or not source['original_P0_and_centrifugal_pressure_increment_retained']):
                raise ValueError('Current analytic core source or domain differs')
            grids=source['ordinary_mixed_profile_grids']
            if set(grids)!={Q,F,V,PD,PI}:raise ValueError('Actual core profile/pressure normalization omitted')
            sourcecounts[name]=0
            for grid in grids.values():
                if set(grid)!=profileindices:raise ValueError('Core mixed4 source index omitted')
                for row in grid.values():interval(row);sourcecounts[name]+=1
            if endpoints(interval(source['Phi_derivatives'][gridkey(0,0)]))[0]<=0:
                raise ValueError('Current positive actual normalized core swirl lost')
            if endpoints(interval(source['rho']))==(0,0):
                if endpoints(interval(source['Phi_derivatives'][gridkey(0,0)]))!=(1,1):raise ValueError('Axis Phi=1 lost')
                for k in range(1,6):exact_zero(source['Phi_derivatives'][gridkey(0,k)])
                for k in range(5):exact_zero(grids[PI][gridkey(0,k)])
            if (physical['datum_enclosure_sha256']!=field.datum_sha or physical[GATE] or physical[COMPOSED]
                    or physical[SOURCE_GATE] or not physical['current_core_axis_source_functional_binding_proved']
                    or physical['current_source_owner']!=field.source_owners['core']['provider']
                    or any(physical[k] for k in SCOPES+OPEN)):
                raise ValueError('Current physical core owner/scope differs')
            native=physical['core_axis_nonsingular_native_map']
            if not native['includes_axis_without_inverse_radius'] or bool(native['axis_only'])!=axis:
                raise ValueError('Nonsingular axis mode omitted or widened')
            spatial=native['cartesian_spatial_multiindices'];time=native['first_fixed_x_physical_time_derivative']
            if set(spatial)!=indices or set(time)!=set(COMPONENTS):raise ValueError('Core spatial4/time1 coverage omitted')
            count=0
            for index,components in spatial.items():
                if set(components)!=set(COMPONENTS):raise ValueError('Core Cartesian component omitted')
                for component,parts in components.items():
                    if set(parts)!=set(BASES[component]):raise ValueError('Core physical source contribution omitted')
                    for label,row in parts.items():
                        interval(row['bracket']);norm=interval(row['absolute_upper'])
                        log=physical['requested_core_spatial_log_bounds'][index][component][label]
                        if bool(row['exactly_zero'])!=(log is None) or endpoints(norm)[0]<0:
                            raise ValueError('Core zero/prefactor bound differs')
                        scale=native['shared_physical_prefactor_bounds'][row['scale_key']]
                        if endpoints(interval(scale['physical_lambda_exponent']))[1]>=0:raise ValueError('Core lambda bound sign lost')
                        if log is None:zero_count+=1
                        else:interval(log)
                        count+=1
            for component,parts in time.items():
                if set(parts)!=set(BASES[component]):raise ValueError('Core fixed-position time contribution omitted')
                for label,row in parts.items():
                    interval(row['bracket']);interval(row['absolute_upper'])
                    log=physical['requested_core_time_log_bounds'][component][label]
                    if bool(row['exactly_zero'])!=(log is None):raise ValueError('Core time zero/prefactor bound differs')
                    if endpoints(interval(row['physical_lambda_exponent']))[1]>=0:raise ValueError('Core time lambda sign lost')
                    if log is None:zero_count+=1
                    else:interval(log)
                    count+=1
            if count!=252:raise ValueError('Core 245 spatial/7 time contributions omitted')
            counts[name]=count
        fresh=field.evaluate('core','0.271','2.3',log_tau='-10')
        if endpoints(fresh['requested_log_tau'])!=(mp.mpf(-10),mp.mpf(-10)):
            raise ValueError('Fresh requested core time lost')
        fresh_axis=field.evaluate('core','0.271',0,axis=True,log_tau='-10')
        if not fresh_axis['core_axis_nonsingular_native_map']['axis_only']:raise ValueError('Fresh current axis acquisition failed')
        for chart,rho,axis in (('core',-1,False),('core','4.01',False),('core',1,True),
                ('axis',0,False),('bridge_first',0,False),('heat_collar',0,True)):
            try:field.evaluate(chart,0,rho,axis=axis)
            except ValueError:pass
            else:raise ValueError('Original current core/axis domain widened')
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_physical_chart_owners_checked=30,
            current_core_source_mixed4_rows_checked=sourcecounts,total_new_core_source_rows_checked=sum(sourcecounts.values()),
            current_core_spatial_and_time_contributions_checked=counts,total_new_core_physical_contributions_checked=sum(counts.values()),
            thirty_regular_source_contributions_checked=raw['retained_twenty_nine_physical_evidence']['regular_contributions']+counts['whole_core'],
            supplemental_axis_contributions_checked=counts['whole_axis'],exact_zero_core_contributions_checked=zero_count,
            actual_current_core_call_bindings=actual_current_core_call_bindings(),
            retained_canonical_core_theorem=THEOREM,reused_independent_nonsingular_core_fixture=raw['reused_independent_nonsingular_core_fixture'],
            retained_twenty_nine_physical_evidence=raw['retained_twenty_nine_physical_evidence'],
            fresh_unsaved_core_axis_and_requested_time_checked=True,
            **{GATE:True,COMPOSED:True,SOURCE_GATE:True},current_core_axis_physical_owner_installed=True,
            core_inner_annulus_interfaces_certified=False,full_pulse_C4_installed=False,**{UNIFORM:False},
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current30 core/axis: 300 new source rows, 1008 physical contributions; prior29 retained',flush=True)
    return result


if __name__=='__main__':run()
