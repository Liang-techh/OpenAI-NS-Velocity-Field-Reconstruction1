"""Focused current three-bridge Cartesian/time admission; retain current30."""
import ast
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_bridge_physical_assembly import (
    CurrentBridgePhysicalAssembly,BRIDGES,CHARTS,SOURCE_VIEWS,GATE,COMPOSED,COVERAGE,
    UNIFORM,SCOPES,OPEN,HERE,PREFIX,sha)
from lei_ren_part1_paper_compliant_global_physical_assembly import INDICES,COMPONENTS,UR,UT,UZ,P
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_bridge_physical_assembly.json'
PARTS={'ux':{UR,UT},'uy':{UR,UT},'uz':{UZ},'p':{P}}


def actual_bridge_call_bindings():
    tree=ast.parse((HERE/(PREFIX+'current_bridge_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentBridgeSourceOverlay')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    expression='self.bridge.evaluate(Z,coordinate,BRIDGES[chart])'
    wanted=ast.dump(ast.parse(expression,mode='eval').body)
    if sum(isinstance(n,ast.Call) and ast.dump(n)==wanted for n in ast.walk(fn))!=1:
        raise ValueError('Original coordinate-labelled actual bridge call missing')
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentBridgePhysicalAssembly')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    expression='BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)'
    wanted=ast.dump(ast.parse(expression,mode='eval').body)
    if sum(isinstance(n,ast.Call) and ast.dump(n)==wanted for n in ast.walk(fn))!=1:
        raise ValueError('Original full physical mapper call missing')
    return dict(same_current_bridge_provider_used_for_packet_and_radius=True,
        actual_coordinate_labelled_bridge_call_AST_bound=True,
        actual_original_full_physical_evaluate_AST_bound=True,passed=True)


def run(base=None):
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentBridgePhysicalAssembly(require_checked=False,base=base)
        if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='whole_current_bridge_physical_maps'}:
            raise ValueError('Fresh current33 physical manifest differs')
        if (tuple(raw['current_physical_chart_owners'])!=CHARTS or len(CHARTS)!=33
                or not all(raw['current_provider_graph_identity'].values())
                or not raw['original_bridge_mapper_source_bindings']['passed']
                or raw[GATE] or raw[COMPOSED] or raw[COVERAGE] or raw[UNIFORM]
                or raw['full_pulse_C4_installed'] or raw['core_inner_annulus_interfaces_certified']
                or any(raw[k] for k in SCOPES+OPEN)):
            raise ValueError('Current33 source coverage/scope differs')
        if not all(v['passed'] for v in raw['reused_independent_unchanged_physical_fixtures'].values()):
            raise ValueError('Unchanged physical Cartesian/time/width fixtures omitted')
        c=MPIntervalContext();c.dps=240
        def finite(row):
            value=read_interval(c,row)
            if not all(mp.isfinite(v) for v in endpoints(value)):raise ValueError('Nonfinite physical bridge row')
            return value
        counts={};term_count=0;zero_count=0
        def row_check(row):
            nonlocal term_count,zero_count
            if (row['source_row_mode']!='uncapped_factored_rows'
                    or not row['original_source_factors_combined_before_enclosure']
                    or not row['global_source_factors_combined_before_final_physical_bound']
                    or not row['positive_source_exponentials_not_materialized']
                    or endpoints(finite(row['physical_lambda_exponent']))[1]>=0
                    or bool(row['exact_zero'])!=(not row['terms'])
                    or bool(row['exact_zero'])!=(row['log_absolute_upper'] is None)):
                raise ValueError('Original factored physical row/lambda semantics differ')
            if row['exact_zero']:zero_count+=1
            else:finite(row['log_absolute_upper'])
            for term in row['terms']:
                if len(term['source_log_exponents'])!=7:raise ValueError('Physical source basis omitted')
                for power in term['source_log_exponents']:
                    value=mp.make_mpf(tuple(power['exact_mpf_tuple']))
                    if not mp.isfinite(value):raise ValueError('Nonfinite formal source power')
                finite(term['signed_coefficient']);finite(term['log_absolute_upper']);term_count+=1
        indices={'x'+str(i)+'_y'+str(j)+'_z'+str(k) for i,j,k in INDICES}
        views=raw['whole_current_bridge_physical_maps']
        if set(views)!=set(SOURCE_VIEWS):raise ValueError('Original bridge domain/endpoints omitted')
        for name,packet in views.items():
            phase,coordinate=SOURCE_VIEWS[name];chart='bridge_'+phase
            if (packet['chart']!=chart or packet[GATE] or packet[COMPOSED] or packet[COVERAGE]
                    or not packet['current_actual_bridge_source_used']
                    or packet['current_source_owner']!=field.source_owners[chart]['provider']
                    or packet['datum_enclosure_sha256']!=field.datum_sha
                    or packet[UNIFORM] or packet['full_pulse_C4_installed']
                    or any(packet[k] for k in SCOPES+OPEN)):
                raise ValueError('Current physical bridge owner/scope differs')
            if (packet['global_source_row_mode']!='uncapped_factored_rows'
                    or not packet['normalization_not_differentiated_twice']
                    or not packet['moving_cylindrical_basis_differentiated']
                    or bool(packet['microscope_phase_to_logR_applied_before_source_bound'])!=(phase!='macro')
                    or endpoints(finite(packet['Z']))!=(-1,1)):
                raise ValueError('Original bridge Cartesian basis/coordinate/domain omitted')
            source=packet['original_radius_source']
            if (not source['microscopic_radius_variation_retained_formally']
                    or not source['absolute_radius_not_rounded_as_source']
                    or not source['numeric_logR_bound_uses_positive_hb_enclosure']
                    or source['formal_hb_radius_correlation_evaluated']):
                raise ValueError('Radius enclosure promoted to correlated exact point')
            if name=='R100_exit':
                lower,upper=endpoints(finite(packet['source_logR_enclosure']))
                if not lower<=mp.log(100)<=upper:raise ValueError('Original exact R100 radius lost')
            spatial=packet['physical_spatial_cartesian_mixed4'];time=packet['first_fixed_x_physical_time_derivative']
            if set(spatial)!=indices or set(time)!=set(COMPONENTS):raise ValueError('Spatial4/time1 indices omitted')
            count=0
            for components in list(spatial.values())+[time]:
                if set(components)!=set(COMPONENTS):raise ValueError('Cartesian component omitted')
                for component,parts in components.items():
                    if set(parts)!=PARTS[component]:raise ValueError('Physical velocity/pressure source contribution omitted')
                    for row in parts.values():row_check(row);count+=1
            if count!=216:raise ValueError('Bridge210 spatial/6 time contributions omitted')
            counts[name]=count
        for chart,coordinate in (('bridge_first','.4'),('bridge_second','1.6'),('bridge_macro','.37')):
            fresh=field.evaluate(chart,'0.271',coordinate,log_tau='-10',theta='.3')
            if endpoints(fresh['requested_log_tau'])!=(-10,-10) or not fresh['current_actual_bridge_source_used']:
                raise ValueError('Fresh current bridge/time request acquisition failed')
        for chart,coordinate,axis in (('bridge_first',-1,False),('bridge_second',0,False),('bridge_macro','1.1',False),
                ('bridge_first',0,True),('axis',0,False),('core',1,True)):
            try:field.evaluate(chart,0,coordinate,axis=axis)
            except ValueError:pass
            else:raise ValueError('Original current33 domain or nonsingular axis widened')
        old=raw['retained_current30_physical_evidence']
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_physical_chart_owners_checked=33,
            current_bridge_spatial_and_time_contributions_checked=counts,
            total_new_bridge_physical_contributions_checked=sum(counts.values()),
            thirty_three_regular_source_contributions_checked=old['regular_contributions']+3*216,
            retained_current30_physical_evidence=old,current_bridge_signed_source_terms_checked=term_count,
            current_bridge_exact_zero_source_contributions_checked=zero_count,
            actual_bridge_call_bindings=actual_bridge_call_bindings(),
            reused_independent_unchanged_physical_fixtures=raw['reused_independent_unchanged_physical_fixtures'],
            fresh_unsaved_Z_all_three_charts_and_requested_time_checked=True,
            **{GATE:True,COMPOSED:True,COVERAGE:True},current_thirty_source_physical_ownership_certified=True,
            core_inner_annulus_interfaces_certified=False,full_pulse_C4_installed=False,**{UNIFORM:False},
            **dict.fromkeys(SCOPES,False),**{k:False for k in OPEN if k!=COVERAGE},all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current33 physical source maps:1944 new bridge contributions; points/global/interfaces/recursion remain open',flush=True)
    return result


if __name__=='__main__':run()
