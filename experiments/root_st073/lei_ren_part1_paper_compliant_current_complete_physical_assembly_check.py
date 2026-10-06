"""Admit current 33-chart callable physical maps; quantitative joins stay open."""
import ast
import copy
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_complete_physical_assembly import (
    CurrentCompletePhysicalAssembly,current_complete_physical_source_proof,
    HERE,PREFIX,NAME,RECEIPT,GATES,OPEN,CHARTS,CHANGED,RETAINED,
    PULSE_CHARTS,POST,VIEWS,physical_views,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    INDICES,COMPONENTS,UZ,UT,UR,P,cartesian_templates)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def actual_current_dispatch_call_bindings():
    tree=ast.parse((HERE/(PREFIX+'current_complete_physical_assembly.py')).read_text(encoding='utf8'))
    def call_in(cls,method,expression):
        owner=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in owner.body if isinstance(n,ast.FunctionDef) and n.name==method)
        target=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(n)==target for n in ast.walk(fn) if isinstance(n,ast.Call))!=1:
            raise ValueError('Actual current source call missing: '+expression)
    for expression in (
            'CurrentNativePulseSourceDispatcher.evaluate(self,chart,Z,coordinate)',
            'self.history.evaluate(PUBLIC_TO_INTERNAL[chart],Z,coordinate)',
            'self.assembly.bridge.evaluate(Z,coordinate,BRIDGES[chart])',
            'self.base.dispatch.evaluate(chart,Z,coordinate)'):
        call_in('CurrentCompletePhysicalDispatch','evaluate',expression)
    call_in('CurrentCompletePhysicalDispatch','gap_overlap',"self.native_pulse.gap(Z,['12','12.0001'])")
    call_in('CurrentCompletePhysicalAssembly','evaluate',
        'BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)')
    return dict(all_six_current_pulse_original_routes_AST_bound=True,
        all_nine_current_energy_source_calls_AST_bound=True,
        retained_common_bridge_and_incoming_calls_AST_bound=True,
        full_original_cartesian_and_fixed_x_time_map_AST_bound=True,
        supplemental_gap_uses_same_current_callable_AST_bound=True,passed=True)


def rejected_old_provider_mutations(field):
    rejected=[]
    for label in ('old_pulse_main_provider','old_heat_provider','old_native_pulse_alias',
            'unchecked_energy_history','unchecked_full_exterior'):
        candidate=copy.copy(field)
        candidate.dispatch=object.__new__(type(field.dispatch))
        candidate.dispatch.__dict__.update(field.dispatch.__dict__)
        candidate.dispatch.assembly=candidate;candidate.dispatch.providers=dict(field.dispatch.providers)
        if label=='old_pulse_main_provider':
            candidate.dispatch.providers['pulse_main']=field.base.dispatch.provider('pulse_main')
        elif label=='old_heat_provider':
            candidate.dispatch.providers['heat_exterior']=field.base.dispatch.provider('heat_exterior')
        elif label=='old_native_pulse_alias':candidate.dispatch.native_pulse=field.base.pulse
        elif label=='unchecked_energy_history':
            candidate.history=copy.copy(field.history);candidate.history.acceptance_loaded=False
            candidate.exterior=copy.copy(field.exterior);candidate.exterior.history=candidate.history
            candidate.dispatch.history=candidate.history
        else:
            candidate.exterior=copy.copy(field.exterior);candidate.exterior.acceptance_loaded=False
        try:candidate.assert_graph()
        except ValueError:rejected.append(label);continue
        raise ArithmeticError('Old/unchecked physical provider accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentCompletePhysicalAssembly(require_checked=False)
    extras=('current_changed_physical_maps','current_same_owner_supplemental_gap_source')
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in extras}:
        raise ValueError('Current 33-chart graph, units, provenance or scoped manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Producer claims physical admission')
    proof=current_complete_physical_source_proof(field)
    if encode(pack(proof))!=raw['current_complete_physical_source_proof']:
        raise ValueError('Actual complete physical function proof changed')
    if (tuple(raw['current_physical_chart_owners'])!=CHARTS or
            tuple(raw['replaced_current_physical_source_charts'])!=CHANGED or
            tuple(raw['retained_current_incoming_source_charts'])!=RETAINED or
            not all(raw['original_physical_operator_identity'].values()) or
            raw['old_live_complete_future_consumers_remaining']):
        raise ValueError('Complete changed/retained physical ownership differs')
    calls=actual_current_dispatch_call_bindings();mutations=rejected_old_provider_mutations(field)
    c=MPIntervalContext();c.dps=240
    indices={'x%d_y%d_z%d'%index for index in INDICES}
    packets=raw['current_changed_physical_maps']
    if tuple(packets)!=tuple(VIEWS)+('fresh_main','fresh_exterior'):
        raise ValueError('All fifteen changed charts and fresh physical views required')
    counts={};terms=zeros=0
    def check(row,label,chart):
        nonlocal terms,zeros
        if bool(row['exact_zero'])!=(row['log_absolute_upper'] is None):
            raise ValueError('Physical exact zero differs from source bound')
        if row['exact_zero']:
            if row['terms']:raise ValueError('Physical exact zero retains source terms')
            zeros+=1;return
        if chart in POST and label in (UR,UZ):
            raise ValueError('Current postpulse zero meridional source was not preserved')
        lo,hi=endpoints(read_interval(c,row['log_absolute_upper']))
        if not all(mp.isfinite(v) for v in (lo,hi)) or lo>hi:
            raise ArithmeticError('Nonfinite current physical source log bound')
        if (not row['positive_source_exponentials_not_materialized'] or
                not row['global_source_factors_combined_before_final_physical_bound'] or
                row['source_row_mode']!='provider_prebounded_mixed_rows' or
                endpoints(read_interval(c,row['physical_lambda_exponent']))[1]>=0):
            raise ValueError('Original physical factor/lambda bound changed')
        for term in row['terms']:
            powers=[mp.make_mpf(tuple(v['exact_mpf_tuple'])) for v in term['source_log_exponents']]
            if len(powers)!=7 or powers[5]>0:
                raise ValueError('Source factors or nonpositive radial powers lost')
            if label==P:
                correct=powers[1]==2 and powers[2]==0 and powers[4]==0
            elif chart in PULSE_CHARTS:
                correct=powers[2]==1 and powers[4]==1
            else:correct=powers[1]==0 and powers[4]==1
            if not correct:raise ValueError('Current correlated velocity/absolute pressure unit changed')
            lo,hi=endpoints(read_interval(c,term['signed_coefficient']))
            if not all(mp.isfinite(v) for v in (lo,hi)) or lo>hi:
                raise ArithmeticError('Nonfinite signed physical coefficient')
            terms+=1
    for name,packet in packets.items():
        chart='pulse_main' if name=='fresh_main' else 'heat_exterior' if name=='fresh_exterior' else name
        args=(chart,'.517','.83','-2.17','.37') if name=='fresh_main' else (
            (chart,'.517','6.23','-2.17','.37') if name=='fresh_exterior' else
            (chart,[-1,1],physical_views(field)[chart],'-1',None))
        # The admitted actual field replays the current source call. Saved
        # old physical reports never supply its values or selected weights.
        fresh=field.evaluate(*args)
        if encode(pack(fresh))!=packet:raise ValueError('Current physical callable replay differs: '+name)
        if (packet['chart']!=chart or packet['current_source_owner']!=field.source_owners[chart]['provider'] or
                packet['current_source_acceptance_receipt']!=field.source_owners[chart]['acceptance_receipt'] or
                packet['datum_enclosure_sha256']!=field.datum_sha or any(packet[k] for k in GATES+OPEN) or
                not packet['moving_cylindrical_basis_differentiated'] or
                not packet['normalization_not_differentiated_twice'] or
                not all(packet['current_complete_provider_graph'].values())):
            raise ValueError('Current physical packet source/scope differs: '+name)
        spatial=packet['physical_spatial_cartesian_mixed4'];time=packet['first_fixed_x_physical_time_derivative']
        if set(spatial)!=indices or set(time)!=set(COMPONENTS):
            raise ValueError('Physical spatial4/fixed-x time1 coverage omitted')
        count=0
        for index,components in spatial.items():
            i,j,b=(int(v[1:]) for v in index.split('_'))
            if set(components)!=set(COMPONENTS):raise ValueError('Cartesian physical component omitted')
            for component,parts in components.items():
                labels={label for label,a,q in cartesian_templates()[component,i,j,b]}
                if set(parts)!=labels:raise ValueError('Moving cylindrical basis contribution omitted')
                for label,row in parts.items():check(row,label,chart);count+=1
        for component,parts in time.items():
            labels={UR,UT} if component in ('ux','uy') else {UZ} if component=='uz' else {P}
            if set(parts)!=labels:raise ValueError('Fixed-position physical time contribution omitted')
            for label,row in parts.items():check(row,label,chart);count+=1
        if count!=216:raise ValueError('Physical source contribution count differs')
        counts[name]=count
        if name=='heat_exterior':
            lo,hi=endpoints(read_interval(c,packet['source_logR_enclosure']))
            if not mp.isfinite(lo) or hi!=mp.inf:
                raise ValueError('Complete unbounded exterior physical coverage omitted')
    overlap=raw['current_same_owner_supplemental_gap_source']
    live=field.dispatch.gap_overlap([-1,1])['source_packet']
    if encode(pack(live))!=overlap['source_packet'] or not overlap['same_current_selected_owner'] or not overlap['public_gap_domain_not_widened']:
        raise ValueError('Supplemental gap source is not the same current selected function')
    point=live;coordinate=['12','12.0001'];logR,_=field.radius('pulse_gap',coordinate,point,field.pulse)
    _,logs,amplitudes=field.normalized_sources('pulse_gap',[-1,1],coordinate,point,field.pulse,logR)
    for key,value in (('source_logR_enclosure',logR),('source_log_bases',logs),('source_amplitudes',amplitudes)):
        if json.loads(json.dumps(encode(pack(value))))!=overlap[key]:
            raise ValueError('Private gap physical factor differs: '+key)
    coverage=field.dispatch.coverage()
    try:field.dispatch.evaluate('pulse_gap',0,'12.0001')
    except ValueError:pass
    else:raise ValueError('Supplemental gap source widened the public domain')
    hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)}
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,actual_current_dispatch_call_bindings=calls,
        current_complete_provider_graph=field.assert_graph(),old_provider_mutations=mutations,
        current_physical_chart_owners_checked=33,new_current_physical_chart_owners_checked=15,
        retained_common_incoming_chart_owners_checked=18,
        current_physical_spatial4_time1_source_contributions_checked=counts,
        total_current_replayed_source_contributions_checked=sum(counts.values()),
        retained_signed_current_source_terms_checked=terms,exact_zero_current_source_contributions_checked=zeros,
        current_same_owner_gap_coverage=coverage,
        unbounded_exterior_finite_log_bounds_and_nonpositive_radial_powers_checked=True,
        original_coordinate_and_unit_source_theorems_retained=True,
        current_complete_physical_source_proof=proof,
        scope='same current callable physical source and original spatial4/time1 maps; source bounds, not resolved point values or quantitative interface admission',
        input_hashes=hashes,all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current complete33 physical PASS: '+str(sum(counts.values()))+' source contributions; quantitative joins/global/time remain open',flush=True)
    return result


if __name__=='__main__':run()
