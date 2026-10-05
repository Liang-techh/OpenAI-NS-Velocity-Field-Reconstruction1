"""Admit changed actual bridge parents with unchanged mixed4 operators."""
import ast
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import (
    CurrentActualBridgeMixedC4,VIEWS,GATE,SCOPES,REPLAY_PROOF,THEOREM,
    HERE,PREFIX,sha,current_coordinate_replay)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_actual_bridge_mixed_C4.json'


def parent_source_bindings():
    assignments=assignment_source_bindings('current_actual_bridge_mixed_C4','current_parents',{
        'incoming':'self.upstream.packet(Z,value,chart)',
        'inp':"self.upstream.prepare(Z)['inputs']",
        'phi':"native.IntervalTaylor(c,incoming['actual_phi_axial5'])",
        'own':"named_moments({name:native.IntervalTaylor(c,row) for name,row in incoming['actual_own_six_moments_axial5'].items()})"})
    known={expression:True for expression in (
        'self.upstream.comparison.macro(Z,value)','self.upstream.comparison.micro(Z,value)',
        "self.upstream.comparison_micro_cover(self.upstream.prepare(Z)['data'],lo,hi)")}
    tree=ast.parse((HERE/(PREFIX+'current_actual_bridge_mixed_C4.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentActualBridgeMixedC4')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='current_parents')
    for expression in known:
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(isinstance(n,ast.Call) and ast.dump(n)==wanted for n in ast.walk(fn))!=1:
            raise ValueError('Known comparison point/whole-chart acquisition changed')
    expected={'Uz_actual_axial5_coefficients':"incoming['actual_raw_V_axial5']",
        'actual_moment_shape_axial5_coefficients':'coefficient_lists(own)',
        'actual_Q_axial4_coefficients':"incoming['actual_radial_Q_axial4']",
        'pressure_axis_axial5_coefficients':"incoming['pressure_axis_axial5']",
        'pressure_increment_true_axial5_divided_by_R_F0_squared':"incoming['actual_pressure_increment_axial5_divided_by_R_F0_squared']"}
    for name,expression in expected.items():
        rows=[kw.value for node in ast.walk(fn) if isinstance(node,ast.Call)
            for kw in node.keywords if kw.arg==name]
        if len(rows)!=1 or ast.dump(rows[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Current actual bridge parent publication changed: '+name)
    return dict(assignments=assignments,known_comparison_calls=known,actual_history_publication=list(expected),
        actual_and_comparison_own_moments_acquired_separately=True,
        actual_axial5_and_comparison_axial6_preserved=True,passed=True)


def retained_feedback_bindings(field,raw):
    """Bind the changed ledger's actual parents to already admitted feedback."""
    name=PREFIX+'actual_bridge_integrals_check.json'
    old=json.loads((HERE/name).read_bytes());_verify_hashes(old)
    if (raw['input_hashes'].get(name)!=sha(name) or not old['all_passed']
            or not old['actual_own_six_moment_feedback_enclosures_available']
            or not old['full_finite_width_signed_bridge_integral_enclosures_available']
            or old['datum_enclosure_sha256']!=field.datum_sha):
        raise ValueError('Admitted actual nonlinear history evidence missing')
    packet=assignment_source_bindings('actual_bridge_integrals','packet',{
        'own':'own_moment_feedback(phi0,V0,initial,theta,delta_phi,history_delta_V)',
        'pressure':"dress(own['actual']['C'],p['inputs']['F0_squared_ratios'])"})
    macro=assignment_source_bindings('actual_bridge_integrals','packet',{
        'own':'self.macro_moments(p,q,inlet,delta_phi,history_delta_V)',
        'inlet':"self.packet(Z,2,'second',root)['actual_own_six_moments_axial5']"})
    micro=assignment_source_bindings('actual_bridge_integrals','own_moment_feedback',{
        'changes':'dict(H=2*delta_phi,M=delta_V,K=2*(phi0*delta_V+V0*delta_phi+delta_phi*delta_V),A=2*V0*delta_V+square(delta_V),B=2*phi0*delta_phi+square(delta_phi),C=2*phi0*delta_phi+square(delta_phi))',
        'actual[name]':'baseline[name]+feedback[name]'})
    quadratic=assignment_source_bindings('actual_bridge_integrals','macro_moments',{
        'quadratic[name]':"{'H':zero,'M':zero,'K':2*delta_phi*delta_V,'A':square(delta_V),'B':square(delta_phi),'C':square(delta_phi)}[name]*W",
        'actual[name]':'baseline[name]+feedback[name]'})
    return dict(retained_actual_history_receipt=name,receipt_sha256=sha(name),
        original_packet_micro_and_pressure_dressing=packet,original_packet_macro_history=macro,
        original_micro_nonlinear_feedback=micro,original_macro_quadratic_feedback=quadratic,
        retained_independent_original_scalar_integral_comparisons=old['independent_original_scalar_integral_comparisons'],
        retained_signed_actual_history_coefficients=old['finite_actual_field_and_own_moment_coefficients_checked'],
        derivative_admission_preserves_unselected_point_history_scope=True,passed=True)


def run(owner=None):
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentActualBridgeMixedC4(owner=owner,require_checked=False)
        manifest={k:v for k,v in raw.items() if k not in (
            'whole_current_bridge_mixed4','same_core_rectangular_extension_proofs',
            'factored_positive_width_amplitude_cap_proofs')}
        if encode(pack(field.manifest()))!=manifest:raise ValueError('Fresh current bridge graph/manifest differs')
        _,proof=current_coordinate_replay()
        feedback=retained_feedback_bindings(field,raw)
        if proof!=REPLAY_PROOF or raw['actual_parent_acquisition_replay']!=proof:
            raise ValueError('Changed bridge operator replay differs')
        if (not all(raw['current_provider_graph_identity'].values()) or raw[GATE]
                or raw['actual_bridge_mixed4_feedback_installed']
                or raw['current_three_bridge_charts_source_ownership_installed']
                or raw['core_bridge_and_R100_local_functional_joins_certified']
                or raw['R100_functional_join_to_existing_switch_installed']
                or raw['all_33_current_source_charts_physical_spatial4_time1_mapped']
                or any(raw[k] for k in SCOPES)):
            raise ValueError('Current bridge graph/admission scope differs')
        for fixture in raw['reused_unchanged_operator_fixtures'].values():
            if not fixture['passed']:raise ValueError('Unchanged bridge derivative evidence omitted')
        c=MPIntervalContext();c.dps=240
        def finite(row):
            value=read_interval(c,row)
            if not all(mp.isfinite(q) for q in endpoints(value)):raise ValueError('Nonfinite current bridge row')
            return value
        views=raw['whole_current_bridge_mixed4'];counts={};ledger_count=0;terms=0
        if set(views)!=set(VIEWS):raise ValueError('Original bridge domains/endpoints omitted')
        for name,packet in views.items():
            chart,coordinate=VIEWS[name];coordinate=coordinate if isinstance(coordinate,list) else [coordinate,coordinate]
            if (packet['chart']!=chart or endpoints(finite(packet['coordinate']))!=tuple(mp.mpf(q) for q in coordinate)
                    or endpoints(finite(packet['Z']))!=(-1,1) or packet[GATE]
                    or packet['actual_bridge_mixed4_feedback_installed']
                    or not packet['current_actual_bridge_mixed4_feedback_proved']
                    or not packet['current_coordinate_labelled_finite_width_history_used']
                    or packet['datum_enclosure_sha256']!=field.datum_sha
                    or any(packet[k] for k in SCOPES)):
                raise ValueError('Current bridge source/domain/scope differs')
            parent=packet['actual_parent_axial5_packet']
            original=parent['current_actual_integral_and_own_feedback_source']
            if parent['comparison_moments_substituted'] or parent['exact_coordinate_labelled_source']!='upstream.packet(Z,value,chart)':
                raise ValueError('Actual histories replaced by known comparison')
            for key in ('F_actual_over_F0_axial5_coefficients','Uz_actual_axial5_coefficients','pressure_axis_axial5_coefficients'):
                if len(parent[key])!=6:raise ValueError('Actual bridge axial5 source omitted')
                for row in parent[key]:finite(row)
            if set(original['actual_own_six_moments_axial5'])!=set('HMKABC'):
                raise ValueError('Actual six moment histories omitted')
            for row in original['actual_own_six_moments_axial5'].values():
                if len(row)!=6:raise ValueError('Own actual history axial5 omitted')
                for value in row:finite(value)
            if len(original['actual_radial_Q_axial4'])!=5:raise ValueError('Actual radial recovery axial4 omitted')
            known=packet['comparison_parent_axial6_packet']
            if len(known['phi'])!=7 or len(known['V'])!=7:raise ValueError('Known comparison axial6 omitted')
            if chart=='macro' and name=='R100_exit':
                if not original['coordinate_R100_endpoint'] or endpoints(finite(packet['R_enclosure_only']))!=(100,100):
                    raise ValueError('Exact labelled R100 interface lost')
            prefix='y' if chart=='macro' else 'phase'
            groups=('physical_velocity_pressure_'+prefix+'_Z_mixed4','physical_five_primitive_'+prefix+'_Z_mixed4')
            rows={row['physical_row']:row for row in packet['final_factored_physical_row_ledgers']}
            coordinate_prefix='y' if chart=='macro' else 's'
            indices={coordinate_prefix+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
            if len(rows)!=135:raise ValueError('Current full mixed4 ledger coverage omitted')
            count=0;logs=[finite(row) for row in packet['factored_source_log_bases']]
            for group in groups:
                expected=4 if group==groups[0] else 5
                if len(packet[group])!=expected:raise ValueError('Velocity/pressure/primitive component omitted')
                for component,grid in packet[group].items():
                    if set(grid)!=indices:raise ValueError('Mixed4 graded derivative indices omitted')
                    for index,bound in grid.items():
                        finite(bound);order=int(index.split('_')[0][1:])
                        ledger=rows[component+'/'+index];term_logs=[]
                        for term in ledger['terms']:
                            coefficient=finite(term['final_ordinary_coefficient']);size=max(abs(v) for v in endpoints(coefficient))
                            if not size:continue
                            exponents=list(term['source_exponents']);exponents[0]-=order if chart!='macro' else 0
                            term_logs.append(sum((logs[i]*p for i,p in enumerate(exponents) if p),c.mpf(0))+c.ln(c.mpf(size)))
                            terms+=1
                        logindex='y'+index[1:]
                        published=packet['physical_logR_Z_mixed4_log_bound_ledger'][group][component][logindex]
                        if bool(published['exact_zero'])!=(not term_logs) or not published['formal_width_powers_cancelled_before_log_bound']:
                            raise ValueError('Original source-width cancellation omitted')
                        if term_logs:
                            required=max(endpoints(row)[1] for row in term_logs)+endpoints(c.ln(len(term_logs)))[1]
                            actual=endpoints(finite(published['log_absolute_upper']))[1]
                            tolerance=max(mp.mpf(1),abs(required))*mp.mpf('1e-220')
                            if actual+tolerance<required:raise ValueError('Current logR derivative ledger upper underestimated')
                        elif published['log_absolute_upper'] is not None:raise ValueError('Exact source zero lost')
                        count+=1;ledger_count+=1
            counts[name]=count
        # One fresh unsaved Z tests direct source acquisition, rather than a
        # chart JSON lookup. Original operator fixtures are retained by hash.
        fresh={chart:field.evaluate('0.271',value,chart) for chart,value in (
            ('first','0.4'),('second','1.6'),('macro','0.37'))}
        for chart,packet in fresh.items():
            if packet['chart']!=chart or not packet['current_coordinate_labelled_finite_width_history_used']:
                raise ValueError('Fresh current bridge source call failed')
        for chart,value in (('first',-1),('first','1.1'),('second',0),('second',3),('macro','1.1'),('axis',0)):
            try:field.evaluate(0,value,chart)
            except ValueError:pass
            else:raise ValueError('Original bridge chart domain widened')
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_actual_bridge_mixed4_rows_checked=counts,
            total_current_bridge_mixed4_rows_checked=sum(counts.values()),
            current_final_factored_logR_rows_checked=ledger_count,current_nonzero_source_terms_checked=terms,
            actual_parent_source_bindings=parent_source_bindings(),
            actual_nonlinear_history_and_pressure_source_bindings=feedback,
            unchanged_original_derivative_and_physical_replay=REPLAY_PROOF,
            reused_unchanged_operator_fixtures=raw['reused_unchanged_operator_fixtures'],retained_operator_receipt=THEOREM,
            fresh_unsaved_Z_each_chart_checked=True,current_three_bridge_charts_source_ownership_installed=True,
            **{GATE:True},actual_bridge_mixed4_feedback_installed=True,
            core_bridge_and_R100_local_functional_joins_certified=False,
            R100_functional_join_to_existing_switch_installed=False,
            all_33_current_source_charts_physical_spatial4_time1_mapped=False,
            **dict.fromkeys(SCOPES,False),all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current actual three-chart bridge mixed4:1215 source/ledger rows; full points/interfaces remain open',flush=True)
    return result


if __name__=='__main__':run()
