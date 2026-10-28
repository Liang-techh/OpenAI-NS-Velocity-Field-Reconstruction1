"""Bind conditional main/gap moment algebra to saved source roles.

No old operator, source-value evaluator, quadrature, graph traversal or
ancestor checker is invoked. Source function semantics, substitution capture
and ordinary-Z transport remain explicit hypotheses, not certified facts.
"""
from __future__ import annotations
import base64
import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

P='lei_ren_part1_paper_compliant_current_original_'
RP=P+'Rp_native_C5_pulse_source.json.gz'
KERNEL=P+'n1_pulse_kernel_functions.json.gz'
C3=P+'n1_pulse_pressure_C3_functions.json.gz'
PRESS_CHECK=P+'n1_pulse_pressure_functions_check.json'
KERNEL_CHECK=P+'n1_pulse_kernel_functions_check.json'
C3_CHECK=P+'n1_pulse_pressure_C3_functions_check.json'
RP_CHECK=P+'Rp_native_C5_pulse_source_check.json'
PIN={RP:'442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
 KERNEL:'9584201b517fe3058ea6d1a169ddb7d319776eb97479b36024af4bffd4973b75',
 C3:'bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d',
 PRESS_CHECK:'caff410405d7c20018df89d5c98e7c69d4bb1271719919ff54eb12aface12749',
 KERNEL_CHECK:'47343a839590c3cb08c8d726ef225415e15f4fb12a5f7cf888049ea1f20cffc0',
 C3_CHECK:'cbac4685fc9031a2ed5c5a7ef624b438e1a4ce3f32b85b569512aecb1d38feb8',
 RP_CHECK:'320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15'}
SEAM_REPORT='experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_m15_source_bound_seams_v2.json'
SEAM_REPORT_SHA='375b7a5d5135c587478cb088d9b90c236fbd95f1d72cc501385a0b03a03318a7'
SEAM_RESULT_SHA='d30e0d78ce4542587156194952b15865b583aa4da8f8570d8354c64e29651a0f'
FIELDS=('R_squared_q','R_squared_q_squared','R_squared_U_squared','R_q_y')


def _hashed(path, expected):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('input hash mismatch: '+path.name)
    return raw


def _gzip(folder,name):
    return json.loads(gzip.decompress(_hashed(folder/name,PIN[name])))


def _report(root):
    report=json.loads(_hashed(root/SEAM_REPORT,SEAM_REPORT_SHA))
    raw=base64.b64decode(report['artifacts']['tmp/frame_m15_source_bound_seams_v2/result.json']['base64'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=SEAM_RESULT_SHA:raise ValueError('saved seam result hash mismatch')
    data=json.loads(raw)
    if data.get('status')!='returned' or data['result'].get('schema')!='frame_m15_source_bound_seams_v2':raise ValueError('unexpected saved seam')
    return data['result']


def recover(root: Path):
    folder=root/'experiments/root_st073'
    checks={name:json.loads(_hashed(folder/name,PIN[name])) for name in (PRESS_CHECK,KERNEL_CHECK,C3_CHECK,RP_CHECK)}
    if any(record.get('all_passed') is not True for record in checks.values()):raise ValueError('source receipt has not passed')
    for name in (PRESS_CHECK,C3_CHECK,RP_CHECK):
        if checks[name]['input_hashes'].get(KERNEL)!=PIN[KERNEL]:raise ValueError('kernel pin differs: '+name)
    press=checks[PRESS_CHECK]['checks']
    selected=press['selected_rows'];actual=press['actual_source']
    for key in ('pointwise_selected_C_u_plus_ap_v_row_identities_exact','finite_mu_second_row_and_explicit_absolute_determinant_proved','original_two_moment_joins_are_function_identities'):
        if selected.get(key) is not True:raise ValueError('missing selected-row premise: '+key)
    for key in ('actual_selected_C1_handles_replayed_from_original_selector','actual_selected_C3_value_and_first_Z_handles_preserved','gap_nonzero_pressure_history_terms_retained'):
        if actual.get(key) is not True:raise ValueError('missing saved actual-source premise: '+key)
    rp=_gzip(folder,RP);kernel=_gzip(folder,KERNEL)
    if rp['source_family']!=kernel['source_family']:raise ValueError('source-family metadata mismatch')
    if any(record.get('source_family') != rp['source_family'] for record in checks.values()):raise ValueError('source receipt family differs')
    if rp['input_hashes'].get(KERNEL)!=PIN[KERNEL] or rp['input_hashes'].get(C3)!=PIN[C3]:raise ValueError('Rp source input pins differ')
    frame=rp['Rp_native_C5_pulse_source_bindings'];h=frame['inputs'][C3];w=kernel['exact_pulse_kernel_functions']
    mapping=frame['input_substitution'];rewrite=frame['rewritten_nodes']
    old_native=h['native']['m1'];new_native=frame['current_frame']['m1']
    if len(old_native)<3 or len(new_native)<3:raise ValueError('missing native m1 rows through Z2')
    for old,new in zip(old_native[:3],new_native[:3]):
        if old in (0,1):
            if old != new:raise ValueError('global constant native row changed')
        elif mapping.get(str(old))!=new:raise ValueError('saved m1 input row mapping mismatch')
    old_ap=h['selected_pulse']['ap'];new_ap=frame['selected_pulse']['ap']
    if len(old_ap)<3 or len(new_ap)<3:raise ValueError('missing selected ap rows through Z2')
    for old,new in zip(old_ap[:3],new_ap[:3]):
        if rewrite.get(str(old))!=new:raise ValueError('saved ap rewrite row mapping mismatch')
    seam=_report(root);occurrences=seam['occurrences'];roles=[]
    for endpoint,chart,index in (('main.upper','main',1),('gap.lower','gap',0)):
        for field in FIELDS:
            for order in (0,1):
                found=[x for x in occurrences if x['endpoint']==endpoint and x['field']==field and x['z_order']==order]
                if len(found)!=1:raise ValueError('saved seam role count mismatch')
                original=h['endpoints'][chart][index][field][order]
                if found[0]['original_ref']!=original:raise ValueError('saved seam original role mismatch')
                roles.append({'endpoint':endpoint,'field':field,'ordinary_Z_order':order,
                    'original_ref':original,'source_ref':found[0]['source_ref']})
    # This module is a new independent symbolic leaf; importing/calling it
    # does not invoke any earlier builder, checker, or source evaluator.
    path=root/'src/openai_ns_reconstruction/frame_main_gap_moment_algebra.py'
    spec=importlib.util.spec_from_file_location('frame_main_gap_moment_algebra_new_leaf',path)
    algebra=importlib.util.module_from_spec(spec);sys.modules[spec.name]=algebra;spec.loader.exec_module(algebra)
    proof=algebra.recover(root)
    return {'schema':'frame_main_gap_source_proof_v1',
      'status':'conditional_algebra_bound_to_saved_source_roles',
      'algebra':proof,'algebra_passed':proof.get('passed') is True,
      'source_roles':{'kernel':{k:w[k] for k in ('mu','lambda_rows','first_row','divided_row','scaled_main_rows','incoming_factors','log_common_prefactor','log_end_scale','original_gp','original_gp_variable')},
          'original_native_m1':old_native,'current_native_m1':new_native,
          'original_selected_ap':old_ap,'current_selected_ap':new_ap,
          'original_selected_C1':h['selected_pulse']['C1'],
          'original_selected_C2':h['selected_pulse']['C2'],
          'current_geometry':{k:frame['geometry'][k] for k in ('mu','delta','Pstar','Rp')},
          'main_gap_endpoint_occurrences':roles},
      'saved_receipt_scope':{'selected_rows':selected,'actual_source':actual,
          'direct_physical_fixture_only':press['direct_physical'].get('independent_mathematical_fixture') is True},
      'source_raw_hashes':{RP:PIN[RP],KERNEL:PIN[KERNEL],'seam_report':SEAM_REPORT_SHA,'seam_result':SEAM_RESULT_SHA},
      'proof_conditions':[
          'actual source definitions implement the symbolic W,H,E and first-row premises without changing physical scaling',
          'the original gp and gp_y vanish at xi=11 for all Z',
          'moment identities hold as C2 functions in the same Z coordinate',
          'the chart radii and theta functions coincide at xi=11 and radius is Z-independent',
          'inlet substitution is capture-free, preserves source algebra and commutes with ordinary Z derivatives',
          'L=1-delta*Z^2 is nonzero on the admitted domain'],
      'symbolic_source_inputs_are_actual_values':False,'source_function_semantics_authenticated_by_this_leaf':False,
      'inlet_substitution_capture_free_proved':False,'ordinary_Z_transport_proved':False,
      'actual_frame_main_gap_endpoint_identity_proved':False,'physical_seam_identity_proved':False,
      'source_values_computed':False,'integrals_evaluated':False,'numeric_values_computed':False,
      'old_builders_called':False,'old_operators_replayed':False,'old_error_terms_removed':False,
      'extra_ledger_terms_added':False,'admissible_for_numeric_use':False,'full_goal_complete':False}
