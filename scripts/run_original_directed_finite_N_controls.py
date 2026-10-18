"""New finite-N propagation; saved inspection never executes math providers."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys
import time

from mpmath.libmp import mpf_cmp


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'experiments/root_st073'
LOOP = 'lei_ren_part1_paper_compliant_current_original_directed_loop_controls.json'
LOOP_SHA = 'adefdae708754a2572eae15a90662274bf17e88a0d72d5db9befbc9e2858b252'
NATIVE = 'lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json'
NATIVE_SHA = '10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a'
NAME = 'lei_ren_part1_paper_compliant_current_original_directed_finite_N_controls.json'
REPORT = BASE/NAME
PENDING = BASE/(NAME+'.pending')
MODULE = 'src/openai_ns_reconstruction/original_directed_finite_N_controls.py'
KERNEL = 'src/openai_ns_reconstruction/original_directed_finite_N_kernel.py'
RUNNER = 'scripts/run_original_directed_finite_N_controls.py'
AUDIT = 'original_directed_finite_N_incoming_audit.json'
FIXTURE = 'original_directed_finite_N_kernel_fixture.json'
REVIEWED = {MODULE:'2fadf55155639e1d3a7403885c2aca5c5740d02df6f84b61cb8bd5ef94d2293e',
            KERNEL:'a7a89dd2b9ca28e7ef7e4f412a61a0e8496d777633bfb86ccebbe31b3f4818fc',
            'experiments/root_st073/'+AUDIT:'41c59b0872848935e20b7b015ef537fe54a66881045b3afe19789e326b1822e4',
            'experiments/root_st073/'+FIXTURE:'93f05c761807ce235ff892973b2caaf7223ee2416140bae5988fb007f4d79bdc'}
CHARTS = ('long_reshape','reference','restoration','post_restore')
SIX = ('q','q_Z','q_squared','U_q_Z','U_squared','U_U_Z')
FALSE_SCOPE = ('physical_y_bounds','full_pressure_revalidated','global_five_moment_closure',
               'full_goal_complete','legacy_packets_recertified','inherited_global_guard_admitted')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):
            digest.update(block)
    return digest.hexdigest()


def input_path(name):
    path = Path(name)
    require(not path.is_absolute(),'Absolute source path')
    if '..' in path.parts:
        require(name == '../../src/openai_ns_reconstruction/schedule_pressure.py','Source traversal')
        return (BASE/path).resolve()
    return ROOT/path if path.parts[0] == 'experiments' else BASE/path.name


def parents():
    loop_raw, native_raw = (BASE/LOOP).read_bytes(), (BASE/NATIVE).read_bytes()
    require(len(loop_raw)==1752702 and hashlib.sha256(loop_raw).hexdigest()==LOOP_SHA,'Loop parent changed')
    require(len(native_raw)==42177499 and hashlib.sha256(native_raw).hexdigest()==NATIVE_SHA,'Native parent changed')
    loop, native = json.loads(loop_raw), json.loads(native_raw)
    require(loop['scope']['fixed_target_phase_Z0_through_Z3_verified'] is True,'New loop premise missing')
    require(loop['parent_report_sha256']==NATIVE_SHA,'Loop/native binding changed')
    for name,digest in loop['input_hashes'].items():
        require(sha(input_path(name))==digest,'Input changed: '+name)
    for name,digest in loop['implementation_hashes'].items():
        path = Path(name)
        require(not path.is_absolute() and '..' not in path.parts,'Unsafe implementation path')
        require(sha(ROOT/path)==digest,'Implementation changed: '+name)
    for path,digest in REVIEWED.items():
        require(digest!='REVIEW_PENDING' and sha(ROOT/path)==digest,'New source review pending/changed: '+path)
    for path in (MODULE,KERNEL,RUNNER):
        raw=(ROOT/path).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'),
                'New Python must use LF and a single EOF newline: '+path)
    return loop,native,loop_raw,native_raw


def endpoint(value):
    require(isinstance(value,(list,tuple)) and len(value)==4
            and all(type(x)is int for x in value),'Malformed tuple')
    sign,mantissa,exponent,bits = value
    require(sign in (0,1) and mantissa>=0 and bits==mantissa.bit_length(),'Nonfinite/malformed endpoint')
    require((mantissa==0 and tuple(value)==(0,0,0,0)) or (mantissa>0 and mantissa&1),'Noncanonical endpoint')
    return tuple(value)


def record(value):
    require(isinstance(value,dict) and value.get('directed')is True,'Directed record required')
    require(mpf_cmp(endpoint(value['lower_exact_mpf_tuple']),endpoint(value['upper_exact_mpf_tuple']))<=0,
            'Reversed interval')
    return 1


def records(tree):
    if isinstance(tree,dict):
        if 'lower_exact_mpf_tuple'in tree or 'upper_exact_mpf_tuple'in tree:
            return record(tree)
        return sum(records(value) for value in tree.values())
    if isinstance(tree,list):
        return sum(records(value) for value in tree)
    return 0


def rows(values,length):
    require(isinstance(values,list) and len(values)==length,'Row count changed')
    for k,row in enumerate(values):
        require(row['ordinary_Z_order']==k,'Ordinary row order changed')
        if 'value'in row:
            record(row['value'])
        else:
            require(type(row['exact_zero'])is bool,'Exact-zero marker missing')
            require(row['log_absolute_upper_bound']is None if row['exact_zero']
                    else isinstance(row['log_absolute_upper_bound'],dict),'Bad log row')
            if not row['exact_zero']:
                record(row['log_absolute_upper_bound'])


def below(values,threshold):
    return all(row['exact_zero'] or mpf_cmp(endpoint(row['log_absolute_upper_bound']['upper_exact_mpf_tuple']),
                endpoint(threshold['lower_exact_mpf_tuple']))<0 for row in values)


def cap_record(row):
    return None if row.get('exact_zero') is True else row['log_absolute_upper_bound']


def same_endpoints(left,right):
    if left is None or right is None:
        return left is right
    return all(endpoint(left[key])==endpoint(right[key]) for key in
               ('lower_exact_mpf_tuple','upper_exact_mpf_tuple'))


def guards(loop):
    path = 'src/openai_ns_reconstruction/original_restoration_finite_N_I1_numeric.py'
    require(sha(ROOT/path)==loop['implementation_hashes'][path],'Source formula file changed')
    tree = ast.parse((ROOT/path).read_bytes(),filename=path)
    names = ('_recover_delta_Q_logs','_difference_logs','_compose_D_error_logs','_boundary_difference_logs')
    functions = {node.name:node for node in tree.body if isinstance(node,ast.FunctionDef)}
    return {name:dict(path=path,source_sha256=sha(ROOT/path),
        AST_sha256=hashlib.sha256(ast.dump(functions[name],include_attributes=False).encode()).hexdigest(),
        source_executed=False) for name in names}


def validate_incoming_audit(loop,native,audit):
    evidence=audit['hash_ast_row_audit']
    for item in evidence['files']:
        path=item['path']
        expected=loop['implementation_hashes'].get(path,
            loop['input_hashes'].get(path,loop['input_hashes'].get(Path(path).name)))
        if path=='experiments/root_st073/'+LOOP:
            expected=LOOP_SHA
        require(expected==item['sha256']==sha(ROOT/path)
                and (ROOT/path).stat().st_size==item['bytes'],'Incoming evidence pin changed: '+path)
    count=0
    for item in evidence['ast_entries']:
        tree=ast.parse((ROOT/item['path']).read_bytes(),filename=item['path'])
        for name,digest in item['functions'].items():
            node=tree
            for part in name.split('.'):
                matches=[child for child in node.body
                         if isinstance(child,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)) and child.name==part]
                require(len(matches)==1,'Ambiguous incoming AST path: '+name)
                node=matches[0]
            require(hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest()==digest,
                    'Incoming AST evidence changed: '+item['path']+':'+name)
            count+=1
    require(count==73,'Incoming AST frontier count changed')
    incoming=native['actual_switch_post_finite_N_I1_numeric']['charts']['post_switch']['physical_error_log_bounds']['delta_M']
    digest=hashlib.sha256(json.dumps(incoming,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    require(digest==audit['saved_delta_M']['rows_canonical_json_sha256'],'Incoming exact-row digest changed')
    require(audit['scope']['producer_builder_owner_checker_selfcheck_recurrence_quadrature_executed'] is False,
            'Incoming audit replayed mathematics')
    return dict(files=len(evidence['files']),AST_functions=count,exact_incoming_rows=4,
                global_historical_certification=False)


def physical_geometry(native):
    result={}
    for name in CHARTS:
        if name in CHARTS[:2]:
            path='actual_long_reference_background_I1_numeric.charts.'+name
            baseline=native['actual_long_reference_background_I1_numeric']['charts'][name]
            require(baseline['W_formula']=='W=R_upper^2; radius not materialized','Long/reference physical W changed')
            require(mpf_cmp(endpoint(baseline['length']['lower_exact_mpf_tuple']),(0,0,0,0))>0,'Nonpositive chart length')
            domain=dict(length=baseline['length'])
        else:
            path='actual_restoration_post_jet_repair.corrected_background.charts.'+name
            baseline=native['actual_restoration_post_jet_repair']['corrected_background']['charts'][name]
            formula='Wrestore=Rrestore^2; radius not materialized' if name=='restoration' else 'Wpost=Rm^2; radius not materialized'
            require(baseline['W_formula']==formula and baseline['coordinate_domain']==['0','1']
                    and baseline['physical_measure']=='R^2*dcoordinate','Restoration physical geometry changed')
            domain=dict(coordinate_domain=baseline['coordinate_domain'],physical_measure=baseline['physical_measure'])
        result[name]=dict(source_native_report_sha256=NATIVE_SHA,source_path=path,
            saved_W_formula=baseline['W_formula'],saved_normalization=baseline['normalization'],
            saved_log_R_upper=baseline['log_R_upper'],saved_log_W=baseline['log_W'],domain=domain,
            source_geometry_reused_without_owner_replay=True,
            measure_identity='R^2 dy; R=R_upper*exp(y-y_upper)',
            normalized_mass_bound='integral exp(2*(y-y_upper))dy <= 1/2',
            endpoint_pair_bound='(R_left^2+R_right^2)/W <= 2',
            radius_or_W_materialized=False)
    return result


def validate(loop,native,report):
    require(report['parent_report']==LOOP and report['parent_report_sha256']==LOOP_SHA
            and report['native_report']==NATIVE and report['native_report_sha256']==NATIVE_SHA,'Parent binding changed')
    inputs = dict(loop['input_hashes'])
    inputs.update({LOOP:LOOP_SHA,AUDIT:REVIEWED['experiments/root_st073/'+AUDIT],
                   FIXTURE:REVIEWED['experiments/root_st073/'+FIXTURE]})
    implementations = dict(loop['implementation_hashes'])
    implementations.update({path:sha(ROOT/path) for path in (MODULE,KERNEL,RUNNER)})
    require(report['input_hashes']==inputs and report['implementation_hashes']==implementations,'Pin maps changed')
    require(report['source_guards']==guards(loop),'Source formula guards changed')
    require(report['physical_geometry']==physical_geometry(native),'Physical geometry source receipt changed')
    require(report['incoming_audit']==json.loads((BASE/AUDIT).read_bytes()),'Incoming audit changed')
    require(report['incoming_audit_check']==validate_incoming_audit(loop,native,report['incoming_audit']),
            'Incoming frontier receipt changed')
    require(report['fresh_fixture']==json.loads((BASE/FIXTURE).read_bytes()),'Fresh fixture changed')
    fixture=report['fresh_fixture']
    require(fixture['module_sha256']==REVIEWED[KERNEL]
            and all(fixture['fixture'][key] is True for key in
                    ('pressure_Q_native','pressure_Q_correction','boundary_difference_endpoint_pair',
                     'D_difference','explicit_None_exact_zero'))
            and all(value is True for value in fixture['fixture']['six_differences'].values()),
            'Fresh rational algebra checks failed or are not source-bound')
    require(report['historical_payload_retained_by_parent_hash']is True,'Preservation link missing')
    require(all(report['scope'][key]is False for key in FALSE_SCOPE),'Downstream scope overclaimed')
    packet = report['actual_directed_finite_N_controls']
    require(set(packet['charts'])==set(CHARTS),'Wrong chart set')
    expected_incoming = native['actual_switch_post_finite_N_I1_numeric']['charts']['post_switch']['physical_error_log_bounds']['delta_M']
    require(packet['incoming_delta_M_source']['rows']==expected_incoming,'Inherited memory reset or changed')
    all_ok = True
    for name in CHARTS:
        chart = packet['charts'][name]
        if name in CHARTS[:2]:
            baseline=native['actual_long_reference_background_I1_numeric']['charts'][name]
            source=native['actual_native_directed_controls_v2']['raw_interval_tree']['long_reference']['charts'][name]
            logPE=source['inverse_Utheta_detail']['log_Utheta_upper']
        else:
            baseline=native['actual_restoration_post_jet_repair']['corrected_background']['charts'][name]
            source=native['actual_native_directed_controls_v2']['raw_interval_tree']['restoration']['charts'][name]
            logPE=source['scales']['log_Utheta_upper']
        centres=chart['baseline_centres']
        for key in ('six_integrals','coefficients','B_over_W','D','S1_over_W'):
            require(centres[key]==baseline['arithmetic'][key],'Signed background changed: '+name+'.'+key)
        require(same_endpoints(chart['physical_error_log_bounds']['log_Pstar_E_upper'],logPE),
                'Physical amplitude changed: '+name)
        rows(chart['incoming_delta_M'],4)
        for original,packed in zip(expected_incoming,chart['incoming_delta_M']):
            require(same_endpoints(cap_record(original),cap_record(packed)),
                    'Nonreset memory edge changed: '+name)
        for field,count in (('delta_U',4),('delta_M',4),('delta_Q',3)):
            rows(chart['physical_error_log_bounds'][field],count)
        outgoing=chart['physical_error_log_bounds']['delta_M']
        for incoming,local,result in zip(chart['incoming_delta_M'],chart['physical_error_log_bounds']['delta_U'],outgoing):
            candidates=[cap_record(row) for row in (incoming,local) if cap_record(row) is not None]
            result=cap_record(result)
            if candidates:
                maximum=endpoint(candidates[0]['upper_exact_mpf_tuple'])
                for candidate in candidates[1:]:
                    value=endpoint(candidate['upper_exact_mpf_tuple'])
                    if mpf_cmp(value,maximum)>0:
                        maximum=value
                require(result is not None and endpoint(result['lower_exact_mpf_tuple'])==maximum
                        and endpoint(result['upper_exact_mpf_tuple'])==maximum,'Memory max cap changed: '+name)
            else:
                require(result is None,'Zero memory was replaced')
        expected_incoming=outgoing
        for key in SIX:
            rows(chart['six_I1_difference_error_log_bounds'][key],2)
        for key in ('B_difference_error_log_bounds','D_difference_error_log_bounds','S1_difference_error_log_bounds'):
            rows(chart[key],2)
        threshold = chart['I1_error_threshold']
        record(threshold)
        ok = all(below(chart['six_I1_difference_error_log_bounds'][key],threshold) for key in SIX)
        ok = ok and below(chart['B_difference_error_log_bounds'],threshold)
        require(chart['target_2_power_minus1024_met']is ok,'Threshold claim does not match exact tuples')
        require(chart['physical_error_log_bounds']['logN']==native['actual_native_directed_controls_v2']['raw_interval_tree']['shared']['logN'],
                'Current N changed')
        require(chart['normalized_by']==dict(long_reshape='Wlong',reference='Wref',restoration='Wrestore',post_restore='Wpost')[name],
                'Wrong own-W normalization')
        if ok:
            require(all(chart[key]is not None for key in ('I1_full','B_full','D_full','S1_full')),'Missing full rows')
            for key in SIX:
                rows(chart['I1_full']['six_integrals'][key],2)
            for key in ('B_full','D_full','S1_full'):
                rows(chart[key]['rows'],2)
            for key in ('I1_full','B_full'):
                radius=chart[key]['symmetric_radius_added_once']
                record(radius)
                require(endpoint(radius['lower_exact_mpf_tuple'])==endpoint(radius['upper_exact_mpf_tuple'])==(0,1,-1024,1),
                        'One-time finite-N radius changed')
            require(chart['D_full']['separate_D_error_not_added'] is True
                    and chart['S1_full']['no_extra_radius'] is True,'D/S1 radius double counted')
        else:
            require(all(chart[key]is None for key in ('I1_full','B_full','D_full','S1_full')),'Unsupported tight full enclosure')
        all_ok = all_ok and ok
    require(packet['all_normalized_I1_and_B_errors_below_2_power_minus1024']is all_ok,'Aggregate threshold changed')
    require(report['execution_audit']['accepted_ancestors_replayed']is False,'Accepted math replay')
    return records(packet)


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write-report',action='store_true')
    mode.add_argument('--inspect-report',action='store_true')
    args = parser.parse_args()
    loop,native,loop_raw,native_raw = parents()
    if args.write_report:
        require(not REPORT.exists() and not PENDING.exists(),'Saved/pending data exists; mathematical replay refused')
        sys.path.insert(0,str(ROOT/'src'))
        from mpmath.ctx_iv import MPIntervalContext
        from openai_ns_reconstruction import original_directed_finite_N_controls as provider
        c=MPIntervalContext()
        c.prec=1600
        started=time.monotonic()
        packet=provider.compute(c,loop_raw,native_raw)
        inputs=dict(loop['input_hashes'])
        inputs.update({LOOP:LOOP_SHA,AUDIT:REVIEWED['experiments/root_st073/'+AUDIT],
                       FIXTURE:REVIEWED['experiments/root_st073/'+FIXTURE]})
        implementations=dict(loop['implementation_hashes'])
        implementations.update({path:sha(ROOT/path) for path in (MODULE,KERNEL,RUNNER)})
        report=dict(stage='current_original_directed_finite_N_controls',parent_report=LOOP,
            parent_report_sha256=LOOP_SHA,native_report=NATIVE,native_report_sha256=NATIVE_SHA,
            input_hashes=inputs,implementation_hashes=implementations,source_guards=guards(loop),
            physical_geometry=physical_geometry(native),
            incoming_audit=json.loads((BASE/AUDIT).read_bytes()),fresh_fixture=json.loads((BASE/FIXTURE).read_bytes()),
            incoming_audit_check=validate_incoming_audit(loop,native,json.loads((BASE/AUDIT).read_bytes())),
            actual_directed_finite_N_controls=packet,historical_payload_retained_by_parent_hash=True,
            scope=dict(new_four_chart_error_transport_computed=True,
                current_N_I1_and_B_target_met=packet['all_normalized_I1_and_B_errors_below_2_power_minus1024'],
                incoming_cap_retained_as_inherited_premise=True,full_inherited_cap_revalidation=False,
                **{key:False for key in FALSE_SCOPE}),
            execution_audit=dict(accepted_ancestors_replayed=False,only_new_error_inputs_propagated=True),
            execution_seconds=time.monotonic()-started)
        PENDING.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
        count=validate(loop,native,report)
        PENDING.rename(REPORT)
    else:
        report=json.loads(REPORT.read_bytes())
        count=validate(loop,native,report)
    print(json.dumps(dict(passed=True,report=str(REPORT),sha256=sha(REPORT),bytes=REPORT.stat().st_size,
        canonical_records=count,input_pins=len(report['input_hashes']),implementation_pins=len(report['implementation_hashes']),
        all_normalized_I1_and_B_errors_below_2_power_minus1024=report['actual_directed_finite_N_controls']['all_normalized_I1_and_B_errors_below_2_power_minus1024'],
        scope=report['scope'],execution_seconds=report['execution_seconds'])))


if __name__=='__main__':
    main()
