"""Check connected local density product roots after frozen diagnostics."""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
RUNNER='scripts/run_original_Rm_local_m_C2_products.py'
OLD_RUNNER='scripts/run_original_Rm_m_C2_C5_source_bridge.py'
OLD_KERNEL='src/openai_ns_reconstruction/original_Rm_local_m_C2_regularity_access.py'
KERNEL='src/openai_ns_reconstruction/original_Rm_local_m_C2_products.py'
OLD_RUNNER_SHA='18308a501dc688f41ce78ad2465afda0217dd4a51eac24c3d7e350cc5c8fdfe9'
OLD_KERNEL_SHA='ce45480452c9825268d40345040a8690b41405eff8ec4afd3d249251bc6f779d'
FAILED_SHA='3df7b1d3b04af01ad10914e356bdf85473d612148ee863e0664e8f95968a6455'
CK=ROOT/'tmp/original_Rm_local_m_C2_products'
SOURCE,ENTRY,RESULT=(CK/name for name in ('new_source.json','new_entered.json','new_result.json'))


def require(condition,message):
    if not condition: raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec)
    sys.modules[name]=value
    spec.loader.exec_module(value)
    return value


def setup():
    require(sha(ROOT/OLD_RUNNER)==OLD_RUNNER_SHA and sha(ROOT/OLD_KERNEL)==OLD_KERNEL_SHA,'frozen failed code changed')
    old=module('frozen_m_source_bridge_data_reader',ROOT/OLD_RUNNER)
    parent,nodes,values,prior=old.setup('local')
    previous=ROOT/'tmp/original_Rm_local_m_C2_access'
    failed_path=previous/'new_result.json'
    require(sha(failed_path)==FAILED_SHA,'failed result changed')
    failed=json.loads(failed_path.read_bytes())
    require(failed['passed'] is False and not failed.get('failed_attempt_preserved')
        and {g.get('handle') for g in failed['unresolved_gates'] if g.get('handle')}=={121079,197899},
        'not the frozen product-root diagnostic')
    previous_source,previous_entry=previous/'new_source.json',previous/'new_entered.json'
    prior_packet,entry=json.loads(previous_source.read_bytes()),json.loads(previous_entry.read_bytes())
    require(failed['source_checkpoint_sha256']==entry['source_sha256']==sha(previous_source)
        and failed['entry_checkpoint_sha256']==sha(previous_entry)
        and entry['kernel_sha256']==OLD_KERNEL_SHA and prior_packet['input_digest']==old.digest(values),
        'previous diagnostic checkpoint chain changed')
    previous_runner='scripts/run_original_Rm_local_m_C2_access.py'
    require(sha(ROOT/previous_runner)=='51e8ca52c99aa3a98b849f99c4bd340da0ac4a8f67ac858acd28aece3d41d71d'
        and sha(ROOT/KERNEL)=='2eda8c1ecd21d7275b9b838f25f01825b135906d7d976ce6eee8e41d17439926',
        'frozen previous runner or new product proof changed')
    for name in (RUNNER,KERNEL):
        raw=(ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'),'LF/EOF changed')
        ast.parse(raw,filename=name)
    packet=dict(schema='Rm_local_m_C2_corrected_source_inputs_v1',prior_source=prior,
        failed_result_sha256=FAILED_SHA,previous_source_sha256=sha(previous_source),previous_entry_sha256=sha(previous_entry),input_digest=old.digest(values),graph_digest=old.GRAPH_SHA,
        implementation_hashes={name:sha(ROOT/name) for name in (RUNNER,KERNEL,OLD_RUNNER,OLD_KERNEL,previous_runner)},
        accepted_mathematical_operators_replayed=False,failed_local_source_not_replayed=True,full_goal_complete=False)
    return old,parent,nodes,values,packet


def prepare():
    old,*_,packet=setup()
    require(not any(p.exists() for p in (SOURCE,ENTRY,RESULT)),'corrected attempt exists')
    old.write_new(SOURCE,packet)
    return dict(passed=True,source_sha256=sha(SOURCE),input_digest=packet['input_digest'])


def run():
    old,_,nodes,values,packet=setup()
    require(json.loads(SOURCE.read_bytes())==packet and not ENTRY.exists() and not RESULT.exists(),'corrected attempt entered or changed')
    kernel=module('new_corrected_actual_local_m_C2',ROOT/KERNEL)
    old.write_new(ENTRY,dict(source_sha256=sha(SOURCE),input_digest=packet['input_digest'],
        kernel_sha256=sha(ROOT/KERNEL),entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started,result=time.perf_counter(),None
    try:
        result=kernel.derive(nodes,**values)
        require(old.digest(nodes)==old.GRAPH_SHA and old.digest(values)==packet['input_digest'],'source input mutated')
        old.validate_result(result,dict(schema_expected=kernel.SCHEMA))
    except Exception as error:
        result=dict(passed=False,failed_attempt_preserved=True,error_type=type(error).__name__,error=str(error),
            traceback=''.join(traceback.format_exception(error)),rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started,source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTRY),runner_full_goal_complete=False)
    old.write_new(RESULT,result)
    return dict(passed=result['passed'],sha256=sha(RESULT),seconds=result['execution_seconds'],
        error=result.get('error'),unresolved_gates=result.get('unresolved_gates',[]))


def saved():
    old,parent,_,_,packet=setup()
    result,entry=json.loads(RESULT.read_bytes()),json.loads(ENTRY.read_bytes())
    require(json.loads(SOURCE.read_bytes())==packet and result['source_checkpoint_sha256']==entry['source_sha256']==sha(SOURCE)
        and result['entry_checkpoint_sha256']==sha(ENTRY) and entry['input_digest']==packet['input_digest']
        and entry['kernel_sha256']==sha(ROOT/KERNEL),'saved corrected chain changed')
    if not result.get('failed_attempt_preserved'):
        old.validate_result(result,dict(schema_expected='original_Rm_local_m_C2_regularity_source_v1'))
    return parent,packet,result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','inspect'))
    action=parser.parse_args().action
    print(json.dumps(dict(passed=True,saved_only=True,source_passed=saved()[2]['passed'],full_goal_complete=False)
        if action=='inspect' else globals()[action]()))
