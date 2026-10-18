"""Single corrected local source call; consume the frozen attempt as data."""
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
RUNNER='scripts/run_original_Rm_local_m_C2_access.py'
OLD_RUNNER='scripts/run_original_Rm_m_C2_C5_source_bridge.py'
OLD_KERNEL='src/openai_ns_reconstruction/original_Rm_local_m_C2_regularity.py'
KERNEL='src/openai_ns_reconstruction/original_Rm_local_m_C2_regularity_access.py'
OLD_RUNNER_SHA='18308a501dc688f41ce78ad2465afda0217dd4a51eac24c3d7e350cc5c8fdfe9'
OLD_KERNEL_SHA='18ad986aba59f91a770757edbdab80539adb41060be1a41c38cd9d5baffcc65f'
FAILED_SHA='d55b48b06206819bd47d45376d4e0abcd1bf4d84704e4f4514db016d3c5c7b3d'
CK=ROOT/'tmp/original_Rm_local_m_C2_access'
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
    failed_path=old.paths('local')[2]
    require(sha(failed_path)==FAILED_SHA,'failed result changed')
    failed=json.loads(failed_path.read_bytes())
    require(failed['passed'] is False and failed['failed_attempt_preserved'] and failed['error_type']=='TypeError'
        and 'argument after ** must be a mapping, not list' in failed['error'],'not the saved call-shape failure')
    _,verified,stored=old.saved_component('local')
    require(verified==prior and stored==failed,'failed checkpoint chain changed')
    before=(ROOT/OLD_KERNEL).read_text()
    expected=before.replace('M_chain_ok = all(_expect_row(nodes, handle, operation, gates, **fields)\n                         for handle, operation, fields in expected["M_chain"])',
        'M_chain_ok = all(_expect_row(nodes, handle, operation, gates, **(fields if isinstance(fields, Mapping) else {"arguments": fields}))\n                         for handle, operation, fields in expected["M_chain"])')
    require(expected!=before and (ROOT/KERNEL).read_text()==expected,'corrected module changed beyond the M chain call')
    for name in (RUNNER,KERNEL):
        raw=(ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'),'LF/EOF changed')
        ast.parse(raw,filename=name)
    packet=dict(schema='Rm_local_m_C2_corrected_source_inputs_v1',prior_source=prior,
        failed_result_sha256=FAILED_SHA,input_digest=old.digest(values),graph_digest=old.GRAPH_SHA,
        implementation_hashes={name:sha(ROOT/name) for name in (RUNNER,KERNEL,OLD_RUNNER,OLD_KERNEL)},
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
