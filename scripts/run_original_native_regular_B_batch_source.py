"""Freeze and enter one regular-B source batch exactly once.

The prepare action is the only action that reads the kernel's literal
configuration.  The batch action consumes that frozen packet, authenticates
the raw inputs, writes its entry checkpoint, imports the new kernel, and
invokes ``derive`` once.  Inspecting a saved batch only checks checkpoint
metadata and the append-only diagnostic journal.
"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import gzip
import hashlib
import json
import os
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
RUNNER = 'scripts/run_original_native_regular_B_batch_source.py'
KERNEL = 'src/openai_ns_reconstruction/original_native_regular_B_batch_source.py'
PARENT = (
    'experiments/root_st073/'
    'lei_ren_part1_paper_compliant_current_original_native_J3_endpoint_transport_v3.json'
)
PARENT_SHA = 'c798cc171350236483e4e5ee4e35ba2bcc2dc4d6b714a491d926eb149175f7d8'
PARENT_COMMIT = '64ce196b8b57cdbd9494333d8f7903c8d6e7c592'
CK = ROOT / 'tmp/original_native_regular_B_batch_source'
COMPONENT = 'regular_B_batch'
INPUT_SCHEMA = 'native_regular_B_batch_inputs_v1'
DIAGNOSTIC_SCHEMA = 'native_regular_B_batch_stage_diagnostic_v1'
INHERITED_SOURCE_COUNT = 1867
FINAL_SOURCE_COUNT = 1870
FAILED_ARTIFACT_COUNT = 26
OLD_CHECKPOINT_COUNT = 140
PARENT_CURRENT_CHECKPOINT_COUNT = 4
NEW_CHECKPOINT_COUNT = 4
FINAL_CHECKPOINT_COUNT = 148

# These are raw-only source families.  The names are intentionally narrow:
# ordinary first-switch source definitions remain available to the new leaf.
RAW_ONLY_MARKERS = (
    'resume_first_switch_2026_10_17',
    'first_switch_history',
    'true_first_switch_history',
    'rm_loop_density_identity',
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
        allow_nan=False).encode('utf-8')).hexdigest()


def source_path(name):
    candidate = ROOT / name
    return candidate if candidate.exists() else ROOT / 'experiments/root_st073' / name


def _raw_only(name):
    normalized = name.replace('-', '_').lower()
    return any(marker in normalized for marker in RAW_ONLY_MARKERS)


def _resolve_literal(node, assignments):
    """Resolve the deliberately small literal grammar used by the kernel."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        require(node.id in assignments, 'Unknown literal name: ' + node.id)
        return _resolve_literal(assignments[node.id], assignments)
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        values = [_resolve_literal(item, assignments) for item in node.elts]
        return values if not isinstance(node, ast.Tuple) else tuple(values)
    if isinstance(node, ast.Dict):
        return {
            _resolve_literal(key, assignments): _resolve_literal(value, assignments)
            for key, value in zip(node.keys, node.values)
        }
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return _resolve_literal(node.left, assignments) + _resolve_literal(node.right, assignments)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _resolve_literal(node.operand, assignments)
        return -value if isinstance(node.op, ast.USub) else value
    raise ValueError('Nonliteral kernel configuration: ' + ast.dump(node))


def literal_configuration():
    """Read only top-level literal assignments from the new kernel.

    This function is called by ``prepare`` only.  It never imports or executes
    the kernel.
    """
    tree = ast.parse((ROOT / KERNEL).read_text(encoding='utf-8'))
    assignments = {}
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name):
                    assignments[target.id] = statement.value
        elif isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name):
            assignments[statement.target.id] = statement.value
    names = ('SCHEMA', 'RECORD_NAMES', 'RECORD_HASHES', 'STAGE_FLAGS', 'FALSE_FLAGS')
    return {name: _resolve_literal(assignments[name], assignments) for name in names}


def _read_record(name, expected, decode=True):
    """Authenticate raw bytes before optional gzip or text decoding."""
    path = source_path(name)
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    require(actual == expected, 'Changed raw source: ' + name)
    record = {
        'raw_sha256': actual,
        'raw_authenticated_before_decode': True,
        'is_gzip': name.endswith('.gz'),
        'source_path': path.relative_to(ROOT).as_posix(),
        'raw_only': not decode,
    }
    if not decode:
        # The forbidden families are intentionally represented by their raw
        # authentication only.  No gzip stream or source text is loaded.
        record['decoded_sha256'] = None
        return record
    decoded = gzip.decompress(raw) if name.endswith('.gz') else raw
    record.update(
        source_text=decoded.decode('utf-8'),
        decoded_sha256=hashlib.sha256(decoded).hexdigest(),
    )
    return record


def _read_raw_hash(name, expected):
    """Authenticate a preserved artifact without decoding its bytes."""
    path = source_path(name)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    require(actual == expected, 'Changed preserved raw artifact: ' + name)
    return actual


def _write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def paths():
    return tuple(CK / ('new_regular_B_batch_source_' + kind + '.json')
                 for kind in ('source', 'entered', 'result'))


def diagnostic_path():
    return CK / 'new_regular_B_batch_source_diagnostics.jsonl'


def _parent_scope(parent):
    require(parent.get('schema') == 'original_native_J3_endpoint_transport_v3_stage_v1',
            'Parent report schema changed')
    require(parent.get('source_only') is True
            and parent.get('accepted_identity') is True
            and parent.get('full_goal_complete') is False,
            'Parent report is not the accepted source-only J3 publication')
    pins = parent.get('source_pins')
    require(isinstance(pins, dict) and len(pins) == INHERITED_SOURCE_COUNT,
            'Parent source inventory changed')
    require(all(isinstance(name, str) and isinstance(value, str)
                and len(value) == 64 for name, value in pins.items()),
            'Parent source inventory hashes changed')
    require(parent.get('input_source_pin_count') == INHERITED_SOURCE_COUNT,
            'Parent input source count changed')
    component_result = parent.get('component_result', {})
    require(isinstance(component_result, dict)
            and component_result.get('passed') is True
            and component_result.get('source_only') is True
            and component_result.get('numeric_values_computed') is False,
            'Parent component result scope changed')
    require(parent.get('ordinary_Z_orders') == [0, 1]
            and parent.get('ordinary_Z_not_temporal_n') is True
            and parent.get('owned_Z_cell') == ['0', '1/2'],
            'Parent ordinary-Z scope changed')
    remaining_nodes = parent.get('remaining_B_identity_nodes')
    require(isinstance(remaining_nodes, list) and remaining_nodes
            and all(isinstance(node, int) for node in remaining_nodes)
            and parent.get('all_twelve_remaining_B_identities_open') is True,
            'Parent remaining-B scope changed')
    require(parent.get('old_proof_functions_called') == 0
            and parent.get('old_operators_replayed') is False,
            'Parent old-consumer scope changed')
    failed = parent.get('preserved_failed_artifact_hashes_union')
    require(isinstance(failed, dict) and len(failed) == FAILED_ARTIFACT_COUNT
            and all(isinstance(name, str) and isinstance(value, str)
                    and len(value) == 64 for name, value in failed.items()),
            'Parent failed-artifact inventory changed')
    return pins


def _parent_checkpoint_inventory(parent):
    """Authenticate the parent 140 old plus 4 current checkpoints distinctly."""
    old_inventory = parent.get('old_checkpoint_inventory')
    require(isinstance(old_inventory, dict), 'Parent old checkpoint inventory missing')
    old_hashes = old_inventory.get('authenticated_checkpoint_hashes')
    current_hashes = parent.get('checkpoint_hashes')
    require(isinstance(old_hashes, dict) and len(old_hashes) == OLD_CHECKPOINT_COUNT,
            'Parent old checkpoint count changed')
    require(isinstance(current_hashes, dict)
            and len(current_hashes) == PARENT_CURRENT_CHECKPOINT_COUNT,
            'Parent current checkpoint count changed')
    require(all(isinstance(name, str) and isinstance(value, str)
                and len(value) == 64 for name, value in
                {**old_hashes, **current_hashes}.items()),
            'Parent checkpoint hashes changed')
    require(set(old_hashes).isdisjoint(current_hashes)
            and len(set(old_hashes) | set(current_hashes)) ==
            OLD_CHECKPOINT_COUNT + PARENT_CURRENT_CHECKPOINT_COUNT,
            'Parent old/current checkpoint sets are not distinct')
    for name, expected in {**old_hashes, **current_hashes}.items():
        _read_raw_hash(name, expected)
    return {
        'authenticated_checkpoint_hashes': dict(old_hashes),
        'parent_current_checkpoint_hashes': dict(current_hashes),
        'authenticated_old_count': OLD_CHECKPOINT_COUNT,
        'parent_current_count': PARENT_CURRENT_CHECKPOINT_COUNT,
        'distinct_parent_count': OLD_CHECKPOINT_COUNT + PARENT_CURRENT_CHECKPOINT_COUNT,
        'raw_sha_only': True,
        'history_decode_forbidden': True,
    }


def _static_scope(config):
    require(config['SCHEMA'] == 'original_native_regular_B_batch_source',
            'Kernel schema changed')
    names = list(config['RECORD_NAMES'])
    hashes = dict(config['RECORD_HASHES'])
    require(names and len(names) == len(set(names))
            and all(isinstance(name, str) for name in names),
            'Kernel record names are not a unique literal list')
    require(set(hashes) == set(names), 'Kernel record names/hashes mismatch')
    require(all(isinstance(value, str) and len(value) == 64 for value in hashes.values()),
            'Kernel record hashes are not SHA-256 literals')
    require(isinstance(config['STAGE_FLAGS'], (list, tuple))
            and bool(config['STAGE_FLAGS'])
            and all(isinstance(flag, str) for flag in config['STAGE_FLAGS'])
            and isinstance(config['FALSE_FLAGS'], (list, tuple))
            and all(isinstance(flag, str) for flag in config['FALSE_FLAGS']),
            'Kernel stage flags are not literal names')
    require(len(config['STAGE_FLAGS']) == len(set(config['STAGE_FLAGS']))
            and len(config['FALSE_FLAGS']) == len(set(config['FALSE_FLAGS'])),
            'Kernel stage flags are duplicated')
    require(set(config['STAGE_FLAGS']).isdisjoint(config['FALSE_FLAGS']),
            'Kernel stage and false flags overlap')
    return names, hashes


def _prepare_packet():
    """Build the immutable source packet; this is the prepare-only path."""
    config = literal_configuration()
    names, expected_hashes = _static_scope(config)
    parent_record = _read_record(PARENT, PARENT_SHA)
    parent = json.loads(parent_record['source_text'])
    inherited = _parent_scope(parent)
    require(PARENT not in inherited
            and RUNNER not in inherited and KERNEL not in inherited,
            'Parent source union already contains a new publication')

    current_code_hashes = {RUNNER: sha(ROOT / RUNNER), KERNEL: sha(ROOT / KERNEL)}
    final_union = dict(inherited)
    final_union[PARENT] = PARENT_SHA
    final_union.update(current_code_hashes)
    require(len(final_union) == FINAL_SOURCE_COUNT,
            'Final source union must be parent report plus code2')

    allowed = set(final_union)
    require(set(names) <= allowed, 'Kernel record escapes the frozen source union')
    require(PARENT in names, 'Parent report missing from kernel records')
    # Authenticate every inherited source pin before reading any configured
    # record.  This is a raw-byte inventory pass; it never decodes a source.
    for name, expected in inherited.items():
        _read_raw_hash(name, expected)
    for name, expected in expected_hashes.items():
        if name in inherited:
            require(expected == inherited[name], 'Inherited record hash changed: ' + name)
        elif name == PARENT:
            require(expected == PARENT_SHA, 'Parent report hash changed')
        elif name in current_code_hashes:
            require(expected == current_code_hashes[name], 'Current code hash changed: ' + name)

    records = {}
    for name in names:
        records[name] = parent_record if name == PARENT else _read_record(
            name, expected_hashes[name], decode=not _raw_only(name))
    failed = parent['preserved_failed_artifact_hashes_union']
    for name, expected in failed.items():
        _read_raw_hash(name, expected)
    checkpoint_inventory = _parent_checkpoint_inventory(parent)
    input_digest = digest(records)
    decoded_hashes = {name: item['decoded_sha256'] for name, item in records.items()}
    raw_only_records = [name for name, item in records.items() if item.get('raw_only')]
    static_receipt = {
        'ast_configuration_read': True,
        'kernel_imported': False,
        'derive_calls': 0,
        'old_producer_calls': 0,
        'old_operator_calls': 0,
        'old_proof_calls': 0,
        'numeric_values_computed': False,
        'runtime_model_metadata': None,
    }
    return dict(
        schema=INPUT_SCHEMA,
        component=COMPONENT,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA,
        parent_commit=PARENT_COMMIT,
        input_digest=input_digest,
        record_names=names,
        record_hashes=expected_hashes,
        record_raw_sha256={name: item['raw_sha256'] for name, item in records.items()},
        record_decoded_sha256=decoded_hashes,
        decoded_record_hashes=decoded_hashes,
        record_source_paths={name: item['source_path'] for name, item in records.items()},
        raw_only_records=raw_only_records,
        inherited_source_hashes=dict(inherited),
        inherited_source_pin_count=INHERITED_SOURCE_COUNT,
        current_code_hashes=current_code_hashes,
        implementation_hashes=current_code_hashes,
        source_pin_manifest={
            'inherited_source_hashes': dict(inherited),
            'inherited_source_pin_count': INHERITED_SOURCE_COUNT,
            'parent_report_hashes': {PARENT: PARENT_SHA},
            'parent_commit': PARENT_COMMIT,
            'current_code_hashes': dict(current_code_hashes),
            'record_raw_sha256': {
                name: item['raw_sha256'] for name, item in records.items()
            },
            'record_decoded_sha256': dict(decoded_hashes),
            'input_digest': input_digest,
            'new_pins': [PARENT, RUNNER, KERNEL],
            'new_pin_count': 3,
            'final_union_hashes': final_union,
            'final_union_count': FINAL_SOURCE_COUNT,
            'parent_report_separate_publication': True,
        },
        preserved_failed_artifact_hashes=dict(failed),
        preserved_failed_artifact_count=FAILED_ARTIFACT_COUNT,
        parent_checkpoint_inventory=checkpoint_inventory,
        new_checkpoint_count=NEW_CHECKPOINT_COUNT,
        final_checkpoint_count=FINAL_CHECKPOINT_COUNT,
        schema_expected=config['SCHEMA'],
        stage_flags=list(config['STAGE_FLAGS']),
        false_flags=list(config['FALSE_FLAGS']),
        ordinary_Z_orders=[0, 1],
        ordinary_Z_not_temporal_n=True,
        identity_domain_scope='common domain of original source definitions',
        source_only=True,
        ancestor_operators_replayed=False,
        old_producer_function_calls=0,
        old_operator_calls=0,
        old_proof_functions_called=0,
        numeric_values_computed=False,
        full_goal_complete=False,
        static_receipt=static_receipt,
        diagnostics_policy='append-only JSONL; snapshot, flush, and fsync before callback return',
        diagnostics_checkpoint_additional_to_source_entered_result=True,
    )


def prepare():
    source, entered, target = paths()
    journal = diagnostic_path()
    require(not any(path.exists() for path in (source, entered, target, journal)),
            'Batch checkpoints already exist; never rerun')
    packet = _prepare_packet()
    _write_new(source, packet)
    return {
        'passed': True,
        'prepared_only': True,
        'source_sha256': sha(source),
        'code_sha256': dict(packet['current_code_hashes']),
        'input_digest': packet['input_digest'],
        'source_pin_count': packet['source_pin_manifest']['final_union_count'],
        'static_receipt': packet['static_receipt'],
    }


def _load_packet(source):
    packet = json.loads(source.read_bytes())
    require(packet.get('schema') == INPUT_SCHEMA
            and packet.get('component') == COMPONENT,
            'Frozen source packet schema changed')
    require(packet.get('source_pin_manifest', {}).get('final_union_count') == FINAL_SOURCE_COUNT,
            'Frozen source union count changed')
    require(packet.get('new_checkpoint_count') == NEW_CHECKPOINT_COUNT,
            'Frozen checkpoint count changed')
    manifest = packet.get('source_pin_manifest', {})
    require(set(packet.get('current_code_hashes', {})) == {RUNNER, KERNEL}
            and packet.get('implementation_hashes') == packet['current_code_hashes'],
            'Frozen implementation inventory changed')
    require(manifest.get('input_digest') == packet.get('input_digest')
            and manifest.get('record_raw_sha256') == packet.get('record_raw_sha256')
            and manifest.get('record_decoded_sha256') == packet.get('record_decoded_sha256'),
            'Frozen source pin manifest changed')
    for name, expected in packet['current_code_hashes'].items():
        require(name in (RUNNER, KERNEL) and sha(ROOT / name) == expected,
                'Current implementation changed after prepare: ' + name)
    return packet


def _load_frozen_records(packet):
    """Re-authenticate packet inputs without consulting kernel literals."""
    inherited = packet['inherited_source_hashes']
    for name, expected in inherited.items():
        _read_raw_hash(name, expected)
    for name, expected in packet['preserved_failed_artifact_hashes'].items():
        _read_raw_hash(name, expected)
    records = {}
    raw_only = set(packet.get('raw_only_records', []))
    require(raw_only == {name for name in packet['record_names'] if _raw_only(name)},
            'Frozen raw-only source policy changed')
    for name in packet['record_names']:
        records[name] = _read_record(
            name, packet['record_hashes'][name], decode=name not in raw_only)
        require(records[name]['raw_sha256'] == packet['record_raw_sha256'][name],
                'Frozen raw record changed: ' + name)
        require(records[name]['decoded_sha256'] == packet['record_decoded_sha256'][name],
                'Frozen decoded record changed: ' + name)
    require(digest(records) == packet['input_digest'], 'Frozen input digest changed')
    return records


def _validate_result(result, packet):
    require(isinstance(result, dict)
            and result.get('schema') == packet['schema_expected']
            and result.get('component') == COMPONENT
            and result.get('source_only') is True
            and type(result.get('passed')) is bool,
            'Kernel result schema or component changed')
    require(result.get('input_digest') == packet['input_digest'],
            'Kernel result input digest changed')
    require(result.get('ordinary_Z_orders') == [0, 1]
            and result.get('ordinary_Z_not_temporal_n') is True,
            'Kernel ordinary-Z scope changed')
    for flag in packet['false_flags']:
        require(result.get(flag) is False, 'Unadmitted scope changed: ' + flag)
    stage_values = [result.get(flag) for flag in packet['stage_flags']]
    require(all(type(value) is bool for value in stage_values),
            'Kernel stage flags are not booleans')
    declared = result.get('stage_flags')
    if declared is not None:
        require(isinstance(declared, dict)
                and all(declared.get(flag) is result.get(flag)
                        for flag in packet['stage_flags']),
                'Kernel declared stage flag map changed')
    all_stages = all(stage_values)
    require((not all_stages) or result['passed'],
            'All declared kernel stages cannot be true on a failed result')
    gates = result.get('unresolved_gates')
    require(isinstance(gates, list), 'Kernel unresolved gates missing')
    if result['passed']:
        require(all_stages and not gates, 'Passed result lacks declared stage completion')
    else:
        require(not all_stages or bool(gates), 'Failed result lacks a failed stage or gate')
    scope = result.get('B_identity_scope', {})
    if isinstance(scope, dict):
        for key in ('newly_proved_B_identity_nodes', 'prospective_B_identity_nodes',
                    'unresolved_B_identity_nodes'):
            values = scope.get(key)
            if values is not None:
                require(isinstance(values, list)
                        and all(isinstance(node, int) and 524 <= node <= 535
                                for node in values),
                        'B identity scope escaped nodes 524 through 535')
    for key in ('old_proof_functions_called', 'old_operator_calls',
                'old_operator_function_calls', 'source_producer_function_calls',
                'old_producer_function_calls'):
        if key in result:
            require(result[key] == 0, 'Old producer/operator/proof call count changed: ' + key)
    require(result.get('old_operators_replayed', False) is False
            and result.get('numeric_values_computed') is False,
            'Unadmitted replay or numeric evaluation')
    audit = result.get('source_inspection_audit')
    if isinstance(audit, dict):
        require(audit.get('new_source_consumer_calls') == 1
                and audit.get('accepted_ancestor_operator_calls') == 0
                and audit.get('old_proof_functions_called', 0) == 0,
                'Source consumer call audit changed')
    require(result.get('runtime_model_metadata') is None,
            'Runtime model metadata must remain null without a tool result')
    remaining = result.get('remaining_construction_gates')
    require(isinstance(remaining, list) and bool(remaining),
            'Open construction gates were dropped')


def _empty_journal(path):
    with path.open('xb') as stream:
        stream.flush()
        os.fsync(stream.fileno())


def batch():
    source, entered, target = paths()
    journal = diagnostic_path()
    require(source.exists() and not entered.exists() and not target.exists()
            and not journal.exists(),
            'Frozen batch is absent or already entered')
    packet = _load_packet(source)
    _write_new(entered, {
        'schema': 'native_regular_B_batch_entered_v1',
        'component': COMPONENT,
        'source_sha256': sha(source),
        'input_digest': packet['input_digest'],
        'implementation_hashes': packet['current_code_hashes'],
        'kernel_sha256': packet['current_code_hashes'][KERNEL],
        'entered_at_utc': datetime.now(timezone.utc).isoformat(),
        'kernel_imported_before_entry': False,
        'derive_calls_before_entry': 0,
    })
    _empty_journal(journal)
    events = []

    def persist_diagnostic(stage, payload):
        require(isinstance(stage, str) and bool(stage), 'Invalid diagnostic stage')
        frozen = json.loads(json.dumps(
            payload, sort_keys=True, ensure_ascii=True, allow_nan=False))
        event = {
            'schema': DIAGNOSTIC_SCHEMA,
            'index': len(events),
            'stage': stage,
            'input_digest': packet['input_digest'],
            'payload': frozen,
        }
        raw = (json.dumps(event, sort_keys=True, ensure_ascii=True,
                          allow_nan=False) + '\n').encode('utf-8')
        with journal.open('ab') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        events.append(event)

    started = time.perf_counter()
    caught = False
    records = None
    try:
        records = _load_frozen_records(packet)
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            'new_native_regular_B_batch_source', ROOT / KERNEL)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result = module.derive(
            records, expected_input_digest=packet['input_digest'],
            diagnostic_callback=persist_diagnostic)
    except Exception as error:
        caught = True
        result = {
            'schema': packet['schema_expected'],
            'component': COMPONENT,
            'passed': False,
            'source_only': True,
            'failed_attempt_preserved': True,
            'error_type': type(error).__name__,
            'error': str(error),
            'traceback': ''.join(traceback.format_exception(error)),
            'input_digest': packet['input_digest'],
        }
    if not isinstance(result, dict):
        caught = True
        result = {
            'schema': packet['schema_expected'],
            'component': COMPONENT,
            'passed': False,
            'source_only': True,
            'failed_attempt_preserved': True,
            'error_type': 'TypeError',
            'error': 'derive returned a non-mapping result',
            'traceback': '',
            'input_digest': packet['input_digest'],
        }
    result.update(
        execution_seconds=time.perf_counter() - started,
        source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered),
        runner_caught_exception=caught,
        runner_full_goal_complete=False,
        runtime_model_metadata=None,
        stage_diagnostics=list(events),
        diagnostic_events=list(events),
        diagnostic_checkpoint_path=journal.relative_to(ROOT).as_posix(),
        diagnostic_checkpoint_sha256=sha(journal),
        diagnostic_events_saved_before_validation=len(events),
        full_result_saved_before_validation=True,
    )
    # This is the one and only result write.  Validation below cannot replace
    # or repair an already saved result.
    _write_new(target, result)
    if records is not None:
        require(digest(records) == packet['input_digest'],
                'Input mutated during regular-B consumer')
    if not caught:
        _validate_result(result, packet)
    return {
        'passed': result['passed'],
        'result_sha256': sha(target),
        'code_sha256': dict(packet['current_code_hashes']),
        'seconds': result['execution_seconds'],
        'error': result.get('error'),
        'diagnostic_events': len(events),
    }


def inspect_batch():
    """Inspect saved metadata only; never decode inputs or call the kernel."""
    source, entered, target = paths()
    journal = diagnostic_path()
    require(all(path.exists() for path in (source, entered, target, journal)),
            'Saved batch checkpoints are incomplete')
    packet = json.loads(source.read_bytes())
    entry = json.loads(entered.read_bytes())
    result = json.loads(target.read_bytes())
    require(packet.get('schema') == INPUT_SCHEMA
            and packet.get('component') == COMPONENT,
            'Saved packet metadata changed')
    require(entry.get('source_sha256') == sha(source)
            and entry.get('input_digest') == packet.get('input_digest')
            and entry.get('implementation_hashes') == packet.get('current_code_hashes')
            and entry.get('kernel_imported_before_entry') is False,
            'Saved entry metadata changed')
    require(result.get('source_checkpoint_sha256') == sha(source)
            and result.get('entry_checkpoint_sha256') == sha(entered)
            and result.get('input_digest') == packet.get('input_digest')
            and result.get('diagnostic_checkpoint_path') ==
            journal.relative_to(ROOT).as_posix()
            and result.get('diagnostic_checkpoint_sha256') == sha(journal),
            'Saved result checkpoint links changed')
    lines = journal.read_text(encoding='utf-8').splitlines()
    events = [json.loads(line) for line in lines]
    require(all(event.get('schema') == DIAGNOSTIC_SCHEMA
                and event.get('input_digest') == packet.get('input_digest')
                and event.get('index') == index
                for index, event in enumerate(events)),
            'Saved diagnostic journal metadata changed')
    require(result.get('stage_diagnostics') == events
            and result.get('diagnostic_events') == events
            and result.get('full_result_saved_before_validation') is True,
            'Saved diagnostic/result snapshots changed')
    if result.get('runner_caught_exception'):
        require(result.get('passed') is False
                and result.get('failed_attempt_preserved') is True
                and isinstance(result.get('traceback'), str),
                'Saved exception envelope changed')
    return {
        'passed': True,
        'saved_only': True,
        'component_passed': result.get('passed'),
        'full_goal_complete': False,
        'code_sha256': packet.get('current_code_hashes'),
        'static_receipt': packet.get('static_receipt'),
        'diagnostic_events': len(events),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare_batch', 'batch', 'inspect_batch'))
    action = parser.parse_args().action
    output = prepare() if action == 'prepare_batch' else (
        batch() if action == 'batch' else inspect_batch())
    print(json.dumps(output, ensure_ascii=True, allow_nan=False))
