"""Read the single saved substitution-variable descriptor, without traversal."""
import gzip
import hashlib
import json

ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz'
ARCHIVE_SHA = '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
PRIOR = 'tmp/frame_replacement_free_names_v2/result.json'
PRIOR_SHA = 'c21c27b02a04f6f909d5cd531053021c2ec50d1e62d4f7c4cc54cffc859ef928'


def recover(root):
    saved_raw = (root / PRIOR).read_bytes()
    if hashlib.sha256(saved_raw).hexdigest() != PRIOR_SHA:
        raise ValueError('v2 result hash changed')
    saved = json.loads(saved_raw)
    unknown = saved['result']['replacement_analysis']['unknown']
    if {(x['node'], x.get('variable')) for x in unknown} != {(265191, '263344'), (265204, '263344')}:
        raise ValueError('v2 diagnostic targets differ')
    raw = (root / ARCHIVE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != ARCHIVE_SHA:
        raise ValueError('raw gzip hash mismatch before decode')
    nodes = json.loads(gzip.decompress(raw))['exact_graph_nodes']
    rows = []
    for index in (265191, 265204):
        row = nodes[index]
        if row.get('operation') != 'function_substitution' or row.get('variable') != 263344:
            raise ValueError('saved substitution descriptor differs')
        rows.append({'node': index, 'operation': row['operation'], 'variable': row['variable']})
    descriptor = nodes[263344]
    if not isinstance(descriptor, dict):
        raise ValueError('descriptor is not a node object')
    # Copy only this one row. No referenced operand or closure is inspected.
    return {'schema': 'frame_binder_descriptor_v1', 'status': 'descriptor_read_without_closure',
            'source_archive_sha256': ARCHIVE_SHA, 'prior_v2_result_sha256': PRIOR_SHA,
            'substitution_rows': rows, 'descriptor_node': 263344, 'descriptor': descriptor,
            'descriptor_operands_traversed': False, 'geometry_closure_traversed': False,
            'source_values_computed': False, 'ordinary_Z_transport_proved': False,
            'actual_capture_free_proved': False, 'actual_capture_proved': False,
            'physical_endpoint_proved': False, 'runtime_model_metadata': None}
