"""Independent current-native scalar and actual selected incoming source check."""
import copy
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_pulse_native_numeric_identity as current


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Changed original native source was admitted')


@current.source_precision
def run(owner,raw):
    began=time.monotonic();owner.assert_graph()
    proof=owner.source_identity()
    assert proof==owner.proof==raw['source_identity'] and proof['passed']
    assert proof['exact_native_scalar_function_identities']==4
    assert proof['actual_buffer_and_constructor_identity_counts']==[7,12]
    assert proof['exact_true_and_numeric_factor_recipe_identities']==2
    assert proof['graph_replay_created_no_new_nodes']
    p=proof['production_Taylor_proof']
    assert p['power_inlet_exact_Taylor_rows']==24
    assert p['selected_incoming_exact_Taylor_rows']==18
    assert p['selected_and_refined_once_scaled_row_identities']==24
    assert p['energy_Z4_Z6_not_discarded'] and p['actual_production_assignments_executed']
    # These are live source/callback corruption checks; the original objects
    # are restored in finally blocks, and no ancestor artifacts are edited.
    saved=owner.selected.constants['C_E']
    try:
        owner.selected.constants['C_E']=saved*2
        rejects_scalar=rejected(owner.assert_graph)
    finally:owner.selected.constants['C_E']=saved
    callback=owner.selected.axial4.incoming
    try:
        owner.selected.axial4.incoming=lambda Z:None
        rejects_callback=rejected(owner.assert_graph)
    finally:del owner.selected.axial4.incoming
    assert owner.selected.axial4.incoming==callback
    selection=owner.selected.axial4.select
    try:
        owner.selected.axial4.select=lambda Z:None
        rejects_selection=rejected(owner.assert_graph)
    finally:del owner.selected.axial4.select
    assert owner.selected.axial4.select==selection
    owner.assert_graph()
    pins=dict(raw['input_hashes']);pins[current.NAME]=current.sha(current.NAME)
    for name,digest in pins.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,source_identity=proof,
        mutated_live_scalar_rejected=rejects_scalar,mutated_live_callback_rejected=rejects_callback,
        mutated_live_selection_callback_rejected=rejects_selection,
        source_identity_not_inferred_from_interval_overlap=True,
        exact_energy_Z4_Z6_and_Taylor_factorials_retained=True,
        pressure_error_bounds_not_transferred_by_this_guard=True,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),input_hashes=pins,
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('TRUE_PULSE_NATIVE_NUMERIC_SOURCE_JOIN_CHECKED',flush=True)
    return receipt
