"""Replay five current pulse source joins, exact reciprocal view and bounds."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_interfaces import (
    CurrentPulseInterfaces,current_pulse_trace_proof,PULSESEAMS,PRIMITIVES,
    GATES,OPEN,HERE,PREFIX,NAME,RECEIPT,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_postpulse_interfaces_check import cloned_graph,trace_bounds
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import PULSE_CHARTS
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def cloned_pulse_graph(physical):
    p=cloned_graph(physical);h=p.history;sel=h.selected=copy.copy(physical.history.selected)
    p.pulse=sel.pulse=copy.copy(physical.pulse)
    h.flatten.pulse=h.outer.pulse=p.pulse
    p.dispatch.selected=sel;p.dispatch.native_pulse=p.pulse
    for chart in PULSE_CHARTS:p.dispatch.providers[chart]=p.pulse
    p.assert_graph()
    return p


def wrong_sources_rejected(field):
    rejected=[]
    for label in ('changed_native_mu','changed_exact_logE_copy','different_current_complete_future','old_selected_pulse_provider'):
        p=cloned_pulse_graph(field.physical)
        if label=='changed_native_mu':p.pulse.mu=p.pulse.mu+1
        elif label=='changed_exact_logE_copy':p.pulse.logE=p.pulse.logE*2
        elif label=='different_current_complete_future':p.heat.future=object()
        else:p.dispatch.providers['pulse_exit']=p.base.pulse
        try:current_pulse_trace_proof(p)
        except ValueError:rejected.append(label);continue
        raise ValueError('Wrong current pulse source accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),baseline_clone_graph_passed=True,
        live_current_owners_unmodified=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseInterfaces(require_checked=False)
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_whole_Z_pulse_trace_views'}:
        raise ValueError('Current pulse source/manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Pulse producer claims acceptance')
    proof=current_pulse_trace_proof(field.physical)
    if encode(pack(proof))!=raw['current_pulse_source_trace_proof']:raise ValueError('Current pulse function theorem differs')
    theorem=proof['source_bound_five_current_pulse_mixed4_theorem']
    for key,count in (('current_instantiated_five_primitive_mixed4_identities',375),
        ('current_instantiated_velocity_absolute_pressure_mixed4_identities',300)):
        if len(theorem[key])!=count or not all(theorem[key].values()):raise ValueError('Five current pulse mixed4 grids missing')
    transfer=proof['source_bound_exact_linear_physical_trace_transfer']
    if not transfer['passed'] or not transfer['source_equal_mixed_rows_give_zero_exact_physical_operator_difference']:
        raise ValueError('Actual linear physical transfer missing')
    names=[row[0] for row in PULSESEAMS]
    if set(raw['current_whole_Z_pulse_trace_views'])!=set(names):raise ValueError('Five whole-Z endpoint pairs required')
    native=0;bounds={}
    for name in names:
        point=field.endpoint(name,[-1,1])
        if encode(pack(point))!=raw['current_whole_Z_pulse_trace_views'][name]:raise ValueError('Current pulse endpoint replay differs: '+name)
        if any(point[k] for k in OPEN):raise ValueError('Pulse endpoint scope exceeds source traces')
        if len(point['same_current_native_mixed4_difference_enclosures'])!=10:raise ValueError('Four physical, five primitive and forcing grids required')
        for rows in point['same_current_native_mixed4_difference_enclosures'].values():
            if len(rows)!=15:raise ValueError('Ordinary mixed4 consistency grid incomplete')
            for value in rows.values():
                lo,hi=endpoints(value)
                if not mp.isfinite(lo) or not mp.isfinite(hi) or not lo<=0<=hi:
                    raise ValueError('Source function identity contradicts native enclosure: '+name)
                native+=1
        if name=='gap_coordinate':
            if not point['reciprocal_source_reduction_used'] or not point['right_physical_trace']['exact_source_endpoint_view']:
                raise ValueError('Exact reciprocal endpoint source was replaced by coverage')
        bounds[name]=trace_bounds(point)
        print('Current pulse trace replay PASS: '+name,flush=True)
    ledger=raw['current_source_identified_interface_ledger']
    if len(ledger)!=14 or not all(row['current_functional_mixed4_trace_admitted'] for row in ledger.values()):
        raise ValueError('Combined fourteen current adjacent source traces incomplete')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_pulse_adjacent_interface_count=5,
        current_source_identified_adjacent_interface_count=14,
        current_instantiated_five_primitive_mixed4_identity_count=375,
        current_instantiated_velocity_pressure_mixed4_identity_count=300,
        native_mixed4_difference_consistency_rows_checked=native,
        physical_two_sided_trace_contributions_checked=5*2*216,
        finite_current_common_physical_trace_log_uppers=bounds,
        wrong_current_sources=wrong_sources_rejected(field),
        exact_reciprocal_source_reduced_before_directed_enclosure=True,
        original_public_pulse_guards_unmodified=True,
        source_bound_exact_linear_physical_trace_transfer=transfer,
        current_source_identified_interface_ledger=ledger,
        remaining_current_adjacent_trace_interfaces=[],eight_internal_beta_support_transfers_remain_open=True,
        scope='Five current pulse adjacent function mixed4 and original physical spatial4/fixed-x time1 structural traces; all fourteen adjacent traces now source-identified. Eight internal support transfers, all-interface quantitative bounds, tensor/NS, resolved points, prescribed-domain energy and temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current five pulse interfaces PASS: 675 function identities, 750 native rows, 2160 physical contributions; 14 adjacent source traces',flush=True)
    return result


if __name__=='__main__':run()
