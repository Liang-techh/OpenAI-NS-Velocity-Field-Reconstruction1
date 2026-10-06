"""Replay current eight postpulse source identities and bounded trace maps."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import (
    CurrentPostpulseInterfaces,current_postpulse_trace_proof,POSTSEAMS,GATES,OPEN,
    HERE,PREFIX,NAME,RECEIPT,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_complete_physical_assembly import PUBLIC_TO_INTERNAL,METHODS
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def cloned_graph(physical):
    """Definition-only clone preserves all baseline links, never mutates live owners."""
    p=copy.copy(physical);h=p.history=copy.copy(physical.history)
    for key in ('flatten','outer','steep','heat'):
        setattr(h,key,copy.copy(getattr(physical.history,key)));setattr(p,key,getattr(h,key))
    h.outer.flatten=h.flatten;h.steep.outer=h.outer;h.heat.outer=h.outer;h.heat.steep=h.steep
    p.exterior=copy.copy(physical.exterior);p.exterior.history=h;p.exterior.heat=h.heat
    # The dispatch forwards missing attributes, so copy.copy's protocol
    # probes would recurse before base is installed on an empty instance.
    p.dispatch=object.__new__(type(physical.dispatch))
    p.dispatch.__dict__=dict(physical.dispatch.__dict__)
    p.dispatch.assembly=p;p.dispatch.history=h
    p.dispatch.providers=dict(physical.dispatch.providers)
    for chart,internal in PUBLIC_TO_INTERNAL.items():p.dispatch.providers[chart]=getattr(h,METHODS[internal][0])
    p.source_assembly=p.dispatch;p.assert_graph()
    return p


def wrong_sources_rejected(field):
    rejected=[]
    for label in ('different_complete_future','different_absolute_pressure_datum','changed_phase_length','different_heat_Ev2_source'):
        p=cloned_graph(field.physical)
        if label=='different_complete_future':p.heat.future=object()
        elif label=='different_absolute_pressure_datum':
            p.flatten.inlet=copy.copy(p.flatten.inlet);p.flatten.inlet.datum=object()
        elif label=='changed_phase_length':p.steep.Ts=p.steep.Ts+1
        else:p.heat.Ev2=p.heat.Ev2+1
        try:current_postpulse_trace_proof(p)
        except ValueError:rejected.append(label);continue
        raise ValueError('Changed current postpulse source accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),baseline_clone_graph_passed=True,
        live_current_owners_unmodified=True,passed=True)


def trace_bounds(point):
    a=point['left_physical_trace'];b=point['right_physical_trace'];bounds={}
    lo,hi=endpoints(a['source_logR_enclosure']-b['source_logR_enclosure'])
    if not mp.isfinite(lo) or not mp.isfinite(hi) or not lo<=0<=hi:
        raise ValueError('Source-proved common radius contradicts its two enclosures')
    for category in ('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative'):
        def leaves(node,path=()):
            if 'exact_zero' in node:yield path,node
            else:
                for key,value in node.items():yield from leaves(value,path+(key,))
        left=dict(leaves(a[category]));right=dict(leaves(b[category]))
        if set(left)!=set(right):raise ValueError('Original physical contribution layout differs')
        for key,row in left.items():
            other=right[key]
            if endpoints(row['physical_lambda_exponent'])!=endpoints(other['physical_lambda_exponent']):
                raise ValueError('Physical derivative time powers differ')
            if bool(row['exact_zero'])!=bool(other['exact_zero']):raise ValueError('Physical source zero trace differs')
            if row['exact_zero']:upper=None
            else:
                values=[endpoints(v['log_absolute_upper'])[1] for v in (row,other)]
                if not all(mp.isfinite(v) for v in values):raise ValueError('Nonfinite postpulse physical trace bound')
                upper=max(values)
            bounds[category+':'+':'.join(key)]=upper
    if len(bounds)!=216:raise ValueError('Postpulse spatial4/time1 trace coverage incomplete')
    return bounds


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPostpulseInterfaces(require_checked=False)
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_whole_Z_postpulse_trace_views'}:
        raise ValueError('Current postpulse source/manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Postpulse producer claims acceptance')
    proof=current_postpulse_trace_proof(field.physical)
    if encode(pack(proof))!=raw['current_postpulse_source_trace_proof']:
        raise ValueError('Current complete-energy functional theorem differs')
    if len(proof['current_instantiated_seam_mixed4_identities'])!=480 or not all(proof['current_instantiated_seam_mixed4_identities'].values()):
        raise ValueError('Eight source-identified mixed4 seam grids missing')
    transfer=proof['source_bound_exact_linear_physical_trace_transfer']
    if not transfer['passed'] or not transfer['source_equal_mixed_rows_give_zero_exact_physical_operator_difference']:
        raise ValueError('Source-bound exact original physical operator transfer missing')
    names=[row[0] for row in POSTSEAMS]
    if set(raw['current_whole_Z_postpulse_trace_views'])!=set(names):raise ValueError('Eight whole-Z endpoint pairs required')
    native=0;bounds={}
    for name in names:
        point=field.endpoint(name,[-1,1])
        if encode(pack(point))!=raw['current_whole_Z_postpulse_trace_views'][name]:raise ValueError('Postpulse endpoint replay differs: '+name)
        if any(point[k] for k in OPEN):raise ValueError('Postpulse endpoint scope exceeds traces')
        if len(point['same_current_native_mixed4_difference_enclosures'])!=6:raise ValueError('Velocity/pressure/angular/energy trace required')
        for rows in point['same_current_native_mixed4_difference_enclosures'].values():
            if len(rows)!=15:raise ValueError('Ordinary mixed4 consistency grid incomplete')
            for value in rows.values():
                lo,hi=endpoints(value)
                if not mp.isfinite(lo) or not mp.isfinite(hi) or not lo<=0<=hi:
                    raise ValueError('Source identity contradicts native endpoint enclosure: '+name)
                native+=1
        bounds[name]=trace_bounds(point)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_postpulse_interface_count=8,
        current_source_identified_adjacent_interface_count=9,
        current_instantiated_primitive_mixed4_identity_count=480,
        native_mixed4_difference_consistency_rows_checked=native,
        physical_two_sided_trace_contributions_checked=8*2*216,
        finite_current_common_physical_trace_log_uppers=bounds,
        wrong_current_sources=wrong_sources_rejected(field),source_identity_precedes_numeric_consistency=True,
        source_bound_exact_linear_physical_trace_transfer=transfer,
        current_source_identified_interface_ledger=raw['current_source_identified_interface_ledger'],
        remaining_current_adjacent_trace_interfaces=raw['remaining_current_adjacent_trace_interfaces'],
        scope='Eight current postpulse ordinary mixed4 and source-bound original physical spatial4/fixed-x time1 structural traces for abs(Z)<1 and positive tau; limiting bounds cover Z=+/-1, finite bounds at the selected positive time sector. Earlier pulse/support, all-interface global bounds, tensor/NS, points and temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current eight postpulse interfaces PASS: 480 source identities, 720 native rows, 3456 physical contributions',flush=True)
    return result


if __name__=='__main__':run()
